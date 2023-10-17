import awkward as ak
import numba as nb
import numpy as np
import pandas as pd
import uproot
import coffea.nanoevents.methods.vector as vector
import os, psutil

from BTVNanoCommissioning.utils.AK4_parameters import correction_config as config

###############
#  HLT table  #
###############
# HLTs for BTA (qcd) and BTA_ttbar workflows
BTA_HLT = [
    "PFJet40",
    "PFJet60",
    "PFJet80",
    "PFJet140",
    "PFJet200",
    "PFJet260",
    "PFJet320",
    "PFJet400",
    "PFJet450",
    "PFJet500",
    "PFJet550",
    "AK8PFJet40",
    "AK8PFJet60",
    "AK8PFJet80",
    "AK8PFJet140",
    "AK8PFJet200",
    "AK8PFJet260",
    "AK8PFJet320",
    "AK8PFJet400",
    "AK8PFJet450",
    "AK8PFJet500",
    "AK8PFJet550",
    "PFHT180",
    "PFHT250",
    "PFHT370",
    "PFHT430",
    "PFHT510",
    "PFHT590",
    "PFHT680",
    "PFHT780",
    "PFHT890",
    "PFHT1050",
    "BTagMu_AK4DiJet20_Mu5",
    "BTagMu_AK4DiJet40_Mu5",
    "BTagMu_AK4DiJet70_Mu5",
    "BTagMu_AK4DiJet110_Mu5",
    "BTagMu_AK4DiJet170_Mu5",
    "BTagMu_AK4Jet300_Mu5",
    "BTagMu_AK8DiJet170_Mu5",
    "BTagMu_AK8Jet300_Mu5",
]
BTA_ttbar_HLT_chns = [
    ("Ele32_WPTight_Gsf", "e"),
    ("IsoMu24", "m"),
    ("Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ", "em"),
    ("Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL", "em"),
]

################
#  Mass table  #
################
# mass table from https://github.com/scikit-hep/particle/blob/master/src/particle/data/particle2022.csv and https://gitlab.cern.ch/lhcb-conddb/DDDB/-/blob/master/param/ParticleTable.txt
df_main = pd.read_csv(
    "src/BTVNanoCommissioning/helpers/particle2022.csv", delimiter=",", skiprows=1
)
df_back = pd.read_csv(
    "src/BTVNanoCommissioning/helpers/ParticleTable.csv", delimiter=",", skiprows=4
)
df_main, df_back = df_main.astype({"ID": int}), df_back.astype({"PDGID": int})
main = dict(zip(df_main.ID, df_main.Mass / 1000.0))
backup = dict(zip(df_back.PDGID, df_back["MASS(GeV)"]))
hadron_mass_table = {**main, **{k: v for k, v in backup.items() if k not in main}}


###############
#  Functions  #
###############
def is_from_GSP(GenPart):
    QGP = ak.zeros_like(GenPart.genPartIdxMother)
    QGP = (
        (GenPart.genPartIdxMother >= 0)
        & (GenPart.genPartIdxMother2 == -1)
        & (abs(GenPart.parent.pdgId) == 21)
        & (
            (abs(GenPart.parent.pdgId) == abs(GenPart.pdgId))
            | abs(GenPart.parent.pdgId)
            == 21
        )
        & ~(
            (abs(GenPart.pdgId) == abs(GenPart.parent.pdgId))
            & (GenPart.parent.status >= 21)
            & (GenPart.parent.status <= 29)
        )
    )
    rest_QGP = ak.mask(GenPart, ~QGP)
    restGenPart = rest_QGP.parent
    while ak.any(ak.is_none(restGenPart.parent.pdgId, axis=-1) == False):
        mask_forbad = (
            (restGenPart.genPartIdxMother >= 0)
            & (restGenPart.genPartIdxMother2 == -1)
            & (
                (abs(restGenPart.parent.pdgId) == abs(restGenPart.pdgId))
                | abs(restGenPart.parent.pdgId)
                == 21
            )
            & ~(
                (abs(restGenPart.pdgId) == abs(GenPart.parent.pdgId))
                & (restGenPart.parent.status >= 21)
                & (restGenPart.parent.status <= 29)
            )
        )
        mask_forbad = ak.fill_none(mask_forbad, False)
        QGP = (mask_forbad & ak.fill_none((abs(restGenPart.pdgId == 21)), False)) | QGP

        rest_QGP = ak.mask(restGenPart, (~QGP) & (mask_forbad))
        restGenPart = rest_QGP.parent
        if ak.all(restGenPart.parent.pdgId == True):
            break

    return ak.fill_none(QGP, False)


@nb.njit
def to_bitwise_trigger(pass_trig, builder):
    for it in pass_trig:
        # group by every 32 bits
        builder.begin_list()
        for bitidx in range(len(it) // 32 + 1):
            trig = 0
            start = bitidx * 32
            end = min((bitidx + 1) * 32, len(it))
            for i, b in enumerate(it[start:end]):
                trig += b << i
            builder.integer(trig)
        builder.end_list()

    return builder


@nb.vectorize([nb.float64(nb.int64)], forceobj=True)
def get_hadron_mass(hadron_pids):
    return hadron_mass_table[abs(hadron_pids)]


def cumsum(array):
    layout = array.layout
    layout.content
    scan = ak.Array(
        ak.layout.ListOffsetArray64(
            layout.offsets, ak.layout.NumpyArray(np.cumsum(layout.content))
        )
    )
    cumsum_array = ak.fill_none(scan - ak.firsts(scan) + ak.firsts(array), [], axis=0)
    return cumsum_array


def calc_ip_vector(obj, dxy, dz, is_3d=False):
    """Calculate the 2D or 3D impact parameter vector, given the track obj (with 4-mom),
    and its dxy and dz, taking the standard definition from NanoAOD"""

    # 2D impact parameter
    pvec = ak.zip(
        {
            "x": obj.px,
            "y": obj.py,
            "z": obj.pz,
        },
        behavior=vector.behavior,
        with_name="ThreeVector",
    )
    zvec = ak.zip(
        {
            "x": ak.zeros_like(dxy),
            "y": ak.zeros_like(dxy),
            "z": ak.zeros_like(dxy) + 1,
        },
        behavior=vector.behavior,
        with_name="ThreeVector",
    )
    # 2D impact parameter vector: (-py, px) / pt * dxy
    ipvec_2d = zvec.cross(pvec) * dxy / obj.pt

    if is_3d == False:
        return ipvec_2d

    # Then calculate the 3D impact parameter vector
    # first, extend ipvec_2d to 3D space
    ipvec_2d_ext = ak.zip(
        {
            "x": ipvec_2d.x,
            "y": ipvec_2d.y,
            "z": dz,
        },
        behavior=vector.behavior,
        with_name="ThreeVector",
    )
    # then, get the closest distance to the track on 3D geometry
    ipvec_3d = ipvec_2d_ext - ipvec_2d_ext.dot(pvec) / pvec.p2 * pvec
    return ipvec_3d


###############
#   Classes   #
###############


class JPCalibHandler(object):
    def __init__(self, campaign, isRealData, dataset):
        r"""
        A tool for calculating the track probability and jet probability
            campaign: campaign name
            isRealData: whether the dataset is real data
            dataset: dataset name from events.metadata["dataset"]
        """
        if isRealData:
            for key in config[campaign]["JPCalib"]:
                if key in dataset:
                    filename = config[campaign]["JPCalib"][key]
                    break
            else:
                raise ValueError(f"No JPCalib file found for dataset {dataset}")
        else:
            filename = config[campaign]["JPCalib"]["MC"]

        # print(f'Using JPCalib file {filename}')

        templates = uproot.open(
            f"src/BTVNanoCommissioning/data/JPCalib/{campaign}/{filename}"
        )
        self.ipsig_histo_val = np.array(
            [templates[f"histoCat{i}"].values() for i in range(10)]
        )
        self.ipsig_histo_tot = np.sum(self.ipsig_histo_val, axis=1)
        self.values_cumsum = np.cumsum(self.ipsig_histo_val[:, ::-1], axis=1)[:, ::-1]
        self.edges = templates["histoCat0"].axes[0].edges()

    def flatten(self, array):
        r"""
        Get the fully flattened array and its layout for each layer
        """
        layouts = []
        array_fl = array
        while str(ak.type(array_fl)).count("*") > 1:
            layouts.append(ak.num(array_fl))
            array_fl = ak.flatten(array_fl)
        return array_fl, layouts

    def unflatten(self, array_fl, layouts):
        r"""
        Recover a flattened array using the original layouts
        """
        array = array_fl
        for layout in layouts[::-1]:
            array = ak.unflatten(array, layout)
        return array

    def calc_track_proba(self, ipsig: ak.Array, cat: ak.Array):
        r"""
        Calculate the track probability from the integral of the track IPsig templates, given the IPsig and category.
        Reference code: https://github.com/cms-sw/cmssw/blob/CMSSW_13_0_X/RecoBTag/TrackProbability/src/HistogramProbabilityEstimator.cc
            ipsig: IP significance array
            cat: category array (0-9)
        """

        if ak.any(cat < 0) or ak.any(cat > 9):
            raise ValueError("Category out of range [0, 9]")

        # get the fully flattened array of the input while storing its layouts for later recovery
        ipsig_fl, layouts = self.flatten(ipsig)
        cat_fl = ak.flatten(cat, axis=None)

        # index of the IPsig bins
        ipsig_fl = abs(ipsig_fl)
        ipsig_fl_index = np.minimum(
            np.searchsorted(self.edges, ipsig_fl), self.ipsig_histo_val.shape[1] - 1
        )

        # retrieve the cumsum value (\int_{ipsig}^{inf} p(ipsig') d(ipsig')) from the correct template
        ipsig_cumsum_fl = self.values_cumsum[cat_fl, ipsig_fl_index]

        # calculate the track probability as (\int_{ipsig}^{inf} ..) / (\int_{0}^{inf} ..) * sign(IPsig)
        proba_fl = (ipsig_cumsum_fl / self.ipsig_histo_tot[cat_fl]) * np.sign(ipsig_fl)

        # recover the original layout
        proba = self.unflatten(proba_fl, layouts)
        return proba

    def calc_jet_proba(self, proba):
        # Calculate jet probability (JP)
        # according to jetProbability func in https://github.com/cms-sw/cmssw/blob/CMSSW_13_0_X/RecoBTag/ImpactParameter/interface/TemplatedJetProbabilityComputer.h

        # minium proba = 0.5%
        proba = np.maximum(proba, 0.005)  # dim: (evt, jet, trk)

        ntrk = ak.num(proba, axis=-1)  # dim: (evt, jet), the number of tracks in a jet
        prodproba_log = ak.sum(
            np.log(proba), axis=-1
        )  # dim: (evt, jet), the log(Π(proba)) of all tracks in a jet
        prodproba_log_m_log = ak.where(
            (ak.num(proba, axis=-1) >= 2) & (prodproba_log < 0),
            np.log(-prodproba_log),
            0,
        )  # log(-logΠ), if >=2 tracks in a jet

        # now calculating Σ_tr{0..N-1} ((-logΠ)^tr / tr!)
        trk_index = ak.local_index(proba)
        fact_array = ak.concatenate(
            [
                [1.0],
                np.arange(1, max(5, ak.max(trk_index) + 1), dtype=np.float64).cumprod(),
            ]
        )  # construct a factorial array
        trk_index_fl, _layouts = self.flatten(trk_index)
        lfact = self.unflatten(
            fact_array[trk_index_fl], _layouts
        )  # dim: (evt, jet, trk), nested factorial array given the track index

        prob = ak.sum(
            np.exp(trk_index * prodproba_log_m_log - np.log(lfact)), axis=-1
        )  # dim: (evt, jet), Σ_tr{0..N-1} ((-logΠ)^tr / tr!)

        prob_jet = np.minimum(
            np.exp(np.maximum(np.log(np.maximum(prob, 1e-30)) + prodproba_log, -30.0)),
            1.0,
        )  # dim: (evt, jet), calculating Π * Σ_tr{0..N-1} ((-logΠ)^tr / tr!)

        prob_jet = ak.where(prodproba_log < 0, prob_jet, 1.0)
        prob_jet = np.maximum(prob_jet, 1e-30)

        return prob_jet
