import importlib.resources
import cloudpickle, gzip, contextlib
import copy, os, re

import numpy as np
import awkward as ak
import uproot
from coffea.lookup_tools import extractor, txt_converters, rochester_lookup
from coffea.lumi_tools import LumiMask

from coffea.jetmet_tools.CorrectedMETFactory import corrected_polar_met

from coffea.analysis_tools import Weights

from coffea.btag_tools import BTagScaleFactor
import correctionlib
from BTVNanoCommissioning.helpers.func import update, _compile_jec_, _load_jmefactory
from BTVNanoCommissioning.helpers.cTagSFReader import getSF

from BTVNanoCommissioning.utils.AK4_parameters import correction_config as config


def load_SF(year, campaign, syst=False):
    """
    Load scale factors (SF) for a given year and campaign.

    This function reads scale factors from the specified campaign configuration and returns them in a suitable format.
    It handles different types of scale factors, such as pileup weights, and checks for the existence of files in
    the jsonpog-integration directory or custom files.

    Parameters:
    year (str): The year for which to load the scale factors.
    campaign (str): The name of the campaign for which to load the scale factors.
    syst (bool, optional): A flag to indicate whether to load systematic variations. Default is False.

    Returns:
    dict: A dictionary containing the scale factors, where keys are the relevant identifiers and values are the scale factors.

    Raises:
    FileNotFoundError: If the specified file does not exist.
    ValueError: If the file content is not in the expected format.
    KeyError: If the specified campaign or year is not found in the configuration.
    """
    # read the configuration file to get the correct SFs
    correct_map = {"campaign": campaign}

    for SF in config[campaign].keys():
        if SF == "lumiMask":
            continue
        ## pileup weight
        if SF == "PU":
            ## Check whether files in jsonpog-integration exist
            if os.path.exists(
                f"/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/LUM/{year}_{campaign}"
            ):
                correct_map["PU"] = correctionlib.CorrectionSet.from_file(
                    f"/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/LUM/{year}_{campaign}/puWeights.json.gz"
                )
            ## Otherwise custom files
            else:
                _pu_path = f"BTVNanoCommissioning.data.PU.{campaign}"
                with importlib.resources.path(
                    _pu_path, config[campaign]["PU"]
                ) as filename:
                    if str(filename).endswith(".pkl.gz"):
                        with gzip.open(filename) as fin:
                            correct_map["PU"] = cloudpickle.load(fin)[
                                "2017_pileupweight"
                            ]
                    elif str(filename).endswith(".json.gz"):
                        correct_map["PU"] = correctionlib.CorrectionSet.from_file(
                            str(filename)
                        )
                    elif str(filename).endswith(".histo.root"):
                        ext = extractor()
                        ext.add_weight_sets([f"* * {filename}"])
                        ext.finalize()
                        correct_map["PU"] = ext.make_evaluator()

        ## btag weight
        elif SF == "BTV":
            if "btag" in config[campaign]["BTV"].keys() and config[campaign]["BTV"][
                "btag"
            ].endswith(".json.gz"):
                correct_map["btag"] = correctionlib.CorrectionSet.from_file(
                    importlib.resources.path(
                        f"BTVNanoCommissioning.data.BTV.{year}_{campaign}", filename
                    )
                )
            if "ctag" in config[campaign]["BTV"].keys() and config[campaign]["BTV"][
                "ctag"
            ].endswith(".json.gz"):
                correct_map["btag"] = correctionlib.CorrectionSet.from_file(
                    importlib.resources.path(
                        f"BTVNanoCommissioning.data.BTV.{year}_{campaign}", filename
                    )
                )
            if os.path.exists(
                f"/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/{year}_{campaign}"
            ):
                correct_map["btag"] = correctionlib.CorrectionSet.from_file(
                    f"/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/{year}_{campaign}/btagging.json.gz"
                )
                correct_map["ctag"] = correctionlib.CorrectionSet.from_file(
                    f"/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/BTV/{year}_{campaign}/ctagging.json.gz"
                )
            else:
                correct_map["btag"] = {}
                correct_map["ctag"] = {}
                correct_map["btv_cfg"] = config[campaign]["BTV"]
                _btag_path = f"BTVNanoCommissioning.data.BTV.{year}_{campaign}"
                for tagger in config[campaign]["BTV"]:
                    with importlib.resources.path(
                        _btag_path, config[campaign]["BTV"][tagger]
                    ) as filename:
                        if "B" in tagger:
                            if filename.endswith(".json.gz"):
                                correct_map["btag"] = (
                                    correctionlib.CorrectionSet.from_file(filename)
                                )
                            else:
                                correct_map["btag"][tagger] = BTagScaleFactor(
                                    filename,
                                    BTagScaleFactor.RESHAPE,
                                    methods="iterativefit,iterativefit,iterativefit",
                                )
                        else:
                            if filename.endswith(".json.gz"):
                                correct_map["ctag"] = (
                                    correctionlib.CorrectionSet.from_file(filename)
                                )
                            else:
                                correct_map["ctag"][tagger] = BTagScaleFactor(
                                    filename,
                                    BTagScaleFactor.RESHAPE,
                                    methods="iterativefit,iterativefit,iterativefit",
                                )

        ## lepton SFs
        elif SF == "LSF":
            correct_map["MUO_cfg"] = {
                mu: f
                for mu, f in config[campaign]["LSF"].items()
                if "mu" in mu and "_json" not in mu
            }
            correct_map["EGM_cfg"] = {
                e: f
                for e, f in config[campaign]["LSF"].items()
                if "ele" in e and "_json" not in e
            }
            ## Muon
            if os.path.exists(
                f"/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/{year}_{campaign}"
            ):
                correct_map["MUO"] = correctionlib.CorrectionSet.from_file(
                    f"/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/MUO/{year}_{campaign}/muon_Z.json.gz"
                )
            if os.path.exists(
                f"/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/{year}_{campaign}"
            ):
                correct_map["EGM"] = correctionlib.CorrectionSet.from_file(
                    f"/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/EGM/{year}_{campaign}/electron.json.gz"
                )
            if any(
                np.char.find(np.array(list(config[campaign]["LSF"].keys())), "mu_json")
                != -1
            ):
                correct_map["MUO"] = correctionlib.CorrectionSet.from_file(
                    f"src/BTVNanoCommissioning/data/LSF/{year}_{campaign}/{config[campaign]['LSF']['mu_json']}"
                )
            if any(
                np.char.find(np.array(list(config[campaign]["LSF"].keys())), "ele_json")
                != -1
            ):
                correct_map["EGM"] = correctionlib.CorrectionSet.from_file(
                    f"src/BTVNanoCommissioning/data/LSF/{year}_{campaign}/{config[campaign]['LSF']['ele_json']}"
                )

            ### Check if any custom corrections needed
            # FIXME: (some low pT muons not supported in jsonpog-integration at the moment)
            if (
                "histo.json" in "\t".join(list(config[campaign]["LSF"].values()))
                or "histo.txt" in "\t".join(list(config[campaign]["LSF"].values()))
                or "histo.root" in "\t".join(list(config[campaign]["LSF"].values()))
            ):
                _mu_path = f"BTVNanoCommissioning.data.LSF.{campaign}"
                ext = extractor()
                with contextlib.ExitStack() as stack:
                    inputs, real_paths = [
                        k
                        for k in correct_map["MUO_cfg"].keys()
                        if "histo.json" in correct_map["MUO_cfg"][k]
                        or "histo.txt" in correct_map["MUO_cfg"][k]
                        or "histo.root" in correct_map["MUO_cfg"][k]
                    ], [
                        stack.enter_context(importlib.resources.path(_mu_path, f))
                        for f in correct_map["MUO_cfg"].values()
                        if ".json" in f or ".txt" in f or ".root" in f
                    ]

                    inputs = [
                        i.split(" ")[0] + " *" if "_low" in i else i for i in inputs
                    ]

                    ext.add_weight_sets(
                        [
                            f"{paths} {file}"
                            for paths, file in zip(inputs, real_paths)
                            if "histo.json" in str(file)
                            or "histo.txt" in str(file)
                            or "histo.root" in str(file)
                        ]
                    )
                    if syst:
                        ext.add_weight_sets(
                            paths.split(" ")[0]
                            + "_error "
                            + paths.split(" ")[1]
                            + "_error "
                            + file
                            for paths, file in zip(inputs, real_paths)
                            if ".root" in str(file)
                        )
                ext.finalize()
                correct_map["MUO_custom"] = ext.make_evaluator()

                _ele_path = f"BTVNanoCommissioning.data.LSF.{campaign}"
                ext = extractor()
                with contextlib.ExitStack() as stack:
                    inputs, real_paths = [
                        k
                        for k in correct_map["EGM_cfg"].keys()
                        if "histo.json" in correct_map["EGM_cfg"][k]
                        or "histo.txt" in correct_map["EGM_cfg"][k]
                        or "histo.root" in correct_map["EGM_cfg"][k]
                    ], [
                        stack.enter_context(importlib.resources.path(_ele_path, f))
                        for f in correct_map["EGM_cfg"].values()
                        if "histo.json" in f or ".txt" in f or ".root" in f
                    ]
                    ext.add_weight_sets(
                        [
                            f"{paths} {file}"
                            for paths, file in zip(inputs, real_paths)
                            if "histo.json" in str(file)
                            or "histo.txt" in str(file)
                            or "histo.root" in str(file)
                        ]
                    )
                    if syst:
                        ext.add_weight_sets(
                            paths.split(" ")[0]
                            + "_error "
                            + paths.split(" ")[1]
                            + "_error "
                            + file
                            for paths, file in zip(inputs, real_paths)
                            if ".root" in str(file)
                        )
                ext.finalize()
                correct_map["EGM_custom"] = ext.make_evaluator()

        ## rochester muon momentum correction
        elif SF == "roccor":
            if "2016postVFP_UL" == campaign:
                filename = "RoccoR2016bUL.txt"
            elif "2016preVFP_UL" in campaign:
                filename = "RoccoR2016aUL.txt"
            elif "2017_UL" in campaign:
                filename = "RoccoR2017UL.txt"
            if "2018_UL" in campaign:
                filename = "RoccoR2018UL.txt"

            full_path = "src/BTVNanoCommissioning/data/LSF/roccor/" + filename
            rochester_data = txt_converters.convert_rochester_file(
                full_path, loaduncs=True
            )
            correct_map["roccor"] = rochester_lookup.rochester_lookup(rochester_data)

        ## JME corrections
        elif SF == "JME":
            if "name" in config[campaign]["JME"].keys():

                if not os.path.exists(
                    f"src/BTVNanoCommissioning/data/JME/{year}_{campaign}/jec_compiled_{config[campaign]['JME']['name']}.pkl.gz"
                ):
                    _compile_jec_(
                        year,
                        campaign,
                        config[campaign]["JME"],
                        f"jec_compiled_{config[campaign]['JME']['name']}",
                    )

                correct_map["JME"] = _load_jmefactory(
                    year,
                    campaign,
                    f"jec_compiled_{config[campaign]['JME']['name']}.pkl.gz",
                )
            elif os.path.exists(
                f"/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/{year}_{campaign}/jet_jerc.json.gz"
            ):
                correct_map["JME"] = correctionlib.CorrectionSet.from_file(
                    f"/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/{year}_{campaign}/jet_jerc.json.gz"
                )
                correct_map["JME_cfg"] = config[campaign]["JME"]
                for dataset in correct_map["JME_cfg"].keys():
                    if (
                        np.all(
                            np.char.find(
                                np.array(list(correct_map["JME"].keys())),
                                correct_map["JME_cfg"][dataset],
                            )
                        )
                        == -1
                    ):
                        raise (
                            f"{dataset} has no JEC map : {correct_map['JME_cfg'][dataset]} available"
                        )
        elif SF == "JMAR":
            if os.path.exists(
                f"/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/{year}_{campaign}/jmar.json.gz"
            ):
                correct_map["JMAR_cfg"] = {
                    j: f for j, f in config[campaign]["JMAR"].items()
                }
                correct_map["JMAR"] = correctionlib.CorrectionSet.from_file(
                    f"/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/JME/{year}_{campaign}/jmar.json.gz"
                )
        elif SF == "jetveto":
            ext = extractor()
            with contextlib.ExitStack() as stack:
                ext.add_weight_sets(
                    [
                        f"{run} {stack.enter_context(importlib.resources.path(f'BTVNanoCommissioning.data.JME.{year}_{campaign}',file))}"
                        for run, file in config[campaign]["jetveto"].items()
                    ]
                )

            ext.finalize()
            correct_map["jetveto_cfg"] = {
                j: f for j, f in config[campaign]["jetveto"].items()
            }
            correct_map["jetveto"] = ext.make_evaluator()

    return correct_map


def load_lumi(campaign):
    """
    Load luminosity mask for a given campaign.

    This function reads the luminosity mask file for the specified campaign and returns a `LumiMask` object.

    Parameters:
    campaign (str): The name of the campaign for which to load the luminosity mask.

    Returns:
    LumiMask: An object representing the luminosity mask for the specified campaign.

    Raises:
    KeyError: If the specified campaign is not found in the configuration.
    FileNotFoundError: If the luminosity mask file does not exist.
    """

    _lumi_path = "BTVNanoCommissioning.data.lumiMasks"
    with importlib.resources.path(_lumi_path, config[campaign]["lumiMask"]) as filename:
        return LumiMask(filename)


# wrapped up common shifts
"""
BTVNanoCommissioning.utils.correction

This module provides functions to handle corrections and uncertainties for scaling factors and weight variables in the BTVNanoCommissioning framework.

Features:
- Scaling factors and weight variables:
    - Corrections are wrapped with `weights` from the `coffea.analysis_tools.Weights` class along the event axis.
    - Standardized correction files are handled using `correctionlib`.
    - Additional methods for applying corrections can be found in the `coffea` documentation.
    - Uncertainties are added using up/down variations, except for btag SFs which use `add_multivariation`.

Example for Scaling Factors (SFs):
```python
## Initialization, add EGM map from correctionlib
correction_map["EGM"] = correctionlib.CorrectionSet.from_file(
                f"src/BTVNanoCommissioning/jsonpog-integration/POG/EGM/{campaign}/electron.json.gz"
            )
## Initialization, add EGM map from custom file by extractor 
 ext = extractor()
 ext.add_weight_sets(["eleID EGamma2D {filename}.root"])
 ext.finalize()
 correction_map["EGM"] = ext.make_evaluator()

 # evaluation depends on file types...ignore here!
 ## add SFs & uncertainties to weight function
weights.add(sf.split(" ")[0], sfs_alle, sfs_alle_up, sfs_alle_down) 
```

"""


def common_shifts(self, events):
    """
    Apply common shifts to a events(mostly affect energy resolution/scale of objects).

    This function applies common shifts to the input DataFrame based on the specified shift type.
    It modifies the DataFrame in place to reflect the systematic variations/dedicated corrections.
    This includes JERC corrections, rochester corrections.


    - Scale/Resolution Corrections:
    Construct a shift list with tuples of (obj_dict, shift_name).
    These corrections are applied independently on all objects by updating the contents of the branch.
    Normally done before selection to apply updated objects.
    Uncertainties are handled by updating object collections with up/down variations.
    Example for Shift List:
    ```python
    # nominal correction
    shift = [({"Jet": jets, "MET": met, "Muon" : muon}, None)]
    # add variations
    shifts += [
                    (
                        {
                            "Jet": jets.JES_Total.up,
                            "MET": met.JES_Total.up,
                        },
                        "JESUp",
                    )]
    shifts += [
                    (
                        {
                            "Jet": jets.JES_Total.down,
                            "MET": met.JES_Total.down,
                        },
                        "JESDown",
                    )]
    ``
    Different treatment for weights and scale/resolution shifts is necessary to ensure accurate corrections and uncertainties are applied to the data.

    Parameters:
    self (dict): The configuration dictionary from SF_map containing the scale factors and other settings.
    events (events): The input events containing the data to be shifted.

    Returns:
    pandas.DataFrame: The DataFrame with the applied systematic shifts.
    """

    isRealData = not hasattr(events, "genWeight")
    dataset = events.metadata["dataset"]
    shifts = []
    if "JME" in self.SF_map.keys():
        syst_JERC = self.isSyst
        if self.isSyst == "JERC_split":
            syst_JERC = "split"
        shifts = JME_shifts(
            shifts, self.SF_map, events, self._campaign, isRealData, syst_JERC
        )
    else:
        if int(self._year) < 2020:
            shifts = [
                ({"Jet": events.Jet, "MET": events.MET, "Muon": events.Muon}, None)
            ]
        else:
            shifts = [
                (
                    {
                        "Jet": events.Jet,
                        "MET": events.PuppiMET,
                        "Muon": events.Muon,
                    },
                    None,
                )
            ]
    if "roccor" in self.SF_map.keys():
        shifts = Roccor_shifts(shifts, self.SF_map, events, isRealData, False)
    else:
        shifts[0][0]["Muon"] = events.Muon
    return shifts


##JEC
# FIXME: would be nicer if we can move to correctionlib in the future together with factory and workable


def add_jec_variables(jets, event_rho):
    jets["pt_raw"] = (1 - jets.rawFactor) * jets.pt
    jets["mass_raw"] = (1 - jets.rawFactor) * jets.mass
    if hasattr(jets, "genJetIdxG"):
        jets["pt_gen"] = ak.values_astype(
            ak.fill_none(jets.matched_gen.pt, 0), np.float32
        )
    else:
        jets["pt_gen"] = ak.zeros_like(jets.pt)
    jets["event_rho"] = ak.broadcast_arrays(event_rho, jets.pt)[0]
    return jets


## Jet Veto
def jetveto(jets, correct_map):
    """
    Apply a veto to jets based on predefined transverse momentum (pt) and pseudorapidity (eta) thresholds.

    This function filters out jets that do not meet the predefined pt and eta criteria. It also utilizes a correction map
    to apply additional corrections or selections to the jets.

    Parameters:
    jets (iterable): A collection of jet objects or dictionaries containing jet properties.
    correct_map (dict): A dictionary containing correction factors or additional selection criteria for the jets.

    Returns:
    jets: A jets of jets that pass the predefined pt and eta criteria and any additional criteria from the correction map.

    Raises:
    TypeError: If the jets parameter is not an iterable.
    KeyError: If the jet objects do not contain the required 'pt' or 'eta' properties.
    """
    return ak.where(
        correct_map["jetveto"][list(correct_map["jetveto"].keys())[0]](
            jets.phi, jets.eta
        )
        > 0,
        ak.ones_like(jets.eta),
        ak.zeros_like(jets.eta),
    )


## JERC
def JME_shifts(
    shifts,
    correct_map,
    events,
    campaign,
    isRealData,
    systematic=False,
    exclude_jetveto=False,
):
    """
    Apply Jet Energy Corrections (JEC) and Jet Energy Resolutions (JER) shifts to events.

    This function applies JEC and JER shifts to the jets in the events based on the provided correction map and campaign.
    It handles both real data and simulated data, and can optionally apply systematic variations and exclude jet vetoes.

    Parameters:
    shifts (list): A list of shift types to apply (e.g., 'up', 'down').
    correct_map (dict): A dictionary containing correction factors and settings for JEC and JER.
    events (awkward.Array): An array of events containing jet information.
    campaign (str): The name of the campaign for which to apply the corrections.
    isRealData (bool): A flag indicating whether the data is real or simulated.
    systematic (bool, optional): A flag to indicate whether to apply systematic variations. Default is False.
    exclude_jetveto (bool, optional): A flag to indicate whether to exclude jet vetoes. Default is False.

    Returns:
    awkward.Array: The events array with applied JEC and JER shifts.

    Raises:
    KeyError: If required keys are missing in the correct_map.
    ValueError: If the campaign is not recognized or supported.
    """
    dataset = events.metadata["dataset"]
    jecname = ""
    # https://cms-jerc.web.cern.ch/JECUncertaintySources/, currently no recommendation of reduced/ full split sources
    syst_list = [
        i.split("_")[3]
        for i in correct_map["JME"].keys()
        if "MC" in i and "L1" not in i and "L2" not in i and "L3" not in i
    ]
    if "JME" in correct_map.keys():
        ## correctionlib
        if "JME_cfg" in correct_map.keys():
            if isRealData:
                jecname = [
                    v
                    for k, v in correct_map["JME_cfg"].items()
                    if k in events.metadata["dataset"]
                ]
                if len(jecname) > 1:
                    raise ("Multiple uncertainties match to this era")
                else:
                    jecname = jecname[0] + "_DATA"
            else:
                jecname = correct_map["JME_cfg"]["MC"] + "_MC"
            corr = correct_map["JME"].compound[f"{jecname}_L1L2L3Res_AK4PFPuppi"]
            nocorrjet = events.Jet
            nocorrjet["pt_raw"] = (1 - nocorrjet["rawFactor"]) * nocorrjet["pt"]
            nocorrjet["mass_raw"] = (1 - nocorrjet["rawFactor"]) * nocorrjet["mass"]
            nocorrjet["rho"] = ak.broadcast_arrays(
                events.fixedGridRhoFastjetAll, nocorrjet.pt
            )[0]
            j, nj = ak.flatten(nocorrjet), ak.num(nocorrjet)
            values = [
                np.array(
                    j[
                        inputs.name.replace("Jet", "")
                        .replace("Pt", "pt")
                        .replace("Phi", "phi")
                        .replace("Eta", "eta")
                        .replace("Mass", "mass")
                        .replace("Rho", "rho")
                        .replace("A", "area")
                    ]
                )
                for inputs in corr.inputs
            ]
            flatCorrFactor = corr.evaluate(*values)
            corrFactor = ak.unflatten(flatCorrFactor, nj)
            jets = copy.copy(nocorrjet)
            jets["pt_orig"] = ak.values_astype(nocorrjet["pt"], np.float32)
            jets["pt"] = ak.values_astype(nocorrjet["pt_raw"] * corrFactor, np.float32)
            jets["mass"] = ak.values_astype(
                nocorrjet["mass_raw"] * corrFactor, np.float32
            )

            # MET correction, from MET correct factory
            # https://github.com/CoffeaTeam/coffea/blob/d7d02634a8d268b130a4d71f76d8eba6e6e27b96/coffea/jetmet_tools/CorrectedMETFactory.py#L105

            nocorrmet = (
                events.PuppiMET if "22" in campaign or "23" in campaign else events.MET
            )
            form = ak.forms.RecordForm(
                {
                    "pt": nocorrmet.pt.layout.form,
                    "phi": nocorrmet.phi.layout.form,
                }
            )

            met = copy.copy(nocorrmet)
            metinfo = [nocorrmet.pt, nocorrmet.phi, jets.pt, jets.phi, jets.pt_raw]
            met["pt"], met["phi"] = (
                ak.values_astype(corrected_polar_met(*metinfo).pt, np.float32),
                ak.values_astype(corrected_polar_met(*metinfo).phi, np.float32),
            )
            met["orig_pt"], met["orig_phi"] = nocorrmet["pt"], nocorrmet["phi"]
            if systematic != False:

                if systematic != "JERC_split":
                    jesuncmap = correct_map["JME"][f"{jecname}_Total_AK4PFPuppi"]

                    jesunc = ak.unflatten(jesuncmap.evaluate(j.eta, j.pt), nj)

                    jets["JES_Total"] = ak.zip(
                        {
                            "up": ak.copy(jets),
                            "down": ak.copy(jets),
                        }
                    )
                    met["JES_Total"] = ak.zip(
                        {
                            "up": ak.copy(met),
                            "down": ak.copy(met),
                        }
                    )
                    for var in ["up", "down"]:
                        fac = 1.0 if var == "up" else -1.0

                        jets["JES_Total"][var]["pt"] = jets["pt"] * (
                            corrFactor + fac * jesunc
                        )
                        jets["JES_Total"][var]["mass"] = jets["mass"] * (
                            corrFactor + fac * jesunc
                        )

                        met["JES_Total"][var]["pt"] = corrected_polar_met(
                            nocorrmet.pt,
                            nocorrmet.phi,
                            jets.JES_Total[var]["pt"],
                            jets.phi,
                            jets.pt_raw,
                        ).pt
                        met["JES_Total"][var]["phi"] = corrected_polar_met(
                            nocorrmet.pt,
                            nocorrmet.phi,
                            jets.JES_Total[var]["pt"],
                            jets.phi,
                            jets.pt_raw,
                        ).phi

        else:
            if isRealData:
                if "2016preVFP_UL" == campaign:
                    if "2016B" in dataset or "2016C" in dataset or "2016D" in dataset:
                        jecname = "BCD"
                    elif "2016E" in dataset or "2016F" in dataset:
                        jecname = "EF"
                elif "2016postVFP_UL" == campaign:
                    jecname = "FGH"
                elif campaign == "Rereco17_94X":
                    jecname = ""
                elif campaign == "Summer23":
                    if "v4" in dataset:
                        jecname = "Cv4"
                    else:
                        jecname = "Cv123"
                elif re.search(r"[Rr]un20\d{2}([A-Z])", dataset):
                    jecname = re.search(r"[Rr]un20\d{2}([A-Z])", dataset).group(1)
                else:
                    print("No valid jec name")
                    raise NameError
                jecname = "data" + jecname
            else:
                jecname = "MC"

            jets = correct_map["JME"]["jet_factory"][jecname].build(
                add_jec_variables(events.Jet, events.fixedGridRhoFastjetAll),
                lazy_cache=events.caches[0],
            )
            met = correct_map["JME"]["met_factory"].build(events.PuppiMET, jets, {})
            ## systematics
        if not isRealData:
            if systematic != False:
                if systematic == "split":
                    for jes in met.fields:
                        if "JES" not in jes or "Total" in jes:
                            continue

                        shifts += [
                            (
                                {
                                    "Jet": jets[jes]["up"],
                                    "MET": met[jes]["up"],
                                },
                                f"{jes}Up",
                            ),
                            (
                                {
                                    "Jet": jets[jes]["down"],
                                    "MET": met[jes]["down"],
                                },
                                f"{jes}Down",
                            ),
                        ]

                else:
                    if "JES_Total" in jets.fields:
                        shifts += [
                            (
                                {
                                    "Jet": jets.JES_Total.up,
                                    "MET": met.JES_Total.up,
                                },
                                "JESUp",
                            ),
                            (
                                {
                                    "Jet": jets.JES_Total.down,
                                    "MET": met.JES_Total.down,
                                },
                                "JESDown",
                            ),
                        ]
                    if "MET_UnclusteredEnergy" in met.fields:
                        shifts += [
                            (
                                {
                                    "Jet": jets,
                                    "MET": met.MET_UnclusteredEnergy.up,
                                },
                                "UESUp",
                            ),
                            (
                                {
                                    "Jet": jets,
                                    "MET": met.MET_UnclusteredEnergy.down,
                                },
                                "UESDown",
                            ),
                        ]
                    if "JER" in jets.fields:
                        shifts += [
                            (
                                {
                                    "Jet": jets.JER.up,
                                    "MET": met.JER.up,
                                },
                                "JERUp",
                            ),
                            (
                                {
                                    "Jet": jets.JER.down,
                                    "MET": met.JER.down,
                                },
                                "JERDown",
                            ),
                        ]

    else:
        met = events.PuppiMET
        jets = events.Jet
    # perform jet veto
    if "jetveto" in correct_map.keys():
        jets = update(jets, {"veto": jetveto(jets, correct_map)})

        if "Summer22" in campaign:
            jets = jets[jets.veto != 1]
    shifts += [({"Jet": jets, "MET": met}, None)]
    return shifts


## Muon Rochester correction
def Roccor_shifts(shifts, correct_map, events, isRealData, systematic=False):
    """
    Apply Rochester corrections (Roccor) shifts to muons in events.

    This function applies Rochester corrections to the muons in the events based on the provided correction map and campaign.
    It handles both real data and simulated data, and can optionally apply systematic variations.

    Parameters:
    shifts (list): A list of shift types to apply (e.g., 'up', 'down').
    correct_map (dict): A dictionary containing correction factors and settings for Rochester corrections.
    events (awkward.Array): An array of events containing muon information.
    campaign (str): The name of the campaign for which to apply the corrections.
    isRealData (bool): A flag indicating whether the data is real or simulated.
    systematic (bool, optional): A flag to indicate whether to apply systematic variations. Default is False.

    Returns:
    awkward.Array: The events array with applied Rochester corrections.

    Raises:
    KeyError: If required keys are missing in the correct_map.
    ValueError: If the campaign is not recognized or supported.
    """
    mu = events.Muon
    if isRealData:
        SF = correct_map["roccor"].kScaleDT(
            events.Muon.charge, events.Muon.pt, events.Muon.eta, events.Muon.phi
        )

    else:
        hasgen = ~np.isnan(ak.fill_none(events.Muon.matched_gen.pt, np.nan))
        mc_kspread = correct_map["roccor"].kSpreadMC(
            events.Muon.charge[hasgen],
            events.Muon.pt[hasgen],
            events.Muon.eta[hasgen],
            events.Muon.phi[hasgen],
            events.Muon.matched_gen.pt[hasgen],
        )
        mc_rand = np.random.rand(len(ak.flatten(events.Muon.pt, axis=1)))
        mc_rand = ak.unflatten(mc_rand, ak.num(events.Muon.pt))
        mc_ksmear = correct_map["roccor"].kSmearMC(
            events.Muon.charge[~hasgen],
            events.Muon.pt[~hasgen],
            events.Muon.eta[~hasgen],
            events.Muon.phi[~hasgen],
            events.Muon.nTrackerLayers[~hasgen],
            mc_rand[~hasgen],
        )
        SF = np.array(ak.flatten(ak.ones_like(events.Muon.pt)))
        hasgen_flat = np.array(ak.flatten(hasgen))
        SF[hasgen_flat] = np.array(ak.flatten(mc_kspread))
        SF[~hasgen_flat] = np.array(ak.flatten(mc_ksmear))
        SF = ak.unflatten(SF, ak.num(events.Muon.pt))

    mu["pt"] = SF * events.Muon.pt
    # add rochester correction to shift
    for i in range(len(shifts)):
        shifts[i][0]["Muon"] = mu

    if systematic:
        if isRealData:
            err = correct_map["roccor"].kScaleDTerror(
                events.Muon.charge, events.Muon.pt, events.Muon.eta, events.Muon.phi
            )
        else:
            mc_errspread = correct_map["roccor"].kSpreadMCerror(
                events.Muon.charge[hasgen],
                events.Muon.pt[hasgen],
                events.Muon.eta[hasgen],
                events.Muon.phi[hasgen],
                events.Muon.matched_gen.pt[hasgen],
            )
            mc_errsmear = correct_map["roccor"].kSmearMCerror(
                events.Muon.charge[~hasgen],
                events.Muon.pt[~hasgen],
                events.Muon.eta[~hasgen],
                events.Muon.phi[~hasgen],
                events.Muon.nTrackerLayers[~hasgen],
                mc_rand[~hasgen],
            )
            err = np.array(ak.flatten(ak.ones_like(events.Muon.pt)))
            err[hasgen_flat] = np.array(ak.flatten(mc_errspread))
            err[~hasgen_flat] = np.array(ak.flatten(mc_errsmear))
            err = ak.unflatten(err, ak.num(events.Muon.pt))
        muup, mudown = events.Muon, events.Muon
        muup["pt"] = (SF + err) * events.Muon.pt
        mudown["pt"] = (SF - err) * events.Muon.pt
        shifts += [
            (
                {"Jet": shifts[0][0]["Jet"], "MET": shifts[0][0]["MET"], "Muon": muup},
                "RoccorUp",
            )
        ]
        shifts += [
            (
                {
                    "Jet": shifts[0][0]["Jet"],
                    "MET": shifts[0][0]["MET"],
                    "Muon": mudown,
                },
                "RoccorDown",
            )
        ]
    return shifts


def puwei(nPU, correct_map, weights, syst=False):
    """
    <<<<<<< HEAD
        Return pileup weight
        Parameters
        ----------
        nPU: ak.Array
        correct_map : dict
        weights : coffea.analysis_tool.weights
        syst: "split", "weight_only"
    =======
        Apply pileup weights to events based on the number of primary vertices (nPU).

        This function applies pileup weights to the events using the provided correction map and weights.
        It can optionally apply systematic variations.

        Parameters:
        nPU (awkward.Array(int)): The number of primary vertices in the event.
        correct_map (dict): A dictionary containing correction factors and settings for pileup weights.
        weights (): A dictionary to store the calculated weights.
        syst (bool, optional): A flag to indicate whether to apply systematic variations. Default is False.

        Returns:
        None: The function modifies the weights dictionary in place.

        Raises:
        KeyError: If required keys are missing in the correct_map.
        ValueError: If the nPU value is not recognized or supported.
    >>>>>>> doc
    """
    if "correctionlib" in str(type(correct_map["PU"])):
        if syst:
            return weights.add(
                "puweight",
                correct_map["PU"][list(correct_map["PU"].keys())[0]].evaluate(
                    nPU, "nominal"
                ),
                correct_map["PU"][list(correct_map["PU"].keys())[0]].evaluate(
                    nPU, "up"
                ),
                correct_map["PU"][list(correct_map["PU"].keys())[0]].evaluate(
                    nPU, "down"
                ),
            )
        else:
            return weights.add(
                "puweight",
                correct_map["PU"][list(correct_map["PU"].keys())[0]].evaluate(
                    nPU, "nominal"
                ),
            )
    else:
        if syst:
            weights.add(
                "puweight",
                correct_map["PU"]["PU"](nPU),
                correct_map["PU"]["PUup"](nPU),
                correct_map["PU"]["PUdown"](nPU),
            )
        else:
            weights.add("puweight", correct_map["PU"]["PU"](nPU))


def btagSFs(jet, correct_map, weights, SFtype, syst=False):
    """
    Apply b-tagging scale factors (SFs) to a single jet.

    This function applies b-tagging scale factors to the given jet based on the provided correction map, weights, and scale factor type.
    It can optionally apply systematic variations.

    Parameters:
    jet (dict): A dictionary containing the properties of the jet.
    correct_map (dict): A dictionary containing correction factors and settings for b-tagging scale factors.x
    weights (coffea.weight.Weight): An instance of coffea's Weight class to store the calculated weights.
    SFtype (str): The type of scale factor to apply , only shape-based C, B are supported.
    syst (bool, optional): A flag to indicate whether to apply systematic variations. Default is False.

    Returns:
    None: The function modifies the weights instance in place.

    Raises:
    KeyError: If required keys are missing in the correct_map.
    ValueError: If the SFtype is not recognized or supported.
    """
    if SFtype.endswith("C"):
        systlist = [
            "Extrap",
            "Interp",
            "LHEScaleWeight_muF",
            "LHEScaleWeight_muR",
            "PSWeightFSR",
            "PSWeightISR",
            "PUWeight",
            "Stat",
            "XSec_BRUnc_DYJets_b",
            "XSec_BRUnc_DYJets_c",
            "XSec_BRUnc_WJets_c",
            "jer",
            "jesTotal",
        ]
    elif SFtype.endswith("B"):
        systlist = [
            "hf",
            "lf",
            "cferr1",
            "cferr2",
            "hfstat1",
            "hfstat2",
            "lfstats1",
            "lfstats2",
        ]
    sfs_up_all, sfs_down_all = {}, {}
    alljet = jet if jet.ndim > 1 else ak.singletons(jet)
    for i, sys in enumerate(systlist):
        sfs, sfs_down, sfs_up = (
            np.ones_like(alljet[:, 0].pt),
            np.ones_like(alljet[:, 0].pt),
            np.ones_like(alljet[:, 0].pt),
        )
        for nj in range(ak.num(alljet.pt)[0]):
            jet = alljet[:, nj]
            masknone = ak.is_none(jet.pt)
            jet.btagDeepFlavCvL = ak.fill_none(jet.btagDeepFlavCvL, 0.0)
            jet.btagDeepFlavCvB = ak.fill_none(jet.btagDeepFlavCvB, 0.0)
            jet.btagDeepCvL = ak.fill_none(jet.btagDeepCvL, 0.0)
            jet.btagDeepCvB = ak.fill_none(jet.btagDeepCvB, 0.0)
            jet.hadronFlavour = ak.fill_none(jet.hadronFlavour, 0)
            if "correctionlib" in str(type(correct_map["ctag"])):
                if SFtype == "DeepJetC":
                    tmp_sfs = np.where(
                        masknone,
                        1.0,
                        correct_map["ctag"]["deepJet_shape"].evaluate(
                            "central",
                            jet.hadronFlavour,
                            jet.btagDeepFlavCvL,
                            jet.btagDeepFlavCvB,
                        ),
                    )
                    if syst:
                        tmp_sfs_up = np.where(
                            masknone,
                            1.0,
                            correct_map["ctag"]["deepJet_shape"].evaluate(
                                f"up_{systlist[i]}",
                                jet.hadronFlavour,
                                jet.btagDeepFlavCvL,
                                jet.btagDeepFlavCvB,
                            ),
                        )
                        tmp_sfs_down = np.where(
                            masknone,
                            1.0,
                            correct_map["ctag"]["deepJet_shape"].evaluate(
                                f"down_{systlist[i]}",
                                jet.hadronFlavour,
                                jet.btagDeepFlavCvL,
                                jet.btagDeepFlavCvB,
                            ),
                        )
                if SFtype == "DeepCSVC":
                    tmp_sfs = np.where(
                        masknone,
                        1.0,
                        correct_map["ctag"]["deepCSV_shape"].evaluate(
                            "central",
                            jet.hadronFlavour,
                            jet.btagDeepCvL,
                            jet.btagDeepCvB,
                        ),
                    )
                    tmp_sfs_up = np.where(
                        masknone,
                        1.0,
                        correct_map["ctag"]["deepCSV_shape"].evaluate(
                            f"up_{systlist[i]}",
                            jet.hadronFlavour,
                            jet.btagDeepCvL,
                            jet.btagDeepCvB,
                        ),
                    )
                    tmp_sfs_down = np.where(
                        masknone,
                        1.0,
                        correct_map["ctag"]["deepCSV_shape"].evaluate(
                            f"down_{systlist[i]}",
                            jet.hadronFlavour,
                            jet.btagDeepCvL,
                            jet.btagDeepCvB,
                        ),
                    )
            if "correctionlib" in str(type(correct_map["btag"])):
                if SFtype == "DeepJetB":
                    tmp_sfs = np.where(
                        masknone,
                        1.0,
                        correct_map["btag"]["deepJet_shape"].evaluate(
                            "central",
                            jet.hadronFlavour,
                            jet.btagDeepFlavCvL,
                            jet.btagDeepFlavCvB,
                        ),
                    )
                    if syst:
                        tmp_sfs_up = np.where(
                            masknone,
                            1.0,
                            correct_map["btag"]["deepJet_shape"].evaluate(
                                f"up_{systlist[i]}",
                                jet.hadronFlavour,
                                jet.btagDeepFlavCvL,
                                jet.btagDeepFlavCvB,
                            ),
                        )
                        tmp_sfs_down = np.where(
                            masknone,
                            1.0,
                            correct_map["btag"]["deepJet_shape"].evaluate(
                                f"down_{systlist[i]}",
                                jet.hadronFlavour,
                                jet.btagDeepFlavCvL,
                                jet.btagDeepFlavCvB,
                            ),
                        )
                if SFtype == "DeepCSVB":
                    tmp_sfs = np.where(
                        masknone,
                        1.0,
                        correct_map["btag"]["deepCSV_shape"].evaluate(
                            "central",
                            jet.hadronFlavour,
                            jet.btagDeepCvL,
                            jet.btagDeepCvB,
                        ),
                    )
                    tmp_sfs_up = np.where(
                        masknone,
                        1.0,
                        correct_map["btag"]["deepCSV_shape"].evaluate(
                            f"up_{systlist[i]}",
                            jet.hadronFlavour,
                            jet.btagDeepCvL,
                            jet.btagDeepCvB,
                        ),
                    )
                    tmp_sfs_down = np.where(
                        masknone,
                        1.0,
                        correct_map["btag"]["deepCSV_shape"].evaluate(
                            f"down_{systlist[i]}",
                            jet.hadronFlavour,
                            jet.btagDeepCvL,
                            jet.btagDeepCvB,
                        ),
                    )

            sfs = sfs * tmp_sfs
            if syst:
                sfs_up = sfs_up * tmp_sfs_up
                sfs_down = sfs_down * tmp_sfs_down

        if i == 0 and syst == False:
            weights.add(SFtype, sfs)
            break
        else:
            sfs_up_all[sys] = sfs_up
            sfs_down_all[sys] = sfs_down
    if syst == True:
        weights.add_multivariation(
            SFtype,
            sfs,
            systlist,
            np.array(list(sfs_up_all.values())),
            np.array(list(sfs_down_all.values())),
        )
    return weights


### Lepton SFs
def eleSFs(ele, correct_map, weights, syst=True, isHLT=False):
    allele = ele if ele.ndim > 1 else ak.singletons(ele)

    for sf in correct_map["EGM_cfg"].keys():
        ## Only apply SFs for lepton pass HLT filter
        if not isHLT and "HLT" in sf:
            continue
        sf_type = sf[: sf.find(" ")]
        if "low" in sf or "high" in sf:
            continue
        for nele in range(ak.num(allele.pt)[0]):
            ele = allele[:, nele]
            ele_eta = ak.fill_none(ele.eta, -2.5)
            ele_pt = ak.fill_none(ele.pt, 20)
            mask = ele.pt > 20.0
            masknone = ak.is_none(ele.pt)
            sfs_alle, sfs_alle_up, sfs_alle_down = (
                np.ones_like(allele[:, 0].pt),
                np.ones_like(allele[:, 0].pt),
                np.ones_like(allele[:, 0].pt),
            )

            if "correctionlib" in str(type(correct_map["EGM"])):
                ## Reco SF -splitted pT
                if "Reco" in sf:
                    ## phi is used in Summer23
                    ele_pt = np.clip(ele.pt, 20.1, 74.9)
                    ele_pt_low = np.where(ele.pt >= 20.0, 19.9, ele.pt)
                    ele_pt_high = np.clip(ele.pt, 75.0, 500.0)
                    if "Summer23" in correct_map["campaign"]:
                        sfs_low = np.where(
                            (ele.pt <= 20.0) & (~masknone),
                            correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                sf.split(" ")[1],
                                "sf",
                                "RecoBelow20",
                                ele_eta,
                                ele_pt_low,
                                ele.phi,
                            ),
                            1.0,
                        )
                        sfs_high = np.where(
                            (ele.pt > 75.0) & (~masknone),
                            correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                sf.split(" ")[1],
                                "sf",
                                "RecoAbove75",
                                ele_eta,
                                ele_pt_high,
                                ele.phi,
                            ),
                            sfs_low,
                        )
                        sfs = np.where(
                            (ele.pt > 20.0) & (ele.pt <= 75.0) & (~masknone),
                            correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                sf.split(" ")[1],
                                "sf",
                                "Reco20to75",
                                ele_eta,
                                ele_pt,
                                ele.phi,
                            ),
                            sfs_high,
                        )

                        sfs = np.where(masknone, 1.0, sfs)

                        if syst:
                            sfs_up_low = np.where(
                                (ele.pt <= 20.0) & ~masknone,
                                correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                    sf.split(" ")[1],
                                    "sfup",
                                    "RecoBelow20",
                                    ele_eta,
                                    ele_pt_low,
                                    ele.phi,
                                ),
                                0.0,
                            )
                            sfs_down_low = np.where(
                                (ele.pt <= 20.0) & ~masknone,
                                correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                    sf.split(" ")[1],
                                    "sfdown",
                                    "RecoBelow20",
                                    ele_eta,
                                    ele_pt_low,
                                    ele.phi,
                                ),
                                0.0,
                            )
                            sfs_up_high = np.where(
                                (ele.pt > 20.0) & (ele.pt <= 75.0) & ~masknone,
                                correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                    sf.split(" ")[1],
                                    "sfup",
                                    "RecoAbove75",
                                    ele_eta,
                                    ele_pt_high,
                                    ele.phi,
                                ),
                                sfs_up_low,
                            )
                            sfs_down_high = np.where(
                                (ele.pt > 20.0) & (ele.pt <= 75.0) & ~masknone,
                                correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                    sf.split(" ")[1],
                                    "sfdown",
                                    "RecoAbove75",
                                    ele_eta,
                                    ele_pt_high,
                                    ele.phi,
                                ),
                                sfs_down_low,
                            )
                            sfs_up = np.where(
                                (ele.pt > 20.0) & (ele.pt <= 75.0) & ~masknone,
                                correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                    sf.split(" ")[1],
                                    "sfup",
                                    "Reco20to75",
                                    ele_eta,
                                    ele_pt,
                                    ele.phi,
                                ),
                                sfs_up_high,
                            )
                            sfs_down = np.where(
                                (ele.pt > 20.0) & (ele.pt <= 75.0) & ~masknone,
                                correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                    sf.split(" ")[1],
                                    "sfdown",
                                    "Reco20to75",
                                    ele_eta,
                                    ele_pt,
                                    ele.phi,
                                ),
                                sfs_down_high,
                            )
                            sfs_up, sfs_down = np.where(
                                masknone, 1.0, sfs_up
                            ), np.where(masknone, 1.0, sfs_down)

                    else:

                        sfs_low = np.where(
                            (ele.pt <= 20.0) & (~masknone),
                            correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                sf.split(" ")[1],
                                "sf",
                                "RecoBelow20",
                                ele_eta,
                                ele_pt_low,
                            ),
                            1.0,
                        )
                        sfs_high = np.where(
                            (ele.pt > 75.0) & (~masknone),
                            correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                sf.split(" ")[1],
                                "sf",
                                "RecoAbove75",
                                ele_eta,
                                ele_pt_high,
                            ),
                            sfs_low,
                        )
                        sfs = np.where(
                            (ele.pt > 20.0) & (ele.pt <= 75.0) & (~masknone),
                            correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                sf.split(" ")[1], "sf", "Reco20to75", ele_eta, ele_pt
                            ),
                            sfs_high,
                        )

                        sfs = np.where(masknone, 1.0, sfs)

                        if syst:
                            sfs_up_low = np.where(
                                (ele.pt <= 20.0) & ~masknone,
                                correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                    sf.split(" ")[1],
                                    "sfup",
                                    "RecoBelow20",
                                    ele_eta,
                                    ele_pt_low,
                                ),
                                0.0,
                            )
                            sfs_down_low = np.where(
                                (ele.pt <= 20.0) & ~masknone,
                                correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                    sf.split(" ")[1],
                                    "sfdown",
                                    "RecoBelow20",
                                    ele_eta,
                                    ele_pt_low,
                                ),
                                0.0,
                            )
                            sfs_up_high = np.where(
                                (ele.pt > 20.0) & (ele.pt <= 75.0) & ~masknone,
                                correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                    sf.split(" ")[1],
                                    "sfup",
                                    "RecoAbove75",
                                    ele_eta,
                                    ele_pt_high,
                                ),
                                sfs_up_low,
                            )
                            sfs_down_high = np.where(
                                (ele.pt > 20.0) & (ele.pt <= 75.0) & ~masknone,
                                correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                    sf.split(" ")[1],
                                    "sfdown",
                                    "RecoAbove75",
                                    ele_eta,
                                    ele_pt_high,
                                ),
                                sfs_down_low,
                            )
                            sfs_up = np.where(
                                (ele.pt > 20.0) & (ele.pt <= 75.0) & ~masknone,
                                correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                    sf.split(" ")[1],
                                    "sfup",
                                    "Reco20to75",
                                    ele_eta,
                                    ele_pt,
                                ),
                                sfs_up_high,
                            )
                            sfs_down = np.where(
                                (ele.pt > 20.0) & (ele.pt <= 75.0) & ~masknone,
                                correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                    sf.split(" ")[1],
                                    "sfdown",
                                    "Reco20to75",
                                    ele_eta,
                                    ele_pt,
                                ),
                                sfs_down_high,
                            )

                            sfs_up, sfs_down = np.where(
                                masknone, 1.0, sfs_up
                            ), np.where(masknone, 1.0, sfs_down)
                ## Other files
                else:
                    if "Summer23" in correct_map["campaign"]:
                        sfs = np.where(
                            masknone,
                            1.0,
                            correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                sf.split(" ")[1],
                                "sf",
                                correct_map["EGM_cfg"][sf],
                                ele_eta,
                                ele_pt,
                                ele.phi,
                            ),
                        )

                        if syst:
                            sfs_up = np.where(
                                masknone,
                                1.0,
                                correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                    sf.split(" ")[1],
                                    "sfup",
                                    correct_map["EGM_cfg"][sf],
                                    ele_eta,
                                    ele_pt,
                                    ele.phi,
                                ),
                            )
                            sfs_down = np.where(
                                masknone,
                                1.0,
                                correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                    sf.split(" ")[1],
                                    "sfdown",
                                    correct_map["EGM_cfg"][sf],
                                    ele_eta,
                                    ele_pt,
                                    ele.phi,
                                ),
                            )
                    else:
                        sfs = np.where(
                            masknone,
                            1.0,
                            correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                sf.split(" ")[1],
                                "sf",
                                correct_map["EGM_cfg"][sf],
                                ele_eta,
                                ele_pt,
                            ),
                        )

                        if syst:
                            sfs_up = np.where(
                                masknone,
                                1.0,
                                correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                    sf.split(" ")[1],
                                    "sfup",
                                    correct_map["EGM_cfg"][sf],
                                    ele_eta,
                                    ele_pt,
                                ),
                            )
                            sfs_down = np.where(
                                masknone,
                                1.0,
                                correct_map["EGM"][sf.split(" ")[2]].evaluate(
                                    sf.split(" ")[1],
                                    "sfdown",
                                    correct_map["EGM_cfg"][sf],
                                    ele_eta,
                                    ele_pt,
                                ),
                            )

            else:
                if "ele_Trig" in sf:
                    sfs = np.where(
                        masknone, 1.0, correct_map["EGM_custom"][sf_type](ele_pt)
                    )
                    if syst:
                        sfs_up = np.where(
                            masknone,
                            1.0,
                            correct_map["EGM_custom"][sf_type](ele_pt)
                            + correct_map["EGM_custom"][f"{sf_type}_error"](ele_pt),
                        )
                        sfs_down = np.where(
                            masknone,
                            1.0,
                            correct_map["EGM_custom"][sf_type](ele_pt)
                            - correct_map["EGM_custom"][f"{sf_type}_error"](ele_pt),
                        )
                elif "ele" in sf:
                    sfs = np.where(
                        masknone,
                        1.0,
                        correct_map["EGM_custom"][sf_type](ele_eta, ele_pt),
                    )
                    if syst:
                        sfs_up = np.where(
                            masknone,
                            1.0,
                            correct_map["EGM_custom"][sf_type](ele_eta, ele_pt)
                            + correct_map["EGM_custom"][f"{sf_type}_error"](
                                ele_eta, ele_pt
                            ),
                        )
                        sfs_down = np.where(
                            masknone,
                            1.0,
                            correct_map["EGM_custom"][sf_type](ele_eta, ele_pt)
                            - correct_map["EGM_custom"][f"{sf_type}_error"](
                                ele_eta, ele_pt
                            ),
                        )
            sfs_alle = sfs_alle * sfs
            if syst:
                sfs_alle_down = sfs_alle_down * sfs_down
                sfs_alle_up = sfs_alle_up * sfs_up

        if syst:
            weights.add(sf.split(" ")[0], sfs_alle, sfs_alle_up, sfs_alle_down)
        else:
            weights.add(sf.split(" ")[0], sfs_alle)
    return weights


def muSFs(mu, correct_map, weights, syst=False, isHLT=False):
    allmu = mu if mu.ndim > 1 else ak.singletons(mu)
    for sf in correct_map["MUO_cfg"].keys():
        ## Only apply SFs for lepton pass HLT filter
        if not isHLT and "HLT" in sf:
            continue
        if "low" in sf:
            continue
        sfs_allmu, sfs_allmu_up, sfs_allmu_down = (
            np.ones_like(allmu[:, 0].pt),
            np.ones_like(allmu[:, 0].pt),
            np.ones_like(allmu[:, 0].pt),
        )
        sf_type = sf[: sf.find(" ")]
        for nmu in range(ak.num(allmu.pt)[0]):
            mu = allmu[:, nmu]
            masknone = ak.is_none(mu.pt)

            mu_pt = np.clip(mu.pt, 15.0, 199.9)
            mu_eta = np.clip(np.abs(mu.eta), 0.0, 2.4)
            mask = mu_pt > 30
            sfs = 1.0
            if "correctionlib" in str(type(correct_map["MUO"])):
                sfs = np.where(
                    masknone,
                    1.0,
                    correct_map["MUO"][correct_map["MUO_cfg"][sf]].evaluate(
                        mu_eta, mu_pt, "nominal"
                    ),
                )
                if syst:

                    sf_unc = np.where(
                        masknone,
                        0.0,
                        correct_map["MUO"][correct_map["MUO_cfg"][sf]].evaluate(
                            mu_eta, mu_pt, "syst"
                        ),
                    )
                    sfs_up, sfs_down = 1.0 + sf_unc, 1.0 - sf_unc
            else:
                if "mu" in sf:
                    sfs = np.where(
                        masknone, 1.0, correct_map["MUO_cfg"][sf_type](mu_eta, mu_pt)
                    )
                    if syst:
                        sfs_up = np.where(
                            masknone,
                            1.0,
                            sfs
                            + correct_map["MUO_custom"][f"{sf_type}_error"](
                                mu_eta, mu_pt
                            ),
                        )
                        sf_down = np.where(
                            masknone,
                            1.0,
                            sfs
                            - correct_map["MUO_custom"][f"{sf_type}_error"](
                                mu_eta, mu_pt
                            ),
                        )

            sfs_allmu = sfs_allmu * sfs
            if syst:
                sfs_allmu_down = sfs_allmu_down * sfs_down
                sfs_allmu_up = sfs_allmu_up * sfs_up
        if syst:
            weights.add(sf.split(" ")[0], sfs_allmu, sfs_allmu_up, sfs_allmu_down)
        else:
            weights.add(sf.split(" ")[0], sfs_allmu)
    return weights


def jmar_sf(jet, correct_map, weights, syst=False):
    alljet = jet if jet.ndim > 1 else ak.singletons(jet)
    for sf in correct_map["JMAR_cfg"].keys():
        sfs_all, sfs_all_up, sfs_all_down = (
            np.ones_like(alljet[:, 0].pt),
            np.ones_like(alljet[:, 0].pt),
            np.ones_like(alljet[:, 0].pt),
        )
        for njet in range(ak.num(alljet.pt)[0]):
            jet = alljet[:, njet]
            jet_eta = ak.fill_none(jet.eta, -2.5)
            jet_pt = ak.fill_none(jet.pt, 20)
            masknone = ak.is_none(jet.pt)
            # PU Jet ID applied only jet pT<50GeV
            if sf == "PUJetID_eff":
                jet_pt = np.where(jet.pt > 50, 20, jet.pt)
                masknone = (ak.is_none(jet.pt)) & (jet.pt > 50)
            if "correctionlib" in str(type(correct_map["JMAR"])):
                sfs = np.where(
                    masknone,
                    1.0,
                    correct_map["JMAR"][sf].evaluate(
                        jet_eta,
                        jet_pt,
                        "nom",
                        correct_map["JMAR_cfg"][sf],
                    ),
                )
                if syst:
                    sfs_up = np.where(
                        masknone,
                        1.0,
                        correct_map["JMAR"][sf].evaluate(
                            jet_eta,
                            jet_pt,
                            "up",
                            correct_map["JMAR_cfg"][sf],
                        ),
                    )
                    sfs_down = np.where(
                        masknone,
                        1.0,
                        correct_map["JMAR"][sf].evaluate(
                            jet_eta,
                            jet_pt,
                            "down",
                            correct_map["JMAR_cfg"][sf],
                        ),
                    )
                    sfs_all_up = sfs_up * sfs_all_up
                    sfs_all_down = sfs_down * sfs_all_down
            sfs_all = sfs * sfs_all
        if syst:
            weights.add(sf, sfs_all, sfs_all_up, sfs_all_down)
        else:
            weights.add(sf, sfs_all)


def add_pdf_weight(weights, pdf_weights):
    nom = np.ones(len(weights.weight()))
    up = np.ones(len(weights.weight()))
    down = np.ones(len(weights.weight()))

    # NNPDF31_nnlo_hessian_pdfas
    # https://lhapdfsets.web.cern.ch/current/NNPDF31_nnlo_hessian_pdfas/NNPDF31_nnlo_hessian_pdfas.info
    if pdf_weights is not None and "306000 - 306102" in pdf_weights.__doc__:
        # Hessian PDF weights
        # Eq. 21 of https://arxiv.org/pdf/1510.03865v1.pdf
        arg = pdf_weights[:, 1:-2] - np.ones((len(weights.weight()), 100))
        summed = ak.sum(np.square(arg), axis=1)
        pdf_unc = np.sqrt((1.0 / 99.0) * summed)
        weights.add("PDF_weight", nom, pdf_unc + nom)

        # alpha_S weights
        # Eq. 27 of same ref
        as_unc = 0.5 * (pdf_weights[:, 102] - pdf_weights[:, 101])
        weights.add("aS_weight", nom, as_unc + nom)

        # PDF + alpha_S weights
        # Eq. 28 of same ref
        pdfas_unc = np.sqrt(np.square(pdf_unc) + np.square(as_unc))
        weights.add("PDFaS_weight", nom, pdfas_unc + nom)

    else:
        weights.add("aS_weight", nom, up, down)
        weights.add("PDF_weight", nom, up, down)
        weights.add("PDFaS_weight", nom, up, down)


# https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting#TOP_PAG_corrections_based_on_the
def top_pT_sf_formula(pt):
    return 0.103 * np.exp(-0.0118 * pt) - 0.000134 * pt + 0.973


def top_pT_reweighting(gen):
    #     """
    #     Apply this SF only to TTbar datasets! Updated to latest suggestion
    #     Documentation:
    #         - https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting
    #         - https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting#TOP_PAG_corrections_based_on_the
    #     """
    top = gen[(gen.pdgId == 6) & gen.hasFlags(["isLastCopy"])]
    anti_top = gen[(gen.pdgId == -6) & gen.hasFlags(["isLastCopy"])]
    return np.sqrt(
        top_pT_sf_formula(ak.flatten(top.pt, axis=-1))
        * top_pT_sf_formula(ak.flatten(anti_top.pt, axis=-1))
    )


# Jennet adds PS weights
# https://github.com/andrzejnovak/boostedhiggs/blob/master/boostedhiggs/corrections.py#L88-L108
def add_ps_weight(weights, ps_weights):
    nom = np.ones(len(weights.weight()))
    up_isr = np.ones(len(weights.weight()))
    down_isr = np.ones(len(weights.weight()))
    up_fsr = np.ones(len(weights.weight()))
    down_fsr = np.ones(len(weights.weight()))

    if ps_weights is not None:
        if len(ps_weights[0]) == 4:
            up_isr = ps_weights[:, 0]
            down_isr = ps_weights[:, 2]
            up_fsr = ps_weights[:, 1]
            down_fsr = ps_weights[:, 3]
        # else:
        #   warnings.warn(f"PS weight vector has length {len(ps_weights[0])}")

    weights.add("UEPS_ISR", nom, up_isr, down_isr)
    weights.add("UEPS_FSR", nom, up_fsr, down_fsr)


def add_scalevar_7pt(weights, lhe_weights):
    nom = np.ones(len(weights.weight()))

    if len(lhe_weights) > 0:
        if len(lhe_weights[0]) == 9:
            up = np.maximum.reduce(
                [
                    lhe_weights[:, 0],
                    lhe_weights[:, 1],
                    lhe_weights[:, 3],
                    lhe_weights[:, 5],
                    lhe_weights[:, 7],
                    lhe_weights[:, 8],
                ]
            )
            down = np.minimum.reduce(
                [
                    lhe_weights[:, 0],
                    lhe_weights[:, 1],
                    lhe_weights[:, 3],
                    lhe_weights[:, 5],
                    lhe_weights[:, 7],
                    lhe_weights[:, 8],
                ]
            )
        elif len(lhe_weights[0]) > 1:
            print("Scale variation vector has length ", len(lhe_weights[0]))
    else:
        up = np.ones(len(weights.weight()))
        down = np.ones(len(weights.weight()))

    weights.add("scalevar_7pt", nom, up, down)


def add_scalevar_3pt(weights, lhe_weights):
    nom = np.ones(len(weights.weight()))

    if len(lhe_weights) > 0:
        if len(lhe_weights[0]) == 9:
            up = np.maximum(lhe_weights[:, 0], lhe_weights[:, 8])
            down = np.minimum(lhe_weights[:, 0], lhe_weights[:, 8])
        elif len(lhe_weights[0]) > 1:
            print("Scale variation vector has length ", len(lhe_weights[0]))
    else:
        up = np.ones(len(weights.weight()))
        down = np.ones(len(weights.weight()))

    weights.add("scalevar_3pt", nom, up, down)


# JP calibration utility
class JPCalibHandler(object):
    def __init__(self, year, campaign, isRealData, dataset, isSyst=False):
        """
        A tool for calculating the track probability and jet probability
            campaign: campaign name
            isRealData: whether the dataset is real data
            dataset: dataset name from events.metadata["dataset"]
        """
        if "JPCalib" not in config[campaign].keys():
            templates = uproot.open(
                "src/BTVNanoCommissioning/data/JPCalib/Summer22Run3/calibeHistoWrite_MC2022_NANO130X_v2.root"
            )
        else:
            if isRealData:
                if isSyst is not False:
                    filename = config[campaign]["JPCalib"]["MC"]
                else:
                    filename = "default"
                    for key in config[campaign]["JPCalib"]:
                        if key in dataset:
                            filename = config[campaign]["JPCalib"][key]
                            break
                    if filename == "default":
                        raise ValueError(f"No JPCalib file found for dataset {dataset}")
            else:
                filename = config[campaign]["JPCalib"]["MC"]

            templates = uproot.open(
                f"src/BTVNanoCommissioning/data/JPCalib/{year}_{campaign}/{filename}"
            )
        self.ipsig_histo_val = np.array(
            [templates[f"histoCat{i}"].values() for i in range(10)]
        )
        self.ipsig_histo_tot = np.sum(self.ipsig_histo_val, axis=1)
        self.values_cumsum = np.cumsum(self.ipsig_histo_val[:, ::-1], axis=1)[:, ::-1]
        self.edges = templates["histoCat0"].axes[0].edges()

    def flatten(self, array):
        """
        Get the fully flattened array and its layout for each layer
        """
        layouts = []
        array_fl = array
        while str(ak.type(array_fl)).count("*") > 1:
            layouts.append(ak.num(array_fl))
            array_fl = ak.flatten(array_fl)
        return array_fl, layouts

    def unflatten(self, array_fl, layouts):
        """
        Recover a flattened array using the original layouts
        """

        array = array_fl
        for layout in layouts[::-1]:
            array = ak.unflatten(array, layout)
        return array

    def calc_track_proba(self, ipsig: ak.Array, cat: ak.Array):
        """
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


def weight_manager(pruned_ev, SF_map, isSyst):
    weights = Weights(len(pruned_ev), storeIndividual=True)
    if len(SF_map.keys()) == 0:
        if "genWeight" in pruned_ev.fields:
            weights.add("genweight", pruned_ev.genWeight)
        return weights
    if "hadronFlavour" in pruned_ev.Jet.fields:

        syst_wei = True if isSyst != False else False
        if "PU" in SF_map.keys():
            puwei(
                pruned_ev.Pileup.nTrueInt,
                SF_map,
                weights,
                syst_wei,
            )
        if "MUO" in SF_map.keys() and "SelMuon" in pruned_ev.fields:
            muSFs(pruned_ev.SelMuon, SF_map, weights, syst_wei, False)
        if "EGM" in SF_map.keys() and "SelElectron" in pruned_ev.fields:
            eleSFs(pruned_ev.SelElectron, SF_map, weights, syst_wei, False)
        if "BTV" in SF_map.keys() and "SelJet" in pruned_ev.fields:
            btagSFs(pruned_ev.SelJet, SF_map, weights, "DeepJetC", syst_wei)
            btagSFs(pruned_ev.SelJet, SF_map, weights, "DeepJetB", syst_wei)
            btagSFs(pruned_ev.SelJet, SF_map, weights, "DeepCSVB", syst_wei)
            btagSFs(pruned_ev.SelJet, SF_map, weights, "DeepCSVC", syst_wei)

    return weights
