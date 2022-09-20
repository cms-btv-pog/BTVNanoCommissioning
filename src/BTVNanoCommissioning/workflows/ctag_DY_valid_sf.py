import gzip
import pickle, os, sys, mplhep as hep, numpy as np
import collections

from matplotlib.pyplot import jet

import coffea
from coffea import processor
import awkward as ak
from coffea.analysis_tools import Weights
from BTVNanoCommissioning.utils.correction import (
    lumiMasks,
    eleSFs,
    muSFs,
    load_pu,
    load_BTV,
    load_jetfactory,
    add_jec_variables,
)
from BTVNanoCommissioning.utils.AK4_parameters import correction_config
from BTVNanoCommissioning.helpers.func import flatten
from BTVNanoCommissioning.helpers.cTagSFReader import getSF
from BTVNanoCommissioning.utils.histogrammer import histogrammer


class NanoProcessor(processor.ProcessorABC):
    # Define histograms

    def num(ar):
        return ak.num(ak.fill_none(ar[~ak.is_none(ar)], 0), axis=0)

    def __init__(self, year="2017", campaign="Rereco17_94X", isCorr=True, isJERC=False):
        self._year = year
        self._campaign = campaign
        ## Load corrections
        self._deepjetc_sf = load_BTV(
            self._campaign, correction_config[self._campaign]["BTV"], "DeepJetC"
        )
        self._deepjetb_sf = load_BTV(
            self._campaign, correction_config[self._campaign]["BTV"], "DeepJetB"
        )
        self._deepcsvc_sf = load_BTV(
            self._campaign, correction_config[self._campaign]["BTV"], "DeepCSVC"
        )
        self._deepcsvb_sf = load_BTV(
            self._campaign, correction_config[self._campaign]["BTV"], "DeepCSVB"
        )

        self._pu = load_pu(self._campaign, correction_config[self._campaign]["PU"])
        self.isJERC = isJERC
        self.isCorr = isCorr
        if isJERC:
            self._jet_factory = load_jetfactory(
                self._campaign, correction_config[self._campaign]["JME"]
            )
        _hist_event_dict = histogrammer("ctag_DY_sf")
        self.make_output = lambda: {
            "sumw": processor.defaultdict_accumulator(float),
            **_hist_event_dict,
        }

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.make_output()
        dataset = events.metadata["dataset"]
        isRealData = not hasattr(events, "genWeight")
        ## Define the CvL, CvB
        if not hasattr(events, "btagDeepFlavCvL"):
            events.Jet["btagDeepFlavCvL"] = np.maximum(
                np.minimum(
                    np.where(
                        (
                            (
                                events.Jet.btagDeepFlavC
                                / (1.0 - events.Jet.btagDeepFlavB)
                            )
                            > 0
                        )
                        & (events.Jet.pt > 15),
                        (events.Jet.btagDeepFlavC / (1.0 - events.Jet.btagDeepFlavB)),
                        -1,
                    ),
                    0.999999,
                ),
                -1,
            )
            events.Jet["btagDeepFlavCvB"] = np.maximum(
                np.minimum(
                    np.where(
                        (
                            (
                                events.Jet.btagDeepFlavC
                                / (events.Jet.btagDeepFlavC + events.Jet.btagDeepFlavB)
                            )
                            > 0
                        )
                        & (events.Jet.pt > 15),
                        (
                            events.Jet.btagDeepFlavC
                            / (events.Jet.btagDeepFlavC + events.Jet.btagDeepFlavB)
                        ),
                        -1,
                    ),
                    0.999999,
                ),
                -1,
            )
            events.Jet["btagDeepCvL"] = np.maximum(
                np.minimum(
                    np.where(
                        (events.Jet.btagDeepC > 0) & (events.Jet.pt > 15),
                        (events.Jet.btagDeepC / (1.0 - events.Jet.btagDeepB)),
                        -1,
                    ),
                    0.999999,
                ),
                -1,
            )
            events.Jet["btagDeepCvB"] = np.maximum(
                np.minimum(
                    np.where(
                        (events.Jet.btagDeepC > 0) & (events.Jet.pt > 15),
                        (
                            events.Jet.btagDeepC
                            / (events.Jet.btagDeepC + events.Jet.btagDeepB)
                        ),
                        -1,
                    ),
                    0.999999,
                ),
                -1,
            )
        if isRealData:

            output["sumw"] = len(events)
        else:
            output["sumw"] = ak.sum(events.genWeight)
            if self.isJERC:
                events.Jet = self._jet_factory["mc"].build(
                    add_jec_variables(events.Jet, events.fixedGridRhoFastjetAll),
                    lazy_cache=events.caches[0],
                )

        req_lumi = np.ones(len(events), dtype="bool")
        if isRealData:
            req_lumi = lumiMasks[self._year](events.run, events.luminosityBlock)
        weights = Weights(len(events), storeIndividual=True)
        if not isRealData:
            weights.add("genweight", events.genWeight)
            if self.isCorr:
                weights.add(
                    "puweight",
                    self._pu[f"{self._year}_pileupweight"](events.Pileup.nPU),
                )
        ##############
        # Trigger level
        mu_triggers = [
            # "HLT_IsoMu24",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
        ]

        trig_arrs = [events.HLT[_trig.strip("HLT_")] for _trig in mu_triggers]
        req_trig = np.zeros(len(events), dtype="bool")
        for t in trig_arrs:
            req_trig = req_trig | t

        ############
        # Event level
        ## Muon cuts
        dilep_mu = events.Muon[
            (events.Muon.pt > 12)
            & (abs(events.Muon.eta) < 2.4)
            & (events.Muon.tightId > 0.5)
            & (events.Muon.pfRelIso04_all <= 0.15)
        ]
        dilep_ele = events.Electron[
            (events.Electron.pt > 15)
            & (
                (abs(events.Electron.eta) < 1.4442)
                | (
                    (abs(events.Electron.eta) < 2.5)
                    & (abs(events.Electron.eta) > 1.566)
                )
            )
            & (events.Electron.mvaFall17V2Iso_WP80 > 0.5)
        ]
        # req_dilepveto = (ak.count(dilep_mu.pt,axis=1)+ak.count(dilep_ele.pt,axis=1)==2)
        pos_dilep = dilep_mu[dilep_mu.charge > 0]
        neg_dilep = dilep_mu[dilep_mu.charge < 0]
        req_dilep = (
            (ak.num(pos_dilep.pt) >= 1)
            & (ak.num(neg_dilep.pt) >= 1)
            & (ak.num(dilep_mu.charge) >= 2)
            & (ak.num(dilep_ele.charge) < 2)
        )
        pos_dilep = ak.pad_none(pos_dilep, 1, axis=1)
        neg_dilep = ak.pad_none(neg_dilep, 1, axis=1)
        # dilepton mass

        dilep_mass = pos_dilep[:, 0] + neg_dilep[:, 0]
        req_dilepmass = (
            (dilep_mass.mass > 81) & (dilep_mass.mass < 101) & (dilep_mass.pt > 15)
        )

        ## Jet cuts
        event_jet = events.Jet[
            (events.Jet.pt > 20)
            & (abs(events.Jet.eta) <= 2.5)
            & ((events.Jet.puId >= 7) & (events.Jet.pt < 50))
            & (events.Jet.jetId >= 3)
            & (ak.all(events.Jet.metric_table(pos_dilep[:, 0]) > 0.4, axis=2))
            & (ak.all(events.Jet.metric_table(neg_dilep[:, 0]) > 0.4, axis=2))
        ]

        req_jets = ak.num(event_jet.puId) >= 1
        event_jet = ak.pad_none(event_jet, 1, axis=1)
        event_level = req_lumi & req_trig & req_dilep & req_dilepmass & req_jets
        if len(event_level) > 0:
            event_level = ak.fill_none(event_level, False)
        # Selected
        selev = events[event_level]

        #########

        ## Hard Muon
        smu = selev.Muon[
            (selev.Muon.pt > 12)
            & (abs(selev.Muon.eta) < 2.4)
            & (selev.Muon.tightId > 0.5)
            & (selev.Muon.pfRelIso04_all <= 0.15)
        ]
        sposmu = selev.Muon[
            (selev.Muon.pt > 12)
            & (abs(selev.Muon.eta) < 2.4)
            & (selev.Muon.tightId > 0.5)
            & (selev.Muon.pfRelIso04_all <= 0.15)
            & (selev.Muon.charge > 0)
        ]
        sposmu = sposmu[:, 0]
        snegmu = selev.Muon[
            (selev.Muon.pt > 12)
            & (abs(selev.Muon.eta) < 2.4)
            & (selev.Muon.tightId > 0.5)
            & (selev.Muon.pfRelIso04_all <= 0.15)
            & (selev.Muon.charge < 0)
        ]
        snegmu = snegmu[:, 0]

        if not isRealData and self.isCorr:
            weights.add(
                "lep1sf",
                np.where(
                    event_level,
                    muSFs(
                        ak.firsts(
                            events.Muon[
                                (events.Muon.pt > 12)
                                & (abs(events.Muon.eta) < 2.4)
                                & (events.Muon.tightId > 0.5)
                                & (events.Muon.pfRelIso04_all <= 0.15)
                                & (events.Muon.charge < 0)
                            ]
                        ),
                        self._campaign,
                        correction_config[self._campaign]["LSF"],
                    ),
                    1.0,
                ),
            )
            weights.add(
                "lep2sf",
                np.where(
                    event_level,
                    muSFs(
                        ak.firsts(
                            events.Muon[
                                (events.Muon.pt > 12)
                                & (abs(events.Muon.eta) < 2.4)
                                & (events.Muon.tightId > 0.5)
                                & (events.Muon.pfRelIso04_all <= 0.15)
                                & (events.Muon.charge > 0)
                            ]
                        ),
                        self._campaign,
                        correction_config[self._campaign]["LSF"],
                    ),
                    1.0,
                ),
            )
        sz = sposmu + snegmu

        ## Jets

        sjets = selev.Jet[
            (selev.Jet.pt > 20)
            & (abs(selev.Jet.eta) <= 2.5)
            & ((selev.Jet.puId >= 7) & (selev.Jet.pt < 50))
            & (selev.Jet.jetId >= 3)
            & (ak.all(selev.Jet.metric_table(sposmu) > 0.4, axis=2))
            & (ak.all(selev.Jet.metric_table(snegmu) > 0.4, axis=2))
            & (selev.Jet.muEF < 0.8)
        ]
        njet = ak.count(sjets.pt, axis=1)

        def flatten(ar):  # flatten awkward into a 1d array to hist
            return ak.flatten(ar, axis=None)

        if isRealData:
            genflavor = ak.zeros_like(sjets.pt)
        else:
            par_flav = (sjets.partonFlavour == 0) & (sjets.hadronFlavour == 0)
            genflavor = sjets.hadronFlavour + 1 * par_flav
            jetsfs_c = collections.defaultdict(dict)
            jetsfs_b = collections.defaultdict(dict)
            csvsfs_c = collections.defaultdict(dict)
            csvsfs_b = collections.defaultdict(dict)

            ## for each jet
            jetsfs_c[0]["SF"] = getSF(
                sjets[:, 0].hadronFlavour,
                sjets[:, 0].btagDeepFlavCvL,
                sjets[:, 0].btagDeepFlavCvB,
                self._deepjetc_sf,
            )
            jetsfs_c[0]["SFup"] = getSF(
                sjets[:, 0].hadronFlavour,
                sjets[:, 0].btagDeepFlavCvL,
                sjets[:, 0].btagDeepFlavCvB,
                self._deepjetc_sf,
                "TotalUncUp",
            )
            jetsfs_c[0]["SFdn"] = getSF(
                sjets[:, 0].hadronFlavour,
                sjets[:, 0].btagDeepFlavCvL,
                sjets[:, 0].btagDeepFlavCvB,
                self._deepjetc_sf,
                "TotalUncDown",
            )
            jetsfs_b[0]["SF"] = self._deepjetb_sf.eval(
                "central",
                sjets[:, 0].hadronFlavour,
                abs(sjets[:, 0].eta),
                sjets[:, 0].pt,
                discr=sjets[:, 0].btagDeepFlavB,
            )
            jetsfs_b[0]["SFup"] = self._deepjetb_sf.eval(
                "up_jes",
                sjets[:, 0].hadronFlavour,
                abs(sjets[:, 0].eta),
                sjets[:, 0].pt,
                discr=sjets[:, 0].btagDeepFlavB,
            )
            jetsfs_b[0]["SFdn"] = self._deepjetb_sf.eval(
                "down_jes",
                sjets[:, 0].hadronFlavour,
                abs(sjets[:, 0].eta),
                sjets[:, 0].pt,
                discr=sjets[:, 0].btagDeepFlavB,
            )
            csvsfs_c[0]["SF"] = getSF(
                sjets[:, 0].hadronFlavour,
                sjets[:, 0].btagDeepCvL,
                sjets[:, 0].btagDeepCvB,
                self._deepcsvc_sf,
            )
            csvsfs_c[0]["SFup"] = getSF(
                sjets[:, 0].hadronFlavour,
                sjets[:, 0].btagDeepCvL,
                sjets[:, 0].btagDeepCvB,
                self._deepcsvc_sf,
                "TotalUncUp",
            )
            csvsfs_c[0]["SFdn"] = getSF(
                sjets[:, 0].hadronFlavour,
                sjets[:, 0].btagDeepCvL,
                sjets[:, 0].btagDeepCvB,
                self._deepcsvc_sf,
                "TotalUncDown",
            )
            csvsfs_b[0]["SFup"] = self._deepcsvb_sf.eval(
                "up_jes",
                sjets[:, 0].hadronFlavour,
                abs(sjets[:, 0].eta),
                sjets[:, 0].pt,
                discr=sjets[:, 0].btagDeepB,
            )
            csvsfs_b[0]["SF"] = self._deepcsvb_sf.eval(
                "central",
                sjets[:, 0].hadronFlavour,
                abs(sjets[:, 0].eta),
                sjets[:, 0].pt,
                discr=sjets[:, 0].btagDeepB,
            )
            csvsfs_b[0]["SFdn"] = self._deepcsvb_sf.eval(
                "down_jes",
                sjets[:, 0].hadronFlavour,
                abs(sjets[:, 0].eta),
                sjets[:, 0].pt,
                discr=sjets[:, 0].btagDeepB,
            )

            disc_list = {
                "btagDeepB": csvsfs_b,
                "btagDeepC": csvsfs_b,
                "btagDeepFlavB": jetsfs_b,
                "btagDeepFlavC": jetsfs_b,
                "btagDeepCvL": csvsfs_c,
                "btagDeepCvB": csvsfs_c,
                "btagDeepFlavCvL": jetsfs_c,
                "btagDeepFlavCvB": jetsfs_c,
            }
        sel_jet = sjets[:, 0]
        for histname, h in output.items():
            if "Deep" in histname and "btag" not in histname:
                h.fill(
                    flatten(genflavor),
                    flatten(sjets[histname]),
                    weight=flatten(
                        ak.broadcast_arrays(weights.weight()[event_level], sjets["pt"])[
                            0
                        ]
                    ),
                )
            elif "posl_" in histname and histname.replace("posl_", "") in sposmu.fields:
                h.fill(
                    flatten(sposmu[histname.replace("posl_", "")]),
                    weight=weights.weight()[event_level],
                )
            elif "negl_" in histname and histname.replace("negl_", "") in snegmu.fields:
                h.fill(
                    flatten(snegmu[histname.replace("negl_", "")]),
                    weight=weights.weight()[event_level],
                )

            elif "jet_" in histname:
                h.fill(
                    genflavor[:, 0],
                    sel_jet[histname.replace("jet_", "")],
                    weight=weights.weight()[event_level],
                )
            elif "btagDeep" in histname and "0" in histname:

                h.fill(
                    flav=genflavor[:, 0],
                    syst="noSF",
                    discr=np.where(
                        sel_jet[histname.replace("_0", "")] < 0,
                        -0.2,
                        sel_jet[histname.replace("_0", "")],
                    ),
                    weight=weights.weight()[event_level],
                )
                if not isRealData and self.isCorr:
                    for syst in disc_list[histname.replace("_0", "")][0].keys():
                        h.fill(
                            flav=genflavor[:, 0],
                            syst=syst,
                            discr=np.where(
                                sel_jet[histname.replace("_0", "")] < 0,
                                -0.2,
                                sel_jet[histname.replace("_0", "")],
                            ),
                            weight=weights.weight()[event_level]
                            * disc_list[histname.replace("_0", "")][0][syst],
                        )
        output["njet"].fill(njet, weight=weights.weight()[event_level])
        output["dr_mumu"].fill(
            dr=snegmu.delta_r(sposmu), weight=weights.weight()[event_level]
        )
        output["z_pt"].fill(flatten(sz.pt), weight=weights.weight()[event_level])
        output["z_eta"].fill(flatten(sz.eta), weight=weights.weight()[event_level])
        output["z_phi"].fill(flatten(sz.phi), weight=weights.weight()[event_level])
        output["z_mass"].fill(flatten(sz.mass), weight=weights.weight()[event_level])
        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator
