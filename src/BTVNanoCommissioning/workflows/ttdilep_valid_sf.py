import pickle, os, sys, mplhep as hep, numpy as np
import collections

from matplotlib.pyplot import jet

import coffea

from coffea import processor
import awkward as ak
import gc
import os, psutil
from BTVNanoCommissioning.utils.AK4_parameters import correction_config
from BTVNanoCommissioning.utils.correction import (
    lumiMasks,
    eleSFs,
    muSFs,
    load_pu,
    load_BTV,
    load_jetfactory,
    add_jec_variables,
)
from coffea.analysis_tools import Weights
from BTVNanoCommissioning.helpers.func import flatten
from BTVNanoCommissioning.helpers.cTagSFReader import getSF
from BTVNanoCommissioning.utils.histogrammer import histogrammer


class NanoProcessor(processor.ProcessorABC):
    # Define histograms
    def __init__(self, year="2017", campaign="Rereco17_94X", isCorr=True, isJERC=False):
        self._year = year
        self._campaign = campaign
        self.isCorr = isCorr
        self.isJERC = isJERC
        ## Load corrections
        if isCorr:
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
        if isJERC:
            self._jet_factory = load_jetfactory(
                self._campaign, correction_config[self._campaign]["JME"]
            )
        _hist_event_dict = histogrammer("ttdilep_sf")
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
        # if isRealData:
        #     req_lumi = lumiMasks[self._year](events.run, events.luminosityBlock)
        weights = Weights(len(events), storeIndividual=True)
        if not isRealData:
            weights.add("genweight", events.genWeight)
            if self.isCorr:
                weights.add(
                    "puweight",
                    self._pu[f"{self._year}_pileupweight"](events.Pileup.nPU),
                )
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
        ##############
        # Trigger level
        triggers = [
            # "HLT_IsoMu24",
            "HLT_IsoMu24",
        ]

        trig_arrs = [events.HLT[_trig.strip("HLT_")] for _trig in triggers]
        req_trig = np.zeros(len(events), dtype="bool")
        for t in trig_arrs:
            req_trig = req_trig | t

        ############
        # Event level

        ## Muon cuts
        # muon twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
        events.Muon = events.Muon[
            (events.Muon.pt > 30)
            & (abs(events.Muon.eta < 2.4))
            & (events.Muon.tightId)
            & (events.Muon.pfRelIso04_all < 0.12)
        ]
        events.Muon = ak.pad_none(events.Muon, 1, axis=1)
        req_muon = ak.count(events.Muon.pt, axis=1) == 1

        # ## Electron cuts
        # # electron twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
        events.Electron = events.Electron[
            (events.Electron.pt > 30)
            & (abs(events.Electron.eta) < 2.4)
            & (events.Electron.cutBased > 3)
        ]
        events.Electron = ak.pad_none(events.Electron, 1, axis=1)
        req_ele = ak.count(events.Electron.pt, axis=1) == 1
        ## Jet cuts

        req_opposite_charge = (
            events.Electron[:, 0:1].charge * events.Muon[:, 0:1].charge
        ) == -1
        req_opposite_charge = ak.fill_none(req_opposite_charge, False)
        req_opposite_charge = ak.flatten(req_opposite_charge)
        event_jet = events.Jet[
            (events.Jet.pt > 25)
            & (abs(events.Jet.eta) <= 2.4)
            # & (events.Jet.puId > 0) # commented out due to run3 condition
            & (events.Jet.jetId > 5)
            & (ak.all(events.Jet.metric_table(events.Muon) > 0.4, axis=2))
            & (ak.all(events.Jet.metric_table(events.Electron) > 0.4, axis=2))
        ]
        req_jets = ak.num(event_jet.pt) >= 2
        event_level = (
            req_trig & req_lumi & req_muon & req_ele & req_jets & req_opposite_charge
        )
        if len(event_level) > 0:
            event_level = ak.fill_none(event_level, False)
        # Selected
        selev = events[event_level]

        #########

        # Per muon
        mu_eta = abs(selev.Muon.eta) < 2.4
        mu_pt = selev.Muon.pt > 30
        mu_idiso = (selev.Muon.tightId > 0.5) & (selev.Muon.pfRelIso04_all < 0.12)
        mu_level = mu_eta & mu_pt & mu_idiso
        smu = selev.Muon[mu_level]
        smu = ak.pad_none(smu, 1, axis=1)

        # Per Electron
        el_eta = abs(selev.Electron.eta) < 2.4
        el_pt = selev.Electron.pt > 30
        el_idiso = selev.Electron.cutBased > 3
        el_level = el_eta & el_pt & el_idiso
        sel = selev.Electron[el_level]
        sel = ak.pad_none(sel, 1, axis=1)

        # Per jet : https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID
        jet_eta = abs(selev.Jet.eta) <= 2.4
        jet_pt = selev.Jet.pt > 25
        # jet_pu = (selev.Jet.puId > 0) & (selev.Jet.jetId > 5) # commented out for run3
        jet_pu = selev.Jet.jetId > 5
        jet_dr = ak.all(selev.Jet.metric_table(smu) > 0.4, axis=2) & ak.all(
            selev.Jet.metric_table(sel) > 0.4, axis=2
        )
        sel = sel[:, 0]
        smu = smu[:, 0]
        jet_level = jet_pu & jet_eta & jet_pt & jet_dr

        sjets = selev.Jet[jet_level]
        sel_jets = sjets
        sjets = sjets[:, :2]
        # b-tag twiki : https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
        bjet_disc = selev.Jet.btagDeepB > 0.4941
        bjet_level = jet_level & bjet_disc

        sbjets = selev.Jet[bjet_level]
        if not isRealData:
            stbjets = sbjets[sbjets.hadronFlavour == 5]
        if not isRealData and self.isCorr:
            weights.add(
                "lep1sf",
                np.where(
                    event_level,
                    muSFs(
                        ak.firsts(
                            events.Muon[
                                (events.Muon.pt > 30)
                                & (abs(events.Muon.eta < 2.4))
                                & (events.Muon.tightId)
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
                                (events.Muon.pt > 30)
                                & (abs(events.Muon.eta < 2.4))
                                & (events.Muon.tightId)
                            ]
                        ),
                        self._campaign,
                        correction_config[self._campaign]["LSF"],
                    ),
                    1.0,
                ),
            )
        if isRealData:
            genflavor = ak.zeros_like(sjets.pt)
        else:
            par_flav = (sjets.partonFlavour == 0) & (sjets.hadronFlavour == 0)
            genflavor = sjets.hadronFlavour + 1 * par_flav
            jetsfs_c = collections.defaultdict(dict)
            jetsfs_b = collections.defaultdict(dict)
            csvsfs_c = collections.defaultdict(dict)
            csvsfs_b = collections.defaultdict(dict)
            if self.isCorr:
                for i in range(2):
                    jetsfs_c[i]["SF"] = getSF(
                        sjets[:, i].hadronFlavour,
                        sjets[:, i].btagDeepFlavCvL,
                        sjets[:, i].btagDeepFlavCvB,
                        self._deepjetc_sf,
                    )
                    jetsfs_c[i]["SFup"] = getSF(
                        sjets[:, i].hadronFlavour,
                        sjets[:, i].btagDeepFlavCvL,
                        sjets[:, i].btagDeepFlavCvB,
                        self._deepjetc_sf,
                        "TotalUncUp",
                    )
                    jetsfs_c[i]["SFdn"] = getSF(
                        sjets[:, i].hadronFlavour,
                        sjets[:, i].btagDeepFlavCvL,
                        sjets[:, i].btagDeepFlavCvB,
                        self._deepjetc_sf,
                        "TotalUncDown",
                    )
                    jetsfs_b[i]["SF"] = self._deepjetb_sf.eval(
                        "central",
                        sjets[:, i].hadronFlavour,
                        abs(sjets[:, i].eta),
                        sjets[:, i].pt,
                        discr=sjets[:, i].btagDeepFlavB,
                    )
                    jetsfs_b[i]["SFup"] = self._deepjetb_sf.eval(
                        "up_jes",
                        sjets[:, i].hadronFlavour,
                        abs(sjets[:, i].eta),
                        sjets[:, i].pt,
                        discr=sjets[:, i].btagDeepFlavB,
                    )
                    jetsfs_b[i]["SFdn"] = self._deepjetb_sf.eval(
                        "down_jes",
                        sjets[:, i].hadronFlavour,
                        abs(sjets[:, i].eta),
                        sjets[:, i].pt,
                        discr=sjets[:, i].btagDeepFlavB,
                    )
                    csvsfs_c[i]["SF"] = getSF(
                        sjets[:, i].hadronFlavour,
                        sjets[:, i].btagDeepCvL,
                        sjets[:, i].btagDeepCvB,
                        self._deepcsvc_sf,
                    )
                    csvsfs_c[i]["SFup"] = getSF(
                        sjets[:, i].hadronFlavour,
                        sjets[:, i].btagDeepCvL,
                        sjets[:, i].btagDeepCvB,
                        self._deepcsvc_sf,
                        "TotalUncUp",
                    )
                    csvsfs_c[i]["SFdn"] = getSF(
                        sjets[:, i].hadronFlavour,
                        sjets[:, i].btagDeepCvL,
                        sjets[:, i].btagDeepCvB,
                        self._deepcsvc_sf,
                        "TotalUncDown",
                    )
                    csvsfs_b[i]["SFup"] = self._deepcsvb_sf.eval(
                        "up_jes",
                        sjets[:, i].hadronFlavour,
                        abs(sjets[:, i].eta),
                        sjets[:, i].pt,
                        discr=sjets[:, i].btagDeepB,
                    )
                    csvsfs_b[i]["SF"] = self._deepcsvb_sf.eval(
                        "central",
                        sjets[:, i].hadronFlavour,
                        abs(sjets[:, i].eta),
                        sjets[:, i].pt,
                        discr=sjets[:, i].btagDeepB,
                    )
                    csvsfs_b[i]["SFdn"] = self._deepcsvb_sf.eval(
                        "down_jes",
                        sjets[:, i].hadronFlavour,
                        abs(sjets[:, i].eta),
                        sjets[:, i].pt,
                        discr=sjets[:, i].btagDeepB,
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
            elif "btagDeep" in histname:
                for i in range(2):
                    sel_jet = sjets[:, i]
                    if str(i) in histname:
                        h.fill(
                            flav=flatten(genflavor[:, i]),
                            syst="noSF",
                            discr=flatten(
                                np.where(
                                    sel_jet[histname.replace(f"_{i}", "")] < 0,
                                    -0.2,
                                    sel_jet[histname.replace(f"_{i}", "")],
                                )
                            ),
                            weight=weights.weight()[event_level],
                        )
                        if not isRealData and self.isCorr:
                            for syst in disc_list[histname.replace(f"_{i}", "")][
                                i
                            ].keys():
                                h.fill(
                                    flav=flatten(genflavor[:, i]),
                                    syst=syst,
                                    discr=flatten(
                                        np.where(
                                            sel_jet[histname.replace(f"_{i}", "")] < 0,
                                            -0.2,
                                            sel_jet[histname.replace(f"_{i}", "")],
                                        )
                                    ),
                                    weight=weights.weight()[event_level]
                                    * disc_list[histname.replace(f"_{i}", "")][i][syst],
                                )
            elif "mu_" in histname and histname.replace("mu_", "") in smu.fields:

                h.fill(
                    flatten(smu[histname.replace("mu_", "")]),
                    weight=weights.weight()[event_level],
                )
            elif "ele_" in histname and histname.replace("ele_", "") in sel.fields:
                h.fill(
                    flatten(sel[histname.replace("ele_", "")]),
                    weight=weights.weight()[event_level],
                )
            elif "jet" in histname and "dr" not in histname and "njet" != histname:
                for i in range(2):
                    sel_jet = sjets[:, i]
                    if str(i) in histname:
                        h.fill(
                            flatten(genflavor[:, i]),
                            flatten(sel_jet[histname.replace(f"jet{i}_", "")]),
                            weight=weights.weight()[event_level],
                        )

        for i in range(2):
            output[f"dr_mujet{i}"].fill(
                flav=flatten(genflavor[:, i]),
                dr=flatten(smu.delta_r(sjets[:, i])),
                weight=weights.weight()[event_level],
            )
        output["njet"].fill(
            ak.count(sjets.pt, axis=1), weight=weights.weight()[event_level]
        )
        print("end : ", psutil.Process(os.getpid()).memory_info().rss / 1024**2, "MB")

        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator
