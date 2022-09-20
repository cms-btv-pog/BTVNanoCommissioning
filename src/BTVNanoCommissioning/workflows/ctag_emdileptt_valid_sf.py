import pickle, os, sys, mplhep as hep, numpy as np
import collections
from coffea import processor
import awkward as ak
from coffea.analysis_tools import Weights
import gc
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
        self.isCorr = isCorr
        self.isJERC = isJERC
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
        _hist_event_dict = histogrammer("emctag_ttdilep_sf")
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
        if isRealData:
            req_lumi = lumiMasks[self._year](events.run, events.luminosityBlock)
        weights = Weights(len(events), storeIndividual=True)
        if not hasattr(events, "DeepCSV_FlavCvL"):
            events.Jet["btagDeepFlavCvL"] = np.minimum(
                np.maximum(
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
                    1,
                ),
                -1,
            )
            events.Jet["btagDeepFlavCvB"] = np.minimum(
                np.maximum(
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
                    1,
                ),
                -1,
            )
            events.Jet["btagDeepCvL"] = np.minimum(
                np.maximum(
                    np.where(
                        (events.Jet.btagDeepC > 0) & (events.Jet.pt > 15),
                        (events.Jet.btagDeepC / (1.0 - events.Jet.btagDeepB)),
                        -1,
                    ),
                    1,
                ),
                -1,
            )
            events.Jet["btagDeepCvB"] = np.minimum(
                np.maximum(
                    np.where(
                        (events.Jet.btagDeepC > 0) & (events.Jet.pt > 15),
                        (
                            events.Jet.btagDeepC
                            / (events.Jet.btagDeepC + events.Jet.btagDeepB)
                        ),
                        -1,
                    ),
                    1,
                ),
                -1,
            )
        if not isRealData:
            weights.add("genweight", events.genWeight)
            if self.isCorr:
                weights.add(
                    "puweight",
                    self._pu[f"{self._year}_pileupweight"](events.Pileup.nPU),
                )
        ##############
        # Trigger level
        triggers = [
            "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
            "Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
        ]
        trig_arrs = [events.HLT[_trig] for _trig in triggers]
        req_trig_ele = np.zeros(len(events), dtype="bool")
        for t in trig_arrs:
            req_trig_ele = req_trig_ele | t
        trig_arrsm = [events.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL]
        req_trig_mu = np.zeros(len(events), dtype="bool")
        for t in trig_arrsm:
            req_trig_mu = req_trig_mu | t

        ############
        # Event level
        ## Muon cuts
        # muon twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
        iso_muon_mu = events.Muon[
            (events.Muon.pt > 14)
            & (abs(events.Muon.eta) < 2.4)
            & (events.Muon.tightId > 0.5)
            & (events.Muon.pfRelIso04_all <= 0.15)
        ]
        iso_ele_mu = events.Electron[
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
        iso_muon_ele = events.Muon[
            (events.Muon.pt > 14)
            & (abs(events.Muon.eta) < 2.4)
            & (events.Muon.tightId > 0.5)
            & (events.Muon.pfRelIso04_all <= 0.15)
        ]
        iso_ele_ele = events.Electron[
            (events.Electron.pt > 27)
            & (
                (abs(events.Electron.eta) < 1.4442)
                | (
                    (abs(events.Electron.eta) < 2.5)
                    & (abs(events.Electron.eta) > 1.566)
                )
            )
            & (events.Electron.mvaFall17V2Iso_WP80 > 0.5)
        ]

        iso_muon_ele = ak.pad_none(iso_muon_ele, 1)
        iso_ele_ele = ak.pad_none(iso_ele_ele, 1)
        iso_muon_mu = ak.pad_none(iso_muon_mu, 1)
        iso_ele_mu = ak.pad_none(iso_ele_mu, 1)
        req_ele = (ak.count(iso_muon_ele.pt, axis=1) == 1) & (
            ak.count(iso_ele_ele.pt, axis=1) == 1
        )
        req_mu = (ak.count(iso_muon_mu.pt, axis=1) == 1) & (
            ak.count(iso_ele_mu.pt, axis=1) == 1
        )
        req_MET = events.METFixEE2017.pt > 40

        ## Jet cuts
        event_jet = events.Jet[
            (events.Jet.pt > 20)
            & (abs(events.Jet.eta) <= 2.5)
            & (events.Jet.jetId >= 5)
            & ((events.Jet.pt > 50) | (events.Jet.puId >= 7))
        ]
        req_jets = ak.count(event_jet.pt, axis=1) >= 2

        ## Soft Muon cuts

        soft_muon = events.Muon[
            (events.Muon.pt < 25)
            & (abs(events.Muon.eta) < 2.4)
            & (events.Muon.tightId > 0.5)
            & (events.Muon.pfRelIso04_all > 0.2)
        ]
        req_softmu = ak.count(soft_muon.pt, axis=1) >= 1

        mu_jet = event_jet[
            (ak.all(event_jet.metric_table(soft_muon) <= 0.4, axis=2))
            & ((event_jet.muonIdx1 != -1) | (event_jet.muonIdx2 != -1))
        ]
        req_mujet = ak.count(mu_jet.pt, axis=1) >= 1

        # dilepton mass

        dilep_mass_ele = iso_muon_ele[:, 0] + iso_ele_ele[:, 0]
        req_dilepmass_ele = (dilep_mass_ele.mass > 12.0) & (
            (dilep_mass_ele.mass < 75) | (dilep_mass_ele.mass > 105)
        )
        dilep_mass_mu = iso_muon_mu[:, 0] + iso_ele_mu[:, 0]
        req_dilepmass_mu = (dilep_mass_mu.mass > 12.0) & (
            (dilep_mass_mu.mass < 75) | (dilep_mass_mu.mass > 105)
        )

        event_level = (
            req_lumi
            & req_MET
            & req_jets
            & req_softmu
            & req_mujet
            & (
                (req_trig_ele & req_dilepmass_ele & req_ele)
                | (req_trig_mu & req_dilepmass_mu & req_mu)
            )
        )

        if len(event_level) > 0:
            event_level = ak.fill_none(event_level, False)
        # Selected
        selev = events[event_level]

        #########

        ## Hard Muon
        shmu = selev.Muon
        shele = selev.Electron

        shmu = selev.Muon[
            (
                (selev.Muon.pt > 25)
                & (abs(selev.Muon.eta) < 2.4)
                & (selev.Muon.tightId > 0.5)
                & (selev.Muon.pfRelIso04_all <= 0.15)
                & (selev.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL == True)
            )
            | (
                (selev.Muon.pt > 14)
                & (abs(selev.Muon.eta) < 2.4)
                & (selev.Muon.tightId > 0.5)
                & (selev.Muon.pfRelIso04_all <= 0.15)
                & (selev.HLT.Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL == True)
            )
        ]
        shele = selev.Electron[
            (
                (selev.Electron.pt > 15)
                & (
                    (abs(selev.Electron.eta) < 1.4442)
                    | (
                        (abs(selev.Electron.eta) < 2.5)
                        & (abs(selev.Electron.eta) > 1.566)
                    )
                )
                & (selev.Electron.mvaFall17V2Iso_WP80 > 0.5)
                & (selev.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL == True)
            )
            | (
                (selev.Electron.pt > 27)
                & (
                    (abs(selev.Electron.eta) < 1.4442)
                    | (
                        (abs(selev.Electron.eta) < 2.5)
                        & (abs(selev.Electron.eta) > 1.566)
                    )
                )
                & (selev.Electron.mvaFall17V2Iso_WP80 > 0.5)
                & (selev.HLT.Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL == True)
            )
        ]

        ## Soft Muon
        ssmu = selev.Muon[
            (selev.Muon.pt < 25)
            & (abs(selev.Muon.eta) < 2.4)
            & (selev.Muon.tightId > 0.5)
            & (selev.Muon.pfRelIso04_all > 0.2)
        ]
        softmu0 = ssmu[:, 0]
        sz = shmu[:, 0] + shele[:, 0]
        if not isRealData and self.isCorr:
            weights.add(
                "lep1sf",
                np.where(
                    events.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL
                    & ~events.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,
                    muSFs(
                        ak.firsts(
                            events.Muon[
                                (events.Muon.pt > 25)
                                & (abs(events.Muon.eta) < 2.4)
                                & (events.Muon.tightId > 0.5)
                                & (events.Muon.pfRelIso04_all <= 0.15)
                            ]
                        ),
                        self._campaign,
                        correction_config[self._campaign]["LSF"],
                    ),
                    muSFs(
                        ak.firsts(
                            events.Muon[
                                (events.Muon.pt > 14)
                                & (abs(events.Muon.eta) < 2.4)
                                & (events.Muon.tightId > 0.5)
                                & (events.Muon.pfRelIso04_all <= 0.15)
                            ]
                        ),
                        self._campaign,
                        correction_config[self._campaign]["LSF"],
                    ),
                ),
            )
            weights.add(
                "lep2sf",
                np.where(
                    events.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL
                    & ~events.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,
                    eleSFs(
                        ak.firsts(
                            events.Electron[
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
                        ),
                        self._campaign,
                        correction_config[self._campaign]["LSF"],
                    ),
                    eleSFs(
                        ak.firsts(
                            events.Electron[
                                (events.Electron.pt > 27)
                                & (
                                    (abs(events.Electron.eta) < 1.4442)
                                    | (
                                        (abs(events.Electron.eta) < 2.5)
                                        & (abs(events.Electron.eta) > 1.566)
                                    )
                                )
                                & (events.Electron.mvaFall17V2Iso_WP80 > 0.5)
                            ]
                        ),
                        self._campaign,
                        correction_config[self._campaign]["LSF"],
                    ),
                ),
            )

        isomu0 = shmu[:, 0]
        isomu1 = shele[:, 0]
        ## Jets

        sjets = selev.Jet[
            (selev.Jet.pt > 20)
            & (abs(selev.Jet.eta) <= 2.5)
            & (selev.Jet.jetId >= 5)
            & ((selev.Jet.pt > 50) | (selev.Jet.puId >= 7))
        ]

        ## Muon Jet
        smuon_jet = sjets[
            (ak.all(sjets.metric_table(ssmu) <= 0.4, axis=2))
            & ((sjets.muonIdx1 != -1) | (sjets.muonIdx2 != -1))
        ]
        smuon_jet = smuon_jet[:, 0]

        njet = ak.count(sjets.pt, axis=1)

        def flatten(ar):  # flatten awkward into a 1d array to hist
            return ak.flatten(ar, axis=None)

        if isRealData:
            genflavor = ak.zeros_like(sjets.pt)
            smflav = ak.zeros_like(smuon_jet.pt)
        else:
            par_flav = (sjets.partonFlavour == 0) & (sjets.hadronFlavour == 0)
            genflavor = sjets.hadronFlavour + 1 * par_flav
            smflav = (
                1 * (smuon_jet.partonFlavour == 0)
                & (smuon_jet.hadronFlavour == 0) + smuon_jet.hadronFlavour
            )
            jetsfs_c = collections.defaultdict(dict)
            jetsfs_b = collections.defaultdict(dict)
            csvsfs_c = collections.defaultdict(dict)
            csvsfs_b = collections.defaultdict(dict)
            jetsfs_c[0]["SF"] = getSF(
                smuon_jet.hadronFlavour,
                smuon_jet.btagDeepFlavCvL,
                smuon_jet.btagDeepFlavCvB,
                self._deepjetc_sf,
            )
            jetsfs_c[0]["SFup"] = getSF(
                smuon_jet.hadronFlavour,
                smuon_jet.btagDeepFlavCvL,
                smuon_jet.btagDeepFlavCvB,
                self._deepjetc_sf,
                "TotalUncUp",
            )
            jetsfs_c[0]["SFdn"] = getSF(
                smuon_jet.hadronFlavour,
                smuon_jet.btagDeepFlavCvL,
                smuon_jet.btagDeepFlavCvB,
                self._deepjetc_sf,
                "TotalUncDown",
            )
            jetsfs_b[0]["SF"] = self._deepjetb_sf.eval(
                "central",
                smuon_jet.hadronFlavour,
                abs(smuon_jet.eta),
                smuon_jet.pt,
                discr=smuon_jet.btagDeepFlavB,
            )
            jetsfs_b[0]["SFup"] = self._deepjetb_sf.eval(
                "up_jes",
                smuon_jet.hadronFlavour,
                abs(smuon_jet.eta),
                smuon_jet.pt,
                discr=smuon_jet.btagDeepFlavB,
            )
            jetsfs_b[0]["SFdn"] = self._deepjetb_sf.eval(
                "down_jes",
                smuon_jet.hadronFlavour,
                abs(smuon_jet.eta),
                smuon_jet.pt,
                discr=smuon_jet.btagDeepFlavB,
            )
            csvsfs_c[0]["SF"] = getSF(
                smuon_jet.hadronFlavour,
                smuon_jet.btagDeepFlavCvL,
                smuon_jet.btagDeepFlavCvB,
                self._deepcsvc_sf,
            )
            csvsfs_c[0]["SFup"] = getSF(
                smuon_jet.hadronFlavour,
                smuon_jet.btagDeepCvL,
                smuon_jet.btagDeepCvB,
                self._deepcsvc_sf,
                "TotalUncUp",
            )
            csvsfs_c[0]["SFdn"] = getSF(
                smuon_jet.hadronFlavour,
                smuon_jet.btagDeepCvL,
                smuon_jet.btagDeepCvB,
                self._deepcsvc_sf,
                "TotalUncDown",
            )
            csvsfs_b[0]["SFup"] = self._deepcsvb_sf.eval(
                "up_jes",
                smuon_jet.hadronFlavour,
                abs(smuon_jet.eta),
                smuon_jet.pt,
                discr=smuon_jet.btagDeepB,
            )
            csvsfs_b[0]["SF"] = self._deepcsvb_sf.eval(
                "central",
                smuon_jet.hadronFlavour,
                abs(smuon_jet.eta),
                smuon_jet.pt,
                discr=smuon_jet.btagDeepB,
            )
            csvsfs_b[0]["SFdn"] = self._deepcsvb_sf.eval(
                "down_jes",
                smuon_jet.hadronFlavour,
                abs(smuon_jet.eta),
                smuon_jet.pt,
                discr=smuon_jet.btagDeepB,
            )
            if all(i > 1 for i in njet):
                jetsfs_c[1]["SF"] = getSF(
                    sjets[:, 1].hadronFlavour,
                    sjets[:, 1].btagDeepFlavCvL,
                    sjets[:, 1].btagDeepFlavCvB,
                    self._deepjetc_sf,
                )
                csvsfs_c[1]["SF"] = getSF(
                    smuon_jet.hadronFlavour,
                    smuon_jet.btagDeepCvL,
                    smuon_jet.btagDeepCvB,
                    self._deepcsvc_sf,
                )
                jetsfs_c[1]["SFup"] = getSF(
                    sjets[:, 1].hadronFlavour,
                    sjets[:, 1].btagDeepFlavCvL,
                    sjets[:, 1].btagDeepFlavCvB,
                    self._deepjetc_sf,
                    "TotalUncUp",
                )
                csvsfs_c[1]["SFup"] = getSF(
                    smuon_jet.hadronFlavour,
                    smuon_jet.btagDeepCvL,
                    smuon_jet.btagDeepCvB,
                    self._deepcsvc_sf,
                    "TotalUncUp",
                )
                jetsfs_c[1]["SFdn"] = getSF(
                    sjets[:, 1].hadronFlavour,
                    sjets[:, 1].btagDeepFlavCvL,
                    sjets[:, 1].btagDeepFlavCvB,
                    self._deepjetc_sf,
                    "TotalUncDown",
                )
                csvsfs_c[1]["SFdn"] = getSF(
                    smuon_jet.hadronFlavour,
                    smuon_jet.btagDeepCvL,
                    smuon_jet.btagDeepCvB,
                    self._deepcsvc_sf,
                    "TotalUncDown",
                )
                jetsfs_b[1]["SF"] = self._deepjetb_sf.eval(
                    "central",
                    sjets[:, 1].hadronFlavour,
                    abs(sjets[:, 1].eta),
                    sjets[:, 1].pt,
                    discr=sjets[:, 1].btagDeepFlavB,
                )
                csvsfs_b[1]["SF"] = self._deepcsvb_sf.eval(
                    "central",
                    sjets[:, 1].hadronFlavour,
                    abs(sjets[:, 1].eta),
                    sjets[:, 1].pt,
                    discr=sjets[:, 1].btagDeepB,
                )
                jetsfs_b[1]["SFup"] = self._deepjetb_sf.eval(
                    "up_jes",
                    sjets[:, 1].hadronFlavour,
                    abs(sjets[:, 1].eta),
                    sjets[:, 1].pt,
                    discr=sjets[:, 1].btagDeepFlavB,
                )
                csvsfs_b[1]["SFup"] = self._deepcsvb_sf.eval(
                    "up_jes",
                    sjets[:, 1].hadronFlavour,
                    abs(sjets[:, 1].eta),
                    sjets[:, 1].pt,
                    discr=sjets[:, 1].btagDeepB,
                )
                jetsfs_b[1]["SFdn"] = self._deepjetb_sf.eval(
                    "down_jes",
                    sjets[:, 1].hadronFlavour,
                    abs(sjets[:, 1].eta),
                    sjets[:, 1].pt,
                    discr=sjets[:, 1].btagDeepFlavB,
                )
                csvsfs_b[1]["SFdn"] = self._deepcsvb_sf.eval(
                    "down_jes",
                    sjets[:, 1].hadronFlavour,
                    abs(sjets[:, 1].eta),
                    sjets[:, 1].pt,
                    discr=sjets[:, 1].btagDeepB,
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
            elif "jet_" in histname and "mu" not in histname:
                h.fill(
                    flatten(genflavor),
                    flatten(sjets[histname.replace("jet_", "")]),
                    weight=flatten(
                        ak.broadcast_arrays(weights.weight()[event_level], sjets["pt"])[
                            0
                        ]
                    ),
                )
            elif "hl_" in histname and histname.replace("hl_", "") in isomu0.fields:

                h.fill(
                    flatten(isomu0[histname.replace("hl_", "")]),
                    weight=weights.weight()[event_level],
                )
            elif "sl_" in histname and histname.replace("sl_", "") in isomu1.fields:
                h.fill(
                    flatten(isomu1[histname.replace("sl_", "")]),
                    weight=weights.weight()[event_level],
                )
            elif "soft_l" in histname and not "ptratio" in histname:
                h.fill(
                    flatten(softmu0[histname.replace("soft_l_", "")]),
                    weight=weights.weight()[event_level],
                )
            elif "MET_" in histname:
                h.fill(
                    flatten(selev.METFixEE2017[histname.replace("MET_", "")]),
                    weight=weights.weight()[event_level],
                )
            elif "lmujet_" in histname:
                h.fill(
                    smflav,
                    flatten(smuon_jet[histname.replace("lmujet_", "")]),
                    weight=weights.weight()[event_level],
                )
            elif "btagDeep" in histname and "0" in histname:
                h.fill(
                    flav=smflav,
                    syst="noSF",
                    discr=np.where(
                        smuon_jet[histname.replace("_0", "")] < 0,
                        -0.2,
                        smuon_jet[histname.replace("_0", "")],
                    ),
                    weight=weights.weight()[event_level],
                )
                if not isRealData and self.isCorr:
                    for syst in disc_list[histname.replace("_0", "")][0].keys():
                        h.fill(
                            flav=smflav,
                            syst=syst,
                            discr=np.where(
                                smuon_jet[histname.replace("_0", "")] < 0,
                                -0.2,
                                smuon_jet[histname.replace("_0", "")],
                            ),
                            weight=weights.weight()[event_level]
                            * disc_list[histname.replace("_0", "")][0][syst],
                        )
            elif (
                "btagDeep" in histname and "1" in histname and all(i > 1 for i in njet)
            ):
                sljets = sjets[:, 1]
                h.fill(
                    flav=genflavor[:, 1],
                    syst="noSF",
                    discr=np.where(
                        sljets[histname.replace("_1", "")] < 0,
                        -0.2,
                        sljets[histname.replace("_1", "")],
                    ),
                    weight=weights.weight()[event_level],
                )
                if not isRealData and self.isCorr:
                    for syst in disc_list[histname.replace("_1", "")][1].keys():
                        h.fill(
                            flav=genflavor[:, 1],
                            syst=syst,
                            discr=np.where(
                                sljets[histname.replace("_1", "")] < 0,
                                -0.2,
                                sljets[histname.replace("_1", "")],
                            ),
                            weight=weights.weight()[event_level]
                            * disc_list[histname.replace("_1", "")][1][syst],
                        )
        output["njet"].fill(njet, weight=weights.weight()[event_level])
        output["hl_ptratio"].fill(
            flav=genflavor[:, 0],
            ratio=isomu0.pt / sjets[:, 0].pt,
            weight=weights.weight()[event_level],
        )
        output["sl_ptratio"].fill(
            flav=genflavor[:, 0],
            ratio=isomu1.pt / sjets[:, 0].pt,
            weight=weights.weight()[event_level],
        )
        output["soft_l_ptratio"].fill(
            flav=smflav,
            ratio=softmu0.pt / smuon_jet.pt,
            weight=weights.weight()[event_level],
        )
        output["dr_lmujetsmu"].fill(
            flav=smflav,
            dr=smuon_jet.delta_r(softmu0),
            weight=weights.weight()[event_level],
        )
        output["dr_lmujethmu"].fill(
            flav=smflav,
            dr=smuon_jet.delta_r(isomu0),
            weight=weights.weight()[event_level],
        )
        output["dr_lmusmu"].fill(
            dr=isomu0.delta_r(softmu0),
            weight=weights.weight()[event_level],
        )
        output["z_pt"].fill(flatten(sz.pt), weight=weights.weight()[event_level])
        output["z_eta"].fill(flatten(sz.eta), weight=weights.weight()[event_level])
        output["z_phi"].fill(flatten(sz.phi), weight=weights.weight()[event_level])
        output["z_mass"].fill(flatten(sz.mass), weight=weights.weight()[event_level])
        del disc_list
        gc.collect()

        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator
