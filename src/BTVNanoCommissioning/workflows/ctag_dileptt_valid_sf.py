import numpy as np
import collections
import re
import coffea
from coffea import processor
import awkward as ak
from coffea.analysis_tools import Weights

from BTVNanoCommissioning.utils.correction import (
    lumiMasks,
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
        _hist_event_dict = histogrammer("ctag_ttdilep_sf")
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
            # "HLT_IsoMu24",
            # "HLT_IsoMu27",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8"
        ]

        trig_arrs = [events.HLT[_trig.strip("HLT_")] for _trig in triggers]
        req_trig = np.zeros(len(events), dtype="bool")
        for t in trig_arrs:
            req_trig = req_trig | t

        ############
        # Event level
        ## Muon cuts
        # muon twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
        iso_muon = events.Muon[
            (events.Muon.pt > 12)
            & (abs(events.Muon.eta) < 2.4)
            & (events.Muon.tightId > 0.5)
            & (events.Muon.pfRelIso04_all <= 0.15)
        ]

        iso_muon = ak.pad_none(iso_muon, 2)

        req_muon = (ak.count(iso_muon.pt, axis=1) == 2) & (iso_muon[:, 0].pt > 20)
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
        req_dilepveto = ak.count(dilep_ele.pt, axis=1) == 0

        req_MET = events.MET.pt > 40

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

        dilep_mass = iso_muon[:, 0] + iso_muon[:, 1]
        req_dilepmass = (dilep_mass.mass > 12.0) & (
            (dilep_mass.mass < 75) | (dilep_mass.mass > 105)
        )

        event_level = (
            req_trig
            & req_lumi
            & req_muon
            & req_dilepveto
            & req_dilepmass
            & req_MET
            & req_jets
            & req_softmu
            & req_mujet
        )

        if len(event_level) > 0:
            event_level = ak.fill_none(event_level, False)
        # Selected
        selev = events[event_level]
        #########

        ## Hard Muon
        shmu = selev.Muon[
            (selev.Muon.pt > 12)
            & (abs(selev.Muon.eta) < 2.4)
            & (selev.Muon.tightId > 0.5)
            & (selev.Muon.pfRelIso04_all <= 0.15)
        ]

        ## Soft Muon
        ssmu = selev.Muon[
            (selev.Muon.pt < 25)
            & (abs(selev.Muon.eta) < 2.4)
            & (selev.Muon.tightId > 0.5)
            & (selev.Muon.pfRelIso04_all > 0.2)
        ]
        softmu0 = ssmu[:, 0]
        sz = shmu[:, 0] + shmu[:, 1]

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
                                (events.Muon.pt < 25)
                                & (abs(events.Muon.eta) < 2.4)
                                & (events.Muon.tightId > 0.5)
                                & (events.Muon.pfRelIso04_all > 0.2)
                            ]
                        ),
                        self._campaign,
                        correction_config[self._campaign]["LSF"],
                    ),
                    1.0,
                ),
            )
        isomu0 = shmu[:, 0]
        isomu1 = shmu[:, 1]

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
            smpu = (smuon_jet.partonFlavour == 0) & (smuon_jet.hadronFlavour == 0)
            smflav = 1 * smpu + smuon_jet.hadronFlavour
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
                smuon_jet.btagDeepCvL,
                smuon_jet.btagDeepCvB,
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
                jetsfs_c[1]["SFup"] = getSF(
                    sjets[:, 1].hadronFlavour,
                    sjets[:, 1].btagDeepFlavCvL,
                    sjets[:, 1].btagDeepFlavCvB,
                    self._deepjetc_sf,
                    "TotalUncUp",
                )
                jetsfs_c[1]["SFdn"] = getSF(
                    smuon_jet.hadronFlavour,
                    smuon_jet.btagDeepFlavCvL,
                    smuon_jet.btagDeepFlavCvB,
                    self._deepjetc_sf,
                    "TotalUncDown",
                )
                jetsfs_b[1]["SF"] = self._deepjetb_sf.eval(
                    "central",
                    sjets[:, 1].hadronFlavour,
                    abs(sjets[:, 1].eta),
                    sjets[:, 1].pt,
                    discr=sjets[:, 1].btagDeepFlavB,
                )
                jetsfs_b[1]["SFup"] = self._deepjetb_sf.eval(
                    "up_jes",
                    sjets[:, 1].hadronFlavour,
                    abs(sjets[:, 1].eta),
                    sjets[:, 1].pt,
                    discr=sjets[:, 1].btagDeepFlavB,
                )
                jetsfs_b[1]["SFdn"] = self._deepjetb_sf.eval(
                    "down_jes",
                    sjets[:, 1].hadronFlavour,
                    abs(sjets[:, 1].eta),
                    sjets[:, 1].pt,
                    discr=sjets[:, 1].btagDeepFlavB,
                )
                csvsfs_c[1]["SF"] = getSF(
                    smuon_jet.hadronFlavour,
                    smuon_jet.btagDeepCvL,
                    smuon_jet.btagDeepCvB,
                    self._deepcsvc_sf,
                )
                csvsfs_c[1]["SFup"] = getSF(
                    smuon_jet.hadronFlavour,
                    smuon_jet.btagDeepCvL,
                    smuon_jet.btagDeepCvB,
                    self._deepcsvc_sf,
                    "TotalUncUp",
                )
                csvsfs_c[1]["SFdn"] = getSF(
                    smuon_jet.hadronFlavour,
                    smuon_jet.btagDeepCvL,
                    smuon_jet.btagDeepCvB,
                    self._deepcsvc_sf,
                    "TotalUncDown",
                )
                csvsfs_b[1]["SF"] = self._deepcsvb_sf.eval(
                    "central",
                    sjets[:, 1].hadronFlavour,
                    abs(sjets[:, 1].eta),
                    sjets[:, 1].pt,
                    discr=sjets[:, 1].btagDeepB,
                )
                csvsfs_b[1]["SFup"] = self._deepcsvb_sf.eval(
                    "up_jes",
                    sjets[:, 1].hadronFlavour,
                    abs(sjets[:, 1].eta),
                    sjets[:, 1].pt,
                    discr=sjets[:, 1].btagDeepB,
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
            elif "hl_" in histname and not "ptratio" in histname:
                h.fill(
                    flatten(isomu0[histname.replace("hl_", "")]),
                    weight=weights.weight()[event_level],
                )
            elif "sl_" in histname and not "ptratio" in histname:
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
        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator
