import collections, numpy as np, awkward as ak

from coffea import processor
from coffea.analysis_tools import Weights
from BTVNanoCommissioning.utils.correction import (
    load_lumi,
    eleSFs,
    muSFs,
    load_pu,
    load_BTV,
    load_jmefactory,
)
from BTVNanoCommissioning.utils.AK4_parameters import correction_config
from BTVNanoCommissioning.helpers.func import flatten, update
from BTVNanoCommissioning.helpers.update_branch import missing_branch, add_jec
from BTVNanoCommissioning.helpers.cTagSFReader import getSF
from BTVNanoCommissioning.utils.histogrammer import histogrammer
from BTVNanoCommissioning.utils.selection import (
    jet_id,
    mu_idiso,
    ele_mvatightid,
    softmu_mask,
)


class NanoProcessor(processor.ProcessorABC):
    def __init__(self, year="2017", campaign="Rereco17_94X", isCorr=True, isJERC=False):
        self._year = year
        self._campaign = campaign
        ## Load corrections
        self.isCorr = isCorr
        self.isJERC = isJERC
        self.lumiMask = load_lumi(correction_config[self._campaign]["lumiMask"])
        if isCorr:
            if "BTV" in correction_config[self._campaign].keys():
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
            if "PU" in correction_config[self._campaign].keys():
                self._pu = load_pu(
                    self._campaign, correction_config[self._campaign]["PU"]
                )
        if isJERC:
            self._jet_factory = load_jmefactory(
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
        events = missing_branch(events)
        weights = Weights(len(events), storeIndividual=True)

        if self.isJERC:
            add_jec(events, self._campaign, self._jet_factory)
        if isRealData:
            output["sumw"] = len(events)
        else:
            output["sumw"] = ak.sum(events.genWeight)
        ####################
        #    Selections    #
        ####################
        ## Lumimask
        req_lumi = np.ones(len(events), dtype="bool")
        if isRealData:
            req_lumi = self.lumiMask(events.run, events.luminosityBlock)

        ## HLT
        trigger_he = [
            "Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
            "Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
        ]
        trigger_hm = [
            "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
            "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
        ]

        checkHLT = ak.Array(
            [hasattr(events.HLT, _trig) for _trig in trigger_he + trigger_hm]
        )
        if ak.all(checkHLT == False):
            raise ValueError(
                "HLT paths:", trigger_he + trigger_hm, " are all invalid in", dataset
            )
        elif ak.any(checkHLT == False):
            print(
                np.array(trigger_he + trigger_hm)[~checkHLT], " not exist in", dataset
            )
        trig_arr_ele = [
            events.HLT[_trig] for _trig in trigger_he if hasattr(events.HLT, _trig)
        ]
        req_trig_ele = np.zeros(len(events), dtype="bool")
        for t in trig_arr_ele:
            req_trig_ele = req_trig_ele | t
        trig_arr_mu = [
            events.HLT[_trig] for _trig in trigger_hm if hasattr(events.HLT, _trig)
        ]
        req_trig_mu = np.zeros(len(events), dtype="bool")
        for t in trig_arr_mu:
            req_trig_mu = req_trig_mu | t

        ## Muon cuts
        iso_muon_mu = events.Muon[
            (events.Muon.pt > 14) & mu_idiso(events, self._campaign)
        ]
        iso_muon_ele = events.Muon[
            (events.Muon.pt > 14) & mu_idiso(events, self._campaign)
        ]

        ## Electron cuts
        iso_ele_ele = events.Electron[
            (events.Electron.pt > 27) & ele_mvatightid(events, self._campaign)
        ]
        iso_ele_mu = events.Electron[
            (events.Electron.pt > 15) & ele_mvatightid(events, self._campaign)
        ]

        ## cross leptons
        req_ele = (ak.count(iso_muon_ele.pt, axis=1) == 1) & (
            ak.count(iso_ele_ele.pt, axis=1) == 1
        )
        req_mu = (ak.count(iso_muon_mu.pt, axis=1) == 1) & (
            ak.count(iso_ele_mu.pt, axis=1) == 1
        )
        iso_ele = ak.concatenate([iso_ele_mu, iso_ele_ele], axis=1)
        iso_mu = ak.concatenate([iso_muon_mu, iso_muon_ele], axis=1)
        iso_ele = ak.pad_none(iso_ele, 1)
        iso_mu = ak.pad_none(iso_mu, 1)

        ## Jet cuts
        event_jet = events.Jet[jet_id(events, self._campaign)]
        req_jets = ak.count(event_jet.pt, axis=1) >= 2

        ## Soft Muon cuts
        soft_muon = events.Muon[softmu_mask(events, self._campaign)]
        req_softmu = ak.count(soft_muon.pt, axis=1) >= 1

        ## Muon jet cuts
        mu_jet = event_jet[
            (ak.all(event_jet.metric_table(soft_muon) <= 0.4, axis=2))
            & ((event_jet.muonIdx1 != -1) | (event_jet.muonIdx2 != -1))
        ]
        req_mujet = ak.count(mu_jet.pt, axis=1) >= 1

        ## store jet index for PFCands, create mask on the jet index
        jetindx = ak.mask(
            ak.local_index(events.Jet.pt),
            (
                jet_id(events, self._campaign)
                & (ak.all(events.Jet.metric_table(soft_muon) <= 0.4, axis=2))
                & ((events.Jet.muonIdx1 != -1) | (events.Jet.muonIdx2 != -1))
            )
            == 1,
        )
        jetindx = ak.pad_none(jetindx, 1)
        jetindx = jetindx[:, 0]

        ## Other cuts
        req_dilepmass = ((iso_mu[:, 0] + iso_ele[:, 0]).mass > 12.0) & (
            ((iso_mu[:, 0] + iso_ele[:, 0]).mass < 75)
            | ((iso_mu[:, 0] + iso_ele[:, 0]).mass > 105)
        )

        if "Run3" not in self._campaign:
            MET = ak.zip(
                {
                    "pt": events.MET.pt,
                    "eta": ak.zeros_like(events.MET.pt),
                    "phi": events.MET.phi,
                    "mass": ak.zeros_like(events.MET.pt),
                },
                with_name="PtEtaPhiMLorentzVector",
            )
        else:
            MET = ak.zip(
                {
                    "pt": events.PuppiMET.pt,
                    "eta": ak.zeros_like(events.PuppiMET.pt),
                    "phi": events.PuppiMET.phi,
                    "mass": ak.zeros_like(events.PuppiMET.pt),
                },
                with_name="PtEtaPhiMLorentzVector",
            )
        req_MET = MET.pt > 40

        event_level = (
            req_lumi
            & req_MET
            & req_jets
            & req_softmu
            & req_mujet
            & req_dilepmass
            & ((req_trig_ele & req_ele) | (req_trig_mu & req_mu))
        )
        event_level = ak.fill_none(event_level, False)

        ####################
        # Selected objects #
        ####################
        shmu = iso_mu[event_level]
        shele = iso_ele[event_level]
        ssmu = soft_muon[event_level]
        softmu0 = ssmu[:, 0]
        sz = shmu[:, 0] + shele[:, 0]
        isomu0 = shmu[:, 0]
        isomu1 = shele[:, 0]
        sjets = event_jet[event_level]
        smuon_jet = mu_jet[event_level]
        smuon_jet = smuon_jet[:, 0]
        smet = MET[event_level]
        njet = ak.count(sjets.pt, axis=1)
        # Find the PFCands associate with selected jets. Search from jetindex->JetPFCands->PFCand
        if self._campaign != "Rereco17_94X":
            spfcands = events[event_level].PFCands[
                events[event_level]
                .JetPFCands[
                    events[event_level].JetPFCands.jetIdx == jetindx[event_level]
                ]
                .pFCandsIdx
            ]

        ####################
        # Weight & Geninfo #
        ####################
        if not isRealData:
            weights.add("genweight", events.genWeight)
        if not isRealData and self.isCorr:
            if "PU" in correction_config[self._campaign].keys():
                if self._campaign == "Rereco17_94X":
                    puname = f"{self._year}_pileupweight"
                else:
                    puname = "PU"
                weights.add("puweight", self._pu[puname](events.Pileup.nTrueInt))
            if "LSF" in correction_config[self._campaign].keys():
                weights.add(
                    "lep1sf",
                    np.where(
                        event_level,
                        muSFs(
                            ak.firsts(iso_mu),
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
                        eleSFs(
                            ak.firsts(iso_ele),
                            self._campaign,
                            correction_config[self._campaign]["LSF"],
                        ),
                        1.0,
                    ),
                )

        if isRealData:
            genflavor = ak.zeros_like(sjets.pt)
            smflav = ak.zeros_like(smuon_jet.pt)
        else:
            par_flav = (sjets.partonFlavour == 0) & (sjets.hadronFlavour == 0)
            genflavor = sjets.hadronFlavour + 1 * par_flav
            smpu = (smuon_jet.partonFlavour == 0) & (smuon_jet.hadronFlavour == 0)
            smflav = 1 * smpu + smuon_jet.hadronFlavour
            if self.isCorr and "BTV" in correction_config[self._campaign].keys():
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

        ####################
        #  Fill histogram  #
        ####################
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
            elif "PFCands" in histname and self._campaign != "Rereco17_94X":
                h.fill(
                    flatten(ak.broadcast_arrays(smflav, spfcands["pt"])[0]),
                    flatten(spfcands[histname.replace("PFCands_", "")]),
                    weight=flatten(
                        ak.broadcast_arrays(
                            weights.weight()[event_level], spfcands["pt"]
                        )[0]
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
                    smflav,
                    flatten(softmu0[histname.replace("soft_l_", "")]),
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
                if (
                    not isRealData
                    and self.isCorr
                    and "BTV" in correction_config[self._campaign].keys()
                    and "_b" not in histname
                    and "_bb" not in histname
                    and "_lepb" not in histname
                ):
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
                if (
                    not isRealData
                    and self.isCorr
                    and "BTV" in correction_config[self._campaign].keys()
                    and "_b" not in histname
                    and "_bb" not in histname
                    and "_lepb" not in histname
                ):
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
            dr=isomu0.delta_r(softmu0), weight=weights.weight()[event_level]
        )
        output["z_pt"].fill(flatten(sz.pt), weight=weights.weight()[event_level])
        output["z_eta"].fill(flatten(sz.eta), weight=weights.weight()[event_level])
        output["z_phi"].fill(flatten(sz.phi), weight=weights.weight()[event_level])
        output["z_mass"].fill(flatten(sz.mass), weight=weights.weight()[event_level])
        output["MET_pt"].fill(flatten(smet.pt), weight=weights.weight()[event_level])
        output["MET_phi"].fill(flatten(smet.phi), weight=weights.weight()[event_level])
        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator
