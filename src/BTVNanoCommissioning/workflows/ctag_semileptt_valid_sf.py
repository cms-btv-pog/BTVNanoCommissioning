import collections, numpy as np, awkward as ak


from coffea import processor
from coffea.analysis_tools import Weights

from BTVNanoCommissioning.utils.correction import (
    load_lumi,
    load_SF,
    eleSFs,
    muSFs,
    puwei,
    btagSFs,
    load_jmefactory,
)
from BTVNanoCommissioning.helpers.func import flatten, update
from BTVNanoCommissioning.helpers.update_branch import missing_branch, add_jec


from BTVNanoCommissioning.utils.histogrammer import histogrammer
from BTVNanoCommissioning.utils.selection import (
    jet_id,
    mu_idiso,
    ele_mvatightid,
    softmu_mask,
)


class NanoProcessor(processor.ProcessorABC):
    def __init__(
        self,
        year="2017",
        campaign="Rereco17_94X",
        isCorr=True,
        isJERC=False,
        isSyst=False,
    ):
        self._year = year
        self._campaign = campaign
        self.isCorr = isCorr
        self.isJERC = isJERC
        self.isSyst = isSyst
        self.lumiMask = load_lumi(self._campaign)
        ## Load corrections
        if isCorr:
            self.SF_map = load_SF(self._campaign)
        if isJERC:
            self._jet_factory = load_jmefactory(self._campaign)

        ## Load histogram
        _hist_event_dict = histogrammer("ctag_ttsemilep_sf")
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

        if self.isJERC:
            events = add_jec(events, self._campaign, self._jet_factory)
        if isRealData:
            output["sumw"] = len(events)
        else:
            output["sumw"] = ak.sum(events.genWeight)

        ###############
        # Event level #
        ###############
        ## Lumimask
        req_lumi = np.ones(len(events), dtype="bool")
        if isRealData:
            req_lumi = self.lumiMask(events.run, events.luminosityBlock)

        ## HLT
        triggers = ["IsoMu24", "IsoMu27"]
        checkHLT = ak.Array([hasattr(events.HLT, _trig) for _trig in triggers])
        if ak.all(checkHLT == False):
            raise ValueError("HLT paths:", triggers, " are all invalid in", dataset)
        elif ak.any(checkHLT == False):
            print(np.array(triggers)[~checkHLT], " not exist in", dataset)
        trig_arrs = [
            events.HLT[_trig] for _trig in triggers if hasattr(events.HLT, _trig)
        ]
        req_trig = np.zeros(len(events), dtype="bool")
        for t in trig_arrs:
            req_trig = req_trig | t

        ## Muon cuts
        # muon twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
        iso_muon = events.Muon[(events.Muon.pt > 30) & mu_idiso(events, self._campaign)]
        req_muon = ak.count(iso_muon.pt, axis=1) == 1
        iso_muon = ak.pad_none(iso_muon, 1, axis=1)
        iso_muon = iso_muon[:, 0]
        iso_muindx = ak.mask(
            ak.local_index(events.Muon.pt),
            ((events.Muon.pt > 30) & mu_idiso(events, self._campaign)) == 1,
        )
        iso_muindx = ak.pad_none(iso_muindx, 1)
        iso_muindx = iso_muindx[:, 0]

        ## Jet cuts
        event_jet = events.Jet[
            jet_id(events, self._campaign)
            & (
                ak.all(
                    events.Jet.metric_table(iso_muon) > 0.5, axis=2, mask_identity=True
                )
            )
        ]
        req_jets = ak.num(event_jet.pt) >= 4

        ## Soft Muon cuts

        soft_muon = events.Muon[softmu_mask(events, self._campaign)]
        req_softmu = ak.count(soft_muon.pt, axis=1) >= 1
        soft_muon = ak.pad_none(soft_muon, 1, axis=1)

        ## Muon-jet cuts
        mu_jet = event_jet[
            (
                ak.all(
                    event_jet.metric_table(soft_muon) <= 0.4, axis=2, mask_identity=True
                )
            )
            & ((event_jet.muonIdx1 != -1) | (event_jet.muonIdx2 != -1))
            & ((event_jet.muEF + event_jet.neEmEF) < 0.7)
        ]
        req_mujet = ak.num(mu_jet.pt, axis=1) >= 1
        mu_jet = ak.pad_none(mu_jet, 1, axis=1)

        ## store jet index for PFCands, create mask on the jet index
        jetindx = ak.mask(
            ak.local_index(events.Jet.pt),
            (
                jet_id(events, self._campaign)
                & (
                    ak.all(
                        events.Jet.metric_table(iso_muon) > 0.5,
                        axis=2,
                        mask_identity=True,
                    )
                )
                & ((events.Jet.muEF + events.Jet.neEmEF) < 0.7)
                & (
                    ak.all(
                        events.Jet.metric_table(iso_muon) > 0.5,
                        axis=2,
                        mask_identity=True,
                    )
                )
                & (
                    ak.all(
                        events.Jet.metric_table(soft_muon) <= 0.4,
                        axis=2,
                        mask_identity=True,
                    )
                )
                & ((events.Jet.muonIdx1 != -1) | (events.Jet.muonIdx2 != -1))
            )
            == 1,
        )
        jetindx = ak.pad_none(jetindx, 1)
        jetindx = jetindx[:, 0]

        ## Other cuts
        req_pTratio = (soft_muon[:, 0].pt / mu_jet[:, 0].pt) < 0.4

        req_QCDveto = (
            (iso_muon.pfRelIso04_all < 0.05)
            & (abs(iso_muon.dz) < 0.01)
            & (abs(iso_muon.dxy) < 0.002)
            & (iso_muon.sip3d < 2)
            & (
                iso_muon.pt
                / ak.firsts(
                    events.Jet[
                        (events.Jet.muonIdx1 == iso_muindx)
                        | ((events.Jet.muonIdx2 == iso_muindx))
                    ].pt
                )
                > 0.75
            )
        )

        dilep_mu = events.Muon[(events.Muon.pt > 12) & mu_idiso(events, self._campaign)]
        dilep_ele = events.Electron[
            (events.Electron.pt > 15) & ele_mvatightid(events, self._campaign)
        ]
        req_dilepveto = (
            ak.count(dilep_mu.pt, axis=1) + ak.count(dilep_ele.pt, axis=1) != 2
        )

        dilep_mass = iso_muon + soft_muon[:, 0]
        req_dilepmass = (dilep_mass.mass > 12.0) & (
            (dilep_mass.mass < 80) | (dilep_mass.mass > 100)
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
        Wmass = MET + iso_muon
        req_Wmass = Wmass.mass > 55

        event_level = (
            req_trig
            & req_lumi
            & req_muon
            & req_jets
            & req_softmu
            & req_dilepmass
            & req_mujet
            & req_Wmass
            & req_dilepveto
            & req_QCDveto
            & req_pTratio
        )

        event_level = ak.fill_none(event_level, False)
        if len(events[event_level]) == 0:
            return {dataset: output}

        ####################
        # Selected objects #
        ####################
        shmu = iso_muon[event_level]
        sjets = event_jet[event_level]
        ssmu = soft_muon[event_level]
        smet = MET[event_level]
        smuon_jet = mu_jet[event_level]
        nsoftmu = ak.count(ssmu.pt, axis=1)
        nmujet = ak.count(smuon_jet.pt, axis=1)
        smuon_jet = smuon_jet[:, 0]
        ssmu = ssmu[:, 0]
        sz = shmu + ssmu
        sw = shmu + smet
        osss = shmu.charge * ssmu.charge * -1
        njet = ak.count(sjets.pt, axis=1)
        # Find the PFCands associate with selected jets. Search from jetindex->JetPFCands->PFCand
        if "PFCands" in events.fields:
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
        weights = Weights(len(events[event_level]), storeIndividual=True)
        if not isRealData:
            weights.add("genweight", events[event_level].genWeight)
        if not isRealData and self.isCorr:
            if "PU" in self.SF_map.keys():
                weights.add(
                    "puweight", puwei(self.SF_map, events[event_level].Pileup.nTrueInt)
                )
            if "MUO" in self.SF_map.keys() or "EGM" in self.SF_map.keys():
                weights.add("lep1sf", muSFs(shmu, self.SF_map, True))
                weights.add("lep2sf", muSFs(ssmu, self.SF_map))

        if isRealData:
            genflavor = ak.zeros_like(sjets.pt)
            smflav = ak.zeros_like(smuon_jet.pt)
        else:
            genflavor = sjets.hadronFlavour + 1 * (
                (sjets.partonFlavour == 0) & (sjets.hadronFlavour == 0)
            )
            smflav = smuon_jet.hadronFlavour + 1 * (
                (smuon_jet.partonFlavour == 0) & (smuon_jet.hadronFlavour == 0)
            )
            jetsfs_c = collections.defaultdict(dict)
            jetsfs_b = collections.defaultdict(dict)
            csvsfs_c = collections.defaultdict(dict)
            csvsfs_b = collections.defaultdict(dict)
            if self.isCorr and (
                "btag" in self.SF_map.keys() or "ctag" in self.SF_map.keys()
            ):
                jetsfs_c[0]["SF"] = btagSFs(smuon_jet, self.SF_map, "DeepJetC")
                jetsfs_b[0]["SF"] = btagSFs(smuon_jet, self.SF_map, "DeepJetB")
                csvsfs_c[0]["SF"] = btagSFs(smuon_jet, self.SF_map, "DeepCSVC")
                csvsfs_b[0]["SF"] = btagSFs(smuon_jet, self.SF_map, "DeepCSVB")

                if self.isSyst:
                    for syst in [
                        "hf",
                        "lf",
                        "cferr1",
                        "cferr2",
                        "hfstat1",
                        "hfstat2",
                        "lfstats1",
                        "lfstats2",
                    ]:
                        jetsfs_c[0][f"SF_{syst}_up"] = btagSFs(
                            smuon_jet, self.SF_map, "DeepJetC", f"up_{syst}"
                        )
                        jetsfs_c[0][f"SF_{syst}_dn"] = btagSFs(
                            smuon_jet, self.SF_map, "DeepJetC", f"down_{syst}"
                        )
                        csvsfs_c[0][f"SF_{syst}_up"] = btagSFs(
                            smuon_jet, self.SF_map, "DeepCSVC", f"up_{syst}"
                        )
                        csvsfs_c[0][f"SF_{syst}_dn"] = btagSFs(
                            smuon_jet, self.SF_map, "DeepCSVC", f"down_{syst}"
                        )
                    csvsfs_b[0][f"SF_{syst}_up"] = btagSFs(
                        smuon_jet, self.SF_map, "DeepCSVB", f"up"
                    )
                    csvsfs_b[0][f"SF_{syst}_dn"] = btagSFs(
                        smuon_jet, self.SF_map, "DeepCSVB", f"down"
                    )
                    jetsfs_b[0][f"SF_{syst}_up"] = btagSFs(
                        smuon_jet, self.SF_map, "DeepJetB", f"up"
                    )
                    jetsfs_b[0][f"SF_{syst}_dn"] = btagSFs(
                        smuon_jet, self.SF_map, "DeepJetB", f"down"
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
            if (
                "Deep" in histname
                and "btag" not in histname
                and histname in events.Jet.fields
            ):
                h.fill(
                    flatten(genflavor),
                    flatten(sjets[histname]),
                    weight=flatten(
                        ak.broadcast_arrays(weights.weight(), sjets["pt"])[0]
                    ),
                )
            elif (
                "PFCands" in events.fields
                and "PFCands" in histname
                and histname.split("_")[1] in events.PFCands.fields
            ):
                h.fill(
                    flatten(ak.broadcast_arrays(smflav, spfcands["pt"])[0]),
                    flatten(spfcands[histname.replace("PFCands_", "")]),
                    weight=flatten(
                        ak.broadcast_arrays(weights.weight(), spfcands["pt"])[0]
                    ),
                )
            elif "jet_" in histname and "mu" not in histname:
                h.fill(
                    flatten(genflavor),
                    flatten(sjets[histname.replace("jet_", "")]),
                    weight=flatten(
                        ak.broadcast_arrays(weights.weight(), sjets["pt"])[0]
                    ),
                )
            elif "hl_" in histname and histname.replace("hl_", "") in shmu.fields:
                h.fill(
                    flatten(shmu[histname.replace("hl_", "")]),
                    weight=weights.weight(),
                )
            elif (
                "soft_l" in histname and histname.replace("soft_l_", "") in ssmu.fields
            ):
                h.fill(
                    smflav,
                    flatten(ssmu[histname.replace("soft_l_", "")]),
                    weight=weights.weight(),
                )
            elif "mujet_" in histname:
                h.fill(
                    smflav,
                    flatten(smuon_jet[histname.replace("mujet_", "")]),
                    weight=weights.weight(),
                )
            elif ["z_mass", "z_pt", "z_eta", "z_phi"] == histname:
                h.fill(
                    flatten(sz[histname.replace("z_", "")]),
                    weight=weights.weight(),
                )
            elif (
                "btagDeep" in histname
                and "0" in histname
                and histname.replace("_0", "") in events.Jet.fields
            ):
                h.fill(
                    flav=smflav,
                    syst="noSF",
                    discr=np.where(
                        smuon_jet[histname.replace("_0", "")] < 0,
                        -0.2,
                        smuon_jet[histname.replace("_0", "")],
                    ),
                    weight=weights.weight(),
                )
                if (
                    not isRealData
                    and self.isCorr
                    and "btag" in self.SF_map.keys()
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
                            weight=weights.weight()
                            * disc_list[histname.replace("_0", "")][0][syst],
                        )
            elif (
                "btagDeep" in histname and "1" in histname and all(i > 1 for i in njet)
            ) and histname.replace("_1", "") in events.Jet.fields:
                sljets = sjets[:, 1]
                h.fill(
                    flav=genflavor[:, 1],
                    syst="noSF",
                    discr=np.where(
                        sljets[histname.replace("_1", "")] < 0,
                        -0.2,
                        sljets[histname.replace("_1", "")],
                    ),
                    weight=weights.weight(),
                )
                if (
                    not isRealData
                    and self.isCorr
                    and "btag" in self.SF_map.keys()
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
                            weight=weights.weight()
                            * disc_list[histname.replace("_1", "")][1][syst],
                        )
        output["njet"].fill(njet, weight=weights.weight())
        output["nmujet"].fill(nmujet, weight=weights.weight())
        output["nsoftmu"].fill(nsoftmu, weight=weights.weight())
        output["hl_ptratio"].fill(
            flav=genflavor[:, 0],
            ratio=shmu.pt / sjets[:, 0].pt,
            weight=weights.weight(),
        )
        output["soft_l_ptratio"].fill(
            flav=smflav,
            ratio=ssmu.pt / smuon_jet.pt,
            weight=weights.weight(),
        )
        output["dr_lmujetsmu"].fill(
            flav=smflav,
            dr=smuon_jet.delta_r(ssmu),
            weight=weights.weight(),
        )
        output["dr_lmujethmu"].fill(
            flav=smflav,
            dr=smuon_jet.delta_r(shmu),
            weight=weights.weight(),
        )
        output["dr_lmusmu"].fill(dr=shmu.delta_r(ssmu), weight=weights.weight())
        output["z_pt"].fill(flatten(sz.pt), weight=weights.weight())
        output["z_eta"].fill(flatten(sz.eta), weight=weights.weight())
        output["z_phi"].fill(flatten(sz.phi), weight=weights.weight())
        output["z_mass"].fill(flatten(sz.mass), weight=weights.weight())
        output["w_pt"].fill(flatten(sw.pt), weight=weights.weight())
        output["w_eta"].fill(flatten(sw.eta), weight=weights.weight())
        output["w_phi"].fill(flatten(sw.phi), weight=weights.weight())
        output["w_mass"].fill(flatten(sw.mass), weight=weights.weight())
        output["MET_pt"].fill(flatten(smet.pt), weight=weights.weight())
        output["MET_phi"].fill(flatten(smet.phi), weight=weights.weight())
        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator
