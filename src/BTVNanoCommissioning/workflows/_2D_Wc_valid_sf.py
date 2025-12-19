import os
import collections, awkward as ak, numpy as np
import uproot
from coffea import processor
from coffea.analysis_tools import Weights

from BTVNanoCommissioning.utils.correction import (
    load_lumi,
    load_SF,
    weight_manager,
    common_shifts,
)
from BTVNanoCommissioning.helpers.func import update, dump_lumi, PFCand_link, flatten
from BTVNanoCommissioning.helpers.update_branch import missing_branch
from BTVNanoCommissioning.utils.histogramming.histogrammer import (
    histogrammer,
    histo_writer,
)
from BTVNanoCommissioning.utils.array_writer import array_writer
from BTVNanoCommissioning.utils.selection import (
    HLT_helper,
    jet_id,
    mu_idiso,
    ele_mvatightid,
    MET_filters,
    softmu_mask,
    btag_wp,
    btag_wp_dict,
    calculate_new_discriminators,
    get_wp_2D,
)


class NanoProcessor(processor.ProcessorABC):
    def __init__(
        self,
        year="2022",
        campaign="Summer22Run3",
        name="",
        isSyst=False,
        isArray=False,
        noHist=False,
        chunksize=75000,
        selectionModifier="WcM",
    ):
        self._year = year
        self._campaign = campaign
        self.name = name
        self.isSyst = isSyst
        self.isArray = isArray
        self.noHist = noHist
        self.lumiMask = load_lumi(self._campaign)
        self.chunksize = chunksize
        ## Load corrections
        self.SF_map = load_SF(self._year, self._campaign)
        self.selMod = selectionModifier

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        events = missing_branch(events)
        vetoed_events, shifts = common_shifts(self, events)

        return processor.accumulate(
            self.process_shift(update(vetoed_events, collections), name)
            for collections, name in shifts
        )

    def process_shift(self, events, shift_name):
        dataset = events.metadata["dataset"]
        isRealData = not hasattr(events, "genWeight")

        isMu = False
        isEle = False
        ### selections from Spandan
        if "WcM" in self.selMod or "semittM" in self.selMod:
            triggers = ["IsoMu27", "IsoMu24"]
            isMu = True
            dxySigcut = 1.0
            muNeEmSum = 0.7
            muonpTratioCut = 0.4
            ### remove muNeEmSum for cutbased
            if "cutbased" in self.selMod:
                muNeEmSum = 1.0
            ### remove muNeEmSum and loosen muon EF for noMuVeto
            if "noMuVeto" in self.selMod:
                muNeEmSum = 1.0
                muonpTratioCut = 0.6
            isolepdz, isolepdxy, isolepsip3d = 0.01, 0.002, 2
        elif "WcE" in self.selMod or "semittE" in self.selMod:
            triggers = ["Ele32_WPTight_Gsf_L1DoubleEG"]
            isEle = True
            dxySigcut = 0.0
            muNeEmSum = 1.0
            muonpTratioCut = 0.6  # 0.6
            isolepdz, isolepdxy, isolepsip3d = 0.02, 0.01, 2.5
        else:
            raise ValueError(self.selMod, "is not a valid selection modifier.")

        #histoname = {
        #    "WcM": "ctag_Wc_sf",
        #    "WcE": "ectag_Wc_sf",
        #    "WcM_2D" : "ctag_Wc_sf_2D",
        #    "WcE_2D" : "ectag_Wc_sf_2D"
        #}
        output = {}
        if not self.noHist:
            output = histogrammer(
                events.Jet.fields,
                obj_list=["hl", "soft_l", "MET", "dilep", "mujet"],
                hist_collections=["common", "fourvec", "Wc"],
                include_nmujet=True,
                include_nsoftmu=True,
                include_osss=True,
                cutbased=("cutbased" in self.selMod),
                year=self._year,
                campaign=self._campaign,
                include_discriminators_2D=True if "2D" in self.selMod else False,
            )

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
        # only dump for nominal case
        if shift_name is None:
            output = dump_lumi(events[req_lumi], output)

        ## HLT
        req_trig = HLT_helper(events, triggers)

        ## Lepton cuts
        if isMu:
            # muon twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
            iso_lep = events.Muon[
                (events.Muon.pt > 30) & mu_idiso(events, self._campaign)
            ]
        elif isEle:
            iso_lep = events.Electron[
                (events.Electron.pt > 34) & ele_mvatightid(events, self._campaign)
            ]
        req_lep = ak.count(iso_lep.pt, axis=1) == 1

        iso_lep = ak.pad_none(iso_lep, 1, axis=1)
        iso_lep = iso_lep[:, 0]
        #sub_lep = events.Muon[
        #        (events.Muon.pt < 30) &  (abs(events.Muon.eta) < 2.4) & (events.Muon.tightId > 0.5)
        #    ]
        #sub_lep = ak.pad_none(sub_lep, 1, axis=1)
        #sub_lep = sub_lep[:, 0]

        """
        if isMu:
            iso_lepindx = ak.mask(
                ak.local_index(events.Muon.pt),
                ((events.Muon.pt > 30) & mu_idiso(events, self._campaign)) == 1,
            )
        elif isEle:
            iso_lepindx = ak.mask(
                ak.local_index(events.Electron.pt),
                ((events.Electron.pt > 34) & ele_mvatightid(events, self._campaign))
                == 1,
            )
        iso_lepindx = ak.pad_none(iso_lepindx, 1)
        iso_lepindx = iso_lepindx[:, 0]
        """

        ## Jet cuts
        jet_sel = ak.fill_none(
            jet_id(events, self._campaign, min_pt=25)
            & (ak.all(events.Jet.metric_table(iso_lep) > 0.5, axis=2)),
            False,
            axis=-1,
        )
        event_jet = events.Jet[jet_sel]
        nseljet = ak.count(event_jet.pt, axis=1)
        #req_base_jets = (nseljet >= 1) & (nseljet <= 3)
        req_jets = (nseljet == 1)

        ## Soft Muon cuts
        soft_muon = events.Muon[
                (events.Muon.pt < 25)
                & (abs(events.Muon.eta) < 2.4)
                & (events.Muon.tightId > 0.5)
                & (events.Muon.pfRelIso04_all > 0.25)
        ]
        soft_muon_tight = soft_muon[
            (abs(soft_muon.dxy / soft_muon.dxyErr) > dxySigcut)
            & (soft_muon.jetIdx != -1)
        ]
        req_softmu = ak.count(soft_muon.pt, axis=1) >= 1
        req_softmu_tight = ak.count(soft_muon_tight.pt, axis=1) >= 1
        req_softmu_pt = ak.count(soft_muon[soft_muon.pt > 5.0].pt, axis=1)  >= 1
        mujetsel = ak.fill_none(
            (
                (ak.all(event_jet.metric_table(soft_muon) <= 0.4, axis=2))
                & ((event_jet.muonIdx1 != -1) | (event_jet.muonIdx2 != -1))
                & ((event_jet.muEF + event_jet.neEmEF) < muNeEmSum)
                & (event_jet.pt > 20)
                & ((event_jet.pt / event_jet.E) > 0.03)
            ),
            False,
            axis=-1,
        )
        mujetsel2 = ak.fill_none(
            (
                ((events.Jet.muEF + events.Jet.neEmEF) < muNeEmSum)
                & (
                    ak.all(
                        events.Jet.metric_table(soft_muon) <= 0.4,
                        axis=2,
                        mask_identity=True,
                    )
                )
                & ((events.Jet.muonIdx1 != -1) | (events.Jet.muonIdx2 != -1))
            ),
            False,
            axis=-1,
        )
        soft_muon = ak.pad_none(soft_muon, 1, axis=1)
        soft_muon["dxySig"] = soft_muon.dxy / soft_muon.dxyErr

        ## Muon-jet cuts
        event_jet["isMuonJet"] = mujetsel
        mu_jet = event_jet[mujetsel]
        otherjets = event_jet[~mujetsel]
        req_mujet = ak.num(mu_jet.pt, axis=1) >= 1
        mu_jet = ak.pad_none(mu_jet, 1, axis=1)

        # jet energy fraction cuts
        muNeEmSum_sel = ak.fill_none(
            (event_jet.muEF + event_jet.neEmEF) < muNeEmSum,
            False,
            axis=-1,
        )
        req_muNeEmSum = ak.num(event_jet[muNeEmSum_sel].pt, axis=1) == 1
        muEF_sel = ak.fill_none(
            (event_jet.muEF < 0.5),
            False,
            axis=-1,
        )
        req_muEF = ak.num(event_jet[muEF_sel].pt, axis=1) == 1


        ## store jet index for PFCands, create mask on the jet index
        jet_selpf = (jet_sel) & (mujetsel2)
        jetindx = ak.mask(ak.local_index(events.Jet.pt), jet_selpf == True)
        jetindx = ak.pad_none(jetindx, 1)
        jetindx = jetindx[:, 0]

        # Other cuts
        # Ratio between soft muon pt and muon jet pt
        #req_pTratio = (soft_muon[:, 0].pt / mu_jet[:, 0].pt) < muonpTratioCut
        idx = np.where(iso_lep.jetIdx == -1, 0, iso_lep.jetIdx)
        ## Additional cut to reject QCD events,used in BTV-20-001
        # req_QCDveto = (
        #     (iso_lep.pfRelIso04_all < 0.05)
        # & (abs(iso_lep.dz) < isolepdz)
        # & (abs(iso_lep.dxy) < isolepdxy)
        # & (iso_lep.sip3d < isolepsip3d)
        # & (
        #     iso_lep.pt
        #     / ak.firsts(
        #         events.Jet[
        #             (events.Jet.muonIdx1 == iso_lepindx)
        #             | ((events.Jet.muonIdx2 == iso_lepindx))
        #         ].pt
        #     )
        #     > 0.75
        # )
        # )

        ## Dilepton veto
        dilep_mu = events.Muon[(events.Muon.pt > 12) & mu_idiso(events, self._campaign)]
        dilep_ele = events.Electron[
            (events.Electron.pt > 15) & ele_mvatightid(events, self._campaign)
        ]
        req_dilepveto = (
            ak.count(dilep_mu.pt, axis=1) + ak.count(dilep_ele.pt, axis=1) != 2
        )
    
        dilep_mass = iso_lep + soft_muon[:, 0]
        if isMu and "noMuVeto" not in self.selMod:
            req_dilepmass = (dilep_mass.mass > 12.0) & (
                (dilep_mass.mass < 80) | (dilep_mass.mass > 100)
            )
        else:
            req_dilepmass = iso_lep.pt > 0

        ## MET and Transverse W mass
        iso_lep_trans = ak.zip(
            {
                "pt": iso_lep.pt,
                "eta": ak.zeros_like(iso_lep.pt),
                "phi": iso_lep.phi,
                "mass": iso_lep.mass,
            },
            with_name="PtEtaPhiMLorentzVector",
        )

        MET = ak.zip(
            {
                "pt": events.PuppiMET.pt,
                "eta": ak.zeros_like(events.PuppiMET.pt),
                "phi": events.PuppiMET.phi,
                "mass": ak.zeros_like(events.PuppiMET.pt),
            },
            with_name="PtEtaPhiMLorentzVector",
        )

        req_metfilter = MET_filters(events, self._campaign)
        event_jet_0 = ak.pad_none(event_jet, 1, axis=1)
        event_jet_0 = event_jet_0[:, 0]
        req_met_pt = MET.pt > 30.0
        jetmet_dphi = abs(event_jet_0.delta_phi(MET))
        req_jetmet_dphi = jetmet_dphi > 1.0
        metTrkmet_dphi = abs(MET.delta_phi(events.TrkMET))
        req_metTrkmet_dphi = (metTrkmet_dphi < 1.0) #& (metTrkmet_dphi >= 0.5)

        # W cuts
        wmasscut = 40 #before 55
        wmasscut_tight = 55
        wptcut = 30
        if "semitt" in self.selMod:
            wmasscut = 0
        Wcand = MET + iso_lep_trans  # transverse mass
        Wmass = Wcand.mass
        Wpt = Wcand.pt
        req_WpT = Wpt > wptcut
        req_mtw = Wmass > wmasscut
        req_mtw_max120 =  (Wmass < 120)
        req_mtw_min55 = (Wmass > wmasscut_tight)
        # req for pt ratio of leading jet and W pt
        req_pTratio = (soft_muon[:, 0].pt / event_jet_0.pt) < muonpTratioCut
        jetw_ptratio = event_jet_0.pt / Wpt
        req_jWpTratio = (jetw_ptratio < 2) & (jetw_ptratio > 0.5)
        #mask on distance between W_phi and jet_phi abs(deltaPhi(ak4jets[0].phi, Vboson.phi) > 2 (if tight sel)
        jetw_dphi = abs(event_jet_0.delta_phi(Wcand))
        req_dphi_Wjet = jetw_dphi > 2
        # delta Phi between leading jet and lepton
        jetl_dphi = abs(event_jet_0.delta_phi(iso_lep))
        req_jetl_dphi = jetl_dphi < 2.

        # ==This is the manual calculation for transverse mass==
        """
        dphi = iso_lep.phi-events.PuppiMET.phi
        dphi = np.where(dphi<np.pi,dphi+2*np.pi,dphi)
        dphi = np.where(dphi>np.pi,dphi-2*np.pi,dphi)
        trans = np.sqrt(2*iso_lep.pt*events.PuppiMET.pt*(1-np.cos(dphi)))
        """

        event_level = (
            req_trig
            & req_lumi
            & req_lep
            & req_jets
            & req_softmu
            & req_dilepmass
            & req_mtw
            & req_metfilter
            & req_met_pt
            & req_jetmet_dphi
            & req_metTrkmet_dphi
            & req_WpT
            & req_jWpTratio
            & req_dphi_Wjet
            & req_jetl_dphi
            & req_muNeEmSum
            & req_muEF
            & req_softmu_pt
            & req_mtw_max120
        )

        event_level = ak.fill_none(event_level, False)
        if len(events[event_level]) == 0:
            if self.isArray:
                array_writer(
                    self,
                    events[event_level],
                    events,
                    None,
                    ["nominal"],
                    dataset,
                    isRealData,
                    empty=True,
                )
            return {dataset: output}

        ####################
        # Selected objects #
        ####################

        shmu = iso_lep[event_level]
        wm = Wmass[event_level]
        wp = Wpt[event_level]
        sjets = event_jet[event_level]
        ssmu = soft_muon[event_level]
        smet = MET[event_level]
        smuon_jet = mu_jet[event_level]
        sotherjets = otherjets[event_level]
        sdilep = dilep_mass[event_level]
        nsoftmu = ak.count(ssmu.pt, axis=1)
        #nmujet = ak.count(smuon_jet.pt, axis=1)
        #smuon_jet = smuon_jet[ak.num(smuon_jet, axis=1) > 0]
        #smuon_jet = ak.pad_none(smuon_jet, 1, axis=1)  # Ensure at least one entry
        #`smuon_jet = smuon_jet[~ak.is_none(smuon_jet, axis=1)]  # Remove None entrie        smuon_jet = smuon_jet[:, 0]
        smuon_jet = sjets[:, 0]
        ssmu = ssmu[:, 0]
        sz = shmu + ssmu
        sw = shmu + smet

        osss = shmu.charge * ssmu.charge * -1
        ossswrite = shmu.charge * ssmu.charge * -1
  
        njet = ak.count(sjets.pt, axis=1)
        # Find the PFCands associate with selected jets. Search from jetindex->JetPFCands->PFCand
        if "PFCands" in events.fields:
            spfcands = PFCand_link(events, event_level, jetindx)

        # Keep the structure of events and pruned the object size
        pruned_ev = events[event_level]
        pruned_ev["SelJet"] = sjets
        #lead_jet = sjets[:, 0]
        #pruned_ev["LeadingJet"] = lead_jet
        if self.selMod.endswith("M") or self.selMod == "WcM_2D":
            pruned_ev["SelMuon"] = shmu
        else:
            pruned_ev["SelElectron"] = shmu
        pruned_ev["MuonJet"] = smuon_jet #leading selected jet fullfilling all event-level cuts
        pruned_ev["SoftMuon"] = ssmu
        pruned_ev["OtherJets"] = sotherjets
        pruned_ev["MET_pt"] = smet.pt
        if "Wc" in self.selMod:
            pruned_ev["osss"] = osss
        else:
            pruned_ev["osss"] = ossswrite
        pruned_ev["njet"] = njet
        pruned_ev["w_mass"] = wm
        pruned_ev["w_pt"] = wp
        pruned_ev["dilep_mass"] = sdilep.mass
        pruned_ev["dilep_pt"] = sdilep.pt
        pruned_ev["jetw_dphi"] = jetw_dphi[event_level]
        pruned_ev["jetmet_dphi"] = jetmet_dphi[event_level]
        pruned_ev["metTrkmet_dphi"] = metTrkmet_dphi[event_level]
        pruned_ev["jetl_dphi"] = jetl_dphi[event_level]
        if "PFCands" in events.fields:
            pruned_ev.PFCands = spfcands
        # Add custom variables

        pruned_ev["dr_mujet_softmu"] = ssmu.delta_r(smuon_jet)
        pruned_ev["dr_mujet_lep1"] = shmu.delta_r(smuon_jet)
        pruned_ev["dr_lep1_softmu"] = shmu.delta_r(ssmu)
        pruned_ev["soft_l_ptratio"] = ssmu.pt / smuon_jet.pt
        pruned_ev["l1_ptratio"] = shmu.pt / smuon_jet.pt
        pruned_ev["MuonJet_beta"] = smuon_jet.pt / smuon_jet.E
        pruned_ev["MuonJet_muneuEF"] = smuon_jet.muEF + smuon_jet.neEmEF
        pruned_ev["jetw_ptratio"] = jetw_ptratio[event_level]

        if "hadronFlavour" in pruned_ev.SelJet.fields:
            isRealData = False
            genflavor = ak.values_astype(
                pruned_ev.SelJet.hadronFlavour
                + 1
                * (
                    (pruned_ev.SelJet.partonFlavour == 0)
                    & (pruned_ev.SelJet.hadronFlavour == 0)
                ),
                int,
            )
            if "MuonJet" in pruned_ev.fields:
                smflav = ak.values_astype(
                    1
                    * (
                        (pruned_ev.MuonJet.partonFlavour == 0)
                        & (pruned_ev.MuonJet.hadronFlavour == 0)
                    )
                    + pruned_ev.MuonJet.hadronFlavour,
                    int,
                )
        else:
            isRealData = True
            genflavor = ak.ones_like(pruned_ev.SelJet.pt, dtype=int)
            if "MuonJet" in pruned_ev.fields:
                smflav = ak.ones_like(pruned_ev.MuonJet.pt, dtype=int)

        if "2D" in self.selMod:
            nj=1
            for i in range(nj):
                btagUParTAK4HFvLF, btagUParTAK4BvC = calculate_new_discriminators(pruned_ev.MuonJet)
                wp2D = ak.Array([get_wp_2D(btagUParTAK4HFvLF[i], btagUParTAK4BvC[i], self._year, self._campaign, "UParTAK4") for i in range(len(btagUParTAK4HFvLF))])
                pruned_ev[f"btagUParTAK4HFvLF_{i}"] = btagUParTAK4HFvLF
                pruned_ev[f"btagUParTAK4BvC_{i}"] = btagUParTAK4BvC
                pruned_ev[f"btagUParTAK4HFvLFt_{i}"] = ak.Array(np.where(btagUParTAK4HFvLF > 0.0, 1.0 - (1.0 - btagUParTAK4HFvLF)**0.5, -1.0))
                pruned_ev[f"btagUParTAK4BvCt_{i}"] = ak.Array(np.where(btagUParTAK4BvC > 0.0, 1.0 - (1.0 - btagUParTAK4BvC)**0.5, -1.0))
                pruned_ev[f"btagUParTAK42D_{i}"] = wp2D
                jet_pt_bins = btag_wp_dict[self._year + "_" + self._campaign]["UParTAK4"]["2D"]["jet_pt_bins"]
                for jet_pt_bin in jet_pt_bins:
                    pruned_ev[f"btagUParTAK42D_pt{jet_pt_bin[0]}to{jet_pt_bin[1]}_{i}"] = [wp2D[ijet] if pt is not None and jet_pt_bin[0] < pt and pt < jet_pt_bin[1] else None for ijet, pt in enumerate(pruned_ev.MuonJet.pt)]
        
        ### Store masks for additional cuts
        #pruned_ev["mask_njets"]     = req_jets[event_level]
        pruned_ev["mask_mujet"]     = req_mujet[event_level]
        pruned_ev["mask_muNeEmSum"] = req_muNeEmSum[event_level]
        pruned_ev["mask_muEF"]     = req_muEF[event_level]
        pruned_ev["mask_jetl_dphi"] = req_jetl_dphi[event_level]
        pruned_ev["mask_dilepveto"] = req_dilepveto[event_level]
        pruned_ev["mask_smj_ptratio"]   = req_pTratio[event_level]
        pruned_ev["mask_met_pt"]    = req_met_pt[event_level]
        pruned_ev["mask_jetmet_dphi"]  = req_jetmet_dphi[event_level]
        pruned_ev["mask_metTrkmet_dphi"]  = req_metTrkmet_dphi[event_level]
        pruned_ev["req_mtw_max120"] = req_mtw_max120[event_level]
        pruned_ev["req_mtw_min55"]  = req_mtw_min55[event_level]
        pruned_ev["mask_wpt"]       = req_WpT[event_level]
        pruned_ev["mask_softmu_tight"] = req_softmu_tight[event_level]
        pruned_ev["mask_softmu_pt"] = req_softmu_pt[event_level]
        pruned_ev["mask_jWpTratio"] = req_jWpTratio[event_level]
        pruned_ev["mask_dphi_Wjet"] = req_dphi_Wjet[event_level]


        ####################
        # Weight & Geninfo #
        ####################

        weights = weight_manager(pruned_ev, self.SF_map, self.isSyst)

        if shift_name is None:
            systematics = ["nominal"] + list(weights.variations)
        else:
            systematics = [shift_name]
        exclude_btv = [
            "DeepCSVC",
            "DeepCSVB",
            "DeepJetB",
            "DeepJetC",
        ]  # exclude b-tag SFs for btag inputs
        ####################
        #  Fill histogram  #
        ####################
        for syst in systematics:
            if self.isSyst == False and syst != "nominal":
                break
            if self.noHist:
                break
            weight = (
                weights.weight()
                if syst == "nominal" or syst == shift_name
                else weights.weight(modifier=syst)
            )
            for histname, h in output.items():
                if (
                    "Deep" in histname
                    and "btag" not in histname
                    and histname in events.Jet.fields
                ):
                    h.fill(
                        syst,
                        flatten(genflavor),
                        flatten(ak.broadcast_arrays(osss, sjets["pt"])[0]),
                        flatten(sjets[histname]),
                        weight=flatten(
                            ak.broadcast_arrays(
                                weights.partial_weight(exclude=exclude_btv), sjets["pt"]
                            )[0]
                        ),
                    )
                elif (
                    "PFCands" in events.fields
                    and "PFCands" in histname
                    and histname.split("_")[1] in events.PFCands.fields
                ):
                    h.fill(
                        syst,
                        flatten(ak.broadcast_arrays(smflav, spfcands["pt"])[0]),
                        flatten(ak.broadcast_arrays(osss, spfcands["pt"])[0]),
                        flatten(spfcands[histname.replace("PFCands_", "")]),
                        weight=flatten(
                            ak.broadcast_arrays(
                                weights.partial_weight(exclude=exclude_btv),
                                spfcands["pt"],
                            )[0]
                        ),
                    )
                elif "jet_" in histname and "mu" not in histname:
                    h.fill(
                        syst,
                        flatten(genflavor),
                        flatten(ak.broadcast_arrays(osss, sjets["pt"])[0]),
                        flatten(sjets[histname.replace("jet_", "")]),
                        weight=flatten(ak.broadcast_arrays(weight, sjets["pt"])[0]),
                    )
                elif "hl_" in histname and histname.replace("hl_", "") in shmu.fields:
                    h.fill(
                        syst,
                        osss,
                        flatten(shmu[histname.replace("hl_", "")]),
                        weight=weight,
                    )
                elif (
                    "soft_l" in histname
                    and histname.replace("soft_l_", "") in ssmu.fields
                ):
                    h.fill(
                        syst,
                        smflav,
                        osss,
                        flatten(ssmu[histname.replace("soft_l_", "")]),
                        weight=weight,
                    )
                elif (
                    "mujet_" in histname
                    and histname.replace("mujet_", "") in smuon_jet.fields
                ):
                    h.fill(
                        syst,
                        smflav,
                        osss,
                        flatten(smuon_jet[histname.replace("mujet_", "")]),
                        weight=weight,
                    )
                elif "btag" in histname and "Trans" not in histname:
                    if (
                        "BvC" not in histname
                        and "HFvLF" not in histname
                        and "2D" not in histname
                    ):
                        for i in range(2):
                            if not histname.endswith(str(i)) or histname.replace(f"_{i}", "") not in smuon_jet.fields:
                                continue
                            h.fill(
                                syst="noSF",
                                flav=smflav,
                                osss=osss,
                                discr=np.where(
                                    smuon_jet[histname.replace(f"_{i}", "")] < 0,
                                    -0.2,
                                    smuon_jet[histname.replace(f"_{i}", "")],
                                ),
                                weight=weights.partial_weight(exclude=exclude_btv),
                            )
                            if not isRealData and "btag" in self.SF_map.keys():
                                h.fill(
                                    syst=syst,
                                    flav=smflav,
                                    osss=osss,
                                    discr=np.where(
                                        smuon_jet[histname.replace(f"_{i}", "")] < 0,
                                        -0.2,
                                        smuon_jet[histname.replace(f"_{i}", "")],
                                    ),
                                    weight=weight,
                                )
                    else:
                        discr_to_plot = pruned_ev[histname]
                        flav_to_plot = smflav
                        osss_to_plot = osss
                        weight_to_plot = weights.partial_weight(exclude=exclude_btv)
                        discr_mask = [False if x is None else True for x in discr_to_plot]
                        if ak.any(np.invert(discr_mask)):
                            flav_to_plot = flav_to_plot[discr_mask]
                            osss_to_plot = osss_to_plot[discr_mask]
                            discr_to_plot = discr_to_plot[discr_mask]
                            weight_to_plot = weight[discr_mask]
                        h.fill(
                            syst=syst,
                            flav=flav_to_plot,
                            osss=osss_to_plot,
                            discr=discr_to_plot,
                            weight=weight_to_plot
                        )
                #elif "btag" in histname and "Trans" in histname:
                #    if histname not in smuon_jet:
                #        continue
                #    for i in range(2):
                #        histname = histname.replace("Trans", "").replace(f"_{i}", "")
                #        h.fill(
                #            syst="noSF",
                #            flav=smflav,
                #            osss=osss,
                #            discr=1.0 / np.tanh(smuon_jet[histname]),
                #            weight=weights.partial_weight(exclude=exclude_btv),
                #        )

            output["njet"].fill(syst, osss, njet, weight=weight)
            #output["nmujet"].fill(syst, osss, nmujet, weight=weight)
            output["nsoftmu"].fill(syst, osss, nsoftmu, weight=weight)
            output["softlpt"].fill(syst, smflav, osss, ssmu.pt, weight=weight)
            output["hl_ptratio"].fill(
                syst,
                genflavor[:, 0],
                osss=osss,
                ratio=shmu.pt / sjets[:, 0].pt,
                weight=weight,
            )
            output["soft_l_ptratio"].fill(
                syst,
                flav=smflav,
                osss=osss,
                ratio=ssmu.pt / smuon_jet.pt,
                weight=weight,
            )
            output["dr_lmujetsmu"].fill(
                syst,
                flav=smflav,
                osss=osss,
                dr=smuon_jet.delta_r(ssmu),
                weight=weight,
            )
            output["dr_lmujethmu"].fill(
                syst,
                flav=smflav,
                osss=osss,
                dr=smuon_jet.delta_r(shmu),
                weight=weight,
            )
            output["dr_hmusmu"].fill(
                syst,
                osss=osss,
                dr=shmu.delta_r(ssmu),
                weight=weight,
            )
            output["mujet_muneuEF"].fill(
                syst,
                smflav,
                osss=osss,
                muneuEF=flatten(pruned_ev["MuonJet_muneuEF"]),
                weight=weight,
            )
            output["dilep_pt"].fill(syst, osss, flatten(sz.pt), weight=weight)
            output["dilep_eta"].fill(syst, osss, flatten(sz.eta), weight=weight)
            output["dilep_phi"].fill(syst, osss, flatten(sz.phi), weight=weight)
            output["dilep_mass"].fill(syst, osss, flatten(sz.mass), weight=weight)
            output["w_pt"].fill(syst, osss, flatten(sw.pt), weight=weight)
            output["w_phi"].fill(syst, osss, flatten(sw.phi), weight=weight)
            output["w_mass"].fill(syst, osss, flatten(sw.mass), weight=weight)
            output["jetw_dphi"].fill(syst, smflav, osss, flatten(jetw_dphi[event_level]), weight=weight)
            output["jetw_ptratio"].fill(syst, smflav, osss, flatten(jetw_ptratio[event_level]), weight=weight)
            output["MET_pt"].fill(syst, osss, flatten(smet.pt), weight=weight)
            output["MET_phi"].fill(syst, osss, flatten(smet.phi), weight=weight)
            output["jetmet_dphi"].fill(syst, smflav, osss, flatten(jetmet_dphi[event_level]), weight=weight)
            output["metTrkmet_dphi"].fill(syst, osss, flatten(metTrkmet_dphi[event_level]), weight=weight)
            output["jetl_dphi"].fill(syst, smflav, osss, flatten(jetl_dphi[event_level]), weight=weight)
        #######################
        #  Create root files  #
        #######################
        if self.isArray:
            array_writer(
                self, pruned_ev, events, weights, systematics, dataset, isRealData
            )

        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator