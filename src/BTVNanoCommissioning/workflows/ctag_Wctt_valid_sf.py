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
    histo_writter,
)
from BTVNanoCommissioning.utils.array_writer import array_writer
from BTVNanoCommissioning.utils.selection import (
    HLT_helper,
    jet_id,
    mu_idiso,
    ele_mvatightid,
    softmu_mask,
    btag_wp,
    btag_wp_dict,
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
        jet_sel = ak.fill_none(
            jet_id(events, self._campaign)
            & (ak.all(events.Jet.metric_table(iso_lep) > 0.5, axis=2)),
            False,
            axis=-1,
        )
        iso_lep = ak.pad_none(iso_lep, 1, axis=1)
        iso_lep = iso_lep[:, 0]
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

        ## Jet cuts
        if "DeepJet_nsv" in events.Jet.fields:
            jet_sel = jet_sel & (events.Jet.DeepJet_nsv > 0)
        event_jet = events.Jet[jet_sel]
        nseljet = ak.count(event_jet.pt, axis=1)
        if "Wc" in self.selMod:
            req_jets = (nseljet >= 1) & (nseljet <= 3)
        else:
            req_jets = nseljet >= 4

        ## Soft Muon cuts
        soft_muon = events.Muon[
            softmu_mask(events, self._campaign)
            & (abs(events.Muon.dxy / events.Muon.dxyErr) > dxySigcut)
        ]
        req_softmu = ak.count(soft_muon.pt, axis=1) >= 1
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

        ## store jet index for PFCands, create mask on the jet index
        jet_selpf = (jet_sel) & (mujetsel2)
        if "DeepJet_nsv" in events.Jet.fields:
            jet_selpf = jet_selpf & (events.Jet.DeepJet_nsv > 0)
        jetindx = ak.mask(ak.local_index(events.Jet.pt), jet_selpf == True)
        jetindx = ak.pad_none(jetindx, 1)
        jetindx = jetindx[:, 0]

        # Other cuts
        req_pTratio = (soft_muon[:, 0].pt / mu_jet[:, 0].pt) < muonpTratioCut
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
        iso_lep_trans = ak.zip(
            {
                "pt": iso_lep.pt,
                "eta": ak.zeros_like(iso_lep.pt),
                "phi": iso_lep.phi,
                "mass": iso_lep.mass,
            },
            with_name="PtEtaPhiMLorentzVector",
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

        wmasscut = 55
        if "semitt" in self.selMod:
            wmasscut = 0
        Wcand = MET + iso_lep_trans  # transverse mass
        Wmass = Wcand.mass
        Wpt = Wcand.pt
        req_mtw = Wmass > wmasscut

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
            & req_mujet
            & req_mtw
            & req_dilepveto
            & req_pTratio
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
        nmujet = ak.count(smuon_jet.pt, axis=1)
        smuon_jet = smuon_jet[:, 0]
        ssmu = ssmu[:, 0]
        sz = shmu + ssmu
        sw = shmu + smet

        osss = 1
        ossswrite = shmu.charge * ssmu.charge * -1
        smuon_jet_passc = {}
        c_algos = []
        c_wps = []
        if "cutbased_Wc" in self.selMod:
            osss = shmu.charge * ssmu.charge * -1
            c_algos = btag_wp_dict[self._year + "_" + self._campaign].keys()
            for c_algo in c_algos:
                smuon_jet_passc[c_algo] = {}
                c_wps = btag_wp_dict[self._year + "_" + self._campaign][c_algo][
                    "c"
                ].keys()
                for c_wp in c_wps:
                    if not "No" in c_wp:
                        smuon_jet_passc[c_algo][c_wp] = btag_wp(
                            smuon_jet, self._year, self._campaign, c_algo, "c", c_wp
                        )
        njet = ak.count(sjets.pt, axis=1)
        # Find the PFCands associate with selected jets. Search from jetindex->JetPFCands->PFCand
        if "PFCands" in events.fields:
            spfcands = PFCand_link(events, event_level, jetindx)

        # Keep the structure of events and pruned the object size
        pruned_ev = events[event_level]
        pruned_ev["SelJet"] = sjets
        if self.selMod.endswith("M"):
            pruned_ev["SelMuon"] = shmu
        else:
            pruned_ev["SelElectron"] = shmu
        pruned_ev["MuonJet"] = smuon_jet
        pruned_ev["SoftMuon"] = ssmu
        pruned_ev["OtherJets"] = sotherjets
        if "Wc" in self.selMod:
            pruned_ev["osss"] = osss
        else:
            pruned_ev["osss"] = ossswrite
        pruned_ev["njet"] = njet
        pruned_ev["W_transmass"] = wm
        pruned_ev["W_pt"] = wp
        pruned_ev["dilep_mass"] = sdilep.mass
        pruned_ev["dilep_pt"] = sdilep.pt
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
                    for i in range(2):
                        if (
                            str(i) not in histname
                            or histname.replace(f"_{i}", "") not in events.Jet.fields
                        ):
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
                elif "btag" in histname and "Trans" in histname:
                    if histname not in smuon_jet:
                        continue
                    for i in range(2):
                        histname = histname.replace("Trans", "").replace(f"_{i}", "")
                        h.fill(
                            syst="noSF",
                            flav=smflav,
                            osss=osss,
                            discr=1.0 / np.tanh(smuon_jet[histname]),
                            weight=weights.partial_weight(exclude=exclude_btv),
                        )

            output["njet"].fill(syst, osss, njet, weight=weight)
            output["nmujet"].fill(syst, osss, nmujet, weight=weight)
            output["nsoftmu"].fill(syst, osss, nsoftmu, weight=weight)
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
            output["dilep_pt"].fill(syst, osss, flatten(sz.pt), weight=weight)
            output["dilep_eta"].fill(syst, osss, flatten(sz.eta), weight=weight)
            output["dilep_phi"].fill(syst, osss, flatten(sz.phi), weight=weight)
            output["dilep_mass"].fill(syst, osss, flatten(sz.mass), weight=weight)
            output["w_pt"].fill(syst, osss, flatten(sw.pt), weight=weight)
            output["w_phi"].fill(syst, osss, flatten(sw.phi), weight=weight)
            output["w_mass"].fill(syst, osss, flatten(sw.mass), weight=weight)
            output["MET_pt"].fill(syst, osss, flatten(smet.pt), weight=weight)
            output["MET_phi"].fill(syst, osss, flatten(smet.phi), weight=weight)
            if "cutbased_Wc" in self.selMod:
                for c_algo in c_algos:
                    for c_wp in c_wps:
                        if not "No" in c_wp:
                            output[f"mujet_pt_{c_algo}{c_wp}"].fill(
                                syst,
                                smflav[smuon_jet_passc[c_algo][c_wp]],
                                osss[smuon_jet_passc[c_algo][c_wp]],
                                flatten(smuon_jet[smuon_jet_passc[c_algo][c_wp]].pt),
                                weight=weight[smuon_jet_passc[c_algo][c_wp]],
                            )
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
