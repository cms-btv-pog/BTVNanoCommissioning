import os
import collections, awkward as ak, numpy as np
import uproot
from coffea import processor
from coffea.analysis_tools import Weights

from BTVNanoCommissioning.utils.correction import (
    load_lumi,
    load_SF,
    muSFs,
    eleSFs,
    puwei,
    btagSFs,
    JME_shifts,
    Roccor_shifts,
)
from BTVNanoCommissioning.helpers.func import (
    flatten,
    update,
    dump_lumi,
)
from BTVNanoCommissioning.helpers.update_branch import missing_branch
from BTVNanoCommissioning.utils.histogrammer import histogrammer
from BTVNanoCommissioning.utils.array_writer import array_writer
from BTVNanoCommissioning.utils.selection import (
    jet_id,
    mu_idiso,
    ele_mvatightid,
    softmu_mask,
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
        self.SF_map = load_SF(self._campaign)
        self.selMod = selectionModifier

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        isRealData = not hasattr(events, "genWeight")
        dataset = events.metadata["dataset"]
        events = missing_branch(events)
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

        return processor.accumulate(
            self.process_shift(update(events, collections), name)
            for collections, name in shifts
        )

    def process_shift(self, events, shift_name):
        dataset = events.metadata["dataset"]
        isRealData = not hasattr(events, "genWeight")

        isMu = False
        isEle = False
        if "WcM" in self.selMod or "semittM" in self.selMod:
            triggers = ["IsoMu27", "IsoMu24"]
            isMu = True
            dxySigcut = 1.0
            muNeEmSum = 0.7
            muonpTratioCut = 0.4
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

        histoname = {
            "WcM": "ctag_Wc_sf",
            "WcE": "ectag_Wc_sf",
            "semittM": "ctag_Wc_sf",  # same histogram representation as W+c
            "semittE": "ectag_Wc_sf",  # same histogram representation as W+c
        }
        _hist_event_dict = (
            {"": None} if self.noHist else histogrammer(events, histoname[self.selMod])
        )

        output = {
            "sumw": processor.defaultdict_accumulator(float),
            **_hist_event_dict,
        }

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
        if isMu:
            req_dilepmass = (dilep_mass.mass > 12.0) & (
                (dilep_mass.mass < 80) | (dilep_mass.mass > 100)
            )
        elif isEle:
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
                    "nominal",
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
        if "Wc" in self.selMod:
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
            genflavor = sjets.hadronFlavour + 1 * (
                (sjets.partonFlavour == 0) & (sjets.hadronFlavour == 0)
            )
            smflav = smuon_jet.hadronFlavour + 1 * (
                (smuon_jet.partonFlavour == 0) & (smuon_jet.hadronFlavour == 0)
            )
            if len(self.SF_map.keys()) > 0:
                syst_wei = True if self.isSyst != False else False
                if "PU" in self.SF_map.keys():
                    puwei(
                        events[event_level].Pileup.nTrueInt,
                        self.SF_map,
                        weights,
                        syst_wei,
                    )
                if isMu and "MUO" in self.SF_map.keys():
                    muSFs(shmu, self.SF_map, weights, syst_wei, False)
                if isEle and "EGM" in self.SF_map.keys():
                    eleSFs(shmu, self.SF_map, weights, syst_wei, False)
                if "BTV" in self.SF_map.keys():
                    btagSFs(smuon_jet, self.SF_map, weights, "DeepJetC", syst_wei)
                    btagSFs(smuon_jet, self.SF_map, weights, "DeepJetB", syst_wei)
                    btagSFs(smuon_jet, self.SF_map, weights, "DeepCSVB", syst_wei)
                    btagSFs(smuon_jet, self.SF_map, weights, "DeepCSVC", syst_wei)

        else:
            genflavor = ak.zeros_like(sjets.pt, dtype=int)
            smflav = ak.zeros_like(smuon_jet.pt, dtype=int)

        # Systematics information
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
                elif "mujet_" in histname:
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
            output["dr_lmusmu"].fill(
                syst,
                osss=osss,
                dr=shmu.delta_r(ssmu),
                weight=weight,
            )
            output["z_pt"].fill(syst, osss, flatten(sz.pt), weight=weight)
            output["z_eta"].fill(syst, osss, flatten(sz.eta), weight=weight)
            output["z_phi"].fill(syst, osss, flatten(sz.phi), weight=weight)
            output["z_mass"].fill(syst, osss, flatten(sz.mass), weight=weight)
            output["w_pt"].fill(syst, osss, flatten(sw.pt), weight=weight)
            output["w_eta"].fill(syst, osss, flatten(sw.eta), weight=weight)
            output["w_phi"].fill(syst, osss, flatten(sw.phi), weight=weight)
            output["w_mass"].fill(syst, osss, flatten(sw.mass), weight=weight)
            output["MET_pt"].fill(syst, osss, flatten(smet.pt), weight=weight)
            output["MET_phi"].fill(syst, osss, flatten(smet.phi), weight=weight)
            output["npvs"].fill(
                syst,
                events[event_level].PV.npvs,
                weight=weight,
            )
            if not isRealData:
                output["pu"].fill(
                    syst,
                    ak.values_astype(events[event_level].Pileup.nTrueInt, int),
                    weight=weight,
                )
        #######################
        #  Create root files  #
        #######################
        if self.isArray:
            # Keep the structure of events and pruned the object size
            pruned_ev = events[event_level]
            pruned_ev["SelJet"] = sjets
            pruned_ev["Muon"] = shmu
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
            if not isRealData:
                pruned_ev["weight"] = weights.weight()
                for ind_wei in weights.weightStatistics.keys():
                    pruned_ev[f"{ind_wei}_weight"] = weights.partial_weight(
                        include=[ind_wei]
                    )

            pruned_ev["dr_mujet_softmu"] = ssmu.delta_r(smuon_jet)
            pruned_ev["dr_mujet_lep1"] = shmu.delta_r(smuon_jet)
            pruned_ev["dr_lep1_softmu"] = shmu.delta_r(ssmu)
            pruned_ev["soft_l_ptratio"] = ssmu.pt / smuon_jet.pt
            pruned_ev["l1_ptratio"] = shmu.pt / smuon_jet.pt
            pruned_ev["MuonJet_beta"] = smuon_jet.pt / smuon_jet.E
            pruned_ev["MuonJet_muneuEF"] = smuon_jet.muEF + smuon_jet.neEmEF

            array_writer(self, pruned_ev, events, systematics[0], dataset, isRealData)

        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator
