import gc, collections, pickle, os, sys, numpy as np, awkward as ak
import os
import uproot

from coffea import processor
from coffea.analysis_tools import Weights

from BTVNanoCommissioning.utils.correction import (
    load_lumi,
    load_SF,
    eleSFs,
    puwei,
    btagSFs,
    JME_shifts,
    Roccor_shifts,
)
from BTVNanoCommissioning.helpers.func import (
    flatten,
    update,
    uproot_writeable,
    dump_lumi,
)
from BTVNanoCommissioning.helpers.update_branch import missing_branch
from BTVNanoCommissioning.utils.histogrammer import histogrammer
from BTVNanoCommissioning.utils.selection import (
    jet_id,
    mu_idiso,
    ele_mvatightid,
    softmu_mask,
)


class NanoProcessor(processor.ProcessorABC):
    # Define histograms
    def __init__(
        self,
        year="2022",
        campaign="Summer22Run3",
        name="",
        isSyst=False,
        isArray=False,
        noHist=False,
        chunksize=75000,
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
            if int(self._year) > 2020:
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
        _hist_event_dict = (
            {"": None} if self.noHist else histogrammer(events, "ectag_ttsemilep_sf")
        )
        if _hist_event_dict == None:
            _hist_event_dict[""]: None
        output = {
            "sumw": processor.defaultdict_accumulator(float),
            **_hist_event_dict,
        }

        if shift_name is None:
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
        triggers = ["Ele32_WPTight_Gsf_L1DoubleEG"]
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

        ## Electron cuts
        iso_ele = events.Electron[
            (events.Electron.pt > 34) & ele_mvatightid(events, self._campaign)
        ]
        req_ele = ak.count(iso_ele.pt, axis=1) == 1
        iso_ele = ak.pad_none(iso_ele, 1, axis=1)
        iso_ele = iso_ele[:, 0]
        iso_eindx = ak.mask(
            ak.local_index(events.Electron.pt),
            ((events.Electron.pt > 34) & ele_mvatightid(events, self._campaign)) == 1,
        )
        iso_eindx = ak.pad_none(iso_eindx, 1)
        iso_eindx = iso_eindx[:, 0]
        ## Jet cuts
        event_jet = events.Jet[
            ak.fill_none(
                jet_id(events, self._campaign)
                & (
                    ak.all(
                        events.Jet.metric_table(iso_ele) > 0.5,
                        axis=2,
                        mask_identity=True,
                    )
                ),
                False,
                axis=-1,
            )
        ]
        req_jets = ak.num(event_jet.pt) >= 4

        ## Soft Muon cuts
        soft_muon = events.Muon[softmu_mask(events, self._campaign)]
        req_softmu = ak.count(soft_muon.pt, axis=1) >= 1
        soft_muon = ak.pad_none(soft_muon, 1, axis=1)

        ## Muon-jet cuts
        mu_jet = event_jet[
            ak.fill_none(
                (
                    ak.all(
                        event_jet.metric_table(soft_muon) <= 0.4,
                        axis=2,
                        mask_identity=True,
                    )
                )
                & ((event_jet.muonIdx1 != -1) | (event_jet.muonIdx2 != -1)),
                False,
                axis=-1,
            )
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
                        events.Jet.metric_table(iso_ele) > 0.5,
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

        # Other cuts
        req_pTratio = (soft_muon[:, 0].pt / mu_jet[:, 0].pt) < 0.6

        req_QCDveto = (
            (iso_ele.pfRelIso03_all < 0.05)
            & (abs(iso_ele.dz) < 0.02)
            & (abs(iso_ele.dxy) < 0.01)
            & (iso_ele.sip3d < 2.5)
            & (
                iso_ele.pt
                / ak.firsts(
                    events.Jet[
                        (events.Jet.electronIdx1 == iso_eindx)
                        | ((events.Jet.electronIdx2 == iso_eindx))
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

        MET = ak.zip(
            {
                "pt": events.MET.pt,
                "eta": ak.zeros_like(events.MET.pt),
                "phi": events.MET.phi,
                "mass": ak.zeros_like(events.MET.pt),
            },
            with_name="PtEtaPhiMLorentzVector",
        )
        Wmass = MET + iso_ele
        req_Wmass = Wmass.mass > 55

        event_level = (
            req_trig
            & req_lumi
            & req_ele
            & req_jets
            & req_softmu
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
        shmu = iso_ele[event_level]
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
        njet = ak.count(sjets.pt, axis=1)
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
            smflav = (
                1 * ((smuon_jet.partonFlavour == 0) & (smuon_jet.hadronFlavour == 0))
                + smuon_jet.hadronFlavour
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
                if "EGM" in self.SF_map.keys():
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
            "DeepJetB",
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
                        flatten(sjets[histname]),
                        weight=flatten(
                            ak.broadcast_arrays(
                                weights.partial_weight(exclude=exclude_btv), sjets["pt"]
                            )[0]
                        ),
                    )
                elif "jet_" in histname and "mu" not in histname:
                    h.fill(
                        syst,
                        flatten(genflavor),
                        flatten(sjets[histname.replace("jet_", "")]),
                        weight=flatten(
                            ak.broadcast_arrays(
                                weights.partial_weight(exclude=exclude_btv), sjets["pt"]
                            )[0]
                        ),
                    )
                elif "hl_" in histname and histname.replace("hl_", "") in shmu.fields:
                    h.fill(
                        syst,
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
                        flatten(ssmu[histname.replace("soft_l_", "")]),
                        weight=weight,
                    )
                elif "mujet_" in histname:
                    h.fill(
                        syst,
                        smflav,
                        flatten(smuon_jet[histname.replace("mujet_", "")]),
                        weight=weight,
                    )
                elif (
                    "PFCands" in events.fields
                    and "PFCands" in histname
                    and histname.split("_")[1] in events.PFCands.fields
                ):
                    h.fill(
                        syst,
                        flatten(ak.broadcast_arrays(smflav, spfcands["pt"])[0]),
                        flatten(spfcands[histname.replace("PFCands_", "")]),
                        weight=flatten(
                            ak.broadcast_arrays(
                                weights.partial_weight(exclude=exclude_btv),
                                spfcands["pt"],
                            )[0]
                        ),
                    )
                elif "btag" in histname:
                    for i in range(2):
                        if (
                            str(i) not in histname
                            or histname.replace(f"_{i}", "") not in events.Jet.fields
                        ):
                            continue
                        h.fill(
                            syst="noSF",
                            flav=smflav,
                            discr=smuon_jet[histname.replace(f"_{i}", "")],
                            weight=weights.partial_weight(exclude=exclude_btv),
                        )
                        if not isRealData and "btag" in self.SF_map.keys():
                            h.fill(
                                syst=syst,
                                flav=smflav,
                                discr=smuon_jet[histname.replace(f"_{i}", "")],
                                weight=weight,
                            )
            output["njet"].fill(syst, njet, weight=weight)
            output["nmujet"].fill(syst, nmujet, weight=weight)
            output["nsoftmu"].fill(syst, nsoftmu, weight=weight)
            output["hl_ptratio"].fill(
                syst,
                genflavor[:, 0],
                ratio=shmu.pt / sjets[:, 0].pt,
                weight=weight,
            )
            output["soft_l_ptratio"].fill(
                syst,
                flav=smflav,
                ratio=ssmu.pt / smuon_jet.pt,
                weight=weight,
            )
            output["dr_lmujetsmu"].fill(
                syst,
                flav=smflav,
                dr=smuon_jet.delta_r(ssmu),
                weight=weight,
            )
            output["dr_lmujethmu"].fill(
                syst,
                flav=smflav,
                dr=smuon_jet.delta_r(shmu),
                weight=weight,
            )
            output["dr_lmusmu"].fill(syst, dr=shmu.delta_r(ssmu), weight=weight)
            output["z_pt"].fill(syst, flatten(sz.pt), weight=weight)
            output["z_eta"].fill(syst, flatten(sz.eta), weight=weight)
            output["z_phi"].fill(syst, flatten(sz.phi), weight=weight)
            output["z_mass"].fill(syst, flatten(sz.mass), weight=weight)
            output["w_pt"].fill(syst, flatten(sw.pt), weight=weight)
            output["w_eta"].fill(syst, flatten(sw.eta), weight=weight)
            output["w_phi"].fill(syst, flatten(sw.phi), weight=weight)
            output["w_mass"].fill(syst, flatten(sw.mass), weight=weight)
            output["MET_pt"].fill(syst, flatten(smet.pt), weight=weight)
            output["MET_phi"].fill(syst, flatten(smet.phi), weight=weight)
            output["npvs"].fill(
                syst,
                events[event_level].PV.npvs,
                weight=weight,
            )
            if not isRealData:
                output["pu"].fill(
                    syst,
                    events[event_level].Pileup.nTrueInt,
                    weight=weight,
                )
        #######################
        #  Create root files  #
        #######################
        if self.isArray:
            # Keep the structure of events and pruned the object size
            pruned_ev = events[event_level]
            pruned_ev.Jet = sjets
            pruned_ev.Muon = shmu
            pruned_ev["MuonJet"] = smuon_jet
            pruned_ev["SoftMuon"] = ssmu
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

            # Create a list of variables want to store. For objects from the PFNano file, specify as {object}_{variable}, wildcard option only accepted at the end of the string
            out_branch = np.setdiff1d(
                np.array(pruned_ev.fields), np.array(events.fields)
            )
            out_branch = np.delete(
                out_branch,
                np.where(
                    (out_branch == "SoftMuon")
                    | (out_branch == "MuonJet")
                    | (out_branch == "dilep")
                ),
            )

            for kin in ["pt", "eta", "phi", "mass", "pfRelIso04_all", "dxy", "dz"]:
                for obj in [
                    "Muon",
                    "Jet",
                    "SoftMuon",
                    "MuonJet",
                    "dilep",
                    "charge",
                    "MET",
                ]:
                    if "MET" in obj and ("pt" != kin or "phi" != kin):
                        continue
                    if (obj != "Muon" and obj != "SoftMuon") and (
                        "pfRelIso04_all" == kin or "d" in kin
                    ):
                        continue
                    out_branch = np.append(out_branch, [f"{obj}_{kin}"])
            out_branch = np.append(
                out_branch, ["Jet_btagDeep*", "Jet_DeepJet*", "PFCands_*", "SV_*"]
            )
            # write to root files
            os.system(f"mkdir -p {self.name}/{dataset}")
            with uproot.recreate(
                f"{self.name}/{dataset}/f{events.metadata['filename'].split('_')[-1].replace('.root','')}_{systematics[0]}_{int(events.metadata['entrystop']/self.chunksize)}.root"
            ) as fout:
                fout["Events"] = uproot_writeable(pruned_ev, include=out_branch)
        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator
