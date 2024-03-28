import collections, awkward as ak, numpy as np
import os
import uproot
from coffea import processor
from coffea.analysis_tools import Weights

# functions to load SFs, corrections
from BTVNanoCommissioning.utils.correction import (
    load_lumi,
    load_SF,
    eleSFs,
    muSFs,
    puwei,
    btagSFs,
    JME_shifts,
    Roccor_shifts,
)

# user helper function
from BTVNanoCommissioning.helpers.func import (
    flatten,
    update,
    uproot_writeable,
    dump_lumi,
)
from BTVNanoCommissioning.helpers.update_branch import missing_branch

## load histograms & selctions for this workflow
from BTVNanoCommissioning.utils.histogrammer import histogrammer
from BTVNanoCommissioning.utils.selection import (
    jet_id,
    mu_idiso,
    ele_cuttightid,
    btag_wp,
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

    ## Apply corrections on momentum/mass on MET, Jet, Muon
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
        ## Create histograms
        _hist_event_dict = (
            {"": None} if self.noHist else histogrammer(events, "ttdilep_sf")
        )
        if _hist_event_dict == None:
            _hist_event_dict[""]
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
        triggers = [
            "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
            "Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
            "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
        ]
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
        events.Muon = events.Muon[
            (events.Muon.pt > 30) & mu_idiso(events, self._campaign)
        ]
        events.Muon = ak.pad_none(events.Muon, 1, axis=1)
        req_muon = ak.count(events.Muon.pt, axis=1) == 1

        ## Electron cuts
        # electron twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
        events.Electron = events.Electron[
            (events.Electron.pt > 30) & ele_cuttightid(events, self._campaign)
        ]
        events.Electron = ak.pad_none(events.Electron, 1, axis=1)
        req_ele = ak.count(events.Electron.pt, axis=1) == 1

        ## Jet cuts
        event_jet = events.Jet[
            ak.fill_none(
                jet_id(events, self._campaign)
                & (
                    ak.all(
                        events.Jet.metric_table(events.Muon) > 0.4,
                        axis=2,
                        mask_identity=True,
                    )
                )
                & (
                    ak.all(
                        events.Jet.metric_table(events.Electron) > 0.4,
                        axis=2,
                        mask_identity=True,
                    )
                ),
                False,
            )
        ]
        req_jets = ak.num(event_jet.pt) >= 2

        ## Other cuts
        req_opposite_charge = (
            events.Electron[:, 0:1].charge * events.Muon[:, 0:1].charge
        ) == -1
        req_opposite_charge = ak.fill_none(req_opposite_charge, False)
        req_opposite_charge = ak.flatten(req_opposite_charge)

        ## store jet index for PFCands, create mask on the jet index
        jetindx = ak.mask(
            ak.local_index(events.Jet.pt),
            (
                jet_id(events, self._campaign)
                & (
                    ak.all(
                        events.Jet.metric_table(events.Muon) > 0.4,
                        axis=2,
                        mask_identity=True,
                    )
                )
                & (
                    ak.all(
                        events.Jet.metric_table(events.Electron) > 0.4,
                        axis=2,
                        mask_identity=True,
                    )
                )
            )
            == 1,
        )
        jetindx = ak.pad_none(jetindx, 2)
        jetindx = jetindx[:, :2]

        event_level = (
            req_trig & req_lumi & req_muon & req_ele & req_jets & req_opposite_charge
        )
        event_level = ak.fill_none(event_level, False)
        if len(events[event_level]) == 0:
            return {dataset: output}

        ####################
        # Selected objects #
        ####################
        smu = events.Muon[event_level]
        smu = smu[:, 0]
        sel = events.Electron[event_level]
        sel = sel[:, 0]
        sjets = event_jet[event_level]
        nseljet = ak.count(sjets.pt, axis=1)
        sjets = sjets[:, :2]
        # Find the PFCands associate with selected jets. Search from jetindex->JetPFCands->PFCand
        if "PFCands" in events.fields:
            jetindx0 = jetindx[:, 0]
            jetindx1 = jetindx[:, 1]
            spfcands = collections.defaultdict(dict)
            spfcands[0] = events[event_level].PFCands[
                events[event_level]
                .JetPFCands[
                    events[event_level].JetPFCands.jetIdx == jetindx0[event_level]
                ]
                .pFCandsIdx
            ]
            spfcands[1] = events[event_level].PFCands[
                events[event_level]
                .JetPFCands[
                    events[event_level].JetPFCands.jetIdx == jetindx1[event_level]
                ]
                .pFCandsIdx
            ]
        ####################
        # Weight & Geninfo #
        ####################
        weights = Weights(len(events[event_level]), storeIndividual=True)
        if not isRealData:
            weights.add("genweight", events[event_level].genWeight)
            par_flav = (sjets.partonFlavour == 0) & (sjets.hadronFlavour == 0)
            genflavor = ak.values_astype(sjets.hadronFlavour + 1 * par_flav, int)
            if len(self.SF_map.keys()) > 0:
                syst_wei = True if self.isSyst != False else False
                if "PU" in self.SF_map.keys():
                    puwei(
                        events[event_level].Pileup.nTrueInt,
                        self.SF_map,
                        weights,
                        syst_wei,
                    )
                if "MUO" in self.SF_map.keys():
                    muSFs(smu, self.SF_map, weights, syst_wei, False)
                if "EGM" in self.SF_map.keys():
                    eleSFs(sel, self.SF_map, weights, syst_wei, False)
                if "BTV" in self.SF_map.keys():
                    btagSFs(sjets, self.SF_map, weights, "DeepJetC", syst_wei)
                    btagSFs(sjets, self.SF_map, weights, "DeepJetB", syst_wei)
                    btagSFs(sjets, self.SF_map, weights, "DeepCSVB", syst_wei)
                    btagSFs(sjets, self.SF_map, weights, "DeepCSVC", syst_wei)
        else:
            genflavor = ak.zeros_like(sjets.pt, dtype=int)

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
                    if syst != "nominal":
                        continue
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
                elif (
                    "PFCands" in events.fields
                    and "PFCands" in histname
                    and histname.split("_")[1] in events.PFCands.fields
                ):
                    if syst != "nominal":
                        continue
                    for i in range(2):
                        h.fill(
                            syst,
                            flatten(
                                ak.broadcast_arrays(
                                    genflavor[:, i],
                                    spfcands[i]["pt"],
                                )[0]
                            ),
                            flatten(spfcands[i][histname.replace("PFCands_", "")]),
                            weight=flatten(
                                ak.broadcast_arrays(
                                    weights.partial_weight(exclude=exclude_btv),
                                    spfcands[i]["pt"],
                                )[0]
                            ),
                        )

                elif "btag" in histname:
                    for i in range(2):
                        sel_jet = sjets[:, i]
                        if (
                            str(i) in histname
                            and histname.replace(f"_{i}", "") in events.Jet.fields
                        ):
                            h.fill(
                                "noSF",
                                flav=flatten(genflavor[:, i]),
                                discr=flatten(sel_jet[histname.replace(f"_{i}", "")]),
                                weight=weights.partial_weight(exclude=exclude_btv),
                            )
                            if not isRealData and "btag" in self.SF_map.keys():
                                h.fill(
                                    syst=syst,
                                    flav=flatten(
                                        ak.values_astype(
                                            genflavor[:, i],
                                            np.uint8,
                                        )
                                    ),
                                    discr=flatten(
                                        sel_jet[histname.replace(f"_{i}", "")]
                                    ),
                                    weight=weight,
                                )
                elif "mu_" in histname and histname.replace("mu_", "") in smu.fields:
                    h.fill(
                        syst,
                        flatten(smu[histname.replace("mu_", "")]),
                        weight=weight,
                    )
                elif "ele_" in histname and histname.replace("ele_", "") in sel.fields:
                    h.fill(
                        syst,
                        flatten(sel[histname.replace("ele_", "")]),
                        weight=weight,
                    )
                elif "jet" in histname and "dr" not in histname and "njet" != histname:
                    for i in range(2):
                        sel_jet = sjets[:, i]
                        if str(i) in histname:
                            h.fill(
                                syst,
                                flatten(genflavor[:, i]),
                                flatten(sel_jet[histname.replace(f"jet{i}_", "")]),
                                weight=weight,
                            )

            for i in range(2):
                output[f"dr_mujet{i}"].fill(
                    syst,
                    flav=flatten(genflavor[:, i]),
                    dr=flatten(smu.delta_r(sjets[:, i])),
                    weight=weight,
                )
            output["njet"].fill(syst, nseljet, weight=weight)
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
            pruned_ev.Electron = sel
            pruned_ev.Muon = smu
            if "PFCands" in events.fields:
                pruned_ev.PFCands = spfcands
            # Add custom variables
            if not isRealData:
                pruned_ev["weight"] = weights.weight()
                for ind_wei in weights.weightStatistics.keys():
                    pruned_ev[f"{ind_wei}_weight"] = weights.partial_weight(
                        include=[ind_wei]
                    )

            pruned_ev["dr_mujet0"] = smu.delta_r(sjets[:, 0])
            pruned_ev["dr_mujet1"] = smu.delta_r(sjets[:, 1])

            # Create a list of variables want to store. For objects from the PFNano file, specify as {object}_{variable}, wildcard option only accepted at the end of the string
            out_branch = np.setdiff1d(
                np.array(pruned_ev.fields), np.array(events.fields)
            )
            for kin in ["pt", "eta", "phi", "mass", "dz", "dxy"]:
                for obj in ["Jet", "Electron", "Muon"]:
                    if obj == "Jet" and "d" in kin:
                        continue
                    out_branch = np.append(out_branch, [f"{obj}_{kin}"])
            out_branch = np.append(
                out_branch,
                [
                    "Jet_btagDeep*",
                    "Jet_DeepJet*",
                    "PFCands_*",
                    "Electron_pfRelIso03_all",
                    "Muon_pfRelIso03_all",
                ],
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
