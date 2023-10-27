import numpy as np, awkward as ak
import os
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
            syst_JERC = True if self.isSyst != None else False
            if self.isSyst == "JERC_split":
                syst_JERC = "split"
            shifts = JME_shifts(
                shifts, self.SF_map, events, self._campaign, isRealData, syst_JERC
            )
        else:
            if "Run3" not in self._campaign:
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
        isRealData = not hasattr(events, "genWeight")
        dataset = events.metadata["dataset"]
        _hist_event_dict = {"": None} if self.noHist else histogrammer(events, "ttcom")

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
        triggers = [
            "Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
            "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
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
        events.Muon = events.Muon[(events.Muon.pt > 30) & (abs(events.Muon.eta) <= 2.4)]
        req_muon = ak.count(events.Muon.pt, axis=1) == 1
        events.Muon = ak.pad_none(events.Muon, 1)

        ## Electron cuts
        # electron twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
        events.Electron = events.Electron[
            (events.Electron.pt > 30) & (abs(events.Electron.eta) <= 2.4)
        ]
        req_ele = ak.count(events.Electron.pt, axis=1) == 1
        events.Electron = ak.pad_none(events.Electron, 1)

        ## Jet cuts
        events.Jet = events.Jet[jet_id(events, self._campaign)]
        req_jets = ak.count(events.Jet.pt, axis=1) >= 2

        ## Other cuts
        req_opposite_charge = (
            events.Electron[:, 0].charge * events.Muon[:, 0].charge
        ) == -1
        event_level = ak.fill_none(
            req_trig & req_jets & req_ele & req_muon & req_opposite_charge, False
        )
        if len(events[event_level]) == 0:
            return {dataset: output}

        ####################
        # Selected objects #
        ####################
        sele = events.Electron[event_level]
        sele = sele[:, 0]
        smu = events.Muon[event_level]
        smu = smu[:, 0]
        sjets = events.Jet[event_level]
        sjets = sjets[:, :2]

        ####################
        # Weight & Geninfo #
        ####################
        weights = Weights(len(events[event_level]), storeIndividual=True)
        if not isRealData:
            weights.add("genweight", events[event_level].genWeight)
            par_flav = (sjets.partonFlavour == 0) & (sjets.hadronFlavour == 0)
            genflavor = sjets.hadronFlavour + 1 * par_flav
            genweiev = ak.flatten(ak.broadcast_arrays(weights.weight(), sjets["pt"])[0])
            if len(self.SF_map.keys()) > 0:
                syst_wei = True if self.isSyst != None else False
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
                    eleSFs(sele, self.SF_map, weights, syst_wei, False)
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
            if self.isSyst == None and syst != "nominal":
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
                elif "jet" in histname:
                    for i in range(2):
                        if histname.replace(f"jet{i}_", "") not in sjets.fields:
                            continue
                        jet = sjets[:, i]
                        h.fill(
                            syst,
                            flatten(genflavor[:, i]),
                            flatten(jet[histname.replace(f"jet{i}_", "")]),
                            weight=weight,
                        )
                elif "mu" in histname and histname.replace("mu_", "") in smu.fields:
                    h.fill(
                        syst,
                        flatten(smu[histname.replace("mu_", "")]),
                        weight=weight,
                    )
                elif "ele" in histname and histname.replace("ele_", "") in sele.fields:
                    h.fill(
                        syst,
                        flatten(sele[histname.replace("ele_", "")]),
                        weight=weight,
                    )

            output["dr_mujet0"].fill(
                syst,
                flatten(genflavor[:, 0]),
                flatten(sjets[:, 0].delta_r(smu)),
                weight=weight,
            )
            output["dr_mujet1"].fill(
                syst,
                flatten(genflavor[:, 1]),
                flatten(sjets[:, 1].delta_r(smu)),
                weight=weight,
            )

        #######################
        #  Create root files  #
        #######################
        if self.isArray:
            # Keep the structure of events and pruned the object size
            pruned_ev = events[event_level]
            pruned_ev.Jet = sjets
            pruned_ev.Electron = sele
            pruned_ev.Muon = smu
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
            for kin in ["pt", "eta", "phi", "mass"]:
                for obj in ["Jet", "Electron", "Muon"]:
                    out_branch = np.append(out_branch, [f"{obj}_{kin}"])
            out_branch = np.append(
                out_branch, ["Jet_btagDeep*", "Jet_DeepJet*", "PFCands_*"]
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
