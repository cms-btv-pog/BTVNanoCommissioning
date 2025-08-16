import collections, numpy as np, awkward as ak, os
from coffea import processor
from coffea.analysis_tools import Weights
from BTVNanoCommissioning.helpers.func import flatten, update, dump_lumi
from BTVNanoCommissioning.utils.histogrammer import histogrammer, histo_writter
from BTVNanoCommissioning.utils.array_writer import array_writer
from BTVNanoCommissioning.helpers.update_branch import missing_branch
from BTVNanoCommissioning.utils.correction import (
    load_lumi,
    load_SF,
    weight_manager,
    common_shifts,
)
from BTVNanoCommissioning.utils.selection import jet_cut, HLT_helper

import correctionlib


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
        addsel=False,
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

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        events = missing_branch(events)
        shifts = common_shifts(self, events)

        return processor.accumulate(
            self.process_shift(update(events, collections), name)
            for collections, name in shifts
        )

    def process_shift(self, events, shift_name):
        isRealData = not hasattr(events, "genWeight")
        dataset = events.metadata["dataset"]
        output = {} if self.noHist else histogrammer(events, "QCD")

        if isRealData:
            output["sumw"] = len(events)
        else:
            output["sumw"] = ak.sum(events.genWeight)

        ####################
        #    Selections    #
        ####################
        ## HLT
        # triggers = {trigger1 : [ptmin, ptmax], ...}
        triggers = {
            "PFJet40": [45, 80],
            "PFJet60": [80, 110],
            "PFJet80": [110, 180],
            "PFJet140": [180, 250],
            "PFJet200": [250, 310],
            "PFJet260": [310, 380],
            "PFJet320": [380, 460],
            "PFJet400": [460, 520],
            "PFJet450": [520, 580],
            "PFJet500": [580, 1e7],
        }
        # This has to be in ascending order, so that the prescale weight of the last passed trigger is kept

        jet_mask = jet_cut(events, self._campaign, ptmin=20)
        event_jet = events.Jet[jet_mask]
        req_jets = ak.count(event_jet.pt, axis=1) >= 1

        req_trig = np.zeros(len(events), dtype="bool")
        trigbools = {}
        for trigger, ptrange in triggers.items():
            ptmin = ptrange[0]
            ptmax = ptrange[1]
            # Require *leading jet* to be in the pT range of the trigger
            thistrigreq = (
                (HLT_helper(events, [trigger]))
                & (ak.fill_none(ak.firsts(event_jet.pt) >= ptmin, False))
                & (ak.fill_none(ak.firsts(event_jet.pt) < ptmax, False))
            )
            trigbools[trigger] = thistrigreq
            req_trig = (req_trig) | (thistrigreq)

        req_lumi = np.ones(len(events), dtype="bool")
        if isRealData:
            req_lumi = self.lumiMask(events.run, events.luminosityBlock)
        if shift_name is None:
            output = dump_lumi(events[req_lumi], output)

        event_level = ak.fill_none(req_lumi & req_trig & req_jets, False)
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
        pruned_ev = events[event_level]
        pruned_ev["SelJet"] = event_jet[event_level]
        pruned_ev["njet"] = ak.count(event_jet[event_level].pt, axis=1)
        for trigger in triggers:
            pruned_ev[trigger] = trigbools[trigger][event_level]

        ###FIXME: https://gitlab.cern.ch/cms-btv-coordination/tasks/-/issues/188
        if "JetSVs" in events.fields:
            valid_events = (ak.count(pruned_ev.Jet.pt, axis=1) > 0) & (
                ak.count(pruned_ev.JetSVs.pt, axis=1) > 0
            )
            filtered_events = pruned_ev[valid_events]

            # Pad pruned_ev.JetSVs.jetIdx to match the length of pruned_ev.Jet
            filtered_events.JetSVs.jetIdx = ak.pad_none(
                filtered_events.JetSVs.jetIdx, len(filtered_events.Jet), clip=True
            )

            # Print the initial state of pruned_ev.Jet and pruned_ev.JetSVs.jetIdx
            # print("pruned_ev.Jet", filtered_events.Jet)
            # print("pruned_ev.JetSVs.jetIdx", filtered_events.JetSVs.jetIdx)

            # Filter events where the number of Jet and JetSVs are the same
            # equal_length_events = ak.count(filtered_events.Jet.pt, axis=1) == ak.count(filtered_events.JetSVs.pt, axis=1)
            # same_length_events = filtered_events[equal_length_events]
            # print("Number of events with equal number of Jet and JetSVs:", len(same_length_events))
            # print("Same length pruned_ev.Jet", same_length_events.Jet)
            # print("Same length pruned_ev.JetSVs.jetIdx", same_length_events.JetSVs.jetIdx)

            # Count and print the number of events where JetSVs is longer than Jet
            # longer_jetSVs_events = ak.count(events.JetSVs.pt, axis=1) > ak.count(events.Jet.pt, axis=1)
            # longer_jetSVs_events = longer_jetSVs_events & (ak.count(events.Jet.pt, axis=1) == 0)
            # num_longer_jetSVs_events = ak.sum(longer_jetSVs_events)
            # print("Number of events where JetSVs is longer than Jet:", num_longer_jetSVs_events)
            # print("Longer SV pruned_ev.Jet", events.Jet[longer_jetSVs_events])
            # print("Longer SV pruned_ev.JetSVs.jetIdx", events.JetSVs.jetIdx[longer_jetSVs_events])

            # Ensure that all indices are within the valid range for pruned_ev.Jet
            valid_indices = (filtered_events.JetSVs.jetIdx >= 0) & (
                filtered_events.JetSVs.jetIdx < ak.num(filtered_events.Jet)
            )
            if not np.all(valid_indices):
                print(
                    "Warning: Some indices in pruned_ev.JetSVs.jetIdx are out of range."
                )
                # Filter out invalid indices
                filtered_events.JetSVs = filtered_events.JetSVs[valid_indices]

            # Check if pruned_ev.JetSVs.jetIdx is empty after filtering
            if len(filtered_events.JetSVs.jetIdx) == 0:
                print("Warning: pruned_ev.JetSVs.jetIdx is empty after filtering.")
                matched_JetSVs = ak.Array([])
                lj_matched_JetSVs = ak.Array([])
                lj_SVs = ak.Array([])
                nJetSVs = ak.Array([])
            else:
                # Proceed with the assignment
                matched_JetSVs = filtered_events.Jet[filtered_events.JetSVs.jetIdx]
                lj_matched_JetSVs = matched_JetSVs[filtered_events.JetSVs.jetIdx == 0]
                lj_SVs = filtered_events.JetSVs[filtered_events.JetSVs.jetIdx == 0]
                nJetSVs = ak.count(lj_SVs.pt, axis=1)

            # Print the final state of the variables
            # print("matched_JetSVs:", matched_JetSVs)
            # print("lj_matched_JetSVs:", lj_matched_JetSVs)
            # print("lj_SVs:", lj_SVs)
            # print("nJetSVs:", nJetSVs)

        ####################
        #     Output       #
        ####################
        # Configure SFs
        weights = weight_manager(pruned_ev, self.SF_map, self.isSyst)
        if isRealData:
            if self._year == "2022":
                run_num = "355374_362760"
            elif self._year == "2023":
                run_num = "366727_370790"

            # if 369869 in pruned_ev.run: continue
            psweight = np.zeros(len(pruned_ev))
            for trigger, trigbool in trigbools.items():
                psfile = f"src/BTVNanoCommissioning/data/Prescales/ps_weight_{trigger}_run{run_num}.json"
                if not os.path.isfile(psfile):
                    raise NotImplementedError(
                        f"Prescale weights not available for {trigger} in {self._year}. Please run `scripts/dump_prescale.py`."
                    )
                pseval = correctionlib.CorrectionSet.from_file(psfile)
                thispsweight = pseval["prescaleWeight"].evaluate(
                    pruned_ev.run,
                    f"HLT_{trigger}",
                    ak.values_astype(pruned_ev.luminosityBlock, np.float32),
                )
                psweight = ak.where(trigbool[event_level], thispsweight, psweight)
            weights.add("psweight", psweight)

        if "JetSVs" in events.fields:
            if len(lj_matched_JetSVs) > 0:
                lj_matched_JetSVs_genflav = ak.zeros_like(
                    lj_matched_JetSVs.pt, dtype=int
                )
            else:
                lj_matched_JetSVs_genflav = ak.Array([])

        ####################
        #     Output       #
        ####################
        # Configure systematics
        if shift_name is None:
            systematics = ["nominal"] + list(weights.variations)
        else:
            systematics = [shift_name]

        # Configure histograms
        if not self.noHist:
            output = histo_writter(
                pruned_ev, output, weights, systematics, self.isSyst, self.SF_map
            )
        # Output arrays
        if self.isArray:
            array_writer(
                self,
                pruned_ev,
                events,
                weights,
                systematics,
                dataset,
                isRealData,
                kinOnly=[],
            )

        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator
