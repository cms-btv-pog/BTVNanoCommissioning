import collections, numpy as np, awkward as ak
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
        output = {"": None} if self.noHist else histogrammer(events, "QCD")

        if isRealData:
            output["sumw"] = len(events)
        else:
            output["sumw"] = ak.sum(events.genWeight)

        ####################
        #    Selections    #
        ####################
        ## HLT
        triggers = [
            "PFJet140",
        ]
        req_trig = HLT_helper(events, triggers)
        req_lumi = np.ones(len(events), dtype="bool")
        if isRealData:
            req_lumi = self.lumiMask(events.run, events.luminosityBlock)
        if shift_name is None:
            output = dump_lumi(events[req_lumi], output)
        ## Jet cuts
        event_jet = events.Jet[jet_cut(events, self._campaign)]
        req_jets = ak.count(event_jet.pt, axis=1) >= 1

        event_level = ak.fill_none(req_lumi & req_trig & req_jets, False)
        if len(events[event_level]) == 0:
            return {dataset: output}

        ####################
        # Selected objects #
        ####################
        pruned_ev = events[event_level]
        pruned_ev["SelJet"] = event_jet[event_level][:, 0]
        pruned_ev["njet"] = ak.count(event_jet[event_level].pt, axis=1)

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
            pseval = correctionlib.CorrectionSet.from_file(
                f"src/BTVNanoCommissioning/data/Prescales/ps_weight_{triggers[0]}_run{run_num}.json"
            )
            # if 369869 in pruned_ev.run: continue
            psweight = pseval["prescaleWeight"].evaluate(
                pruned_ev.run,
                f"HLT_{triggers[0]}",
                ak.values_astype(pruned_ev.luminosityBlock, np.float32),
            )
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
                self, pruned_ev, events, weights, systematics, dataset, isRealData
            )

        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator
