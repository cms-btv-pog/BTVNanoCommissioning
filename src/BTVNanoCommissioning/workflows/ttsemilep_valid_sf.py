import collections, gc
import os
import uproot
import numpy as np, awkward as ak

from coffea import processor
from coffea.analysis_tools import Weights
from BTVNanoCommissioning.utils.correction import (
    load_lumi,
    load_SF,
    weight_manager,
    common_shifts,
)

from BTVNanoCommissioning.helpers.func import update, dump_lumi, PFCand_link
from BTVNanoCommissioning.helpers.update_branch import missing_branch
from BTVNanoCommissioning.utils.histogrammer import histogrammer, histo_writter
from BTVNanoCommissioning.utils.array_writer import array_writer
from BTVNanoCommissioning.utils.selection import (
    HLT_helper,
    jet_id,
    btag_mu_idiso,
    MET_filters,
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
        selectionModifier="tt_semilep",
    ):
        self._year = year
        self._campaign = campaign
        self.name = name
        self.isSyst = isSyst
        self.isArray = isArray
        self.noHist = noHist
        self.lumiMask = load_lumi(self._campaign)
        self.chunksize = chunksize
        ### Added selection for ttbar semileptonic
        self.ttaddsel = selectionModifier
        ## Load corrections
        self.SF_map = load_SF(self._year, self._campaign)

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
        histname = (
            "ttsemilep_sf" if self.ttaddsel != "c_tt_semilep" else "c_ttsemilep_sf"
        )
        output = {} if self.noHist else histogrammer(events, histname)

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
        triggers = ["IsoMu24"]
        req_trig = HLT_helper(events, triggers)

        ## Muon cuts
        # muon twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
        event_muon = events.Muon[
            (events.Muon.pt > 30) & btag_mu_idiso(events, self._campaign)
        ]
        # event_muon = ak.pad_none(events.Muon, 1, axis=1)
        req_muon = ak.count(event_muon.pt, axis=1) == 1

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
                ),
                False,
                axis=-1,
            )
        ]
        req_jets = ak.num(event_jet.pt) >= 4

        if self.ttaddsel == "c_tt_semilep":

            # Sort all jets by pt in descending order
            sorted_jets = event_jet[ak.argsort(event_jet.pt, axis=-1, ascending=False)]

            # Pad the sorted jets to ensure there are at least two jets in each event
            padded_jets = ak.pad_none(sorted_jets, 2, axis=1)

            # Determine the leading and subleading jets
            leading_jet = padded_jets[:, 0]
            subleading_jet = padded_jets[:, 1]

            c_selection = (
                event_jet.btagDeepFlavC > 0.6
            )  # & (event_jet.pt == leading_jet.pt)
            ### plot # c jets , DeepFlavC of the c jets
            c_selection = ak.fill_none(c_selection, False)
            req_c_jets = ak.num(event_jet[c_selection].pt) >= 1
            req_c_jets = ak.fill_none(req_c_jets, False)
            c_jets = event_jet[c_selection]
            single_c_jet = ak.firsts(c_jets)
            remaining_jets = event_jet[~c_selection]
            nob_selection = remaining_jets.btagDeepFlavB < 0.5

            nob_candidate_jets = remaining_jets[nob_selection]

            # Create combinations of single_c_jet and nob_candidate_jets
            combinations = ak.cartesian(
                {"c_jet": single_c_jet, "nob_jet": nob_candidate_jets}, axis=1
            )

            # Calculate the invariant mass of each pair using the .mass attribute
            masses = (combinations["c_jet"] + combinations["nob_jet"]).mass

            # Create a mask for pairs with mass between 55 and 110
            mass_mask = (masses > 50) & (masses < 110)

            # Apply the mask to get the valid combinations
            valid_combinations = combinations[mass_mask]
            # Create an event-level mask
            event_mask = ak.any(mass_mask, axis=1)

            req_w_mass = ak.fill_none(event_mask, False)

            nob_w_selection = (
                event_jet.btagDeepFlavB < 0.5
            )  # & (event_jet.pt == subleading_jet.pt)
            nob_w_selection = ak.fill_none(nob_w_selection, False)
            w_cand_b_jet = ak.firsts(event_jet[nob_w_selection])

            w_mass_lead_sublead = (single_c_jet + w_cand_b_jet).mass

            w_mask = (w_mass_lead_sublead > 50) & (w_mass_lead_sublead < 110)
            w_mask = ak.fill_none(w_mask, False)

            # req_w_jets = ak.fill_none(req_w_jets, False)

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
            )
            == 1,
        )
        jetindx = ak.pad_none(jetindx, 4)
        jetindx = jetindx[:, :4]

        ## other cuts

        MET = ak.zip(
            {
                "pt": events.MET.pt,
                "eta": ak.zeros_like(events.MET.pt),
                "phi": events.MET.phi,
                "mass": ak.zeros_like(events.MET.pt),
            },
            with_name="PtEtaPhiMLorentzVector",
        )

        req_MET = MET.pt > 50

        req_metfilter = MET_filters(events, self._campaign)
        event_level = ak.fill_none(
            req_trig & req_jets & req_muon & req_MET & req_lumi & req_metfilter, False
        )
        # Combine masks with the original conditions
        if self.ttaddsel == "c_tt_semilep":
            event_level = ak.fill_none(event_level & req_w_mass & req_c_jets, False)
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
        pruned_ev["SelJet"] = (
            event_jet[event_level][:, :4]
            if self.ttaddsel != "c_tt_semilep"
            else c_jets[event_level][:, 0]
        )
        pruned_ev["SelMuon"] = event_muon[event_level][:, 0]
        pruned_ev["njet"] = ak.count(event_jet[event_level].pt, axis=1)
        if "PFCands" in events.fields:
            pruned_ev.PFCands = PFCand_link(events, event_level, jetindx)
        if self.ttaddsel != "c_tt_semilep":
            for i in range(4):
                pruned_ev[f"dr_mujet{i}"] = pruned_ev.SelMuon.delta_r(
                    pruned_ev.SelJet[:, i]
                )
        else:
            pruned_ev[f"dr_mujet0"] = pruned_ev.SelMuon.delta_r(pruned_ev.SelJet)

        ####################
        #     Output       #
        ####################
        # Configure SFs
        weights = weight_manager(pruned_ev, self.SF_map, self.isSyst)
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
