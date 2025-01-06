import collections, awkward as ak, numpy as np
import os
import uproot
from coffea import processor
from coffea.analysis_tools import Weights

# functions to load SFs, corrections
from BTVNanoCommissioning.utils.correction import (
    load_lumi,
    load_SF,
    common_shifts,
    weight_manager,
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
from BTVNanoCommissioning.utils.histogrammer import histogrammer, histo_writter
from BTVNanoCommissioning.utils.array_writer import array_writer
from BTVNanoCommissioning.utils.selection import (
    HLT_helper,
    jet_id,
    mu_idiso,
    ele_cuttightid,
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
        self.SF_map = load_SF(self._year, self._campaign)

    @property
    def accumulator(self):
        return self._accumulator

    ## Apply corrections on momentum/mass on MET, Jet, Muon
    def process(self, events):
        events = missing_branch(events)
        shifts = common_shifts(self, events)

        return processor.accumulate(
            self.process_shift(update(events, collections), name)
            for collections, name in shifts
        )

    ## Processed events per-chunk, made selections, filled histogram, stored root files
    def process_shift(self, events, shift_name):
        dataset = events.metadata["dataset"]
        isRealData = not hasattr(events, "genWeight")
        ######################
        #  Create histogram  # : Get the histogram dict from `histogrammer`
        ######################
        _hist_event_dict = (
            {"": None}
            if self.noHist
            else histogrammer(events, "example")  # this is the place to modify
        )

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
        ##====> start here, make your customize modification
        ## HLT
        triggers = [
            "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
            "Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
            "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
        ]
        req_trig = HLT_helper(events, triggers)

        ##### Add some selections
        ## Muon cuts
        # muon twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2

        muon_sel = (events.Muon.pt > 15) & (mu_idiso(events, self._campaign))
        event_mu = events.Muon[muon_sel]
        req_muon = ak.num(event_mu.pt) == 1

        # Electron cut
        ele_sel = (events.Electron.pt > 15) & (ele_cuttightid(events, self._campaign))
        event_e = events.Electron[ele_sel]
        req_ele = ak.num(event_e.pt) == 1

        req_leadlep_pt = ak.any(event_e.pt > 25, axis=-1) | ak.any(
            event_mu.pt > 25, axis=-1
        )

        ## Jet cuts
        jet_sel = (events.Jet.pt > 30) & (jet_id(events, self._campaign))
        event_jet = events.Jet[jet_sel]
        req_jet = ak.num(event_jet) >= 1
        ## Other cuts

        ## Apply all selections
        event_level = (
            req_trig & req_lumi & req_jet & req_muon & req_ele & req_leadlep_pt
        )

        ##<==== finish selection
        event_level = ak.fill_none(event_level, False)
        # Skip empty events -
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
        ##===>  Ntuplization  : store custom information
        ####################
        # Selected objects # : Pruned objects with reduced event_level
        ####################
        # Keep the structure of events and pruned the object size
        pruned_ev = events[event_level]
        pruned_ev["SelJet"] = event_jet[event_level][:, 0]
        pruned_ev["SelMuon"] = event_mu[event_level][:, 0]
        pruned_ev["SelElectron"] = event_e[event_level][:, 0]
        pruned_ev["mujet_ptratio"] = pruned_ev.Muon.pt / pruned_ev.SelJet.pt
        pruned_ev["mujet_dr"] = pruned_ev.Muon.delta_r(pruned_ev.SelJet)

        ## <========= end: store custom objects
        ####################
        #     Output       #
        ####################
        # Configure SFs - read pruned objects from the pruned_ev and apply SFs and call the systematics
        weights = weight_manager(pruned_ev, self.SF_map, self.isSyst)
        # Configure systematics shifts
        if shift_name is None:
            systematics = ["nominal"] + list(
                weights.variations
            )  # nominal + weight variation systematics
        else:
            systematics = [shift_name]  # JES/JER systematics

        # Fill the weight to output arrys
        if not isRealData:
            pruned_ev["weight"] = weights.weight()
            for ind_wei in weights.weightStatistics.keys():
                pruned_ev[f"{ind_wei}_weight"] = weights.partial_weight(
                    include=[ind_wei]
                )
        # Configure histograms- fill the histograms with pruned objects
        if not self.noHist:
            output = histo_writter(
                pruned_ev, output, weights, systematics, self.isSyst, self.SF_map
            )
        # Output arrays - store the pruned objects in the output arrays
        if self.isArray:
            array_writer(self, pruned_ev, events, systematics[0], dataset, isRealData)

        return {dataset: output}

    ## post process, return the accumulator, compressed
    def postprocess(self, accumulator):
        return accumulator
