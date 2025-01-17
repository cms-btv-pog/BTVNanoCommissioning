# https://git.rwth-aachen.de/3pia/cms_analyses/jet_tagging_sf/-/blob/master/recipes/selection3.py?ref_type=heads

# import collections
import awkward as ak
import numpy as np

# import os
# import uproot
from coffea import processor

# user helper function
from BTVNanoCommissioning.helpers.func import (
    # uproot_writeable,
    dump_lumi,
    update,
    PFCand_link,
)
from BTVNanoCommissioning.helpers.update_branch import missing_branch
from BTVNanoCommissioning.utils.array_writer import array_writer

from BTVNanoCommissioning.utils.selection import HLT_helper

# functions to load SFs, corrections
from BTVNanoCommissioning.utils.correction import (
    common_shifts,
    weight_manager,
    load_lumi,
    load_SF,
)

## load histograms & selctions for this workflow
from BTVNanoCommissioning.utils.histogrammer import histogrammer, histo_writter
from BTVNanoCommissioning.utils.selection import ele_mvatightid, jet_id, mu_idiso


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
        self.SF_map = load_SF(year=self._year, campaign=self._campaign)

    @property
    def accumulator(self):
        return self._accumulator

    ## Apply corrections on momentum/mass on MET, Jet, Muon
    def process(self, events):
        events = missing_branch(events)
        shifts = common_shifts(self, events)

        return processor.accumulate(
            self.process_shift(update(events, collections), name) for collections, name in shifts
        )

    ## Processed events per-chunk, made selections, filled histogram, stored root files
    def process_shift(self, events, shift_name):
        """Selection following
        https://git.rwth-aachen.de/3pia/cms_analyses/jet_tagging_sf/-/blob/master/recipes/selection3.py?ref_type=heads
        """
        dataset = events.metadata["dataset"]
        isRealData = not hasattr(events, "genWeight")
        ######################
        #  Create histogram  #
        ######################
        _hist_event_dict = (
            {"": None} if self.noHist else histogrammer(events, "btag_ttbar_sf")
        )
        if _hist_event_dict is None:
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
        # TODO: check trigger paths
        triggers = [
            "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
            "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
        ]
        req_trig = HLT_helper(events, triggers)

        ##### Add some selections
        ## Muon cuts
        # muon twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2

        # TODO: compare muon/electron selection
        muon_sel = (events.Muon.pt > 15) & (mu_idiso(events, self._campaign))
        event_mu = events.Muon[muon_sel]
        req_muon = ak.num(event_mu.pt) == 1
        event_mu = ak.pad_none(event_mu, 1, axis=1)
        events.Muon = event_mu


        # Electron cut
        ele_sel = (events.Electron.pt > 15) & (ele_mvatightid(events, self._campaign))
        event_e = events.Electron[ele_sel]
        req_ele = ak.num(event_e.pt) == 1
        event_e = ak.pad_none(event_e, 1, axis=1)
        events.Electron = event_e


        req_leadlep_pt = ak.any(event_e.pt > 25, axis=-1) | ak.any(
            event_mu.pt > 25, axis=-1
        )

        ## Jet cuts
        # TODO: check jet selection
        jet_sel = (events.Jet.pt > 20.0) & (jet_id(events, self._campaign))
        event_jet = events.Jet[jet_sel]
        req_jet = ak.num(event_jet) == 2
        # event_jet = ak.pad_none(event_jet, 2, axis=1)

        ## store jet index for PFCands, create mask on the jet index
        jetindx = ak.mask(ak.local_index(events.Jet.pt), jet_sel)
        jetindx = ak.pad_none(jetindx, 2)
        jetindx = jetindx[:, :2]
        ## Other cuts

        ## Apply all selections
        event_level = (
            req_trig & req_lumi & req_jet & req_muon & req_ele & req_leadlep_pt
        )
        event_level = ak.fill_none(event_level, False)
        # Skip empty events
        if len(events[event_level]) == 0:
            if self.isArray:
                array_writer(
                    self,
                    events[event_level],
                    events,
                    "nominal",
                    dataset,
                    isRealData,
                )
            return {dataset: output}

        ####################
        # Selected objects # : Pruned objects with reduced event_level
        ####################
        pruned_ev = events[event_level]
        pruned_ev["SelMuon"] = event_mu[event_level][:, 0]
        # TODO: 2 jets?
        pruned_ev["SelJet"] = event_jet[event_level][:, :2]
        pruned_ev["SelElectron"] = event_e[event_level][:, 0]
        pruned_ev["njet"] = ak.count(event_jet[event_level].pt, axis=1)
        pruned_ev["MET"] = events.MET[event_level]
        if "PFCands" in events.fields:
            pruned_ev.PFCands = PFCand_link(events, event_level, jetindx)
        
        ####################
        #     Output       #
        ####################
        # Configure SFs
        weights = weight_manager(pruned_ev, self.SF_map, self.isSyst)       

        # configure systematics
        if shift_name is None:
            systematics = ["nominal"] + list(weights.variations)
        else:
            systematics = [shift_name]

        if not isRealData:
            pruned_ev["weight"] = weights.weight()
            for ind_wei in weights.weightStatistics.keys():
                pruned_ev[f"{ind_wei}_weight"] = weights.partial_weight(include=[ind_wei])

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
                systematics[0],
                dataset,
                isRealData,
            )

        return {dataset: output}

    ## post process, return the accumulator, compressed
    def postprocess(self, accumulator):
        return accumulator
