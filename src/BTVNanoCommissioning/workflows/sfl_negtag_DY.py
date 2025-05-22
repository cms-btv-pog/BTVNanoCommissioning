import awkward as ak
import numpy as np
import os
import uproot
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
    mu_idiso,
    ele_mvatightid,
    MET_filters,
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
        # Load corrections
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
        dataset = events.metadata["dataset"]
        isRealData = not hasattr(events, "genWeight")

        # TODO: check which triggers to use
        triggers = ["Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8", "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"]

        histname = "sfl_negtag_DY"
        output = {} if self.noHist else histogrammer(events, histname)

        if isRealData:
            output["sumw"] = len(events)
        else:
            output["sumw"] = ak.sum(events.genWeight)

        ####################
        #    Selections    #
        ####################

        # Lumimask
        req_lumi = np.ones(len(events), dtype="bool")
        if isRealData:
            req_lumi = self.lumiMask(events.run, events.luminosityBlock)
        # only dump for nominal case
        if shift_name is None:
            output = dump_lumi(events[req_lumi], output)

        # HLT
        req_trig = HLT_helper(events, triggers)
        req_metfilter = MET_filters(events, self._campaign)

        # Muon cuts
        # TODO: check if muon object selection cuts could be optimised
        muons = events.Muon[(events.Muon.pt > 12) & mu_idiso(events, self._campaign)]

        # Electron cuts
        # TODO: check if electron object selection cuts could be optimised
        electrons = events.Electron[
            (events.Electron.pt > 15) & ele_mvatightid(events, self._campaign)
        ]

        # TODO add requirement to have at least two electrons or at least two muons
        req_dilep = True

        # TODO: implement logic to find lepton pairs from muons and electrons 

        # TODO: implement dilepton reconstruction
        # # dilepton
        # dilep_mass = lepton0[:, 0] + lepton1[:, 0]
        # req_dilepmass = (
        #     (dilep_mass.mass > 81) & (dilep_mass.mass < 101) & (dilep_mass.pt > 15)
        # )

        # TODO: add jet selection
        jets = events.Jet[
            (events.Jet.pt > 20) # TODO: add more selection cuts (e.g. check utils/selection.py for jet_id(...))
        ]

        req_jets = ak.num(jets.pt) >= 1

        # event level selection
        event_level = ak.fill_none(
            req_lumi & req_trig & req_jets, # & req_dilep & req_dilepmass,
            False,
        )
        # return empty arrays if no events are surviving the selection
        # TODO: implement empty histograms
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
        # TODO: add selected objects, such as electron pair / muon pair from Z and leading jet in event


        # Only store events which survive event selection ("pruned events")
        pruned_ev = events[event_level]
        # jets
        pruned_ev["SelJet"] = jets[event_level][:, :1]
        pruned_ev["njet"] = ak.count(jets[event_level].pt, axis=1)

        # muons
        # pruned_ev["posl"] = sposmu
        # pruned_ev["negl"] = snegmu
        # pruned_ev["SelMuon"] = smu

        # electrons
        # pruned_ev["posl"] = sposmu
        # pruned_ev["negl"] = snegmu
        # pruned_ev["SelElectron"] = smu

        # dilepton
        # pruned_ev["dilep"] = sz
        # pruned_ev["dilep", "pt"] = pruned_ev.dilep.pt
        # pruned_ev["dilep", "eta"] = pruned_ev.dilep.eta
        # pruned_ev["dilep", "phi"] = pruned_ev.dilep.phi
        # pruned_ev["dilep", "mass"] = pruned_ev.dilep.mass

        # Find the PFCands associate with selected jets. Search from jetindex->JetPFCands->PFCand
        if "PFCands" in events.fields:
            pruned_ev["PFCands"] = PFCand_link(events, event_level, jetindx)

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
        # TODO: add histograms to utils/histogrammer.py
        if not self.noHist:
            output = histo_writter(
                pruned_ev, output, weights, systematics, self.isSyst, self.SF_map
            )
        # Output arrays
        # TODO: add output arrays to utils/array_writer.py
        if self.isArray:
            array_writer(
                self, pruned_ev, events, weights, systematics, dataset, isRealData
            )

        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator
