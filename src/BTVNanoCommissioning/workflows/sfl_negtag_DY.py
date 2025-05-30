import awkward as ak
import numpy as np
import os  # noqa: F401
import uproot  # noqa: F401
from coffea import processor
from coffea.analysis_tools import Weights  # noqa: F401

from BTVNanoCommissioning.utils.correction import (
    load_lumi,
    load_SF,
    weight_manager,
    common_shifts,
)

from BTVNanoCommissioning.helpers.func import update, dump_lumi
from BTVNanoCommissioning.helpers.update_branch import missing_branch
from BTVNanoCommissioning.utils.histogrammer import histogrammer, histo_writter
from BTVNanoCommissioning.utils.array_writer import array_writer
from BTVNanoCommissioning.utils.selection import (
    HLT_helper,
    jet_id,  # noqa: F401
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
        triggers = [
            "Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
            "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
        ]

        # print("Available HLT paths:\n", events.HLT.fields)

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
        muons = events.Muon[
            (events.Muon.pt > 12)
            & mu_idiso(events, self._campaign)
            & (abs(events.Muon.eta) < 2.5)
        ]
        muons_tag = ak.with_field(
            muons, ak.full_like(muons.pt, 13, dtype=int), "flavor"
        )
        # Electron cuts
        # TODO: check if electron object selection cuts could be optimised
        electrons = events.Electron[
            (events.Electron.pt > 15)
            & ele_mvatightid(events, self._campaign)
            & (abs(events.Electron.eta) < 2.5)
        ]

        electrons_tag = ak.with_field(
            electrons, ak.full_like(electrons.pt, 11, dtype=int), "flavor"
        )

        # TODO add requirement to have at least two electrons or at least two muons
        nmu = ak.num(muons, axis=1)
        nel = ak.num(electrons, axis=1)
        req_dilep = (nmu >= 2) | (nel >= 2)

        # TODO: implement logic to find lepton pairs from muons and electrons
        leptons = ak.concatenate([muons_tag, electrons_tag], axis=1)
        dilepton = ak.combinations(leptons, 2, fields=["l1", "l2"])
        dilepton = dilepton[
            (dilepton.l1.charge != dilepton.l2.charge)
            & (dilepton.l1.flavor == dilepton.l2.flavor)
        ]
        closest_to_Z = ak.argmin(
            abs((dilepton.l1 + dilepton.l2).mass - 91),
            axis=1,
            keepdims=True,
        )
        closest_to_Z = closest_to_Z[
            ~ak.is_none(closest_to_Z, axis=-1)
        ]  # equivalent to ak.drop_none()
        dilepton = dilepton[closest_to_Z]
        leptons = ak.concatenate([dilepton.l1, dilepton.l2], axis=1)  # update leptons
        # TODO: implement dilepton reconstruction
        # # dilepton
        dilep_system = dilepton.l1 + dilepton.l2
        req_dilepmass = ak.any(
            (dilep_system.mass > 81)
            & (dilep_system.mass < 101)
            & (dilep_system.pt > 15),
            axis=1,
        )

        # TODO: add jet selection
        overlap_l1 = ak.all(
            events.Jet.metric_table(dilepton.l1) > 0.4, axis=2, mask_identity=True
        )
        overlap_l2 = ak.all(
            events.Jet.metric_table(dilepton.l2) > 0.4, axis=2, mask_identity=True
        )
        jets = events.Jet[
            (events.Jet.pt > 20)
            & (abs(events.Jet.eta) < 2.5)
            & overlap_l1
            & overlap_l2
            & jet_id(events, self._campaign)
            # TODO: add more selection cuts (e.g. check utils/selection.py for jet_id(...))
        ]

        req_jets = ak.num(jets.pt) >= 1

        # event level selection
        event_level = ak.fill_none(
            req_lumi & req_trig & req_jets & req_dilep & req_dilepmass & req_metfilter,
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
        positive_lepton = ak.where(dilepton.l1.charge > 0, dilepton.l1, dilepton.l2)
        negative_lepton = ak.where(dilepton.l1.charge < 0, dilepton.l1, dilepton.l2)
        pruned_dilep = dilepton.l1[event_level] + dilepton.l2[event_level]
        smu = leptons[(dilepton.l1.flavor == 13) & event_level]
        sel = leptons[(dilepton.l1.flavor == 11) & event_level]
        # print("dilepton structure", dilepton)
        # print("dilepton structure", dilepton.l1)
        # print("dilepton structure", dilepton.l1.flavor)
        # print("leptons", leptons)
        # print("leptons", leptons.flavor)
        # Only store events which survive event selection ("pruned events")
        pruned_ev = events[event_level]
        print("events", len(events))
        print("pruned_events", len(pruned_ev))
        # print("pruned events field", pruned_ev.fields)
        # print("somthing to prove this file been executed")
        # jets
        pruned_ev["SelJet"] = jets[event_level][:, :1]
        pruned_ev["njet"] = ak.count(jets[event_level].pt, axis=1)
        pruned_ev["dilep"] = pruned_dilep
        pruned_ev["dilep", "pt"] = pruned_dilep.pt
        pruned_ev["dilep", "eta"] = pruned_dilep.eta
        pruned_ev["dilep", "phi"] = pruned_dilep.phi
        pruned_ev["dilep", "mass"] = pruned_dilep.mass
        pruned_ev["posl"] = positive_lepton[event_level]
        pruned_ev["negl"] = negative_lepton[event_level]
        # pruned_ev["SelMuon"] = smu
        # pruned_ev["SelMuon", "pt"] = smu.pt
        # pruned_ev["SelMuon", "eta"] = smu.eta
        # pruned_ev["SelMuon", "phi"] = smu.phi
        # pruned_ev["SelMuon", "charge"] = smu.charge
        # pruned_ev["SelElectron"] = sel
        # pruned_ev["SelElectron", "pt"] = sel.pt
        # pruned_ev["SelElectron", "eta"] = sel.eta
        # pruned_ev["SelElectron", "phi"] = sel.phi
        # pruned_ev["SelElectron", "charge"] = sel.charge

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
        # if "PFCands" in events.fields:
        #     pruned_ev["PFCands"] = PFCand_link(events, event_level, jetindx)

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
