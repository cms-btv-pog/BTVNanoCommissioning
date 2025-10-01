import awkward as ak
import numpy as np
import os
import uproot
from coffea import processor
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
    softmu_mask,
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
        vetoed_events, shifts = common_shifts(self, events)

        return processor.accumulate(
            self.process_shift(update(vetoed_events, collections), name)
            for collections, name in shifts
        )

    def process_shift(self, events, shift_name):
        dataset = events.metadata["dataset"]
        isRealData = not hasattr(events, "genWeight")
        output = {} if self.noHist else histogrammer(events, "emctag_ttdilep_sf")

        if shift_name is None:
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
        trigger_he = [
            "Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
            "Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
        ]
        trigger_hm = [
            "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
            "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
        ]

        req_trig_ele = HLT_helper(events, trigger_he)
        req_trig_mu = HLT_helper(events, trigger_hm)

        # Muon cuts
        iso_muon_mu = events.Muon[
            (events.Muon.pt > 24) & mu_idiso(events, self._campaign)
        ]
        iso_muon_ele = events.Muon[
            (events.Muon.pt > 14) & mu_idiso(events, self._campaign)
        ]

        # Electron cuts
        iso_ele_ele = events.Electron[
            (events.Electron.pt > 27) & ele_mvatightid(events, self._campaign)
        ]
        iso_ele_mu = events.Electron[
            (events.Electron.pt > 15) & ele_mvatightid(events, self._campaign)
        ]

        # cross leptons
        req_ele = (ak.count(iso_muon_ele.pt, axis=1) == 1) & (
            ak.count(iso_ele_ele.pt, axis=1) == 1
        )
        req_mu = (ak.count(iso_muon_mu.pt, axis=1) == 1) & (
            ak.count(iso_ele_mu.pt, axis=1) == 1
        )
        iso_ele = ak.concatenate([iso_ele_mu, iso_ele_ele], axis=1)
        iso_mu = ak.concatenate([iso_muon_mu, iso_muon_ele], axis=1)
        iso_ele = ak.pad_none(iso_ele, 1)
        iso_mu = ak.pad_none(iso_mu, 1)

        # Jet cuts
        req_ele_iso = ak.all(
            events.Jet.metric_table(iso_ele) > 0.4, axis=2, mask_identity=True
        )
        req_mu_iso = ak.all(
            events.Jet.metric_table(iso_mu) > 0.4, axis=2, mask_identity=True
        )
        jetsel = ak.fill_none(
            jet_id(events, self._campaign) & req_ele_iso & req_mu_iso,
            False,
            axis=-1,
        )
        event_jet = events.Jet[jetsel]
        req_jets = ak.count(event_jet.pt, axis=1) >= 2

        # Soft Muon cuts
        soft_muon = events.Muon[softmu_mask(events, self._campaign)]
        req_softmu = ak.count(soft_muon.pt, axis=1) >= 1

        # Muon jet cuts
        smu_iso = ak.all(
            events.Jet.metric_table(soft_muon) > 0.4, axis=2, mask_identity=True
        )
        mujetsel = ak.fill_none(
            smu_iso & ((events.Jet.muonIdx1 != -1) | (events.Jet.muonIdx2 != -1)),
            False,
            axis=-1,
        )
        mu_jet = events.Jet[mujetsel & jetsel]
        req_mujet = ak.count(mu_jet.pt, axis=1) >= 1

        # store jet index for PFCands, create mask on the jet index
        jetindx = ak.mask(ak.local_index(events.Jet.pt), mujetsel)
        jetindx = ak.pad_none(jetindx, 1)
        jetindx = jetindx[:, 0]

        # Other cuts
        req_dilepmass = ((iso_mu[:, 0] + iso_ele[:, 0]).mass > 12.0) & (
            ((iso_mu[:, 0] + iso_ele[:, 0]).mass < 75)
            | ((iso_mu[:, 0] + iso_ele[:, 0]).mass > 105)
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
        req_MET = MET.pt > 40

        event_level = (
            req_lumi
            & req_MET
            & req_jets
            & req_softmu
            & req_mujet
            & req_dilepmass
            & ((req_trig_ele & req_ele) | (req_trig_mu & req_mu))
        )
        event_level = ak.fill_none(event_level, False)
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
        # Keep the structure of events and pruned the object size
        pruned_ev = events[event_level]
        pruned_ev["SelJet"] = event_jet[event_level]
        pruned_ev["SelMuon"] = iso_mu[event_level][:, 0]
        pruned_ev["SelElectron"] = iso_ele[event_level][:, 0]
        pruned_ev["MuonJet"] = mu_jet[event_level][:, 0]
        pruned_ev["SoftMuon"] = soft_muon[event_level][:, 0]
        pruned_ev["dilep"] = pruned_ev.SelMuon + pruned_ev.SelElectron
        pruned_ev["dilep", "pt"] = pruned_ev.dilep.pt
        pruned_ev["dilep", "eta"] = pruned_ev.dilep.eta
        pruned_ev["dilep", "phi"] = pruned_ev.dilep.phi
        pruned_ev["dilep", "mass"] = pruned_ev.dilep.mass
        if "PFCands" in events.fields:
            pruned_ev["PFCands"] = PFCand_link(events, event_level, jetindx)
        # Add custom variables
        pruned_ev["njet"] = ak.count(event_jet[event_level].pt, axis=1)
        pruned_ev["dr_mujet_softmu"] = pruned_ev.SoftMuon.delta_r(pruned_ev.MuonJet)
        pruned_ev["dr_mujet_lep1"] = pruned_ev.SelMuon.delta_r(pruned_ev.MuonJet)
        pruned_ev["dr_mujet_lep2"] = pruned_ev.SelElectron.delta_r(pruned_ev.MuonJet)
        pruned_ev["dr_lep1_softmu"] = pruned_ev.SelMuon.delta_r(pruned_ev.SoftMuon)
        pruned_ev["soft_l_ptratio"] = pruned_ev.SoftMuon.pt / pruned_ev.MuonJet.pt
        pruned_ev["l1_ptratio"] = pruned_ev.SelMuon.pt / pruned_ev.MuonJet.pt
        pruned_ev["l2_ptratio"] = pruned_ev.SelElectron.pt / pruned_ev.MuonJet.pt
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
