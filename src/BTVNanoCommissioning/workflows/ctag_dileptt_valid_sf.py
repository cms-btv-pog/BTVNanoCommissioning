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
        selectionModifier="DilepTTEM",
    ):
        self._year = year
        self._campaign = campaign
        self.name = name
        self.isSyst = isSyst
        self.isArray = isArray
        self.noHist = noHist
        self.lumiMask = load_lumi(self._campaign)
        self.chunksize = chunksize
        self.selMod = selectionModifier
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

        ####################
        #    Selections    #
        ####################
        # Lumimask
        req_lumi = np.ones(len(events), dtype="bool")
        if isRealData:
            req_lumi = self.lumiMask(events.run, events.luminosityBlock)

        # HLT
        if self.selMod == "dilepttM":
            triggers = ["Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8"]
        elif self.selMod == "dilepttE":
            triggers = ["Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"]
        req_trig = HLT_helper(events, triggers)

        # Lepton selections
        if self.selMod == "dilepttM":
            # dilepton selections
            iso_muon = events.Muon[
                (events.Muon.pt > 12) & mu_idiso(events, self._campaign)
            ]
            iso_lep = ak.pad_none(iso_muon, 2)
            req_lep = (ak.count(iso_lep.pt, axis=1) == 2) & (iso_lep[:, 0].pt > 20)
            # veto other flavors
            dilep_ele = events.Electron[
                (events.Electron.pt > 15) & ele_mvatightid(events, self._campaign)
            ]
            req_dilepveto = ak.count(dilep_ele.pt, axis=1) == 0
        elif self.selMod == "dilepttE":
            # dilepton selections
            iso_ele = events.Electron[
                (events.Electron.pt > 25) & ele_mvatightid(events, self._campaign)
            ]
            iso_lep = ak.pad_none(iso_ele, 2)
            req_lep = (ak.count(iso_lep.pt, axis=1) == 2) & (iso_lep[:, 0].pt > 27)
            # veto other flavors
            dilep_mu = events.Muon[
                (events.Muon.pt > 12) & mu_idiso(events, self._campaign)
            ]
            req_dilepveto = ak.count(dilep_mu.pt, axis=1) == 0

        # veto Z events
        dilep_mass = iso_lep[:, 0] + iso_lep[:, 1]
        req_dilepmass = (dilep_mass.mass > 12.0) & (
            (dilep_mass.mass < 75) | (dilep_mass.mass > 105)
        )
        # Jet cuts
        lep_iso = ak.all(
            events.Jet.metric_table(iso_lep) > 0.4, axis=2, mask_identity=True
        )
        event_jet = events.Jet[
            ak.fill_none(jet_id(events, self._campaign) & lep_iso, False, axis=-1)
        ]
        req_jets = ak.count(event_jet.pt, axis=1) >= 2

        # Soft Muon cuts
        soft_muon = events.Muon[softmu_mask(events, self._campaign)]
        req_softmu = ak.count(soft_muon.pt, axis=1) >= 1

        # Muon jet cuts
        mu_jet = events.Jet[
            ak.fill_none(
                (
                    ak.all(
                        event_jet.metric_table(soft_muon) <= 0.4,
                        axis=2,
                        mask_identity=True,
                    )
                )
                & ((event_jet.muonIdx1 != -1) | (event_jet.muonIdx2 != -1)),
                False,
                axis=-1,
            )
        ]

        req_mujet = ak.count(mu_jet.pt, axis=1) >= 1

        # store jet index for PFCands, create mask on the jet index
        softmu_iso = ak.all(
            events.Jet.metric_table(soft_muon) <= 0.4, axis=2, mask_identity=True
        )
        jetindx = ak.mask(
            ak.local_index(events.Jet.pt),
            (
                jet_id(events, self._campaign)
                & softmu_iso
                & ((events.Jet.muonIdx1 != -1) | (events.Jet.muonIdx2 != -1))
            )
            == 1,
        )
        jetindx = ak.pad_none(jetindx, 1)
        jetindx = jetindx[:, 0]

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
            req_trig
            & req_lumi
            & req_lep
            & req_dilepveto
            & req_dilepmass
            & req_MET
            & req_jets
            & req_softmu
            & req_mujet
        )
        event_level = ak.fill_none(event_level, False)
        histname = {
            "dilepttM": "ctag_ttdilep_sf",
            "dilepttE": "ectag_ttdilep_sf",
        }
        output = {} if self.noHist else histogrammer(events, histname[self.selMod])

        if shift_name is None:
            if isRealData:
                output["sumw"] = len(events)
            else:
                output["sumw"] = ak.sum(events.genWeight)
        # only dump for nominal case
        if shift_name is None:
            output = dump_lumi(events[req_lumi], output)
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
        shlep = iso_lep[event_level]

        # Keep the structure of events and pruned the object size
        pruned_ev = events[event_level]
        pruned_ev["SelJet"] = event_jet[event_level]

        if self.selMod == "dilepttM":
            pruned_ev["SelMuon"] = shlep[:, :2]
            pruned_ev["hl"] = shlep[:, 0]
            pruned_ev["sl"] = shlep[:, 1]
            pruned_ev["dilep"] = shlep[:, 0] + shlep[:, 1]
        if self.selMod == "dilepttE":
            pruned_ev["SelElectron"] = shlep[:, :2]
            pruned_ev["hl"] = shlep[:, 0]
            pruned_ev["sl"] = shlep[:, 1]
            pruned_ev["dilep"] = shlep[:, 0] + shlep[:, 1]

        pruned_ev["MuonJet"] = mu_jet[event_level][:, 0]
        pruned_ev["SoftMuon"] = soft_muon[event_level][:, 0]
        pruned_ev["njet"] = ak.count(event_jet[event_level].pt, axis=1)
        pruned_ev["dilep", "pt"] = pruned_ev.dilep.pt
        pruned_ev["dilep", "eta"] = pruned_ev.dilep.eta
        pruned_ev["dilep", "phi"] = pruned_ev.dilep.phi
        pruned_ev["dilep", "mass"] = pruned_ev.dilep.mass
        if "PFCands" in events.fields:
            pruned_ev.PFCands = PFCand_link(events, event_level, jetindx)

        pruned_ev["dr_mujet_softmu"] = pruned_ev.SoftMuon.delta_r(pruned_ev.MuonJet)
        pruned_ev["soft_l_ptratio"] = pruned_ev.SoftMuon.pt / pruned_ev.MuonJet.pt

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
