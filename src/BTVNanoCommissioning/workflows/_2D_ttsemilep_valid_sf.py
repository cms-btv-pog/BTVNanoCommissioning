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
from BTVNanoCommissioning.utils.histogramming.histogrammer import (
    histogrammer,
    histo_writer,
)
from BTVNanoCommissioning.utils.array_writer import array_writer
from BTVNanoCommissioning.utils.selection import (
    HLT_helper,
    jet_id,
    # mu_idiso,
    mu_promptmvaid,
    # ele_mvatightid,
    ele_promptmvaid,
    MET_filters,
    btag_wp,
)


class NanoProcessor(processor.ProcessorABC):
    ## Define histograms
    def __init__(
        self,
        year="2024",
        campaign="Summer24",
        name="",
        isSyst=False,
        isArray=False,
        noHist=False,
        chunksize=75000,
        selectionModifier="semittE",
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
        ## Add selection for ttbar semileptonic
        self.selMod = selectionModifier

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

        objs = ["MET", "jet0", "jet1", "jet2", "jet3"]
        if self.selMod == "semittE":
            objs.append("ele")
            chn = "ele"
            triggers = ["Ele30_WPTight_Gsf"]
        elif self.selMod == "semittM":
            objs.append("mu")
            chn = "mu"
            triggers = ["IsoMu24"]
        else:
            raise ValueError(self.selMod, "is not a valid selection modifier.")

        ## Create histograms
        output = {}
        if not self.noHist:
            output = histogrammer(
                events.Jet.fields,
                obj_list=objs,
                hist_collections=["common", "fourvec", "ttsemilep_2D"],
                channel=chn,
                njet=4,
            )

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
        req_trig = HLT_helper(events, triggers)

        ## Lepton cuts
        if self.selMod == "semittE":
            event_iso_lep = events.Electron[
                (events.Electron.pt > 32) & ele_promptmvaid(events, self._campaign)
            ]
        elif self.selMod == "semittM":
            event_iso_lep = events.Muon[
                (events.Muon.pt > 26) & mu_promptmvaid(events, self._campaign)
            ]
            event_soft_mu = events.Muon[
                (events.Muon.pt > 5)
                & (abs(events.Muon.eta) < 2.5)
                & (events.Muon.tightId > 0.5)
                & (events.Muon.pfRelIso04_all > 0.25)
            ]
        req_lep = ak.count(event_iso_lep.pt, axis=1) == 1

        # DY -> mumu veto
        event_iso_lep = ak.pad_none(event_iso_lep, 1, axis=1)
        if self.selMod == "semittE":
            req_DYveto = np.ones(len(events), dtype="bool")
        elif self.selMod == "semittM":
            event_di_mu = event_soft_mu + event_iso_lep[:, 0]
            req_DYveto = (event_di_mu.mass > 12) & (
                (event_di_mu.mass < 81) | (event_di_mu.mass > 101)
            )
            req_DYveto = ak.all(req_DYveto, axis=1)

        ## Jet cuts
        jet_sel = ak.fill_none(
            jet_id(events, self._campaign, min_pt=25)
            & (
                ak.all(
                    events.Jet.metric_table(event_iso_lep) > 0.4,
                    axis=2,
                    mask_identity=True,
                )
            ),
            False,
            axis=-1,
        )
        event_jet = events.Jet[jet_sel]
        n_jet = ak.count(event_jet.pt, axis=1)
        # req_jets = (n_jet >= 3) & (n_jet <= 4)
        req_jets = n_jet >= 4

        # b-tagged jets requirement
        mask_bjets = btag_wp(
            event_jet, self._year, self._campaign, "UParTAK4", "b", "M"
        )
        mask_cjets = btag_wp(
            event_jet, self._year, self._campaign, "UParTAK4", "c", "M"
        )
        n_bjets = ak.count(event_jet[mask_bjets].pt, axis=1)
        n_cjets = ak.count(event_jet[mask_cjets].pt, axis=1)
        n_hfjets = n_bjets + n_cjets
        # req_b_jets = (n_bjets >= 1)
        req_b_jets = (n_bjets >= 1) & (n_hfjets >= 3)

        ## Other cuts

        ## MET cuts
        event_MET = ak.zip(
            {
                "pt": events.PuppiMET.pt,
                "eta": ak.zeros_like(events.PuppiMET.pt),
                "phi": events.PuppiMET.phi,
                "mass": ak.zeros_like(events.PuppiMET.pt),
            },
            with_name="PtEtaPhiMLorentzVector",
        )
        req_MET = event_MET.pt > 20.0
        req_metfilter = MET_filters(events, self._campaign)

        # Cut on tranverse W mass
        # Calculation taken from https://github.com/cms-btv-pog/BTVNanoCommissioning/blob/9f56298727c9c05df49b5701a0e92f8bf43d610e/src/BTVNanoCommissioning/workflows/ctag_Wctt_valid_sf.py#L299-L305
        event_dphi = event_iso_lep[:, 0].phi - events.PuppiMET.phi
        event_dphi = np.where(event_dphi < np.pi, event_dphi + 2 * np.pi, event_dphi)
        event_dphi = np.where(event_dphi > np.pi, event_dphi - 2 * np.pi, event_dphi)
        event_Wmt = np.sqrt(
            2 * event_iso_lep[:, 0].pt * events.PuppiMET.pt * (1 - np.cos(event_dphi))
        )
        req_Wmt = event_Wmt > 40

        event_level = ak.fill_none(
            req_lumi & req_trig & req_lep
            # & req_DYveto
            & req_jets
            # & req_b_jets
            & req_MET & req_metfilter,
            # & req_Wmt,
            False,
        )
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
        pruned_ev["SelJet"] = event_jet[event_level][:, :4]
        if self.selMod == "semittE":
            pruned_ev["SelElectron"] = event_iso_lep[event_level][:, 0]
        elif self.selMod == "semittM":
            pruned_ev["SelMuon"] = event_iso_lep[event_level][:, 0]
        pruned_ev["njet"] = ak.count(event_jet[event_level].pt, axis=1)
        b_jet_mask = btag_wp(
            event_jet[event_level], self._year, self._campaign, "UParTAK4", "b", "M"
        )
        c_jet_mask = btag_wp(
            event_jet[event_level], self._year, self._campaign, "UParTAK4", "c", "M"
        )
        pruned_ev["nbjet"] = ak.count(event_jet[event_level].pt[b_jet_mask], axis=1)
        pruned_ev["ncjet"] = ak.count(event_jet[event_level].pt[c_jet_mask], axis=1)
        pruned_ev["MET"] = event_MET[event_level]
        pruned_ev["w_mt"] = event_Wmt[event_level]

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
            output = histo_writer(
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
