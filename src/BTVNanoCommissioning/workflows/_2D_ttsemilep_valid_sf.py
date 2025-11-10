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
from BTVNanoCommissioning.utils.histogramming.histogrammer import histogrammer, histo_writer
from BTVNanoCommissioning.utils.array_writer import array_writer
from BTVNanoCommissioning.utils.selection import (
    HLT_helper,
    jet_id,
    mu_idiso,
    ele_mvatightid,
    MET_filters,
    calculate_new_discriminators,
    get_wp_2D,
    btag_wp_dict,
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
        selectionModifier="semittM",
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
            triggers = ["Ele30_WPTight_Gsf"] # so far only 2024 triggers
        elif self.selMod == "semittM":
            objs.append("mu")
            chn = "mu"
            triggers = ["IsoMu24"] # so far only 2024 triggers
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
                include_discriminators_2D=True,
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
                (events.Electron.pt > 32) & ele_mvatightid(events, self._campaign)
            ]
        elif self.selMod == "semittM":
            # muon twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
            event_iso_lep = events.Muon[
                (events.Muon.pt > 30) & mu_idiso(events, self._campaign)
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
            req_DYveto = (event_di_mu.mass > 12) & ((event_di_mu.mass < 81) | (event_di_mu.mass > 101))
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
        req_jets = (ak.num(event_jet.pt) >= 3) & (ak.num(event_jet.pt) <= 4)
        event_jet = ak.pad_none(event_jet, 4, axis=1)

        # b-tagged jets requirement
        WP = btag_wp_dict[self._year + "_" + self._campaign]["UParTAK4"]["2D"]
        btagUParTAK4HFvLF1, btagUParTAK4BvC1 = calculate_new_discriminators(event_jet[:, 0])
        btagUParTAK4HFvLF2, btagUParTAK4BvC2 = calculate_new_discriminators(event_jet[:, 1])
        wp2D_1 = ak.Array([get_wp_2D(btagUParTAK4HFvLF1[i], btagUParTAK4BvC1[i], self._year, self._campaign, "UParTAK4") for i in range(len(btagUParTAK4HFvLF1))])
        wp2D_2 = ak.Array([get_wp_2D(btagUParTAK4HFvLF2[i], btagUParTAK4BvC2[i], self._year, self._campaign, "UParTAK4") for i in range(len(btagUParTAK4HFvLF2))])
        req_b_jets = ((wp2D_1 >= WP["mapping"]["B2"]) & (wp2D_1 <= WP["mapping"]["B4"])) | ((wp2D_2 >= WP["mapping"]["B2"]) & (wp2D_2 <= WP["mapping"]["B4"]))

        ## Store jet index for PFCands, create mask on the jet index
        jetindx = ak.mask(ak.local_index(events.Jet.pt), jet_sel==True)
        jetindx = ak.pad_none(jetindx, 1)
        jetindx = jetindx[:, 0]

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
        req_MET = event_MET.pt > 50
        req_metfilter = MET_filters(events, self._campaign)

        ## Cut on tranverse W mass
        # Calculation taken from https://github.com/cms-btv-pog/BTVNanoCommissioning/blob/9f56298727c9c05df49b5701a0e92f8bf43d610e/src/BTVNanoCommissioning/workflows/ctag_Wctt_valid_sf.py#L299-L305
        event_dphi = event_iso_lep[:, 0].phi - events.PuppiMET.phi
        event_dphi = np.where(event_dphi < np.pi, event_dphi + 2 * np.pi, event_dphi)
        event_dphi = np.where(event_dphi > np.pi, event_dphi - 2 * np.pi, event_dphi)
        event_Wmt = np.sqrt(2 * event_iso_lep[:, 0].pt * events.PuppiMET.pt * (1 - np.cos(event_dphi)))
        req_Wmt = event_Wmt > 40

        event_level = ak.fill_none(
            req_lumi
            & req_trig
            & req_lep
            & req_DYveto
            & req_jets
            & req_b_jets
            & req_MET
            & req_metfilter
            & req_Wmt,
            False
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
        if "PFCands" in events.fields:
            pruned_ev.PFCands = PFCand_link(events, event_level, jetindx)
        pruned_ev["MET"] = event_MET[event_level]
        pruned_ev["w_mt"] = event_Wmt[event_level]

        pruned_ev[f"nbjet"] = ak.zeros_like(pruned_ev["SelJet"][:, 0].pt)
        pruned_ev[f"ncjet"] = ak.zeros_like(pruned_ev["SelJet"][:, 0].pt)
        for i in range(4):
            ith_jets = pruned_ev.SelJet[:, i]
            btagUParTAK4HFvLF, btagUParTAK4BvC = calculate_new_discriminators(ith_jets)
            wp2D = ak.Array([get_wp_2D(btagUParTAK4HFvLF[i], btagUParTAK4BvC[i], self._year, self._campaign, "UParTAK4") for i in range(len(btagUParTAK4HFvLF))])
            nbjet = np.where(((wp2D >= WP["mapping"]["B2"]) & (wp2D <= WP["mapping"]["B4"])), 1, 0)
            ncjet = np.where(((wp2D >= WP["mapping"]["C1"]) & (wp2D <= WP["mapping"]["C4"])), 1, 0)
            if i == 3:
                nbjet = ak.fill_none(nbjet, 0)
                ncjet = ak.fill_none(ncjet, 0)
            pruned_ev[f"nbjet"] = pruned_ev[f"nbjet"] + nbjet
            pruned_ev[f"ncjet"] = pruned_ev[f"ncjet"] + ncjet
            pruned_ev[f"btagUParTAK4HFvLF_{i}"] = btagUParTAK4HFvLF
            pruned_ev[f"btagUParTAK4BvC_{i}"] = btagUParTAK4BvC
            pruned_ev[f"btagUParTAK4HFvLFt_{i}"] = ak.Array(np.where(btagUParTAK4HFvLF >= 0.0, 1.0 - (1.0 - btagUParTAK4HFvLF)**0.5, -1.0))
            pruned_ev[f"btagUParTAK4BvCt_{i}"] = ak.Array(np.where(btagUParTAK4BvC >= 0.0, 1.0 - (1.0 - btagUParTAK4BvC)**0.5, -1.0))
            pruned_ev[f"btagUParTAK42D_{i}"] = wp2D
            jet_pt_bins = WP["jet_pt_bins"]
            for jet_pt_bin in jet_pt_bins:
                 pruned_ev[f"btagUParTAK42D_pt{jet_pt_bin[0]}to{jet_pt_bin[1]}_{i}"] = [wp2D[ijet] if pt is not None and jet_pt_bin[0] < pt and pt < jet_pt_bin[1] else None for ijet, pt in enumerate(ith_jets.pt)]
        pruned_ev.nbjet = ak.values_astype(pruned_ev.nbjet, "int64")
        pruned_ev.ncjet = ak.values_astype(pruned_ev.ncjet, "int64")

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
