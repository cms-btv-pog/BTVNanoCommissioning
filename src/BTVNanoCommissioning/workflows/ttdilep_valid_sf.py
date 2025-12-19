import collections, awkward as ak, numpy as np
import os
import uproot
from coffea import processor
from coffea.analysis_tools import Weights

# functions to load SFs, corrections
from BTVNanoCommissioning.utils.correction import (
    load_lumi,
    load_SF,
    weight_manager,
    common_shifts,
)

# user helper function
from BTVNanoCommissioning.helpers.func import flatten, update, dump_lumi, PFCand_link
from BTVNanoCommissioning.helpers.update_branch import missing_branch

## load histograms & selctions for this workflow
from BTVNanoCommissioning.utils.histogramming.histogrammer import (
    histogrammer,
    histo_writer,
)
from BTVNanoCommissioning.utils.array_writer import array_writer
from BTVNanoCommissioning.utils.selection import (
    HLT_helper,
    jet_id,
    mu_idiso,
    ele_cuttightid,
    ele_mvatightid,
    ele_promptmvaid,
    MET_filters,
    calculate_new_discriminators,
    get_wp_2D,
    btag_wp_dict,
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
        selectionModifier="tt_dilep",
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
        ## Load corrections
        self.SF_map = load_SF(self._year, self._campaign)

    @property
    def accumulator(self):
        return self._accumulator

    ## Apply corrections on momentum/mass on MET, Jet, Muon
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
        ## Create histograms
        objs = ["mu", "ele", "jet0", "jet1"]
        if self.selMod == "ttdilep_sf_2D":
            objs.append("MET")
        output = {}
        if not self.noHist:
            output = histogrammer(
                events.Jet.fields,
                obj_list=objs,
                hist_collections=["common", "fourvec", "ttdilep"],
                njet=2,
                include_discriminators_2D=True if "2D" in self.selMod else False,
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
        triggers = [
            "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
            "Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
            "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
        ]
        req_trig = HLT_helper(events, triggers)

        ## Muon cuts
        # muon twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
        events.Muon = events.Muon[
            (events.Muon.pt > 30) & mu_idiso(events, self._campaign)
        ]
        events.Muon = ak.pad_none(events.Muon, 1, axis=1)
        req_muon = ak.count(events.Muon.pt, axis=1) == 1

        ## Electron cuts
        # electron twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
        if self.selMod == "ttdilep_sf_2D":
            events.Electron = events.Electron[
                (events.Electron.pt > 30) & ele_mvatightid(events, self._campaign)
            ]
        else:
            events.Electron = events.Electron[
                (events.Electron.pt > 30) & ele_cuttightid(events, self._campaign)
            ]
        events.Electron = ak.pad_none(events.Electron, 1, axis=1)
        req_ele = ak.count(events.Electron.pt, axis=1) == 1

        ## Jet cuts
        jetsel = ak.fill_none(
            jet_id(events, self._campaign, min_pt=25)
            & (
                ak.all(
                    events.Jet.metric_table(events.Muon) > 0.4,
                    axis=2,
                    mask_identity=True,
                )
            )
            & (
                ak.all(
                    events.Jet.metric_table(events.Electron) > 0.4,
                    axis=2,
                    mask_identity=True,
                )
            ),
            False,
        )
        event_jet = events.Jet[jetsel]
        req_jets = ak.num(event_jet.pt) >= 2

        ## Other cuts
        req_opposite_charge = (
            events.Electron[:, 0:1].charge * events.Muon[:, 0:1].charge
        ) == -1
        req_opposite_charge = ak.fill_none(req_opposite_charge, False)
        req_opposite_charge = ak.flatten(req_opposite_charge)

        ## store jet index for PFCands, create mask on the jet index
        jetindx = ak.mask(ak.local_index(events.Jet.pt), jetsel)
        jetindx = ak.pad_none(jetindx, 2)
        jetindx = jetindx[:, :2]

        ## MET
        if self.selMod == "ttdilep_sf_2D":
            event_MET = ak.zip(
                {
                    "pt": events.PuppiMET.pt,
                    "eta": ak.zeros_like(events.PuppiMET.pt),
                    "phi": events.PuppiMET.phi,
                    "mass": ak.zeros_like(events.PuppiMET.pt),
                },
                with_name="PtEtaPhiMLorentzVector",
            )
        req_metfilter = MET_filters(events, self._campaign)

        event_level = (
            req_trig & req_lumi & req_muon & req_ele & req_jets & req_opposite_charge & req_metfilter
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
        pruned_ev["SelJet"] = event_jet[event_level][:, :2]
        pruned_ev["SelMuon"] = events.Muon[event_level][:, 0]
        pruned_ev["SelElectron"] = events.Electron[event_level][:, 0]
        pruned_ev["dr_mujet0"] = pruned_ev.Muon.delta_r(pruned_ev.Jet[:, 0])
        pruned_ev["dr_mujet1"] = pruned_ev.Muon.delta_r(pruned_ev.Jet[:, 1])
        pruned_ev["njet"] = ak.count(event_jet[event_level].pt, axis=1)
        # Find the PFCands associate with selected jets. Search from jetindex->JetPFCands->PFCand
        if "PFCands" in events.fields:
            pruned_ev.PFCands = PFCand_link(events, event_level, jetindx)

        if self.selMod == "ttdilep_sf_2D":
            nj = 2
            pruned_ev["MET"] = event_MET[event_level]
            for i in range(nj):
                btagUParTAK4HFvLF, btagUParTAK4BvC = calculate_new_discriminators(
                    pruned_ev.SelJet[:, i]
                )
                wp2D = ak.Array(
                    [
                        get_wp_2D(
                            btagUParTAK4HFvLF[i],
                            btagUParTAK4BvC[i],
                            self._year,
                            self._campaign,
                            "UParTAK4",
                        )
                        for i in range(len(btagUParTAK4HFvLF))
                    ]
                )
                pruned_ev[f"btagUParTAK4HFvLF_{i}"] = btagUParTAK4HFvLF
                pruned_ev[f"btagUParTAK4BvC_{i}"] = btagUParTAK4BvC
                pruned_ev[f"btagUParTAK4HFvLFt_{i}"] = ak.Array(
                    np.where(
                        btagUParTAK4HFvLF > 0.0,
                        1.0 - (1.0 - btagUParTAK4HFvLF) ** 0.5,
                        -1.0,
                    )
                )
                pruned_ev[f"btagUParTAK4BvCt_{i}"] = ak.Array(
                    np.where(
                        btagUParTAK4BvC > 0.0,
                        1.0 - (1.0 - btagUParTAK4BvC) ** 0.5,
                        -1.0,
                    )
                )
                pruned_ev[f"btagUParTAK42D_{i}"] = wp2D
                jet_pt_bins = btag_wp_dict[self._year + "_" + self._campaign][
                    "UParTAK4"
                ]["2D"]["jet_pt_bins"]
                for jet_pt_bin in jet_pt_bins:
                    pruned_ev[
                        f"btagUParTAK42D_pt{jet_pt_bin[0]}to{jet_pt_bin[1]}_{i}"
                    ] = [
                        (
                            wp2D[ijet]
                            if pt is not None
                            and jet_pt_bin[0] < pt
                            and pt < jet_pt_bin[1]
                            else None
                        )
                        for ijet, pt in enumerate(pruned_ev.SelJet[:, i].pt)
                    ]

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
