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
        selectionModifier="mu",  # possible choices: mu, el to select muon/electron channel
    ):
        self._year = year
        self._campaign = campaign
        self.name = name
        self.isSyst = isSyst
        self.isArray = isArray
        self.noHist = noHist
        self.lumiMask = load_lumi(self._campaign)
        self.chunksize = chunksize
        ### Added selection to select muon/electron channel
        if selectionModifier not in ["el", "mu"]:
            raise ValueError(f"Invalid selectionModifier: {selectionModifier}")
        self.channel = selectionModifier
        ## Load corrections
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
        output = {} if self.noHist else histogrammer(events, "sf_ttsemilep_tnp")

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
        if self.channel == "mu":
            # Iso trigger: IsoMu24_v.*
            # High pt trigger: Mu50_v.*
            triggers = ["IsoMu24", "Mu50"]
            isMu = True
        elif self.channel == "el":
            # Iso trigger: Ele32_WPTight_Gsf_v.*
            # High pt trigger: Ele50_CaloIdVT_GsfTrkIdT_PFJet165_v.*
            triggers = ["Ele32_WPTight_Gsf", "Ele50_CaloIdVT_GsfTrkIdT_PFJet165"]
        else:
            raise ValueError(self.channel, "is not a valid selection modifier.")
        req_trig = HLT_helper(events, triggers)

        ## Lepton selection
        veto_mu = events.Muon[
            (events.Muon.pt > 15) & (abs(events.Muon.eta) < 2.4) & (mu_iso(events, self._campaign))
        ]
        veto_el = events.Electron[
            (events.Electron.pt > 15) & (abs(events.Muon.eta) < 2.4) & (mu_iso(events, self._campaign))
        ]
        if self.channel == "mu":
            # muon twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
            # select muons with pt > 30 GeV, abs(eta) < 2.4, tight ID
            mu_mask = (
                (events.Muon.pt > 30)
                & (abs(events.Muon.eta) < 2.4)
                & mu_idiso(events, self._campaign)  # per-muon boolean mask expected
            )
            iso_lep = events.Muon[mu_mask]
            req_lepveto = (ak.num(veto_mu, axis=1) == 1) & (ak.num(veto_el, axis=1) == 0)

        elif self.channel == "el":
            # select electrons with pt > 30 GeV, abs(eta) < 2.4, tight ID
            el_mask = (
                (events.Electron.pt > 30)
                & (abs(events.Electron.eta) < 2.4)
                & ele_mvatightid(events, self._campaign)  # per-electron boolean mask expected
            )
            iso_lep = events.Electron[el_mask]
            req_lepveto = (ak.num(veto_mu, axis=1) == 0) & (ak.num(veto_el, axis=1) == 1)

        # exactly one selected lepton
        req_lep = ak.num(iso_lep, axis=1) == 1


        ## Jet selection
        event_jet = events.Jet[
            ak.fill_none(
                jet_id(events, self._campaign)
                & (
                    ak.all(
                        events.Jet.metric_table(events.Muon) > 0.4,
                        axis=2,
                        mask_identity=True,
                    )
                ) & (
                    ak.all(
                        events.Jet.metric_table(events.Electron) > 0.4,
                        axis=2,
                        mask_identity=True,
                    )
                ),
                False,
                axis=-1,
            )
        ]
        req_jets = ak.num(event_jet.pt) == 4

        # at least one medium tagged jet ("tag")

        # ... add stuff here ...

        # binnings
        # bins_toppt = {0.0, 40.0, 80.0, 120.0, 160.0, 200.0, 250.0, 300.0, 350.0, 400.0, 500., 1500};
        # bins_topy = {0.0, 0.2, 0.4, 0.6, 0.8, 1., 1.2, 1.4, 1.6, 1.8, 2.0, 2.5};
        # bins_ttm = {250., 380., 460., 540., 620., 700., 800., 900., 1000, 1200, 1400, 1600, 3500};
        # bins_ttpt = {0.0, 50.0, 100.0, 150., 200.0, 300.0, 400.0, 500., 600., 800.0, 1200.0};
        # bins_tty = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.4};
        # bins_ttphi = {0.0, 40., 80., 120., 140., 150., 160., 170., 180.};




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
