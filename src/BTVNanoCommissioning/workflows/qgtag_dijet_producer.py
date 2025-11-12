import awkward as ak
import numpy as np
import correctionlib
import os
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
from BTVNanoCommissioning.utils.histogramming.histogrammer import (
    histogrammer,
    histo_writter,
)
from BTVNanoCommissioning.utils.histogramming.histograms.qgtag import qg_writer
from BTVNanoCommissioning.utils.array_writer import array_writer
from BTVNanoCommissioning.utils.selection import (
    HLT_helper,
    jet_id,
    mu_idiso,
    ele_cuttightid,
    MET_filters,
)


class NanoProcessor(processor.ProcessorABC):
    def __init__(
        self,
        year="2022",
        campaign="Summer22",
        name="",
        isSyst=False,
        isArray=False,
        noHist=False,
        chunksize=75000,
        selectionModifier="DiPFJetAve",
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
        self.selectionModifier = selectionModifier

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

    ## Processed events per-chunk, made selections, filled histogram, stored root files
    def process_shift(self, events, shift_name):
        dataset = events.metadata["dataset"]
        isRealData = not hasattr(events, "genWeight")
        ######################
        #  Create histogram  # : Get the histogram dict from `histogrammer`
        ######################

        output = {}

        if not self.noHist:
            output = histogrammer(
                jet_fields=events.Jet.fields,
                obj_list=[],
                hist_collections=["qgtag"],
                axes_collections=["qgtag"],
                is_dijet=True,
            )

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
        if self.selectionModifier == "ZB":
            triggers = {
                "ZeroBias": [15, 80],
            }
            ptmin, ptmax = 15, 80
        elif self.selectionModifier == "PFJet500":
            triggers = {
                "PFJet500": [575, 1e7],
            }
            ptmin, ptmax = 575, 1e7
        elif self.selectionModifier == "PFJet":
            triggers = {
                "PFJet40": [60, 86],
                "PFJet60": [86, 110],
                "PFJet80": [110, 170],
                "PFJet140": [170, 236],
                "PFJet200": [236, 305],
                "PFJet260": [305, 375],
                "PFJet320": [375, 460],
                "PFJet400": [460, 575],
                # "PFJet450": [470, 530],
                "PFJet500": [575, 1e7],
                # "PFJet550": [600, 1e7],
            }
            # if int(self._year) > 2022:
                # triggers["PFJet80"] = [110, 140]
                # triggers["PFJet110"] = [140, 180]
                # triggers["PFJet140"] = [180, 220]
            ptmin, ptmax = 60, 1e7
        elif self.selectionModifier == "DiPFJetAve":
            # These act on the minimum pT of the two leading jets
            # Should test effect with average pT of the two leading jets
            triggers = {
                "DiPFJetAve40": [40, 60],
                "DiPFJetAve60": [60, 80],
                "DiPFJetAve80": [80, 140],
                "DiPFJetAve140": [140, 200],
                "DiPFJetAve200": [200, 260],
                "DiPFJetAve260": [260, 320],
                "DiPFJetAve320": [320, 400],
                "DiPFJetAve400": [400, 500],
                "DiPFJetAve500": [500, 1e7],
            }
            ptmin, ptmax = 40, 1e7
        elif self.selectionModifier == "DiPFJet_HF":
            triggers = {
                "DiPFJetAve60_HFJEC": [60, 80],
                "DiPFJetAve80_HFJEC": [80, 100],
                "DiPFJetAve100_HFJEC": [100, 160],
                "DiPFJetAve160_HFJEC": [160, 220],
                "DiPFJetAve220_HFJEC": [220, 300],
                "DiPFJetAve300_HFJEC": [300, 1e7],
            }
            if int(self._year) > 2022:
                triggers["DiPFJetAve220_HFJEC"] = [220, 260]
                triggers["DiPFJetAve260_HFJEC"] = [260, 300]
                triggers["DiPFJetAve300_HFJEC"] = [300, 1e7]
            ptmin, ptmax = 60, 1e7
        else:
            raise ValueError(
                self.selectionModifier, "is not a valid selection modifier."
            )

        # Redefine triggers to match the jet selection
        events.Jet = ak.pad_none(events.Jet, 2)
        if isRealData:
            for trg in triggers:
                events.HLT[trg] = (
                    events.HLT[trg]
                    & ((events.Jet[:, 0].pt >= triggers[trg][0]))
                    & ((events.Jet[:, 0].pt < triggers[trg][1]))
                )

        req_metfilter = MET_filters(events, self._campaign)

        event_level = req_lumi & req_metfilter

        ##### Add some selections
        ## Jet cuts
        jet_sel = jet_id(events, self._campaign, max_eta=4.7, min_pt=20)

        if self._year == "2016":
            jet_puid = events.Jet.puId >= 1
        elif self._year in ["2017", "2018"]:
            jet_puid = events.Jet.puId >= 4
        else:
            jet_puid = ak.ones_like(jet_sel)

        jet_sel = jet_sel & jet_puid
        event_jet = ak.mask(events.Jet, jet_sel)
        event_jet = ak.pad_none(event_jet, 3)
        events.Jet = ak.pad_none(
            events.Jet, 3
        )  # Make sure that the shape is consistent

        req_leadjet = event_jet[:, 0].pt < ptmax
        req_jet = ak.count(event_jet.pt, axis=1) > 1
        req_dphi = abs(event_jet[:, 0].delta_phi(event_jet[:, 1])) > 2.7
        req_subjet = ak.where(
            ak.count(event_jet.pt, axis=1) > 2,
            event_jet[:, 2].pt / (0.5 * (event_jet[:, 0] + event_jet[:, 1])).pt < 0.15,
            ak.ones_like(req_jet, dtype=bool),
        )
        req_bal = event_jet[:, 1].pt / event_jet[:, 0].pt > 0.7

        req_trig = HLT_helper(events, list(triggers.keys()))
        # req_trig = np.zeros_like(len(events), dtype=bool)
        # trigbools = {}
        # for trigger, pt_range in triggers.items():
            # rmin, rmax = pt_range
            # if "DiPFJetAve" in self.selectionModifier:
                # thistrigreq = (
                    # (HLT_helper(events, [trigger]))
                    # # Average
                    # # & (ak.fill_none(0.5 * (event_jet[:, 0].pt + event_jet[:, 1].pt) >= rmin, False))
                    # # & (ak.fill_none(0.5 * (event_jet[:, 0].pt + event_jet[:, 1].pt) < rmax, False))
                    # # Minimum
                    # & (ak.fill_none(event_jet[:, 1].pt >= rmin, False))
                    # & (ak.fill_none(event_jet[:, 1].pt < rmax, False))
                # )
                # trigbools[trigger] = thistrigreq
                # req_trig = req_trig | thistrigreq
            # else:
                # thistrigreq = (
                    # (HLT_helper(events, [trigger]))
                    # # & (ak.fill_none(event_jet[:, 0].pt >= rmin, False))
                    # # & (ak.fill_none(event_jet[:, 0].pt < rmax, False))
                    # & (event_jet[:, 0].pt >= rmin)
                    # & (event_jet[:, 0].pt < rmax)
                # )
                # trigbools[trigger] = thistrigreq
                # req_trig = req_trig | thistrigreq

        event_level = event_level & req_jet & req_dphi & req_subjet & req_trig & req_bal & req_leadjet

        ## MC only: require gen vertex to be close to reco vertex
        # if "GenVtx_z" in events.fields:
            # req_vtx = np.abs(events.GenVtx_z - events.PV_z) < 0.2
        # else:
            # req_vtx = ak.ones_like(events.run, dtype=bool)

        # event_level = event_level & req_vtx

        ##<==== finish selection
        event_level = ak.fill_none(event_level, False)
        # Skip empty events -
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

        ##===>  Ntuplization  : store custom information
        ####################
        # Selected objects # : Pruned objects with reduced event_level
        ####################
        # Keep the structure of events and pruned the object size
        pruned_ev = events[event_level]

        # Central jet, Forward jet, Random jet, Selected two jets
        pruned_ev["CenJet"] = ak.where(
            np.abs(pruned_ev.Jet[:, 0].eta) < np.abs(pruned_ev.Jet[:, 1].eta),
            pruned_ev.Jet[:, 0],
            pruned_ev.Jet[:, 1],
        )
        pruned_ev["FwdJet"] = ak.where(
            np.abs(pruned_ev.Jet[:, 0].eta) > np.abs(pruned_ev.Jet[:, 1].eta),
            pruned_ev.Jet[:, 0],
            pruned_ev.Jet[:, 1],
        )
        # pruned_ev["RndJet"] = pruned_ev.Jet[
        # :, np.random.randint(0, 2)#, size=len(pruned_ev))
        # ]
        pruned_ev["RndJet"] = ak.where(
            np.random.randint(0, 2, size=len(pruned_ev)) == 0,
            pruned_ev.Jet[:, 0],
            pruned_ev.Jet[:, 1],
        )
        pruned_ev["SelJet"] = pruned_ev.Jet[:, :2]
        pruned_ev["njet"] = ak.count(pruned_ev.Jet.pt, axis=1)

        ## <========= end: store custom objects

        ####################
        #     Output       #
        ####################
        # Configure SFs
        weights = weight_manager(pruned_ev, self.SF_map, self.isSyst)
        if isRealData:
            if self._year == "2022":
                run_num = "355374_362760"
            elif self._year == "2023":
                run_num = "366727_370790"
            elif self._year == "2024":
                run_num = "378985_386951"

            pruned_ev["psweight"] = np.zeros(len(pruned_ev))
            trglist = sorted(list(triggers.keys()), reverse=True, key=lambda x: triggers[x][0])
            
            for trigger in triggers:
                psfile = f"src/BTVNanoCommissioning/data/Prescales/ps_weight_{trigger}_run{run_num}.json"
                if not os.path.isfile(psfile):
                    raise NotImplementedError(
                        f"Prescale weights not available for {trigger} in {self._year}. Please run `scripts/dump_prescale.py`."
                    )
                pseval = correctionlib.CorrectionSet.from_file(psfile)
                thispsweight = pseval["prescaleWeight"].evaluate(
                    pruned_ev.run,
                    f"HLT_{trigger}",
                    ak.values_astype(pruned_ev.luminosityBlock, np.float32),
                )
                pruned_ev["psweight"] = ak.where(
                    (pruned_ev.HLT[trigger]) & (pruned_ev["psweight"] == 0), thispsweight, pruned_ev["psweight"]
                )
            weights.add("psweight", pruned_ev["psweight"])

        # Configure systematics
        if shift_name is None:
            systematics = ["nominal"] + list(weights.variations)
        else:
            systematics = [shift_name]

        # Configure histograms
        if not self.noHist:
            output = qg_writer(
                pruned_ev, output, weights, systematics, self.isSyst, self.SF_map
            )
        # Output arrays
        if self.isArray:
            othersData = [
                "SV_*",
                "PV_npvs",
                "PV_npvsGood",
                "Rho_*",
                "SoftMuon_dxySig",
                "Muon_sip3d",
                "run",
                "luminosityBlock",
            ]
            for trigger in triggers:
                othersData.append(f"HLT_{trigger}")
            array_writer(
                self,
                pruned_ev,
                events,
                weights,
                systematics,
                dataset,
                isRealData,
                othersData=othersData,
            )

        return {dataset: output}

    ## post process, return the accumulator, compressed
    def postprocess(self, accumulator):
        return accumulator
