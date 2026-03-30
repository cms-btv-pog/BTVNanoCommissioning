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
    reweighting,
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
        sumws = reweighting(events, self.isSyst)
        vetoed_events, shifts = common_shifts(self, events)

        return processor.accumulate(
            self.process_shift(update(vetoed_events, collections), sumws, name)
            for collections, name in shifts
        )

    ## Processed events per-chunk, made selections, filled histogram, stored root files
    def process_shift(self, events, sumws, shift_name):
        dataset = events.metadata["dataset"]
        isRealData = not hasattr(events, "genWeight")

        ####################
        #    Selections    #
        ####################
        ## Lumimask
        req_lumi = np.ones(len(events), dtype="bool")
        if isRealData:
            req_lumi = self.lumiMask(events.run, events.luminosityBlock)

        ## HLT
        if self.selectionModifier == "ZB":
            triggers = {
                "ZeroBias": [15, 60],
            }
            ptmin, ptmax = 15, 80
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
                "PFJet500": [575, 1e7],
            }
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
        else:
            raise ValueError(
                self.selectionModifier, "is not a valid selection modifier."
            )

        # Redefine triggers to match the jet selection
        events.Jet = ak.pad_none(events.Jet, 2)

        # Sort the jets by pt
        events.Jet = events.Jet[ak.argsort(events.Jet.pt, axis=1, ascending=False)]

        for trg in triggers:
            events.HLT[trg] = (
                events.HLT[trg]
                & ((events.Jet[:, 0].pt >= triggers[trg][0]))
                & ((events.Jet[:, 0].pt < triggers[trg][1]))
            )

        req_trig = HLT_helper(events, list(triggers.keys()))

        req_metfilter = MET_filters(events, self._campaign)

        event_level = req_lumi & req_metfilter

        ##### Add some selections
        ## Jet cuts
        jet_sel = jet_id(events, self._campaign, max_eta=5.0, min_pt=20)

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
        req_bal = np.abs(1.0 - event_jet[:, 1].pt / event_jet[:, 0].pt) < 0.3

        event_level = (
            event_level
            & req_jet
            & req_dphi
            & req_subjet
            & req_trig
            & req_bal
            & req_leadjet
        )

        ## MC only: require gen vertex to be close to reco vertex
        if "GenVtx_z" in events.fields:
            req_vtx = np.abs(events.GenVtx_z - events.PV_z) < 0.2
        else:
            req_vtx = ak.ones_like(events.run, dtype=bool)

        event_level = event_level & req_vtx

        ##<==== finish selection

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
        # only dump for nominal case
        if shift_name is None:
            output = dump_lumi(events[req_lumi], output)

        if shift_name is None:
            output["sumw"] = sumws["sumw"]
            if not isRealData and self.isSyst:
                if "LHEPdfWeight" in events.fields:
                    output["PDF_sumwUp"] = sumws["PDF_sumwUp"]
                    output["PDF_sumwDown"] = sumws["PDF_sumwDown"]
                    output["aS_sumwUp"] = sumws["aS_sumwUp"]
                    output["aS_sumwDown"] = sumws["aS_sumwDown"]
                    output["PDFaS_sumwUp"] = sumws["PDFaS_sumwUp"]
                    output["PDFaS_sumwDown"] = sumws["PDFaS_sumwDown"]
                if "LHEScaleWeight" in events.fields:
                    output["muR_sumwUp"] = sumws["muR_sumwUp"]
                    output["muR_sumwDown"] = sumws["muR_sumwDown"]
                    output["muF_sumwUp"] = sumws["muF_sumwUp"]
                    output["muF_sumwDown"] = sumws["muF_sumwDown"]
                if "PSWeight" in events.fields:
                    if len(events.PSWeight[0]) == 4:
                        output["ISR_sumwUp"] = sumws["ISR_sumwUp"]
                        output["ISR_sumwDown"] = sumws["ISR_sumwDown"]
                        output["FSR_sumwUp"] = sumws["FSR_sumwUp"]
                        output["FSR_sumwDown"] = sumws["FSR_sumwDown"]

        event_level = ak.fill_none(event_level, False)
        if shift_name is None:
            output = dump_lumi(events[req_lumi], output)

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
        pruned_ev["RndJet"] = ak.where(
            np.random.randint(0, 2, size=len(pruned_ev)) == 0,
            pruned_ev.Jet[:, 0],
            pruned_ev.Jet[:, 1],
        )
        pruned_ev["LeadJet"] = ak.where(
            pruned_ev.Jet[:, 0].pt > pruned_ev.Jet[:, 1].pt,
            pruned_ev.Jet[:, 0],
            pruned_ev.Jet[:, 1],
        )
        pruned_ev["SubleadJet"] = ak.where(
            pruned_ev.Jet[:, 0].pt < pruned_ev.Jet[:, 1].pt,
            pruned_ev.Jet[:, 0],
            pruned_ev.Jet[:, 1],
        )
        pruned_ev["SubleadJet"] = ak.where(
            pruned_ev.Jet[:, 0].pt < pruned_ev.Jet[:, 1].pt,
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
            elif self._year == "2025":
                run_num = "391658_398860"
            else:
                raise NotImplementedError(
                    f"Prescale weights not available for data in {self._year}."
                )

            pruned_ev["psweight"] = np.zeros(len(pruned_ev))
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
                    (pruned_ev.HLT[trigger]) & (pruned_ev["psweight"] == 0),
                    thispsweight,
                    pruned_ev["psweight"],
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
