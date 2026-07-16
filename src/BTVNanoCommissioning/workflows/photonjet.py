import awkward as ak
import os
import numpy as np
import correctionlib
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
        selectionModifier="",
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

    ## Apply corrections on momentum/mass on MET, Jet, Muon
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
        if self._year == "2022" or self._year == "2023":
            triggers = {
                "Photon20_HoverELoose": [20, 30],
                "Photon30EB_TightID_TightIso": [30, 50],
                "Photon50EB_TightID_TightIso": [50, 75],
                "Photon75EB_TightID_TightIso": [75, 90],
                "Photon90EB_TightID_TightIso": [90, 110],
                "Photon110EB_TightID_TightIso": [110, 200],
                "Photon200": [200, 9999],
            }
        elif self._year == "2024":
            triggers = {
                "Photon20_HoverELoose": [20, 30],
                "Photon30EB_TightID_TightIso": [30, 50],
                "Photon50EB_TightID_TightIso": [50, 110],
                "Photon110EB_TightID_TightIso": [110, 200],
                "Photon200": [200, 9999],
            }
        elif self._year == "2025":
            triggers = {
                "Photon20_HoverELoose": [20, 30],
                "Photon30EB_TightID_TightIso": [30, 40],
                "Photon40EB_TightID_TightIso": [40, 110],
                "Photon110EB_TightID_TightIso": [110, 200],
                "Photon200": [200, 9999],
            }
        else:
            raise ValueError(self._year, "is not a valid selection modifier.")

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

        ## Photon cuts
        photon_sel = (
            (events.Photon.cutBased == 3)
            & (events.Photon.hoe < 0.02148)
            & (events.Photon.r9 > 0.94)
            & (events.Photon.r9 < 1.0)
            & (np.abs(events.Photon.eta) < 1.3)
        )

        # Index with the selection rather than ak.mask: masking leaves the failing
        # candidates in place as None, so slot 0 remains the leading *unselected*
        # candidate. The photon is itself clustered as the leading jet in ~86% of
        # gamma+jet events, so the cuts below were being evaluated on that
        # photon-jet, which then dropped the event via None propagation.
        event_ph = ak.pad_none(events.Photon[photon_sel], 1)
        event_jet = ak.pad_none(events.Jet[jet_sel], 1)

        req_photon = ak.count(event_ph.pt, axis=1) > 0
        req_jet = ak.count(event_jet.pt, axis=1) > 0

        # Validate the paths against the sample (raises if none of them exist).
        # The returned OR is unused: each path is pt-binned individually below.
        HLT_helper(events, list(triggers.keys()))

        # NB: `events.HLT[trg] = ...` writes into a temporary copy of the HLT
        # record and is silently discarded, leaving the pt binning with no effect.
        # Keep the binned decisions in a local dict instead.
        trig_pass = {}
        for trg in triggers:
            if not hasattr(events.HLT, trg):
                continue
            trig_pass[trg] = ak.fill_none(
                events.HLT[trg]
                & (event_ph[:, 0].pt >= triggers[trg][0])
                & (event_ph[:, 0].pt < triggers[trg][1]),
                False,
            )

        req_trig = np.zeros(len(events), dtype="bool")
        for trg_pass in trig_pass.values():
            req_trig = req_trig | trg_pass

        req_dphi = np.abs(event_jet[:, 0].delta_phi(event_ph[:, 0])) > 2.7
        req_scale = np.abs(1.0 - event_jet[:, 0].pt / event_ph[:, 0].pt) < 0.3

        event_level = (
            event_level & req_photon & req_jet & req_dphi & req_scale & req_trig
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
                is_dijet=False,
            )

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

        # Take the leading candidate from the *selected* collections, so the
        # stored objects are the ones the cuts above were evaluated on.
        pruned_sel_jet = event_jet[event_level]
        pruned_ev["Tag"] = event_ph[event_level][:, 0]
        pruned_ev["Tag", "pt"] = pruned_ev["Tag"].pt
        pruned_ev["Tag", "eta"] = pruned_ev["Tag"].eta
        pruned_ev["Tag", "phi"] = pruned_ev["Tag"].phi
        pruned_ev["SelJet"] = pruned_sel_jet[:, 0]

        pruned_ev["njet"] = ak.count(pruned_sel_jet.pt, axis=1)

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
                raise ValueError(self._year, "is not supported for prescale weights.")

            pruned_ev["psweight"] = np.zeros(len(pruned_ev))
            for trigger in trig_pass:
                # Check if the prescale weight file exists for the given trigger and year
                psfile = f"src/BTVNanoCommissioning/data/Prescales/ps_weight_{trigger}_year{self._year}.json"
                if not os.path.isfile(psfile):
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
                # Use the pt-binned decision: the raw HLT bits overlap, so the
                # lowest (most prescaled) path would otherwise claim every event.
                pruned_ev["psweight"] = ak.where(
                    (trig_pass[trigger][event_level]) & (pruned_ev["psweight"] == 0),
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
            for trigger in trig_pass:
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
