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
from BTVNanoCommissioning.utils.histogrammer import histogrammer, histo_writter
from BTVNanoCommissioning.utils.array_writer import array_writer
from BTVNanoCommissioning.utils.selection import (
    HLT_helper,
    jet_id,
    mu_idiso,
    ele_cuttightid,
    MET_filters,
)

import hist as Hist
hists = {
    "jet0_pt": Hist.Hist(
        Hist.axis.StrCategory([], name="syst", growth=True),
        Hist.axis.IntCategory([0, 1, 4, 5, 6], name="flav", label="Genflavour"),
        Hist.axis.Regular(100, 0, 200, name="pt", underflow=False, overflow=False),
        Hist.storage.Weight(),
    ),
    "njet": Hist.Hist(
        Hist.axis.StrCategory([], name="syst", growth=True),
        Hist.axis.Regular(10, 0, 10, name="njet", underflow=False, overflow=False),
        Hist.storage.Weight(),
    ),
}


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
        output = (
            {}
            if self.noHist
            else hists
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
        if self._year == "2022" or self._year == "2023":
            triggers = {
                "Photon20_HoverELoose": [20, 30],
                "Photon30EB_TightID_TightIso": [30, 50],
                "Photon50": [50, 75],
                "Photon75": [75, 90],
                "Photon90": [90, 110],
                "Photon110EB_TightID_TightIso": [110, 200],
                "Photon200": [200, 9999],
            }
        elif self._year == "2024" or self._year == "2025":
            triggers = {
                "Photon30EB_TightID_TightIso": [30, 40],
                "Photon40EB_TightID_TightIso": [40, 50],
                "Photon50EB_TightID_TightIso": [50, 55],
                "Photon55EB_TightID_TightIso": [55, 75],
                "Photon75EB_TightID_TightIso": [75, 90],
                "Photon90EB_TightID_TightIso": [90, 110],
                "Photon110EB_TightID_TightIso": [110, 200],
                "Photon200": [200, 9999],
            }
        else:
            raise ValueError(
                self._year, "is not a valid selection modifier."
            )

        req_metfilter = MET_filters(events, self._campaign)

        event_level = req_lumi & req_metfilter

        ##### Add some selections
        ## Jet cuts
        jet_sel = (
            (events.Jet.pt >= 15)
            & (abs(events.Jet.eta) < 5.13)
            & (events.Jet.jetId >= 4)
        )

        ## Photon cuts
        photon_sel = (
            (events.Photon.cutBased == 3)
            & (events.Photon.hoe < 0.02148)
            & (events.Photon.r9 > 0.94)
            & (events.Photon.r9 < 1.0)
        )

        event_ph = ak.mask(events.Photon, photon_sel)
        event_ph = ak.pad_none(event_ph, 1)
        events.Photon = ak.pad_none(events.Photon, 1)

        req_photon = ak.count(event_ph.pt, axis=1) > 0

        req_trig = np.zeros(len(events), dtype="bool")
        trigbools = {}
        for trigger, ptrange in triggers.items():
            ptmin = ptrange[0]
            ptmax = ptrange[1]
            # Require *leading photon* to be in the pT range of the trigger
            thistrigreq = (
                (HLT_helper(events, [trigger]))
                & (ak.fill_none(ak.firsts(event_ph.pt) >= ptmin, False))
                & (ak.fill_none(ak.firsts(event_ph.pt) < ptmax, False))
            )
            trigbools[trigger] = thistrigreq
            req_trig = (req_trig) | (thistrigreq)

        if self._year == "2016":
            jet_puid = events.Jet.puId >= 1
        elif self._year in ["2017", "2018"]:
            jet_puid = events.Jet.puId >= 4
        else:
            jet_puid = ak.ones_like(jet_sel)

        jet_sel = jet_sel & jet_puid
        event_jet = ak.mask(events.Jet, jet_sel)
        event_jet = ak.pad_none(event_jet, 1)
        events.Jet = ak.pad_none(
            events.Jet, 1
        )  # Make sure that the shape is consistent

        req_jet = ak.count(event_jet.pt, axis=1) > 0
        req_dphi = np.abs(event_jet[:, 0].delta_phi(event_ph[:, 0])) > 2.7
        req_scale = np.abs(1.0 - event_jet[:, 0].pt / event_ph[:, 0].pt) < 0.7

        event_level = event_level & req_jet & req_dphi & req_scale & req_trig

        ## MC only: require gen vertex to be close to reco vertex
        if "GenVtx_z" in events.fields:
            req_vtx = np.abs(events.GenVtx_z - events.PV_z) < 0.2
        else:
            req_vtx = ak.ones_like(events.run, dtype=bool)

        event_level = event_level & req_vtx

        ##<==== finish selection
        event_level = ak.fill_none(event_level, False)
        # Skip empty events -
        if len(events[event_level]) == 0:
            if self.isArray:
                array_writer(
                    self,
                    events[event_level],
                    events,
                    "nominal",
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

        pruned_ev["Tag"] = pruned_ev.Photon[:, 0]
        pruned_ev["Tag", "pt"] = pruned_ev["Tag"].pt
        pruned_ev["Tag", "eta"] = pruned_ev["Tag"].eta
        pruned_ev["Tag", "phi"] = pruned_ev["Tag"].phi
        pruned_ev["Tag", "mass"] = pruned_ev["Tag"].mass
        pruned_ev["SelJet"] = pruned_ev.Jet[:, 0]

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

            psweight = np.zeros(len(pruned_ev))
            for trigger, trigbool in trigbools.items():
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
                psweight = ak.where(trigbool[event_level], thispsweight, psweight)
            weights.add("psweight", psweight)

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
