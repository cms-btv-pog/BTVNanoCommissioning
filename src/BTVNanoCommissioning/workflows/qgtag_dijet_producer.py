import awkward as ak
import numpy as np
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

    ## Apply corrections on momentum/mass on MET, Jet, Muon
    def process(self, events):
        events = missing_branch(events)
        shifts = common_shifts(self, events)

        return processor.accumulate(
            self.process_shift(update(events, collections), name)
            for collections, name in shifts
        )

    ## Processed events per-chunk, made selections, filled histogram, stored root files
    def process_shift(self, events, shift_name):
        dataset = events.metadata["dataset"]
        isRealData = not hasattr(events, "genWeight")
        ######################
        #  Create histogram  # : Get the histogram dict from `histogrammer`
        ######################
        _hist_event_dict = (
            {}
            if self.noHist
            else histogrammer(events, "qgtag_dijet")  # this is the place to modify
        )

        output = {
            "sumw": processor.defaultdict_accumulator(float),
            **_hist_event_dict,
        }
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
            triggers = [
                "ZeroBias",
            ]
        elif self.selectionModifier == "PFJet":
            triggers = [
                "PFJet40",
                "PFJet60",
                "PFJet80",
                "PFJet110",
                "PFJet140",
                "PFJet200",
                "PFJet260",
                "PFJet320",
                "PFJet400",
                "PFJet450",
                "PFJet500",
                "PFJet550",
            ]
        elif self.selectionModifier == "DiPFJetAve":
            triggers = [
                "DiPFJetAve40",
                "DiPFJetAve60",
                "DiPFJetAve80",
                "DiPFJetAve140",
                "DiPFJetAve200",
                "DiPFJetAve260",
                "DiPFJetAve320",
                "DiPFJetAve400",
                "DiPFJetAve500",
            ]
        else:
            raise ValueError(
                self.selectionModifier, "is not a valid selection modifier."
            )

        req_trig = HLT_helper(events, triggers)
        req_metfilter = MET_filters(events, self._campaign)

        event_level = req_trig & req_lumi & req_metfilter

        ##### Add some selections
        ## Jet cuts
        jet_sel = (
            (events.Jet.pt >= 15)
            & (abs(events.Jet.eta) < 5.31)
            & (events.Jet.jetId >= 4)
        )

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

        req_jet = ak.count(event_jet.pt, axis=1) > 1
        req_dphi = abs(event_jet[:, 0].delta_phi(event_jet[:, 1])) > 2.7
        req_subjet = ak.where(
            ak.count(event_jet.pt, axis=1) > 2,
            event_jet[:, 2].pt / (0.5 * (event_jet[:, 0] + event_jet[:, 1])).pt < 1.0,
            ak.ones_like(req_jet),
        )
        req_tagjet = ak.where(
            req_jet,
            ak.fill_none(np.abs(event_jet.eta[:, 0]) < 1.3, False)
            | ak.fill_none(np.abs(event_jet.eta[:, 1]) < 1.3, False),
            ak.zeros_like(req_jet),
        )

        event_level = event_level & req_jet & req_dphi & req_subjet # & req_tagjet

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

        # pruned_ev["Tag"] = ak.where(
            # ak.fill_none(np.abs(pruned_ev.Jet.eta)[:, 0] < 1.3, False)
            # & ak.fill_none(np.abs(pruned_ev.Jet.eta)[:, 1] < 1.3, False),
            # pruned_ev.Jet[
                # :, np.random.randint(0, 2)
            # ],  # Is this correct? Is a single value returned for the whole program or for each event?
            # ak.where(
                # np.abs(pruned_ev.Jet.eta)[:, 0] < 1.3,
                # pruned_ev.Jet[:, 0],
                # pruned_ev.Jet[:, 1],
            # ),
        # )
        pruned_ev["Tag"] = pruned_ev.Jet[:, np.random.randint(0, 2)]
        pruned_ev["Tag", "pt"] = pruned_ev["Tag"].pt
        pruned_ev["Tag", "eta"] = pruned_ev["Tag"].eta
        pruned_ev["Tag", "phi"] = pruned_ev["Tag"].phi
        pruned_ev["Tag", "mass"] = pruned_ev["Tag"].mass
        pruned_ev["SelJet"] = ak.where(
            pruned_ev.Jet[:, 0].pt == pruned_ev.Tag.pt,
            pruned_ev.Jet[:, 1],
            pruned_ev.Jet[:, 0],
        )
        pruned_ev["njet"] = ak.count(pruned_ev.Jet.pt, axis=1)

        ## <========= end: store custom objects

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

        print({dataset: output})

        return {dataset: output}

    ## post process, return the accumulator, compressed
    def postprocess(self, accumulator):
        return accumulator
