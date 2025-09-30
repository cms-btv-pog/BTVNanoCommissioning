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
            triggers = {
                "ZeroBias": [15, 45],
            }
        elif self.selectionModifier == "PFJet":
            triggers = {
                "PFJet40": [45, 80],
                "PFJet60": [80, 110],
                "PFJet80": [110, 160],
                "PFJet110": [160, 180],
                "PFJet140": [180, 220],
                "PFJet200": [220, 300],
                "PFJet260": [300, 340],
                "PFJet320": [340, 420],
                "PFJet400": [420, 470],
                "PFJet450": [470, 530],
                "PFJet500": [530, 600],
                "PFJet550": [600, 1e7],
            }
        elif self.selectionModifier == "ZBpPFJet":
            triggers = {
                "ZeroBias": [15, 45],
                "PFJet40": [45, 80],
                "PFJet60": [80, 110],
                "PFJet80": [110, 160],
                "PFJet110": [160, 180],
                "PFJet140": [180, 220],
                "PFJet200": [220, 300],
                "PFJet260": [300, 340],
                "PFJet320": [340, 420],
                "PFJet400": [420, 470],
                "PFJet450": [470, 530],
                "PFJet500": [530, 600],
                "PFJet550": [600, 1e7],
            }
        elif self.selectionModifier == "DiPFJetAve":
            # These act on the average pT of the two leading jets
            # Should test effect with minimum pT of the two leading jets
            triggers = {
                "HLT_DiPFJetAve40": [45, 80],
                "HLT_DiPFJetAve60": [80, 110],
                "HLT_DiPFJetAve80": [110, 170],
                "HLT_DiPFJetAve140": [170, 230],
                "HLT_DiPFJetAve200": [230, 295],
                "HLT_DiPFJetAve260": [295, 360],
                "HLT_DiPFJetAve320": [360, 440],
                "HLT_DiPFJetAve400": [440, 540],
                "HLT_DiPFJetAve500": [540, 1e7],
            }
        elif self.selectionModifier == "DiPFJet_HF":
            triggers = {
                "HLT_DiPFJetAve60_HFJEC": [80, 100],
                "HLT_DiPFJetAve80_HFJEC": [100, 110],
                "HLT_DiPFJetAve100_HFJEC": [110, 180],
                "HLT_DiPFJetAve160_HFJEC": [180, 235],
                "HLT_DiPFJetAve220_HFJEC": [235, 280],
                "HLT_DiPFJetAve260_HFJEC": [280, 320],
                "HLT_DiPFJetAve300_HFJEC": [320, 295],
            }

        else:
            raise ValueError(
                self.selectionModifier, "is not a valid selection modifier."
            )

        # req_trig = HLT_helper(events, triggers)

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

        req_trig = np.zeros_like(len(events), dtype=bool)
        trigbools = {}
        for trigger, pt_range in triggers.items():
            rmin, rmax = pt_range
            if "DiPFJetAve" in self.selectionModifier:
                thistrigreq = (
                    (HLT_helper(events, [trigger])) 
                    & (ak.fill_none(0.5 * (event_jet[:, 0].pt + event_jet[:, 1].pt) >= rmin, False))
                    & (ak.fill_none(0.5 * (event_jet[:, 0].pt + event_jet[:, 1].pt) < rmax, False))
                )
                trigbools[trigger] = thistrigreq
                req_trig = req_trig | thistrigreq


        event_level = event_level & req_jet & req_dphi & req_subjet & req_trig

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

        # Central jet, Forward jet, Random jet, Selected two jets
        pruned_ev["CntJet"] = ak.where(
            np.abs(pruned_ev.Jet[:, 0].eta) < np.abs(pruned_ev.Jet[:, 1].eta),
            pruned_ev.Jet[:, 0],
            pruned_ev.Jet[:, 1],
        )
        pruned_ev["FwdJet"] = ak.where(
            np.abs(pruned_ev.Jet[:, 0].eta) > np.abs(pruned_ev.Jet[:, 1].eta),
            pruned_ev.Jet[:, 0],
            pruned_ev.Jet[:, 1],
        )
        pruned_ev["RndJet"] = pruned_ev.Jet[:, np.random.randint(0, 2, size=len(pruned_ev))]
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

            # if 369869 in pruned_ev.run: continue
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

        print({dataset: output})

        return {dataset: output}

    ## post process, return the accumulator, compressed
    def postprocess(self, accumulator):
        return accumulator
