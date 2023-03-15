import collections, numpy as np, awkward as ak
import collections

from coffea import processor
from coffea.analysis_tools import Weights

from BTVNanoCommissioning.helpers.func import flatten, update
from BTVNanoCommissioning.helpers.definitions import definitions
from BTVNanoCommissioning.helpers.update_branch import missing_branch
from BTVNanoCommissioning.utils.correction import load_lumi
from BTVNanoCommissioning.utils.histogrammer import histogrammer
from BTVNanoCommissioning.utils.selection import jet_id, mu_idiso, ele_mvatightid


class NanoProcessor(processor.ProcessorABC):
    # Define histograms
    def __init__(self, year="2017", campaign="Rereco17_94X"):
        self._year = year
        self._campaign = campaign
        _hist_event_dict = histogrammer("validation")
        self.make_output = lambda: {
            "sumw": processor.defaultdict_accumulator(float),
            **_hist_event_dict,
        }
        self.lumiMask = load_lumi(self._campaign)

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.make_output()
        dataset = events.metadata["dataset"]
        isRealData = not hasattr(events, "genWeight")
        events = missing_branch(events)
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

        ## HLT
        triggers = ["IsoMu24"]
        checkHLT = ak.Array([hasattr(events.HLT, _trig) for _trig in triggers])
        if ak.all(checkHLT == False):
            raise ValueError("HLT paths:", triggers, " are all invalid in", dataset)
        elif ak.any(checkHLT == False):
            print(np.array(triggers)[~checkHLT], " not exist in", dataset)
        trig_arrs = [
            events.HLT[_trig] for _trig in triggers if hasattr(events.HLT, _trig)
        ]
        req_trig = np.zeros(len(events), dtype="bool")
        for t in trig_arrs:
            req_trig = req_trig | t

        ## Jet cuts
        event_jet = events.Jet[
            jet_id(events, self._campaign)
            & (ak.all(events.Jet.metric_table(events.Muon) > 0.4, axis=2))
            & (ak.all(events.Jet.metric_table(events.Electron) > 0.4, axis=2))
        ]
        req_jets = ak.num(event_jet.pt) >= 2

        ## store jet index for PFCands, create mask on the jet index
        jetindx = ak.mask(
            ak.local_index(events.Jet.pt),
            (
                jet_id(events, self._campaign)
                & (ak.all(events.Jet.metric_table(events.Muon) > 0.4, axis=2))
                & (ak.all(events.Jet.metric_table(events.Electron) > 0.4, axis=2))
            )
            == 1,
        )
        jetindx = ak.pad_none(jetindx, 2)
        jetindx = jetindx[:, :2]

        event_level = ak.fill_none(req_jets & req_lumi, False)
        weights = Weights(len(events[event_level]), storeIndividual=True)
        if not isRealData:
            weights.add("genweight", events[event_level].genWeight)
        ####################
        # Selected objects #
        ####################
        sjets = event_jet[event_level]
        sjets = sjets[:, :2]
        # Find the PFCands associate with selected jets. Search from jetindex->JetPFCands->PFCand
        if "PFCands" in events.fields:
            jetindx0 = jetindx[:, 0]
            jetindx1 = jetindx[:, 1]
            spfcands = collections.defaultdict(dict)
            spfcands[0] = events[event_level].PFCands[
                events[event_level]
                .JetPFCands[
                    events[event_level].JetPFCands.jetIdx == jetindx0[event_level]
                ]
                .pFCandsIdx
            ]
            spfcands[1] = events[event_level].PFCands[
                events[event_level]
                .JetPFCands[
                    events[event_level].JetPFCands.jetIdx == jetindx1[event_level]
                ]
                .pFCandsIdx
            ]

        ####################
        # Weight & Geninfo #
        ####################
        weights = Weights(len(events[event_level]), storeIndividual=True)
        if isRealData:
            genflavor = ak.zeros_like(sjets.pt)
        else:
            par_flav = (sjets.partonFlavour == 0) & (sjets.hadronFlavour == 0)
            genflavor = sjets.hadronFlavour + 1 * par_flav
            genweiev = ak.flatten(ak.broadcast_arrays(weights.weight(), sjets["pt"])[0])
        ####################
        #  Fill histogram  #
        ####################
        for histname, h in output.items():
            if (
                "Deep" in histname
                and "btag" not in histname
                and histname in events.Jet.fields
            ):
                h.fill(
                    flatten(genflavor),
                    flatten(sjets[histname]),
                    weight=flatten(
                        ak.broadcast_arrays(weights.weight(), sjets["pt"])[0]
                    ),
                )
            elif (
                "PFCands" in events.fields
                and "PFCands" in histname
                and histname.split("_")[1] in events.PFCands.fields
            ):
                for i in range(2):
                    h.fill(
                        flatten(
                            ak.broadcast_arrays(genflavor[:, i], spfcands[i]["pt"])[0]
                        ),
                        flatten(spfcands[i][histname.replace("PFCands_", "")]),
                        weight=flatten(
                            ak.broadcast_arrays(weights.weight(), spfcands[i]["pt"])[0]
                        ),
                    )
            elif "jet" in histname and histname.replace("jet0_", "") in sjets.fields:
                jet = sjets[:, 0]
                h.fill(
                    flatten(genflavor[:, 0]),
                    flatten(jet[histname.replace(f"jet0_", "")]),
                    weight=weights.weight(),
                )
            elif "jet" in histname and histname.replace("jet1_", "") in sjets.fields:
                jet = sjets[:, 1]
                h.fill(
                    flatten(genflavor[:, 1]),
                    flatten(jet[histname.replace(f"jet1_", "")]),
                    weight=weights.weight(),
                )
        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator
