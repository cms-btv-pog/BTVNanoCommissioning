import gzip
import pickle, os, sys, mplhep as hep, numpy as np, awkward as ak
import collections

import coffea
from coffea import processor
from coffea.analysis_tools import Weights
from BTVNanoCommissioning.helpers.func import flatten, update

from BTVNanoCommissioning.helpers.definitions import definitions
from BTVNanoCommissioning.utils.AK4_parameters import correction_config
from BTVNanoCommissioning.utils.correction import lumiMasks, load_pu
from BTVNanoCommissioning.utils.histogrammer import histogrammer


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
        self._pu = load_pu(self._campaign, correction_config[self._campaign]["PU"])

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.make_output()
        dataset = events.metadata["dataset"]
        isRealData = not hasattr(events, "genWeight")

        if isRealData:
            output["sumw"] = len(events)
        else:
            output["sumw"] = ak.sum(events.genWeight)

        req_lumi = np.ones(len(events), dtype="bool")
        if isRealData:
            req_lumi = lumiMasks[self._year](events.run, events.luminosityBlock)
        weights = Weights(len(events), storeIndividual=True)
        if not isRealData:
            weights.add("genweight", events.genWeight)
        if not hasattr(events, "btagDeepFlavCvL"):
            events.Jet["btagDeepFlavCvL"] = np.maximum(
                np.minimum(
                    np.where(
                        (
                            (
                                events.Jet.btagDeepFlavC
                                / (1.0 - events.Jet.btagDeepFlavB)
                            )
                            > 0
                        )
                        & (events.Jet.pt > 15),
                        (events.Jet.btagDeepFlavC / (1.0 - events.Jet.btagDeepFlavB)),
                        -1,
                    ),
                    0.999999,
                ),
                -1,
            )
            events.Jet["btagDeepFlavCvB"] = np.maximum(
                np.minimum(
                    np.where(
                        (
                            (
                                events.Jet.btagDeepFlavC
                                / (events.Jet.btagDeepFlavC + events.Jet.btagDeepFlavB)
                            )
                            > 0
                        )
                        & (events.Jet.pt > 15),
                        (
                            events.Jet.btagDeepFlavC
                            / (events.Jet.btagDeepFlavC + events.Jet.btagDeepFlavB)
                        ),
                        -1,
                    ),
                    0.999999,
                ),
                -1,
            )
            events.Jet["btagDeepCvL"] = np.maximum(
                np.minimum(
                    np.where(
                        (events.Jet.btagDeepC > 0) & (events.Jet.pt > 15),
                        (events.Jet.btagDeepC / (1.0 - events.Jet.btagDeepB)),
                        -1,
                    ),
                    0.999999,
                ),
                -1,
            )
            events.Jet["btagDeepCvB"] = np.maximum(
                np.minimum(
                    np.where(
                        (events.Jet.btagDeepC > 0) & (events.Jet.pt > 15),
                        (
                            events.Jet.btagDeepC
                            / (events.Jet.btagDeepC + events.Jet.btagDeepB)
                        ),
                        -1,
                    ),
                    0.999999,
                ),
                -1,
            )
        ##############
        # Trigger level
        triggers = [
            "HLT_IsoMu24",
        ]

        trig_arrs = [events.HLT[_trig.strip("HLT_")] for _trig in triggers]
        req_trig = np.zeros(len(events), dtype="bool")
        for t in trig_arrs:
            req_trig = req_trig | t

        ############
        # Event level
        event_jet = events.Jet[
            (events.Jet.pt > 25)
            & (abs(events.Jet.eta) <= 2.4)
            & (events.Jet.puId > 0)
            & (events.Jet.jetId > 5)
            & (ak.all(events.Jet.metric_table(events.Muon) > 0.4, axis=2))
            & (ak.all(events.Jet.metric_table(events.Electron) > 0.4, axis=2))
        ]
        req_jets = ak.num(event_jet.puId) >= 2
        event_level = req_jets
        if len(event_level) > 0:
            event_level = ak.fill_none(event_level, False)
        # Selected
        selev = events[event_level]
        #########
        # Per jet : https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID
        jet_eta = abs(selev.Jet.eta) <= 2.4
        jet_pt = selev.Jet.pt > 25
        jet_pu = (selev.Jet.puId > 0) & (selev.Jet.jetId > 5)
        jet_level = jet_pu & jet_eta & jet_pt
        sjets = selev.Jet[jet_level]
        sel_jets = sjets
        sjets = sjets[:, :2]

        if isRealData:
            genflavor = ak.zeros_like(sjets.pt)
        else:
            par_flav = (sjets.partonFlavour == 0) & (sjets.hadronFlavour == 0)
            genflavor = sjets.hadronFlavour + 1 * par_flav
            genweiev = ak.flatten(
                ak.broadcast_arrays(weights.weight()[event_level], sjets["pt"])[0]
            )

        for histname, h in output.items():
            if "Deep" in histname and "btag" not in histname:
                h.fill(
                    flatten(genflavor),
                    flatten(sjets[histname]),
                    weight=flatten(
                        ak.broadcast_arrays(weights.weight()[event_level], sjets["pt"])[
                            0
                        ]
                    ),
                )
            elif "jet" in histname and histname.replace("jet0_", "") in sjets.fields:
                jet = sjets[:, 0]
                h.fill(
                    flatten(genflavor[:, 0]),
                    flatten(jet[histname.replace(f"jet0_", "")]),
                    weight=weights.weight()[event_level],
                )
            elif "jet" in histname and histname.replace("jet1_", "") in sjets.fields:
                jet = sjets[:, 1]
                h.fill(
                    flatten(genflavor[:, 1]),
                    flatten(jet[histname.replace(f"jet1_", "")]),
                    weight=weights.weight()[event_level],
                )
        return output

    def postprocess(self, accumulator):
        return accumulator
