import coffea
from coffea import processor
from coffea.analysis_tools import Weights
import numpy as np
import awkward as ak
from BTVNanoCommissioning.helpers.func import flatten, update

from BTVNanoCommissioning.utils.correction import lumiMasks
from BTVNanoCommissioning.utils.AK4_parameters import correction_config
from BTVNanoCommissioning.helpers.definitions import definitions
from BTVNanoCommissioning.utils.histogrammer import histogrammer


class NanoProcessor(processor.ProcessorABC):
    # Define histograms
    def __init__(self, year="2017", campaign="Rereco17_94X"):
        self._year = year
        self._campaign = campaign
        self._year = year
        self._campaign = campaign
        _hist_event_dict = histogrammer("ttcom")
        self.make_output = lambda: {
            "sumw": processor.defaultdict_accumulator(float),
            **_hist_event_dict,
        }

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.make_output()
        isRealData = not hasattr(events, "genWeight")

        dataset = events.metadata["dataset"]
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
            "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
            "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
        ]

        trig_arrs = [events.HLT[_trig.strip("HLT_")] for _trig in triggers]
        req_trig = np.zeros(len(events), dtype="bool")
        for t in trig_arrs:
            req_trig = req_trig | t

        ############
        # Event level

        ## Muon cuts
        # muon twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
        events.Muon = events.Muon[(events.Muon.pt > 30) & (abs(events.Muon.eta) <= 2.4)]
        req_muon = ak.count(events.Muon.pt, axis=1) == 1
        events.Muon = ak.pad_none(events.Muon, 1)

        ## Electron cuts
        # electron twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
        events.Electron = events.Electron[
            (events.Electron.pt > 30) & (abs(events.Electron.eta) <= 2.4)
        ]
        req_ele = ak.count(events.Electron.pt, axis=1) == 1
        events.Electron = ak.pad_none(events.Electron, 1)

        ## Jet cuts
        events.Jet = events.Jet[
            (abs(events.Jet.eta) <= 2.4)
            & (events.Jet.pt > 25)
            & (((events.Jet.puId > 6) & (events.Jet.pt < 50)) | (events.Jet.pt > 50))
            & (events.Jet.jetId >= 2)
        ]
        req_jets = ak.count(events.Jet.pt, axis=1) >= 2

        req_opposite_charge = (
            events.Electron[:, 0].charge * events.Muon[:, 0].charge
        ) == -1
        event_level = ak.fill_none(
            req_trig & req_jets & req_ele & req_muon & req_opposite_charge, False
        )

        selev = events[event_level]

        #########
        # Per electron
        el_eta = abs(selev.Electron.eta) <= 2.4
        el_pt = selev.Electron.pt > 30
        el_level = el_eta & el_pt

        # Per muon
        mu_eta = abs(selev.Muon.eta) <= 2.4
        mu_pt = selev.Muon.pt > 30
        mu_level = mu_eta & mu_pt

        # Per jet
        jet_eta = abs(selev.Jet.eta) <= 2.4
        jet_pt = selev.Jet.pt > 25
        jet_pu = ((selev.Jet.puId > 6) & (selev.Jet.pt < 50)) | (selev.Jet.pt > 50)
        jet_id = selev.Jet.jetId >= 2
        jet_level = jet_pu & jet_eta & jet_pt & jet_id

        sele = selev.Electron[el_level]
        sele = sele[:, 0]
        smu = selev.Muon[mu_level]
        smu = smu[:, 0]
        sjets = selev.Jet[jet_level]
        sjets = sjets[:, :2]

        if isRealData:
            genflavor = ak.zeros_like(sjets.pt)
        else:
            par_flav = (sjets.partonFlavour == 0) & (sjets.hadronFlavour == 0)
            genflavor = sjets.hadronFlavour + 1 * par_flav
            genweiev = ak.flatten(
                ak.broadcast_arrays(weights.weight()[event_level], sjets["pt"])[0]
            )
        # Fill histograms dynamically
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

            elif "mu" in histname and histname.replace("mu_", "") in smu.fields:

                h.fill(
                    flatten(smu[histname.replace("mu_", "")]),
                    weight=weights.weight()[event_level],
                )
            elif "ele" in histname and histname.replace("ele_", "") in sele.fields:
                h.fill(
                    flatten(sele[histname.replace("ele_", "")]),
                    weight=weights.weight()[event_level],
                )

        output["dr_mujet0"].fill(
            flatten(genflavor[:, 0]),
            flatten(sjets[:, 0].delta_r(smu)),
            weight=weights.weight()[event_level],
        )
        output["dr_mujet1"].fill(
            flatten(genflavor[:, 1]),
            flatten(sjets[:, 1].delta_r(smu)),
            weight=weights.weight()[event_level],
        )

        return output

    def postprocess(self, accumulator):
        return accumulator
