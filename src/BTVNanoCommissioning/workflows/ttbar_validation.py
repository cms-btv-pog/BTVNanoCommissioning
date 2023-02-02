import numpy as np, awkward as ak

from coffea import processor
from coffea.analysis_tools import Weights

from BTVNanoCommissioning.helpers.func import flatten, update
from BTVNanoCommissioning.utils.correction import load_lumi
from BTVNanoCommissioning.utils.AK4_parameters import correction_config
from BTVNanoCommissioning.helpers.definitions import definitions
from BTVNanoCommissioning.utils.histogrammer import histogrammer
from BTVNanoCommissioning.utils.selection import jet_id
from BTVNanoCommissioning.helpers.update_branch import missing_branch


class NanoProcessor(processor.ProcessorABC):
    # Define histograms
    def __init__(self, year="2017", campaign="Rereco17_94X"):
        self._year = year
        self._campaign = campaign
        self._year = year
        self.lumiMask = load_lumi(correction_config[self._campaign]["lumiMask"])
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
        weights = Weights(len(events), storeIndividual=True)
        if not isRealData:
            weights.add("genweight", events.genWeight)

        ## HLT
        triggers = [
            "Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
            "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
        ]
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
        events.Jet = events.Jet[jet_id(events, self._campaign)]
        req_jets = ak.count(events.Jet.pt, axis=1) >= 2

        ## Other cuts
        req_opposite_charge = (
            events.Electron[:, 0].charge * events.Muon[:, 0].charge
        ) == -1
        event_level = ak.fill_none(
            req_trig & req_jets & req_ele & req_muon & req_opposite_charge, False
        )

        ####################
        # Selected objects #
        ####################
        sele = events.Electron[event_level]
        sele = sele[:, 0]
        smu = events.Muon[event_level]
        smu = smu[:, 0]
        sjets = events.Jet[event_level]
        sjets = sjets[:, :2]

        ####################
        # Weight & Geninfo #
        ####################
        if isRealData:
            genflavor = ak.zeros_like(sjets.pt)
        else:
            par_flav = (sjets.partonFlavour == 0) & (sjets.hadronFlavour == 0)
            genflavor = sjets.hadronFlavour + 1 * par_flav
            genweiev = ak.flatten(
                ak.broadcast_arrays(weights.weight()[event_level], sjets["pt"])[0]
            )
        ####################
        #  Fill histogram  #
        ####################
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
