import pickle, os, sys, mplhep as hep, numpy as np
import collections
from coffea import hist, processor
import awkward as ak
from coffea.analysis_tools import Weights
import gc

from BTVNanoCommissioning.utils.correction import (
    lumiMasks,
    load_pu,
    load_BTV,
    load_jetfactory,
    add_jec_variables,
)
from BTVNanoCommissioning.helpers.definitions import definitions
from BTVNanoCommissioning.utils.AK4_parameters import correction_config


class NanoProcessor(processor.ProcessorABC):
    # Define histograms
    def __init__(self, year="2017", campaign="Rereco17_94X"):
        self._year = year
        self._campaign = campaign
        # Define axes
        # Should read axes from NanoAOD config
        dataset_axis = hist.Cat("dataset", "Primary dataset")
        flav_axis = hist.Bin("flav", r"Genflavour", [0, 1, 4, 5, 6])
        cutflow_axis = hist.Cat("cut", "Cut")

        # Events
        nmu_axis = hist.Bin("nmu", r"N muons", [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        njet_axis = hist.Bin("njet", r"N jets", [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        nbjet_axis = hist.Bin("nbjet", r"N b-jets", [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        lmupt_axis = hist.Bin("lmupt", r"Muon pt", 45, 20, 200)

        # Jet
        jet_pt_axis = hist.Bin("pt", r"Jet $p_{T}$ [GeV]", 100, 0, 400)
        jet_rawpt_axis = hist.Bin("rawpt", r"Jet $p_{T}$ [GeV]", 100, 20, 400)
        jet_rawfactor_axis = hist.Bin("rawfactor", r"raw factor", 50, 0, 10)

        jet_eta_axis = hist.Bin("eta", r"Jet $\eta$", 60, -3, 3)
        jet_phi_axis = hist.Bin("phi", r"Jet $\phi$", 60, -3, 3)
        jet_mass_axis = hist.Bin("mass", r"Jet $m$ [GeV]", 100, 0, 50)
        jet_dr_axis = hist.Bin("dr", r"Jet $\Delta$R(l,j)", 50, 0, 5)
        ljpt_axis = hist.Bin("ljpt", r"Leading jet $p_{T}$ [GeV]", 100, 20, 400)
        sljpt_axis = hist.Bin("sljpt", r"Subleading jet $p_{T}$ [GeV]", 100, 20, 400)
        ljrawpt_axis = hist.Bin("ljrawpt", r"Leading jet $p_{T}$ [GeV]", 100, 20, 400)
        sljrawpt_axis = hist.Bin(
            "sljrawpt", r"Subleading jet $p_{T}$ [GeV]", 100, 20, 400
        )
        ljdr_axis = hist.Bin("ljdr", "Leading jet $\Delta$R(l,j)", 50, 0, 5)
        sljdr_axis = hist.Bin("sljdr", "Subleading jet $\Delta$R(l,j)", 50, 0, 5)

        # Define similar axes dynamically
        disc_list = [
            "btagCMVA",
            "btagCSVV2",
            "btagDeepB",
            "btagDeepC",
            "btagDeepFlavB",
            "btagDeepFlavC",
            "btagDeepCvB",
            "btagDeepCvL",
            "btagDeepFlavCvB",
            "btagDeepFlavCvL",
        ]
        discSF_list = [
            "btagDeepBSF",
            "btagDeepCSF",
            "btagDeepFlavBSF",
            "btagDeepFlavCSF",
        ]
        ddx_list = ["btagDDBvLV2", "btagDDCvBV2", "btagDDCvLV2"]
        btag_axes = []
        for d in disc_list:
            btag_axes.append(hist.Bin(d, d, 30, -0.2, 1))

        syst_axis = hist.Cat("syst", ["noSF", "SF", "SFup", "SFdn"])

        btag_axes = []
        for d in disc_list:
            btag_axes.append(hist.Bin(d, d, 30, -0.2, 1))
        btagDeepaxes = []
        bininfo = definitions()
        for d in bininfo.keys():
            ranges = bininfo[d]["manual_ranges"]
            binning = bininfo[d]["bins"]
            if ranges[1] is None:
                ranges[1] = 0.0
            if ranges[0] is None:
                ranges[0] = -0.5
            btagDeepaxes.append(hist.Bin(d, d, binning, ranges[0], ranges[1]))
        # Define histograms from axes
        _hist_jet_dict = {
            "pt": hist.Hist("Counts", dataset_axis, flav_axis, jet_pt_axis),
            "rawpt": hist.Hist("Counts", dataset_axis, flav_axis, jet_rawpt_axis),
            "eta": hist.Hist("Counts", dataset_axis, flav_axis, jet_eta_axis),
            "phi": hist.Hist("Counts", dataset_axis, flav_axis, jet_phi_axis),
            "mass": hist.Hist("Counts", dataset_axis, flav_axis, jet_mass_axis),
            "rawpt": hist.Hist(
                "Counts",
                dataset_axis,
                flav_axis,
                jet_pt_axis,
                jet_eta_axis,
                jet_rawpt_axis,
            ),
            "rawfactor": hist.Hist(
                "Counts",
                dataset_axis,
                flav_axis,
                jet_pt_axis,
                jet_eta_axis,
                jet_rawfactor_axis,
            )
            # 'dr': hist.Hist("Counts", dataset_axis, flav_axis,jet_dr_axis)
        }
        _hist_btagDeepdict = {
            "pt": hist.Hist("Counts", dataset_axis, flav_axis, jet_pt_axis),
            "rawpt": hist.Hist("Counts", dataset_axis, flav_axis, jet_rawpt_axis),
            "eta": hist.Hist("Counts", dataset_axis, flav_axis, jet_eta_axis),
            "phi": hist.Hist("Counts", dataset_axis, flav_axis, jet_phi_axis),
            "mass": hist.Hist("Counts", dataset_axis, flav_axis, jet_mass_axis),
            "rawpt": hist.Hist(
                "Counts",
                dataset_axis,
                flav_axis,
                jet_pt_axis,
                jet_eta_axis,
                jet_rawpt_axis,
            ),
            # 'rawfactor': hist.Hist("Counts",dataset_axis, flav_axis,jet_pt_axis, jet_eta_axis, jet_rawfactor_axis)
        }

        # Generate some histograms dynamically
        for disc, axis in zip(disc_list, btag_axes):
            _hist_btagDeepdict[disc] = hist.Hist(
                "Counts", dataset_axis, flav_axis, axis
            )
        for deepcsv, axises in zip(bininfo.keys(), btagDeepaxes):
            _hist_btagDeepdict[deepcsv] = hist.Hist(
                "Counts", dataset_axis, flav_axis, axises
            )
        _hist_event_dict = {
            "njet": hist.Hist("Counts", dataset_axis, njet_axis),
            "nbjet": hist.Hist("Counts", dataset_axis, nbjet_axis),
            "nmu": hist.Hist("Counts", dataset_axis, nmu_axis),
            "lmupt": hist.Hist("Counts", dataset_axis, lmupt_axis),
            "ljpt": hist.Hist("Counts", dataset_axis, flav_axis, ljpt_axis),
            "sljpt": hist.Hist("Counts", dataset_axis, flav_axis, sljpt_axis),
            "ljrawpt": hist.Hist("Counts", dataset_axis, flav_axis, ljrawpt_axis),
            "sljrawpt": hist.Hist("Counts", dataset_axis, flav_axis, sljrawpt_axis),
            "ljdr": hist.Hist("Counts", dataset_axis, flav_axis, ljdr_axis),
            "sljdr": hist.Hist("Counts", dataset_axis, flav_axis, sljdr_axis),
        }
        self.jet_hists = list(_hist_jet_dict.keys())
        self.btagDeephists = list(_hist_btagDeepdict.keys())
        # self.deepcsvSF_hists = list(_hist_deepcsvSF_dict.keys())
        # self.deepddx_hists = list(_hist_deepddx_dict.keys())
        self.event_hists = list(_hist_event_dict.keys())
        _hist_dict = {**_hist_btagDeepdict, **_hist_event_dict}
        self._accumulator = processor.dict_accumulator(_hist_dict)
        self._accumulator["sumw"] = processor.defaultdict_accumulator(float)
        ## Load corrections
        self._jet_factory = load_jetfactory(
            self._campaign, correction_config[self._campaign]["JME"]
        )
        (
            self._deepcsvb_sf,
            self._deepcsvc_sf,
            self._deepjetb_sf,
            self._deepjetc_sf,
        ) = load_BTV(self._campaign, correction_config[self._campaign]["BTV"])
        self._pu = load_pu(self._campaign, correction_config[self._campaign]["PU"])

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.accumulator.identity()
        dataset = events.metadata["dataset"]
        isRealData = not hasattr(events, "genWeight")

        if isRealData:
            output["sumw"][dataset] += 1.0
        else:
            output["sumw"][dataset] += ak.sum(events.genWeight)
        req_lumi = np.ones(len(events), dtype="bool")
        if isRealData:
            req_lumi = lumiMasks[self._year](events.run, events.luminosityBlock)
        weights = Weights(len(events), storeIndividual=True)
        ## Define the CvL, CvB
        if not hasattr(events, "btagDeepFlavCvL"):
            events.Jet["btagDeepFlavCvL"] = np.minimum(
                np.maximum(
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
                    1,
                ),
                -1,
            )
            events.Jet["btagDeepFlavCvB"] = np.minimum(
                np.maximum(
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
                    1,
                ),
                -1,
            )
            events.Jet["btagDeepCvL"] = np.minimum(
                np.maximum(
                    np.where(
                        (events.Jet.btagDeepC > 0) & (events.Jet.pt > 15),
                        (events.Jet.btagDeepC / (1.0 - events.Jet.btagDeepB)),
                        -1,
                    ),
                    1,
                ),
                -1,
            )
            events.Jet["btagDeepCvB"] = np.minimum(
                np.maximum(
                    np.where(
                        (events.Jet.btagDeepC > 0) & (events.Jet.pt > 15),
                        (
                            events.Jet.btagDeepC
                            / (events.Jet.btagDeepC + events.Jet.btagDeepB)
                        ),
                        -1,
                    ),
                    1,
                ),
                -1,
            )
        if not isRealData:
            weights.add("genweight", events.genWeight)
            weights.add(
                "puweight", self._pu[f"{self._year}_pileupweight"](events.Pileup.nPU)
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

        ## Muon cuts
        # muon twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
        events.Muon = events.Muon[
            (events.Muon.pt > 30) & (abs(events.Muon.eta < 2.4)) & (events.Muon.tightId)
        ]
        events.Muon = ak.pad_none(events.Muon, 1, axis=1)
        req_muon = ak.count(events.Muon.pt, axis=1) == 1

        # ## Electron cuts
        # # electron twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
        events.Electron = events.Electron[
            (events.Electron.pt > 30)
            & (abs(events.Electron.eta) < 2.4)
            & (events.Electron.cutBased > 3)
        ]
        events.Electron = ak.pad_none(events.Electron, 1, axis=1)
        req_ele = ak.count(events.Electron.pt, axis=1) == 1
        ## Jet cuts

        req_opposite_charge = (
            events.Electron[:, 0].charge * events.Muon[:, 0].charge
        ) == -1

        if not isRealData:
            corrected_jets = self._jet_factory["mc"].build(
                add_jec_variables(events.Jet, events.fixedGridRhoFastjetAll),
                lazy_cache=events.caches[0],
            )
        else:
            corrected_jets = self._jet_factory["data"].build(
                add_jec_variables(events.Jet, events.fixedGridRhoFastjetAll),
                lazy_cache=events.caches[0],
            )

        event_jet = corrected_jets[
            ((corrected_jets.pt * (1 - corrected_jets.rawFactor)) > 50)
            & (abs(corrected_jets.eta) <= 2.4)
            & (corrected_jets.puId > 0)
            & (corrected_jets.jetId > 5)
            & (corrected_jets.btagDeepB > 0.0)
            & (corrected_jets.btagDeepB < 1.0)
            & (corrected_jets.btagDeepC > 0.0)
            & (corrected_jets.btagDeepC < 1.0)
            & (corrected_jets.btagDeepFlavB > 0.0)
            & (corrected_jets.btagDeepFlavB < 1.0)
            & (corrected_jets.btagDeepFlavC > 0.0)
            & (corrected_jets.btagDeepFlavC < 1.0)
        ]
        req_jets = ak.num(event_jet.puId) >= 2

        event_level = (
            req_trig & req_lumi & req_muon & req_ele & req_jets & req_opposite_charge
        )
        if len(event_level) > 0:
            event_level = ak.fill_none(event_level, False)
        # Selected
        selev = events[event_level]

        #########

        # Per muon
        mu_eta = abs(selev.Muon.eta) < 2.4
        mu_pt = selev.Muon.pt > 30
        mu_idiso = (selev.Muon.tightId > 0.5) & (selev.Muon.pfRelIso04_all < 0.12)
        mu_level = mu_eta & mu_pt & mu_idiso

        smu = selev.Muon[mu_level]
        # Per jet : https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID
        jet_eta = abs(corrected_jets[event_level].eta) <= 2.4
        jet_pt = (
            corrected_jets[event_level].pt * (1 - corrected_jets[event_level].rawFactor)
        ) > 50
        jet_pu = (corrected_jets[event_level].puId > 0) & (
            corrected_jets[event_level].jetId > 5
        )
        jet_clean = (
            (corrected_jets[event_level].btagDeepB > 0.0)
            & (corrected_jets[event_level].btagDeepB < 1.0)
            & (corrected_jets[event_level].btagDeepC > 0.0)
            & (corrected_jets[event_level].btagDeepC < 1.0)
            & (corrected_jets[event_level].btagDeepFlavB > 0.0)
            & (corrected_jets[event_level].btagDeepFlavB < 1.0)
            & (corrected_jets[event_level].btagDeepFlavC > 0.0)
            & (corrected_jets[event_level].btagDeepFlavC < 1.0)
        )

        sjets = corrected_jets[event_level]

        # sbjets = selev.Jet[bjet_level]

        if isRealData:
            genflavor = ak.zeros_like(sjets.pt)
        else:
            par_flav = (sjets.partonFlavour == 0) & (sjets.hadronFlavour == 0)
            genflavor = sjets.hadronFlavour + 1 * par_flav
        ## Fill histograms dynamically
        for histname, h in output.items():
            if histname in self.btagDeephists:
                if isRealData:
                    if histname == "rawpt":
                        h.fill(
                            dataset=dataset,
                            flav=ak.flatten(genflavor),
                            pt=ak.flatten(sjets.pt),
                            eta=ak.flatten(sjets.eta),
                            rawpt=ak.flatten(sjets.pt * (1 - sjets.rawFactor)),
                        )
                    elif "btagDeep" in histname:
                        fields = {
                            l: np.where(sjets[l] < 0, -0.2, sjets[l])
                            for l in h.fields
                            if l in dir(sjets)
                        }
                        h.fill(dataset=dataset, flav=genflavor, **fields)

                else:
                    genweiev = ak.flatten(
                        ak.broadcast_arrays(weights.weight()[event_level], sjets["pt"])[
                            0
                        ]
                    )
                    if histname == "rawpt":
                        h.fill(
                            dataset=dataset,
                            flav=ak.flatten(genflavor),
                            pt=ak.flatten(sjets.pt),
                            eta=ak.flatten(sjets.eta),
                            rawpt=ak.flatten(sjets.pt * (1 - sjets.rawFactor)),
                            weight=genweiev,
                        )
                    elif "btagDeep" in histname:

                        fields = {
                            l: ak.flatten(np.where(sjets[l] < 0, -0.2, sjets[l]))
                            for l in h.fields
                            if l in dir(sjets)
                        }
                        h.fill(
                            dataset=dataset,
                            flav=ak.flatten(genflavor),
                            **fields,
                            weight=genweiev,
                        )

        def flatten(ar):  # flatten awkward into a 1d array to hist
            return ak.flatten(ar, axis=None)

        def num(ar):
            return ak.num(ak.fill_none(ar[~ak.is_none(ar)], 0), axis=0)

        if isRealData:
            output["njet"].fill(dataset=dataset, njet=flatten(ak.num(sjets)))
            # output['nmu'].fill(dataset=dataset,   nmu=flatten(ak.num(smu)))
            # output['lmupt'].fill(dataset=dataset, lmupt=flatten((smu.pt)))
            output["ljpt"].fill(dataset=dataset, flav=0, ljpt=flatten(sjets[:, 0].pt))
            output["sljpt"].fill(dataset=dataset, flav=0, sljpt=flatten(sjets[:, 1].pt))
            output["ljrawpt"].fill(
                dataset=dataset,
                flav=0,
                ljrawpt=flatten(sjets[:, 0].pt * (1 - sjets[:, 0].rawFactor)),
            )
            output["sljrawpt"].fill(
                dataset=dataset,
                flav=0,
                sljrawpt=flatten(sjets[:, 1].pt * (1 - sjets[:, 1].rawFactor)),
            )
        else:

            output["njet"].fill(
                dataset=dataset,
                njet=flatten(ak.num(sjets)),
                weight=weights.weight()[event_level],
            )
            # output['nmu'].fill(dataset=dataset,   nmu=flatten(ak.num(smu)),weight=weights.weight()[event_level])
            # output['lmupt'].fill(dataset=dataset, lmupt=flatten((smu.pt)),weight=weights.weight()[event_level])
            output["ljpt"].fill(
                dataset=dataset,
                flav=flatten(sjets[:, 0].hadronFlavour),
                ljpt=flatten(sjets[:, 0].pt),
                weight=weights.weight()[event_level],
            )
            output["sljpt"].fill(
                dataset=dataset,
                flav=flatten(sjets[:, 1].hadronFlavour),
                sljpt=flatten(sjets[:, 1].pt),
                weight=weights.weight()[event_level],
            )
            output["ljrawpt"].fill(
                dataset=dataset,
                flav=flatten(sjets[:, 0].hadronFlavour),
                ljrawpt=flatten(sjets[:, 0].pt * (1 - sjets[:, 0].rawFactor)),
                weight=weights.weight()[event_level],
            )
            output["sljrawpt"].fill(
                dataset=dataset,
                flav=flatten(sjets[:, 1].hadronFlavour),
                sljrawpt=flatten(sjets[:, 1].pt * (1 - sjets[:, 1].rawFactor)),
                weight=weights.weight()[event_level],
            )

        return output

    def postprocess(self, accumulator):
        return accumulator
