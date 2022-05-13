import pickle, os, sys, mplhep as hep, numpy as np
import collections

from matplotlib.pyplot import jet

import coffea
from coffea import hist, processor
import awkward as ak
import gc

from BTVNanoCommissioning.utils.AK4_parameters import correction_config
from BTVNanoCommissioning.utils.correction import lumiMasks, muSFs, load_pu, load_BTV
from BTVNanoCommissioning.helpers.definitions import definitions
from coffea.analysis_tools import Weights
from BTVNanoCommissioning.helpers.cTagSFReader import getSF


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
        ntbjet_axis = hist.Bin(
            "ntbjet", r"N b-jets", [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        )
        lmupt_axis = hist.Bin("lmupt", r"Muon pt", 40, 0, 200)

        # Jet
        jet_pt_axis = hist.Bin("pt", r"Jet $p_{T}$ [GeV]", 50, 0, 500)
        jet_ptwide_axis = hist.Bin(
            "ptwide",
            r"Jet $p_{T}$ [GeV]",
            [25, 30, 40, 60, 80, 100, 150, 200, 300, 500],
        )
        jet_rawpt_axis = hist.Bin("rawpt", r"Jet $p_{T}$ [GeV]", 50, 0, 500)
        jet_eta_axis = hist.Bin("eta", r"Jet $\eta$", 25, -2.5, 2.5)
        jet_etawide_axis = hist.Bin(
            "etawide", r"Jet $\eta$", [-2.5, -2.0, -1.5, -0.5, 0.0, 0.5, 1.5, 2.0, 2.5]
        )
        jet_phi_axis = hist.Bin("phi", r"Jet $\phi$", 30, -3, 3)
        jet_mass_axis = hist.Bin("mass", r"Jet $m$ [GeV]", 50, 0, 500)
        jet_dr_axis = hist.Bin("dr", r"Jet $\Delta$R(l,j)", 20, 0, 5)
        ljpt_axis = hist.Bin("ljpt", r"Leading jet $p_{T}$ [GeV]", 50, 0, 500)
        sljpt_axis = hist.Bin("sljpt", r"Subleading jet $p_{T}$ [GeV]", 50, 0, 500)
        ljeta_axis = hist.Bin("ljeta", r"Leading Jet $\eta$", 25, -2.5, 2.5)
        sljeta_axis = hist.Bin("sljeta", r"Leading Jet $\eta$", 25, -2.5, 2.5)

        ljdr_axis = hist.Bin("ljdr", "Leading jet $\Delta$R(l,j)", 20, 0, 5)
        sljdr_axis = hist.Bin("sljdr", "Subleading jet $\Delta$R(l,j)", 20, 0, 5)
        syst_axis = hist.Cat("syst", ["noSF", "SF", "SFup", "SFdn"])

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
        # Define similar axes dynamically
        disc_list = [
            "btagDeepB",
            "btagDeepC",
            "btagDeepFlavB",
            "btagDeepFlavC",
            "btagDeepCvL",
            "btagDeepCvB",
            "btagDeepFlavCvL",
            "btagDeepFlavCvB",
        ]
        btag_axes = []
        for d in disc_list:
            btag_axes.append(hist.Bin("%s" % (d), "%s" % (d), 30, -0.2, 1))

        _hist_sf_dict = {}
        _hist_btagDeepdict = {}
        for disc, axis in zip(disc_list, btag_axes):
            for i in range(2):
                _hist_sf_dict["%s_%d" % (disc, i)] = hist.Hist(
                    "Counts", dataset_axis, flav_axis, syst_axis, axis
                )
        for deepcsv, axises in zip(bininfo.keys(), btagDeepaxes):
            _hist_btagDeepdict["%s" % (deepcsv)] = hist.Hist(
                "Counts", dataset_axis, flav_axis, axises
            )

        _hist_event_dict = {
            "njet": hist.Hist("Counts", dataset_axis, njet_axis),
            "nbjet": hist.Hist("Counts", dataset_axis, nbjet_axis),
            "nbjet_up": hist.Hist("Counts", dataset_axis, nbjet_axis),
            "nbjet_dn": hist.Hist("Counts", dataset_axis, nbjet_axis),
            "ntbjet": hist.Hist("Counts", dataset_axis, ntbjet_axis),
            "nmu": hist.Hist("Counts", dataset_axis, nmu_axis),
            "lmupt": hist.Hist("Counts", dataset_axis, lmupt_axis),
            "ljpt": hist.Hist("Counts", dataset_axis, flav_axis, ljpt_axis),
            "sljpt": hist.Hist("Counts", dataset_axis, flav_axis, sljpt_axis),
            "ljeta": hist.Hist("Counts", dataset_axis, flav_axis, ljeta_axis),
            "sljeta": hist.Hist("Counts", dataset_axis, flav_axis, sljeta_axis),
            "ljdr": hist.Hist("Counts", dataset_axis, flav_axis, ljdr_axis),
            "sljdr": hist.Hist("Counts", dataset_axis, flav_axis, sljdr_axis),
        }

        self.sf_dict = list(_hist_sf_dict.keys())
        self.event_hists = list(_hist_event_dict.keys())
        self.btagDeephists = list(_hist_btagDeepdict.keys())
        _hist_dict = {**_hist_sf_dict, **_hist_event_dict, **_hist_btagDeepdict}
        # }
        self._accumulator = processor.dict_accumulator(_hist_dict)
        self._accumulator["sumw"] = processor.defaultdict_accumulator(float)
        ## Load corrections
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
        if not isRealData:
            weights.add("genweight", events.genWeight)
            weights.add(
                "puweight", self._pu[f"{self._year}_pileupweight"](events.Pileup.nPU)
            )
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
        ##############
        # Trigger level
        triggers = [
            # "HLT_IsoMu24",
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
        # req_opposite_charge = ak.fill_none(req_opposite_charge,False)
        event_jet = events.Jet[
            (events.Jet.pt > 25)
            & (abs(events.Jet.eta) <= 2.4)
            & (events.Jet.puId > 0)
            & (events.Jet.jetId > 5)
            & (events.Jet.btagDeepB > 0.0)
            & (events.Jet.btagDeepB < 1.0)
            & (events.Jet.btagDeepC > 0.0)
            & (events.Jet.btagDeepC < 1.0)
            & (events.Jet.btagDeepFlavB > 0.0)
            & (events.Jet.btagDeepFlavB < 1.0)
            & (events.Jet.btagDeepFlavC > 0.0)
            & (events.Jet.btagDeepFlavC < 1.0)
            & (ak.all(events.Jet.metric_table(events.Muon) > 0.4, axis=2))
            & (ak.all(events.Jet.metric_table(events.Electron) > 0.4, axis=2))
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

        # Per Electron
        el_eta = abs(selev.Electron.eta) < 2.4
        el_pt = selev.Electron.pt > 30
        el_idiso = selev.Electron.cutBased > 3
        el_level = el_eta & el_pt & el_idiso
        sel = selev.Electron[el_level]
        # Per jet : https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID
        jet_eta = abs(selev.Jet.eta) <= 2.4
        jet_pt = selev.Jet.pt > 25
        jet_pu = (selev.Jet.puId > 0) & (selev.Jet.jetId > 5)
        jet_dr = ak.all(selev.Jet.metric_table(smu) > 0.4, axis=2) & ak.all(
            selev.Jet.metric_table(sel) > 0.4, axis=2
        )
        jet_clean = (
            (selev.Jet.btagDeepB > 0.0)
            & (selev.Jet.btagDeepB < 1.0)
            & (selev.Jet.btagDeepC > 0.0)
            & (selev.Jet.btagDeepC < 1.0)
            & (selev.Jet.btagDeepFlavB > 0.0)
            & (selev.Jet.btagDeepFlavB < 1.0)
            & (selev.Jet.btagDeepFlavC > 0.0)
            & (selev.Jet.btagDeepFlavC < 1.0)
        )

        jet_level = jet_pu & jet_eta & jet_pt & jet_clean & jet_dr
        sjets = selev.Jet[jet_level]
        sel_jets = sjets
        sjets = sjets[:, :2]
        # b-tag twiki : https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
        bjet_disc = selev.Jet.btagDeepB > 0.4941
        bjet_level = jet_level & bjet_disc

        sbjets = selev.Jet[bjet_level]
        if not isRealData:
            stbjets = sbjets[sbjets.hadronFlavour == 5]
        if not isRealData:
            weights.add(
                "lep1sf",
                np.where(
                    event_level,
                    muSFs(
                        ak.firsts(
                            events.Muon[
                                (events.Muon.pt > 30)
                                & (abs(events.Muon.eta < 2.4))
                                & (events.Muon.tightId)
                            ]
                        ),
                        self._campaign,
                        correction_config[self._campaign]["LSF"],
                    ),
                    1.0,
                ),
            )
            weights.add(
                "lep2sf",
                np.where(
                    event_level,
                    muSFs(
                        ak.firsts(
                            events.Muon[
                                (events.Muon.pt > 30)
                                & (abs(events.Muon.eta < 2.4))
                                & (events.Muon.tightId)
                            ]
                        ),
                        self._campaign,
                        correction_config[self._campaign]["LSF"],
                    ),
                    1.0,
                ),
            )
        if isRealData:
            genflavor = ak.zeros_like(sjets.pt)
        else:
            par_flav = (sjets.partonFlavour == 0) & (sjets.hadronFlavour == 0)
            genflavor = sjets.hadronFlavour + 1 * par_flav
            jetsfs_c = collections.defaultdict(dict)
            jetsfs_b = collections.defaultdict(dict)
            csvsfs_c = collections.defaultdict(dict)
            csvsfs_b = collections.defaultdict(dict)
            for i in range(2):
                jetsfs_c[i]["SF"] = getSF(
                    sjets[:, i].hadronFlavour,
                    sjets[:, i].btagDeepFlavCvL,
                    sjets[:, i].btagDeepFlavCvB,
                    self._deepjetc_sf,
                )
                jetsfs_c[i]["SFup"] = getSF(
                    sjets[:, i].hadronFlavour,
                    sjets[:, i].btagDeepFlavCvL,
                    sjets[:, i].btagDeepFlavCvB,
                    self._deepjetc_sf,
                    "TotalUncUp",
                )
                jetsfs_c[i]["SFdn"] = getSF(
                    sjets[:, i].hadronFlavour,
                    sjets[:, i].btagDeepFlavCvL,
                    sjets[:, i].btagDeepFlavCvB,
                    self._deepjetc_sf,
                    "TotalUncDown",
                )
                jetsfs_b[i]["SF"] = self._deepjetb_sf.eval(
                    "central",
                    sjets[:, i].hadronFlavour,
                    abs(sjets[:, i].eta),
                    sjets[:, i].pt,
                    discr=sjets[:, i].btagDeepFlavB,
                )
                jetsfs_b[i]["SFup"] = self._deepjetb_sf.eval(
                    "up_jes",
                    sjets[:, i].hadronFlavour,
                    abs(sjets[:, i].eta),
                    sjets[:, i].pt,
                    discr=sjets[:, i].btagDeepFlavB,
                )
                jetsfs_b[i]["SFdn"] = self._deepjetb_sf.eval(
                    "down_jes",
                    sjets[:, i].hadronFlavour,
                    abs(sjets[:, i].eta),
                    sjets[:, i].pt,
                    discr=sjets[:, i].btagDeepFlavB,
                )
                csvsfs_c[i]["SF"] = getSF(
                    sjets[:, i].hadronFlavour,
                    sjets[:, i].btagDeepCvL,
                    sjets[:, i].btagDeepCvB,
                    self._deepcsvc_sf,
                )
                csvsfs_c[i]["SFup"] = getSF(
                    sjets[:, i].hadronFlavour,
                    sjets[:, i].btagDeepCvL,
                    sjets[:, i].btagDeepCvB,
                    self._deepcsvc_sf,
                    "TotalUncUp",
                )
                csvsfs_c[i]["SFdn"] = getSF(
                    sjets[:, i].hadronFlavour,
                    sjets[:, i].btagDeepCvL,
                    sjets[:, i].btagDeepCvB,
                    self._deepcsvc_sf,
                    "TotalUncDown",
                )
                csvsfs_b[i]["SFup"] = self._deepcsvb_sf.eval(
                    "up_jes",
                    sjets[:, i].hadronFlavour,
                    abs(sjets[:, i].eta),
                    sjets[:, i].pt,
                    discr=sjets[:, i].btagDeepB,
                )
                csvsfs_b[i]["SF"] = self._deepcsvb_sf.eval(
                    "central",
                    sjets[:, i].hadronFlavour,
                    abs(sjets[:, i].eta),
                    sjets[:, i].pt,
                    discr=sjets[:, i].btagDeepB,
                )
                csvsfs_b[i]["SFdn"] = self._deepcsvb_sf.eval(
                    "down_jes",
                    sjets[:, i].hadronFlavour,
                    abs(sjets[:, i].eta),
                    sjets[:, i].pt,
                    discr=sjets[:, i].btagDeepB,
                )
        disc_list = {
            "btagDeepB": csvsfs_b,
            "btagDeepC": csvsfs_b,
            "btagDeepFlavB": jetsfs_b,
            "btagDeepFlavC": jetsfs_b,
            "btagDeepCvL": csvsfs_c,
            "btagDeepCvB": csvsfs_c,
            "btagDeepFlavCvL": jetsfs_c,
            "btagDeepFlavCvB": jetsfs_c,
        }
        for histname, h in output.items():
            if histname in self.btagDeephists:
                if isRealData:
                    fields = {
                        l: ak.flatten(sjets[l], axis=None)
                        for l in h.fields
                        if l in dir(sjets)
                    }
                    h.fill(dataset=dataset, flav=5, **fields)
                else:
                    genweiev = ak.flatten(
                        ak.broadcast_arrays(weights.weight()[event_level], sjets["pt"])[
                            0
                        ]
                    )
                    fields = {
                        l: ak.flatten(sjets[histname])
                        for l in h.fields
                        if l in dir(sjets)
                    }
                    h.fill(
                        dataset=dataset,
                        flav=ak.flatten(genflavor),
                        **fields,
                        weight=genweiev,
                    )
            elif "btagDeep" in histname:
                i = int(histname[-1])
                sel_jet = sjets[:, i]
                fields = {
                    l: np.where(sel_jet[l] < 0, -0.2, sel_jet[l])
                    for l in h.fields
                    if l in dir(sel_jet)
                }
                if isRealData:
                    h.fill(dataset=dataset, flav=genflavor[:, i], syst="noSF", **fields)
                else:
                    h.fill(
                        dataset=dataset,
                        flav=genflavor[:, i],
                        syst="noSF",
                        **fields,
                        weight=weights.weight()[event_level],
                    )
                    for syst in disc_list[histname.replace(f"_{i}", "")][i].keys():
                        h.fill(
                            dataset=dataset,
                            flav=genflavor[:, i],
                            syst=syst,
                            **fields,
                            weight=weights.weight()[event_level]
                            * disc_list[histname.replace(f"_{i}", "")][i][syst],
                        )

        def flatten(ar):  # flatten awkward into a 1d array to hist
            return ak.flatten(ar, axis=None)

        def num(ar):
            return ak.num(ak.fill_none(ar[~ak.is_none(ar)], 0), axis=0)

        if isRealData:
            output["njet"].fill(dataset=dataset, njet=flatten(ak.num(sel_jets)))
            output["nbjet"].fill(dataset=dataset, nbjet=flatten(ak.num(sbjets)))
            output["ntbjet"].fill(dataset=dataset, ntbjet=flatten(ak.num(sbjets)))
            output["ljpt"].fill(dataset=dataset, flav=0, ljpt=flatten(sjets[:, 0].pt))
            output["sljpt"].fill(dataset=dataset, flav=0, sljpt=flatten(sjets[:, 1].pt))
            output["ljeta"].fill(
                dataset=dataset, flav=0, ljeta=flatten(sjets[:, 0].eta)
            )
            output["sljeta"].fill(
                dataset=dataset, flav=0, sljeta=flatten(sjets[:, 1].eta)
            )
        else:
            output["njet"].fill(
                dataset=dataset,
                njet=flatten(ak.num(sel_jets)),
                weight=weights.weight()[event_level],
            )
            output["nbjet"].fill(
                dataset=dataset,
                nbjet=flatten(ak.num(sbjets)),
                weight=weights.weight()[event_level]
                * csvsfs_b[0]["SF"]
                * csvsfs_b[1]["SF"],
            )
            output["nbjet_up"].fill(
                dataset=dataset,
                nbjet=flatten(ak.num(sbjets)),
                weight=weights.weight()[event_level]
                * csvsfs_b[0]["SFup"]
                * csvsfs_b[1]["SFup"],
            )
            output["nbjet_dn"].fill(
                dataset=dataset,
                nbjet=flatten(ak.num(sbjets)),
                weight=weights.weight()[event_level]
                * csvsfs_b[0]["SFdn"]
                * csvsfs_b[1]["SFdn"],
            )
            output["ntbjet"].fill(
                dataset=dataset,
                ntbjet=flatten(ak.num(stbjets)),
                weight=weights.weight()[event_level],
            )
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
            output["ljeta"].fill(
                dataset=dataset,
                flav=flatten(sjets[:, 0].hadronFlavour),
                ljeta=flatten(sjets[:, 0].eta),
                weight=weights.weight()[event_level],
            )
            output["sljeta"].fill(
                dataset=dataset,
                flav=flatten(sjets[:, 1].hadronFlavour),
                sljeta=flatten(sjets[:, 1].eta),
                weight=weights.weight()[event_level],
            )

        return output

    def postprocess(self, accumulator):
        return accumulator
