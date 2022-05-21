import gzip
import pickle, os, sys, mplhep as hep, numpy as np
import collections

from matplotlib.pyplot import jet

from coffea import hist, processor
import awkward as ak
from coffea.analysis_tools import Weights

from BTVNanoCommissioning.utils.correction import (
    lumiMasks,
    eleSFs,
    muSFs,
    load_pu,
    load_BTV,
)
from BTVNanoCommissioning.helpers.definitions import definitions
from BTVNanoCommissioning.helpers.cTagSFReader import getSF
from BTVNanoCommissioning.utils.AK4_parameters import correction_config


class NanoProcessor(processor.ProcessorABC):
    # Define histograms

    def num(ar):
        return ak.num(ak.fill_none(ar[~ak.is_none(ar)], 0), axis=0)

    def __init__(self, year="2017", campaign="Rereco17_94X"):
        self._year = year
        self._campaign = campaign
        # Define axes
        # Should read axes from NanoAOD config
        dataset_axis = hist.Cat("dataset", "Primary dataset")
        flav_axis = hist.Bin("flav", r"Genflavour", [0, 1, 4, 5, 6])

        # Jet
        jet_pt_axis = hist.Bin("pt", r"Jet $p_{T}$ [GeV]", 50, 0, 500)
        jet_eta_axis = hist.Bin("eta", r"Jet $\eta$", 25, -2.5, 2.5)
        jet_phi_axis = hist.Bin("phi", r"Jet $\phi$", 30, -3, 3)
        jet_mass_axis = hist.Bin("mass", r"Jet $m$ [GeV]", 50, 0, 500)

        # Events
        njet_axis = hist.Bin("njet", r"N jets", [0, 1, 2, 3])
        ## Muons
        hard_lpt_axis = hist.Bin("pt", r"Hard lepton pt", 40, 0, 200)
        hard_leta_axis = hist.Bin("eta", r"Lepton $\eta$", 25, -2.5, 2.5)
        hard_lphi_axis = hist.Bin("phi", r"Lepton $\phi$", 30, -3, 3)
        hard_liso_axis = hist.Bin("pfRelIso03_all", r"Hard Muon Rel. Iso", 40, 0, 4.0)
        soft_lpt_axis = hist.Bin("pt", r"Soft lepton pt", 40, 0, 40)
        soft_leta_axis = hist.Bin("eta", r"Lepton $\eta$", 25, -2.5, 2.5)
        soft_lphi_axis = hist.Bin("phi", r"Lepton $\phi$", 30, -3, 3)
        soft_liso_axis = hist.Bin("pfRelIso03_all", r"Soft Muon Rel. Iso", 40, 0, 4.0)
        l_dxy_axis = hist.Bin("dxy", r"dxy", 20, 0, 0.002)
        l_dz_axis = hist.Bin("dz", r"dz", 20, 0, 0.01)
        l_sip3d_axis = hist.Bin("dz", r"dz", 20, 0, 0.2)

        ## Z
        zmass_axis = hist.Bin("zmass", r"Z Mass", 25, 50, 100)
        zpt_axis = hist.Bin("zpt", r"Z $p_{T}$", 25, 0, 100)
        zeta_axis = hist.Bin("zeta", r"Z $\eta$", 25, -2.5, 2.5)
        zphi_axis = hist.Bin("zphi", r"Z $\phi$", 30, -3, 3)
        drmumu_axis = hist.Bin(
            "drmumu", r"$\Delta$R($\mu_{soft}$,$\mu_{hard}$)", 25, 0, 5
        )

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
        _hist_btagDeepdict = {
            "pt": hist.Hist("Counts", dataset_axis, flav_axis, jet_pt_axis),
            "eta": hist.Hist("Counts", dataset_axis, flav_axis, jet_eta_axis),
            "phi": hist.Hist("Counts", dataset_axis, flav_axis, jet_phi_axis),
            "mass": hist.Hist("Counts", dataset_axis, flav_axis, jet_mass_axis),
        }
        for disc, axis in zip(disc_list, btag_axes):
            for i in range(1):
                _hist_sf_dict["%s_%d" % (disc, i)] = hist.Hist(
                    "Counts", dataset_axis, flav_axis, syst_axis, axis
                )
        for deepcsv, axises in zip(bininfo.keys(), btagDeepaxes):
            _hist_btagDeepdict["%s" % (deepcsv)] = hist.Hist(
                "Counts", dataset_axis, flav_axis, axises
            )

        _hist_event_dict = {
            "njet": hist.Hist("Counts", dataset_axis, njet_axis),
            "pos_pt": hist.Hist("Counts", dataset_axis, hard_lpt_axis),
            "pos_eta": hist.Hist("Counts", dataset_axis, hard_leta_axis),
            "pos_phi": hist.Hist("Counts", dataset_axis, hard_lphi_axis),
            "pos_pfRelIso04_all": hist.Hist("Counts", dataset_axis, hard_liso_axis),
            "neg_pt": hist.Hist("Counts", dataset_axis, soft_lpt_axis),
            "neg_eta": hist.Hist("Counts", dataset_axis, soft_leta_axis),
            "neg_phi": hist.Hist("Counts", dataset_axis, soft_lphi_axis),
            "neg_pfRelIso04_all": hist.Hist("Counts", dataset_axis, soft_liso_axis),
            "zmass": hist.Hist("Counts", dataset_axis, zmass_axis),
            "zpt": hist.Hist("Counts", dataset_axis, zpt_axis),
            "zeta": hist.Hist("Counts", dataset_axis, zeta_axis),
            "zphi": hist.Hist("Counts", dataset_axis, zphi_axis),
            "drmumu": hist.Hist("Counts", dataset_axis, drmumu_axis),
        }

        self.sf_hists = list(_hist_sf_dict.keys())
        self.event_hists = list(_hist_event_dict.keys())
        self.btagDeephists = list(_hist_btagDeepdict.keys())
        _hist_dict = {**_hist_sf_dict, **_hist_event_dict, **_hist_btagDeepdict}
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
        mu_triggers = [
            # "HLT_IsoMu24",
            "Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
        ]

        trig_arrs = [events.HLT[_trig] for _trig in mu_triggers]
        req_trig = np.zeros(len(events), dtype="bool")
        for t in trig_arrs:
            req_trig = req_trig | t

        ############
        # Event level
        ## Muon cuts
        dilep_mu = events.Muon[
            (events.Muon.pt > 12)
            & (abs(events.Muon.eta) < 2.4)
            & (events.Muon.tightId > 0.5)
            & (events.Muon.pfRelIso04_all <= 0.15)
        ]
        dilep_ele = events.Electron[
            (events.Electron.pt > 15)
            & (
                (abs(events.Electron.eta) < 1.4442)
                | (
                    (abs(events.Electron.eta) < 2.5)
                    & (abs(events.Electron.eta) > 1.566)
                )
            )
            & (events.Electron.mvaFall17V2Iso_WP80 > 0.5)
        ]
        pos_dilep = dilep_ele[dilep_ele.charge > 0]
        neg_dilep = dilep_ele[dilep_ele.charge < 0]
        req_dilep = (
            (ak.num(pos_dilep.pt) >= 1)
            & (ak.num(neg_dilep.pt) >= 1)
            & (ak.num(dilep_ele.charge) >= 2)
            & (ak.num(dilep_mu.charge) == 0)
        )
        pos_dilep = ak.pad_none(pos_dilep, 1, axis=1)
        neg_dilep = ak.pad_none(neg_dilep, 1, axis=1)
        # dilepton mass

        dilep_mass = pos_dilep[:, 0] + neg_dilep[:, 0]
        req_dilepmass = (
            (dilep_mass.mass > 81)
            & (dilep_mass.mass < 101)
            & (dilep_mass.pt > 15)
            & ((pos_dilep[:, 0].pt > 27) | (neg_dilep[:, 0].pt > 27))
        )

        ## Jet cuts
        event_jet = events.Jet[
            (events.Jet.pt > 20)
            & (abs(events.Jet.eta) <= 2.5)
            & ((events.Jet.puId >= 7) & (events.Jet.pt < 50))
            & (events.Jet.jetId >= 3)
            & (events.Jet.btagDeepB > 0.0)
            & (events.Jet.btagDeepB < 1.0)
            & (events.Jet.btagDeepC > 0.0)
            & (events.Jet.btagDeepC < 1.0)
            & (events.Jet.btagDeepFlavB > 0.0)
            & (events.Jet.btagDeepFlavB < 1.0)
            & (events.Jet.btagDeepFlavC > 0.0)
            & (events.Jet.btagDeepFlavC < 1.0)
            & (ak.all(events.Jet.metric_table(pos_dilep[:, 0]) > 0.4, axis=2))
            & (ak.all(events.Jet.metric_table(neg_dilep[:, 0]) > 0.4, axis=2))
        ]

        req_jets = ak.num(event_jet.puId) >= 1
        event_jet = ak.pad_none(event_jet, 1, axis=1)
        event_level = req_lumi & req_trig & req_dilep & req_dilepmass & req_jets
        if len(event_level) > 0:
            event_level = ak.fill_none(event_level, False)
        # Selected
        selev = events[event_level]

        #########

        ## Hard Muon

        sposmu = selev.Electron[
            (selev.Electron.pt > 15)
            & (
                (abs(selev.Electron.eta) < 1.4442)
                | ((abs(selev.Electron.eta) < 2.5) & (abs(selev.Electron.eta) > 1.566))
            )
            & (selev.Electron.mvaFall17V2Iso_WP80 > 0.5)
            & (selev.Electron.charge > 0)
        ]
        sposmu = sposmu[:, 0]
        snegmu = selev.Electron[
            (selev.Electron.pt > 15)
            & (
                (abs(selev.Electron.eta) < 1.4442)
                | ((abs(selev.Electron.eta) < 2.5) & (abs(selev.Electron.eta) > 1.566))
            )
            & (selev.Electron.mvaFall17V2Iso_WP80 > 0.5)
            & (selev.Electron.charge < 0)
        ]
        snegmu = snegmu[:, 0]

        if not isRealData:
            weights.add(
                "lep1sf",
                np.where(
                    event_level,
                    muSFs(
                        ak.firsts(
                            events.Muon[
                                (events.Muon.pt > 12)
                                & (abs(events.Muon.eta) < 2.4)
                                & (events.Muon.tightId > 0.5)
                                & (events.Muon.pfRelIso04_all <= 0.15)
                                & (events.Muon.charge < 0)
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
                                (events.Muon.pt > 12)
                                & (abs(events.Muon.eta) < 2.4)
                                & (events.Muon.tightId > 0.5)
                                & (events.Muon.pfRelIso04_all <= 0.15)
                                & (events.Muon.charge > 0)
                            ]
                        ),
                        self._campaign,
                        correction_config[self._campaign]["LSF"],
                    ),
                    1.0,
                ),
            )
        sz = sposmu + snegmu

        ## Jets

        sjets = selev.Jet[
            (selev.Jet.pt > 20)
            & (abs(selev.Jet.eta) <= 2.5)
            & ((selev.Jet.puId >= 7) & (selev.Jet.pt < 50))
            & (selev.Jet.jetId >= 3)
            & (selev.Jet.btagDeepFlavC < 1.0)
            & (ak.all(selev.Jet.metric_table(sposmu) > 0.4, axis=2))
            & (ak.all(selev.Jet.metric_table(snegmu) > 0.4, axis=2))
            & (selev.Jet.muEF < 0.8)
        ]

        njet = ak.count(sjets.pt, axis=1)

        def flatten(ar):  # flatten awkward into a 1d array to hist
            return ak.flatten(ar, axis=None)

        if isRealData:
            genflavor = ak.zeros_like(sjets.pt)
        else:
            par_flav = (sjets.partonFlavour == 0) & (sjets.hadronFlavour == 0)
            genflavor = sjets.hadronFlavour + 1 * par_flav
            jetsfs_c = collections.defaultdict(dict)
            jetsfs_b = collections.defaultdict(dict)
            csvsfs_c = collections.defaultdict(dict)
            csvsfs_b = collections.defaultdict(dict)

            ## for each jet
            jetsfs_c[0]["SF"] = getSF(
                sjets[:, 0].hadronFlavour,
                sjets[:, 0].btagDeepFlavCvL,
                sjets[:, 0].btagDeepFlavCvB,
                self._deepjetc_sf,
            )
            jetsfs_c[0]["SFup"] = getSF(
                sjets[:, 0].hadronFlavour,
                sjets[:, 0].btagDeepFlavCvL,
                sjets[:, 0].btagDeepFlavCvB,
                self._deepjetc_sf,
                "TotalUncUp",
            )
            jetsfs_c[0]["SFdn"] = getSF(
                sjets[:, 0].hadronFlavour,
                sjets[:, 0].btagDeepFlavCvL,
                sjets[:, 0].btagDeepFlavCvB,
                self._deepjetc_sf,
                "TotalUncDown",
            )
            jetsfs_b[0]["SF"] = self._deepjetb_sf.eval(
                "central",
                sjets[:, 0].hadronFlavour,
                abs(sjets[:, 0].eta),
                sjets[:, 0].pt,
                discr=sjets[:, 0].btagDeepFlavB,
            )
            jetsfs_b[0]["SFup"] = self._deepjetb_sf.eval(
                "up_jes",
                sjets[:, 0].hadronFlavour,
                abs(sjets[:, 0].eta),
                sjets[:, 0].pt,
                discr=sjets[:, 0].btagDeepFlavB,
            )
            jetsfs_b[0]["SFdn"] = self._deepjetb_sf.eval(
                "down_jes",
                sjets[:, 0].hadronFlavour,
                abs(sjets[:, 0].eta),
                sjets[:, 0].pt,
                discr=sjets[:, 0].btagDeepFlavB,
            )
            csvsfs_c[0]["SF"] = getSF(
                sjets[:, 0].hadronFlavour,
                sjets[:, 0].btagDeepCvL,
                sjets[:, 0].btagDeepCvB,
                self._deepcsvc_sf,
            )
            csvsfs_c[0]["SFup"] = getSF(
                sjets[:, 0].hadronFlavour,
                sjets[:, 0].btagDeepCvL,
                sjets[:, 0].btagDeepCvB,
                self._deepcsvc_sf,
                "TotalUncUp",
            )
            csvsfs_c[0]["SFdn"] = getSF(
                sjets[:, 0].hadronFlavour,
                sjets[:, 0].btagDeepCvL,
                sjets[:, 0].btagDeepCvB,
                self._deepcsvc_sf,
                "TotalUncDown",
            )
            csvsfs_b[0]["SFup"] = self._deepcsvb_sf.eval(
                "up_jes",
                sjets[:, 0].hadronFlavour,
                abs(sjets[:, 0].eta),
                sjets[:, 0].pt,
                discr=sjets[:, 0].btagDeepB,
            )
            csvsfs_b[0]["SF"] = self._deepcsvb_sf.eval(
                "central",
                sjets[:, 0].hadronFlavour,
                abs(sjets[:, 0].eta),
                sjets[:, 0].pt,
                discr=sjets[:, 0].btagDeepB,
            )
            csvsfs_b[0]["SFdn"] = self._deepcsvb_sf.eval(
                "down_jes",
                sjets[:, 0].hadronFlavour,
                abs(sjets[:, 0].eta),
                sjets[:, 0].pt,
                discr=sjets[:, 0].btagDeepB,
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
                fields = {
                    l: ak.flatten(sjets[histname]) for l in h.fields if l in dir(sjets)
                }
                if isRealData:
                    h.fill(dataset=dataset, flav=5, **fields)
                else:
                    h.fill(
                        dataset=dataset,
                        flav=ak.flatten(genflavor),
                        **fields,
                        weight=weights.weight()[event_level],
                    )
            elif "pos_" in histname:
                fields = {l: sposmu[l] for l in h.fields if l in dir(sposmu)}
                if isRealData:
                    h.fill(dataset=dataset, **fields)
                else:
                    h.fill(
                        dataset=dataset, **fields, weight=weights.weight()[event_level]
                    )
            elif "neg_" in snegmu:
                fields = {l: snegmu[l] for l in h.fields if l in dir(snegmu)}
                if isRealData:
                    h.fill(dataset=dataset, **fields)
                else:
                    h.fill(
                        dataset=dataset, **fields, weight=weights.weight()[event_level]
                    )

            elif ["zmass", "zpt", "zeta", "zphi"] == histname:
                fields = {l: sz[l] for l in h.fields if l in dir(sz)}
                if isRealData:
                    h.fill(dataset=dataset, **fields)
                else:
                    h.fill(
                        dataset=dataset, **fields, weight=weights.weight()[event_level]
                    )
            elif "btagDeep" in histname and "0" in histname:
                sel_jet = sjets[:, 0]
                fields = {
                    l: np.where(sel_jet[l] < 0, -0.2, sel_jet[l])
                    for l in h.fields
                    if l in dir(sel_jet)
                }
                if isRealData:
                    h.fill(dataset=dataset, flav=5, syst="noSF", **fields)
                else:
                    h.fill(
                        dataset=dataset,
                        flav=genflavor[:, 0],
                        syst="noSF",
                        **fields,
                        weight=weights.weight()[event_level],
                    )
                    for syst in disc_list[histname.replace("_0", "")][0].keys():
                        h.fill(
                            dataset=dataset,
                            flav=genflavor[:, 0],
                            syst=syst,
                            **fields,
                            weight=weights.weight()[event_level]
                            * disc_list[histname.replace("_0", "")][0][syst],
                        )
        if not isRealData:
            output["drmumu"].fill(
                dataset=dataset,
                drmumu=snegmu.delta_r(sposmu),
                weight=weights.weight()[event_level],
            )
        else:
            output["drmumu"].fill(dataset=dataset, drmumu=snegmu.delta_r(sposmu))

        return output

    def postprocess(self, accumulator):
        return accumulator
