import gzip
import pickle, os, sys, mplhep as hep, numpy as np
import collections

from matplotlib.pyplot import jet

import coffea
from coffea import hist, processor
import awkward as ak
from coffea.analysis_tools import Weights
from BTVNanoCommissioning.utils.correction import lumiMasks, eleSFs, load_pu, load_BTV
from BTVNanoCommissioning.helpers.definitions import definitions
import gc
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
        charge_axis = hist.Bin("char", r"Charge", [-2, 0, 2])

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
        hard_lpt_ratio_axis = hist.Bin(
            "hlptratio", r"Hard $\mu$/Jet $p_{T}$ [GeV]", 40, 0, 1
        )
        soft_lpt_axis = hist.Bin("pt", r"Soft lepton pt", 40, 0, 40)
        soft_leta_axis = hist.Bin("eta", r"Lepton $\eta$", 25, -2.5, 2.5)
        soft_lphi_axis = hist.Bin("phi", r"Lepton $\phi$", 30, -3, 3)
        soft_liso_axis = hist.Bin("pfRelIso04_all", r"Soft Muon Rel. Iso", 40, 0, 4.0)
        soft_lpt_ratio_axis = hist.Bin(
            "slptratio", r"Soft $\mu$/Jet $p_{T}$ [GeV]", 40, 0, 1
        )
        l_dxy_axis = hist.Bin("dxy", r"dxy", 20, 0, 0.002)
        l_dz_axis = hist.Bin("dz", r"dz", 20, 0, 0.01)
        l_sip3d_axis = hist.Bin("dz", r"dz", 20, 0, 0.2)

        ## Z/W
        zmass_axis = hist.Bin("zmass", r"Z Mass", 25, 50, 100)
        zpt_axis = hist.Bin("zpt", r"Z $p_{T}$", 25, 0, 100)
        zeta_axis = hist.Bin("zeta", r"Z $\eta$", 25, -2.5, 2.5)
        zphi_axis = hist.Bin("zphi", r"Z $\phi$", 30, -3, 3)
        drmumu_axis = hist.Bin(
            "drmumu", r"$\Delta$R($\mu_{soft}$,$\mu_{hard}$)", 25, 0, 5
        )
        wmass_axis = hist.Bin("wmass", r"W mass", 25, 50, 100)
        wpt_axis = hist.Bin("wpt", r"W $p_{T}$", 25, 0, 100)
        weta_axis = hist.Bin("weta", r"W $\eta$", 25, -2.5, 2.5)
        wphi_axis = hist.Bin("wphi", r"W $\phi$", 30, -3, 3)

        ## MET
        met_axis = hist.Bin("pt", r"MET $p_{T}$", 50, 0, 500)
        metphi_axis = hist.Bin("phi", r"met $\phi$", 30, -3, 3)

        ## Muon jets
        mujet_pt_axis = hist.Bin("pt", r"Jet $p_{T}$ [GeV]", 50, 0, 500)
        mujet_eta_axis = hist.Bin("eta", r"Jet $\eta$", 25, -2.5, 2.5)
        mujet_phi_axis = hist.Bin("phi", r"Jet $\phi$", 30, -3, 3)
        mujet_mass_axis = hist.Bin("mass", r"Jet $m$ [GeV]", 50, 0, 500)
        dr_mujetsoftmu_axis = hist.Bin(
            "drjet_smu", r"$\Delta$R($\mu_{soft}$,j)", 25, 0, 5
        )
        dr_mujethardmu_axis = hist.Bin(
            "drjet_hmu", r"$\Delta$R($\mu_{hard}$,j)", 25, 0, 5
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
            "hl_pt": hist.Hist("Counts", dataset_axis, hard_lpt_axis),
            "hl_eta": hist.Hist("Counts", dataset_axis, hard_leta_axis),
            "hl_phi": hist.Hist("Counts", dataset_axis, hard_lphi_axis),
            "hl_pfRelIso04_all": hist.Hist("Counts", dataset_axis, hard_liso_axis),
            "sl_pt": hist.Hist("Counts", dataset_axis, soft_lpt_axis),
            "sl_eta": hist.Hist("Counts", dataset_axis, soft_leta_axis),
            "sl_phi": hist.Hist("Counts", dataset_axis, soft_lphi_axis),
            "sl_pfRelIso04_all": hist.Hist("Counts", dataset_axis, soft_liso_axis),
            "sl_dxy": hist.Hist("Counts", dataset_axis, l_dxy_axis),
            "sl_dz": hist.Hist("Counts", dataset_axis, l_dz_axis),
            "sl_sip3d": hist.Hist("Counts", dataset_axis, l_sip3d_axis),
            "hlptratio": hist.Hist(
                "Counts", dataset_axis, flav_axis, hard_lpt_ratio_axis
            ),
            "slptratio": hist.Hist(
                "Counts", dataset_axis, flav_axis, soft_lpt_ratio_axis
            ),
            "zmass": hist.Hist("Counts", dataset_axis, zmass_axis),
            "zpt": hist.Hist("Counts", dataset_axis, zpt_axis),
            "zeta": hist.Hist("Counts", dataset_axis, zeta_axis),
            "zphi": hist.Hist("Counts", dataset_axis, zphi_axis),
            "wmass": hist.Hist("Counts", dataset_axis, wmass_axis),
            "wpt": hist.Hist("Counts", dataset_axis, wpt_axis),
            "weta": hist.Hist("Counts", dataset_axis, weta_axis),
            "wphi": hist.Hist("Counts", dataset_axis, wphi_axis),
            "drmumu": hist.Hist("Counts", dataset_axis, drmumu_axis),
            "metpt": hist.Hist("Counts", dataset_axis, met_axis),
            "metphi": hist.Hist("Counts", dataset_axis, metphi_axis),
            "drjet_smu": hist.Hist(
                "Counts", dataset_axis, flav_axis, dr_mujetsoftmu_axis
            ),
            "drjet_hmu": hist.Hist(
                "Counts", dataset_axis, flav_axis, dr_mujethardmu_axis
            ),
        }

        self.sf_hists = list(_hist_sf_dict.keys())
        self.event_hists = list(_hist_event_dict.keys())
        self.btagDeephists = list(_hist_btagDeepdict.keys())
        _hist_dict = {**_hist_sf_dict, **_hist_event_dict, **_hist_btagDeepdict}
        # ,**_hist_btagDeepdict}
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
            "HLT_Ele32_WPTight_Gsf_L1DoubleEG",
        ]

        trig_arrs = [events.HLT[_trig.strip("HLT_")] for _trig in triggers]
        req_trig = np.zeros(len(events), dtype="bool")
        for t in trig_arrs:
            req_trig = req_trig | t

        ############
        # Event level
        ## Electron cuts
        iso_ele = events.Electron[
            (events.Electron.pt > 34)
            & (abs(events.Electron.eta) < 2.5)
            & (events.Electron.mvaFall17V2Iso_WP80 > 0.5)
        ]

        req_ele = ak.count(iso_ele.pt, axis=1) == 1
        iso_ele = ak.pad_none(iso_ele, 1, axis=1)
        iso_ele = iso_ele[:, 0]

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
        req_dilepveto = (
            ak.count(dilep_mu.pt, axis=1) + ak.count(dilep_ele.pt, axis=1) != 2
        )

        MET = ak.zip(
            {
                "pt": events.METFixEE2017.pt,
                "eta": ak.zeros_like(events.METFixEE2017.pt),
                "phi": events.METFixEE2017.phi,
                "mass": ak.zeros_like(events.METFixEE2017.pt),
            },
            with_name="PtEtaPhiMLorentzVector",
        )

        Wmass = MET + iso_ele
        req_Wmass = Wmass.mass > 55

        ## Jet cuts
        event_jet = events.Jet[
            (events.Jet.pt > 20)
            & (abs(events.Jet.eta) <= 2.5)
            & (((events.Jet.puId >= 7) & (events.Jet.pt < 50)) | (events.Jet.pt >= 50))
            & (events.Jet.jetId >= 3)
            & (events.Jet.btagDeepB > 0.0)
            & (events.Jet.btagDeepB < 1.0)
            & (events.Jet.btagDeepC > 0.0)
            & (events.Jet.btagDeepC < 1.0)
            & (events.Jet.btagDeepFlavB > 0.0)
            & (events.Jet.btagDeepFlavB < 1.0)
            & (events.Jet.btagDeepFlavC > 0.0)
            & (events.Jet.btagDeepFlavC < 1.0)
            & (ak.all(events.Jet.metric_table(iso_ele) > 0.5, axis=2))
        ]

        req_jets = ak.num(event_jet.puId) >= 4

        ## Soft Muon cuts

        soft_muon = events.Muon[
            (events.Muon.pt < 25)
            & (abs(events.Muon.eta) < 2.4)
            & (events.Muon.tightId > 0.5)
            & (events.Muon.pfRelIso04_all > 0.2)
            & (events.Muon.jetIdx != -1)
        ]
        req_softmu = ak.count(soft_muon.pt, axis=1) >= 1

        mu_jet = event_jet[
            (ak.all(event_jet.metric_table(soft_muon) <= 0.4, axis=2))
            & ((event_jet.muonIdx1 != -1) | (event_jet.muonIdx2 != -1))
        ]
        req_mujet = ak.num(mu_jet.pt, axis=1) >= 1
        # mu_jet = ak.fill_none(mu_jet.pt,0)
        mu_jet = ak.pad_none(mu_jet, 1, axis=1)
        soft_muon = ak.pad_none(soft_muon, 1, axis=1)
        # pT ratio
        req_pTratio = (soft_muon[:, 0].pt / mu_jet[:, 0].pt) < 0.6

        req_QCDveto = (
            (iso_ele.pfRelIso03_all < 0.05)
            & (abs(iso_ele.dz) < 0.01)
            & (abs(iso_ele.dxy) < 0.002)
            & (iso_ele.ip3d < 0.2)
            & (
                (iso_ele.pt / mu_jet[:, 0].pt < 0.0)
                | (iso_ele.pt / mu_jet[:, 0].pt > 0.75)
            )
        )
        event_level = (
            req_trig
            & req_lumi
            & req_ele
            & req_jets
            & req_softmu
            & req_mujet
            & req_Wmass
            & req_dilepveto
            & req_QCDveto
            & req_pTratio
        )
        if len(event_level) > 0:
            event_level = ak.fill_none(event_level, False)
        # Selected
        selev = events[event_level]

        #########

        ## Hard Muon
        shmu = selev.Electron[
            (selev.Electron.pt > 34)
            & (abs(selev.Electron.eta) < 2.5)
            & (selev.Electron.mvaFall17V2Iso_WP80 > 0.5)
        ]

        shmu = shmu[:, 0]

        sjets = selev.Jet[
            (selev.Jet.pt > 20)
            & (abs(selev.Jet.eta) <= 2.5)
            & (((selev.Jet.puId >= 7) & (selev.Jet.pt < 50)) | (selev.Jet.pt >= 50))
            & (selev.Jet.jetId >= 3)
            & (selev.Jet.btagDeepB > 0.0)
            & (selev.Jet.btagDeepB < 1.0)
            & (selev.Jet.btagDeepC > 0.0)
            & (selev.Jet.btagDeepC < 1.0)
            & (selev.Jet.btagDeepFlavB > 0.0)
            & (selev.Jet.btagDeepFlavB < 1.0)
            & (selev.Jet.btagDeepFlavC > 0.0)
            & (selev.Jet.btagDeepFlavC < 1.0)
            & (ak.all(selev.Jet.metric_table(shmu) > 0.5, axis=2))
        ]
        ## Soft Muon
        ssmu = selev.Muon[
            (selev.Muon.pt < 25)
            & (abs(selev.Muon.eta) < 2.4)
            & (selev.Muon.tightId > 0.5)
            & (selev.Muon.pfRelIso04_all > 0.2)
            & (selev.Muon.jetIdx != -1)
        ]
        # ssmu=ssmu[:,0]

        ## MET
        smet = ak.zip(
            {
                "pt": selev.METFixEE2017.pt,
                "eta": ak.zeros_like(selev.METFixEE2017.pt),
                "phi": selev.METFixEE2017.phi,
                "mass": ak.zeros_like(selev.METFixEE2017.pt),
            },
            with_name="PtEtaPhiMLorentzVector",
        )

        ## Muon Jet
        # print(sjets.pt)
        smuon_jet = sjets[
            (ak.all(sjets.metric_table(ssmu) <= 0.4, axis=2))
            & ((sjets.muonIdx1 != -1) | (sjets.muonIdx2 != -1))
        ]
        # print(smuon_jet.pt)
        smuon_jet = smuon_jet[:, 0]
        ssmu = ssmu[:, 0]
        sz = shmu + ssmu
        sw = shmu + smet
        if not isRealData:

            weights.add(
                "lep1sf",
                np.where(
                    event_level,
                    eleSFs(
                        ak.firsts(
                            events.Electron[
                                (events.Electron.pt > 34)
                                & (abs(events.Electron.eta) < 2.5)
                                & (events.Electron.mvaFall17V2Iso_WP80 > 0.5)
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
                    eleSFs(
                        ak.firsts(
                            events.Electron[
                                (events.Electron.pt > 34)
                                & (abs(events.Electron.eta) < 2.5)
                                & (events.Electron.mvaFall17V2Iso_WP80 > 0.5)
                            ]
                        ),
                        self._campaign,
                        correction_config[self._campaign]["LSF"],
                    ),
                    1.0,
                ),
            )
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

            jetsfs_c[0]["SF"] = getSF(
                smuon_jet.hadronFlavour,
                smuon_jet.btagDeepFlavCvL,
                smuon_jet.btagDeepFlavCvB,
                self._deepjetc_sf,
            )
            jetsfs_c[0]["SFup"] = getSF(
                smuon_jet.hadronFlavour,
                smuon_jet.btagDeepFlavCvL,
                smuon_jet.btagDeepFlavCvB,
                self._deepjetc_sf,
                "TotalUncUp",
            )
            jetsfs_c[0]["SFdn"] = getSF(
                smuon_jet.hadronFlavour,
                smuon_jet.btagDeepFlavCvL,
                smuon_jet.btagDeepFlavCvB,
                self._deepjetc_sf,
                "TotalUncDown",
            )
            jetsfs_b[0]["SF"] = self._deepjetb_sf.eval(
                "central",
                smuon_jet.hadronFlavour,
                abs(smuon_jet.eta),
                smuon_jet.pt,
                discr=smuon_jet.btagDeepFlavB,
            )
            jetsfs_b[0]["SFup"] = self._deepjetb_sf.eval(
                "up_jes",
                smuon_jet.hadronFlavour,
                abs(smuon_jet.eta),
                smuon_jet.pt,
                discr=smuon_jet.btagDeepFlavB,
            )
            jetsfs_b[0]["SFdn"] = self._deepjetb_sf.eval(
                "down_jes",
                smuon_jet.hadronFlavour,
                abs(smuon_jet.eta),
                smuon_jet.pt,
                discr=smuon_jet.btagDeepFlavB,
            )
            csvsfs_c[0]["SF"] = getSF(
                smuon_jet.hadronFlavour,
                smuon_jet.btagDeepCvL,
                smuon_jet.btagDeepCvB,
                self._deepcsvc_sf,
            )
            csvsfs_c[0]["SFup"] = getSF(
                smuon_jet.hadronFlavour,
                smuon_jet.btagDeepCvL,
                smuon_jet.btagDeepCvB,
                self._deepcsvc_sf,
                "TotalUncUp",
            )
            csvsfs_c[0]["SFdn"] = getSF(
                smuon_jet.hadronFlavour,
                smuon_jet.btagDeepCvL,
                smuon_jet.btagDeepCvB,
                self._deepcsvc_sf,
                "TotalUncDown",
            )
            csvsfs_b[0]["SFup"] = self._deepcsvb_sf.eval(
                "up_jes",
                smuon_jet.hadronFlavour,
                abs(smuon_jet.eta),
                smuon_jet.pt,
                discr=smuon_jet.btagDeepB,
            )
            csvsfs_b[0]["SF"] = self._deepcsvb_sf.eval(
                "central",
                smuon_jet.hadronFlavour,
                abs(smuon_jet.eta),
                smuon_jet.pt,
                discr=smuon_jet.btagDeepB,
            )
            csvsfs_b[0]["SFdn"] = self._deepcsvb_sf.eval(
                "down_jes",
                smuon_jet.hadronFlavour,
                abs(smuon_jet.eta),
                smuon_jet.pt,
                discr=smuon_jet.btagDeepB,
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
            smpu = (smuon_jet.partonFlavour == 0) & (smuon_jet.hadronFlavour == 0)
            genflavor = 1 * smpu + smuon_jet.hadronFlavour
            if histname in self.btagDeephists:
                fields = {
                    l: smuon_jet[histname] for l in h.fields if l in dir(smuon_jet)
                }
                if isRealData:
                    h.fill(dataset=dataset, flav=5, **fields, weight=osss)
                else:
                    h.fill(
                        dataset=dataset,
                        flav=genflavor,
                        **fields,
                        weight=weights.weight()[event_level],
                    )
            elif "hl_" in histname:
                fields = {l: shmu[l] for l in h.fields if l in dir(shmu)}
                if isRealData:
                    h.fill(dataset=dataset, **fields)
                else:
                    h.fill(
                        dataset=dataset, **fields, weight=weights.weight()[event_level]
                    )
            elif "sl_" in histname:
                fields = {l: ssmu[l] for l in h.fields if l in dir(ssmu)}
                if isRealData:
                    h.fill(dataset=dataset, **fields)
                else:
                    h.fill(
                        dataset=dataset, **fields, weight=weights.weight()[event_level]
                    )
            elif "met" in histname:
                fields = {
                    l: selev.METFixEE2017[l]
                    for l in h.fields
                    if l in dir(selev.METFixEE2017)
                }
                if isRealData:
                    h.fill(dataset=dataset, **fields)
                else:
                    h.fill(
                        dataset=dataset, **fields, weight=weights.weight()[event_level]
                    )
            elif "mujet_" in histname:
                fields = {l: smuon_jet[l] for l in h.fields if l in dir(smuon_jet)}
                if isRealData:
                    h.fill(dataset=dataset, flav=5, **fields)
                else:
                    h.fill(
                        dataset=dataset,
                        flav=5,
                        **fields,
                        weight=weights.weight()[event_level],
                    )
            elif ["zmass", "zpt", "zeta", "zphi"] == histname:
                fields = {l: sz[l] for l in h.fields if l in dir(sz)}
                if isRealData:
                    h.fill(dataset=dataset, **fields)
                else:
                    h.fill(
                        dataset=dataset, **fields, weight=weights.weight()[event_level]
                    )
            elif ["wmass", "wpt", "weta", "wphi"] == histname:
                fields = {l: sw[l] for l in h.fields if l in dir(sw)}
                if isRealData:
                    h.fill(dataset=dataset, **fields)
                else:
                    h.fill(
                        dataset=dataset, **fields, weight=weights.weight()[event_level]
                    )
            elif "btagDeep" in histname and "0" in histname:
                fields = {
                    l: np.where(smuon_jet[l] < 0, -0.2, smuon_jet[l])
                    for l in h.fields
                    if l in dir(smuon_jet)
                }
                if isRealData:
                    h.fill(dataset=dataset, flav=5, syst="noSF", **fields)
                else:
                    h.fill(
                        dataset=dataset,
                        flav=genflavor,
                        syst="noSF",
                        **fields,
                        weight=weights.weight()[event_level],
                    )
                    for syst in disc_list[histname.replace("_0", "")][0].keys():
                        h.fill(
                            dataset=dataset,
                            flav=genflavor,
                            syst=syst,
                            **fields,
                            weight=weights.weight()[event_level]
                            * disc_list[histname.replace("_0", "")][0][syst],
                        )

        if not isRealData:
            output["drmumu"].fill(
                dataset=dataset,
                drmumu=ssmu.delta_r(shmu),
                weight=weights.weight()[event_level],
            )
            output["hlptratio"].fill(
                dataset=dataset,
                flav=genflavor,
                hlptratio=shmu.pt / smuon_jet.pt,
                weight=weights.weight()[event_level],
            )
            output["slptratio"].fill(
                dataset=dataset,
                flav=genflavor,
                slptratio=ssmu.pt / smuon_jet.pt,
                weight=weights.weight()[event_level],
            )
            output["drjet_hmu"].fill(
                dataset=dataset,
                flav=genflavor,
                drjet_hmu=smuon_jet.delta_r(shmu),
                weight=weights.weight()[event_level],
            )
            output["drjet_smu"].fill(
                dataset=dataset,
                flav=genflavor,
                drjet_smu=smuon_jet.delta_r(ssmu),
                weight=weights.weight()[event_level],
            )
        else:
            output["drmumu"].fill(dataset=dataset, drmumu=ssmu.delta_r(shmu))
            output["hlptratio"].fill(
                dataset=dataset, flav=5, hlptratio=shmu.pt / smuon_jet.pt
            )
            output["slptratio"].fill(
                dataset=dataset, flav=5, slptratio=ssmu.pt / smuon_jet.pt
            )
            output["drjet_hmu"].fill(
                dataset=dataset, flav=5, drjet_hmu=smuon_jet.delta_r(shmu)
            )
            output["drjet_smu"].fill(
                dataset=dataset, flav=5, drjet_smu=smuon_jet.delta_r(ssmu)
            )

        return output

    def postprocess(self, accumulator):
        return accumulator
