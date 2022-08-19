import numpy as np
import collections


from coffea import hist, processor
import awkward as ak
from coffea.analysis_tools import Weights

from BTVNanoCommissioning.utils.correction import (
    lumiMasks,
    muSFs,
    load_pu,
    load_BTV,
)
from BTVNanoCommissioning.utils.AK4_parameters import correction_config
from BTVNanoCommissioning.helpers.definitions import definitions
from BTVNanoCommissioning.helpers.cTagSFReader import getSF


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
        lpt_axis = hist.Bin("pt", r"leading Lepton pt", 40, 0, 200)
        leta_axis = hist.Bin("eta", r"leading Lepton $\eta$", 25, -2.5, 2.5)
        lphi_axis = hist.Bin("phi", r"l eadingLepton $\phi$", 30, -3, 3)
        liso_axis = hist.Bin("pfRelIso04_all", r"leading Muon Rel. Iso", 40, 0, 4.0)
        lpt_ratio_axis = hist.Bin(
            "lptratio", r"leading $\mu$/Jet $p_{T}$ [GeV]", 40, 0, 1
        )
        slpt_axis = hist.Bin("pt", r"subleading Lepton pt", 40, 0, 40)
        sleta_axis = hist.Bin("eta", r"subleading Lepton $\eta$", 25, -2.5, 2.5)
        slphi_axis = hist.Bin("phi", r"subleading Lepton $\phi$", 30, -3, 3)
        sliso_axis = hist.Bin(
            "pfRelIso04_all", r"subleading Soft Muon Rel. Iso", 40, 0, 4.0
        )
        slpt_ratio_axis = hist.Bin(
            "slptratio", r"subleading Soft $\mu$/Jet $p_{T}$ [GeV]", 40, 0, 1
        )

        ## Soft Muons
        softpt_axis = hist.Bin("pt", r"Soft lepton pt", 25, 0, 25)
        softeta_axis = hist.Bin("eta", r"Soft Lepton $\eta$", 25, -2.5, 2.5)
        softphi_axis = hist.Bin("phi", r"Soft Lepton $\phi$", 30, -3, 3)
        softiso_axis = hist.Bin("pfRelIso04_all", r"Soft Muon Rel. Iso", 40, 0, 4.0)
        softpt_ratio_axis = hist.Bin(
            "softptratio", r"Soft $\mu$/Jet $p_{T}$ [GeV]", 40, 0, 1
        )

        ## Z/W
        zmass_axis = hist.Bin("zmass", r"Z Mass", 25, 50, 100)
        zpt_axis = hist.Bin("zpt", r"Z $p_{T}$", 25, 0, 100)
        zeta_axis = hist.Bin("zeta", r"Z $\eta$", 25, -2.5, 2.5)
        zphi_axis = hist.Bin("zphi", r"Z $\phi$", 30, -3, 3)

        ## MET
        met_axis = hist.Bin("pt", r"MET $p_{T}$", 50, 0, 500)
        metphi_axis = hist.Bin("phi", r"met $\phi$", 30, -3, 3)

        ## Muon jets
        lmujet_pt_axis = hist.Bin("pt", r"Jet $p_{T}$ [GeV]", 50, 0, 500)
        lmujet_eta_axis = hist.Bin("eta", r"Jet $\eta$", 25, -2.5, 2.5)
        lmujet_phi_axis = hist.Bin("phi", r"Jet $\phi$", 30, -3, 3)
        lmujet_mass_axis = hist.Bin("mass", r"Jet $m$ [GeV]", 50, 0, 500)
        slmujet_pt_axis = hist.Bin("pt", r"Jet $p_{T}$ [GeV]", 50, 0, 500)
        slmujet_eta_axis = hist.Bin("eta", r"Jet $\eta$", 25, -2.5, 2.5)
        slmujet_phi_axis = hist.Bin("phi", r"Jet $\phi$", 30, -3, 3)
        slmujet_mass_axis = hist.Bin("mass", r"Jet $m$ [GeV]", 50, 0, 500)
        dr_lmujetsmu_axis = hist.Bin(
            "dr_lmujetsmu", r"$\Delta$R($\mu_{soft}$,j)", 25, 0, 5
        )
        dr_slmujetsmu_axis = hist.Bin(
            "dr_slmujetsmu", r"$\Delta$R($\mu_{soft}$,j)", 25, 0, 5
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
            for i in range(3):
                _hist_sf_dict["%s_%d" % (disc, i)] = hist.Hist(
                    "Counts", dataset_axis, flav_axis, syst_axis, axis
                )
        for deepcsv, axises in zip(bininfo.keys(), btagDeepaxes):
            _hist_btagDeepdict["%s" % (deepcsv)] = hist.Hist(
                "Counts", dataset_axis, flav_axis, axises
            )

        _hist_event_dict = {
            "njet": hist.Hist("Counts", dataset_axis, njet_axis),
            "hl_pt": hist.Hist("Counts", dataset_axis, lpt_axis),
            "hl_eta": hist.Hist("Counts", dataset_axis, leta_axis),
            "hl_phi": hist.Hist("Counts", dataset_axis, lphi_axis),
            "hl_pfRelIso04_all": hist.Hist("Counts", dataset_axis, liso_axis),
            "sl_pt": hist.Hist("Counts", dataset_axis, slpt_axis),
            "sl_eta": hist.Hist("Counts", dataset_axis, sleta_axis),
            "sl_phi": hist.Hist("Counts", dataset_axis, slphi_axis),
            "sl_pfRelIso04_all": hist.Hist("Counts", dataset_axis, sliso_axis),
            "hlptratio": hist.Hist("Counts", dataset_axis, flav_axis, lpt_ratio_axis),
            "slptratio": hist.Hist("Counts", dataset_axis, flav_axis, slpt_ratio_axis),
            "soft_lpt": hist.Hist("Counts", dataset_axis, softpt_axis),
            "soft_leta": hist.Hist("Counts", dataset_axis, softeta_axis),
            "soft_lphi": hist.Hist("Counts", dataset_axis, softphi_axis),
            "soft_liso": hist.Hist("Counts", dataset_axis, softiso_axis),
            "softlptratio": hist.Hist(
                "Counts", dataset_axis, flav_axis, softpt_ratio_axis
            ),
            "zmass": hist.Hist("Counts", dataset_axis, zmass_axis),
            "zpt": hist.Hist("Counts", dataset_axis, zpt_axis),
            "zeta": hist.Hist("Counts", dataset_axis, zeta_axis),
            "zphi": hist.Hist("Counts", dataset_axis, zphi_axis),
            "metpt": hist.Hist("Counts", dataset_axis, met_axis),
            "metphi": hist.Hist("Counts", dataset_axis, metphi_axis),
            "lmujet_pt": hist.Hist("Counts", dataset_axis, flav_axis, lmujet_pt_axis),
            "lmujet_eta": hist.Hist("Counts", dataset_axis, flav_axis, lmujet_eta_axis),
            "lmujet_phi": hist.Hist("Counts", dataset_axis, flav_axis, lmujet_phi_axis),
            "lmujet_mass": hist.Hist(
                "Counts", dataset_axis, flav_axis, lmujet_mass_axis
            ),
            "slmujet_pt": hist.Hist("Counts", dataset_axis, flav_axis, slmujet_pt_axis),
            "slmujet_eta": hist.Hist(
                "Counts", dataset_axis, flav_axis, slmujet_eta_axis
            ),
            "slmujet_phi": hist.Hist(
                "Counts", dataset_axis, flav_axis, slmujet_phi_axis
            ),
            "slmujet_mass": hist.Hist(
                "Counts", dataset_axis, flav_axis, slmujet_mass_axis
            ),
            "dr_lmujetsmu": hist.Hist(
                "Counts", dataset_axis, flav_axis, dr_lmujetsmu_axis
            ),
            "dr_slmujetsmu": hist.Hist(
                "Counts", dataset_axis, flav_axis, dr_slmujetsmu_axis
            ),
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
        triggers = [
            # "HLT_IsoMu24",
            # "HLT_IsoMu27",
            "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8"
        ]

        trig_arrs = [events.HLT[_trig.strip("HLT_")] for _trig in triggers]
        req_trig = np.zeros(len(events), dtype="bool")
        for t in trig_arrs:
            req_trig = req_trig | t

        ############
        # Event level
        ## Muon cuts
        # muon twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
        iso_muon = events.Muon[
            (events.Muon.pt > 12)
            & (abs(events.Muon.eta) < 2.4)
            & (events.Muon.tightId > 0.5)
            & (events.Muon.pfRelIso04_all <= 0.15)
        ]

        iso_muon = ak.pad_none(iso_muon, 2)

        req_muon = (ak.count(iso_muon.pt, axis=1) == 2) & (iso_muon[:, 0].pt > 20)
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
        req_dilepveto = ak.count(dilep_ele.pt, axis=1) == 0

        req_MET = events.MET.pt > 40

        ## Jet cuts
        event_jet = events.Jet[
            (events.Jet.pt > 20)
            & (abs(events.Jet.eta) <= 2.5)
            & (events.Jet.jetId >= 5)
            & ((events.Jet.pt > 50) | (events.Jet.puId >= 7))
        ]
        req_jets = ak.count(event_jet.pt, axis=1) >= 2

        ## Soft Muon cuts

        soft_muon = events.Muon[
            (events.Muon.pt < 25)
            & (abs(events.Muon.eta) < 2.4)
            & (events.Muon.tightId > 0.5)
            & (events.Muon.pfRelIso04_all > 0.2)
        ]
        req_softmu = ak.count(soft_muon.pt, axis=1) >= 1

        mu_jet = event_jet[
            (ak.all(event_jet.metric_table(soft_muon) <= 0.4, axis=2))
            & ((event_jet.muonIdx1 != -1) | (event_jet.muonIdx2 != -1))
        ]
        req_mujet = ak.count(mu_jet.pt, axis=1) >= 1

        # dilepton mass

        dilep_mass = iso_muon[:, 0] + iso_muon[:, 1]
        req_dilepmass = (dilep_mass.mass > 12.0) & (
            (dilep_mass.mass < 75) | (dilep_mass.mass > 105)
        )

        event_level = (
            req_trig
            & req_lumi
            & req_muon
            & req_dilepveto
            & req_dilepmass
            & req_MET
            & req_jets
            & req_softmu
            & req_mujet
        )

        if len(event_level) > 0:
            event_level = ak.fill_none(event_level, False)
        # Selected
        selev = events[event_level]
        #########

        ## Hard Muon
        shmu = selev.Muon[
            (selev.Muon.pt > 12)
            & (abs(selev.Muon.eta) < 2.4)
            & (selev.Muon.tightId > 0.5)
            & (selev.Muon.pfRelIso04_all <= 0.15)
        ]

        # print(ak.num(shmu)

        # shmu=shmu[:,:1]
        ## Soft Muon
        ssmu = selev.Muon[
            (selev.Muon.pt < 25)
            & (abs(selev.Muon.eta) < 2.4)
            & (selev.Muon.tightId > 0.5)
            & (selev.Muon.pfRelIso04_all > 0.2)
        ]
        softmu0 = ssmu[:, 0]
        sz = shmu[:, 0] + shmu[:, 1]
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
                                (events.Muon.pt < 25)
                                & (abs(events.Muon.eta) < 2.4)
                                & (events.Muon.tightId > 0.5)
                                & (events.Muon.pfRelIso04_all > 0.2)
                            ]
                        ),
                        self._campaign,
                        correction_config[self._campaign]["LSF"],
                    ),
                    1.0,
                ),
            )
        isomu0 = shmu[:, 0]
        isomu1 = shmu[:, 1]

        ## Jets

        sjets = selev.Jet[
            (selev.Jet.pt > 20)
            & (abs(selev.Jet.eta) <= 2.5)
            & (selev.Jet.jetId >= 5)
            & ((selev.Jet.pt > 50) | (selev.Jet.puId >= 7))
        ]

        ## Muon Jet
        smuon_jet = sjets[
            (ak.all(sjets.metric_table(ssmu) <= 0.4, axis=2))
            & ((sjets.muonIdx1 != -1) | (sjets.muonIdx2 != -1))
        ]
        smuon_jet = smuon_jet[:, 0]

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
            if all(i > 1 for i in njet):
                jetsfs_c[1]["SF"] = getSF(
                    sjets[:, 1].hadronFlavour,
                    sjets[:, 1].btagDeepFlavCvL,
                    sjets[:, 1].btagDeepFlavCvB,
                    self._deepjetc_sf,
                )
                jetsfs_c[1]["SFup"] = getSF(
                    sjets[:, 1].hadronFlavour,
                    sjets[:, 1].btagDeepFlavCvL,
                    sjets[:, 1].btagDeepFlavCvB,
                    self._deepjetc_sf,
                    "TotalUncUp",
                )
                jetsfs_c[1]["SFdn"] = getSF(
                    smuon_jet.hadronFlavour,
                    smuon_jet.btagDeepFlavCvL,
                    smuon_jet.btagDeepFlavCvB,
                    self._deepjetc_sf,
                    "TotalUncDown",
                )
                jetsfs_b[1]["SF"] = self._deepjetb_sf.eval(
                    "central",
                    sjets[:, 1].hadronFlavour,
                    abs(sjets[:, 1].eta),
                    sjets[:, 1].pt,
                    discr=sjets[:, 1].btagDeepFlavB,
                )
                jetsfs_b[1]["SFup"] = self._deepjetb_sf.eval(
                    "up_jes",
                    sjets[:, 1].hadronFlavour,
                    abs(sjets[:, 1].eta),
                    sjets[:, 1].pt,
                    discr=sjets[:, 1].btagDeepFlavB,
                )
                jetsfs_b[1]["SFdn"] = self._deepjetb_sf.eval(
                    "down_jes",
                    sjets[:, 1].hadronFlavour,
                    abs(sjets[:, 1].eta),
                    sjets[:, 1].pt,
                    discr=sjets[:, 1].btagDeepFlavB,
                )
                csvsfs_c[1]["SF"] = getSF(
                    smuon_jet.hadronFlavour,
                    smuon_jet.btagDeepCvL,
                    smuon_jet.btagDeepCvB,
                    self._deepcsvc_sf,
                )
                csvsfs_c[1]["SFup"] = getSF(
                    smuon_jet.hadronFlavour,
                    smuon_jet.btagDeepCvL,
                    smuon_jet.btagDeepCvB,
                    self._deepcsvc_sf,
                    "TotalUncUp",
                )
                csvsfs_c[1]["SFdn"] = getSF(
                    smuon_jet.hadronFlavour,
                    smuon_jet.btagDeepCvL,
                    smuon_jet.btagDeepCvB,
                    self._deepcsvc_sf,
                    "TotalUncDown",
                )
                csvsfs_b[1]["SF"] = self._deepcsvb_sf.eval(
                    "central",
                    sjets[:, 1].hadronFlavour,
                    abs(sjets[:, 1].eta),
                    sjets[:, 1].pt,
                    discr=sjets[:, 1].btagDeepB,
                )
                csvsfs_b[1]["SFup"] = self._deepcsvb_sf.eval(
                    "up_jes",
                    sjets[:, 1].hadronFlavour,
                    abs(sjets[:, 1].eta),
                    sjets[:, 1].pt,
                    discr=sjets[:, 1].btagDeepB,
                )
                csvsfs_b[1]["SFdn"] = self._deepcsvb_sf.eval(
                    "down_jes",
                    sjets[:, 1].hadronFlavour,
                    abs(sjets[:, 1].eta),
                    sjets[:, 1].pt,
                    discr=sjets[:, 1].btagDeepB,
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
            if not isRealData:
                smpu = (smuon_jet.partonFlavour == 0) & (smuon_jet.hadronFlavour == 0)
                smflav = 1 * smpu + smuon_jet.hadronFlavour
            if histname in self.btagDeephists:
                fields = {
                    l: ak.flatten(sjets[histname]) for l in h.fields if l in dir(sjets)
                }
                wei = ak.flatten(
                    ak.broadcast_arrays(weights.weight()[event_level], sjets["pt"])[0]
                )
                if isRealData:
                    h.fill(dataset=dataset, flav=5, **fields)
                else:
                    h.fill(
                        dataset=dataset,
                        flav=ak.flatten(sjets.hadronFlavour),
                        **fields,
                        weight=wei,
                    )
            elif "hl_" in histname:
                fields = {l: isomu0[l] for l in h.fields if l in dir(isomu0)}
                if isRealData:
                    h.fill(dataset=dataset, **fields)
                else:
                    h.fill(
                        dataset=dataset, **fields, weight=weights.weight()[event_level]
                    )
            elif "sl_" in histname:
                fields = {l: isomu1[l] for l in h.fields if l in dir(isomu1)}
                if isRealData:
                    h.fill(dataset=dataset, **fields)
                else:
                    h.fill(
                        dataset=dataset, **fields, weight=weights.weight()[event_level]
                    )
            elif "soft_l" in histname:
                fields = {l: softmu0[l] for l in h.fields if l in dir(softmu0)}
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
            elif "lmujet_" in histname:
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
            elif "btagDeep" in histname and "0" in histname:
                fields = {
                    l: np.where(smuon_jet[l] < 0, -0.2, smuon_jet[l])
                    for l in h.fields
                    if l in dir(smuon_jet)
                }
                if isRealData:
                    h.fill(dataset=dataset, flav=5, syst="noSF", **fields)
                else:
                    smpu = (smuon_jet.partonFlavour == 0) & (
                        smuon_jet.hadronFlavour == 0
                    )
                    smflav = 1 * smpu + smuon_jet.hadronFlavour
                    h.fill(
                        dataset=dataset,
                        flav=smflav,
                        syst="noSF",
                        **fields,
                        weight=weights.weight()[event_level],
                    )
                    for syst in disc_list[histname.replace("_0", "")][0].keys():
                        h.fill(
                            dataset=dataset,
                            flav=smflav,
                            syst=syst,
                            **fields,
                            weight=weights.weight()[event_level]
                            * disc_list[histname.replace("_0", "")][0][syst],
                        )
            elif (
                "btagDeep" in histname and "1" in histname and all(i > 1 for i in njet)
            ):
                sljets = sjets[:, 1]
                fields = {
                    l: np.where(sljets[l] < 0, -0.2, sljets[l])
                    for l in h.fields
                    if l in dir(sljets)
                }
                if isRealData:
                    h.fill(dataset=dataset, flav=5, syst="noSF", **fields)
                else:
                    h.fill(
                        dataset=dataset,
                        flav=genflavor[:, 1],
                        syst="noSF",
                        **fields,
                        weight=weights.weight()[event_level],
                    )
                    for syst in disc_list[histname.replace("_1", "")][1].keys():
                        h.fill(
                            dataset=dataset,
                            flav=genflavor[:, 1],
                            syst=syst,
                            **fields,
                            weight=weights.weight()[event_level]
                            * disc_list[histname.replace("_1", "")][1][syst],
                        )

        if not isRealData:
            output["hlptratio"].fill(
                dataset=dataset,
                flav=sjets[:, 0].hadronFlavour,
                lptratio=isomu0.pt / sjets[:, 0].pt,
                weight=weights.weight()[event_level],
            )
            output["softlptratio"].fill(
                dataset=dataset,
                flav=smflav,
                softptratio=softmu0.pt / smuon_jet.pt,
                weight=weights.weight()[event_level],
            )
            output["dr_lmujetsmu"].fill(
                dataset=dataset,
                flav=smflav,
                dr_lmujetsmu=smuon_jet.delta_r(softmu0),
                weight=weights.weight()[event_level],
            )
        else:
            output["hlptratio"].fill(
                dataset=dataset, flav=5, lptratio=isomu0.pt / smuon_jet.pt
            )
            output["softlptratio"].fill(
                dataset=dataset, flav=5, softptratio=ssmu[:, 0].pt / smuon_jet.pt
            )
            output["dr_lmujetsmu"].fill(
                dataset=dataset, flav=5, dr_lmujetsmu=smuon_jet.delta_r(ssmu[:, 0])
            )
        return output

    def postprocess(self, accumulator):
        return accumulator
