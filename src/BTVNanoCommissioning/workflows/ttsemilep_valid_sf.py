import collections
from coffea import hist, processor
import numpy as np
import awkward as ak
from coffea.analysis_tools import Weights
from BTVNanoCommissioning.utils.correction import lumiMasks, muSFs, load_pu, load_BTV
import gc
from BTVNanoCommissioning.helpers.definitions import definitions
from BTVNanoCommissioning.helpers.cTagSFReader import getSF
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
        ntbjet_axis = hist.Bin(
            "ntbjet", r"N b-jets", [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        )
        lmupt_axis = hist.Bin("lmupt", r"Muon pt", 45, 20, 200)
        met_axis = hist.Bin("met", r"Missing ET", 50, 0, 500)
        npv_axis = hist.Bin("npv", r"N-primary vertices", 40, 0, 80)
        nsv_axis = hist.Bin("nsv", r"N-secondary vertices", 10, 0, 20)

        # Jet
        jet_pt_axis = hist.Bin("pt", r"Jet $p_{T}$ [GeV]", 50, 0, 500)
        jet_ptwide_axis = hist.Bin(
            "ptwide",
            r"Jet $p_{T}$ [GeV]",
            [25, 30, 40, 60, 80, 100, 150, 200, 300, 500],
        )
        jet_eta_axis = hist.Bin("eta", r"Jet $\eta$", 25, -2.5, 2.5)
        jet_etawide_axis = hist.Bin(
            "etawide", r"Jet $\eta$", [-2.5, -2.0, -1.5, -0.5, 0.0, 0.5, 1.5, 2.0, 2.5]
        )
        jet_phi_axis = hist.Bin("phi", r"Jet $\phi$", 30, -3, 3)
        jet_mass_axis = hist.Bin("mass", r"Jet $m$ [GeV]", 50, 0, 500)
        jet_dr_axis = hist.Bin("dr", r"Jet $\Delta$R(l,j)", 20, 0, 5)
        ljeta_axis = hist.Bin("ljeta", r"Leading Jet $\eta$", 25, -2.5, 2.5)
        sljeta_axis = hist.Bin("sljeta", r"Leading Jet $\eta$", 25, -2.5, 2.5)
        ssljeta_axis = hist.Bin("ssljeta", r"Leading Jet $\eta$", 25, -2.5, 2.5)
        sssljeta_axis = hist.Bin("sssljeta", r"Leading Jet $\eta$", 25, -2.5, 2.5)
        ljpt_axis = hist.Bin("ljpt", r"Leading jet $p_{T}$ [GeV]", 50, 0, 500)
        sljpt_axis = hist.Bin("sljpt", r"Subleading jet $p_{T}$ [GeV]", 50, 0, 500)
        ssljpt_axis = hist.Bin("ssljpt", r"Subleading jet $p_{T}$ [GeV]", 50, 0, 500)
        sssljpt_axis = hist.Bin("sssljpt", r"Subleading jet $p_{T}$ [GeV]", 50, 0, 500)
        ljdr_axis = hist.Bin("ljdr", "Leading jet $\Delta$R(l,j)", 20, 0, 5)
        sljdr_axis = hist.Bin("sljdr", "Subleading jet $\Delta$R(l,j)", 20, 0, 5)
        ljdr_axis = hist.Bin("ljdr", "Leading jet $\Delta$R(l,j)", 25, 0, 5)
        sljdr_axis = hist.Bin("sljdr", "Subleading jet $\Delta$R(l,j)", 25, 0, 5)
        ssljdr_axis = hist.Bin("ssljdr", "ssSubleading jet $\Delta$R(l,j)", 25, 0, 5)
        sssljdr_axis = hist.Bin("sssljdr", "sssSubleading jet $\Delta$R(l,j)", 25, 0, 5)

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
        for disc, axis in zip(disc_list, btag_axes):
            for i in range(4):
                _hist_sf_dict["%s_%d" % (disc, i)] = hist.Hist(
                    "Counts", dataset_axis, flav_axis, syst_axis, axis
                )

        for deepcsv, axises in zip(bininfo.keys(), btagDeepaxes):
            _hist_btagDeepdict[deepcsv] = hist.Hist(
                "Counts", dataset_axis, flav_axis, axises
            )

        # Generate some histograms dynamically
        _hist_event_dict = {
            "njet": hist.Hist("Counts", dataset_axis, njet_axis),
            "nbjet": hist.Hist("Counts", dataset_axis, nbjet_axis),
            "nbjet_up": hist.Hist("Counts", dataset_axis, nbjet_axis),
            "nbjet_dn": hist.Hist("Counts", dataset_axis, nbjet_axis),
            "ntbjet": hist.Hist("Counts", dataset_axis, ntbjet_axis),
            "nmu": hist.Hist("Counts", dataset_axis, nmu_axis),
            "lmupt": hist.Hist("Counts", dataset_axis, lmupt_axis),
            "ljeta": hist.Hist("Counts", dataset_axis, flav_axis, ljeta_axis),
            "sljeta": hist.Hist("Counts", dataset_axis, flav_axis, sljeta_axis),
            "ssljeta": hist.Hist("Counts", dataset_axis, flav_axis, ssljeta_axis),
            "sssljeta": hist.Hist("Counts", dataset_axis, flav_axis, sssljeta_axis),
            "ljpt": hist.Hist("Counts", dataset_axis, flav_axis, ljpt_axis),
            "sljpt": hist.Hist("Counts", dataset_axis, flav_axis, sljpt_axis),
            "ssljpt": hist.Hist("Counts", dataset_axis, flav_axis, ssljpt_axis),
            "sssljpt": hist.Hist("Counts", dataset_axis, flav_axis, sssljpt_axis),
            "ljdr": hist.Hist("Counts", dataset_axis, flav_axis, ljdr_axis),
            "sljdr": hist.Hist("Counts", dataset_axis, flav_axis, sljdr_axis),
            "ssljdr": hist.Hist("Counts", dataset_axis, flav_axis, ssljdr_axis),
            "sssljdr": hist.Hist("Counts", dataset_axis, flav_axis, sssljdr_axis),
            "met": hist.Hist("Counts", dataset_axis, met_axis),
        }
        self.sf_dict = list(_hist_sf_dict.keys())
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
            (events.Muon.pt > 30)
            & (abs(events.Muon.eta) < 2.4)
            & (events.Muon.tightId > 0.5)
            & (events.Muon.pfRelIso04_all < 0.12)
        ]  #

        event_muon = ak.pad_none(events.Muon, 1, axis=1)

        req_muon = ak.count(event_muon.pt, axis=1) == 1
        ## Jet cuts

        event_jet = events.Jet[
            (events.Jet.pt > 25)
            & (abs(events.Jet.eta) <= 2.4)
            & (events.Jet.puId > 0)
            & (events.Jet.jetId > 0)
            & (ak.all(events.Jet.metric_table(events.Muon) > 0.4, axis=2))
        ]

        req_jets = ak.num(event_jet) >= 4
        if hasattr(events, "METFixEE2017"):
            req_MET = events.METFixEE2017.pt > 50
        else:
            req_MET = events.MET.pt > 50

        event_level = req_trig & req_jets & req_muon & req_MET & req_lumi
        # Selected
        selev = events[event_level]
        if hasattr(events, "METFixEE2017"):
            MET = selev.METFixEE2017.pt
        else:
            MET = selev.MET.pt
        #########

        # Per muon
        mu_eta = abs(selev.Muon.eta) < 2.4
        mu_pt = selev.Muon.pt > 30
        mu_idiso = (selev.Muon.tightId > 0.5) & (selev.Muon.pfRelIso04_all < 0.12)
        mu_level = mu_eta & mu_pt & mu_idiso

        smu = selev.Muon[mu_level]
        # Per jet : https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID
        jet_eta = abs(selev.Jet.eta) <= 2.4
        jet_pt = selev.Jet.pt > 25
        jet_pu = (selev.Jet.puId > 0) & (selev.Jet.jetId > 0)
        jet_dr = ak.all(selev.Jet.metric_table(smu) > 0.4, axis=2)

        jet_level = jet_pu & jet_eta & jet_pt & jet_dr

        # b-tag twiki : https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
        bjet_disc = selev.Jet.btagDeepB > 0.4941  # L=0.0494, M=0.2770, T=0.7264
        bjet_level = jet_level & bjet_disc

        sjets = selev.Jet[jet_level]
        sel_jets = sjets
        sjets = sjets[:, :4]
        sbjets = selev.Jet[bjet_level]
        if not isRealData:
            weights.add(
                "lep1sf",
                np.where(
                    event_level,
                    muSFs(
                        ak.firsts(
                            events.Muon[
                                (events.Muon.pt > 30)
                                & (abs(events.Muon.eta) < 2.4)
                                & (events.Muon.tightId > 0.5)
                                & (events.Muon.pfRelIso04_all < 0.12)
                            ]
                        ),
                        self._campaign,
                        correction_config[self._campaign]["LSF"],
                    ),
                    1.0,
                ),
            )
        if not isRealData:
            stbjets = sbjets[sbjets.hadronFlavour == 5]
        if isRealData:
            genflavor = ak.zeros_like(sjets.pt)
        else:
            par_flav = (sjets.partonFlavour == 0) & (sjets.hadronFlavour == 0)
            genflavor = sjets.hadronFlavour + 1 * par_flav
            jetsfs_c = collections.defaultdict(dict)
            jetsfs_b = collections.defaultdict(dict)
            csvsfs_c = collections.defaultdict(dict)
            csvsfs_b = collections.defaultdict(dict)
            for i in range(4):
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
            output["nmu"].fill(dataset=dataset, nmu=flatten(ak.num(smu)))
            output["lmupt"].fill(dataset=dataset, lmupt=flatten((smu.pt)))
            output["ljpt"].fill(dataset=dataset, flav=0, ljpt=flatten(sjets[:, 0].pt))
            output["sljpt"].fill(dataset=dataset, flav=0, sljpt=flatten(sjets[:, 1].pt))
            output["ssljpt"].fill(
                dataset=dataset, flav=0, ssljpt=flatten(sjets[:, 2].pt)
            )
            output["sssljpt"].fill(
                dataset=dataset, flav=0, sssljpt=flatten(sjets[:, 3].pt)
            )
            output["ljeta"].fill(
                dataset=dataset, flav=0, ljeta=flatten(sjets[:, 0].eta)
            )
            output["sljeta"].fill(
                dataset=dataset, flav=0, sljeta=flatten(sjets[:, 1].eta)
            )
            output["ssljeta"].fill(
                dataset=dataset, flav=0, ssljeta=flatten(sjets[:, 2].eta)
            )
            output["sssljeta"].fill(
                dataset=dataset, flav=0, sssljeta=flatten(sjets[:, 3].eta)
            )
            output["ljdr"].fill(
                dataset=dataset, flav=0, ljdr=flatten(sjets[:, 0].delta_r(smu))
            )
            output["sljdr"].fill(
                dataset=dataset, flav=0, sljdr=flatten(sjets[:, 1].delta_r(smu))
            )
            output["ssljdr"].fill(
                dataset=dataset, flav=0, ssljdr=flatten(sjets[:, 2].delta_r(smu))
            )
            output["sssljdr"].fill(
                dataset=dataset, flav=0, sssljdr=flatten(sjets[:, 3].delta_r(smu))
            )
            output["met"].fill(dataset=dataset, met=flatten((MET)))
        else:
            output["njet"].fill(
                dataset=dataset,
                njet=flatten(ak.num(sel_jets)),
                weight=weights.weight()[event_level],
            )
            output["ntbjet"].fill(
                dataset=dataset,
                ntbjet=flatten(ak.num(stbjets)),
                weight=weights.weight()[event_level],
            )
            output["nmu"].fill(
                dataset=dataset,
                nmu=flatten(ak.num(smu)),
                weight=weights.weight()[event_level],
            )
            output["lmupt"].fill(
                dataset=dataset,
                lmupt=flatten((smu.pt)),
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
            output["ssljpt"].fill(
                dataset=dataset,
                flav=flatten(sjets[:, 2].hadronFlavour),
                ssljpt=flatten(sjets[:, 2].pt),
                weight=weights.weight()[event_level],
            )
            output["sssljpt"].fill(
                dataset=dataset,
                flav=flatten(sjets[:, 3].hadronFlavour),
                sssljpt=flatten(sjets[:, 3].pt),
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
            output["ssljeta"].fill(
                dataset=dataset,
                flav=flatten(sjets[:, 2].hadronFlavour),
                ssljeta=flatten(sjets[:, 2].eta),
                weight=weights.weight()[event_level],
            )
            output["sssljeta"].fill(
                dataset=dataset,
                flav=flatten(sjets[:, 3].hadronFlavour),
                sssljeta=flatten(sjets[:, 3].eta),
                weight=weights.weight()[event_level],
            )

            output["ljdr"].fill(
                dataset=dataset,
                flav=flatten(sjets[:, 0].hadronFlavour),
                ljdr=flatten(sjets[:, 0].delta_r(smu)),
                weight=weights.weight()[event_level],
            )
            output["sljdr"].fill(
                dataset=dataset,
                flav=flatten(sjets[:, 1].hadronFlavour),
                sljdr=flatten(sjets[:, 1].delta_r(smu)),
                weight=weights.weight()[event_level],
            )
            output["ssljdr"].fill(
                dataset=dataset,
                flav=flatten(sjets[:, 2].hadronFlavour),
                ssljdr=flatten(sjets[:, 2].delta_r(smu)),
                weight=weights.weight()[event_level],
            )
            output["sssljdr"].fill(
                dataset=dataset,
                flav=flatten(sjets[:, 3].hadronFlavour),
                sssljdr=flatten(sjets[:, 3].delta_r(smu)),
                weight=weights.weight()[event_level],
            )
            output["met"].fill(
                dataset=dataset,
                met=flatten((MET)),
                weight=weights.weight()[event_level],
            )

        return output

    def postprocess(self, accumulator):
        return accumulator
