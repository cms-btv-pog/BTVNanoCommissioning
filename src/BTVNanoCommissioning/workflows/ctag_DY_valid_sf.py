import collections, numpy as np, awkward as ak

from coffea import processor
from coffea.analysis_tools import Weights

from BTVNanoCommissioning.utils.correction import (
    load_lumi,
    eleSFs,
    muSFs,
    load_pu,
    load_BTV,
    load_jmefactory,
)
from BTVNanoCommissioning.utils.AK4_parameters import correction_config
from BTVNanoCommissioning.helpers.func import flatten, update
from BTVNanoCommissioning.helpers.update_branch import missing_branch, add_jec
from BTVNanoCommissioning.helpers.cTagSFReader import getSF
from BTVNanoCommissioning.utils.histogrammer import histogrammer
from BTVNanoCommissioning.utils.selection import jet_id, mu_idiso, ele_mvatightid


class NanoProcessor(processor.ProcessorABC):
    def __init__(self, year="2017", campaign="Rereco17_94X", isCorr=True, isJERC=False):
        self._year = year
        self._campaign = campaign
        self.isJERC = isJERC
        self.isCorr = isCorr
        self.lumiMask = load_lumi(correction_config[self._campaign]["lumiMask"])

        ## Load corrections
        if isCorr:
            if "BTV" in correction_config[self._campaign].keys():
                self._deepjetc_sf = load_BTV(
                    self._campaign, correction_config[self._campaign]["BTV"], "DeepJetC"
                )
                self._deepjetb_sf = load_BTV(
                    self._campaign, correction_config[self._campaign]["BTV"], "DeepJetB"
                )
                self._deepcsvc_sf = load_BTV(
                    self._campaign, correction_config[self._campaign]["BTV"], "DeepCSVC"
                )
                self._deepcsvb_sf = load_BTV(
                    self._campaign, correction_config[self._campaign]["BTV"], "DeepCSVB"
                )
            if "PU" in correction_config[self._campaign].keys():
                self._pu = load_pu(
                    self._campaign, correction_config[self._campaign]["PU"]
                )

        if isJERC:
            self._jet_factory = load_jmefactory(
                self._campaign, correction_config[self._campaign]["JME"]
            )
        _hist_event_dict = histogrammer("ctag_DY_sf")
        self.make_output = lambda: {
            "sumw": processor.defaultdict_accumulator(float),
            **_hist_event_dict,
        }

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.make_output()
        dataset = events.metadata["dataset"]
        isRealData = not hasattr(events, "genWeight")
        events = missing_branch(events)
        weights = Weights(len(events), storeIndividual=True)
        if self.isJERC:
            add_jec(events, self._campaign, self._jet_factory)
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
        triggers = ["Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8"]
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
        dilep_mu = events.Muon[(events.Muon.pt > 12) & mu_idiso(events, self._campaign)]
        ## Electron cuts
        dilep_ele = events.Electron[
            (events.Electron.pt > 15) & ele_mvatightid(events, self._campaign)
        ]
        ## dilepton
        pos_dilep = dilep_mu[dilep_mu.charge > 0]
        neg_dilep = dilep_mu[dilep_mu.charge < 0]
        req_dilep = (
            (ak.num(pos_dilep.pt) >= 1)
            & (ak.num(neg_dilep.pt) >= 1)
            & (ak.num(dilep_mu.charge) >= 2)
            & (ak.num(dilep_ele.charge) == 0)
        )
        pos_dilep = ak.pad_none(pos_dilep, 1, axis=1)
        neg_dilep = ak.pad_none(neg_dilep, 1, axis=1)

        dilep_mass = pos_dilep[:, 0] + neg_dilep[:, 0]
        req_dilepmass = (
            (dilep_mass.mass > 81) & (dilep_mass.mass < 101) & (dilep_mass.pt > 15)
        )

        ## Jet cuts
        event_jet = events.Jet[
            jet_id(events, self._campaign)
            & (ak.all(events.Jet.metric_table(pos_dilep[:, 0]) > 0.4, axis=2))
            & (ak.all(events.Jet.metric_table(neg_dilep[:, 0]) > 0.4, axis=2))
        ]
        req_jets = ak.num(event_jet.pt) >= 1
        event_jet = ak.pad_none(event_jet, 1, axis=1)

        ## store jet index for PFCands, create mask on the jet index
        jetindx = ak.mask(
            ak.local_index(events.Jet.pt),
            (
                jet_id(events, self._campaign)
                & (ak.all(events.Jet.metric_table(pos_dilep[:, 0]) > 0.4, axis=2))
                & (ak.all(events.Jet.metric_table(neg_dilep[:, 0]) > 0.4, axis=2))
            )
            == 1,
        )
        jetindx = ak.pad_none(jetindx, 1)
        jetindx = jetindx[:, 0]

        event_level = ak.fill_none(
            req_lumi & req_trig & req_dilep & req_dilepmass & req_jets, False
        )

        ####################
        # Selected objects #
        ####################
        sposmu = pos_dilep[event_level]
        sposmu = sposmu[:, 0]
        snegmu = neg_dilep[event_level]
        snegmu = snegmu[:, 0]
        sz = sposmu + snegmu
        sjets = event_jet[event_level]
        njet = ak.count(sjets.pt, axis=1)
        # Find the PFCands associate with selected jets. Search from jetindex->JetPFCands->PFCand
        if self._campaign != "Rereco17_94X":
            spfcands = events[event_level].PFCands[
                events[event_level]
                .JetPFCands[
                    events[event_level].JetPFCands.jetIdx == jetindx[event_level]
                ]
                .pFCandsIdx
            ]
        sel_jet = sjets[:, 0]
        ####################
        # Weight & Geninfo #
        ####################
        if not isRealData:
            weights.add("genweight", events.genWeight)
        if not isRealData and self.isCorr:
            if "PU" in correction_config[self._campaign].keys():
                if self._campaign == "Rereco17_94X":
                    puname = f"{self._year}_pileupweight"
                else:
                    puname = "PU"
                weights.add("puweight", self._pu[puname](events.Pileup.nTrueInt))
            if "LSF" in correction_config[self._campaign].keys():
                weights.add(
                    "lep1sf",
                    np.where(
                        event_level,
                        muSFs(
                            ak.firsts(pos_dilep),
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
                            ak.firsts(neg_dilep),
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
            if self.isCorr and "BTV" in correction_config[self._campaign].keys():
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
            elif "PFCands" in histname and self._campaign != "Rereco17_94X":
                h.fill(
                    flatten(ak.broadcast_arrays(genflavor[:, 0], spfcands["pt"])[0]),
                    flatten(spfcands[histname.replace("PFCands_", "")]),
                    weight=flatten(
                        ak.broadcast_arrays(
                            weights.weight()[event_level], spfcands["pt"]
                        )[0]
                    ),
                )
            elif "posl_" in histname and histname.replace("posl_", "") in sposmu.fields:
                h.fill(
                    flatten(sposmu[histname.replace("posl_", "")]),
                    weight=weights.weight()[event_level],
                )
            elif "negl_" in histname and histname.replace("negl_", "") in snegmu.fields:
                h.fill(
                    flatten(snegmu[histname.replace("negl_", "")]),
                    weight=weights.weight()[event_level],
                )

            elif "jet_" in histname:
                h.fill(
                    genflavor[:, 0],
                    sel_jet[histname.replace("jet_", "")],
                    weight=weights.weight()[event_level],
                )
            elif "btagDeep" in histname and "0" in histname:

                h.fill(
                    flav=genflavor[:, 0],
                    syst="noSF",
                    discr=np.where(
                        sel_jet[histname.replace("_0", "")] < 0,
                        -0.2,
                        sel_jet[histname.replace("_0", "")],
                    ),
                    weight=weights.weight()[event_level],
                )
                if (
                    not isRealData
                    and self.isCorr
                    and "BTV" in correction_config[self._campaign].keys()
                    and "_b" not in histname
                    and "_bb" not in histname
                    and "_lepb" not in histname
                ):
                    for syst in disc_list[histname.replace("_0", "")][0].keys():
                        h.fill(
                            flav=genflavor[:, 0],
                            syst=syst,
                            discr=np.where(
                                sel_jet[histname.replace("_0", "")] < 0,
                                -0.2,
                                sel_jet[histname.replace("_0", "")],
                            ),
                            weight=weights.weight()[event_level]
                            * disc_list[histname.replace("_0", "")][0][syst],
                        )
        output["njet"].fill(njet, weight=weights.weight()[event_level])
        output["dr_mumu"].fill(
            dr=snegmu.delta_r(sposmu), weight=weights.weight()[event_level]
        )
        output["z_pt"].fill(flatten(sz.pt), weight=weights.weight()[event_level])
        output["z_eta"].fill(flatten(sz.eta), weight=weights.weight()[event_level])
        output["z_phi"].fill(flatten(sz.phi), weight=weights.weight()[event_level])
        output["z_mass"].fill(flatten(sz.mass), weight=weights.weight()[event_level])
        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator
