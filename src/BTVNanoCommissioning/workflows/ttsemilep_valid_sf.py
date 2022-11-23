import collections, gc
import numpy as np, awkward as ak

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
from BTVNanoCommissioning.helpers.func import flatten, update
from BTVNanoCommissioning.helpers.update_branch import missing_branch, add_jec
from BTVNanoCommissioning.helpers.cTagSFReader import getSF
from BTVNanoCommissioning.utils.AK4_parameters import correction_config
from BTVNanoCommissioning.utils.histogrammer import histogrammer
from BTVNanoCommissioning.utils.selection import jet_id, btag_mu_idiso, ele_mvatightid


class NanoProcessor(processor.ProcessorABC):
    # Define histograms
    def __init__(self, year="2017", campaign="Rereco17_94X", isCorr=True, isJERC=False):
        self._year = year
        self._campaign = campaign
        self.isCorr = isCorr
        self.isJERC = isJERC
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
        _hist_event_dict = histogrammer("ttsemilep_sf")
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
        triggers = ["IsoMu24"]
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
        events.Muon = events.Muon[
            (events.Muon.pt > 30) & btag_mu_idiso(events, self._campaign)
        ]
        event_muon = ak.pad_none(events.Muon, 1, axis=1)
        req_muon = ak.count(event_muon.pt, axis=1) == 1

        ## Jet cuts
        event_jet = events.Jet[
            jet_id(events, self._campaign)
            & (ak.all(events.Jet.metric_table(events.Muon) > 0.4, axis=2))
        ]
        req_jets = ak.num(event_jet.pt) >= 4

        ## store jet index for PFCands, create mask on the jet index
        jetindx = ak.mask(
            ak.local_index(events.Jet.pt),
            (
                jet_id(events, self._campaign)
                & (ak.all(events.Jet.metric_table(events.Muon) > 0.4, axis=2))
            )
            == 1,
        )
        jetindx = ak.pad_none(jetindx, 4)
        jetindx = jetindx[:, :4]

        ## other cuts
        if "Run3" not in self._campaign:
            MET = ak.zip(
                {
                    "pt": events.MET.pt,
                    "eta": ak.zeros_like(events.MET.pt),
                    "phi": events.MET.phi,
                    "mass": ak.zeros_like(events.MET.pt),
                },
                with_name="PtEtaPhiMLorentzVector",
            )
        else:
            MET = ak.zip(
                {
                    "pt": events.PuppiMET.pt,
                    "eta": ak.zeros_like(events.PuppiMET.pt),
                    "phi": events.PuppiMET.phi,
                    "mass": ak.zeros_like(events.PuppiMET.pt),
                },
                with_name="PtEtaPhiMLorentzVector",
            )
        req_MET = MET.pt > 50

        event_level = ak.fill_none(
            req_trig & req_jets & req_muon & req_MET & req_lumi, False
        )

        ####################
        # Selected objects #
        ####################
        smu = event_muon[event_level]
        sjets = event_jet[event_level]
        sel_jets = sjets
        sjets = sjets[:, :4]
        smet = MET[event_level]
        # Find the PFCands associate with selected jets. Search from jetindex->JetPFCands->PFCand
        if self._campaign != "Rereco17_94X":
            jetindexall = collections.defaultdict(dict)
            spfcands = collections.defaultdict(dict)
            for i in range(4):
                jetindexall[i] = jetindx[:, i]
                spfcands[i] = events[event_level].PFCands[
                    events[event_level]
                    .JetPFCands[
                        events[event_level].JetPFCands.jetIdx
                        == jetindexall[i][event_level]
                    ]
                    .pFCandsIdx
                ]

            # spfcands[1] = events[event_level].PFCands[events[event_level].JetPFCands[events[event_level].JetPFCands.jetIdx==jetindx0[event_level]].pFCandsIdx]

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
                            ak.firsts(event_muon),
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

            if self.isCorr and "BTV" in correction_config[self._campaign].keys():
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
                for i in range(4):
                    h.fill(
                        flatten(
                            ak.broadcast_arrays(genflavor[:, i], spfcands[i]["pt"])[0]
                        ),
                        flatten(spfcands[i][histname.replace("PFCands_", "")]),
                        weight=flatten(
                            ak.broadcast_arrays(
                                weights.weight()[event_level], spfcands[i]["pt"]
                            )[0]
                        ),
                    )
            elif "btagDeep" in histname:
                for i in range(4):
                    sel_jet = sjets[:, i]
                    if str(i) in histname:
                        h.fill(
                            flav=genflavor[:, i],
                            syst="noSF",
                            discr=np.where(
                                sel_jet[histname.replace(f"_{i}", "")] < 0,
                                -0.2,
                                sel_jet[histname.replace(f"_{i}", "")],
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
                            for syst in disc_list[histname.replace(f"_{i}", "")][
                                i
                            ].keys():
                                h.fill(
                                    flav=genflavor[:, i],
                                    syst=syst,
                                    discr=np.where(
                                        sel_jet[histname.replace(f"_{i}", "")] < 0,
                                        -0.2,
                                        sel_jet[histname.replace(f"_{i}", "")],
                                    ),
                                    weight=weights.weight()[event_level]
                                    * disc_list[histname.replace(f"_{i}", "")][i][syst],
                                )
            elif "mu_" in histname and histname.replace("mu_", "") in smu.fields:
                h.fill(
                    flatten(smu[histname.replace("mu_", "")]),
                    weight=weights.weight()[event_level],
                )
            elif "jet" in histname and "dr" not in histname and "njet" != histname:
                for i in range(4):
                    sel_jet = sjets[:, i]
                    if str(i) in histname:
                        h.fill(
                            flatten(genflavor[:, i]),
                            flatten(sel_jet[histname.replace(f"jet{i}_", "")]),
                            weight=weights.weight()[event_level],
                        )

        for i in range(4):
            output[f"dr_mujet{i}"].fill(
                flav=flatten(genflavor[:, i]),
                dr=flatten(smu.delta_r(sjets[:, i])),
                weight=weights.weight()[event_level],
            )
        output["njet"].fill(
            ak.count(sjets.pt, axis=1), weight=weights.weight()[event_level]
        )
        output["MET_pt"].fill(flatten(smet.pt), weight=weights.weight()[event_level])
        output["MET_phi"].fill(flatten(smet.phi), weight=weights.weight()[event_level])
        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator
