import collections, gc
import numpy as np, awkward as ak

from coffea import processor
from coffea.analysis_tools import Weights

from BTVNanoCommissioning.utils.correction import (
    load_lumi,
    load_SF,
    muSFs,
    puwei,
    btagSFs,
    load_jmefactory,
)
from BTVNanoCommissioning.helpers.func import flatten
from BTVNanoCommissioning.helpers.update_branch import missing_branch, add_jec


from BTVNanoCommissioning.utils.histogrammer import histogrammer
from BTVNanoCommissioning.utils.selection import jet_id, btag_mu_idiso


class NanoProcessor(processor.ProcessorABC):
    # Define histograms
    def __init__(
        self,
        year="2017",
        campaign="Rereco17_94X",
        isCorr=True,
        isJERC=False,
        isSyst=False,
    ):
        self._year = year
        self._campaign = campaign
        self.isCorr = isCorr
        self.isJERC = isJERC
        self.isSyst = isSyst
        self.lumiMask = load_lumi(self._campaign)
        ## Load corrections
        if isCorr:
            self.SF_map = load_SF(self._campaign)
        if isJERC:
            self._jet_factory = load_jmefactory(self._campaign)
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

        if self.isJERC:
            events = add_jec(events, self._campaign, self._jet_factory)
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
            & (
                ak.all(
                    events.Jet.metric_table(events.Muon) > 0.4,
                    axis=2,
                    mask_identity=True,
                )
            )
        ]
        req_jets = ak.num(event_jet.pt) >= 4

        ## store jet index for PFCands, create mask on the jet index
        jetindx = ak.mask(
            ak.local_index(events.Jet.pt),
            (
                jet_id(events, self._campaign)
                & (
                    ak.all(
                        events.Jet.metric_table(events.Muon) > 0.4,
                        axis=2,
                        mask_identity=True,
                    )
                )
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
        if len(events[event_level]) == 0:
            return {dataset: output}
        ####################
        # Selected objects #
        ####################
        smu = event_muon[event_level]
        sjets = event_jet[event_level]
        nseljet = ak.count(sjets.pt, axis=1)
        sjets = sjets[:, :4]
        smet = MET[event_level]
        # Find the PFCands associate with selected jets. Search from jetindex->JetPFCands->PFCand
        if "PFCands" in events.fields:
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

        ####################
        # Weight & Geninfo #
        ####################
        weights = Weights(len(events[event_level]), storeIndividual=True)
        if not isRealData:
            weights.add("genweight", events[event_level].genWeight)
        if not isRealData and self.isCorr:
            if "PU" in self.SF_map.keys():
                weights.add(
                    "puweight", puwei(self.SF_map, events[event_level].Pileup.nTrueInt)
                )
            if "MUO" in self.SF_map.keys() or "EGM" in self.SF_map.keys():
                weights.add("lep1sf", muSFs(smu[:, 0], self.SF_map, True))

        if isRealData:
            genflavor = ak.zeros_like(sjets.pt)
        else:
            par_flav = (sjets.partonFlavour == 0) & (sjets.hadronFlavour == 0)
            genflavor = sjets.hadronFlavour + 1 * par_flav
            jetsfs_c = collections.defaultdict(dict)
            jetsfs_b = collections.defaultdict(dict)
            csvsfs_c = collections.defaultdict(dict)
            csvsfs_b = collections.defaultdict(dict)

            if self.isCorr and (
                "btag" in self.SF_map.keys() or "ctag" in self.SF_map.keys()
            ):
                for i in range(4):
                    jetsfs_c[i]["SF"] = btagSFs(sjets[:, i], self.SF_map, "DeepJetC")
                    jetsfs_b[i]["SF"] = btagSFs(sjets[:, i], self.SF_map, "DeepJetB")
                    csvsfs_c[i]["SF"] = btagSFs(sjets[:, i], self.SF_map, "DeepCSVC")
                    csvsfs_b[i]["SF"] = btagSFs(sjets[:, i], self.SF_map, "DeepCSVB")
                    if self.isSyst:
                        for syst in [
                            "hf",
                            "lf",
                            "cferr1",
                            "cferr2",
                            "hfstat1",
                            "hfstat2",
                            "lfstats1",
                            "lfstats2",
                        ]:
                            jetsfs_c[i][f"SF_{syst}_up"] = btagSFs(
                                sjets[:, i], self.SF_map, "DeepJetC", f"up_{syst}"
                            )
                            jetsfs_c[i][f"SF_{syst}_dn"] = btagSFs(
                                sjets[:, i], self.SF_map, "DeepJetC", f"down_{syst}"
                            )
                            csvsfs_c[i][f"SF_{syst}_up"] = btagSFs(
                                sjets[:, i], self.SF_map, "DeepCSVC", f"up_{syst}"
                            )
                            csvsfs_c[i][f"SF_{syst}_dn"] = btagSFs(
                                sjets[:, i], self.SF_map, "DeepCSVC", f"down_{syst}"
                            )
                        csvsfs_b[i][f"SF_{syst}_up"] = btagSFs(
                            sjets[:, i], self.SF_map, "DeepCSVB", f"up"
                        )
                        csvsfs_b[i][f"SF_{syst}_dn"] = btagSFs(
                            sjets[:, i], self.SF_map, "DeepCSVB", f"down"
                        )
                        jetsfs_b[i][f"SF_{syst}_up"] = btagSFs(
                            sjets[:, i], self.SF_map, "DeepJetB", f"up"
                        )
                        jetsfs_b[i][f"SF_{syst}_dn"] = btagSFs(
                            sjets[:, i], self.SF_map, "DeepJetB", f"down"
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
            if (
                "Deep" in histname
                and "btag" not in histname
                and histname in events.Jet.fields
            ):
                h.fill(
                    flatten(genflavor),
                    flatten(sjets[histname]),
                    weight=flatten(
                        ak.broadcast_arrays(weights.weight(), sjets["pt"])[0]
                    ),
                )
            elif (
                "PFCands" in histname
                and "PFCands" in events.fields
                and histname.split("_")[1] in events.PFCands.fields
            ):
                for i in range(4):
                    h.fill(
                        flatten(
                            ak.broadcast_arrays(genflavor[:, i], spfcands[i]["pt"])[0]
                        ),
                        flatten(spfcands[i][histname.replace("PFCands_", "")]),
                        weight=flatten(
                            ak.broadcast_arrays(weights.weight(), spfcands[i]["pt"])[0]
                        ),
                    )
            elif "btagDeep" in histname:
                for i in range(4):
                    sel_jet = sjets[:, i]
                    if (
                        str(i) in histname
                        and histname.replace(f"_{i}", "") in events.Jet.fields
                    ):
                        h.fill(
                            flav=genflavor[:, i],
                            syst="noSF",
                            discr=np.where(
                                sel_jet[histname.replace(f"_{i}", "")] < 0,
                                -0.2,
                                sel_jet[histname.replace(f"_{i}", "")],
                            ),
                            weight=weights.weight(),
                        )
                        if (
                            not isRealData
                            and self.isCorr
                            and "btag" in self.SF_map.keys()
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
                                    weight=weights.weight()
                                    * disc_list[histname.replace(f"_{i}", "")][i][syst],
                                )
            elif "mu_" in histname and histname.replace("mu_", "") in smu.fields:
                h.fill(
                    flatten(smu[histname.replace("mu_", "")]),
                    weight=weights.weight(),
                )
            elif "jet" in histname and "dr" not in histname and "njet" != histname:
                for i in range(4):
                    sel_jet = sjets[:, i]
                    if str(i) in histname:
                        h.fill(
                            flatten(genflavor[:, i]),
                            flatten(sel_jet[histname.replace(f"jet{i}_", "")]),
                            weight=weights.weight(),
                        )

        for i in range(4):
            output[f"dr_mujet{i}"].fill(
                flav=flatten(genflavor[:, i]),
                dr=flatten(smu.delta_r(sjets[:, i])),
                weight=weights.weight(),
            )
        output["njet"].fill(nseljet, weight=weights.weight())
        output["MET_pt"].fill(flatten(smet.pt), weight=weights.weight())
        output["MET_phi"].fill(flatten(smet.phi), weight=weights.weight())
        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator
