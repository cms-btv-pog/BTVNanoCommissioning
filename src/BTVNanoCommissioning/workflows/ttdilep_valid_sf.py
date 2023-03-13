import collections, numpy as np, awkward as ak

from coffea import processor
from coffea.analysis_tools import Weights


from BTVNanoCommissioning.utils.correction import (
    load_lumi,
    load_SF,
    eleSFs,
    muSFs,
    puwei,
    btagSFs,
    load_jmefactory,
)
from BTVNanoCommissioning.helpers.func import flatten
from BTVNanoCommissioning.helpers.update_branch import missing_branch, add_jec

from BTVNanoCommissioning.utils.histogrammer import histogrammer
from BTVNanoCommissioning.utils.selection import jet_id, mu_idiso, ele_cuttightid


class NanoProcessor(processor.ProcessorABC):
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
        _hist_event_dict = histogrammer("ttdilep_sf")
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
        triggers = [
            "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
            "Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
            "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
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
        events.Muon = events.Muon[
            (events.Muon.pt > 30) & mu_idiso(events, self._campaign)
        ]
        events.Muon = ak.pad_none(events.Muon, 1, axis=1)
        req_muon = ak.count(events.Muon.pt, axis=1) == 1

        ## Electron cuts
        # electron twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
        events.Electron = events.Electron[
            (events.Electron.pt > 30) & ele_cuttightid(events, self._campaign)
        ]
        events.Electron = ak.pad_none(events.Electron, 1, axis=1)
        req_ele = ak.count(events.Electron.pt, axis=1) == 1

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
            & (
                ak.all(
                    events.Jet.metric_table(events.Electron) > 0.4,
                    axis=2,
                    mask_identity=True,
                )
            )
        ]
        req_jets = ak.num(event_jet.pt) >= 2

        ## Other cuts
        req_opposite_charge = (
            events.Electron[:, 0:1].charge * events.Muon[:, 0:1].charge
        ) == -1
        req_opposite_charge = ak.fill_none(req_opposite_charge, False)
        req_opposite_charge = ak.flatten(req_opposite_charge)

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
                & (
                    ak.all(
                        events.Jet.metric_table(events.Electron) > 0.4,
                        axis=2,
                        mask_identity=True,
                    )
                )
            )
            == 1,
        )
        jetindx = ak.pad_none(jetindx, 2)
        jetindx = jetindx[:, :2]

        event_level = (
            req_trig & req_lumi & req_muon & req_ele & req_jets & req_opposite_charge
        )
        event_level = ak.fill_none(event_level, False)
        if len(events[event_level]) == 0:
            return {dataset: output}

        ####################
        # Selected objects #
        ####################
        smu = events.Muon[event_level]
        sel = events.Electron[event_level]
        sjets = event_jet[event_level]
        nseljet = ak.count(sjets.pt, axis=1)
        sjets = sjets[:, :2]
        # Find the PFCands associate with selected jets. Search from jetindex->JetPFCands->PFCand
        if "PFCands" in events.fields:
            jetindx0 = jetindx[:, 0]
            jetindx1 = jetindx[:, 1]
            spfcands = collections.defaultdict(dict)
            spfcands[0] = events[event_level].PFCands[
                events[event_level]
                .JetPFCands[
                    events[event_level].JetPFCands.jetIdx == jetindx0[event_level]
                ]
                .pFCandsIdx
            ]
            spfcands[1] = events[event_level].PFCands[
                events[event_level]
                .JetPFCands[
                    events[event_level].JetPFCands.jetIdx == jetindx1[event_level]
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
                weights.add("lep2sf", eleSFs(sel[:, 0], self.SF_map))

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
                for i in range(2):
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
                "PFCands" in events.fields
                and "PFCands" in histname
                and histname.split("_")[1] in events.PFCands.fields
            ):
                for i in range(2):
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
                for i in range(2):
                    sel_jet = sjets[:, i]
                    if (
                        str(i) in histname
                        and histname.replace(f"_{i}", "") in events.Jet.fields
                    ):
                        h.fill(
                            flav=flatten(genflavor[:, i]),
                            syst="noSF",
                            discr=flatten(
                                np.where(
                                    sel_jet[histname.replace(f"_{i}", "")] < 0,
                                    -0.2,
                                    sel_jet[histname.replace(f"_{i}", "")],
                                )
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
                                    flav=flatten(genflavor[:, i]),
                                    syst=syst,
                                    discr=flatten(
                                        np.where(
                                            sel_jet[histname.replace(f"_{i}", "")] < 0,
                                            -0.2,
                                            sel_jet[histname.replace(f"_{i}", "")],
                                        )
                                    ),
                                    weight=weights.weight()
                                    * disc_list[histname.replace(f"_{i}", "")][i][syst],
                                )
            elif "mu_" in histname and histname.replace("mu_", "") in smu.fields:
                h.fill(
                    flatten(smu[histname.replace("mu_", "")]),
                    weight=weights.weight(),
                )
            elif "ele_" in histname and histname.replace("ele_", "") in sel.fields:
                h.fill(
                    flatten(sel[histname.replace("ele_", "")]),
                    weight=weights.weight(),
                )
            elif "jet" in histname and "dr" not in histname and "njet" != histname:
                for i in range(2):
                    sel_jet = sjets[:, i]
                    if str(i) in histname:
                        h.fill(
                            flatten(genflavor[:, i]),
                            flatten(sel_jet[histname.replace(f"jet{i}_", "")]),
                            weight=weights.weight(),
                        )

        for i in range(2):
            output[f"dr_mujet{i}"].fill(
                flav=flatten(genflavor[:, i]),
                dr=flatten(smu.delta_r(sjets[:, i])),
                weight=weights.weight(),
            )
        output["njet"].fill(nseljet, weight=weights.weight())

        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator
