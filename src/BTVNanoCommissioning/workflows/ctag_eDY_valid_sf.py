import collections, numpy as np, awkward as ak

from coffea import processor
from coffea.analysis_tools import Weights

from BTVNanoCommissioning.utils.correction import (
    load_lumi,
    load_SF,
    eleSFs,
    puwei,
    btagSFs,
    load_jmefactory,
)
from BTVNanoCommissioning.helpers.func import flatten
from BTVNanoCommissioning.helpers.update_branch import missing_branch, add_jec


from BTVNanoCommissioning.utils.histogrammer import histogrammer
from BTVNanoCommissioning.utils.selection import jet_id, mu_idiso, ele_mvatightid


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
        _hist_event_dict = histogrammer("ectag_DY_sf")
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
        triggers = ["Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"]
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

        dilep_mass = pos_dilep[:, 0] + neg_dilep[:, 0]
        req_dilepmass = (
            (dilep_mass.mass > 81)
            & (dilep_mass.mass < 101)
            & (dilep_mass.pt > 15)
            & ((pos_dilep[:, 0].pt > 27) | (neg_dilep[:, 0].pt > 27))
        )

        ## Jet cuts
        event_jet = events.Jet[
            jet_id(events, self._campaign)
            & (
                ak.all(
                    events.Jet.metric_table(pos_dilep[:, 0]) > 0.4,
                    axis=2,
                    mask_identity=True,
                )
            )
            & (
                ak.all(
                    events.Jet.metric_table(neg_dilep[:, 0]) > 0.4,
                    axis=2,
                    mask_identity=True,
                )
            )
        ]
        req_jets = ak.num(event_jet.pt) >= 1

        ## store jet index for PFCands, create mask on the jet index
        jetindx = ak.mask(
            ak.local_index(events.Jet.pt),
            (
                jet_id(events, self._campaign)
                & (
                    ak.all(
                        events.Jet.metric_table(pos_dilep[:, 0]) > 0.4,
                        axis=2,
                        mask_identity=True,
                    )
                )
                & (
                    ak.all(
                        events.Jet.metric_table(neg_dilep[:, 0]) > 0.4,
                        axis=2,
                        mask_identity=True,
                    )
                )
            )
            == 1,
        )
        jetindx = ak.pad_none(jetindx, 1)
        jetindx = jetindx[:, 0]

        event_level = ak.fill_none(
            req_lumi & req_trig & req_dilep & req_dilepmass & req_jets, False
        )
        if len(events[event_level]) == 0:
            return {dataset: output}

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
        if "PFCands" in events.fields:
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
        weights = Weights(len(events[event_level]), storeIndividual=True)
        if not isRealData:
            weights.add("genweight", events[event_level].genWeight)
        if not isRealData and self.isCorr:
            if "PU" in self.SF_map.keys():
                weights.add(
                    "puweight", puwei(self.SF_map, events[event_level].Pileup.nTrueInt)
                )
            if "MUO" in self.SF_map.keys() or "EGM" in self.SF_map.keys():
                weights.add("lep1sf", eleSFs(sposmu, self.SF_map, True))
                weights.add("lep2sf", eleSFs(snegmu, self.SF_map, True))

        if isRealData:
            genflavor = ak.zeros_like(sjets.pt)
        else:
            par_flav = (sjets.partonFlavour == 0) & (sjets.hadronFlavour == 0)
            genflavor = sjets.hadronFlavour + 1 * par_flav
            if self.isCorr and (
                "btag" in self.SF_map.keys() or "ctag" in self.SF_map.keys()
            ):
                jetsfs_c = collections.defaultdict(dict)
                jetsfs_b = collections.defaultdict(dict)
                csvsfs_c = collections.defaultdict(dict)
                csvsfs_b = collections.defaultdict(dict)

                if self.isCorr and (
                    "btag" in self.SF_map.keys() or "ctag" in self.SF_map.keys()
                ):
                    jetsfs_c[0]["SF"] = btagSFs(sjets[:, 0], self.SF_map, "DeepJetC")
                    jetsfs_b[0]["SF"] = btagSFs(sjets[:, 0], self.SF_map, "DeepJetB")
                    csvsfs_c[0]["SF"] = btagSFs(sjets[:, 0], self.SF_map, "DeepCSVC")
                    csvsfs_b[0]["SF"] = btagSFs(sjets[:, 0], self.SF_map, "DeepCSVB")
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
                            jetsfs_c[0][f"SF_{syst}_up"] = btagSFs(
                                sjets[:, 0], self.SF_map, "DeepJetC", f"up_{syst}"
                            )
                            jetsfs_c[0][f"SF_{syst}_dn"] = btagSFs(
                                sjets[:, 0], self.SF_map, "DeepJetC", f"down_{syst}"
                            )
                            csvsfs_c[0][f"SF_{syst}_up"] = btagSFs(
                                sjets[:, 0], self.SF_map, "DeepCSVC", f"up_{syst}"
                            )
                            csvsfs_c[0][f"SF_{syst}_dn"] = btagSFs(
                                sjets[:, 0], self.SF_map, "DeepCSVC", f"down_{syst}"
                            )
                        csvsfs_b[0][f"SF_{syst}_up"] = btagSFs(
                            sjets[:, 0], self.SF_map, "DeepCSVB", f"up"
                        )
                        csvsfs_b[0][f"SF_{syst}_dn"] = btagSFs(
                            sjets[:, 0], self.SF_map, "DeepCSVB", f"down"
                        )
                        jetsfs_b[0][f"SF_{syst}_up"] = btagSFs(
                            sjets[:, 0], self.SF_map, "DeepJetB", f"up"
                        )
                        jetsfs_b[0][f"SF_{syst}_dn"] = btagSFs(
                            sjets[:, 0], self.SF_map, "DeepJetB", f"down"
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
                and histname.split("_")[0] in events.PFCands.fields
            ):
                h.fill(
                    flatten(ak.broadcast_arrays(genflavor[:, 0], spfcands["pt"])[0]),
                    flatten(spfcands[histname.replace("PFCands_", "")]),
                    weight=flatten(
                        ak.broadcast_arrays(weights.weight(), spfcands["pt"])[0]
                    ),
                )
            elif "posl_" in histname and histname.replace("posl_", "") in sposmu.fields:
                h.fill(
                    flatten(sposmu[histname.replace("posl_", "")]),
                    weight=weights.weight(),
                )
            elif "negl_" in histname and histname.replace("negl_", "") in snegmu.fields:
                h.fill(
                    flatten(snegmu[histname.replace("negl_", "")]),
                    weight=weights.weight(),
                )
            elif "jet_" in histname:
                h.fill(
                    genflavor[:, 0],
                    sel_jet[histname.replace("jet_", "")],
                    weight=weights.weight(),
                )
            elif (
                "btagDeep" in histname
                and "0" in histname
                and histname.replace("_0", "") in events.Jet.fields
            ):
                h.fill(
                    flav=genflavor[:, 0],
                    syst="noSF",
                    discr=np.where(
                        sel_jet[histname.replace("_0", "")] < 0,
                        -0.2,
                        sel_jet[histname.replace("_0", "")],
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
                    for syst in disc_list[histname.replace("_0", "")][0].keys():
                        h.fill(
                            flav=genflavor[:, 0],
                            syst=syst,
                            discr=np.where(
                                sel_jet[histname.replace("_0", "")] < 0,
                                -0.2,
                                sel_jet[histname.replace("_0", "")],
                            ),
                            weight=weights.weight()
                            * disc_list[histname.replace("_0", "")][0][syst],
                        )
        output["njet"].fill(njet, weight=weights.weight())
        output["dr_mumu"].fill(dr=snegmu.delta_r(sposmu), weight=weights.weight())
        output["z_pt"].fill(flatten(sz.pt), weight=weights.weight())
        output["z_eta"].fill(flatten(sz.eta), weight=weights.weight())
        output["z_phi"].fill(flatten(sz.phi), weight=weights.weight())
        output["z_mass"].fill(flatten(sz.mass), weight=weights.weight())
        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator
