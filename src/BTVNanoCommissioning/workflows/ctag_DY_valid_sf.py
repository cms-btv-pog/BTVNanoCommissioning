import collections, awkward as ak, numpy as np
import os
import uproot
from coffea import processor
from coffea.analysis_tools import Weights

from BTVNanoCommissioning.utils.correction import (
    load_lumi,
    load_SF,
    muSFs,
    eleSFs,
    puwei,
    btagSFs,
    JME_shifts,
    Roccor_shifts,
)

from BTVNanoCommissioning.helpers.func import (
    flatten,
    update,
    dump_lumi,
)
from BTVNanoCommissioning.helpers.update_branch import missing_branch
from BTVNanoCommissioning.utils.histogrammer import histogrammer
from BTVNanoCommissioning.utils.array_writer import array_writer
from BTVNanoCommissioning.utils.selection import (
    jet_id,
    mu_idiso,
    ele_mvatightid,
    MET_filters,
)


class NanoProcessor(processor.ProcessorABC):
    def __init__(
        self,
        year="2022",
        campaign="Summer22Run3",
        name="",
        isSyst=False,
        isArray=False,
        noHist=False,
        chunksize=75000,
        selectionModifier="DYM",
    ):
        self._year = year
        self._campaign = campaign
        self.name = name
        self.isSyst = isSyst
        self.isArray = isArray
        self.noHist = noHist
        self.lumiMask = load_lumi(self._campaign)
        self.chunksize = chunksize
        self.selMod = selectionModifier
        ## Load corrections
        self.SF_map = load_SF(self._campaign)

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        isRealData = not hasattr(events, "genWeight")
        dataset = events.metadata["dataset"]
        events = missing_branch(events)
        shifts = []
        if "JME" in self.SF_map.keys():
            syst_JERC = self.isSyst
            if self.isSyst == "JERC_split":
                syst_JERC = "split"
            shifts = JME_shifts(
                shifts, self.SF_map, events, self._campaign, isRealData, syst_JERC
            )
        else:
            if int(self._year) < 2020:
                shifts = [
                    ({"Jet": events.Jet, "MET": events.MET, "Muon": events.Muon}, None)
                ]
            else:
                shifts = [
                    (
                        {
                            "Jet": events.Jet,
                            "MET": events.PuppiMET,
                            "Muon": events.Muon,
                        },
                        None,
                    )
                ]
        if "roccor" in self.SF_map.keys():
            shifts = Roccor_shifts(shifts, self.SF_map, events, isRealData, False)
        else:
            shifts[0][0]["Muon"] = events.Muon

        return processor.accumulate(
            self.process_shift(update(events, collections), name)
            for collections, name in shifts
        )

    def process_shift(self, events, shift_name):
        dataset = events.metadata["dataset"]
        isRealData = not hasattr(events, "genWeight")

        isMu = False
        isEle = False
        if "DYM" in self.selMod:
            triggers = ["Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8"]
            isMu = True
        elif "DYE" in self.selMod:
            triggers = ["Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"]
            isEle = True
        else:
            raise ValueError(self.selMod, "is not a valid selection modifier.")

        histname = {"DYM": "ctag_DY_sf", "DYE": "ectag_DY_sf"}
        _hist_event_dict = (
            {"": None} if self.noHist else histogrammer(events, histname[self.selMod])
        )
        output = {
            "sumw": processor.defaultdict_accumulator(float),
            **_hist_event_dict,
        }

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
        # only dump for nominal case
        if shift_name is None:
            output = dump_lumi(events[req_lumi], output)

        ## HLT
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

        req_metfilter = MET_filters(events, self._campaign)

        ## Muon cuts
        dilep_mu = events.Muon[(events.Muon.pt > 12) & mu_idiso(events, self._campaign)]
        ## Electron cuts
        dilep_ele = events.Electron[
            (events.Electron.pt > 15) & ele_mvatightid(events, self._campaign)
        ]
        if isMu:
            thisdilep = dilep_mu
            otherdilep = dilep_ele
        else:
            thisdilep = dilep_ele
            otherdilep = dilep_mu
        ## dilepton
        pos_dilep = thisdilep[thisdilep.charge > 0]
        neg_dilep = thisdilep[thisdilep.charge < 0]
        req_dilep = ak.fill_none(
            (
                (ak.num(pos_dilep.pt) >= 1)
                & (ak.num(neg_dilep.pt) >= 1)
                & (ak.num(thisdilep.charge) >= 2)
                & (ak.num(otherdilep.charge) == 0)
            ),
            False,
            axis=-1,
        )

        jet_sel = ak.fill_none(
            jet_id(events, self._campaign)
            & (
                ak.all(
                    events.Jet.metric_table(pos_dilep) > 0.4,
                    axis=2,
                    mask_identity=True,
                )
            )
            & (
                ak.all(
                    events.Jet.metric_table(neg_dilep) > 0.4,
                    axis=2,
                    mask_identity=True,
                )
            ),
            False,
            axis=-1,
        )

        pos_dilep = ak.pad_none(pos_dilep, 1, axis=1)
        neg_dilep = ak.pad_none(neg_dilep, 1, axis=1)

        dilep_mass = pos_dilep[:, 0] + neg_dilep[:, 0]
        req_dilepmass = (
            (dilep_mass.mass > 81) & (dilep_mass.mass < 101) & (dilep_mass.pt > 15)
        )

        ## Jet cuts
        event_jet = events.Jet[
            ak.fill_none(
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
                ),
                False,
                axis=-1,
            )
        ]
        req_jets = ak.count(event_jet.pt, axis=1) >= 1
        # event_jet = ak.pad_none(event_jet, 1, axis=1)

        ## store jet index for PFCands, create mask on the jet index
        jetindx = ak.mask(
            ak.local_index(events.Jet.pt),
            jet_sel == 1,
        )
        jetindx = ak.pad_none(jetindx, 1)
        jetindx = jetindx[:, 0]

        event_level = ak.fill_none(
            req_lumi & req_trig & req_dilep & req_dilepmass & req_jets & req_metfilter,
            False,
        )
        if len(events[event_level]) == 0:
            if self.isArray:
                array_writer(
                    self,
                    events[event_level],
                    events,
                    "nominal",
                    dataset,
                    isRealData,
                    empty=True,
                )
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
        sel_jet = sjets[:, 0]
        njet = ak.count(sjets.pt, axis=1)
        sel_mu = ak.concatenate([sposmu, snegmu])
        smu = ak.zip(
            {
                b: ak.Array(np.reshape(sel_mu[b].to_numpy(), (len(sposmu[b]), 2)))
                for b in sposmu.fields
            }
        )

        # Find the PFCands associate with selected jets. Search from jetindex->JetPFCands->PFCand
        if "PFCands" in events.fields:
            spfcands = events[event_level].PFCands[
                events[event_level]
                .JetPFCands[
                    events[event_level].JetPFCands.jetIdx == jetindx[event_level]
                ]
                .pFCandsIdx
            ]

        ####################
        # Weight & Geninfo #
        ####################
        weights = Weights(len(events[event_level]), storeIndividual=True)
        if not isRealData:
            weights.add("genweight", events[event_level].genWeight)
            par_flav = (sel_jet.partonFlavour == 0) & (sel_jet.hadronFlavour == 0)
            genflavor = sel_jet.hadronFlavour + 1 * par_flav
            if len(self.SF_map.keys()) > 0:
                syst_wei = True if self.isSyst != False else False
                if "PU" in self.SF_map.keys():
                    puwei(
                        events[event_level].Pileup.nTrueInt,
                        self.SF_map,
                        weights,
                        syst_wei,
                    )
                if isMu and "MUO" in self.SF_map.keys():
                    muSFs(smu, self.SF_map, weights, syst_wei, False)
                if isEle and "EGM" in self.SF_map.keys():
                    eleSFs(smu, self.SF_map, weights, syst_wei, False)
                if "BTV" in self.SF_map.keys():
                    btagSFs(sel_jet, self.SF_map, weights, "DeepJetC", syst_wei)
                    btagSFs(sel_jet, self.SF_map, weights, "DeepJetB", syst_wei)
                    btagSFs(sel_jet, self.SF_map, weights, "DeepCSVB", syst_wei)
                    btagSFs(sel_jet, self.SF_map, weights, "DeepCSVC", syst_wei)
        else:
            genflavor = ak.zeros_like(sel_jet.pt, dtype=int)

        # Systematics information
        if shift_name is None:
            systematics = ["nominal"] + list(weights.variations)
        else:
            systematics = [shift_name]
        exclude_btv = [
            "DeepCSVC",
            "DeepCSVB",
            "DeepJetB",
            "DeepJetC",
        ]  # exclude b-tag SFs for btag inputs

        ####################
        #  Fill histogram  #
        ####################
        for syst in systematics:
            if self.isSyst == False and syst != "nominal":
                break
            if self.noHist:
                break
            weight = (
                weights.weight()
                if syst == "nominal" or syst == shift_name
                else weights.weight(modifier=syst)
            )
            for histname, h in output.items():
                if (
                    "Deep" in histname
                    and "btag" not in histname
                    and histname in events.Jet.fields
                ):
                    h.fill(
                        syst,
                        genflavor,
                        sel_jet[histname],
                        weight=weights.partial_weight(exclude=exclude_btv),
                    )
                elif (
                    "PFCands" in events.fields
                    and "PFCands" in histname
                    and histname.split("_")[1] in events.PFCands.fields
                ):
                    h.fill(
                        syst,
                        flatten(ak.broadcast_arrays(genflavor, spfcands["pt"])[0]),
                        flatten(spfcands[histname.replace("PFCands_", "")]),
                        weight=flatten(
                            ak.broadcast_arrays(
                                weights.partial_weight(exclude=exclude_btv),
                                spfcands["pt"],
                            )[0]
                        ),
                    )
                elif (
                    "posl_" in histname
                    and histname.replace("posl_", "") in sposmu.fields
                ):
                    h.fill(
                        syst,
                        flatten(sposmu[histname.replace("posl_", "")]),
                        weight=weight,
                    )
                elif (
                    "negl_" in histname
                    and histname.replace("negl_", "") in snegmu.fields
                ):
                    h.fill(
                        syst,
                        flatten(snegmu[histname.replace("negl_", "")]),
                        weight=weight,
                    )

                elif "jet_" in histname:
                    h.fill(
                        syst,
                        genflavor,
                        sel_jet[histname.replace("jet_", "")],
                        weight=weight,
                    )
                elif (
                    "btag" in histname
                    and "0" in histname
                    and histname.replace("_0", "") in events.Jet.fields
                ):
                    h.fill(
                        syst="noSF",
                        flav=genflavor,
                        discr=np.where(
                            sel_jet[histname.replace("_0", "")] < 0,
                            -0.2,
                            sel_jet[histname.replace("_0", "")],
                        ),
                        weight=weights.partial_weight(exclude=exclude_btv),
                    )
                    if not isRealData and "btag" in self.SF_map.keys():
                        h.fill(
                            syst=syst,
                            flav=genflavor,
                            discr=np.where(
                                sel_jet[histname.replace("_0", "")] < 0,
                                -0.2,
                                sel_jet[histname.replace("_0", "")],
                            ),
                            weight=weight,
                        )
            output["njet"].fill(syst, njet, weight=weight)
            output["dr_mumu"].fill(syst, snegmu.delta_r(sposmu), weight=weight)
            output["z_pt"].fill(syst, flatten(sz.pt), weight=weight)
            output["z_eta"].fill(syst, flatten(sz.eta), weight=weight)
            output["z_phi"].fill(syst, flatten(sz.phi), weight=weight)
            output["z_mass"].fill(syst, flatten(sz.mass), weight=weight)
            output["npvs"].fill(
                syst,
                events[event_level].PV.npvs,
                weight=weight,
            )
            if not isRealData:
                output["pu"].fill(
                    syst,
                    events[event_level].Pileup.nTrueInt,
                    weight=weight,
                )
        #######################
        #  Create root files  #
        #######################
        if self.isArray:
            # Keep the structure of events and pruned the object size
            pruned_ev = events[event_level]
            pruned_ev["SelJet"] = sjets
            pruned_ev["LeadJet"] = sel_jet
            if isMu:
                pruned_ev["MuonPlus"] = sposmu
                pruned_ev["MuonMinus"] = snegmu
                kinOnly = ["Muon", "MuonPlus", "MuonMinus", "Jet"]
            else:
                pruned_ev["ElectronPlus"] = sposmu
                pruned_ev["ElectronMinus"] = snegmu
                kinOnly = ["Electron", "ElectronPlus", "ElectronMinus", "Jet"]
            pruned_ev["dilep"] = sposmu + snegmu
            pruned_ev["dilep", "pt"] = pruned_ev.dilep.pt
            pruned_ev["dilep", "eta"] = pruned_ev.dilep.eta
            pruned_ev["dilep", "phi"] = pruned_ev.dilep.phi
            pruned_ev["dilep", "mass"] = pruned_ev.dilep.mass
            pruned_ev["njet"] = njet
            if "PFCands" in events.fields:
                pruned_ev.PFCands = spfcands
            # Add custom variables
            if not isRealData:
                pruned_ev["weight"] = weights.weight()
                for ind_wei in weights.weightStatistics.keys():
                    pruned_ev[f"{ind_wei}_weight"] = weights.partial_weight(
                        include=[ind_wei]
                    )

            pruned_ev["dr_mu1jet"] = sposmu.delta_r(sel_jet)
            pruned_ev["dr_mu2jet"] = snegmu.delta_r(sel_jet)

            array_writer(
                self,
                pruned_ev,
                events,
                systematics[0],
                dataset,
                isRealData,
                kinOnly=kinOnly,
                remove=kinOnly,
            )

        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator
