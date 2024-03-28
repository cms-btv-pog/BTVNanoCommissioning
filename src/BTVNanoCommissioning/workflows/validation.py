import collections, awkward as ak, numpy as np
import os
import uproot


from coffea import processor
from coffea.analysis_tools import Weights

from BTVNanoCommissioning.helpers.func import flatten, update
from BTVNanoCommissioning.helpers.definitions import definitions
from BTVNanoCommissioning.helpers.update_branch import missing_branch
from BTVNanoCommissioning.utils.correction import (
    load_lumi,
    load_SF,
    puwei,
    btagSFs,
    JME_shifts,
    Roccor_shifts,
)

from BTVNanoCommissioning.helpers.func import (
    flatten,
    update,
    uproot_writeable,
    dump_lumi,
)
from BTVNanoCommissioning.helpers.update_branch import missing_branch

from BTVNanoCommissioning.utils.histogrammer import histogrammer
from BTVNanoCommissioning.utils.selection import (
    jet_id,
    mu_idiso,
    ele_mvatightid,
    btag_wp,
    btag_wp_dict,
)


class NanoProcessor(processor.ProcessorABC):
    # Define histograms
    def __init__(
        self,
        year="2022",
        campaign="Summer22Run3",
        name="",
        isSyst=False,
        isArray=False,
        noHist=False,
        chunksize=75000,
    ):
        self._year = year
        self._campaign = campaign
        self.name = name
        self.isSyst = False
        self.isArray = isArray
        self.noHist = noHist
        self.lumiMask = load_lumi(self._campaign)
        self.chunksize = chunksize
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
            if int(self._year) > 2020:
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
        _hist_event_dict = (
            {"": None} if self.noHist else histogrammer(events, "validation")
        )

        output = {
            "sumw": processor.defaultdict_accumulator(float),
            **_hist_event_dict,
        }

        if shift_name is None:
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

        ## Jet cuts
        event_jet = events.Jet[
            jet_id(events, self._campaign)
            & (ak.all(events.Jet.metric_table(events.Muon) > 0.4, axis=2))
            & (ak.all(events.Jet.metric_table(events.Electron) > 0.4, axis=2))
        ]
        req_jets = ak.num(event_jet.pt) >= 2

        ## store jet index for PFCands, create mask on the jet index
        jetindx = ak.mask(
            ak.local_index(events.Jet.pt),
            (
                jet_id(events, self._campaign)
                & (ak.all(events.Jet.metric_table(events.Muon) > 0.4, axis=2))
                & (ak.all(events.Jet.metric_table(events.Electron) > 0.4, axis=2))
            )
            == 1,
        )
        jetindx = ak.pad_none(jetindx, 2)
        jetindx = jetindx[:, :2]

        event_level = ak.fill_none(req_jets & req_lumi, False)
        if len(events[event_level]) == 0:
            return {dataset: output}
        ####################
        # Selected objects #
        ####################
        sjets = event_jet[event_level]
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
            par_flav = (sjets.partonFlavour == 0) & (sjets.hadronFlavour == 0)
            genflavor = ak.values_astype(sjets.hadronFlavour + 1 * par_flav, int)
            genweiev = ak.flatten(
                ak.broadcast_arrays(events[event_level].genWeight, sjets["pt"])[0]
            )
            if len(self.SF_map.keys()) > 0:
                syst_wei = True if self.isSyst == True else False
                if "PU" in self.SF_map.keys():
                    puwei(
                        events[event_level].Pileup.nTrueInt,
                        self.SF_map,
                        weights,
                        syst_wei,
                    )
                if "BTV" in self.SF_map.keys():
                    btagSFs(sjets, self.SF_map, weights, "DeepJetC", syst_wei)
                    btagSFs(sjets, self.SF_map, weights, "DeepJetB", syst_wei)
                    btagSFs(sjets, self.SF_map, weights, "DeepCSVB", syst_wei)
                    btagSFs(sjets, self.SF_map, weights, "DeepCSVC", syst_wei)
        else:
            genflavor = ak.zeros_like(sjets.pt, dtype=int)

        # Systematics information
        if shift_name is None:
            systematics = ["nominal"] + list(weights.variations)
        else:
            systematics = [shift_name]
        exclude_btv = [
            "DeepCSVC",
            "DeepCSVB",
            "DeepJetB",
            "DeepJetB",
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
                        flatten(genflavor),
                        flatten(sjets[histname]),
                        weight=flatten(
                            ak.broadcast_arrays(
                                weights.partial_weight(exclude=exclude_btv), sjets["pt"]
                            )[0]
                        ),
                    )
                elif "WP" in histname:
                    jet = sjets[:, 0]

                    for tagger in btag_wp_dict[self._campaign].keys():
                        if "bjet" in histname:
                            for wp in btag_wp_dict[self._campaign][tagger]["b"].keys():
                                wp_weight = weight[
                                    btag_wp(
                                        jet,
                                        self._campaign,
                                        tagger,
                                        "b",
                                        wp,
                                    )
                                    & (jet.hadronFlavour == 5)
                                ]
                                wp_jet = jet[
                                    btag_wp(
                                        jet,
                                        self._campaign,
                                        tagger,
                                        "b",
                                        wp,
                                    )
                                    & (jet.hadronFlavour == 5)
                                ]

                                if "discr" in histname:
                                    h.fill(
                                        wp,
                                        tagger,
                                        wp_jet[f"btag{tagger}B"],
                                        weight=wp_weight,
                                    )
                                else:
                                    h.fill(
                                        wp,
                                        tagger,
                                        wp_jet[histname.replace("bjet_WP_", "")],
                                        weight=wp_weight,
                                    )
                        elif "cjet" in histname:
                            for wp in btag_wp_dict[self._campaign][tagger]["c"].keys():
                                wp_weight = weight[
                                    btag_wp(
                                        jet,
                                        self._campaign,
                                        tagger,
                                        "c",
                                        wp,
                                    )
                                    & (jet.hadronFlavour == 4)
                                ]
                                wp_jet = jet[
                                    btag_wp(
                                        jet,
                                        self._campaign,
                                        tagger,
                                        "c",
                                        wp,
                                    )
                                    & (jet.hadronFlavour == 4)
                                ]

                                if "discr" in histname:
                                    h.fill(
                                        wp,
                                        tagger,
                                        wp_jet[f"btag{tagger}CvL"],
                                        wp_jet[f"btag{tagger}CvB"],
                                        weight=wp_weight,
                                    )
                                else:
                                    h.fill(
                                        wp,
                                        tagger,
                                        wp_jet[histname.replace("cjet_WP_", "")],
                                        weight=wp_weight,
                                    )
                elif "jet" in histname:
                    for i in range(2):
                        if histname.replace(f"jet{i}_", "") not in sjets.fields:
                            continue
                        jet = sjets[:, i]
                        h.fill(
                            syst,
                            flatten(genflavor[:, i]),
                            flatten(jet[histname.replace(f"jet{i}_", "")]),
                            weight=weight,
                        )
                elif "btag" in histname:
                    for i in range(2):
                        if histname.replace(f"_{i}", "") not in sjets.fields:
                            continue
                        # print(histname.replace(f"_{i}", ""))
                        jet = sjets[:, i]
                        h.fill(
                            "noSF",
                            flav=flatten(genflavor[:, i]),
                            discr=jet[histname.replace(f"_{i}", "")],
                            weight=weight,
                        )
        #######################
        #  Create root files  #
        #######################
        if self.isArray:
            # Keep the structure of events and pruned the object size
            pruned_ev = events[event_level]
            pruned_ev.Jet = sjets
            # Add custom variables
            if not isRealData:
                pruned_ev["weight"] = weights.weight()
                for ind_wei in weights.weightStatistics.keys():
                    pruned_ev[f"{ind_wei}_weight"] = weights.partial_weight(
                        include=[ind_wei]
                    )

            # Create a list of variables want to store. For objects from the PFNano file, specify as {object}_{variable}, wildcard option only accepted at the end of the string
            out_branch = np.setdiff1d(
                np.array(pruned_ev.fields), np.array(events.fields)
            )
            for kin in ["pt", "eta", "phi", "mass"]:
                for obj in ["Jet"]:
                    out_branch = np.append(out_branch, [f"{obj}_{kin}"])
            out_branch = np.append(
                out_branch, ["Jet_btagDeep*", "Jet_DeepJet*", "PFCands_*"]
            )
            # write to root files
            os.system(f"mkdir -p {self.name}/{dataset}")
            with uproot.recreate(
                f"{self.name}/{dataset}/f{events.metadata['filename'].split('_')[-1].replace('.root','')}_{systematics[0]}_{int(events.metadata['entrystop']/self.chunksize)}.root"
            ) as fout:
                fout["Events"] = uproot_writeable(pruned_ev, include=out_branch)
        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator
