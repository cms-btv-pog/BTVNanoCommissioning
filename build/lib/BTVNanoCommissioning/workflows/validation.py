import collections, awkward as ak, numpy as np
import os
import uproot


from coffea import processor
from coffea.analysis_tools import Weights

from BTVNanoCommissioning.utils.correction import (
    load_lumi,
    load_SF,
    weight_manager,
    common_shifts,
)

from BTVNanoCommissioning.helpers.func import update, dump_lumi, PFCand_link, flatten
from BTVNanoCommissioning.helpers.update_branch import missing_branch
from BTVNanoCommissioning.utils.histogramming.histogrammer import (
    histogrammer,
    histo_writter,
)
from BTVNanoCommissioning.utils.array_writer import array_writer
from BTVNanoCommissioning.utils.selection import (
    HLT_helper,
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
        self.SF_map = load_SF(self._year, self._campaign)

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        events = missing_branch(events)
        vetoed_events, shifts = common_shifts(self, events)

        return processor.accumulate(
            self.process_shift(update(vetoed_events, collections), name)
            for collections, name in shifts
        )

    def process_shift(self, events, shift_name):
        dataset = events.metadata["dataset"]
        isRealData = not hasattr(events, "genWeight")
        output = {}
        if not self.noHist:
            output = histogrammer(
                events.Jet.fields,
                obj_list=["jet0", "jet1"],
                hist_collections=["common", "fourvec", "validation"],
            )

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
        req_trig = HLT_helper(events, triggers)

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
        # Keep the structure of events and pruned the object size
        pruned_ev = events[event_level]
        pruned_ev["SelJet"] = event_jet[event_level][:, :2]

        ####################
        #     Output       #
        ####################
        # Configure SFs
        weights = weight_manager(pruned_ev, self.SF_map, self.isSyst)
        # Configure systematics
        if shift_name is None:
            systematics = ["nominal"] + list(weights.variations)
        else:
            systematics = [shift_name]

        # Configure histograms
        if not self.noHist:
            exclude_btv = [
                "DeepCSVC",
                "DeepCSVB",
                "DeepJetB",
                "DeepJetC",
            ]  # exclude b-tag SFs for btag inputs
        if isRealData:
            genflavor = ak.zeros_like(pruned_ev.SelJet.pt, axis=-1)
        else:
            genflavor = ak.values_astype(
                pruned_ev.SelJet.hadronFlavour
                + 1 * (pruned_ev.SelJet.partonFlavour == 0)
                & (pruned_ev.SelJet.hadronFlavour == 0),
                int,
            )
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
                        flatten(pruned_ev.SelJet[histname]),
                        weight=flatten(
                            ak.broadcast_arrays(
                                weights.partial_weight(exclude=exclude_btv),
                                pruned_ev.SelJet["pt"],
                            )[0]
                        ),
                    )
                elif "WP" in histname:
                    jet = pruned_ev.SelJet[:, 0]

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
                        if (
                            histname.replace(f"jet{i}_", "")
                            not in pruned_ev.SelJet.fields
                        ):
                            continue
                        jet = pruned_ev.SelJet[:, i]
                        h.fill(
                            syst,
                            flatten(genflavor[:, i]),
                            flatten(jet[histname.replace(f"jet{i}_", "")]),
                            weight=weight,
                        )
                elif "btag" in histname:
                    for i in range(2):
                        if histname.replace(f"_{i}", "") not in pruned_ev.SelJet.fields:
                            continue
                        jet = pruned_ev.SelJet[:, i]
                        h.fill(
                            "noSF",
                            flav=flatten(genflavor[:, i]),
                            discr=jet[histname.replace(f"_{i}", "")],
                            weight=weight,
                        )

            # Add custom variables

        # Output arrays
        if self.isArray:

            array_writer(
                self, pruned_ev, events, weights, systematics, dataset, isRealData
            )
        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator
