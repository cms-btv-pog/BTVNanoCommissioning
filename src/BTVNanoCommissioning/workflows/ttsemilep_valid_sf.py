import collections, gc
import os
import uproot
import numpy as np, awkward as ak

from coffea import processor
from coffea.analysis_tools import Weights

from BTVNanoCommissioning.utils.correction import (
    load_lumi,
    load_SF,
    muSFs,
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
from BTVNanoCommissioning.utils.array_writer import array_writer
from BTVNanoCommissioning.utils.selection import jet_id, btag_mu_idiso, MET_filters
import hist


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
        selectionModifier="tt_semilep",
    ):
        self._year = year
        self._campaign = campaign
        self.name = name
        self.isSyst = isSyst
        self.isArray = isArray
        self.noHist = noHist
        self.lumiMask = load_lumi(self._campaign)
        self.chunksize = chunksize
        ### Added selection for ttbar semileptonic
        self.ttaddsel = selectionModifier
        ## Load corrections
        self.SF_map = load_SF(self._campaign)

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        isRealData = not hasattr(events, "genWeight")
        dataset = events.metadata["dataset"]

        attributes = dir(events.Jet)

        # Filter attributes that contain 'btagDeepFlav' in their name
        btag_attributes = [attr for attr in attributes if "btagDeepFlav" in attr]

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
        _hist_event_dict = (
            {"": None} if self.noHist else histogrammer(events, "ttsemilep_sf")
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

        ## Muon cuts
        # muon twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
        event_muon = events.Muon[
            (events.Muon.pt > 30) & btag_mu_idiso(events, self._campaign)
        ]
        # event_muon = ak.pad_none(events.Muon, 1, axis=1)
        req_muon = ak.count(event_muon.pt, axis=1) == 1

        ## Jet cuts
        event_jet = events.Jet[
            ak.fill_none(
                jet_id(events, self._campaign)
                & (
                    ak.all(
                        events.Jet.metric_table(events.Muon) > 0.4,
                        axis=2,
                        mask_identity=True,
                    )
                ),
                False,
                axis=-1,
            )
        ]
        req_jets = ak.num(event_jet.pt) >= 4

        if self.ttaddsel == "c_tt_semilep":
            print("ttbar semileptonic selection")
            # Sort jets by pt in descending order

            # Sort all jets by pt in descending order
            sorted_jets = event_jet[ak.argsort(event_jet.pt, axis=-1, ascending=False)]

            # Pad the sorted jets to ensure there are at least two jets in each event
            padded_jets = ak.pad_none(sorted_jets, 2, axis=1)

            # Determine the leading and subleading jets
            leading_jet = padded_jets[:, 0]
            subleading_jet = padded_jets[:, 1]

            c_selection = (
                event_jet.btagDeepFlavC > 0.6
            )  # & (event_jet.pt == leading_jet.pt)
            ### plot # c jets , DeepFlavC of the c jets
            c_selection = ak.fill_none(c_selection, False)
            req_c_jets = ak.num(event_jet[c_selection].pt) >= 1
            req_c_jets = ak.fill_none(req_c_jets, False)
            c_jets = event_jet[c_selection]
            single_c_jet = ak.firsts(c_jets)
            remaining_jets = event_jet[~c_selection]
            nob_selection = remaining_jets.btagDeepFlavB < 0.5

            nob_candidate_jets = remaining_jets[nob_selection]

            # Create combinations of single_c_jet and nob_candidate_jets
            combinations = ak.cartesian(
                {"c_jet": single_c_jet, "nob_jet": nob_candidate_jets}, axis=1
            )

            # Calculate the invariant mass of each pair using the .mass attribute
            masses = (combinations["c_jet"] + combinations["nob_jet"]).mass

            # Create a mask for pairs with mass between 55 and 110
            mass_mask = (masses > 50) & (masses < 110)

            # Apply the mask to get the valid combinations
            valid_combinations = combinations[mass_mask]
            # Create an event-level mask
            event_mask = ak.any(mass_mask, axis=1)

            req_w_mass = ak.fill_none(event_mask, False)

            nob_w_selection = (
                event_jet.btagDeepFlavB < 0.5
            )  # & (event_jet.pt == subleading_jet.pt)
            nob_w_selection = ak.fill_none(nob_w_selection, False)
            w_cand_b_jet = ak.firsts(event_jet[nob_w_selection])

            w_mass_lead_sublead = (single_c_jet + w_cand_b_jet).mass

            w_mask = (w_mass_lead_sublead > 50) & (w_mass_lead_sublead < 110)
            w_mask = ak.fill_none(w_mask, False)

            # req_w_jets = ak.fill_none(req_w_jets, False)

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

        MET = ak.zip(
            {
                "pt": events.MET.pt,
                "eta": ak.zeros_like(events.MET.pt),
                "phi": events.MET.phi,
                "mass": ak.zeros_like(events.MET.pt),
            },
            with_name="PtEtaPhiMLorentzVector",
        )

        req_MET = MET.pt > 50

        req_metfilter = MET_filters(events, self._campaign)
        event_level = ak.fill_none(
            req_trig & req_jets & req_muon & req_MET & req_lumi & req_metfilter, False
        )
        if self.ttaddsel == "c_tt_semilep":

            # Combine masks with the original conditions
            event_level = ak.fill_none(event_level & req_w_mass & req_c_jets, False)
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
        smu = event_muon[event_level]
        smu = smu[:, 0]
        sjets = event_jet[event_level]

        # Filter jets with btagDeepFlavC > 0.6
        ### With such a cut we reach c-purity of 70 %
        jet_c_sel_hist = sjets.btagDeepFlavC > 0.6
        if self.ttaddsel == "c_tt_semilep":
            scjets = sjets[jet_c_sel_hist]

        # Ensure the filtered jets are padded to maintain a consistent shape
        # scjets = ak.pad_none(scjets, 4, axis=1)

        # Calculate the number of selected jets
        if self.ttaddsel == "c_tt_semilep":
            ncseljet = ak.count(scjets.pt, axis=1)
        nseljet = ak.count(sjets.pt, axis=1)

        # Limit the number of jets to 4
        sjets = sjets[:, :4]

        if self.ttaddsel == "c_tt_semilep":
            scjets = scjets[:, 0]

        sw = sjets[:, 0] + sjets[:, 1]
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
            par_flav = (sjets.partonFlavour == 0) & (sjets.hadronFlavour == 0)
            if self.ttaddsel == "c_tt_semilep":
                par_flav_c = (scjets.partonFlavour == 0) & (scjets.hadronFlavour == 0)
            genflavor = ak.values_astype(sjets.hadronFlavour + 1 * par_flav, int)
            if self.ttaddsel == "c_tt_semilep":
                genflavor_c = ak.values_astype(
                    scjets.hadronFlavour + 1 * par_flav_c, int
                )
            if len(self.SF_map.keys()) > 0:
                syst_wei = True if self.isSyst != False else False
                if "PU" in self.SF_map.keys():
                    puwei(
                        events[event_level].Pileup.nTrueInt,
                        self.SF_map,
                        weights,
                        syst_wei,
                    )
                if "MUO" in self.SF_map.keys():
                    muSFs(smu, self.SF_map, weights, syst_wei, False)
                if "BTV" in self.SF_map.keys():
                    btagSFs(sjets, self.SF_map, weights, "DeepJetC", syst_wei)
                    btagSFs(sjets, self.SF_map, weights, "DeepJetB", syst_wei)
                    btagSFs(sjets, self.SF_map, weights, "DeepCSVB", syst_wei)
                    btagSFs(sjets, self.SF_map, weights, "DeepCSVC", syst_wei)
                    if self.ttaddsel == "c_tt_semilep":
                        btagSFs(scjets, self.SF_map, weights, "DeepJetC", syst_wei)
                        btagSFs(scjets, self.SF_map, weights, "DeepJetB", syst_wei)
                        btagSFs(scjets, self.SF_map, weights, "DeepCSVB", syst_wei)
                        btagSFs(scjets, self.SF_map, weights, "DeepCSVC", syst_wei)
        else:
            genflavor = ak.zeros_like(sjets.pt, dtype=int)
            if self.ttaddsel == "c_tt_semilep":
                genflavor_c = ak.zeros_like(scjets.pt, dtype=int)

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
            if not self.ttaddsel == "c_tt_semilep":
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
                                    weights.partial_weight(exclude=exclude_btv),
                                    sjets["pt"],
                                )[0]
                            ),
                        )
                    elif (
                        "PFCands" in histname
                        and "PFCands" in events.fields
                        and histname.split("_")[1] in events.PFCands.fields
                    ):
                        for i in range(4):
                            h.fill(
                                syst,
                                flatten(
                                    ak.broadcast_arrays(
                                        genflavor[:, i],
                                        spfcands[i]["pt"],
                                    )[0]
                                ),
                                flatten(spfcands[i][histname.replace("PFCands_", "")]),
                                weight=flatten(
                                    ak.broadcast_arrays(
                                        weights.partial_weight(exclude=exclude_btv),
                                        spfcands[i]["pt"],
                                    )[0]
                                ),
                            )
                    elif "btag" in histname:
                        for i in range(4):
                            sel_jet = sjets[:, i]
                            if (
                                str(i) in histname
                                and histname.replace(f"_{i}", "") in events.Jet.fields
                            ):
                                h.fill(
                                    syst="noSF",
                                    flav=genflavor[:, i],
                                    discr=sel_jet[histname.replace(f"_{i}", "")],
                                    weight=weight,
                                )
                                if (
                                    not isRealData
                                    and "btag" in self.SF_map.keys()
                                    and "_b" not in histname
                                    and "_bb" not in histname
                                    and "_lepb" not in histname
                                ):
                                    h.fill(
                                        syst=syst,
                                        flav=genflavor[:, i],
                                        discr=sel_jet[histname.replace(f"_{i}", "")],
                                        weight=weight,
                                    )
                    elif (
                        "mu_" in histname and histname.replace("mu_", "") in smu.fields
                    ):
                        h.fill(
                            syst,
                            flatten(smu[histname.replace("mu_", "")]),
                            weight=weight,
                        )
                    elif (
                        "jet" in histname
                        and "dr" not in histname
                        and "njet" != histname
                    ):
                        for i in range(4):
                            sel_jet = sjets[:, i]
                            if str(i) in histname:
                                h.fill(
                                    syst,
                                    flatten(genflavor[:, i]),
                                    flatten(sel_jet[histname.replace(f"jet{i}_", "")]),
                                    weight=weight,
                                )

                for i in range(4):
                    output[f"dr_mujet{i}"].fill(
                        syst,
                        flav=flatten(genflavor[:, i]),
                        dr=flatten(smu.delta_r(sjets[:, i])),
                        weight=weight,
                    )
                output["njet"].fill(syst, nseljet, weight=weight)
                output["MET_pt"].fill(syst, flatten(smet.pt), weight=weight)
                output["MET_phi"].fill(syst, flatten(smet.phi), weight=weight)
                output["npvs"].fill(
                    syst,
                    events[event_level].PV.npvs,
                    weight=weight,
                )
                output["w_mass"].fill(
                    syst,
                    flatten(sw.mass),
                    weight=weight,
                )
            else:
                for histname, h in output.items():
                    if (
                        "Deep" in histname
                        and "btag" not in histname
                        and histname in events.Jet.fields
                    ):
                        h.fill(
                            syst,
                            flatten(genflavor_c),
                            flatten(scjets[histname]),
                            weight=flatten(
                                ak.broadcast_arrays(
                                    weights.partial_weight(exclude=exclude_btv),
                                    scjets["pt"],
                                )[0]
                            ),
                        )
                    elif (
                        "PFCands" in histname
                        and "PFCands" in events.fields
                        and histname.split("_")[1] in events.PFCands.fields
                    ):
                        h.fill(
                            syst,
                            flatten(
                                ak.broadcast_arrays(
                                    genflavor_c,
                                    spfcands[0]["pt"],
                                )[0]
                            ),
                            flatten(spfcands[0][histname.replace("PFCands_", "")]),
                            weight=flatten(
                                ak.broadcast_arrays(
                                    weights.partial_weight(exclude=exclude_btv),
                                    spfcands[0]["pt"],
                                )[0]
                            ),
                        )
                    elif "btag" in histname:
                        sel_jet = scjets
                        if histname.replace(f"c_", "") in events.Jet.fields:
                            h.fill(
                                syst="noSF",
                                flav=genflavor_c,
                                discr=sel_jet[histname.replace(f"c_", "")],
                                weight=weight,
                            )
                            if (
                                not isRealData
                                and "btag" in self.SF_map.keys()
                                and "_b" not in histname
                                and "_bb" not in histname
                                and "_lepb" not in histname
                            ):
                                h.fill(
                                    syst=syst,
                                    flav=genflavor_c,
                                    discr=sel_jet[histname.replace(f"c_", "")],
                                    weight=weight,
                                )
                    elif (
                        "mu_" in histname and histname.replace("mu_", "") in smu.fields
                    ):
                        mu_array = flatten(smu[histname.replace("mu_", "")])
                        h.fill(
                            syst,
                            flatten(smu[histname.replace("mu_", "")]),
                            weight=weight,
                        )
                    elif (
                        "jet" in histname
                        and "dr" not in histname
                        and "njet" != histname
                    ):
                        if "c" in histname:
                            jet_array = flatten(scjets[histname.replace(f"cjet_", "")])

                            h.fill(
                                syst,
                                flatten(genflavor_c),
                                flatten(scjets[histname.replace(f"cjet_", "")]),
                                weight=weight,
                            )

                output[f"dr_cjet"].fill(
                    syst,
                    flav=flatten(genflavor_c),
                    dr=flatten(smu.delta_r(scjets)),
                    weight=weight,
                )
                output["njet"].fill(syst, ncseljet, weight=weight)
                output["MET_pt"].fill(syst, flatten(smet.pt), weight=weight)
                output["MET_phi"].fill(syst, flatten(smet.phi), weight=weight)
                output["npvs"].fill(
                    syst,
                    events[event_level].PV.npvs,
                    weight=weight,
                )
                output["w_mass"].fill(
                    syst,
                    flatten(sw.mass),
                    weight=weight,
                )
        #######################
        #  Create root files  #
        #######################
        if self.isArray:
            # Keep the structure of events and pruned the object size
            pruned_ev = events[event_level]
            pruned_ev["SelJet"] = sjets
            pruned_ev["Muon"] = smu
            if "PFCands" in events.fields:
                pruned_ev.PFCands = spfcands
            # Add custom variables
            if not isRealData:
                pruned_ev["weight"] = weights.weight()
                for ind_wei in weights.weightStatistics.keys():
                    pruned_ev[f"{ind_wei}_weight"] = weights.partial_weight(
                        include=[ind_wei]
                    )

            for i in range(4):
                pruned_ev[f"dr_mujet{i}"] = smu.delta_r(sjets[:, i])
            array_writer(self, pruned_ev, events, systematics[0], dataset, isRealData)
        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator
