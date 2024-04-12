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
from BTVNanoCommissioning.utils.selection import jet_id, btag_mu_idiso


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
        self.isSyst = isSyst
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

        event_level = ak.fill_none(
            req_trig & req_jets & req_muon & req_MET & req_lumi, False
        )
        if len(events[event_level]) == 0:
            return {dataset: output}
        ####################
        # Selected objects #
        ####################
        smu = event_muon[event_level]
        smu = smu[:, 0]
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
            par_flav = (sjets.partonFlavour == 0) & (sjets.hadronFlavour == 0)
            genflavor = ak.values_astype(sjets.hadronFlavour + 1 * par_flav, int)
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
                elif "mu_" in histname and histname.replace("mu_", "") in smu.fields:
                    h.fill(
                        syst,
                        flatten(smu[histname.replace("mu_", "")]),
                        weight=weight,
                    )
                elif "jet" in histname and "dr" not in histname and "njet" != histname:
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
            pruned_ev.Jet = sjets
            pruned_ev.Muon = smu
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

            # Create a list of variables want to store. For objects from the PFNano file, specify as {object}_{variable}, wildcard option only accepted at the end of the string
            out_branch = np.setdiff1d(
                np.array(pruned_ev.fields), np.array(events.fields)
            )
            for kin in ["pt", "eta", "phi", "mass"]:  # ,"pfRelIso04_all","dxy", "dz"]:
                # for obj in ["Jet",  "Muon", "MET","PuppiMET"]:
                for obj in ["Muon"]:
                    if "MET" in obj and ("pt" != kin or "phi" != kin):
                        continue
                    if (obj != "Muon") and ("pfRelIso04_all" == kin or "d" in kin):
                        continue
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
