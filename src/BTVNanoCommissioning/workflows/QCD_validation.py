import collections, numpy as np, awkward as ak
from coffea import processor
from coffea.analysis_tools import Weights
from BTVNanoCommissioning.utils.selection import jet_cut
from BTVNanoCommissioning.helpers.func import flatten, update
from BTVNanoCommissioning.utils.histogrammer import histogrammer
from BTVNanoCommissioning.helpers.update_branch import missing_branch, add_jec
from BTVNanoCommissioning.utils.correction import (
    load_SF,
    load_lumi,
    muSFs,
    eleSFs,
    puwei,
    btagSFs,
    JME_shifts,
    Roccor_shifts,
    jetveto,
)
import correctionlib


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
        if "JME" in self.SF_map.keys() or "jetveto" in self.SF_map.keys():
            syst_JERC = self.isSyst
            if self.isSyst == "JERC_split":
                syst_JERC = "split"
            shifts = JME_shifts(
                shifts, self.SF_map, events, self._campaign, isRealData, False, True
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
        isRealData = not hasattr(events, "genWeight")
        dataset = events.metadata["dataset"]
        _hist_event_dict = {"": None} if self.noHist else histogrammer(events, "QCD")
        if _hist_event_dict == None:
            _hist_event_dict[""]
        output = {
            "sumw": processor.defaultdict_accumulator(float),
            **_hist_event_dict,
        }
        if "JME" in self.SF_map.keys() or "jetveto" in self.SF_map.keys():
            events.Jet = update(events.Jet, {"veto": jetveto(events, self.SF_map)})

        if isRealData:
            output["sumw"] = len(events)
        else:
            output["sumw"] = ak.sum(events.genWeight)

        ####################
        #    Selections    #
        ####################
        ## HLT
        triggers = [
            "PFJet140",
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

        ## Jet cuts
        events.Jet = events.Jet[
            jet_cut(events, self._campaign) & (events.Jet.veto != 1)
        ]
        req_jets = ak.count(events.Jet.pt, axis=1) >= 1

        event_level = ak.fill_none(req_trig & req_jets, False)
        if len(events[event_level]) == 0:
            return {dataset: output}

        ####################
        # Selected objects #
        ####################
        sjets = events.Jet[event_level]
        njet = ak.count(sjets.pt, axis=1)
        ###############
        # Selected SV #
        ###############
        selev = events[event_level]
        matched_JetSVs = selev.Jet[selev.JetSVs.jetIdx]
        lj_matched_JetSVs = matched_JetSVs[selev.JetSVs.jetIdx == 0]
        lj_SVs = selev.JetSVs[selev.JetSVs.jetIdx == 0]
        nJetSVs = ak.count(lj_SVs.pt, axis=1)

        ####################
        # Weight & Geninfo #
        ####################
        weights = Weights(len(selev), storeIndividual=True)
        if not isRealData:
            weights.add("genweight", selev.genWeight)
            par_flav = (sjets.partonFlavour == 0) & (sjets.hadronFlavour == 0)
            genflavor = sjets.hadronFlavour + 1 * par_flav
            # genweiev = ak.flatten(
            # ak.broadcast_arrays(weights.weight()[event_level], sjets["pt"])[0]
            # )
            lj_matched_JetSVs_par_flav = (lj_matched_JetSVs.partonFlavour == 0) & (
                lj_matched_JetSVs.hadronFlavour == 0
            )
            lj_matched_JetSVs_genflav = (
                lj_matched_JetSVs.hadronFlavour + 1 * lj_matched_JetSVs_par_flav
            )
            syst_wei = True if self.isSyst != None else False
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

        if isRealData:
            if self._year == "2022":
                run_num = "355374_362760"
            elif self._year == "2023":
                run_num = "366727_370790"
            pseval = correctionlib.CorrectionSet.from_file(
                f"src/BTVNanoCommissioning/data/Prescales/ps_weight_{triggers[0]}_run{run_num}.json"
            )
            print(triggers[0])
            psweight = pseval["prescaleWeight"].evaluate(
                selev.run,
                f"HLT_{triggers[0]}",
                ak.values_astype(selev.luminosityBlock, np.float32),
            )
            weights.add("psweight", psweight)
            genflavor = ak.zeros_like(sjets.pt)
            lj_matched_JetSVs_genflav = ak.zeros_like(lj_matched_JetSVs.pt)

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
                    "DeepCSV" in histname
                    and "btag" not in histname
                    and histname in events.Jet.fields
                ):
                    h.fill(
                        syst,
                        flatten(genflavor),
                        flatten(sjets[histname]),
                        weight=flatten(ak.broadcast_arrays(weight, sjets["pt"])[0]),
                    )
                elif (
                    "jet" in histname and histname.replace("jet0_", "") in sjets.fields
                ):
                    h.fill(
                        syst,
                        flatten(genflavor[:, 0]),
                        flatten(sjets[:, 0][histname.replace(f"jet0_", "")]),
                        weight=weight,
                    )
                elif "JetSVs_" in histname:
                    h.fill(
                        syst,
                        flatten(lj_matched_JetSVs_genflav),
                        flatten(lj_SVs[histname.replace("JetSVs_", "")]),
                        weight=flatten(
                            ak.broadcast_arrays(weight, lj_matched_JetSVs["pt"])[0]
                        ),
                    )
                elif (
                    "btag" in histname
                    and "0" in histname
                    and histname.replace("_0", "") in events.Jet.fields
                ):
                    sel_jet = sjets[:, 0]
                    if syst == "nominal":
                        h.fill(
                            syst="noSF",
                            flav=genflavor[:, 0],
                            discr=np.where(
                                sel_jet[histname.replace("_0", "")] < 0,
                                -0.2,
                                sel_jet[histname.replace("_0", "")],
                            ),
                            weight=weight,
                        )
                    if not isRealData and "btag" in self.SF_map.keys():
                        h.fill(
                            syst=syst,
                            flav=genflavor[:, 0],
                            discr=np.where(
                                sel_jet[histname.replace("_0", "")] < 0,
                                -0.2,
                                sel_jet[histname.replace("_0", "")],
                            ),
                            weight=weight,
                        )

            output["njet"].fill(syst, njet, weight=weight)
            output["nJetSVs"].fill(syst, nJetSVs, weight=weight)

            output["dr_SVjet0"].fill(
                syst,
                flatten(lj_matched_JetSVs_genflav),
                flatten(abs(lj_SVs.deltaR) - 0.1),
                weight=flatten(ak.broadcast_arrays(weight, lj_matched_JetSVs["pt"])[0]),
            )
            output["npvs"].fill(syst, flatten(selev.PV.npvsGood), weight=weight)
            if not isRealData:
                if syst == "nominal":
                    output["pu"].fill(
                        "noPU",
                        flatten(selev.Pileup.nTrueInt),
                        weight=np.ones_like(weight),
                    )
                    output["npvs"].fill(
                        "noPU",
                        flatten(selev.Pileup.nTrueInt),
                        weight=np.ones_like(weight),
                    )
                    output["pu"].fill(syst, flatten(selev.PV.npvsGood), weight=weight)

        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator
