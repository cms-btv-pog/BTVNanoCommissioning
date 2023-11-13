import collections, numpy as np, awkward as ak
from coffea import processor
import correctionlib
from coffea.analysis_tools import Weights
"""
from BTVNanoCommissioning.utils.correction import load_lumi
from BTVNanoCommissioning.utils.AK4_parameters import correction_config
"""
from BTVNanoCommissioning.helpers.definitions import definitions,SV_definitions
from BTVNanoCommissioning.utils.selection import jet_cut
from BTVNanoCommissioning.helpers.func import flatten, update
from BTVNanoCommissioning.utils.histogrammer import histogrammer
from BTVNanoCommissioning.helpers.update_branch import missing_branch, add_jec
from BTVNanoCommissioning.utils.correction import (
    load_SF,
    muSFs,
    eleSFs,
    puwei,
    btagSFs,
    JME_shifts,
    Roccor_shifts,
)

class NanoProcessor(processor.ProcessorABC):
    # Define histograms
    def __init__(
        self, 
        year="2018", 
        campaign="2018_UL", 
        isCorr=False, 
        isJERC=False,
        isSyst=False, 
        isArray=False, 
        noHist=False,
        chunksize=75000,
    ):
        self._year = year
        self._campaign = campaign
        self.isCorr = isCorr
        self.isJERC = isJERC
        self.isSyst = isSyst
        self.isArray = isArray
        self.noHist = noHist
        self.chunksize = chunksize
        ## Load corrections  
        #if isCorr:
        self.SF_map = load_SF(self._campaign,True)
        #if isJERC:
            #self._jet_factory = load_jmefactory(self._campaign)
         

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        #output = self.make_output()
        isRealData = not hasattr(events, "genWeight")
        dataset = events.metadata["dataset"]
        events = missing_branch(events)
        shifts = []
        if "JME" in self.SF_map.keys() and self.isJERC:
            syst_JERC = True if self.isSyst != None else False
            if self.isSyst == "JERC_split":
                syst_JERC = "split"
                shifts = JME_shifts(
                    shifts, self.SF_map, events, self._campaign, isRealData, syst_JERC
            )
        else:
            shifts = [
                ({"Jet": events.Jet, "MET": events.MET, "Muon": events.Muon}, None)
            ]

        return processor.accumulate(
            self.process_shift(update(events, collections), name)
            for collections, name in shifts
        )
    """
        if self.isJERC:
            events = add_jec(events, self._campaign, self._jet_factory)
        if isRealData:
            output["sumw"] = len(events)
        else:
            output["sumw"] = ak.sum(events.genWeight)
      """      
            
    def process_shift(self, events, shift_name):
        isRealData = not hasattr(events, "genWeight")
        dataset = events.metadata["dataset"]
        _hist_event_dict = {"": None} if self.noHist else histogrammer("QCD")
        if _hist_event_dict == None:
            _hist_event_dict[""]
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
        events.Jet = events.Jet[jet_cut(events, self._campaign)]
        req_jets = ak.count(events.Jet.pt, axis=1) >= 1

        event_level = ak.fill_none(
            req_trig & req_jets, False
        )
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
        matched_JetSVs=  selev.Jet[selev.JetSVs.jetIdx]
        lj_matched_JetSVs = matched_JetSVs[selev.JetSVs.jetIdx==0]
        lj_SVs = selev.JetSVs[selev.JetSVs.jetIdx==0]
        nJetSVs = ak.count(lj_SVs.pt, axis=1)

        ####################
        # Weight & Geninfo #
        ####################
        weights = Weights(len(events), storeIndividual=True)
        if not isRealData:
            weights.add("genweight", events.genWeight)
            par_flav = (sjets.partonFlavour == 0) & (sjets.hadronFlavour == 0)
            genflavor = sjets.hadronFlavour + 1 * par_flav
            genweiev = ak.flatten(
            ak.broadcast_arrays(weights.weight()[event_level], sjets["pt"])[0]
            )
            lj_matched_JetSVs_par_flav = (lj_matched_JetSVs.partonFlavour == 0) & (lj_matched_JetSVs.hadronFlavour == 0)
            lj_matched_JetSVs_genflav = lj_matched_JetSVs.hadronFlavour + 1 * lj_matched_JetSVs_par_flav
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
            pseval = correctionlib.CorrectionSet.from_file(f"src/BTVNanoCommissioning/data/Prescales/ps_weight_JSON_PFJet140.json")
            psweight = pseval['prescaleWeight'].evaluate(events.run,"PFJet140",ak.values_astype(events.luminosityBlock, np.float32))
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
            if self.isSyst == None and syst != "nominal":
                break
            if self.noHist:
                break
            weight = (
                weights.weight()
                if syst == "nominal" or syst == shift_name
                else weights.weight(modifier=syst)
            )
            for histname, h in output.items():
                if "DeepCSV" in histname and "btag" not in histname:
                    h.fill(
                        flatten(genflavor),
                        flatten(sjets[histname]),
                        weight=flatten(
                            ak.broadcast_arrays(weights.weight()[event_level], sjets["pt"])[
                                0
                            ]
                        ),
                    )
                elif "jet" in histname and histname.replace("jet0_", "") in sjets.fields: 
                    h.fill(
                        syst,
                        flatten(genflavor[:, 0]),
                        flatten(sjets[:, 0][histname.replace(f"jet0_", "")]),
                        weight=weights.weight()[event_level],
                    )
                elif "JetSVs_" in histname:
                    h.fill(
                        flatten(lj_matched_JetSVs_genflav),
                        flatten(lj_SVs[histname.replace("JetSVs_", "")]),
                        weight=flatten(ak.broadcast_arrays(weights.weight()[event_level],lj_matched_JetSVs["pt"])[0]),
                    )
                elif "btagDeep" in histname:
                    for i in range(1):
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
                                weight=weights.weight()[event_level],
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
                                        weight=weights.weight()[event_level]
                                        * disc_list[histname.replace(f"_{i}", "")][i][syst],
                                    )            
        #output["njet"].fill(njet, weight=weights.weight()[event_level])
        #output["nJetSVs"].fill(nJetSVs, weight=weights.weight()[event_level]) 
        
        #output["dr_SVjet0"].fill(flatten(lj_matched_JetSVs_genflav),
        #                         flatten(abs(lj_SVs.deltaR)-0.1),
        #                         weight=flatten(ak.broadcast_arrays(weights.weight()[event_level],lj_matched_JetSVs["pt"])[0]),)
        #output["unw_PV_npvsGood"].fill(flatten(events.PV.npvsGood))
        #output["weighted_PV_npvsGood"].fill(flatten(events.PV.npvsGood),weight=weights.weight())

        #if not isRealData:
            #output["genWeight"].fill(flatten(events.genWeight))
            #output["unw_Pileup_nTrueInt"].fill(flatten(events.Pileup.nTrueInt))
            #output["weighted_Pileup_nTrueInt"].fill(flatten(events.Pileup.nTrueInt),weight=weights.weight())

        return {dataset: output}

    def postprocess(self, accumulator):
        return accumulator
