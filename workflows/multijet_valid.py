import coffea
from coffea import hist, processor
import numpy as np
import awkward as ak


class NanoProcessor(processor.ProcessorABC):
    # Define histograms
    def __init__(self):        
        # Define axes
        # Should read axes from NanoAOD config
        dataset_axis = hist.Cat("dataset", "Primary dataset")
        flav_axis = hist.Bin("flav", r"Genflavour",[0,4,5,6])
        cutflow_axis   = hist.Cat("cut",   "Cut")

        # Events
        njet_axis  = hist.Bin("njet",  r"N jets",      [0,1,2,3,4,5,6,7,8,9,10])
        nbjet_axis = hist.Bin("nbjet", r"N b-jets",    [0,1,2,3,4,5,6,7,8,9,10])

        # Jet
        jet_pt_axis   = hist.Bin("pt",   r"Jet $p_{T}$ [GeV]", 100, 20, 400)
        jet_eta_axis  = hist.Bin("eta",  r"Jet $\eta$", 60, -3, 3)
        jet_phi_axis  = hist.Bin("phi",  r"Jet $\phi$", 60, -3, 3)
        jet_mass_axis = hist.Bin("mass", r"Jet $m$ [GeV]", 100, 0, 50)
        fatjet_pt_axis   = hist.Bin("fatjet_pt",   r"FatJet $p_{T}$ [GeV]", 50, 350, 1000)
        fatjet_eta_axis  = hist.Bin("fatjet_eta",  r"FatJet $\eta$", 60, -3, 3)
        fatjet_phi_axis  = hist.Bin("fatjet_phi",  r"FatJet $\phi$", 60, -3, 3)
        fatjet_mass_axis = hist.Bin("fatjet_mass", r"FatJet $m$ [GeV]", 100, 0, 1000)
        subjet_pt_axis   = hist.Bin("subjet_pt",   r"subJet $p_{T}$ [GeV]", 50, 350, 1000)
        subjet_eta_axis  = hist.Bin("subjet_eta",  r"subJet $\eta$", 60, -3, 3)
        subjet_phi_axis  = hist.Bin("subjet_phi",  r"subJet $\phi$", 60, -3, 3)
        subjet_mass_axis = hist.Bin("subjet_mass", r"subJet $m$ [GeV]", 100, 0, 1000)
        ljpt_axis     = hist.Bin("ljpt", r"Leading jet $p_{T}$ [GeV]", 100, 20, 400)
        sljpt_axis     = hist.Bin("sljpt", r"Subleading jet $p_{T}$ [GeV]", 100, 20, 400)

        # Define similar axes dynamically
        disc_list = ["btagCMVA", "btagCSVV2", 'btagDeepB', 'btagDeepC', 'btagDeepFlavB', 'btagDeepFlavC','btagDeepB_bb']
        ddx_list = ["btagDDBvLV2","btagDDCvBV2","btagDDCvLV2"]
        btag_axes = []
        for d in disc_list:
            btag_axes.append(hist.Bin(d, d, 50, 0, 1))     
        deepddx_list = ["DDX_jetNTracks","DDX_jetNSecondaryVertices","DDX_tau1_trackEtaRel_0","DDX_tau1_trackEtaRel_1","DDX_tau1_trackEtaRel_2","DDX_tau2_trackEtaRel_0","DDX_tau2_trackEtaRel_1","DDX_tau2_trackEtaRel_3","DDX_tau1_flightDistance2dSig","DDX_tau2_flightDistance2dSig","DDX_tau1_vertexDeltaR","DDX_tau1_vertexEnergyRatio","DDX_tau2_vertexEnergyRatio","DDX_tau1_vertexMass","DDX_tau2_vertexMass","DDX_trackSip2dSigAboveBottom_0","DDX_trackSip2dSigAboveBottom_1","DDX_trackSip2dSigAboveCharm","DDX_trackSip3dSig_0","DDX_tau1_trackSip3dSig_0","DDX_tau1_trackSip3dSig_1","DDX_trackSip3dSig_1","DDX_tau2_trackSip3dSig_0","DDX_tau2_trackSip3dSig_1"]
        deepcsv_list = [
        "DeepCSV_trackDecayLenVal_0", "DeepCSV_trackDecayLenVal_1", "DeepCSV_trackDecayLenVal_2", "DeepCSV_trackDecayLenVal_3", "DeepCSV_trackDecayLenVal_4", "DeepCSV_trackDecayLenVal_5", 
        "DeepCSV_trackDeltaR_0", "DeepCSV_trackDeltaR_1", "DeepCSV_trackDeltaR_2", "DeepCSV_trackDeltaR_3", "DeepCSV_trackDeltaR_4", "DeepCSV_trackDeltaR_5",
        "DeepCSV_trackEtaRel_0","DeepCSV_trackEtaRel_1","DeepCSV_trackEtaRel_2","DeepCSV_trackEtaRel_3", 	
        "DeepCSV_trackJetDistVal_0","DeepCSV_trackJetDistVal_1","DeepCSV_trackJetDistVal_2","DeepCSV_trackJetDistVal_3","DeepCSV_trackJetDistVal_4","DeepCSV_trackJetDistVal_5", 
        "DeepCSV_trackPtRatio_0","DeepCSV_trackPtRatio_1","DeepCSV_trackPtRatio_2","DeepCSV_trackPtRatio_3","DeepCSV_trackPtRatio_4","DeepCSV_trackPtRatio_5", 
        "DeepCSV_trackPtRel_0", "DeepCSV_trackPtRel_1","DeepCSV_trackPtRel_2","DeepCSV_trackPtRel_3","DeepCSV_trackPtRel_4","DeepCSV_trackPtRel_5",
        "DeepCSV_trackSip3dSig_0","DeepCSV_trackSip3dSig_1","DeepCSV_trackSip3dSig_2","DeepCSV_trackSip3dSig_3","DeepCSV_trackSip3dSig_4","DeepCSV_trackSip3dSig_5",
        "DeepCSV_trackSip2dSig_0","DeepCSV_trackSip2dSig_1","DeepCSV_trackSip2dSig_2","DeepCSV_trackSip2dSig_3","DeepCSV_trackSip2dSig_4","DeepCSV_trackSip2dSig_5",
        "DeepCSV_trackSip2dValAboveCharm","DeepCSV_trackSip2dSigAboveCharm","DeepCSV_trackSip3dValAboveCharm","DeepCSV_trackSip3dSigAboveCharm",
        "DeepCSV_vertexCategory","DeepCSV_vertexEnergyRatio", "DeepCSV_vertexJetDeltaR","DeepCSV_vertexMass", 
        "DeepCSV_flightDistance2dVal","DeepCSV_flightDistance2dSig","DeepCSV_flightDistance3dVal","DeepCSV_flightDistance3dSig","DeepCSV_trackJetPt", 
        "DeepCSV_jetNSecondaryVertices","DeepCSV_jetNSelectedTracks","DeepCSV_jetNTracksEtaRel","DeepCSV_trackSumJetEtRatio","DeepCSV_trackSumJetDeltaR","DeepCSV_vertexNTracks"]
    
        deepcsv_axes = []
        for d in deepcsv_list:
            if "jetN" in d:
                deepcsv_axes.append(hist.Bin(d, d, 16, 0, 15))
            elif "vertexN" in d:
                deepcsv_axes.append(hist.Bin(d, d, 16, 0, 15))
            elif "vertexCategory" in d:
                deepcsv_axes.append(hist.Bin(d, d, 3, 0, 2))
            elif "EtaRel" in d:
                deepcsv_axes.append(hist.Bin(d, d, 50, 0, 10))
            elif "ValAboveCharm" in d:
                deepcsv_axes.append(hist.Bin(d, d, 50, 0, 0.3))
            elif "EnergyRatio" in d:
                deepcsv_axes.append(hist.Bin(d, d, 50, 0, 1))
            elif "trackJetDistVal" in d:
                deepcsv_axes.append(hist.Bin(d, d, 50, -0.1, 0.1))
            elif "trackPtRatio" in d:
                deepcsv_axes.append(hist.Bin(d, d, 50, 0, 0.3))
            elif "DeltaR" in d: 
                deepcsv_axes.append(hist.Bin(d, d, 50, 0, 0.3))
            elif "jetNSecondaryVertices" in d:
                deepcsv_axes.append(hist.Bin(d, d, 5, 0, 5))
            elif "vertexCategory" in d:
                deepcsv_axes.append(hist.Bin(d, d, 5, 0, 5))
            elif "jetNSelectedTracks" in d:
                deepcsv_axes.append(hist.Bin(d, d, 5, 0, 5))
            elif "vertexNTracks" in d:
                deepcsv_axes.append(hist.Bin(d, d, 5, 0, 5))
            else:
                deepcsv_axes.append(hist.Bin(d, d, 50, 0, 5.))
        deepddx_axes = []
        for d in deepddx_list:
            if "NTracks" in d:
                deepddx_axes.append(hist.Bin(d, d, 10, 0, 10))
            elif "NSecondaryVertices" in d:
                deepddx_axes.append(hist.Bin(d, d, 10, 0, 10))
            elif "EtaRel" in d:
                deepddx_axes.append(hist.Bin(d, d, 50, 0, 10))
            elif "DeltaR" in d: 
                deepddx_axes.append(hist.Bin(d, d, 50, 0, 0.3)) 
            elif "ValAboveCharm" in d:
                deepddx_axes.append(hist.Bin(d, d, 50, 0, 0.3))
            elif "EnergyRatio" in d:
                deepddx_axes.append(hist.Bin(d, d, 50, 0, 1))
            else :
                deepddx_axes.append(hist.Bin(d, d, 50, 0, 5.))
        
        # Define histograms from axes
        _hist_jet_dict = {
                'jet_pt'  : hist.Hist("Counts", dataset_axis, flav_axis,jet_pt_axis),
                'jet_eta' : hist.Hist("Counts", dataset_axis, flav_axis,jet_eta_axis),
                'jet_phi' : hist.Hist("Counts", dataset_axis, flav_axis,jet_phi_axis),
                'jet_mass': hist.Hist("Counts", dataset_axis, flav_axis,jet_mass_axis)
            }
        _hist_deepcsv_dict = {
                'pt'  : hist.Hist("Counts", dataset_axis, flav_axis,jet_pt_axis),
                'eta' : hist.Hist("Counts", dataset_axis, flav_axis,jet_eta_axis),
                'phi' : hist.Hist("Counts", dataset_axis, flav_axis,jet_phi_axis),
                'mass': hist.Hist("Counts", dataset_axis, flav_axis,jet_mass_axis),
            }
        _hist_deepddx_dict = {
                'pt'  : hist.Hist("Counts", dataset_axis, flav_axis,jet_pt_axis),
                'eta' : hist.Hist("Counts", dataset_axis, flav_axis,jet_eta_axis),
                'phi' : hist.Hist("Counts", dataset_axis, flav_axis,jet_phi_axis),
                'mass': hist.Hist("Counts", dataset_axis, flav_axis,jet_mass_axis),
            }
       
        # Generate some histograms dynamically
        
        for disc, axis in zip(ddx_list, btag_axes):
            _hist_deepddx_dict[disc] = hist.Hist("Counts", dataset_axis,  flav_axis, axis)
        for deepcsv, axises in zip(deepcsv_list, deepcsv_axes):
            _hist_deepcsv_dict[deepcsv] = hist.Hist("Counts", dataset_axis,flav_axis,  axises)
        _hist_event_dict = {
                'njet'  : hist.Hist("Counts", dataset_axis, njet_axis),
                'nbjet' : hist.Hist("Counts", dataset_axis, nbjet_axis),
                # 'ljpt'  : hist.Hist("Counts", dataset_axis, ljpt_axis),
                # 'sljpt'  : hist.Hist("Counts", dataset_axis, sljpt_axis),
            }
        self.jet_hists = list(_hist_jet_dict.keys())
        self.deepcsv_hists = list(_hist_deepcsv_dict.keys())
        self.deepddx_hists = list(_hist_deepddx_dict.keys())
        self.event_hists = list(_hist_event_dict.keys())
        _hist_dict = {**_hist_deepcsv_dict,**_hist_event_dict}
        self._accumulator = processor.dict_accumulator(_hist_dict)
        self._accumulator['sumw'] = processor.defaultdict_accumulator(float)


    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.accumulator.identity()
        # print(output.items())
        dataset = events.metadata['dataset']
        # if(len(events)==0):continue
        isRealData = not hasattr(events, "genWeight")
        if(isRealData):output['sumw'][dataset] += 1.
        else:output['sumw'][dataset] += ak.sum(events.genWeight)
        
        ##############
        # Trigger level
        triggers = [
        # "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
        # "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",   
        "HLT_PFJet40" 
        ]
        
        trig_arrs = [events.HLT[_trig.strip("HLT_")] for _trig in triggers]
        req_trig = np.zeros(len(events), dtype='bool')
        for t in trig_arrs:
            req_trig = req_trig | t

        ############
        # Event level        
        
        ## Jet cuts
        events.Jet = events.Jet[(events.Jet.pt > 25) & (abs(events.Jet.eta) <= 2.5)]
        req_jets = (ak.count(events.Jet.pt, axis=1) >= 1)    
        
        
        
        event_level = req_trig & req_jets
        
        # Selected
        selev = events[event_level]    
        
        #########
        
        

        # Per jet : https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID
        jet_eta    = (abs(selev.Jet.eta) <= 2.4)
        jet_pt     = selev.Jet.pt > 40
        jet_pu     = selev.Jet.puId > 6
        jet_level  = jet_pu & jet_eta & jet_pt
        
        # b-tag twiki : https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
        bjet_disc  = selev.Jet.btagDeepB > 0.7264 # L=0.0494, M=0.2770, T=0.7264
        bjet_level = jet_level & bjet_disc
        
        sjets  = selev.Jet[jet_level]
        sbjets = selev.Jet[bjet_level]
        sjets_b = sjets
        sjets_c = sjets
        sjets_l = sjets
        if(not isRealData):
            sjets_b = selev.Jet[(jet_level)&(abs(selev.Jet.hadronFlavour)==5)]
            sjets_c = selev.Jet[(jet_level)&(abs(selev.Jet.hadronFlavour)==4)]
            sjets_l = selev.Jet[(jet_level)&(abs(selev.Jet.hadronFlavour)<4)]
        
        if isRealData :
            genflavor = ak.zeros_like(sjets.pt)
        else:
            genflavor = selev.Jet.hadronFlavour
            #weight = prepro_fcn_to_match_other_inputs(events.genWeight)

        for histname, h in output.items():
            if(isRealData):
                if (histname in self.jet_hists) or (histname in self.deepcsv_hists):
                        fields = {l: ak.flatten(sjets_b[l.replace('jet_','')], axis=None) for l in h.fields if l.replace('jet_','') in dir(sjets_b)}
                        h.fill(dataset=dataset,flav=5, **fields)
            else:
                if (histname in self.jet_hists) or (histname in self.deepcsv_hists):
                    fields = {l: ak.flatten(ak.fill_none(sjets[histname],np.nan)) for l in h.fields if l in dir(sjets)}
                    genweiev=ak.flatten(ak.broadcast_arrays(selev.genWeight,ak.fill_none(sjets[histname],np.nan))[0])
                    h.fill(dataset=dataset,flav=ak.flatten(ak.fill_none(genflavor,np.nan)), **fields,weight=genweiev)


        def flatten(ar): # flatten awkward into a 1d array to hist
            return ak.flatten(ar, axis=None)

        def num(ar):
            return ak.num(ak.fill_none(ar[~ak.is_none(ar)], 0), axis=0)

        output['njet'].fill(dataset=dataset,  njet=flatten(ak.num(sjets)))
        output['nbjet'].fill(dataset=dataset, nbjet=flatten(ak.num(sbjets)))
        # output['ljpt'].fill(dataset=dataset,  ljpt=flatten(selev.Jet[:, 0].pt))
        # output['sljpt'].fill(dataset=dataset,  sljpt=flatten(selev.Jet[:, 1].pt))

        return output

    def postprocess(self, accumulator):
        return accumulator
