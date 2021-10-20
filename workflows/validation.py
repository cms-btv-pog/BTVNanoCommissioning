import gzip
import pickle, os, sys, mplhep as hep, numpy as np

from matplotlib.pyplot import jet

import coffea
from coffea import hist, processor
import awkward as ak
from coffea.analysis_tools import Weights
from coffea.lumi_tools import LumiMask
from coffea.btag_tools import BTagScaleFactor

import gc

class NanoProcessor(processor.ProcessorABC):
    # Define histograms
    def __init__(self):        
        # Define axes
        # Should read axes from NanoAOD config
        dataset_axis = hist.Cat("dataset", "Primary dataset")
        flav_axis = hist.Bin("flav", r"Genflavour",[0,1,4,5,6])
        cutflow_axis   = hist.Cat("cut",   "Cut")

        # Events
        nmu_axis   = hist.Bin("nmu",   r"N muons",     [0,1,2,3,4,5,6,7,8,9,10])
        njet_axis  = hist.Bin("njet",  r"N jets",      [0,1,2,3,4,5,6,7,8,9,10])
        nbjet_axis = hist.Bin("nbjet", r"N b-jets",    [0,1,2,3,4,5,6,7,8,9,10])       
        ntbjet_axis = hist.Bin("ntbjet", r"N b-jets",    [0,1,2,3,4,5,6,7,8,9,10])       
        lmupt_axis   = hist.Bin("lmupt",   r"Muon pt", 40, 0, 200)


        # Jet
        jet_pt_axis   = hist.Bin("pt",   r"Jet $p_{T}$ [GeV]", 50, 0, 500)
        jet_ptwide_axis   = hist.Bin("ptwide",   r"Jet $p_{T}$ [GeV]", [25, 30, 40, 60, 80, 100,150,200,300,500])
        jet_rawpt_axis   = hist.Bin("rawpt",   r"Jet $p_{T}$ [GeV]", 50, 0, 500)
        jet_eta_axis  = hist.Bin("eta",  r"Jet $\eta$", 25, -2.5, 2.5)
        jet_etawide_axis  = hist.Bin("etawide",  r"Jet $\eta$", [-2.5,-2.0,-1.5,-0.5,0.,0.5,1.5,2.0,2.5])
        jet_phi_axis  = hist.Bin("phi",  r"Jet $\phi$", 30, -3, 3)
        jet_mass_axis = hist.Bin("mass", r"Jet $m$ [GeV]", 50, 0, 500)
        jet_dr_axis = hist.Bin("dr", r"Jet $\Delta$R(l,j)", 20, 0, 5)
        ljpt_axis     = hist.Bin("ljpt", r"Leading jet $p_{T}$ [GeV]", 50, 0, 500)
        sljpt_axis     = hist.Bin("sljpt", r"Subleading jet $p_{T}$ [GeV]", 50, 0, 500)
        ljeta_axis = hist.Bin("ljeta",  r"Leading Jet $\eta$", 25, -2.5, 2.5)
        sljeta_axis = hist.Bin("sljeta",  r"Leading Jet $\eta$", 25, -2.5, 2.5)
        
        ljdr_axis = hist.Bin("ljdr", "Leading jet $\Delta$R(l,j)", 20,0,5)
        sljdr_axis = hist.Bin("sljdr", "Subleading jet $\Delta$R(l,j)", 20,0,5)
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
            if "flightDistance2dSig" in d or "flightDistance3dSig" in d :
                deepcsv_axes.append(hist.Bin(d, d, 101, -0.1, 100))
            elif "flightDistance2dVal" in d :
                deepcsv_axes.append(hist.Bin(d, d, 27, -0.1, 2.6))
            elif "flightDistance3dVal" in d :
                deepcsv_axes.append(hist.Bin(d, d, 51, -0.1, 5.))
            elif "trackDecayLenVal" in d :
                deepcsv_axes.append(hist.Bin(d, d, 22, -0.1, 1.))
            elif "trackPtRatio" in d or "DeltaR" in d:
                deepcsv_axes.append(hist.Bin(d, d, 30, -0.001, 0.301))
            elif "trackEtaRel" in d:
                deepcsv_axes.append(hist.Bin(d, d, 30, 0, 9))
            elif "trackJetDistVal" in d :
                deepcsv_axes.append(hist.Bin(d, d, 35, -0.08,0.0025))
            elif "trackJetPt" in d : 
                deepcsv_axes.append(hist.Bin(d, d, 50, 0.,250.))
            elif "trackPtRel" in d:
                deepcsv_axes.append(hist.Bin(d, d, 32, -0.1, 3.1))
            elif "trackSip2dSigAboveCharm" in d:
                deepcsv_axes.append(hist.Bin(d, d, 22, -5.5, 5.5))
            elif "trackSip2dSig_0" in d:
                deepcsv_axes.append(hist.Bin(d, d, 21, -5, 16))
            elif "trackSip2dSig_1" in d:
                deepcsv_axes.append(hist.Bin(d, d, 18, -5, 13))
            elif "trackSip2dSig_2" in d:
                deepcsv_axes.append(hist.Bin(d, d, 16, -6, 10))
            elif "trackSip2dSig_3" in d:
                deepcsv_axes.append(hist.Bin(d, d, 26, -6, 7))
            elif "trackSip2dSig_4" in d:
                deepcsv_axes.append(hist.Bin(d, d, 22, -6.5, 4.5))
            elif "trackSip2dSig_5" in d:
                deepcsv_axes.append(hist.Bin(d, d, 18, -7, 2))
            elif "trackSip2dValAboveCharm" in d or "trackSip3dValAboveCharm " in d:
                deepcsv_axes.append(hist.Bin(d, d, 24, -0.06, 0.06))
            elif "trackSip3dSigAboveCharm" in d :
                deepcsv_axes.append(hist.Bin(d, d, 26, -6.5, 6.5)) 
            elif "trackSip3dSig" in d : 
                deepcsv_axes.append(hist.Bin(d, d, 25, -25, 50)) 
            elif "trackSumJetEtRatio" in d:
                deepcsv_axes.append(hist.Bin(d, d, 25, 0., 1.4)) 
            elif "vertexCategory" in d:
                deepcsv_axes.append(hist.Bin(d, d, 32, -0.6,2.6)) 
            elif "vertexEnergyRatio" in d:
                deepcsv_axes.append(hist.Bin(d, d, 25, 0,2.5)) 
            elif "vertexMass" in d:
                deepcsv_axes.append(hist.Bin(d, d, 20, 0,20))
            else:
                deepcsv_axes.append(hist.Bin(d, d, 25, -0.5,0.))
        # Define similar axes dynamically
        disc_list = ['btagDeepB_b',  'btagDeepFlavB', 'btagDeepCvB','btagDeepCvL','btagDeepFlavCvB','btagDeepFlavCvL']
        syst_list = ['','SF','_up','_dn']
        varlist=[]
        btag_axes = []
        for d in disc_list:
            btag_axes.append(hist.Bin("%s" %d, "%s" %(d), 30, -0.1, 1.1)) 
            # for s in syst_list:
            #     btag_axes.append(hist.Bin("%s%s" %(d,s), "%s%s" %(d,s), 50, 0, 1))  
            #     varlist.append("%s%s" %(d,s))
                
        _hist_sf_dict={}   
        _hist_deepcsv_dict={}
        for disc, axis in zip(disc_list, btag_axes):
            _hist_deepcsv_dict["%s" %(disc)] = hist.Hist("Counts", dataset_axis,flav_axis,  axis)
            # for i in range(2):
                # _hist_sf_dict["%s_%d" %(disc,i)] = hist.Hist("Counts", dataset_axis, flav_axis, jet_etawide_axis, jet_ptwide_axis,axis)
        for deepcsv, axises in zip(deepcsv_list, deepcsv_axes):
             _hist_deepcsv_dict["%s" %(deepcsv)] = hist.Hist("Counts", dataset_axis,flav_axis,  axises)
        # for deepcsv, axises in zip(deepcsv_list, deepcsv_axes):
            # _hist_deepcsv_dict[deepcsv] = hist.Hist("Counts", dataset_axis,flav_axis,  axises)
        _hist_event_dict = {
                'njet'  : hist.Hist("Counts", dataset_axis,  njet_axis),
                'nbjet' : hist.Hist("Counts", dataset_axis,  nbjet_axis),
                'nbjet_up' : hist.Hist("Counts", dataset_axis,  nbjet_axis),
                'nbjet_dn' : hist.Hist("Counts", dataset_axis,  nbjet_axis),
                'ntbjet' : hist.Hist("Counts", dataset_axis,  ntbjet_axis),
                'nmu'   : hist.Hist("Counts", dataset_axis,  nmu_axis),
                'lmupt' : hist.Hist("Counts", dataset_axis,  lmupt_axis),
                'ljpt'  : hist.Hist("Counts", dataset_axis, flav_axis, ljpt_axis),
                'sljpt'  : hist.Hist("Counts", dataset_axis, flav_axis, sljpt_axis),
                'ljeta' : hist.Hist("Counts", dataset_axis, flav_axis, ljeta_axis),
                'sljeta' : hist.Hist("Counts", dataset_axis, flav_axis, sljeta_axis),
                'ljdr'  : hist.Hist("Counts", dataset_axis, flav_axis, ljdr_axis),
                'sljdr'  : hist.Hist("Counts", dataset_axis, flav_axis, sljdr_axis)
            }

        # self.sf_dict = list(_hist_sf_dict.keys())
        self.event_hists = list(_hist_event_dict.keys())
        self.deepcsv_hists = list(_hist_deepcsv_dict.keys())
        print(_hist_deepcsv_dict.keys())
        _hist_dict = {**_hist_event_dict,**_hist_deepcsv_dict}
        #}
        self._accumulator = processor.dict_accumulator(_hist_dict)
        self._accumulator['sumw'] = processor.defaultdict_accumulator(float)


    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.accumulator.identity()
        dataset = events.metadata['dataset']
        isRealData = not hasattr(events, "genWeight")
        
        if(isRealData):output['sumw'][dataset] += 1.
        else:output['sumw'][dataset] += ak.sum(events.genWeight)
        req_lumi=np.ones(len(events), dtype='bool')

        weights = Weights(len(events), storeIndividual=True)
        if not isRealData:
            weights.add('genweight',events.genWeight)

        ##############
        # Trigger level
        triggers = [
        # "HLT_IsoMu24",
        "HLT_IsoMu24",
        ]
        
        trig_arrs = [events.HLT[_trig.strip("HLT_")] for _trig in triggers]
        req_trig = np.zeros(len(events), dtype='bool')
        for t in trig_arrs:
            req_trig = req_trig | t

        ############
        # Event level
        
        ## Muon cuts
        # muon twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
        events.Muon = events.Muon[(events.Muon.pt > 30) & (abs(events.Muon.eta < 2.4))& (events.Muon.tightId)]  
        events.Muon = ak.pad_none(events.Muon, 1, axis=1) 
        req_muon =(ak.count(events.Muon.pt, axis=1) == 1)
        
        # ## Electron cuts
        # # electron twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
        events.Electron = events.Electron[(events.Electron.pt > 30) & (abs(events.Electron.eta) < 2.4)&(events.Electron.cutBased>3)]
        events.Electron = ak.pad_none(events.Electron, 1, axis=1) 
        req_ele = (ak.count(events.Electron.pt, axis=1) == 1)
        ## Jet cuts 
        
        req_opposite_charge = (events.Electron[:, 0].charge * events.Muon[:, 0].charge) == -1
        # req_opposite_charge = ak.fill_none(req_opposite_charge,False)
        event_jet = events.Jet[(events.Jet.pt > 25) & (abs(events.Jet.eta) <= 2.4)&(events.Jet.puId > 0) &(events.Jet.jetId>5)&(ak.all(events.Jet.metric_table(events.Muon) > 0.4, axis=2))&(ak.all(events.Jet.metric_table(events.Electron) > 0.4, axis=2))]
        req_jets = (ak.num(event_jet.puId) >= 2)
        event_level= req_jets
        # event_level=req_trig & req_lumi & req_muon & req_ele & req_jets & req_opposite_charge
        if(len(event_level)>0):event_level = ak.fill_none(event_level,False)        
        # Selected
        selev = events[event_level]    
        
        #########
        
        # Per muon
        # mu_eta   = (abs(selev.Muon.eta) < 2.4)
        # mu_pt    = selev.Muon.pt > 30
        # mu_idiso = (selev.Muon.tightId > .5)&(selev.Muon.pfRelIso04_all<0.12)
        # mu_level = mu_eta & mu_pt & mu_idiso
        # smu=selev.Muon[mu_level]

        # #Per Electron
        # el_eta = (abs(selev.Electron.eta)<2.4)
        # el_pt =  selev.Electron.pt > 30
        # el_idiso = selev.Electron.cutBased>3
        # el_level = el_eta & el_pt & el_idiso
        # sel = selev.Electron[el_level] 
        # Per jet : https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID
        jet_eta    = (abs(selev.Jet.eta) <= 2.4)
        jet_pt     = selev.Jet.pt > 25 
        jet_pu     = (selev.Jet.puId > 0) &(selev.Jet.jetId>5)
        # jet_dr     = (ak.all(selev.Jet.metric_table(smu) > 0.4, axis=2) & ak.all(selev.Jet.metric_table(sel) > 0.4, axis=2) )
        # jet_clean  = (selev.Jet.btagDeepB>0.) & (selev.Jet.btagDeepB<1.) & (selev.Jet.btagDeepC>0.) & (selev.Jet.btagDeepC<1.) & (selev.Jet.btagDeepFlavB>0.) & (selev.Jet.btagDeepFlavB<1.) & (selev.Jet.btagDeepFlavC>0.) & (selev.Jet.btagDeepFlavC<1.)
        
        jet_level  = jet_pu & jet_eta & jet_pt # & jet_dr
        sjets  = selev.Jet[jet_level]
        sel_jets = sjets
        sjets = sjets[:,:2]
        # b-tag twiki : https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
        bjet_disc  = selev.Jet.btagDeepB > 0.4941 
        bjet_level = jet_level & bjet_disc

        sbjets = selev.Jet[bjet_level]
        if not isRealData: stbjets = sbjets[sbjets.hadronFlavour==5]
        
        if isRealData :
            genflavor = ak.zeros_like(sjets.pt)
        else:
            par_flav = (sjets.partonFlavour == 0 ) & (sjets.hadronFlavour==0)
            # genflavor = sjets.hadronFlavour + 1*par_flav 
            genflavor=sjets.hadronFlavour
            genweiev=ak.flatten(ak.broadcast_arrays(weights.weight()[event_level],sjets['pt'])[0])
            
            
        for histname, h in output.items():
            if histname in self.deepcsv_hists:
                if(isRealData):
                    fields = {l: ak.flatten(sjets[l], axis=None) for l in h.fields if l in dir(sjets)}
                    h.fill(dataset=dataset,flav=5, **fields)
                else:
                    fields = {l: ak.flatten(sjets[histname]) for l in h.fields if l in dir(sjets)}
                    genweiev=ak.flatten(ak.broadcast_arrays(weights.weight()[event_level],sjets['pt'])[0])
                    h.fill(dataset=dataset,flav=ak.flatten(genflavor), **fields,weight=genweiev)
       
            
        # if not isRealData:
        #     output['btagDeepFlavB_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepFlavB=sjets[:,0].btagDeepFlavB,weight=weights.weight()[event_level])
        #     output['btagDeepB_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepB=sjets[:,0].btagDeepB_b,weight=weights.weight()[event_level])
        #     output['deepcsv_CvB_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepcsv_CvB=sjets[:,0].btagDeepCvB,weight=weights.weight()[event_level])
        #     output['deepcsv_CvL_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepcsv_CvL=sjets[:,0].btagDeepCvL,weight=weights.weight()[event_level])
        #     output['deepflav_CvB_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepflav_CvB=sjets[:,0].btagDeepFlavCvB,weight=weights.weight()[event_level])
        #     output['deepflav_CvL_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepflav_CvL=sjets[:,0].btagDeepFlavCvL,weight=weights.weight()[event_level])
        #     output['btagDeepFlavB_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepFlavB=sjets[:,1].btagDeepFlavB,weight=weights.weight()[event_level])
        #     output['btagDeepB_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepB=sjets[:,1].btagDeepB_b,weight=weights.weight()[event_level])
        #     output['deepcsv_CvB_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepcsv_CvB=sjets[:,1].btagDeepCvB,weight=weights.weight()[event_level])
        #     output['deepcsv_CvL_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepcsv_CvL=sjets[:,1].btagDeepCvL,weight=weights.weight()[event_level])
        #     output['deepflav_CvB_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepflav_CvB=sjets[:,1].btagDeepFlavCvB,weight=weights.weight()[event_level])
        #     output['deepflav_CvL_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepflav_CvL=sjets[:,1].btagDeepFlavCvL,weight=weights.weight()[event_level])            
        # else:
        #     output['btagDeepFlavB_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepFlavB=sjets[:,0].btagDeepFlavB)
        #     output['btagDeepB_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepB=sjets[:,0].btagDeepB_b)
        #     output['deepcsv_CvB_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepcsv_CvB=sjets[:,0].btagDeepCvB)
        #     output['deepcsv_CvL_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepcsv_CvL=sjets[:,0].btagDeepCvL)
        #     output['deepflav_CvB_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepflav_CvB=sjets[:,0].btagDeepFlavCvB)
        #     output['deepflav_CvL_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepflav_CvL=sjets[:,0].btagDeepFlavCvL)
        #     output['btagDeepFlavB_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepFlavB=sjets[:,1].btagDeepFlavB)
        #     output['btagDeepB_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepB=sjets[:,1].btagDeepB_b)
        #     output['deepcsv_CvB_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepcsv_CvB=sjets[:,1].btagDeepCvB)
        #     output['deepcsv_CvL_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepcsv_CvL=sjets[:,1].btagDeepCvL)
        #     output['deepflav_CvB_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepflav_CvB=sjets[:,1].btagDeepFlavCvB)
        #     output['deepflav_CvL_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepflav_CvL=sjets[:,1].btagDeepFlavCvL)

        def flatten(ar): # flatten awkward into a 1d array to hist
            return ak.flatten(ar, axis=None)

        def num(ar):
            return ak.num(ak.fill_none(ar[~ak.is_none(ar)], 0), axis=0)
        if(isRealData):
            # output['njet'].fill(dataset=dataset,  njet=flatten(ak.num(sel_jets)))
            # output['nbjet'].fill(dataset=dataset,  nbjet=flatten(ak.num(sbjets)))
            # output['ntbjet'].fill(dataset=dataset,  ntbjet=flatten(ak.num(sbjets)))
            output['ljpt'].fill(dataset=dataset, flav=0, ljpt=flatten(sjets[:,0].pt))
            output['sljpt'].fill(dataset=dataset, flav=0, sljpt=flatten(sjets[:,1].pt))
            output['ljeta'].fill(dataset=dataset, flav=0, ljeta=flatten(sjets[:,0].eta))
            output['sljeta'].fill(dataset=dataset, flav=0, sljeta=flatten(sjets[:,1].eta))
        else:
            output['njet'].fill(dataset=dataset,  njet=flatten(ak.num(sel_jets)),weight=weights.weight()[event_level])
            # output['nbjet'].fill(dataset=dataset,  nbjet=flatten(ak.num(sbjets)),weight=weights.weight()[event_level])
            # output['nbjet_up'].fill(dataset=dataset,  nbjet=flatten(ak.num(sbjets)),weight=weights.weight()[event_level]*csvsfs_b_up_lj*csvsfs_b_up_sj)
            # output['nbjet_dn'].fill(dataset=dataset,  nbjet=flatten(ak.num(sbjets)),weight=weights.weight()[event_level]*csvsfs_b_dn_lj*csvsfs_b_dn_sj)
            output['ntbjet'].fill(dataset=dataset,  ntbjet=flatten(ak.num(stbjets)),weight=weights.weight()[event_level])
            output['ljpt'].fill(dataset=dataset, flav=flatten(sjets[:,0].hadronFlavour), ljpt=flatten(sjets[:,0].pt),weight=weights.weight()[event_level])
            output['sljpt'].fill(dataset=dataset, flav=flatten(sjets[:,1].hadronFlavour), sljpt=flatten(sjets[:,1].pt),weight=weights.weight()[event_level])
            output['ljeta'].fill(dataset=dataset, flav=flatten(sjets[:,0].hadronFlavour), ljeta=flatten(sjets[:,0].eta),weight=weights.weight()[event_level])
            output['sljeta'].fill(dataset=dataset, flav=flatten(sjets[:,1].hadronFlavour), sljeta=flatten(sjets[:,1].eta),weight=weights.weight()[event_level])
        gc.collect()
        schedule.every(20).minutes.do(dosomething)

        return output

    def postprocess(self, accumulator):
        return accumulator
