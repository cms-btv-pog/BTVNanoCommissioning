import gzip
import pickle, os, sys, mplhep as hep, numpy as np
import re

from matplotlib.pyplot import jet

import coffea
from coffea import hist, processor
import awkward as ak
from coffea.analysis_tools import Weights

from BTVNanoCommissioning.utils.correction import lumiMasks, compiled, eleSFs,muSFs,deepcsvb_sf,deepcsvc_sf,deepjetb_sf,deepjetc_sf
from BTVNanoCommissioning.helpers.definitions import definitions
from BTVNanoCommissioning.helpers.cTagSFReader import getSF

class NanoProcessor(processor.ProcessorABC):
    # Define histograms
    

    def num(ar):
        return ak.num(ak.fill_none(ar[~ak.is_none(ar)], 0), axis=0)
    def __init__(self):        
        # Define axes
        # Should read axes from NanoAOD config
        dataset_axis = hist.Cat("dataset", "Primary dataset")
        flav_axis = hist.Bin("flav", r"Genflavour",[0,1,4,5,6])

        # Jet
        jet_pt_axis   = hist.Bin("pt",   r"Jet $p_{T}$ [GeV]", 50, 0, 500)
        jet_eta_axis  = hist.Bin("eta",  r"Jet $\eta$", 25, -2.5, 2.5)
        jet_phi_axis  = hist.Bin("phi",  r"Jet $\phi$", 30, -3, 3)
        jet_mass_axis = hist.Bin("mass", r"Jet $m$ [GeV]", 50, 0, 500)

        # Events

        njet_axis  = hist.Bin("njet",  r"N jets",      [0,1,2,3])
        ## Muons
        lpt_axis = hist.Bin("pt",   r"leading Lepton pt", 40, 0, 200)
        leta_axis  = hist.Bin("eta",  r"leading Lepton $\eta$", 25, -2.5, 2.5)
        lphi_axis  = hist.Bin("phi",  r"l eadingLepton $\phi$", 30, -3, 3)
        liso_axis = hist.Bin("pfRelIso03_all", r"leading Muon Rel. Iso", 40,0,4.)
        lpt_ratio_axis = hist.Bin("lptratio", r"leading $\mu$/Jet $p_{T}$ [GeV]",40,0,1)
        slpt_axis   = hist.Bin("pt",   r"subleading Lepton pt", 40, 0, 40)
        sleta_axis  = hist.Bin("eta",  r"subleading Lepton $\eta$", 25, -2.5, 2.5)
        slphi_axis  = hist.Bin("phi",  r"subleading Lepton $\phi$", 30, -3, 3)
        sliso_axis = hist.Bin("pfRelIso03_all", r"subleading Soft Muon Rel. Iso", 40,0,4.)
        slpt_ratio_axis = hist.Bin("slptratio", r"subleading Soft $\mu$/Jet $p_{T}$ [GeV]",40,0,1)

        ## Soft Muons
        softpt_axis = hist.Bin("pt",   r"Soft lepton pt", 25, 0, 25)
        softeta_axis  = hist.Bin("eta",  r"Soft Lepton $\eta$", 25, -2.5, 2.5)
        softphi_axis  = hist.Bin("phi",  r"Soft Lepton $\phi$", 30, -3, 3)
        softiso_axis = hist.Bin("pfRelIso04_all", r"Soft Muon Rel. Iso", 40,0,4.)
        softpt_ratio_axis = hist.Bin("softptratio", r"Soft $\mu$/Jet $p_{T}$ [GeV]",40,0,1)
        
        

        ## Z/W
        zmass_axis = hist.Bin("zmass", r"Z Mass", 25, 50,100)
        zpt_axis = hist.Bin("zpt",r"Z $p_{T}$", 25,0,100)
        zeta_axis  = hist.Bin("zeta",  r"Z $\eta$", 25, -2.5, 2.5)
        zphi_axis  = hist.Bin("zphi",  r"Z $\phi$", 30, -3, 3)
        drmumu_axis = hist.Bin("drmumu", r"$\Delta$R($\mu_{soft}$,$\mu_{hard}$)", 25,0,5)
        
        ## MET
        met_axis = hist.Bin("pt", r"MET $p_{T}$", 50, 0,500)
        metphi_axis = hist.Bin("phi",  r"met $\phi$", 30, -3, 3)

        ## Muon jets
        lmujet_pt_axis   = hist.Bin("pt",   r"Jet $p_{T}$ [GeV]", 50, 0, 500)
        lmujet_eta_axis  = hist.Bin("eta",  r"Jet $\eta$", 25, -2.5, 2.5)
        lmujet_phi_axis  = hist.Bin("phi",  r"Jet $\phi$", 30, -3, 3)
        lmujet_mass_axis = hist.Bin("mass", r"Jet $m$ [GeV]", 50, 0, 500)
        slmujet_pt_axis   = hist.Bin("pt",   r"Jet $p_{T}$ [GeV]", 50, 0, 500)
        slmujet_eta_axis  = hist.Bin("eta",  r"Jet $\eta$", 25, -2.5, 2.5)
        slmujet_phi_axis  = hist.Bin("phi",  r"Jet $\phi$", 30, -3, 3)
        slmujet_mass_axis = hist.Bin("mass", r"Jet $m$ [GeV]", 50, 0, 500)
        dr_lmujetsmu_axis = hist.Bin("dr_lmujetsmu", r"$\Delta$R($\mu_{soft}$,j)",25,0,5)
        dr_slmujetsmu_axis = hist.Bin("dr_slmujetsmu", r"$\Delta$R($\mu_{soft}$,j)",25,0,5)

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
        input_names,manual_ranges,bins = definitions()
        bininfo = dict(zip(input_names,zip(bins,manual_ranges)))
        for d in deepcsv_list:
            binning, ranges = bininfo["Jet_%s"%d]
            if ranges[1] is None : ranges[1] = 0.
            if ranges[0] is None : ranges[0] = -0.5
            deepcsv_axes.append(hist.Bin(d,d,binning,ranges[0],ranges[1]))
        # Define similar axes dynamically
        disc_list = ['btagDeepB', 'btagDeepC', 'btagDeepFlavB', 'btagDeepFlavC','deepcsv_CvL','deepcsv_CvB','deepflav_CvL','deepflav_CvB']
        syst_list = ['','SF','_up','_dn']
        varlist=[]
        btag_axes = []
        for d in disc_list:
            for s in syst_list:
                btag_axes.append(hist.Bin("%s%s" %(d,s), "%s%s" %(d,s), 51, -0.2, 1))  
                varlist.append("%s%s" %(d,s))
                
        _hist_sf_dict={}   
        _hist_deepcsv_dict={
            'pt'  : hist.Hist("Counts", dataset_axis, flav_axis, jet_pt_axis),
            'eta' : hist.Hist("Counts", dataset_axis, flav_axis,jet_eta_axis),
            'phi' : hist.Hist("Counts", dataset_axis, flav_axis,jet_phi_axis),
            'mass': hist.Hist("Counts", dataset_axis, flav_axis,jet_mass_axis)
        }
        for disc, axis in zip(varlist, btag_axes):
            for i in range(3):
                _hist_sf_dict["%s_%d" %(disc,i)] = hist.Hist("Counts", dataset_axis, flav_axis,axis)
        for deepcsv, axises in zip(deepcsv_list, deepcsv_axes):
             _hist_deepcsv_dict["%s" %(deepcsv)] = hist.Hist("Counts", dataset_axis,flav_axis, axises)
        
        _hist_event_dict = {
            'njet' : hist.Hist("Counts", dataset_axis, njet_axis),
            'hl_pt' : hist.Hist("Counts", dataset_axis, lpt_axis),
            'hl_eta' : hist.Hist("Counts", dataset_axis, leta_axis),
            'hl_phi' : hist.Hist("Counts", dataset_axis, lphi_axis),
            'hl_pfRelIso03_all' : hist.Hist("Counts", dataset_axis, liso_axis),
            'sl_pt' : hist.Hist("Counts", dataset_axis, slpt_axis),
            'sl_eta' : hist.Hist("Counts", dataset_axis, sleta_axis),
            'sl_phi' : hist.Hist("Counts", dataset_axis, slphi_axis),
            'sl_pfRelIso03_all' : hist.Hist("Counts", dataset_axis, sliso_axis),
            'hlptratio' :  hist.Hist("Counts", dataset_axis, flav_axis,lpt_ratio_axis),
            'slptratio': hist.Hist("Counts", dataset_axis, flav_axis,slpt_ratio_axis),
            'soft_lpt':hist.Hist("Counts",dataset_axis,softpt_axis),
            'soft_leta':hist.Hist("Counts",dataset_axis,softeta_axis),
            'soft_lphi':hist.Hist("Counts",dataset_axis,softphi_axis),
            'soft_liso' : hist.Hist("Counts", dataset_axis, softiso_axis),
            'softlptratio': hist.Hist("Counts", dataset_axis, flav_axis,softpt_ratio_axis),
            'zmass': hist.Hist("Counts",dataset_axis,zmass_axis),
            'zpt': hist.Hist("Counts",dataset_axis,zpt_axis),
            'zeta': hist.Hist("Counts",dataset_axis,zeta_axis),
            'zphi': hist.Hist("Counts",dataset_axis,zphi_axis),
            'metpt' : hist.Hist("Counts", dataset_axis, met_axis), 
            'metphi' : hist.Hist("Counts", dataset_axis, metphi_axis),
            'lmujet_pt': hist.Hist("Counts",dataset_axis,flav_axis,lmujet_pt_axis),
            'lmujet_eta': hist.Hist("Counts",dataset_axis,flav_axis,lmujet_eta_axis),
            'lmujet_phi': hist.Hist("Counts",dataset_axis,flav_axis,lmujet_phi_axis),
            'lmujet_mass': hist.Hist("Counts",dataset_axis,flav_axis,lmujet_mass_axis),
            'slmujet_pt': hist.Hist("Counts",dataset_axis,flav_axis,slmujet_pt_axis),
            'slmujet_eta': hist.Hist("Counts",dataset_axis,flav_axis,slmujet_eta_axis),
            'slmujet_phi': hist.Hist("Counts",dataset_axis,flav_axis,slmujet_phi_axis),
            'slmujet_mass': hist.Hist("Counts",dataset_axis,flav_axis,slmujet_mass_axis),
            'dr_lmujetsmu':hist.Hist("Counts",dataset_axis,flav_axis,dr_lmujetsmu_axis),
            'dr_slmujetsmu':hist.Hist("Counts",dataset_axis,flav_axis,dr_slmujetsmu_axis),
            }

        self.sf_hists = list(_hist_sf_dict.keys())
        self.event_hists = list(_hist_event_dict.keys())
        self.deepcsv_hists = list(_hist_deepcsv_dict.keys())
        _hist_dict = {**_hist_sf_dict,**_hist_event_dict,**_hist_deepcsv_dict}
        #,**_hist_deepcsv_dict}
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
        if(isRealData): req_lumi=lumiMasks['2017'](events.run, events.luminosityBlock)
        weights = Weights(len(events), storeIndividual=True)
        if not hasattr(events,"btagDeepFlavCvL"): 
            events.Jet['btagDeepFlavCvL'] = np.where(((events.Jet.btagDeepFlavC/(1.-events.Jet.btagDeepFlavB))>0)&(events.Jet.pt>15),(events.Jet.btagDeepFlavC/(1.-events.Jet.btagDeepFlavB)),-1)
            events.Jet['btagDeepFlavCvB'] = np.where(((events.Jet.btagDeepFlavC/(events.Jet.btagDeepFlavC+events.Jet.btagDeepFlavB))>0)&(events.Jet.pt>15),(events.Jet.btagDeepFlavC/(events.Jet.btagDeepFlavC+events.Jet.btagDeepFlavB)),-1)
            events.Jet['btagDeepCvL'] = np.where((events.Jet.btagDeepC>0)&(events.Jet.pt>15),(events.Jet.btagDeepC/(1.-events.Jet.btagDeepB)),-1)
            events.Jet['btagDeepCvB'] = np.where((events.Jet.btagDeepC>0)&(events.Jet.pt>15),(events.Jet.btagDeepC/(events.Jet.btagDeepC+events.Jet.btagDeepB)),-1)
        if not isRealData:
            weights.add('genweight',events.genWeight)
            weights.add('puweight', compiled['2017_pileupweight'](events.Pileup.nPU))
        ##############
        # Trigger level
        triggers = [
        'Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL',
        'Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL',
        ]
        trig_arrs = [events.HLT[_trig] for _trig in triggers]
        req_trig_ele = np.zeros(len(events), dtype='bool')
        for t in trig_arrs:
            req_trig_ele = req_trig_ele | t
        trig_arrsm = [events.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL]
        req_trig_mu = np.zeros(len(events), dtype='bool')
        for t in trig_arrsm: req_trig_mu = req_trig_mu | t

        ############
        # Event level
        ## Muon cuts
        # muon twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
        # iso_muon = events.Muon[(events.Muon.pt > 12) &(abs(events.Muon.eta) < 2.4) & (events.Muon.tightId > .5)&(events.Muon.pfRelIso04_all<=0.15)]
        iso_muon_mu = events.Muon[(events.Muon.pt > 14) &(abs(events.Muon.eta) < 2.4) & (events.Muon.tightId > .5)&(events.Muon.pfRelIso04_all<=0.15)]
        iso_ele_mu = events.Electron[(events.Electron.pt>15)&((abs(events.Electron.eta) < 1.4442)|((abs(events.Electron.eta) < 2.5)&(abs(events.Electron.eta) >1.566)))& (events.Electron.mvaFall17V2Iso_WP80 > .5)]
        iso_muon_ele = events.Muon[(events.Muon.pt > 14) &(abs(events.Muon.eta) < 2.4) & (events.Muon.tightId > .5)&(events.Muon.pfRelIso04_all<=0.15)]
        iso_ele_ele = events.Electron[(events.Electron.pt>27)&((abs(events.Electron.eta) < 1.4442)|((abs(events.Electron.eta) < 2.5)&(abs(events.Electron.eta) >1.566)))& (events.Electron.mvaFall17V2Iso_WP80 > .5)]
        
     
        
        iso_muon_ele=ak.pad_none(iso_muon_ele,1)
        iso_ele_ele=ak.pad_none(iso_ele_ele,1)
        iso_muon_mu =ak.pad_none(iso_muon_mu,1)
        iso_ele_mu =ak.pad_none(iso_ele_mu,1)
        req_ele = (ak.count(iso_muon_ele.pt,axis=1)==1) & (ak.count(iso_ele_ele.pt,axis=1)==1)
        req_mu = (ak.count(iso_muon_mu.pt,axis=1)==1) & (ak.count(iso_ele_mu.pt,axis=1)==1)

        # req_muon =(ak.count(iso_muon.pt, axis=1) == 2) & (iso_muon[:,0].pt>20)
        # dilep_ele = events.Electron[(events.Electron.pt>15)&((abs(events.Electron.eta) < 1.4442)|((abs(events.Electron.eta) < 2.5)&(abs(events.Electron.eta) >1.566)))& (events.Electron.mvaFall17V2Iso_WP80 > .5)]
        # req_dilepveto = (ak.count(dilep_ele.pt,axis=1)==0)

        req_MET = (events.METFixEE2017.pt>40)
        
        ## Jet cuts 
        event_jet =  events.Jet[(events.Jet.pt > 20) & (abs(events.Jet.eta) <= 2.5)&(events.Jet.jetId>=5)& ((events.Jet.pt>50)|(events.Jet.puId>=7))]
        req_jets = (ak.count(event_jet.pt, axis=1) >=2)
        
        ## Soft Muon cuts 
        
        soft_muon = events.Muon[(events.Muon.pt < 25) & (abs(events.Muon.eta) < 2.4) & (events.Muon.tightId > .5)&(events.Muon.pfRelIso04_all>0.2)]
        req_softmu=(ak.count(soft_muon.pt, axis=1)>=1)
        
        
        mu_jet = event_jet[(ak.all(event_jet.metric_table(soft_muon) <= 0.4, axis=2))&((event_jet.muonIdx1!=-1)|(event_jet.muonIdx2!=-1))]
        req_mujet = (ak.count(mu_jet.pt,axis=1)>=1)
       
        #dilepton mass
        
        dilep_mass_ele = iso_muon_ele[:,0]+iso_ele_ele[:,0]
        req_dilepmass_ele = ((dilep_mass_ele.mass>12.)&((dilep_mass_ele.mass<75)| (dilep_mass_ele.mass>105)))
        dilep_mass_mu = iso_muon_mu[:,0]+iso_ele_mu[:,0]
        req_dilepmass_mu = ((dilep_mass_mu.mass>12.)&((dilep_mass_mu.mass<75)| (dilep_mass_mu.mass>105)))

        event_level =  req_lumi & req_MET & req_jets & req_softmu&req_mujet & ((req_trig_ele& req_dilepmass_ele&req_ele)|(req_trig_mu& req_dilepmass_mu & req_mu))
        
        if(len(event_level)>0):event_level = ak.fill_none(event_level,False)        
        # Selected
        selev = events[event_level]    
        
        #########
        
        ## Hard Muon
        shmu =  selev.Muon
        shele = selev.Electron
        'Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL',
        'Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL',
        'Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ',
        'Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ',

        # if selev.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL:
        shmu =  selev.Muon[((selev.Muon.pt > 25) &(abs(selev.Muon.eta) < 2.4) & (selev.Muon.tightId > .5)&(selev.Muon.pfRelIso04_all<=0.15)&(selev.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL==True))|((selev.Muon.pt > 14) &(abs(selev.Muon.eta) < 2.4) & (selev.Muon.tightId > .5)&(selev.Muon.pfRelIso04_all<=0.15)&(selev.HLT.Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL==True))]
        shele = selev.Electron[((selev.Electron.pt>15)&((abs(selev.Electron.eta) < 1.4442)|((abs(selev.Electron.eta) < 2.5)&(abs(selev.Electron.eta) >1.566)))& (selev.Electron.mvaFall17V2Iso_WP80 > .5)&(selev.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL==True))|((selev.Electron.pt>27)&((abs(selev.Electron.eta) < 1.4442)|((abs(selev.Electron.eta) < 2.5)&(abs(selev.Electron.eta) >1.566)))& (selev.Electron.mvaFall17V2Iso_WP80 > .5)&(selev.HLT.Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL==True))]

        ## Soft Muon
        ssmu = selev.Muon[(selev.Muon.pt < 25) & (abs(selev.Muon.eta) < 2.4) & (selev.Muon.tightId > .5)&(selev.Muon.pfRelIso04_all>0.2)]
        softmu0=ssmu[:,0]
        sz=shmu[:,0]+shele[:,0]
        if not isRealData:
                weights.add('lep1sf',np.where(events.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL&~events.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,muSFs(ak.firsts(events.Muon[(events.Muon.pt>25)&(abs(events.Muon.eta) < 2.4)& (events.Muon.tightId > .5)&(events.Muon.pfRelIso04_all<=0.15)])),muSFs(ak.firsts(events.Muon[(events.Muon.pt>14)&(abs(events.Muon.eta) < 2.4)& (events.Muon.tightId > .5)&(events.Muon.pfRelIso04_all<=0.15)]))))
                weights.add('lep2sf',np.where(events.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL&~events.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL,eleSFs(ak.firsts(events.Electron[(events.Electron.pt>15)&((abs(events.Electron.eta) < 1.4442)|((abs(events.Electron.eta) < 2.5)&(abs(events.Electron.eta) >1.566)))& (events.Electron.mvaFall17V2Iso_WP80 > .5)])),eleSFs(ak.firsts(events.Electron[(events.Electron.pt>27)&((abs(events.Electron.eta) < 1.4442)|((abs(events.Electron.eta) < 2.5)&(abs(events.Electron.eta) >1.566)))& (events.Electron.mvaFall17V2Iso_WP80 > .5)]))))
            
        isomu0 = shmu[:,0] 
        isomu1 = shele[:,0]
        # softmu0 = ssmu[:,0]
        # softmu1 = ssmu[:,1]
        ## Jets
        
        sjets = selev.Jet[(selev.Jet.pt > 20) & (abs(selev.Jet.eta) <= 2.5)&(selev.Jet.jetId>=5)& ((selev.Jet.pt>50)|(selev.Jet.puId>=7))]
        
       
        ## Muon Jet 
        smuon_jet = sjets[(ak.all(sjets.metric_table(ssmu) <= 0.4, axis=2))&((sjets.muonIdx1!=-1)|(sjets.muonIdx2!=-1))]
        smuon_jet = smuon_jet[:,0]
        
        njet= ak.count(sjets.pt,axis=1)   
        def flatten(ar): # flatten awkward into a 1d array to hist
            return ak.flatten(ar, axis=None)
        # FIXME - adapt getSF to parse awkard - move `ak.to_numpy` down
        # FIXME - define reused variables aka `sjets[:,1].hadronFlavour` -> `_jet_hf`
        if isRealData :
            genflavor = ak.zeros_like(sjets.pt)
        else:
            par_flav = (sjets.partonFlavour == 0 ) & (sjets.hadronFlavour==0)
            genflavor = sjets.hadronFlavour + 1*par_flav 
            # genweiev=ak.flatten(ak.broadcast_arrays(weights.weight()[event_level],sjets['pt'])[0])
            jetsfs_c_lj = getSF(ak.to_numpy(smuon_jet.hadronFlavour),ak.to_numpy(smuon_jet.btagDeepFlavCvL),ak.to_numpy(smuon_jet.btagDeepFlavCvB),deepjetc_sf)
            jetsfs_c_up_lj = getSF(ak.to_numpy(smuon_jet.hadronFlavour),ak.to_numpy(smuon_jet.btagDeepFlavCvL),ak.to_numpy(smuon_jet.btagDeepFlavCvB),deepjetc_sf,"TotalUncUp")
            jetsfs_c_dn_lj = getSF(ak.to_numpy(smuon_jet.hadronFlavour),ak.to_numpy(smuon_jet.btagDeepFlavCvL),ak.to_numpy(smuon_jet.btagDeepFlavCvB),deepjetc_sf,"TotalUncDown")
            jetsfs_b_lj = deepjetb_sf.eval('central',smuon_jet.hadronFlavour,abs(smuon_jet.eta),smuon_jet.pt,discr=smuon_jet.btagDeepFlavB)
            jetsfs_b_up_lj = deepjetb_sf.eval('up_jes',smuon_jet.hadronFlavour,abs(smuon_jet.eta),smuon_jet.pt,discr=smuon_jet.btagDeepFlavB)
            jetsfs_b_dn_lj = deepjetb_sf.eval('down_jes',smuon_jet.hadronFlavour,abs(smuon_jet.eta),smuon_jet.pt,discr=smuon_jet.btagDeepFlavB)
            csvsfs_c_lj = getSF(ak.to_numpy(smuon_jet.hadronFlavour),ak.to_numpy(smuon_jet.btagDeepFlavCvL),ak.to_numpy(smuon_jet.btagDeepFlavCvB),deepcsvc_sf)
            csvsfs_c_up_lj = getSF(ak.to_numpy(smuon_jet.hadronFlavour),ak.to_numpy(smuon_jet.btagDeepCvL),ak.to_numpy(smuon_jet.btagDeepCvB),deepcsvc_sf,"TotalUncUp")
            csvsfs_c_dn_lj = getSF(ak.to_numpy(smuon_jet.hadronFlavour),ak.to_numpy(smuon_jet.btagDeepCvL),ak.to_numpy(smuon_jet.btagDeepCvB),deepcsvc_sf,"TotalUncDown")
            csvsfs_b_up_lj = deepcsvb_sf.eval('up_jes',smuon_jet.hadronFlavour,abs(smuon_jet.eta),smuon_jet.pt,discr=smuon_jet.btagDeepB)
            csvsfs_b_lj = deepcsvb_sf.eval('central',smuon_jet.hadronFlavour,abs(smuon_jet.eta),smuon_jet.pt,discr=smuon_jet.btagDeepB)
            csvsfs_b_dn_lj = deepcsvb_sf.eval('down_jes',smuon_jet.hadronFlavour,abs(smuon_jet.eta),smuon_jet.pt,discr=smuon_jet.btagDeepB)
            if(all(i>1 for i in njet)):
                jetsfs_c_sj = getSF(ak.to_numpy(sjets[:,1].hadronFlavour),ak.to_numpy(sjets[:,1].btagDeepFlavCvL),ak.to_numpy(sjets[:,1].btagDeepFlavCvB),deepjetc_sf)
                csvsfs_c_sj = getSF(ak.to_numpy(sjets[:,1].hadronFlavour),ak.to_numpy(sjets[:,1].btagDeepC/(1.-sjets[:,1].btagDeepB)),ak.to_numpy(sjets[:,1].btagDeepC/(sjets[:,1].btagDeepC+sjets[:,1].btagDeepB)),deepcsvc_sf)
                jetsfs_c_up_sj = getSF(ak.to_numpy(sjets[:,1].hadronFlavour),ak.to_numpy(sjets[:,1].btagDeepFlavCvL),ak.to_numpy(sjets[:,1].btagDeepFlavCvB),deepjetc_sf,"TotalUncUp")
                csvsfs_c_up_sj = getSF(ak.to_numpy(sjets[:,1].hadronFlavour),ak.to_numpy(sjets[:,1].btagDeepC/(1.-sjets[:,1].btagDeepB)),ak.to_numpy(sjets[:,1].btagDeepC/(sjets[:,1].btagDeepC+sjets[:,1].btagDeepB)),deepcsvc_sf,"TotalUncUp")
                jetsfs_c_dn_sj = getSF(ak.to_numpy(sjets[:,1].hadronFlavour),ak.to_numpy(sjets[:,1].btagDeepFlavCvL),ak.to_numpy(sjets[:,1].btagDeepFlavCvB),deepjetc_sf,"TotalUncDown")
                csvsfs_c_dn_sj = getSF(ak.to_numpy(sjets[:,1].hadronFlavour),ak.to_numpy(sjets[:,1].btagDeepC/(1.-sjets[:,1].btagDeepB)),ak.to_numpy(sjets[:,1].btagDeepC/(sjets[:,1].btagDeepC+sjets[:,1].btagDeepB)),deepcsvc_sf,"TotalUncDown")
                jetsfs_b_sj = deepjetb_sf.eval('central',sjets[:,1].hadronFlavour,abs(sjets[:,1].eta),sjets[:,1].pt,discr=sjets[:,1].btagDeepFlavB)
                csvsfs_b_sj = deepcsvb_sf.eval('central',sjets[:,1].hadronFlavour,abs(sjets[:,1].eta),sjets[:,1].pt,discr=sjets[:,1].btagDeepB)
                jetsfs_b_up_sj = deepjetb_sf.eval('up_jes',sjets[:,1].hadronFlavour,abs(sjets[:,1].eta),sjets[:,1].pt,discr=sjets[:,1].btagDeepFlavB)
                csvsfs_b_up_sj = deepcsvb_sf.eval('up_jes',sjets[:,1].hadronFlavour,abs(sjets[:,1].eta),sjets[:,1].pt,discr=sjets[:,1].btagDeepB)
                jetsfs_b_dn_sj = deepjetb_sf.eval('down_jes',sjets[:,1].hadronFlavour,abs(sjets[:,1].eta),sjets[:,1].pt,discr=sjets[:,1].btagDeepFlavB)
                csvsfs_b_dn_sj = deepcsvb_sf.eval('down_jes',sjets[:,1].hadronFlavour,abs(sjets[:,1].eta),sjets[:,1].pt,discr=sjets[:,1].btagDeepB)
            
            
        for histname, h in output.items():
            if 'sumw' in histname:print('sumw')
            # elif histname in self.deepcsv_hists:
            #     fields = {l: ak.flatten(sjets[histname]) for l in h.fields if l in dir(sjets)}
            #     if(isRealData):
            #         h.fill(dataset=dataset,flav=5, **fields)
            #     else:                    
            #         genweiev=ak.flatten(ak.broadcast_arrays(weights.weight()[event_level],sjets['pt'])[0])
            #         h.fill(dataset=dataset,flav=ak.flatten(sjets.hadronFlavour), **fields,weight=genweiev)
            elif 'hl_' in histname:
                fields = {l: isomu0[l] for l in h.fields if l in dir(isomu0)}
                if(isRealData): h.fill(dataset=dataset,**fields)
                else: h.fill(dataset=dataset,**fields,weight=weights.weight()[event_level])
            elif 'sl_' in histname:
                fields = {l: isomu1[l] for l in h.fields if l in dir(isomu1)}
                if(isRealData): h.fill(dataset=dataset,**fields)
                else: h.fill(dataset=dataset,**fields,weight=weights.weight()[event_level])
            elif 'soft_l' in histname:
                fields = {l: softmu0[l] for l in h.fields if l in dir(softmu0)}
                if(isRealData): h.fill(dataset=dataset,**fields)
                else: h.fill(dataset=dataset,**fields,weight=weights.weight()[event_level])
            # elif 'soft_sl' in histname:
            #     fields = {l: softmu1[l] for l in h.fields if l in dir(softmu1)}
            #     if(isRealData): h.fill(dataset=dataset,**fields)
            #     else: h.fill(dataset=dataset,**fields,weight=weights.weight()[event_level])
            elif 'met' in histname: 
                fields = {l: selev.METFixEE2017[l] for l in h.fields if l in dir(selev.METFixEE2017)}
                if(isRealData): h.fill(dataset=dataset,**fields)
                else: h.fill(dataset=dataset,**fields,weight=weights.weight()[event_level])
            elif 'lmujet_' in histname:
                fields = {l: smuon_jet[l] for l in h.fields if l in dir(smuon_jet)}
                if isRealData: h.fill(dataset=dataset,flav=5,  **fields)
                else :h.fill(dataset=dataset,flav=5,  **fields,weight = weights.weight()[event_level])
            # elif 'slmujet_' in histname:
            #     fields = {l: sjets[:,1][l] for l in h.fields if l in dir(sjets[:,1])}
            #     if isRealData: h.fill(dataset=dataset,flav=5,  **fields)
            #     else :h.fill(dataset=dataset,flav=5,  **fields,weight = weights.weight()[event_level])
         
        if not isRealData:
            ###Fill no SFs
            smpu= (smuon_jet.partonFlavour == 0 ) & (smuon_jet.hadronFlavour==0)
            smflav= 1*smpu+smuon_jet.hadronFlavour
            output['zmass'].fill(dataset=dataset,zmass=sz.mass,weight=weights.weight()[event_level])
            output['zpt'].fill(dataset=dataset,zpt=sz.pt,weight=weights.weight()[event_level])
            output['zeta'].fill(dataset=dataset,zeta=sz.eta,weight=weights.weight()[event_level])
            output['zphi'].fill(dataset=dataset,zphi=sz.phi,weight=weights.weight()[event_level])
            output['hlptratio'].fill(dataset=dataset,flav=sjets[:,0].hadronFlavour,lptratio=isomu0.pt/sjets[:,0].pt,weight=weights.weight()[event_level])
            # output['slptratio'].fill(dataset=dataset,flav=smflav[:,1],slptratio=isomu1.pt/sjets[:,1].pt,weight=weights.weight()[event_level])
            output['softlptratio'].fill(dataset=dataset,flav=smflav,softptratio=softmu0.pt/smuon_jet.pt,weight=weights.weight()[event_level])
            # output['softslptratio'].fill(dataset=dataset,flav=smflav[:,1],softptratio=softmu0.pt/sjets[:,1].pt,weight=weights.weight()[event_level])
            
            output['dr_lmujetsmu'].fill(dataset=dataset,flav=smflav,dr_lmujetsmu=smuon_jet.delta_r(softmu0),weight=weights.weight()[event_level])
            # output['dr_slmujetsmu'].fill(dataset=dataset,flav=smflav[:,1],dr_slmujetsmu=sjets[:,1].delta_r(softmu0),weight=weights.weight()[event_level])

            ## discri
            output['btagDeepFlavB_0'].fill(dataset=dataset,flav=smflav, btagDeepFlavB=smuon_jet.btagDeepFlavB,weight=weights.weight()[event_level])
            output['btagDeepFlavC_0'].fill(dataset=dataset,flav=smflav, btagDeepFlavC=smuon_jet.btagDeepFlavC,weight=weights.weight()[event_level])
            output['btagDeepB_0'].fill(dataset=dataset,flav=smflav, btagDeepB=smuon_jet.btagDeepB,weight=weights.weight()[event_level])
            output['btagDeepC_0'].fill(dataset=dataset,flav=smflav, btagDeepC=smuon_jet.btagDeepC,weight=weights.weight()[event_level])
            output['deepcsv_CvB_0'].fill(dataset=dataset,flav=smflav, deepcsv_CvB=np.where(smuon_jet.btagDeepCvB<0,-0.2,smuon_jet.btagDeepCvB),weight=weights.weight()[event_level])
            output['deepcsv_CvL_0'].fill(dataset=dataset,flav=smflav, deepcsv_CvL=np.where(smuon_jet.btagDeepCvL<0,-0.2,smuon_jet.btagDeepCvL),weight=weights.weight()[event_level])
            output['deepflav_CvB_0'].fill(dataset=dataset,flav=smflav, deepflav_CvB=np.where(smuon_jet.btagDeepFlavCvB<0,-0.2,smuon_jet.btagDeepFlavCvB),weight=weights.weight()[event_level])
            output['deepflav_CvL_0'].fill(dataset=dataset,flav=smflav, deepflav_CvL=np.where(smuon_jet.btagDeepFlavCvL<0,-0.2,smuon_jet.btagDeepFlavCvL),weight=weights.weight()[event_level])
            output['btagDeepFlavBSF_0'].fill(dataset=dataset,flav=smflav, btagDeepFlavBSF=smuon_jet.btagDeepFlavB,weight=(weights.weight()[event_level]*jetsfs_b_lj))      
            output['btagDeepFlavCSF_0'].fill(dataset=dataset,flav=smflav, btagDeepFlavCSF=smuon_jet.btagDeepFlavC,weight=weights.weight()[event_level]*jetsfs_c_lj)
            output['btagDeepBSF_0'].fill(dataset=dataset,flav=smflav, btagDeepBSF=smuon_jet.btagDeepB,weight=weights.weight()[event_level]*csvsfs_b_lj)
            output['btagDeepCSF_0'].fill(dataset=dataset,flav=smflav, btagDeepCSF=smuon_jet.btagDeepC,weight=weights.weight()[event_level]*csvsfs_c_lj)
            output['deepcsv_CvBSF_0'].fill(dataset=dataset,flav=smflav, deepcsv_CvBSF=np.where(smuon_jet.btagDeepCvB<0,-0.2,smuon_jet.btagDeepCvB),weight=weights.weight()[event_level]*csvsfs_c_lj)
            output['deepcsv_CvLSF_0'].fill(dataset=dataset,flav=smflav, deepcsv_CvLSF=np.where(smuon_jet.btagDeepCvL<0,-0.2,smuon_jet.btagDeepCvL),weight=weights.weight()[event_level]*csvsfs_c_lj)
            output['deepflav_CvBSF_0'].fill(dataset=dataset,flav=smflav, deepflav_CvBSF=np.where(smuon_jet.btagDeepFlavCvB<0,-0.2,smuon_jet.btagDeepFlavCvB),weight=weights.weight()[event_level]*jetsfs_c_lj)
            output['deepflav_CvLSF_0'].fill(dataset=dataset,flav=smflav, deepflav_CvLSF=np.where(smuon_jet.btagDeepFlavCvL<0,-0.2,smuon_jet.btagDeepFlavCvL),weight=weights.weight()[event_level]*jetsfs_c_lj)
            output['btagDeepFlavB_up_0'].fill(dataset=dataset,flav=smflav, btagDeepFlavB_up=smuon_jet.btagDeepFlavB,weight=weights.weight()[event_level]*jetsfs_b_up_lj)
            output['btagDeepFlavC_up_0'].fill(dataset=dataset,flav=smflav, btagDeepFlavC_up=smuon_jet.btagDeepFlavC,weight=weights.weight()[event_level]*jetsfs_c_up_lj)
            output['btagDeepB_up_0'].fill(dataset=dataset,flav=smflav, btagDeepB_up=smuon_jet.btagDeepB,weight=weights.weight()[event_level]*csvsfs_b_up_lj)
            output['btagDeepC_up_0'].fill(dataset=dataset,flav=smflav, btagDeepC_up=smuon_jet.btagDeepC,weight=weights.weight()[event_level]*csvsfs_c_up_lj)
            output['deepcsv_CvB_up_0'].fill(dataset=dataset,flav=smflav, deepcsv_CvB_up=np.where(smuon_jet.btagDeepCvB<0,-0.2,smuon_jet.btagDeepCvB),weight=weights.weight()[event_level]*csvsfs_c_up_lj)
            output['deepcsv_CvL_up_0'].fill(dataset=dataset,flav=smflav, deepcsv_CvL_up=np.where(smuon_jet.btagDeepCvL<0,-0.2,smuon_jet.btagDeepCvL),weight=weights.weight()[event_level]*csvsfs_c_up_lj)
            output['deepflav_CvB_up_0'].fill(dataset=dataset,flav=smflav, deepflav_CvB_up=np.where(smuon_jet.btagDeepFlavCvB<0,-0.2,smuon_jet.btagDeepFlavCvB),weight=weights.weight()[event_level]*jetsfs_c_up_lj)
            output['deepflav_CvL_up_0'].fill(dataset=dataset,flav=smflav, deepflav_CvL_up=np.where(smuon_jet.btagDeepFlavCvL<0,-0.2,smuon_jet.btagDeepFlavCvL),weight=weights.weight()[event_level]*jetsfs_c_up_lj)
            output['btagDeepFlavB_dn_0'].fill(dataset=dataset,flav=smflav, btagDeepFlavB_dn=smuon_jet.btagDeepFlavB,weight=weights.weight()[event_level]*jetsfs_b_dn_lj)
            output['btagDeepFlavC_dn_0'].fill(dataset=dataset,flav=smflav, btagDeepFlavC_dn=smuon_jet.btagDeepFlavC,weight=weights.weight()[event_level]*jetsfs_c_dn_lj)
            output['btagDeepB_dn_0'].fill(dataset=dataset,flav=smflav, btagDeepB_dn=smuon_jet.btagDeepB,weight=weights.weight()[event_level]*csvsfs_b_dn_lj)
            output['btagDeepC_dn_0'].fill(dataset=dataset,flav=smflav, btagDeepC_dn=smuon_jet.btagDeepC,weight=weights.weight()[event_level]*csvsfs_c_dn_lj)
            output['deepcsv_CvB_dn_0'].fill(dataset=dataset,flav=smflav, deepcsv_CvB_dn=np.where(smuon_jet.btagDeepCvB<0,-0.2,smuon_jet.btagDeepCvB),weight=weights.weight()[event_level]*csvsfs_c_dn_lj)
            output['deepcsv_CvL_dn_0'].fill(dataset=dataset,flav=smflav, deepcsv_CvL_dn=np.where(smuon_jet.btagDeepCvL<0,-0.2,smuon_jet.btagDeepCvL),weight=weights.weight()[event_level]*csvsfs_c_dn_lj)
            output['deepflav_CvB_dn_0'].fill(dataset=dataset,flav=smflav, deepflav_CvB_dn=np.where(smuon_jet.btagDeepFlavCvB<0,-0.2,smuon_jet.btagDeepFlavCvB),weight=weights.weight()[event_level]*jetsfs_c_dn_lj)
            output['deepflav_CvL_dn_0'].fill(dataset=dataset,flav=smflav, deepflav_CvL_dn=np.where(smuon_jet.btagDeepFlavCvL<0,-0.2,smuon_jet.btagDeepFlavCvL),weight=weights.weight()[event_level]*jetsfs_c_dn_lj)
            # print(smuon_jet,ak.count(smuon_jet.pt,axis=1))
            if(all(i>1 for i in njet)):
                output['btagDeepFlavB_1'].fill(dataset=dataset,flav=genflavor[:,1], btagDeepFlavB=sjets[:,1].btagDeepFlavB,weight=weights.weight()[event_level])
                output['btagDeepFlavC_1'].fill(dataset=dataset,flav=genflavor[:,1], btagDeepFlavC=sjets[:,1].btagDeepFlavC,weight=weights.weight()[event_level])
                output['btagDeepB_1'].fill(dataset=dataset,flav=genflavor[:,1], btagDeepB=sjets[:,1].btagDeepB,weight=weights.weight()[event_level])
                output['btagDeepC_1'].fill(dataset=dataset,flav=genflavor[:,1], btagDeepC=sjets[:,1].btagDeepC,weight=weights.weight()[event_level])
                output['deepcsv_CvB_1'].fill(dataset=dataset,flav=genflavor[:,1], deepcsv_CvB=np.where(sjets[:,1].btagDeepCvB<0,-0.2,sjets[:,1].btagDeepCvB),weight=weights.weight()[event_level])
                output['deepcsv_CvL_1'].fill(dataset=dataset,flav=genflavor[:,1], deepcsv_CvL=np.where(sjets[:,1].btagDeepCvL<0,-0.2,sjets[:,1].btagDeepCvL),weight=weights.weight()[event_level])
                output['deepflav_CvB_1'].fill(dataset=dataset,flav=genflavor[:,1], deepflav_CvB=np.where(sjets[:,1].btagDeepFlavCvB<0,-0.2,sjets[:,1].btagDeepFlavCvB),weight=weights.weight()[event_level])
                output['deepflav_CvL_1'].fill(dataset=dataset,flav=genflavor[:,1], deepflav_CvL=np.where(sjets[:,1].btagDeepFlavCvL<0,-0.2,sjets[:,1].btagDeepFlavCvL),weight=weights.weight()[event_level])
                output['btagDeepFlavBSF_1'].fill(dataset=dataset,flav=genflavor[:,1], btagDeepFlavBSF=sjets[:,1].btagDeepFlavB,weight=weights.weight()[event_level]*jetsfs_b_sj)
                output['btagDeepFlavCSF_1'].fill(dataset=dataset,flav=genflavor[:,1], btagDeepFlavCSF=sjets[:,1].btagDeepFlavC,weight=weights.weight()[event_level]*jetsfs_c_sj)
                output['btagDeepBSF_1'].fill(dataset=dataset,flav=genflavor[:,1], btagDeepBSF=sjets[:,1].btagDeepB,weight=weights.weight()[event_level]*csvsfs_b_sj)
                output['btagDeepCSF_1'].fill(dataset=dataset,flav=genflavor[:,1], btagDeepCSF=sjets[:,1].btagDeepC,weight=weights.weight()[event_level]*csvsfs_c_sj)
                output['deepcsv_CvBSF_1'].fill(dataset=dataset,flav=genflavor[:,1], deepcsv_CvBSF=np.where(sjets[:,1].btagDeepCvB<0,-0.2,sjets[:,1].btagDeepCvB),weight=weights.weight()[event_level]*csvsfs_c_sj)
                output['deepcsv_CvLSF_1'].fill(dataset=dataset,flav=genflavor[:,1], deepcsv_CvLSF=np.where(sjets[:,1].btagDeepCvL<0,-0.2,sjets[:,1].btagDeepCvL),weight=weights.weight()[event_level]*csvsfs_c_sj)
                output['deepflav_CvBSF_1'].fill(dataset=dataset,flav=genflavor[:,1], deepflav_CvBSF=np.where(sjets[:,1].btagDeepFlavCvB<0,-0.2,sjets[:,1].btagDeepFlavCvB),weight=weights.weight()[event_level]*jetsfs_c_sj)
                output['deepflav_CvLSF_1'].fill(dataset=dataset,flav=genflavor[:,1], deepflav_CvLSF=np.where(sjets[:,1].btagDeepFlavCvL<0,-0.2,sjets[:,1].btagDeepFlavCvL),weight=weights.weight()[event_level]*jetsfs_c_sj)
                output['btagDeepFlavB_up_1'].fill(dataset=dataset,flav=genflavor[:,1], btagDeepFlavB_up=sjets[:,1].btagDeepFlavB,weight=weights.weight()[event_level]*jetsfs_b_up_sj)
                output['btagDeepFlavC_up_1'].fill(dataset=dataset,flav=genflavor[:,1], btagDeepFlavC_up=sjets[:,1].btagDeepFlavC,weight=weights.weight()[event_level]*jetsfs_c_up_sj)
                output['btagDeepB_up_1'].fill(dataset=dataset,flav=genflavor[:,1], btagDeepB_up=sjets[:,1].btagDeepB,weight=weights.weight()[event_level]*csvsfs_b_up_sj)
                output['btagDeepC_up_1'].fill(dataset=dataset,flav=genflavor[:,1], btagDeepC_up=sjets[:,1].btagDeepC,weight=weights.weight()[event_level]*csvsfs_c_up_sj)
                output['deepcsv_CvB_up_1'].fill(dataset=dataset,flav=genflavor[:,1], deepcsv_CvB_up=np.where(sjets[:,1].btagDeepCvB<0,-0.2,sjets[:,1].btagDeepCvB),weight=weights.weight()[event_level]*csvsfs_c_up_sj)
                output['deepcsv_CvL_up_1'].fill(dataset=dataset,flav=genflavor[:,1], deepcsv_CvL_up=np.where(sjets[:,1].btagDeepCvL<0,-0.2,sjets[:,1].btagDeepCvL),weight=weights.weight()[event_level]*csvsfs_c_up_sj)
                output['deepflav_CvB_up_1'].fill(dataset=dataset,flav=genflavor[:,1], deepflav_CvB_up=np.where(sjets[:,1].btagDeepFlavCvB<0,-0.2,sjets[:,1].btagDeepFlavCvB),weight=weights.weight()[event_level]*jetsfs_c_up_sj)
                output['deepflav_CvL_up_1'].fill(dataset=dataset,flav=genflavor[:,1], deepflav_CvL_up=np.where(sjets[:,1].btagDeepFlavCvL<0,-0.2,sjets[:,1].btagDeepFlavCvL),weight=weights.weight()[event_level]*jetsfs_c_up_sj)
                output['btagDeepFlavB_dn_1'].fill(dataset=dataset,flav=genflavor[:,1], btagDeepFlavB_dn=sjets[:,1].btagDeepFlavB,weight=weights.weight()[event_level]*jetsfs_b_dn_sj)
                output['btagDeepFlavC_dn_1'].fill(dataset=dataset,flav=genflavor[:,1], btagDeepFlavC_dn=sjets[:,1].btagDeepFlavC,weight=weights.weight()[event_level]*jetsfs_c_dn_sj)
                output['btagDeepB_dn_1'].fill(dataset=dataset,flav=genflavor[:,1], btagDeepB_dn=sjets[:,1].btagDeepB,weight=weights.weight()[event_level]*csvsfs_b_dn_sj)
                output['btagDeepC_dn_1'].fill(dataset=dataset,flav=genflavor[:,1], btagDeepC_dn=sjets[:,1].btagDeepC,weight=weights.weight()[event_level]*csvsfs_c_dn_sj)
                output['deepcsv_CvB_dn_1'].fill(dataset=dataset,flav=genflavor[:,1], deepcsv_CvB_dn=np.where(sjets[:,1].btagDeepCvB<0,-0.2,sjets[:,1].btagDeepCvB),weight=weights.weight()[event_level]*csvsfs_c_dn_sj)
                output['deepcsv_CvL_dn_1'].fill(dataset=dataset,flav=genflavor[:,1], deepcsv_CvL_dn=np.where(sjets[:,1].btagDeepCvL<0,-0.2,sjets[:,1].btagDeepCvL),weight=weights.weight()[event_level]*csvsfs_c_dn_sj)
                output['deepflav_CvB_dn_1'].fill(dataset=dataset,flav=genflavor[:,1], deepflav_CvB_dn=np.where(sjets[:,1].btagDeepFlavCvB<0,-0.2,sjets[:,1].btagDeepFlavCvB),weight=weights.weight()[event_level]*jetsfs_c_dn_sj)
                output['deepflav_CvL_dn_1'].fill(dataset=dataset,flav=genflavor[:,1], deepflav_CvL_dn=np.where(sjets[:,1].btagDeepFlavCvL<0,-0.2,sjets[:,1].btagDeepFlavCvL),weight=weights.weight()[event_level]*jetsfs_c_dn_sj)
        else:
            ###Fill no SFs
            output['zmass'].fill(dataset=dataset,zmass=sz.mass)
            output['zpt'].fill(dataset=dataset,zpt=sz.pt,weight=weights.weight()[event_level])
            output['zeta'].fill(dataset=dataset,zeta=sz.eta,weight=weights.weight()[event_level])
            output['zphi'].fill(dataset=dataset,zphi=sz.phi,weight=weights.weight()[event_level])
            output['hlptratio'].fill(dataset=dataset,flav=5,lptratio=isomu0.pt/smuon_jet.pt)
            # output['slptratio'].fill(dataset=dataset,flav=5,slptratio=isomu1.pt/sjets[:,1].pt)
            output['softlptratio'].fill(dataset=dataset,flav=5,softptratio=ssmu[:,0].pt/smuon_jet.pt)
            # output['softslptratio'].fill(dataset=dataset,flav=5,softptratio=ssmu[:,1].pt/sjets[:,1].pt)
            
            output['dr_lmujetsmu'].fill(dataset=dataset,flav=5,dr_lmujetsmu=smuon_jet.delta_r(ssmu[:,0]))
            # output['dr_slmujetsmu'].fill(dataset=dataset,flav=5,dr_slmujetsmu=sjets[:,1].delta_r(ssmu[:,1]))
                
                ## discri
            output['btagDeepFlavB_0'].fill(dataset=dataset,flav=5, btagDeepFlavB=smuon_jet.btagDeepFlavB)
            output['btagDeepFlavC_0'].fill(dataset=dataset,flav=5, btagDeepFlavC=smuon_jet.btagDeepFlavC)
            output['btagDeepB_0'].fill(dataset=dataset,flav=5, btagDeepB=smuon_jet.btagDeepB)
            output['btagDeepC_0'].fill(dataset=dataset,flav=5, btagDeepC=smuon_jet.btagDeepC)
            output['deepcsv_CvB_0'].fill(dataset=dataset,flav=5, deepcsv_CvB=np.where(smuon_jet.btagDeepCvB<0,-0.2,smuon_jet.btagDeepCvB))
            output['deepcsv_CvL_0'].fill(dataset=dataset,flav=5, deepcsv_CvL=np.where(smuon_jet.btagDeepCvL<0,-0.2,smuon_jet.btagDeepCvL))
            output['deepflav_CvB_0'].fill(dataset=dataset,flav=5, deepflav_CvB=np.where(smuon_jet.btagDeepFlavCvB<0,-0.2,smuon_jet.btagDeepFlavCvB))
            output['deepflav_CvL_0'].fill(dataset=dataset,flav=5, deepflav_CvL=np.where(smuon_jet.btagDeepFlavCvL<0,-0.2,smuon_jet.btagDeepFlavCvL))
            output['btagDeepFlavBSF_0'].fill(dataset=dataset,flav=5, btagDeepFlavBSF=smuon_jet.btagDeepFlavB)
            output['btagDeepFlavCSF_0'].fill(dataset=dataset,flav=5, btagDeepFlavCSF=smuon_jet.btagDeepFlavC)
            output['btagDeepBSF_0'].fill(dataset=dataset,flav=5, btagDeepBSF=smuon_jet.btagDeepB)
            output['btagDeepCSF_0'].fill(dataset=dataset,flav=5, btagDeepCSF=smuon_jet.btagDeepC)
            output['deepcsv_CvBSF_0'].fill(dataset=dataset,flav=5, deepcsv_CvBSF=np.where(smuon_jet.btagDeepCvB<0,-0.2,smuon_jet.btagDeepCvB))
            output['deepcsv_CvLSF_0'].fill(dataset=dataset,flav=5, deepcsv_CvLSF=np.where(smuon_jet.btagDeepCvL<0,-0.2,smuon_jet.btagDeepCvL))
            output['deepflav_CvBSF_0'].fill(dataset=dataset,flav=5, deepflav_CvBSF=np.where(smuon_jet.btagDeepFlavCvB<0,-0.2,smuon_jet.btagDeepFlavCvB))
            output['deepflav_CvLSF_0'].fill(dataset=dataset,flav=5, deepflav_CvLSF=np.where(smuon_jet.btagDeepFlavCvL<0,-0.2,smuon_jet.btagDeepFlavCvL))
            output['btagDeepFlavB_up_0'].fill(dataset=dataset,flav=5, btagDeepFlavB_up=smuon_jet.btagDeepFlavB)
            output['btagDeepFlavC_up_0'].fill(dataset=dataset,flav=5, btagDeepFlavC_up=smuon_jet.btagDeepFlavC)
            output['btagDeepB_up_0'].fill(dataset=dataset,flav=5, btagDeepB_up=smuon_jet.btagDeepB)
            output['btagDeepC_up_0'].fill(dataset=dataset,flav=5, btagDeepC_up=smuon_jet.btagDeepC)
            output['deepcsv_CvB_up_0'].fill(dataset=dataset,flav=5, deepcsv_CvB_up=np.where(smuon_jet.btagDeepCvB<0,-0.2,smuon_jet.btagDeepCvB))
            output['deepcsv_CvL_up_0'].fill(dataset=dataset,flav=5, deepcsv_CvL_up=np.where(smuon_jet.btagDeepCvL<0,-0.2,smuon_jet.btagDeepCvL))
            output['deepflav_CvB_up_0'].fill(dataset=dataset,flav=5, deepflav_CvB_up=np.where(smuon_jet.btagDeepFlavCvB<0,-0.2,smuon_jet.btagDeepFlavCvB))
            output['deepflav_CvL_up_0'].fill(dataset=dataset,flav=5, deepflav_CvL_up=np.where(smuon_jet.btagDeepFlavCvL<0,-0.2,smuon_jet.btagDeepFlavCvL))
            output['btagDeepFlavB_dn_0'].fill(dataset=dataset,flav=5, btagDeepFlavB_dn=smuon_jet.btagDeepFlavB)
            output['btagDeepFlavC_dn_0'].fill(dataset=dataset,flav=5, btagDeepFlavC_dn=smuon_jet.btagDeepFlavC)
            output['btagDeepB_dn_0'].fill(dataset=dataset,flav=5, btagDeepB_dn=smuon_jet.btagDeepB)
            output['btagDeepC_dn_0'].fill(dataset=dataset,flav=5, btagDeepC_dn=smuon_jet.btagDeepC)
            output['deepcsv_CvB_dn_0'].fill(dataset=dataset,flav=5, deepcsv_CvB_dn=np.where(smuon_jet.btagDeepCvB<0,-0.2,smuon_jet.btagDeepCvB))
            output['deepcsv_CvL_dn_0'].fill(dataset=dataset,flav=5, deepcsv_CvL_dn=np.where(smuon_jet.btagDeepCvL<0,-0.2,smuon_jet.btagDeepCvL))
            output['deepflav_CvB_dn_0'].fill(dataset=dataset,flav=5, deepflav_CvB_dn=np.where(smuon_jet.btagDeepFlavCvB<0,-0.2,smuon_jet.btagDeepFlavCvB))
            output['deepflav_CvL_dn_0'].fill(dataset=dataset,flav=5, deepflav_CvL_dn=np.where(smuon_jet.btagDeepFlavCvL<0,-0.2,smuon_jet.btagDeepFlavCvL))
            if(all(i>1 for i in njet)):
                output['btagDeepFlavB_1'].fill(dataset=dataset,flav=5, btagDeepFlavB=sjets[:,1].btagDeepFlavB)
                output['btagDeepFlavC_1'].fill(dataset=dataset,flav=5, btagDeepFlavC=sjets[:,1].btagDeepFlavC)
                output['btagDeepB_1'].fill(dataset=dataset,flav=5, btagDeepB=sjets[:,1].btagDeepB)
                output['btagDeepC_1'].fill(dataset=dataset,flav=5, btagDeepC=sjets[:,1].btagDeepC)
                output['deepcsv_CvB_1'].fill(dataset=dataset,flav=5, deepcsv_CvB=np.where(sjets[:,1].btagDeepCvB<0,-0.2,sjets[:,1].btagDeepCvB))
                output['deepcsv_CvL_1'].fill(dataset=dataset,flav=5, deepcsv_CvL=np.where(sjets[:,1].btagDeepCvL<0,-0.2,sjets[:,1].btagDeepCvL))
                output['deepflav_CvB_1'].fill(dataset=dataset,flav=5, deepflav_CvB=np.where(sjets[:,1].btagDeepFlavCvB<0,-0.2,sjets[:,1].btagDeepFlavCvB))
                output['deepflav_CvL_1'].fill(dataset=dataset,flav=5, deepflav_CvL=np.where(sjets[:,1].btagDeepFlavCvL<0,-0.2,sjets[:,1].btagDeepFlavCvL))
                output['btagDeepFlavBSF_1'].fill(dataset=dataset,flav=5, btagDeepFlavBSF=sjets[:,1].btagDeepFlavB)
                output['btagDeepFlavCSF_1'].fill(dataset=dataset,flav=5, btagDeepFlavCSF=sjets[:,1].btagDeepFlavC)
                output['btagDeepBSF_1'].fill(dataset=dataset,flav=5, btagDeepBSF=sjets[:,1].btagDeepB)
                output['btagDeepCSF_1'].fill(dataset=dataset,flav=5, btagDeepCSF=sjets[:,1].btagDeepC)
                output['deepcsv_CvBSF_1'].fill(dataset=dataset,flav=5, deepcsv_CvBSF=np.where(sjets[:,1].btagDeepCvB<0,-0.2,sjets[:,1].btagDeepCvB))
                output['deepcsv_CvLSF_1'].fill(dataset=dataset,flav=5, deepcsv_CvLSF=np.where(sjets[:,1].btagDeepCvL<0,-0.2,sjets[:,1].btagDeepCvL))
                output['deepflav_CvBSF_1'].fill(dataset=dataset,flav=5, deepflav_CvBSF=np.where(sjets[:,1].btagDeepFlavCvB<0,-0.2,sjets[:,1].btagDeepFlavCvB))
                output['deepflav_CvLSF_1'].fill(dataset=dataset,flav=5, deepflav_CvLSF=np.where(sjets[:,1].btagDeepFlavCvL<0,-0.2,sjets[:,1].btagDeepFlavCvL))
                output['btagDeepFlavB_up_1'].fill(dataset=dataset,flav=5, btagDeepFlavB_up=sjets[:,1].btagDeepFlavB)
                output['btagDeepFlavC_up_1'].fill(dataset=dataset,flav=5, btagDeepFlavC_up=sjets[:,1].btagDeepFlavC)
                output['btagDeepB_up_1'].fill(dataset=dataset,flav=5, btagDeepB_up=sjets[:,1].btagDeepB)
                output['btagDeepC_up_1'].fill(dataset=dataset,flav=5, btagDeepC_up=sjets[:,1].btagDeepC)
                output['deepcsv_CvB_up_1'].fill(dataset=dataset,flav=5, deepcsv_CvB_up=np.where(sjets[:,1].btagDeepCvB<0,-0.2,sjets[:,1].btagDeepCvB))
                output['deepcsv_CvL_up_1'].fill(dataset=dataset,flav=5, deepcsv_CvL_up=np.where(sjets[:,1].btagDeepCvL<0,-0.2,sjets[:,1].btagDeepCvL))
                output['deepflav_CvB_up_1'].fill(dataset=dataset,flav=5, deepflav_CvB_up=np.where(sjets[:,1].btagDeepFlavCvB<0,-0.2,sjets[:,1].btagDeepFlavCvB))
                output['deepflav_CvL_up_1'].fill(dataset=dataset,flav=5, deepflav_CvL_up=np.where(sjets[:,1].btagDeepFlavCvL<0,-0.2,sjets[:,1].btagDeepFlavCvL))
                output['btagDeepFlavB_dn_1'].fill(dataset=dataset,flav=5, btagDeepFlavB_dn=sjets[:,1].btagDeepFlavB)
                output['btagDeepFlavC_dn_1'].fill(dataset=dataset,flav=5, btagDeepFlavC_dn=sjets[:,1].btagDeepFlavC)
                output['btagDeepB_dn_1'].fill(dataset=dataset,flav=5, btagDeepB_dn=sjets[:,1].btagDeepB)
                output['btagDeepC_dn_1'].fill(dataset=dataset,flav=5, btagDeepC_dn=sjets[:,1].btagDeepC)
                output['deepcsv_CvB_dn_1'].fill(dataset=dataset,flav=5, deepcsv_CvB_dn=np.where(sjets[:,1].btagDeepCvB<0,-0.2,sjets[:,1].btagDeepCvB))
                output['deepcsv_CvL_dn_1'].fill(dataset=dataset,flav=5, deepcsv_CvL_dn=np.where(sjets[:,1].btagDeepCvL<0,-0.2,sjets[:,1].btagDeepCvL))
                output['deepflav_CvB_dn_1'].fill(dataset=dataset,flav=5, deepflav_CvB_dn=np.where(sjets[:,1].btagDeepFlavCvB<0,-0.2,sjets[:,1].btagDeepFlavCvB))
                output['deepflav_CvL_dn_1'].fill(dataset=dataset,flav=5, deepflav_CvL_dn=np.where(sjets[:,1].btagDeepFlavCvL<0,-0.2,sjets[:,1].btagDeepFlavCvL))
            
        
        
         
        # gc.collect()
        # schedule.every(20).minutes.do(dosomething)

        return output

    def postprocess(self, accumulator):
        return accumulator
