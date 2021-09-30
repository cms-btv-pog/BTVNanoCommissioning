import gzip
import pickle, os, sys, mplhep as hep, numpy as np

from matplotlib.pyplot import jet

import coffea
from coffea import hist, processor
import awkward as ak
from coffea.analysis_tools import Weights
from coffea.lumi_tools import LumiMask
from coffea.btag_tools import BTagScaleFactor
from cTagSFReader import *
import gc
# from pympler import tracker,summary,muppy
# import objgraph
# import schedule

deepcsvb_sf = BTagScaleFactor("data/DeepCSV_94XSF_V5_B_F.csv",BTagScaleFactor.RESHAPE,methods='iterativefit,iterativefit,iterativefit')
deepcsvc_sf = "data/DeepCSV_ctagSF_MiniAOD94X_2017_pTincl_v3_2_interp.root"
deepjetb_sf = BTagScaleFactor("data/DeepFlavour_94XSF_V4_B_F.csv",BTagScaleFactor.RESHAPE,methods='iterativefit,iterativefit,iterativefit')
deepjetc_sf = "data/DeepJet_ctagSF_MiniAOD94X_2017_pTincl_v3_2_interp.root"
##import corrections, masks
with gzip.open("data/corrections.pkl.gz") as fin:
    compiled = pickle.load(fin)
def build_lumimask(filename):
    return LumiMask(filename)
lumiMasks = {
    '2016': build_lumimask('data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'),
    '2017': build_lumimask('data/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt'),
    '2018': build_lumimask('data/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'),
}





class NanoProcessor(processor.ProcessorABC):
    # Define histograms
    def flatten(ar): # flatten awkward into a 1d array to hist
            return ak.flatten(ar, axis=None)

    def num(ar):
        return ak.num(ak.fill_none(ar[~ak.is_none(ar)], 0), axis=0)
    def __init__(self):        
        # Define axes
        # Should read axes from NanoAOD config
        dataset_axis = hist.Cat("dataset", "Primary dataset")
        flav_axis = hist.Bin("flav", r"Genflavour",[0,1,4,5,6])
        charge_axis = hist.Bin("char", r"Charge", [1,-1])
        cutflow_axis   = hist.Cat("cut",   "Cut")

        # Events
        njet_axis  = hist.Bin("njet",  r"N jets",      [0,1,2,3])
        hard_lpt_axis = hist.Bin("hard_lpt",   r"Hard lepton pt", 40, 0, 200)
        soft_lpt_axis   = hist.Bin("soft_lpt",   r"Soft lepton pt", 40, 0, 40)
        l_eta_axis  = hist.Bin("eta",  r"Lepton $\eta$", 25, -2.5, 2.5)
        l_phi_axis  = hist.Bin("phi",  r"Lepton $\phi$", 30, -3, 3)
        l_iso_axis = hist.Bin("pfRelIso04_all", r"Rel. Iso", 40,0,4.)
        l_dxy_axis = hist.Bin("dxy", r"dxy", 20,0,0.002)    
        l_dz_axis = hist.Bin("dz", r"dz", 20,0,0.01)    
        l_sip3d_axis = hist.Bin("dz", r"dz", 20,0,0.2)
        soft_lpt_ratio_axis = hist.Bin("soft_lpt_ratio", r"Soft $\mu$/Jet $p_{T}$ [GeV]",40,0,1)

        # Jet
        jet_pt_axis   = hist.Bin("pt",   r"Jet $p_{T}$ [GeV]", 50, 0, 500)
        jet_eta_axis  = hist.Bin("eta",  r"Jet $\eta$", 25, -2.5, 2.5)
        jet_phi_axis  = hist.Bin("phi",  r"Jet $\phi$", 30, -3, 3)
        jet_mass_axis = hist.Bin("mass", r"Jet $m$ [GeV]", 50, 0, 500)
        # jet_dr_axis = hist.Bin("dr", r"Jet $\Delta$R(l,j)", 20, 0, 5)
        ljpt_axis     = hist.Bin("ljpt", r"Leading jet $p_{T}$ [GeV]", 50, 0, 500)
        sljpt_axis     = hist.Bin("sljpt", r"Subleading jet $p_{T}$ [GeV]", 50, 0, 500)
        ljeta_axis = hist.Bin("ljeta",  r"Leading Jet $\eta$", 25, -2.5, 2.5)
        sljeta_axis = hist.Bin("sljeta",  r"Leading Jet $\eta$", 25, -2.5, 2.5)
        
        # ljdr_axis = hist.Bin("ljdr", "Leading jet $\Delta$R(l,j)", 20,0,5)
        # sljdr_axis = hist.Bin("sljdr", "Subleading jet $\Delta$R(l,j)", 20,0,5)
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
                deepcsv_axes.append(hist.Bin(d, d, 20, 0, 5))
            elif "ValAboveCharm" in d:
                deepcsv_axes.append(hist.Bin(d, d, 30, 0, 0.3))
            elif "EnergyRatio" in d:
                deepcsv_axes.append(hist.Bin(d, d, 20, 0, 1))
            elif "trackJetDistVal" in d:
                deepcsv_axes.append(hist.Bin(d, d, 20, 0., 0.1))
            elif "trackPtRatio" in d:
                deepcsv_axes.append(hist.Bin(d, d, 30, 0, 0.3))
            elif "DeltaR" in d: 
                deepcsv_axes.append(hist.Bin(d, d, 30, 0, 0.3))
            elif "jetNSecondaryVertices" in d:
                deepcsv_axes.append(hist.Bin(d, d, 5, 0, 5))
            elif "vertexCategory" in d:
                deepcsv_axes.append(hist.Bin(d, d, 5, 0, 5))
            elif "jetNSelectedTracks" in d:
                deepcsv_axes.append(hist.Bin(d, d, 5, 0, 5))
            elif "vertexNTracks" in d:
                deepcsv_axes.append(hist.Bin(d, d, 5, 0, 5))
            else:
                deepcsv_axes.append(hist.Bin(d, d, 25, 0, 5.))
        # Define similar axes dynamically
        disc_list = ['btagDeepB', 'btagDeepC', 'btagDeepFlavB', 'btagDeepFlavC','deepcsv_CvL','deepcsv_CvB','deepflav_CvL','deepflav_CvB']
        syst_list = ['','SF','_up','_dn']
        varlist=[]
        btag_axes = []
        for d in disc_list:
            for s in syst_list:
                btag_axes.append(hist.Bin("%s%s" %(d,s), "%s%s" %(d,s), 50, 0, 1))  
                varlist.append("%s%s" %(d,s))
                
        _hist_sf_dict={}   
        _hist_deepcsv_dict={}
        for disc, axis in zip(varlist, btag_axes):
            for i in range(2):
                _hist_sf_dict["%s_%d" %(disc,i)] = hist.Hist("Counts", dataset_axis, flav_axis, jet_etawide_axis, jet_ptwide_axis,axis)
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

        self.sf_dict = list(_hist_sf_dict.keys())
        self.event_hists = list(_hist_event_dict.keys())
        self.deepcsv_hists = list(_hist_deepcsv_dict.keys())
        _hist_dict = {**_hist_sf_dict,**_hist_event_dict}
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
        if not isRealData:
            weights.add('genweight',events.genWeight)
            weights.add('puweight', compiled['2017_pileupweight'](events.Pileup.nPU))
        ##############
        # Trigger level
        triggers = [
        # "HLT_IsoMu24",
        "HLT_IsoMu27",
        ]
        
        trig_arrs = [events.HLT[_trig.strip("HLT_")] for _trig in triggers]
        req_trig = np.zeros(len(events), dtype='bool')
        for t in trig_arrs:
            req_trig = req_trig | t

        ############
        # Event level
        ## Muon cuts
        # muon twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
        iso_muon = events.Muon[(events.Muon.pt > 30) & (abs(events.Muon.eta) < 2.4) & (events.Muon.tightId > .5)&(events.Muon.pfRelIso04_all<=0.15)]
        iso_muon = ak.pad_none(iso_muon,1, axis=1)
        req_muon =(ak.count(iso_muon.pt, axis=1) == 1)


        dilep_mu = events.Muon[(events.Muon.pt>12)&(abs(events.Muon.eta) < 2.4)& (events.Muon.tightId > .5)&(events.Muon.pfRelIso04_all<=0.15)]
        dilep_ele = events.Electron[(events.Electron.pt>15)&((abs(events.Electron.eta) < 1.4442)|((abs(events.Electron.eta) < 2.5)&(abs(events.Electron.eta) >1.566)))& (events.Electron.mvaFall17V2Iso_WP80 > .5)]
        req_dilepveto = (ak.count(dilep_mu.pt,axis=1)+ak.count(dilep_ele.pt,axis=1)==2)

        # MET = events.METFixEE2017
        
        MET =ak.zip({
            "pt": events.METFixEE2017.pt,
            "eta": ak.zeros_like(events.METFixEE2017.pt),
            "phi": events.METFixEE2017.phi,
            "mass":ak.zeros_like(events.METFixEE2017.pt),
            }, with_name="PtEtaPhiMLorentzVector")

        Wmass = MET+iso_muon[:,0]
        req_Wmass = (Wmass.mass>55)
        ## Jet cuts 
        event_jet =  events.Jet[(events.Jet.pt > 20) & (abs(events.Jet.eta) <= 2.5)&((events.Jet.puId >=7)&( events.Jet.pt<50)) &(events.Jet.jetId>=3)&(events.Jet.btagDeepB>0.) & (events.Jet.btagDeepB<1.) & (events.Jet.btagDeepC>0.) & (events.Jet.btagDeepC<1.) & (events.Jet.btagDeepFlavB>0.) & (events.Jet.btagDeepFlavB<1.) & (events.Jet.btagDeepFlavC>0.) & (events.Jet.btagDeepFlavC<1.)& (ak.all(events.Jet.metric_table(iso_muon[:,0]) > 0.5, axis=2))&((events.Jet.muEF+events.Jet.neEmEF)<0.7)]

        req_jets = ((ak.num(event_jet.puId) >=1) & (ak.num(event_jet.puId) <=3))
        
        ## Soft Muon cuts 
        
        soft_muon = events.Muon[(events.Muon.pt < 25) & (abs(events.Muon.eta) < 2.4) & (events.Muon.tightId > .5)&(events.Muon.pfRelIso04_all>0.2)&(events.Muon.jetIdx!=-1)]
        req_softmu=(ak.count(soft_muon.pt, axis=1) >= 1)
        soft_muon= ak.pad_none(soft_muon,1,axis=1)
        
        mu_jet = event_jet[(ak.all(event_jet.metric_table(soft_muon) <= 0.4, axis=2))&((event_jet.muonIdx1!=-1)|(event_jet.muonIdx2!=-1))]
        req_mujet = (ak.count(mu_jet.pt,axis=1)>=1)
        # mu_jet = ak.fill_none(mu_jet.pt,0)
        mu_jet = ak.pad_none(mu_jet,1, axis=1)

        # pT ratio
        req_pTratio = ((soft_muon[:,0].pt/mu_jet[:,0].pt)<0.4)
        #dilepton mass
        req_dilepmass = np.zeros(len(events), dtype='bool')
        req_dilepmass = req_muon&req_softmu
        dilep_mass = iso_muon[:,0]+soft_muon[:,0]
        req_dilepmass = ((dilep_mass.mass>12.)&((dilep_mass.mass<80)| (dilep_mass.mass>100)))

        req_QCDveto = (iso_muon[:,0].pfRelIso04_all<0.05) & (abs(iso_muon[:,0].dz)<0.01) & (abs(iso_muon[:,0].dxy)<0.002) & (iso_muon[:,0].ip3d < 0.2) & ((iso_muon[:,0].pt/mu_jet[:,0].pt<0.)|(iso_muon[:,0].pt/mu_jet[:,0].pt>0.75))
        event_level = req_trig & req_lumi & req_muon &  req_jets & req_softmu   & req_dilepmass &req_mujet & req_Wmass & req_dilepveto & req_QCDveto
        if(len(event_level)>0):event_level = ak.fill_none(event_level,False)        
        # Selected
        selev = events[event_level]    
        
        #########
        
        ## Hard Muon
        shmu = selev.Muon[(selev.Muon.pt > 30) & (abs(selev.Muon.eta) < 2.4) & (selev.Muon.tightId > .5)&(selev.Muon.pfRelIso04_all<=0.15)]
        ## Soft Muon
        ssmu = selev.Muon[(selev.Muon.pt < 25) & (abs(selev.Muon.eta) < 2.4) & (selev.Muon.tightId > .5)&(selev.Muon.pfRelIso04_all>0.2)&(selev.Muon.jetIdx!=-1)]
        ## Jets
        sjets = selev.Jet[(selev.Jet.pt > 20) & (abs(selev.Jet.eta) <= 2.5)&((selev.Jet.puId >=7)&( selev.Jet.pt<50)) &(selev.Jet.jetId>=3)&(selev.Jet.btagDeepB>0.) & (selev.Jet.btagDeepB<1.) & (selev.Jet.btagDeepC>0.) & (selev.Jet.btagDeepC<1.) & (selev.Jet.btagDeepFlavB>0.) & (selev.Jet.btagDeepFlavB<1.) & (selev.Jet.btagDeepFlavC>0.) & (selev.Jet.btagDeepFlavC<1.)& (ak.all(selev.Jet.metric_table(shmu[:,0]) > 0.5, axis=2))&((selev.Jet.muEF+selev.Jet.neEmEF)<0.7)]
        
        ## Muon Jet 
        smuon_jet = sjets[(ak.all(sjets.metric_table(ssmu) <= 0.4, axis=2))&((sjets.muonIdx1!=-1)|(sjets.muonIdx2!=-1))]

    
        
        if isRealData :
            genflavor = ak.zeros_like(sjets.pt)
        else:
            par_flav = (sjets.partonFlavour == 0 ) & (sjets.hadronFlavour==0)
            genflavor = sjets.hadronFlavour + 1*par_flav 
            # genweiev=ak.flatten(ak.broadcast_arrays(weights.weight()[event_level],sjets['pt'])[0])
            jetsfs_c_lj = getSF(ak.to_numpy(sjets[:,0].hadronFlavour),ak.to_numpy(sjets[:,0].btagDeepFlavC/(1.-sjets[:,0].btagDeepFlavB)),ak.to_numpy(sjets[:,0].btagDeepFlavC/(sjets[:,0].btagDeepFlavC+sjets[:,0].btagDeepFlavB)),deepjetc_sf)
            jetsfs_c_up_lj = getSF(ak.to_numpy(sjets[:,0].hadronFlavour),ak.to_numpy(sjets[:,0].btagDeepFlavC/(1.-sjets[:,0].btagDeepFlavB)),ak.to_numpy(sjets[:,0].btagDeepFlavC/(sjets[:,0].btagDeepFlavC+sjets[:,0].btagDeepFlavB)),deepjetc_sf,"TotalUncUp")
            jetsfs_c_dn_lj = getSF(ak.to_numpy(sjets[:,0].hadronFlavour),ak.to_numpy(sjets[:,0].btagDeepFlavC/(1.-sjets[:,0].btagDeepFlavB)),ak.to_numpy(sjets[:,0].btagDeepFlavC/(sjets[:,0].btagDeepFlavC+sjets[:,0].btagDeepFlavB)),deepjetc_sf,"TotalUncDown")
            jetsfs_b_lj = deepjetb_sf.eval('central',sjets[:,0].hadronFlavour,abs(sjets[:,0].eta),sjets[:,0].pt,discr=sjets[:,0].btagDeepFlavB)
            jetsfs_b_up_lj = deepjetb_sf.eval('up_jes',sjets[:,0].hadronFlavour,abs(sjets[:,0].eta),sjets[:,0].pt,discr=sjets[:,0].btagDeepFlavB)
            jetsfs_b_dn_lj = deepjetb_sf.eval('down_jes',sjets[:,0].hadronFlavour,abs(sjets[:,0].eta),sjets[:,0].pt,discr=sjets[:,0].btagDeepFlavB)
            csvsfs_c_lj = getSF(ak.to_numpy(sjets[:,0].hadronFlavour),ak.to_numpy(sjets[:,0].btagDeepFlavC/(1.-sjets[:,0].btagDeepFlavB)),ak.to_numpy(sjets[:,0].btagDeepFlavC/(sjets[:,0].btagDeepFlavC+sjets[:,0].btagDeepFlavB)),deepcsvc_sf)
            csvsfs_c_up_lj = getSF(ak.to_numpy(sjets[:,0].hadronFlavour),ak.to_numpy(sjets[:,0].btagDeepFlavC/(1.-sjets[:,0].btagDeepFlavB)),ak.to_numpy(sjets[:,0].btagDeepFlavC/(sjets[:,0].btagDeepFlavC+sjets[:,0].btagDeepFlavB)),deepcsvc_sf,"TotalUncUp")
            csvsfs_c_dn_lj = getSF(ak.to_numpy(sjets[:,0].hadronFlavour),ak.to_numpy(sjets[:,0].btagDeepFlavC/(1.-sjets[:,0].btagDeepFlavB)),ak.to_numpy(sjets[:,0].btagDeepFlavC/(sjets[:,0].btagDeepFlavC+sjets[:,0].btagDeepFlavB)),deepcsvc_sf,"TotalUncDown")
            csvsfs_b_up_lj = deepcsvb_sf.eval('up_jes',sjets[:,0].hadronFlavour,abs(sjets[:,0].eta),sjets[:,0].pt,discr=sjets[:,0].btagDeepB)
            csvsfs_b_lj = deepcsvb_sf.eval('central',sjets[:,0].hadronFlavour,abs(sjets[:,0].eta),sjets[:,0].pt,discr=sjets[:,0].btagDeepB)
            csvsfs_b_dn_lj = deepcsvb_sf.eval('down_jes',sjets[:,0].hadronFlavour,abs(sjets[:,0].eta),sjets[:,0].pt,discr=sjets[:,0].btagDeepB)
            jetsfs_c_sj = getSF(ak.to_numpy(sjets[:,1].hadronFlavour),ak.to_numpy(sjets[:,1].btagDeepFlavC/(1.-sjets[:,1].btagDeepFlavB)),ak.to_numpy(sjets[:,1].btagDeepFlavC/(sjets[:,1].btagDeepFlavC+sjets[:,1].btagDeepFlavB)),deepjetc_sf)
            csvsfs_c_sj = getSF(ak.to_numpy(sjets[:,1].hadronFlavour),ak.to_numpy(sjets[:,1].btagDeepFlavC/(1.-sjets[:,1].btagDeepFlavB)),ak.to_numpy(sjets[:,1].btagDeepFlavC/(sjets[:,1].btagDeepFlavC+sjets[:,1].btagDeepFlavB)),deepcsvc_sf)
            jetsfs_c_up_sj = getSF(ak.to_numpy(sjets[:,1].hadronFlavour),ak.to_numpy(sjets[:,1].btagDeepFlavC/(1.-sjets[:,1].btagDeepFlavB)),ak.to_numpy(sjets[:,1].btagDeepFlavC/(sjets[:,1].btagDeepFlavC+sjets[:,1].btagDeepFlavB)),deepjetc_sf,"TotalUncUp")
            csvsfs_c_up_sj = getSF(ak.to_numpy(sjets[:,1].hadronFlavour),ak.to_numpy(sjets[:,1].btagDeepFlavC/(1.-sjets[:,1].btagDeepFlavB)),ak.to_numpy(sjets[:,1].btagDeepFlavC/(sjets[:,1].btagDeepFlavC+sjets[:,1].btagDeepFlavB)),deepcsvc_sf,"TotalUncUp")
            jetsfs_c_dn_sj = getSF(ak.to_numpy(sjets[:,1].hadronFlavour),ak.to_numpy(sjets[:,1].btagDeepFlavC/(1.-sjets[:,1].btagDeepFlavB)),ak.to_numpy(sjets[:,1].btagDeepFlavC/(sjets[:,1].btagDeepFlavC+sjets[:,1].btagDeepFlavB)),deepjetc_sf,"TotalUncDown")
            csvsfs_c_dn_sj = getSF(ak.to_numpy(sjets[:,1].hadronFlavour),ak.to_numpy(sjets[:,1].btagDeepFlavC/(1.-sjets[:,1].btagDeepFlavB)),ak.to_numpy(sjets[:,1].btagDeepFlavC/(sjets[:,1].btagDeepFlavC+sjets[:,1].btagDeepFlavB)),deepcsvc_sf,"TotalUncDown")
            jetsfs_b_sj = deepjetb_sf.eval('central',sjets[:,1].hadronFlavour,abs(sjets[:,1].eta),sjets[:,1].pt,discr=sjets[:,1].btagDeepFlavB)
            csvsfs_b_sj = deepcsvb_sf.eval('central',sjets[:,1].hadronFlavour,abs(sjets[:,1].eta),sjets[:,1].pt,discr=sjets[:,1].btagDeepB)
            jetsfs_b_up_sj = deepjetb_sf.eval('up_jes',sjets[:,1].hadronFlavour,abs(sjets[:,1].eta),sjets[:,1].pt,discr=sjets[:,1].btagDeepFlavB)
            csvsfs_b_up_sj = deepcsvb_sf.eval('up_jes',sjets[:,1].hadronFlavour,abs(sjets[:,1].eta),sjets[:,1].pt,discr=sjets[:,1].btagDeepB)
            jetsfs_b_dn_sj = deepjetb_sf.eval('down_jes',sjets[:,1].hadronFlavour,abs(sjets[:,1].eta),sjets[:,1].pt,discr=sjets[:,1].btagDeepFlavB)
            csvsfs_b_dn_sj = deepcsvb_sf.eval('down_jes',sjets[:,1].hadronFlavour,abs(sjets[:,1].eta),sjets[:,1].pt,discr=sjets[:,1].btagDeepB)
            
        for histname, h in output.items():
            if histname in self.deepcsv_hists:
                if(isRealData):
                    fields = {l: ak.flatten(sjets[l], axis=None) for l in h.fields if l in dir(sjets)}
                    h.fill(dataset=dataset,flav=5, **fields)
                else:
                    fields = {l: ak.flatten(sjets[histname]) for l in h.fields if l in dir(sjets)}
                    genweiev=ak.flatten(ak.broadcast_arrays(weights.weight()[event_level],sjets['pt'])[0])
                    h.fill(dataset=dataset,flav=ak.flatten(genflavor), **fields,weight=genweiev)
       
            
        if not isRealData:
            ###Fill no SFs
            # output['btagDeepFlavB_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepFlavB=sjets[:,0].btagDeepFlavB,weight=weights.weight()[event_level])
            # output['btagDeepFlavC_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepFlavC=sjets[:,0].btagDeepFlavC,weight=weights.weight()[event_level])
            # output['btagDeepB_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepB=sjets[:,0].btagDeepB,weight=weights.weight()[event_level])
            # output['btagDeepC_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepC=sjets[:,0].btagDeepC,weight=weights.weight()[event_level])
            # output['deepcsv_CvB_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepcsv_CvB=sjets[:,0].btagDeepC/(1.-sjets[:,0].btagDeepB),weight=weights.weight()[event_level])
            # output['deepcsv_CvL_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepcsv_CvL=sjets[:,0].btagDeepC/(sjets[:,0].btagDeepC+sjets[:,0].btagDeepB),weight=weights.weight()[event_level])
            # output['deepflav_CvB_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepflav_CvB=sjets[:,0].btagDeepFlavC/(1.-sjets[:,0].btagDeepFlavB),weight=weights.weight()[event_level])
            # output['deepflav_CvL_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepflav_CvL=sjets[:,0].btagDeepFlavC/(sjets[:,0].btagDeepFlavC+sjets[:,0].btagDeepFlavB),weight=weights.weight()[event_level])
            # output['btagDeepFlavB_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepFlavB=sjets[:,1].btagDeepFlavB,weight=weights.weight()[event_level])
            # output['btagDeepFlavC_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepFlavC=sjets[:,1].btagDeepFlavC,weight=weights.weight()[event_level])
            # output['btagDeepB_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepB=sjets[:,1].btagDeepB,weight=weights.weight()[event_level])
            # output['btagDeepC_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepC=sjets[:,1].btagDeepC,weight=weights.weight()[event_level])
            # output['deepcsv_CvB_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepcsv_CvB=sjets[:,1].btagDeepC/(1.-sjets[:,1].btagDeepB),weight=weights.weight()[event_level])
            # output['deepcsv_CvL_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepcsv_CvL=sjets[:,1].btagDeepC/(sjets[:,1].btagDeepC+sjets[:,1].btagDeepB),weight=weights.weight()[event_level])
            # output['deepflav_CvB_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepflav_CvB=sjets[:,1].btagDeepFlavC/(1.-sjets[:,1].btagDeepFlavB),weight=weights.weight()[event_level])
            # output['deepflav_CvL_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepflav_CvL=sjets[:,1].btagDeepFlavC/(sjets[:,1].btagDeepFlavC+sjets[:,1].btagDeepFlavB),weight=weights.weight()[event_level])
            # ##SFs 
            # print(ak.type(jetsfs_b_lj),jetsfs_b_lj)
            # output['btagDeepFlavBSF_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepFlavBSF=sjets[:,0].btagDeepFlavB,weight=(weights.weight()[event_level]*jetsfs_b_lj))      
            # output['btagDeepFlavCSF_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepFlavCSF=sjets[:,0].btagDeepFlavC,weight=weights.weight()[event_level]*jetsfs_c_lj)
            # output['btagDeepBSF_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepBSF=sjets[:,0].btagDeepB,weight=weights.weight()[event_level]*csvsfs_b_lj)
            # output['btagDeepCSF_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepCSF=sjets[:,0].btagDeepC,weight=weights.weight()[event_level]*csvsfs_c_lj)
            # output['deepcsv_CvBSF_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepcsv_CvBSF=sjets[:,0].btagDeepC/(1.-sjets[:,0].btagDeepB),weight=weights.weight()[event_level]*csvsfs_c_lj)
            # output['deepcsv_CvLSF_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepcsv_CvLSF=sjets[:,0].btagDeepC/(sjets[:,0].btagDeepC+sjets[:,0].btagDeepB),weight=weights.weight()[event_level]*csvsfs_c_lj)
            # output['deepflav_CvBSF_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepflav_CvBSF=sjets[:,0].btagDeepFlavC/(1.-sjets[:,0].btagDeepFlavB),weight=weights.weight()[event_level]*jetsfs_c_lj)
            # output['deepflav_CvLSF_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepflav_CvLSF=sjets[:,0].btagDeepFlavC/(sjets[:,0].btagDeepFlavC+sjets[:,0].btagDeepFlavB),weight=weights.weight()[event_level]*jetsfs_c_lj)
            # output['btagDeepFlavBSF_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepFlavBSF=sjets[:,1].btagDeepFlavB,weight=weights.weight()[event_level]*jetsfs_b_sj)
            # output['btagDeepFlavCSF_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepFlavCSF=sjets[:,1].btagDeepFlavC,weight=weights.weight()[event_level]*jetsfs_c_sj)
            # output['btagDeepBSF_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepBSF=sjets[:,1].btagDeepB,weight=weights.weight()[event_level]*csvsfs_b_sj)
            # output['btagDeepCSF_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepCSF=sjets[:,1].btagDeepC,weight=weights.weight()[event_level]*csvsfs_c_sj)
            # output['deepcsv_CvBSF_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepcsv_CvBSF=sjets[:,1].btagDeepC/(1.-sjets[:,1].btagDeepB),weight=weights.weight()[event_level]*csvsfs_c_sj)
            # output['deepcsv_CvLSF_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepcsv_CvLSF=sjets[:,1].btagDeepC/(sjets[:,1].btagDeepC+sjets[:,1].btagDeepB),weight=weights.weight()[event_level]*csvsfs_c_sj)
            # output['deepflav_CvBSF_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepflav_CvBSF=sjets[:,1].btagDeepFlavC/(1.-sjets[:,1].btagDeepFlavB),weight=weights.weight()[event_level]*jetsfs_c_sj)
            # output['deepflav_CvLSF_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepflav_CvLSF=sjets[:,1].btagDeepFlavC/(sjets[:,1].btagDeepFlavC+sjets[:,1].btagDeepFlavB),weight=weights.weight()[event_level]*jetsfs_c_sj)
            # ###SFs up
            # output['btagDeepFlavB_up_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepFlavB_up=sjets[:,0].btagDeepFlavB,weight=weights.weight()[event_level]*jetsfs_b_up_lj)
            # output['btagDeepFlavC_up_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepFlavC_up=sjets[:,0].btagDeepFlavC,weight=weights.weight()[event_level]*jetsfs_c_up_lj)
            # output['btagDeepB_up_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepB_up=sjets[:,0].btagDeepB,weight=weights.weight()[event_level]*csvsfs_b_up_lj)
            # output['btagDeepC_up_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepC_up=sjets[:,0].btagDeepC,weight=weights.weight()[event_level]*csvsfs_c_up_lj)
            # output['deepcsv_CvB_up_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepcsv_CvB_up=sjets[:,0].btagDeepC/(1.-sjets[:,0].btagDeepB),weight=weights.weight()[event_level]*csvsfs_c_up_lj)
            # output['deepcsv_CvL_up_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepcsv_CvL_up=sjets[:,0].btagDeepC/(sjets[:,0].btagDeepC+sjets[:,0].btagDeepB),weight=weights.weight()[event_level]*csvsfs_c_up_lj)
            # output['deepflav_CvB_up_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepflav_CvB_up=sjets[:,0].btagDeepFlavC/(1.-sjets[:,0].btagDeepFlavB),weight=weights.weight()[event_level]*jetsfs_c_up_lj)
            # output['deepflav_CvL_up_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepflav_CvL_up=sjets[:,0].btagDeepFlavC/(sjets[:,0].btagDeepFlavC+sjets[:,0].btagDeepFlavB),weight=weights.weight()[event_level]*jetsfs_c_up_lj)
            # output['btagDeepFlavB_up_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepFlavB_up=sjets[:,1].btagDeepFlavB,weight=weights.weight()[event_level]*jetsfs_b_up_sj)
            # output['btagDeepFlavC_up_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepFlavC_up=sjets[:,1].btagDeepFlavC,weight=weights.weight()[event_level]*jetsfs_c_up_sj)
            # output['btagDeepB_up_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepB_up=sjets[:,1].btagDeepB,weight=weights.weight()[event_level]*csvsfs_b_up_sj)
            # output['btagDeepC_up_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepC_up=sjets[:,1].btagDeepC,weight=weights.weight()[event_level]*csvsfs_c_up_sj)
            # output['deepcsv_CvB_up_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepcsv_CvB_up=sjets[:,1].btagDeepC/(1.-sjets[:,1].btagDeepB),weight=weights.weight()[event_level]*csvsfs_c_up_sj)
            # output['deepcsv_CvL_up_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepcsv_CvL_up=sjets[:,1].btagDeepC/(sjets[:,1].btagDeepC+sjets[:,1].btagDeepB),weight=weights.weight()[event_level]*csvsfs_c_up_sj)
            # output['deepflav_CvB_up_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepflav_CvB_up=sjets[:,1].btagDeepFlavC/(1.-sjets[:,1].btagDeepFlavB),weight=weights.weight()[event_level]*jetsfs_c_up_sj)
            # output['deepflav_CvL_up_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepflav_CvL_up=sjets[:,1].btagDeepFlavC/(sjets[:,1].btagDeepFlavC+sjets[:,1].btagDeepFlavB),weight=weights.weight()[event_level]*jetsfs_c_up_sj)
            # ###SFs down
            # output['btagDeepFlavB_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepFlavB_dn=sjets[:,0].btagDeepFlavB,weight=weights.weight()[event_level]*jetsfs_b_dn_lj)
            # output['btagDeepFlavC_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepFlavC_dn=sjets[:,0].btagDeepFlavC,weight=weights.weight()[event_level]*jetsfs_c_dn_lj)
            # output['btagDeepB_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepB_dn=sjets[:,0].btagDeepB,weight=weights.weight()[event_level]*csvsfs_b_dn_lj)
            # output['btagDeepC_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepC_dn=sjets[:,0].btagDeepC,weight=weights.weight()[event_level]*csvsfs_c_dn_lj)
            # output['deepcsv_CvB_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepcsv_CvB_dn=sjets[:,0].btagDeepC/(1.-sjets[:,0].btagDeepB),weight=weights.weight()[event_level]*csvsfs_c_dn_lj)
            # output['deepcsv_CvL_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepcsv_CvL_dn=sjets[:,0].btagDeepC/(sjets[:,0].btagDeepC+sjets[:,0].btagDeepB),weight=weights.weight()[event_level]*csvsfs_c_dn_lj)
            # output['deepflav_CvB_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepflav_CvB_dn=sjets[:,0].btagDeepFlavC/(1.-sjets[:,0].btagDeepFlavB),weight=weights.weight()[event_level]*jetsfs_c_dn_lj)
            # output['deepflav_CvL_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepflav_CvL_dn=sjets[:,0].btagDeepFlavC/(sjets[:,0].btagDeepFlavC+sjets[:,0].btagDeepFlavB),weight=weights.weight()[event_level]*jetsfs_c_dn_lj)
            # output['btagDeepFlavB_dn_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepFlavB_dn=sjets[:,1].btagDeepFlavB,weight=weights.weight()[event_level]*jetsfs_b_dn_sj)
            # output['btagDeepFlavC_dn_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepFlavC_dn=sjets[:,1].btagDeepFlavC,weight=weights.weight()[event_level]*jetsfs_c_dn_sj)
            # output['btagDeepB_dn_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepB_dn=sjets[:,1].btagDeepB,weight=weights.weight()[event_level]*csvsfs_b_dn_sj)
            # output['btagDeepC_dn_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepC_dn=sjets[:,1].btagDeepC,weight=weights.weight()[event_level]*csvsfs_c_dn_sj)
            # output['deepcsv_CvB_dn_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepcsv_CvB_dn=sjets[:,1].btagDeepC/(1.-sjets[:,1].btagDeepB),weight=weights.weight()[event_level]*csvsfs_c_dn_sj)
            # output['deepcsv_CvL_dn_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepcsv_CvL_dn=sjets[:,1].btagDeepC/(sjets[:,1].btagDeepC+sjets[:,1].btagDeepB),weight=weights.weight()[event_level]*csvsfs_c_dn_sj)
            # output['deepflav_CvB_dn_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepflav_CvB_dn=sjets[:,1].btagDeepFlavC/(1.-sjets[:,1].btagDeepFlavB),weight=weights.weight()[event_level]*jetsfs_c_dn_sj)
            # output['deepflav_CvL_dn_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepflav_CvL_dn=sjets[:,1].btagDeepFlavC/(sjets[:,1].btagDeepFlavC+sjets[:,1].btagDeepFlavB),weight=weights.weight()[event_level]*jetsfs_c_dn_sj)
            
        else:
            # output['btagDeepFlavB_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepFlavB=sjets[:,0].btagDeepFlavB)
            # output['btagDeepFlavC_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepFlavC=sjets[:,0].btagDeepFlavC)
            # output['btagDeepB_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepB=sjets[:,0].btagDeepB)
            # output['btagDeepC_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepC=sjets[:,0].btagDeepC)
            # output['deepcsv_CvB_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepcsv_CvB=sjets[:,0].btagDeepC/(1.-sjets[:,0].btagDeepB))
            # output['deepcsv_CvL_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepcsv_CvL=sjets[:,0].btagDeepC/(sjets[:,0].btagDeepC+sjets[:,0].btagDeepB))
            # output['deepflav_CvB_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepflav_CvB=sjets[:,0].btagDeepFlavC/(1.-sjets[:,0].btagDeepFlavB))
            # output['deepflav_CvL_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepflav_CvL=sjets[:,0].btagDeepFlavC/(sjets[:,0].btagDeepFlavC+sjets[:,0].btagDeepFlavB))
            # output['btagDeepFlavB_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepFlavB=sjets[:,1].btagDeepFlavB)
            # output['btagDeepFlavC_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepFlavC=sjets[:,1].btagDeepFlavC)
            # output['btagDeepB_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepB=sjets[:,1].btagDeepB)
            # output['btagDeepC_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepC=sjets[:,1].btagDeepC)
            # output['deepcsv_CvB_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepcsv_CvB=sjets[:,1].btagDeepC/(1.-sjets[:,1].btagDeepB))
            # output['deepcsv_CvL_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepcsv_CvL=sjets[:,1].btagDeepC/(sjets[:,1].btagDeepC+sjets[:,1].btagDeepB))
            # output['deepflav_CvB_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepflav_CvB=sjets[:,1].btagDeepFlavC/(1.-sjets[:,1].btagDeepFlavB))
            # output['deepflav_CvL_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepflav_CvL=sjets[:,1].btagDeepFlavC/(sjets[:,1].btagDeepFlavC+sjets[:,1].btagDeepFlavB))
            # ##SFs 
            # # print(ak.type(jetsfs_b_lj),jetsfs_b_lj)
            # output['btagDeepFlavBSF_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepFlavBSF=sjets[:,0].btagDeepFlavB)
            # output['btagDeepFlavCSF_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepFlavCSF=sjets[:,0].btagDeepFlavC)
            # output['btagDeepBSF_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepBSF=sjets[:,0].btagDeepB)
            # output['btagDeepCSF_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepCSF=sjets[:,0].btagDeepC)
            # output['deepcsv_CvBSF_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepcsv_CvBSF=sjets[:,0].btagDeepC/(1.-sjets[:,0].btagDeepB))
            # output['deepcsv_CvLSF_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepcsv_CvLSF=sjets[:,0].btagDeepC/(sjets[:,0].btagDeepC+sjets[:,0].btagDeepB))
            # output['deepflav_CvBSF_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepflav_CvBSF=sjets[:,0].btagDeepFlavC/(1.-sjets[:,0].btagDeepFlavB))
            # output['deepflav_CvLSF_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepflav_CvLSF=sjets[:,0].btagDeepFlavC/(sjets[:,0].btagDeepFlavC+sjets[:,0].btagDeepFlavB))
            # output['btagDeepFlavBSF_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepFlavBSF=sjets[:,1].btagDeepFlavB)
            # output['btagDeepFlavCSF_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepFlavCSF=sjets[:,1].btagDeepFlavC)
            # output['btagDeepBSF_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepBSF=sjets[:,1].btagDeepB)
            # output['btagDeepCSF_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepCSF=sjets[:,1].btagDeepC)
            # output['deepcsv_CvBSF_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepcsv_CvBSF=sjets[:,1].btagDeepC/(1.-sjets[:,1].btagDeepB))
            # output['deepcsv_CvLSF_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepcsv_CvLSF=sjets[:,1].btagDeepC/(sjets[:,1].btagDeepC+sjets[:,1].btagDeepB))
            # output['deepflav_CvBSF_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepflav_CvBSF=sjets[:,1].btagDeepFlavC/(1.-sjets[:,1].btagDeepFlavB))
            # output['deepflav_CvLSF_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepflav_CvLSF=sjets[:,1].btagDeepFlavC/(sjets[:,1].btagDeepFlavC+sjets[:,1].btagDeepFlavB))
            # ###SFs up
            # output['btagDeepFlavB_up_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepFlavB_up=sjets[:,0].btagDeepFlavB)
            # output['btagDeepFlavC_up_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepFlavC_up=sjets[:,0].btagDeepFlavC)
            # output['btagDeepB_up_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepB_up=sjets[:,0].btagDeepB)
            # output['btagDeepC_up_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepC_up=sjets[:,0].btagDeepC)
            # output['deepcsv_CvB_up_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepcsv_CvB_up=sjets[:,0].btagDeepC/(1.-sjets[:,0].btagDeepB))
            # output['deepcsv_CvL_up_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepcsv_CvL_up=sjets[:,0].btagDeepC/(sjets[:,0].btagDeepC+sjets[:,0].btagDeepB))
            # output['deepflav_CvB_up_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepflav_CvB_up=sjets[:,0].btagDeepFlavC/(1.-sjets[:,0].btagDeepFlavB))
            # output['deepflav_CvL_up_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepflav_CvL_up=sjets[:,0].btagDeepFlavC/(sjets[:,0].btagDeepFlavC+sjets[:,0].btagDeepFlavB))
            # output['btagDeepFlavB_up_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepFlavB_up=sjets[:,1].btagDeepFlavB)
            # output['btagDeepFlavC_up_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepFlavC_up=sjets[:,1].btagDeepFlavC)
            # output['btagDeepB_up_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepB_up=sjets[:,1].btagDeepB)
            # output['btagDeepC_up_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepC_up=sjets[:,1].btagDeepC)
            # output['deepcsv_CvB_up_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepcsv_CvB_up=sjets[:,1].btagDeepC/(1.-sjets[:,1].btagDeepB))
            # output['deepcsv_CvL_up_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepcsv_CvL_up=sjets[:,1].btagDeepC/(sjets[:,1].btagDeepC+sjets[:,1].btagDeepB))
            # output['deepflav_CvB_up_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepflav_CvB_up=sjets[:,1].btagDeepFlavC/(1.-sjets[:,1].btagDeepFlavB))
            # output['deepflav_CvL_up_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepflav_CvL_up=sjets[:,1].btagDeepFlavC/(sjets[:,1].btagDeepFlavC+sjets[:,1].btagDeepFlavB))
            # # ###SFs down
            # output['btagDeepFlavB_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepFlavB_dn=sjets[:,0].btagDeepFlavB)
            # output['btagDeepFlavC_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepFlavC_dn=sjets[:,0].btagDeepFlavC)
            # output['btagDeepB_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepB_dn=sjets[:,0].btagDeepB)
            # output['btagDeepC_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, btagDeepC_dn=sjets[:,0].btagDeepC)
            # output['deepcsv_CvB_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepcsv_CvB_dn=sjets[:,0].btagDeepC/(1.-sjets[:,0].btagDeepB))
            # output['deepcsv_CvL_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepcsv_CvL_dn=sjets[:,0].btagDeepC/(sjets[:,0].btagDeepC+sjets[:,0].btagDeepB))
            # output['deepflav_CvB_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepflav_CvB_dn=sjets[:,0].btagDeepFlavC/(1.-sjets[:,0].btagDeepFlavB))
            # output['deepflav_CvL_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], etawide=sjets[:,0].eta,ptwide=sjets[:,0].pt, deepflav_CvL_dn=sjets[:,0].btagDeepFlavC/(sjets[:,0].btagDeepFlavC+sjets[:,0].btagDeepFlavB))
            # output['btagDeepFlavB_dn_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepFlavB_dn=sjets[:,1].btagDeepFlavB)
            # output['btagDeepFlavC_dn_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepFlavC_dn=sjets[:,1].btagDeepFlavC)
            # output['btagDeepB_dn_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepB_dn=sjets[:,1].btagDeepB)
            # output['btagDeepC_dn_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, btagDeepC_dn=sjets[:,1].btagDeepC)
            # output['deepcsv_CvB_dn_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepcsv_CvB_dn=sjets[:,1].btagDeepC/(1.-sjets[:,1].btagDeepB))
            # output['deepcsv_CvL_dn_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepcsv_CvL_dn=sjets[:,1].btagDeepC/(sjets[:,1].btagDeepC+sjets[:,1].btagDeepB))
            # output['deepflav_CvB_dn_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepflav_CvB_dn=sjets[:,1].btagDeepFlavC/(1.-sjets[:,1].btagDeepFlavB))
            # output['deepflav_CvL_dn_1'].fill(dataset=dataset,flav=genflavor[:,1], etawide=sjets[:,1].eta,ptwide=sjets[:,1].pt, deepflav_CvL_dn=sjets[:,1].btagDeepFlavC/(sjets[:,1].btagDeepFlavC+sjets[:,1].btagDeepFlavB))

        '''
        if(isRealData):
            output['njet'].fill(dataset=dataset,  njet=flatten(ak.num(sel_jets)))
            output['nbjet'].fill(dataset=dataset,  nbjet=flatten(ak.num(sbjets)))
            output['ntbjet'].fill(dataset=dataset,  ntbjet=flatten(ak.num(sbjets)))
            output['ljpt'].fill(dataset=dataset, flav=0, ljpt=flatten(sjets[:,0].pt))
            output['sljpt'].fill(dataset=dataset, flav=0, sljpt=flatten(sjets[:,1].pt))
            output['ljeta'].fill(dataset=dataset, flav=0, ljeta=flatten(sjets[:,0].eta))
            output['sljeta'].fill(dataset=dataset, flav=0, sljeta=flatten(sjets[:,1].eta))
        else:
            output['njet'].fill(dataset=dataset,  njet=flatten(ak.num(sel_jets)),weight=weights.weight()[event_level])
            output['nbjet'].fill(dataset=dataset,  nbjet=flatten(ak.num(sbjets)),weight=weights.weight()[event_level]*csvsfs_b_lj*csvsfs_b_sj)
            output['nbjet_up'].fill(dataset=dataset,  nbjet=flatten(ak.num(sbjets)),weight=weights.weight()[event_level]*csvsfs_b_up_lj*csvsfs_b_up_sj)
            output['nbjet_dn'].fill(dataset=dataset,  nbjet=flatten(ak.num(sbjets)),weight=weights.weight()[event_level]*csvsfs_b_dn_lj*csvsfs_b_dn_sj)
            output['ntbjet'].fill(dataset=dataset,  ntbjet=flatten(ak.num(stbjets)),weight=weights.weight()[event_level])
            output['ljpt'].fill(dataset=dataset, flav=flatten(sjets[:,0].hadronFlavour), ljpt=flatten(sjets[:,0].pt),weight=weights.weight()[event_level])
            output['sljpt'].fill(dataset=dataset, flav=flatten(sjets[:,1].hadronFlavour), sljpt=flatten(sjets[:,1].pt),weight=weights.weight()[event_level])
            output['ljeta'].fill(dataset=dataset, flav=flatten(sjets[:,0].hadronFlavour), ljeta=flatten(sjets[:,0].eta),weight=weights.weight()[event_level])
            output['sljeta'].fill(dataset=dataset, flav=flatten(sjets[:,1].hadronFlavour), sljeta=flatten(sjets[:,1].eta),weight=weights.weight()[event_level])
        '''
        gc.collect()
        # schedule.every(20).minutes.do(dosomething)

        return output

    def postprocess(self, accumulator):
        return accumulator
