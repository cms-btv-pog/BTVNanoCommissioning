import gzip
import pickle
import coffea
from coffea import hist, processor
import numpy as np
import awkward as ak
from coffea.analysis_tools import Weights
from BTVNanoCommissioning.utils.correction import lumiMasks, compiled, eleSFs,muSFs,deepcsvb_sf,deepcsvc_sf,deepjetb_sf,deepjetc_sf
import gc
from BTVNanoCommissioning.helpers.definitions import definitions
from BTVNanoCommissioning.helpers.cTagSFReader import getSF

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
        lmupt_axis   = hist.Bin("lmupt",   r"Muon pt", 45, 20, 200)
        met_axis = hist.Bin("met",   r"Missing ET", 50, 0, 500)
        npv_axis = hist.Bin("npv", r"N-primary vertices",  40,  0,80) 
        nsv_axis = hist.Bin("nsv", r"N-secondary vertices",  10,  0,20) 

        # Jet
        jet_pt_axis   = hist.Bin("pt",   r"Jet $p_{T}$ [GeV]", 50, 0, 500)
        jet_ptwide_axis   = hist.Bin("ptwide",   r"Jet $p_{T}$ [GeV]", [25, 30, 40, 60, 80, 100,150,200,300,500])
        jet_eta_axis  = hist.Bin("eta",  r"Jet $\eta$", 25, -2.5, 2.5)
        jet_etawide_axis  = hist.Bin("etawide",  r"Jet $\eta$", [-2.5,-2.0,-1.5,-0.5,0.,0.5,1.5,2.0,2.5])
        jet_phi_axis  = hist.Bin("phi",  r"Jet $\phi$", 30, -3, 3)
        jet_mass_axis = hist.Bin("mass", r"Jet $m$ [GeV]", 50, 0, 500)
        jet_dr_axis = hist.Bin("dr", r"Jet $\Delta$R(l,j)", 20, 0, 5)
        ljeta_axis = hist.Bin("ljeta",  r"Leading Jet $\eta$", 25, -2.5, 2.5)
        sljeta_axis = hist.Bin("sljeta",  r"Leading Jet $\eta$", 25, -2.5, 2.5)
        ssljeta_axis = hist.Bin("ssljeta",  r"Leading Jet $\eta$", 25, -2.5, 2.5)
        sssljeta_axis = hist.Bin("sssljeta",  r"Leading Jet $\eta$", 25, -2.5, 2.5)
        ljpt_axis     = hist.Bin("ljpt", r"Leading jet $p_{T}$ [GeV]", 50, 0, 500)
        sljpt_axis     = hist.Bin("sljpt", r"Subleading jet $p_{T}$ [GeV]", 50, 0, 500)        
        ssljpt_axis     = hist.Bin("ssljpt", r"Subleading jet $p_{T}$ [GeV]", 50, 0, 500)        
        sssljpt_axis     = hist.Bin("sssljpt", r"Subleading jet $p_{T}$ [GeV]", 50, 0, 500)        
        ljdr_axis = hist.Bin("ljdr", "Leading jet $\Delta$R(l,j)", 20,0,5)
        sljdr_axis = hist.Bin("sljdr", "Subleading jet $\Delta$R(l,j)", 20,0,5)
        ljdr_axis = hist.Bin("ljdr", "Leading jet $\Delta$R(l,j)", 25,0,5)
        sljdr_axis = hist.Bin("sljdr", "Subleading jet $\Delta$R(l,j)", 25,0,5)
        ssljdr_axis = hist.Bin("ssljdr", "ssSubleading jet $\Delta$R(l,j)", 25,0,5)
        sssljdr_axis = hist.Bin("sssljdr", "sssSubleading jet $\Delta$R(l,j)", 25,0,5)

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
        _hist_deepcsv_dict={}
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
        for disc, axis in zip(varlist, btag_axes):
            for i in range(4):
                # _hist_sf_dict["%s_%d" %(disc,i)] = hist.Hist("Counts", dataset_axis, flav_axis, jet_etawide_axis, jet_ptwide_axis,axis)
                _hist_sf_dict["%s_%d" %(disc,i)] = hist.Hist("Counts", dataset_axis, flav_axis,axis)
        for deepcsv, axises in zip(deepcsv_list, deepcsv_axes):
             _hist_deepcsv_dict[deepcsv] = hist.Hist("Counts", dataset_axis,flav_axis,  axises)
        
    
        # Generate some histograms dynamically
        _hist_event_dict = {
                'njet'  : hist.Hist("Counts", dataset_axis,  njet_axis),
                'nbjet' : hist.Hist("Counts", dataset_axis,  nbjet_axis),
                'nbjet_up' : hist.Hist("Counts", dataset_axis,  nbjet_axis),
                'nbjet_dn' : hist.Hist("Counts", dataset_axis,  nbjet_axis),
                'ntbjet' : hist.Hist("Counts", dataset_axis,  ntbjet_axis),
                'nmu'   : hist.Hist("Counts", dataset_axis,  nmu_axis),
                'lmupt' : hist.Hist("Counts", dataset_axis,  lmupt_axis),
                'ljeta' : hist.Hist("Counts", dataset_axis, flav_axis, ljeta_axis),
                'sljeta' : hist.Hist("Counts", dataset_axis, flav_axis, sljeta_axis),
                'ssljeta' : hist.Hist("Counts", dataset_axis, flav_axis, ssljeta_axis),
                'sssljeta' : hist.Hist("Counts", dataset_axis, flav_axis, sssljeta_axis),
                'ljpt'  : hist.Hist("Counts", dataset_axis, flav_axis, ljpt_axis),
                'sljpt'  : hist.Hist("Counts", dataset_axis, flav_axis, sljpt_axis),
                'ssljpt'  : hist.Hist("Counts", dataset_axis, flav_axis, ssljpt_axis),
                'sssljpt'  : hist.Hist("Counts", dataset_axis, flav_axis, sssljpt_axis),
                'ljdr'  : hist.Hist("Counts", dataset_axis, flav_axis, ljdr_axis),
                'sljdr'  : hist.Hist("Counts", dataset_axis, flav_axis, sljdr_axis),
                'ssljdr'  : hist.Hist("Counts", dataset_axis, flav_axis, ssljdr_axis),
                'sssljdr'  : hist.Hist("Counts", dataset_axis, flav_axis, sssljdr_axis),
                'met' : hist.Hist("Counts", dataset_axis, met_axis)
            }
        self.sf_dict = list(_hist_sf_dict.keys())
        self.event_hists = list(_hist_event_dict.keys())
        self.deepcsv_hists = list(_hist_deepcsv_dict.keys())
        _hist_dict = {**_hist_sf_dict,**_hist_event_dict,**_hist_deepcsv_dict}

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
        events.Muon = events.Muon[(events.Muon.pt > 30) & (abs(events.Muon.eta) < 2.4) & (events.Muon.tightId > .5)&(events.Muon.pfRelIso04_all<0.12)] #
        

        event_muon = ak.pad_none(events.Muon, 1, axis=1) 
        
        req_muon =(ak.count(event_muon.pt, axis=1) == 1)
        ## Jet cuts 
        
        event_jet = events.Jet[(events.Jet.pt> 25) & (abs(events.Jet.eta) <= 2.4)&(events.Jet.puId > 0)&(events.Jet.jetId > 0)&(ak.all(events.Jet.metric_table(events.Muon) > 0.4, axis=2))]
        # print(events.Jet.muonIdx1)                                                                               
        
        req_jets = (ak.num(event_jet)>=4)
        if  hasattr(events, "METFixEE2017"):
            req_MET = events.METFixEE2017.pt>50
        else: req_MET = events.MET.pt>50
        
        event_level = req_trig  & req_jets & req_muon&req_MET &req_lumi
        #event_level = req_trig  & req_jets & req_muon&req_MET
        # Selected
        selev = events[event_level]    
        if  hasattr(events, "METFixEE2017"): 
            MET=selev.METFixEE2017.pt
        else:MET = selev.MET.pt
        #########
        
        # Per muon
        mu_eta   = (abs(selev.Muon.eta) < 2.4)
        mu_pt    = selev.Muon.pt > 30
        mu_idiso = (selev.Muon.tightId > .5)&(selev.Muon.pfRelIso04_all<0.12)
        mu_level = mu_eta & mu_pt & mu_idiso
        
        smu=selev.Muon[mu_level]
        # Per jet : https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID
        jet_eta    = (abs(selev.Jet.eta) <= 2.4)
        jet_pt     = selev.Jet.pt > 25 
        jet_pu     = (selev.Jet.puId > 0) &(selev.Jet.jetId>0)
        jet_dr     = (ak.all(selev.Jet.metric_table(smu) > 0.4, axis=2))

        
        jet_level  = jet_pu & jet_eta & jet_pt & jet_dr
        
        # b-tag twiki : https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
        bjet_disc  = selev.Jet.btagDeepB > 0.4941 # L=0.0494, M=0.2770, T=0.7264
        bjet_level = jet_level & bjet_disc
        
        
        sjets  = selev.Jet[jet_level]
        sel_jets = sjets
        sjets = sjets[:,:4]
        sbjets = selev.Jet[bjet_level]
        
        if not isRealData: stbjets = sbjets[sbjets.hadronFlavour==5]
        if isRealData :
            genflavor = ak.zeros_like(sjets.pt)
        else:
            par_flav = (sjets.partonFlavour == 0 ) & (sjets.hadronFlavour==0)
            genflavor = sjets.hadronFlavour + 1*par_flav 
            jetsfs_c_lj = getSF(ak.to_numpy(sjets[:,0].hadronFlavour),ak.to_numpy(sjets[:,0].btagDeepFlavCvL),ak.to_numpy(sjets[:,0].btagDeepFlavCvB),deepjetc_sf)
            jetsfs_c_sj = getSF(ak.to_numpy(sjets[:,1].hadronFlavour),ak.to_numpy(sjets[:,1].btagDeepFlavCvL),ak.to_numpy(sjets[:,1].btagDeepFlavCvB),deepjetc_sf)
            jetsfs_c_tj = getSF(ak.to_numpy(sjets[:,2].hadronFlavour),ak.to_numpy(sjets[:,2].btagDeepFlavCvL),ak.to_numpy(sjets[:,2].btagDeepFlavCvB),deepjetc_sf)
            jetsfs_c_fj = getSF(ak.to_numpy(sjets[:,3].hadronFlavour),ak.to_numpy(sjets[:,3].btagDeepFlavCvL),ak.to_numpy(sjets[:,3].btagDeepFlavCvB),deepjetc_sf)
            csvsfs_c_lj = getSF(ak.to_numpy(sjets[:,0].hadronFlavour),ak.to_numpy(sjets[:,0].btagDeepCvL),ak.to_numpy(sjets[:,0].btagDeepCvB),deepcsvc_sf)
            csvsfs_c_sj = getSF(ak.to_numpy(sjets[:,1].hadronFlavour),ak.to_numpy(sjets[:,1].btagDeepCvL),ak.to_numpy(sjets[:,1].btagDeepCvB),deepcsvc_sf)
            csvsfs_c_tj = getSF(ak.to_numpy(sjets[:,2].hadronFlavour),ak.to_numpy(sjets[:,2].btagDeepCvL),ak.to_numpy(sjets[:,2].btagDeepCvB),deepcsvc_sf)
            csvsfs_c_fj = getSF(ak.to_numpy(sjets[:,3].hadronFlavour),ak.to_numpy(sjets[:,3].btagDeepCvL),ak.to_numpy(sjets[:,3].btagDeepCvB),deepcsvc_sf)
            jetsfs_c_up_lj = getSF(ak.to_numpy(sjets[:,0].hadronFlavour),ak.to_numpy(sjets[:,0].btagDeepFlavCvL),ak.to_numpy(sjets[:,0].btagDeepFlavCvB),deepjetc_sf,"TotalUncUp")
            jetsfs_c_up_sj = getSF(ak.to_numpy(sjets[:,1].hadronFlavour),ak.to_numpy(sjets[:,1].btagDeepFlavCvL),ak.to_numpy(sjets[:,1].btagDeepFlavCvB),deepjetc_sf,"TotalUncUp")
            jetsfs_c_up_tj = getSF(ak.to_numpy(sjets[:,2].hadronFlavour),ak.to_numpy(sjets[:,2].btagDeepFlavCvL),ak.to_numpy(sjets[:,2].btagDeepFlavCvB),deepjetc_sf,"TotalUncUp")
            jetsfs_c_up_fj = getSF(ak.to_numpy(sjets[:,3].hadronFlavour),ak.to_numpy(sjets[:,3].btagDeepFlavCvL),ak.to_numpy(sjets[:,3].btagDeepFlavCvB),deepjetc_sf,"TotalUncUp")
            csvsfs_c_up_lj = getSF(ak.to_numpy(sjets[:,0].hadronFlavour),ak.to_numpy(sjets[:,0].btagDeepCvL),ak.to_numpy(sjets[:,0].btagDeepCvB),deepcsvc_sf,"TotalUncUp")
            csvsfs_c_up_sj = getSF(ak.to_numpy(sjets[:,1].hadronFlavour),ak.to_numpy(sjets[:,1].btagDeepCvL),ak.to_numpy(sjets[:,1].btagDeepCvB),deepcsvc_sf,"TotalUncUp")
            csvsfs_c_up_tj = getSF(ak.to_numpy(sjets[:,2].hadronFlavour),ak.to_numpy(sjets[:,2].btagDeepCvL),ak.to_numpy(sjets[:,2].btagDeepCvB),deepcsvc_sf,"TotalUncUp")
            csvsfs_c_up_fj = getSF(ak.to_numpy(sjets[:,3].hadronFlavour),ak.to_numpy(sjets[:,3].btagDeepCvL),ak.to_numpy(sjets[:,3].btagDeepCvB),deepcsvc_sf,"TotalUncUp")
            jetsfs_c_dn_lj = getSF(ak.to_numpy(sjets[:,0].hadronFlavour),ak.to_numpy(sjets[:,0].btagDeepFlavCvL),ak.to_numpy(sjets[:,0].btagDeepFlavCvB),deepjetc_sf,"TotalUncDown")
            jetsfs_c_dn_sj = getSF(ak.to_numpy(sjets[:,1].hadronFlavour),ak.to_numpy(sjets[:,1].btagDeepFlavCvL),ak.to_numpy(sjets[:,1].btagDeepFlavCvB),deepjetc_sf,"TotalUncDown")
            jetsfs_c_dn_tj = getSF(ak.to_numpy(sjets[:,2].hadronFlavour),ak.to_numpy(sjets[:,2].btagDeepFlavCvL),ak.to_numpy(sjets[:,2].btagDeepFlavCvB),deepjetc_sf,"TotalUncDown")
            jetsfs_c_dn_fj = getSF(ak.to_numpy(sjets[:,3].hadronFlavour),ak.to_numpy(sjets[:,3].btagDeepFlavCvL),ak.to_numpy(sjets[:,3].btagDeepFlavCvB),deepjetc_sf,"TotalUncDown")
            csvsfs_c_dn_lj = getSF(ak.to_numpy(sjets[:,0].hadronFlavour),ak.to_numpy(sjets[:,0].btagDeepCvL),ak.to_numpy(sjets[:,0].btagDeepCvB),deepcsvc_sf,"TotalUncDown")
            csvsfs_c_dn_sj = getSF(ak.to_numpy(sjets[:,1].hadronFlavour),ak.to_numpy(sjets[:,1].btagDeepCvL),ak.to_numpy(sjets[:,1].btagDeepCvB),deepcsvc_sf,"TotalUncDown")
            csvsfs_c_dn_tj = getSF(ak.to_numpy(sjets[:,2].hadronFlavour),ak.to_numpy(sjets[:,2].btagDeepCvL),ak.to_numpy(sjets[:,2].btagDeepCvB),deepcsvc_sf,"TotalUncDown")
            csvsfs_c_dn_fj = getSF(ak.to_numpy(sjets[:,3].hadronFlavour),ak.to_numpy(sjets[:,3].btagDeepCvL),ak.to_numpy(sjets[:,3].btagDeepCvB),deepcsvc_sf,"TotalUncDown")
            jetsfs_b_lj = ak.to_numpy(deepjetb_sf.eval('central',sjets[:,0].hadronFlavour,abs(sjets[:,0].eta),sjets[:,0].pt,discr=sjets[:,0].btagDeepFlavB))
            jetsfs_b_sj = ak.to_numpy(deepjetb_sf.eval('central',sjets[:,1].hadronFlavour,abs(sjets[:,1].eta),sjets[:,1].pt,discr=sjets[:,1].btagDeepFlavB))
            jetsfs_b_tj = ak.to_numpy(deepjetb_sf.eval('central',sjets[:,2].hadronFlavour,abs(sjets[:,2].eta),sjets[:,2].pt,discr=sjets[:,2].btagDeepFlavB))
            jetsfs_b_fj = ak.to_numpy(deepjetb_sf.eval('central',sjets[:,3].hadronFlavour,abs(sjets[:,3].eta),sjets[:,3].pt,discr=sjets[:,3].btagDeepFlavB))
            csvsfs_b_lj = ak.to_numpy(deepcsvb_sf.eval('central',sjets[:,0].hadronFlavour,abs(sjets[:,0].eta),sjets[:,0].pt,discr=sjets[:,0].btagDeepB))
            csvsfs_b_sj = ak.to_numpy(deepcsvb_sf.eval('central',sjets[:,1].hadronFlavour,abs(sjets[:,1].eta),sjets[:,1].pt,discr=sjets[:,1].btagDeepB))
            csvsfs_b_tj = ak.to_numpy(deepcsvb_sf.eval('central',sjets[:,2].hadronFlavour,abs(sjets[:,2].eta),sjets[:,2].pt,discr=sjets[:,2].btagDeepB))
            csvsfs_b_fj = ak.to_numpy(deepcsvb_sf.eval('central',sjets[:,3].hadronFlavour,abs(sjets[:,3].eta),sjets[:,3].pt,discr=sjets[:,3].btagDeepB))
            jetsfs_b_up_lj = ak.to_numpy(deepjetb_sf.eval('up_jes',sjets[:,0].hadronFlavour,abs(sjets[:,0].eta),sjets[:,0].pt,discr=sjets[:,0].btagDeepFlavB))
            jetsfs_b_up_sj = ak.to_numpy(deepjetb_sf.eval('up_jes',sjets[:,1].hadronFlavour,abs(sjets[:,1].eta),sjets[:,1].pt,discr=sjets[:,1].btagDeepFlavB))
            jetsfs_b_up_tj = ak.to_numpy(deepjetb_sf.eval('up_jes',sjets[:,2].hadronFlavour,abs(sjets[:,2].eta),sjets[:,2].pt,discr=sjets[:,2].btagDeepFlavB))
            jetsfs_b_up_fj = ak.to_numpy(deepjetb_sf.eval('up_jes',sjets[:,3].hadronFlavour,abs(sjets[:,3].eta),sjets[:,3].pt,discr=sjets[:,3].btagDeepFlavB))
            csvsfs_b_up_lj = ak.to_numpy(deepcsvb_sf.eval('up_jes',sjets[:,0].hadronFlavour,abs(sjets[:,0].eta),sjets[:,0].pt,discr=sjets[:,0].btagDeepB))
            csvsfs_b_up_sj = ak.to_numpy(deepcsvb_sf.eval('up_jes',sjets[:,1].hadronFlavour,abs(sjets[:,1].eta),sjets[:,1].pt,discr=sjets[:,1].btagDeepB))
            csvsfs_b_up_tj = ak.to_numpy(deepcsvb_sf.eval('up_jes',sjets[:,2].hadronFlavour,abs(sjets[:,2].eta),sjets[:,2].pt,discr=sjets[:,2].btagDeepB))
            csvsfs_b_up_fj = ak.to_numpy(deepcsvb_sf.eval('up_jes',sjets[:,3].hadronFlavour,abs(sjets[:,3].eta),sjets[:,3].pt,discr=sjets[:,3].btagDeepB))
            jetsfs_b_dn_lj = ak.to_numpy(deepjetb_sf.eval('down_jes',sjets[:,0].hadronFlavour,abs(sjets[:,0].eta),sjets[:,0].pt,discr=sjets[:,0].btagDeepFlavB))
            jetsfs_b_dn_sj = ak.to_numpy(deepjetb_sf.eval('down_jes',sjets[:,1].hadronFlavour,abs(sjets[:,1].eta),sjets[:,1].pt,discr=sjets[:,1].btagDeepFlavB))
            jetsfs_b_dn_tj = ak.to_numpy(deepjetb_sf.eval('down_jes',sjets[:,2].hadronFlavour,abs(sjets[:,2].eta),sjets[:,2].pt,discr=sjets[:,2].btagDeepFlavB))
            jetsfs_b_dn_fj = ak.to_numpy(deepjetb_sf.eval('down_jes',sjets[:,3].hadronFlavour,abs(sjets[:,3].eta),sjets[:,3].pt,discr=sjets[:,3].btagDeepFlavB))
            csvsfs_b_dn_lj = ak.to_numpy(deepcsvb_sf.eval('down_jes',sjets[:,0].hadronFlavour,abs(sjets[:,0].eta),sjets[:,0].pt,discr=sjets[:,0].btagDeepB))
            csvsfs_b_dn_sj = ak.to_numpy(deepcsvb_sf.eval('down_jes',sjets[:,1].hadronFlavour,abs(sjets[:,1].eta),sjets[:,1].pt,discr=sjets[:,1].btagDeepB))
            csvsfs_b_dn_tj = ak.to_numpy(deepcsvb_sf.eval('down_jes',sjets[:,2].hadronFlavour,abs(sjets[:,2].eta),sjets[:,2].pt,discr=sjets[:,2].btagDeepB))
            csvsfs_b_dn_fj = ak.to_numpy(deepcsvb_sf.eval('down_jes',sjets[:,3].hadronFlavour,abs(sjets[:,3].eta),sjets[:,3].pt,discr=sjets[:,3].btagDeepB))

            
                
            #print(ak.type(selev))
            #print(ak.type(selev.PV.npvs))
            # # Fill histograms dynamically  
        # for histname, h in output.items():
        #     if histname in self.deepcsv_hists:
        #         if(isRealData):
        #             fields = {l: ak.flatten(sjets[l], axis=None) for l in h.fields if l in dir(sjets)}
        #             # h.fill(dataset=dataset,flav=5,npv=ak.flatten(ak.broadcast_arrays(selev.PV.npvs,sjets['pt'])[0]), **fields)
        #             h.fill(dataset=dataset,flav=5, **fields)
        #         else:
        #             fields = {l: ak.flatten(sjets[histname]) for l in h.fields if l in dir(sjets)}
        #             genweiev=ak.flatten(ak.broadcast_arrays(weights.weight()[event_level],sjets['pt'])[0])
        #             h.fill(dataset=dataset,flav=ak.flatten(genflavor),**fields,weight=genweiev)
       
            
        
              
        if not isRealData:
            output['btagDeepFlavB_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepFlavB=sjets[:,0].btagDeepFlavB,weight=weights.weight()[event_level])
            output['btagDeepFlavC_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepFlavC=sjets[:,0].btagDeepFlavC,weight=weights.weight()[event_level])
            output['btagDeepB_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepB=sjets[:,0].btagDeepB,weight=weights.weight()[event_level])
            output['btagDeepC_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepC=sjets[:,0].btagDeepC,weight=weights.weight()[event_level])
            output['deepcsv_CvB_0'].fill(dataset=dataset,flav=genflavor[:,0], deepcsv_CvB=np.where(sjets[:,0].btagDeepCvB<0,-0.2,sjets[:,0].btagDeepCvB),weight=weights.weight()[event_level])
            output['deepcsv_CvL_0'].fill(dataset=dataset,flav=genflavor[:,0], deepcsv_CvL=np.where(sjets[:,0].btagDeepCvL<0,-0.2,sjets[:,0].btagDeepCvL),weight=weights.weight()[event_level])
            output['deepflav_CvB_0'].fill(dataset=dataset,flav=genflavor[:,0], deepflav_CvB=np.where(sjets[:,0].btagDeepFlavCvB<0,-0.2,sjets[:,0].btagDeepFlavCvB),weight=weights.weight()[event_level])
            output['deepflav_CvL_0'].fill(dataset=dataset,flav=genflavor[:,0], deepflav_CvL=np.where(sjets[:,0].btagDeepFlavCvL<0,-0.2,sjets[:,0].btagDeepFlavCvL),weight=weights.weight()[event_level])
            output['btagDeepFlavB_1'].fill(dataset=dataset,flav=genflavor[:,1],  btagDeepFlavB=sjets[:,1].btagDeepFlavB,weight=weights.weight()[event_level])
            output['btagDeepFlavC_1'].fill(dataset=dataset,flav=genflavor[:,1],  btagDeepFlavC=sjets[:,1].btagDeepFlavC,weight=weights.weight()[event_level])
            output['btagDeepB_1'].fill(dataset=dataset,flav=genflavor[:,1],  btagDeepB=sjets[:,1].btagDeepB,weight=weights.weight()[event_level])
            output['btagDeepC_1'].fill(dataset=dataset,flav=genflavor[:,1],  btagDeepC=sjets[:,1].btagDeepC,weight=weights.weight()[event_level])
            output['deepcsv_CvB_1'].fill(dataset=dataset,flav=genflavor[:,1],  deepcsv_CvB=np.where(sjets[:,1].btagDeepCvB<0,-0.2,sjets[:,1].btagDeepCvB),weight=weights.weight()[event_level])
            output['deepcsv_CvL_1'].fill(dataset=dataset,flav=genflavor[:,1],  deepcsv_CvL=np.where(sjets[:,1].btagDeepCvL<0,-0.2,sjets[:,1].btagDeepCvL),weight=weights.weight()[event_level])
            output['deepflav_CvB_1'].fill(dataset=dataset,flav=genflavor[:,1],  deepflav_CvB=np.where(sjets[:,1].btagDeepFlavCvB<0,-0.2,sjets[:,1].btagDeepFlavCvB),weight=weights.weight()[event_level])
            output['deepflav_CvL_1'].fill(dataset=dataset,flav=genflavor[:,1],  deepflav_CvL=np.where(sjets[:,1].btagDeepFlavCvL<0,-0.2,sjets[:,1].btagDeepFlavCvL),weight=weights.weight()[event_level])
            output['btagDeepFlavB_2'].fill(dataset=dataset,flav=genflavor[:,2],  btagDeepFlavB=sjets[:,2].btagDeepFlavB,weight=weights.weight()[event_level])
            output['btagDeepFlavC_2'].fill(dataset=dataset,flav=genflavor[:,2],  btagDeepFlavC=sjets[:,2].btagDeepFlavC,weight=weights.weight()[event_level])
            output['btagDeepB_2'].fill(dataset=dataset,flav=genflavor[:,2],  btagDeepB=sjets[:,2].btagDeepB,weight=weights.weight()[event_level])
            output['btagDeepC_2'].fill(dataset=dataset,flav=genflavor[:,2],  btagDeepC=sjets[:,2].btagDeepC,weight=weights.weight()[event_level])
            output['deepcsv_CvB_2'].fill(dataset=dataset,flav=genflavor[:,2],  deepcsv_CvB=np.where(sjets[:,2].btagDeepCvB<0,-0.2,sjets[:,2].btagDeepCvB),weight=weights.weight()[event_level])
            output['deepcsv_CvL_2'].fill(dataset=dataset,flav=genflavor[:,2],  deepcsv_CvL=np.where(sjets[:,2].btagDeepCvL<0,-0.2,sjets[:,2].btagDeepCvL),weight=weights.weight()[event_level])
            output['deepflav_CvB_2'].fill(dataset=dataset,flav=genflavor[:,2],  deepflav_CvB=np.where(sjets[:,2].btagDeepFlavCvB<0,-0.2,sjets[:,2].btagDeepFlavCvB),weight=weights.weight()[event_level])
            output['deepflav_CvL_2'].fill(dataset=dataset,flav=genflavor[:,2],  deepflav_CvL=np.where(sjets[:,2].btagDeepFlavCvL<0,-0.2,sjets[:,2].btagDeepFlavCvL),weight=weights.weight()[event_level])
            output['btagDeepFlavB_3'].fill(dataset=dataset,flav=genflavor[:,3],  btagDeepFlavB=sjets[:,3].btagDeepFlavB,weight=weights.weight()[event_level])
            output['btagDeepFlavC_3'].fill(dataset=dataset,flav=genflavor[:,3],  btagDeepFlavC=sjets[:,3].btagDeepFlavC,weight=weights.weight()[event_level])
            output['btagDeepB_3'].fill(dataset=dataset,flav=genflavor[:,3],  btagDeepB=sjets[:,3].btagDeepB,weight=weights.weight()[event_level])
            output['btagDeepC_3'].fill(dataset=dataset,flav=genflavor[:,3],  btagDeepC=sjets[:,3].btagDeepC,weight=weights.weight()[event_level])
            output['deepcsv_CvB_3'].fill(dataset=dataset,flav=genflavor[:,3],  deepcsv_CvB=np.where(sjets[:,3].btagDeepCvB<0,-0.2,sjets[:,3].btagDeepCvB),weight=weights.weight()[event_level])
            output['deepcsv_CvL_3'].fill(dataset=dataset,flav=genflavor[:,3],  deepcsv_CvL=np.where(sjets[:,3].btagDeepCvL<0,-0.2,sjets[:,3].btagDeepCvL),weight=weights.weight()[event_level])
            output['deepflav_CvB_3'].fill(dataset=dataset,flav=genflavor[:,3],  deepflav_CvB=np.where(sjets[:,3].btagDeepFlavCvB<0,-0.2,sjets[:,3].btagDeepFlavCvB),weight=weights.weight()[event_level])
            output['deepflav_CvL_3'].fill(dataset=dataset,flav=genflavor[:,3],  deepflav_CvL=np.where(sjets[:,3].btagDeepFlavCvL<0,-0.2,sjets[:,3].btagDeepFlavCvL),weight=weights.weight()[event_level])
            
            ##SFs 
            output['btagDeepFlavBSF_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepFlavBSF=sjets[:,0].btagDeepFlavB,weight=(weights.weight()[event_level]*jetsfs_b_lj))      
            output['btagDeepFlavCSF_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepFlavCSF=sjets[:,0].btagDeepFlavC,weight=weights.weight()[event_level]*jetsfs_c_lj)
            output['btagDeepBSF_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepBSF=sjets[:,0].btagDeepB,weight=weights.weight()[event_level]*csvsfs_b_lj)
            output['btagDeepCSF_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepCSF=sjets[:,0].btagDeepC,weight=weights.weight()[event_level]*csvsfs_c_lj)
            output['deepcsv_CvBSF_0'].fill(dataset=dataset,flav=genflavor[:,0], deepcsv_CvBSF=np.where(sjets[:,0].btagDeepCvB<0,-0.2,sjets[:,0].btagDeepCvB),weight=weights.weight()[event_level]*csvsfs_c_lj)
            output['deepcsv_CvLSF_0'].fill(dataset=dataset,flav=genflavor[:,0], deepcsv_CvLSF=np.where(sjets[:,0].btagDeepCvL<0,-0.2,sjets[:,0].btagDeepCvL),weight=weights.weight()[event_level]*csvsfs_c_lj)
            output['deepflav_CvBSF_0'].fill(dataset=dataset,flav=genflavor[:,0], deepflav_CvBSF=np.where(sjets[:,0].btagDeepFlavCvB<0,-0.2,sjets[:,0].btagDeepFlavCvB),weight=weights.weight()[event_level]*jetsfs_c_lj)
            output['deepflav_CvLSF_0'].fill(dataset=dataset,flav=genflavor[:,0], deepflav_CvLSF=np.where(sjets[:,0].btagDeepFlavCvL<0,-0.2,sjets[:,0].btagDeepFlavCvL),weight=weights.weight()[event_level]*jetsfs_c_lj)
            output['btagDeepFlavBSF_1'].fill(dataset=dataset,flav=genflavor[:,1],  btagDeepFlavBSF=sjets[:,1].btagDeepFlavB,weight=weights.weight()[event_level]*jetsfs_b_sj)
            output['btagDeepFlavCSF_1'].fill(dataset=dataset,flav=genflavor[:,1],  btagDeepFlavCSF=sjets[:,1].btagDeepFlavC,weight=weights.weight()[event_level]*jetsfs_c_sj)
            output['btagDeepBSF_1'].fill(dataset=dataset,flav=genflavor[:,1],  btagDeepBSF=sjets[:,1].btagDeepB,weight=weights.weight()[event_level]*csvsfs_b_sj)
            output['btagDeepCSF_1'].fill(dataset=dataset,flav=genflavor[:,1],  btagDeepCSF=sjets[:,1].btagDeepC,weight=weights.weight()[event_level]*csvsfs_c_sj)
            output['deepcsv_CvBSF_1'].fill(dataset=dataset,flav=genflavor[:,1],  deepcsv_CvBSF=np.where(sjets[:,1].btagDeepCvB<0,-0.2,sjets[:,1].btagDeepCvB),weight=weights.weight()[event_level]*csvsfs_c_sj)
            output['deepcsv_CvLSF_1'].fill(dataset=dataset,flav=genflavor[:,1],  deepcsv_CvLSF=np.where(sjets[:,1].btagDeepCvL<0,-0.2,sjets[:,1].btagDeepCvL),weight=weights.weight()[event_level]*csvsfs_c_sj)
            output['deepflav_CvBSF_1'].fill(dataset=dataset,flav=genflavor[:,1],  deepflav_CvBSF=np.where(sjets[:,1].btagDeepFlavCvB<0,-0.2,sjets[:,1].btagDeepFlavCvB),weight=weights.weight()[event_level]*jetsfs_c_sj)
            output['deepflav_CvLSF_1'].fill(dataset=dataset,flav=genflavor[:,1],  deepflav_CvLSF=np.where(sjets[:,1].btagDeepFlavCvL<0,-0.2,sjets[:,1].btagDeepFlavCvL),weight=weights.weight()[event_level]*jetsfs_c_sj)
            output['btagDeepFlavBSF_2'].fill(dataset=dataset,flav=genflavor[:,2],  btagDeepFlavBSF=sjets[:,2].btagDeepFlavB,weight=(weights.weight()[event_level]*jetsfs_b_tj))      
            output['btagDeepFlavCSF_2'].fill(dataset=dataset,flav=genflavor[:,2],  btagDeepFlavCSF=sjets[:,2].btagDeepFlavC,weight=weights.weight()[event_level]*jetsfs_c_tj)
            output['btagDeepBSF_2'].fill(dataset=dataset,flav=genflavor[:,2],  btagDeepBSF=sjets[:,2].btagDeepB,weight=weights.weight()[event_level]*csvsfs_b_tj)
            output['btagDeepCSF_2'].fill(dataset=dataset,flav=genflavor[:,2],  btagDeepCSF=sjets[:,2].btagDeepC,weight=weights.weight()[event_level]*csvsfs_c_tj)
            output['deepcsv_CvBSF_2'].fill(dataset=dataset,flav=genflavor[:,2],  deepcsv_CvBSF=np.where(sjets[:,2].btagDeepCvB<0,-0.2,sjets[:,2].btagDeepCvB),weight=weights.weight()[event_level]*csvsfs_c_tj)
            output['deepcsv_CvLSF_2'].fill(dataset=dataset,flav=genflavor[:,2],  deepcsv_CvLSF=np.where(sjets[:,2].btagDeepCvL<0,-0.2,sjets[:,2].btagDeepCvL),weight=weights.weight()[event_level]*csvsfs_c_tj)
            output['deepflav_CvBSF_2'].fill(dataset=dataset,flav=genflavor[:,2],  deepflav_CvBSF=np.where(sjets[:,2].btagDeepFlavCvB<0,-0.2,sjets[:,2].btagDeepFlavCvB),weight=weights.weight()[event_level]*jetsfs_c_tj)
            output['deepflav_CvLSF_2'].fill(dataset=dataset,flav=genflavor[:,2],  deepflav_CvLSF=np.where(sjets[:,2].btagDeepFlavCvL<0,-0.2,sjets[:,2].btagDeepFlavCvL),weight=weights.weight()[event_level]*jetsfs_c_tj)
            output['btagDeepFlavBSF_3'].fill(dataset=dataset,flav=genflavor[:,3],  btagDeepFlavBSF=sjets[:,3].btagDeepFlavB,weight=(weights.weight()[event_level]*jetsfs_b_fj))      
            output['btagDeepFlavCSF_3'].fill(dataset=dataset,flav=genflavor[:,3],  btagDeepFlavCSF=sjets[:,3].btagDeepFlavC,weight=weights.weight()[event_level]*jetsfs_c_fj)
            output['btagDeepBSF_3'].fill(dataset=dataset,flav=genflavor[:,3],  btagDeepBSF=sjets[:,3].btagDeepB,weight=weights.weight()[event_level]*csvsfs_b_fj)
            output['btagDeepCSF_3'].fill(dataset=dataset,flav=genflavor[:,3],  btagDeepCSF=sjets[:,3].btagDeepC,weight=weights.weight()[event_level]*csvsfs_c_fj)
            output['deepcsv_CvBSF_3'].fill(dataset=dataset,flav=genflavor[:,3],  deepcsv_CvBSF=np.where(sjets[:,3].btagDeepCvB<0,-0.2,sjets[:,3].btagDeepCvB),weight=weights.weight()[event_level]*csvsfs_c_fj)
            output['deepcsv_CvLSF_3'].fill(dataset=dataset,flav=genflavor[:,3],  deepcsv_CvLSF=np.where(sjets[:,3].btagDeepCvL<0,-0.2,sjets[:,3].btagDeepCvL),weight=weights.weight()[event_level]*csvsfs_c_fj)
            output['deepflav_CvBSF_3'].fill(dataset=dataset,flav=genflavor[:,3],  deepflav_CvBSF=np.where(sjets[:,3].btagDeepFlavCvB<0,-0.2,sjets[:,3].btagDeepFlavCvB),weight=weights.weight()[event_level]*jetsfs_c_fj)
            output['deepflav_CvLSF_3'].fill(dataset=dataset,flav=genflavor[:,3],  deepflav_CvLSF=np.where(sjets[:,3].btagDeepFlavCvL<0,-0.2,sjets[:,3].btagDeepFlavCvL),weight=weights.weight()[event_level]*jetsfs_c_fj)
            ###SFs up
            output['btagDeepFlavB_up_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepFlavB_up=sjets[:,0].btagDeepFlavB,weight=weights.weight()[event_level]*jetsfs_b_up_lj)
            output['btagDeepFlavC_up_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepFlavC_up=sjets[:,0].btagDeepFlavC,weight=weights.weight()[event_level]*jetsfs_c_up_lj)
            output['btagDeepB_up_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepB_up=sjets[:,0].btagDeepB,weight=weights.weight()[event_level]*csvsfs_b_up_lj)
            output['btagDeepC_up_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepC_up=sjets[:,0].btagDeepC,weight=weights.weight()[event_level]*csvsfs_c_up_lj)
            output['deepcsv_CvB_up_0'].fill(dataset=dataset,flav=genflavor[:,0], deepcsv_CvB_up=np.where(sjets[:,0].btagDeepCvB<0,-0.2,sjets[:,0].btagDeepCvB),weight=weights.weight()[event_level]*csvsfs_c_up_lj)
            output['deepcsv_CvL_up_0'].fill(dataset=dataset,flav=genflavor[:,0], deepcsv_CvL_up=np.where(sjets[:,0].btagDeepCvL<0,-0.2,sjets[:,0].btagDeepCvL),weight=weights.weight()[event_level]*csvsfs_c_up_lj)
            output['deepflav_CvB_up_0'].fill(dataset=dataset,flav=genflavor[:,0], deepflav_CvB_up=np.where(sjets[:,0].btagDeepFlavCvB<0,-0.2,sjets[:,0].btagDeepFlavCvB),weight=weights.weight()[event_level]*jetsfs_c_up_lj)
            output['deepflav_CvL_up_0'].fill(dataset=dataset,flav=genflavor[:,0], deepflav_CvL_up=np.where(sjets[:,0].btagDeepFlavCvL<0,-0.2,sjets[:,0].btagDeepFlavCvL),weight=weights.weight()[event_level]*jetsfs_c_up_lj)
            output['btagDeepFlavB_up_1'].fill(dataset=dataset,flav=genflavor[:,1],  btagDeepFlavB_up=sjets[:,1].btagDeepFlavB,weight=weights.weight()[event_level]*jetsfs_b_up_sj)
            output['btagDeepFlavC_up_1'].fill(dataset=dataset,flav=genflavor[:,1],  btagDeepFlavC_up=sjets[:,1].btagDeepFlavC,weight=weights.weight()[event_level]*jetsfs_c_up_sj)
            output['btagDeepB_up_1'].fill(dataset=dataset,flav=genflavor[:,1],  btagDeepB_up=sjets[:,1].btagDeepB,weight=weights.weight()[event_level]*csvsfs_b_up_sj)
            output['btagDeepC_up_1'].fill(dataset=dataset,flav=genflavor[:,1],  btagDeepC_up=sjets[:,1].btagDeepC,weight=weights.weight()[event_level]*csvsfs_c_up_sj)
            output['deepcsv_CvB_up_1'].fill(dataset=dataset,flav=genflavor[:,1],  deepcsv_CvB_up=np.where(sjets[:,1].btagDeepCvB<0,-0.2,sjets[:,1].btagDeepCvB),weight=weights.weight()[event_level]*csvsfs_c_up_sj)
            output['deepcsv_CvL_up_1'].fill(dataset=dataset,flav=genflavor[:,1],   deepcsv_CvL_up=np.where(sjets[:,1].btagDeepCvL<0,-0.2,sjets[:,1].btagDeepCvL),weight=weights.weight()[event_level]*csvsfs_c_up_sj)
            output['deepflav_CvB_up_1'].fill(dataset=dataset,flav=genflavor[:,1],  deepflav_CvB_up=np.where(sjets[:,1].btagDeepFlavCvB<0,-0.2,sjets[:,1].btagDeepFlavCvB),weight=weights.weight()[event_level]*jetsfs_c_up_sj)
            output['deepflav_CvL_up_1'].fill(dataset=dataset,flav=genflavor[:,1],  deepflav_CvL_up=np.where(sjets[:,1].btagDeepFlavCvL<0,-0.2,sjets[:,1].btagDeepFlavCvL),weight=weights.weight()[event_level]*jetsfs_c_up_sj)
            output['btagDeepFlavB_up_2'].fill(dataset=dataset,flav=genflavor[:,2],  btagDeepFlavB_up=sjets[:,2].btagDeepFlavB,weight=weights.weight()[event_level]*jetsfs_b_up_tj)
            output['btagDeepFlavC_up_2'].fill(dataset=dataset,flav=genflavor[:,2],  btagDeepFlavC_up=sjets[:,2].btagDeepFlavC,weight=weights.weight()[event_level]*jetsfs_c_up_tj)
            output['btagDeepB_up_2'].fill(dataset=dataset,flav=genflavor[:,2],  btagDeepB_up=sjets[:,2].btagDeepB,weight=weights.weight()[event_level]*csvsfs_b_up_tj)
            output['btagDeepC_up_2'].fill(dataset=dataset,flav=genflavor[:,2],  btagDeepC_up=sjets[:,2].btagDeepC,weight=weights.weight()[event_level]*csvsfs_c_up_tj)
            output['deepcsv_CvB_up_2'].fill(dataset=dataset,flav=genflavor[:,2],  deepcsv_CvB_up=np.where(sjets[:,2].btagDeepCvB<0,-0.2,sjets[:,2].btagDeepCvB),weight=weights.weight()[event_level]*csvsfs_c_up_tj)
            output['deepcsv_CvL_up_2'].fill(dataset=dataset,flav=genflavor[:,2],   deepcsv_CvL_up=np.where(sjets[:,2].btagDeepCvL<0,-0.2,sjets[:,2].btagDeepCvL),weight=weights.weight()[event_level]*csvsfs_c_up_tj)
            output['deepflav_CvB_up_2'].fill(dataset=dataset,flav=genflavor[:,2],  deepflav_CvB_up=np.where(sjets[:,2].btagDeepFlavCvB<0,-0.2,sjets[:,2].btagDeepFlavCvB),weight=weights.weight()[event_level]*jetsfs_c_up_tj)
            output['deepflav_CvL_up_2'].fill(dataset=dataset,flav=genflavor[:,2],  deepflav_CvL_up=np.where(sjets[:,2].btagDeepFlavCvL<0,-0.2,sjets[:,2].btagDeepFlavCvL),weight=weights.weight()[event_level]*jetsfs_c_up_tj)
            output['btagDeepFlavB_up_3'].fill(dataset=dataset,flav=genflavor[:,3],  btagDeepFlavB_up=sjets[:,3].btagDeepFlavB,weight=weights.weight()[event_level]*jetsfs_b_up_fj)
            output['btagDeepFlavC_up_3'].fill(dataset=dataset,flav=genflavor[:,3],  btagDeepFlavC_up=sjets[:,3].btagDeepFlavC,weight=weights.weight()[event_level]*jetsfs_c_up_fj)
            output['btagDeepB_up_3'].fill(dataset=dataset,flav=genflavor[:,3],  btagDeepB_up=sjets[:,3].btagDeepB,weight=weights.weight()[event_level]*csvsfs_b_up_fj)
            output['btagDeepC_up_3'].fill(dataset=dataset,flav=genflavor[:,3],  btagDeepC_up=sjets[:,3].btagDeepC,weight=weights.weight()[event_level]*csvsfs_c_up_fj)
            output['deepcsv_CvB_up_3'].fill(dataset=dataset,flav=genflavor[:,3],  deepcsv_CvB_up=np.where(sjets[:,3].btagDeepCvB<0,-0.2,sjets[:,3].btagDeepCvB),weight=weights.weight()[event_level]*csvsfs_c_up_fj)
            output['deepcsv_CvL_up_3'].fill(dataset=dataset,flav=genflavor[:,3],   deepcsv_CvL_up=np.where(sjets[:,3].btagDeepCvL<0,-0.2,sjets[:,3].btagDeepCvL),weight=weights.weight()[event_level]*csvsfs_c_up_fj)
            output['deepflav_CvB_up_3'].fill(dataset=dataset,flav=genflavor[:,3],  deepflav_CvB_up=np.where(sjets[:,3].btagDeepFlavCvB<0,-0.2,sjets[:,3].btagDeepFlavCvB),weight=weights.weight()[event_level]*jetsfs_c_up_fj)
            output['deepflav_CvL_up_3'].fill(dataset=dataset,flav=genflavor[:,3],  deepflav_CvL_up=np.where(sjets[:,3].btagDeepFlavCvL<0,-0.2,sjets[:,3].btagDeepFlavCvL),weight=weights.weight()[event_level]*jetsfs_c_up_fj)            
            ###SFs down
            output['btagDeepFlavB_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepFlavB_dn=sjets[:,0].btagDeepFlavB,weight=weights.weight()[event_level]*jetsfs_b_dn_lj)
            output['btagDeepFlavC_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepFlavC_dn=sjets[:,0].btagDeepFlavC,weight=weights.weight()[event_level]*jetsfs_c_dn_lj)
            output['btagDeepB_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepB_dn=sjets[:,0].btagDeepB,weight=weights.weight()[event_level]*csvsfs_b_dn_lj)
            output['btagDeepC_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepC_dn=sjets[:,0].btagDeepC,weight=weights.weight()[event_level]*csvsfs_c_dn_lj)
            output['deepcsv_CvB_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], deepcsv_CvB_dn=np.where(sjets[:,0].btagDeepCvB<0,-0.2,sjets[:,0].btagDeepCvB),weight=weights.weight()[event_level]*csvsfs_c_dn_lj)
            output['deepcsv_CvL_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], deepcsv_CvL_dn=np.where(sjets[:,0].btagDeepCvL<0,-0.2,sjets[:,0].btagDeepCvL),weight=weights.weight()[event_level]*csvsfs_c_dn_lj)
            output['deepflav_CvB_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], deepflav_CvB_dn=np.where(sjets[:,0].btagDeepFlavCvB<0,-0.2,sjets[:,0].btagDeepFlavCvB),weight=weights.weight()[event_level]*jetsfs_c_dn_lj)
            output['deepflav_CvL_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], deepflav_CvL_dn=np.where(sjets[:,0].btagDeepFlavCvL<0,-0.2,sjets[:,0].btagDeepFlavCvL),weight=weights.weight()[event_level]*jetsfs_c_dn_lj)
            output['btagDeepFlavB_dn_1'].fill(dataset=dataset,flav=genflavor[:,1],  btagDeepFlavB_dn=sjets[:,1].btagDeepFlavB,weight=weights.weight()[event_level]*jetsfs_b_dn_sj)
            output['btagDeepFlavC_dn_1'].fill(dataset=dataset,flav=genflavor[:,1],  btagDeepFlavC_dn=sjets[:,1].btagDeepFlavC,weight=weights.weight()[event_level]*jetsfs_c_dn_sj)
            output['btagDeepB_dn_1'].fill(dataset=dataset,flav=genflavor[:,1],  btagDeepB_dn=sjets[:,1].btagDeepB,weight=weights.weight()[event_level]*csvsfs_b_dn_sj)
            output['btagDeepC_dn_1'].fill(dataset=dataset,flav=genflavor[:,1],  btagDeepC_dn=sjets[:,1].btagDeepC,weight=weights.weight()[event_level]*csvsfs_c_dn_sj)
            output['deepcsv_CvB_dn_1'].fill(dataset=dataset,flav=genflavor[:,1],  deepcsv_CvB_dn=np.where(sjets[:,1].btagDeepCvB<0,-0.2,sjets[:,1].btagDeepCvB),weight=weights.weight()[event_level]*csvsfs_c_dn_sj)
            output['deepcsv_CvL_dn_1'].fill(dataset=dataset,flav=genflavor[:,1],   deepcsv_CvL_dn=np.where(sjets[:,1].btagDeepCvL<0,-0.2,sjets[:,1].btagDeepCvL),weight=weights.weight()[event_level]*csvsfs_c_dn_sj)
            output['deepflav_CvB_dn_1'].fill(dataset=dataset,flav=genflavor[:,1],  deepflav_CvB_dn=np.where(sjets[:,1].btagDeepFlavCvB<0,-0.2,sjets[:,1].btagDeepFlavCvB),weight=weights.weight()[event_level]*jetsfs_c_dn_sj)
            output['deepflav_CvL_dn_1'].fill(dataset=dataset,flav=genflavor[:,1],  deepflav_CvL_dn=np.where(sjets[:,1].btagDeepFlavCvL<0,-0.2,sjets[:,1].btagDeepFlavCvL),weight=weights.weight()[event_level]*jetsfs_c_dn_sj)
            output['btagDeepFlavB_dn_2'].fill(dataset=dataset,flav=genflavor[:,2],  btagDeepFlavB_dn=sjets[:,2].btagDeepFlavB,weight=weights.weight()[event_level]*jetsfs_b_dn_tj)
            output['btagDeepFlavC_dn_2'].fill(dataset=dataset,flav=genflavor[:,2],  btagDeepFlavC_dn=sjets[:,2].btagDeepFlavC,weight=weights.weight()[event_level]*jetsfs_c_dn_tj)
            output['btagDeepB_dn_2'].fill(dataset=dataset,flav=genflavor[:,2],  btagDeepB_dn=sjets[:,2].btagDeepB,weight=weights.weight()[event_level]*csvsfs_b_dn_tj)
            output['btagDeepC_dn_2'].fill(dataset=dataset,flav=genflavor[:,2],  btagDeepC_dn=sjets[:,2].btagDeepC,weight=weights.weight()[event_level]*csvsfs_c_dn_tj)
            output['deepcsv_CvB_dn_2'].fill(dataset=dataset,flav=genflavor[:,2],  deepcsv_CvB_dn=np.where(sjets[:,2].btagDeepCvB<0,-0.2,sjets[:,2].btagDeepCvB),weight=weights.weight()[event_level]*csvsfs_c_dn_tj)
            output['deepcsv_CvL_dn_2'].fill(dataset=dataset,flav=genflavor[:,2],   deepcsv_CvL_dn=np.where(sjets[:,2].btagDeepCvL<0,-0.2,sjets[:,2].btagDeepCvL),weight=weights.weight()[event_level]*csvsfs_c_dn_tj)
            output['deepflav_CvB_dn_2'].fill(dataset=dataset,flav=genflavor[:,2],  deepflav_CvB_dn=np.where(sjets[:,2].btagDeepFlavCvB<0,-0.2,sjets[:,2].btagDeepFlavCvB),weight=weights.weight()[event_level]*jetsfs_c_dn_tj)
            output['deepflav_CvL_dn_2'].fill(dataset=dataset,flav=genflavor[:,2],  deepflav_CvL_dn=np.where(sjets[:,2].btagDeepFlavCvL<0,-0.2,sjets[:,2].btagDeepFlavCvL),weight=weights.weight()[event_level]*jetsfs_c_dn_tj)
            output['btagDeepFlavB_dn_3'].fill(dataset=dataset,flav=genflavor[:,3],  btagDeepFlavB_dn=sjets[:,3].btagDeepFlavB,weight=weights.weight()[event_level]*jetsfs_b_dn_fj)
            output['btagDeepFlavC_dn_3'].fill(dataset=dataset,flav=genflavor[:,3],  btagDeepFlavC_dn=sjets[:,3].btagDeepFlavC,weight=weights.weight()[event_level]*jetsfs_c_dn_fj)
            output['btagDeepB_dn_3'].fill(dataset=dataset,flav=genflavor[:,3],  btagDeepB_dn=sjets[:,3].btagDeepB,weight=weights.weight()[event_level]*csvsfs_b_dn_fj)
            output['btagDeepC_dn_3'].fill(dataset=dataset,flav=genflavor[:,3],  btagDeepC_dn=sjets[:,3].btagDeepC,weight=weights.weight()[event_level]*csvsfs_c_dn_fj)
            output['deepcsv_CvB_dn_3'].fill(dataset=dataset,flav=genflavor[:,3],  deepcsv_CvB_dn=np.where(sjets[:,3].btagDeepCvB<0,-0.2,sjets[:,3].btagDeepCvB),weight=weights.weight()[event_level]*csvsfs_c_dn_fj)
            output['deepcsv_CvL_dn_3'].fill(dataset=dataset,flav=genflavor[:,3],   deepcsv_CvL_dn=np.where(sjets[:,3].btagDeepCvL<0,-0.2,sjets[:,3].btagDeepCvL),weight=weights.weight()[event_level]*csvsfs_c_dn_fj)
            output['deepflav_CvB_dn_3'].fill(dataset=dataset,flav=genflavor[:,3],  deepflav_CvB_dn=np.where(sjets[:,3].btagDeepFlavCvB<0,-0.2,sjets[:,3].btagDeepFlavCvB),weight=weights.weight()[event_level]*jetsfs_c_dn_fj)
            output['deepflav_CvL_dn_3'].fill(dataset=dataset,flav=genflavor[:,3],  deepflav_CvL_dn=np.where(sjets[:,3].btagDeepFlavCvL<0,-0.2,sjets[:,3].btagDeepFlavCvL),weight=weights.weight()[event_level]*jetsfs_c_dn_fj)
        else:
            output['btagDeepFlavB_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepFlavB=sjets[:,0].btagDeepFlavB)
            output['btagDeepFlavC_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepFlavC=sjets[:,0].btagDeepFlavC)
            output['btagDeepB_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepB=sjets[:,0].btagDeepB)
            output['btagDeepC_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepC=sjets[:,0].btagDeepC)
            output['deepcsv_CvB_0'].fill(dataset=dataset,flav=genflavor[:,0], deepcsv_CvB=np.where(sjets[:,0].btagDeepCvB<0,-0.2,sjets[:,0].btagDeepCvB))
            output['deepcsv_CvL_0'].fill(dataset=dataset,flav=genflavor[:,0], deepcsv_CvL=np.where(sjets[:,0].btagDeepCvL<0,-0.2,sjets[:,0].btagDeepCvL))
            output['deepflav_CvB_0'].fill(dataset=dataset,flav=genflavor[:,0], deepflav_CvB=np.where(sjets[:,0].btagDeepFlavCvB<0,-0.2,sjets[:,0].btagDeepFlavCvB))
            output['deepflav_CvL_0'].fill(dataset=dataset,flav=genflavor[:,0], deepflav_CvL=np.where(sjets[:,0].btagDeepFlavCvL<0,-0.2,sjets[:,0].btagDeepFlavCvL))
            output['btagDeepFlavB_1'].fill(dataset=dataset,flav=genflavor[:,1],  btagDeepFlavB=sjets[:,1].btagDeepFlavB)
            output['btagDeepFlavC_1'].fill(dataset=dataset,flav=genflavor[:,1],  btagDeepFlavC=sjets[:,1].btagDeepFlavC)
            output['btagDeepB_1'].fill(dataset=dataset,flav=genflavor[:,1],  btagDeepB=sjets[:,1].btagDeepB)
            output['btagDeepC_1'].fill(dataset=dataset,flav=genflavor[:,1],  btagDeepC=sjets[:,1].btagDeepC)
            output['deepcsv_CvB_1'].fill(dataset=dataset,flav=genflavor[:,1],  deepcsv_CvB=np.where(sjets[:,1].btagDeepCvB<0,-0.2,sjets[:,1].btagDeepCvB))
            output['deepcsv_CvL_1'].fill(dataset=dataset,flav=genflavor[:,1],  deepcsv_CvL=np.where(sjets[:,1].btagDeepCvL<0,-0.2,sjets[:,1].btagDeepCvL))
            output['deepflav_CvB_1'].fill(dataset=dataset,flav=genflavor[:,1],  deepflav_CvB=np.where(sjets[:,1].btagDeepFlavCvB<0,-0.2,sjets[:,1].btagDeepFlavCvB))
            output['deepflav_CvL_1'].fill(dataset=dataset,flav=genflavor[:,1],  deepflav_CvL=np.where(sjets[:,1].btagDeepFlavCvL<0,-0.2,sjets[:,1].btagDeepFlavCvL))
            output['btagDeepFlavB_2'].fill(dataset=dataset,flav=genflavor[:,2],  btagDeepFlavB=sjets[:,2].btagDeepFlavB)
            output['btagDeepFlavC_2'].fill(dataset=dataset,flav=genflavor[:,2],  btagDeepFlavC=sjets[:,2].btagDeepFlavC)
            output['btagDeepB_2'].fill(dataset=dataset,flav=genflavor[:,2],  btagDeepB=sjets[:,2].btagDeepB)
            output['btagDeepC_2'].fill(dataset=dataset,flav=genflavor[:,2],  btagDeepC=sjets[:,2].btagDeepC)
            output['deepcsv_CvB_2'].fill(dataset=dataset,flav=genflavor[:,2],  deepcsv_CvB=np.where(sjets[:,2].btagDeepCvB<0,-0.2,sjets[:,2].btagDeepCvB))
            output['deepcsv_CvL_2'].fill(dataset=dataset,flav=genflavor[:,2],  deepcsv_CvL=np.where(sjets[:,2].btagDeepCvL<0,-0.2,sjets[:,2].btagDeepCvL))
            output['deepflav_CvB_2'].fill(dataset=dataset,flav=genflavor[:,2],  deepflav_CvB=np.where(sjets[:,2].btagDeepFlavCvB<0,-0.2,sjets[:,2].btagDeepFlavCvB))
            output['deepflav_CvL_2'].fill(dataset=dataset,flav=genflavor[:,2],  deepflav_CvL=np.where(sjets[:,2].btagDeepFlavCvL<0,-0.2,sjets[:,2].btagDeepFlavCvL))
            output['btagDeepFlavB_3'].fill(dataset=dataset,flav=genflavor[:,3],  btagDeepFlavB=sjets[:,3].btagDeepFlavB)
            output['btagDeepFlavC_3'].fill(dataset=dataset,flav=genflavor[:,3],  btagDeepFlavC=sjets[:,3].btagDeepFlavC)
            output['btagDeepB_3'].fill(dataset=dataset,flav=genflavor[:,3],  btagDeepB=sjets[:,3].btagDeepB)
            output['btagDeepC_3'].fill(dataset=dataset,flav=genflavor[:,3],  btagDeepC=sjets[:,3].btagDeepC)
            output['deepcsv_CvB_3'].fill(dataset=dataset,flav=genflavor[:,3],  deepcsv_CvB=np.where(sjets[:,3].btagDeepCvB<0,-0.2,sjets[:,3].btagDeepCvB))
            output['deepcsv_CvL_3'].fill(dataset=dataset,flav=genflavor[:,3],  deepcsv_CvL=np.where(sjets[:,3].btagDeepCvL<0,-0.2,sjets[:,3].btagDeepCvL))
            output['deepflav_CvB_3'].fill(dataset=dataset,flav=genflavor[:,3],  deepflav_CvB=np.where(sjets[:,3].btagDeepFlavCvB<0,-0.2,sjets[:,3].btagDeepFlavCvB))
            output['deepflav_CvL_3'].fill(dataset=dataset,flav=genflavor[:,3],  deepflav_CvL=np.where(sjets[:,3].btagDeepFlavCvL<0,-0.2,sjets[:,3].btagDeepFlavCvL))
            output['btagDeepFlavBSF_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepFlavBSF=sjets[:,0].btagDeepFlavB)      
            output['btagDeepFlavCSF_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepFlavCSF=sjets[:,0].btagDeepFlavC)
            output['btagDeepBSF_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepBSF=sjets[:,0].btagDeepB)
            output['btagDeepCSF_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepCSF=sjets[:,0].btagDeepC)
            output['deepcsv_CvBSF_0'].fill(dataset=dataset,flav=genflavor[:,0], deepcsv_CvBSF=np.where(sjets[:,0].btagDeepCvB<0,-0.2,sjets[:,0].btagDeepCvB))
            output['deepcsv_CvLSF_0'].fill(dataset=dataset,flav=genflavor[:,0], deepcsv_CvLSF=np.where(sjets[:,0].btagDeepCvL<0,-0.2,sjets[:,0].btagDeepCvL))
            output['deepflav_CvBSF_0'].fill(dataset=dataset,flav=genflavor[:,0], deepflav_CvBSF=np.where(sjets[:,0].btagDeepFlavCvB<0,-0.2,sjets[:,0].btagDeepFlavCvB))
            output['deepflav_CvLSF_0'].fill(dataset=dataset,flav=genflavor[:,0], deepflav_CvLSF=np.where(sjets[:,0].btagDeepFlavCvL<0,-0.2,sjets[:,0].btagDeepFlavCvL))
            output['btagDeepFlavBSF_1'].fill(dataset=dataset,flav=genflavor[:,1],  btagDeepFlavBSF=sjets[:,1].btagDeepFlavB)
            output['btagDeepFlavCSF_1'].fill(dataset=dataset,flav=genflavor[:,1],  btagDeepFlavCSF=sjets[:,1].btagDeepFlavC)
            output['btagDeepBSF_1'].fill(dataset=dataset,flav=genflavor[:,1],  btagDeepBSF=sjets[:,1].btagDeepB)
            output['btagDeepCSF_1'].fill(dataset=dataset,flav=genflavor[:,1],  btagDeepCSF=sjets[:,1].btagDeepC)
            output['deepcsv_CvBSF_1'].fill(dataset=dataset,flav=genflavor[:,1],  deepcsv_CvBSF=np.where(sjets[:,1].btagDeepCvB<0,-0.2,sjets[:,1].btagDeepCvB))
            output['deepcsv_CvLSF_1'].fill(dataset=dataset,flav=genflavor[:,1],  deepcsv_CvLSF=np.where(sjets[:,1].btagDeepCvL<0,-0.2,sjets[:,1].btagDeepCvL))
            output['deepflav_CvBSF_1'].fill(dataset=dataset,flav=genflavor[:,1],  deepflav_CvBSF=np.where(sjets[:,1].btagDeepFlavCvB<0,-0.2,sjets[:,1].btagDeepFlavCvB))
            output['deepflav_CvLSF_1'].fill(dataset=dataset,flav=genflavor[:,1],  deepflav_CvLSF=np.where(sjets[:,1].btagDeepFlavCvL<0,-0.2,sjets[:,1].btagDeepFlavCvL))
            output['btagDeepFlavBSF_2'].fill(dataset=dataset,flav=genflavor[:,2],  btagDeepFlavBSF=sjets[:,2].btagDeepFlavB)
            output['btagDeepFlavCSF_2'].fill(dataset=dataset,flav=genflavor[:,2],  btagDeepFlavCSF=sjets[:,2].btagDeepFlavC)
            output['btagDeepBSF_2'].fill(dataset=dataset,flav=genflavor[:,2],  btagDeepBSF=sjets[:,2].btagDeepB)
            output['btagDeepCSF_2'].fill(dataset=dataset,flav=genflavor[:,2],  btagDeepCSF=sjets[:,2].btagDeepC)
            output['deepcsv_CvBSF_2'].fill(dataset=dataset,flav=genflavor[:,2],  deepcsv_CvBSF=np.where(sjets[:,2].btagDeepCvB<0,-0.2,sjets[:,2].btagDeepCvB))
            output['deepcsv_CvLSF_2'].fill(dataset=dataset,flav=genflavor[:,2],  deepcsv_CvLSF=np.where(sjets[:,2].btagDeepCvL<0,-0.2,sjets[:,2].btagDeepCvL))
            output['deepflav_CvBSF_2'].fill(dataset=dataset,flav=genflavor[:,2],  deepflav_CvBSF=np.where(sjets[:,2].btagDeepFlavCvB<0,-0.2,sjets[:,2].btagDeepFlavCvB))
            output['deepflav_CvLSF_2'].fill(dataset=dataset,flav=genflavor[:,2],  deepflav_CvLSF=np.where(sjets[:,2].btagDeepFlavCvL<0,-0.2,sjets[:,2].btagDeepFlavCvL))
            output['btagDeepFlavBSF_3'].fill(dataset=dataset,flav=genflavor[:,3],  btagDeepFlavBSF=sjets[:,3].btagDeepFlavB)
            output['btagDeepFlavCSF_3'].fill(dataset=dataset,flav=genflavor[:,3],  btagDeepFlavCSF=sjets[:,3].btagDeepFlavC)
            output['btagDeepBSF_3'].fill(dataset=dataset,flav=genflavor[:,3],  btagDeepBSF=sjets[:,3].btagDeepB)
            output['btagDeepCSF_3'].fill(dataset=dataset,flav=genflavor[:,3],  btagDeepCSF=sjets[:,3].btagDeepC)
            output['deepcsv_CvBSF_3'].fill(dataset=dataset,flav=genflavor[:,3],  deepcsv_CvBSF=np.where(sjets[:,3].btagDeepCvB<0,-0.2,sjets[:,3].btagDeepCvB))
            output['deepcsv_CvLSF_3'].fill(dataset=dataset,flav=genflavor[:,3],  deepcsv_CvLSF=np.where(sjets[:,3].btagDeepCvL<0,-0.2,sjets[:,3].btagDeepCvL))
            output['deepflav_CvBSF_3'].fill(dataset=dataset,flav=genflavor[:,3],  deepflav_CvBSF=np.where(sjets[:,3].btagDeepFlavCvB<0,-0.2,sjets[:,3].btagDeepFlavCvB))
            output['deepflav_CvLSF_3'].fill(dataset=dataset,flav=genflavor[:,3],  deepflav_CvLSF=np.where(sjets[:,3].btagDeepFlavCvL<0,-0.2,sjets[:,3].btagDeepFlavCvL))
        def flatten(ar): # flatten awkward into a 1d array to hist
            return ak.flatten(ar, axis=None)

        def num(ar):
            return ak.num(ak.fill_none(ar[~ak.is_none(ar)], 0), axis=0)
        if(isRealData):
            output['njet'].fill(dataset=dataset,  njet=flatten(ak.num(sel_jets)))
            output['nbjet'].fill(dataset=dataset,  nbjet=flatten(ak.num(sbjets)))
            output['ntbjet'].fill(dataset=dataset,  ntbjet=flatten(ak.num(sbjets)))
            output['nmu'].fill(dataset=dataset,   nmu=flatten(ak.num(smu)))
            output['lmupt'].fill(dataset=dataset, lmupt=flatten((smu.pt)))
            output['ljpt'].fill(dataset=dataset, flav=0, ljpt=flatten(sjets[:,0].pt))
            output['sljpt'].fill(dataset=dataset, flav=0, sljpt=flatten(sjets[:,1].pt))
            output['ssljpt'].fill(dataset=dataset, flav=0, ssljpt=flatten(sjets[:,2].pt))
            output['sssljpt'].fill(dataset=dataset, flav=0, sssljpt=flatten(sjets[:,3].pt))
            output['ljeta'].fill(dataset=dataset, flav=0, ljeta=flatten(sjets[:,0].eta))
            output['sljeta'].fill(dataset=dataset, flav=0, sljeta=flatten(sjets[:,1].eta))
            output['ssljeta'].fill(dataset=dataset, flav=0, ssljeta=flatten(sjets[:,2].eta))
            output['sssljeta'].fill(dataset=dataset, flav=0, sssljeta=flatten(sjets[:,3].eta))
            output['ljdr'].fill(dataset=dataset, flav=0, ljdr=flatten(sjets[:,0].delta_r(smu)))
            output['sljdr'].fill(dataset=dataset, flav=0, sljdr=flatten(sjets[:,1].delta_r(smu)))
            output['ssljdr'].fill(dataset=dataset, flav=0, ssljdr=flatten(sjets[:,2].delta_r(smu)))
            output['sssljdr'].fill(dataset=dataset, flav=0, sssljdr=flatten(sjets[:,3].delta_r(smu)))
            output['met'].fill(dataset=dataset, met=flatten((MET)))
        else:
            output['njet'].fill(dataset=dataset,  njet=flatten(ak.num(sel_jets)),weight=weights.weight()[event_level])
            output['nbjet'].fill(dataset=dataset,  nbjet=flatten(ak.num(sbjets)),weight=weights.weight()[event_level]*csvsfs_b_lj*csvsfs_b_sj*csvsfs_b_tj*csvsfs_b_fj)
            output['nbjet_up'].fill(dataset=dataset,  nbjet=flatten(ak.num(sbjets)),weight=weights.weight()[event_level]*csvsfs_b_up_lj*csvsfs_b_up_sj*csvsfs_b_up_tj*csvsfs_b_up_fj)
            output['nbjet_dn'].fill(dataset=dataset,  nbjet=flatten(ak.num(sbjets)),weight=weights.weight()[event_level]*csvsfs_b_dn_lj*csvsfs_b_dn_sj*csvsfs_b_dn_tj*csvsfs_b_dn_fj)
            output['ntbjet'].fill(dataset=dataset,  ntbjet=flatten(ak.num(stbjets)),weight=weights.weight()[event_level])
            output['nmu'].fill(dataset=dataset,   nmu=flatten(ak.num(smu)),weight=weights.weight()[event_level])
            output['lmupt'].fill(dataset=dataset, lmupt=flatten((smu.pt)),weight=weights.weight()[event_level])
            output['ljpt'].fill(dataset=dataset, flav=flatten(sjets[:,0].hadronFlavour), ljpt=flatten(sjets[:,0].pt),weight=weights.weight()[event_level])
            output['sljpt'].fill(dataset=dataset, flav=flatten(sjets[:,1].hadronFlavour), sljpt=flatten(sjets[:,1].pt),weight=weights.weight()[event_level])
            output['ssljpt'].fill(dataset=dataset, flav=flatten(sjets[:,2].hadronFlavour), ssljpt=flatten(sjets[:,2].pt),weight=weights.weight()[event_level])
            output['sssljpt'].fill(dataset=dataset, flav=flatten(sjets[:,3].hadronFlavour), sssljpt=flatten(sjets[:,3].pt),weight=weights.weight()[event_level])
            output['ljeta'].fill(dataset=dataset, flav=flatten(sjets[:,0].hadronFlavour), ljeta=flatten(sjets[:,0].eta),weight=weights.weight()[event_level])
            output['sljeta'].fill(dataset=dataset, flav=flatten(sjets[:,1].hadronFlavour), sljeta=flatten(sjets[:,1].eta),weight=weights.weight()[event_level])
            output['ssljeta'].fill(dataset=dataset, flav=flatten(sjets[:,2].hadronFlavour), ssljeta=flatten(sjets[:,2].eta),weight=weights.weight()[event_level])
            output['sssljeta'].fill(dataset=dataset, flav=flatten(sjets[:,3].hadronFlavour), sssljeta=flatten(sjets[:,3].eta),weight=weights.weight()[event_level])
            
            output['ljdr'].fill(dataset=dataset, flav=flatten(sjets[:,0].hadronFlavour), ljdr=flatten(sjets[:,0].delta_r(smu)),weight=weights.weight()[event_level])
            output['sljdr'].fill(dataset=dataset, flav=flatten(sjets[:,1].hadronFlavour), sljdr=flatten(sjets[:,1].delta_r(smu)),weight=weights.weight()[event_level])
            output['ssljdr'].fill(dataset=dataset, flav=flatten(sjets[:,2].hadronFlavour), ssljdr=flatten(sjets[:,2].delta_r(smu)),weight=weights.weight()[event_level])
            output['sssljdr'].fill(dataset=dataset, flav=flatten(sjets[:,3].hadronFlavour), sssljdr=flatten(sjets[:,3].delta_r(smu)),weight=weights.weight()[event_level])
            output['met'].fill(dataset=dataset, met=flatten((MET)),weight=weights.weight()[event_level])
        gc.collect()
        return output

    def postprocess(self, accumulator):
        return accumulator
