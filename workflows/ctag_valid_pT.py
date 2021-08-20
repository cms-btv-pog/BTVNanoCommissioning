import gzip
import pickle
import coffea
from coffea import hist, processor
import numpy as np
import awkward as ak
from coffea.analysis_tools import Weights
from coffea.lumi_tools import LumiMask
from coffea.btag_tools import BTagScaleFactor
from cTagSFReader import *
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


        # Jet
        jet_pt_axis   = hist.Bin("pt",   r"Jet $p_{T}$ [GeV]", 100, 0, 500)
        jet_rawpt_axis   = hist.Bin("rawpt",   r"Jet $p_{T}$ [GeV]", 100, 0, 500)
        jet_eta_axis  = hist.Bin("eta",  r"Jet $\eta$", 60, -3, 3)
        jet_phi_axis  = hist.Bin("phi",  r"Jet $\phi$", 60, -3, 3)
        jet_mass_axis = hist.Bin("mass", r"Jet $m$ [GeV]", 100, 0, 50)
        jet_dr_axis = hist.Bin("dr", r"Jet $\Delta$R(l,j)", 50, 0, 5)
        ljpt_axis     = hist.Bin("ljpt", r"Leading jet $p_{T}$ [GeV]", 100, 20, 400)
        sljpt_axis     = hist.Bin("sljpt", r"Subleading jet $p_{T}$ [GeV]", 100, 20, 400)        
        ssljpt_axis     = hist.Bin("ssljpt", r"3rd jet $p_{T}$ [GeV]", 100, 20, 400)
        sssljpt_axis     = hist.Bin("sssljpt", r"4th jet $p_{T}$ [GeV]", 100, 20, 400) 
        ljrawpt_axis     = hist.Bin("ljrawpt", r"Leading jet $p_{T}$ [GeV]", 100, 20, 400)
        sljrawpt_axis     = hist.Bin("sljrawpt", r"Subleading jet $p_{T}$ [GeV]", 100, 20, 400)        
        ssljrawpt_axis     = hist.Bin("ssljrawpt", r"3rd jet $p_{T}$ [GeV]", 100, 20, 400)
        sssljrawpt_axis     = hist.Bin("sssljrawpt", r"4th jet $p_{T}$ [GeV]", 100, 20, 400) 
        ljdr_axis = hist.Bin("ljdr", "Leading jet $\Delta$R(l,j)", 50,0,5)
        sljdr_axis = hist.Bin("sljdr", "Subleading jet $\Delta$R(l,j)", 50,0,5)
        ssljdr_axis = hist.Bin("ssljdr", "3rd jet $\Delta$R(l,j)", 50,0,5)
        sssljdr_axis = hist.Bin("sssljdr", "4th jet $\Delta$R(l,j)", 50,0,5)

        # Define similar axes dynamically
        disc_list = ["btagCMVA", "btagCSVV2", 'btagDeepB', 'btagDeepC', 'btagDeepFlavB', 'btagDeepFlavC']
        discSF_list = ['btagDeepBSF', 'btagDeepCSF', 'btagDeepFlavBSF', 'btagDeepFlavCSF','btagDeepB', 'btagDeepC', 'btagDeepFlavB', 'btagDeepFlavC']
        ddx_list = ["btagDDBvLV2","btagDDCvBV2","btagDDCvLV2"]
        btag_axes = []
        for d in disc_list:
            btag_axes.append(hist.Bin(d, d, 50, 0, 1))     
        deepb_axis = hist.Bin("btagDeepBSF", "btagDeepBSF", 50, 0, 1)     
        deepc_axis = hist.Bin("btagDeepCSF", "btagDeepCSF", 50, 0, 1) 
        deepcsv_CvL_axis = hist.Bin("deepcsv_CvL", "deepcsv_CvL", 50, 0, 1)  
        deepcsv_BvC_axis = hist.Bin("deepcsv_BvC", "deepcsv_BvC", 50, 0, 1)  
        deepflav_CvL_axis = hist.Bin("deepflav_CvL", "deepflav_CvL", 50, 0, 1)  
        deepflav_BvC_axis = hist.Bin("deepflav_BvC", "deepflav_BvC", 50, 0, 1)
        deepcsv_CvL_up_axis = hist.Bin("deepcsv_CvL_up", "deepcsv_CvL_up", 50, 0, 1)  
        deepcsv_BvC_up_axis = hist.Bin("deepcsv_BvC_up", "deepcsv_BvC_up", 50, 0, 1)  
        deepflav_CvL_up_axis = hist.Bin("deepflav_CvL_up", "deepflav_CvL_up", 50, 0, 1)  
        deepflav_BvC_up_axis = hist.Bin("deepflav_BvC_up", "deepflav_BvC_up", 50, 0, 1)
        deepcsv_CvL_dn_axis = hist.Bin("deepcsv_CvL_dn", "deepcsv_CvL_dn", 50, 0, 1)  
        deepcsv_BvC_dn_axis = hist.Bin("deepcsv_BvC_dn", "deepcsv_BvC_dn", 50, 0, 1)  
        deepflav_CvL_dn_axis = hist.Bin("deepflav_CvL_dn", "deepflav_CvL_dn", 50, 0, 1)  
        deepflav_BvC_dn_axis = hist.Bin("deepflav_BvC_dn", "deepflav_BvC_dn", 50, 0, 1)  
        deepjb_axis = hist.Bin("btagDeepFlavBSF", "btagDeepFlavBSF", 50, 0, 1)     
        deepjc_axis = hist.Bin("btagDeepFlavCSF", "btagDeepFlavCSF", 50, 0, 1)  
        deepb_axis_up = hist.Bin("btagDeepBSF_up", "btagDeepBSF_up", 50, 0, 1)     
        deepc_axis_up = hist.Bin("btagDeepCSF_up", "btagDeepCSF_up", 50, 0, 1)  
        deepjb_axis_up = hist.Bin("btagDeepFlavBSF_up", "btagDeepFlavBSF_up", 50, 0, 1)     
        deepjc_axis_up = hist.Bin("btagDeepFlavCSF_up", "btagDeepFlavCSF_up", 50, 0, 1)
        deepb_axis_down = hist.Bin("btagDeepBSF_down", "btagDeepBSF_down", 50, 0, 1)     
        deepc_axis_down = hist.Bin("btagDeepCSF_down", "btagDeepCSF_down", 50, 0, 1)  
        deepjb_axis_down = hist.Bin("btagDeepFlavBSF_down", "btagDeepFlavBSF_down", 50, 0, 1)     
        deepjc_axis_down = hist.Bin("btagDeepFlavCSF_down", "btagDeepFlavCSF_down", 50, 0, 1)     
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
        
        
        # Define histograms from axes
        _hist_jet_dict = {
                'pt'  : hist.Hist("Counts", dataset_axis, flav_axis,jet_pt_axis),
                'rawpt'  : hist.Hist("Counts", dataset_axis, flav_axis,jet_rawpt_axis),
                'eta' : hist.Hist("Counts", dataset_axis, flav_axis,jet_eta_axis),
                'phi' : hist.Hist("Counts", dataset_axis, flav_axis,jet_phi_axis),
                'mass': hist.Hist("Counts", dataset_axis, flav_axis,jet_mass_axis),
                # 'dr': hist.Hist("Counts", dataset_axis, flav_axis,jet_dr_axis)
            }
        _hist_deepcsv_dict = {
                'pt'  : hist.Hist("Counts", dataset_axis, flav_axis,jet_pt_axis),
                'rawpt'  : hist.Hist("Counts", dataset_axis, flav_axis,jet_rawpt_axis),
                'eta' : hist.Hist("Counts", dataset_axis, flav_axis,jet_eta_axis),
                'phi' : hist.Hist("Counts", dataset_axis, flav_axis,jet_phi_axis),
                'mass': hist.Hist("Counts", dataset_axis, flav_axis,jet_mass_axis),

            }
        
        _hist_deepcsvSF_dict={
            'btagDeepBSF' : hist.Hist("Counts", dataset_axis, flav_axis,jet_eta_axis, jet_pt_axis, deepb_axis),
            # 'btagDeepCSF' : hist.Hist("Counts", dataset_axis, flav_axis,jet_eta_axis, jet_pt_axis, deepc_axis),
            # 'btagDeepBvC' : hist.Hist("Counts", dataset_axis, flav_axis,jet_eta_axis, jet_pt_axis, deepcsv_BvC_axis),
            # 'btagDeepCvL' : hist.Hist("Counts", dataset_axis, flav_axis,jet_eta_axis, jet_pt_axis, deepcsv_CvL_axis),
            'btagDeepFlavBSF' : hist.Hist("Counts", dataset_axis, flav_axis, jet_eta_axis, jet_pt_axis,deepjb_axis),
            # 'btagDeepFlavCSF' : hist.Hist("Counts", dataset_axis, flav_axis,jet_eta_axis, jet_pt_axis, deepjc_axis),
            # 'btagDeepFlavBvC' : hist.Hist("Counts", dataset_axis, flav_axis,jet_eta_axis, jet_pt_axis, deepflav_BvC_axis),
            # 'btagDeepFlavCvL' : hist.Hist("Counts", dataset_axis, flav_axis,jet_eta_axis, jet_pt_axis, deepflav_CvL_axis),
            'btagDeepBSF_up' : hist.Hist("Counts", dataset_axis, flav_axis,jet_eta_axis, jet_pt_axis, deepb_axis_up),
            # 'btagDeepCSF_up' : hist.Hist("Counts", dataset_axis, flav_axis,jet_eta_axis, jet_pt_axis, deepc_axis_up),
            # 'btagDeepBvC_up' : hist.Hist("Counts", dataset_axis, flav_axis,jet_eta_axis, jet_pt_axis, deepcsv_BvC_up_axis),
            # 'btagDeepCvL_up' : hist.Hist("Counts", dataset_axis, flav_axis,jet_eta_axis, jet_pt_axis, deepcsv_CvL_up_axis),
            'btagDeepFlavBSF_up' : hist.Hist("Counts", dataset_axis, flav_axis, jet_eta_axis, jet_pt_axis,deepjb_axis_up),
            # 'btagDeepFlavCSF_up' : hist.Hist("Counts", dataset_axis, flav_axis,jet_eta_axis, jet_pt_axis, deepjc_axis_up),
            # 'btagDeepFlavBvC_up' : hist.Hist("Counts", dataset_axis, flav_axis,jet_eta_axis, jet_pt_axis, deepflav_BvC_up_axis),
            # 'btagDeepFlavCvL_up' : hist.Hist("Counts", dataset_axis, flav_axis,jet_eta_axis, jet_pt_axis, deepflav_CvL_up_axis),
            'btagDeepBSF_down' : hist.Hist("Counts", dataset_axis, flav_axis,jet_eta_axis, jet_pt_axis, deepb_axis_down),
            # 'btagDeepCSF_down' : hist.Hist("Counts", dataset_axis, flav_axis,jet_eta_axis, jet_pt_axis, deepc_axis_down),
            # 'btagDeepBvC_down' : hist.Hist("Counts", dataset_axis, flav_axis,jet_eta_axis, jet_pt_axis, deepcsv_BvC_dn_axis),
            # 'btagDeepCvL_down' : hist.Hist("Counts", dataset_axis, flav_axis,jet_eta_axis, jet_pt_axis, deepcsv_CvL_dn_axis),
            'btagDeepFlavBSF_down' : hist.Hist("Counts", dataset_axis, flav_axis, jet_eta_axis, jet_pt_axis,deepjb_axis_down),
            # 'btagDeepFlavCSF_down' : hist.Hist("Counts", dataset_axis, flav_axis,jet_eta_axis, jet_pt_axis, deepjc_axis_down),
            # 'btagDeepFlavBvC_down' : hist.Hist("Counts", dataset_axis, flav_axis,jet_eta_axis, jet_pt_axis, deepflav_BvC_dn_axis),
            # 'btagDeepFlavCvL_down' : hist.Hist("Counts", dataset_axis, flav_axis,jet_eta_axis, jet_pt_axis, deepflav_CvL_dn_axis)       
        }
        # Generate some histograms dynamically
        for disc, axis in zip(disc_list, btag_axes):
            _hist_deepcsvSF_dict[disc] = hist.Hist("Counts", dataset_axis, flav_axis,jet_pt_axis,jet_eta_axis, axis)
        for disc, axis in zip(disc_list, btag_axes):
            _hist_deepcsv_dict[disc] = hist.Hist("Counts", dataset_axis, flav_axis,jet_pt_axis,jet_eta_axis, axis)
        for deepcsv, axises in zip(deepcsv_list, deepcsv_axes):
             _hist_deepcsv_dict[deepcsv] = hist.Hist("Counts", dataset_axis,flav_axis,  axises)
        _hist_event_dict = {
                'njet'  : hist.Hist("Counts", dataset_axis,  njet_axis),
                'nbjet' : hist.Hist("Counts", dataset_axis,  nbjet_axis),
                'ntbjet' : hist.Hist("Counts", dataset_axis,  ntbjet_axis),
                'nmu'   : hist.Hist("Counts", dataset_axis,  nmu_axis),
                'lmupt' : hist.Hist("Counts", dataset_axis,  lmupt_axis),
                'ljpt'  : hist.Hist("Counts", dataset_axis, flav_axis, ljpt_axis),
                'sljpt'  : hist.Hist("Counts", dataset_axis, flav_axis, sljpt_axis),
                'ssljpt'  : hist.Hist("Counts", dataset_axis, flav_axis, ssljpt_axis),
                'sssljpt'  : hist.Hist("Counts", dataset_axis, flav_axis, sssljpt_axis),
                'ljrawpt'  : hist.Hist("Counts", dataset_axis, flav_axis, ljrawpt_axis),
                'sljrawpt'  : hist.Hist("Counts", dataset_axis, flav_axis, sljrawpt_axis),
                'ssljrawpt'  : hist.Hist("Counts", dataset_axis, flav_axis, ssljrawpt_axis),
                'sssljrawpt'  : hist.Hist("Counts", dataset_axis, flav_axis, sssljrawpt_axis),
                'ljdr'  : hist.Hist("Counts", dataset_axis, flav_axis, ljdr_axis),
                'sljdr'  : hist.Hist("Counts", dataset_axis, flav_axis, sljdr_axis),
                'ssljdr'  : hist.Hist("Counts", dataset_axis, flav_axis, ssljdr_axis),
                'sssljdr'  : hist.Hist("Counts", dataset_axis, flav_axis, sssljdr_axis),
                'met' : hist.Hist("Counts", dataset_axis, met_axis)
            }
        self.jet_hists = list(_hist_jet_dict.keys())
        # self.deepcsv_hists = list(_hist_deepcsv_dict.keys())
        self.deepcsvSF_hists = list(_hist_deepcsvSF_dict.keys())
        # self.deepddx_hists = list(_hist_deepddx_dict.keys())
        self.event_hists = list(_hist_event_dict.keys())
        # _hist_dict = {**_hist_deepcsv_dict,**_hist_event_dict,**_hist_deepcsvSF_dict}
        _hist_dict = {**_hist_event_dict,**_hist_deepcsvSF_dict}
        self._accumulator = processor.dict_accumulator(_hist_dict)
        self._accumulator['sumw'] = processor.defaultdict_accumulator(float)


    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.accumulator.identity()
        # print(output.items())
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
        
        # print(ak.all(events.Jet.metric_table(events.Muon) > 0.4, axis=2))
        event_muon = ak.pad_none(events.Muon, 1, axis=1) 
        
        req_muon =(ak.count(event_muon.pt, axis=1) == 1)
        ## Jet cuts 
        
        event_jet = events.Jet[(events.Jet.pt> 25) & (abs(events.Jet.eta) <= 2.4)&(events.Jet.puId > 0)&(events.Jet.jetId > 0)&(ak.all(events.Jet.metric_table(events.Muon) > 0.4, axis=2))]
        # print(events.Jet.muonIdx1)                                                                               
        
        req_jets = (ak.num(event_jet)>=4)
        req_MET = events.METFixEE2017.pt>50
        
        event_level = req_trig  & req_jets & req_muon&req_MET &req_lumi
        #event_level = req_trig  & req_jets & req_muon&req_MET
        # Selected
        selev = events[event_level]    
        
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
        sbjets = selev.Jet[bjet_level]
        
        # sjets.btag_CvL_flav = sjets.btagDeepFlavC/(1.-sjets.btagDeepFlavB)
        # sjets.btag_CvB_flav = sjets.btagDeepFlavC/(sjets.btagDeepFlavC+sjets.btagDeepFlavB)
        # sjets.btag_CvL_csv = sjets.btagDeepC/(1.-sjets.btagDeepB)
        # sjets.btag_CvB_csv = sjets.btagDeepC/(sjets.btagDeepC+sjets.btagDeepB)
        if not isRealData: stbjets = sbjets[sbjets.hadronFlavour==5]
        if isRealData :
            genflavor = ak.zeros_like(sjets.pt)
        else:
            par_flav = (sjets.partonFlavour == 0 ) & (sjets.hadronFlavour==0)
            genflavor = sjets.hadronFlavour + 1*par_flav 
            genweiev=ak.flatten(ak.broadcast_arrays(weights.weight()[event_level],sjets['pt'])[0])
            # jetsfs_c = getSF(ak.to_numpy(ak.flatten(sjets.hadronFlavour)),ak.to_numpy(ak.flatten(sjets.btagDeepFlavC/(1.-sjets.btagDeepFlavB))),ak.to_numpy(ak.flatten(sjets.btagDeepFlavC/(sjets.btagDeepFlavC+sjets.btagDeepFlavB))),deepjetc_sf)*genweiev
            # csvsfs_c = getSF(ak.to_numpy(ak.flatten(sjets.hadronFlavour)),ak.to_numpy(ak.flatten(sjets.btagDeepC/(1.-sjets.btagDeepB))),ak.to_numpy(ak.flatten(sjets.btagDeepC/(sjets.btagDeepC+sjets.btagDeepB))),deepcsvc_sf)*genweiev
            # jetsfs_cup = getSF(ak.to_numpy(ak.flatten(sjets.hadronFlavour)),ak.to_numpy(ak.flatten(sjets.btagDeepFlavC/(1.-sjets.btagDeepFlavB))),ak.to_numpy(ak.flatten(sjets.btagDeepFlavC/(sjets.btagDeepFlavC+sjets.btagDeepFlavB))),deepjetc_sf,"TotalUncUp")*genweiev
            # csvsfs_cup = getSF(ak.to_numpy(ak.flatten(sjets.hadronFlavour)),ak.to_numpy(ak.flatten(sjets.btagDeepC/(1.-sjets.btagDeepB))),ak.to_numpy(ak.flatten(sjets.btagDeepC/(sjets.btagDeepC+sjets.btagDeepB))),deepcsvc_sf,"TotalUncUp")*genweiev
            # jetsfs_cdn = getSF(ak.to_numpy(ak.flatten(sjets.hadronFlavour)),ak.to_numpy(ak.flatten(sjets.btagDeepFlavC/(1.-sjets.btagDeepFlavB))),ak.to_numpy(ak.flatten(sjets.btagDeepFlavC/(sjets.btagDeepFlavC+sjets.btagDeepFlavB))),deepjetc_sf,"TotalUncDown")*genweiev
            # csvsfs_cdn = getSF(ak.to_numpy(ak.flatten(sjets.hadronFlavour)),ak.to_numpy(ak.flatten(sjets.btagDeepC/(1.-sjets.btagDeepB))),ak.to_numpy(ak.flatten(sjets.btagDeepC/(sjets.btagDeepC+sjets.btagDeepB))),deepcsvc_sf,"TotalUncDown")*genweiev
            jetsfs_b=ak.flatten(ak.broadcast_arrays(weights.weight()[event_level]*ak.fill_none(ak.prod(deepjetb_sf.eval('central',sjets.hadronFlavour,abs(sjets.eta),sjets.pt,discr=sjets.btagDeepFlavB),axis=-1),1.),sjets['pt'])[0])
            csvsfs_b=ak.flatten(ak.broadcast_arrays(weights.weight()[event_level]*ak.fill_none(ak.prod(deepcsvb_sf.eval('central',sjets.hadronFlavour,abs(sjets.eta),sjets.pt,discr=sjets.btagDeepB),axis=-1),1.),sjets['pt'])[0])
            csvsfs_bup=ak.flatten(ak.broadcast_arrays(weights.weight()[event_level]*ak.fill_none(ak.prod(deepcsvb_sf.eval('up_jes',sjets.hadronFlavour,abs(sjets.eta),sjets.pt,discr=sjets.btagDeepB),axis=-1),1.),sjets['pt'])[0])
            jetsfs_bup=ak.flatten(ak.broadcast_arrays(weights.weight()[event_level]*ak.fill_none(ak.prod(deepjetb_sf.eval('up_jes',sjets.hadronFlavour,abs(sjets.eta),sjets.pt,discr=sjets.btagDeepFlavB),axis=-1),1.),sjets['pt'])[0])
            csvsfs_bdn=ak.flatten(ak.broadcast_arrays(weights.weight()[event_level]*ak.fill_none(ak.prod(deepcsvb_sf.eval('down_jes',sjets.hadronFlavour,abs(sjets.eta),sjets.pt,discr=sjets.btagDeepB),axis=-1),1.),sjets['pt'])[0])
            jetsfs_bdn=ak.flatten(ak.broadcast_arrays(weights.weight()[event_level]*ak.fill_none(ak.prod(deepjetb_sf.eval('down_jes',sjets.hadronFlavour,abs(sjets.eta),sjets.pt,discr=sjets.btagDeepFlavB),axis=-1),1.),sjets['pt'])[0])
            
                
            
        # # Fill histograms dynamically  
        # for histname, h in output.items():
        #     if histname in self.deepcsv_hists:
        #         if(isRealData):
        #             if histname == 'rawpt':
        #                 h.fill(dataset=dataset,flav=ak.flatten(genflavor), rawpt =ak.flatten(sjets.pt))
        #             else:
        #                 fields = {l: ak.flatten(sjets[l.replace('jet_','')], axis=None) for l in h.fields if l.replace('jet_','') in dir(sjets)}
        #                 h.fill(dataset=dataset,flav=5, **fields)
        #         else:
        #             if histname == 'rawpt':
        #                 genweiev=ak.flatten(ak.broadcast_arrays(weights.weight()[event_level],sjets['pt'])[0])
        #                 h.fill(dataset=dataset,flav=ak.flatten(genflavor), rawpt =ak.flatten(sjets.pt*(1-sjets.rawFactor)) ,weight=genweiev)
        #             else:
        #                 fields = {l: ak.flatten(sjets[histname]) for l in h.fields if l in dir(sjets)}
        #                 genweiev=ak.flatten(ak.broadcast_arrays(weights.weight()[event_level],sjets['pt'])[0])
        #                 h.fill(dataset=dataset,flav=ak.flatten(genflavor), **fields,weight=genweiev)
            
       
        if not isRealData:
            ###Fill no SFs
            output['btagDeepFlavB'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepFlavB=ak.flatten(sjets.btagDeepFlavB),weight=genweiev)
            # output['btagDeepFlavC'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepFlavC=ak.flatten(sjets.btagDeepFlavC),weight=genweiev)
            output['btagDeepB'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepB=ak.flatten(sjets.btagDeepB),weight=genweiev)
            # output['btagDeepC'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepC=ak.flatten(sjets.btagDeepC),weight=genweiev)
            
            ###Fill SFs-nominal
            output['btagDeepFlavBSF'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepFlavBSF=ak.flatten(sjets.btagDeepFlavB),weight=jetsfs_b)
            output['btagDeepBSF'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepBSF=ak.flatten(sjets.btagDeepB),weight=csvsfs_b)
            output['btagDeepFlavBSF_up'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepFlavBSF_up=ak.flatten(sjets.btagDeepFlavB),weight=jetsfs_bup)
            output['btagDeepBSF_up'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepBSF_up=ak.flatten(sjets.btagDeepB),weight=csvsfs_bup)
            output['btagDeepFlavBSF_down'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepFlavBSF_down=ak.flatten(sjets.btagDeepFlavB),weight=jetsfs_bdn)
            output['btagDeepBSF_down'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepBSF_down=ak.flatten(sjets.btagDeepB),weight=csvsfs_bdn)
            # output['btagDeepFlavCSF'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepFlavCSF=ak.flatten(sjets.btagDeepFlavC),weight=jetsfs_c)
            # output['btagDeepFlavBvC'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt),deepflav_BvC=ak.flatten(sjets.btag_CvB_flav),weight=jetsfs_c)
            # output['btagDeepFlavCvL'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt),deepflav_CvL=ak.flatten(sjets.btag_CvL_flav),weight=jetsfs_c)
            # output['btagDeepCSF'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepCSF=ak.flatten(sjets.btagDeepC),weight=csvsfs_c)
            # output['btagDeepBvC'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt),deepcsv_BvC=ak.flatten(sjets.btag_CvB_csv),weight=csvsfs_c)
            # output['btagDeepCvL'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt),deepcsv_CvL=ak.flatten(sjets.btag_CvL_csv),weight=csvsfs_c)
            # output['btagDeepCSF_up'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepCSF_up=ak.flatten(sjets.btagDeepC),weight=csvsfs_cup)
            # output['btagDeepBvC_up'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt),deepcsv_BvC_up= ak.flatten(sjets.btag_CvB_csv),weight=csvsfs_cup)
            # output['btagDeepCvL_up'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt),deepcsv_CvL_up= ak.flatten(sjets.btag_CvL_csv),weight=csvsfs_cup)
            # output['btagDeepCSF_down'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepCSF_down=ak.flatten(sjets.btagDeepC),weight=csvsfs_cdn)
            # output['btagDeepBvC_down'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt),deepcsv_BvC_dn=ak.flatten(sjets.btag_CvB_csv),weight=csvsfs_cdn)
            # output['btagDeepCvL_down'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt),deepcsv_CvL_dn=ak.flatten(sjets.btag_CvL_csv),weight=csvsfs_cdn)
            # output['btagDeepFlavCSF_up'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepFlavCSF_up=ak.flatten(sjets.btagDeepFlavC),weight=jetsfs_cup)
            # output['btagDeepFlavBvC_up'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), deepflav_BvC_up=ak.flatten(sjets.btag_CvB_flav),weight=jetsfs_cup)
            # output['btagDeepFlavCvL_up'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), deepflav_CvL_up=ak.flatten(sjets.btag_CvL_flav),weight=jetsfs_cup)
            # output['btagDeepFlavCSF_down'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepFlavCSF_down=ak.flatten(sjets.btagDeepFlavC),weight=jetsfs_cdn)
            # output['btagDeepFlavBvC_down'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt),deepflav_BvC_dn=ak.flatten(sjets.btag_CvB_flav),weight=jetsfs_cdn)
            # output['btagDeepFlavCvL_down'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt),deepflav_CvL_dn=ak.flatten(sjets.btag_CvL_flav),weight=jetsfs_cdn)
        else:
            output['btagDeepFlavB'].fill(dataset=dataset,flav=0, eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepFlavB=ak.flatten(sjets.btagDeepFlavB))
            # output['btagDeepFlavC'].fill(dataset=dataset,flav=0, eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepFlavC=ak.flatten(sjets.btagDeepFlavC))
            output['btagDeepB'].fill(dataset=dataset,flav=0, eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepB=ak.flatten(sjets.btagDeepB))
            # output['btagDeepC'].fill(dataset=dataset,flav=0, eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepC=ak.flatten(sjets.btagDeepC))
            ###Fill SFs-nominal
            output['btagDeepFlavBSF'].fill(dataset=dataset,flav=0, eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepFlavBSF=ak.flatten(sjets.btagDeepFlavB))
            output['btagDeepBSF'].fill(dataset=dataset,flav=0, eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepBSF=ak.flatten(sjets.btagDeepB))
            output['btagDeepFlavBSF_up'].fill(dataset=dataset,flav=0, eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepFlavBSF_up=ak.flatten(sjets.btagDeepFlavB))
            output['btagDeepBSF_up'].fill(dataset=dataset,flav=0, eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepBSF_up=ak.flatten(sjets.btagDeepB))
            output['btagDeepFlavBSF_down'].fill(dataset=dataset,flav=0, eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepFlavBSF_down=ak.flatten(sjets.btagDeepFlavB))
            output['btagDeepBSF_down'].fill(dataset=dataset,flav=0, eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepBSF_down=ak.flatten(sjets.btagDeepB))
            # output['btagDeepFlavCSF'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepFlavCSF=ak.flatten(sjets.btagDeepFlavC))
            # output['btagDeepFlavBvC'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt),deepflav_BvC=ak.flatten(sjets.btag_CvB_flav))
            # output['btagDeepFlavCvL'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt),deepflav_CvL=ak.flatten(sjets.btag_CvL_flav))
            # output['btagDeepCSF'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepCSF=ak.flatten(sjets.btagDeepC))
            # output['btagDeepBvC'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt),deepcsv_BvC=ak.flatten(sjets.btag_CvB_csv))
            # output['btagDeepCvL'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt),deepcsv_CvL=ak.flatten(sjets.btag_CvL_csv))
            # output['btagDeepCSF_up'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepCSF_up=ak.flatten(sjets.btagDeepC))
            # output['btagDeepBvC_up'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt),deepcsv_BvC_up= ak.flatten(sjets.btag_CvB_csv))
            # output['btagDeepCvL_up'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt),deepcsv_CvL_up= ak.flatten(sjets.btag_CvL_csv))
            # output['btagDeepCSF_down'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepCSF_down=ak.flatten(sjets.btagDeepC))
            # output['btagDeepBvC_down'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt),deepcsv_BvC_dn=ak.flatten(sjets.btag_CvB_csv))
            # output['btagDeepCvL_down'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt),deepcsv_CvL_dn=ak.flatten(sjets.btag_CvL_csv))
            # output['btagDeepFlavCSF_up'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepFlavCSF_up=ak.flatten(sjets.btagDeepFlavC))
            # output['btagDeepFlavBvC_up'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), deepflav_BvC_up=ak.flatten(sjets.btag_CvB_flav))
            # output['btagDeepFlavCvL_up'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), deepflav_CvL_up=ak.flatten(sjets.btag_CvL_flav))
            # output['btagDeepFlavCSF_down'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt), btagDeepFlavCSF_down=ak.flatten(sjets.btagDeepFlavC))
            # output['btagDeepFlavBvC_down'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt),deepflav_BvC_dn=ak.flatten(sjets.btag_CvB_flav))
            # output['btagDeepFlavCvL_down'].fill(dataset=dataset,flav=ak.flatten(genflavor), eta=ak.flatten(sjets.eta),pt=ak.flatten(sjets.pt),deepflav_CvL_dn=ak.flatten(sjets.btag_CvL_flav))

        def flatten(ar): # flatten awkward into a 1d array to hist
            return ak.flatten(ar, axis=None)

        def num(ar):
            return ak.num(ak.fill_none(ar[~ak.is_none(ar)], 0), axis=0)
        if(isRealData):
            output['njet'].fill(dataset=dataset,  njet=flatten(ak.num(sjets)))
            output['nbjet'].fill(dataset=dataset,  nbjet=flatten(ak.num(sbjets)))
            output['ntbjet'].fill(dataset=dataset,  ntbjet=flatten(ak.num(sbjets)))
            output['nmu'].fill(dataset=dataset,   nmu=flatten(ak.num(smu)))
            output['lmupt'].fill(dataset=dataset, lmupt=flatten((smu.pt)))
            output['ljpt'].fill(dataset=dataset, flav=0, ljpt=flatten(sjets[:,0].pt))
            output['sljpt'].fill(dataset=dataset, flav=0, sljpt=flatten(sjets[:,1].pt))
            output['ssljpt'].fill(dataset=dataset, flav=0, ssljpt=flatten(sjets[:,2].pt))
            output['sssljpt'].fill(dataset=dataset, flav=0, sssljpt=flatten(sjets[:,3].pt))
            output['ljrawpt'].fill(dataset=dataset, flav=0, ljrawpt=flatten(sjets[:,0].pt))
            output['sljrawpt'].fill(dataset=dataset, flav=0, sljrawpt=flatten(sjets[:,1].pt))
            output['ssljrawpt'].fill(dataset=dataset, flav=0, ssljrawpt=flatten(sjets[:,2].pt))
            output['sssljrawpt'].fill(dataset=dataset, flav=0, sssljrawpt=flatten(sjets[:,3].pt))
            output['ljdr'].fill(dataset=dataset, flav=0, ljdr=flatten(sjets[:,0].delta_r(smu)))
            output['sljdr'].fill(dataset=dataset, flav=0, sljdr=flatten(sjets[:,1].delta_r(smu)))
            output['ssljdr'].fill(dataset=dataset, flav=0, ssljdr=flatten(sjets[:,2].delta_r(smu)))
            output['sssljdr'].fill(dataset=dataset, flav=0, sssljdr=flatten(sjets[:,3].delta_r(smu)))
            output['met'].fill(dataset=dataset, met=flatten((selev.METFixEE2017.pt)))
        else:
            output['njet'].fill(dataset=dataset,  njet=flatten(ak.num(sjets)),weight=weights.weight()[event_level])
            output['nbjet'].fill(dataset=dataset,  nbjet=flatten(ak.num(sbjets)),weight=weights.weight()[event_level])
            output['ntbjet'].fill(dataset=dataset,  ntbjet=flatten(ak.num(stbjets)),weight=weights.weight()[event_level])
            output['nmu'].fill(dataset=dataset,   nmu=flatten(ak.num(smu)),weight=weights.weight()[event_level])
            output['lmupt'].fill(dataset=dataset, lmupt=flatten((smu.pt)),weight=weights.weight()[event_level])
            output['ljpt'].fill(dataset=dataset, flav=flatten(sjets[:,0].hadronFlavour), ljpt=flatten(sjets[:,0].pt),weight=weights.weight()[event_level])
            output['sljpt'].fill(dataset=dataset, flav=flatten(sjets[:,1].hadronFlavour), sljpt=flatten(sjets[:,1].pt),weight=weights.weight()[event_level])
            output['ssljpt'].fill(dataset=dataset, flav=flatten(sjets[:,2].hadronFlavour), ssljpt=flatten(sjets[:,2].pt),weight=weights.weight()[event_level])
            output['sssljpt'].fill(dataset=dataset, flav=flatten(sjets[:,3].hadronFlavour), sssljpt=flatten(sjets[:,3].pt),weight=weights.weight()[event_level])
            output['ljrawpt'].fill(dataset=dataset, flav=flatten(sjets[:,0].hadronFlavour), ljrawpt=flatten(sjets[:,0].pt*(1-sjets[:,0].rawFactor)),weight=weights.weight()[event_level])
            output['sljrawpt'].fill(dataset=dataset, flav=flatten(sjets[:,1].hadronFlavour), sljrawpt=flatten(sjets[:,1].pt*(1-sjets[:,1].rawFactor)),weight=weights.weight()[event_level])
            output['ssljrawpt'].fill(dataset=dataset, flav=flatten(sjets[:,2].hadronFlavour), ssljrawpt=flatten(sjets[:,2].pt*(1-sjets[:,2].rawFactor)),weight=weights.weight()[event_level])
            output['sssljrawpt'].fill(dataset=dataset, flav=flatten(sjets[:,3].hadronFlavour), sssljrawpt=flatten(sjets[:,3].pt*(1-sjets[:,3].rawFactor)),weight=weights.weight()[event_level])
            output['ljdr'].fill(dataset=dataset, flav=flatten(sjets[:,0].hadronFlavour), ljdr=flatten(sjets[:,0].delta_r(smu)),weight=weights.weight()[event_level])
            output['sljdr'].fill(dataset=dataset, flav=flatten(sjets[:,1].hadronFlavour), sljdr=flatten(sjets[:,1].delta_r(smu)),weight=weights.weight()[event_level])
            output['ssljdr'].fill(dataset=dataset, flav=flatten(sjets[:,2].hadronFlavour), ssljdr=flatten(sjets[:,2].delta_r(smu)),weight=weights.weight()[event_level])
            output['sssljdr'].fill(dataset=dataset, flav=flatten(sjets[:,3].hadronFlavour), sssljdr=flatten(sjets[:,3].delta_r(smu)),weight=weights.weight()[event_level])
            output['met'].fill(dataset=dataset, met=flatten((selev.METFixEE2017.pt)),weight=weights.weight()[event_level])
        return output

    def postprocess(self, accumulator):
        return accumulator
