import gzip
import pickle, os, sys, mplhep as hep, numpy as np

import coffea
from coffea import hist, processor
import awkward as ak
from coffea.analysis_tools import Weights

# from BTVNanoCommissioning.utils.correction import lumiMasks, compiled, eleSFs,muSFs,deepcsvb_sf,deepcsvc_sf,deepjetb_sf,deepjetc_sfmuSFs
# from BTVNanoCommissioning.helpers.definitions import definitions
# from BTVNanoCommissioning.helpers.cTagSFReader import getSF
from coffea.lumi_tools import LumiMask
from coffea.btag_tools import BTagScaleFactor

class NanoProcessor(processor.ProcessorABC):
    # Define histograms   
    def num(ar):
        return ak.num(ak.fill_none(ar[~ak.is_none(ar)], 0), axis=0)
    def __init__(self):#,year='2017',campaign='Rereco17_94X',jec=True):        
        # Define axes
        # Should read axes from NanoAOD config
        dataset_axis = hist.Cat("dataset", "Primary dataset")
        flav_axis = hist.Bin("flav", r"Genflavour",[0,1,4,5,6])
        charge_axis = hist.Bin("char", r"Charge", [-2,0,2])

        # Jet
        jet_pt_axis   = hist.Bin("pt",   r"Jet $p_{T}$ [GeV]", 50, 0, 500)
        jet_eta_axis  = hist.Bin("eta",  r"Jet $\eta$", 25, -2.5, 2.5)
        jet_phi_axis  = hist.Bin("phi",  r"Jet $\phi$", 30, -3, 3)
        jet_mass_axis = hist.Bin("mass", r"Jet $m$ [GeV]", 50, 0, 500)

        # Events

        njet_axis  = hist.Bin("njet",  r"N jets",      [0,1,2,3])
        ## Muons
        hard_lpt_axis = hist.Bin("pt",   r"Hard lepton pt", 40, 0, 200)
        
        hard_leta_axis  = hist.Bin("eta",  r"Lepton $\eta$", 25, -2.5, 2.5)
        hard_lphi_axis  = hist.Bin("phi",  r"Lepton $\phi$", 30, -3, 3)
        hard_liso_axis = hist.Bin("pfRelIso04_all", r"Hard Muon Rel. Iso", 40,0,4.)
        hard_lpt_ratio_axis = hist.Bin("hlptratio", r"Hard $\mu$/Jet $p_{T}$ [GeV]",40,0,1)
        soft_lpt_axis   = hist.Bin("pt",   r"Soft lepton pt", 40, 0, 40)
        soft_leta_axis  = hist.Bin("eta",  r"Lepton $\eta$", 25, -2.5, 2.5)
        soft_lphi_axis  = hist.Bin("phi",  r"Lepton $\phi$", 30, -3, 3)
        soft_liso_axis = hist.Bin("pfRelIso04_all", r"Soft Muon Rel. Iso", 40,0,4.)
        soft_lpt_ratio_axis = hist.Bin("slptratio", r"Soft $\mu$/Jet $p_{T}$ [GeV]",40,0,1)
        l_dxy_axis = hist.Bin("dxy", r"dxy", 20,0,0.002)    
        l_dz_axis = hist.Bin("dz", r"dz", 20,0,0.01)    
        l_sip3d_axis = hist.Bin("dz", r"dz", 20,0,0.2)
        

        ## Z/W
        zmass_axis = hist.Bin("zmass", r"Z Mass", 25, 50,100)
        zpt_axis = hist.Bin("zpt",r"Z $p_{T}$", 25,0,100)
        zeta_axis  = hist.Bin("zeta",  r"Z $\eta$", 25, -2.5, 2.5)
        zphi_axis  = hist.Bin("zphi",  r"Z $\phi$", 30, -3, 3)
        drmumu_axis = hist.Bin("drmumu", r"$\Delta$R($\mu_{soft}$,$\mu_{hard}$)", 25,0,5)
        wmass_axis = hist.Bin("wmass", r"W mass", 25,50,100)
        wpt_axis = hist.Bin("wpt",r"W $p_{T}$", 25,0,100)
        weta_axis  = hist.Bin("weta",  r"W $\eta$", 25, -2.5, 2.5)
        wphi_axis  = hist.Bin("wphi",  r"W $\phi$", 30, -3, 3)
        
        ## MET
        met_axis = hist.Bin("pt", r"MET $p_{T}$", 50, 0,500)
        metphi_axis = hist.Bin("phi",  r"met $\phi$", 30, -3, 3)

        ## Muon jets
        mujet_pt_axis   = hist.Bin("pt",   r"Jet $p_{T}$ [GeV]", 50, 0, 500)
        mujet_eta_axis  = hist.Bin("eta",  r"Jet $\eta$", 25, -2.5, 2.5)
        mujet_phi_axis  = hist.Bin("phi",  r"Jet $\phi$", 30, -3, 3)
        mujet_mass_axis = hist.Bin("mass", r"Jet $m$ [GeV]", 50, 0, 500)
        dr_mujetsoftmu_axis = hist.Bin("drjet_smu", r"$\Delta$R($\mu_{soft}$,j)",25,0,5)
        dr_mujethardmu_axis = hist.Bin("drjet_hmu", r"$\Delta$R($\mu_{hard}$,j)",25,0,5)

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
        # input_names,manual_ranges,bins = definitions()
        # bininfo = dict(zip(input_names,zip(bins,manual_ranges)))
        # for d in deepcsv_list:
        #     binning, ranges = bininfo["Jet_%s"%d]
        #     if ranges[1] is None : ranges[1] = 0.
        #     if ranges[0] is None : ranges[0] = -0.5
        #     deepcsv_axes.append(hist.Bin(d,d,binning,ranges[0],ranges[1]))
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
            'pt'  : hist.Hist("Counts", dataset_axis, flav_axis, charge_axis,jet_pt_axis),
            'eta' : hist.Hist("Counts", dataset_axis, flav_axis,charge_axis,jet_eta_axis),
            'phi' : hist.Hist("Counts", dataset_axis, flav_axis,charge_axis,jet_phi_axis),
            'mass': hist.Hist("Counts", dataset_axis, flav_axis,charge_axis,jet_mass_axis)
        }
        for disc, axis in zip(varlist, btag_axes):
            for i in range(1):
                _hist_sf_dict["%s_%d" %(disc,i)] = hist.Hist("Counts", dataset_axis, flav_axis,charge_axis,axis)
        for deepcsv, axises in zip(deepcsv_list, deepcsv_axes):
             _hist_deepcsv_dict["%s" %(deepcsv)] = hist.Hist("Counts", dataset_axis,flav_axis,charge_axis, axises)
        
        _hist_event_dict = {
            'njet' : hist.Hist("Counts", dataset_axis, njet_axis),
            'hl_pt' : hist.Hist("Counts", dataset_axis, charge_axis ,hard_lpt_axis),
            'hl_eta' : hist.Hist("Counts", dataset_axis, charge_axis ,hard_leta_axis),
            'hl_phi' : hist.Hist("Counts", dataset_axis, charge_axis ,hard_lphi_axis),
            'hl_pfRelIso04_all' : hist.Hist("Counts", dataset_axis, charge_axis ,hard_liso_axis),
            'sl_pt' : hist.Hist("Counts", dataset_axis, charge_axis ,soft_lpt_axis),
            'sl_eta' : hist.Hist("Counts", dataset_axis, charge_axis ,soft_leta_axis),
            'sl_phi' : hist.Hist("Counts", dataset_axis, charge_axis ,soft_lphi_axis),
            'sl_pfRelIso04_all' : hist.Hist("Counts", dataset_axis, charge_axis ,soft_liso_axis),
            'sl_dxy' : hist.Hist("Counts", dataset_axis, charge_axis, l_dxy_axis),
            'sl_dz' : hist.Hist("Counts", dataset_axis, charge_axis, l_dz_axis),
            'sl_sip3d' : hist.Hist("Counts", dataset_axis,charge_axis, l_sip3d_axis),
            'hlptratio': hist.Hist("Counts", dataset_axis, flav_axis,charge_axis,hard_lpt_ratio_axis),
            'slptratio': hist.Hist("Counts", dataset_axis, flav_axis,charge_axis,soft_lpt_ratio_axis),
            'zmass': hist.Hist("Counts",dataset_axis,charge_axis,zmass_axis),
            'zpt': hist.Hist("Counts",dataset_axis,charge_axis,zpt_axis),
            'zeta': hist.Hist("Counts",dataset_axis,charge_axis,zeta_axis),
            'zphi': hist.Hist("Counts",dataset_axis,charge_axis,zphi_axis),
            'wmass': hist.Hist("Counts",dataset_axis,charge_axis,wmass_axis),
            'wpt': hist.Hist("Counts",dataset_axis,charge_axis,wpt_axis),
            'weta': hist.Hist("Counts",dataset_axis,charge_axis,weta_axis),
            'wphi': hist.Hist("Counts",dataset_axis,charge_axis,wphi_axis),
            'drmumu':hist.Hist("Counts",dataset_axis,charge_axis,drmumu_axis),
            'metpt' : hist.Hist("Counts", dataset_axis, met_axis), 
            'metphi' : hist.Hist("Counts", dataset_axis, metphi_axis),
            'drjet_smu':hist.Hist("Counts",dataset_axis,flav_axis,charge_axis,dr_mujetsoftmu_axis),
            'drjet_hmu':hist.Hist("Counts",dataset_axis,flav_axis,charge_axis,dr_mujethardmu_axis),
            }
        with gzip.open("data/PU/Rereco17_94X/94XPUwei_corrections.pkl.gz") as fin:self.compiled = pickle.load(fin)
        self.deepcsvb_sf = BTagScaleFactor("data/BTV/Rereco17_94X/DeepCSV_94XSF_V5_B_F.csv",BTagScaleFactor.RESHAPE,methods='iterativefit,iterativefit,iterativefit')
        self.deepcsvc_sf = "data/BTV/Rereco17_94X/DeepCSV_ctagSF_MiniAOD94X_2017_pTincl_v3_2_interp.root"
        self.deepjetb_sf = BTagScaleFactor("data/BTV/Rereco17_94X/DeepFlavour_94XSF_V4_B_F.csv",BTagScaleFactor.RESHAPE,methods='iterativefit,iterativefit,iterativefit')
        self.deepjetc_sf = "data/BTV/Rereco17_94X/DeepJet_ctagSF_MiniAOD94X_2017_pTincl_v3_2_interp.root"
        self.lumiMasks = {
            '2016': LumiMask('data/lumiMasks/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'),
            '2017': LumiMask('data/lumiMasks/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt'),
            '2018': LumiMask('data/lumiMasks/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'),
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
        # ## Define the CvL, CvB 
        # if not hasattr(events,"btagDeepFlavCvL"): 
        #     events.Jet['btagDeepFlavCvL'] = np.where(((events.Jet.btagDeepFlavC/(1.-events.Jet.btagDeepFlavB))>0)&(events.Jet.pt>15),(events.Jet.btagDeepFlavC/(1.-events.Jet.btagDeepFlavB)),-1)
        #     events.Jet['btagDeepFlavCvB'] = np.where(((events.Jet.btagDeepFlavC/(events.Jet.btagDeepFlavC+events.Jet.btagDeepFlavB))>0)&(events.Jet.pt>15),(events.Jet.btagDeepFlavC/(events.Jet.btagDeepFlavC+events.Jet.btagDeepFlavB)),-1)
        #     events.Jet['btagDeepCvL'] = np.where((events.Jet.btagDeepC>0)&(events.Jet.pt>15),(events.Jet.btagDeepC/(1.-events.Jet.btagDeepB)),-1)
        #     events.Jet['btagDeepCvB'] = np.where((events.Jet.btagDeepC>0)&(events.Jet.pt>15),(events.Jet.btagDeepC/(events.Jet.btagDeepC+events.Jet.btagDeepB)),-1)
        
        if(isRealData):output['sumw'][dataset] += len(events)
        else:output['sumw'][dataset] += ak.sum(events.genWeight)
        # req_lumi=np.ones(len(events), dtype='bool')
        # if(isRealData): req_lumi=self.lumiMasks['2017'](events.run, events.luminosityBlock)
        # weights = Weights(len(events), storeIndividual=True)
        # if not isRealData:
        #     weights.add('genweight',events.genWeight)
        #     weights.add('puweight', self.compiled['2017_pileupweight'](events.Pileup.nPU))
        # ##############
        # # Trigger level
        # triggers = [
        # # "HLT_IsoMu24",
        # "HLT_IsoMu27",
        # ]
        
        # trig_arrs = [events.HLT[_trig.strip("HLT_")] for _trig in triggers]
        # req_trig = np.zeros(len(events), dtype='bool')
        # for t in trig_arrs:
        #     req_trig = req_trig | t

        # ############
        # # Event level
        # ## Muon cuts
        # # muon twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
        # iso_muon = events.Muon[(events.Muon.pt > 30) & (abs(events.Muon.eta) < 2.4) & (events.Muon.tightId > .5)&(events.Muon.pfRelIso04_all<=0.15)]
        
        # req_muon =(ak.count(iso_muon.pt, axis=1) == 1)
        # iso_muon = ak.pad_none(iso_muon,1, axis=1)
        # iso_muon=iso_muon[:,0]

        # dilep_mu = events.Muon[(events.Muon.pt>12)&(abs(events.Muon.eta) < 2.4)& (events.Muon.tightId > .5)&(events.Muon.pfRelIso04_all<=0.15)]
        # dilep_ele = events.Electron[(events.Electron.pt>15)&((abs(events.Electron.eta) < 1.4442)|((abs(events.Electron.eta) < 2.5)&(abs(events.Electron.eta) >1.566)))& (events.Electron.mvaFall17V2Iso_WP80 > .5)]
        # req_dilepveto = (ak.count(dilep_mu.pt,axis=1)+ak.count(dilep_ele.pt,axis=1)!=2)

        # # MET = events.METFixEE2017
        
        # MET =ak.zip({
        #     "pt": events.METFixEE2017.pt,
        #     "eta": ak.zeros_like(events.METFixEE2017.pt),
        #     "phi": events.METFixEE2017.phi,
        #     "mass":ak.zeros_like(events.METFixEE2017.pt),
        #     }, with_name="PtEtaPhiMLorentzVector")

        # Wmass = MET+iso_muon
        # req_Wmass = (Wmass.mass>55)
        
        # ## Jet cuts 
        # event_jet =  events.Jet[(events.Jet.pt > 20) & (abs(events.Jet.eta) <= 2.5)&(((events.Jet.puId >=7)&( events.Jet.pt<50))|(events.Jet.pt>=50)) &(events.Jet.jetId>=3)&(events.Jet.btagDeepB>0.) & (events.Jet.btagDeepB<1.) & (events.Jet.btagDeepC>0.) & (events.Jet.btagDeepC<1.) & (events.Jet.btagDeepFlavB>0.) & (events.Jet.btagDeepFlavB<1.) & (events.Jet.btagDeepFlavC>0.) & (events.Jet.btagDeepFlavC<1.)& (ak.all(events.Jet.metric_table(iso_muon) > 0.5, axis=2))&((events.Jet.muEF+events.Jet.neEmEF)<0.7)]

        # req_jets = ((ak.num(event_jet.puId) >=1) & (ak.num(event_jet.puId) <=3))
        
        # ## Soft Muon cuts 
        
        # soft_muon = events.Muon[(events.Muon.pt < 25) & (abs(events.Muon.eta) < 2.4) & (events.Muon.tightId > .5)&(events.Muon.pfRelIso04_all>0.2)&(events.Muon.jetIdx!=-1)]
        # req_softmu=(ak.count(soft_muon.pt, axis=1) >= 1)
        
        
        # mu_jet = event_jet[(ak.all(event_jet.metric_table(soft_muon) <= 0.4, axis=2))&((event_jet.muonIdx1!=-1)|(event_jet.muonIdx2!=-1))]
        # req_mujet = (ak.num(mu_jet.pt,axis=1)>=1)
        # # mu_jet = ak.fill_none(mu_jet.pt,0)
        # mu_jet = ak.pad_none(mu_jet,1, axis=1)
        # soft_muon= ak.pad_none(soft_muon,1,axis=1)
        # # pT ratio
        # req_pTratio = ((soft_muon[:,0].pt/mu_jet[:,0].pt)<0.4)
        # #dilepton mass
        # # req_dilepmass = np.zeros(len(events), dtype='bool')
        # # req_dilepmass = req_muon&req_softmu
        # dilep_mass = iso_muon+soft_muon[:,0]
        # req_dilepmass = ((dilep_mass.mass>12.)&((dilep_mass.mass<80)| (dilep_mass.mass>100)))

        # req_QCDveto = (iso_muon.pfRelIso04_all<0.05) & (abs(iso_muon.dz)<0.01) & (abs(iso_muon.dxy)<0.002) & (iso_muon.ip3d < 0.2) & ((iso_muon.pt/mu_jet[:,0].pt<0.)|(iso_muon.pt/mu_jet[:,0].pt>0.75))
        # event_level = req_trig & req_lumi & req_muon &  req_jets & req_softmu   & req_dilepmass &req_mujet & req_Wmass & req_dilepveto & req_QCDveto & req_pTratio
        # if(len(event_level)>0):event_level = ak.fill_none(event_level,False)        
        # # Selected
        # selev = events[event_level]    
        
        # #########
        
        # ## Hard Muon
        # shmu = selev.Muon[(selev.Muon.pt > 30) & (abs(selev.Muon.eta) < 2.4) & (selev.Muon.tightId > .5)&(selev.Muon.pfRelIso04_all<=0.15)]
       
        # shmu=shmu[:,0]
        
        # sjets = selev.Jet[(selev.Jet.pt > 20) & (abs(selev.Jet.eta) <= 2.5)&(((selev.Jet.puId >=7)&( selev.Jet.pt<50))|(selev.Jet.pt>=50)) &(selev.Jet.jetId>=3)&(selev.Jet.btagDeepB>0.) & (selev.Jet.btagDeepB<1.) & (selev.Jet.btagDeepC>0.) & (selev.Jet.btagDeepC<1.) & (selev.Jet.btagDeepFlavB>0.) & (selev.Jet.btagDeepFlavB<1.) & (selev.Jet.btagDeepFlavC>0.) & (selev.Jet.btagDeepFlavC<1.)& (ak.all(selev.Jet.metric_table(shmu) > 0.5, axis=2))&((selev.Jet.muEF+selev.Jet.neEmEF)<0.7)]
        # ## Soft Muon
        # ssmu = selev.Muon[(selev.Muon.pt < 25) & (abs(selev.Muon.eta) < 2.4) & (selev.Muon.tightId > .5)&(selev.Muon.pfRelIso04_all>0.2)&(selev.Muon.jetIdx!=-1)]
        # # ssmu=ssmu[:,0]
        # ## MET
        # smet =ak.zip({
        #     "pt": selev.METFixEE2017.pt,
        #     "eta": ak.zeros_like(selev.METFixEE2017.pt),
        #     "phi": selev.METFixEE2017.phi,
        #     "mass":ak.zeros_like(selev.METFixEE2017.pt),
        #     }, with_name="PtEtaPhiMLorentzVector")
       
        # ## Jets
        
        
        
       
        # ## Muon Jet 
        # # print(sjets.pt)
        # smuon_jet = sjets[(ak.all(sjets.metric_table(ssmu) <= 0.4, axis=2))&((sjets.muonIdx1!=-1)|(sjets.muonIdx2!=-1))]
        # # print(smuon_jet.pt)
        # smuon_jet = smuon_jet[:,0]
        # ssmu = ssmu[:,0]
        # sz=shmu+ssmu
        # sw=shmu+smet
        # osss=shmu.charge*ssmu.charge*-1
        # # if not isRealData:
        # #     weights.add('lep1sf',np.where(event_level,muSFs(ak.firsts(events.Muon[(events.Muon.pt>12)&(abs(events.Muon.eta) < 2.4)& (events.Muon.tightId > .5)&(events.Muon.pfRelIso04_all<=0.15)])),1.))
        # #     weights.add('lep2sf',np.where(event_level,muSFs(ak.firsts(events.Muon[(events.Muon.pt < 25) & (abs(events.Muon.eta) < 2.4) & (events.Muon.tightId > .5)&(events.Muon.pfRelIso04_all>0.2)])),1.))
        # njet= ak.count(sjets.pt,axis=1)   
        # def flatten(ar): # flatten awkward into a 1d array to hist
        #     return ak.flatten(ar, axis=None)
        # if isRealData :
        #     genflavor = ak.zeros_like(sjets.pt)
        # else:
        #     par_flav = (smuon_jet.partonFlavour == 0 ) & (smuon_jet.hadronFlavour==0)
        #     genflavor = smuon_jet.hadronFlavour + 1*par_flav 
        #     # jetsfs_c_lj = getSF(ak.to_numpy(smuon_jet.hadronFlavour),ak.to_numpy(smuon_jet.btagDeepFlavCvL),ak.to_numpy(smuon_jet.btagDeepFlavCvB),self.deepjetc_sf)
        #     # jetsfs_c_up_lj = getSF(ak.to_numpy(smuon_jet.hadronFlavour),ak.to_numpy(smuon_jet.btagDeepFlavCvL),ak.to_numpy(smuon_jet.btagDeepFlavCvB),self.deepjetc_sf,"TotalUncUp")
        #     # jetsfs_c_dn_lj = getSF(ak.to_numpy(smuon_jet.hadronFlavour),ak.to_numpy(smuon_jet.btagDeepFlavCvL),ak.to_numpy(smuon_jet.btagDeepFlavCvB),self.deepjetc_sf,"TotalUncDown")
        #     jetsfs_b_lj = self.deepjetb_sf.eval('central',smuon_jet.hadronFlavour,abs(smuon_jet.eta),smuon_jet.pt,discr=smuon_jet.btagDeepFlavB)
        #     jetsfs_b_up_lj = self.deepjetb_sf.eval('up_jes',smuon_jet.hadronFlavour,abs(smuon_jet.eta),smuon_jet.pt,discr=smuon_jet.btagDeepFlavB)
        #     jetsfs_b_dn_lj = self.deepjetb_sf.eval('down_jes',smuon_jet.hadronFlavour,abs(smuon_jet.eta),smuon_jet.pt,discr=smuon_jet.btagDeepFlavB)
        #     # csvsfs_c_lj = getSF(ak.to_numpy(smuon_jet.hadronFlavour),ak.to_numpy(smuon_jet.btagDeepCvL),ak.to_numpy(smuon_jet.btagDeepCvB),self.deepcsvc_sf)
        #     # csvsfs_c_up_lj = getSF(ak.to_numpy(smuon_jet.hadronFlavour),ak.to_numpy(smuon_jet.btagDeepCvL),ak.to_numpy(smuon_jet.btagDeepCvB),self.deepcsvc_sf,"TotalUncUp")
        #     # csvsfs_c_dn_lj = getSF(ak.to_numpy(smuon_jet.hadronFlavour),ak.to_numpy(ak.to_numpy(smuon_jet.btagDeepCvL)),ak.to_numpy(smuon_jet.btagDeepCvB),self.deepcsvc_sf,"TotalUncDown")
        #     csvsfs_b_up_lj = self.deepcsvb_sf.eval('up_jes',smuon_jet.hadronFlavour,abs(smuon_jet.eta),smuon_jet.pt,discr=smuon_jet.btagDeepB)
        #     csvsfs_b_lj = self.deepcsvb_sf.eval('central',smuon_jet.hadronFlavour,abs(smuon_jet.eta),smuon_jet.pt,discr=smuon_jet.btagDeepB)
        #     csvsfs_b_dn_lj = self.deepcsvb_sf.eval('down_jes',smuon_jet.hadronFlavour,abs(smuon_jet.eta),smuon_jet.pt,discr=smuon_jet.btagDeepB)
            
            
        # for histname, h in output.items():
        #     # if histname in self.deepcsv_hists:
        #     #     #ch = ak.flatten(ak.broadcast_arrays(osss,sjets['pt'])[0])
        #     #     fields = {l: smuon_jet[histname] for l in h.fields if l in dir(smuon_jet)}
        #     #     if(isRealData):
        #     #         h.fill(dataset=dataset,flav=5, char= osss,**fields,weight=osss)
        #     #     else:                    
        #     #         h.fill(dataset=dataset,flav=genflavor,char=osss, **fields,weight=weights.weight()[event_level]*osss)
        #     if 'hl_' in histname:
        #         fields = {l: shmu[l] for l in h.fields if l in dir(shmu)}
        #         if(isRealData): h.fill(dataset=dataset,char=osss,**fields,weight=osss)
        #         else: h.fill(dataset=dataset,char=osss,**fields,weight=weights.weight()[event_level]*osss)
        #     elif 'sl_' in histname:
        #         fields = {l: ssmu[l] for l in h.fields if l in dir(ssmu)}
        #         if(isRealData): h.fill(dataset=dataset,char=osss,**fields,weight=osss)
        #         else: h.fill(dataset=dataset,char=osss,**fields,weight=weights.weight()[event_level]*osss)
        #     elif 'met' in histname: 
        #         fields = {l: selev.METFixEE2017[l] for l in h.fields if l in dir(selev.METFixEE2017)}
        #         if(isRealData): h.fill(dataset=dataset,**fields,weight=osss)
        #         else: h.fill(dataset=dataset,**fields,weight=weights.weight()[event_level]*osss)
            
        # if not isRealData:
        #     ###Fill no SFs
        #     output['zmass'].fill(dataset=dataset,char=osss,zmass=sz.mass,weight=weights.weight()[event_level]*osss)
        #     output['zpt'].fill(dataset=dataset,char=osss,zpt=sz.pt,weight=weights.weight()[event_level]*osss)
        #     output['zeta'].fill(dataset=dataset,char=osss,zeta=sz.eta,weight=weights.weight()[event_level]*osss)
        #     output['zphi'].fill(dataset=dataset,char=osss,zphi=sz.phi,weight=weights.weight()[event_level]*osss)
        #     output['wmass'].fill(dataset=dataset,char=osss,wmass=sw.mass,weight=weights.weight()[event_level]*osss)
        #     output['wpt'].fill(dataset=dataset,char=osss,wpt=sw.pt,weight=weights.weight()[event_level]*osss)
        #     output['weta'].fill(dataset=dataset,char=osss,weta=sw.eta,weight=weights.weight()[event_level]*osss)
        #     output['wphi'].fill(dataset=dataset,char=osss,wphi=sw.phi,weight=weights.weight()[event_level]*osss)
        #     output['drmumu'].fill(dataset=dataset,char=osss,drmumu=ssmu.delta_r(shmu),weight=weights.weight()[event_level]*osss)
        #     output['hlptratio'].fill(dataset=dataset,flav=genflavor,char=osss,hlptratio=shmu.pt/smuon_jet.pt,weight=weights.weight()[event_level]*osss)
        #     output['slptratio'].fill(dataset=dataset,flav=genflavor,char=osss,slptratio=ssmu.pt/smuon_jet.pt,weight=weights.weight()[event_level]*osss)
        #     output['drjet_hmu'].fill(dataset=dataset,flav=genflavor,char=osss,drjet_hmu=smuon_jet.delta_r(shmu),weight=weights.weight()[event_level]*osss)
        #     output['drjet_smu'].fill(dataset=dataset,flav=genflavor,char=osss,drjet_smu=smuon_jet.delta_r(ssmu),weight=weights.weight()[event_level]*osss)

        #     ## discri
            
        #     output['btagDeepFlavB_0'].fill(dataset=dataset,flav=genflavor,char=osss, btagDeepFlavB=smuon_jet.btagDeepFlavB,weight=weights.weight()[event_level]*osss)
        #     output['btagDeepFlavC_0'].fill(dataset=dataset,flav=genflavor,char=osss, btagDeepFlavC=smuon_jet.btagDeepFlavC,weight=weights.weight()[event_level]*osss)
        #     output['btagDeepB_0'].fill(dataset=dataset,flav=genflavor,char=osss, btagDeepB=smuon_jet.btagDeepB,weight=weights.weight()[event_level]*osss)
        #     output['btagDeepC_0'].fill(dataset=dataset,flav=genflavor,char=osss, btagDeepC=smuon_jet.btagDeepC,weight=weights.weight()[event_level]*osss)
        #     output['deepcsv_CvB_0'].fill(dataset=dataset,flav=genflavor,char=osss, deepcsv_CvB=np.where(smuon_jet.btagDeepCvB<0,-0.2,smuon_jet.btagDeepCvB),weight=weights.weight()[event_level]*osss)
        #     output['deepcsv_CvL_0'].fill(dataset=dataset,flav=genflavor,char=osss, deepcsv_CvL=np.where(smuon_jet.btagDeepCvL<0,-0.2,smuon_jet.btagDeepCvL),weight=weights.weight()[event_level]*osss)
        #     output['deepflav_CvB_0'].fill(dataset=dataset,flav=genflavor,char=osss, deepflav_CvB=np.where(smuon_jet.btagDeepFlavCvB<0,-0.2,smuon_jet.btagDeepFlavCvB),weight=weights.weight()[event_level]*osss)
        #     output['deepflav_CvL_0'].fill(dataset=dataset,flav=genflavor,char=osss, deepflav_CvL=np.where(smuon_jet.btagDeepFlavCvL<0,-0.2,smuon_jet.btagDeepFlavCvL),weight=weights.weight()[event_level]*osss)
        #     # output['btagDeepFlavBSF_0'].fill(dataset=dataset,flav=genflavor,char=osss, btagDeepFlavBSF=smuon_jet.btagDeepFlavB,weight=(weights.weight()[event_level]*jetsfs_b_lj))      
        #     # output['btagDeepFlavCSF_0'].fill(dataset=dataset,flav=genflavor,char=osss, btagDeepFlavCSF=smuon_jet.btagDeepFlavC,weight=weights.weight()[event_level]*osss*jetsfs_c_lj)
        #     # output['btagDeepBSF_0'].fill(dataset=dataset,flav=genflavor,char=osss, btagDeepBSF=smuon_jet.btagDeepB,weight=weights.weight()[event_level]*osss*csvsfs_b_lj)
        #     # output['btagDeepCSF_0'].fill(dataset=dataset,flav=genflavor,char=osss, btagDeepCSF=smuon_jet.btagDeepC,weight=weights.weight()[event_level]*osss*csvsfs_c_lj)
        #     # output['deepcsv_CvBSF_0'].fill(dataset=dataset,flav=genflavor,char=osss, deepcsv_CvBSF=np.where(smuon_jet.btagDeepCvB<0,-0.2,smuon_jet.btagDeepCvB),weight=weights.weight()[event_level]*osss*csvsfs_c_lj)
        #     # output['deepcsv_CvLSF_0'].fill(dataset=dataset,flav=genflavor,char=osss, deepcsv_CvLSF=np.where(smuon_jet.btagDeepCvL<0,-0.2,smuon_jet.btagDeepCvL),weight=weights.weight()[event_level]*osss*csvsfs_c_lj)
        #     # output['deepflav_CvBSF_0'].fill(dataset=dataset,flav=genflavor,char=osss, deepflav_CvBSF=np.where(smuon_jet.btagDeepFlavCvB<0,-0.2,smuon_jet.btagDeepFlavCvB),weight=weights.weight()[event_level]*osss*jetsfs_c_lj)
        #     # output['deepflav_CvLSF_0'].fill(dataset=dataset,flav=genflavor,char=osss, deepflav_CvLSF=np.where(smuon_jet.btagDeepFlavCvL<0,-0.2,smuon_jet.btagDeepFlavCvL),weight=weights.weight()[event_level]*osss*jetsfs_c_lj)
        #     # output['btagDeepFlavB_up_0'].fill(dataset=dataset,flav=genflavor,char=osss, btagDeepFlavB_up=smuon_jet.btagDeepFlavB,weight=weights.weight()[event_level]*osss*jetsfs_b_up_lj)
        #     # output['btagDeepFlavC_up_0'].fill(dataset=dataset,flav=genflavor,char=osss, btagDeepFlavC_up=smuon_jet.btagDeepFlavC,weight=weights.weight()[event_level]*osss*jetsfs_c_up_lj)
        #     # output['btagDeepB_up_0'].fill(dataset=dataset,flav=genflavor,char=osss, btagDeepB_up=smuon_jet.btagDeepB,weight=weights.weight()[event_level]*osss*csvsfs_b_up_lj)
        #     # output['btagDeepC_up_0'].fill(dataset=dataset,flav=genflavor,char=osss, btagDeepC_up=smuon_jet.btagDeepC,weight=weights.weight()[event_level]*osss*csvsfs_c_up_lj)
        #     # output['deepcsv_CvB_up_0'].fill(dataset=dataset,flav=genflavor,char=osss, deepcsv_CvB_up=np.where(smuon_jet.btagDeepCvB<0,-0.2,smuon_jet.btagDeepCvB),weight=weights.weight()[event_level]*osss*csvsfs_c_up_lj)
        #     # output['deepcsv_CvL_up_0'].fill(dataset=dataset,flav=genflavor,char=osss, deepcsv_CvL_up=np.where(smuon_jet.btagDeepCvL<0,-0.2,smuon_jet.btagDeepCvL),weight=weights.weight()[event_level]*osss*csvsfs_c_up_lj)
        #     # output['deepflav_CvB_up_0'].fill(dataset=dataset,flav=genflavor,char=osss, deepflav_CvB_up=np.where(smuon_jet.btagDeepFlavCvB<0,-0.2,smuon_jet.btagDeepFlavCvB),weight=weights.weight()[event_level]*osss*jetsfs_c_up_lj)
        #     # output['deepflav_CvL_up_0'].fill(dataset=dataset,flav=genflavor,char=osss, deepflav_CvL_up=np.where(smuon_jet.btagDeepFlavCvL<0,-0.2,smuon_jet.btagDeepFlavCvL),weight=weights.weight()[event_level]*osss*jetsfs_c_up_lj)
        #     # output['btagDeepFlavB_dn_0'].fill(dataset=dataset,flav=genflavor,char=osss, btagDeepFlavB_dn=smuon_jet.btagDeepFlavB,weight=weights.weight()[event_level]*osss*jetsfs_b_dn_lj)
        #     # output['btagDeepFlavC_dn_0'].fill(dataset=dataset,flav=genflavor,char=osss, btagDeepFlavC_dn=smuon_jet.btagDeepFlavC,weight=weights.weight()[event_level]*osss*jetsfs_c_dn_lj)
        #     # output['btagDeepB_dn_0'].fill(dataset=dataset,flav=genflavor,char=osss, btagDeepB_dn=smuon_jet.btagDeepB,weight=weights.weight()[event_level]*osss*csvsfs_b_dn_lj)
        #     # output['btagDeepC_dn_0'].fill(dataset=dataset,flav=genflavor,char=osss, btagDeepC_dn=smuon_jet.btagDeepC,weight=weights.weight()[event_level]*osss*csvsfs_c_dn_lj)
        #     # output['deepcsv_CvB_dn_0'].fill(dataset=dataset,flav=genflavor,char=osss, deepcsv_CvB_dn=np.where(smuon_jet.btagDeepCvB<0,-0.2,smuon_jet.btagDeepCvB),weight=weights.weight()[event_level]*osss*csvsfs_c_dn_lj)
        #     # output['deepcsv_CvL_dn_0'].fill(dataset=dataset,flav=genflavor,char=osss, deepcsv_CvL_dn=np.where(smuon_jet.btagDeepCvL<0,-0.2,smuon_jet.btagDeepCvL),weight=weights.weight()[event_level]*osss*csvsfs_c_dn_lj)
        #     # output['deepflav_CvB_dn_0'].fill(dataset=dataset,flav=genflavor,char=osss, deepflav_CvB_dn=np.where(smuon_jet.btagDeepFlavCvB<0,-0.2,smuon_jet.btagDeepFlavCvB),weight=weights.weight()[event_level]*osss*jetsfs_c_dn_lj)
        #     # output['deepflav_CvL_dn_0'].fill(dataset=dataset,flav=genflavor,char=osss, deepflav_CvL_dn=np.where(smuon_jet.btagDeepFlavCvL<0,-0.2,smuon_jet.btagDeepFlavCvL),weight=weights.weight()[event_level]*osss*jetsfs_c_dn_lj)
            
        # else:
        #     ###Fill no SFs
        #     output['zmass'].fill(dataset=dataset,char=osss,zmass=sz.mass,weight=osss)
        #     output['zpt'].fill(dataset=dataset,char=osss,zpt=sz.pt,weight=osss)
        #     output['zeta'].fill(dataset=dataset,char=osss,zeta=sz.eta,weight=osss)
        #     output['zphi'].fill(dataset=dataset,char=osss,zphi=sz.phi,weight=osss)
        #     output['wmass'].fill(dataset=dataset,char=osss,wmass=sw.mass,weight=osss)
        #     output['wpt'].fill(dataset=dataset,char=osss,wpt=sw.pt,weight=osss)
        #     output['weta'].fill(dataset=dataset,char=osss,weta=sw.eta,weight=osss)
        #     output['wphi'].fill(dataset=dataset,char=osss,wphi=sw.phi,weight=osss)
        #     output['drmumu'].fill(dataset=dataset,char=osss,drmumu=ssmu.delta_r(shmu),weight=osss)
        #     output['hlptratio'].fill(dataset=dataset,flav=5,char=osss,hlptratio=shmu.pt/smuon_jet.pt,weight=osss)
        #     output['slptratio'].fill(dataset=dataset,flav=5,char=osss,slptratio=ssmu.pt/smuon_jet.pt,weight=osss)
        #     output['drjet_hmu'].fill(dataset=dataset,flav=5,char=osss,drjet_hmu=smuon_jet.delta_r(shmu),weight=osss)
        #     output['drjet_smu'].fill(dataset=dataset,flav=5,char=osss,drjet_smu=smuon_jet.delta_r(ssmu),weight=osss)
                
        #         ## discri
        #     output['btagDeepFlavB_0'].fill(dataset=dataset,flav=5,char=osss, btagDeepFlavB=smuon_jet.btagDeepFlavB,weight=osss)
        #     output['btagDeepFlavC_0'].fill(dataset=dataset,flav=5,char=osss, btagDeepFlavC=smuon_jet.btagDeepFlavC,weight=osss)
        #     output['btagDeepB_0'].fill(dataset=dataset,flav=5,char=osss, btagDeepB=smuon_jet.btagDeepB,weight=osss)
        #     output['btagDeepC_0'].fill(dataset=dataset,flav=5,char=osss, btagDeepC=smuon_jet.btagDeepC,weight=osss)
        #     output['deepcsv_CvB_0'].fill(dataset=dataset,flav=5,char=osss, deepcsv_CvB=np.where(smuon_jet.btagDeepCvB<0,-0.2,smuon_jet.btagDeepCvB),weight=osss)
        #     output['deepcsv_CvL_0'].fill(dataset=dataset,flav=5,char=osss, deepcsv_CvL=np.where(smuon_jet.btagDeepCvL<0,-0.2,smuon_jet.btagDeepCvL),weight=osss)
        #     output['deepflav_CvB_0'].fill(dataset=dataset,flav=5,char=osss, deepflav_CvB=np.where(smuon_jet.btagDeepFlavCvB<0,-0.2,smuon_jet.btagDeepFlavCvB),weight=osss)
        #     output['deepflav_CvL_0'].fill(dataset=dataset,flav=5,char=osss, deepflav_CvL=np.where(smuon_jet.btagDeepFlavCvL<0,-0.2,smuon_jet.btagDeepFlavCvL),weight=osss)
        #     output['btagDeepFlavBSF_0'].fill(dataset=dataset,flav=5,char=osss, btagDeepFlavBSF=smuon_jet.btagDeepFlavB,weight=osss)
        #     output['btagDeepFlavCSF_0'].fill(dataset=dataset,flav=5,char=osss, btagDeepFlavCSF=smuon_jet.btagDeepFlavC,weight=osss)
        #     output['btagDeepBSF_0'].fill(dataset=dataset,flav=5,char=osss, btagDeepBSF=smuon_jet.btagDeepB,weight=osss)
        #     output['btagDeepCSF_0'].fill(dataset=dataset,flav=5,char=osss, btagDeepCSF=smuon_jet.btagDeepC,weight=osss)
        #     output['deepcsv_CvBSF_0'].fill(dataset=dataset,flav=5,char=osss, deepcsv_CvBSF=np.where(smuon_jet.btagDeepCvB<0,-0.2,smuon_jet.btagDeepCvB),weight=osss)
        #     output['deepcsv_CvLSF_0'].fill(dataset=dataset,flav=5,char=osss, deepcsv_CvLSF=np.where(smuon_jet.btagDeepCvL<0,-0.2,smuon_jet.btagDeepCvL),weight=osss)
        #     output['deepflav_CvBSF_0'].fill(dataset=dataset,flav=5,char=osss, deepflav_CvBSF=np.where(smuon_jet.btagDeepFlavCvB<0,-0.2,smuon_jet.btagDeepFlavCvB),weight=osss)
        #     output['deepflav_CvLSF_0'].fill(dataset=dataset,flav=5,char=osss, deepflav_CvLSF=np.where(smuon_jet.btagDeepFlavCvL<0,-0.2,smuon_jet.btagDeepFlavCvL),weight=osss)
        #     output['btagDeepFlavB_up_0'].fill(dataset=dataset,flav=5,char=osss, btagDeepFlavB_up=smuon_jet.btagDeepFlavB,weight=osss)
        #     output['btagDeepFlavC_up_0'].fill(dataset=dataset,flav=5,char=osss, btagDeepFlavC_up=smuon_jet.btagDeepFlavC,weight=osss)
        #     output['btagDeepB_up_0'].fill(dataset=dataset,flav=5,char=osss, btagDeepB_up=smuon_jet.btagDeepB,weight=osss)
        #     output['btagDeepC_up_0'].fill(dataset=dataset,flav=5,char=osss, btagDeepC_up=smuon_jet.btagDeepC,weight=osss)
        #     output['deepcsv_CvB_up_0'].fill(dataset=dataset,flav=5,char=osss, deepcsv_CvB_up=np.where(smuon_jet.btagDeepCvB<0,-0.2,smuon_jet.btagDeepCvB),weight=osss)
        #     output['deepcsv_CvL_up_0'].fill(dataset=dataset,flav=5,char=osss, deepcsv_CvL_up=np.where(smuon_jet.btagDeepCvL<0,-0.2,smuon_jet.btagDeepCvL),weight=osss)
        #     output['deepflav_CvB_up_0'].fill(dataset=dataset,flav=5,char=osss, deepflav_CvB_up=np.where(smuon_jet.btagDeepFlavCvB<0,-0.2,smuon_jet.btagDeepFlavCvB),weight=osss)
        #     output['deepflav_CvL_up_0'].fill(dataset=dataset,flav=5,char=osss, deepflav_CvL_up=np.where(smuon_jet.btagDeepFlavCvL<0,-0.2,smuon_jet.btagDeepFlavCvL),weight=osss)
        #     output['btagDeepFlavB_dn_0'].fill(dataset=dataset,flav=5,char=osss, btagDeepFlavB_dn=smuon_jet.btagDeepFlavB,weight=osss)
        #     output['btagDeepFlavC_dn_0'].fill(dataset=dataset,flav=5,char=osss, btagDeepFlavC_dn=smuon_jet.btagDeepFlavC,weight=osss)
        #     output['btagDeepB_dn_0'].fill(dataset=dataset,flav=5,char=osss, btagDeepB_dn=smuon_jet.btagDeepB,weight=osss)
        #     output['btagDeepC_dn_0'].fill(dataset=dataset,flav=5,char=osss, btagDeepC_dn=smuon_jet.btagDeepC,weight=osss)
        #     output['deepcsv_CvB_dn_0'].fill(dataset=dataset,flav=5,char=osss, deepcsv_CvB_dn=np.where(smuon_jet.btagDeepCvB<0,-0.2,smuon_jet.btagDeepCvB),weight=osss)
        #     output['deepcsv_CvL_dn_0'].fill(dataset=dataset,flav=5,char=osss, deepcsv_CvL_dn=np.where(smuon_jet.btagDeepCvL<0,-0.2,smuon_jet.btagDeepCvL),weight=osss)
        #     output['deepflav_CvB_dn_0'].fill(dataset=dataset,flav=5,char=osss, deepflav_CvB_dn=np.where(smuon_jet.btagDeepFlavCvB<0,-0.2,smuon_jet.btagDeepFlavCvB),weight=osss)
        #     output['deepflav_CvL_dn_0'].fill(dataset=dataset,flav=5,char=osss, deepflav_CvL_dn=np.where(smuon_jet.btagDeepFlavCvL<0,-0.2,smuon_jet.btagDeepFlavCvL),weight=osss)
        # print(output)    
        return output

    def postprocess(self, accumulator):
        return accumulator
