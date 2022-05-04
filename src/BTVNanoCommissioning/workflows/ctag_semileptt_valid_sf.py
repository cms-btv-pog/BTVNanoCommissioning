import gzip
import pickle, os, sys, mplhep as hep, numpy as np

from matplotlib.pyplot import jet

import coffea
from coffea import hist, processor
import awkward as ak
from coffea.analysis_tools import Weights
import gc
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
        input_names,manual_ranges,bins = definitions()
        bininfo = dict(zip(input_names,zip(bins,manual_ranges)))
        for d in deepcsv_list:
            binning, ranges = bininfo["Jet_%s"%d]
            if ranges[1] is None : ranges[1] = 0.
            if ranges[0] is None : ranges[0] = -0.5
            print(binning,ranges[0],ranges[1])
            deepcsv_axes.append(hist.Bin(d,d,binning,ranges[0],ranges[1]))
        # Define similar axes dynamically
        disc_list = ['btagDeepB', 'btagDeepC', 'btagDeepFlavB', 'btagDeepFlavC','deepcsv_CvL','deepcsv_CvB','deepflav_CvL','deepflav_CvB']
        syst_list = ['','SF','_up','_dn']
        varlist=[]
        btag_axes = []
        for d in disc_list:
            for s in syst_list:
                btag_axes.append(hist.Bin("%s%s" %(d,s), "%s%s" %(d,s), 30, -0.2, 1))  
                varlist.append("%s%s" %(d,s))
                
        _hist_sf_dict={}   
        _hist_deepcsv_dict={
            'pt'  : hist.Hist("Counts", dataset_axis, flav_axis,jet_pt_axis),
            'eta' : hist.Hist("Counts", dataset_axis, flav_axis,jet_eta_axis),
            'phi' : hist.Hist("Counts", dataset_axis, flav_axis,jet_phi_axis),
            'mass': hist.Hist("Counts", dataset_axis, flav_axis,jet_mass_axis)
        }
        for disc, axis in zip(varlist, btag_axes):
            for i in range(1):
                _hist_sf_dict["%s_%d" %(disc,i)] = hist.Hist("Counts", dataset_axis, flav_axis,axis)
        for deepcsv, axises in zip(deepcsv_list, deepcsv_axes):
             _hist_deepcsv_dict["%s" %(deepcsv)] = hist.Hist("Counts", dataset_axis,flav_axis, axises)
        
        _hist_event_dict = {
            'njet' : hist.Hist("Counts", dataset_axis, njet_axis),
            'hl_pt' : hist.Hist("Counts", dataset_axis, hard_lpt_axis),
            'hl_eta' : hist.Hist("Counts", dataset_axis, hard_leta_axis),
            'hl_phi' : hist.Hist("Counts", dataset_axis, hard_lphi_axis),
            'hl_pfRelIso04_all' : hist.Hist("Counts", dataset_axis, hard_liso_axis),
            'sl_pt' : hist.Hist("Counts", dataset_axis, soft_lpt_axis),
            'sl_eta' : hist.Hist("Counts", dataset_axis, soft_leta_axis),
            'sl_phi' : hist.Hist("Counts", dataset_axis, soft_lphi_axis),
            'sl_pfRelIso04_all' : hist.Hist("Counts", dataset_axis, soft_liso_axis),
            'sl_dxy' : hist.Hist("Counts", dataset_axis,  l_dxy_axis),
            'sl_dz' : hist.Hist("Counts", dataset_axis,  l_dz_axis),
            'sl_sip3d' : hist.Hist("Counts", dataset_axis, l_sip3d_axis),
            'hlptratio': hist.Hist("Counts", dataset_axis, flav_axis,hard_lpt_ratio_axis),
            'slptratio': hist.Hist("Counts", dataset_axis, flav_axis,soft_lpt_ratio_axis),
            'zmass': hist.Hist("Counts",dataset_axis,zmass_axis),
            'zpt': hist.Hist("Counts",dataset_axis,zpt_axis),
            'zeta': hist.Hist("Counts",dataset_axis,zeta_axis),
            'zphi': hist.Hist("Counts",dataset_axis,zphi_axis),
            'wmass': hist.Hist("Counts",dataset_axis,wmass_axis),
            'wpt': hist.Hist("Counts",dataset_axis,wpt_axis),
            'weta': hist.Hist("Counts",dataset_axis,weta_axis),
            'wphi': hist.Hist("Counts",dataset_axis,wphi_axis),
            'drmumu':hist.Hist("Counts",dataset_axis,drmumu_axis),
            'metpt' : hist.Hist("Counts", dataset_axis, met_axis), 
            'metphi' : hist.Hist("Counts", dataset_axis, metphi_axis),
            'mujet_pt': hist.Hist("Counts",dataset_axis,flav_axis,mujet_pt_axis),
            'mujet_eta': hist.Hist("Counts",dataset_axis,flav_axis,mujet_eta_axis),
            'mujet_phi': hist.Hist("Counts",dataset_axis,flav_axis,mujet_phi_axis),
            'mujet_mass': hist.Hist("Counts",dataset_axis,flav_axis,mujet_mass_axis),
            'drjet_smu':hist.Hist("Counts",dataset_axis,flav_axis,dr_mujetsoftmu_axis),
            'drjet_hmu':hist.Hist("Counts",dataset_axis,flav_axis,dr_mujethardmu_axis),
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
        
        req_muon =(ak.count(iso_muon.pt, axis=1) == 1)
        iso_muon = ak.pad_none(iso_muon,1, axis=1)
        iso_muon=iso_muon[:,0]

        dilep_mu = events.Muon[(events.Muon.pt>12)&(abs(events.Muon.eta) < 2.4)& (events.Muon.tightId > .5)&(events.Muon.pfRelIso04_all<=0.15)]
        dilep_ele = events.Electron[(events.Electron.pt>15)&((abs(events.Electron.eta) < 1.4442)|((abs(events.Electron.eta) < 2.5)&(abs(events.Electron.eta) >1.566)))& (events.Electron.mvaFall17V2Iso_WP80 > .5)]
        req_dilepveto = (ak.count(dilep_mu.pt,axis=1)+ak.count(dilep_ele.pt,axis=1)!=2)

        # MET = events.METFixEE2017
        
        MET =ak.zip({
            "pt": events.METFixEE2017.pt,
            "eta": ak.zeros_like(events.METFixEE2017.pt),
            "phi": events.METFixEE2017.phi,
            "mass":ak.zeros_like(events.METFixEE2017.pt),
            }, with_name="PtEtaPhiMLorentzVector")

        Wmass = MET+iso_muon
        req_Wmass = (Wmass.mass>55)
        
        ## Jet cuts 
        event_jet =  events.Jet[(events.Jet.pt > 20) & (abs(events.Jet.eta) <= 2.5)&(((events.Jet.puId >=7)&( events.Jet.pt<50))|(events.Jet.pt>=50)) &(events.Jet.jetId>=3)&(events.Jet.btagDeepB>0.) & (events.Jet.btagDeepB<1.) & (events.Jet.btagDeepC>0.) & (events.Jet.btagDeepC<1.) & (events.Jet.btagDeepFlavB>0.) & (events.Jet.btagDeepFlavB<1.) & (events.Jet.btagDeepFlavC>0.) & (events.Jet.btagDeepFlavC<1.)& (ak.all(events.Jet.metric_table(iso_muon) > 0.5, axis=2))&((events.Jet.muEF+events.Jet.neEmEF)<0.7)]

        req_jets = (ak.num(event_jet.puId) >=4)
        
        ## Soft Muon cuts 
        
        soft_muon = events.Muon[(events.Muon.pt < 25) & (abs(events.Muon.eta) < 2.4) & (events.Muon.tightId > .5)&(events.Muon.pfRelIso04_all>0.2)&(events.Muon.jetIdx!=-1)]
        req_softmu=(ak.count(soft_muon.pt, axis=1) >= 1)
        
        
        mu_jet = event_jet[(ak.all(event_jet.metric_table(soft_muon) <= 0.4, axis=2))&((event_jet.muonIdx1!=-1)|(event_jet.muonIdx2!=-1))]
        req_mujet = (ak.num(mu_jet.pt,axis=1)>=1)
        # mu_jet = ak.fill_none(mu_jet.pt,0)
        mu_jet = ak.pad_none(mu_jet,1, axis=1)
        soft_muon= ak.pad_none(soft_muon,1,axis=1)
        # pT ratio
        req_pTratio = ((soft_muon[:,0].pt/mu_jet[:,0].pt)<0.4)
        #dilepton mass
        # req_dilepmass = np.zeros(len(events), dtype='bool')
        # req_dilepmass = req_muon&req_softmu
        dilep_mass = iso_muon+soft_muon[:,0]
        req_dilepmass = ((dilep_mass.mass>12.)&((dilep_mass.mass<80)| (dilep_mass.mass>100)))

        req_QCDveto = (iso_muon.pfRelIso04_all<0.05) & (abs(iso_muon.dz)<0.01) & (abs(iso_muon.dxy)<0.002) & (iso_muon.ip3d < 0.2) & ((iso_muon.pt/mu_jet[:,0].pt<0.)|(iso_muon.pt/mu_jet[:,0].pt>0.75))
        event_level = req_trig & req_lumi & req_muon &  req_jets & req_softmu   & req_dilepmass &req_mujet & req_Wmass & req_dilepveto & req_QCDveto & req_pTratio
        if(len(event_level)>0):event_level = ak.fill_none(event_level,False)        
        # Selected
        selev = events[event_level]    
        #########
        
        ## Hard Muon
        shmu = selev.Muon[(selev.Muon.pt > 30) & (abs(selev.Muon.eta) < 2.4) & (selev.Muon.tightId > .5)&(selev.Muon.pfRelIso04_all<=0.15)]
       
        shmu=shmu[:,0]
        
        sjets = selev.Jet[(selev.Jet.pt > 20) & (abs(selev.Jet.eta) <= 2.5)&(((selev.Jet.puId >=7)&( selev.Jet.pt<50))|(selev.Jet.pt>=50)) &(selev.Jet.jetId>=3)&(selev.Jet.btagDeepB>0.) & (selev.Jet.btagDeepB<1.) & (selev.Jet.btagDeepC>0.) & (selev.Jet.btagDeepC<1.) & (selev.Jet.btagDeepFlavB>0.) & (selev.Jet.btagDeepFlavB<1.) & (selev.Jet.btagDeepFlavC>0.) & (selev.Jet.btagDeepFlavC<1.)& (ak.all(selev.Jet.metric_table(shmu) > 0.5, axis=2))&((selev.Jet.muEF+selev.Jet.neEmEF)<0.7)]
        ## Soft Muon
        ssmu = selev.Muon[(selev.Muon.pt < 25) & (abs(selev.Muon.eta) < 2.4) & (selev.Muon.tightId > .5)&(selev.Muon.pfRelIso04_all>0.2)&(selev.Muon.jetIdx!=-1)]
        # ssmu=ssmu[:,0]
       
        ## MET
        smet =ak.zip({
            "pt": selev.METFixEE2017.pt,
            "eta": ak.zeros_like(selev.METFixEE2017.pt),
            "phi": selev.METFixEE2017.phi,
            "mass":ak.zeros_like(selev.METFixEE2017.pt),
            }, with_name="PtEtaPhiMLorentzVector")
       
        ## Jets
        
        
        
       
        ## Muon Jet 
        # print(sjets.pt)
        smuon_jet = sjets[(ak.all(sjets.metric_table(ssmu) <= 0.4, axis=2))&((sjets.muonIdx1!=-1)|(sjets.muonIdx2!=-1))]
        # print(smuon_jet.pt)
        smuon_jet = smuon_jet[:,0]
        ssmu = ssmu[:,0]
        sz=shmu+ssmu
        sw=shmu+smet
        osss=shmu.charge*ssmu.charge*-1
        njet= ak.count(sjets.pt,axis=1)   
        def flatten(ar): # flatten awkward into a 1d array to hist
            return ak.flatten(ar, axis=None)
        if isRealData :
            genflavor = ak.zeros_like(sjets.pt)
        else:
            par_flav = (sjets.partonFlavour == 0 ) & (sjets.hadronFlavour==0)
            genflavor = sjets.hadronFlavour + 1*par_flav 
            # genweiev=ak.flatten(ak.broadcast_arrays(weights.weight()[event_level],sjets['pt'])[0])
            jetsfs_c_lj = getSF(ak.to_numpy(smuon_jet.hadronFlavour),ak.to_numpy(smuon_jet.btagDeepFlavC/(1.-smuon_jet.btagDeepFlavB)),ak.to_numpy(smuon_jet.btagDeepFlavC/(smuon_jet.btagDeepFlavC+smuon_jet.btagDeepFlavB)),deepjetc_sf)
            jetsfs_c_up_lj = getSF(ak.to_numpy(smuon_jet.hadronFlavour),ak.to_numpy(smuon_jet.btagDeepFlavC/(1.-smuon_jet.btagDeepFlavB)),ak.to_numpy(smuon_jet.btagDeepFlavC/(smuon_jet.btagDeepFlavC+smuon_jet.btagDeepFlavB)),deepjetc_sf,"TotalUncUp")
            jetsfs_c_dn_lj = getSF(ak.to_numpy(smuon_jet.hadronFlavour),ak.to_numpy(smuon_jet.btagDeepFlavC/(1.-smuon_jet.btagDeepFlavB)),ak.to_numpy(smuon_jet.btagDeepFlavC/(smuon_jet.btagDeepFlavC+smuon_jet.btagDeepFlavB)),deepjetc_sf,"TotalUncDown")
            jetsfs_b_lj = deepjetb_sf.eval('central',smuon_jet.hadronFlavour,abs(smuon_jet.eta),smuon_jet.pt,discr=smuon_jet.btagDeepFlavB)
            jetsfs_b_up_lj = deepjetb_sf.eval('up_jes',smuon_jet.hadronFlavour,abs(smuon_jet.eta),smuon_jet.pt,discr=smuon_jet.btagDeepFlavB)
            jetsfs_b_dn_lj = deepjetb_sf.eval('down_jes',smuon_jet.hadronFlavour,abs(smuon_jet.eta),smuon_jet.pt,discr=smuon_jet.btagDeepFlavB)
            csvsfs_c_lj = getSF(ak.to_numpy(smuon_jet.hadronFlavour),ak.to_numpy(smuon_jet.btagDeepC/(1.-smuon_jet.btagDeepB)),ak.to_numpy(smuon_jet.btagDeepC/(smuon_jet.btagDeepC+smuon_jet.btagDeepB)),deepcsvc_sf)
            csvsfs_c_up_lj = getSF(ak.to_numpy(smuon_jet.hadronFlavour),ak.to_numpy(smuon_jet.btagDeepC/(1.-smuon_jet.btagDeepB)),ak.to_numpy(smuon_jet.btagDeepC/(smuon_jet.btagDeepC+smuon_jet.btagDeepB)),deepcsvc_sf,"TotalUncUp")
            csvsfs_c_dn_lj = getSF(ak.to_numpy(smuon_jet.hadronFlavour),ak.to_numpy(smuon_jet.btagDeepC/(1.-smuon_jet.btagDeepB)),ak.to_numpy(smuon_jet.btagDeepC/(smuon_jet.btagDeepC+smuon_jet.btagDeepB)),deepcsvc_sf,"TotalUncDown")
            csvsfs_b_up_lj = deepcsvb_sf.eval('up_jes',smuon_jet.hadronFlavour,abs(smuon_jet.eta),smuon_jet.pt,discr=smuon_jet.btagDeepB)
            csvsfs_b_lj = deepcsvb_sf.eval('central',smuon_jet.hadronFlavour,abs(smuon_jet.eta),smuon_jet.pt,discr=smuon_jet.btagDeepB)
            csvsfs_b_dn_lj = deepcsvb_sf.eval('down_jes',smuon_jet.hadronFlavour,abs(smuon_jet.eta),smuon_jet.pt,discr=smuon_jet.btagDeepB)
            
            
        for histname, h in output.items():
            if 'sumw' in histname:print('sumw')
            elif histname in self.deepcsv_hists:
                ch = ak.flatten(ak.broadcast_arrays(osss,sjets['pt'])[0])
                fields = {l: ak.flatten(sjets[histname]) for l in h.fields if l in dir(sjets)}
                if(isRealData):
                    h.fill(dataset=dataset,flav=5, **fields,weight=ch)
                else:                    
                    genweiev=ak.flatten(ak.broadcast_arrays(weights.weight()[event_level],sjets['pt'])[0])
                    h.fill(dataset=dataset,flav=ak.flatten(genflavor),**fields,weight=genweiev*ch)
            elif 'hl_' in histname:
                fields = {l: shmu[l] for l in h.fields if l in dir(shmu)}
                if(isRealData): h.fill(dataset=dataset,**fields)
                else: h.fill(dataset=dataset,**fields,weight=weights.weight()[event_level])
            elif 'sl_' in histname:
                fields = {l: ssmu[l] for l in h.fields if l in dir(ssmu)}
                if(isRealData): h.fill(dataset=dataset,**fields)
                else: h.fill(dataset=dataset,**fields,weight=weights.weight()[event_level])
            elif 'met' in histname: 
                fields = {l: selev.METFixEE2017[l] for l in h.fields if l in dir(selev.METFixEE2017)}
                if(isRealData): h.fill(dataset=dataset,**fields)
                else: h.fill(dataset=dataset,**fields,weight=weights.weight()[event_level])
            elif 'mujet_' in histname:
                fields = {l: smuon_jet[l] for l in h.fields if l in dir(smuon_jet)}
                if isRealData: h.fill(dataset=dataset,flav=5,  **fields)
                else :h.fill(dataset=dataset,flav=5,  **fields,weight = weights.weight()[event_level])
        
        if not isRealData:
            ###Fill no SFs
            smpu= (smuon_jet.partonFlavour == 0 ) & (smuon_jet.hadronFlavour==0)
            genflavor= 1*smpu+smuon_jet.hadronFlavour
            output['zmass'].fill(dataset=dataset,zmass=sz.mass,weight=weights.weight()[event_level])
            output['zpt'].fill(dataset=dataset,zpt=sz.pt,weight=weights.weight()[event_level])
            output['zeta'].fill(dataset=dataset,zeta=sz.eta,weight=weights.weight()[event_level])
            output['zphi'].fill(dataset=dataset,zphi=sz.phi,weight=weights.weight()[event_level])
            output['wmass'].fill(dataset=dataset,wmass=sw.mass,weight=weights.weight()[event_level])
            output['wpt'].fill(dataset=dataset,wpt=sw.pt,weight=weights.weight()[event_level])
            output['weta'].fill(dataset=dataset,weta=sw.eta,weight=weights.weight()[event_level])
            output['wphi'].fill(dataset=dataset,wphi=sw.phi,weight=weights.weight()[event_level])
            output['drmumu'].fill(dataset=dataset,drmumu=ssmu.delta_r(shmu),weight=weights.weight()[event_level])
            output['hlptratio'].fill(dataset=dataset,flav=genflavor,hlptratio=shmu.pt/smuon_jet.pt,weight=weights.weight()[event_level])
            output['slptratio'].fill(dataset=dataset,flav=genflavor,slptratio=ssmu.pt/smuon_jet.pt,weight=weights.weight()[event_level])
            output['drjet_hmu'].fill(dataset=dataset,flav=genflavor,drjet_hmu=smuon_jet.delta_r(shmu),weight=weights.weight()[event_level])
            output['drjet_smu'].fill(dataset=dataset,flav=genflavor,drjet_smu=smuon_jet.delta_r(ssmu),weight=weights.weight()[event_level])

            ## discri
            
            output['btagDeepFlavB_0'].fill(dataset=dataset,flav=genflavor, btagDeepFlavB=smuon_jet.btagDeepFlavB,weight=weights.weight()[event_level])
            output['btagDeepFlavC_0'].fill(dataset=dataset,flav=genflavor, btagDeepFlavC=smuon_jet.btagDeepFlavC,weight=weights.weight()[event_level])
            output['btagDeepB_0'].fill(dataset=dataset,flav=genflavor, btagDeepB=smuon_jet.btagDeepB,weight=weights.weight()[event_level])
            output['btagDeepC_0'].fill(dataset=dataset,flav=genflavor, btagDeepC=smuon_jet.btagDeepC,weight=weights.weight()[event_level])
            output['deepcsv_CvB_0'].fill(dataset=dataset,flav=genflavor, deepcsv_CvB=smuon_jet.btagDeepC/(1.-smuon_jet.btagDeepB),weight=weights.weight()[event_level])
            output['deepcsv_CvL_0'].fill(dataset=dataset,flav=genflavor, deepcsv_CvL=smuon_jet.btagDeepC/(smuon_jet.btagDeepC+smuon_jet.btagDeepB),weight=weights.weight()[event_level])
            output['deepflav_CvB_0'].fill(dataset=dataset,flav=genflavor, deepflav_CvB=smuon_jet.btagDeepFlavC/(1.-smuon_jet.btagDeepFlavB),weight=weights.weight()[event_level])
            output['deepflav_CvL_0'].fill(dataset=dataset,flav=genflavor, deepflav_CvL=smuon_jet.btagDeepFlavC/(smuon_jet.btagDeepFlavC+smuon_jet.btagDeepFlavB),weight=weights.weight()[event_level])
            output['btagDeepFlavBSF_0'].fill(dataset=dataset,flav=genflavor, btagDeepFlavBSF=smuon_jet.btagDeepFlavB,weight=(weights.weight()[event_level]*jetsfs_b_lj))      
            output['btagDeepFlavCSF_0'].fill(dataset=dataset,flav=genflavor, btagDeepFlavCSF=smuon_jet.btagDeepFlavC,weight=weights.weight()[event_level]*jetsfs_c_lj)
            output['btagDeepBSF_0'].fill(dataset=dataset,flav=genflavor, btagDeepBSF=smuon_jet.btagDeepB,weight=weights.weight()[event_level]*csvsfs_b_lj)
            output['btagDeepCSF_0'].fill(dataset=dataset,flav=genflavor, btagDeepCSF=smuon_jet.btagDeepC,weight=weights.weight()[event_level]*csvsfs_c_lj)
            output['deepcsv_CvBSF_0'].fill(dataset=dataset,flav=genflavor, deepcsv_CvBSF=smuon_jet.btagDeepC/(1.-smuon_jet.btagDeepB),weight=weights.weight()[event_level]*csvsfs_c_lj)
            output['deepcsv_CvLSF_0'].fill(dataset=dataset,flav=genflavor, deepcsv_CvLSF=smuon_jet.btagDeepC/(smuon_jet.btagDeepC+smuon_jet.btagDeepB),weight=weights.weight()[event_level]*csvsfs_c_lj)
            output['deepflav_CvBSF_0'].fill(dataset=dataset,flav=genflavor, deepflav_CvBSF=smuon_jet.btagDeepFlavC/(1.-smuon_jet.btagDeepFlavB),weight=weights.weight()[event_level]*jetsfs_c_lj)
            output['deepflav_CvLSF_0'].fill(dataset=dataset,flav=genflavor, deepflav_CvLSF=smuon_jet.btagDeepFlavC/(smuon_jet.btagDeepFlavC+smuon_jet.btagDeepFlavB),weight=weights.weight()[event_level]*jetsfs_c_lj)
            output['btagDeepFlavB_up_0'].fill(dataset=dataset,flav=genflavor, btagDeepFlavB_up=smuon_jet.btagDeepFlavB,weight=weights.weight()[event_level]*jetsfs_b_up_lj)
            output['btagDeepFlavC_up_0'].fill(dataset=dataset,flav=genflavor, btagDeepFlavC_up=smuon_jet.btagDeepFlavC,weight=weights.weight()[event_level]*jetsfs_c_up_lj)
            output['btagDeepB_up_0'].fill(dataset=dataset,flav=genflavor, btagDeepB_up=smuon_jet.btagDeepB,weight=weights.weight()[event_level]*csvsfs_b_up_lj)
            output['btagDeepC_up_0'].fill(dataset=dataset,flav=genflavor, btagDeepC_up=smuon_jet.btagDeepC,weight=weights.weight()[event_level]*csvsfs_c_up_lj)
            output['deepcsv_CvB_up_0'].fill(dataset=dataset,flav=genflavor, deepcsv_CvB_up=smuon_jet.btagDeepC/(1.-smuon_jet.btagDeepB),weight=weights.weight()[event_level]*csvsfs_c_up_lj)
            output['deepcsv_CvL_up_0'].fill(dataset=dataset,flav=genflavor, deepcsv_CvL_up=smuon_jet.btagDeepC/(smuon_jet.btagDeepC+smuon_jet.btagDeepB),weight=weights.weight()[event_level]*csvsfs_c_up_lj)
            output['deepflav_CvB_up_0'].fill(dataset=dataset,flav=genflavor, deepflav_CvB_up=smuon_jet.btagDeepFlavC/(1.-smuon_jet.btagDeepFlavB),weight=weights.weight()[event_level]*jetsfs_c_up_lj)
            output['deepflav_CvL_up_0'].fill(dataset=dataset,flav=genflavor, deepflav_CvL_up=smuon_jet.btagDeepFlavC/(smuon_jet.btagDeepFlavC+smuon_jet.btagDeepFlavB),weight=weights.weight()[event_level]*jetsfs_c_up_lj)
            output['btagDeepFlavB_dn_0'].fill(dataset=dataset,flav=genflavor, btagDeepFlavB_dn=smuon_jet.btagDeepFlavB,weight=weights.weight()[event_level]*jetsfs_b_dn_lj)
            output['btagDeepFlavC_dn_0'].fill(dataset=dataset,flav=genflavor, btagDeepFlavC_dn=smuon_jet.btagDeepFlavC,weight=weights.weight()[event_level]*jetsfs_c_dn_lj)
            output['btagDeepB_dn_0'].fill(dataset=dataset,flav=genflavor, btagDeepB_dn=smuon_jet.btagDeepB,weight=weights.weight()[event_level]*csvsfs_b_dn_lj)
            output['btagDeepC_dn_0'].fill(dataset=dataset,flav=genflavor, btagDeepC_dn=smuon_jet.btagDeepC,weight=weights.weight()[event_level]*csvsfs_c_dn_lj)
            output['deepcsv_CvB_dn_0'].fill(dataset=dataset,flav=genflavor, deepcsv_CvB_dn=smuon_jet.btagDeepC/(1.-smuon_jet.btagDeepB),weight=weights.weight()[event_level]*csvsfs_c_dn_lj)
            output['deepcsv_CvL_dn_0'].fill(dataset=dataset,flav=genflavor, deepcsv_CvL_dn=smuon_jet.btagDeepC/(smuon_jet.btagDeepC+smuon_jet.btagDeepB),weight=weights.weight()[event_level]*csvsfs_c_dn_lj)
            output['deepflav_CvB_dn_0'].fill(dataset=dataset,flav=genflavor, deepflav_CvB_dn=smuon_jet.btagDeepFlavC/(1.-smuon_jet.btagDeepFlavB),weight=weights.weight()[event_level]*jetsfs_c_dn_lj)
            output['deepflav_CvL_dn_0'].fill(dataset=dataset,flav=genflavor, deepflav_CvL_dn=smuon_jet.btagDeepFlavC/(smuon_jet.btagDeepFlavC+smuon_jet.btagDeepFlavB),weight=weights.weight()[event_level]*jetsfs_c_dn_lj)
            
        else:
            ###Fill no SFs
            output['zmass'].fill(dataset=dataset,zmass=sz.mass)
            output['zpt'].fill(dataset=dataset,zpt=sz.pt)
            output['zeta'].fill(dataset=dataset,zeta=sz.eta)
            output['zphi'].fill(dataset=dataset,zphi=sz.phi)
            output['wmass'].fill(dataset=dataset,wmass=sw.mass)
            output['wpt'].fill(dataset=dataset,wpt=sw.pt)
            output['weta'].fill(dataset=dataset,weta=sw.eta)
            output['wphi'].fill(dataset=dataset,wphi=sw.phi)
            output['drmumu'].fill(dataset=dataset,drmumu=ssmu.delta_r(shmu))
            output['hlptratio'].fill(dataset=dataset,flav=5,hlptratio=shmu.pt/smuon_jet.pt)
            output['slptratio'].fill(dataset=dataset,flav=5,slptratio=ssmu.pt/smuon_jet.pt)
            output['drjet_hmu'].fill(dataset=dataset,flav=5,drjet_hmu=smuon_jet.delta_r(shmu))
            output['drjet_smu'].fill(dataset=dataset,flav=5,drjet_smu=smuon_jet.delta_r(ssmu))
                
                ## discri
            output['btagDeepFlavB_0'].fill(dataset=dataset,flav=5, btagDeepFlavB=smuon_jet.btagDeepFlavB)
            output['btagDeepFlavC_0'].fill(dataset=dataset,flav=5, btagDeepFlavC=smuon_jet.btagDeepFlavC)
            output['btagDeepB_0'].fill(dataset=dataset,flav=5, btagDeepB=smuon_jet.btagDeepB)
            output['btagDeepC_0'].fill(dataset=dataset,flav=5, btagDeepC=smuon_jet.btagDeepC)
            output['deepcsv_CvB_0'].fill(dataset=dataset,flav=5, deepcsv_CvB=smuon_jet.btagDeepC/(1.-smuon_jet.btagDeepB))
            output['deepcsv_CvL_0'].fill(dataset=dataset,flav=5, deepcsv_CvL=smuon_jet.btagDeepC/(smuon_jet.btagDeepC+smuon_jet.btagDeepB))
            output['deepflav_CvB_0'].fill(dataset=dataset,flav=5, deepflav_CvB=smuon_jet.btagDeepFlavC/(1.-smuon_jet.btagDeepFlavB))
            output['deepflav_CvL_0'].fill(dataset=dataset,flav=5, deepflav_CvL=smuon_jet.btagDeepFlavC/(smuon_jet.btagDeepFlavC+smuon_jet.btagDeepFlavB))
            output['btagDeepFlavBSF_0'].fill(dataset=dataset,flav=5, btagDeepFlavBSF=smuon_jet.btagDeepFlavB)
            output['btagDeepFlavCSF_0'].fill(dataset=dataset,flav=5, btagDeepFlavCSF=smuon_jet.btagDeepFlavC)
            output['btagDeepBSF_0'].fill(dataset=dataset,flav=5, btagDeepBSF=smuon_jet.btagDeepB)
            output['btagDeepCSF_0'].fill(dataset=dataset,flav=5, btagDeepCSF=smuon_jet.btagDeepC)
            output['deepcsv_CvBSF_0'].fill(dataset=dataset,flav=5, deepcsv_CvBSF=smuon_jet.btagDeepC/(1.-smuon_jet.btagDeepB))
            output['deepcsv_CvLSF_0'].fill(dataset=dataset,flav=5, deepcsv_CvLSF=smuon_jet.btagDeepC/(smuon_jet.btagDeepC+smuon_jet.btagDeepB))
            output['deepflav_CvBSF_0'].fill(dataset=dataset,flav=5, deepflav_CvBSF=smuon_jet.btagDeepFlavC/(1.-smuon_jet.btagDeepFlavB))
            output['deepflav_CvLSF_0'].fill(dataset=dataset,flav=5, deepflav_CvLSF=smuon_jet.btagDeepFlavC/(smuon_jet.btagDeepFlavC+smuon_jet.btagDeepFlavB))
            output['btagDeepFlavB_up_0'].fill(dataset=dataset,flav=5, btagDeepFlavB_up=smuon_jet.btagDeepFlavB)
            output['btagDeepFlavC_up_0'].fill(dataset=dataset,flav=5, btagDeepFlavC_up=smuon_jet.btagDeepFlavC)
            output['btagDeepB_up_0'].fill(dataset=dataset,flav=5, btagDeepB_up=smuon_jet.btagDeepB)
            output['btagDeepC_up_0'].fill(dataset=dataset,flav=5, btagDeepC_up=smuon_jet.btagDeepC)
            output['deepcsv_CvB_up_0'].fill(dataset=dataset,flav=5, deepcsv_CvB_up=smuon_jet.btagDeepC/(1.-smuon_jet.btagDeepB))
            output['deepcsv_CvL_up_0'].fill(dataset=dataset,flav=5, deepcsv_CvL_up=smuon_jet.btagDeepC/(smuon_jet.btagDeepC+smuon_jet.btagDeepB))
            output['deepflav_CvB_up_0'].fill(dataset=dataset,flav=5, deepflav_CvB_up=smuon_jet.btagDeepFlavC/(1.-smuon_jet.btagDeepFlavB))
            output['deepflav_CvL_up_0'].fill(dataset=dataset,flav=5, deepflav_CvL_up=smuon_jet.btagDeepFlavC/(smuon_jet.btagDeepFlavC+smuon_jet.btagDeepFlavB))
            output['btagDeepFlavB_dn_0'].fill(dataset=dataset,flav=5, btagDeepFlavB_dn=smuon_jet.btagDeepFlavB)
            output['btagDeepFlavC_dn_0'].fill(dataset=dataset,flav=5, btagDeepFlavC_dn=smuon_jet.btagDeepFlavC)
            output['btagDeepB_dn_0'].fill(dataset=dataset,flav=5, btagDeepB_dn=smuon_jet.btagDeepB)
            output['btagDeepC_dn_0'].fill(dataset=dataset,flav=5, btagDeepC_dn=smuon_jet.btagDeepC)
            output['deepcsv_CvB_dn_0'].fill(dataset=dataset,flav=5, deepcsv_CvB_dn=smuon_jet.btagDeepC/(1.-smuon_jet.btagDeepB))
            output['deepcsv_CvL_dn_0'].fill(dataset=dataset,flav=5, deepcsv_CvL_dn=smuon_jet.btagDeepC/(smuon_jet.btagDeepC+smuon_jet.btagDeepB))
            output['deepflav_CvB_dn_0'].fill(dataset=dataset,flav=5, deepflav_CvB_dn=smuon_jet.btagDeepFlavC/(1.-smuon_jet.btagDeepFlavB))
            output['deepflav_CvL_dn_0'].fill(dataset=dataset,flav=5, deepflav_CvL_dn=smuon_jet.btagDeepFlavC/(smuon_jet.btagDeepFlavC+smuon_jet.btagDeepFlavB))
            
        
         
        gc.collect()

        return output

    def postprocess(self, accumulator):
        return accumulator
