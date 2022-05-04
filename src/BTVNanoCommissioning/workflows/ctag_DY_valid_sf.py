import gzip
import pickle, os, sys, mplhep as hep, numpy as np

from matplotlib.pyplot import jet

import coffea
from coffea import hist, processor
import awkward as ak
from coffea.analysis_tools import Weights

from BTVNanoCommissioning.helpers.definitions import definitions
from BTVNanoCommissioning.utils.correction import lumiMasks, compiled, eleSFs,muSFs,deepcsvb_sf,deepcsvc_sf,deepjetb_sf,deepjetc_sf



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
        soft_lpt_axis   = hist.Bin("pt",   r"Soft lepton pt", 40, 0, 40)
        soft_leta_axis  = hist.Bin("eta",  r"Lepton $\eta$", 25, -2.5, 2.5)
        soft_lphi_axis  = hist.Bin("phi",  r"Lepton $\phi$", 30, -3, 3)
        soft_liso_axis = hist.Bin("pfRelIso04_all", r"Soft Muon Rel. Iso", 40,0,4.)
        l_dxy_axis = hist.Bin("dxy", r"dxy", 20,0,0.002)    
        l_dz_axis = hist.Bin("dz", r"dz", 20,0,0.01)    
        l_sip3d_axis = hist.Bin("dz", r"dz", 20,0,0.2)
        

        ## Z
        zmass_axis = hist.Bin("zmass", r"Z Mass", 25, 50,100)
        zpt_axis = hist.Bin("zpt",r"Z $p_{T}$", 25,0,100)
        zeta_axis  = hist.Bin("zeta",  r"Z $\eta$", 25, -2.5, 2.5)
        zphi_axis  = hist.Bin("zphi",  r"Z $\phi$", 30, -3, 3)
        drmumu_axis = hist.Bin("drmumu", r"$\Delta$R($\mu_{soft}$,$\mu_{hard}$)", 25,0,5)

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
                btag_axes.append(hist.Bin("%s%s" %(d,s), "%s%s" %(d,s), 51, -0.2, 1))  
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
            'pos_pt' : hist.Hist("Counts", dataset_axis ,hard_lpt_axis),
            'pos_eta' : hist.Hist("Counts", dataset_axis ,hard_leta_axis),
            'pos_phi' : hist.Hist("Counts", dataset_axis ,hard_lphi_axis),
            'pos_pfRelIso04_all' : hist.Hist("Counts", dataset_axis ,hard_liso_axis),
            'neg_pt' : hist.Hist("Counts", dataset_axis ,soft_lpt_axis),
            'neg_eta' : hist.Hist("Counts", dataset_axis ,soft_leta_axis),
            'neg_phi' : hist.Hist("Counts", dataset_axis ,soft_lphi_axis),
            'neg_pfRelIso04_all' : hist.Hist("Counts", dataset_axis ,soft_liso_axis),
            'zmass': hist.Hist("Counts",dataset_axis,zmass_axis),
            'zpt': hist.Hist("Counts",dataset_axis,zpt_axis),
            'zeta': hist.Hist("Counts",dataset_axis,zeta_axis),
            'zphi': hist.Hist("Counts",dataset_axis,zphi_axis),
            'drmumu':hist.Hist("Counts",dataset_axis,drmumu_axis),
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
        ## Define the CvL, CvB 
        if not hasattr(events,"btagDeepFlavCvL"): 
            events.Jet['btagDeepFlavCvL'] = np.where(((events.Jet.btagDeepFlavC/(1.-events.Jet.btagDeepFlavB))>0)&(events.Jet.pt>15),(events.Jet.btagDeepFlavC/(1.-events.Jet.btagDeepFlavB)),-1)
            events.Jet['btagDeepFlavCvB'] = np.where(((events.Jet.btagDeepFlavC/(events.Jet.btagDeepFlavC+events.Jet.btagDeepFlavB))>0)&(events.Jet.pt>15),(events.Jet.btagDeepFlavC/(events.Jet.btagDeepFlavC+events.Jet.btagDeepFlavB)),-1)
            events.Jet['btagDeepCvL'] = np.where((events.Jet.btagDeepC>0)&(events.Jet.pt>15),(events.Jet.btagDeepC/(1.-events.Jet.btagDeepB)),-1)
            events.Jet['btagDeepCvB'] = np.where((events.Jet.btagDeepC>0)&(events.Jet.pt>15),(events.Jet.btagDeepC/(events.Jet.btagDeepC+events.Jet.btagDeepB)),-1)
        if(isRealData):output['sumw'][dataset] += len(events)
        else:output['sumw'][dataset] += ak.sum(events.genWeight)
        req_lumi=np.ones(len(events), dtype='bool')
        if(isRealData): req_lumi=lumiMasks['2017'](events.run, events.luminosityBlock)
        weights = Weights(len(events), storeIndividual=True)
        if not isRealData:
            weights.add('genweight',events.genWeight)
            weights.add('puweight', compiled['2017_pileupweight'](events.Pileup.nPU))
        ##############
        # Trigger level
        mu_triggers = [
        # "HLT_IsoMu24",
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
        ]
        
        trig_arrs = [events.HLT[_trig.strip("HLT_")] for _trig in mu_triggers]
        req_trig = np.zeros(len(events), dtype='bool')
        for t in trig_arrs:
            req_trig = req_trig | t

        ############
        # Event level
        ## Muon cuts
        dilep_mu = events.Muon[(events.Muon.pt>12)&(abs(events.Muon.eta) < 2.4)& (events.Muon.tightId > .5)&(events.Muon.pfRelIso04_all<=0.15)]
        dilep_ele = events.Electron[(events.Electron.pt>15)&((abs(events.Electron.eta) < 1.4442)|((abs(events.Electron.eta) < 2.5)&(abs(events.Electron.eta) >1.566)))& (events.Electron.mvaFall17V2Iso_WP80 > .5)]
        # req_dilepveto = (ak.count(dilep_mu.pt,axis=1)+ak.count(dilep_ele.pt,axis=1)==2)
        pos_dilep = dilep_mu[dilep_mu.charge>0]
        neg_dilep = dilep_mu[dilep_mu.charge<0]
        req_dilep = (ak.num(pos_dilep.pt)>=1) & (ak.num(neg_dilep.pt)>=1) & (ak.num(dilep_mu.charge)>=2) & (ak.num(dilep_ele.charge)<2)
        pos_dilep = ak.pad_none(pos_dilep,1, axis=1)
        neg_dilep = ak.pad_none(neg_dilep,1, axis=1)
        #dilepton mass
        
        dilep_mass = pos_dilep[:,0]+neg_dilep[:,0]
        req_dilepmass = (dilep_mass.mass>81) & (dilep_mass.mass<101) & (dilep_mass.pt>15)

        ## Jet cuts 
        event_jet =  events.Jet[(events.Jet.pt > 20) & (abs(events.Jet.eta) <= 2.5)&((events.Jet.puId >=7)&( events.Jet.pt<50)) &(events.Jet.jetId>=3)&(events.Jet.btagDeepB>0.) & (events.Jet.btagDeepB<1.) & (events.Jet.btagDeepC>0.) & (events.Jet.btagDeepC<1.) & (events.Jet.btagDeepFlavB>0.) & (events.Jet.btagDeepFlavB<1.) & (events.Jet.btagDeepFlavC>0.) & (events.Jet.btagDeepFlavC<1.)& (ak.all(events.Jet.metric_table(pos_dilep[:,0]) > 0.4, axis=2))&(ak.all(events.Jet.metric_table(neg_dilep[:,0]) > 0.4, axis=2))]

        req_jets = (ak.num(event_jet.puId) >=1)
        event_jet = ak.pad_none(event_jet,1, axis=1)
        event_level= req_lumi & req_trig & req_dilep & req_dilepmass & req_jets
        if(len(event_level)>0):event_level = ak.fill_none(event_level,False)        
        # Selected
        selev = events[event_level]    
        
        #########
        
        ## Hard Muon
        smu = selev.Muon[(selev.Muon.pt>12)&(abs(selev.Muon.eta) < 2.4)& (selev.Muon.tightId > .5)&(selev.Muon.pfRelIso04_all<=0.15)]
        # mu_pairs = smu.distincts()
        # print(mu_pairs.i0.charge,mu_pairs.i1.charge)
        sposmu = selev.Muon[(selev.Muon.pt>12)&(abs(selev.Muon.eta) < 2.4)& (selev.Muon.tightId > .5)&(selev.Muon.pfRelIso04_all<=0.15)&(selev.Muon.charge>0)]
        sposmu = sposmu[:,0]
        snegmu = selev.Muon[(selev.Muon.pt>12)&(abs(selev.Muon.eta) < 2.4)& (selev.Muon.tightId > .5)&(selev.Muon.pfRelIso04_all<=0.15)&(selev.Muon.charge<0)]
        snegmu = snegmu[:,0]
        
        if not isRealData:
            # lepsf = muSFs(sposmu)*muSFs(snegmu)
            weights.add('lep1sf',np.where(event_level,muSFs(ak.firsts(events.Muon[(events.Muon.pt>12)&(abs(events.Muon.eta) < 2.4)& (events.Muon.tightId > .5)&(events.Muon.pfRelIso04_all<=0.15)&(events.Muon.charge<0)])),1.))
            weights.add('lep2sf',np.where(event_level,muSFs(ak.firsts(events.Muon[(events.Muon.pt>12)&(abs(events.Muon.eta) < 2.4)& (events.Muon.tightId > .5)&(events.Muon.pfRelIso04_all<=0.15)&(events.Muon.charge>0)])),1.))
            # weights.add('lep2sf',muSFs(sposmu))
        sz=sposmu+snegmu
        
        ## Jets
        
        sjets = selev.Jet[(selev.Jet.pt > 20) & (abs(selev.Jet.eta) <= 2.5)&((selev.Jet.puId >=7)&( selev.Jet.pt<50)) &(selev.Jet.jetId>=3)&(selev.Jet.btagDeepB>0.) & (selev.Jet.btagDeepB<1.) & (selev.Jet.btagDeepC>0.) & (selev.Jet.btagDeepC<1.) & (selev.Jet.btagDeepFlavB>0.) & (selev.Jet.btagDeepFlavB<1.) & (selev.Jet.btagDeepFlavC>0.) & (selev.Jet.btagDeepFlavC<1.)& (ak.all(selev.Jet.metric_table(sposmu) > 0.4, axis=2))&(selev.Jet.btagDeepFlavC<1.)& (ak.all(selev.Jet.metric_table(snegmu) > 0.4, axis=2))&(selev.Jet.muEF<0.8)]
    
        njet= ak.count(sjets.pt,axis=1)   
        def flatten(ar): # flatten awkward into a 1d array to hist
            return ak.flatten(ar, axis=None)
        if isRealData :
            genflavor = ak.zeros_like(sjets.pt)
        else:
            par_flav = (sjets.partonFlavour == 0 ) & (sjets.hadronFlavour==0)
            genflavor = sjets.hadronFlavour + 1*par_flav 
            # genweiev=ak.flatten(ak.broadcast_arrays(weights.weight()[event_level],sjets['pt'])[0])
            jetsfs_c=jetsfs_c_up=jetsfs_c_dn=jetsfs_b=jetsfs_b_up=jetsfs_b_dn=csvsfs_c=csvsfs_c_up=csvsfs_c_dn=csvsfs_b=csvsfs_b_up=csvsfs_b_dn=1.
            
            ## for each jet
            jetsfs_c_lj = getSF(ak.to_numpy(sjets[:,0].hadronFlavour),ak.to_numpy(sjets[:,0].btagDeepFlavCvL),ak.to_numpy(sjets[:,0].btagDeepFlavCvB),deepjetc_sf)
            jetsfs_c_up_lj = getSF(ak.to_numpy(sjets[:,0].hadronFlavour),ak.to_numpy(sjets[:,0].btagDeepFlavCvL),ak.to_numpy(sjets[:,0].btagDeepFlavCvB),deepjetc_sf,"TotalUncUp")
            jetsfs_c_dn_lj = getSF(ak.to_numpy(sjets[:,0].hadronFlavour),ak.to_numpy(sjets[:,0].btagDeepFlavCvL),ak.to_numpy(sjets[:,0].btagDeepFlavCvB),deepjetc_sf,"TotalUncDown")
            jetsfs_b_lj = deepjetb_sf.eval('central',sjets[:,0].hadronFlavour,abs(sjets[:,0].eta),sjets[:,0].pt,discr=sjets[:,0].btagDeepFlavB)
            jetsfs_b_up_lj = deepjetb_sf.eval('up_jes',sjets[:,0].hadronFlavour,abs(sjets[:,0].eta),sjets[:,0].pt,discr=sjets[:,0].btagDeepFlavB)
            jetsfs_b_dn_lj = deepjetb_sf.eval('down_jes',sjets[:,0].hadronFlavour,abs(sjets[:,0].eta),sjets[:,0].pt,discr=sjets[:,0].btagDeepFlavB)
            csvsfs_c_lj = getSF(ak.to_numpy(sjets[:,0].hadronFlavour),ak.to_numpy(sjets[:,0].btagDeepCvL),ak.to_numpy(sjets[:,0].btagDeepCvB),deepcsvc_sf)
            csvsfs_c_up_lj = getSF(ak.to_numpy(sjets[:,0].hadronFlavour),ak.to_numpy(np.where(sjets[:,0].btagDeepCvL<0,-0.2,sjets[:,0].btagDeepCvL)),ak.to_numpy(sjets[:,0].btagDeepCvB),deepcsvc_sf,"TotalUncUp")
            csvsfs_c_dn_lj = getSF(ak.to_numpy(sjets[:,0].hadronFlavour),ak.to_numpy(np.where(sjets[:,0].btagDeepCvL<0,-0.2,sjets[:,0].btagDeepCvL)),ak.to_numpy(sjets[:,0].btagDeepCvB),deepcsvc_sf,"TotalUncDown")
            csvsfs_b_up_lj = deepcsvb_sf.eval('up_jes',sjets[:,0].hadronFlavour,abs(sjets[:,0].eta),sjets[:,0].pt,discr=sjets[:,0].btagDeepB)
            csvsfs_b_lj = deepcsvb_sf.eval('central',sjets[:,0].hadronFlavour,abs(sjets[:,0].eta),sjets[:,0].pt,discr=sjets[:,0].btagDeepB)
            csvsfs_b_dn_lj = deepcsvb_sf.eval('down_jes',sjets[:,0].hadronFlavour,abs(sjets[:,0].eta),sjets[:,0].pt,discr=sjets[:,0].btagDeepB)
            
            
        for histname, h in output.items():
            if 'sumw' in histname:print('sumw')

            elif histname in self.deepcsv_hists:
                 fields = {l: ak.flatten(sjets[histname]) for l in h.fields if l in dir(sjets)}
            #     if(isRealData):
            #         h.fill(dataset=dataset,flav=5, **fields)
            #     else:                    
            #         genweiev=ak.flatten(ak.broadcast_arrays(weights.weight()[event_level],sjets['pt'])[0])
            #         h.fill(dataset=dataset,flav=ak.flatten(genflavor), **fields,weight=genweiev)
            elif 'pos_' in histname:
                fields = {l: sposmu[l] for l in h.fields if l in dir(sposmu)}
                if(isRealData): h.fill(dataset=dataset,**fields)
                else: h.fill(dataset=dataset,**fields,weight=weights.weight()[event_level])
            elif 'neg_' in snegmu:
                fields = {l: snegmu[l] for l in h.fields if l in dir(snegmu)}
                if(isRealData): h.fill(dataset=dataset,**fields)
                else: h.fill(dataset=dataset,**fields,weight=weights.weight()[event_level])
            
         
        if not isRealData:
            ###Fill no SFs
            output['zmass'].fill(dataset=dataset,zmass=sz.mass,weight=weights.weight()[event_level])
            output['zpt'].fill(dataset=dataset,zpt=sz.pt,weight=weights.weight()[event_level])
            output['zeta'].fill(dataset=dataset,zeta=sz.eta,weight=weights.weight()[event_level])
            output['zphi'].fill(dataset=dataset,zphi=sz.phi,weight=weights.weight()[event_level])
            output['drmumu'].fill(dataset=dataset,drmumu=snegmu.delta_r(sposmu))
           
            ## discri
            output['btagDeepFlavB_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepFlavB=sjets[:,0].btagDeepFlavB,weight=weights.weight()[event_level])
            output['btagDeepFlavC_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepFlavC=sjets[:,0].btagDeepFlavC,weight=weights.weight()[event_level])
            output['btagDeepB_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepB=sjets[:,0].btagDeepB,weight=weights.weight()[event_level])
            output['btagDeepC_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepC=sjets[:,0].btagDeepC,weight=weights.weight()[event_level])
            output['deepcsv_CvB_0'].fill(dataset=dataset,flav=genflavor[:,0], deepcsv_CvB=np.where(sjets[:,0].btagDeepCvB<0,-0.2,sjets[:,0].btagDeepCvB),weight=weights.weight()[event_level])
            output['deepcsv_CvL_0'].fill(dataset=dataset,flav=genflavor[:,0], deepcsv_CvL=np.where(sjets[:,0].btagDeepCvL<0,-0.2,sjets[:,0].btagDeepCvL),weight=weights.weight()[event_level])
            output['deepflav_CvB_0'].fill(dataset=dataset,flav=genflavor[:,0], deepflav_CvB=np.where(sjets[:,0].btagDeepFlavCvB<0,-0.2,sjets[:,0].btagDeepFlavCvB),weight=weights.weight()[event_level])
            output['deepflav_CvL_0'].fill(dataset=dataset,flav=genflavor[:,0], deepflav_CvL=np.where(sjets[:,0].btagDeepFlavCvL<0,-0.2,sjets[:,0].btagDeepFlavCvL),weight=weights.weight()[event_level])
            output['btagDeepFlavBSF_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepFlavBSF=sjets[:,0].btagDeepFlavB,weight=(weights.weight()[event_level]*jetsfs_b_lj))      
            output['btagDeepFlavCSF_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepFlavCSF=sjets[:,0].btagDeepFlavC,weight=weights.weight()[event_level]*jetsfs_c_lj)
            output['btagDeepBSF_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepBSF=sjets[:,0].btagDeepB,weight=weights.weight()[event_level]*csvsfs_b_lj)
            output['btagDeepCSF_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepCSF=sjets[:,0].btagDeepC,weight=weights.weight()[event_level]*csvsfs_c_lj)
            output['deepcsv_CvBSF_0'].fill(dataset=dataset,flav=genflavor[:,0], deepcsv_CvBSF=np.where(sjets[:,0].btagDeepCvB<0,-0.2,sjets[:,0].btagDeepCvB),weight=weights.weight()[event_level]*csvsfs_c_lj)
            output['deepcsv_CvLSF_0'].fill(dataset=dataset,flav=genflavor[:,0], deepcsv_CvLSF=np.where(sjets[:,0].btagDeepCvL<0,-0.2,sjets[:,0].btagDeepCvL),weight=weights.weight()[event_level]*csvsfs_c_lj)
            output['deepflav_CvBSF_0'].fill(dataset=dataset,flav=genflavor[:,0], deepflav_CvBSF=np.where(sjets[:,0].btagDeepFlavCvB<0,-0.2,sjets[:,0].btagDeepFlavCvB),weight=weights.weight()[event_level]*jetsfs_c_lj)
            output['deepflav_CvLSF_0'].fill(dataset=dataset,flav=genflavor[:,0], deepflav_CvLSF=np.where(sjets[:,0].btagDeepFlavCvL<0,-0.2,sjets[:,0].btagDeepFlavCvL),weight=weights.weight()[event_level]*jetsfs_c_lj)
            output['btagDeepFlavB_up_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepFlavB_up=sjets[:,0].btagDeepFlavB,weight=weights.weight()[event_level]*jetsfs_b_up_lj)
            output['btagDeepFlavC_up_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepFlavC_up=sjets[:,0].btagDeepFlavC,weight=weights.weight()[event_level]*jetsfs_c_up_lj)
            output['btagDeepB_up_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepB_up=sjets[:,0].btagDeepB,weight=weights.weight()[event_level]*csvsfs_b_up_lj)
            output['btagDeepC_up_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepC_up=sjets[:,0].btagDeepC,weight=weights.weight()[event_level]*csvsfs_c_up_lj)
            output['deepcsv_CvB_up_0'].fill(dataset=dataset,flav=genflavor[:,0], deepcsv_CvB_up=np.where(sjets[:,0].btagDeepCvB<0,-0.2,sjets[:,0].btagDeepCvB),weight=weights.weight()[event_level]*csvsfs_c_up_lj)
            output['deepcsv_CvL_up_0'].fill(dataset=dataset,flav=genflavor[:,0], deepcsv_CvL_up=np.where(sjets[:,0].btagDeepCvL<0,-0.2,sjets[:,0].btagDeepCvL),weight=weights.weight()[event_level]*csvsfs_c_up_lj)
            output['deepflav_CvB_up_0'].fill(dataset=dataset,flav=genflavor[:,0], deepflav_CvB_up=np.where(sjets[:,0].btagDeepFlavCvB<0,-0.2,sjets[:,0].btagDeepFlavCvB),weight=weights.weight()[event_level]*jetsfs_c_up_lj)
            output['deepflav_CvL_up_0'].fill(dataset=dataset,flav=genflavor[:,0], deepflav_CvL_up=np.where(sjets[:,0].btagDeepFlavCvL<0,-0.2,sjets[:,0].btagDeepFlavCvL),weight=weights.weight()[event_level]*jetsfs_c_up_lj)
            output['btagDeepFlavB_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepFlavB_dn=sjets[:,0].btagDeepFlavB,weight=weights.weight()[event_level]*jetsfs_b_dn_lj)
            output['btagDeepFlavC_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepFlavC_dn=sjets[:,0].btagDeepFlavC,weight=weights.weight()[event_level]*jetsfs_c_dn_lj)
            output['btagDeepB_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepB_dn=sjets[:,0].btagDeepB,weight=weights.weight()[event_level]*csvsfs_b_dn_lj)
            output['btagDeepC_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], btagDeepC_dn=sjets[:,0].btagDeepC,weight=weights.weight()[event_level]*csvsfs_c_dn_lj)
            output['deepcsv_CvB_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], deepcsv_CvB_dn=np.where(sjets[:,0].btagDeepCvB<0,-0.2,sjets[:,0].btagDeepCvB),weight=weights.weight()[event_level]*csvsfs_c_dn_lj)
            output['deepcsv_CvL_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], deepcsv_CvL_dn=np.where(sjets[:,0].btagDeepCvL<0,-0.2,sjets[:,0].btagDeepCvL),weight=weights.weight()[event_level]*csvsfs_c_dn_lj)
            output['deepflav_CvB_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], deepflav_CvB_dn=np.where(sjets[:,0].btagDeepFlavCvB<0,-0.2,sjets[:,0].btagDeepFlavCvB),weight=weights.weight()[event_level]*jetsfs_c_dn_lj)
            output['deepflav_CvL_dn_0'].fill(dataset=dataset,flav=genflavor[:,0], deepflav_CvL_dn=np.where(sjets[:,0].btagDeepFlavCvL<0,-0.2,sjets[:,0].btagDeepFlavCvL),weight=weights.weight()[event_level]*jetsfs_c_dn_lj)
            # print(sjets,ak.count(sjets.pt,axis=1))
            
        else:
        
            output['zmass'].fill(dataset=dataset,zmass=sz.mass)
            output['zpt'].fill(dataset=dataset,zpt=sz.pt)
            output['zeta'].fill(dataset=dataset,zeta=sz.eta)
            output['zphi'].fill(dataset=dataset,zphi=sz.phi)
            output['drmumu'].fill(dataset=dataset,drmumu=snegmu.delta_r(sposmu))                
                ## discri
            output['btagDeepFlavB_0'].fill(dataset=dataset,flav=5, btagDeepFlavB=sjets[:,0].btagDeepFlavB)
            output['btagDeepFlavC_0'].fill(dataset=dataset,flav=5, btagDeepFlavC=sjets[:,0].btagDeepFlavC)
            output['btagDeepB_0'].fill(dataset=dataset,flav=5, btagDeepB=sjets[:,0].btagDeepB)
            output['btagDeepC_0'].fill(dataset=dataset,flav=5, btagDeepC=sjets[:,0].btagDeepC)
            output['deepcsv_CvB_0'].fill(dataset=dataset,flav=5, deepcsv_CvB=np.where(sjets[:,0].btagDeepCvB<0,-0.2,sjets[:,0].btagDeepCvB))
            output['deepcsv_CvL_0'].fill(dataset=dataset,flav=5, deepcsv_CvL=np.where(sjets[:,0].btagDeepCvL<0,-0.2,sjets[:,0].btagDeepCvL))
            output['deepflav_CvB_0'].fill(dataset=dataset,flav=5, deepflav_CvB=np.where(sjets[:,0].btagDeepFlavCvB<0,-0.2,sjets[:,0].btagDeepFlavCvB))
            output['deepflav_CvL_0'].fill(dataset=dataset,flav=5, deepflav_CvL=np.where(sjets[:,0].btagDeepFlavCvL<0,-0.2,sjets[:,0].btagDeepFlavCvL))
            output['btagDeepFlavBSF_0'].fill(dataset=dataset,flav=5, btagDeepFlavBSF=sjets[:,0].btagDeepFlavB)
            output['btagDeepFlavCSF_0'].fill(dataset=dataset,flav=5, btagDeepFlavCSF=sjets[:,0].btagDeepFlavC)
            output['btagDeepBSF_0'].fill(dataset=dataset,flav=5, btagDeepBSF=sjets[:,0].btagDeepB)
            output['btagDeepCSF_0'].fill(dataset=dataset,flav=5, btagDeepCSF=sjets[:,0].btagDeepC)
            output['deepcsv_CvBSF_0'].fill(dataset=dataset,flav=5, deepcsv_CvBSF=np.where(sjets[:,0].btagDeepCvB<0,-0.2,sjets[:,0].btagDeepCvB))
            output['deepcsv_CvLSF_0'].fill(dataset=dataset,flav=5, deepcsv_CvLSF=np.where(sjets[:,0].btagDeepCvL<0,-0.2,sjets[:,0].btagDeepCvL))
            output['deepflav_CvBSF_0'].fill(dataset=dataset,flav=5, deepflav_CvBSF=np.where(sjets[:,0].btagDeepFlavCvB<0,-0.2,sjets[:,0].btagDeepFlavCvB))
            output['deepflav_CvLSF_0'].fill(dataset=dataset,flav=5, deepflav_CvLSF=np.where(sjets[:,0].btagDeepFlavCvL<0,-0.2,sjets[:,0].btagDeepFlavCvL))
            output['btagDeepFlavB_up_0'].fill(dataset=dataset,flav=5, btagDeepFlavB_up=sjets[:,0].btagDeepFlavB)
            output['btagDeepFlavC_up_0'].fill(dataset=dataset,flav=5, btagDeepFlavC_up=sjets[:,0].btagDeepFlavC)
            output['btagDeepB_up_0'].fill(dataset=dataset,flav=5, btagDeepB_up=sjets[:,0].btagDeepB)
            output['btagDeepC_up_0'].fill(dataset=dataset,flav=5, btagDeepC_up=sjets[:,0].btagDeepC)
            output['deepcsv_CvB_up_0'].fill(dataset=dataset,flav=5, deepcsv_CvB_up=np.where(sjets[:,0].btagDeepCvB<0,-0.2,sjets[:,0].btagDeepCvB))
            output['deepcsv_CvL_up_0'].fill(dataset=dataset,flav=5, deepcsv_CvL_up=np.where(sjets[:,0].btagDeepCvL<0,-0.2,sjets[:,0].btagDeepCvL))
            output['deepflav_CvB_up_0'].fill(dataset=dataset,flav=5, deepflav_CvB_up=np.where(sjets[:,0].btagDeepFlavCvB<0,-0.2,sjets[:,0].btagDeepFlavCvB))
            output['deepflav_CvL_up_0'].fill(dataset=dataset,flav=5, deepflav_CvL_up=np.where(sjets[:,0].btagDeepFlavCvL<0,-0.2,sjets[:,0].btagDeepFlavCvL))
            output['btagDeepFlavB_dn_0'].fill(dataset=dataset,flav=5, btagDeepFlavB_dn=sjets[:,0].btagDeepFlavB)
            output['btagDeepFlavC_dn_0'].fill(dataset=dataset,flav=5, btagDeepFlavC_dn=sjets[:,0].btagDeepFlavC)
            output['btagDeepB_dn_0'].fill(dataset=dataset,flav=5, btagDeepB_dn=sjets[:,0].btagDeepB)
            output['btagDeepC_dn_0'].fill(dataset=dataset,flav=5, btagDeepC_dn=sjets[:,0].btagDeepC)
            output['deepcsv_CvB_dn_0'].fill(dataset=dataset,flav=5, deepcsv_CvB_dn=np.where(sjets[:,0].btagDeepCvB<0,-0.2,sjets[:,0].btagDeepCvB))
            output['deepcsv_CvL_dn_0'].fill(dataset=dataset,flav=5, deepcsv_CvL_dn=np.where(sjets[:,0].btagDeepCvL<0,-0.2,sjets[:,0].btagDeepCvL))
            output['deepflav_CvB_dn_0'].fill(dataset=dataset,flav=5, deepflav_CvB_dn=np.where(sjets[:,0].btagDeepFlavCvB<0,-0.2,sjets[:,0].btagDeepFlavCvB))
            output['deepflav_CvL_dn_0'].fill(dataset=dataset,flav=5, deepflav_CvL_dn=np.where(sjets[:,0].btagDeepFlavCvL<0,-0.2,sjets[:,0].btagDeepFlavCvL))
            
        
        
         
        gc.collect()
        # schedule.every(20).minutes.do(dosomething)

        return output

    def postprocess(self, accumulator):
        return accumulator
