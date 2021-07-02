import gzip
import pickle
import coffea
from coffea import hist, processor
import numpy as np
import awkward as ak
from coffea.analysis_tools import Weights
from coffea.lumi_tools import LumiMask
from coffea.btag_tools import BTagScaleFactor
deepcsv_sf = BTagScaleFactor("data/DeepCSV_94XSF_V5_B_F.csv",BTagScaleFactor.RESHAPE,methods='iterativefit,iterativefit,iterativefit')
deepjet_sf = BTagScaleFactor("data/DeepFlavour_94XSF_V4_B_F.csv",BTagScaleFactor.RESHAPE,methods='iterativefit,iterativefit,iterativefit')

with gzip.open("corrections.pkl.gz") as fin:
    compiled = pickle.load(fin)
def build_lumimask(filename):
    return LumiMask(filename)
lumiMasks = {
    '2016': build_lumimask('Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'),
    '2017': build_lumimask('Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt'),
    '2018': build_lumimask('Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'),
}


class NanoProcessor(processor.ProcessorABC):
    # Define histograms
    def __init__(self):        
        # Define axes
        # Should read axes from NanoAOD config
        dataset_axis = hist.Cat("dataset", "Primary dataset")
        flav_axis = hist.Bin("flav", r"Genflavour",[0,1,4,5,6])
        cutflow_axis   = hist.Cat("cut",   "Cut")


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
                'pt'  : hist.Hist("Counts", dataset_axis, flav_axis,jet_pt_axis),
                'eta' : hist.Hist("Counts", dataset_axis, flav_axis,jet_eta_axis),
                'phi' : hist.Hist("Counts", dataset_axis, flav_axis,jet_phi_axis),
                'mass': hist.Hist("Counts", dataset_axis, flav_axis,jet_mass_axis)
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
        
        for disc, axis in zip(disc_list, btag_axes):
            _hist_deepcsv_dict[disc] = hist.Hist("Counts", dataset_axis,  flav_axis, axis)
        for deepcsv, axises in zip(deepcsv_list, deepcsv_axes):
            _hist_deepcsv_dict[deepcsv] = hist.Hist("Counts", dataset_axis,flav_axis,  axises)
        _hist_event_dict = {
                'njet'  : hist.Hist("Counts", dataset_axis, njet_axis),
                'nbjet' : hist.Hist("Counts", dataset_axis, nbjet_axis),
            'ljpt'  : hist.Hist("Counts", dataset_axis, flav_axis,ljpt_axis),
            'sljpt'  : hist.Hist("Counts", dataset_axis, flav_axis,sljpt_axis),
            }
        self.jet_hists = list(_hist_jet_dict.keys())
        self.deepcsv_hists = list(_hist_deepcsv_dict.keys())
#        self.deepddx_hists = list(_hist_deepddx_dict.keys())
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
        req_lumi=np.ones(len(events), dtype='bool')
        if(isRealData): req_lumi=lumiMasks['2017'](events.run, events.luminosityBlock)
        ##############
        # Trigger level
        triggers = [
            "HLT_IsoMu24"
        #"HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
        #"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
        # "HLT_PFJet40" 
        ]
        
        trig_arrs = [events.HLT[_trig.strip("HLT_")] for _trig in triggers]
        req_trig = np.zeros(len(events), dtype='bool')
        for t in trig_arrs:
            req_trig = req_trig | t

        ############
        ## Muon cuts
        # muon twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
        events.Muon = events.Muon[(events.Muon.pt > 30) & (abs(events.Muon.eta < 2.4))& (events.Muon.tightId)]  
        events.Muon = ak.pad_none(events.Muon, 1, axis=1) 
        req_muon =(ak.count(events.Muon.pt, axis=1) == 1)
        #event_muon = ak.pad_none(events.Muon, 1, axis=1)
        ## Electron cuts
        # electron twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
        events.Electron = events.Electron[(events.Electron.pt > 30) & (abs(events.Electron.eta) < 2.4)&(events.Electron.cutBased>3)]
        #event_electron = ak.pad_none(events.Electron, 1, axis=1)
        events.Electron = ak.pad_none(events.Electron, 1, axis=1) 
        req_ele = (ak.count(events.Electron.pt, axis=1) == 1)

        ## Jet cuts
        events.Jet = events.Jet[(events.Jet.pt > 25) & (abs(events.Jet.eta) <= 2.5)&(events.Jet.puId > 5) & (events.Jet.jetId>0)]
        #&(ak.all(events.Jet.metric_table(events.Muon) > 0.4, axis=2))&(ak.all(events.Jet.metric_table(events.Electron) > 0.4, axis=2))]
        req_jets = (ak.count(events.Jet.pt, axis=1) >= 2)    
        
        
        req_opposite_charge = events.Electron[:, 0].charge * events.Muon[:, 0].charge == -1
        event_dilep = events.Electron[:, 0]+events.Muon[:, 0]
        req_dilep_mass = event_dilep.mass >= 12        
        event_level = req_trig & req_muon & req_ele & req_opposite_charge & req_jets & req_lumi #& req_dilep_mass
        
        # Selected
        selev = events[event_level]    
        
        #########
        smu = selev.Muon[(selev.Muon.pt > 30) & (abs(selev.Muon.eta < 2.4))& (selev.Muon.tightId)]
        sele = selev.Electron[(selev.Electron.pt > 30) & (abs(selev.Electron.eta) < 2.4)&(selev.Electron.cutBased>3)]
        
        # Per jet : https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID
        jet_eta    = (abs(selev.Jet.eta) <= 2.4)
        jet_pt     = selev.Jet.pt > 25
        jet_pu     = (selev.Jet.puId > 5)&(selev.Jet.jetId>0)
        jet_dr = (ak.all(selev.Jet.metric_table(smu) > 0.4, axis=2)) & (ak.all(selev.Jet.metric_table(sele) > 0.4, axis=2)) 
        jet_level  = jet_pu & jet_eta & jet_pt #& jet_dr
        
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
        
        else:
            par_flav = (sjets.partonFlavour == 0 ) & (sjets.hadronFlavour==0)
            genflavor = sjets.hadronFlavour + 1*par_flav 
            weights.add('genweight',selev.genWeight)
            weights.add('puweight', compiled['2017_pileupweight'](selev.Pileup.nPU))
            
        # Fill histograms dynamically  
        for histname, h in output.items():
            if histname in self.deepcsv_hists:
                if(isRealData):
                    if histname == 'rawpt':
                        h.fill(dataset=dataset,flav=ak.flatten(genflavor), rawpt =ak.flatten(sjets.pt))
                    else:
                        fields = {l: ak.flatten(sjets[l.replace('jet_','')], axis=None) for l in h.fields if l.replace('jet_','') in dir(sjets)}
                        h.fill(dataset=dataset,flav=5, **fields)
                else:
                    if histname == 'rawpt':
                        genweiev=ak.flatten(ak.broadcast_arrays(weights.weight(),sjets['pt'])[0])
                        h.fill(dataset=dataset,flav=ak.flatten(genflavor), rawpt =ak.flatten(sjets.pt*(1-sjets.rawFactor)) ,weight=genweiev)
                    else:
                        fields = {l: ak.flatten(sjets[histname]) for l in h.fields if l in dir(sjets)}
                        genweiev=ak.flatten(ak.broadcast_arrays(weights.weight(),sjets['pt'])[0])
                        h.fill(dataset=dataset,flav=ak.flatten(genflavor), **fields,weight=genweiev)
            if histname in self.deepcsvSF_hists:
                fields = {l: ak.flatten(sjets[l.replace('SF','')], axis=None) for l in h.fields if l.replace('SF','') in dir(sjets)}
                if(isRealData):h.fill(dataset=dataset,flav=5, **fields)
        if not isRealData:
            if 'Flav' in histname:
                jetsfs_b=ak.fill_none(ak.prod(deepjet_sf.eval('central',sjets.hadronFlavour,abs(sjets.eta),sjets.pt,discr=sjets.btagDeepFlavB),axis=-1),1.)
                flavsf_b=ak.flatten(ak.broadcast_arrays(weights.weight()*jetsfs_b,sjets['pt'])[0])
                output['btagDeepFlavBSF'].fill(dataset=dataset,flav=ak.flatten(genflavor), btagDeepFlavBSF=ak.flatten(sjets.btagDeepFlavB),weight=flavsf_b)
                jetsfs_c=ak.fill_none(ak.prod(deepjet_sf.eval('central',sjets.hadronFlavour,abs(sjets.eta),sjets.pt,discr=sjets.btagDeepFlavC),axis=-1),1.)
                flavsf_c=ak.flatten(ak.broadcast_arrays(weights.weight()*jetsfs_c,sjets['pt'])[0])
                output['btagDeepFlavCSF'].fill(dataset=dataset,flav=ak.flatten(genflavor), btagDeepFlavCSF=ak.flatten(sjets.btagDeepFlavC),weight=flavsf_c)
            else:
                csvsfs_b=ak.fill_none(ak.prod(deepcsv_sf.eval('central',sjets.hadronFlavour,abs(sjets.eta),sjets.pt,discr=sjets.btagDeepB),axis=-1),1.)
                deepcsvsf_b=ak.flatten(ak.broadcast_arrays(weights.weight()*csvsfs_b,sjets['pt'])[0])
                output['btagDeepBSF'].fill(dataset=dataset,flav=ak.flatten(genflavor), btagDeepBSF=ak.flatten(sjets.btagDeepB),weight=deepcsvsf_b)
                csvsfs_c=ak.fill_none(ak.prod(deepcsv_sf.eval('central',sjets.hadronFlavour,abs(sjets.eta),sjets.pt,discr=sjets.btagDeepC),axis=-1),1.)
                deepcsvsf_c=ak.flatten(ak.broadcast_arrays(weights.weight()*csvsfs_c,sjets['pt'])[0])
                output['btagDeepCSF'].fill(dataset=dataset,flav=ak.flatten(genflavor), btagDeepCSF=ak.flatten(sjets.btagDeepC),weight=deepcsvsf_c)



        def flatten(ar): # flatten awkward into a 1d array to hist
            return ak.flatten(ar, axis=None)

        def num(ar):
            return ak.num(ak.fill_none(ar[~ak.is_none(ar)], 0), axis=0)
        if(isRealData):
            output['njet'].fill(dataset=dataset,  njet=flatten(ak.num(sjets)))
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
            
            output['njet'].fill(dataset=dataset,  njet=flatten(ak.num(sjets)),weight=weights.weight())
            output['nmu'].fill(dataset=dataset,   nmu=flatten(ak.num(smu)),weight=weights.weight())
            output['lmupt'].fill(dataset=dataset, lmupt=flatten((smu.pt)),weight=weights.weight())
            output['ljpt'].fill(dataset=dataset, flav=flatten(sjets[:,0].hadronFlavour), ljpt=flatten(sjets[:,0].pt),weight=weights.weight())
            output['sljpt'].fill(dataset=dataset, flav=flatten(sjets[:,1].hadronFlavour), sljpt=flatten(sjets[:,1].pt),weight=weights.weight())
            output['ssljpt'].fill(dataset=dataset, flav=flatten(sjets[:,2].hadronFlavour), ssljpt=flatten(sjets[:,2].pt),weight=weights.weight())
            output['sssljpt'].fill(dataset=dataset, flav=flatten(sjets[:,3].hadronFlavour), sssljpt=flatten(sjets[:,3].pt),weight=weights.weight())
            output['ljrawpt'].fill(dataset=dataset, flav=flatten(sjets[:,0].hadronFlavour), ljrawpt=flatten(sjets[:,0].pt*(1-sjets[:,0].rawFactor)),weight=weights.weight())
            output['sljrawpt'].fill(dataset=dataset, flav=flatten(sjets[:,1].hadronFlavour), sljrawpt=flatten(sjets[:,1].pt*(1-sjets[:,1].rawFactor)),weight=weights.weight())
            output['ssljrawpt'].fill(dataset=dataset, flav=flatten(sjets[:,2].hadronFlavour), ssljrawpt=flatten(sjets[:,2].pt*(1-sjets[:,2].rawFactor)),weight=weights.weight())
            output['sssljrawpt'].fill(dataset=dataset, flav=flatten(sjets[:,3].hadronFlavour), sssljrawpt=flatten(sjets[:,3].pt*(1-sjets[:,3].rawFactor)),weight=weights.weight())
            output['ljdr'].fill(dataset=dataset, flav=flatten(sjets[:,0].hadronFlavour), ljdr=flatten(sjets[:,0].delta_r(smu)),weight=weights.weight())
            output['sljdr'].fill(dataset=dataset, flav=flatten(sjets[:,1].hadronFlavour), sljdr=flatten(sjets[:,1].delta_r(smu)),weight=weights.weight())
            output['ssljdr'].fill(dataset=dataset, flav=flatten(sjets[:,2].hadronFlavour), ssljdr=flatten(sjets[:,2].delta_r(smu)),weight=weights.weight())
            output['sssljdr'].fill(dataset=dataset, flav=flatten(sjets[:,3].hadronFlavour), sssljdr=flatten(sjets[:,3].delta_r(smu)),weight=weights.weight())
            output['met'].fill(dataset=dataset, met=flatten((selev.METFixEE2017.pt)),weight=weights.weight())
        return output

    def postprocess(self, accumulator):
        return accumulator
