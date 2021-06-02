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
        nel_axis   = hist.Bin("nel",   r"N electrons", [0,1,2,3,4,5,6,7,8,9,10])
        nmu_axis   = hist.Bin("nmu",   r"N muons",     [0,1,2,3,4,5,6,7,8,9,10])
        njet_axis  = hist.Bin("njet",  r"N jets",      [0,1,2,3,4,5,6,7,8,9,10])
        nbjet_axis = hist.Bin("nbjet", r"N b-jets",    [0,1,2,3,4,5,6,7,8,9,10])

        # Electron
        el_pt_axis   = hist.Bin("el_pt",    r"Electron $p_{T}$ [GeV]", 100, 20, 400)
        el_eta_axis  = hist.Bin("el_eta",   r"Electron $\eta$", 60, -3, 3)
        el_phi_axis  = hist.Bin("el_phi",   r"Electron $\phi$", 60, -3, 3)
        lelpt_axis   = hist.Bin("lelpt", r"Leading electron $p_{T}$ [GeV]", 100, 20, 200)
        
        # Muons
        mu_pt_axis   = hist.Bin("mu_pt",    r"Muon $p_{T}$ [GeV]", 100, 20, 400)
        mu_eta_axis  = hist.Bin("mu_eta",   r"Muon $\eta$", 60, -3, 3)
        mu_phi_axis  = hist.Bin("mu_phi",   r"Muon $\phi$", 60, -3, 3)
        lmupt_axis   = hist.Bin("lmupt", r"Leading muon $p_{T}$ [GeV]", 100, 20, 200)
            

        # Jet
        jet_pt_axis   = hist.Bin("jet_pt",   r"Jet $p_{T}$ [GeV]", 100, 20, 400)
        jet_eta_axis  = hist.Bin("jet_eta",  r"Jet $\eta$", 60, -3, 3)
        jet_phi_axis  = hist.Bin("jet_phi",  r"Jet $\phi$", 60, -3, 3)
        jet_mass_axis = hist.Bin("jet_mass", r"Jet $m$ [GeV]", 100, 0, 50)
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
        _hist_fatjet_dict = {
                'fatjet_pt'  : hist.Hist("Counts", dataset_axis, flav_axis,fatjet_pt_axis),
                'fatjet_eta' : hist.Hist("Counts", dataset_axis, flav_axis,fatjet_eta_axis),
                'fatjet_phi' : hist.Hist("Counts", dataset_axis, flav_axis,fatjet_phi_axis),
                'fatjet_mass': hist.Hist("Counts", dataset_axis, flav_axis,fatjet_mass_axis)
            
            }
        _hist_subjet_dict = {
                'subjet_pt'  : hist.Hist("Counts", dataset_axis, flav_axis,subjet_pt_axis),
                'subjet_eta' : hist.Hist("Counts", dataset_axis, flav_axis,subjet_eta_axis),
                'subjet_phi' : hist.Hist("Counts", dataset_axis, flav_axis,subjet_phi_axis),
                'subjet_mass': hist.Hist("Counts", dataset_axis, flav_axis,subjet_mass_axis)
                
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
            _hist_deepcsv_dict[disc] = hist.Hist("Counts", dataset_axis, flav_axis, axis)
        for disc, axis in zip(ddx_list, btag_axes):
            _hist_deepddx_dict[disc] = hist.Hist("Counts", dataset_axis, flav_axis, axis)
        for deepcsv, axises in zip(deepcsv_list, deepcsv_axes):
            _hist_deepcsv_dict[deepcsv] = hist.Hist("Counts", dataset_axis,flav_axis,  axises)
        for deepddx, axiss in zip(deepddx_list, deepddx_axes):
            _hist_deepddx_dict[deepddx] = hist.Hist("Counts", dataset_axis,flav_axis, axiss)
        _hist_event_dict = {
                'njet'  : hist.Hist("Counts", dataset_axis, njet_axis),
                'nbjet' : hist.Hist("Counts", dataset_axis, nbjet_axis),
                'nel'   : hist.Hist("Counts", dataset_axis, nel_axis),
                'nmu'   : hist.Hist("Counts", dataset_axis, nmu_axis),
                'lelpt' : hist.Hist("Counts", dataset_axis, lelpt_axis),
                'lmupt' : hist.Hist("Counts", dataset_axis, lmupt_axis),
                'ljpt'  : hist.Hist("Counts", dataset_axis, ljpt_axis),
                'sljpt'  : hist.Hist("Counts", dataset_axis, sljpt_axis),
            }
        self.jet_hists = list(_hist_jet_dict.keys())
        self.fatjet_hists = list(_hist_fatjet_dict.keys())
        self.subjet_hists = list(_hist_subjet_dict.keys())
        self.deepcsv_hists = list(_hist_deepcsv_dict.keys())
        self.deepddx_hists = list(_hist_deepddx_dict.keys())
        self.event_hists = list(_hist_event_dict.keys())
        _hist_dict = {**_hist_jet_dict,**_hist_deepcsv_dict,**_hist_deepddx_dict, **_hist_fatjet_dict, **_hist_subjet_dict,**_hist_event_dict}
        self._accumulator = processor.dict_accumulator(_hist_dict)
        self._accumulator['sumw'] = processor.defaultdict_accumulator(float)


    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        print(events)
        output = self.accumulator.identity()
        # print(output.items())
        dataset = events.metadata['dataset']
        isRealData = not hasattr(events, "genWeight")
        if(isRealData):output['sumw'][dataset] += ak.sum(1.)
        else:output['sumw'][dataset] += ak.sum(events.genWeight)
        
        ##############
        # Trigger level
        triggers = [
        "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
        "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",    
        ]
        
        trig_arrs = [events.HLT[_trig.strip("HLT_")] for _trig in triggers]
        req_trig = np.zeros(len(events), dtype='bool')
        for t in trig_arrs:
            req_trig = req_trig | t

        ############
        # Event level
        
        ## Muon cuts
        # muon twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
        events.Muon = events.Muon[(events.Muon.pt > 30) & (abs(events.Muon.eta < 2.4))] # & (events.Muon.tightId > .5)
        events.Muon = ak.pad_none(events.Muon, 1, axis=1) 
        req_muon =(ak.count(events.Muon.pt, axis=1) == 1)
        
        ## Electron cuts
        # electron twiki: https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
        events.Electron = events.Electron[(events.Electron.pt > 30) & (abs(events.Electron.eta) < 2.4)]
        events.Electron = ak.pad_none(events.Electron, 1, axis=1) 
        req_ele = (ak.count(events.Electron.pt, axis=1) == 1)
        
        ## Jet cuts
        events.Jet = events.Jet[(events.Jet.pt > 25) & (abs(events.Jet.eta) <= 2.5)]
        req_jets = (ak.count(events.Jet.pt, axis=1) >= 2)    
        
        
        req_opposite_charge = events.Electron[:, 0].charge * events.Muon[:, 0].charge == -1
        
        event_level = req_trig & req_muon & req_ele & req_opposite_charge & req_jets
        
        # Selected
        selev = events[event_level]    
        
        #########
        
        # Per electron
        el_eta   = (abs(selev.Electron.eta) <= 2.4)
        el_pt    = selev.Electron.pt > 30
        el_level = el_eta & el_pt
        
        # Per muon
        mu_eta   = (abs(selev.Muon.eta) <= 2.4)
        mu_pt    = selev.Muon.pt > 30
        mu_level = mu_eta & mu_pt
        

        # Per jet : https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJetID
        jet_eta    = (abs(selev.Jet.eta) <= 2.4)
        jet_pt     = selev.Jet.pt > 25 
        jet_pu     = selev.Jet.puId > 6
        jet_level  = jet_pu & jet_eta & jet_pt
        
        # b-tag twiki : https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
        bjet_disc  = selev.Jet.btagDeepB > 0.7264 # L=0.0494, M=0.2770, T=0.7264
        bjet_level = jet_level & bjet_disc
        

        fatjet_eta = (abs(selev.FatJet.eta)<=2.5)
        fatjet_pt  = selev.FatJet.pt>350.
        fatjet_ID  = selev.FatJet.jetId!=0
        fatjet_prunemass = selev.FatJet.mass<1e6
        fatjet_tau21=selev.FatJet.tau2/selev.FatJet.tau1
        fatjet_tau21cut = (fatjet_tau21<=0.6)&(fatjet_tau21>=0.0)
        
        fatjet_level = fatjet_eta&fatjet_pt&fatjet_ID&fatjet_prunemass
        subjet_pt = selev.SubJet.pt>0.
        subjet_level = subjet_pt
       

        sel    = selev.Electron[el_level]
        smu    = selev.Muon[mu_level]
        sjets  = selev.Jet[jet_level]
        sjets_b = selev.Jet[(jet_level)&(abs(selev.Jet.hadronFlavour)==5)]
        sjets_c = selev.Jet[(jet_level)&(abs(selev.Jet.hadronFlavour)==4)]
        sjets_l = selev.Jet[(jet_level)&(abs(selev.Jet.hadronFlavour)<4)]
        sbjets = selev.Jet[bjet_level]
        sfatjets = selev.FatJet[fatjet_level]
        sfatjets_b = selev.FatJet[(fatjet_level)&(abs(selev.FatJet.hadronFlavour)==5)]
        sfatjets_c = selev.FatJet[(fatjet_level)&(abs(selev.FatJet.hadronFlavour)==4)]
        sfatjets_l = selev.FatJet[(fatjet_level)&(abs(selev.FatJet.hadronFlavour)<4)]
        ssubjets = selev.SubJet[subjet_level]
        ssubjets_b = selev.SubJet[(subjet_level)&(abs(selev.SubJet.hadronFlavour)==5)]
        ssubjets_c = selev.SubJet[(subjet_level)&(abs(selev.SubJet.hadronFlavour)==4)]
        ssubjets_l = selev.SubJet[(subjet_level)&(abs(selev.SubJet.hadronFlavour)<4)]
        
        # output['pt'].fill(dataset=dataset, pt=selev.Jet.pt.flatten())
        # Fill histograms dynamically  
        for histname, h in output.items():
            if (histname in self.jet_hists) or (histname in self.deepcsv_hists):
                fields = {l: ak.flatten(sjets_b[l.replace('jet_','')], axis=None) for l in h.fields if l.replace('jet_','') in dir(sjets_b)}
                h.fill(dataset=dataset,flav=5, **fields)
                fields2 = {l: ak.flatten(sjets_c[l.replace('jet_','')], axis=None) for l in h.fields if l.replace('jet_','') in dir(sjets_c)}
                h.fill(dataset=dataset,flav=4, **fields2)
                fields3 = {l: ak.flatten(sjets_l[l.replace('jet_','')], axis=None) for l in h.fields if l.replace('jet_','') in dir(sjets_l)}
                h.fill(dataset=dataset,flav=0, **fields3)
            elif (histname in self.deepddx_hists): 
                fields_ddx = {l: ak.flatten(sfatjets_b[l.replace('jet_','')], axis=None) for l in h.fields if l.replace('jet_','') in dir(sfatjets_b)}
                h.fill(dataset=dataset,flav=5, **fields_ddx)
                fields_ddx2 = {l: ak.flatten(sfatjets_c[l.replace('jet_','')], axis=None) for l in h.fields if l.replace('jet_','') in dir(sfatjets_c)}
                h.fill(dataset=dataset,flav=4, **fields_ddx2)
                fields_ddx3 = {l: ak.flatten(sfatjets_l[l.replace('jet_','')], axis=None) for l in h.fields if l.replace('jet_','') in dir(sfatjets_l)}
                h.fill(dataset=dataset,flav=0, **fields_ddx3)
            elif (histname in self.subjet_hists):
                for l in h.fields:
                    if l.replace('subjet_','') in dir(ssubjets_b):
                        fields_s = {l: ak.flatten(ssubjets_b[l.replace('subjet_','')], axis=None)}
                        h.fill(dataset=dataset, flav=5, **fields_s)
                    if l.replace('subjet_','') in dir(ssubjets_c):
                        fields_s1 = {l: ak.flatten(ssubjets_c[l.replace('subjet_','')], axis=None)}
                        h.fill(dataset=dataset, flav=4, **fields_s1)
                    if l.replace('subjet_','') in dir(ssubjets_l):
                        fields_s2 = {l: ak.flatten(ssubjets_l[l.replace('subjet_','')], axis=None)}
                        h.fill(dataset=dataset, flav=0, **fields_s2)
            elif (histname in self.fatjet_hists):
                for l in h.fields:
                    if l.replace('fatjet_','') in dir(sfatjets_b):
                        fields_f = {l: ak.flatten(sfatjets_b[l.replace('fatjet_','')], axis=None)}
                        h.fill(dataset=dataset,flav=5, **fields_f)
                    if l.replace('fatjet_','') in dir(sfatjets_c):
                        fields_f1 = {l: ak.flatten(sfatjets_c[l.replace('fatjet_','')], axis=None)}
                        h.fill(dataset=dataset,flav=4, **fields_f1)
                    if l.replace('fatjet_','') in dir(sfatjets_l):
                        fields_f2 = {l: ak.flatten(sfatjets_l[l.replace('fatjet_','')], axis=None)}
                        h.fill(dataset=dataset,flav=0, **fields_f2)


        def flatten(ar): # flatten awkward into a 1d array to hist
            return ak.flatten(ar, axis=None)

        def num(ar):
            return ak.num(ak.fill_none(ar[~ak.is_none(ar)], 0), axis=0)

        output['njet'].fill(dataset=dataset,  njet=flatten(ak.num(sjets)))
        output['nbjet'].fill(dataset=dataset, nbjet=flatten(ak.num(sbjets)))
        output['nel'].fill(dataset=dataset,   nel=flatten(ak.num(sel)))
        output['nmu'].fill(dataset=dataset,   nmu=flatten(ak.num(smu)))

        output['lelpt'].fill(dataset=dataset, lelpt=flatten(selev.Electron[:, 0].pt))
        output['lmupt'].fill(dataset=dataset, lmupt=flatten(selev.Muon[:, 0].pt))
        output['ljpt'].fill(dataset=dataset,  ljpt=flatten(selev.Jet[:, 0].pt))
        output['sljpt'].fill(dataset=dataset,  sljpt=flatten(selev.Jet[:, 1].pt))

        return output

    def postprocess(self, accumulator):
        return accumulator
