import coffea
from coffea import hist, processor
import numpy as np
import awkward1 as ak


class NanoProcessor(processor.ProcessorABC):
    # Define histograms
    def __init__(self):        
        # Define axes
        # Should read axes from NanoAOD config
        dataset_axis = hist.Cat("dataset", "Primary dataset")
        cutflow_axis   = hist.Cat("cut",   "Cut")

        # Events
        nel_axis   = hist.Bin("nel",   r"N electrons", [0,1,2,3,4,5,6,7,8,9,10])
        nmu_axis   = hist.Bin("nmu",   r"N muons",     [0,1,2,3,4,5,6,7,8,9,10])
        njet_axis  = hist.Bin("njet",  r"N jets",      [0,1,2,3,4,5,6,7,8,9,10])
        nbjet_axis = hist.Bin("nbjet", r"N b-jets",    [0,1,2,3,4,5,6,7,8,9,10])

        # Electron
        el_pt_axis   = hist.Bin("pt",    r"Electron $p_{T}$ [GeV]", 100, 20, 400)
        el_eta_axis  = hist.Bin("eta",   r"$\eta$", 60, -3, 3)
        el_phi_axis  = hist.Bin("phi",   r"$\phi$", 60, -3, 3)
        lelpt_axis   = hist.Bin("lelpt", r"Leading electron $p_{T}$ [GeV]", 100, 20, 200)
        
        # Muons
        mu_pt_axis   = hist.Bin("pt",    r"Muon $p_{T}$ [GeV]", 100, 20, 400)
        mu_eta_axis  = hist.Bin("eta",   r"$\eta$", 60, -3, 3)
        mu_phi_axis  = hist.Bin("phi",   r"$\phi$", 60, -3, 3)
        lmupt_axis   = hist.Bin("lmupt", r"Leading muon $p_{T}$ [GeV]", 100, 20, 200)
        
        # Jet
        jet_pt_axis   = hist.Bin("pt",   r"Jet $p_{T}$ [GeV]", 100, 20, 400)
        jet_eta_axis  = hist.Bin("eta",  r"$\eta$", 60, -3, 3)
        jet_phi_axis  = hist.Bin("phi",  r"$\phi$", 60, -3, 3)
        jet_mass_axis = hist.Bin("mass", r"Jet $m$ [GeV]", 100, 0, 50)
        ljpt_axis     = hist.Bin("ljpt", r"Leading jet $p_{T}$ [GeV]", 100, 20, 400)
        sljpt_axis     = hist.Bin("sljpt", r"Subleading jet $p_{T}$ [GeV]", 100, 20, 400)
        
        # Define similar axes dynamically
        disc_list = ["btagCMVA", "btagCSVV2", 'btagDeepB', 'btagDeepC', 'btagDeepFlavB', 'btagDeepFlavC',]
        btag_axes = []
        for d in disc_list:
            btag_axes.append(hist.Bin(d, d, 50, 0, 1))        
        
        deepcsv_list = ["DeepCSV_trackDecayLenVal_0", "DeepCSV_trackDecayLenVal_1", "DeepCSV_trackDecayLenVal_2", "DeepCSV_trackDecayLenVal_3", "DeepCSV_trackDecayLenVal_4", "DeepCSV_trackDecayLenVal_5", 
	"DeepCSV_trackDeltaR_0", "DeepCSV_trackDeltaR_1", "DeepCSV_trackDeltaR_2", "DeepCSV_trackDeltaR_3", "DeepCSV_trackDeltaR_4", "DeepCSV_trackDeltaR_5",
	"DeepCSV_trackEtaRel_0","DeepCSV_trackEtaRel_1","DeepCSV_trackEtaRel_2","DeepCSV_trackEtaRel_3", 	
	"DeepCSV_trackJetDistVal_0","DeepCSV_trackJetDistVal_1","DeepCSV_trackJetDistVal_2","DeepCSV_trackJetDistVal_3","DeepCSV_trackJetDistVal_4","DeepCSV_trackJetDistVal_5", 
	"DeepCSV_trackJetPt", "DeepCSV_trackPtRatio_0","DeepCSV_trackPtRatio_1","DeepCSV_trackPtRatio_2","DeepCSV_trackPtRatio_3","DeepCSV_trackPtRatio_4","DeepCSV_trackPtRatio_5", 
	"DeepCSV_trackPtRel_0", "DeepCSV_trackPtRel_1","DeepCSV_trackPtRel_2","DeepCSV_trackPtRel_3","DeepCSV_trackPtRel_4","DeepCSV_trackPtRel_5",
	"DeepCSV_vertexCategory","DeepCSV_vertexEnergyRatio", "DeepCSV_vertexJetDeltaR","DeepCSV_vertexMass", "DeepCSV_jetNSecondaryVertices","DeepCSV_jetNSelectedTracks","DeepCSV_jetNTracksEtaRel","DeepCSV_vertexNTracks"]
        deepcsv_axes = []
        for d in deepcsv_list:
            if "DeltaR" in d: 
                deepcsv_axes.append(hist.Bin(d, d, 50, 0, 0.3))
            elif "trackPtRatio" in d:
                deepcsv_axes.append(hist.Bin(d, d, 50, 0, 0.3))
            elif "vertexMass" in d:
                deepcsv_axes.append(hist.Bin(d, d, 50, 0, 5))
            elif "trackEtaRel" in d:
                deepcsv_axes.append(hist.Bin(d, d, 50, 0, 6))
            elif "EnergyRatio" in d:
                deepcsv_axes.append(hist.Bin(d, d, 50, 0, 1))
            elif "trackJetDistVal" in d:
                deepcsv_axes.append(hist.Bin(d, d, 50, -0.1, 0))
            elif "vertexCategory" in d:
                deepcsv_axes.append(hist.Bin(d, d, 3, 0, 2))
            elif "jetN" in d:
                deepcsv_axes.append(hist.Bin(d, d, 16, 0, 15))
            elif "vertexN" in d:
                deepcsv_axes.append(hist.Bin(d, d, 16, 0, 15))
            else:
                deepcsv_axes.append(hist.Bin(d, d, 50, 0, 3))

        # Define histograms from axes
        _hist_jet_dict = {
                'pt'  : hist.Hist("Counts", dataset_axis, jet_pt_axis),
                'eta' : hist.Hist("Counts", dataset_axis, jet_eta_axis),
                'phi' : hist.Hist("Counts", dataset_axis, jet_phi_axis),
                'mass': hist.Hist("Counts", dataset_axis, jet_mass_axis),
            }
        _hist_deepcsv_dict = {
                'pt'  : hist.Hist("Counts", dataset_axis, jet_pt_axis),
                'eta' : hist.Hist("Counts", dataset_axis, jet_eta_axis),
                'phi' : hist.Hist("Counts", dataset_axis, jet_phi_axis),
                'mass': hist.Hist("Counts", dataset_axis, jet_mass_axis),
            }
 
        # Generate some histograms dynamically
        for disc, axis in zip(disc_list, btag_axes):
            _hist_jet_dict[disc] = hist.Hist("Counts", dataset_axis, axis)
        for deepcsv, axis in zip(deepcsv_list, deepcsv_axes):
            _hist_deepcsv_dict[deepcsv] = hist.Hist("Counts", dataset_axis, axis)
        
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
        self.deepcsv_hists = list(_hist_deepcsv_dict.keys())
        self.event_hists = list(_hist_event_dict.keys())
    
        _hist_dict = {**_hist_jet_dict, **_hist_deepcsv_dict, **_hist_event_dict}
        self._accumulator = processor.dict_accumulator(_hist_dict)


    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.accumulator.identity()

        dataset = events.metadata['dataset']
        
        ##############
        # Trigger level
        triggers = [
        #"HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
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
        
        # Per jet
        jet_eta    = (abs(selev.Jet.eta) <= 2.4)
        jet_pt     = selev.Jet.pt > 25
        jet_pu     = selev.Jet.puId > 6
        jet_level  = jet_pu & jet_eta & jet_pt
        
        # b-tag twiki : https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
        bjet_disc  = selev.Jet.btagDeepB > 0.7264 # L=0.0494, M=0.2770, T=0.7264
        bjet_level = jet_level & bjet_disc
        

        sel    = selev.Electron[el_level]
        smu    = selev.Muon[mu_level]
        sjets  = selev.Jet[jet_level]
        sbjets = selev.Jet[bjet_level]
        
        # output['pt'].fill(dataset=dataset, pt=selev.Jet.pt.flatten())
        # Fill histograms dynamically  
        for histname, h in output.items():
            if (histname not in self.jet_hists) and (histname not in self.deepcsv_hists): continue
            # Get valid fields perhistogram to fill
            fields = {k: ak.flatten(sjets[k], axis=None) for k in h.fields if k in dir(sjets)}
            h.fill(dataset=dataset, **fields)


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
