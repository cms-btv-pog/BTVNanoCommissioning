import gzip
import pickle
import coffea
from coffea import hist, processor
import numpy as np
import awkward as ak
from coffea.analysis_tools import Weights,PackedSelection
from coffea.lumi_tools import LumiMask
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
        flav_axis = hist.Bin("flav", r"Genflavour",[0,4,5,6])
        cutflow_axis   = hist.Cat("region",   r"Cut", ['musel','jetsel','metsel','trigger'])

        # Events
        nmu_axis   = hist.Bin("nmu",   r"N muons",     [0,1,2,3,4,5,6,7,8,9,10])
        njet_axis  = hist.Bin("njet",  r"N jets",      [0,1,2,3,4,5,6,7,8,9,10])
        lmupt_axis   = hist.Bin("lmupt",   r"Muon pt", 45, 20, 200)
        met_axis = hist.Bin("met",   r"Missing ET", 50, 0, 500)


        # Jet
        jet_pt_axis   = hist.Bin("pt",   r"Jet $p_{T}$ [GeV]", 100, 20, 400)
        jet_eta_axis  = hist.Bin("eta",  r"Jet $\eta$", 60, -3, 3)
        jet_phi_axis  = hist.Bin("phi",  r"Jet $\phi$", 60, -3, 3)
        jet_mass_axis = hist.Bin("mass", r"Jet $m$ [GeV]", 100, 0, 50)
        jet_dr_axis = hist.Bin("dr", r"Jet $\Delta$R(l,j)", 50, 0, 5)
        ljpt_axis     = hist.Bin("ljpt", r"Leading jet $p_{T}$ [GeV]", 100, 20, 400)
        ljdr_axis = hist.Bin("ljdr", "Leading jet $\Delta$R(l,j)", 50,0,5)

        

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


        # Define histograms from axes
        _hist_jet_dict = {
                'pt'  : hist.Hist("Counts", dataset_axis, cutflow_axis, flav_axis,jet_pt_axis),
                'eta' : hist.Hist("Counts", dataset_axis, cutflow_axis, flav_axis,jet_eta_axis),
                'phi' : hist.Hist("Counts", dataset_axis, cutflow_axis, flav_axis,jet_phi_axis),
                'mass': hist.Hist("Counts", dataset_axis, cutflow_axis, flav_axis,jet_mass_axis)
            }
        _hist_deepcsv_dict = {
                'pt'  : hist.Hist("Counts", dataset_axis, cutflow_axis, flav_axis,jet_pt_axis),
                'eta' : hist.Hist("Counts", dataset_axis, cutflow_axis, flav_axis,jet_eta_axis),
                'phi' : hist.Hist("Counts", dataset_axis, cutflow_axis, flav_axis,jet_phi_axis),
                'mass': hist.Hist("Counts", dataset_axis, cutflow_axis, flav_axis,jet_mass_axis),
            }


        # Generate some histograms dynamically
        for disc, axis in zip(disc_list, btag_axes):
            _hist_deepcsv_dict[disc] = hist.Hist("Counts", dataset_axis, cutflow_axis, flav_axis, axis)
        for deepcsv, axises in zip(deepcsv_list, deepcsv_axes):
             _hist_deepcsv_dict[deepcsv] = hist.Hist("Counts", dataset_axis, cutflow_axis, flav_axis,  axises)
        _hist_event_dict = {
                'njet'  : hist.Hist("Counts", dataset_axis, cutflow_axis, flav_axis, njet_axis),
                'nbjet' : hist.Hist("Counts", dataset_axis, cutflow_axis, flav_axis, nbjet_axis),
                'nmu'   : hist.Hist("Counts", dataset_axis, cutflow_axis,  nmu_axis),
                'lmupt' : hist.Hist("Counts", dataset_axis, cutflow_axis,  lmupt_axis),
                'ljpt'  : hist.Hist("Counts", dataset_axis, cutflow_axis, flav_axis, ljpt_axis),
                'sljpt'  : hist.Hist("Counts", dataset_axis, cutflow_axis, flav_axis, sljpt_axis),
                'ssljpt'  : hist.Hist("Counts", dataset_axis, cutflow_axis, flav_axis, ssljpt_axis),
                'sssljpt'  : hist.Hist("Counts", dataset_axis, cutflow_axis, flav_axis, sssljpt_axis),
                'ljdr'  : hist.Hist("Counts", dataset_axis, cutflow_axis, flav_axis, ljdr_axis),
                'sljdr'  : hist.Hist("Counts", dataset_axis, cutflow_axis, flav_axis, sljdr_axis),
                'ssljdr'  : hist.Hist("Counts", dataset_axis, cutflow_axis, flav_axis, ssljdr_axis),
                'sssljdr'  : hist.Hist("Counts", dataset_axis, cutflow_axis, flav_axis, sssljdr_axis),
                'met' : hist.Hist("Counts", dataset_axis, cutflow_axis, met_axis)
            }
        self.jet_hists = list(_hist_jet_dict.keys())
        self.deepcsv_hists = list(_hist_deepcsv_dict.keys())
        # self.deepddx_hists = list(_hist_deepddx_dict.keys())
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
        isRealData = not hasattr(events, "genWeight")
        weights = Weights(len(events), storeIndividual=True)
        selection = PackedSelection()
        if(isRealData):output['sumw'][dataset] += 1.
        else:output['sumw'][dataset] += ak.sum(events.genWeight)
        req_lumi=np.ones(len(events), dtype='bool')
        if(isRealData): req_lumi=lumiMasks['2017'](events.run, events.luminosityBlock)
        selection.add('lumimask',req_lumi)
        ##############
        # Trigger level
        triggers = [
        "HLT_IsoMu24",
        ]

        trig_arrs = [events.HLT[_trig.strip("HLT_")] for _trig in triggers]
        req_trig = np.zeros(len(events), dtype='bool')
        for t in trig_arrs:
            req_trig = req_trig | t
        selection.add('trigger',(req_trig)&(req_lumi))
        ############
        # Event level

        ## Muon cuts
        # muon twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
        event_muon = events.Muon[(events.Muon.pt > 30) & (abs(events.Muon.eta) < 2.4) & (events.Muon.tightId > .5)&(events.Muon.pfRelIso04_all<0.12)] #
        
        # print(ak.all(events.Jet.metric_table(events.Muon) > 0.4, axis=2))
        # event_muon = ak.pad_none(events.Muon, 1, axis=1) 
        
        req_muon =(ak.count(event_muon.pt, axis=1) == 1)
        
        # smu=events.Muon[mu_level]
        selection.add('musel',(req_trig)&(req_muon)&(req_lumi))
        # selection.add('musel',(req_trig)&(req_muon))

        ## Jet cuts

        event_jet = events.Jet[(events.Jet.pt > 25) & (abs(events.Jet.eta) <= 2.4)]
        # print(events.Jet.muonIdx1)
        req_jetkin = (ak.num(event_jet)>=4)
        selection.add('jetkin',(req_lumi)&(req_trig)&(req_muon)&(req_jetkin))
        event_jet_sel = events.Jet[(events.Jet.pt > 25) & (abs(events.Jet.eta) <= 2.4)&(ak.all(events.Jet.metric_table(events.Muon) > 0.4, axis=2))]
        req_jetdr = (ak.num(event_jet_sel)>=4)
        selection.add('jetdr',(req_lumi)&(req_trig)&(req_muon)&(req_jetkin)&(req_jetdr))
        event_jet_sel = events.Jet[(events.Jet.pt > 25) & (abs(events.Jet.eta) <= 2.4)&(ak.all(events.Jet.metric_table(events.Muon) > 0.4, axis=2))&(events.Jet.puId > 0)]
        req_jetpuid = (ak.num(event_jet_sel)>=4)
        selection.add('jetpuid',(req_lumi)&(req_trig)&(req_muon)&(req_jetkin)&(req_jetdr)&(req_jetpuid))
        event_jet_sel = events.Jet[(events.Jet.pt > 25) & (abs(events.Jet.eta) <= 2.4)&(ak.all(events.Jet.metric_table(events.Muon) > 0.4, axis=2))&(events.Jet.puId > 0)& (events.Jet.jetId > 0)]
        req_jetjetid = (ak.num(event_jet_sel)>=4)
        selection.add('jetjetid',(req_lumi)&(req_trig)&(req_muon)&(req_jetkin)&(req_jetdr)&(req_jetpuid)&(req_jetjetid))
        
        req_MET = events.METFixEE2017.pt>50
        selection.add('metsel',(req_lumi)&(req_trig)&(req_muon)&(req_jetkin)&(req_jetdr)&(req_jetpuid)&(req_jetjetid)&(req_MET))


        regions=['lumimask','trigger', 'musel', 'jetkin', 'jetdr', 'jetpuid', 'jetjetid','metsel']
        def flatten(ar): # flatten awkward into a 1d array to hist
            return ak.flatten(ar, axis=None)

        def num(ar):
            return ak.num(ak.fill_none(ar[~ak.is_none(ar)], 0), axis=0)
        def norm(val, cut):
            if cut is None:
                ar = ak.to_numpy(ak.fill_none(val, np.nan))
                return ar
            else:
                ar = ak.to_numpy(ak.fill_none(val[cut], np.nan))
                return ar

        if not isRealData :
            # genflavor = sjets.hadronFlavour
            weights.add('genweight',events.genWeight)
            weights.add('puweight',compiled['2017_pileupweight'](events.Pileup.nPU))
        for sel in regions:
        # # Fill histograms dynamically
        #     print(sel)
            # print(selection.names)
            
            cut = selection.all(*regions)
            for histname, h in output.items():
                if(isRealData):
                    if (histname in self.jet_hists) or (histname in self.deepcsv_hists):
                        fields = {l: ak.flatten(ak.fill_none(event_jet[cut][l],np.nan)) for l in h.fields if l in dir(event_jet)}
                        h.fill(dataset=dataset,region=sel,flav=5, **fields)
                else:
                    if (histname in self.jet_hists) or (histname in self.deepcsv_hists):
                        fields = {l: ak.flatten(ak.fill_none(event_jet[cut][l],np.nan)) for l in h.fields if l in dir(event_jet)}
                        # print(ak.type(weights.weight()[cut]))
                        # print(ak.type(sjets[cut]))
                        genweiev=ak.flatten(ak.broadcast_arrays(weights.weight()[cut],event_jet[cut][histname])[0])
                        h.fill(dataset=dataset,region=sel,flav=ak.flatten(event_jet[cut].hadronFlavour), **fields,weight=genweiev)





            if(isRealData):
                output['njet'].fill(dataset=dataset, region=sel, flav=0, njet=flatten(ak.num(event_jet[cut])))
                output['nmu'].fill(dataset=dataset,   region= sel,nmu=flatten(ak.num(event_muon[cut])))
                output['lmupt'].fill(dataset=dataset, region= sel,lmupt=flatten((event_muon[cut].pt)))
                output['ljpt'].fill(dataset=dataset, region = sel, flav=0, ljpt=flatten(event_jet[cut][:,0].pt))
                print(ak.num(event_jet))
                output['sljpt'].fill(dataset=dataset, region= sel,flav=0, sljpt=flatten(event_jet[cut][:,1].pt))
                output['ssljpt'].fill(dataset=dataset, region= sel,flav=0, ssljpt=flatten(event_jet[cut][:,2].pt))
                output['sssljpt'].fill(dataset=dataset,region= sel, flav=0, sssljpt=flatten(event_jet[cut][:,3].pt))
                output['ljdr'].fill(dataset=dataset, region= sel,flav=0, ljdr=flatten(event_jet[cut][:,0].delta_r(event_muon[cut])))
                output['sljdr'].fill(dataset=dataset, region= sel,flav=0, sljdr=flatten(event_jet[cut][:,1].delta_r(event_muon[cut])))
                output['ssljdr'].fill(dataset=dataset, region= sel,flav=0, ssljdr=flatten(event_jet[cut][:,2].delta_r(event_muon[cut])))
                output['sssljdr'].fill(dataset=dataset,region= sel, flav=0, sssljdr=flatten(event_jet[cut][:,3].delta_r(event_muon[cut])))
            else:
                
                output['njet'].fill(dataset=dataset,region=sel, flav=flatten(event_jet[cut][:,0].hadronFlavour), njet=flatten(ak.num(event_jet[cut])),weight=weights.weight()[cut])
                output['nmu'].fill(dataset=dataset,  region=sel, nmu=flatten(ak.num(event_muon[cut])),weight=weights.weight()[cut])
                output['lmupt'].fill(dataset=dataset, region=sel,lmupt=flatten((event_muon[cut].pt)),weight=weights.weight()[cut])
                output['ljpt'].fill(dataset=dataset, region=sel,flav=flatten(event_jet[cut][:,0].hadronFlavour), ljpt=flatten(event_jet[cut][:,0].pt),weight=weights.weight()[cut])
                output['sljpt'].fill(dataset=dataset, region=sel,flav=flatten(event_jet[cut][:,1].hadronFlavour), sljpt=flatten(event_jet[cut][:,1].pt),weight=weights.weight()[cut])
                output['ssljpt'].fill(dataset=dataset, region=sel,flav=flatten(event_jet[cut][:,2].hadronFlavour), ssljpt=flatten(event_jet[cut][:,2].pt),weight=weights.weight()[cut])
                output['sssljpt'].fill(dataset=dataset, region=sel,flav=flatten(event_jet[cut][:,3].hadronFlavour), sssljpt=flatten(event_jet[cut][:,3].pt),weight=weights.weight()[cut])
                output['ljdr'].fill(dataset=dataset, region=sel,flav=flatten(event_jet[cut][:,0].hadronFlavour), ljdr=flatten(event_jet[cut][:,0].delta_r(event_muon[cut])),weight=weights.weight()[cut])
                output['sljdr'].fill(dataset=dataset, region=sel,flav=flatten(event_jet[cut][:,1].hadronFlavour), sljdr=flatten(event_jet[cut][:,1].delta_r(event_muon[cut])),weight=weights.weight()[cut])
                output['ssljdr'].fill(dataset=dataset, region=sel,flav=flatten(event_jet[cut][:,2].hadronFlavour), ssljdr=flatten(event_jet[cut][:,2].delta_r(event_muon[cut])),weight=weights.weight()[cut])
                output['sssljdr'].fill(dataset=dataset, region=sel,flav=flatten(event_jet[cut][:,3].hadronFlavour), sssljdr=flatten(event_jet[cut][:,3].delta_r(event_muon[cut])),weight=weights.weight()[cut])
        return output

    def postprocess(self, accumulator):
        return accumulator
