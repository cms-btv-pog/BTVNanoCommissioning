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
        njet_axis  = hist.Bin("njet",  r"N jets",      [0,1,2,3,4,5,6,7,8,9,10])
        nbjet_axis = hist.Bin("nbjet", r"N b-jets",    [0,1,2,3,4,5,6,7,8,9,10])
        met_axis = hist.Bin("met",   r"Missing ET", 50, 0, 500)


        # Jet
        jet_pt_axis   = hist.Bin("pt",   r"Jet $p_{T}$ [GeV]", 100, 20, 400)
        jet_eta_axis  = hist.Bin("eta",  r"Jet $\eta$", 60, -3, 3)
        jet_phi_axis  = hist.Bin("phi",  r"Jet $\phi$", 60, -3, 3)
        jet_mass_axis = hist.Bin("mass", r"Jet $m$ [GeV]", 100, 0, 50)
        ljpt_axis     = hist.Bin("ljpt", r"Leading jet $p_{T}$ [GeV]", 100, 20, 400)

        

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
                'ljpt'  : hist.Hist("Counts", dataset_axis, cutflow_axis, flav_axis, ljpt_axis),
                # 'met' : hist.Hist("Counts", dataset_axis, cutflow_axis, met_axis)
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
        #triggers = [
        #"HLT_PFJet40",
        #]

        #trig_arrs = [events.HLT[_trig.strip("HLT_")] for _trig in triggers]
        req_hlt40 = np.zeros(len(events), dtype='bool')
        req_hlt60 = np.zeros(len(events), dtype='bool')
        req_hlt80 = np.zeros(len(events), dtype='bool')
        req_hlt40 = req_hlt40|events.HLT['PFJet40']
        req_hlt60 = req_hlt60|events.HLT['PFJet60']
        req_hlt80 = req_hlt80|events.HLT['PFJet80']
        #for t in trig_arrs:
         #   req_trig = req_trig | t
        selection.add('hlt40',(req_hlt40)&(req_lumi))
        selection.add('hlt60',(req_hlt60)&(req_lumi))
        selection.add('hlt80',(req_hlt80)&(req_lumi))
        ############

        ## Jet cuts
        req_jetkin40 = (ak.num((events.Jet.pt > 50) & (abs(events.Jet.eta) <= 2.4)&(events.Jet.puId > 0)& (events.Jet.jetId > 0))>=1)
        selection.add('jetkin40',req_jetkin40)
        req_jetkin60 = (ak.num((events.Jet.pt > 70) & (abs(events.Jet.eta) <= 2.4)&(events.Jet.puId > 0)& (events.Jet.jetId > 0))>=1)
        selection.add('jetkin60',req_jetkin60)
        req_jetkin80 = (ak.num((events.Jet.pt > 100) & (abs(events.Jet.eta) <= 2.4)&(events.Jet.puId > 0)& (events.Jet.jetId > 0))>=1)
        selection.add('jetkin80',req_jetkin80)

        
        
        regions={
            "PFJet40" : ['lumimask','hlt40','jetkin40'],
            "PFJet60" : ['lumimask','hlt60','jetkin60'],
            "PFJet80" : ['lumimask','hlt80','jetkin80']
        }

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
        #for sel in regions:
        for region, cuts in regions.items():
        # # Fill histograms dynamically
        #     print(sel)
            # print(selection.names)
            selections = regions[region]
            cut = selection.all(*selections)
            for histname, h in output.items():
                if(isRealData):
                    if (histname in self.jet_hists) or (histname in self.deepcsv_hists):
                        fields = {l: ak.flatten(ak.fill_none(events.Jet[cut][l],np.nan)) for l in h.fields if l in dir(events.Jet)}
                        h.fill(dataset=dataset,region=region,flav=5, **fields)
                else:
                    if (histname in self.jet_hists) or (histname in self.deepcsv_hists):
                        fields = {l: ak.flatten(ak.fill_none(events.Jet[cut][l],np.nan)) for l in h.fields if l in dir(events.Jet)}
                        # print(ak.type(weights.weight()[cut]))
                        # print(ak.type(sjets[cut]))
                        genweiev=ak.flatten(ak.broadcast_arrays(weights.weight()[cut],events.Jet[cut][histname])[0])
                        h.fill(dataset=dataset,region=region,flav=ak.flatten(events.Jet[cut].hadronFlavour), **fields,weight=genweiev)





            if(isRealData):
                output['njet'].fill(dataset=dataset, region=region, flav=0, njet=flatten(ak.num(events.Jet[cut])))
                print(events.Jet[cut][:,0].pt)
                # output['ljpt'].fill(dataset=dataset, region =region, flav=0, ljpt=flatten(events.Jet[cut][:,0].pt))
                # print(ak.num(events.Jet))
            else:
                
                output['njet'].fill(dataset=dataset,region=region, flav=flatten(events.Jet[cut][:,0].hadronFlavour), njet=flatten(ak.num(events.Jet[cut])),weight=weights.weight()[cut])
                # output['ljpt'].fill(dataset=dataset, region=region, flav=events.Jet)
                output['ljpt'].fill(dataset=dataset, region=region,flav=events.Jet[cut][:,0].hadronFlavour, ljpt=events.Jet[cut][:,0].pt,weight=weights.weight()[cut])
              
        return output

    def postprocess(self, accumulator):
        return accumulator
