import gc
import numpy as np
import awkward as ak
import gzip
import pickle

from coffea.lookup_tools.lookup_base import lookup_base

from coffea.lookup_tools import extractor


from coffea.lumi_tools import LumiMask
from coffea.btag_tools import BTagScaleFactor
from helpers.cTagSFReader import *
def build_lumimask(filename):
    return LumiMask(filename)


lumiMasks = {
    '2016': build_lumimask('data/lumiMasks/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'),
    '2017': build_lumimask('data/lumiMasks/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt'),
    '2018': build_lumimask('data/lumiMasks/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'),
}




'''with gzip.open("data/compiled_jec.pkl.gzt") as fin:
    jmestuff = cloudpickle.load(fin)

jet_factory = jmestuff["jet_factory"]
def add_jec_variables(jets, event_rho):
    print("========> begin jec var <==========")
    #objgraph.show_growth()
    jets["pt_raw"] = (1 - jets.rawFactor)*jets.pt
    jets["mass_raw"] = (1 - jets.rawFactor)*jets.mass
    jets["pt_gen"] = ak.values_astype(ak.fill_none(jets.matched_gen.pt, 0), np.float32)
    jets["event_rho"] = ak.broadcast_arrays(event_rho, jets.pt)[0]

    gc.collect()
    print("====>after gc<=======")
    #objgraph.show_growth()
    return jets'''
<<<<<<< HEAD
<<<<<<< HEAD
with gzip.open("data/Rereco17_94X/94XPUwei_corrections.pkl.gz") as fin:
    compiled = pickle.load(fin)

deepcsvb_sf = BTagScaleFactor("data/Rereco17_94X/BTag/DeepCSV_94XSF_V5_B_F.csv",BTagScaleFactor.RESHAPE,methods='iterativefit,iterativefit,iterativefit')
deepcsvc_sf = "data/Rereco17_94X/BTag/DeepCSV_ctagSF_MiniAOD94X_2017_pTincl_v3_2_interp.root"
deepjetb_sf = BTagScaleFactor("data/Rereco17_94X/BTag/DeepFlavour_94XSF_V4_B_F.csv",BTagScaleFactor.RESHAPE,methods='iterativefit,iterativefit,iterativefit')
deepjetc_sf = "data/Rereco17_94X/BTag/DeepJet_ctagSF_MiniAOD94X_2017_pTincl_v3_2_interp.root"

### Lepton SFs
ext = extractor()
ext.add_weight_sets(["ele_Trig TrigSF data/Rereco17_94X/Lepton/Ele32_L1DoubleEG_TrigSF_vhcc.histo.root"])
ext.add_weight_sets(["ele_ID EGamma_SF2D data/Rereco17_94X/Lepton/ElectronIDSF_94X_MVA80WP.histo.root"])
ext.add_weight_sets(["ele_Rereco EGamma_SF2D data/Rereco17_94X/Lepton/ElectronRecoSF_94X.histo.root"])
# ext.add_weight_sets(["mu_TrigSF IsoMu27_PtEtaBins/pt_abseta_ratio data/Rereco17_94X/singleMuonTrig.histo.root"])
ext.add_weight_sets(["mu_ID NUM_TightID_DEN_genTracks_pt_abseta data/Rereco17_94X/Lepton/RunBCDEF_SF_ID.histo.root"])
ext.add_weight_sets(["mu_ID_low NUM_TightID_DEN_genTracks_pt_abseta data/Rereco17_94X/Lepton/RunBCDEF_SF_MuID_lowpT.histo.root"])
ext.add_weight_sets(["mu_Iso NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta data/Rereco17_94X/Lepton/RunBCDEF_SF_ISO.histo.root"])
=======
=======
>>>>>>> 55c5cc4594c4800d0b7e32d1443ee6a4241a0d36
with gzip.open("data/Rereco17_94X/corrections.pkl.gz") as fin:
    compiled = pickle.load(fin)

deepcsvb_sf = BTagScaleFactor("data/Rereco17_94X/DeepCSV_94XSF_V5_B_F.csv",BTagScaleFactor.RESHAPE,methods='iterativefit,iterativefit,iterativefit')
deepcsvc_sf = "data/Rereco17_94X/DeepCSV_ctagSF_MiniAOD94X_2017_pTincl_v3_2_interp.root"
deepjetb_sf = BTagScaleFactor("data/Rereco17_94X/DeepFlavour_94XSF_V4_B_F.csv",BTagScaleFactor.RESHAPE,methods='iterativefit,iterativefit,iterativefit')
deepjetc_sf = "data/Rereco17_94X/DeepJet_ctagSF_MiniAOD94X_2017_pTincl_v3_2_interp.root"

### Lepton SFs
ext = extractor()
ext.add_weight_sets(["ele_Trig TrigSF data/Rereco17_94X/Ele32_L1DoubleEG_TrigSF_vhcc.histo.root"])
ext.add_weight_sets(["ele_ID EGamma_SF2D data/Rereco17_94X/ElectronIDSF_94X_MVA80WP.histo.root"])
ext.add_weight_sets(["ele_Rereco EGamma_SF2D data/Rereco17_94X/ElectronRecoSF_94X.histo.root"])
# ext.add_weight_sets(["mu_TrigSF IsoMu27_PtEtaBins/pt_abseta_ratio data/Rereco17_94X/singleMuonTrig.histo.root"])
ext.add_weight_sets(["mu_ID NUM_TightID_DEN_genTracks_pt_abseta data/Rereco17_94X/RunBCDEF_SF_ID.histo.root"])
ext.add_weight_sets(["mu_ID_low NUM_TightID_DEN_genTracks_pt_abseta data/Rereco17_94X/RunBCDEF_SF_MuID_lowpT.histo.root"])
ext.add_weight_sets(["mu_Iso NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta data/Rereco17_94X/RunBCDEF_SF_ISO.histo.root"])
<<<<<<< HEAD
>>>>>>> 55c5cc4594c4800d0b7e32d1443ee6a4241a0d36
=======
>>>>>>> 55c5cc4594c4800d0b7e32d1443ee6a4241a0d36
ext.finalize()
evaluator = ext.make_evaluator()

def eleSFs(ele):
    weight = evaluator["ele_Trig"](ele.eta,ele.pt)*evaluator["ele_ID"](ele.eta,ele.pt)*evaluator["ele_Rereco"](ele.eta,ele.pt)

    return weight

def muSFs(mu):
    
    # weight 
    mu_eta=ak.fill_none(mu.eta,0.)
    mu_pt= ak.fill_none(mu.pt,0.)
    weight = evaluator["mu_ID"](mu_eta,mu_pt)*evaluator["mu_ID_low"](mu_eta,mu_pt)*evaluator["mu_Iso"](mu_eta,mu_pt)

    return weight
