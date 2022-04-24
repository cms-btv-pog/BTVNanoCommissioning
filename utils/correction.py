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
from coffea.jetmet_tools import JECStack, CorrectedJetsFactory
from coffea.lookup_tools import extractor

lumiMasks = {
    '2016': LumiMask('data/lumiMasks/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'),
    '2017': LumiMask('data/lumiMasks/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt'),
    '2018': LumiMask('data/lumiMasks/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'),
}
##JEC 
with gzip.open("data/JME/Rereco17_94X/jec_compiled.pkl.gz") as fin:
    jmestuff = pickle.load(fin)
jet_factory = jmestuff["jet_factory"]

def add_jec_variables(jets, event_rho):
    jets["pt_raw"] = (1 - jets.rawFactor)*jets.pt
    jets["mass_raw"] = (1 - jets.rawFactor)*jets.mass
    if hasattr(jets, "genJetIdxG"):jets["pt_gen"] = ak.values_astype(ak.fill_none(jets.matched_gen.pt, 0), np.float32)
    jets["event_rho"] = ak.broadcast_arrays(event_rho, jets.pt)[0]
    return jets
## PU weight
with gzip.open("data/PU/Rereco17_94X/94XPUwei_corrections.pkl.gz") as fin:
    compiled = pickle.load(fin)

## BTag SFs
deepcsvb_sf = BTagScaleFactor("data/BTV/Rereco17_94X/DeepCSV_94XSF_V5_B_F.csv",BTagScaleFactor.RESHAPE,methods='iterativefit,iterativefit,iterativefit')
deepcsvc_sf = "data/BTV/Rereco17_94X/DeepCSV_ctagSF_MiniAOD94X_2017_pTincl_v3_2_interp.root"
deepjetb_sf = BTagScaleFactor("data/BTV/Rereco17_94X/DeepFlavour_94XSF_V4_B_F.csv",BTagScaleFactor.RESHAPE,methods='iterativefit,iterativefit,iterativefit')
deepjetc_sf = "data/BTV/Rereco17_94X/DeepJet_ctagSF_MiniAOD94X_2017_pTincl_v3_2_interp.root"

### Lepton SFs
ext = extractor()
ext.add_weight_sets(["ele_Trig TrigSF data/LSF/Rereco17_94X/Ele32_L1DoubleEG_TrigSF_vhcc.histo.root"])
ext.add_weight_sets(["ele_ID EGamma_SF2D data/LSF/Rereco17_94X/ElectronIDSF_94X_MVA80WP.histo.root"])
ext.add_weight_sets(["ele_Rereco EGamma_SF2D data/LSF/Rereco17_94X/ElectronRecoSF_94X.histo.root"])
ext.add_weight_sets(["mu_ID NUM_TightID_DEN_genTracks_pt_abseta data/LSF/Rereco17_94X/RunBCDEF_SF_ID.histo.root"])
ext.add_weight_sets(["mu_ID_low NUM_TightID_DEN_genTracks_pt_abseta data/LSF/Rereco17_94X/RunBCDEF_SF_MuID_lowpT.histo.root"])
ext.add_weight_sets(["mu_Iso NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta data/LSF/Rereco17_94X/RunBCDEF_SF_ISO.histo.root"])

ext.finalize()
evaluator = ext.make_evaluator()

def eleSFs(ele):
    ele_eta = ak.fill_none(ele.eta,0.)
    ele_pt = ak.fill_none(ele.pt,0.)

    weight =  evaluator["ele_ID"](ele_eta,ele_pt)*evaluator["ele_Rereco"](ele_eta,ele_pt)  
    return weight

def muSFs(mu):
    
    # weight 
    mu_eta=ak.fill_none(mu.eta,0.)
    mu_pt= ak.fill_none(mu.pt,0.)
    weight = evaluator["mu_ID"](mu_eta,mu_pt)*evaluator["mu_ID_low"](mu_eta,mu_pt)*evaluator["mu_Iso"](mu_eta,mu_pt)

    return weight
