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
# with gzip.open("data/compiled_jec.pkl.gzt") as fin:
#     jmestuff = cloudpickle.load(fin)
ext = extractor()
ext.add_weight_sets([
    "* * data/Rereco17_94X/JEC_JERSF/MC/Fall17_17Nov2017_V32_MC_L1FastJet_AK4PFchs.jec.txt",
    "* * data/Rereco17_94X/JEC_JERSF/MC/Fall17_17Nov2017_V32_MC_L2Relative_AK4PFchs.jec.txt",
    "* * data/Rereco17_94X/JEC_JERSF/MC/Fall17_17Nov2017_V32_MC_L3Absolute_AK4PFchs.jec.txt",
    "* * data/Rereco17_94X/JEC_JERSF/MC/Fall17_17Nov2017_V32_MC_L2L3Residual_AK4PFchs.jec.txt",
])

ext.finalize()
jec_stack_names = [
                    "Fall17_17Nov2017_V32_MC_L1FastJet_AK4PFchs",
                   "Fall17_17Nov2017_V32_MC_L2Relative_AK4PFchs",
                   "Fall17_17Nov2017_V32_MC_L3Absolute_AK4PFchs",
                   "Fall17_17Nov2017_V32_MC_L2L3Residual_AK4PFchs",
                   ]
evaluator = ext.make_evaluator()
jec_inputs = {name: evaluator[name] for name in jec_stack_names}
jec_stack = JECStack(jec_inputs)
name_map = jec_stack.blank_name_map
name_map['JetPt'] = 'pt'
name_map['JetMass'] = 'mass'
name_map['JetEta'] = 'eta'
name_map['JetA'] = 'area'

def add_jec_variables(events,jets, event_rho):
    jets["pt_raw"] = (1 - jets.rawFactor)*jets.pt
    jets["mass_raw"] = (1 - jets.rawFactor)*jets.mass
    jets["pt_gen"] = ak.values_astype(ak.fill_none(jets.matched_gen.pt, 0), np.float32)
    jets["event_rho"] = ak.broadcast_arrays(event_rho, jets.pt)[0]
    jet_factory = CorrectedJetsFactory(name_map, jec_stack)
    events_cache = events.caches[0]
    corrected_jets = jet_factory.build(jets, lazy_cache=events_cache)
    gc.collect()
    return corrected_jets
## PU weight
with gzip.open("data/Rereco17_94X/94XPUwei_corrections.pkl.gz") as fin:
    compiled = pickle.load(fin)

## BTag SFs
deepcsvb_sf = BTagScaleFactor("data/Rereco17_94X/BTag/DeepCSV_94XSF_V5_B_F.csv",BTagScaleFactor.RESHAPE,methods='iterativefit,iterativefit,iterativefit')
deepcsvc_sf = "data/Rereco17_94X/BTag/DeepCSV_ctagSF_MiniAOD94X_2017_pTincl_v3_2_interp.root"
deepjetb_sf = BTagScaleFactor("data/Rereco17_94X/BTag/DeepFlavour_94XSF_V4_B_F.csv",BTagScaleFactor.RESHAPE,methods='iterativefit,iterativefit,iterativefit')
deepjetc_sf = "data/Rereco17_94X/BTag/DeepJet_ctagSF_MiniAOD94X_2017_pTincl_v3_2_interp.root"

### Lepton SFs
ext = extractor()
ext.add_weight_sets(["ele_Trig TrigSF data/Rereco17_94X/Lepton/Ele32_L1DoubleEG_TrigSF_vhcc.histo.root"])
ext.add_weight_sets(["ele_ID EGamma_SF2D data/Rereco17_94X/Lepton/ElectronIDSF_94X_MVA80WP.histo.root"])
ext.add_weight_sets(["ele_Rereco EGamma_SF2D data/Rereco17_94X/Lepton/ElectronRecoSF_94X.histo.root"])
ext.add_weight_sets(["mu_ID NUM_TightID_DEN_genTracks_pt_abseta data/Rereco17_94X/Lepton/RunBCDEF_SF_ID.histo.root"])
ext.add_weight_sets(["mu_ID_low NUM_TightID_DEN_genTracks_pt_abseta data/Rereco17_94X/Lepton/RunBCDEF_SF_MuID_lowpT.histo.root"])
ext.add_weight_sets(["mu_Iso NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta data/Rereco17_94X/Lepton/RunBCDEF_SF_ISO.histo.root"])

ext.finalize()
evaluator = ext.make_evaluator()

def eleSFs(ele):
    ele_eta = ak.fill_none(ele.eta,0.)
    ele_pt = ak.fill_none(ele.pt,0.)
    print(ak.type(ele_eta))
    #weight = evaluator["ele_Trig"](ele_eta,ele_pt)*evaluator["ele_ID"](ele_eta,ele_pt)*evaluator["ele_Rereco"](ele_eta,ele_pt)
    weight =  evaluator["ele_ID"](ele_eta,ele_pt)*evaluator["ele_Rereco"](ele_eta,ele_pt)  
    return weight

def muSFs(mu):
    
    # weight 
    mu_eta=ak.fill_none(mu.eta,0.)
    mu_pt= ak.fill_none(mu.pt,0.)
    weight = evaluator["mu_ID"](mu_eta,mu_pt)*evaluator["mu_ID_low"](mu_eta,mu_pt)*evaluator["mu_Iso"](mu_eta,mu_pt)

    return weight
