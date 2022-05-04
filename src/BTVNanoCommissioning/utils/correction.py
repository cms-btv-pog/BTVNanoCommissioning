import importlib
import gzip
import pickle
import contextlib
import cloudpickle

import numpy as np
import awkward as ak
from coffea.lookup_tools.lookup_base import lookup_base
from coffea.lookup_tools import extractor


from coffea.lumi_tools import LumiMask
from coffea.btag_tools import BTagScaleFactor
from coffea.jetmet_tools import JECStack, CorrectedJetsFactory
from coffea.lookup_tools import extractor
# from BTVNanoCommissioning.helpers.cTagSFReader import getSF

# FIXME make jsons campaing configurable
lumiMasks = {}
_lumi_path = "BTVNanoCommissioning.data.lumiMasks"
with importlib.resources.path(_lumi_path, "Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt") as filename:
    lumiMasks["2016"] = LumiMask(filename)
with importlib.resources.path(_lumi_path, "Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt") as filename:
    lumiMasks["2017"] = LumiMask(filename)
with importlib.resources.path(_lumi_path, "Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt") as filename:
    lumiMasks["2018"] = LumiMask(filename)

##JEC
with importlib.resources.path("BTVNanoCommissioning.data.JME.Rereco17_94X", "jec_compiled.pkl.gz") as path:
    with gzip.open(path) as fin:
        jmestuff = cloudpickle.load(fin)

# with gzip.open("data/JME/Rereco17_94X/jec_compiled.pkl.gz") as fin:
#     jmestuff = pickle.load(fin)
jet_factory = jmestuff["jet_factory"]

def add_jec_variables(jets, event_rho):
    jets["pt_raw"] = (1 - jets.rawFactor)*jets.pt
    jets["mass_raw"] = (1 - jets.rawFactor)*jets.mass
    if hasattr(jets, "genJetIdxG"):jets["pt_gen"] = ak.values_astype(ak.fill_none(jets.matched_gen.pt, 0), np.float32)
    jets["event_rho"] = ak.broadcast_arrays(event_rho, jets.pt)[0]
    return jets

## PU weight
with importlib.resources.path("BTVNanoCommissioning.data.PU.Rereco17_94X", "94XPUwei_corrections.pkl.gz") as path:
    with gzip.open(path) as fin:
        compiled = cloudpickle.load(fin)

## BTag SFs
_btag_path = "BTVNanoCommissioning.data.BTV.Rereco17_94X"
with importlib.resources.path(_btag_path, "DeepCSV_94XSF_V5_B_F.csv") as filename:
    deepcsvb_sf = BTagScaleFactor(filename, BTagScaleFactor.RESHAPE,
                                  methods='iterativefit,iterativefit,iterativefit')
with importlib.resources.path(_btag_path, "DeepFlavour_94XSF_V4_B_F.csv") as filename:
    deepjetb_sf = BTagScaleFactor(filename, BTagScaleFactor.RESHAPE,
                                  methods='iterativefit,iterativefit,iterativefit')

deepcsvc_sf = "DeepCSV_ctagSF_MiniAOD94X_2017_pTincl_v3_2_interp.root"
deepjetc_sf = "DeepJet_ctagSF_MiniAOD94X_2017_pTincl_v3_2_interp.root"

### Lepton SFs
ext = extractor()
ele_sf_mapping = {
    "ele_Trig TrigSF": "Ele32_L1DoubleEG_TrigSF_vhcc.histo.root",
    "ele_ID EGamma_SF2D": "ElectronIDSF_94X_MVA80WP.histo.root",
    "ele_Rereco EGamma_SF2D": "ElectronRecoSF_94X.histo.root",
    "mu_ID NUM_TightID_DEN_genTracks_pt_abseta": "RunBCDEF_SF_ID.histo.root",
    "mu_ID_low NUM_TightID_DEN_genTracks_pt_abseta": "RunBCDEF_SF_MuID_lowpT.histo.root",
    "mu_Iso NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta": "RunBCDEF_SF_ISO.histo.root",
}

with contextlib.ExitStack() as stack:
    # this would work even in zipballs but since extractor keys on file extension and
    # importlib make a random tempfile, it won't work. coffea needs to enable specifying the type manually
    # for now we run this whole module as $ python -m boostedhiggs.build_jec boostedhiggs/data/jec_compiled.pkl.gz
    # so the compiled value can be loaded using the importlib tool in corrections.py
    _ele_path = "BTVNanoCommissioning.data.LSF.Rereco17_94X"
    real_paths = [stack.enter_context(importlib.resources.path(_ele_path, f)) for f in ele_sf_mapping.values()]
    ext.add_weight_sets([f"{path} {file}" for path, file in zip(ele_sf_mapping.keys(), real_paths)])
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
