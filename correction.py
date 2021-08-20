import gc
import numpy as np
import awkward as ak
import gzip
import pickle
import cloudpickle
import importlib.resources
import objgraph

from coffea.lookup_tools.lookup_base import lookup_base
from coffea import lookup_tools
from coffea import util

from coffea.lumi_tools import LumiMask

def build_lumimask(filename):
    return LumiMask(filename)


lumiMasks = {
    '2016': build_lumimask('data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'),
    '2017': build_lumimask('data/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt'),
    '2018': build_lumimask('data/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'),
}




with gzip.open("data/compiled_jec.pkl.gzt") as fin:
    jmestuff = cloudpickle.load(fin)

jet_factory = jmestuff["jet_factory"]
def add_jec_variables(jets, event_rho):
    print("========> begin jec var <==========")
    objgraph.show_growth()
    jets["pt_raw"] = (1 - jets.rawFactor)*jets.pt
    jets["mass_raw"] = (1 - jets.rawFactor)*jets.mass
    jets["pt_gen"] = ak.values_astype(ak.fill_none(jets.matched_gen.pt, 0), np.float32)
    jets["event_rho"] = ak.broadcast_arrays(event_rho, jets.pt)[0]

    gc.collect()
    print("====>after gc<=======")
    objgraph.show_growth()
    return jets
with gzip.open("data/corrections.pkl.gz") as fin:
    compiled = pickle.load(fin)
