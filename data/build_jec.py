import importlib.resources
import contextlib
from coffea.lookup_tools import extractor
from coffea.jetmet_tools import JECStack, CorrectedJetsFactory, CorrectedMETFactory

jec_name_map = {
    'JetPt': 'pt',
    'JetMass': 'mass',
    'JetEta': 'eta',
    'JetA': 'area',
    'ptGenJet': 'pt_gen',
    'ptRaw': 'pt_raw',
    'massRaw': 'mass_raw',
    'Rho': 'event_rho',
    'METpt': 'pt',
    'METphi': 'phi',
    'JetPhi': 'phi',
    'UnClusteredEnergyDeltaX': 'MetUnclustEnUpDeltaX',
    'UnClusteredEnergyDeltaY': 'MetUnclustEnUpDeltaY',
}


def jet_factory_factory(files):
    ext = extractor()
    with contextlib.ExitStack() as stack:
        # this would work even in zipballs but since extractor keys on file extension and
        # importlib make a random tempfile, it won't work. coffea needs to enable specifying the type manually
        # for now we run this whole module as $ python -m boostedhiggs.build_jec boostedhiggs/data/jec_compiled.pkl.gzt
        # so the compiled value can be loaded using the importlib tool in corrections.py
        #real_files = [stack.enter_context(importlib.resources.path("data","Summer19UL17_JRV2_MC_SF_AK4PFchs.txt")) for f in files]

        real_files = [stack.enter_context(open(f)) for f in files]

    
        ext.add_weight_sets([f"* * {f}" for f in files])
        ext.finalize()

    jec_stack = JECStack(ext.make_evaluator())
    return CorrectedJetsFactory(jec_name_map, jec_stack)


jet_factory = {

    "2017mc": jet_factory_factory(
        files=[
            # https://raw.githubusercontent.com/cms-jet/JECDatabase/master/textFiles/Summer19UL17_V5_MC/Summer19UL17_V5_MC_L1FastJet_AK4PFchs.txt
            "Summer19UL17_V5_MC_L1FastJet_AK4PFchs.jec.txt",
            # https://github.com/cms-jet/JECDatabase/blob/master/textFiles/Summer19UL17_V5_MC/Summer19UL17_V5_MC_L2Relative_AK4PFchs.txt
            "Summer19UL17_V5_MC_L2Relative_AK4PFchs.jec.txt",
            # https://github.com/cms-jet/JECDatabase/blob/master/textFiles/Summer19UL17_V5_MC/RegroupedV2_Summer19UL17_V5_MC_UncertaintySources_AK4PFchs.txt
            "RegroupedV2_Summer19UL17_V5_MC_UncertaintySources_AK4PFchs.junc.txt",
            # https://github.com/cms-jet/JECDatabase/blob/master/textFiles/Summer19UL17_V5_MC/Summer19UL17_V5_MC_Uncertainty_AK4PFchs.txt
            "Summer19UL17_V5_MC_Uncertainty_AK4PFchs.junc.txt",
            # https://github.com/cms-jet/JRDatabase/raw/master/textFiles/Fall17_V3b_MC/Fall17_V3b_MC_PtResolution_AK4PFchs.txt
            "Summer19UL17_JRV2_MC_PtResolution_AK4PFchs.jr.txt",
            "Summer19UL17_JRV2_MC_SF_AK4PFchs.jersf.txt"
        ]
    ),
"2017data": jet_factory_factory(
        files=[
            "Summer19UL17_RunB_V5_DATA_L1FastJet_AK4PFchs.jec.txt",
            "Summer19UL17_RunB_V5_DATA_L2Relative_AK4PFchs.jec.txt",
            "Summer19UL17_RunB_V5_DATA_L3Absolute_AK4PFchs.jec.txt",
            "Summer19UL17_RunB_V5_DATA_Uncertainty_AK4PFchs.junc.txt",
            "Summer19UL17_JRV2_DATA_PtResolution_AK4PFchs.jr.txt",
            "Summer19UL17_JRV2_DATA_SF_AK4PFchs.jersf.txt",
        ])
}







if __name__ == "__main__":
    import sys
    import gzip
    # jme stuff not pickleable in coffea
    import cloudpickle

    with gzip.open(sys.argv[-1], "wb") as fout:

        cloudpickle.dump(
            {
                "jet_factory": jet_factory
            },
            fout
        )