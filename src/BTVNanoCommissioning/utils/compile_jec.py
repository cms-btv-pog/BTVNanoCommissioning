import importlib.resources
import contextlib
from coffea.lookup_tools import extractor
from coffea.jetmet_tools import JECStack, CorrectedJetsFactory

jec_name_map = {
    "JetPt": "pt",
    "JetMass": "mass",
    "JetEta": "eta",
    "JetA": "area",
    "ptGenJet": "pt_gen",
    "ptRaw": "pt_raw",
    "massRaw": "mass_raw",
    "Rho": "event_rho",
}


def jet_factory_factory(files):
    ext = extractor()
    ext.add_weight_sets([f"* * {file}" for file in files])
    ext.finalize()
    jec_stack = JECStack(ext.make_evaluator())

    return CorrectedJetsFactory(jec_name_map, jec_stack)


def jet_factories(campaign):
    jet_factory = {
        "Rereco17_94X": {
            "mc": jet_factory_factory(
                files=[
                    "data/JME/Rereco17_94X/MC/Fall17_V3b_MC_PtResolution_AK4PFchs.jr.txt",
                    "data/JME/Rereco17_94X/MC/Fall17_V3b_MC_SF_AK4PFchs.jersf.txt",
                    "data/JME/Rereco17_94X/MC/Fall17_17Nov2017_V32_MC_L1FastJet_AK4PFchs.jec.txt",
                    "data/JME/Rereco17_94X/MC/Fall17_17Nov2017_V32_MC_L2Relative_AK4PFchs.jec.txt",
                    "data/JME/Rereco17_94X/MC/Fall17_17Nov2017_V32_MC_L3Absolute_AK4PFchs.jec.txt",
                    "data/JME/Rereco17_94X/MC/Fall17_17Nov2017_V32_MC_L2L3Residual_AK4PFchs.jec.txt",
                ]
            ),
            "data": jet_factory_factory(
                files=[
                    "data/JME/Rereco17_94X/RunDE/Fall17_17Nov2017DE_V32_DATA_L1FastJet_AK4PFchs.jec.txt",
                    "data/JME/Rereco17_94X/RunDE/Fall17_17Nov2017DE_V32_DATA_L2Relative_AK4PFchs.jec.txt",
                    "data/JME/Rereco17_94X/RunDE/Fall17_17Nov2017DE_V32_DATA_L3Absolute_AK4PFchs.jec.txt",
                ]
            ),
        },
        "UL17_106X": {
            "mc": jet_factory_factory(
                files=[
                    "data/JME/UL17_106X/MC/Summer19UL17_JRV2_MC_PtResolution_AK4PFchs.jr.txt",
                    "data/JME/UL17_106X/MC/Summer19UL17_JRV2_MC_SF_AK4PFchs.jersf.txt",
                    "data/JME/UL17_106X/MC/Summer19UL17_V5_MC_L1FastJet_AK4PFchs.jec.txt",
                    "data/JME/UL17_106X/MC/Summer19UL17_V5_MC_L2Relative_AK4PFchs.jec.txt",
                    "data/JME/UL17_106X/MC/Summer19UL17_V5_MC_L3Absolute_AK4PFchs.jec.txt",
                    "data/JME/UL17_106X/MC/Summer19UL17_V5_MC_L2L3Residual_AK4PFchs.jec.txt",
                ]
            ),
            "dataB": jet_factory_factory(
                files=[
                    "data/JME/UL17_106X/RunB/Summer19UL17_RunB_V5_DATA_L1FastJet_AK4PFchs.jec.txt",
                    "data/JME/UL17_106X/RunB/Summer19UL17_RunB_V5_DATA_L2Relative_AK4PFchs.jec.txt",
                    "data/JME/UL17_106X/RunB/Summer19UL17_RunB_V5_DATA_L3Absolute_AK4PFchs.jec.txt",
                ]
            ),
            "dataC": jet_factory_factory(
                files=[
                    "data/JME/UL17_106X/RunC/Summer19UL17_RunC_V5_DATA_L1FastJet_AK4PFchs.jec.txt",
                    "data/JME/UL17_106X/RunC/Summer19UL17_RunC_V5_DATA_L2Relative_AK4PFchs.jec.txt",
                    "data/JME/UL17_106X/RunC/Summer19UL17_RunC_V5_DATA_L3Absolute_AK4PFchs.jec.txt",
                ]
            ),
            "dataD": jet_factory_factory(
                files=[
                    "data/JME/UL17_106X/RunD/Summer19UL17_RunD_V5_DATA_L1FastJet_AK4PFchs.jec.txt",
                    "data/JME/UL17_106X/RunD/Summer19UL17_RunD_V5_DATA_L2Relative_AK4PFchs.jec.txt",
                    "data/JME/UL17_106X/RunD/Summer19UL17_RunD_V5_DATA_L3Absolute_AK4PFchs.jec.txt",
                ]
            ),
            "dataE": jet_factory_factory(
                files=[
                    "data/JME/UL17_106X/RunE/Summer19UL17_RunE_V5_DATA_L1FastJet_AK4PFchs.jec.txt",
                    "data/JME/UL17_106X/RunE/Summer19UL17_RunE_V5_DATA_L2Relative_AK4PFchs.jec.txt",
                    "data/JME/UL17_106X/RunE/Summer19UL17_RunE_V5_DATA_L3Absolute_AK4PFchs.jec.txt",
                ]
            ),
            "dataF": jet_factory_factory(
                files=[
                    "data/JME/UL17_106X/RunF/Summer19UL17_RunF_V5_DATA_L1FastJet_AK4PFchs.jec.txt",
                    "data/JME/UL17_106X/RunF/Summer19UL17_RunF_V5_DATA_L2Relative_AK4PFchs.jec.txt",
                    "data/JME/UL17_106X/RunF/Summer19UL17_RunF_V5_DATA_L3Absolute_AK4PFchs.jec.txt",
                ]
            ),
        },
    }
    return jet_factory[campaign]


if __name__ == "__main__":
    import sys
    import gzip

    # jme stuff not pickleable in coffea
    import cloudpickle

    campaign = sys.argv[-2]
    with gzip.open(sys.argv[-1], "wb") as fout:
        cloudpickle.dump(
            {
                "jet_factory": jet_factories(campaign),
            },
            fout,
        )
