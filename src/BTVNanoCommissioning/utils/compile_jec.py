import importlib.resources
import contextlib
from coffea.lookup_tools import extractor
from coffea.jetmet_tools import JECStack, CorrectedJetsFactory, CorrectedMETFactory

jec_name_map = {
    "JetPt": "pt",
    "JetMass": "mass",
    "JetEta": "eta",
    "JetA": "area",
    "ptGenJet": "pt_gen",
    "ptRaw": "pt_raw",
    "massRaw": "mass_raw",
    "Rho": "event_rho",
    "METpt": "pt",
    "METphi": "phi",
    "JetPhi": "phi",
    "UnClusteredEnergyDeltaX": "MetUnclustEnUpDeltaX",
    "UnClusteredEnergyDeltaY": "MetUnclustEnUpDeltaY",
}


def jet_factory_factory(files):
    ext = extractor()
    ext.add_weight_sets([f"* * {file}" for file in files])
    ext.finalize()
    jec_stack = JECStack(ext.make_evaluator())

    return CorrectedJetsFactory(jec_name_map, jec_stack)


def jet_factories(campaign):
    jet_factory = {
        # "Rereco17_94X": {
        #     "mc": jet_factory_factory(
        #         files=[
        #             ""data/JME/Rereco17_94X/MC/Fall17_V3b_MC_PtResolution_AK4PFchs.jr.txt",",
        #             ""data/JME/Rereco17_94X/MC/Fall17_V3b_MC_SF_AK4PFchs.jersf.txt",",
        #             ""data/JME/Rereco17_94X/MC/Fall17_17Nov2017_V32_MC_L1FastJet_AK4PFchs.jec.txt",",
        #             ""data/JME/Rereco17_94X/MC/Fall17_17Nov2017_V32_MC_L2Relative_AK4PFchs.jec.txt",",
        #             ""data/JME/Rereco17_94X/MC/Fall17_17Nov2017_V32_MC_L3Absolute_AK4PFchs.jec.txt",",
        #             ""data/JME/Rereco17_94X/MC/Fall17_17Nov2017_V32_MC_L2L3Residual_AK4PFchs.jec.txt",",
        #         ]
        #     ),
        #     "data": jet_factory_factory(
        #         files=[
        #             ""data/JME/Rereco17_94X/RunDE/Fall17_17Nov2017DE_V32_DATA_L1FastJet_AK4PFchs.jec.txt",",
        #             ""data/JME/Rereco17_94X/RunDE/Fall17_17Nov2017DE_V32_DATA_L2Relative_AK4PFchs.jec.txt",",
        #             ""data/JME/Rereco17_94X/RunDE/Fall17_17Nov2017DE_V32_DATA_L3Absolute_AK4PFchs.jec.txt",",
        #         ]
        #     ),
        # },
        "UL16_preAPV": {
            "mc": jet_factory_factory(
                files=[
                    "data/JME/UL16_preAPV/Summer19UL16APV_V7_MC_L1FastJet_AK4PFchs.jec.txt",
                    "data/JME/UL16_preAPV/Summer19UL16APV_V7_MC_L2Relative_AK4PFchs.jec.txt",
                    "data/JME/UL16_preAPV/Summer19UL16APV_V7_MC_L3Absolute_AK4PFchs.jec.txt",
                    "data/JME/UL16_preAPV/Summer20UL16APV_JRV3_DATA_PtResolution_AK4PFchs.jr.txt",
                    "data/JME/UL16_preAPV/Summer20UL16APV_JRV3_MC_SF_AK4PFchs.jersf.txt",
                    "data/JME/UL16_preAPV/RegroupedV2_Summer19UL16APV_V7_MC_UncertaintySources_AK4PFchs.junc.txt",
                    "data/JME/UL16_preAPV/Summer19UL16APV_V7_MC_Uncertainty_AK4PFchs.junc.txt",
                ]
            ),
            "dataBCD": jet_factory_factory(
                files=[
                    "data/JME/UL16_preAPV/Summer19UL16APV_RunBCD_V7_DATA_L1FastJet_AK4PFchs.jec.txt",
                    "data/JME/UL16_preAPV/Summer19UL16APV_RunBCD_V7_DATA_L2L3Residual_AK4PFchs.jec.txt",
                    "data/JME/UL16_preAPV/Summer19UL16APV_RunBCD_V7_DATA_L2Relative_AK4PFchs.jec.txt",
                    "data/JME/UL16_preAPV/Summer19UL16APV_RunBCD_V7_DATA_L3Absolute_AK4PFchs.jec.txt",
                    "data/JME/UL16_preAPV/Summer20UL16APV_JRV3_DATA_SF_AK4PFchs.jersf.txt",
                    "data/JME/UL16_preAPV/Summer20UL16APV_JRV3_DATA_PtResolution_AK4PFchs.jr.txt",
                    # "data/JME/UL16_preAPV/Summer19UL16APV_RunBCD_V7_DATA_UncertaintySources_AK4PFchs.junc.txt",
                    # "data/JME/UL16_preAPV/Summer19UL16APV_RunBCD_V7_DATA_Uncertainty_AK4PFchs.junc.txt",
                ]
            ),
            "dataEF": jet_factory_factory(
                files=[
                    "data/JME/UL16_preAPV/Summer19UL16APV_RunEF_V7_DATA_L1FastJet_AK4PFchs.jec.txt",
                    "data/JME/UL16_preAPV/Summer19UL16APV_RunEF_V7_DATA_L2L3Residual_AK4PFchs.jec.txt",
                    "data/JME/UL16_preAPV/Summer19UL16APV_RunEF_V7_DATA_L2Relative_AK4PFchs.jec.txt",
                    "data/JME/UL16_preAPV/Summer19UL16APV_RunEF_V7_DATA_L3Absolute_AK4PFchs.jec.txt",
                    "data/JME/UL16_preAPV/Summer20UL16APV_JRV3_DATA_SF_AK4PFchs.jersf.txt",
                    "data/JME/UL16_preAPV/Summer20UL16APV_JRV3_DATA_PtResolution_AK4PFchs.jr.txt",
                    # "data/JME/UL16/Summer19UL16APV_RunEF_V7_DATA_Uncertainty_AK4PFchs.junc.txt",
                    # "data/JME/UL16/Summer19UL16APV_RunEF_V7_DATA_UncertaintySources_AK4PFchs.junc.txt",
                ]
            ),
        },
        "UL16_postAPV": {
            "mc": jet_factory_factory(
                files=[
                    "data/JME/UL16_postAPV/RegroupedV2_Summer19UL16_V7_MC_UncertaintySources_AK4PFchs.junc.txt",
                    "data/JME/UL16_postAPV/Summer19UL16_V7_MC_L1FastJet_AK4PFchs.jec.txt",
                    "data/JME/UL16_postAPV/Summer19UL16_V7_MC_L2Relative_AK4PFchs.jec.txt",
                    "data/JME/UL16_postAPV/Summer19UL16_V7_MC_L3Absolute_AK4PFchs.jec.txt",
                    "data/JME/UL16_postAPV/Summer19UL16_V7_MC_Uncertainty_AK4PFchs.junc.txt",
                    "data/JME/UL16_postAPV/Summer20UL16_JRV3_MC_PtResolution_AK4PFchs.jr.txt",
                    "data/JME/UL16_postAPV/Summer20UL16_JRV3_MC_SF_AK4PFchs.jersf.txt",
                ]
            ),
            "dataFGH": jet_factory_factory(
                files=[
                    "data/JME/UL16_postAPV/Summer20UL16_JRV3_DATA_PtResolution_AK4PFchs.jr.txt",
                    "data/JME/UL16_postAPV/Summer20UL16_JRV3_DATA_SF_AK4PFchs.jersf.txt",
                    "data/JME/UL16_postAPV/Summer19UL16_RunFGH_V7_DATA_L1FastJet_AK4PFchs.jec.txt",
                    "data/JME/UL16_postAPV/Summer19UL16_RunFGH_V7_DATA_L2L3Residual_AK4PFchs.jec.txt",
                    "data/JME/UL16_postAPV/Summer19UL16_RunFGH_V7_DATA_L2Relative_AK4PFchs.jec.txt",
                    "data/JME/UL16_postAPV/Summer19UL16_RunFGH_V7_DATA_L3Absolute_AK4PFchs.jec.txt",
                    # "data/JME/UL16/Summer19UL16_RunFGH_V7_DATA_Uncertainty_AK4PFchs.junc.txt",
                    # "data/JME/UL16/Summer19UL16_RunFGH_V7_DATA_UncertaintySources_AK4PFchs.junc.txt",
                ]
            ),
        },
        "UL17": {
            "mc": jet_factory_factory(
                files=[
                    "data/JME/UL17/RegroupedV2_Summer19UL17_V5_MC_UncertaintySources_AK4PFchs.junc.txt",
                    "data/JME/UL17/Summer19UL17_JRV2_MC_PtResolution_AK4PFchs.jr.txt",
                    "data/JME/UL17/Summer19UL17_JRV2_MC_SF_AK4PFchs.jersf.txt",
                    "data/JME/UL17/Summer19UL17_V5_MC_L1FastJet_AK4PFchs.jec.txt",
                    "data/JME/UL17/Summer19UL17_V5_MC_L2Relative_AK4PFchs.jec.txt",
                    "data/JME/UL17/Summer19UL17_V5_MC_L3Absolute_AK4PFchs.jec.txt",
                    "data/JME/UL17/Summer19UL17_V5_MC_Uncertainty_AK4PFchs.junc.txt",
                ]
            ),
            # "dataB": jet_factory_factory(
            #     files=[
            #         "data/JME/UL17/Summer19UL17_RunB_V5_DATA_L1FastJet_AK4PFchs.jec.txt",
            #         "data/JME/UL17/Summer19UL17_RunB_V5_DATA_L2L3Residual_AK4PFchs.jec.txt",
            #         "data/JME/UL17/Summer19UL17_RunB_V5_DATA_L2Relative_AK4PFchs.jec.txt",
            #         "data/JME/UL17/Summer19UL17_RunB_V5_DATA_L3Absolute_AK4PFchs.jec.txt",
            #         # "data/JME/UL17/Summer19UL17_RunB_V5_DATA_Uncertainty_AK4PFchs.junc.txt",
            #         # "data/JME/UL17/Summer19UL17_RunB_V5_DATA_UncertaintySources_AK4PFchs.junc.txt",
            #     ]
            # ),
            # "dataC": jet_factory_factory(
            #     files=[
            #         "data/JME/UL17/Summer19UL17_RunC_V5_DATA_L1FastJet_AK4PFchs.jec.txt",
            #         "data/JME/UL17/Summer19UL17_RunC_V5_DATA_L2L3Residual_AK4PFchs.jec.txt",
            #         "data/JME/UL17/Summer19UL17_RunC_V5_DATA_L2Relative_AK4PFchs.jec.txt",
            #         "data/JME/UL17/Summer19UL17_RunC_V5_DATA_L3Absolute_AK4PFchs.jec.txt",
            #         # "data/JME/UL17/Summer19UL17_RunC_V5_DATA_Uncertainty_AK4PFchs.junc.txt",
            #         # "data/JME/UL17/Summer19UL17_RunC_V5_DATA_UncertaintySources_AK4PFchs.junc.txt",
            #     ]
            # ),
            # "dataD": jet_factory_factory(
            #     files=[
            #         "data/JME/UL17/Summer19UL17_RunD_V5_DATA_L1FastJet_AK4PFchs.jec.txt",
            #         "data/JME/UL17/Summer19UL17_RunD_V5_DATA_L2L3Residual_AK4PFchs.jec.txt",
            #         "data/JME/UL17/Summer19UL17_RunD_V5_DATA_L2Relative_AK4PFchs.jec.txt",
            #         "data/JME/UL17/Summer19UL17_RunD_V5_DATA_L3Absolute_AK4PFchs.jec.txt",
            #         # "data/JME/UL17/Summer19UL17_RunD_V5_DATA_Uncertainty_AK4PFchs.junc.txt",
            #         # "data/JME/UL17/Summer19UL17_RunD_V5_DATA_UncertaintySources_AK4PFchs.junc.txt",
            #     ]
            # ),
            # "dataE": jet_factory_factory(
            #     files=[
            #         "data/JME/UL17/Summer19UL17_RunE_V5_DATA_L1FastJet_AK4PFchs.jec.txt",
            #         "data/JME/UL17/Summer19UL17_RunE_V5_DATA_L2L3Residual_AK4PFchs.jec.txt",
            #         "data/JME/UL17/Summer19UL17_RunE_V5_DATA_L2Relative_AK4PFchs.jec.txt",
            #         "data/JME/UL17/Summer19UL17_RunE_V5_DATA_L3Absolute_AK4PFchs.jec.txt",
            #         # "data/JME/UL17/Summer19UL17_RunE_V5_DATA_Uncertainty_AK4PFchs.junc.txt",
            #         # "data/JME/UL17/Summer19UL17_RunE_V5_DATA_UncertaintySources_AK4PFchs.junc.txt",
            #     ]
            # ),
            # "dataF": jet_factory_factory(
            #     files=[
            #         "data/JME/UL17/Summer19UL17_RunF_V5_DATA_L1FastJet_AK4PFchs.jec.txt",
            #         "data/JME/UL17/Summer19UL17_RunF_V5_DATA_L2L3Residual_AK4PFchs.jec.txt",
            #         "data/JME/UL17/Summer19UL17_RunF_V5_DATA_L2Relative_AK4PFchs.jec.txt",
            #         "data/JME/UL17/Summer19UL17_RunF_V5_DATA_L3Absolute_AK4PFchs.jec.txt",
            #         # "data/JME/UL17/Summer19UL17_RunF_V5_DATA_Uncertainty_AK4PFchs.junc.txt",
            #         # "data/JME/UL17/Summer19UL17_RunF_V5_DATA_UncertaintySources_AK4PFchs.junc.txt",
            #     ]
            # ),
        },
        "UL18": {
            "mc": jet_factory_factory(
                files=[
                    "data/JME/UL18/RegroupedV2_Summer19UL18_V5_MC_UncertaintySources_AK4PFchs.junc.txt",
                    "data/JME/UL18/Summer19UL18_JRV2_MC_PtResolution_AK4PFchs.jr.txt",
                    "data/JME/UL18/Summer19UL18_JRV2_MC_SF_AK4PFchs.jersf.txt",
                    "data/JME/UL18/Summer19UL18_V5_MC_L1FastJet_AK4PFchs.jec.txt",
                    "data/JME/UL18/Summer19UL18_V5_MC_L2L3Residual_AK4PFchs.jec.txt",
                    "data/JME/UL18/Summer19UL18_V5_MC_L2Relative_AK4PFchs.jec.txt",
                    "data/JME/UL18/Summer19UL18_V5_MC_L3Absolute_AK4PFchs.jec.txt",
                    "data/JME/UL18/Summer19UL18_V5_MC_Uncertainty_AK4PFchs.junc.txt",
                ]
            ),
            "dataA": jet_factory_factory(
                files=[
                    "data/JME/UL18/Summer19UL18_RunA_V5_DATA_L1FastJet_AK4PFchs.jec.txt",
                    "data/JME/UL18/Summer19UL18_RunA_V5_DATA_L2L3Residual_AK4PFchs.jec.txt",
                    "data/JME/UL18/Summer19UL18_RunA_V5_DATA_L2Relative_AK4PFchs.jec.txt",
                    "data/JME/UL18/Summer19UL18_RunA_V5_DATA_L3Absolute_AK4PFchs.jec.txt",
                    # "data/JME/UL18/Summer19UL18_RunA_V5_DATA_Uncertainty_AK4PFchs.junc.txt",
                    # "data/JME/UL18/Summer19UL18_RunA_V5_DATA_UncertaintySources_AK4PFchs.junc.txt",
                ]
            ),
            "dataB": jet_factory_factory(
                files=[
                    "data/JME/UL18/Summer19UL18_RunB_V5_DATA_L1FastJet_AK4PFchs.jec.txt",
                    "data/JME/UL18/Summer19UL18_RunB_V5_DATA_L2L3Residual_AK4PFchs.jec.txt",
                    "data/JME/UL18/Summer19UL18_RunB_V5_DATA_L2Relative_AK4PFchs.jec.txt",
                    "data/JME/UL18/Summer19UL18_RunB_V5_DATA_L3Absolute_AK4PFchs.jec.txt",
                    # "data/JME/UL18/Summer19UL18_RunB_V5_DATA_Uncertainty_AK4PFchs.junc.txt",
                    # "data/JME/UL18/Summer19UL18_RunB_V5_DATA_UncertaintySources_AK4PFchs.junc.txt",
                ]
            ),
            "dataC": jet_factory_factory(
                files=[
                    "data/JME/UL18/Summer19UL18_RunC_V5_DATA_L1FastJet_AK4PFchs.jec.txt",
                    "data/JME/UL18/Summer19UL18_RunC_V5_DATA_L2L3Residual_AK4PFchs.jec.txt",
                    "data/JME/UL18/Summer19UL18_RunC_V5_DATA_L2Relative_AK4PFchs.jec.txt",
                    "data/JME/UL18/Summer19UL18_RunC_V5_DATA_L3Absolute_AK4PFchs.jec.txt",
                    # "data/JME/UL18/Summer19UL18_RunC_V5_DATA_Uncertainty_AK4PFchs.junc.txt",
                    # "data/JME/UL18/Summer19UL18_RunC_V5_DATA_UncertaintySources_AK4PFchs.junc.txt",
                ]
            ),
            "dataD": jet_factory_factory(
                files=[
                    "data/JME/UL18/Summer19UL18_RunD_V5_DATA_L1FastJet_AK4PFchs.jec.txt",
                    "data/JME/UL18/Summer19UL18_RunD_V5_DATA_L2L3Residual_AK4PFchs.jec.txt",
                    "data/JME/UL18/Summer19UL18_RunD_V5_DATA_L2Relative_AK4PFchs.jec.txt",
                    "data/JME/UL18/Summer19UL18_RunD_V5_DATA_L3Absolute_AK4PFchs.jec.txt",
                    # "data/JME/UL18/Summer19UL18_RunD_V5_DATA_Uncertainty_AK4PFchs.junc.txt",
                    # "data/JME/UL18/Summer19UL18_RunD_V5_DATA_UncertaintySources_AK4PFchs.junc.txt",
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
                "met_factory": CorrectedMETFactory(jec_name_map),
            },
            fout,
        )
