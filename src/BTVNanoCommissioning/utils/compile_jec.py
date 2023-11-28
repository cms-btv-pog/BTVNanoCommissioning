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
        "Rereco17_94X": {
            "mc": jet_factory_factory(
                files=[
                    "src/BTVNanoCommissioning/data/JME/Rereco17_94X/MC/Fall17_V3b_MC_PtResolution_AK4PFchs.jr.txt",
                    "src/BTVNanoCommissioning/data/JME/Rereco17_94X/MC/Fall17_V3b_MC_SF_AK4PFchs.jersf.txt",
                    "src/BTVNanoCommissioning/data/JME/Rereco17_94X/MC/Fall17_17Nov2017_V32_MC_L1FastJet_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Rereco17_94X/MC/Fall17_17Nov2017_V32_MC_L2Relative_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Rereco17_94X/MC/Fall17_17Nov2017_V32_MC_L3Absolute_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Rereco17_94X/MC/Fall17_17Nov2017_V32_MC_L2L3Residual_AK4PFchs.jec.txt",
                ]
            ),
            "data": jet_factory_factory(
                files=[
                    "src/BTVNanoCommissioning/data/JME/Rereco17_94X/RunDE/Fall17_17Nov2017DE_V32_DATA_L1FastJet_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Rereco17_94X/RunDE/Fall17_17Nov2017DE_V32_DATA_L2Relative_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Rereco17_94X/RunDE/Fall17_17Nov2017DE_V32_DATA_L3Absolute_AK4PFchs.jec.txt",
                ]
            ),
        },
        "2016preVFP_UL": {
            "mc": jet_factory_factory(
                files=[
                    "src/BTVNanoCommissioning/data/JME/2016preVFP_UL/Summer19UL16APV_V7_MC_L1FastJet_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2016preVFP_UL/Summer19UL16APV_V7_MC_L2Relative_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2016preVFP_UL/Summer19UL16APV_V7_MC_L3Absolute_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2016preVFP_UL/Summer20UL16APV_JRV3_DATA_PtResolution_AK4PFchs.jr.txt",
                    "src/BTVNanoCommissioning/data/JME/2016preVFP_UL/Summer20UL16APV_JRV3_MC_SF_AK4PFchs.jersf.txt",
                    "src/BTVNanoCommissioning/data/JME/2016preVFP_UL/RegroupedV2_Summer19UL16APV_V7_MC_UncertaintySources_AK4PFchs.junc.txt",
                    "src/BTVNanoCommissioning/data/JME/2016preVFP_UL/Summer19UL16APV_V7_MC_Uncertainty_AK4PFchs.junc.txt",
                ]
            ),
            "dataBCD": jet_factory_factory(
                files=[
                    "src/BTVNanoCommissioning/data/JME/2016preVFP_UL/Summer19UL16APV_RunBCD_V7_DATA_L1FastJet_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2016preVFP_UL/Summer19UL16APV_RunBCD_V7_DATA_L2L3Residual_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2016preVFP_UL/Summer19UL16APV_RunBCD_V7_DATA_L2Relative_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2016preVFP_UL/Summer19UL16APV_RunBCD_V7_DATA_L3Absolute_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2016preVFP_UL/Summer20UL16APV_JRV3_DATA_SF_AK4PFchs.jersf.txt",
                    "src/BTVNanoCommissioning/data/JME/2016preVFP_UL/Summer20UL16APV_JRV3_DATA_PtResolution_AK4PFchs.jr.txt",
                    "src/BTVNanoCommissioning/data/JME/2016preVFP_UL/Summer19UL16APV_RunBCD_V7_DATA_UncertaintySources_AK4PFchs.junc.txt",
                    "src/BTVNanoCommissioning/data/JME/2016preVFP_UL/Summer19UL16APV_RunBCD_V7_DATA_Uncertainty_AK4PFchs.junc.txt",
                ]
            ),
            "dataEF": jet_factory_factory(
                files=[
                    "src/BTVNanoCommissioning/data/JME/2016preVFP_UL/Summer19UL16APV_RunEF_V7_DATA_L1FastJet_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2016preVFP_UL/Summer19UL16APV_RunEF_V7_DATA_L2L3Residual_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2016preVFP_UL/Summer19UL16APV_RunEF_V7_DATA_L2Relative_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2016preVFP_UL/Summer19UL16APV_RunEF_V7_DATA_L3Absolute_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2016preVFP_UL/Summer20UL16APV_JRV3_DATA_SF_AK4PFchs.jersf.txt",
                    "src/BTVNanoCommissioning/data/JME/2016preVFP_UL/Summer20UL16APV_JRV3_DATA_PtResolution_AK4PFchs.jr.txt",
                    "src/BTVNanoCommissioning/data/JME/2016preVFP_UL/Summer19UL16APV_RunEF_V7_DATA_Uncertainty_AK4PFchs.junc.txt",
                    "src/BTVNanoCommissioning/data/JME/2016preVFP_UL/Summer19UL16APV_RunEF_V7_DATA_UncertaintySources_AK4PFchs.junc.txt",
                ]
            ),
        },
        "2016postVFP_UL": {
            "mc": jet_factory_factory(
                files=[
                    "src/BTVNanoCommissioning/data/JME/2016postVFP_UL/RegroupedV2_Summer19UL16_V7_MC_UncertaintySources_AK4PFchs.junc.txt",
                    "src/BTVNanoCommissioning/data/JME/2016postVFP_UL/Summer19UL16_V7_MC_L1FastJet_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2016postVFP_UL/Summer19UL16_V7_MC_L2Relative_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2016postVFP_UL/Summer19UL16_V7_MC_L3Absolute_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2016postVFP_UL/Summer19UL16_V7_MC_Uncertainty_AK4PFchs.junc.txt",
                    "src/BTVNanoCommissioning/data/JME/2016postVFP_UL/Summer20UL16_JRV3_MC_PtResolution_AK4PFchs.jr.txt",
                    "src/BTVNanoCommissioning/data/JME/2016postVFP_UL/Summer20UL16_JRV3_MC_SF_AK4PFchs.jersf.txt",
                ]
            ),
            "dataFGH": jet_factory_factory(
                files=[
                    "src/BTVNanoCommissioning/data/JME/2016postVFP_UL/Summer20UL16_JRV3_DATA_PtResolution_AK4PFchs.jr.txt",
                    "src/BTVNanoCommissioning/data/JME/2016postVFP_UL/Summer20UL16_JRV3_DATA_SF_AK4PFchs.jersf.txt",
                    "src/BTVNanoCommissioning/data/JME/2016postVFP_UL/Summer19UL16_RunFGH_V7_DATA_L1FastJet_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2016postVFP_UL/Summer19UL16_RunFGH_V7_DATA_L2L3Residual_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2016postVFP_UL/Summer19UL16_RunFGH_V7_DATA_L2Relative_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2016postVFP_UL/Summer19UL16_RunFGH_V7_DATA_L3Absolute_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2016postVFP_UL/Summer19UL16_RunFGH_V7_DATA_Uncertainty_AK4PFchs.junc.txt",
                    "src/BTVNanoCommissioning/data/JME/2016postVFP_UL/Summer19UL16_RunFGH_V7_DATA_UncertaintySources_AK4PFchs.junc.txt",
                ]
            ),
        },
        "2017_UL": {
            "mc": jet_factory_factory(
                files=[
                    "src/BTVNanoCommissioning/data/JME/2017_UL/RegroupedV2_Summer19UL17_V5_MC_UncertaintySources_AK4PFchs.junc.txt",
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_JRV2_MC_PtResolution_AK4PFchs.jr.txt",
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_JRV2_MC_SF_AK4PFchs.jersf.txt",
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_V5_MC_L1FastJet_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_V5_MC_L2Relative_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_V5_MC_L3Absolute_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_V5_MC_Uncertainty_AK4PFchs.junc.txt",
                ]
            ),
            "dataB": jet_factory_factory(
                files=[
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_RunB_V5_DATA_L1FastJet_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_RunB_V5_DATA_L2L3Residual_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_RunB_V5_DATA_L2Relative_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_RunB_V5_DATA_L3Absolute_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_RunB_V5_DATA_Uncertainty_AK4PFchs.junc.txt",
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_RunB_V5_DATA_UncertaintySources_AK4PFchs.junc.txt",
                ]
            ),
            "dataC": jet_factory_factory(
                files=[
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_RunC_V5_DATA_L1FastJet_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_RunC_V5_DATA_L2L3Residual_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_RunC_V5_DATA_L2Relative_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_RunC_V5_DATA_L3Absolute_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_RunC_V5_DATA_Uncertainty_AK4PFchs.junc.txt",
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_RunC_V5_DATA_UncertaintySources_AK4PFchs.junc.txt",
                ]
            ),
            "dataD": jet_factory_factory(
                files=[
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_RunD_V5_DATA_L1FastJet_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_RunD_V5_DATA_L2L3Residual_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_RunD_V5_DATA_L2Relative_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_RunD_V5_DATA_L3Absolute_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_RunD_V5_DATA_Uncertainty_AK4PFchs.junc.txt",
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_RunD_V5_DATA_UncertaintySources_AK4PFchs.junc.txt",
                ]
            ),
            "dataE": jet_factory_factory(
                files=[
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_RunE_V5_DATA_L1FastJet_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_RunE_V5_DATA_L2L3Residual_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_RunE_V5_DATA_L2Relative_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_RunE_V5_DATA_L3Absolute_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_RunE_V5_DATA_Uncertainty_AK4PFchs.junc.txt",
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_RunE_V5_DATA_UncertaintySources_AK4PFchs.junc.txt",
                ]
            ),
            "dataF": jet_factory_factory(
                files=[
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_RunF_V5_DATA_L1FastJet_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_RunF_V5_DATA_L2L3Residual_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_RunF_V5_DATA_L2Relative_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_RunF_V5_DATA_L3Absolute_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_RunF_V5_DATA_Uncertainty_AK4PFchs.junc.txt",
                    "src/BTVNanoCommissioning/data/JME/2017_UL/Summer19UL17_RunF_V5_DATA_UncertaintySources_AK4PFchs.junc.txt",
                ]
            ),
        },
        "2018_UL": {
            "mc": jet_factory_factory(
                files=[
                    "src/BTVNanoCommissioning/data/JME/2018_UL/RegroupedV2_Summer19UL18_V5_MC_UncertaintySources_AK4PFchs.junc.txt",
                    "src/BTVNanoCommissioning/data/JME/2018_UL/Summer19UL18_JRV2_MC_PtResolution_AK4PFchs.jr.txt",
                    "src/BTVNanoCommissioning/data/JME/2018_UL/Summer19UL18_JRV2_MC_SF_AK4PFchs.jersf.txt",
                    "src/BTVNanoCommissioning/data/JME/2018_UL/Summer19UL18_V5_MC_L1FastJet_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2018_UL/Summer19UL18_V5_MC_L2L3Residual_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2018_UL/Summer19UL18_V5_MC_L2Relative_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2018_UL/Summer19UL18_V5_MC_L3Absolute_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2018_UL/Summer19UL18_V5_MC_Uncertainty_AK4PFchs.junc.txt",
                ]
            ),
            "dataA": jet_factory_factory(
                files=[
                    "src/BTVNanoCommissioning/data/JME/2018_UL/Summer19UL18_RunA_V5_DATA_L1FastJet_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2018_UL/Summer19UL18_RunA_V5_DATA_L2L3Residual_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2018_UL/Summer19UL18_RunA_V5_DATA_L2Relative_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2018_UL/Summer19UL18_RunA_V5_DATA_L3Absolute_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2018_UL/Summer19UL18_RunA_V5_DATA_Uncertainty_AK4PFchs.junc.txt",
                    "src/BTVNanoCommissioning/data/JME/2018_UL/Summer19UL18_RunA_V5_DATA_UncertaintySources_AK4PFchs.junc.txt",
                ]
            ),
            "dataB": jet_factory_factory(
                files=[
                    "src/BTVNanoCommissioning/data/JME/2018_UL/Summer19UL18_RunB_V5_DATA_L1FastJet_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2018_UL/Summer19UL18_RunB_V5_DATA_L2L3Residual_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2018_UL/Summer19UL18_RunB_V5_DATA_L2Relative_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2018_UL/Summer19UL18_RunB_V5_DATA_L3Absolute_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2018_UL/Summer19UL18_RunB_V5_DATA_Uncertainty_AK4PFchs.junc.txt",
                    "src/BTVNanoCommissioning/data/JME/2018_UL/Summer19UL18_RunB_V5_DATA_UncertaintySources_AK4PFchs.junc.txt",
                ]
            ),
            "dataC": jet_factory_factory(
                files=[
                    "src/BTVNanoCommissioning/data/JME/2018_UL/Summer19UL18_RunC_V5_DATA_L1FastJet_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2018_UL/Summer19UL18_RunC_V5_DATA_L2L3Residual_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2018_UL/Summer19UL18_RunC_V5_DATA_L2Relative_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2018_UL/Summer19UL18_RunC_V5_DATA_L3Absolute_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2018_UL/Summer19UL18_RunC_V5_DATA_Uncertainty_AK4PFchs.junc.txt",
                    "src/BTVNanoCommissioning/data/JME/2018_UL/Summer19UL18_RunC_V5_DATA_UncertaintySources_AK4PFchs.junc.txt",
                ]
            ),
            "dataD": jet_factory_factory(
                files=[
                    "src/BTVNanoCommissioning/data/JME/2018_UL/Summer19UL18_RunD_V5_DATA_L1FastJet_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2018_UL/Summer19UL18_RunD_V5_DATA_L2L3Residual_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2018_UL/Summer19UL18_RunD_V5_DATA_L2Relative_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2018_UL/Summer19UL18_RunD_V5_DATA_L3Absolute_AK4PFchs.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/2018_UL/Summer19UL18_RunD_V5_DATA_Uncertainty_AK4PFchs.junc.txt",
                    "src/BTVNanoCommissioning/data/JME/2018_UL/Summer19UL18_RunD_V5_DATA_UncertaintySources_AK4PFchs.junc.txt",
                ]
            ),
        },
        "Winter22Run3": {
            "mc": jet_factory_factory(
                files=[
                    "src/BTVNanoCommissioning/data/JME/Winter22Run3/Winter22Run3_V1_MC_PtResolution_AK4PFPuppi.jr.txt",
                    "src/BTVNanoCommissioning/data/JME/Winter22Run3/Winter22Run3_V1_MC_SF_AK4PFPuppi.jersf.txt",
                    "src/BTVNanoCommissioning/data/JME/Winter22Run3/Winter22Run3_V2_MC_L1FastJet_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Winter22Run3/Winter22Run3_V2_MC_L2Relative_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Winter22Run3/Winter22Run3_V2_MC_L2Residual_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Winter22Run3/Winter22Run3_V2_MC_L3Absolute_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Winter22Run3/Winter22Run3_V2_MC_UncertaintySources_AK4PFPuppi.junc.txt",
                    "src/BTVNanoCommissioning/data/JME/Winter22Run3/Winter22Run3_V2_MC_Uncertainty_AK4PFPuppi.junc.txt",
                ]
            ),
            "dataC": jet_factory_factory(
                files=[
                    "src/BTVNanoCommissioning/data/JME/Winter22Run3/Winter22Run3_RunC_V2_DATA_L1FastJet_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Winter22Run3/Winter22Run3_RunC_V2_DATA_L2Relative_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Winter22Run3/Winter22Run3_RunC_V2_DATA_L3Absolute_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Winter22Run3/Winter22Run3_RunC_V2_DATA_L2L3Residual_AK4PFPuppi.jec.txt",
                ]
            ),
            "dataD": jet_factory_factory(
                files=[
                    "src/BTVNanoCommissioning/data/JME/Winter22Run3/Winter22Run3_RunD_V2_DATA_L1FastJet_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Winter22Run3/Winter22Run3_RunD_V2_DATA_L2Relative_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Winter22Run3/Winter22Run3_RunD_V2_DATA_L3Absolute_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Winter22Run3/Winter22Run3_RunD_V2_DATA_L2L3Residual_AK4PFPuppi.jec.txt",
                ]
            ),
        },
        "Summer22Run3": {
            "mc": jet_factory_factory(
                files=[
                    "src/BTVNanoCommissioning/data/JME/Summer22Run3/Summer22_V1_MC_L1FastJet_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Summer22Run3/Summer22_V1_MC_L2Relative_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Summer22Run3/Summer22_V1_MC_L2Residual_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Summer22Run3/Summer22_V1_MC_L3Absolute_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Summer22Run3/Summer22_V1_MC_UncertaintySources_AK4PFPuppi.junc.txt",
                    "src/BTVNanoCommissioning/data/JME/Summer22Run3/Summer22_V1_MC_Uncertainty_AK4PFPuppi.junc.txt",
                ]
            ),
            "dataC": jet_factory_factory(
                files=[
                    "src/BTVNanoCommissioning/data/JME/Summer22Run3/Summer22_RunCD_V1_DATA_L1FastJet_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Summer22Run3/Summer22_RunCD_V1_DATA_L2Relative_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Summer22Run3/Summer22_RunCD_V1_DATA_L3Absolute_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Summer22Run3/Summer22_RunCD_V1_DATA_L2L3Residual_AK4PFPuppi.jec.txt",
                ]
            ),
            "dataD": jet_factory_factory(
                files=[
                    "src/BTVNanoCommissioning/data/JME/Summer22Run3/Summer22_RunCD_V1_DATA_L1FastJet_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Summer22Run3/Summer22_RunCD_V1_DATA_L2Relative_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Summer22Run3/Summer22_RunCD_V1_DATA_L3Absolute_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Summer22Run3/Summer22_RunCD_V1_DATA_L2L3Residual_AK4PFPuppi.jec.txt",
                ]
            ),
        },
        "Summer22EERun3": {
            "mc": jet_factory_factory(
                files=[
                    "src/BTVNanoCommissioning/data/JME/Summer22EERun3/Summer22EEPrompt22_JRV1_MC_SF_AK4PFPuppi.jersf.txt",
                    "src/BTVNanoCommissioning/data/JME/Summer22EERun3/Summer22EEPrompt22_JRV1_MC_PtResolution_AK4PFPuppi.jr.txt",
                    "src/BTVNanoCommissioning/data/JME/Summer22EERun3/Summer22EE_V1_MC_L1FastJet_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Summer22EERun3/Summer22EE_V1_MC_L2Relative_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Summer22EERun3/Summer22EE_V1_MC_L2Residual_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Summer22EERun3/Summer22EE_V1_MC_L3Absolute_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Summer22EERun3/Summer22EE_V1_MC_UncertaintySources_AK4PFPuppi.junc.txt",
                    "src/BTVNanoCommissioning/data/JME/Summer22EERun3/Summer22EE_V1_MC_Uncertainty_AK4PFPuppi.junc.txt",
                ]
            ),
            "dataE": jet_factory_factory(
                files=[
                    "src/BTVNanoCommissioning/data/JME/Summer22EERun3/Summer22EE_RunE_V1_DATA_L1FastJet_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Summer22EERun3/Summer22EE_RunE_V1_DATA_L2Relative_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Summer22EERun3/Summer22EE_RunE_V1_DATA_L3Absolute_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Summer22EERun3/Summer22EE_RunE_V1_DATA_L2L3Residual_AK4PFPuppi.jec.txt",
                ]
            ),
            "dataF": jet_factory_factory(
                files=[
                    "src/BTVNanoCommissioning/data/JME/Summer22EERun3/Summer22EE_RunF_V1_DATA_L1FastJet_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Summer22EERun3/Summer22EE_RunF_V1_DATA_L2Relative_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Summer22EERun3/Summer22EE_RunF_V1_DATA_L3Absolute_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Summer22EERun3/Summer22EE_RunF_V1_DATA_L2L3Residual_AK4PFPuppi.jec.txt",
                ]
            ),
            "dataG": jet_factory_factory(
                files=[
                    "src/BTVNanoCommissioning/data/JME/Summer22EERun3/Summer22EE_RunG_V1_DATA_L1FastJet_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Summer22EERun3/Summer22EE_RunG_V1_DATA_L2Relative_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Summer22EERun3/Summer22EE_RunG_V1_DATA_L3Absolute_AK4PFPuppi.jec.txt",
                    "src/BTVNanoCommissioning/data/JME/Summer22EERun3/Summer22EE_RunG_V1_DATA_L2L3Residual_AK4PFPuppi.jec.txt",
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
    with gzip.open(
        f"src/BTVNanoCommissioning/data/JME/{campaign}/{sys.argv[-1]}.pkl.gz", "wb"
    ) as fout:
        cloudpickle.dump(
            {
                "jet_factory": jet_factories(campaign),
                "met_factory": CorrectedMETFactory(jec_name_map),
            },
            fout,
        )
