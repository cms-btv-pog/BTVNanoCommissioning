correction_config = {
    "Rereco17_94X": {
        "lumiMask": "Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt",
        "PU": "94XPUwei_corrections.pkl.gz",
        "JME": "jec_compiled.pkl.gz",
        "BTV": {
            "DeepCSVB": "DeepCSV_94XSF_V5_B_F.csv",
            "DeepCSVC": "DeepCSV_ctagSF_MiniAOD94X_2017_pTincl_v3_2_interp.root",
            "DeepJetB": "DeepFlavour_94XSF_V4_B_F.csv",
            "DeepJetC": "DeepJet_ctagSF_MiniAOD94X_2017_pTincl_v3_2_interp.root",
        },
        "LSF": {
            "ele_Trig TrigSF": "Ele32_L1DoubleEG_TrigSF_vhcc.histo.root",
            "ele_ID EGamma_SF2D": "ElectronIDSF_94X_MVA80WP.histo.root",
            "ele_Rereco EGamma_SF2D": "ElectronRecoSF_94X.histo.root",
            "mu_ID NUM_TightID_DEN_genTracks_pt_abseta": "RunBCDEF_MuIDSF.histo.root",
            "mu_ID_low NUM_TightID_DEN_genTracks_pt_abseta": "RunBCDEF_MuIDSF_lowpT.histo.root",
            "mu_Iso NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta": "RunBCDEF_MuISOSF.histo.root",
        },
    },
    "2017_UL": {
        "lumiMask": "Cert_294927-306462_13TeV_UL2017_Collisions17_MuonJSON.txt",
        "PU": None,
        "JME": "jec_compiled.pkl.gz",
        "BTV": {"deepCSV": "shape", "deepJet": "shape"},
        "LSF": {
            "ele_ID 2017 UL-Electron-ID-SF": "wp90iso",
            "ele_Reco 2017 UL-Electron-ID-SF": "RecoAbove20",
            "mu_Reco 2017_UL": "NUM_TrackerMuons_DEN_genTracks",
            "mu_HLT 2017_UL": "NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight",
            "mu_ID 2017_UL": "NUM_TightID_DEN_TrackerMuons",
            "mu_Iso 2017_UL": "NUM_TightRelIso_DEN_TightIDandIPCut",
            "mu_ID_low *": "Efficiency_muon_trackerMuon_Run2017_UL_ID.histo.json",
            "mu_Reco_low *": "Efficiency_muon_generalTracks_Run2017_UL_trackerMuon.histo.json",
        },
    },
    "Winter22Run3": {
        "lumiMask": "Cert_Collisions2022_355100_357900_Golden.json",
        "PU": "puweight_Winter22Run3.histo.root",
        "JME": "jec_compiled.pkl.gz",
        "jetveto": {"Run2022CD jetvetomap": "Winter22Run3_RunCD_v1.histo.root"},
    },
    "Summer22": {
        "lumiMask": "Cert_Collisions2022_355100_362760_Golden.json",
        "PU": "puwei_2022_preEE.histo.root",  # new PU files, based on preEE
        "JME": {
            "MC": "Summer22_22Sep2023_V2 Summer22_22Sep2023_JRV1",
            "Run2022C": "Summer22_22Sep2023_RunCD_V2",
            "Run2022D": "Summer22_22Sep2023_RunCD_V2",
        },  # update to latest JEC
        "LSF": {
            "mu_ID": "NUM_TightID_DEN_TrackerMuons",
            "mu_Iso": "NUM_TightPFIso_DEN_TightID",
            "ele_ID 2022Re-recoBCD Electron-ID-SF": "Tight",
            "ele_Reco 2022Re-recoBCD Electron-ID-SF": "",
        },
        "JPCalib": {
            "Run2022C": "calibeHistoWrite_Data2022C_NANO130X_v1.root",
            "Run2022D": "calibeHistoWrite_Data2022D_NANO130X_v1.root",
            "MC": "calibeHistoWrite_MC2022_NANO130X_v2.root",
        },
        "jetveto": {"Summer22_23Sep2023_RunCD_V1": "jetvetomap"},
        "muonSS": "",
        "electronSS": [
            "EGMScale_Compound_Ele_2022preEE",
            "EGMSmearAndSyst_ElePTsplit_2022preEE",
        ],
    },
    "Summer22EE": {
        "lumiMask": "Cert_Collisions2022_355100_362760_Golden.json",
        "PU": "puwei_2022_postEE.histo.root",  # new PU file, post EE
        "JME": {
            "MC": "Summer22EE_22Sep2023_V2 Summer22EE_22Sep2023_JRV1",
            "Run2022E": "Summer22EE_22Sep2023_RunE_V2",
            "Run2022F": "Summer22EE_22Sep2023_RunF_V2",
            "Run2022G": "Summer22EE_22Sep2023_RunG_V2",
        },
        "LSF": {
            "mu_ID": "NUM_TightID_DEN_TrackerMuons",
            "mu_Iso": "NUM_TightPFIso_DEN_TightID",
            "ele_ID 2022Re-recoE+PromptFG Electron-ID-SF": "Tight",
            "ele_Reco 2022Re-recoE+PromptFG Electron-ID-SF": "",
        },
        "jetveto": {"Summer22EE_23Sep2023_RunEFG_V1": "jetvetomap"},
        # use for BTA production, jet probablity
        "JPCalib": {
            "Run2022E": "calibeHistoWrite_Data2022F_NANO130X_v1.root",
            "Run2022F": "calibeHistoWrite_Data2022F_NANO130X_v1.root",
            "Run2022G": "calibeHistoWrite_Data2022G_NANO130X_v1.root",
            "MC": "calibeHistoWrite_MC2022EE_NANO130X_v1.root",
        },
        "muonSS": "",
        "electronSS": [
            "EGMScale_Compound_Ele_2022postEE",
            "EGMSmearAndSyst_ElePTsplit_2022postEE",
        ],
    },
    "Summer23": {
        "lumiMask": "Cert_Collisions2023_366442_370790_Golden.json",
        "PU": "puwei_Summer23.histo.root",
        "JME": {
            "name": "V1_AK4PFPuppi",
            "MC": [
                "Summer23Prompt23_V1_MC_L1FastJet_AK4PFPuppi",
                "Summer23Prompt23_V1_MC_L2Relative_AK4PFPuppi",
                "Summer23Prompt23_V1_MC_L2Residual_AK4PFPuppi",
                "Summer23Prompt23_V1_MC_L3Absolute_AK4PFPuppi",
                "Summer23Prompt23_V1_MC_UncertaintySources_AK4PFPuppi",
                "Summer23Prompt23_V1_MC_Uncertainty_AK4PFPuppi",
                "Summer23Prompt23_JRV1_MC_SF_AK4PFPuppi",
                "Summer23Prompt23_JRV1_MC_PtResolution_AK4PFPuppi",
            ],
            "dataCv123": [
                "Summer23Prompt23_RunCv123_V1_DATA_L1FastJet_AK4PFPuppi",
                "Summer23Prompt23_RunCv123_V1_DATA_L2Relative_AK4PFPuppi",
                "Summer23Prompt23_RunCv123_V1_DATA_L3Absolute_AK4PFPuppi",
                "Summer23Prompt23_RunCv123_V1_DATA_L2L3Residual_AK4PFPuppi",
            ],
            "dataCv4": [
                "Summer23Prompt23_RunCv4_V1_DATA_L1FastJet_AK4PFPuppi",
                "Summer23Prompt23_RunCv4_V1_DATA_L2Relative_AK4PFPuppi",
                "Summer23Prompt23_RunCv4_V1_DATA_L3Absolute_AK4PFPuppi",
                "Summer23Prompt23_RunCv4_V1_DATA_L2L3Residual_AK4PFPuppi",
            ],
        },
        "jetveto": {"Summer23Prompt23_RunC_V1": "jetvetomap"},
        "JPCalib": {
            "Run2023C-22Sep2023_v1": "calibeHistoWrite_Data2023C-22Sep2023_v1.root",
            "Run2023C-22Sep2023_v2": "calibeHistoWrite_Data2023C-22Sep2023_v2.root",
            "Run2023C-22Sep2023_v3": "calibeHistoWrite_Data2023C-22Sep2023_v3.root",
            "Run2023C-22Sep2023_v4": "calibeHistoWrite_Data2023C-22Sep2023_v4.root",
            "MC": "calibeHistoWrite_MC2023_Summer23.root",
        },
        "LSF": {
            "mu_ID": "NUM_TightID_DEN_TrackerMuons",
            "mu_Iso": "NUM_TightPFIso_DEN_TightID",
            "ele_ID 2023PromptC Electron-ID-SF": "Tight",
            "ele_Reco 2023PromptC Electron-ID-SF": "",
        },
        "muonSS": "",
        "electronSS": [
            "EGMScale_Compound_Ele_2023preBPIX",
            "EGMSmearAndSyst_ElePTsplit_2023preBPIX",
        ],
    },
    "Summer23BPix": {
        "lumiMask": "Cert_Collisions2023_366442_370790_Golden.json",
        "PU": "puwei_Summer23BPix.histo.root",
        "JME": {
            "MC": "Summer23BPixPrompt23_V1 Summer23BPixPrompt23_RunD_JRV1",
            "Run2023D": "Summer23BPixPrompt23_RunD_V1",
        },
        # "JME": "jec_compiled.pkl.gz",
        "LSF": {
            "mu_ID": "NUM_TightID_DEN_TrackerMuons",
            "mu_Iso": "NUM_TightPFIso_DEN_TightID",
            "ele_ID 2023PromptD Electron-ID-SF": "Tight",
            "ele_Reco 2023PromptD Electron-ID-SF": "",
        },
        # This is from Mikko https://indico.cern.ch/event/1315421/contributions/5532963/attachments/2697975/4683826/2023_08_16_jetvetomaps_v3.pdf
        "jetveto": {"Summer23BPixPrompt23_RunD_V1": "jetvetomap"},
        "JPCalib": {
            "Run2023D-22Sep2023_v1": "calibeHistoWrite_Data2023D-22Sep2023_v1.root",
            "Run2023D-22Sep2023_v2": "calibeHistoWrite_Data2023D-22Sep2023_v2.root",
            "MC": "calibeHistoWrite_MC2023_Summer23BPix.root",
        },
        "muonSS": "",
        "electronSS": [
            "EGMScale_Compound_Ele_2023postBPIX",
            "EGMSmearAndSyst_ElePTsplit_2023postBPIX",
        ],
    },
    "Summer24": {
        "lumiMask": "Cert_Collisions2024_378981_386951_Golden.json",
        "PU": "PU_weights_Summer24.histo.root",
        "JME": {
            # TODO: JER are a placeholder for now (July 2025)
            "MC": "Summer24Prompt24_V1 Summer23BPixPrompt23_RunD_JRV1",
            "Run2024C": "Summer24Prompt24_V1",
            "Run2024D": "Summer24Prompt24_V1",
            "Run2024E": "Summer24Prompt24_V1",
            "Run2024F": "Summer24Prompt24_V1",
            "Run2024G": "Summer24Prompt24_V1",
            "Run2024H": "Summer24Prompt24_V1",
            "Run2024I": "Summer24Prompt24_V1",
        },
        "jetveto": {"Summer24Prompt24_RunBCDEFGHI_V1": "jetvetomap"},
        "LSF": {
            "mu_ID": "NUM_TightID_DEN_TrackerMuons",
            "mu_Iso": "NUM_TightPFIso_DEN_TightID",
            "ele_Reco 2024 Electron-ID-SF": "",
            "ele_ID 2024 Electron-ID-SF": "wp80iso",
        },
        "muonSS": "",
        "electronSS": ["EGMScale_Compound_Ele_2024", "EGMSmearAndSyst_ElePTsplit_2024"],
    },
    "prompt_dataMC": {"lumiMask": "$PROMPT_DATAMC"},
}
