correction_config = {
    "Rereco17_94X": {
        "DC": "Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt",
        "LUM": "94XPUwei_corrections.pkl.gz",
        "JME": "jec_compiled.pkl.gz",
        "BTV": {
            "DeepCSVB": "DeepCSV_94XSF_V5_B_F.csv",
            "DeepCSVC": "DeepCSV_ctagSF_MiniAOD94X_2017_pTincl_v3_2_interp.root",
            "DeepJetB": "DeepFlavour_94XSF_V4_B_F.csv",
            "DeepJetC": "DeepJet_ctagSF_MiniAOD94X_2017_pTincl_v3_2_interp.root",
        },
        "EGM": {
            "ele_Trig TrigSF": "Ele32_L1DoubleEG_TrigSF_vhcc.histo.root",
            "ele_ID EGamma_SF2D": "ElectronIDSF_94X_MVA80WP.histo.root",
            "ele_Rereco EGamma_SF2D": "ElectronRecoSF_94X.histo.root",
        },
        "MUO": {
            "mu_ID NUM_TightID_DEN_genTracks_pt_abseta": "RunBCDEF_MuIDSF.histo.root",
            "mu_ID_low NUM_TightID_DEN_genTracks_pt_abseta": "RunBCDEF_MuIDSF_lowpT.histo.root",
            "mu_Iso NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta": "RunBCDEF_MuISOSF.histo.root",
        },
    },
    "2017-UL": {
        "DC": "Cert_294927-306462_13TeV_UL2017_Collisions17_MuonJSON.txt",
        "LUM": None,
        "JME": {
            "MC": "Summer20UL17_V1 Summer19UL17_JRV3",
            "Run2017B": "Summer20UL17_RunB_V1",
            "Run2017C": "Summer20UL17_RunC_V1",
            "Run2017D": "Summer20UL17_RunD_V1",
            "Run2017E": "Summer20UL17_RunE_V1",
            "Run2017F": "Summer20UL17_RunF_V1",
        },
        "JME_path": "src/BTVNanoCommissioning/data/JME/2017-UL/jet_jerc.json.gz",
        "BTV": {"deepCSV": "shape", "deepJet": "shape"},
        "EGM": {
            "ele_ID 2017 UL-Electron-ID-SF": "wp90iso",
            "ele_Reco 2017 UL-Electron-ID-SF": "RecoAbove20",
        },
        "MUO": {
            "mu_Reco 2017_UL": "NUM_TrackerMuons_DEN_genTracks",
            "mu_HLT 2017_UL": "NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight",
            "mu_ID 2017_UL": "NUM_TightID_DEN_TrackerMuons",
            "mu_Iso 2017_UL": "NUM_TightRelIso_DEN_TightIDandIPCut",
            "mu_ID_low *": "Efficiency_muon_trackerMuon_Run2017_UL_ID.histo.json",
            "mu_Reco_low *": "Efficiency_muon_generalTracks_Run2017_UL_trackerMuon.histo.json",
        },
    },
    "2016preVFP-UL": {
        "DC": "Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt",
        "LUM": None,
        "JME": {
            "MC": "Summer20UL16APV_V1 Summer20UL16APV_JRV3",
            "Run2016B": "Summer20UL16APV_RunBCD_V1",
            "Run2016C": "Summer20UL16APV_RunBCD_V1",
            "Run2016D": "Summer20UL16APV_RunBCD_V1",
            "Run2016E": "Summer20UL16APV_RunEF_V1",
            "Run2016F": "Summer20UL16APV_RunEF_V1",
        },
        "JME_path": "src/BTVNanoCommissioning/data/JME/2016preVFP-UL/jet_jerc.json.gz",
        "BTV": {"deepCSV": "shape", "deepJet": "shape"},
        "EGM": {
            "ele_ID 2016preVFP UL-Electron-ID-SF": "wp90iso",
            "ele_Reco 2016preVFP UL-Electron-ID-SF": "RecoAbove20",
        },
        "MUO": {
            "mu_Reco 2016preVFP_UL": "NUM_TrackerMuons_DEN_genTracks",
            "mu_ID 2016preVFP_UL": "NUM_TightID_DEN_TrackerMuons",
            "mu_Iso 2016preVFP_UL": "NUM_TightRelIso_DEN_TightIDandIPCut",
        },
    },
    "2016postVFP-UL": {
        "DC": "Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt",
        "LUM": None,
        "JME": {
            "MC": "Summer20UL16_V1 Summer20UL16_JRV3",
            "Run2016F": "Summer20UL16_RunFGH_V1",
            "Run2016G": "Summer20UL16_RunFGH_V1",
            "Run2016H": "Summer20UL16_RunFGH_V1",
        },
        "JME_path": "src/BTVNanoCommissioning/data/JME/2016postVFP-UL/jet_jerc.json.gz",
        "BTV": {"deepCSV": "shape", "deepJet": "shape"},
        "EGM": {
            "ele_ID 2016postVFP UL-Electron-ID-SF": "wp90iso",
            "ele_Reco 2016postVFP UL-Electron-ID-SF": "RecoAbove20",
        },
        "MUO": {
            "mu_Reco 2016postVFP_UL": "NUM_TrackerMuons_DEN_genTracks",
            "mu_ID 2016postVFP_UL": "NUM_TightID_DEN_TrackerMuons",
            "mu_Iso 2016postVFP_UL": "NUM_TightRelIso_DEN_TightIDandIPCut",
        },
    },
    "2018-UL": {
        "DC": "Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt",
        "LUM": None,
        "JME": {
            "MC": "Summer20UL18_V1 Summer19UL18_JRV2",
            "Run2018A": "Summer20UL18_RunA_V1",
            "Run2018B": "Summer20UL18_RunB_V1",
            "Run2018C": "Summer20UL18_RunC_V1",
            "Run2018D": "Summer20UL18_RunD_V1",
        },
        "JME_path": "src/BTVNanoCommissioning/data/JME/2018-UL/jet_jerc.json.gz",
        "BTV": {"deepCSV": "shape", "deepJet": "shape"},
        "EGM": {
            "ele_ID 2018 UL-Electron-ID-SF": "wp90iso",
            "ele_Reco 2018 UL-Electron-ID-SF": "RecoAbove20",
        },
        "MUO": {
            "mu_Reco 2018_UL": "NUM_TrackerMuons_DEN_genTracks",
            "mu_ID 2018_UL": "NUM_TightID_DEN_TrackerMuons",
            "mu_Iso 2018_UL": "NUM_TightRelIso_DEN_TightIDandIPCut",
        },
    },
    "Winter22Run3": {
        "DC": "Cert_Collisions2022_355100_357900_Golden.json",
        "LUM": "puweight_Winter22Run3.histo.root",
        "JME": "jec_compiled.pkl.gz",
        "jetveto": {"Run2022CD jetvetomap": "Winter22Run3_RunCD_v1.histo.root"},
    },
    "Summer22": {
        "DC": "Cert_Collisions2022_355100_362760_Golden.json",
        "LUM": "puwei_2022_preEE.histo.root",  # new PU files, based on preEE
        "JME": {
            "MC": "Summer22_22Sep2023_V3 Summer22_22Sep2023_JRV1",
            "Run2022C": "Summer22_22Sep2023_RunCD_V3",
            "Run2022D": "Summer22_22Sep2023_RunCD_V3",
        },  # update to latest JEC
        "MUO": {
            "mu_ID": "NUM_TightID_DEN_TrackerMuons",
            "mu_Iso": "NUM_TightPFIso_DEN_TightID",
        },
        "EGM": {
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
            "Scale",
            "SmearAndSyst",
        ],
    },
    "Summer22EE": {
        "DC": "Cert_Collisions2022_355100_362760_Golden.json",
        "LUM": "puwei_2022_postEE.histo.root",  # new PU file, post EE
        "JME": {
            "MC": "Summer22EE_22Sep2023_V3 Summer22EE_22Sep2023_JRV1",
            "Run2022E": "Summer22EE_22Sep2023_RunE_V3",
            "Run2022F": "Summer22EE_22Sep2023_RunF_V3",
            "Run2022G": "Summer22EE_22Sep2023_RunG_V3",
        },
        "MUO": {
            "mu_ID": "NUM_TightID_DEN_TrackerMuons",
            "mu_Iso": "NUM_TightPFIso_DEN_TightID",
        },
        "EGM": {
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
            "Scale",
            "SmearAndSyst",
        ],
    },
    "Summer23": {
        "DC": "Cert_Collisions2023_366442_370790_Golden.json",
        "LUM": "puwei_Summer23.histo.root",
        "JME": {
            "MC": "Summer23Prompt23_V2 Summer23Prompt23_RunCv1234_JRV1",
            "Run2023C": "Summer23Prompt23_V2",
        },
        "jetveto": {"Summer23Prompt23_RunC_V1": "jetvetomap"},
        "JPCalib": {
            "Run2023C-22Sep2023_v1": "calibeHistoWrite_Data2023C-22Sep2023_v1.root",
            "Run2023C-22Sep2023_v2": "calibeHistoWrite_Data2023C-22Sep2023_v2.root",
            "Run2023C-22Sep2023_v3": "calibeHistoWrite_Data2023C-22Sep2023_v3.root",
            "Run2023C-22Sep2023_v4": "calibeHistoWrite_Data2023C-22Sep2023_v4.root",
            "MC": "calibeHistoWrite_MC2023_Summer23.root",
        },
        "MUO": {
            "mu_ID": "NUM_TightID_DEN_TrackerMuons",
            "mu_Iso": "NUM_TightPFIso_DEN_TightID",
        },
        "EGM": {
            "ele_ID 2023PromptC Electron-ID-SF": "Tight",
            "ele_Reco 2023PromptC Electron-ID-SF": "",
        },
        "muonSS": "",
        "electronSS": [
            "Scale",
            "SmearAndSyst",
        ],
    },
    "Summer23BPix": {
        "DC": "Cert_Collisions2023_366442_370790_Golden.json",
        "LUM": "puwei_Summer23BPix.histo.root",
        "JME": {
            "MC": "Summer23BPixPrompt23_V3 Summer23BPixPrompt23_RunD_JRV1",
            "Run2023D": "Summer23BPixPrompt23_V3",
        },
        "MUO": {
            "mu_ID": "NUM_TightID_DEN_TrackerMuons",
            "mu_Iso": "NUM_TightPFIso_DEN_TightID",
        },
        "EGM": {
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
            "Scale",
            "SmearAndSyst",
        ],
    },
    "Summer24": {
        "DC": "Cert_Collisions2024_378981_386951_Golden.json",
        "LUM": "puWeights_BCDEFGHI.json.gz",
        "JME": {
            # TODO: JER are a placeholder for now (July 2025)
            "MC": "Summer24Prompt24_V2 Summer23BPixPrompt23_RunD_JRV1",
            "Run2024C": "Summer24Prompt24_V2",
            "Run2024D": "Summer24Prompt24_V2",
            "Run2024E": "Summer24Prompt24_V2",
            "Run2024F": "Summer24Prompt24_V2",
            "Run2024G": "Summer24Prompt24_V2",
            "Run2024H": "Summer24Prompt24_V2",
            "Run2024I": "Summer24Prompt24_V2",
        },
        "jetveto": {"Summer24Prompt24_RunBCDEFGHI_V1": "jetvetomap"},
        "MUO": {
            "mu_ID": "NUM_TightID_DEN_TrackerMuons",
            "mu_Iso": "NUM_TightPFIso_DEN_TightID",
        },
        "EGM": {
            "ele_Reco 2024 Electron-ID-SF": "",
            "ele_ID 2024 Electron-ID-SF": "wp80iso",
        },
        "muonSS": "",
        "electronSS": [
            "Scale",
            "SmearAndSyst",
        ],
    },
    "Winter25": {
        "DC": "Cert_Collisions2025_391658_398860_Golden.json",
        "LUM": "puWeights2025.json.gz",
        "JME": {
            # MC: use Summer24 MC truth JECs from the Summer24 CVMFS era
            # JER: placeholder from Summer23BPix until dedicated 2025 JER is derived.
            "MC": "Summer24Prompt24_V2 Summer23BPixPrompt23_RunD_JRV1",
            "Run2025C": "Winter25Prompt25_V3",
            "Run2025D": "Winter25Prompt25_V3",
            "Run2025E": "Winter25Prompt25_V3",
            "Run2025F": "Winter25Prompt25_V3",
            "Run2025G": "Winter25Prompt25_V3",
        },
        "jetveto": {"Winter25Prompt25_RunCDEFG_V1": "jetvetomap"},
        "MUO": {
            "mu_ID": "NUM_TightID_DEN_TrackerMuons",
            "mu_Iso": "NUM_TightPFIso_DEN_TightID",
        },
        "EGM": {
            "ele_Reco 2024 Electron-ID-SF": "",
            "ele_ID 2024 Electron-ID-SF": "wp80iso",
        },
        # Muon scale & smearing: reuse 2024 from Run3-24CDE...Summer24
        "muonSS": "",
        # Electron scale & smearing: use 2025 SaS from Run3-25Prompt-Summer24
        "electronSS": [
            "Scale",
            "SmearAndSyst",
        ],
        # Per-POG CVMFS path overrides.
        # JME_MC: MC truth JECs must come from the Summer24 era (L2Relative differs).
        # MUO: 2025 SFs available (muon_Z.json.gz under muo_SF25_Z_and_highpt/,
        #       copied to local data/ fallback since latest/ is still empty on CVMFS).
        # muonSS: no 2025 muon scale/smearing yet, reuse Summer24.
        # EGM: no 2025 electron ID SFs yet (electron.json.gz), reuse Summer24.
        # electronSS: 2025 SaS available.
        "cvmfs_override": {
            "JME_MC": "Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15",
            "MUO": "Run3-25Prompt-Summer24-NanoAODv15",
            "muonSS": "Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15",
            "EGM": "Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15",
            "electronSS": "Run3-25Prompt-Summer24-NanoAODv15",
        },
    },
    "prompt_dataMC": {"DC": "$PROMPT_DATAMC"},
}
