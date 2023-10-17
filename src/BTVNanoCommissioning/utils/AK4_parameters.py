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
    },
    "Summer22Run3": {
        "lumiMask": "Cert_Collisions2022_355100_362760_Golden.json",
        "PU": "puweight_Winter22Run3.histo.root",  # use same 69.2mb for Summer22
        "JME": "winter_jec_compiled.pkl.gz",  # not this is from Winter22Run3 since this is not yet finished
        "LSF": {
            "mu_json": "ScaleFactors_Muon_trackerMuons_Z_2022_Prompt_ID_ISO_schemaV2.json",
            "mu_ID": "NUM_TightID_DEN_TrackerMuons",
            "mu_Iso": "NUM_TightPFIso_DEN_TightID",
        },
        "JPCalib": {
            "Run2022C": "calibeHistoWrite_Data2022F_NANO130X_v1.root",  # to be updated
            "Run2022D": "calibeHistoWrite_Data2022F_NANO130X_v1.root",  # to be updated
            "MC": "calibeHistoWrite_MC2022EE_NANO130X_v1.root",  # to be updated
        },
    },
    "Summer22EERun3": {
        "lumiMask": "Cert_Collisions2022_355100_362760_Golden.json",
        "PU": "puweight_Summer22EERun3.histo.root",  # 80mb
        "JME": "jec_compiled.pkl.gz",
        "LSF": {
            "ele_json": "electron.json.gz",
            "mu_json": "ScaleFactors_Muon_trackerMuons_Z_2022EE_Prompt_ID_ISO_schemaV2.json",
            "mu_ID": "NUM_TightID_DEN_TrackerMuons",
            "mu_Iso": "NUM_TightPFIso_DEN_TightID",
            "ele_ID 2022FG 2022FG-Electron-ID-SF": "Tight",
            "ele_Reco_low 2022FG 2022FG-Electron-ID-SF": "RecoBelow20",
            "ele_Reco_med 2022FG 2022FG-Electron-ID-SF": "Reco20to75",
            "ele_Reco_high 2022FG 2022FG-Electron-ID-SF": "RecoAbove75",
        },
        "JPCalib": {
            "Run2022E": "calibeHistoWrite_Data2022F_NANO130X_v1.root",  # to be updated
            "Run2022F": "calibeHistoWrite_Data2022F_NANO130X_v1.root",
            "Run2022G": "calibeHistoWrite_Data2022G_NANO130X_v1.root",
            "MC": "calibeHistoWrite_MC2022EE_NANO130X_v1.root",
        },
    },
}
