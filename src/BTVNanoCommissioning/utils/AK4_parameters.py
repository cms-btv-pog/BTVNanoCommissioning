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
            "ele_ID 2017": "wp90iso",
            "ele_Reco 2017": "RecoAbove20",
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
}
