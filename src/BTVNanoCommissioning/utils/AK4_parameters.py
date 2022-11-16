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
    "UL17": {
        "lumiMask": "Cert_294927-306462_13TeV_UL2017_Collisions17_MuonJSON.txt",
        "PU": "puweight_UL17.histo.root",
        "JME": "mc_compile_jec.pkl.gz",
        "BTV": {"DeepJetC": "DeepJet_ctagSF_Summer20UL17_interp.root"},
        "LSF": {
            "ele_Rereco_above20 EGamma_SF2D": "egammaEffi_ptAbove20.txt_EGM2D_UL2017.histo.root",
            "ele_Rereco_below20 EGamma_SF2D": "egammaEffi_ptBelow20.txt_EGM2D_UL2017.histo.root",
            "ele_ID EGamma_SF2D": "egammaEffi.txt_EGM2D_MVA90iso_UL17.histo.root",
            "mu_ID NUM_TightID_DEN_TrackerMuons_abseta_pt": "Efficiencies_muon_generalTracks_Z_Run2017_UL_ID.histo.root",
            "mu_Iso NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt": "Efficiencies_muon_generalTracks_Z_Run2017_UL_ISO.histo.root",
        },
    },
    "Winter22Run3": {
        "lumiMask": "Cert_Collisions2022_355100_357900_Golden.json",
        "PU": "puweight_Winter22Run3.histo.root",
        "JME": "jec_compiled.pkl.gz",
    },
}
