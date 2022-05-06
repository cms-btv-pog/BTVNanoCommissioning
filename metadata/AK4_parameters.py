correction_config= {
    "Rereco17_94X": {
        "PU":"data/PU/Rereco17_94X/94XPUwei_corrections.pkl.gz",
        "JME" :"data/JME/Rereco17_94X/jec_compiled.pkl.gz",
        "BTV" : {
            "DeepCSVB":"data/BTV/Rereco17_94X/DeepCSV_94XSF_V5_B_F.csv",
            "DeepCSVC":"data/BTV/Rereco17_94X/DeepCSV_ctagSF_MiniAOD94X_2017_pTincl_v3_2_interp.root",
            "DeepJetB":"data/BTV/Rereco17_94X/DeepFlavour_94XSF_V4_B_F.csv",
            "DeepJetC":"data/BTV/Rereco17_94X/DeepJet_ctagSF_MiniAOD94X_2017_pTincl_v3_2_interp.root"
        },
        "LSF":{
            "ele_Trig":"ele_Trig TrigSF data/LSF/Rereco17_94X/Ele32_L1DoubleEG_TrigSF_vhcc.histo.root",
            "ele_ID":"ele_ID EGamma_SF2D data/LSF/Rereco17_94X/ElectronIDSF_94X_MVA80WP.histo.root",
            "ele_Rereco":"ele_Rereco EGamma_SF2D data/LSF/Rereco17_94X/ElectronRecoSF_94X.histo.root",
            "mu_ID":"mu_ID NUM_TightID_DEN_genTracks_pt_abseta data/LSF/Rereco17_94X/RunBCDEF_SF_ID.histo.root",
            "mu_ID_low":"mu_ID_low NUM_TightID_DEN_genTracks_pt_abseta data/LSF/Rereco17_94X/RunBCDEF_SF_MuID_lowpT.histo.root",
            "mu_Iso":"mu_Iso NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta data/LSF/Rereco17_94X/RunBCDEF_SF_ISO.histo.root"
        }
    }
}