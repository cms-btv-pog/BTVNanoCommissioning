from BTVNanoCommissioning.helpers.definitions import definitions, SV_definitions
import hist as Hist


def histogrammer(events, workflow):
    _hist_dict = {}
    ## Common variables
    flav_axis = Hist.axis.IntCategory([0, 1, 4, 5, 6], name="flav", label="Genflavour")
    syst_axis = Hist.axis.StrCategory([], name="syst", growth=True)
    pt_axis = Hist.axis.Regular(60, 0, 300, name="pt", label=" $p_{T}$ [GeV]")
    jpt_axis = Hist.axis.Regular(300, 0, 3000, name="pt", label=" $p_{T}$ [GeV]")
    softlpt_axis = Hist.axis.Regular(25, 0, 25, name="pt", label=" $p_{T}$ [GeV]")
    mass_axis = Hist.axis.Regular(50, 0, 300, name="mass", label=" $p_{T}$ [GeV]")
    eta_axis = Hist.axis.Regular(25, -2.5, 2.5, name="eta", label=" $\eta$")
    phi_axis = Hist.axis.Regular(30, -3, 3, name="phi", label="$\phi$")
    mt_axis = Hist.axis.Regular(30, 0, 300, name="mt", label=" $m_{T}$ [GeV]")
    iso_axis = Hist.axis.Regular(30, 0, 0.05, name="pfRelIso03_all", label="Rel. Iso")
    softliso_axis = Hist.axis.Regular(
        20, 0.2, 6.2, name="pfRelIso03_all", label="Rel. Iso"
    )
    npvs_axis = Hist.axis.Integer(0, 100, name="npv", label="N PVs")
    dr_axis = Hist.axis.Regular(20, 0, 8, name="dr", label="$\Delta$R")
    dr_s_axis = Hist.axis.Regular(20, 0, 0.5, name="dr", label="$\Delta$R")
    dr_SV_axis = Hist.axis.Regular(20, 0, 1.0, name="dr", label="$\Delta$R")
    dxy_axis = Hist.axis.Regular(40, -0.05, 0.05, name="dxy", label="d_{xy}")
    dz_axis = Hist.axis.Regular(40, -0.01, 0.01, name="dz", label="d_{z}")
    qcddxy_axis = Hist.axis.Regular(40, -0.002, 0.002, name="dxy", label="d_{xy}")
    sip3d_axis = Hist.axis.Regular(20, 0, 0.2, name="sip3d", label="SIP 3D")
    ptratio_axis = Hist.axis.Regular(50, 0, 1, name="ratio", label="ratio")
    n_axis = Hist.axis.Integer(0, 10, name="n", label="N obj")
    osss_axis = Hist.axis.IntCategory([1, -1], name="osss", label="OS(+)/SS(-)")
    _hist_dict["npvs"] = Hist.Hist(syst_axis, npvs_axis, Hist.storage.Weight())
    _hist_dict["pu"] = Hist.Hist(syst_axis, npvs_axis, Hist.storage.Weight())
    ### Workflow specific
    if "example" == workflow:
        obj_list = [
            "jet",
            "mu",
        ]  # store basic 4-vector, pt,eta, phi, mass for the object
        _hist_dict[f"dr_mujet"] = Hist.Hist(
            syst_axis, flav_axis, dr_axis, Hist.storage.Weight()
        )  # create cutstomize histogram
    elif "QCD" == workflow:
        obj_list = ["jet0"]
        _hist_dict["dr_SVjet0"] = Hist.Hist(
            syst_axis, flav_axis, dr_SV_axis, Hist.storage.Weight()
        )
        _hist_dict["nJetSVs"] = Hist.Hist(syst_axis, n_axis, Hist.storage.Weight())

    elif "validation" == workflow:
        obj_list = ["jet0", "jet1"]
        _hist_dict["bjet_WP_pt"] = Hist.Hist(
            Hist.axis.StrCategory([], name="WP", growth=True),
            Hist.axis.StrCategory(
                ["DeepFlav", "PNet", "RobustParTAK4"], name="tagger", growth=True
            ),
            pt_axis,
            Hist.storage.Weight(),
        )
        _hist_dict["bjet_WP_eta"] = Hist.Hist(
            Hist.axis.StrCategory([], name="WP", growth=True),
            Hist.axis.StrCategory(
                ["DeepFlav", "PNet", "RobustParTAK4"], name="tagger", growth=True
            ),
            eta_axis,
            Hist.storage.Weight(),
        )
        _hist_dict["bjet_WP_phi"] = Hist.Hist(
            Hist.axis.StrCategory([], name="WP", growth=True),
            Hist.axis.StrCategory(
                ["DeepFlav", "PNet", "RobustParTAK4"], name="tagger", growth=True
            ),
            phi_axis,
            Hist.storage.Weight(),
        )
        _hist_dict["bjet_WP_discr"] = Hist.Hist(
            Hist.axis.StrCategory([], name="WP", growth=True),
            Hist.axis.StrCategory(
                ["DeepFlav", "PNet", "RobustParTAK4"], name="tagger", growth=True
            ),
            Hist.axis.Regular(25, 0, 1, name="B"),
            Hist.storage.Weight(),
        )
        _hist_dict["cjet_WP_pt"] = Hist.Hist(
            Hist.axis.StrCategory([], name="WP", growth=True),
            Hist.axis.StrCategory(
                ["DeepFlav", "PNet", "RobustParTAK4"], name="tagger", growth=True
            ),
            pt_axis,
            Hist.storage.Weight(),
        )
        _hist_dict["cjet_WP_eta"] = Hist.Hist(
            Hist.axis.StrCategory([], name="WP", growth=True),
            Hist.axis.StrCategory(
                ["DeepFlav", "PNet", "RobustParTAK4"], name="tagger", growth=True
            ),
            eta_axis,
            Hist.storage.Weight(),
        )
        _hist_dict["cjet_WP_phi"] = Hist.Hist(
            Hist.axis.StrCategory([], name="WP", growth=True),
            Hist.axis.StrCategory(
                ["DeepFlav", "PNet", "RobustParTAK4"], name="tagger", growth=True
            ),
            phi_axis,
            Hist.storage.Weight(),
        )
        _hist_dict["cjet_WP_discr"] = Hist.Hist(
            Hist.axis.StrCategory([], name="WP", growth=True),
            Hist.axis.StrCategory(
                ["DeepFlav", "PNet", "RobustParTAK4"], name="tagger", growth=True
            ),
            Hist.axis.Regular(25, 0, 1, name="CvL"),
            Hist.axis.Regular(25, 0, 1, name="CvB"),
            Hist.storage.Weight(),
        )

    elif "ttdilep_sf" == workflow:
        obj_list = ["mu", "ele"]
        for i in range(2):
            obj_list.append(f"jet{i}")
            _hist_dict[f"dr_mujet{i}"] = Hist.Hist(
                syst_axis, flav_axis, dr_axis, Hist.storage.Weight()
            )
        for i in ["mu", "ele"]:
            if i == "mu":
                _hist_dict[f"{i}_pfRelIso04_all"] = Hist.Hist(
                    syst_axis, iso_axis, Hist.storage.Weight()
                )
            else:
                _hist_dict[f"{i}_pfRelIso03_all"] = Hist.Hist(
                    syst_axis, iso_axis, Hist.storage.Weight()
                )
            _hist_dict[f"{i}_dxy"] = Hist.Hist(
                syst_axis, dxy_axis, Hist.storage.Weight()
            )
            _hist_dict[f"{i}_dz"] = Hist.Hist(syst_axis, dz_axis, Hist.storage.Weight())
    elif "ttsemilep_sf" == workflow:
        obj_list = ["mu", "MET"]
        for i in range(4):
            obj_list.append(f"jet{i}")
            _hist_dict[f"dr_mujet{i}"] = Hist.Hist(
                syst_axis, flav_axis, dr_axis, Hist.storage.Weight()
            )
        for i in ["mu"]:
            _hist_dict[f"{i}_pfRelIso04_all"] = Hist.Hist(
                syst_axis, iso_axis, Hist.storage.Weight()
            )
            _hist_dict[f"{i}_dxy"] = Hist.Hist(
                syst_axis, dxy_axis, Hist.storage.Weight()
            )
            _hist_dict[f"{i}_dz"] = Hist.Hist(syst_axis, dz_axis, Hist.storage.Weight())
    elif "ctag_ttdilep_sf" in workflow:
        obj_list = ["hl", "sl", "soft_l", "MET", "z", "lmujet"]
        _hist_dict["z_mass"] = Hist.Hist(
            syst_axis,
            Hist.axis.Regular(
                50, 50, 100, name="mass", label=" $m_{\\ell\\ell}$ [GeV]"
            ),
            Hist.storage.Weight(),
        )
        # delta R between soft muon and mu-jet
        _hist_dict["dr_lmujetsmu"] = Hist.Hist(
            syst_axis, flav_axis, dr_s_axis, Hist.storage.Weight()
        )
        # delta R between hard muon and mu-jet
        _hist_dict["dr_lmujethmu"] = Hist.Hist(
            syst_axis, flav_axis, dr_axis, Hist.storage.Weight()
        )
        # delta R between soft muon and hard muon
        _hist_dict["dr_lmusmu"] = Hist.Hist(syst_axis, dr_axis, Hist.storage.Weight())
        for i in ["hl", "sl", "soft_l"]:
            if i == "soft_l":
                _hist_dict[f"soft_l_pfRelIso04_all"] = Hist.Hist(
                    syst_axis, flav_axis, softliso_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{i}_dxy"] = Hist.Hist(
                    syst_axis, flav_axis, dxy_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{i}_dz"] = Hist.Hist(
                    syst_axis, flav_axis, dz_axis, Hist.storage.Weight()
                )
            else:
                if "m" in workflow:
                    _hist_dict[f"{i}_pfRelIso04_all"] = Hist.Hist(
                        syst_axis, iso_axis, Hist.storage.Weight()
                    )
                _hist_dict[f"{i}_dxy"] = Hist.Hist(
                    syst_axis, qcddxy_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{i}_dz"] = Hist.Hist(
                    syst_axis, dz_axis, Hist.storage.Weight()
                )
            # lepton / jet pT ratio
            _hist_dict[f"{i}_ptratio"] = Hist.Hist(
                syst_axis, flav_axis, ptratio_axis, Hist.storage.Weight()
            )
    elif "ctag_ttsemilep_sf" in workflow:
        obj_list = ["hl", "soft_l", "MET", "z", "w", "mujet"]
        _hist_dict["z_mass"] = Hist.Hist(
            syst_axis,
            Hist.axis.Regular(50, 50, 100, name="mass", label="$m_{\\ell\\ell}$ [GeV]"),
            Hist.storage.Weight(),
        )
        _hist_dict["w_mass"] = Hist.Hist(
            syst_axis,
            Hist.axis.Regular(50, 50, 100, name="mass", label="$m_{\\ell\\nu}$ [GeV]"),
            Hist.storage.Weight(),
        )
        # delta R between soft muon and mu-jet
        _hist_dict["dr_lmujetsmu"] = Hist.Hist(
            syst_axis, flav_axis, dr_s_axis, Hist.storage.Weight()
        )
        # delta R between hard muon and mu-jet
        _hist_dict["dr_lmujethmu"] = Hist.Hist(
            syst_axis, flav_axis, dr_axis, Hist.storage.Weight()
        )
        # delta R between hard muon and soft-muon
        _hist_dict["dr_lmusmu"] = Hist.Hist(syst_axis, dr_axis, Hist.storage.Weight())
        for i in ["hl", "soft_l"]:
            if i == "soft_l":
                _hist_dict[f"soft_l_pfRelIso04_all"] = Hist.Hist(
                    syst_axis, flav_axis, softliso_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{i}_dxy"] = Hist.Hist(
                    syst_axis, flav_axis, dxy_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{i}_dz"] = Hist.Hist(
                    syst_axis, flav_axis, dz_axis, Hist.storage.Weight()
                )
            else:
                _hist_dict[f"{i}_pfRelIso04_all"] = Hist.Hist(
                    syst_axis, iso_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{i}_dxy"] = Hist.Hist(
                    syst_axis, dxy_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{i}_dz"] = Hist.Hist(
                    syst_axis, dz_axis, Hist.storage.Weight()
                )
            _hist_dict[f"{i}_ptratio"] = Hist.Hist(
                syst_axis, flav_axis, ptratio_axis, Hist.storage.Weight()
            )
    elif "Wc_sf" in workflow:
        obj_list = ["hl", "soft_l", "MET", "z", "w", "mujet"]
        _hist_dict["SV_charge"] = Hist.Hist(
            syst_axis,
            osss_axis,
            Hist.axis.Regular(20, -10, 10, name="charge", label="SV charge"),
            Hist.storage.Weight(),
        )
        _hist_dict["z_mass"] = Hist.Hist(
            syst_axis,
            osss_axis,
            Hist.axis.Regular(50, 50, 100, name="mass", label="$m_{\\ell\\ell}$ [GeV]"),
            Hist.storage.Weight(),
        )
        _hist_dict["w_mass"] = Hist.Hist(
            syst_axis,
            osss_axis,
            Hist.axis.Regular(50, 50, 100, name="mass", label="$m_{\\ell\\nu}$ [GeV]"),
            Hist.storage.Weight(),
        )
        # delta R between soft muon and mu-jet
        _hist_dict["dr_lmujetsmu"] = Hist.Hist(
            syst_axis, flav_axis, osss_axis, dr_s_axis, Hist.storage.Weight()
        )
        # delta R between hard muon and mu-jet
        _hist_dict["dr_lmujethmu"] = Hist.Hist(
            syst_axis, flav_axis, osss_axis, dr_axis, Hist.storage.Weight()
        )
        # delta R between hard muon and soft-muon
        _hist_dict["dr_lmusmu"] = Hist.Hist(
            syst_axis, osss_axis, dr_axis, Hist.storage.Weight()
        )
        for i in ["hl", "soft_l"]:
            if i == "soft_l":
                _hist_dict[f"soft_l_pfRelIso04_all"] = Hist.Hist(
                    syst_axis,
                    flav_axis,
                    osss_axis,
                    softliso_axis,
                    Hist.storage.Weight(),
                )
                _hist_dict[f"{i}_dxy"] = Hist.Hist(
                    syst_axis, flav_axis, osss_axis, dxy_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{i}_dz"] = Hist.Hist(
                    syst_axis, flav_axis, osss_axis, dz_axis, Hist.storage.Weight()
                )
            else:
                _hist_dict[f"{i}_pfRelIso04_all"] = Hist.Hist(
                    syst_axis, osss_axis, iso_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{i}_dxy"] = Hist.Hist(
                    syst_axis, osss_axis, qcddxy_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{i}_dz"] = Hist.Hist(
                    syst_axis, osss_axis, dz_axis, Hist.storage.Weight()
                )
            _hist_dict[f"{i}_ptratio"] = Hist.Hist(
                syst_axis, flav_axis, osss_axis, ptratio_axis, Hist.storage.Weight()
            )
    elif "DY_sf" in workflow:
        obj_list = ["posl", "negl", "z", "jet"]
        _hist_dict["z_mass"] = Hist.Hist(
            syst_axis,
            Hist.axis.Regular(50, 50, 100, name="mass", label="$m_{\\ell\\ell}$ [GeV]"),
            Hist.storage.Weight(),
        )
        _hist_dict["dr_mumu"] = Hist.Hist(syst_axis, dr_axis, Hist.storage.Weight())
        for i in ["posl", "negl"]:
            if "m" in workflow:
                _hist_dict[f"{i}_pfRelIso04_all"] = Hist.Hist(
                    syst_axis, iso_axis, Hist.storage.Weight()
                )
            _hist_dict[f"{i}_dxy"] = Hist.Hist(
                syst_axis, dxy_axis, Hist.storage.Weight()
            )
            _hist_dict[f"{i}_dz"] = Hist.Hist(syst_axis, dz_axis, Hist.storage.Weight())
        for i in range(2):
            obj_list.append(f"jet{i}")
            _hist_dict[f"dr_mujet{i}"] = Hist.Hist(
                syst_axis, flav_axis, dr_axis, Hist.storage.Weight()
            )
    elif "QCDmuen" in workflow:
        obj_list = ["mu", "MET"]
        for i in range(2):
            obj_list.append(f"jet{i}")
            _hist_dict[f"dr_mujet{i}"] = Hist.Hist(
                syst_axis, flav_axis, dr_axis, Hist.storage.Weight()
            )
            _hist_dict["mu_pfRelIso04_all"] = Hist.Hist(
                syst_axis, iso_axis, Hist.storage.Weight()
            )
            _hist_dict[f"muon_dxy"] = Hist.Hist(
                syst_axis, dxy_axis, Hist.storage.Weight()
            )
            _hist_dict[f"muon_dz"] = Hist.Hist(
                syst_axis, dz_axis, Hist.storage.Weight()
            )
            _hist_dict[f"mu_ptratio"] = Hist.Hist(
                syst_axis, flav_axis, osss_axis, ptratio_axis, Hist.storage.Weight()
            )

            _hist_dict[f"{i}_dz"] = Hist.Hist(syst_axis, dz_axis, Hist.storage.Weight())

    ### Common kinematic variables histogram creation
    if "Wc_sf" not in workflow:
        _hist_dict["njet"] = Hist.Hist(syst_axis, n_axis, Hist.storage.Weight())
        if "ctag_tt" in workflow:
            _hist_dict["nmujet"] = Hist.Hist(syst_axis, n_axis, Hist.storage.Weight())
            _hist_dict["nsoftmu"] = Hist.Hist(syst_axis, n_axis, Hist.storage.Weight())
        for obj in obj_list:
            if "jet" in obj or "soft_l" in obj:
                if obj == "soft_l":
                    _hist_dict["soft_l_pt"] = Hist.Hist(
                        syst_axis, flav_axis, softlpt_axis, Hist.storage.Weight()
                    )
                else:
                    _hist_dict[f"{obj}_pt"] = Hist.Hist(
                        syst_axis, flav_axis, jpt_axis, Hist.storage.Weight()
                    )
                _hist_dict[f"{obj}_eta"] = Hist.Hist(
                    syst_axis, flav_axis, eta_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{obj}_phi"] = Hist.Hist(
                    syst_axis, flav_axis, phi_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{obj}_mass"] = Hist.Hist(
                    syst_axis, flav_axis, mass_axis, Hist.storage.Weight()
                )
            else:
                _hist_dict[f"{obj}_pt"] = Hist.Hist(
                    syst_axis, pt_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{obj}_phi"] = Hist.Hist(
                    syst_axis, phi_axis, Hist.storage.Weight()
                )
                if obj != "MET":
                    _hist_dict[f"{obj}_eta"] = Hist.Hist(
                        syst_axis, eta_axis, Hist.storage.Weight()
                    )
    else:
        _hist_dict["njet"] = Hist.Hist(
            syst_axis, osss_axis, n_axis, Hist.storage.Weight()
        )
        _hist_dict["nmujet"] = Hist.Hist(
            syst_axis, osss_axis, n_axis, Hist.storage.Weight()
        )
        _hist_dict["nsoftmu"] = Hist.Hist(
            syst_axis, osss_axis, n_axis, Hist.storage.Weight()
        )
        for obj in obj_list:
            if "jet" in obj or "soft_l" in obj:
                if obj == "soft_l":
                    _hist_dict["soft_l_pt"] = Hist.Hist(
                        syst_axis,
                        flav_axis,
                        osss_axis,
                        softlpt_axis,
                        Hist.storage.Weight(),
                    )
                else:
                    _hist_dict[f"{obj}_pt"] = Hist.Hist(
                        syst_axis, flav_axis, osss_axis, jpt_axis, Hist.storage.Weight()
                    )
                _hist_dict[f"{obj}_eta"] = Hist.Hist(
                    syst_axis, flav_axis, osss_axis, eta_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{obj}_phi"] = Hist.Hist(
                    syst_axis, flav_axis, osss_axis, phi_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{obj}_mass"] = Hist.Hist(
                    syst_axis, flav_axis, osss_axis, mass_axis, Hist.storage.Weight()
                )
            else:
                _hist_dict[f"{obj}_pt"] = Hist.Hist(
                    syst_axis, osss_axis, pt_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{obj}_phi"] = Hist.Hist(
                    syst_axis, osss_axis, phi_axis, Hist.storage.Weight()
                )
                if obj != "MET":
                    _hist_dict[f"{obj}_eta"] = Hist.Hist(
                        syst_axis, osss_axis, eta_axis, Hist.storage.Weight()
                    )

    ### Btag input variables & PFCands
    bininfo = definitions()
    for d in bininfo.keys():
        if d not in events.Jet.fields:
            continue
        ranges = bininfo[d]["manual_ranges"]
        binning = bininfo[d]["bins"]
        labels = (
            bininfo[d]["displayname"] + " [" + bininfo[d]["inputVar_units"] + "]"
            if bininfo[d]["inputVar_units"] is not None
            else bininfo[d]["displayname"]
        )
        if "Wc_sf" in workflow:
            _hist_dict[d] = Hist.Hist(
                syst_axis,
                flav_axis,
                osss_axis,
                Hist.axis.Regular(binning, ranges[0], ranges[1], name=d, label=labels),
                Hist.storage.Weight(),
            )
        else:
            _hist_dict[d] = Hist.Hist(
                syst_axis,
                flav_axis,
                Hist.axis.Regular(binning, ranges[0], ranges[1], name=d, label=labels),
                Hist.storage.Weight(),
            )
    ### JetSVs variables
    SV_bininfo = SV_definitions()
    for d in SV_bininfo.keys():
        ranges = SV_bininfo[d]["manual_ranges"]
        binning = SV_bininfo[d]["bins"]
        labels = (
            SV_bininfo[d]["displayname"] + " [" + SV_bininfo[d]["inputVar_units"] + "]"
            if SV_bininfo[d]["inputVar_units"] is not None
            else SV_bininfo[d]["displayname"]
        )
        _hist_dict[d] = Hist.Hist(
            syst_axis,
            flav_axis,
            Hist.axis.Regular(binning, ranges[0], ranges[1], name=d, label=labels),
            Hist.storage.Weight(),
        )
    ### discriminators
    disc_list = [
        "btagDeepFlavB",
        "btagDeepFlavC",
        "btagDeepFlavCvL",
        "btagDeepFlavCvB",
        "btagDeepFlavB_b",
        "btagDeepFlavB_bb",
        "btagPNetB",
        "btagPNetCvB",
        "btagPNetCvL",
        "btagPNetCvNotB",
        "btagPNetProbB",
        "btagPNetProbC",
        "btagPNetProbG",
        "btagPNetProbUDS",
        "btagPNetQvG",
        "btagPNetTauVJet",
        "btagRobustParTAK4B",
        "btagRobustParTAK4B_b",
        "btagRobustParTAK4B_bb",
        "btagRobustParTAK4B_lepb",
        "btagRobustParTAK4C",
        "btagRobustParTAK4G",
        "btagRobustParTAK4UDS",
        "btagRobustParTAK4CvB",
        "btagRobustParTAK4CvL",
        "btagRobustParTAK4QG",
        "btagUParTAK4B",
        "btagUParTAK4CvL",
        "btagUParTAK4CvB",
        "btagUParTAK4CvNotB",
        "btagUParTAK4QvG",
        "btagUParTAK4TauVJet",
        ## Negative tagger
        "btagNegDeepFlavB",
        "btagNegDeepFlavB_b",
        "btagNegDeepFlavB_bb",
        "btagNegDeepFlavB_lepb",
        "btagNegDeepFlavC",
        "btagNegDeepFlavCvB",
        "btagNegDeepFlavCvL",
        "btagNegDeepFlavG",
        "btagNegDeepFlavQG",
        "btagNegDeepFlavUDS",
        "btagNegPNetB",
        "btagNegPNetC",
        "btagNegPNetCvB",
        "btagNegPNetCvL",
        "btagNegPNetProbB",
        "btagNegPNetProbC",
        "btagNegPNetProbG",
        "btagNegPNetProbUDS",
        "btagNegRobustParTAK4B",
        "btagNegRobustParTAK4B_b",
        "btagNegRobustParTAK4B_bb",
        "btagNegRobustParTAK4B_lepb",
        "btagNegRobustParTAK4C",
        "btagNegRobustParTAK4CvB",
        "btagNegRobustParTAK4CvL",
        "btagNegRobustParTAK4G",
        "btagNegRobustParTAK4QG",
        "btagNegRobustParTAK4UDS",
        # other prob info
        "PNetRegPtRawCorr",
        "PNetRegPtRawCorrNeutrino",
        "PNetRegPtRawRes",
        "UParTAK4RegPtRawRes",
        "UParTAK4RegPtRawCorrNeutrino",
        "UParTAK4RegPtRawCorr",
        "Bprob",
        "BprobN",
        "ProbaN",
    ]
    for disc in disc_list:
        if disc not in events.Jet.fields and "Trans" not in disc:
            continue
        njet = 1
        if "ttdilep_sf" in workflow:
            njet = 2
        elif "ttsemilep_sf" in workflow:
            njet = 4
        for i in range(njet):
            if "Wc_sf" in workflow:
                if "Trans" in disc:
                    _hist_dict[f"{disc}_{i}"] = Hist.Hist(
                        syst_axis,
                        flav_axis,
                        osss_axis,
                        Hist.axis.Regular(40, 0, 8, name="discr", label=disc),
                        Hist.storage.Weight(),
                    )
                elif "btag" in disc or "ProbaN" == disc:
                    _hist_dict[f"{disc}_{i}"] = Hist.Hist(
                        syst_axis,
                        flav_axis,
                        osss_axis,
                        Hist.axis.Regular(50, 0.0, 1, name="discr", label=disc),
                        Hist.storage.Weight(),
                    )
                elif "Bprob" in disc:
                    _hist_dict[f"{disc}_{i}"] = Hist.Hist(
                        syst_axis,
                        flav_axis,
                        osss_axis,
                        Hist.axis.Regular(50, 0, 10, name="discr", label=disc),
                        Hist.storage.Weight(),
                    )
                elif "Res" in disc:
                    _hist_dict[f"{disc}_{i}"] = Hist.Hist(
                        syst_axis,
                        flav_axis,
                        osss_axis,
                        Hist.axis.Regular(40, 0, 1, name="discr", label=disc),
                        Hist.storage.Weight(),
                    )
                elif "Corr" in disc:
                    _hist_dict[f"{disc}_{i}"] = Hist.Hist(
                        syst_axis,
                        flav_axis,
                        osss_axis,
                        Hist.axis.Regular(40, 0, 2, name="discr", label=disc),
                        Hist.storage.Weight(),
                    )

            else:
                if "Trans" in disc:
                    _hist_dict[f"{disc}_{i}"] = Hist.Hist(
                        syst_axis,
                        flav_axis,
                        Hist.axis.Regular(40, 0, 8, name="discr", label=disc),
                        Hist.storage.Weight(),
                    )
                if "btag" in disc or "ProbaN" == disc:
                    _hist_dict[f"{disc}_{i}"] = Hist.Hist(
                        syst_axis,
                        flav_axis,
                        Hist.axis.Regular(50, 0.0, 1, name="discr", label=disc),
                        Hist.storage.Weight(),
                    )
                elif "Bprob" in disc:
                    _hist_dict[f"{disc}_{i}"] = Hist.Hist(
                        syst_axis,
                        flav_axis,
                        Hist.axis.Regular(50, 0, 10, name="discr", label=disc),
                        Hist.storage.Weight(),
                    )
                elif "Res" == disc:
                    _hist_dict[f"{disc}_{i}"] = Hist.Hist(
                        syst_axis,
                        flav_axis,
                        Hist.axis.Regular(40, 0, 1, name="discr", label=disc),
                        Hist.storage.Weight(),
                    )
                elif "Corr" in disc:
                    _hist_dict[f"{disc}_{i}"] = Hist.Hist(
                        syst_axis,
                        flav_axis,
                        Hist.axis.Regular(40, 0, 2, name="discr", label=disc),
                        Hist.storage.Weight(),
                    )
    return _hist_dict
