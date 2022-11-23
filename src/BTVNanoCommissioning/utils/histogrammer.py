from BTVNanoCommissioning.helpers.definitions import definitions
import hist as Hist


def histogrammer(workflow):
    _hist_dict = {}
    ## Common variables
    flav_axis = Hist.axis.IntCategory([0, 1, 4, 5, 6], name="flav", label="Genflavour")
    syst_axis = Hist.axis.StrCategory([], name="syst", growth=True)
    pt_axis = Hist.axis.Regular(50, 0, 200, name="pt", label=" $p_{T}$ [GeV]")
    jpt_axis = Hist.axis.Regular(50, 0, 300, name="pt", label=" $p_{T}$ [GeV]")
    softlpt_axis = Hist.axis.Regular(25, 0, 25, name="pt", label=" $p_{T}$ [GeV]")
    mass_axis = Hist.axis.Regular(50, 0, 300, name="mass", label=" $p_{T}$ [GeV]")
    eta_axis = Hist.axis.Regular(25, -2.5, 2.5, name="eta", label=" $\eta$")
    phi_axis = Hist.axis.Regular(30, -3, 3, name="phi", label="$\phi$")
    mt_axis = Hist.axis.Regular(30, 0, 300, name="mt", label=" $m_{T}$ [GeV]")
    iso_axis = Hist.axis.Regular(30, 0, 0.05, name="pfRelIso03_all", label="Rel. Iso")
    softliso_axis = Hist.axis.Regular(
        20, 0.2, 6.2, name="pfRelIso03_all", label="Rel. Iso"
    )
    dr_axis = Hist.axis.Regular(20, 0, 8, name="dr", label="$\Delta$R")
    dxy_axis = Hist.axis.Regular(40, -0.05, 0.05, name="dxy", label="d_{xy}")
    dz_axis = Hist.axis.Regular(40, -0.01, 0.01, name="dz", label="d_{z}")
    qcddxy_axis = Hist.axis.Regular(40, -0.002, 0.002, name="dxy", label="d_{xy}")
    sip3d_axis = Hist.axis.Regular(20, 0, 0.2, name="sip3d", label="SIP 3D")
    ptratio_axis = Hist.axis.Regular(50, 0, 1, name="ratio", label="ratio")
    n_axis = Hist.axis.IntCategory(
        [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], name="n", label="N obj"
    )
    osss_axis = Hist.axis.IntCategory([1, -1], name="osss", label="OS(+)/SS(-)")
    ### Workflow specific
    if "validation" == workflow:
        obj_list = ["jet0", "jet1"]
    elif "ttcom" == workflow:
        obj_list = ["mu", "ele"]
        for i in range(2):
            obj_list.append(f"jet{i}")
            _hist_dict[f"dr_mujet{i}"] = Hist.Hist(
                flav_axis, dr_axis, Hist.storage.Weight()
            )
    elif "ttdilep_sf" == workflow:
        obj_list = ["mu", "ele"]
        for i in range(2):
            obj_list.append(f"jet{i}")
            _hist_dict[f"dr_mujet{i}"] = Hist.Hist(
                flav_axis, dr_axis, Hist.storage.Weight()
            )
        for i in ["mu", "ele"]:
            if i == "mu":
                _hist_dict[f"{i}_pfRelIso04_all"] = Hist.Hist(
                    iso_axis, Hist.storage.Weight()
                )
            else:
                _hist_dict[f"{i}_pfRelIso03_all"] = Hist.Hist(
                    iso_axis, Hist.storage.Weight()
                )
            _hist_dict[f"{i}_dxy"] = Hist.Hist(dxy_axis, Hist.storage.Weight())
            _hist_dict[f"{i}_dz"] = Hist.Hist(dz_axis, Hist.storage.Weight())
    elif "ttsemilep_sf" == workflow:
        obj_list = ["mu", "MET"]
        for i in range(4):
            obj_list.append(f"jet{i}")
            _hist_dict[f"dr_mujet{i}"] = Hist.Hist(
                flav_axis, dr_axis, Hist.storage.Weight()
            )
        for i in ["mu"]:
            _hist_dict[f"{i}_pfRelIso04_all"] = Hist.Hist(
                iso_axis, Hist.storage.Weight()
            )
            _hist_dict[f"{i}_dxy"] = Hist.Hist(dxy_axis, Hist.storage.Weight())
            _hist_dict[f"{i}_dz"] = Hist.Hist(dz_axis, Hist.storage.Weight())
    elif "ctag_ttdilep_sf" in workflow:
        obj_list = ["hl", "sl", "soft_l", "MET", "z", "lmujet"]
        _hist_dict["z_mass"] = Hist.Hist(
            Hist.axis.Regular(50, 50, 100, name="mass", label=" $m_Z$ [GeV]"),
            Hist.storage.Weight(),
        )
        # delta R between soft muon and mu-jet
        _hist_dict["dr_lmujetsmu"] = Hist.Hist(
            flav_axis, dr_axis, Hist.storage.Weight()
        )
        # delta R between hard muon and mu-jet
        _hist_dict["dr_lmujethmu"] = Hist.Hist(
            flav_axis, dr_axis, Hist.storage.Weight()
        )
        # delta R between soft muon and hard muon
        _hist_dict["dr_lmusmu"] = Hist.Hist(dr_axis, Hist.storage.Weight())
        for i in ["hl", "sl", "soft_l"]:
            if i == "soft_l":
                _hist_dict[f"soft_l_pfRelIso04_all"] = Hist.Hist(
                    flav_axis, softliso_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{i}_dxy"] = Hist.Hist(
                    flav_axis, dxy_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{i}_dz"] = Hist.Hist(
                    flav_axis, dz_axis, Hist.storage.Weight()
                )
            else:
                _hist_dict[f"{i}_pfRelIso04_all"] = Hist.Hist(
                    iso_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{i}_dxy"] = Hist.Hist(qcddxy_axis, Hist.storage.Weight())
                _hist_dict[f"{i}_dz"] = Hist.Hist(dz_axis, Hist.storage.Weight())
            # lepton / jet pT ratio
            _hist_dict[f"{i}_ptratio"] = Hist.Hist(
                flav_axis, ptratio_axis, Hist.storage.Weight()
            )

    elif "ctag_ttsemilep_sf" in workflow:
        obj_list = ["hl", "soft_l", "MET", "z", "w", "mujet"]
        _hist_dict["z_mass"] = Hist.Hist(
            Hist.axis.Regular(50, 50, 100, name="mass", label="$m_Z$ [GeV]"),
            Hist.storage.Weight(),
        )
        _hist_dict["w_mass"] = Hist.Hist(
            Hist.axis.Regular(50, 50, 100, name="mass", label="$m_W$ [GeV]"),
            Hist.storage.Weight(),
        )
        # delta R between soft muon and mu-jet
        _hist_dict["dr_lmujetsmu"] = Hist.Hist(
            flav_axis, dr_axis, Hist.storage.Weight()
        )
        # delta R between hard muon and mu-jet
        _hist_dict["dr_lmujethmu"] = Hist.Hist(
            flav_axis, dr_axis, Hist.storage.Weight()
        )
        # delta R between hard muon and soft-muon
        _hist_dict["dr_lmusmu"] = Hist.Hist(dr_axis, Hist.storage.Weight())
        for i in ["hl", "soft_l"]:
            if i == "soft_l":
                _hist_dict[f"soft_l_pfRelIso04_all"] = Hist.Hist(
                    flav_axis, softliso_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{i}_dxy"] = Hist.Hist(
                    flav_axis, dxy_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{i}_dz"] = Hist.Hist(
                    flav_axis, dz_axis, Hist.storage.Weight()
                )
            else:
                _hist_dict[f"{i}_pfRelIso04_all"] = Hist.Hist(
                    iso_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{i}_dxy"] = Hist.Hist(dxy_axis, Hist.storage.Weight())
                _hist_dict[f"{i}_dz"] = Hist.Hist(dz_axis, Hist.storage.Weight())
            _hist_dict[f"{i}_ptratio"] = Hist.Hist(
                flav_axis, ptratio_axis, Hist.storage.Weight()
            )

    elif "Wc_sf" in workflow:
        obj_list = ["hl", "soft_l", "MET", "z", "w", "mujet"]
        _hist_dict["SV_charge"] = Hist.Hist(
            osss_axis,
            Hist.axis.Regular(20, -10, 10, name="charge", label="SV charge"),
            Hist.storage.Weight(),
        )
        _hist_dict["z_mass"] = Hist.Hist(
            osss_axis,
            Hist.axis.Regular(50, 50, 100, name="mass", label="$m_Z$ [GeV]"),
            Hist.storage.Weight(),
        )
        _hist_dict["w_mass"] = Hist.Hist(
            osss_axis,
            Hist.axis.Regular(50, 50, 100, name="mass", label="$m_W$ [GeV]"),
            Hist.storage.Weight(),
        )
        # delta R between soft muon and mu-jet
        _hist_dict["dr_lmujetsmu"] = Hist.Hist(
            flav_axis, osss_axis, dr_axis, Hist.storage.Weight()
        )
        # delta R between hard muon and mu-jet
        _hist_dict["dr_lmujethmu"] = Hist.Hist(
            flav_axis, osss_axis, dr_axis, Hist.storage.Weight()
        )
        # delta R between hard muon and soft-muon
        _hist_dict["dr_lmusmu"] = Hist.Hist(osss_axis, dr_axis, Hist.storage.Weight())
        for i in ["hl", "soft_l"]:
            if i == "soft_l":
                _hist_dict[f"soft_l_pfRelIso04_all"] = Hist.Hist(
                    flav_axis, osss_axis, softliso_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{i}_dxy"] = Hist.Hist(
                    flav_axis, osss_axis, dxy_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{i}_dz"] = Hist.Hist(
                    flav_axis, osss_axis, dz_axis, Hist.storage.Weight()
                )
            else:
                _hist_dict[f"{i}_pfRelIso04_all"] = Hist.Hist(
                    osss_axis, iso_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{i}_dxy"] = Hist.Hist(
                    osss_axis, qcddxy_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{i}_dz"] = Hist.Hist(
                    osss_axis, dz_axis, Hist.storage.Weight()
                )
            _hist_dict[f"{i}_ptratio"] = Hist.Hist(
                flav_axis, osss_axis, ptratio_axis, Hist.storage.Weight()
            )

    elif "DY_sf" in workflow:
        obj_list = ["posl", "negl", "z", "jet"]
        _hist_dict["z_mass"] = Hist.Hist(
            Hist.axis.Regular(50, 50, 100, name="mass", label="$m_Z$ [GeV]"),
            Hist.storage.Weight(),
        )
        _hist_dict["dr_mumu"] = Hist.Hist(dr_axis, Hist.storage.Weight())
        for i in ["posl", "negl"]:
            _hist_dict[f"{i}_pfRelIso04_all"] = Hist.Hist(
                iso_axis, Hist.storage.Weight()
            )
            _hist_dict[f"{i}_dxy"] = Hist.Hist(dxy_axis, Hist.storage.Weight())
            _hist_dict[f"{i}_dz"] = Hist.Hist(dz_axis, Hist.storage.Weight())
    ### Common kinematic variables
    if "Wc_sf" not in workflow:
        _hist_dict["njet"] = Hist.Hist(n_axis, Hist.storage.Weight())
        for obj in obj_list:
            if "jet" in obj or "soft_l" in obj:
                if obj == "soft_l":
                    _hist_dict["soft_l_pt"] = Hist.Hist(
                        flav_axis, softlpt_axis, Hist.storage.Weight()
                    )
                else:
                    _hist_dict[f"{obj}_pt"] = Hist.Hist(
                        flav_axis, jpt_axis, Hist.storage.Weight()
                    )
                _hist_dict[f"{obj}_eta"] = Hist.Hist(
                    flav_axis, eta_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{obj}_phi"] = Hist.Hist(
                    flav_axis, phi_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{obj}_mass"] = Hist.Hist(
                    flav_axis, mass_axis, Hist.storage.Weight()
                )
            else:
                _hist_dict[f"{obj}_pt"] = Hist.Hist(pt_axis, Hist.storage.Weight())
                _hist_dict[f"{obj}_phi"] = Hist.Hist(phi_axis, Hist.storage.Weight())
                if obj != "MET":
                    _hist_dict[f"{obj}_eta"] = Hist.Hist(
                        eta_axis, Hist.storage.Weight()
                    )
    else:
        _hist_dict["njet"] = Hist.Hist(osss_axis, n_axis, Hist.storage.Weight())
        for obj in obj_list:
            if "jet" in obj or "soft_l" in obj:
                if obj == "soft_l":
                    _hist_dict["soft_l_pt"] = Hist.Hist(
                        flav_axis, osss_axis, softlpt_axis, Hist.storage.Weight()
                    )
                _hist_dict[f"{obj}_pt"] = Hist.Hist(
                    flav_axis, osss_axis, jpt_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{obj}_eta"] = Hist.Hist(
                    flav_axis, osss_axis, eta_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{obj}_phi"] = Hist.Hist(
                    flav_axis, osss_axis, phi_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{obj}_mass"] = Hist.Hist(
                    flav_axis, osss_axis, mass_axis, Hist.storage.Weight()
                )
            else:
                _hist_dict[f"{obj}_pt"] = Hist.Hist(
                    osss_axis, pt_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{obj}_phi"] = Hist.Hist(
                    osss_axis, phi_axis, Hist.storage.Weight()
                )
                if obj != "MET":
                    _hist_dict[f"{obj}_eta"] = Hist.Hist(
                        osss_axis, eta_axis, Hist.storage.Weight()
                    )

    ### Btag variables
    bininfo = definitions()
    for d in bininfo.keys():
        ranges = bininfo[d]["manual_ranges"]
        binning = bininfo[d]["bins"]
        labels = (
            bininfo[d]["displayname"] + " [" + bininfo[d]["inputVar_units"] + "]"
            if bininfo[d]["inputVar_units"] is not None
            else bininfo[d]["displayname"]
        )
        if "Wc_sf" in workflow:
            _hist_dict[d] = Hist.Hist(
                flav_axis,
                osss_axis,
                Hist.axis.Regular(binning, ranges[0], ranges[1], name=d, label=labels),
                Hist.storage.Weight(),
            )
        else:
            _hist_dict[d] = Hist.Hist(
                flav_axis,
                Hist.axis.Regular(binning, ranges[0], ranges[1], name=d, label=labels),
                Hist.storage.Weight(),
            )
    ### discriminators
    disc_list = [
        "btagDeepB",
        "btagDeepC",
        "btagDeepFlavB",
        "btagDeepFlavC",
        "btagDeepCvL",
        "btagDeepCvB",
        "btagDeepFlavCvL",
        "btagDeepFlavCvB",
        "btagDeepB_b",
        "btagDeepB_bb",
        "btagDeepFlavB_b",
        "btagDeepFlavB_bb",
        "btagDeepFlavB_lepb",
    ]
    for disc in disc_list:
        njet = 1
        if "ttdilep_sf" in workflow:
            njet = 2
        elif "ttsemilep_sf" in workflow:
            njet = 4
        for i in range(njet):
            if "Wc_sf" in workflow:
                _hist_dict[f"{disc}_{i}"] = Hist.Hist(
                    flav_axis,
                    osss_axis,
                    syst_axis,
                    Hist.axis.Regular(30, -0.2, 1, name="discr", label=disc),
                    Hist.storage.Weight(),
                )
            else:
                _hist_dict[f"{disc}_{i}"] = Hist.Hist(
                    flav_axis,
                    syst_axis,
                    Hist.axis.Regular(30, -0.2, 1, name="discr", label=disc),
                    Hist.storage.Weight(),
                )
    return _hist_dict
