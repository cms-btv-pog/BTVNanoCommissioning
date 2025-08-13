from BTVNanoCommissioning.utils.selection import btag_wp_dict
from BTVNanoCommissioning.helpers.definitions import (
    definitions,
    SV_definitions,
    disc_list,
)
import hist as Hist
import awkward as ak
from BTVNanoCommissioning.helpers.func import flatten


def histogrammer(events, workflow, year="2022", campaign="Summer22"):
    """
    Most of workflows require same set of variables. Collect axis, histograms definition in single file
    To contribute: Add additional axis, histogram using [hist](https://hist.readthedocs.io/en/latest/) for dedicated workflow into the `_hist_dict`. For the new histogram, please have the `syst_axis` as first axis, and `Weight` as last axis.

    Parameters:
    events (awkward.Array): The events data to be histogrammed.
    workflow (str): The workflow identifier to determine specific histogramming logic.

    Example:

    ```python
    # axis example
    mass_axis = Hist.axis.Regular(50, 0, 300, name="mass", label=" $p_{T}$ [GeV]")
    # hist example
    _hist_dict["dr_lmusmu"] = Hist.Hist(syst_axis, dr_axis, Hist.storage.Weight())
    ```

    Returns:
    dict: A dictionary containing the defined histograms.
    """

    _hist_dict = {}
    ## Common axis
    flav_axis = Hist.axis.IntCategory([0, 1, 4, 5, 6], name="flav", label="Genflavour")
    syst_axis = Hist.axis.StrCategory([], name="syst", growth=True)
    pt_axis = Hist.axis.Regular(60, 0, 300, name="pt", label=" $p_{T}$ [GeV]")
    jpt_axis = Hist.axis.Regular(300, 0, 3000, name="pt", label=" $p_{T}$ [GeV]")
    softlpt_axis = Hist.axis.Regular(25, 0, 25, name="pt", label=" $p_{T}$ [GeV]")
    mass_axis = Hist.axis.Regular(50, 0, 300, name="mass", label=" mass [GeV]")
    bdt_axis = Hist.axis.Regular(50, 0, 1, name="bdt", label=" BDT discriminant")
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
    ## create histograms for each workflow
    ### Workflow specific
    if "example" == workflow:
        obj_list = [
            "jet",
            "mu",
        ]  # store basic 4-vector, pt,eta, phi, mass for the object
        _hist_dict[f"dr_mujet0"] = Hist.Hist(
            syst_axis, flav_axis, dr_axis, Hist.storage.Weight()
        )  # create cutstomize histogram
    elif "QCD" == workflow:
        obj_list = ["jet0"]
        # FIXME: commented SVJet related histogram until fixing linkinf of BTVNano
        # _hist_dict["dr_SVjet0"] = Hist.Hist(
        #     syst_axis, flav_axis, dr_SV_axis, Hist.storage.Weight()
        # )x
        # _hist_dict["nJetSVs"] = Hist.Hist(syst_axis, n_axis, Hist.storage.Weight())

    elif "QCD_smu" == workflow:
        obj_list = ["mujet", "hl", "soft_l"]
        # _hist_dict["dr_SVjet0"] = Hist.Hist(
        #     syst_axis, flav_axis, dr_SV_axis, Hist.storage.Weight()
        # )
        # _hist_dict["nJetSVs"] = Hist.Hist(syst_axis, n_axis, Hist.storage.Weight())
        _hist_dict["nmujet"] = Hist.Hist(syst_axis, n_axis, Hist.storage.Weight())
        _hist_dict["nssmu"] = Hist.Hist(syst_axis, n_axis, Hist.storage.Weight())
        _hist_dict["dr_lmujetsmu"] = Hist.Hist(
            syst_axis, flav_axis, dr_s_axis, Hist.storage.Weight()
        )
        for i in ["soft_l"]:
            if i == "soft_l":
                _hist_dict[f"soft_l_pfRelIso04_all"] = Hist.Hist(
                    syst_axis,
                    flav_axis,
                    softliso_axis,
                    Hist.storage.Weight(),
                )
                _hist_dict[f"{i}_dxy"] = Hist.Hist(
                    syst_axis, flav_axis, dxy_axis, Hist.storage.Weight()
                )
                _hist_dict[f"{i}_dz"] = Hist.Hist(
                    syst_axis, flav_axis, dz_axis, Hist.storage.Weight()
                )
            _hist_dict[f"{i}_ptratio"] = Hist.Hist(
                syst_axis, flav_axis, ptratio_axis, Hist.storage.Weight()
            )
        _hist_dict["dr_lmusmujetsmu"] = Hist.Hist(
            syst_axis, flav_axis, dr_s_axis, Hist.storage.Weight()
        )

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
        obj_list.append("cjet")
        _hist_dict["dr_cjet"] = Hist.Hist(
            syst_axis, flav_axis, dr_axis, Hist.storage.Weight()
        )
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

    elif "c_ttsemilep_sf" == workflow:
        obj_list = ["mu", "MET"]
        obj_list.append("cjet")
        _hist_dict["dr_cjet"] = Hist.Hist(
            syst_axis, flav_axis, dr_axis, Hist.storage.Weight()
        )
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
        obj_list = ["hl", "sl", "soft_l", "MET", "dilep", "lmujet"]
        _hist_dict["dilep_mass"] = Hist.Hist(
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
        _hist_dict["dr_hmusmu"] = Hist.Hist(syst_axis, dr_axis, Hist.storage.Weight())
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
                else:
                    _hist_dict[f"{i}_pfRelIso03_all"] = Hist.Hist(
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
        obj_list = ["hl", "soft_l", "MET", "dilep", "mujet"]
        _hist_dict["dilep_mass"] = Hist.Hist(
            syst_axis,
            Hist.axis.Regular(50, 50, 100, name="mass", label="$m_{\\ell\\ell}$ [GeV]"),
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
        _hist_dict["dr_hmusmu"] = Hist.Hist(syst_axis, dr_axis, Hist.storage.Weight())
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
        obj_list = ["hl", "soft_l", "MET", "dilep", "mujet"]
        _hist_dict["SV_charge"] = Hist.Hist(
            syst_axis,
            osss_axis,
            Hist.axis.Regular(20, -10, 10, name="charge", label="SV charge"),
            Hist.storage.Weight(),
        )
        _hist_dict["dilep_mass"] = Hist.Hist(
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
        _hist_dict["w_pt"] = Hist.Hist(
            syst_axis,
            osss_axis,
            Hist.axis.Regular(50, 0, 300, name="pT", label="$p_T^{\\ell\\nu}$ [GeV]"),
            Hist.storage.Weight(),
        )
        _hist_dict["w_phi"] = Hist.Hist(
            syst_axis,
            osss_axis,
            phi_axis,
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
        _hist_dict["dr_hmusmu"] = Hist.Hist(
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
        obj_list = ["posl", "negl", "dilep", "jet0"]
        _hist_dict["dilep_mass"] = Hist.Hist(
            syst_axis,
            Hist.axis.Regular(50, 50, 100, name="mass", label="$m_{\\ell\\ell}$ [GeV]"),
            Hist.storage.Weight(),
        )
        _hist_dict["dr_poslnegl"] = Hist.Hist(syst_axis, dr_axis, Hist.storage.Weight())
        _hist_dict["dr_posljet"] = Hist.Hist(syst_axis, dr_axis, Hist.storage.Weight())
        _hist_dict["dr_negljet"] = Hist.Hist(syst_axis, dr_axis, Hist.storage.Weight())
        for i in ["posl", "negl"]:
            if "m" in workflow:
                _hist_dict[f"{i}_pfRelIso04_all"] = Hist.Hist(
                    syst_axis, iso_axis, Hist.storage.Weight()
                )
            _hist_dict[f"{i}_dxy"] = Hist.Hist(
                syst_axis, dxy_axis, Hist.storage.Weight()
            )
            _hist_dict[f"{i}_dz"] = Hist.Hist(syst_axis, dz_axis, Hist.storage.Weight())
            _hist_dict[f"dr_{i}jet"] = Hist.Hist(
                syst_axis, flav_axis, dr_axis, Hist.storage.Weight()
            )
    elif "sf_ttdilep_kin" in workflow:
        obj_list = ["dilep"]

        _hist_dict["kindisc"] = Hist.Hist(
            syst_axis, flav_axis, bdt_axis, Hist.storage.Weight()
        )

        _hist_dict["close_mlj"] = Hist.Hist(
            syst_axis, flav_axis, mass_axis, Hist.storage.Weight()
        )
        _hist_dict["close_deta"] = Hist.Hist(
            syst_axis, flav_axis, eta_axis, Hist.storage.Weight()
        )
        _hist_dict["close_dphi"] = Hist.Hist(
            syst_axis, flav_axis, phi_axis, Hist.storage.Weight()
        )
        _hist_dict["close_ptrel"] = Hist.Hist(
            syst_axis, flav_axis, pt_axis, Hist.storage.Weight()
        )
        _hist_dict["close_lj2ll_deta"] = Hist.Hist(
            syst_axis, flav_axis, eta_axis, Hist.storage.Weight()
        )
        _hist_dict["close_lj2ll_dphi"] = Hist.Hist(
            syst_axis, flav_axis, phi_axis, Hist.storage.Weight()
        )

        _hist_dict["far_mlj"] = Hist.Hist(
            syst_axis, flav_axis, mass_axis, Hist.storage.Weight()
        )
        _hist_dict["far_deta"] = Hist.Hist(
            syst_axis, flav_axis, eta_axis, Hist.storage.Weight()
        )
        _hist_dict["far_dphi"] = Hist.Hist(
            syst_axis, flav_axis, phi_axis, Hist.storage.Weight()
        )
        _hist_dict["far_ptrel"] = Hist.Hist(
            syst_axis, flav_axis, pt_axis, Hist.storage.Weight()
        )
        _hist_dict["far_lj2ll_deta"] = Hist.Hist(
            syst_axis, flav_axis, eta_axis, Hist.storage.Weight()
        )
        _hist_dict["far_lj2ll_dphi"] = Hist.Hist(
            syst_axis, flav_axis, phi_axis, Hist.storage.Weight()
        )

        _hist_dict["j2ll_deta"] = Hist.Hist(
            syst_axis, flav_axis, eta_axis, Hist.storage.Weight()
        )
        _hist_dict["j2ll_dphi"] = Hist.Hist(
            syst_axis, flav_axis, phi_axis, Hist.storage.Weight()
        )

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
                        syst_axis, flav_axis, pt_axis, Hist.storage.Weight()
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
            # mujet pt passing tagger WPs
            if "mujet" in obj:
                if "cutbased" in workflow:
                    for tagger in btag_wp_dict[year + "_" + campaign].keys():
                        for wp in btag_wp_dict[year + "_" + campaign][tagger][
                            "c"
                        ].keys():
                            if not "No" in wp:
                                _hist_dict[f"{obj}_pt_{tagger}{wp}"] = Hist.Hist(
                                    syst_axis,
                                    flav_axis,
                                    osss_axis,
                                    pt_axis,
                                    Hist.storage.Weight(),
                                )

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
                        syst_axis, flav_axis, osss_axis, pt_axis, Hist.storage.Weight()
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
    if "QCD_sf" in workflow:
        _hist_dict[f"{obj}_pt"] = Hist.Hist(
            syst_axis, flav_axis, jpt_axis, Hist.storage.Weight()
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
        if "WP" not in workflow:
            break
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
    ### FIXME: Commented out JetSV distrobution until btvnano is fixed
    # SV_bininfo = SV_definitions()
    # for d in SV_bininfo.keys():
    #     ranges = SV_bininfo[d]["manual_ranges"]
    #     binning = SV_bininfo[d]["bins"]
    #     labels = (
    #         SV_bininfo[d]["displayname"] + " [" + SV_bininfo[d]["inputVar_units"] + "]"
    #         if SV_bininfo[d]["inputVar_units"] is not None
    #         else SV_bininfo[d]["displayname"]
    #     )
    #     _hist_dict[d] = Hist.Hist(
    #         syst_axis,
    #         flav_axis,
    #         Hist.axis.Regular(binning, ranges[0], ranges[1], name=d, label=labels),
    #         Hist.storage.Weight(),
    #     )
    ### discriminators
    for disc in disc_list:
        if disc not in events.Jet.fields:
            continue
        njet = 1
        if "ttdilep_sf" in workflow:
            njet = 2
        elif "ttsemilep_sf" in workflow:
            njet = 4
            if "btag" in disc or "ProbaN" == disc:
                _hist_dict[f"c_{disc}"] = Hist.Hist(
                    syst_axis,
                    flav_axis,
                    Hist.axis.Regular(50, 0.0, 1, name="discr", label=disc),
                    Hist.storage.Weight(),
                )
            elif "Bprob" in disc:
                _hist_dict[f"c_{disc}"] = Hist.Hist(
                    syst_axis,
                    flav_axis,
                    Hist.axis.Regular(50, 0, 10, name="discr", label=disc),
                    Hist.storage.Weight(),
                )
            elif "PNetRegPtRawRes" == disc:
                _hist_dict[f"c_{disc}"] = Hist.Hist(
                    syst_axis,
                    flav_axis,
                    Hist.axis.Regular(40, 0, 1, name="discr", label=disc),
                    Hist.storage.Weight(),
                )
            elif "PNetRegPtRawCorr" in disc:
                _hist_dict[f"c_{disc}"] = Hist.Hist(
                    syst_axis,
                    flav_axis,
                    Hist.axis.Regular(40, 0, 2, name="discr", label=disc),
                    Hist.storage.Weight(),
                )
        for i in range(njet):
            if "Wc_sf" in workflow:
                if "btag" in disc or "ProbaN" == disc:
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
                elif "PNetRegPtRawRes" == disc:
                    _hist_dict[f"{disc}_{i}"] = Hist.Hist(
                        syst_axis,
                        flav_axis,
                        Hist.axis.Regular(40, 0, 1, name="discr", label=disc),
                        Hist.storage.Weight(),
                    )
                elif "PNetRegPtRawCorr" in disc:
                    _hist_dict[f"{disc}_{i}"] = Hist.Hist(
                        syst_axis,
                        flav_axis,
                        Hist.axis.Regular(40, 0, 2, name="discr", label=disc),
                        Hist.storage.Weight(),
                    )
    return _hist_dict


# Filled common histogram
def histo_writter(pruned_ev, output, weights, systematics, isSyst, SF_map):
    """
    Write histograms to the output dictionary based on pruned events and other parameters.

    This function processes the pruned events and writes the histograms to the `output` dictionary. It takes into account the weights, systematics, and scale factors.

    Parameters:
    pruned_ev (coffea.nanoaodevents): The pruned events data to be histogrammed.
    output (dict): The output dictionary where histograms will be stored.
    weights (coffea.analysis_tools.Weights): The weights object for the events.
    systematics (list): A list of systematic variations to be considered.
    isSyst (str,bool): Indicating whether systematic variations are to be applied.
    SF_map (dict): A dictionary containing scale factors for different variables.

    Example:
    ```python
    histo_writter(pruned_ev, output, weights, systematics, isSyst, SF_map)
    ```

    Returns:
    None
    """
    exclude_btv = [
        "DeepCSVC",
        "DeepCSVB",
        "DeepJetB",
        "DeepJetC",
    ]  # exclude b-tag SFs for btag inputs
    # define Jet flavor

    # Reduce the jet to the correct dimension in the plot
    nj = 4 if "jet4" in output.keys() else 2 if "jet2" in output.keys() else 1
    pruned_ev.SelJet = pruned_ev.SelJet if nj == 1 else pruned_ev.SelJet[:, :nj]
    if "var" in str(ak.type(pruned_ev.SelJet.pt)) and nj == 1:
        pruned_ev.SelJet = pruned_ev.SelJet[:, 0]
    if "hadronFlavour" in pruned_ev.SelJet.fields:
        isRealData = False
        genflavor = ak.values_astype(
            pruned_ev.SelJet.hadronFlavour
            + 1
            * (
                (pruned_ev.SelJet.partonFlavour == 0)
                & (pruned_ev.SelJet.hadronFlavour == 0)
            ),
            int,
        )
        if "MuonJet" in pruned_ev.fields:
            smflav = ak.values_astype(
                1
                * (
                    (pruned_ev.MuonJet.partonFlavour == 0)
                    & (pruned_ev.MuonJet.hadronFlavour == 0)
                )
                + pruned_ev.MuonJet.hadronFlavour,
                int,
            )
    else:
        isRealData = True
        genflavor = ak.zeros_like(pruned_ev.SelJet.pt, dtype=int)
        if "MuonJet" in pruned_ev.fields:
            smflav = ak.zeros_like(pruned_ev.MuonJet.pt, dtype=int)
    # Loop over the systematic variations
    for syst in systematics:
        if isSyst == False and syst != "nominal":
            break
        # weight modifications for systematics
        weight = (
            weights.weight()
            if syst == "nominal" or syst not in list(weights.variations)
            else weights.weight(modifier=syst)
        )
        # Loop over the histograms
        for histname, h in output.items():
            # tagger score histograms
            if (
                "Deep" in histname
                and "btag" not in histname
                and histname in pruned_ev.SelJet.fields
            ):

                h.fill(
                    syst,
                    flatten(genflavor),
                    flatten(pruned_ev.SelJet[histname]),
                    weight=flatten(
                        ak.broadcast_arrays(
                            weights.partial_weight(exclude=exclude_btv),
                            pruned_ev.SelJet["pt"],
                        )[0]
                    ),
                )
            # PFcands histograms
            elif (
                "PFCands" in pruned_ev.fields
                and "PFCands" in histname
                and histname.split("_")[1] in pruned_ev.PFCands.fields
            ):
                if "MuonJet" in pruned_ev.fields:
                    h.fill(
                        syst,
                        flatten(
                            ak.broadcast_arrays(smflav, pruned_ev.PFCands["pt"])[0]
                        ),
                        flatten(pruned_ev.PFCands[histname.replace("PFCands_", "")]),
                        weight=flatten(
                            ak.broadcast_arrays(
                                weights.partial_weight(exclude=exclude_btv),
                                pruned_ev.PFCands["pt"],
                            )[0]
                        ),
                    )
                else:
                    h.fill(
                        syst,
                        flatten(
                            ak.broadcast_arrays(
                                genflavor[:, 0], pruned_ev.PFCands["pt"]
                            )[0]
                        ),
                        flatten(pruned_ev.PFCands[histname.replace("PFCands_", "")]),
                        weight=flatten(
                            ak.broadcast_arrays(
                                weights.partial_weight(exclude=exclude_btv),
                                pruned_ev.PFCands["pt"],
                            )[0]
                        ),
                    )
            # leading lepton histograms
            elif (
                "hl_" in histname
                and "hl" in pruned_ev.fields
                and histname.replace("hl_", "") in pruned_ev.hl.fields
            ):

                h.fill(
                    syst,
                    flatten(pruned_ev.hl[histname.replace("hl_", "")]),
                    weight=weight,
                )
            # subleading lepton histograms
            elif (
                "sl_" in histname
                and "sl" in pruned_ev.fields
                and histname.replace("sl_", "") in pruned_ev.sl.fields
            ):
                h.fill(
                    syst,
                    flatten(pruned_ev.sl[histname.replace("sl_", "")]),
                    weight=weight,
                )
            # Selected electron histograms
            elif (
                "ele_" in histname
                and histname.replace("ele_", "") in pruned_ev.SelElectron.fields
                and "Plus" not in pruned_ev.fields
            ):

                h.fill(
                    syst,
                    flatten(pruned_ev.SelElectron[histname.replace("ele_", "")]),
                    weight=weight,
                )
            # Selected muon histograms
            elif (
                "mu_" in histname
                and histname.replace("mu_", "") in pruned_ev.SelMuon.fields
                and "Plus" not in pruned_ev.fields
            ):
                h.fill(
                    syst,
                    flatten(pruned_ev.SelMuon[histname.replace("mu_", "")]),
                    weight=weight,
                )
            # Negatively charged lepton histograms-in DY workflow
            elif (
                "negl_" in histname
                and histname.replace("negl_", "") in pruned_ev.negl.fields
            ):
                h.fill(
                    syst,
                    flatten(pruned_ev.negl[histname.replace("negl_", "")]),
                    weight=weight,
                )
            # Posively charged lepton histograms-in DY workflow
            elif (
                "posl_" in histname
                and histname.replace("posl_", "") in pruned_ev.posl.fields
            ):
                h.fill(
                    syst,
                    flatten(pruned_ev.posl[histname.replace("posl_", "")]),
                    weight=weight,
                )
            # Soft muon histograms
            elif "soft_l" in histname and not "ptratio" in histname:
                h.fill(
                    syst,
                    smflav,
                    flatten(pruned_ev.SoftMuon[histname.replace("soft_l_", "")]),
                    weight=weight,
                )
            elif "njet" == histname:
                output["njet"].fill(syst, pruned_ev.njet, weight=weight)
            # Jet kinmeatics & deltaR between jet and lepton
            elif (
                "jet" in histname and "posl" not in histname and "negl" not in histname
            ):
                for i in range(nj):
                    if f"jet{i}" not in histname:
                        continue
                    if nj == 1:
                        sel_jet, flav = pruned_ev.SelJet, genflavor
                    else:
                        sel_jet, flav = pruned_ev.SelJet[:, i], genflavor[:, i]

                    if str(i) in histname:
                        if "dr_mujet" in histname:
                            h.fill(
                                syst,
                                flatten(flav),
                                flatten(sel_jet.delta_r(pruned_ev.SelMuon)),
                                weight=weight,
                            )
                        else:
                            h.fill(
                                syst,
                                flatten(flav),
                                flatten(sel_jet[histname.replace(f"jet{i}_", "")]),
                                weight=weight,
                            )
            # mu-Jets distribution
            elif "lmujet_" in histname:
                h.fill(
                    syst,
                    smflav,
                    flatten(pruned_ev.MuonJet[histname.replace("lmujet_", "")]),
                    weight=weight,
                )
            # filled discriminants
            elif "btag" in histname or "PNet" in histname:
                # Events with muon jet
                if "MuonJet" in pruned_ev.fields:
                    flavs, seljets = smflav, pruned_ev.MuonJet
                    nj = 1
                else:
                    flavs, seljets = genflavor, pruned_ev.SelJet

                for i in range(nj):
                    if not histname.endswith(str(i)):
                        continue
                    if nj > 1:
                        flav, seljet = flavs[:, i], seljets[:, i]
                    else:
                        flav, seljet = flavs, seljets
                    h.fill(
                        syst=syst,
                        flav=flav,
                        discr=seljet[histname.replace(f"_{i}", "")],
                        weight=weights.partial_weight(exclude=exclude_btv),
                    )

        if "dr_poslnegl" in output.keys():
            # DY histograms
            output["dr_poslnegl"].fill(
                syst, pruned_ev.posl.delta_r(pruned_ev.negl), weight=weight
            )
            output["dr_posljet"].fill(
                syst,
                genflavor,
                pruned_ev.posl.delta_r(pruned_ev.SelJet),
                weight=weight,
            )
            output["dr_negljet"].fill(
                syst,
                genflavor,
                pruned_ev.negl.delta_r(pruned_ev.SelJet),
                weight=weight,
            )
        # Muon enriched jet histograms
        if "MuonJet" in pruned_ev.fields:
            if (
                "hl" not in pruned_ev.fields
                and "SelElectron" in pruned_ev.fields
                and "SelMuon" in pruned_ev.fields
            ):
                pruned_ev["hl"] = ak.zip(
                    {
                        "pt": (
                            pruned_ev.SelMuon.pt
                            if pruned_ev.SelMuon.pt > pruned_ev.SelElectron.pt
                            else pruned_ev.SelElectron.pt
                        ),
                        "eta": (
                            pruned_ev.SelMuon.eta
                            if pruned_ev.SelMuon.pt > pruned_ev.SelElectron.pt
                            else pruned_ev.SelElectron.eta
                        ),
                        "phi": (
                            pruned_ev.SelMuon.phi
                            if pruned_ev.SelMuon.pt > pruned_ev.SelElectron.pt
                            else pruned_ev.SelElectron.phi
                        ),
                        "energy": (
                            pruned_ev.SelMuon.energy
                            if pruned_ev.SelMuon.pt > pruned_ev.SelElectron.pt
                            else pruned_ev.SelElectron.energy
                        ),
                    },
                    with_name="PtEtaPhiECandidate",
                )
            output["soft_l_ptratio"].fill(
                syst,
                flav=smflav,
                ratio=pruned_ev.SoftMuon.pt / pruned_ev.MuonJet.pt,
                weight=weight,
            )
            output["dr_lmujetsmu"].fill(
                syst,
                flav=smflav,
                dr=pruned_ev.MuonJet.delta_r(pruned_ev.SoftMuon),
                weight=weight,
            )
            if "hl" in pruned_ev.fields:
                output["hl_ptratio"].fill(
                    syst,
                    smflav,
                    ratio=pruned_ev.hl.pt / pruned_ev.MuonJet.pt,
                    weight=weight,
                )
            if "sl" in pruned_ev.fields:
                output["sl_ptratio"].fill(
                    syst,
                    smflav,
                    ratio=pruned_ev.sl.pt / pruned_ev.MuonJet.pt,
                    weight=weight,
                )

            if "SelMuon" in pruned_ev.fields and "hl" not in pruned_ev.fields:
                output["dr_lmujethmu"].fill(
                    syst,
                    flav=smflav,
                    dr=pruned_ev.MuonJet.delta_r(pruned_ev.SelMuon),
                    weight=weight,
                )
                output["dr_hmusmu"].fill(
                    syst, pruned_ev.SelMuon.delta_r(pruned_ev.SoftMuon), weight=weight
                )
        # dilepton system histograms: DY workflow
        if "dilep" in pruned_ev.fields:
            output["dilep_pt"].fill(syst, flatten(pruned_ev.dilep.pt), weight=weight)
            output["dilep_eta"].fill(syst, flatten(pruned_ev.dilep.eta), weight=weight)
            output["dilep_phi"].fill(syst, flatten(pruned_ev.dilep.phi), weight=weight)
            output["dilep_mass"].fill(
                syst, flatten(pruned_ev.dilep.mass), weight=weight
            )

        if "MET" in output.keys():
            output["MET_pt"].fill(syst, flatten(pruned_ev.MET.pt), weight=weight)
            output["MET_phi"].fill(syst, flatten(pruned_ev.MET.phi), weight=weight)

        # ttbar dilepton kin workflow
        if "kindisc" in output.keys():
            output["kindisc"].fill(
                syst,
                flatten(pruned_ev.flavour),
                flatten(pruned_ev.kindisc),
                weight=flatten(
                    ak.broadcast_arrays(
                        weights.partial_weight(exclude=exclude_btv), pruned_ev.kindisc
                    )[0]
                ),
            )
            output["close_mlj"].fill(
                syst,
                flatten(pruned_ev.flavour),
                flatten(pruned_ev.close_mlj),
                weight=flatten(ak.broadcast_arrays(weight, pruned_ev.close_mlj)[0]),
            )
            output["close_deta"].fill(
                syst,
                flatten(pruned_ev.flavour),
                flatten(pruned_ev.close_deta),
                weight=flatten(ak.broadcast_arrays(weight, pruned_ev.close_deta)[0]),
            )
            output["close_dphi"].fill(
                syst,
                flatten(pruned_ev.flavour),
                flatten(pruned_ev.close_dphi),
                weight=flatten(ak.broadcast_arrays(weight, pruned_ev.close_dphi)[0]),
            )
            output["close_ptrel"].fill(
                syst,
                flatten(pruned_ev.flavour),
                flatten(pruned_ev.close_ptrel),
                weight=flatten(ak.broadcast_arrays(weight, pruned_ev.close_ptrel)[0]),
            )
            output["close_lj2ll_deta"].fill(
                syst,
                flatten(pruned_ev.flavour),
                flatten(pruned_ev.close_lj2ll_deta),
                weight=flatten(
                    ak.broadcast_arrays(weight, pruned_ev.close_lj2ll_deta)[0]
                ),
            )
            output["close_lj2ll_dphi"].fill(
                syst,
                flatten(pruned_ev.flavour),
                flatten(pruned_ev.close_lj2ll_dphi),
                weight=flatten(
                    ak.broadcast_arrays(weight, pruned_ev.close_lj2ll_dphi)[0]
                ),
            )
            output["far_mlj"].fill(
                syst,
                flatten(pruned_ev.flavour),
                flatten(pruned_ev.far_mlj),
                weight=flatten(ak.broadcast_arrays(weight, pruned_ev.far_mlj)[0]),
            )
            output["far_deta"].fill(
                syst,
                flatten(pruned_ev.flavour),
                flatten(pruned_ev.far_deta),
                weight=flatten(ak.broadcast_arrays(weight, pruned_ev.far_deta)[0]),
            )
            output["far_dphi"].fill(
                syst,
                flatten(pruned_ev.flavour),
                flatten(pruned_ev.far_dphi),
                weight=flatten(ak.broadcast_arrays(weight, pruned_ev.far_dphi)[0]),
            )
            output["far_ptrel"].fill(
                syst,
                flatten(pruned_ev.flavour),
                flatten(pruned_ev.far_ptrel),
                weight=flatten(ak.broadcast_arrays(weight, pruned_ev.far_ptrel)[0]),
            )
            output["far_lj2ll_deta"].fill(
                syst,
                flatten(pruned_ev.flavour),
                flatten(pruned_ev.far_lj2ll_deta),
                weight=flatten(
                    ak.broadcast_arrays(weight, pruned_ev.far_lj2ll_deta)[0]
                ),
            )
            output["far_lj2ll_dphi"].fill(
                syst,
                flatten(pruned_ev.flavour),
                flatten(pruned_ev.far_lj2ll_dphi),
                weight=flatten(
                    ak.broadcast_arrays(weight, pruned_ev.far_lj2ll_dphi)[0]
                ),
            )
            output["j2ll_deta"].fill(
                syst,
                flatten(pruned_ev.flavour),
                flatten(pruned_ev.j2ll_deta),
                weight=flatten(ak.broadcast_arrays(weight, pruned_ev.j2ll_deta)[0]),
            )
            output["j2ll_dphi"].fill(
                syst,
                flatten(pruned_ev.flavour),
                flatten(pruned_ev.j2ll_dphi),
                weight=flatten(ak.broadcast_arrays(weight, pruned_ev.j2ll_dphi)[0]),
            )

    return output
