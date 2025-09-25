import hist as Hist

from BTVNanoCommissioning.utils.selection import btag_wp_dict


def get_histograms(axes, **kwargs):

    year = kwargs.get("year", None)
    if year == None:
        raise ValueError(
            "year is not specified. Please specify the year in histogrammer."
        )

    campaign = kwargs.get("campaign", None)
    if campaign == None:
        raise ValueError(
            "campaign is not specified. Please specify the campaign in histogrammer."
        )

    cutbased = kwargs.get("cutbased", None)
    if cutbased == None:
        raise ValueError(
            "cutbased is not specified. Please specify whether to use a cutbased wf or not."
        )

    obj_list = kwargs.get("obj_list", [])
    if not obj_list:
        raise ValueError(
            "obj_list is empty. Please specify the objects to be histogrammed."
        )

    hists = {}

    hists["SV_charge"] = Hist.Hist(
        axes["syst"],
        axes["osss"],
        Hist.axis.Regular(20, -10, 10, name="charge", label="SV charge"),
        Hist.storage.Weight(),
    )
    hists["dilep_mass"] = Hist.Hist(
        axes["syst"],
        axes["osss"],
        Hist.axis.Regular(50, 50, 100, name="mass", label=" $m_{\\ell\\ell}$ [GeV]"),
        Hist.storage.Weight(),
    )
    hists["w_mass"] = Hist.Hist(
        axes["syst"],
        axes["osss"],
        Hist.axis.Regular(50, 50, 100, name="mass", label=" $m_{\\ell\\nu}$ [GeV]"),
        Hist.storage.Weight(),
    )
    hists["w_pt"] = Hist.Hist(
        axes["syst"],
        axes["osss"],
        Hist.axis.Regular(50, 0, 300, name="pT", label=" $p_T^{\\ell\\nu}$ [GeV]"),
        Hist.storage.Weight(),
    )
    hists["w_phi"] = Hist.Hist(
        axes["syst"],
        axes["osss"],
        Hist.axis.Regular(50, 0, 360, name="phi", label=" $\\phi_{\\ell\\nu}$ [deg]"),
        Hist.storage.Weight(),
    )
    hists["dr_lmujetsmu"] = Hist.Hist(
        axes["syst"], axes["flav"], axes["osss"], axes["dr_s"], Hist.storage.Weight()
    )
    hists["dr_lmujethmu"] = Hist.Hist(
        axes["syst"], axes["flav"], axes["osss"], axes["dr"], Hist.storage.Weight()
    )
    hists["dr_hmusmu"] = Hist.Hist(
        axes["syst"], axes["osss"], axes["dr"], Hist.storage.Weight()
    )
    for i in ["hl", "sl", "soft_l"]:
        if i == "soft_l":
            hists[f"soft_l_pfRelIso04_all"] = Hist.Hist(
                axes["syst"],
                axes["flav"],
                axes["osss"],
                axes["softliso"],
                Hist.storage.Weight(),
            )
            hists[f"{i}_dxy"] = Hist.Hist(
                axes["syst"],
                axes["flav"],
                axes["osss"],
                axes["dxy"],
                Hist.storage.Weight(),
            )
            hists[f"{i}_dz"] = Hist.Hist(
                axes["syst"],
                axes["flav"],
                axes["osss"],
                axes["dz"],
                Hist.storage.Weight(),
            )
        else:
            hists[f"{i}_pfRelIso04_all"] = Hist.Hist(
                axes["syst"], axes["osss"], axes["iso"], Hist.storage.Weight()
            )
            hists[f"{i}_dxy"] = Hist.Hist(
                axes["syst"], axes["osss"], axes["qcddxy"], Hist.storage.Weight()
            )
            hists[f"{i}_dz"] = Hist.Hist(
                axes["syst"], axes["osss"], axes["dz"], Hist.storage.Weight()
            )

        hists[f"{i}_ptratio"] = Hist.Hist(
            axes["syst"],
            axes["flav"],
            axes["osss"],
            axes["ptratio"],
            Hist.storage.Weight(),
        )

    for obj in obj_list:
        if "mujet" in obj and cutbased:
            for tagger in btag_wp_dict[year + "_" + campaign].keys():
                for wp in btag_wp_dict[year + "_" + campaign][tagger]["c"].keys():
                    if not "No" in wp:
                        hists[f"{obj}_pt_{tagger}{wp}"] = Hist.Hist(
                            axes["syst"],
                            axes["flav"],
                            axes["osss"],
                            axes["pt"],
                            Hist.storage.Weight(),
                        )

    return hists
