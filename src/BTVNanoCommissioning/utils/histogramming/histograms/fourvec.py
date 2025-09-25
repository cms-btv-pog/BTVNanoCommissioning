import hist as Hist


def get_histograms(axes, **kwargs):
    hists = {}

    obj_list = kwargs.get("obj_list", [])
    if len(obj_list) == 0:
        raise ValueError("obj_list is not specified when running fourvec workflow.")
    include_osss = kwargs.get("include_osss", False)

    if not include_osss:
        pt_axes = [axes["syst"], axes["flav"], axes["pt"]]
        pt_axes_no_flav = [axes["syst"], axes["pt"]]
        soft_l_pt_axes = [axes["syst"], axes["flav"], axes["softlpt"]]
        eta_axes = [axes["syst"], axes["flav"], axes["eta"]]
        eta_axes_no_flav = [axes["syst"], axes["eta"]]
        phi_axes = [axes["syst"], axes["flav"], axes["phi"]]
        phi_axes_no_flav = [axes["syst"], axes["phi"]]
        mass_axes = [axes["syst"], axes["flav"], axes["mass"]]
    else:
        pt_axes = [axes["syst"], axes["flav"], axes["osss"], axes["pt"]]
        pt_axes_no_flav = [axes["syst"], axes["osss"], axes["pt"]]
        soft_l_pt_axes = [axes["syst"], axes["flav"], axes["osss"], axes["softlpt"]]
        eta_axes = [axes["syst"], axes["flav"], axes["osss"], axes["eta"]]
        eta_axes_no_flav = [axes["syst"], axes["osss"], axes["eta"]]
        phi_axes = [axes["syst"], axes["flav"], axes["osss"], axes["phi"]]
        phi_axes_no_flav = [axes["syst"], axes["osss"], axes["phi"]]
        mass_axes = [axes["syst"], axes["flav"], axes["osss"], axes["mass"]]

    for obj in obj_list:
        if "jet" in obj or "soft_l" in obj:
            if obj == "soft_l":
                hists["softlpt"] = Hist.Hist(*soft_l_pt_axes, Hist.storage.Weight())
            else:
                hists[f"{obj}_pt"] = Hist.Hist(*pt_axes, Hist.storage.Weight())
            hists[f"{obj}_eta"] = Hist.Hist(*eta_axes, Hist.storage.Weight())
            hists[f"{obj}_phi"] = Hist.Hist(*phi_axes, Hist.storage.Weight())
            hists[f"{obj}_mass"] = Hist.Hist(*mass_axes, Hist.storage.Weight())

        else:
            hists[f"{obj}_pt"] = Hist.Hist(*pt_axes_no_flav, Hist.storage.Weight())
            hists[f"{obj}_phi"] = Hist.Hist(*phi_axes_no_flav, Hist.storage.Weight())
            if obj != "MET":
                hists[f"{obj}_eta"] = Hist.Hist(
                    *eta_axes_no_flav, Hist.storage.Weight()
                )

    return hists
