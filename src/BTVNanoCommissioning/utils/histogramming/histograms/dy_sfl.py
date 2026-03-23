import hist as Hist

from BTVNanoCommissioning.utils.selection import btag_wp_dict


def get_histograms(axes, **kwargs):
    hists = {}

    include_m = kwargs.get("include_m", False)
    year = kwargs.get("year", "2024")
    campaign = kwargs.get("campaign", "Summer24")

    hists["dilep_mass"] = Hist.Hist(
        axes["syst"],
        Hist.axis.Regular(50, 50, 100, name="mass", label=" $m_{\\ell\\ell}$ [GeV]"),
        Hist.storage.Weight(),
    )
    hists["dr_poslnegl"] = Hist.Hist(axes["syst"], axes["dr"], Hist.storage.Weight())
    hists["dr_posljet"] = Hist.Hist(
        axes["syst"], axes["flav"], axes["dr"], Hist.storage.Weight()
    )
    hists["dr_negljet"] = Hist.Hist(
        axes["syst"], axes["flav"], axes["dr"], Hist.storage.Weight()
    )

    for i in ["posl", "negl"]:
        if include_m:
            hists[f"{i}_pfRelIso04_all"] = Hist.Hist(
                axes["syst"], axes["iso"], Hist.storage.Weight()
            )
        hists[f"{i}_dxy"] = Hist.Hist(axes["syst"], axes["dxy"], Hist.storage.Weight())
        hists[f"{i}_dz"] = Hist.Hist(axes["syst"], axes["dz"], Hist.storage.Weight())
        hists[f"dr_{i}jet"] = Hist.Hist(
            axes["syst"], axes["flav"], axes["dr"], Hist.storage.Weight()
        )
    for tagger in btag_wp_dict[year + "_" + campaign].keys():
        for stringency in btag_wp_dict[year + "_" + campaign][tagger]["b"].keys():
            if "No" in stringency:
                continue
            hists[f"{tagger}{stringency}_postag_jet_pt"] = Hist.Hist(
                axes["syst"], axes["flav"], axes["ljpt"], Hist.storage.Weight()
            )
            hists[f"{tagger}{stringency}_negtag_jet_pt"] = Hist.Hist(
                axes["syst"], axes["flav"], axes["ljpt"], Hist.storage.Weight()
            )
    hists["jet_pt"] = Hist.Hist(
        axes["syst"], axes["flav"], axes["ljpt"], Hist.storage.Weight()
    )
    return hists
