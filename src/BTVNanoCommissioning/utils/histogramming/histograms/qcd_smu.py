import hist as Hist


def get_histograms(axes, **kwargs):
    hists = {}

    hists["nmujet"] = Hist.Hist(axes["syst"], axes["n"], Hist.storage.Weight())
    hists["nssmu"] = Hist.Hist(axes["syst"], axes["n"], Hist.storage.Weight())
    hists["dr_lmujetsmu"] = Hist.Hist(
        axes["syst"], axes["flav"], axes["dr_s"], Hist.storage.Weight()
    )
    for i in ["soft_l"]:
        if i == "soft_l":
            hists[f"soft_l_pfRelIso04_all"] = Hist.Hist(
                axes["syst"],
                axes["flav"],
                axes["softliso"],
                Hist.storage.Weight(),
            )
            hists[f"{i}_dxy"] = Hist.Hist(
                axes["syst"], axes["flav"], axes["dxy"], Hist.storage.Weight()
            )
            hists[f"{i}_dz"] = Hist.Hist(
                axes["syst"], axes["flav"], axes["dz"], Hist.storage.Weight()
            )
        hists[f"{i}_ptratio"] = Hist.Hist(
            axes["syst"], axes["flav"], axes["ptratio"], Hist.storage.Weight()
        )
    hists["dr_lmusmujetsmu"] = Hist.Hist(
        axes["syst"], axes["flav"], axes["dr_s"], Hist.storage.Weight()
    )

    return hists
