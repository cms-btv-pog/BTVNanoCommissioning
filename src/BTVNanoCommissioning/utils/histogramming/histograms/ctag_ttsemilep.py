import hist as Hist


def get_histograms(axes, **kwargs):
    hists = {}

    hists["dilep_mass"] = Hist.Hist(
        axes["syst"],
        Hist.axis.Regular(50, 50, 100, name="mass", label=" $m_{\\ell\\ell}$ [GeV]"),
        Hist.storage.Weight(),
    )
    # delta R between soft muon and mu-jet
    hists["dr_lmujetsmu"] = Hist.Hist(
        axes["syst"], axes["flav"], axes["dr_s"], Hist.storage.Weight()
    )
    # delta R between hard muon and mu-jet
    hists["dr_lmujethmu"] = Hist.Hist(
        axes["syst"], axes["flav"], axes["dr"], Hist.storage.Weight()
    )
    # delta R between hard muon and soft-muon
    hists["dr_hmusmu"] = Hist.Hist(axes["syst"], axes["dr"], Hist.storage.Weight())
    for i in ["hl", "soft_l"]:
        if i == "soft_l":
            hists[f"soft_l_pfRelIso04_all"] = Hist.Hist(
                axes["syst"], axes["flav"], axes["softliso"], Hist.storage.Weight()
            )
            hists[f"{i}_dxy"] = Hist.Hist(
                axes["syst"], axes["flav"], axes["dxy"], Hist.storage.Weight()
            )
            hists[f"{i}_dz"] = Hist.Hist(
                axes["syst"], axes["flav"], axes["dz"], Hist.storage.Weight()
            )
        else:
            hists[f"{i}_pfRelIso04_all"] = Hist.Hist(
                axes["syst"], axes["iso"], Hist.storage.Weight()
            )
            hists[f"{i}_dxy"] = Hist.Hist(
                axes["syst"], axes["dxy"], Hist.storage.Weight()
            )
            hists[f"{i}_dz"] = Hist.Hist(
                axes["syst"], axes["dz"], Hist.storage.Weight()
            )
        hists[f"{i}_ptratio"] = Hist.Hist(
            axes["syst"], axes["flav"], axes["ptratio"], Hist.storage.Weight()
        )

    return hists
