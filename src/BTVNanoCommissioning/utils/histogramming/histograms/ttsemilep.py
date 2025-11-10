import hist as Hist


def get_histograms(axes, **kwargs):
    hists = {}

    c_ttsemilep = kwargs.get("c_ttsemilep", None)
    if c_ttsemilep == None:
        raise ValueError(
            "c_ttsemilep is not specified when running ttsemilep workflow."
        )

    hists["dr_cjet"] = Hist.Hist(
        axes["syst"], axes["flav"], axes["dr"], Hist.storage.Weight()
    )
    if not c_ttsemilep:
        for i in range(4):
            hists[f"dr_mujet{i}"] = Hist.Hist(
                axes["syst"], axes["flav"], axes["dr"], Hist.storage.Weight()
            )
    for i in ["mu"]:
        hists[f"{i}_pfRelIso04_all"] = Hist.Hist(
            axes["syst"], axes["iso"], Hist.storage.Weight()
        )
        hists[f"{i}_dxy"] = Hist.Hist(axes["syst"], axes["dxy"], Hist.storage.Weight())
        hists[f"{i}_dz"] = Hist.Hist(axes["syst"], axes["dz"], Hist.storage.Weight())

    return hists
