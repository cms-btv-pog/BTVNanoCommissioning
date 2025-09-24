import hist as Hist

def get_histograms(axes, **kwargs):
    hists = {}

    hists["dr_cjet"] = Hist.Hist(
        axes["syst"], axes["flav"], axes["dr"], Hist.storage.Weight()
    )

    # for i in range(4):
    #     hists[f"dr_mujet{i}"] = Hist.Hist(
    #         axes["syst"], axes["flav"], axes["dr"], Hist.storage.Weight()
    #     )

    for i in ["mu"]:
        hists[f"{i}_pfRelIso04_all"] = Hist.Hist(
            axes["syst"], axes["iso"], Hist.storage.Weight()
        )
        hists[f"{i}_dxy"] = hist.Hist(
            axes["syst"], axes["dxy"], Hist.storage.Weight()
        )
        hists[f"{i}_dz"] = hist.Hist(
            axes["syst"], axes["dz"], Hist.storage.Weight()
        )

    return hists