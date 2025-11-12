import hist as Hist


def get_histograms(axes, **kwargs):
    hists = {}

    hists["dilep_ptratio"] = Hist.Hist(
        axes["syst"], axes["flav"], axes["ptratio"], Hist.storage.Weight()
    )

    return hists
