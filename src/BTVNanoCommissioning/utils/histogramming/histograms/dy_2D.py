import hist as Hist


def get_histograms(axes, **kwargs):
    hists = {}

    hists["dilep_ptratio"] = Hist.Hist(
        axes["syst"], axes["flav"], axes["ptratio"], Hist.storage.Weight()
    )
    hists["top_pt"] = Hist.Hist(axes["syst"], axes["jpt"], Hist.storage.Weight())
    hists["antitop_pt"] = Hist.Hist(axes["syst"], axes["jpt"], Hist.storage.Weight())

    return hists
