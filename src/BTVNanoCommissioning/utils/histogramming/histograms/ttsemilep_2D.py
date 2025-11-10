import hist as Hist


def get_histograms(axes, **kwargs):
    hists = {}

    hists[f"nbjet"] = Hist.Hist(axes["syst"], axes["n"], Hist.storage.Weight())
    hists[f"ncjet"] = Hist.Hist(axes["syst"], axes["n"], Hist.storage.Weight())
    hists[f"w_mt"] = Hist.Hist(axes["syst"], axes["mt"], Hist.storage.Weight())

    channel = kwargs.get("channel", "mu")
    hists[f"{channel}_pfRelIso04_all"] = Hist.Hist(axes["syst"], axes["iso"], Hist.storage.Weight())
    hists[f"{channel}_dxy"] = Hist.Hist(axes["syst"], axes["dxy"], Hist.storage.Weight())
    hists[f"{channel}_dz"] = Hist.Hist(axes["syst"], axes["dz"], Hist.storage.Weight())

    return hists
