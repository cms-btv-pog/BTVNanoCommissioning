import hist as Hist


def get_histograms(axes, **kwargs):
    hists = {
        "dr_mujet0": Hist.Hist(
            axes["syst"], axes["flav"], axes["dr"], Hist.storage.Weight()
        )  # create customize histogram
    }
    return hists
