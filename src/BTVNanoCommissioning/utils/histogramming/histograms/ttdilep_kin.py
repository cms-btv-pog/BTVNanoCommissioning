import hist as Hist


def get_histograms(axes, **kwargs):
    hists = {}

    hists["kindisc"] = Hist.Hist(
        axes["syst"], axes["flav"], axes["bdt"], Hist.storage.Weight()
    )
    hists["close_mlj"] = Hist.Hist(
        axes["syst"], axes["flav"], axes["mass"], Hist.storage.Weight()
    )
    hists["close_deta"] = Hist.Hist(
        axes["syst"], axes["flav"], axes["eta"], Hist.storage.Weight()
    )
    hists["close_dphi"] = Hist.Hist(
        axes["syst"], axes["flav"], axes["phi"], Hist.storage.Weight()
    )
    hists["close_ptrel"] = Hist.Hist(
        axes["syst"], axes["flav"], axes["pt"], Hist.storage.Weight()
    )
    hists["close_lj2ll_deta"] = Hist.Hist(
        axes["syst"], axes["flav"], axes["eta"], Hist.storage.Weight()
    )
    hists["close_lj2ll_dphi"] = Hist.Hist(
        axes["syst"], axes["flav"], axes["phi"], Hist.storage.Weight()
    )

    hists["far_mlj"] = Hist.Hist(
        axes["syst"], axes["flav"], axes["mass"], Hist.storage.Weight()
    )
    hists["far_deta"] = Hist.Hist(
        axes["syst"], axes["flav"], axes["eta"], Hist.storage.Weight()
    )
    hists["far_dphi"] = Hist.Hist(
        axes["syst"], axes["flav"], axes["phi"], Hist.storage.Weight()
    )
    hists["far_ptrel"] = Hist.Hist(
        axes["syst"], axes["flav"], axes["pt"], Hist.storage.Weight()
    )
    hists["far_lj2ll_deta"] = Hist.Hist(
        axes["syst"], axes["flav"], axes["eta"], Hist.storage.Weight()
    )
    hists["far_lj2ll_dphi"] = Hist.Hist(
        axes["syst"], axes["flav"], axes["phi"], Hist.storage.Weight()
    )

    hists["j2ll_deta"] = Hist.Hist(
        axes["syst"], axes["flav"], axes["eta"], Hist.storage.Weight()
    )
    hists["j2ll_dphi"] = Hist.Hist(
        axes["syst"], axes["flav"], axes["phi"], Hist.storage.Weight()
    )

    return hists
