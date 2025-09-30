import hist as Hist


def get_histograms(axes, **kwargs):
    hists = {}

    hists["bjet_WP_pt"] = Hist.Hist(
        Hist.axis.StrCategory([], name="WP", growth=True),
        Hist.axis.StrCategory(
            ["DeepFlav", "PNet", "RobustParTAK4"], name="tagger", growth=True
        ),
        axes["pt"],
        Hist.storage.Weight(),
    )
    hists["bjet_WP_eta"] = Hist.Hist(
        Hist.axis.StrCategory([], name="WP", growth=True),
        Hist.axis.StrCategory(
            ["DeepFlav", "PNet", "RobustParTAK4"], name="tagger", growth=True
        ),
        axes["eta"],
        Hist.storage.Weight(),
    )
    hists["bjet_WP_phi"] = Hist.Hist(
        Hist.axis.StrCategory([], name="WP", growth=True),
        Hist.axis.StrCategory(
            ["DeepFlav", "PNet", "RobustParTAK4"], name="tagger", growth=True
        ),
        axes["phi"],
        Hist.storage.Weight(),
    )
    hists["bjet_WP_discr"] = Hist.Hist(
        Hist.axis.StrCategory([], name="WP", growth=True),
        Hist.axis.StrCategory(
            ["DeepFlav", "PNet", "RobustParTAK4"], name="tagger", growth=True
        ),
        Hist.axis.Regular(25, 0, 1, name="B"),
        Hist.storage.Weight(),
    )
    hists["cjet_WP_pt"] = Hist.Hist(
        Hist.axis.StrCategory([], name="WP", growth=True),
        Hist.axis.StrCategory(
            ["DeepFlav", "PNet", "RobustParTAK4"], name="tagger", growth=True
        ),
        axes["pt"],
        Hist.storage.Weight(),
    )
    hists["cjet_WP_eta"] = Hist.Hist(
        Hist.axis.StrCategory([], name="WP", growth=True),
        Hist.axis.StrCategory(
            ["DeepFlav", "PNet", "RobustParTAK4"], name="tagger", growth=True
        ),
        axes["eta"],
        Hist.storage.Weight(),
    )
    hists["cjet_WP_phi"] = Hist.Hist(
        Hist.axis.StrCategory([], name="WP", growth=True),
        Hist.axis.StrCategory(
            ["DeepFlav", "PNet", "RobustParTAK4"], name="tagger", growth=True
        ),
        axes["phi"],
        Hist.storage.Weight(),
    )
    hists["cjet_WP_discr"] = Hist.Hist(
        Hist.axis.StrCategory([], name="WP", growth=True),
        Hist.axis.StrCategory(
            ["DeepFlav", "PNet", "RobustParTAK4"], name="tagger", growth=True
        ),
        Hist.axis.Regular(25, 0, 1, name="CvL"),
        Hist.axis.Regular(25, 0, 1, name="CvB"),
        Hist.storage.Weight(),
    )

    return hists
