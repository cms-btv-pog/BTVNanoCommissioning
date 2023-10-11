import matplotlib.pyplot as plt
import numpy as np
import scipy.stats
import warnings
import hist

errband_opts = {
    "hatch": "////",
    "facecolor": "none",
    "lw": 0,
    "color": "k",
    "alpha": 0.4,
}
markers = [".", "o", "^", "s", "+", "x", "D", "*"]
color_map = {
    "$t\\bar{t}$": "#008695",
    "Single top": "#3969ac",
    "VV": "#7f3c8d",
    "Z+jets": "#e68310",
    "W+jets": "#80ba5a",
    "QCD": "#e73f74",
    "QCD($\\mu$)": "#a5aa99",
    "udsg": "tab:blue",
    "pu": "tab:orange",
    "c": "tab:green",
    "b": "tab:red",
}
sample_mergemap = {
    # ttbar
    "$t\\bar{t}$": [
        "TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
        "TTto4Q_TuneCP5_13p6TeV_powheg-pythia8",
        "TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",
        "TTTo2J1L1Nu_CP5_13p6TeV_powheg-pythia8",
    ],
    # single top
    "Single top": [
        "TbarBQ_t-channel_4FS_CP5_13p6TeV_powheg-madspin-pythia8",
        "TBbarQ_t-channel_4FS_CP5_13p6TeV_powheg-madspin-pythia8",
        "TWminus_DR_AtLeastOneLepton_CP5_13p6TeV_powheg-pythia8",
        "TbarBQ_t-channel_4FS_CP5_13p6TeV_powheg-madspin-pythia8",
        "TbarWplus_DR_AtLeastOneLepton_CP5_13p6TeV_powheg-pythia8",
        # decay
        "TWminusto4Q_TuneCP5_13p6TeV_powheg-pythia8",
        "TWminustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",
        "TWminusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
        "TbarWplusto4Q_TuneCP5_13p6TeV_powheg-pythia8",
        "TbarWplustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",
        "TbarWplusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
        "TBbarQ_t-channel_4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8",
        "TbarBQ_t-channel_4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8",
    ],
    # diboson
    "VV": [
        "WW_TuneCP5_13p6TeV-pythia8",
        "WZ_TuneCP5_13p6TeV-pythia8",
        "ZZ_TuneCP5_13p6TeV-pythia8"
        # decay
        "ZZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8",
        "ZZto2Nu2Q_TuneCP5_13p6TeV_powheg-pythia8",
        "ZZto4L_TuneCP5_13p6TeV_powheg-pythia8",
        "ZZto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
        "WZtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",
        "WZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8",
        "WZto3LNu_TuneCP5_13p6TeV_powheg-pythia8",
        "WWto4Q_TuneCP5_13p6TeV_powheg-pythia8",
        "WWtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",
        "WWto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
    ],
    "Z+jets": [
        "DYto2L-2Jets_MLL-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "DYto2L-2Jets_MLL-50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "DYJetsToLL_M-10to50_TuneCP5_13p6TeV-madgraphMLM-pythia8",
        "DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8",
    ],
    "W+jets": [
        "WJetsToLNu_TuneCP5_13p6TeV-madgraphMLM-pythia8",
        "WtoLNu-4Jets_TuneCP5_13p6TeV_madgraphMLM-pythia8",
        # jet binned sample
        "WtoLNu-4Jets_1J_TuneCP5_13p6TeV_madgraphMLM-pythia8",
        "WtoLNu-4Jets_2J_TuneCP5_13p6TeV_madgraphMLM-pythia8",
        "WtoLNu-4Jets_3J_TuneCP5_13p6TeV_madgraphMLM-pythia8",
    ],
    # QCD
    "QCD": [
        "QCD_PT-15to30_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-30to50_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-50to80_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-50to80_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-80to120_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-120to170_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-170to300_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-170to300_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-300to470_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-300to470_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-470to600_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-470to600_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-600to800_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-600to800_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-800to1000_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-800to1000_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-1000to1400_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-1000to1400_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-1400to1800_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-1400to1800_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-1800to2400_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-1800to2400_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-2400to3200_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-2400to3200_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-3200_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-3200_TuneCP5_13p6TeV_pythia8",
    ],
    # QCD muon enriched
    "QCD($\\mu$)": [
        "QCD_PT-15to20_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-20to30_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-30to50_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-50to80_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-80to120_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-120to170_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-170to300_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-300to470_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-470to600_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-600to800_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-800to1000_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
        "QCD_PT-1000_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8",
    ],
}
### copy functions coffea.hist.plotratio https://github.com/CoffeaTeam/coffea/blob/master/coffea/hist/plot.py to boost-hist
################
## ratio plot ##
################

_coverage1sd = scipy.stats.norm.cdf(1) - scipy.stats.norm.cdf(-1)


def compatible(self, other):
    """Checks if this histogram is compatible with another, i.e. they have identical binning"""
    if len(self.axes) != len(other.axes):
        return False
    if set(self.axes.name) != set(other.axes.name):
        return False
    # if not all(d1 == d2 for d1, d2 in zip(self.dense_axes(), other.dense_axes())):
    #     return False
    return True


def poisson_interval(sumw, sumw2, coverage=_coverage1sd):
    """Frequentist coverage interval for Poisson-distributed observations
    Parameters
    ----------
        sumw : numpy.ndarray
            Sum of weights vector
        sumw2 : numpy.ndarray
            Sum weights squared vector
        coverage : float, optional
            Central coverage interval, defaults to 68%
    Calculates the so-called 'Garwood' interval,
    c.f. https://www.ine.pt/revstat/pdf/rs120203.pdf or
    http://ms.mcmaster.ca/peter/s743/poissonalpha.html
    For weighted data, this approximates the observed count by ``sumw**2/sumw2``, which
    effectively scales the unweighted poisson interval by the average weight.
    This may not be the optimal solution: see https://arxiv.org/pdf/1309.1287.pdf for a proper treatment.
    When a bin is zero, the scale of the nearest nonzero bin is substituted to scale the nominal upper bound.
    If all bins zero, a warning is generated and interval is set to ``sumw``.
    """
    scale = np.empty_like(sumw)
    scale[sumw != 0] = sumw2[sumw != 0] / sumw[sumw != 0]
    if np.sum(sumw == 0) > 0:
        missing = np.where(sumw == 0)
        available = np.nonzero(sumw)
        if len(available[0]) == 0:
            warnings.warn(
                "All sumw are zero!  Cannot compute meaningful error bars",
                RuntimeWarning,
            )
            return np.vstack([sumw, sumw])
        nearest = sum(
            [np.subtract.outer(d, d0) ** 2 for d, d0 in zip(available, missing)]
        ).argmin(axis=0)
        argnearest = tuple(dim[nearest] for dim in available)
        scale[missing] = scale[argnearest]
    counts = sumw / scale
    lo = scale * scipy.stats.chi2.ppf((1 - coverage) / 2, 2 * counts) / 2.0
    hi = scale * scipy.stats.chi2.ppf((1 + coverage) / 2, 2 * (counts + 1)) / 2.0
    interval = np.array([lo, hi])
    interval[interval == np.nan] = 0.0  # chi2.ppf produces nan for counts=0
    return interval


def normal_interval(pw, tw, pw2, tw2, coverage=_coverage1sd):
    """Compute errors based on the expansion of pass/(pass + fail), possibly weighted
    Parameters
    ----------
    pw : np.ndarray
        Numerator, or number of (weighted) successes, vectorized
    tw : np.ndarray
        Denominator or number of (weighted) trials, vectorized
    pw2 : np.ndarray
        Numerator sum of weights squared, vectorized
    tw2 : np.ndarray
        Denominator sum of weights squared, vectorized
    coverage : float, optional
        Central coverage interval, defaults to 68%
    c.f. https://root.cern.ch/doc/master/TEfficiency_8cxx_source.html#l02515
    """

    eff = pw / tw

    variance = (pw2 * (1 - 2 * eff) + tw2 * eff**2) / (tw**2)
    sigma = np.sqrt(variance)

    prob = 0.5 * (1 - coverage)
    delta = np.zeros_like(sigma)
    delta[sigma != 0] = scipy.stats.norm.ppf(prob, scale=sigma[sigma != 0])

    lo = eff - np.minimum(eff + delta, np.ones_like(eff))
    hi = np.maximum(eff - delta, np.zeros_like(eff)) - eff

    return np.array([lo, hi])


def clopper_pearson_interval(num, denom, coverage=_coverage1sd):
    """Compute Clopper-Pearson coverage interval for a binomial distribution
    Parameters
    ----------
        num : numpy.ndarray
            Numerator, or number of successes, vectorized
        denom : numpy.ndarray
            Denominator or number of trials, vectorized
        coverage : float, optional
            Central coverage interval, defaults to 68%
    c.f. http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    """
    if np.any(num > denom):
        raise ValueError(
            "Found numerator larger than denominator while calculating binomial uncertainty"
        )
    lo = scipy.stats.beta.ppf((1 - coverage) / 2, num, denom - num + 1)
    hi = scipy.stats.beta.ppf((1 + coverage) / 2, num + 1, denom - num)
    interval = np.array([lo, hi])
    interval[:, num == 0.0] = 0.0
    interval[1, num == denom] = 1.0
    return interval


## ratioplot function
def plotratio(
    num,
    denom,
    ax=None,
    clear=True,
    flow=None,
    xerr=False,
    error_opts={},
    denom_fill_opts={},
    guide_opts={},
    unc="num",
    label=None,
    ext_denom_error=None,
):
    """Create a ratio plot, dividing two compatible histograms
    Parameters
    ----------
        num : Hist
            Numerator, a single-axis histogram
        denom : Hist
            Denominator, a single-axis histogram
        ax : matplotlib.axes.Axes, optional
            Axes object (if None, one is created)
        clear : bool, optional
            Whether to clear Axes before drawing (if passed); if False, this function will skip drawing the legend
        flow : str, optional {None, "show", "sum"}
            Whether plot the under/overflow bin. If "show", add additional under/overflow bin. If "sum", add the under/overflow bin content to first/last bin.
        xerr: bool, optional
            If true, then error bars are drawn for x-axis to indicate the size of the bin.
        error_opts : dict, optional
            A dictionary of options to pass to the matplotlib
            `ax.errorbar <https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.errorbar.html>`_ call
            internal to this function.  Leave blank for defaults.  Some special options are interpreted by
            this function and not passed to matplotlib: 'emarker' (default: '') specifies the marker type
            to place at cap of the errorbar.
        denom_fill_opts : dict, optional
            A dictionary of options to pass to the matplotlib
            `ax.fill_between <https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.axes.Axes.fill_between.html>`_ call
            internal to this function, filling the denominator uncertainty band.  Leave blank for defaults.
        guide_opts : dict, optional
            A dictionary of options to pass to the matplotlib
            `ax.axhline <https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.axes.Axes.axhline.html>`_ call
            internal to this function, to plot a horizontal guide line at ratio of 1.  Leave blank for defaults.
        unc : str, optional
            Uncertainty calculation option: 'clopper-pearson' interval for efficiencies; 'poisson-ratio' interval
            for ratio of poisson distributions; 'num' poisson interval of numerator scaled by denominator value
            (common for data/mc, for better or worse).
        label : str, optional
            Associate a label to this entry (note: y axis label set by ``num.label``)
        ext_denom_error: list of np.array[error_up,error_down], optional
            External MC errors not stored in the original histogram
    Returns
    -------
        ax : matplotlib.axes.Axes
            A matplotlib `Axes <https://matplotlib.org/3.1.1/api/axes_api.html>`_ object
    """
    if ax is None:
        fig, ax = plt.subplots(1, 1)
    else:
        if not isinstance(ax, plt.Axes):
            raise ValueError("ax must be a matplotlib Axes object")
        if clear:
            ax.clear()
    if not compatible(num, denom):
        raise ValueError(
            "numerator and denominator histograms have incompatible axis definitions"
        )
    if len(num.axes) > 1:
        raise ValueError("plotratio() can only support one-dimensional histograms")
    if error_opts is None and denom_fill_opts is None and guide_opts is None:
        error_opts = {}
        denom_fill_opts = {}

    axis = num.axes[0]

    ax.set_xlabel(axis.label)
    ax.set_ylabel(num.label)
    edges = axis.edges
    if flow == "show":
        edges = np.array(
            [
                edges[0] - (edges[1] - edges[0]) * 3,
                edges[0] - (edges[1] - edges[0]),
                *edges,
                edges[-1] + (edges[1] - edges[0]),
                edges[-1] + (edges[1] - edges[0]) * 3,
            ]
        )
    centers = (edges[1:] + edges[:-1]) / 2
    ranges = (edges[1:] - edges[:-1]) / 2 if xerr else None
    sumw_num, sumw2_num = num.values(), num.variances()
    sumw_denom, sumw2_denom = denom.values(), denom.variances()

    if flow == "show":
        print("Show under/overflow bin ")
        sumw_num, sumw2_num = (
            num.view(flow=True)["value"],
            num.view(flow=True)["variance"],
        )
        sumw_num, sumw2_num = np.insert(sumw_num, -1, 0), np.insert(sumw2_num, -1, 0)
        sumw_num, sumw2_num = np.insert(sumw_num, 1, 0), np.insert(sumw2_num, 1, 0)
        sumw_denom, sumw2_denom = (
            denom.view(flow=True)["value"],
            denom.view(flow=True)["variance"],
        )
        sumw_denom, sumw2_denom = np.insert(sumw_denom, -1, 0), np.insert(
            sumw2_denom, -1, 0
        )
        sumw_denom, sumw2_denom = np.insert(sumw_denom, 1, 0), np.insert(
            sumw2_denom, 1, 0
        )
    elif flow == "sum":
        print("Merge under/overflow bin to first/last bin")

        sumw_num[0], sumw2_num[0] = (
            sumw_num[0] + num.view(flow=True)["value"][0],
            sumw2_num[0] + num.view(flow=True)["value"][0],
        )
        sumw_num[-1], sumw2_num[-1] = (
            sumw_num[-1] + num.view(flow=True)["value"][-1],
            sumw2_num[-1] + num.view(flow=True)["value"][-1],
        )
        sumw_denom[0], sumw2_denom[0] = (
            sumw_denom[0] + denom.view(flow=True)["value"][0],
            sumw2_denom[0] + denom.view(flow=True)["value"][0],
        )
        sumw_denom[-1], sumw2_denom[-1] = (
            sumw_denom[-1] + denom.view(flow=True)["value"][-1],
            sumw2_denom[-1] + denom.view(flow=True)["value"][-1],
        )
    else:
        sumw_num, sumw2_num = num.values(), num.variances()
        sumw_denom, sumw2_denom = denom.values(), denom.variances()
    rsumw = sumw_num / sumw_denom
    if unc == "clopper-pearson":
        rsumw_err = np.abs(clopper_pearson_interval(sumw_num, sumw_denom) - rsumw)
    elif unc == "poisson-ratio":
        # poisson ratio n/m is equivalent to binomial n/(n+m)
        rsumw_err = np.abs(
            clopper_pearson_interval(sumw_num, sumw_num + sumw_denom) - rsumw
        )
    elif unc == "num":
        rsumw_err = np.abs(poisson_interval(rsumw, sumw2_num / sumw_denom**2) - rsumw)
    elif unc == "efficiency":
        rsumw_err = np.abs(
            normal_interval(sumw_num, sumw_denom, sumw2_num, sumw2_denom)
        )
    else:
        raise ValueError("Unrecognized uncertainty option: %r" % unc)

    ## if additional uncertainties
    if ext_denom_error is not None:
        if denom_fill_opts is {}:
            print("suggest to use different style for additional error")
        if np.shape(rsumw_err) != np.shape(ext_denom_error / sumw_denom):
            raise ValueError("Imcompatible error length")
        rsumw_err = np.sqrt(rsumw_err**2 + (ext_denom_error / sumw_denom) ** 2)

    if error_opts is not None:
        opts = {
            "label": label,
            "linestyle": "none",
            "lw": 1,
            "marker": "o",
            "color": "k",
        }
        opts.update(error_opts)
        emarker = opts.pop("emarker", "")
        errbar = ax.errorbar(x=centers, y=rsumw, xerr=ranges, yerr=rsumw_err, **opts)
        plt.setp(errbar[1], "marker", emarker)
    if denom_fill_opts is not None:
        unity = np.ones_like(sumw_denom)
        denom_unc = poisson_interval(unity, sumw2_denom / sumw_denom**2)
        opts = {
            "hatch": "////",
            "facecolor": "none",
            "lw": 0,
            "color": "k",
            "alpha": 0.4,
        }
        if ext_denom_error is not None:
            denom_unc[0] = (
                unity[0]
                - np.sqrt(
                    (denom_unc - unity) ** 2 + (ext_denom_error / sumw_denom) ** 2
                )[0]
            )
            denom_unc[1] = (
                unity[1]
                + np.sqrt(
                    (denom_unc - unity) ** 2 + (ext_denom_error / sumw_denom) ** 2
                )[1]
            )
            opts = denom_fill_opts
        ax.stairs(denom_unc[0], edges=edges, baseline=denom_unc[1], **opts)
    if guide_opts is not None:
        opts = {"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1}
        opts.update(guide_opts)
        if clear is not False:
            ax.axhline(1.0, **opts)

    if clear:
        ax.autoscale(axis="x", tight=True)
        ax.set_ylim(0, None)

    return ax


## Calculate SF error
def SFerror(collated, discr, flow=None):
    err_up = (
        collated["mc"][discr][{"syst": "SFup", "flav": sum}].values()
        - collated["mc"][discr][{"syst": "SF", "flav": sum}].values()
    )
    err_dn = (
        collated["mc"][discr][{"syst": "SFdn", "flav": sum}].values()
        - collated["mc"][discr][{"syst": "SF", "flav": sum}].values()
    )
    if flow is not None:
        err_up = (
            collated["mc"][discr][{"syst": "SFup", "flav": sum}].view(flow=True)[
                "value"
            ]
            - collated["mc"][discr][{"syst": "SF", "flav": sum}].view(flow=True)[
                "value"
            ]
        )
        err_dn = (
            collated["mc"][discr][{"syst": "SFdn", "flav": sum}].view(flow=True)[
                "value"
            ]
            - collated["mc"][discr][{"syst": "SF", "flav": sum}].view(flow=True)[
                "value"
            ]
        )

    if "C" in discr:  ## scale uncertainties for charm tagger by 2
        err_up = np.sqrt(
            np.add(
                np.power(err_up, 2),
                np.power(
                    2
                    * (
                        collated["mc"][discr][{"syst": "SFup", "flav": 4}].values()
                        - collated["mc"][discr][{"syst": "SF", "flav": 4}].values()
                    ),
                    2,
                ),
            )
        )
        err_dn = np.sqrt(
            np.add(
                np.power(err_dn, 2),
                np.power(
                    2
                    * (
                        collated["mc"][discr][{"syst": "SFdn", "flav": 4}].values()
                        - collated["mc"][discr][{"syst": "SF", "flav": 4}].values()
                    ),
                    2,
                ),
            )
        )
        if flow:
            err_up = np.sqrt(
                np.add(
                    np.power(err_up, 2),
                    np.power(
                        2
                        * (
                            collated["mc"][discr][{"syst": "SFup", "flav": 4}].view(
                                flow=True
                            )["value"]
                            - collated["mc"][discr][{"syst": "SF", "flav": 4}].view(
                                flow=True
                            )["value"]
                        ),
                        2,
                    ),
                )
            )
            err_dn = np.sqrt(
                np.add(
                    np.power(err_dn, 2),
                    np.power(
                        2
                        * (
                            collated["mc"][discr][{"syst": "SFdn", "flav": 4}].view(
                                flow=True
                            )["value"]
                            - collated["mc"][discr][{"syst": "SF", "flav": 4}].view(
                                flow=True
                            )["value"]
                        ),
                        2,
                    ),
                )
            )

    return np.array([err_dn, err_up])


def autoranger(hist, flow=None):
    val, axis = hist.values(), hist.axes[-1].edges
    if flow == "show":
        val = hist.view(flow=True)["value"]
        axis = np.array(
            [
                axis[0] - (axis[1] - axis[0]) * 3,
                axis[0] - (axis[1] - axis[0]),
                *axis,
                axis[-1] + (axis[1] - axis[0]),
                axis[-1] + (axis[1] - axis[0]) * 3,
            ]
        )
    for i in range(len(val)):
        if val[i] != 0:
            mins = i
            break
    for i in reversed(range(len(val))):
        if val[i] != 0:
            maxs = i + 1
            break
    return axis[mins], axis[maxs]


def MCerrorband(
    hmc,
    ax=None,
    flow=None,
    label="Stat. unc.",
    fill_opts=None,
    ext_error=None,
    clear=False,
):
    """Create a ratio plot, dividing two compatible histograms
    Parameters
    ----------
        hmc : Hist
            A single-axis histogram
        ax : matplotlib.axes.Axes, optional
            Axes object (if None, one is created)
        flow : str, optional {None, "show", "sum"}
            Whether plot the under/overflow bin. If "show", add additional under/overflow bin. If "sum", add the under/overflow bin content to first/last bin.
        fill_opts : dict, optional
            A dictionary of options to pass to the matplotlib
            `ax.fill_between <https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.axes.Axes.fill_between.html>`_ call
            internal to this function, filling the denominator uncertainty band.  Leave blank for defaults.
        label : str, optional
            Associate a label to this entry (note: y axis label set by ``num.label``)
        ext_error: list of np.array[error_up,error_down], optional
            External MC errors not stored in the original histogram
        clear : bool, optional
            Whether to clear Axes before drawing (if passed); if False, this function will skip drawing the legend

    Returns
    -------
        ax : matplotlib.axes.Axes
            A matplotlib `Axes <https://matplotlib.org/3.1.1/api/axes_api.html>`_ object
    """
    if ax is None:
        fig, ax = plt.subplots(1, 1)
    else:
        if not isinstance(ax, plt.Axes):
            raise ValueError("ax must be a matplotlib Axes object")
        if clear:
            ax.clear()

    if ext_error is None:
        values = hmc.values() + np.sqrt(hmc.variances())
        baseline = hmc.values() - np.sqrt(hmc.variances())
        edges = hmc.axes[0].edges

        if flow == "show":
            edges = np.array(
                [
                    edges[0] - (edges[1] - edges[0]) * 2,
                    *edges,
                    edges[-1] + (edges[1] - edges[0]) * 2,
                ]
            )
            values = hmc.view(flow=True)["value"] + np.sqrt(
                hmc.view(flow=True)["variance"]
            )
            baseline = hmc.view(flow=True)["value"] - np.sqrt(
                hmc.view(flow=True)["variance"]
            )
        if flow == "sum":
            values[0], values[-1] = (
                hmc.view(flow=True)["value"] + np.sqrt(hmc.view(flow=True)["variance"])
            )[0] + values[0], (
                hmc.view(flow=True)["value"] + np.sqrt(hmc.view(flow=True)["variance"])
            )[
                -1
            ] + values[
                -1
            ]
            baseline[0], baseline[-1] = (
                hmc.view(flow=True)["value"] - np.sqrt(hmc.view(flow=True)["variance"])
            )[0] + baseline[0], (
                hmc.view(flow=True)["value"] - np.sqrt(hmc.view(flow=True)["variance"])
            )[
                -1
            ] + baseline[
                -1
            ]

        ax.stairs(
            values=values,
            baseline=baseline,
            edges=edges,
            label=label,
            **errband_opts,
        )
    else:
        values = hmc.values() + ext_error[1]
        baseline = hmc.values() - ext_error[0]
        edges = hmc.axes[-1].edges

        if flow == "show":
            edges = np.array(
                [
                    edges[0] - (edges[1] - edges[0]) * 2,
                    *edges,
                    edges[-1] + (edges[1] - edges[0]) * 2,
                ]
            )
            values = hmc.view(flow=True)["value"] + ext_error[1]
            baseline = hmc.view(flow=True)["value"] - ext_error[0]
        if flow == "sum":
            values[0], values[-1] = (hmc.view(flow=True)["value"] + ext_error[1])[
                0
            ] + values[0], (hmc.view(flow=True)["value"] + ext_error[0])[-1] + values[
                -1
            ]
            baseline[0], baseline[-1] = (hmc.view(flow=True)["value"] - ext_error[1])[
                0
            ] + baseline[0], (hmc.view(flow=True)["value"] - ext_error[0])[
                -1
            ] + baseline[
                -1
            ]

        ax.stairs(
            values=values,
            baseline=baseline,
            edges=edges,
            label=label,
            **fill_opts,
        )


## Very nice implementtion from Kenneth Long:
## https://gist.github.com/kdlong/d697ee691c696724fc656186c25f8814
def rebin_hist(h, axis_name, edges):
    if type(edges) == int:
        return h[{axis_name: hist.rebin(edges)}]

    ax = h.axes[axis_name]
    ax_idx = [a.name for a in h.axes].index(axis_name)
    if not all([np.isclose(x, ax.edges).any() for x in edges]):
        raise ValueError(
            f"Cannot rebin histogram due to incompatible edges for axis '{ax.name}'\n"
            f"Edges of histogram are {ax.edges}, requested rebinning to {edges}"
        )

    # If you rebin to a subset of initial range, keep the overflow and underflow
    overflow = ax.traits.overflow or (
        edges[-1] < ax.edges[-1] and not np.isclose(edges[-1], ax.edges[-1])
    )
    underflow = ax.traits.underflow or (
        edges[0] > ax.edges[0] and not np.isclose(edges[0], ax.edges[0])
    )
    flow = overflow or underflow
    new_ax = hist.axis.Variable(
        edges, name=ax.name, label=ax.label, overflow=overflow, underflow=underflow
    )
    axes = list(h.axes)
    axes[ax_idx] = new_ax

    hnew = hist.Hist(*axes, name=h.name, storage=h._storage_type())

    # Offset from bin edge to avoid numeric issues
    offset = 0.5 * np.min(ax.edges[1:] - ax.edges[:-1])
    edges_eval = edges + offset
    edge_idx = ax.index(edges_eval)
    # Avoid going outside the range, reduceat will add the last index anyway
    if edge_idx[-1] == ax.size + ax.traits.overflow:
        edge_idx = edge_idx[:-1]

    if underflow:
        # Only if the original axis had an underflow should you offset
        if ax.traits.underflow:
            edge_idx += 1
        edge_idx = np.insert(edge_idx, 0, 0)

    # Take is used because reduceat sums i:len(array) for the last entry, in the case
    # where the final bin isn't the same between the initial and rebinned histogram, you
    # want to drop this value. Add tolerance of 1/2 min bin width to avoid numeric issues
    hnew.values(flow=flow)[...] = np.add.reduceat(
        h.values(flow=flow), edge_idx, axis=ax_idx
    ).take(indices=range(new_ax.size + underflow + overflow), axis=ax_idx)
    if hnew._storage_type() == hist.storage.Weight():
        hnew.variances(flow=flow)[...] = np.add.reduceat(
            h.variances(flow=flow), edge_idx, axis=ax_idx
        ).take(indices=range(new_ax.size + underflow + overflow), axis=ax_idx)
    return hnew
