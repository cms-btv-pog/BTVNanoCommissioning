import matplotlib.pyplot as plt
import numpy as np
from hist import intervals
import scipy.stats
import warnings

errband_opts = {
    "hatch": "////",
    "facecolor": "none",
    "lw": 0,
    "color": "k",
    "alpha": 0.4,
}
markers = [".", "o", "^", "s", "+", "x", "D", "*"]

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
    overflow="none",
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
        overflow : str, optional
            If overflow behavior is not 'none', extra bins will be drawn on either end of the nominal
            axis range, to represent the contents of the overflow bins.  See `Hist.sum` documentation
            for a description of the options.
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
    centers = axis.centers
    ranges = (edges[1:] - edges[:-1]) / 2 if xerr else None

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
def SFerror(collated, discr):

    err_up = (
        collated["mc"][discr][{"syst": "SFup", "flav": sum}].values()
        - collated["mc"][discr][{"syst": "SF", "flav": sum}].values()
    )
    err_dn = (
        collated["mc"][discr][{"syst": "SFdn", "flav": sum}].values()
        - collated["mc"][discr][{"syst": "SF", "flav": sum}].values()
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

    return np.array([err_dn, err_up])


def autoranger(hist):
    val, axis = hist.values(), hist.axes[-1].edges
    for i in range(len(val)):
        if val[i] != 0:
            mins = i - 1
            break
    for i in reversed(range(len(val))):
        if val[i] != 0:
            maxs = i + 1
            break
    if mins == -1:
        mins = 0
    #if maxs == len(val):
    #    maxs = maxs - 1
    return axis[mins], axis[maxs]
