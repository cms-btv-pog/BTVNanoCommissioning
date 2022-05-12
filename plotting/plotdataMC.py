import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from coffea.util import load
from coffea.hist import plot
from coffea import hist
import re
import argparse
import sys

data_err_opts = {
    "linestyle": "none",
    "marker": ".",
    "markersize": 10.0,
    "color": "k",
    "elinewidth": 1,
}
from cycler import cycler

sys.path.append("..")
from BTVNanoCommissioning.utils.xs_scaler import scale_xs

ptbin = [25, 30, 40, 60, 80, 100, 150, 200, 300, 500]
etabin = [-2.5, -2.0, -1.5, -0.5, 0.0, 0.5, 1.5, 2.0, 2.5]
notdata = re.compile("(?!data)")
parser = argparse.ArgumentParser(description="hist plotter for commissioning")
parser.add_argument("--lumi", required=True, type=float, help="luminosity in /pb")
parser.add_argument("-c", "--combine", type=bool, help="combined all the jets")
parser.add_argument(
    "-p",
    "--phase",
    required=True,
    choices=[
        "dilep_sf",
        "ttsemilep_sf",
        "ctag_Wc_sf",
        "ctag_DY_sf",
        "ctag_ttsemilep_sf",
        "ctag_ttdilep_sf",
    ],
    dest="phase",
    help="which phase space",
)
parser.add_argument("--log", type=bool, help="log on x axis")
parser.add_argument(
    "--norm",
    default=False,
    type=bool,
    help="Use for reshape SF, scale to same yield as no SFs case",
)
parser.add_argument(
    "-d",
    "--discr_list",
    nargs="+",
    default=[
        "deepcsv_CvL",
        "deepcsv_CvB",
        "deepflav_CvL",
        "deepflav_CvB",
        "btagDeepB",
        "btagDeepC",
        "btagDeepFlavB",
        "btagDeepFlavC",
    ],
    help="discriminators",
)
parser.add_argument("--ext", type=str, default="data", help="addional name")
parser.add_argument("-i", "--input", type=str, default="", help="input coffea files")
arg = parser.parse_args()
datas = re.compile(f"(?={arg.ext})")

output = load(arg.input)
events = output["sumw"]
if arg.phase == "dilep":
    input_txt = "dilepton ttbar"
    nj = 2
elif arg.phase == "ctag":
    input_txt = "semileptonic ttbar"
    nj = 4
else:
    if "Wc" in arg.phase:
        input_txt = "W+c"
    elif "DY" in arg.phase:
        input_txt = "DY+jets"
    elif "ttsemilep" in arg.phase:
        input_txt = "semileptonic ttbar"
    elif "ttdilep" in arg.phase:
        input_txt = "dileptonic ttbar"
    nj = 1
if arg.combine:
    nj = 1
if "njet" in arg.discr_list or "nbjet" in arg.discr_list or "mu" in arg.discr_list:
    nj = 1


for j in range(nj):
    for discr in arg.discr_list:
        if arg.combine:
            hflav_0 = output[f"{discr}_0"]
            hflav_1 = output[f"{discr}_1"]
            hflav = hflav_0 + hflav_1
            if arg.phase == "ctag":
                hflav_2 = output[f"{discr}_2"]
                hflav_3 = output[f"{discr}_3"]
                hflav = hflav_0 + hflav_1 + hflav_2 + hflav_3
        else:
            if "btag" in discr or "CvL" in discr or "CvB" in discr:
                hflav = output[f"{discr}_{j}"]
            else:
                hflav = output[discr]

        if "btag" in discr or "DeepCSV" in discr:
            hflav = hflav.rebin("flav", hist.Bin("flav", "flav", [0, 1, 4, 5, 6]))

        if ("btag" in discr or "CvL" in discr or "CvB" in discr) and arg.norm:
            if "Wc" in arg.phase:
                scale_sf = sum(
                    hflav[notdata]
                    .integrate("dataset")
                    .integrate("syst", "noSF")
                    .integrate("flav")
                    .integrate("char")
                    .values()[()]
                ) / sum(
                    hflav[notdata]
                    .integrate("dataset")
                    .integrate("syst", "SF")
                    .integrate("flav")
                    .integrate("char")
                    .values()[()]
                )

            else:
                scale_sf = sum(
                    hflav[notdata]
                    .integrate("dataset")
                    .integrate("syst", "noSF")
                    .integrate("flav")
                    .values()[()]
                ) / sum(
                    hflav[notdata]
                    .integrate("dataset")
                    .integrate("syst", "SF")
                    .integrate("flav")
                    .values()[()]
                )

        if not arg.norm:
            scale_sf = 1.0
        hflav = scale_xs(hflav, arg.lumi, events)
        if "btag" in discr or "CvL" in discr or "CvB" in discr:
            hflav = scale_xs(hflav, arg.lumi, events)
            if "Wc" in arg.phase:
                print("============")
                print(
                    sum(
                        hflav[notdata]
                        .integrate("dataset")
                        .integrate("syst", "noSF")
                        .integrate("flav")
                        .integrate("char")
                        .values()[()]
                    )
                )
                print(
                    "b:%.3f\%",
                    sum(
                        hflav[notdata]
                        .integrate("dataset")
                        .integrate("syst", "noSF")
                        .integrate("flav", slice(5, 6))
                        .integrate("char")
                        .values()[()]
                    )
                    / sum(
                        hflav[notdata]
                        .integrate("dataset")
                        .integrate("syst", "noSF")
                        .integrate("flav")
                        .integrate("char")
                        .values()[()]
                    ),
                )
                print(
                    "c:%.3f\%",
                    sum(
                        hflav[notdata]
                        .integrate("dataset")
                        .integrate("syst", "noSF")
                        .integrate("flav", slice(4, 5))
                        .integrate("char")
                        .values()[()]
                    )
                    / sum(
                        hflav[notdata]
                        .integrate("dataset")
                        .integrate("syst", "noSF")
                        .integrate("flav")
                        .integrate("char")
                        .values()[()]
                    ),
                )
                print(
                    "l:%.3f\%",
                    sum(
                        hflav[notdata]
                        .integrate("dataset")
                        .integrate("syst", "noSF")
                        .integrate("flav", slice(0, 4))
                        .integrate("char")
                        .values()[()]
                    )
                    / sum(
                        hflav[notdata]
                        .integrate("dataset")
                        .integrate("syst", "noSF")
                        .integrate("flav")
                        .integrate("char")
                        .values()[()]
                    ),
                )
                print("============")

            else:
                print("============")
                print(
                    sum(
                        hflav[notdata]
                        .integrate("dataset")
                        .integrate("syst", "noSF")
                        .integrate("flav")
                        .values()[()]
                    )
                )
                print(
                    "b:%.3f",
                    sum(
                        hflav[notdata]
                        .integrate("dataset")
                        .integrate("syst", "noSF")
                        .integrate("flav", slice(5, 6))
                        .values()[()]
                    )
                    / sum(
                        hflav[notdata]
                        .integrate("dataset")
                        .integrate("syst", "noSF")
                        .integrate("flav")
                        .values()[()]
                    ),
                )
                print(
                    "c:%.3f",
                    sum(
                        hflav[notdata]
                        .integrate("dataset")
                        .integrate("syst", "noSF")
                        .integrate("flav", slice(4, 5))
                        .values()[()]
                    )
                    / sum(
                        hflav[notdata]
                        .integrate("dataset")
                        .integrate("syst", "noSF")
                        .integrate("flav")
                        .values()[()]
                    ),
                )
                print(
                    "l:%.3f",
                    sum(
                        hflav[notdata]
                        .integrate("dataset")
                        .integrate("syst", "noSF")
                        .integrate("flav", slice(0, 4))
                        .values()[()]
                    )
                    / sum(
                        hflav[notdata]
                        .integrate("dataset")
                        .integrate("syst", "noSF")
                        .integrate("flav")
                        .values()[()]
                    ),
                )
                print("============")
        elif "nbjet" in discr:
            hflav = scale_xs(hflav, arg.lumi, events)
        fig, ((ax), (rax)) = plt.subplots(
            2, 1, figsize=(6, 6), gridspec_kw={"height_ratios": (3, 1)}, sharex=True
        )
        fig.subplots_adjust(hspace=0.07)
        if "btag" in discr or "CvL" in discr or "CvB" in discr:
            if "Wc" in arg.phase:
                if "C" in discr:
                    err_up = (
                        hflav[notdata]
                        .integrate("dataset")
                        .integrate("syst", "SFup")
                        .integrate("flav")
                        .integrate("char")
                        .values()[()]
                    )
                    err_dn = (
                        hflav[notdata]
                        .integrate("dataset")
                        .integrate("syst", "SFdn")
                        .integrate("flav")
                        .integrate("char")
                        .values()[()]
                    )
                else:
                    err_up = np.sqrt(
                        np.add(
                            np.power(
                                np.add(
                                    hflav[notdata]
                                    .integrate("dataset")
                                    .integrate("syst", "SFup")
                                    .integrate("flav")
                                    .integrate("char")
                                    .values()[()],
                                    -1.0
                                    * hflav[notdata]
                                    .integrate("dataset")
                                    .integrate("syst", "SF")
                                    .integrate("flav")
                                    .integrate("char")
                                    .values()[()],
                                ),
                                2,
                            ),
                            np.power(
                                2
                                * np.add(
                                    hflav[notdata]
                                    .integrate("dataset")
                                    .integrate("syst", "SFup")
                                    .integrate("flav", slice(5, 6))
                                    .integrate("char")
                                    .values()[()],
                                    -1.0
                                    * hflav[notdata]
                                    .integrate("dataset")
                                    .integrate("flav", slice(5, 6))
                                    .integrate("char")
                                    .values()[()],
                                ),
                                2,
                            ),
                        )
                    )
                    err_dn = np.sqrt(
                        np.add(
                            np.power(
                                np.add(
                                    hflav[notdata]
                                    .integrate("dataset")
                                    .integrate("syst", "SF")
                                    .integrate("flav")
                                    .integrate("char")
                                    .values()[()],
                                    -1.0
                                    * hflav[notdata]
                                    .integrate("dataset")
                                    .integrate("syst", "SFdn")
                                    .integrate("flav")
                                    .integrate("char")
                                    .values()[()],
                                ),
                                2,
                            ),
                            np.power(
                                2
                                * np.add(
                                    hflav[notdata]
                                    .integrate("dataset")
                                    .integrate("flav", slice(5, 6))
                                    .integrate("char")
                                    .values()[()],
                                    -1.0
                                    * hflav[notdata]
                                    .integrate("dataset")
                                    .integrate("syst", "SFdn")
                                    .integrate("flav", slice(5, 6))
                                    .integrate("char")
                                    .values()[()],
                                ),
                                2,
                            ),
                        )
                    )
                data = (
                    hflav[datas]
                    .integrate("dataset")
                    .integrate("syst", "noSF")
                    .integrate("flav")
                    .integrate("char")
                    .values()[()]
                )
                maximum = max(
                    max(
                        (
                            hflav[notdata]
                            .integrate("dataset")
                            .integrate("syst", "SF")
                            .integrate("flav")
                            .integrate("char")
                            .values()[()]
                            + err_up
                        )
                    ),
                    max(data),
                )
            else:
                if "C" in discr:
                    err_up = (
                        hflav[notdata]
                        .integrate("dataset")
                        .integrate("syst", "SFup")
                        .integrate("flav")
                        .values()[()]
                    )
                    err_dn = (
                        hflav[notdata]
                        .integrate("dataset")
                        .integrate("syst", "SFdn")
                        .integrate("flav")
                        .values()[()]
                    )
                else:
                    err_up = np.sqrt(
                        np.add(
                            np.power(
                                np.add(
                                    hflav[notdata]
                                    .integrate("dataset")
                                    .integrate("syst", "SFup")
                                    .integrate("flav")
                                    .values()[()],
                                    -1.0
                                    * hflav[notdata]
                                    .integrate("dataset")
                                    .integrate("syst", "SF")
                                    .integrate("flav")
                                    .values()[()],
                                ),
                                2,
                            ),
                            np.power(
                                2
                                * np.add(
                                    hflav[notdata]
                                    .integrate("dataset")
                                    .integrate("syst", "SFup")
                                    .integrate("flav", slice(5, 6))
                                    .values()[()],
                                    -1.0
                                    * hflav[notdata]
                                    .integrate("dataset")
                                    .integrate("flav", slice(5, 6))
                                    .values()[()],
                                ),
                                2,
                            ),
                        )
                    )
                    err_dn = np.sqrt(
                        np.add(
                            np.power(
                                np.add(
                                    hflav[notdata]
                                    .integrate("dataset")
                                    .integrate("syst", "SF")
                                    .integrate("flav")
                                    .values()[()],
                                    -1.0
                                    * hflav[notdata]
                                    .integrate("dataset")
                                    .integrate("syst", "SFdn")
                                    .integrate("flav")
                                    .values()[()],
                                ),
                                2,
                            ),
                            np.power(
                                2
                                * np.add(
                                    hflav[notdata]
                                    .integrate("dataset")
                                    .integrate("flav", slice(5, 6))
                                    .values()[()],
                                    -1.0
                                    * hflav[notdata]
                                    .integrate("dataset")
                                    .integrate("syst", "SFdn")
                                    .integrate("flav", slice(5, 6))
                                    .values()[()],
                                ),
                                2,
                            ),
                        )
                    )
                data = hflav[datas].integrate("dataset").integrate("flav").values()[()]
                maximum = max(
                    max(
                        (
                            hflav[notdata]
                            .integrate("dataset")
                            .integrate("syst", "SF")
                            .integrate("flav")
                            .values()[()]
                            + err_up
                        )
                    ),
                    max(data),
                )
            ratio_up = np.divide(
                err_up, data, out=np.zeros_like(err_up), where=data != 0
            )
            ratio_dn = np.divide(
                err_dn, data, out=np.zeros_like(err_up), where=data != 0
            )

            if "Wc" in arg.phase:
                ax = plot.plot1d(
                    hflav[notdata].sum("dataset").integrate("char"),
                    overlay="flav",
                    fill_opts={},
                    error_opts=None,
                    ax=ax,
                    stack=True,
                )
                plot.plot1d(
                    hflav[notdata].sum("dataset", sumw2=True).sum("flav").sum("char"),
                    ax=ax,
                    density=False,
                    clear=False,
                    error_opts={
                        "linestyle": "none",
                        "markersize": 0,
                        "elinewidth": 10,
                        "color": "tab:brown",
                        "alpha": 0.3,
                        "yerr": [err_up, err_dn],
                    },
                )
                plot.plot1d(
                    hflav[notdata].sum("dataset").sum("flav").sum("char"),
                    ax=ax,
                    error_opts={
                        "linestyle": "none",
                        "markersize": 0,
                        "elinewidth": 10,
                        "color": "tab:gray",
                        "alpha": 0.3,
                    },
                    clear=False,
                    density=False,
                )

                plot.plot1d(
                    hflav[notdata]
                    .sum("dataset")
                    .integrate("syst", "noSF")
                    .integrate("char")
                    .sum("flav"),
                    error_opts={
                        "linestyle": "none",
                        "marker": "o",
                        "markersize": 5.0,
                        "mfc": "none",
                        "color": "tab:pink",
                        "elinewidth": 1.5,
                    },
                    ax=ax,
                    clear=False,
                )
                plot.plot1d(
                    hflav[datas].sum("dataset").sum("flav").sum("char"),
                    error_opts=data_err_opts,
                    ax=ax,
                    clear=False,
                    density=False,
                )
                ax.legend(
                    ncol=2,
                    loc="upper right",
                    handles=ax.get_legend_handles_labels()[0],
                    labels=[
                        "b",
                        "c",
                        "pileup",
                        "udsg",
                        "SFs Unc.",
                        "stat. Unc.",
                        "w/o SFs",
                        arg.ext,
                    ],
                    fontsize=13,
                )
                ax.set_xlabel(None)

                ax.set_xticklabels(ax.get_xticklabels(), fontsize=0)
                rax = plot.plotratio(
                    num=hflav[datas].sum("dataset").sum("flav").sum("char"),
                    denom=hflav[notdata].sum("dataset").sum("flav").sum("char"),
                    ax=rax,
                    error_opts={
                        "linestyle": "none",
                        "marker": ".",
                        "markersize": 0.0,
                        "color": "k",
                    },
                    denom_fill_opts={"yerr": [ratio_up, ratio_dn]},
                    # denom_fill_opts={},
                    guide_opts={},
                    unc="num",
                    clear=False,
                )
                plot.plotratio(
                    num=hflav[datas].sum("dataset").sum("flav").sum("char"),
                    denom=hflav[notdata].sum("dataset").sum("flav").sum("char"),
                    ax=rax,
                    error_opts=data_err_opts,
                    denom_fill_opts={},
                    guide_opts={},
                    unc="num",
                    clear=False,
                )
                plot.plotratio(
                    num=hflav[datas]
                    .sum("dataset")
                    .integrate("syst", "noSF")
                    .sum("flav")
                    .sum("char"),
                    denom=hflav[notdata]
                    .sum("dataset")
                    .integrate("syst", "noSF")
                    .sum("flav")
                    .sum("char"),
                    ax=rax,
                    error_opts={
                        "linestyle": "none",
                        "marker": "o",
                        "markersize": 5.0,
                        "mfc": "none",
                        "color": "tab:pink",
                        "elinewidth": 1.5,
                    },
                    denom_fill_opts={},
                    guide_opts={},
                    unc="num",
                    clear=False,
                )
            else:
                ax = plot.plot1d(
                    hflav[notdata].sum("dataset"),
                    overlay="flav",
                    fill_opts={},
                    error_opts=None,
                    ax=ax,
                    stack=True,
                )

                plot.plot1d(
                    hflav[notdata].sum("dataset", sumw2=True).sum("flav"),
                    ax=ax,
                    density=False,
                    clear=False,
                    error_opts={
                        "linestyle": "none",
                        "markersize": 0,
                        "elinewidth": 10,
                        "color": "tab:brown",
                        "alpha": 0.3,
                        "yerr": [err_up, err_dn],
                    },
                )
                plot.plot1d(
                    hflav[notdata].sum("dataset").sum("flav"),
                    ax=ax,
                    error_opts={
                        "linestyle": "none",
                        "markersize": 0,
                        "elinewidth": 10,
                        "color": "tab:gray",
                        "alpha": 0.3,
                    },
                    clear=False,
                    density=False,
                )

                plot.plot1d(
                    hflav[notdata].sum("dataset").integrate("syst", "noSF").sum("flav"),
                    error_opts={
                        "linestyle": "none",
                        "marker": "o",
                        "markersize": 5.0,
                        "mfc": "none",
                        "color": "tab:pink",
                        "elinewidth": 1.5,
                    },
                    ax=ax,
                    clear=False,
                )
                plot.plot1d(
                    hflav[datas].sum("dataset").sum("flav"),
                    error_opts=data_err_opts,
                    ax=ax,
                    clear=False,
                    density=False,
                )
                ax.legend(
                    ncol=2,
                    loc="upper right",
                    handles=ax.get_legend_handles_labels()[0],
                    labels=[
                        "b",
                        "c",
                        "pileup",
                        "udsg",
                        "SFs Unc.",
                        "stat. Unc.",
                        "w/o SFs",
                        arg.ext,
                    ],
                    fontsize=13,
                )

                ax.set_xlabel(None)

                ax.set_xticklabels(ax.get_xticklabels(), fontsize=0)
                rax = plot.plotratio(
                    num=hflav[datas].sum("dataset").sum("flav"),
                    denom=hflav[notdata].sum("dataset").sum("flav"),
                    ax=rax,
                    error_opts={
                        "linestyle": "none",
                        "marker": ".",
                        "markersize": 0.0,
                        "color": "k",
                    },
                    denom_fill_opts={"yerr": [ratio_up, ratio_dn]},
                    # denom_fill_opts={},
                    guide_opts={},
                    unc="num",
                    clear=False,
                )
                plot.plotratio(
                    num=hflav[datas].sum("dataset").sum("flav"),
                    denom=hflav[notdata].sum("dataset").sum("flav"),
                    ax=rax,
                    error_opts=data_err_opts,
                    denom_fill_opts={},
                    guide_opts={},
                    unc="num",
                    clear=False,
                )
                plot.plotratio(
                    num=hflav[datas]
                    .sum("dataset")
                    .integrate("syst", "noSF")
                    .sum("flav"),
                    denom=hflav[notdata]
                    .sum("dataset")
                    .integrate("syst", "noSF")
                    .sum("flav"),
                    ax=rax,
                    error_opts={
                        "linestyle": "none",
                        "marker": "o",
                        "markersize": 5.0,
                        "mfc": "none",
                        "color": "tab:pink",
                        "elinewidth": 1.5,
                    },
                    denom_fill_opts={},
                    guide_opts={},
                    unc="num",
                    clear=False,
                )
        elif "nbjet" in discr:

            err_up = np.add(
                hflav[notdata]
                .integrate("dataset")
                .integrate("syst", "SFup")
                .values()[()],
                -1.0 * hflav[notdata].integrate("dataset").values()[()],
            )
            err_dn = np.add(
                hflav[notdata].integrate("dataset").values()[()],
                -1.0
                * hflav[notdata]
                .integrate("dataset")
                .integrate("syst", "SFdn")
                .values()[()],
            )
            data = hflav[datas].integrate("dataset").values()[()]
            if not arg.log:
                maximum = max(
                    max((hflav[notdata].integrate("dataset").values()[()] + err_up)),
                    max(data),
                )
            ratio_up = np.divide(
                err_up, data, out=np.zeros_like(err_up), where=data != 0
            )
            ratio_dn = np.divide(
                err_dn, data, out=np.zeros_like(err_up), where=data != 0
            )

            ax = plot.plot1d(
                hflav[notdata].sum("dataset"), fill_opts={}, error_opts=None, ax=ax
            )
            plot.plot1d(
                hflav[notdata].sum("dataset", sumw2=True),
                ax=ax,
                density=False,
                clear=False,
                error_opts={
                    "linestyle": "none",
                    "markersize": 0,
                    "elinewidth": 10,
                    "color": "tab:brown",
                    "alpha": 0.3,
                    "yerr": [err_up, err_dn],
                },
            )
            plot.plot1d(
                hflav[notdata].sum("dataset"),
                ax=ax,
                error_opts={
                    "linestyle": "none",
                    "markersize": 0,
                    "elinewidth": 10,
                    "color": "tab:gray",
                    "alpha": 0.3,
                },
                clear=False,
                density=False,
            )

            plot.plot1d(
                hflav[datas].sum("dataset"),
                error_opts=data_err_opts,
                ax=ax,
                clear=False,
                density=False,
            )
            ax.legend(
                ncol=2,
                loc="upper right",
                handles=ax.get_legend_handles_labels()[0],
                labels=["MC", "SFs Unc.", "stat. Unc.", arg.ext],
                fontsize=13,
            )
            ax.set_xlabel(None)

            rax = plot.plotratio(
                num=hflav[datas].sum("dataset"),
                denom=hflav[notdata].sum("dataset"),
                ax=rax,
                error_opts={
                    "linestyle": "none",
                    "marker": ".",
                    "markersize": 0.0,
                    "color": "k",
                },
                denom_fill_opts={"yerr": [ratio_up, ratio_dn]},
                # denom_fill_opts={},
                guide_opts={},
                unc="num",
                clear=False,
            )
            plot.plotratio(
                num=hflav[datas].sum("dataset"),
                denom=hflav[notdata].sum("dataset"),
                ax=rax,
                error_opts=data_err_opts,
                denom_fill_opts={},
                guide_opts={},
                unc="num",
                clear=False,
            )

        elif "DeepCSV_" in discr:

            if "Wc" in arg.phase:
                maximum = max(
                    max(
                        (
                            hflav[notdata]
                            .sum("dataset")
                            .sum("char")
                            .sum("flav")
                            .values()[()]
                        )
                    ),
                    max(
                        hflav[datas].sum("dataset").sum("flav").sum("char").values()[()]
                    ),
                )
            else:
                maximum = max(
                    max((hflav[notdata].sum("dataset").sum("flav").values()[()])),
                    max(hflav[datas].sum("dataset").sum("flav").values()[()]),
                )

            # maximum=max(max((hflav[notdata].integrate("dataset").integrate("syst","SF").integrate("flav").sum("c").values()[()])),max(data))
            if "Wc" in arg.phase:
                # ax=plot.plot1d(hflav[notdata].sum("flav").integrate("char"),overlay="dataset",fill_opts={},error_opts=None,ax=ax,stack=True)
                print(hflav[datas].sum("flav").values())
                ax = plot.plot1d(
                    hflav[datas].sum("dataset").sum("flav").integrate("char"),
                    error_opts=data_err_opts,
                    ax=ax,
                    density=False,
                )
                # ax.legend(ncol=2,loc="upper right",handles=ax.get_legend_handles_labels()[0],labels=['b','c','pileup','udsg',arg.ext],fontsize=13)
                ax.legend(fontsize=8)
                ax.set_xlabel(None)

                ax.set_xticklabels(ax.get_xticklabels(), fontsize=0)

                rax = plot.plotratio(
                    num=hflav[datas].sum("dataset").sum("flav").sum("char"),
                    denom=hflav[notdata].sum("dataset").sum("flav").sum("char"),
                    ax=rax,
                    error_opts=data_err_opts,
                    denom_fill_opts={},
                    guide_opts={},
                    unc="num",
                )

            else:
                ax = plot.plot1d(
                    hflav[notdata].sum("dataset"),
                    overlay="flav",
                    fill_opts={},
                    error_opts=None,
                    ax=ax,
                    stack=True,
                )

                plot.plot1d(
                    hflav[notdata].sum("dataset", sumw2=True).sum("flav"),
                    ax=ax,
                    density=False,
                    clear=False,
                    error_opts={
                        "linestyle": "none",
                        "markersize": 0,
                        "elinewidth": 10,
                        "color": "tab:brown",
                        "alpha": 0.3,
                        "yerr": [err_up, err_dn],
                    },
                )
                plot.plot1d(
                    hflav[notdata].sum("dataset").sum("flav"),
                    ax=ax,
                    error_opts={
                        "linestyle": "none",
                        "markersize": 0,
                        "elinewidth": 10,
                        "color": "tab:gray",
                        "alpha": 0.3,
                    },
                    clear=False,
                    density=False,
                )

                plot.plot1d(
                    hflav[notdata].sum("dataset").integrate("syst", "noSF").sum("flav"),
                    error_opts={
                        "linestyle": "none",
                        "marker": "o",
                        "markersize": 5.0,
                        "mfc": "none",
                        "color": "tab:pink",
                        "elinewidth": 1.5,
                    },
                    ax=ax,
                    clear=False,
                )
                plot.plot1d(
                    hflav[datas].sum("dataset").sum("flav"),
                    error_opts=data_err_opts,
                    ax=ax,
                    clear=False,
                    density=False,
                )
                ax.legend(
                    ncol=2,
                    loc="upper right",
                    handles=ax.get_legend_handles_labels()[0],
                    labels=[
                        "b",
                        "c",
                        "pileup",
                        "udsg",
                        "SFs Unc.",
                        "stat. Unc.",
                        "w/o SFs",
                        arg.ext,
                    ],
                    fontsize=13,
                )

                ax.set_xlabel(None)

                # ax.set_xticklabels(ax.get_xticklabels(), fontsize=0)
                rax = plot.plotratio(
                    num=hflav[datas].sum("dataset").sum("flav"),
                    denom=hflav[notdata].sum("dataset").sum("flav"),
                    ax=rax,
                    error_opts={
                        "linestyle": "none",
                        "marker": ".",
                        "markersize": 0.0,
                        "color": "k",
                    },
                    denom_fill_opts={"yerr": [ratio_up, ratio_dn]},
                    # denom_fill_opts={},
                    guide_opts={},
                    unc="num",
                    clear=False,
                )
                plot.plotratio(
                    num=hflav[datas].sum("dataset").sum("flav"),
                    denom=hflav[notdata].sum("dataset").sum("flav"),
                    ax=rax,
                    error_opts=data_err_opts,
                    denom_fill_opts={},
                    guide_opts={},
                    unc="num",
                    clear=False,
                )
                plot.plotratio(
                    num=hflav[datas]
                    .sum("dataset")
                    .integrate("syst", "noSF")
                    .sum("flav"),
                    denom=hflav[notdata]
                    .sum("dataset")
                    .integrate("syst", "noSF")
                    .sum("flav"),
                    ax=rax,
                    error_opts={
                        "linestyle": "none",
                        "marker": "o",
                        "markersize": 5.0,
                        "mfc": "none",
                        "color": "tab:pink",
                        "elinewidth": 1.5,
                    },
                    denom_fill_opts={},
                    guide_opts={},
                    unc="num",
                    clear=False,
                )
        elif "mu" not in discr and "njet" != discr:
            data = hflav[datas].sum("dataset").sum("char").values()[()]
            # print(data)
            # if not arg.log:
            maximum = max(
                max((hflav[notdata].sum("dataset").sum("char").values()[()])), max(data)
            )

            ax = plot.plot1d(
                hflav[notdata].sum("dataset").sum("char"), error_opts=None, ax=ax
            )
            plot.plot1d(
                hflav[datas].sum("dataset").sum("char"),
                error_opts=data_err_opts,
                ax=ax,
                clear=False,
                density=False,
            )
            ax.legend(
                ncol=2,
                loc="upper right",
                handles=ax.get_legend_handles_labels()[0],
                labels=["MC", arg.ext],
                fontsize=13,
            )
            ax.set_xlabel(None)
            # ax.set_xticklabels(ax.get_xticklabels(), fontsize=0)

            rax = plot.plotratio(
                num=hflav[datas].sum("dataset").sum("char"),
                denom=hflav[notdata].sum("dataset").sum("char"),
                ax=rax,
                error_opts=data_err_opts,
                denom_fill_opts={},
                guide_opts={},
                unc="num",
                clear=False,
            )
        else:
            if not arg.log:
                maximum = max(
                    max(
                        (
                            hflav[notdata]
                            .integrate("dataset")
                            .integrate("syst", "SF")
                            .integrate("flav")
                            .values()[()]
                        )
                    ),
                    max(data),
                )

            ax = plot.plot1d(
                hflav[notdata].sum("dataset"),
                overlay="flav",
                fill_opts={},
                error_opts=None,
                ax=ax,
                stack=True,
            )
            plot.plot1d(
                hflav[datas].sum("dataset").sum("flav"),
                error_opts=data_err_opts,
                ax=ax,
                clear=False,
                density=False,
            )
            ax.legend(
                ncol=2,
                loc="upper right",
                handles=ax.get_legend_handles_labels()[0],
                labels=["b", "c", "pileup", "udsg", arg.ext],
                fontsize=13,
            )
            ax.set_xlabel(None)
            # ax.set_xticklabels(ax.get_xticklabels(), fontsize=0)

            rax = plot.plotratio(
                num=hflav[datas].sum("dataset").sum("flav"),
                denom=hflav[notdata].sum("dataset").sum("flav"),
                ax=rax,
                error_opts=data_err_opts,
                denom_fill_opts={},
                guide_opts={},
                unc="num",
                clear=False,
            )
        # if arg.log:
        #     ax.set_ylim(1,maximum*100)
        #     ax.semilogy()
        # else:ax.set_ylim(0,maximum*1.8)

        ax.set_ylabel("Events", fontsize=15)
        rax.set_ylabel("Data/MC", fontsize=15)
        if "CvL" in discr:
            discrs = discr.replace("CvL", "CvB")
        elif "CvB" in discr:
            discrs = discr.replace("CvB", "CvL")
        else:
            discrs = discr
        if arg.combine:
            rax.set_xlabel(discrs, fontsize=15)
        else:
            rax.set_xlabel(f"{discrs}[{j}]", fontsize=15)
        rax.set_ylim(0.5, 1.5)

        at = AnchoredText(
            input_txt + "\n"
            # + "inclusive pT, $\eta$"
            ,
            loc=2,
            prop=dict(size=15),
            frameon=False,
        )
        ax.add_artist(at)
        scale = ""
        if arg.norm:
            scale = "_norm"
        name = "all"
        if not arg.combine:
            name = str(j)
        if arg.log:
            fig.savefig(
                f"{arg.phase}_unc_{discrs}_inclusive{scale}_{arg.ext}_{name}.pdf"
            )
        else:
            fig.savefig(
                f"{arg.phase}_unc_lin_{discrs}_inclusive{scale}_{arg.ext}_{name}.pdf"
            )
