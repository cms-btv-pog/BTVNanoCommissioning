import numpy as np
import argparse, sys, os, arrow, glob
from coffea.util import load
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText

import mplhep as hep
import hist
from hist.intervals import ratio_uncertainty

plt.style.use(hep.style.ROOT)
from BTVNanoCommissioning.utils.xs_scaler import getSumW, collate, scaleSumW

parser = argparse.ArgumentParser(description="hist plotter for commissioning")
parser.add_argument("--lumi", required=True, type=float, help="luminosity in /pb")
parser.add_argument("--com", default="13", type=str, help="sqrt(s) in TeV")
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
parser.add_argument(
    "-i",
    "--input",
    type=str,
    default="",
    help="input coffea files (str), splitted different files with ','. Wildcard option * available as well.",
)
parser.add_argument("--log", action="store_true", help="log on y axis")
parser.add_argument(
    "--norm",
    default=False,
    type=bool,
    help="Use for reshape SF, scale to same yield as no SFs case",
)
parser.add_argument(
    "-v",
    "--variable",
    type=str,
    help="variables to plot, splitted by ,. Wildcard option * available as well. Specifying `all` will run through all variables.",
)
parser.add_argument(
    "--SF", action="store_true", default=False, help="make w/, w/o SF comparisons"
)
parser.add_argument("--ext", type=str, default="data", help="prefix name")
parser.add_argument(
    "--autorebin",
    type=int,
    default=1,
    help="Rebin the plotting variables by merging N bins in case the current binning is too fine for you ",
)
arg = parser.parse_args()
time = arrow.now().format("YY_MM_DD")
if not os.path.isdir(f"plot/BTV/{arg.phase}_{arg.ext}_{time}/"):
    os.makedirs(f"plot/BTV/{arg.phase}_{arg.ext}_{time}/")
if len(arg.input.split(",")) > 1:
    output = {i: load(i) for i in arg.input.split(",")}
    for out in output.keys():
        output[out] = scaleSumW(output[out], arg.lumi, getSumW(output[out]))
elif "*" in arg.input:
    files = glob.glob(arg.input)
    output = {i: load(i) for i in files}
    for out in output.keys():
        output[out] = scaleSumW(output[out], arg.lumi, getSumW(output[out]))
else:
    output = load(arg.input)
    output = scaleSumW(output, arg.lumi, getSumW(output))
mergemap = {}
if not any(".coffea" in o for o in output.keys()):
    mergemap["data"] = [m for m in output.keys() if "Run" in m]
    mergemap["mc"] = [m for m in output.keys() if "Run" not in m]
else:
    datalist = []
    mclist = []
    for f in output.keys():
        datalist.extend([m for m in output[f].keys() if "Run" in m])
        mclist.extend([m for m in output[f].keys() if "Run" not in m])
    mergemap["mc"] = mclist
    mergemap["data"] = datalist
collated = collate(output, mergemap)
if "Wc" in arg.phase:
    input_txt = "W+c"
elif "DY" in arg.phase:
    input_txt = "DY+jets"
elif "semilep" in arg.phase:
    input_txt = "semileptonic ttbar"
    nj = 4
elif "dilep" in arg.phase:
    input_txt = "dileptonic ttbar"
    nj = 2
if (
    "njet" in arg.variable.split(",")
    or "nbjet" in arg.variable.split(",")
    or "mu" in arg.variable.split(",")
):
    nj = 1

if arg.variable == "all":
    var_set = collated["mc"].keys()
elif "*" in arg.variable:
    var_set = [
        var for var in collated["mc"].keys() if arg.variable.replace("*", "") in var
    ]
else:
    var_set = arg.variable.split(",")

for discr in var_set:
    if "sumw" == discr:
        continue
    if arg.autorebin != 1:
        rebin_factor = arg.autorebin
        if len(collated["data"][discr].axes.name) == 3:
            collated["data"][discr] = collated["data"][discr][
                :, :, hist.rebin(rebin_factor)
            ]
            collated["mc"][discr] = collated["mc"][discr][
                :, :, hist.rebin(rebin_factor)
            ]
        elif len(collated["data"][discr].axes.name) == 2:
            collated["data"][discr] = collated["data"][discr][
                :, hist.rebin(rebin_factor)
            ]
            collated["mc"][discr] = collated["mc"][discr][:, hist.rebin(rebin_factor)]
        else:
            collated["data"][discr] = collated["data"][discr][hist.rebin(rebin_factor)]
            collated["mc"][discr] = collated["mc"][discr][hist.rebin(rebin_factor)]

    if (
        "flav" in collated["mc"][discr].axes.name
        and "syst" in collated["mc"][discr].axes.name
        and arg.SF
    ):
        scale_sf = np.sum(
            collated["mc"][discr][{"syst": "SF", "flav": sum}].values()
        ) / np.sum(collated["mc"][discr][{"syst": "noSF", "flav": sum}].values())
    else:
        scale_sf = 1.0
    if (
        "flav" in collated["mc"][discr].axes.name
        and "syst" in collated["mc"][discr].axes.name
        and arg.SF
    ):
        print("============> fraction of each flavor in MC")
        print(
            "b:%.3f",
            np.sum(collated["mc"][discr][{"syst": "SF", "flav": 3}].values())
            / np.sum(collated["mc"][discr][{"syst": "SF", "flav": sum}].values()),
        )
        print(
            "c:%.3f",
            np.sum(collated["mc"][discr][{"syst": "SF", "flav": 2}].values())
            / np.sum(collated["mc"][discr][{"syst": "SF", "flav": sum}].values()),
        )
        print(
            "pu:%.3f",
            np.sum(collated["mc"][discr][{"syst": "SF", "flav": 1}].values())
            / np.sum(collated["mc"][discr][{"syst": "SF", "flav": sum}].values()),
        )
        print(
            "l:%.3f",
            np.sum(collated["mc"][discr][{"syst": "SF", "flav": 0}].values())
            / np.sum(collated["mc"][discr][{"syst": "SF", "flav": sum}].values()),
        )
        print("============")
    fig, ((ax), (rax)) = plt.subplots(
        2, 1, figsize=(10, 10), gridspec_kw={"height_ratios": (3, 1)}, sharex=True
    )
    fig.subplots_adjust(hspace=0.06, top=0.92, bottom=0.1, right=0.97)
    hep.cms.label(
        "Preliminary", data=True, lumi=arg.lumi / 1000.0, com=arg.com, loc=0, ax=ax
    )

    if (
        "flav" in collated["mc"][discr].axes.name
        and "syst" in collated["mc"][discr].axes.name
        and arg.SF
    ):

        err_up = (
            collated["mc"][discr][{"syst": "SFup", "flav": sum}].values()
            - collated["mc"][discr][{"syst": "SF", "flav": sum}].values()
        )
        err_dn = (
            collated["mc"][discr][{"syst": "SF", "flav": sum}].values()
            - collated["mc"][discr][{"syst": "SFdn", "flav": sum}].values()
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
                            collated["mc"][discr][{"syst": "SF", "flav": 4}].values()
                            - collated["mc"][discr][
                                {"syst": "SFdn", "flav": 4}
                            ].values()
                        ),
                        2,
                    ),
                )
            )

        hdata = collated["data"][discr][{"syst": "noSF", "flav": sum}]
        ratio_up = ratio_uncertainty(err_up, hdata.values())
        ratio_dn = ratio_uncertainty(err_dn, hdata.values())

        hep.histplot(
            [collated["mc"][discr][{"syst": "SF", "flav": i}] for i in range(4)],
            stack=True,
            label=["udsg", "pileup", "c", "b"],
            histtype="fill",
            yerr=True,
            ax=ax,
        )
        hep.histplot(
            collated["mc"][discr][{"syst": "noSF", "flav": sum}],
            label=["w/o SF"],
            color="tab:gray",
            width=2,
            yerr=True,
            ax=ax,
        )
        hep.histplot(
            hdata,
            histtype="errorbar",
            color="black",
            label="Data",
            yerr=True,
            ax=ax,
        )
        rax.errorbar(
            x=hdata.axes[0].centers,
            y=hdata.values()
            / collated["mc"][discr][{"syst": "SF", "flav": sum}].values(),
            yerr=ratio_uncertainty(
                hdata.values(),
                collated["mc"][discr][{"syst": "SF", "flav": sum}].values(),
            ),
            color="k",
            linestyle="none",
            marker="o",
            elinewidth=1,
        )
        rax.errorbar(
            x=hdata.axes[0].centers,
            y=hdata.values()
            / collated["mc"][discr][{"syst": "noSF", "flav": sum}].values(),
            yerr=ratio_uncertainty(
                hdata.values(),
                collated["mc"][discr][{"syst": "noSF", "flav": sum}].values(),
            ),
            color="tab:brown",
            linestyle="none",
            marker="o",
            elinewidth=1,
        )
        ### FIXME: errorband calculation
        # stat_denom_unc = ratio_uncertainty(
        #     hdata.values(),
        #     collated["mc"][discr][{"syst": "noSF", "flav": sum}].values(),
        # )
        # ax.fill_between(
        #     hdata.axes.edges,
        #     np.ones(stat_denom_unc[0])
        #     - np.r_[stat_denom_unc[0], stat_denom_unc[0, -1]],
        #     np.ones(stat_denom_unc[0])
        #     + np.r_[stat_denom_unc[1], stat_denom_unc[1, -1]],
        #     {"facecolor": "tab:gray", "linewidth": 0},
        # )
        # ax.fill_between(
        #     hdata.axes.edges,
        #     np.ones(stat_denom_unc[0]) - np.r_[err_dn, err_dn[-1]],
        #     np.ones(stat_denom_unc[0]) + np.r_[err_up, err_up[-1]],
        #     {"facecolor": "tab:brown", "linewidth": 0},
        # )
    elif "syst" in collated["mc"][discr].axes.name and not arg.SF:
        hep.histplot(
            [collated["mc"][discr][{"flav": i, "syst": "noSF"}] for i in range(4)],
            stack=True,
            histtype="fill",
            label=["udsg", "pileup", "c", "b"],
            yerr=True,
            ax=ax,
        )
        hep.histplot(
            collated["mc"][discr][{"flav": sum, "syst": "noSF"}],
            histtype="errorbar",
            label=["Stat. Unc"],
            yerr=True,
            color="gray",
            markersize=0.0,
            ax=ax,
        )
        hep.histplot(
            collated["data"][discr][{"flav": sum, "syst": "noSF"}],
            histtype="errorbar",
            color="black",
            label="Data",
            yerr=True,
            ax=ax,
        )
        rax.errorbar(
            x=collated["mc"][discr][{"flav": sum, "syst": "noSF"}].axes[0].centers,
            y=collated["data"][discr][{"flav": sum, "syst": "noSF"}].values()
            / collated["mc"][discr][{"flav": sum, "syst": "noSF"}].values(),
            yerr=ratio_uncertainty(
                collated["data"][discr][{"flav": sum, "syst": "noSF"}].values(),
                collated["mc"][discr][{"flav": sum, "syst": "noSF"}].values(),
            ),
            color="k",
            linestyle="none",
            marker="o",
            elinewidth=1,
        )
    elif "flav" in collated["mc"][discr].axes.name:
        hep.histplot(
            [collated["mc"][discr][{"flav": i}] for i in range(4)],
            stack=True,
            histtype="fill",
            label=["udsg", "pileup", "c", "b"],
            yerr=True,
            ax=ax,
        )
        hep.histplot(
            collated["mc"][discr][{"flav": sum}],
            histtype="errorbar",
            label=["Stat. Unc."],
            yerr=True,
            color="gray",
            markersize=0.0,
            ax=ax,
        )
        hep.histplot(
            collated["data"][discr][{"flav": sum}],
            histtype="errorbar",
            color="black",
            label="Data",
            yerr=True,
            ax=ax,
        )
        rax.errorbar(
            x=collated["mc"][discr][{"flav": sum}].axes[0].centers,
            y=collated["data"][discr][{"flav": sum}].values()
            / collated["mc"][discr][{"flav": sum}].values(),
            yerr=ratio_uncertainty(
                collated["data"][discr][{"flav": sum}].values(),
                collated["mc"][discr][{"flav": sum}].values(),
            ),
            color="k",
            linestyle="none",
            marker="o",
            elinewidth=1,
        )

    else:
        hep.histplot(
            collated["mc"][discr],
            stack=True,
            label=["MC"],
            yerr=True,
            ax=ax,
        )
        hep.histplot(
            collated["mc"][discr],
            histtype="errorbar",
            label=["Stat. Unc."],
            yerr=True,
            color="gray",
            markersize=0.0,
            ax=ax,
        )
        hep.histplot(
            collated["data"][discr],
            histtype="errorbar",
            color="black",
            label="Data",
            yerr=True,
            ax=ax,
        )
        rax.errorbar(
            x=collated["data"][discr].axes[0].centers,
            y=collated["data"][discr].values() / collated["mc"][discr].values(),
            yerr=ratio_uncertainty(
                collated["data"][discr].values(), collated["mc"][discr].values()
            ),
            color="k",
            linestyle="none",
            marker="o",
            elinewidth=1,
        )
    ax.set_xlabel(None)
    ax.set_ylabel("Events")
    rax.set_ylabel("Data/MC")
    rax.set_xlabel(discr)
    rax.axhline(y=1.0, linestyle="dashed", color="gray")
    ax.legend()
    rax.set_ylim(0.5, 1.5)
    ax.set_ylim(bottom=0.0)
    at = AnchoredText(
        input_txt + "\n" + arg.ext,
        loc=2,
        frameon=False,
    )
    ax.add_artist(at)
    scale = ""
    if arg.norm:
        scale = "_norm"
    name = "all"
    hep.mpl_magic(ax=ax)
    if arg.log:
        ax.set_yscale("log")
        name = "log"
        ax.set_ylim(bottom=0.1)
        hep.mpl_magic(ax=ax)
        fig.savefig(
            f"plot/BTV/{arg.phase}_{arg.ext}_{time}/unc_{discr}_inclusive{scale}_{name}.pdf"
        )
        fig.savefig(
            f"plot/BTV/{arg.phase}_{arg.ext}_{time}/unc_{discr}_inclusive{scale}_{name}.png"
        )
    else:
        fig.savefig(
            f"plot/BTV/{arg.phase}_{arg.ext}_{time}/unc_{discr}_inclusive{scale}_{name}.pdf"
        )
        fig.savefig(
            f"plot/BTV/{arg.phase}_{arg.ext}_{time}/unc_{discr}_inclusive{scale}_{name}.png"
        )
