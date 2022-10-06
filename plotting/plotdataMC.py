import numpy as np
import argparse, sys, os, arrow
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
    type=str,
    help="discriminators",
)
parser.add_argument("--SF", action="store_true", default=False, help="SF comparisons")
parser.add_argument("--ext", type=str, default="data", help="addional name")
parser.add_argument("-o", "--output", type=str, default="", help="input coffea files")
arg = parser.parse_args()
time = arrow.now().format("YY_MM_DD")
if not os.path.isdir(f"plot/BTV/{arg.phase}_{arg.ext}_{time}/"):
    os.makedirs(f"plot/BTV/{arg.phase}_{arg.ext}_{time}/")
if len(arg.output.split(",")) > 1:
    output = {i.replace(".coffea", ""): load(i) for i in arg.output.split(",")}
    for out in output.keys():
        output[out] = scaleSumW(output[out], arg.lumi, getSumW(output[out]))
else:
    output = load(arg.output)
    output = scaleSumW(output, arg.lumi, getSumW(output))
mergemap = {}
if "sumw" in output.keys():
    mergemap["data"] = [
        m for m in output.keys() if "Run" in m or "data" in m or "Data" in m
    ]
    mergemap["mc"] = [
        m
        for m in output.keys()
        if "Run" not in m and "data" not in m and "Data" not in m
    ]
else:
    datalist = []
    mclist = []
    for f in output.keys():
        datalist.extend(
            [m for m in output[f].keys() if "Run" in m or "data" in m or "Data" in m]
        )
        mclist.extend(
            [
                m
                for m in output[f].keys()
                if "Run" not in m and "data" not in m and "Data" not in m
            ]
        )
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
    "njet" in arg.discr_list.split(",")
    or "nbjet" in arg.discr_list.split(",")
    or "mu" in arg.discr_list.split(",")
):
    nj = 1


for discr in arg.discr_list.split(","):
    if (
        "flav" in collated["mc"][discr].axes.name
        and "syst" in collated["mc"][discr].axes.name
    ):
        scale_sf = np.sum(
            collated["mc"][discr][{"syst": "SF", "flav": sum}].values()
        ) / np.sum(collated["mc"][discr][{"syst": "noSF", "flav": sum}].values())
    else:
        scale_sf = 1.0
    if (
        "flav" in collated["mc"][discr].axes.name
        and "syst" in collated["mc"][discr].axes.name
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
    fig.subplots_adjust(hspace=0.05)
    hep.cms.label("Preliminary", data=True, lumi=arg.lumi / 1000.0, loc=0, ax=ax)

    if (
        "flav" in collated["mc"][discr].axes.name
        and "syst" in collated["mc"][discr].axes.name
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
        stat_denom_unc = ratio_uncertainty(
            hdata.values(),
            collated["mc"][discr][{"syst": "noSF", "flav": sum}].values(),
        )
        ax.fill_between(
            hdata.axes.edges,
            np.ones(stat_denom_unc[0])
            - np.r_[stat_denom_unc[0], stat_denom_unc[0, -1]],
            np.ones(stat_denom_unc[0])
            + np.r_[stat_denom_unc[1], stat_denom_unc[1, -1]],
            {"facecolor": "tab:gray", "linewidth": 0},
        )
        ax.fill_between(
            hdata.axes.edges,
            np.ones(stat_denom_unc[0]) - np.r_[err_dn, err_dn[-1]],
            np.ones(stat_denom_unc[0]) + np.r_[err_up, err_up[-1]],
            {"facecolor": "tab:brown", "linewidth": 0},
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

    elif "mu" not in discr and "njet" != discr:
        hep.histplot(
            collated["mc"][discr],
            stack=True,
            label=["MC"],
            yerr=True,
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
    ax.set_xticklabels(ax.get_xticklabels(), fontsize=0)
    ax.set_ylabel("Events")
    rax.set_ylabel("Data/MC")
    rax.set_xlabel(discr)
    ax.legend()
    rax.set_ylim(0.5, 1.5)

    at = AnchoredText(
        input_txt + "\n",
        loc=2,
        frameon=False,
    )
    ax.add_artist(at)
    scale = ""
    if arg.norm:
        scale = "_norm"
    name = "all"
    # ax.set_ylim(0,500)
    # hep.mpl_magic(ax=ax)
    if arg.log:
        fig.savefig(
            f"plot/BTV/{arg.phase}_{arg.ext}_{time}/{arg.phase}_unc_{discr}_inclusive{scale}_{arg.ext}_{name}.pdf"
        )
    else:
        fig.savefig(
            f"plot/BTV/{arg.phase}_{arg.ext}_{time}/{arg.phase}_unc_lin_{discr}_inclusive{scale}_{arg.ext}_{name}.pdf"
        )
