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

markers = [".", "o", "^", "s", "+", "x", "D", "*"]
parser = argparse.ArgumentParser(description="hist plotter for commissioning")
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

parser.add_argument("--ext", type=str, default="data", help="addional name")
parser.add_argument("-o", "--output", required=True, type=str, help="files set")
parser.add_argument("-r", "--ref", required=True, help="referance dataset")
parser.add_argument(
    "-c",
    "--compared",
    required=True,
    type=str,
    help="compared dataset",
)
parser.add_argument("--sepflav", default=False, type=bool, help="seperate flavour")
parser.add_argument("--log", action="store_true", help="log on x axis")
parser.add_argument(
    "-d",
    "--discr_list",
    nargs="+",
    default=[
        "btagDeepCvL",
        "btagDeepCvB",
        "btagDeepFlavCvL",
        "btagDeepFlavCvB",
        "btagDeepB_b",
        "btagDeepFlavB",
    ],
    help="discriminators",
)
parser.add_argument("--ext", type=str, default="data", help="addional name")


args = parser.parse_args()
output = {}
if len(args.output.split(",")) > 1:
    output = {i: load({args.output.split(",")[i]}) for i in args.output.split(",")}
else:
    output = load(args.output)
mergemap = {}
time = arrow.now().format("YY_MM_DD")
if not os.path.isdir(f"plot/BTV/{args.phase}_{args.ext}_{time}/"):
    os.makedirs(f"plot/BTV/{args.phase}_{args.ext}_{time}/")
if "sumw" in output.keys():
    mergemap[args.ref] = [m for m in output.keys() if args.ref == m]
    for c in args.compared.split(","):
        mergemap[c] = [m for m in output.keys() if c == m]
else:
    reflist = []
    comparelist = []
    for f in output.keys():
        reflist.extend([m for m in output[f].keys() if args.ref == m])
        for c in args.compared.split(","):
            comparelist.extend([m for m in output[f].keys() if c == m])
    mergemap[args.ref] = reflist
    mergemap[c] = comparelist
collated = collate(output, mergemap)
### style settings
if "Run" in args.ref or "data" in args.ref or "Data" in args.ref:
    hist_type = "errorbar"
else:
    hist_type = "step"

if "ttdilep" in args.phase:
    input_txt = "dilepton ttbar"
    nj = 2
elif "ttsemilep" in args.phase:
    input_txt = "semileptonic ttbar"
    nj = 4
else:
    if "Wc" in args.phase:
        input_txt = "W+c"
    elif "DY" in args.phase:
        input_txt = "DY+jets"
    nj = 1


for discr in args.discr_list:
    if args.sepflav:  # split into 3 panels for different flavor
        fig, (ax, rax, rax2, rax3) = plt.subplots(
            4,
            1,
            figsize=(8, 8),
            gridspec_kw={"height_ratios": (3, 1, 1, 1)},
            sharex=True,
        )
        fig.subplots_adjust(hspace=0.07)
        ax.set_xlabel(None)
        ax.set_xticklabels(ax.get_xticklabels(), fontsize=0)
        hep.cms.label(
            "Preliminary",
            data=True,
            loc=0,
            ax=ax,
        )
        laxis = {"flav": 0}
        puaxis = {"flav": 1}
        caxis = {"flav": 2}
        baxis = {"flav": 3}

        if "syst" in collated[args.ref][discr].axes.name:
            laxis["syst"] = "noSF"
            puaxis["syst"] = "noSF"
            caxis["syst"] = "noSF"
            baxis["syst"] = "noSF"

        hep.histplot(
            collated[args.ref][discr][laxis] + collated[args.ref][discr][puaxis],
            label=args.ref + "-l",
            color="b",
            histtype=hist_type,
            yerr=True,
            ax=ax,
        )
        hep.histplot(
            collated[args.ref][discr][caxis],
            label=args.ref + "-c",
            color="g",
            histtype=hist_type,
            yerr=True,
            ax=ax,
        )
        hep.histplot(
            collated[args.ref][discr][baxis],
            label=args.ref + "-b",
            yerr=True,
            color="r",
            histtype=hist_type,
            ax=ax,
        )

        for c in args.compared.split(","):

            hep.histplot(
                collated[c][discr][laxis] + collated[c][discr][puaxis],
                label=c + "-l",
                color="b",
                marker=markers[c + 1],
                histtype="errorbar",
                yerr=True,
                ax=ax,
            )
            hep.histplot(
                collated[c][discr][caxis],
                label=c + "-c",
                color="g",
                marker=markers[c + 1],
                histtype="errorbar",
                yerr=True,
                ax=ax,
            )
            hep.histplot(
                collated[c][discr][baxis],
                label=c + "-b",
                yerr=True,
                color="r",
                marker=markers[c + 1],
                histtype="errorbar",
                ax=ax,
            )
        ax.legend(
            ncol=3,
            loc="upper right",
        )

        for c in args.compared.split(","):

            rax.errorbar(
                x=collated[c][discr][laxis].axes[0].centers,
                y=(
                    collated[c][discr][laxis].values()
                    + collated[c][discr][puaxis].values()
                )
                / (
                    collated[args.ref][discr][laxis].values()
                    + collated[args.ref][discr][puaxis].values()
                ),
                yerr=ratio_uncertainty(
                    (
                        collated[c][discr][laxis].values()
                        + collated[c][discr][puaxis].values()
                    ),
                    collated[args.ref][discr][laxis].values()
                    + collated[args.ref][discr][puaxis].values(),
                ),
                color="b",
                linestyle="none",
                marker=markers[c + 1],
            )
            rax2.errorbar(
                x=collated[c][discr][caxis].axes[0].centers,
                y=collated[c][discr][caxis].values()
                / collated[args.ref][discr][caxis].values(),
                yerr=ratio_uncertainty(
                    collated[c][discr][caxis].values(),
                    collated[args.ref][discr][caxis].values(),
                ),
                color="g",
                linestyle="none",
                marker=markers[c + 1],
            )
            rax3.errorbar(
                x=collated[c][discr][baxis].axes[0].centers,
                y=collated[c][discr][baxis].values()
                / collated[args.ref][discr][baxis].values(),
                yerr=ratio_uncertainty(
                    collated[c][discr][baxis].values(),
                    collated[args.ref][discr][baxis].values(),
                ),
                color="r",
                linestyle="none",
                marker=markers[c + 1],
            )

        discrs = discr
        ax.set_xlabel("A.U.")
        rax3.set_xlabel(discrs)
        rax.set_ylabel("udsg-jets")
        rax2.set_ylabel("c-jets")
        rax3.set_ylabel("b-jets")
        rax.set_ylim(0.5, 1.5)
        rax2.set_ylim(0.5, 1.5)
        rax3.set_ylim(0.5, 1.5)
        rax3.set_xlabel(discr)
        ax.legend()
        at = AnchoredText(
            "",
            # + "inclusive pT, $\eta$"
            loc=2,
            prop=dict(size=15),
            frameon=False,
        )
        ax.add_artist(at)
        hep.mpl_magic(ax=ax)
        fig.savefig(
            f"plot/BTV/{args.phase}_{args.ext}_{time}/compare_{args.phase}_inclusive{discrs}.png"
        )
        fig.savefig(
            f"plot/BTV/{args.phase}_{args.ext}_{time}/compare_{args.phase}_inclusive{discrs}.pdf"
        )

    else:
        fig, ((ax), (rax)) = plt.subplots(
            2, 1, figsize=(8, 8), gridspec_kw={"height_ratios": (3, 1)}, sharex=True
        )
        fig.subplots_adjust(hspace=0.07)
        hep.cms.label(
            "Preliminary",
            data=True,
            loc=0,
            ax=ax,
        )
        ax.set_xlabel(None)
        ax.set_xticklabels(ax.get_xticklabels(), fontsize=0)
        allaxis = {}
        if "flav" in collated[args.ref][discr].axes.name:
            allaxis["flav"] = sum
        if "syst" in collated[args.ref][discr].axes.name:
            allaxis["syst"] = sum
        hep.histplot(
            collated[args.ref][discr][allaxis],
            label=args.ref,
            histtype=hist_type,
            yerr=True,
            ax=ax,
        )
        for c in args.compared.split(","):
            hep.histplot(
                collated[c][discr][allaxis],
                label=c,
                histtype=hist_type,
                yerr=True,
                ax=ax,
            )

        for c in args.compared.split(","):
            rax.errorbar(
                x=collated[c][discr][allaxis].axes[0].centers,
                y=collated[c][discr][allaxis].values()
                / collated[args.ref][discr][allaxis].values(),
                yerr=ratio_uncertainty(
                    collated[c][discr][allaxis].values(),
                    collated[args.ref][discr][allaxis].values(),
                ),
                color="k",
                linestyle="none",
                marker="o",
                elinewidth=1,
            )
        rax.set_xlabel(discr)
        ax.set_xlabel(None)
        ax.set_ylabel("Events")
        rax.set_ylabel("Other/Ref")
        ax.legend()
        rax.set_ylim(0.0, 2.0)

        at = AnchoredText(
            "",
            loc=2,
            frameon=False,
        )
        ax.add_artist(at)
        hep.mpl_magic(ax=ax)
        fig.savefig(
            f"plot/BTV/{args.phase}_{args.ext}_{time}/compare_{args.phase}_lin_inclusive{discr}.pdf"
        )
        fig.savefig(
            f"plot/BTV/{args.phase}_{args.ext}_{time}/compare_{args.phase}_lin_inclusive{discr}.png"
        )
