import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import uproot, sys, os, argparse, hist
from coffea.util import load
from BTVNanoCommissioning.utils.plot_utils import plotratio
from matplotlib.offsetbox import AnchoredText

parser = argparse.ArgumentParser(description="Get validation plots(ROC&efficiencies)")
parser.add_argument(
    "-i",
    "--input",
    type=str,
    required=True,
    help="Input coffea file(s)",
)

parser.add_argument(
    "-v",
    "--version",
    type=str,
    required=True,
    help="version name",
)
tagger = ["btagDeepFlav", "btagPNet", "btagRobustParTAK4"]
color = ["tab:blue", "tab:orange", "tab:green"]
if __name__ == "__main__":
    args = parser.parse_args()
    output = load(args.input)
    output = output[list(output.keys())[0]]
    ## b-jet ROC
    colors = dict(zip(tagger, color))
    plt.style.use(hep.style.ROOT)
    fig, ax = plt.subplots()
    for tag in tagger:
        axis = {"syst": "noSF", "flav": 3}  # consider b-jet only
        ROC = [
            np.sum(output[f"{tag}B_0"][axis].values()[i:])
            / np.sum(output[f"{tag}B_0"][axis].values())
            for i in range(len(output[f"{tag}B_0"][axis].values()))
        ]
        ax.plot(
            np.linspace(0, 1, 50), ROC, label=tag.replace("btag", ""), color=colors[tag]
        )
    ax.legend()
    ax.set_ylabel("b-jet efficiency")
    ax.set_xlabel("tagger score")
    ax.set_title(args.version)
    ax.set_ylim(top=1.0)
    ax.set_ylim(0, 1)
    at = AnchoredText("", loc=2, frameon=False)
    ax.add_artist(at)
    hep.mpl_magic(ax=ax)
    fig.savefig(f"{args.version}_bjet_ROC.pdf")
    fig, ax = plt.subplots()
    ## c-jet ROC
    for tag in tagger:
        axis = {"syst": "noSF", "flav": 2}  # consider c-jet only
        ROC_CvL = [
            np.sum(output[f"{tag}CvL_0"][axis].values()[i:])
            / np.sum(output[f"{tag}CvL_0"][axis].values())
            for i in range(len(output[f"{tag}CvL_0"][axis].values()))
        ]
        ROC_CvB = [
            np.sum(output[f"{tag}CvB_0"][axis].values()[i:])
            / np.sum(output[f"{tag}CvB_0"][axis].values())
            for i in range(len(output[f"{tag}CvB_0"][axis].values()))
        ]
        ax.plot(
            np.linspace(0, 1, 50),
            ROC_CvL,
            label=tag.replace("btag", "") + " CvL",
            color=colors[tag],
        )
        ax.plot(
            np.linspace(0, 1, 50),
            ROC_CvB,
            label=tag.replace("btag", "") + " CvB",
            color=colors[tag],
            ls=":",
        )
    ax.legend()
    at = AnchoredText("", loc=2, frameon=False)
    ax.add_artist(at)
    ax.set_ylabel("c-jet efficiency")
    ax.set_xlabel("tagger score")
    ax.set_ylim(top=1.0)
    ax.set_xlim(0, 1)
    hep.mpl_magic(ax=ax)
    ax.set_title(args.version)
    fig.savefig(f"{args.version}_cjet_ROC.pdf")
    ## Efficiency plot of kinmetic variables
    color = ["r", "g", "b", "c", "m"]
    for j in ["b", "c"]:
        for var in ["pt", "eta", "phi"]:
            for t in range(output[f"{j}jet_WP_{var}"].axes["tagger"].size):
                fig, ax = plt.subplots()
                for wp in range(1, output[f"{j}jet_WP_{var}"].axes["WP"].size):
                    label = output[f"{j}jet_WP_{var}"].axes["WP"].value(wp)
                    if wp == 1:
                        plotratio(
                            output[f"{j}jet_WP_{var}"][wp, t, :],
                            output[f"{j}jet_WP_{var}"][0, t, :],
                            denom_fill_opts=None,
                            error_opts={"color": color[wp - 1]},
                            ax=ax,
                            label=label,
                        )  # ,denom_fill_opts=None,label=label,ax=ax)
                    else:
                        plotratio(
                            output[f"{j}jet_WP_{var}"][wp, t, :],
                            output[f"{j}jet_WP_{var}"][0, t, :],
                            denom_fill_opts=None,
                            clear=False,
                            ax=ax,
                            error_opts={"color": color[wp - 1]},
                            label=label,
                        )
                ax.set_xlabel(
                    "jet " + output[f"{j}jet_WP_{var}"][0, t, :].axes[0].label
                )
                ax.set_ylabel(f"{j}-jet efficiency")
                ax.set_title(args.version)
                at = AnchoredText("", loc=2, frameon=False)
                ax.add_artist(at)
                ax.set_title(
                    output[f"{j}jet_WP_{var}"]
                    .axes["tagger"]
                    .value(t)
                    .replace("btag", "")
                )
                ax.legend()
                hep.mpl_magic(ax=ax)
                fig.savefig(
                    f'{args.version}_{output[f"{j}jet_WP_{var}"].axes["tagger"].value(t)}_jet{var}_{j}jet_eff.pdf'
                )
