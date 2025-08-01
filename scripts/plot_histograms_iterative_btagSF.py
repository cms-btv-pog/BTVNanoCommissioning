from pathlib import Path
import argparse
import warnings
import math

import matplotlib.pyplot as plt
from matplotlib.axes import Axes
import mplhep as hep
from coffea.util import load
from matplotlib import rc_context

from BTVNanoCommissioning.helpers.xs_scaler import collate, scaleSumW
from BTVNanoCommissioning.utils.plot_utils import MCerrorband, plotratio

parser = argparse.ArgumentParser(
    description="Plot histograms for iterative b-tagging SFs"
)
parser.add_argument("--suffix", type=str, help="suffix for input histograms", default="")
parser.add_argument(
    "--separate-bins",
    "-s",
    action="store_true",
    help="plot each eta and pt bin separately",
)
parser.add_argument(
    "--lumi",
    help="luminosity in pb^-1",
    type=float,
)
args = parser.parse_args()

lumi = args.lumi  # pb^-1
# need to devide, mplhep uses fb^-1
lumi_label = lumi / 1000  # fb^-1
com = 13.6
suffix = args.suffix
# if suffix != "":
suffix = f"_{suffix}"

# load the coffea file
# mc_path = "/net/data_cms3a-1/fzinn/BTV/btag_sf/nobackup/Summer23BPix/btag_ttbar_sf/hists_MC/hists_MC.coffea"
mc_path = Path(Path.cwd(), f"hists_MC{suffix}", f"hists_MC{suffix}.coffea")
# data_path = "/net/data_cms3a-1/fzinn/BTV/btag_sf/nobackup/Summer23BPix/btag_ttbar_sf/hists_data/hists_data.coffea"
data_path = Path(Path.cwd(), f"hists_data{suffix}", f"hists_data{suffix}.coffea")
print(f"Loading MC from {mc_path}")
print(f"Loading Data from {data_path}")

output = {
    "mc": load(mc_path),
    "data": load(data_path),
}

output = scaleSumW(output, lumi)
mclist = [m for m in output.keys() if "Run" not in m]
datalist = [m for m in output.keys() if "Run" in m]
mergemap = {"data": datalist, "mc": mclist}
collated = collate(output, mergemap)

variables = [
    "iterative_btagDeepFlavB",
    "iterative_btagPNetB",
    "iterative_btagRobustParTAK4B",
    "iterative_btagUParTAK4B",
]

# this would be to BTV Framework plotting style
# flav_labels = {"c": [4], "b": [5], "udsg": [0], "pu": [1]}
# color_list = [color_map[flav] for flav in flav_labels.keys()]

# this is with only b,c,l
# flavor -> index
# {
#     "b": "5",
#     "c": "4",
#     "l": "<4 or >5",
# }
# flav_axis has 0,1,4,5,6
flav_labels = {
    "c": [4],
    "b": [5],
    "l": [0, 1, 6],
}
COLORS = {
    "b": "#5790fc",
    "c": "#f89c20",
    "l": "#964a8b",
}
color_list = [COLORS[flav] for flav in flav_labels.keys()]


def set_yaxis_limit_ratio(ax: Axes) -> None:
    """Set the y-axis limit for the ratio plot.

    Sets the y-axis limit between 0 and 2,
    Takes the current y-axis limit and sets:
    - the lower limit to 0 if it is below 0
    - the upper limit rounded up to the nearest 0.1
        if it is above 2, set to 2

    :param ax: the axes to set the limit for
    :type ax: Axes
    """
    ax.set_ylim(
        max(0, ax.get_ylim()[0]),
        min(2, math.ceil(10 * ax.get_ylim()[1]) / 10),
    )


with rc_context(hep.style.CMS) as cms, warnings.catch_warnings() as catch_warnings:
    warnings.simplefilter("ignore", category=UserWarning)
    for variable in variables:
        try:
            histogram = collated["mc"][variable]
            data_histogram = collated["data"][variable]
        except KeyError as e:
            warnings.warn(f"Variable {variable} not in the files: {e}")
        # syst, flav, eta, pt, region, jet_index, discr

        for region in histogram.axes["region"]:
            flav_axis = histogram.axes["flav"]

            hist_list = [
                histogram[
                    "nominal", list(flav_axis.index(flav)), sum, sum, region, sum, :
                ][sum, :]
                for flav in flav_labels.values()
            ]

            # sum of MC histograms
            summed_hist = histogram["nominal", sum, sum, sum, region, sum, :]
            fig, (ax, ratio_ax) = plt.subplots(
                2,
                1,
                figsize=(10, 10),
                gridspec_kw={"height_ratios": (3, 1)},
                sharex=True,
            )
            fig.subplots_adjust(hspace=0.06, top=0.92, bottom=0.1, right=0.97)
            hep.cms.label(
                "Preliminary", com=com, lumi=lumi_label, ax=ax, loc=0, data=True
            )
            hep.histplot(
                summed_hist,
                label="MC",
                histtype="fill",
                yerr=True,
                ax=ax,
                color="#f89c20",
            )
            hep.histplot(
                data_histogram["nominal", 0, sum, sum, region, sum, :],
                label="Data",
                histtype="errorbar",
                color="black",
                yerr=True,
                ax=ax,
            )
            MCerrorband(summed_hist, ax=ax)
            ax.set_xlim(0, 1)
            ax.set_xlabel(None)
            ax.set_ylabel("Events")
            # plot ratio
            rax = plotratio(
                data_histogram["nominal", 0, sum, sum, region, sum, :],
                summed_hist,
                ax=ratio_ax,
            )
            rax.set_ylabel("Data/MC")
            set_yaxis_limit_ratio(rax)

            # title
            fig.suptitle(f"{region=}")

            save_path = Path("plots_iterative_btagSF", suffix, variable)
            save_path.mkdir(parents=True, exist_ok=True)
            fig.tight_layout()
            fig.savefig(save_path / f"{region}_MC_sum.png")
            ax.set_yscale("log")
            fig.savefig(save_path / f"{region}_MC_sum_log.png")
            plt.close(fig)

            # MC histograms stacked
            # plot sum over all eta and pt bins
            fig, (ax, ratio_ax) = plt.subplots(
                2,
                1,
                figsize=(10, 10),
                gridspec_kw={"height_ratios": (3, 1)},
                sharex=True,
            )
            fig.subplots_adjust(hspace=0.06, top=0.92, bottom=0.1, right=0.97)
            hep.cms.label(
                "Preliminary", com=com, lumi=lumi_label, ax=ax, loc=0, data=True
            )
            hep.histplot(
                hist_list,
                label=list(flav_labels.keys()),
                stack=True,
                histtype="fill",
                yerr=True,
                ax=ax,
                color=color_list,
            )
            hep.histplot(
                data_histogram["nominal", 0, sum, sum, region, sum, :],
                label="Data",
                histtype="errorbar",
                color="black",
                yerr=True,
                ax=ax,
            )
            MCerrorband(histogram["nominal", sum, sum, sum, region, sum, :], ax=ax)
            ax.set_xlim(0, 1)
            ax.set_xlabel(None)
            ax.set_ylabel("Events")
            ax.legend()

            # plot ratio
            rax = plotratio(
                data_histogram["nominal", 0, sum, sum, region, sum, :],
                histogram["nominal", sum, sum, sum, region, sum, :],
                ax=ratio_ax,
            )
            set_yaxis_limit_ratio(rax)
            rax.set_ylabel("Data/MC")

            fig.suptitle(f"{region=}")

            save_path = Path("plots_iterative_btagSF", suffix, variable)
            save_path.mkdir(parents=True, exist_ok=True)
            fig.tight_layout()
            fig.savefig(save_path / f"{region}.png")
            ax.set_yscale("log")
            fig.savefig(save_path / f"{region}_log.png")
            plt.close(fig)

            # MC histograms separately
            fig, ax = plt.subplots(
                1,
                1,
                figsize=(10, 8),
            )
            fig.subplots_adjust(hspace=0.06, top=0.92, bottom=0.1, right=0.97)
            hep.cms.label(
                "Preliminary", com=com, lumi=lumi_label, ax=ax, loc=0, data=True
            )
            hep.histplot(
                hist_list,
                label=list(flav_labels.keys()),
                stack=False,
                histtype="step",
                yerr=False,
                ax=ax,
                color=color_list,
                lw=2,
            )
            # MCerrorband(histogram["nominal", sum, sum, sum, region, sum, :], ax=ax)
            ax.set_xlim(0, 1)
            ax.set_ylim(bottom=max(1e-1, ax.get_ylim()[0]))
            ax.set_xlabel(None)
            ax.set_ylabel("Events")
            ax.legend()

            fig.suptitle(f"{region=}")
            fig.tight_layout()

            fig.savefig(save_path / f"{region}_hists.png")
            ax.set_yscale("log")
            fig.savefig(save_path / f"{region}_hists_log.png")
            plt.close(fig)

            if not args.separate_bins:
                continue

            save_path_split = save_path / "eta_pt_bins"
            save_path_split.mkdir(parents=True, exist_ok=True)
            for eta_index, eta in enumerate(histogram.axes["eta"]):
                eta_bin = f"eta{eta[0]:.1f}to{eta[1]:.1f}"
                for pt_index, pt in enumerate(histogram.axes["pt"]):
                    pt_bin = f"pt{pt[0]:.0f}to{pt[1]:.0f}"

                    # plot for each eta and pt bin
                    fig, (ax, ratio_ax) = plt.subplots(
                        2,
                        1,
                        figsize=(10, 10),
                        gridspec_kw={"height_ratios": (3, 1)},
                        sharex=True,
                    )
                    fig.subplots_adjust(hspace=0.06, top=0.92, bottom=0.1, right=0.97)
                    hep.cms.label(
                        "Preliminary", com=com, lumi=lumi_label, ax=ax, loc=0, data=True
                    )
                    hep.histplot(
                        [
                            histogram[
                                "nominal",
                                list(flav_axis.index(flav)),
                                # index,
                                eta_index,
                                pt_index,
                                region,
                                sum,
                                :,
                            ][sum, :]
                            for flav in flav_labels.values()
                            # for index in range(len(flav_axis) - 1)
                        ],
                        label=flav_labels,
                        stack=True,
                        histtype="fill",
                        yerr=True,
                        ax=ax,
                        color=color_list,
                    )
                    hep.histplot(
                        data_histogram[
                            "nominal", 0, eta_index, pt_index, region, sum, :
                        ],
                        label="Data",
                        histtype="errorbar",
                        color="black",
                        yerr=True,
                        ax=ax,
                    )
                    ax.set_xlim(0, 1)
                    ax.legend()

                    # plot ratio
                    rax = plotratio(
                        data_histogram[
                            "nominal", 0, eta_index, pt_index, region, sum, :
                        ],
                        histogram["nominal", sum, eta_index, pt_index, region, sum, :],
                        ax=ratio_ax,
                    )
                    set_yaxis_limit_ratio(rax)
                    rax.set_ylabel("Data/MC")
                    fig.suptitle(f"{region=}")
                    fig.tight_layout()

                    # save
                    fig.savefig(save_path_split / f"{region}_{eta_bin}_{pt_bin}.png")
                    plt.close(fig)
