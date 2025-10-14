import numpy as np
import argparse
import hist as Hist
import os

import matplotlib.pyplot as plt
import mplhep as hep

from itertools import product
from BTVNanoCommissioning.helpers.xs_scaler import scaleSumW
from BTVNanoCommissioning.utils.plot_utils import plotratio, MCerrorband
from qghelpers import load_inputs, get_slices, stack_hists, plot_roc

plt.style.use(hep.style.ROOT)

color_map = {
    "$t\\bar{t}$": "#008695",
    "Single top": "#3969ac",
    "VV": "#7f3c8d",
    "Z+jets": "#e68310",
    "W+jets": "#80ba5a",
    "QCD": "#e73f74",
    "QCD($\\mu$)": "#a5aa99",
    "Ph+jets (HT-binned)": "#f2ab6d",
    "Ph+jets (HT+PTG-binned)": "tab:brown",
    "udsg": "tab:blue",
    "ud": "tab:blue",
    "s": "tab:olive",
    "g": "tab:purple",
    "pu": "tab:orange",
    "c": "tab:green",
    "b": "tab:red",
    "other": "tab:gray",
}

parser = argparse.ArgumentParser("Plotter for QG workflow outputs")
parser.add_argument(
    "-l", "--lumi", type=float, default=1.0, help="Lumi to scale MC to (in /pb)"
)
parser.add_argument(
    "-i", "--input", type=str, required=True, help="Input files", nargs="+"
)
parser.add_argument(
    "-o", "--output", type=str, required=True, help="Output directory"
)
parser.add_argument(
    "-v", "--var", type=str, default=["all"], help="Variables to plot", nargs="+"
)
parser.add_argument(
    "-s", "--split", type=str, default=["flav"], help="Split by 'flav' or 'sample'", choices=["flav", "sample"], nargs="+"
)

ptr_parser = parser.add_mutually_exclusive_group()
ptr_parser.add_argument(
    "--ptrange", type=str, default=None, help="pT range to select in data coordinates, e.g. '20,50'"
)
ptr_parser.add_argument(
    "--iptrange", type=str, default=None, help="pT range to select in bin coordinates, e.g. '0,20'"
)

etar_parser = parser.add_mutually_exclusive_group()
etar_parser.add_argument(
    "--etarange", type=str, default=None, help="eta range to select in data coordinates, e.g. '0,2.5'",
)
etar_parser.add_argument(
    "--ietarange", type=str, default=None, help="eta range to select in bin coordinates, e.g. '0,5'"
)

args = parser.parse_args()

# Load input files and merge histograms
input_hists, available_vars = load_inputs(args.input)

# Scale MC to lumi
scaled_hists = scaleSumW(input_hists, args.lumi)

vars_to_plot = set()
if "all" not in args.var:
    for var in args.var:
        if var not in available_vars:
            print(f"Variable {var} not found in input files, available: {available_vars}")
            exit(1)
        vars_to_plot.add(var)
else:
    vars_to_plot = available_vars

vars_to_remove = set()
for v in vars_to_plot:
    if "Tag" in v:
        print(f"Variable {v} is a Tag variable, skipping as not implemented")
        vars_to_remove.add(v)
    if "Obj" not in v or "Var" not in v:
        print(f"Variable {v} does not follow naming convention, skipping")
        vars_to_remove.add(v)

vars_to_plot -= vars_to_remove

print(f"Variables to plot: {vars_to_plot}")

# Select pT/eta ranges if given
pt_slices = []
if args.ptrange or args.iptrange:
    pt_slices = get_slices(args.ptrange if args.ptrange else args.iptrange, bincoords=bool(args.iptrange))
    print(f"Selecting pT slices: {pt_slices}")

eta_slices = []
if args.etarange or args.ietarange:
    eta_slices = get_slices(args.etarange if args.etarange else args.ietarange, bincoords=bool(args.ietarange))
    print(f"Selecting eta slices: {eta_slices}")

slices = [{"pt": ps, "eta": es} for ps in (pt_slices if pt_slices else [slice(None)]) for es in (eta_slices if eta_slices else [slice(None)])]

hists = stack_hists(scaled_hists, vars_to_plot, regions=slices)

var_split = list(product(list(hists.keys()), args.split))

if not os.path.exists(args.output):
    os.makedirs(args.output)

for var, split_by in var_split:
    if "_pteta" in var and "Tag" not in var:
        plot_roc(hists[var], systname="nominal", var=var)

    ptrange_str = "$p_{T}$: " + var.split("PT")[1].split("_")[0] if "PT" in var else ""
    etarange_str = "$\\eta$: " + var.split("ETA")[1].split("_")[0] if "ETA" in var else ""
    if ptrange_str:
        ptrange_str = ptrange_str + " GeV"
    if etarange_str:
        etarange_str = etarange_str


    fig, ((ax), (rax)) = plt.subplots(
        2,
        1,
        figsize=(10,10),
        gridspec_kw={"height_ratios": (3, 1)},
        sharex=True,
    )
    # fig.subplots_adjust(hspace=0.06, top=0.92, bottom=0.1, right=0.97)
    hep.cms.label(
        "Private", data=True, lumi=args.lumi / 1000, com=13.6, loc=0, ax=ax
    )

    dhists = [hists[var]["data"][s]["nominal",:] for s in hists[var]["data"]]
    hep.histplot(
        dhists,
        histtype="errorbar",
        color="black",
        # stack=True,
        label="Data",
        yerr=True,
        ax=ax,
        flow="none",
    )

    mc_hists = hists[var]["mc_flav"] if split_by == "flav" else hists[var]["mc_sample"]
    
    hep.histplot(
        [mc_hists[s]["nominal",:] for s in mc_hists],
        histtype="fill",
        stack=True,
        ax=ax,
        label=[k for k in mc_hists],
        sort="y",
        color=[color_map.get(k, None) for k in mc_hists],
        flow="none",
    )

    rax = plotratio(
        sum(dhists),
        sum([mc_hists[s]["nominal",:] for s in mc_hists]),
        ax=rax,
    )

    rax = MCerrorband(
        sum([mc_hists[s]["nominal",:] for s in mc_hists]),
        ax=rax,
    )

    std_legend = ax.legend(ncol=2)
    ax.add_artist(std_legend)

    # ptartist = ax.scatter([], [], marker="$p_{T}$", color="black")
    # etaartist = ax.scatter([], [], marker="$|\\eta|$", color="black")
    # ax.legend(
        # handles=[ptartist, etaartist],
        # labels=[ptrange_str.replace("j", ""), etarange_str.replace("j", "")],
        # loc="lower right",
    # )
    # print(ptrange_str, etarange_str)

    ax.set_xlabel(None)
    ax.set_ylabel("Events")

    ax.text(
        0.95,
        0.2,
        f"{ptrange_str.replace('j', '')}\n{etarange_str.replace('j', '')}",
        horizontalalignment="right",
        verticalalignment="top",
        transform=ax.transAxes,
        fontsize=12,
        bbox=dict(facecolor="white", edgecolor="black", boxstyle="round,pad=0.3"),
    )

    rax.set_ylabel("Data/MC")
    rax.set_ylim(0.5, 1.5)
    rax.hlines([1], xmin=rax.get_xlim()[0], xmax=rax.get_xlim()[1], color="black", ls="--")
    rax.hlines([0.8, 1.2], xmin=rax.get_xlim()[0], xmax=rax.get_xlim()[1], color="black", ls=":")

    ax.set_xlim(rax.get_xlim())
    ax.ticklabel_format(style="sci", scilimits=(-4, 4))
    ax.get_yaxis().get_offset_text().set_position((-0.065, 1.05))

    plt.savefig(f"{args.output}/{var}_{split_by}.png")
    plt.close()

    # Plot the flavour fractions
    if split_by == "flav":
        fig, ax = plt.subplots(figsize=(10,7))
        fig.subplots_adjust(top=0.92, bottom=0.1, right=0.97)
        hep.cms.label(
            "Private", data=False, com=13.6, loc=0, ax=ax
        )

        mc_hists = hists[var]["mc_flav"]
        total_mc = sum([mc_hists[s]["nominal",:] for s in mc_hists]).values()
        frac_hists = {s: mc_hists[s]["nominal",:] * 1/ total_mc for s in mc_hists}

        hep.histplot(
            [frac_hists[s] for s in frac_hists],
            histtype="fill",
            stack=True,
            ax=ax,
            label=[k for k in frac_hists],
            sort="y",
            color=[color_map.get(k, None) for k in frac_hists],
            flow="none",
        )

        ax.legend(ncol=3)
        ax.set_xlabel(var.split("_Var")[1])
        ax.set_ylabel("Flavour fraction")
        ax.set_ylim(0,1.3)
        ax.hlines(1, xmin=0, xmax=1.0, color="black", ls="--")
        if "Varpt" in var:
            ax.set_xlim(30,80)
        elif "Vareta" in var:
            ax.set_xlim(0,2.5)
        elif "_pteta" in var:
            ax.set_xlim(0,1.0)

        plt.savefig(f"{args.output}/{var}_flav_fraction.png")
        plt.close()

    