import numpy as np
import argparse, os, arrow, glob
from coffea.util import load
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import mplhep as hep
import hist

plt.style.use(hep.style.ROOT)
from BTVNanoCommissioning.workflows import workflows
from BTVNanoCommissioning.utils.xs_scaler import collate, scaleSumW
from BTVNanoCommissioning.helpers.definitions import definitions, axes_name
from BTVNanoCommissioning.utils.plot_utils import (
    plotratio,
    SFerror,
    MCerrorband,
    autoranger,
    rebin_hist,
)

bininfo = definitions()
parser = argparse.ArgumentParser(description="hist plotter for commissioning")
parser.add_argument("--lumi", required=True, type=float, help="luminosity in /pb")
parser.add_argument("--com", default="13", type=str, help="sqrt(s) in TeV")
parser.add_argument(
    "-p",
    "--phase",
    required=True,
    choices=list(workflows.keys()),
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
    "--xlabel",
    type=str,
    help="rename the label of variables to plot, splitted by ,.  Wildcard option * NOT available here",
)
parser.add_argument(
    "--ylabel", type=str, default=None, help="Modify y-axis label of plot"
)
parser.add_argument(
    "--SF", action="store_true", default=False, help="make w/, w/o SF comparisons"
)
parser.add_argument(
    "--xrange", type=str, default=None, help="custom x-range, --xrange xmin,xmax"
)
parser.add_argument("--ext", type=str, default="", help="prefix name")
parser.add_argument(
    "--autorebin",
    default=None,
    help="Rebin the plotting variables, input `int` or `list`. int: merge N bins. list of number: rebin edges(non-uniform bin is possible)",
)
parser.add_argument(
    "--flow",
    type=str,
    default=None,
    help="str, optional {None, 'show', 'sum'} Whether plot the under/overflow bin. If 'show', add additional under/overflow bin. If 'sum', add the under/overflow bin content to first/last bin.",
)
parser.add_argument(
    "--splitOSSS",
    type=int,
    default=None,
    help="Only for W+c phase space, split opposite sign(1) and same sign events(-1), if not specified, the combined OS-SS phase space is used",
)


arg = parser.parse_args()
time = arrow.now().format("YY_MM_DD")
if not os.path.isdir(f"plot/BTV/{arg.phase}_{arg.ext}_{time}/"):
    os.makedirs(f"plot/BTV/{arg.phase}_{arg.ext}_{time}/")
if len(arg.input.split(",")) > 1:
    output = {i: load(i) for i in arg.input.split(",")}
elif "*" in arg.input:
    files = glob.glob(arg.input)
    output = {i: load(i) for i in files}
else:
    output = {arg.input: load(arg.input)}
output = scaleSumW(output, arg.lumi)
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
### input text settings
if "Wc" in arg.phase:
    input_txt = "W+c"
    if arg.splitOSSS == 1:
        input_txt = input_txt + " OS"
    elif arg.splitOSSS == -1:
        input_txt = input_txt + " SS"
    else:
        input_txt = input_txt + " OS-SS"
elif "DY" in arg.phase:
    input_txt = "DY+jets"
elif "semilep" in arg.phase:
    input_txt = r"t$\bar{t}$ semileptonic"
    nj = 4
elif "dilep" in arg.phase:
    input_txt = r"t$\bar{t}$ dileptonic"
    nj = 2
if (
    "njet" in arg.variable.split(",")
    or "nbjet" in arg.variable.split(",")
    or "mu" in arg.variable.split(",")
):
    nj = 1
if "emctag" in arg.phase:
    input_txt = input_txt + " (e$\mu$)"
elif "ectag" in arg.phase:
    input_txt = input_txt + " (e)"
elif "ttdilep_sf" == arg.phase:
    input_txt = input_txt + " (e$\mu$)"
else:
    input_txt = input_txt + " ($\mu$)"
if "ctag" in arg.phase and "DY" not in arg.phase:
    input_txt = input_txt + "\nw/ soft-$\mu$"
if arg.variable == "all":
    var_set = collated["mc"].keys()
elif "*" in arg.variable:
    var_set = [
        var for var in collated["mc"].keys() if arg.variable.replace("*", "") in var
    ]
else:
    var_set = arg.variable.split(",")

for index, discr in enumerate(var_set):
    if "sumw" == discr:
        continue
    allaxis = {}
    if "Wc" in arg.phase:
        if arg.splitOSSS is None:  # OS-SS
            collated["mc"][discr] = (
                collated["mc"][discr][{"osss": 0}]
                + collated["mc"][discr][{"osss": 1}] * -1
            )
            collated["data"][discr] = (
                collated["data"][discr][{"osss": 0}]
                + collated["data"][discr][{"osss": 1}] * -1
            )
        elif arg.splitOSSS == 1:
            allaxis["osss"] = 0  # opposite sign
        elif arg.splitOSSS == -1:
            allaxis["osss"] = 1  # same sign
    if "flav" in collated["mc"][discr].axes.name:
        allaxis["flav"] = sum
        SF_axis = allaxis
        noSF_axis = allaxis
    if "syst" in collated["mc"][discr].axes.name:
        allaxis["syst"] = "nominal"
        SF_axis = allaxis
        noSF_axis = allaxis
        systlist = [i for i in range(collated["mc"][discr].axes[0].size)]
        if "noSF" in systlist:
            noSF_axis["syst"] = "noSF"

    do_xerr = False
    if arg.autorebin is not None:
        if arg.autorebin.isdigit():
            rebin = int(arg.autorebin)
        else:
            rebin = np.array([float(i) for i in arg.autorebin.split(",")])
            do_xerr = True
        collated["mc"][discr] = rebin_hist(
            collated["mc"][discr], collated["mc"][discr].axes[-1].name, rebin
        )
        collated["data"][discr] = rebin_hist(
            collated["data"][discr], collated["data"][discr].axes[-1].name, rebin
        )

    if (
        "flav" in collated["mc"][discr].axes.name
        and "syst" in collated["mc"][discr].axes.name
        and arg.SF
    ):
        scale_sf = np.sum(collated["mc"][discr][SF_axis].values()) / np.sum(
            collated["mc"][discr][noSF_axis].values()
        )
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
        splitflav_stack = []
        splitflav_axis = SF_axis
        for i in range(4):
            splitflav_axis["flav"] = i
            splitflav_stack += [[collated["mc"][discr][splitflav_axis]]]

        SF_axis["flav"] = sum
        hep.histplot(
            splitflav_stack,
            stack=True,
            label=["udsg", "pileup", "c", "b"],
            histtype="fill",
            yerr=True,
            ax=ax,
            # flow=arg.flow
        )
        hep.histplot(
            collated["mc"][discr][noSF_axis],
            label=["w/o SF"],
            color="tab:gray",
            width=2,
            yerr=True,
            ax=ax,
            # flow=arg.flow
        )
        hep.histplot(
            collated["data"][discr][noSF_axis],
            histtype="errorbar",
            color="black",
            label="Data",
            yerr=True,
            ax=ax,
            xerr=do_xerr
            # flow=arg.flow
        )
        hmc = collated["mc"][discr][SF_axis]
        MCerrorband(hmc, ax=ax, flow=arg.flow)  # stat. unc. errorband
        SFerror = SFerror(collated, discr, flow=arg.flow)
        other = {"hatch": "\\\\", "lw": 0, "color": "r", "alpha": 0.4}
        MCerrorband(
            hmc,
            ax=ax,
            flow=arg.flow,
            ext_error=SFerror,
            label="SF unc.",
            fill_opts=other,
        )  # stat. unc. errorband
        plotratio(
            collated["data"][discr][noSF_axis],
            collated["mc"][discr][noSF_axis],
            ax=rax,
            flow=arg.flow,
            xerr=do_xerr,
        )
        plotratio(
            collated["data"][discr][noSF_axis],
            hmc,
            ax=rax,
            ext_denom_error=SFerror,
            error_opts={"color": "r", "marker": "v"},
            denom_fill_opts=other,
            clear=False,
            label="SF unc.",
            flow=arg.flow,
            xerr=do_xerr,
        )

    elif (
        "syst" in collated["mc"][discr].axes.name and not arg.SF and "noSF" in systlist
    ):
        splitflav_stack = []
        splitflav_axis = noSF_axis
        for i in range(4):
            splitflav_axis["flav"] = i
            splitflav_stack += [collated["mc"][discr][splitflav_axis]]
        noSF_axis["flav"] = sum
        hep.histplot(
            splitflav_stack,
            stack=True,
            histtype="fill",
            label=["udsg", "pileup", "c", "b"],
            yerr=True,
            ax=ax,
            # flow=arg.flow
        )
        hep.histplot(
            collated["data"][discr][noSF_axis],
            histtype="errorbar",
            color="black",
            label="Data",
            yerr=True,
            ax=ax,
            xerr=do_xerr
            # flow=arg.flow
        )
        hmc = collated["mc"][discr][noSF_axis]
        MCerrorband(hmc, ax=ax, flow=arg.flow)  # stat. unc. errorband
        rax = plotratio(
            collated["data"][discr][noSF_axis], hmc, ax=rax, flow=arg.flow, xerr=do_xerr
        )
    elif "flav" in collated["mc"][discr].axes.name:
        splitflav_stack = []
        splitflav_axis = allaxis
        for i in range(4):
            splitflav_axis["flav"] = i
            splitflav_stack += [collated["mc"][discr][splitflav_axis]]

        allaxis["flav"] = sum
        hep.histplot(
            splitflav_stack,
            stack=True,
            histtype="fill",
            label=["udsg", "pileup", "c", "b"],
            yerr=True,
            ax=ax,
            # flow=arg.flow
        )
        hep.histplot(
            collated["data"][discr][allaxis],
            histtype="errorbar",
            color="black",
            label="Data",
            yerr=True,
            ax=ax,
            xerr=do_xerr
            # flow=arg.flow
        )
        hmc = collated["mc"][discr][allaxis]
        MCerrorband(hmc, ax=ax)
        rax = plotratio(
            collated["data"][discr][allaxis], hmc, ax=rax, flow=arg.flow, xerr=do_xerr
        )
    else:
        hmc = collated["mc"][discr][allaxis]
        hep.histplot(
            hmc,
            color="tab:orange",
            histtype="fill",
            label=["MC"],
            yerr=True,
            ax=ax,
            # flow=arg.flow
        )
        hep.histplot(
            collated["data"][discr][allaxis],
            histtype="errorbar",
            color="black",
            label="Data",
            yerr=True,
            ax=ax,
            xerr=do_xerr
            # flow=arg.flow
        )
        MCerrorband(hmc, ax=ax, flow=arg.flow)  # stat. unc. errorband
        rax = plotratio(
            collated["data"][discr][allaxis], hmc, ax=rax, flow=arg.flow, xerr=do_xerr
        )

    ax.set_xlabel(None)
    ax.set_ylabel("Events")
    rax.set_ylabel("Data/MC")
    ax.ticklabel_format(style="sci", scilimits=(-3, 3))
    ax.get_yaxis().get_offset_text().set_position((-0.065, 1.05))
    # FIXME: add wildcard option for xlabel
    xlabel = (
        arg.xlabel if arg.xlabel is not None else collated["data"][discr].axes[-1].label
    )  # Use label from stored hists
    ## FIXME: Set temporary fix for the x-axis
    if arg.xlabel is not None:
        arg.xlabel.split(",")[index]
    elif "DeepJet" in discr or "DeepCSV" in discr or "PFCands" in discr:
        xlabel = (
            bininfo[discr]["displayname"]
            + " ["
            + bininfo[discr]["inputVar_units"]
            + "]"
            if (bininfo[discr]["inputVar_units"] is not None)
            and (bininfo[discr]["inputVar_units"] == "")
            else bininfo[discr]["displayname"]
        )
    else:
        xlabel = axes_name(discr)

    rax.set_xlabel(xlabel)
    ax.legend()
    rax.set_ylim(0, 2.0)
    ax.set_ylim(bottom=0.0)

    rax.autoscale(True, axis="x", tight=True)
    xmin, xmax = autoranger(
        collated["data"][discr][allaxis] + collated["mc"][discr][allaxis], flow=arg.flow
    )
    if arg.xrange is not None:
        xmin, xmax = float(arg.xrange.split(",")[0]), float(arg.xrange.split(",")[1])
    rax.set_xlim(xmin, xmax)
    at = AnchoredText(input_txt + "\n" + arg.ext, loc=2, frameon=False)
    ax.add_artist(at)
    scale = ""
    if arg.norm:
        scale = "_norm"
    name = "all"
    hep.mpl_magic(ax=ax)
    if arg.log:
        print(
            "creating:",
            f"plot/BTV/{arg.phase}_{arg.ext}_{time}/unc_{discr}_inclusive{scale}_{name}.png",
        )
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
        print(
            "creating:",
            f"plot/BTV/{arg.phase}_{arg.ext}_{time}/unc_{discr}_inclusive{scale}_{name}.png",
        )
        fig.savefig(
            f"plot/BTV/{arg.phase}_{arg.ext}_{time}/unc_{discr}_inclusive{scale}_{name}.pdf"
        )
        fig.savefig(
            f"plot/BTV/{arg.phase}_{arg.ext}_{time}/unc_{discr}_inclusive{scale}_{name}.png"
        )
