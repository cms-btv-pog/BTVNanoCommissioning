import numpy as np
import argparse, os, arrow, glob, re
from coffea.util import load
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import mplhep as hep
import hist

plt.style.use(hep.style.ROOT)
from BTVNanoCommissioning.workflows import workflows
from BTVNanoCommissioning.helpers.xs_scaler import collate, scaleSumW
from BTVNanoCommissioning.helpers.definitions import definitions, axes_name
from BTVNanoCommissioning.utils.plot_utils import (
    plotratio,
    SFerror,
    MCerrorband,
    autoranger,
    rebin_hist,
    sample_mergemap,
    color_map,
)

bininfo = definitions()
parser = argparse.ArgumentParser(description="hist plotter for commissioning")
parser.add_argument("--lumi", required=True, type=float, help="luminosity in /pb")
parser.add_argument("--com", default="13.6", type=str, help="sqrt(s) in TeV")
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
    default="show",
    help="str, optional {None, 'show', 'sum'} Whether plot the under/overflow bin. If 'show', add additional under/overflow bin. If 'sum', add the under/overflow bin content to first/last bin.",
)
parser.add_argument(
    "--split",
    type=str,
    default="flavor",
    choices=["flavor", "sample", "sample_flav"],
    help="Decomposition of MC samples. Default is split to jet flavor, possible to split by group of MC samples. Combination of jetflavor+ sample split is also possible",
)
parser.add_argument(
    "--splitOSSS",
    type=int,
    default=None,
    help="Only for W+c phase space, split opposite sign(1) and same sign events(-1), if not specified, the combined OS-SS phase space is used",
)


args = parser.parse_args()
time = arrow.now().format("YY_MM_DD")
if not os.path.isdir(f"plot/BTV/{args.phase}_{args.ext}_{time}/"):
    os.makedirs(f"plot/BTV/{args.phase}_{args.ext}_{time}/")
if len(args.input.split(",")) > 1:
    output = {i: load(i) for i in args.input.split(",")}
elif "*" in args.input:
    files = glob.glob(args.input)
    output = {i: load(i) for i in files}
else:
    output = {args.input: load(args.input)}
output = scaleSumW(output, args.lumi)
mergemap = {}
## create merge map from sample set/data MC
if not any(".coffea" in o for o in output.keys()):
    if "sample" in args.split:
        mergemap = sample_mergemap
    mergemap["data"] = [m for m in output.keys() if "Run" in m]
    mergemap["mc"] = [m for m in output.keys() if "Run" not in m]
else:
    datalist = []
    mclist = []
    for f in output.keys():
        datalist.extend([m for m in output[f].keys() if "Run" in m])
        mclist.extend([m for m in output[f].keys() if "Run" not in m])
    if "sample" in args.split:
        mergemap = sample_mergemap
    mergemap["mc"] = mclist
    mergemap["data"] = datalist

collated = collate(output, mergemap)

for sample in mergemap.keys():
    if type(collated[sample]) is not dict:
        del collated[sample]
### input text settings
if "Wc" in args.phase:
    input_txt = "W+c"
    if args.splitOSSS == 1:
        input_txt = input_txt + " OS"
    elif args.splitOSSS == -1:
        input_txt = input_txt + " SS"
    else:
        input_txt = input_txt + " OS-SS"
elif "DY" in args.phase:
    input_txt = "DY+jets"
elif "semilep" in args.phase:
    input_txt = r"t$\bar{t}$ semileptonic"
    nj = 4
elif "dilep" in args.phase:
    input_txt = r"t$\bar{t}$ dileptonic"
    nj = 2
if (
    "njet" in args.variable.split(",")
    or "nbjet" in args.variable.split(",")
    or "mu" in args.variable.split(",")
):
    nj = 1
if "emctag" in args.phase:
    input_txt = input_txt + " (e$\mu$)"
elif "ectag" in args.phase:
    input_txt = input_txt + " (e)"
elif "ttdilep_sf" == args.phase:
    input_txt = input_txt + " (e$\mu$)"
else:
    input_txt = input_txt + " ($\mu$)"
if "ctag" in args.phase and "DY" not in args.phase:
    input_txt = input_txt + "\nw/ soft-$\mu$"
if args.variable == "all":
    var_set = [
        var
        for var in collated["mc"].keys()
        if var not in ["fname", "run", "lumi", "sumw"]
    ]
elif "*" in args.variable:
    if args.variable.count("*") > 1:
        var_set = [
            var
            for var in collated["mc"].keys()
            if args.variable.replace("*", "") in var
        ]
    elif args.variable.startswith("*") or args.variable.endswith("*"):
        var_set = [
            var
            for var in collated["mc"].keys()
            if var.startswith(args.variable.replace("*", ""))
            or var.endswith(args.variable.replace("*", ""))
        ]
    else:
        var_set = [
            var
            for var in collated["mc"].keys()
            if re.match(
                f"^{args.variable.split('*')[0]}.*{args.variable.split('*')[1]}$", var
            )
            != None
        ]

else:
    var_set = args.variable.split(",")
for index, discr in enumerate(var_set):
    ## remove empty
    if (
        discr not in collated["mc"].keys()
        or discr not in collated["data"].keys()
        or (collated["mc"][discr].values() == 0).all()
        or (collated["data"][discr].values() == 0).all()
    ):
        print(discr, "not in file or empty")
        continue

    ## axis info
    allaxis = {}
    if "Wc" in args.phase:
        if args.splitOSSS is None:  # OS-SS
            for sample in collated.keys():
                if discr not in collated[sample].keys():
                    continue
                collated[sample][discr] = (
                    collated[sample][discr][{"osss": 0}]
                    + collated[sample][discr][{"osss": 1}] * -1
                )
        elif args.splitOSSS == 1:
            allaxis["osss"] = 0  # opposite sign
        elif args.splitOSSS == -1:
            allaxis["osss"] = 1  # same sign
    if "flav" in collated["mc"][discr].axes.name:
        allaxis["flav"] = sum
        SF_axis = allaxis
        noSF_axis = allaxis
    if "syst" in collated["mc"][discr].axes.name:
        allaxis["syst"] = "nominal"
        SF_axis = allaxis
        noSF_axis = allaxis
        systlist = [
            collated["mc"][discr].axes[0].value(i)
            for i in range(collated["mc"][discr].axes[0].size)
        ]
        if "noSF" in systlist:
            noSF_axis["syst"] = "noSF"

    ## rebin config, add xerr
    do_xerr = False
    if args.autorebin is not None:
        if args.autorebin.isdigit():
            rebin = int(args.autorebin)
        else:
            rebin = np.array([float(i) for i in args.autorebin.split(",")])
            do_xerr = True
        for s in collated.keys():
            collated[s][discr] = rebin_hist(
                collated[s][discr], collated[s][discr].axes[-1].name, rebin
            )

    ## Rescale noSF & SF to same MC yields
    if (
        "flav" in collated["mc"][discr].axes.name
        and "syst" in collated["mc"][discr].axes.name
        and args.SF
    ):
        scale_sf = np.sum(collated["mc"][discr][SF_axis].values()) / np.sum(
            collated["mc"][discr][noSF_axis].values()
        )
    else:
        scale_sf = 1.0
    if (
        "flav" in collated["mc"][discr].axes.name
        and "syst" in collated["mc"][discr].axes.name
        and args.SF
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
        "Preliminary", data=True, lumi=args.lumi / 1000.0, com=args.com, loc=0, ax=ax
    )

    ## w/ & w/o btag SF
    if (
        "flav" in collated["mc"][discr].axes.name
        and "syst" in collated["mc"][discr].axes.name
        and args.SF
        and args.split != "sample"
    ):
        splitflav_stack, labels, color_config = [], [], {}
        color_config["color"], color_config["alpha"] = [], []
        splitflav_axis = SF_axis
        for sample in collated.keys():
            if sample == "data":
                continue
            if args.split == "flavor" and sample != "mc":
                continue
            elif args.split == "sample_flav" and sample == "mc":
                continue
            for i, t in enumerate(["udsg", "pu", "c", "b"]):
                splitflav_axis["flav"] = i
                splitflav_stack.append(collated[sample][discr][splitflav_axis])
                if args.split == "flavor":
                    labels.append(t)
                    color_config["color"].append(color_map[t])
                    color_config["alpha"].append(1.0)
                else:
                    labels.append(f"{sample} ({t})")
                    color_config["color"].append(color_map[sample])
                    color_config["alpha"].append(1.0 - 0.2 * i)

        SF_axis["flav"] = sum
        hep.histplot(
            splitflav_stack,
            stack=True,
            label=labels,
            histtype="fill",
            yerr=True,
            ax=ax,
            flow=args.flow,
            **color_config,
        )
        hep.histplot(
            collated["mc"][discr][noSF_axis],
            label=["w/o SF"],
            color="tab:gray",
            width=2,
            yerr=True,
            ax=ax,
            flow=args.flow,
        )
        hep.histplot(
            collated["data"][discr][noSF_axis],
            histtype="errorbar",
            color="black",
            label="Data",
            yerr=True,
            ax=ax,
            xerr=do_xerr,
            flow=args.flow,
        )
        hmc = collated["mc"][discr][SF_axis]
        MCerrorband(hmc, ax=ax, flow=args.flow)  # stat. unc. errorband
        SFerror = SFerror(collated, discr, flow=args.flow)
        other = {"hatch": "\\\\", "lw": 0, "color": "r", "alpha": 0.4}
        MCerrorband(
            hmc,
            ax=ax,
            flow=args.flow,
            ext_error=SFerror,
            label="SF unc.",
            fill_opts=other,
        )  # stat. unc. errorband
        plotratio(
            collated["data"][discr][noSF_axis],
            collated["mc"][discr][noSF_axis],
            ax=rax,
            flow=args.flow,
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
            flow=args.flow,
            xerr=do_xerr,
        )

    elif (
        "syst" in collated["mc"][discr].axes.name
        and not args.SF
        and "noSF" in systlist
        and args.split != "sample"
    ):
        splitflav_stack, labels, color_config = [], [], {}
        color_config["color"], color_config["alpha"] = [], []
        splitflav_axis = noSF_axis
        for sample in collated.keys():
            if sample == "data":
                continue
            if args.split == "flavor" and sample != "mc":
                continue
            elif args.split == "sample_flav" and sample == "mc":
                continue
            for i, t in enumerate(["udsg", "pu", "c", "b"]):
                splitflav_axis["flav"] = i
                splitflav_stack.append(collated[sample][discr][splitflav_axis])
                if args.split == "flavor":
                    labels.append(t)
                    color_config["color"].append(color_map[t])
                    color_config["alpha"].append(1.0)
                else:
                    labels.append(f"{sample} ({t})")
                    color_config["color"].append(color_map[sample])
                    color_config["alpha"].append(1.0 - 0.2 * i)
        noSF_axis["flav"] = sum
        hep.histplot(
            splitflav_stack,
            stack=True,
            histtype="fill",
            label=labels,
            yerr=True,
            ax=ax,
            flow=args.flow,
            **color_config,
        )
        hep.histplot(
            collated["data"][discr][noSF_axis],
            histtype="errorbar",
            color="black",
            label="Data",
            yerr=True,
            ax=ax,
            xerr=do_xerr,
            flow=args.flow,
        )
        hmc = collated["mc"][discr][noSF_axis]
        MCerrorband(hmc, ax=ax, flow=args.flow)  # stat. unc. errorband
        rax = plotratio(
            collated["data"][discr][noSF_axis],
            hmc,
            ax=rax,
            flow=args.flow,
            xerr=do_xerr,
        )
    elif "flav" in collated["mc"][discr].axes.name and args.split != "sample":
        splitflav_stack, labels, color_config = [], [], {}
        color_config["color"], color_config["alpha"] = [], []
        splitflav_axis = allaxis
        for sample in collated.keys():
            if sample == "data":
                continue
            if args.split == "flavor" and sample != "mc":
                continue
            elif args.split == "sample_flav" and sample == "mc":
                continue
            for i, t in enumerate(["udsg", "pu", "c", "b"]):
                splitflav_axis["flav"] = i
                splitflav_stack.append(collated[sample][discr][splitflav_axis])

                if args.split == "flavor":
                    labels.append(t)
                    color_config["color"].append(color_map[t])
                    color_config["alpha"].append(1.0)
                else:
                    labels.append(f"{sample} ({t})")
                    color_config["color"].append(color_map[sample])

                    color_config["alpha"].append(1.0 - 0.2 * i)
        allaxis["flav"] = sum
        hep.histplot(
            splitflav_stack,
            stack=True,
            histtype="fill",
            label=labels,
            yerr=True,
            ax=ax,
            flow=args.flow,
            **color_config,
        )
        hep.histplot(
            collated["data"][discr][allaxis],
            histtype="errorbar",
            color="black",
            label="Data",
            yerr=True,
            ax=ax,
            xerr=do_xerr,
            flow=args.flow,
        )
        hmc = collated["mc"][discr][allaxis]
        MCerrorband(hmc, ax=ax)
        rax = plotratio(
            collated["data"][discr][allaxis], hmc, ax=rax, flow=args.flow, xerr=do_xerr
        )
    else:
        hmc = collated["mc"][discr][allaxis]
        if "sample" in args.split:
            hep.histplot(
                [
                    collated[s][discr][allaxis]
                    for s in collated.keys()
                    if s != "mc" and s != "data"
                ],
                histtype="fill",
                stack=True,
                label=[s for s in collated.keys() if s != "mc" and s != "data"],
                yerr=True,
                ax=ax,
                color=[
                    color_map[s] for s in collated.keys() if s != "mc" and s != "data"
                ],
                flow=args.flow,
            )
        else:
            hep.histplot(
                hmc,
                color="tab:orange",
                histtype="fill",
                label=["MC"],
                yerr=True,
                ax=ax,
                flow=args.flow,
            )
        hep.histplot(
            collated["data"][discr][allaxis],
            histtype="errorbar",
            color="black",
            label="Data",
            yerr=True,
            ax=ax,
            xerr=do_xerr,
            flow=args.flow,
        )
        MCerrorband(hmc, ax=ax, flow=args.flow)  # stat. unc. errorband
        rax = plotratio(
            collated["data"][discr][allaxis], hmc, ax=rax, flow=args.flow, xerr=do_xerr
        )

    ax.set_xlabel(None)
    ax.set_ylabel("Events")
    rax.set_ylabel("Data/MC")
    ax.ticklabel_format(style="sci", scilimits=(-3, 3))
    ax.get_yaxis().get_offset_text().set_position((-0.065, 1.05))
    # FIXME: add wildcard option for xlabel
    xlabel = (
        args.xlabel
        if args.xlabel is not None
        else collated["data"][discr].axes[-1].label
    )  # Use label from stored hists
    ## FIXME: Set temporary fix for the x-axis
    if args.xlabel is not None:
        args.xlabel.split(",")[index]
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
    if "sample" in args.split:
        ax.legend(ncols=2, prop={"size": 16})
    else:
        ax.legend()
    rax.set_ylim(0.5, 1.5)
    ax.set_ylim(bottom=0.0)

    rax.autoscale(True, axis="x", tight=True)
    xmin, xmax = autoranger(
        collated["data"][discr][allaxis] + collated["mc"][discr][allaxis],
        flow=args.flow,
    )
    if args.xrange is not None:
        xmin, xmax = float(args.xrange.split(",")[0]), float(args.xrange.split(",")[1])
    rax.set_xlim(xmin, xmax)
    at = AnchoredText(input_txt + "\n" + args.ext, loc=2, frameon=False)
    ax.add_artist(at)
    scale = ""
    if args.norm:
        scale = "_norm"
    name = "all"
    hep.mpl_magic(ax=ax)
    if args.log:
        print(
            "creating:",
            f"plot/BTV/{args.phase}_{args.ext}_{time}/unc_{discr}_inclusive{scale}_{name}.png",
        )
        ax.set_yscale("log")
        name = "log"
        ax.set_ylim(bottom=0.1)
        hep.mpl_magic(ax=ax)
        fig.savefig(
            f"plot/BTV/{args.phase}_{args.ext}_{time}/unc_{discr}_inclusive{scale}_{name}.pdf"
        )
        fig.savefig(
            f"plot/BTV/{args.phase}_{args.ext}_{time}/unc_{discr}_inclusive{scale}_{name}.png"
        )
    else:
        print(
            "creating:",
            f"plot/BTV/{args.phase}_{args.ext}_{time}/unc_{discr}_inclusive{scale}_{name}.png",
        )
        fig.savefig(
            f"plot/BTV/{args.phase}_{args.ext}_{time}/unc_{discr}_inclusive{scale}_{name}.pdf"
        )
        fig.savefig(
            f"plot/BTV/{args.phase}_{args.ext}_{time}/unc_{discr}_inclusive{scale}_{name}.png"
        )
