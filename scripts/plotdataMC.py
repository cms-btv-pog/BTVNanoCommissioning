import numpy as np
import fnmatch
import argparse, os, glob, sys
import matplotlib.pyplot as plt
import mplhep as hep
import hist
from coffea.util import load
from matplotlib.offsetbox import AnchoredText

plt.style.use(hep.style.ROOT)

from BTVNanoCommissioning.workflows import workflows
from BTVNanoCommissioning.helpers.xs_scaler import collate, scaleSumW
from BTVNanoCommissioning.helpers.definitions import (
    get_definitions,
    axes_name,
)
from BTVNanoCommissioning.utils.plot_utils import (
    plotratio,
    SFerror,
    MCerrorband,
    autoranger,
    rebin_hist,
    sample_mergemap,
    color_map,
)

bininfo = get_definitions()
SV_bininfo = get_definitions(include_definitions=["SV"])


def get_parser():
    parser = argparse.ArgumentParser(description="Hist plotter for commissioning")
    parser.add_argument("--lumi", required=True, type=float, help="Luminosity in /pb")
    parser.add_argument("--com", default="13.6", type=str, help="Sqrt(s) in TeV")
    parser.add_argument(
        "-p",
        "--phase",
        required=True,
        choices=list(workflows.keys()),
        dest="phase",
        help="Which phase space",
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        default="",
        help="Input coffea files (str), split different files with ','. Wildcard option * available as well.",
    )
    parser.add_argument("--log", "--logy", action="store_true", help="Log on y axis")
    parser.add_argument("--logx", action="store_true", help="Log on x axis")
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
        default="all",
        help="Variables to plot, split by ,. Wildcard option * available as well. Specifying `all` will run through all variables.",
    )
    parser.add_argument(
        "--xlabel",
        type=str,
        help="Rename the label of variables to plot, split by ,.  Wildcard option * NOT available here",
    )
    parser.add_argument(
        "--ylabel", type=str, default=None, help="Modify y-axis label of plot"
    )
    parser.add_argument(
        "--SF", action="store_true", default=False, help="Make with SF comparisons"
    )
    parser.add_argument(
        "--xrange",
        type=str,
        default=None,
        help="Custom x-range. Usage: --xrange xmin,xmax",
    )
    parser.add_argument(
        "--ratio-range",
        type=str,
        default="0.5,1.5",
        help="y-range for ratio plot, default is 0.5 to 1.5. Usage: --ratio-range ymin,ymax",
    )
    parser.add_argument("--ext", type=str, default="", help="Prefix name")
    parser.add_argument(
        "--autorebin",
        default=None,
        help="Rebin the plotting variables, input `int` or `list`. int: merge N bins. list of number: rebin edges(non-uniform bin is possible)",
    )
    parser.add_argument(
        "--flow",
        type=str,
        default="none",
        help="str, optional {'none', 'show', 'sum'} Whether plot the under/overflow bin. If 'show', add additional under/overflow bin. If 'sum', add the under/overflow bin content to first/last bin.",
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
    parser.add_argument(
        "--save-as",
        type=str,
        default="png,pdf",
        help="File format to save the plots, split by ,. Default is png and pdf",
    )
    parser.add_argument(
        "--CMS-label",
        type=str,
        default="Preliminary",
        help="CMS label to put on the plot, default is Preliminary",
    )

    return parser.parse_args()


def main(args):
    inputs = args.input.split(",")

    # Load the inputs
    output = {}
    for inp in inputs:
        if not os.path.exists(inp) and "*" not in inp:
            print(f"{inp} does not exist!")
            return 1

        files = glob.glob(inp)

        for i in files:
            if ".coffea" not in i:
                print(f"{i} is not a valid coffea file")
                return 1
            fdict = load(i)
            for sample in fdict:
                if sample not in output:
                    output[sample] = fdict[sample]
                else:
                    for key in fdict[sample]:
                        if key not in output[sample]:
                            continue
                        output[sample][key] += fdict[sample][key]

    output = scaleSumW(output, args.lumi)

    mergemap = {}
    ## create merge map from sample set/data MC
    if not os.path.isdir(f"plot/{args.phase}_{args.ext}/"):
        os.makedirs(f"plot/{args.phase}_{args.ext}/")
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
    collated = {
        key: value for key, value in collated.items() if isinstance(collated[key], dict)
    }
    print(collated.keys())

    # Which flavor to split, default is udsg+pu+c+b, for QG workflow, split to ud+s+c+b+g+other
    flav_set = ["udsg", "pu", "c", "b"]
    if "qg" in args.phase.lower():
        flav_set = ["ud", "s", "c", "b", "g", "other"]

    ### input text settings
    phase_label_map = {
        "Wc": "W+c",
        "DY": "DY+jets",
        "QCD": "QCD",
        "semilep": r"t$\bar{t}$ semileptonic",
        "dilep": r"t$\bar{t}$ dileptonic",
        "dijet": "QCD dijet",
    }
    phase_nj_map = {"semilep": 4, "dilep": 2}
    osss_suffix = {1: " OS", -1: " SS"}

    input_txt = next(
        (label for key, label in phase_label_map.items() if key in args.phase),
        "placeholder",
    )
    if "Wc" in args.phase:
        input_txt += osss_suffix.get(args.splitOSSS, " OS-SS")
    for key, nj_val in phase_nj_map.items():
        if key in args.phase:
            nj = nj_val
            break

    if (
        "njet" in args.variable.split(",")
        or "nbjet" in args.variable.split(",")
        or "mu" in args.variable.split(",")
    ):
        nj = 1

    if "emctag" in args.phase:
        input_txt += r" (e$\mu$)"
    elif "ectag" in args.phase:
        input_txt += " (e)"
    elif args.phase == "ttdilep_sf":
        input_txt += r" (e$\mu$)"
    elif args.phase == "QCD" or "qg" in args.phase.lower():
        pass
    else:
        input_txt += r" ($\mu$)"

    if "ctag" in args.phase and "DY" not in args.phase:
        input_txt += "\nw/ soft-$\mu$"

    # Find the variables to plot
    if args.variable == "all":
        var_set = [var for var in collated["mc"].keys()]
    else:
        all_vars = collated["mc"].keys()
        var_set = [
            var
            for pattern in args.variable.split(",")
            for var in fnmatch.filter(all_vars, pattern)
        ]

    for index, discr in enumerate(var_set):
        # Skip non-histogram objects
        try:
            if not isinstance(collated["mc"][discr], hist.hist.Hist):
                continue
        except:
            print(f"{discr} not found. Variable must be in", collated["mc"].keys())
            return 1

        # Skip missing and empty hists
        if discr not in collated["mc"].keys() or discr not in collated["data"].keys():
            print(discr, "not in files, skipping")
            continue
        elif (collated["mc"][discr].values() == 0).all() or (
            collated["data"][discr].values() == 0
        ).all():
            print(discr, "is empty, skipping")
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
            # Get list of available systematics
            systlist = [
                collated["mc"][discr]
                .axes[collated["mc"][discr].axes.name.index("syst")]
                .value(i)
                for i in range(
                    collated["mc"][discr]
                    .axes[collated["mc"][discr].axes.name.index("syst")]
                    .size
                )
            ]
            print(f"Available systematics: {systlist}")

            # Choose the appropriate systematic name
            if "nominal" in systlist:
                allaxis["syst"] = "nominal"
            elif "SF" in systlist:  # For ctag_Wc_sf workflow
                allaxis["syst"] = "SF"
            elif len(systlist) > 0:
                # Fallback to first available systematic
                allaxis["syst"] = systlist[0]
                print(
                    f"Warning: 'nominal' not found in systematics, using '{systlist[0]}' instead"
                )

            SF_axis = allaxis.copy()
            noSF_axis = allaxis.copy()
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

        # check if plottable, skip otherwise
        hmc = collated["mc"][discr][allaxis]
        if len(hmc.axes) > 2:
            print(f"Skipping {discr}: after selections still has {len(hmc.axes)} axes")
            continue

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
            args.CMS_label,
            data=True,
            lumi=args.lumi / 1000.0,
            com=args.com,
            loc=0,
            ax=ax,
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
                for i, t in enumerate(flav_set):
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
                sort="y",  # sort by yield
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
                sort="y",  # sort by yield
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
                sort="y",  # sort by yield
            )
            hmc = collated["mc"][discr][SF_axis]
            MCerrorband(hmc, ax=ax, flow=args.flow)  # stat. unc. errorband
            sf_err = SFerror(collated, discr, flow=args.flow)
            other = {"hatch": "\\\\", "lw": 0, "color": "r", "alpha": 0.4}
            MCerrorband(
                hmc,
                ax=ax,
                flow=args.flow,
                ext_error=sf_err,
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
                ext_denom_error=sf_err,
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
                for i, t in enumerate(flav_set):
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
                sort="y",  # sort by yield
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
                sort="y",  # sort by yield
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
                for i, t in enumerate(flav_set):
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
                sort="y",  # sort by yield
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
                sort="y",  # sort by yield
            )
            hmc = collated["mc"][discr][allaxis]
            MCerrorband(hmc, ax=ax)
            rax = plotratio(
                collated["data"][discr][allaxis],
                hmc,
                ax=rax,
                flow=args.flow,
                xerr=do_xerr,
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
                        color_map[s]
                        for s in collated.keys()
                        if s != "mc" and s != "data"
                    ],
                    flow=args.flow,
                    sort="y",  # sort by yield
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
                    sort="y",  # sort by yield
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
                sort="y",  # sort by yield
            )
            MCerrorband(hmc, ax=ax, flow=args.flow)  # stat. unc. errorband
            rax = plotratio(
                collated["data"][discr][allaxis],
                hmc,
                ax=rax,
                flow=args.flow,
                xerr=do_xerr,
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
        elif "JetSVs_" in discr:
            xlabel = (
                SV_bininfo[discr]["displayname"]
                + " ["
                + SV_bininfo[discr]["inputVar_units"]
                + "]"
                if (SV_bininfo[discr]["inputVar_units"] is not None)
                and (SV_bininfo[discr]["inputVar_units"] == "")
                else SV_bininfo[discr]["displayname"]
            )
        else:
            xlabel = axes_name(discr)

        rax.set_xlabel(xlabel)

        if "sample" in args.split:
            ax.legend(ncols=2, prop={"size": 16})
        else:
            ax.legend(
                ncols=(2 if len(flav_set) > 3 else 1),
            )

        rax.set_ylim(
            float(args.ratio_range.split(",")[0]), float(args.ratio_range.split(",")[1])
        )
        ax.set_ylim(bottom=0.0)

        rax.autoscale(True, axis="x", tight=True)
        xmin, xmax = autoranger(
            collated["data"][discr][allaxis] + collated["mc"][discr][allaxis],
            flow=args.flow,
        )

        if args.xrange is not None:
            xmin, xmax = float(args.xrange.split(",")[0]), float(
                args.xrange.split(",")[1]
            )

        rax.set_xlim(xmin, xmax)
        at = AnchoredText(input_txt + "\n" + args.ext, loc=2, frameon=False)
        ax.add_artist(at)

        scale = ""
        if args.norm:
            scale = "_norm"
        name = "all"
        if args.split == "sample":
            name += "_sample"

        try:
            hep.mpl_magic(ax=ax)
        except RuntimeError as e:
            print(f"Warning: {e}")
            print("Using soft_fail=True for legend placement")
            try:
                # Try with soft_fail=True
                hep.mpl_magic(ax=ax, soft_fail=True)
            except Exception as e2:
                print(f"Still failed: {e2}")
                # Continue anyway - the plot will still be usable

        if args.log:
            name += "_log"
            ax.set_yscale("log")
            ax.set_ylim(bottom=0.1)
            hep.mpl_magic(ax=ax)

        if args.logx:
            name += "_logx"
            ax.set_xscale("log")

        # Save the plots
        for filetype in args.save_as.split(","):
            if filetype not in ["png", "pdf"]:
                print(
                    f"Unsupported file type: {filetype}. Only 'png' and 'pdf' are supported."
                )
                continue
            print(
                "creating:",
                f"plot/{args.phase}_{args.ext}/unc_{discr}_inclusive{scale}_{name}.{filetype}",
            )
            fig.savefig(
                f"plot/{args.phase}_{args.ext}/unc_{discr}_inclusive{scale}_{name}.{filetype}"
            )

        plt.close(fig)


if __name__ == "__main__":
    args = get_parser()
    sys.exit(main(args))
