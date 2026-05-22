import numpy as np
import fnmatch
import argparse, os, glob, json, sys
import matplotlib.pyplot as plt
import mplhep as hep
import hist
from matplotlib.offsetbox import AnchoredText
from coffea.util import load
from BTVNanoCommissioning.workflows import workflows
from BTVNanoCommissioning.helpers.definitions import get_definitions, axes_name
from BTVNanoCommissioning.utils.plot_utils import (
    plotratio,
    markers,
    autoranger,
    rebin_hist,
)

plt.style.use(hep.style.ROOT)
custom_palette = [
    "#3f90da",
    "#ffa90e",
    "#bd1f01",
    "#94a4a2",
    "#832db6",
    "#a96b59",
    "#e76300",
    "#b9ac70",
    "#717581",
    "#92dadd",
]

# Apply the color cycle globally
plt.rcParams["axes.prop_cycle"] = plt.cycler(color=custom_palette)
from BTVNanoCommissioning.helpers.xs_scaler import collate

bininfo = get_definitions()


def get_parser():
    parser = argparse.ArgumentParser(
        description="make comparison for different campaigns"
    )
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
        required=True,
        type=str,
        help="input coffea files (str), splitted different files with ','. Wildcard option * available as well.",
    )
    parser.add_argument("-r", "--ref", required=True, help="referance dataset")
    parser.add_argument(
        "-c",
        "--compared",
        required=True,
        type=str,
        help="compared datasets, splitted by ,",
    )
    parser.add_argument(
        "--sepflav", action="store_true", help="seperate flavour(b/c/light)"
    )
    parser.add_argument("--log", "--logy", action="store_true", help="log on y axis")
    parser.add_argument("--logx", action="store_true", help="log on x axis")
    parser.add_argument(
        "--norm",
        action="store_true",
        help="compare shape, normalized yield to reference",
        default=False,
    )
    parser.add_argument(
        "-v",
        "--variable",
        type=str,
        default="all",
        help="variables to plot, splitted by ,. Wildcard option * available as well. Specifying `all` will run through all variables.",
    )
    parser.add_argument(
        "--xrange", type=str, default=None, help="custom x-range, --xrange xmin,xmax"
    )
    parser.add_argument(
        "--ratio-range",
        type=str,
        default="0.5,1.5",
        help="y-range for ratio panels, default is 0.5 to 1.5. Usage: --ratio-range ymin,ymax",
    )
    parser.add_argument(
        "--flow",
        type=str,
        default="none",
        help="str, optional {none, 'show', 'sum'} Whether plot the under/overflow bin. If 'show', add additional under/overflow bin. If 'sum', add the under/overflow bin content to first/last bin.",
    )
    parser.add_argument("--ext", type=str, default="", help="prefix name/btv name tag")
    parser.add_argument("--com", default="13.6", type=str, help="sqrt(s) in TeV")
    parser.add_argument(
        "--mergemap",
        default=None,
        type=str,
        help="Group list of sample(keys in coffea) as reference/compare set as dictionary format. Keys would be the new lables of the group, it can be also a .json file",
    )
    parser.add_argument(
        "--autorebin",
        default=None,
        help="Rebin the plotting variables, input `int` or `list`. int: merge N bins. list of number: rebin edges(non-uniform bin is possible)",
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
        default=None,
        help="Override the CMS label. By default derived from whether the reference is data or MC.",
    )
    return parser.parse_args()


def main(args):
    output = {}
    if len(args.input.split(",")) > 1:
        output = {i: load(i) for i in args.input.split(",")}
    elif "*" in args.input:
        files = glob.glob(args.input)
        output = {i: load(i) for i in files}
    else:
        output = {args.input: load(args.input)}
    sample_list = []
    for f in output.keys():
        sample_list.extend([m for m in output[f].keys()])

    mergemap = {}
    if not os.path.isdir(f"plot/{args.phase}_{args.ext}/"):
        os.makedirs(f"plot/{args.phase}_{args.ext}/")

    if args.mergemap is None:
        # Single coffea file
        if not any(".coffea" in o for o in output.keys()):
            mergemap[args.ref] = [m for m in output.keys() if args.ref == m]
            for c in args.compared.split(","):
                mergemap[c] = [m for m in output.keys() if c == m]
        # Multiple files
        else:
            if len(list(set(sample_list))) < len(sample_list):
                raise Exception("coffea files contain same datasets! Need mergemap")
            else:
                reflist = []
                for f in output.keys():
                    reflist.extend([m for m in output[f].keys() if args.ref == m])
                mergemap[args.ref] = reflist

                for c in args.compared.split(","):
                    comparelist = []
                    for f in output.keys():
                        comparelist.extend([m for m in output[f].keys() if c == m])
                    mergemap[c] = comparelist
    else:
        if "json" in args.mergemap:
            args.mergemap = open(args.mergemap)
            mergemap = json.load(args.mergemap)
        else:
            mergemap = json.loads(args.mergemap)

    # if dataset duplicated
    collated = {}
    if len(list(set(sample_list))) < len(sample_list):
        collated[args.ref] = collate(
            output[mergemap["fname"][args.ref]], mergemap["dataset"]
        )
        for c in args.compared.split(","):
            collated[c] = collate(output[mergemap["fname"][c]], mergemap["dataset"])
    else:
        collated = collate(output, mergemap)

    ### style settings
    if args.CMS_label is not None:
        label = args.CMS_label
    elif "Run" in args.ref:
        label = "Preliminary"
    else:
        label = "Simulation Preliminary"
    hist_type = "errorbar" if "Run" in args.ref else "step"

    ### input text settings
    phase_label_map = {
        "Wc": "W+c",
        "DY": "DY+jets",
        "QCD": "QCD",
        "semilep": r"t$\bar{t}$ semileptonic",
        "dilep": r"t$\bar{t}$ dileptonic",
        "dijet": "QCD dijet",
    }
    osss_suffix = {1: " OS", -1: " SS"}

    input_txt = next(
        (lbl for key, lbl in phase_label_map.items() if key in args.phase),
        "placeholder",
    )
    if "Wc" in args.phase:
        input_txt += osss_suffix.get(args.splitOSSS, " OS-SS")

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

    # Find variables to plot
    all_vars = collated[args.ref].keys()
    if args.variable == "all":
        var_set = list(all_vars)
    else:
        var_set = [
            var
            for pattern in args.variable.split(",")
            for var in fnmatch.filter(all_vars, pattern)
        ]

    ratio_ymin = float(args.ratio_range.split(",")[0])
    ratio_ymax = float(args.ratio_range.split(",")[1])

    for index, discr in enumerate(var_set):
        if not isinstance(collated[args.ref][discr], hist.hist.Hist):
            continue

        # Skip missing or empty hists
        if (
            discr not in collated[args.ref].keys()
            or (collated[args.ref][discr].values() == 0).all()
        ):
            print(discr, "not in file or empty")
            continue
        if any(
            discr not in collated[c].keys() or (collated[c][discr].values() == 0).all()
            for c in args.compared.split(",")
        ):
            print(discr, "not in file or empty in a compared dataset")
            continue

        allaxis = {}
        if "flav" in collated[args.ref][discr].axes.name:
            allaxis["flav"] = sum
        if "syst" in collated[args.ref][discr].axes.name:
            allaxis["syst"] = "nominal"
        if "osss" in collated[args.ref][discr].axes.name:
            if args.splitOSSS is None:  # OS-SS
                collated[args.ref][discr] = (
                    collated[args.ref][discr][{"osss": 0}]
                    + collated[args.ref][discr][{"osss": 1}] * -1
                )
                for c in args.compared.split(","):
                    collated[c][discr] = (
                        collated[c][discr][{"osss": 0}]
                        + collated[c][discr][{"osss": 1}] * -1
                    )
            elif args.splitOSSS == 1:
                allaxis["osss"] = 0  # opposite sign
            elif args.splitOSSS == -1:
                allaxis["osss"] = 1  # same sign

        do_xerr = False
        if args.autorebin is not None:
            if args.autorebin.isdigit():
                rebin = int(args.autorebin)
            else:
                rebin = np.array([float(i) for i in args.autorebin.split(",")])
                do_xerr = True
            collated[args.ref][discr] = rebin_hist(
                collated[args.ref][discr],
                collated[args.ref][discr].axes[-1].name,
                rebin,
            )
            for c in args.compared.split(","):
                collated[c][discr] = rebin_hist(
                    collated[c][discr], collated[c][discr].axes[-1].name, rebin
                )

        if args.xlabel is not None:
            xlabel = args.xlabel.split(",")[index]
        elif "DeepJet" in discr or "DeepCSV" in discr or "PFCands" in discr:
            xlabel = (
                bininfo[discr]["displayname"]
                + " ["
                + bininfo[discr]["inputVar_units"]
                + "]"
                if bininfo[discr]["inputVar_units"] is not None
                else bininfo[discr]["displayname"]
            )
        else:
            xlabel = axes_name(discr)

        text = args.ext
        ref_norm = np.sum(collated[args.ref][discr][allaxis].values())
        if args.norm:
            text = args.ext + "\nNormalized to Ref."
            for c in args.compared.split(","):
                collated[c][discr] = collated[c][discr] * float(
                    ref_norm / np.sum(collated[c][discr][allaxis].values())
                )

        if args.sepflav:  # split into 3 panels for different flavor
            fig, (ax, rax, rax2, rax3) = plt.subplots(
                4,
                1,
                figsize=(8, 8),
                gridspec_kw={"height_ratios": (3, 1, 1, 1)},
                sharex=True,
            )
            fig.subplots_adjust(hspace=0.06, top=0.92, bottom=0.1, right=0.97)
            ax.set_xlabel(None)
            hep.cms.label(label, com=args.com, data=True, loc=0, ax=ax)
            laxis = {"flav": 0}
            puaxis = {"flav": 1}
            caxis = {"flav": 2}
            baxis = {"flav": 3}

            if "syst" in collated[args.ref][discr].axes.name:
                laxis["syst"] = "nominal"
                puaxis["syst"] = "nominal"
                caxis["syst"] = "nominal"
                baxis["syst"] = "nominal"

            hep.histplot(
                collated[args.ref][discr][laxis] + collated[args.ref][discr][puaxis],
                label=args.ref + "-l",
                color="b",
                histtype=hist_type,
                yerr=True,
                ax=ax,
                xerr=do_xerr,
                flow=args.flow,
            )
            hep.histplot(
                collated[args.ref][discr][caxis],
                label=args.ref + "-c",
                color="g",
                histtype=hist_type,
                yerr=True,
                ax=ax,
                xerr=do_xerr,
                flow=args.flow,
            )
            hep.histplot(
                collated[args.ref][discr][baxis],
                label=args.ref + "-b",
                yerr=True,
                color="r",
                histtype=hist_type,
                ax=ax,
                xerr=do_xerr,
                flow=args.flow,
            )

            mindex = 0
            ax.legend(ncol=3, loc=1)
            for c in args.compared.split(","):
                hep.histplot(
                    collated[c][discr][laxis] + collated[c][discr][puaxis],
                    label=c + "-l",
                    color="b",
                    marker=markers[mindex + 1],
                    histtype="errorbar",
                    yerr=True,
                    ax=ax,
                    xerr=do_xerr,
                    flow=args.flow,
                )
                hep.histplot(
                    collated[c][discr][caxis],
                    label=c + "-c",
                    color="g",
                    marker=markers[mindex + 1],
                    histtype="errorbar",
                    yerr=True,
                    ax=ax,
                    xerr=do_xerr,
                    flow=args.flow,
                )
                hep.histplot(
                    collated[c][discr][baxis],
                    label=c + "-b",
                    yerr=True,
                    color="r",
                    marker=markers[mindex + 1],
                    histtype="errorbar",
                    ax=ax,
                    xerr=do_xerr,
                    flow=args.flow,
                )
                # comparison splitted by flavor
                rax = plotratio(
                    collated[c][discr][laxis] + collated[c][discr][puaxis],
                    collated[args.ref][discr][laxis]
                    + collated[args.ref][discr][puaxis],
                    ax=rax,
                    denom_fill_opts=None,
                    error_opts={"color": "b", "marker": markers[mindex + 1]},
                    clear=False,
                    xerr=do_xerr,
                    flow=args.flow,
                )
                rax2 = plotratio(
                    collated[c][discr][caxis],
                    collated[args.ref][discr][caxis],
                    ax=rax2,
                    denom_fill_opts=None,
                    error_opts={"color": "g", "marker": markers[mindex + 1]},
                    clear=False,
                    xerr=do_xerr,
                    flow=args.flow,
                )
                rax3 = plotratio(
                    collated[c][discr][baxis],
                    collated[args.ref][discr][baxis],
                    ax=rax3,
                    denom_fill_opts=None,
                    error_opts={"color": "r", "marker": markers[mindex + 1]},
                    clear=False,
                    xerr=do_xerr,
                    flow=args.flow,
                )
                mindex += 1

            ax.set_xlabel("A.U.")
            rax.set_ylabel("udsg-jets")
            rax2.set_ylabel("c-jets")
            rax3.set_ylabel("b-jets")
            rax.set_ylim(ratio_ymin, ratio_ymax)
            rax2.set_ylim(ratio_ymin, ratio_ymax)
            rax3.set_ylim(ratio_ymin, ratio_ymax)
            rax3.set_xlabel(xlabel)
            ax.legend()

            at = AnchoredText(input_txt + "\n" + text, loc=2, frameon=False)
            ax.add_artist(at)

            if args.log:
                ax.set_yscale("log")
            if args.logx:
                ax.set_xscale("log")
            hep.mpl_magic(ax=ax, soft_fail=True)

            for filetype in args.save_as.split(","):
                fig.savefig(
                    f"plot/{args.phase}_{args.ext}/compare_{args.phase}_inclusive{discr}.{filetype}"
                )

        else:
            fig, ((ax), (rax)) = plt.subplots(
                2,
                1,
                figsize=(12, 12),
                gridspec_kw={"height_ratios": (3, 1)},
                sharex=True,
            )
            fig.subplots_adjust(hspace=0.06, top=0.92, bottom=0.1, right=0.97)
            hep.cms.label(label, com=args.com, data=True, loc=0, ax=ax)
            ax.set_xlabel(None)

            hep.histplot(
                collated[args.ref][discr][allaxis],
                label=args.ref + " (Ref)",
                histtype=hist_type,
                yerr=True,
                ax=ax,
                xerr=do_xerr,
                flow=args.flow,
            )

            for i, c in enumerate(args.compared.split(",")):
                hep.histplot(
                    collated[c][discr][allaxis],
                    label=c,
                    histtype=hist_type,
                    yerr=True,
                    ax=ax,
                    xerr=do_xerr,
                    flow=args.flow,
                )
                plotratio(
                    collated[c][discr][allaxis],
                    collated[args.ref][discr][allaxis],
                    ax=rax,
                    denom_fill_opts=None,
                    error_opts={"color": custom_palette[i + 1]},
                    clear=False,
                    xerr=do_xerr,
                    flow=args.flow,
                )  ## No error band used

            xmin, xmax = autoranger(collated[args.ref][discr][allaxis])
            if args.xrange is not None:
                xmin, xmax = float(args.xrange.split(",")[0]), float(
                    args.xrange.split(",")[1]
                )
            rax.set_xlim(xmin, xmax)
            rax.set_xlabel(xlabel)
            ax.set_xlabel(None)
            ax.set_ylabel(args.ylabel)
            rax.set_ylabel("Other/Ref")
            ax.legend(loc=1)
            rax.axhline(1.0, ls=":")
            rax.set_ylim(ratio_ymin, ratio_ymax)
            ax.set_ylim(bottom=0)

            at = AnchoredText(input_txt + "\n" + text, loc=2, frameon=False)
            ax.add_artist(at)

            logext = ""
            normtext = "_norm" if args.norm else ""
            if args.log:
                ax.set_yscale("log")
                ax.set_ylim(bottom=0.1)
                logext = "_log"
            if args.logx:
                ax.set_xscale("log")
            hep.mpl_magic(ax=ax, soft_fail=True)

            print(
                "creating:",
                f"plot/{args.phase}_{args.ext}/compare_{args.phase}_inclusive{discr}{logext}{normtext}",
            )
            for filetype in args.save_as.split(","):
                fig.savefig(
                    f"plot/{args.phase}_{args.ext}/compare_{args.phase}_inclusive{discr}{logext}{normtext}.{filetype}"
                )

        plt.close(fig)


if __name__ == "__main__":
    args = get_parser()
    sys.exit(main(args))
