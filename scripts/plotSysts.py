from coffea.util import load
from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import glob, json, argparse, math, copy, os, warnings

from BTVNanoCommissioning.utils.plot_utils import (
    sample_mergemap,
    MCerrorband,
    errband_opts,
    plotratio,
    color_map,
)
from BTVNanoCommissioning.helpers.xs_scaler import collate, scaleSumW
from BTVNanoCommissioning.helpers.definitions import axes_name

warnings.filterwarnings("ignore")
hep.style.use("CMS")
hep.style.use("CMS")


parser = argparse.ArgumentParser(description="Make plots with systematics from coffea files.")
parser.add_argument(
    "-i",
    "--input",
    type=str,
    required=True,
    help="Input coffea file(s).",
)
parser.add_argument(
    "-o",
    "--outdir",
    type=str,
    default="plot",
    help="Output directory.",
)
parser.add_argument(
    "-v",
    "--vars",
    type=str,
    required=True,
    help="Variable(s) to plot (histogram name), separated by commas.",
)
parser.add_argument(
    "-w",
    "--workflow",
    type=str,
    required=True,
    choices=[
        "2D_e_DY_sf",
        "2D_mu_DY_sf",
        "2D_e_Wc_sf",
        "2D_mu_Wc_sf",
        "2D_emu_ttdilep_sf",
        "2D_e_ttsemilep_sf",
        "2D_mu_ttsemilep_sf",
    ],
    help="Workflow.",
)
parser.add_argument(
    "--campaign",
    type=str,
    choices=[
        "Run3",
        "2425",
        "Summer22",
        "Summer22EE",
        "Summer23",
        "Summer23BPix",
        "Summer24",
        "Prompt25",
    ],
    required=True,
    help="Campaign.",
)
parser.add_argument(
    "-s",
    "--split",
    type=str,
    default="sample",
    choices=["sample", "flavour"],
    help="Histogram stack splitting.",
)
parser.add_argument(
    "--splitOSSS",
    type=int,
    default=None,
    choices=[None, 1, -1],
    help="Only for W+c phase space, split opposite sign (1) and same sign events (-1). If not specified, the combined OS-SS phase space is used.",
)
parser.add_argument(
    "-e",
    "--exp_systs",
    action="append",
    default=[],
    choices=[
        "puweight",
        "JESRegrouped_Absolute_$CAMPAIGN",
        "JESRegrouped_Absolute",
        "JESRegrouped_BBEC1_$CAMPAIGN",
        "JESRegrouped_BBEC1",
        "JESRegrouped_EC2_$CAMPAIGN",
        "JESRegrouped_EC2",
        "JESRegrouped_FlavorQCD",
        "JESRegrouped_HF_$CAMPAIGN",
        "JESRegrouped_HF",
        "JESRegrouped_RelativeBal",
        "JESRegrouped_RelativeSample_$CAMPAIGN",
        "JESTotal",
        "JEReta0to1p93",
        "JEReta1p93to2p5",
        "JERTotal",
        "UES",
        "ele_Reco",
        "ele_ID",
        "ele_Trig",
        "ElectronScale",
        "ElectronSmear",
        "mu_ID",
        "mu_Iso",
        "mu_Trig",
        "MuonScale",
        "MuonResol",
    ],
    help="Experimental uncertainties to include in the plots.",
)
parser.add_argument(
    "-t",
    "--th_systs",
    action="append",
    default=[],
    choices=[
        "ttbar_weight",
        "PDF_weight",
        "aS_weight",
        "scalevar_muF",
        "scalevar_muR",
        "UEPS_ISR",
        "UEPS_FSR",
    ],
    help="Theory uncertainties to include in the plots.",
)
parser.add_argument(
    "--lumi",
    type=float,
    required=True,
    help="Luminosity in /pb.",
)
parser.add_argument(
    "--mergemap",
    type=str,
    default=sample_mergemap,
    help="Specify mergemap as dict: `{merge1:[dataset1,dataset2]...}`. Also works with a json file containing a dict.",
)
parser.add_argument(
    "--log", action="store_true", help="Save an additional plot with log scale on y axis (default = False)."
)


exp_systs = {
    "puweight": "Pileup",
    "JESRegrouped_Absolute_$CAMPAIGN": "JES Absolute ($CAMPAIGN)",
    "JESRegrouped_Absolute": "JES Absolute",
    "JESRegrouped_BBEC1_$CAMPAIGN": "JES BBEC1 ($CAMPAIGN)",
    "JESRegrouped_BBEC1": "JES BBEC1",
    "JESRegrouped_EC2_$CAMPAIGN": "JES EC2 ($CAMPAIGN)",
    "JESRegrouped_EC2": "JES EC2",
    "JESRegrouped_FlavorQCD": "JES FlavorQCD",
    "JESRegrouped_HF_$CAMPAIGN": "JES HF ($CAMPAIGN)",
    "JESRegrouped_HF": "JES HF",
    "JESRegrouped_RelativeBal": "JES RelativeBal",
    "JESRegrouped_RelativeSample_$CAMPAIGN": "JES RelativeSample ($CAMPAIGN)",
    "JESTotal": "JES Total",
    "JEReta0to1p93": r"JER (0.0 $\leq$ $\eta$ < 1.93)",
    "JEReta1p93to2p5": r"JER (1.93 $\leq$ $\eta$ < 2.5)",
    "JERTotal": "JER Total",
    "UES": "MET Uncl. Energy",
    "ele_Reco": "Electron Reco.",
    "ele_ID": "Electron ID",
    "ele_Trig": "Electron Trigger",
    "ElectronScale": "Electron Scale",
    "ElectronSmear": "Electron Smear",
    "mu_ID": "Muon ID",
    "mu_Iso": "Muon Iso.",
    "mu_Trig": "Muon Trigger",
    "MuonScale": "Muon Scale",
    "MuonResol": "Muon Resol.",
}
th_systs = {
    "ttbar_weight": r"Top $p_T$",
    "PDF_weight": "PDF",
    "aS_weight": r"$\alpha_S$",
    "scalevar_muF": r"$\mu_F$",
    "scalevar_muR": r"$\mu_R$",
    "UEPS_ISR": "PS ISR",
    "UEPS_FSR": "PS FSR",
}
colors10 = {
    "puweight": "#000000",  # black
    "JEReta0to1p93": "#3f90da",
    "JEReta1p93to2p5": "#ffa90e",
    "UES": "#bd1f01",
    "ele_Reco": "#94a4a2",
    "ele_ID": "#832db6",
    "ele_Trig": "#a96b59",
    "ElectronScale": "#e76300",
    "ElectronSmear": "#b9ac70",
    "mu_ID": "#717581",
    "mu_Iso": "#92dadd",
    "mu_Trig": "#d900b2",  # magenta
    "MuonScale": "#05a30f",  # green
    "MuonResol": "#0a00ff",  # blue
}
colorsJES = {
    "JESTotal": "#3f90da",
    "JESRegrouped_Absolute_$CAMPAIGN": "#3f90da",
    "JESRegrouped_Absolute": "#ffa90e",
    "JESRegrouped_BBEC1_$CAMPAIGN": "#bd1f01",
    "JESRegrouped_BBEC1": "#94a4a2",
    "JESRegrouped_EC2_$CAMPAIGN": "#832db6",
    "JESRegrouped_EC2": "#a96b59",
    "JESRegrouped_FlavorQCD": "#e76300",
    "JESRegrouped_HF_$CAMPAIGN": "#b9ac70",
    "JESRegrouped_HF": "#717581",
    "JESRegrouped_RelativeBal": "#92dadd",
    "JESRegrouped_RelativeSample_$CAMPAIGN": "#d900b2",
}
colors8 = {
    "ttbar_weight": "#1845fb",
    "PDF_weight": "#ff5e02",
    "aS_weight": "#c91f16",
    "scalevar_muF": "#c849a9",
    "scalevar_muR": "#adad7d",
    "UEPS_ISR": "#86c8dd",
    "UEPS_FSR": "#578dff",
    # "": "#656364",
}
colors6 = ["#5790fc", "#f89c20", "#e42536", "#964a8b", "#9c9ca1", "#7a21dd"]


def jes_year(campaign):
    if campaign == "Summer22":
        return "2022"
    if campaign == "Summer22EE":
        return "2022EE"
    if campaign == "Summer23":
        return "2023"
    if campaign == "Summer23BPix":
        return "2023BPix"
    if campaign == "Summer24":
        return "2024"
    if campaign == "Prompt25":
        return "2025"
    return campaign


def do_osss(osss, h):  # For W+c workflows
    if osss == 1:  # OS
        h = h[{"osss": 0}]
    elif osss == -1:  # SS
        h = h[{"osss": 1}]
    else:  # OS - SS
        h = h[{"osss": 0}] + h[{"osss": 1}] * -1
    return h


def hide_inf_nan(arr, side):
    if side != "up" and side != "down":
        print("WARNING: incorrect option for parameter `side` in `hide_inf_nan()`!")

    for i in range(len(arr)):
        if math.isnan(arr[i]):
            arr[i] = 100.0 if side == "up" else -100.0  # dummy value
        if not np.isfinite(arr[i]):
            arr[i] = 100.0 if side == "up" else -100.0  # dummy value
    return arr


def get_error_magnitude(nom, up, down):
    ret_up = np.zeros_like(nom)
    ret_down = np.zeros_like(nom)

    for i in range(len(nom)):
        if up[i] > nom[i]:
            ret_up[i] += up[i] - nom[i]
        else:
            ret_down[i] += nom[i] - up[i]
        if down[i] < nom[i]:
            ret_down[i] += nom[i] - down[i]
        else:
            ret_up[i] += down[i] - nom[i]

    return ret_up, ret_down


def plot_uncertainty(ax, h_nom, unc_up, unc_down, syst_name, col):
    ax.hist(
        x=h_nom.axes[0].edges[:-1],
        bins=h_nom.axes[0].edges,
        weights=unc_up,
        label=syst_name,
        histtype="step",
        color=col,
    )
    ax.hist(
        x=h_nom.axes[0].edges[:-1],
        bins=h_nom.axes[0].edges,
        weights=unc_down,
        histtype="step",
        color=col,
        linestyle=":",
    )


def plot_band(ax, unc_up, unc_down, edges, label):
    ax.stairs(
        values=unc_up,
        baseline=unc_down,
        edges=edges,
        label=label,
        lw=0,
        color="#d3d3d3",
        alpha=0.4,
        fill=True,
    )


def plot_systs(
    inputs,
    variable,
    workflow,
    campaign,
    split,
    mergemap,
    osss,
    lumi,
    log,
    outdir,
    scale=1.0,
):
    inputs = scaleSumW(inputs, lumi)
    mergemap["data"] = [m for m in inputs.keys() if "Run" in m]
    collated = collate(inputs, mergemap)
    collated = {
        key: value for key, value in collated.items() if isinstance(collated[key], dict)
    }
    print("Samples:", collated.keys())

    h_mc = collated["mc"][variable]
    h_data = collated["data"][variable]

    if scale != 1.0 and scale != None:
        print("INFO: scaling MC histograms by %.4f." % scale)

    if split == "flavour" and "flav" not in h_mc.axes.name:
        print(f"WARNING: {variable} is not split by flavour -- not making the plot.")
        return

    # For W+c workflows
    if "_Wc_" in workflow:
        h_mc = do_osss(osss, h_mc)
        h_data = do_osss(osss, h_data)

    # Build base axis used for nominal and systematic variations
    axis = {}
    if "flav" in h_mc.axes.name:
        axis["flav"] = sum

    # Build nominal axis
    axis_nom = copy.deepcopy(axis)
    axis_nom["syst"] = "nominal"

    # Get nominal MC and data histograms
    h_mc_nom = h_mc[axis_nom]
    if "top_pt" not in var:  # data is empty for top_pt and antitop_pt
        h_data = h_data[axis_nom]

    # Create figure
    fig, axs = plt.subplots(
        5,  # Y
        1,  # X
        figsize=(10, 16),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1, 1, 1, 1]},
    )
    fig.subplots_adjust(hspace=0.1)  # Adjust vertical space between Axes

    # Plot MC and data
    stack, labels, colors = [], [], []
    if split == "flavour":
        axis_flav = {"syst": "nominal"}
        for i, flav in enumerate(["udsg", "pu", "c", "b"]):
            axis_flav["flav"] = i
            stack.append(h_mc[axis_flav] * scale)
            labels.append(flav)
            colors.append(color_map[flav])
    else:  # Split plots by sample
        for s in collated.keys():
            if s == "mc" or s == "data":
                continue
            if len(collated[s][variable].values()) == 0:
                print(
                    f"WARNING: {variable} histogram for sample {s} is empty! Continuing..."
                )
                continue
            h_sample = collated[s][variable][axis_nom]
            if "_Wc_" in workflow:  # For W+c workflows
                h_sample = do_osss(osss, h_sample)
            stack.append(h_sample * scale)
            labels.append(s)
            colors.append(color_map[s])
    hep.histplot(
        stack,
        histtype="fill",
        stack=True,
        label=labels,
        yerr=False,
        ax=axs[0],
        color=colors,
        sort="y",
        flow="none",
    )
    if "top_pt" not in variable:
        hep.histplot(
            h_data,
            histtype="errorbar",
            color="black",
            label="Data",
            yerr=True,
            ax=axs[0],
            xerr=False,
            flow="none",
        )

    # Add black line to ratio plots
    axs[1].plot(
        [h_mc_nom.axes[0].edges[0], h_mc_nom.axes[0].edges[-1]],
        [1.0, 1.0],
        linestyle=":",
        color="black",
    )
    axs[2].plot(
        [h_mc_nom.axes[0].edges[0], h_mc_nom.axes[0].edges[-1]],
        [1.0, 1.0],
        linestyle=":",
        color="black",
    )
    axs[3].plot(
        [h_mc_nom.axes[0].edges[0], h_mc_nom.axes[0].edges[-1]],
        [1.0, 1.0],
        linestyle=":",
        color="black",
    )
    axs[4].plot(
        [h_mc_nom.axes[0].edges[0], h_mc_nom.axes[0].edges[-1]],
        [1.0, 1.0],
        linestyle=":",
        color="black",
    )

    total_up = np.zeros_like(h_mc_nom.values())
    total_down = np.zeros_like(h_mc_nom.values())

    ###################################
    # Reduced split JES uncertainties #
    ###################################

    jes_total_up = np.zeros_like(h_mc_nom.values())
    jes_total_down = np.zeros_like(h_mc_nom.values())

    for syst in args.exp_systs:
        if "JES" not in syst:
            continue

        syst_campaign = syst.replace("$CAMPAIGN", jes_year(campaign))

        if (
            f"{syst_campaign}Up" not in h_mc.axes["syst"]
            or f"{syst_campaign}Down" not in h_mc.axes["syst"]
        ):
            print(
                f"WARNING: {syst_campaign}Up or {syst_campaign}Down are not available!"
            )
            print(f'         Available systematics are: {h_mc.axes["syst"]}')
            continue

        axis_up = copy.deepcopy(axis)
        axis_up["syst"] = f"{syst_campaign}Up"
        axis_down = copy.deepcopy(axis)
        axis_down["syst"] = f"{syst_campaign}Down"

        h_up = h_mc[axis_up]
        h_down = h_mc[axis_down]

        rat_up = hide_inf_nan(h_up.values() / h_mc_nom.values(), "up")
        rat_down = hide_inf_nan(h_down.values() / h_mc_nom.values(), "down")
        plot_uncertainty(
            axs[2],
            h_mc_nom,
            rat_up,
            rat_down,
            exp_systs[syst].replace("$CAMPAIGN", jes_year(campaign)),
            colorsJES[syst],
        )

        syst_up, syst_down = get_error_magnitude(
            h_mc_nom.values(), h_up.values(), h_down.values()
        )
        jes_total_up += syst_up * syst_up
        jes_total_down += syst_down * syst_down

    total_up += jes_total_up
    total_down += jes_total_down
    jes_total_up = np.sqrt(jes_total_up)
    jes_total_down = np.sqrt(jes_total_down)
    jes_up = hide_inf_nan(1.0 + jes_total_up / h_mc_nom.values(), "up")
    jes_down = hide_inf_nan(1.0 - jes_total_down / h_mc_nom.values(), "down")
    plot_band(axs[2], jes_up, jes_down, h_mc_nom.axes[0].edges, "JES Total")

    ##############################
    # Experimental uncertainties #
    ##############################

    exp_total_up = jes_total_up * jes_total_up
    exp_total_down = jes_total_down * jes_total_down

    for syst in args.exp_systs:
        if "JES" in syst:
            continue

        if (
            f"{syst}Up" not in h_mc.axes["syst"]
            or f"{syst}Down" not in h_mc.axes["syst"]
        ):
            print(f"WARNING: {syst}Up or {syst}Down are not available!")
            print(f'         Available systematics are: {h_mc.axes["syst"]}')
            continue

        axis_up = copy.deepcopy(axis)
        axis_up["syst"] = f"{syst}Up"
        axis_down = copy.deepcopy(axis)
        axis_down["syst"] = f"{syst}Down"

        h_up = h_mc[axis_up]
        h_down = h_mc[axis_down]

        rat_up = hide_inf_nan(h_up.values() / h_mc_nom.values(), "up")
        rat_down = hide_inf_nan(h_down.values() / h_mc_nom.values(), "down")
        plot_uncertainty(
            axs[1], h_mc_nom, rat_up, rat_down, exp_systs[syst], colors10[syst]
        )

        syst_up, syst_down = get_error_magnitude(
            h_mc_nom.values(), h_up.values(), h_down.values()
        )
        exp_total_up += syst_up * syst_up
        exp_total_down += syst_down * syst_down

    total_up += exp_total_up
    total_down += exp_total_down
    exp_total_up = np.sqrt(exp_total_up)
    exp_total_down = np.sqrt(exp_total_down)
    exp_up = hide_inf_nan(1.0 + exp_total_up / h_mc_nom.values(), "up")
    exp_down = hide_inf_nan(1.0 - exp_total_down / h_mc_nom.values(), "down")
    plot_band(axs[1], exp_up, exp_down, h_mc_nom.axes[0].edges, "Exp. Total")
    plot_uncertainty(axs[4], h_mc_nom, exp_up, exp_down, "Experimental", colors6[0])

    #############################
    # Theoretical uncertainties #
    #############################

    th_total_up = np.zeros_like(h_mc_nom.values())
    th_total_down = np.zeros_like(h_mc_nom.values())

    for syst in args.th_systs:
        if (
            f"{syst}Up" not in h_mc.axes["syst"]
            or f"{syst}Down" not in h_mc.axes["syst"]
        ):
            print(f"WARNING: {syst}Up or {syst}Down are not available!")
            print(f'         Available systematics are: {h_mc.axes["syst"]}')
            continue

        axis_up = copy.deepcopy(axis)
        axis_up["syst"] = f"{syst}Up"
        axis_down = copy.deepcopy(axis)
        axis_down["syst"] = f"{syst}Down"

        h_up = h_mc[axis_up]
        h_down = h_mc[axis_down]
        if "PDF" in syst:
            h_up = collated["mc"][f"{variable}_PDF_weightUp"][axis_up]
            h_down = collated["mc"][f"{variable}_PDF_weightDown"][axis_down]
        elif "aS" in syst:
            h_up = collated["mc"][f"{variable}_aS_weightUp"][axis_up]
            h_down = collated["mc"][f"{variable}_aS_weightDown"][axis_down]
        elif "muF" in syst:
            h_up = collated["mc"][f"{variable}_scalevar_muFUp"][axis_up]
            h_down = collated["mc"][f"{variable}_scalevar_muFDown"][axis_down]
        elif "muR" in syst:
            h_up = collated["mc"][f"{variable}_scalevar_muRUp"][axis_up]
            h_down = collated["mc"][f"{variable}_scalevar_muRDown"][axis_down]
        elif "ISR" in syst:
            h_up = collated["mc"][f"{variable}_UEPS_ISRUp"][axis_up]
            h_down = collated["mc"][f"{variable}_UEPS_ISRDown"][axis_down]
        elif "FSR" in syst:
            h_up = collated["mc"][f"{variable}_UEPS_FSRUp"][axis_up]
            h_down = collated["mc"][f"{variable}_UEPS_FSRDown"][axis_down]
        # For W+c workflows
        if "_Wc_" in workflow:
            h_up = do_osss(osss, h_up)
            h_down = do_osss(osss, h_down)

        rat_up = hide_inf_nan(h_up.values() / h_mc_nom.values(), "up")
        rat_down = hide_inf_nan(h_down.values() / h_mc_nom.values(), "down")
        plot_uncertainty(
            axs[3], h_mc_nom, rat_up, rat_down, th_systs[syst], colors8[syst]
        )

        syst_up, syst_down = get_error_magnitude(
            h_mc_nom.values(), h_up.values(), h_down.values()
        )
        th_total_up += syst_up * syst_up
        th_total_down += syst_down * syst_down

    total_up += th_total_up
    total_down += th_total_down
    th_total_up = np.sqrt(th_total_up)
    th_total_down = np.sqrt(th_total_down)
    th_up = hide_inf_nan(1.0 + th_total_up / h_mc_nom.values(), "up")
    th_down = hide_inf_nan(1.0 - th_total_down / h_mc_nom.values(), "down")
    plot_band(axs[3], th_up, th_down, h_mc_nom.axes[0].edges, "Th. Total")
    plot_uncertainty(axs[4], h_mc_nom, th_up, th_down, "Theoretical", colors6[2])

    ###########################
    # Statistical uncertainty #
    ###########################

    stat_total = h_mc_nom.variances()

    total_up += stat_total
    total_down += stat_total
    stat_total = np.sqrt(stat_total)
    stat_up = hide_inf_nan(1.0 + stat_total / h_mc_nom.values(), "up")
    stat_down = hide_inf_nan(1.0 - stat_total / h_mc_nom.values(), "down")
    plot_uncertainty(axs[4], h_mc_nom, stat_up, stat_down, "MC Stat.", colors6[5])

    ########################
    # All MC uncertainties #
    ########################

    total_up = np.sqrt(total_up)
    total_down = np.sqrt(total_down)
    total_up_ratio = hide_inf_nan(1.0 + total_up / h_mc_nom.values(), "up")
    total_down_ratio = hide_inf_nan(1.0 - total_down / h_mc_nom.values(), "down")
    plot_band(axs[4], total_up_ratio, total_down_ratio, h_mc_nom.axes[0].edges, "Total")

    # Add total MC uncertainty bars to MC histogram
    MCerrorband(
        h_mc_nom * scale,
        ax=axs[0],
        ext_error=[total_up, total_down],
        label=r"MC stat. $\oplus$ syst.",
        flow=None,
        fill_opts=errband_opts,
    )
    # Add data/MC to bottom plot
    syst_up = np.sqrt(exp_total_up * exp_total_up + th_total_up * th_total_up)
    syst_down = np.sqrt(exp_total_down * exp_total_down + th_total_down * th_total_down)
    if "top_pt" not in variable:
        plotratio(
            h_data,
            h_mc_nom * scale,
            ax=axs[4],
            flow=None,
            xerr=False,
            ext_denom_error=np.array([syst_up, syst_down]),
            clear=False,
            label="Data/MC",
            denom_fill_opts={"lw": 0},
        )

    #######################
    # Plot customisations #
    #######################

    # Set legend font sizes and columns
    axs[0].legend(fontsize=15, ncols=3)
    axs[1].legend(fontsize=12, ncols=3)
    axs[2].legend(fontsize=12, ncols=3)
    axs[3].legend(fontsize=12, ncols=3)
    axs[4].legend(fontsize=12, ncols=3)

    # Add additional plot info
    input_txt = ""
    if "Wc" in workflow:
        input_txt += "W+c"
        if osss == 1:
            input_txt += " OS"
        elif osss == -1:
            input_txt += " SS"
        else:
            input_txt += " OS-SS"
    elif "DY" in workflow:
        input_txt += "DY+jets"
    elif "ttdilep" in workflow:
        input_txt += r"t$\bar{t}$ dileptonic"
    elif "ttsemilep" in workflow:
        input_txt += r"t$\bar{t}$ semileptonic"
    if "_e_" in workflow:
        if "DY" in workflow:
            input_txt += " (ee)"
        else:
            input_txt += " (e)"
    elif "_mu_" in workflow:
        if "DY" in workflow:
            input_txt += r" ($\mu\mu$)"
        else:
            input_txt += r" ($\mu$)"
    elif "_emu_" in workflow:
        input_txt += r" (e$\mu$)"
    if "2Dbin_pt" in variable or "HFvLF_pt" in variable or "BvC_pt" in variable:
        tmp = variable.strip().split("_")[1][2:]
        tmp = tmp.split("to")
        low_edge = tmp[0]
        up_edge = tmp[1]
        input_txt += "\n" + low_edge + r" < jet $p_T$ < " + up_edge + " GeV"
    if scale != 1.0:
        input_txt += "\n" + r"MC$\times%.4f$" % scale
    at = AnchoredText(
        input_txt,
        loc=2,
        frameon=False,
        prop=dict(fontsize=15),
    )
    axs[0].add_artist(at)

    # Find borders of x axis
    bin_min = None
    bin_max = None
    for i in range(len(h_mc_nom.values())):
        if h_mc_nom.values()[i] != 0.0:
            bin_min = i
            break
    for i in range(len(h_mc_nom.values())):
        if h_mc_nom.values()[len(h_mc_nom.values()) - (i + 1)] != 0.0:
            bin_max = len(h_mc_nom.values()) - i
            break

    # Find maximum y value
    max_y = max(
        max(h_mc_nom.values() + total_up),
        (
            max(h_data.values() + np.sqrt(h_data.variances()))
            if "top_pt" not in variable
            else -100.0
        ),
    )

    # Set x and y axis limits
    axs[4].set_xlim(
        (
            h_mc_nom.axes[0].edges[0]
            if bin_min is None
            else h_mc_nom.axes[0].edges[bin_min]
        ),
        (
            h_mc_nom.axes[0].edges[-1]
            if bin_max is None
            else h_mc_nom.axes[0].edges[bin_max]
        ),
    )
    axs[1].set_ylim(0.9, 1.1)
    axs[2].set_ylim(0.9, 1.1)
    axs[3].set_ylim(0.9, 1.1)
    axs[4].set_ylim(0.5, 1.5)

    # Set labels
    hep.cms.label("Private Work", data=True, lumi=lumi / 1000.0, com=13.6, ax=axs[0])
    if "2Dbin" in variable:
        axs[4].set_xticks(
            [0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5],
            ["L0", "C0", "C1", "C2", "C3", "C4", "B0", "B1", "B2", "B3", "B4"],
        )
    if variable == "npv":
        axs[4].locator_params(axis="x", nbins=10)
    axs[0].set_xlabel("")
    axs[4].set_xlabel(axes_name(variable))
    axs[0].set_ylabel("Events")
    axs[1].set_ylabel("Exp. Unc.")
    axs[2].set_ylabel("JES Unc.")
    axs[3].set_ylabel("Th. Unc.")
    axs[4].set_ylabel("Total Unc.")

    outname = f"{outdir}/systs_{var}_"
    if split == "flavour":
        outname += "flav"
    else:
        outname += "sample"
    if scale != 1.0:
        outname += "_scale"

    axs[0].set_ylim(top=1.35 * max_y)
    # https://cms-analysis.docs.cern.ch/guidelines/plotting/general/#scientific-notation
    # Scientific notation
    axs[0].ticklabel_format(axis="y", style="sci", scilimits=(-3, 3), useMathText=True)
    # Shift multiplier position out
    axs[0].get_yaxis().get_offset_text().set_position((-0.085, 1.05))
    print(f"Saving {outname}.pdf")
    fig.savefig(f"{outname}.pdf", bbox_inches="tight", pad_inches=0.3)
    fig.savefig(f"{outname}.png", bbox_inches="tight", pad_inches=0.3)

    if log:
        axs[0].set_yscale("log")
        axs[0].set_ylim(
            bottom=0.1,
            top=60.0 * max_y,
        )
        print(f"Saving {outname}_log.pdf")
        fig.savefig(f"{outname}_log.pdf", bbox_inches="tight", pad_inches=0.3)
        fig.savefig(f"{outname}_log.png", bbox_inches="tight", pad_inches=0.3)

    return sum(h_data.values()) / sum(h_mc_nom.values())


if __name__ == "__main__":

    args = parser.parse_args()

    # Get list of coffea files
    if len(args.input.split(",")) > 1 and "*" in args.input:
        directories = args.input.split(",")
        files = []
        for directory in directories:
            files.extend(glob.glob(directory))
    elif len(args.input.split(",")) > 1:
        files = args.input.split(",")
    elif "*" in args.input:
        files = glob.glob(args.input)
    else:
        files = [args.input]
    inputs = {i: load(i) for i in files if ".coffea" in i}

    # Get mergemap
    if type(args.mergemap) == str and args.mergemap.endswith(".json"):
        mergemap = json.load(open(args.mergemap))
    else:
        mergemap = args.mergemap
    # Remove extra mergemap key(s)
    del mergemap["DY"]

    datalist = []
    mclist = []
    for f in inputs.keys():
        datalist.extend([m for m in inputs[f].keys() if "Run" in m])
        mclist.extend([m for m in inputs[f].keys() if "Run" not in m])
    mergemap["mc"] = mclist
    mergemap["data"] = datalist

    vars = args.vars.strip().split(",")

    outdir = f"{args.outdir}/{args.campaign}/{args.workflow}_/"
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    for var in vars:
        data_over_mc = plot_systs(
            inputs,
            var,
            args.workflow,
            args.campaign,
            args.split,
            mergemap,
            args.splitOSSS,
            args.lumi,
            args.log,
            outdir,
        )
        if "top_pt" not in var:
            plot_systs(
                inputs,
                var,
                args.workflow,
                args.campaign,
                args.split,
                mergemap,
                args.splitOSSS,
                args.lumi,
                False,
                outdir,
                data_over_mc,
            )
