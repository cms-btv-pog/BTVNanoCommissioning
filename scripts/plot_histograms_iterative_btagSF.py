from pathlib import Path

import matplotlib.pyplot as plt
import mplhep as hep
from coffea.util import load
from matplotlib import rc_context

from BTVNanoCommissioning.helpers.xs_scaler import collate, scaleSumW
from BTVNanoCommissioning.utils.plot_utils import MCerrorband, color_map, plotratio

lumi = 9451  # pb^-1
# need to devide, mplhep uses fb^-1
lumi_label = lumi / 1000  # fb^-1
com = 13.6

# load the coffea file
# mc_path = "/home/home1/institut_3a/zinn/analyses/BTV/BTVNanoCommissioning/hists_btag_ttbar_sf_mumu_MC_Summer23_2023_BTV_Run3_2023_Comm_MINIAODv4_NanoV12/hists_btag_ttbar_sf_mumu_MC_Summer23_2023_BTV_Run3_2023_Comm_MINIAODv4_NanoV12.coffea"
mc_path = "/home/home1/institut_3a/zinn/analyses/BTV/BTVNanoCommissioning/hists_btag_ttbar_sf_MC_Summer23BPix_2023_BTV_Run3_2023_Comm_MINIAODv4_NanoV12/hists_btag_ttbar_sf_MC_Summer23BPix_2023_BTV_Run3_2023_Comm_MINIAODv4_NanoV12.coffea"
mc_output = load(mc_path)

# data_muBTV_path = "/home/home1/institut_3a/zinn/analyses/BTV/BTVNanoCommissioning/hists_btag_ttbar_sf_mumu_data_Summer23_2023_mu_BTV_Run3_2023_Comm_MINIAODv4_NanoV12/hists_btag_ttbar_sf_mumu_data_Summer23_2023_mu_BTV_Run3_2023_Comm_MINIAODv4_NanoV12.coffea"
data_muBTV_path = "/home/home1/institut_3a/zinn/analyses/BTV/BTVNanoCommissioning/hists_btag_ttbar_sf_data_Summer23BPix_2023_mu_BTV_Run3_2023_Comm_MINIAODv4_NanoV12/hists_btag_ttbar_sf_data_Summer23BPix_2023_mu_BTV_Run3_2023_Comm_MINIAODv4_NanoV12.coffea"
data_muBTV_output = load(data_muBTV_path)
data_emBTV = "/home/home1/institut_3a/zinn/analyses/BTV/BTVNanoCommissioning/hists_btag_ttbar_sf_data_Summer23BPix_2023_em_BTV_Run3_2023_Comm_MINIAODv4_NanoV12/hists_btag_ttbar_sf_data_Summer23BPix_2023_em_BTV_Run3_2023_Comm_MINIAODv4_NanoV12.coffea"
data_emBTV_output = load(data_emBTV)

output = {
    "mc": mc_output,
    "data_muBTV": data_muBTV_output,
    # "data_emBTV": data_emBTV_output,
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
]

# this would be to BTV Framework plotting style
# flav_labels = {"c": 4, "b": 5, "udsg": 0, "pu": 1}
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


with rc_context(hep.style.CMS):
    for variable in variables:
        histogram = collated["mc"][variable]
        data_histogram = collated["data"][variable]
        # syst, flav, eta, pt, region, jet_index, discr

        for region in histogram.axes["region"]:
            flav_axis = histogram.axes["flav"]

            hist_list = [
                histogram[
                    "nominal", list(flav_axis.index(flav)), sum, sum, region, sum, :
                ][sum, :]
                for flav in flav_labels.values()
            ]

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
            rax.set_ylabel("Data/MC")

            save_path = Path("plots", variable)
            save_path.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path / f"{region}.png")
            ax.set_yscale("log")
            fig.savefig(save_path / f"{region}_log.png")
            plt.close(fig)

            # plot not stacked
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

            save_path = Path("plots", variable)
            save_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path / f"{region}_hists.png")
            ax.set_yscale("log")
            fig.savefig(save_path / f"{region}_hists_log.png")
            plt.close(fig)

            save_path_split = save_path / "eta_pt_bins"
            save_path_split.mkdir(parents=True, exist_ok=True)
            for eta_index, eta in enumerate(histogram.axes["eta"]):
                eta_bin = f"eta{eta[0]:.1f}to{eta[1]:.1f}"
                for pt_index, pt in enumerate(histogram.axes["pt"]):
                    pt_bin = f"pt{pt[0]:.0f}to{pt[1]:.0f}"

                    # plot for each eta and pt bin
                    fig, ax = plt.subplots()
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
                            ]
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
                        collated["data"][variable][
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
                    fig.tight_layout()

                    # save
                    fig.savefig(save_path_split / f"{region}_{eta_bin}_{pt_bin}.png")
                    plt.close(fig)
