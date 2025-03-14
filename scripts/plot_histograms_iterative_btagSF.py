from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import rc_context
import mplhep as hep
from coffea.util import load

from BTVNanoCommissioning.helpers.xs_scaler import collate, scaleSumW
from BTVNanoCommissioning.utils.plot_utils import MCerrorband, plotratio

lumi = 9451
com = 13.6

# load the coffea file
mc_path = "/home/home1/institut_3a/zinn/analyses/BTV/BTVNanoCommissioning/hists_btag_ttbar_sf_mumu_MC_Summer23BPix_2023_BTV_Run3_2023_Comm_MINIAODv4_NanoV12/hists_btag_ttbar_sf_mumu_MC_Summer23BPix_2023_BTV_Run3_2023_Comm_MINIAODv4_NanoV12.coffea"
mc_output = load(mc_path)

data_path = "/home/home1/institut_3a/zinn/analyses/BTV/BTVNanoCommissioning/hists_btag_ttbar_sf_mumu_data_Summer23BPix_2023_mu_BTV_Run3_2023_Comm_MINIAODv4_NanoV12/hists_btag_ttbar_sf_mumu_data_Summer23BPix_2023_mu_BTV_Run3_2023_Comm_MINIAODv4_NanoV12.coffea"
data_output = load(data_path)

output = {"mc": mc_output, "data": data_output}

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

flav_labels = ("udsg", "pu", "c", "b")

with rc_context(hep.style.CMS):
    for variable in variables:
        histogram = collated["mc"][variable]
        data_histogram = collated["data"][variable]
        # syst, flav, eta, pt, region, jet_index, discr

        for region in histogram.axes["region"]:
            flav_axis = histogram.axes["flav"]

            # plot sum over all eta and pt bins
            fig, (ax, ratio_ax) = plt.subplots(
                2,
                1,
                figsize=(10, 10),
                gridspec_kw={"height_ratios": (3, 1)},
                sharex=True,
            )
            fig.subplots_adjust(hspace=0.06, top=0.92, bottom=0.1, right=0.97)
            hep.cms.label("Preliminary", com=com, lumi=lumi, ax=ax, loc=0, data=True)
            hep.histplot(
                [
                    histogram["nominal", index, sum, sum, region, sum, :]
                    for index in range(len(flav_axis) - 1)
                ],
                label=flav_labels,
                stack=True,
                histtype="fill",
                yerr=True,
                ax=ax,
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

            save_path = Path("plots", variable, f"{region}.png")
            save_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_path)
            ax.set_yscale("log")
            save_path = Path("plots", variable, f"{region}_log.png")
            fig.savefig(save_path)
            plt.close(fig)

            for eta_index, eta in enumerate(histogram.axes["eta"]):
                eta_bin = f"eta{eta[0]:.1f}to{eta[1]:.1f}"
                for pt_index, pt in enumerate(histogram.axes["pt"]):
                    pt_bin = f"pt{pt[0]:.0f}to{pt[1]:.0f}"

                    # plot for each eta and pt bin
                    fig, ax = plt.subplots()
                    hep.cms.label("Preliminary", com=com)
                    hep.histplot(
                        [
                            histogram[
                                "nominal",
                                # flav_axis.index(flav),
                                index,
                                eta_index,
                                pt_index,
                                region,
                                sum,
                                :,
                            ]
                            # for flav in flav_axis
                            for index in range(len(flav_axis) - 1)
                        ],
                        label=flav_labels,
                        stack=True,
                        histtype="fill",
                        yerr=True,
                        ax=ax,
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
                    save_path = Path(
                        "plots", variable, f"{region}_{eta_bin}_{pt_bin}.png"
                    )
                    fig.savefig(save_path)
                    plt.close(fig)
