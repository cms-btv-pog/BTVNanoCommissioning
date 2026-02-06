from .axes.common import axes as common_axes

from .histograms.common import get_histograms as common_hists
from .histograms.ctag_ttdilep import get_histograms as ctag_ttdilep_hists
from .histograms.ctag_ttsemilep import get_histograms as ctag_ttsemilep_hists
from .histograms.dy import get_histograms as dy_hists
from .histograms.dy_sfl import get_histograms as dy_sfl_hists
from .histograms.example import get_histograms as example_hists
from .histograms.fourvec import get_histograms as fourvec_hists
from .histograms.ttdilep_kin import get_histograms as ttdilep_kin_hists
from .histograms.ttdilep import get_histograms as ttdilep_hists
from .histograms.ttsemilep import get_histograms as ttsemilep_hists
from .histograms.qcd import get_histograms as qcd_hists
from .histograms.qcd_smu import get_histograms as qcd_smu_hists
from .histograms.validation import get_histograms as validation_hists
from .histograms.wc import get_histograms as wc_hists


def get_axes_collections(
    axes_list: list = None,
):
    available_collections = {
        "common": common_axes,
    }

    output = {}

    if axes_list is None:
        for ax in available_collections:
            output = output | available_collections[ax]
    else:
        for ax in axes_list:
            if ax not in available_collections:
                raise ValueError(
                    f"Axis {ax} not found in {available_collections.keys()}. Check utils/histogramming/axes"
                )
            output = output | available_collections[ax]

    return output


def get_hist_collections(axes: dict, hist_collections: list, **kwargs):
    available_collections = {
        "example": example_hists,
        "common": common_hists,
        "fourvec": fourvec_hists,
        "QCD": qcd_hists,
        "QCD_smu": qcd_smu_hists,
        "ctag_ttdilep": ctag_ttdilep_hists,
        "ctag_ttsemilep": ctag_ttsemilep_hists,
        "DY": dy_hists,
        "DY_sfl": dy_sfl_hists,
        "ttdilep_kin": ttdilep_kin_hists,
        "ttdilep": ttdilep_hists,
        "ttsemilep": ttsemilep_hists,
        "validation": validation_hists,
        "Wc": wc_hists,
    }

    output = {}

    for h in hist_collections:
        if h not in available_collections:
            raise ValueError(
                f"Histogram collection {h} not found in {available_collections.keys()}. Check utils/histogramming/histograms"
            )
        output = output | available_collections[h](axes, **kwargs)

    return output
