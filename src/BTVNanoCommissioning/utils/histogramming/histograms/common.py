import hist as Hist

from BTVNanoCommissioning.utils.histogramming import hist_helpers
from BTVNanoCommissioning.helpers.definitions import get_definitions, get_discriminators


def get_histograms(axes, **kwargs):
    hists = {}

    include_osss = kwargs.get("include_osss", False)

    include_npv = kwargs.get("include_npv", True)
    include_njet = kwargs.get("include_njet", True)
    include_nmujet = kwargs.get("include_nmujet", False)
    include_nsoftmu = kwargs.get("include_nsoftmu", False)

    include_discriminators = kwargs.get("include_discriminators", True)
    include_common_definitions = kwargs.get("include_common_definitions", True)
    include_sv_variables = kwargs.get("include_sv_variables", False)

    if include_osss:
        n_axes = [axes["syst"], axes["osss"], axes["n"]]
        npv_axes = [axes["syst"], axes["osss"], axes["npv"]]
    else:
        n_axes = [axes["syst"], axes["n"]]
        npv_axes = [axes["syst"], axes["npv"]]

    if include_npv:
        hists["npv"] = Hist.Hist(*npv_axes, Hist.storage.Weight())

    if include_njet:
        hists["njet"] = Hist.Hist(*n_axes, Hist.storage.Weight())

    if include_nmujet:
        hists["nmujet"] = Hist.Hist(*n_axes, Hist.storage.Weight())

    if include_nsoftmu:
        hists["nsoftmu"] = Hist.Hist(*n_axes, Hist.storage.Weight())

    if include_discriminators:
        hists = hists | _get_discriminators(axes, **kwargs)

    if include_common_definitions:
        hists = hists | _get_common_definitions(axes, **kwargs)

    if include_sv_variables:
        hists = hists | _get_sv_variables(axes, **kwargs)

    return hists


def _get_common_definitions(axes, **kwargs):

    include_osss = kwargs.get("include_osss", False)
    jet_fields = kwargs.get("jet_fields", None)

    if jet_fields is None:
        raise ValueError(
            "Tried to get btag input histograms without passing jet_fields (list of keys in the jet collection)"
        )

    hists = {}

    definitions = get_definitions()

    for d in definitions:
        if d not in jet_fields:
            continue

        ranges = definitions[d]["manual_ranges"]
        binning = definitions[d]["bins"]
        labels = (
            definitions[d]["displayname"]
            + " ["
            + definitions[d]["inputVar_units"]
            + "]"
            if definitions[d]["inputVar_units"] is not None
            else definitions[d]["displayname"]
        )

        if include_osss:
            ax = [
                axes["syst"],
                axes["flav"],
                axes["osss"],
                Hist.axis.Regular(binning, ranges[0], ranges[1], name=d, label=labels),
            ]
        else:
            ax = [
                axes["syst"],
                axes["flav"],
                Hist.axis.Regular(binning, ranges[0], ranges[1], name=d, label=labels),
            ]

        hists[d] = Hist.Hist(*ax, Hist.storage.Weight())

    return hists


# NOTE: This was commented out in the old histogrammer
# Disabled by default in get_histograms()
def _get_sv_variables(axes, **kwargs):

    include_osss = kwargs.get("include_osss", False)

    hists = {}

    definitions = hist_helpers.get_definitions(include_defitions=["SV"])

    for d in definitions:
        ranges = definitions[d]["manual_ranges"]
        binning = definitions[d]["bins"]
        labels = (
            definitions[d]["displayname"]
            + " ["
            + definitions[d]["inputVar_units"]
            + "]"
            if definitions[d]["inputVar_units"] is not None
            else definitions[d]["displayname"]
        )

        if include_osss:
            ax = [
                axes["syst"],
                axes["flav"],
                axes["osss"],
                Hist.axis.Regular(binning, ranges[0], ranges[1], name=d, label=labels),
            ]
        else:
            ax = [
                axes["syst"],
                axes["flav"],
                Hist.axis.Regular(binning, ranges[0], ranges[1], name=d, label=labels),
            ]

        hists[d] = Hist.Hist(*ax, Hist.storage.Weight())

    return hists


def _get_discriminators(axes, **kwargs):
    jet_fields = kwargs.get("jet_fields", None)

    if jet_fields is None:
        raise ValueError(
            "Tried to get discriminator histograms without passing jet_fields (list of keys in the jet collection)"
        )

    include_osss = kwargs.get("include_osss", False)
    njet = kwargs.get("njet", 1)
    c_wf = kwargs.get("c_wf", False)

    disc_list = get_discriminators()

    hists = {}

    common_axes = [axes["syst"], axes["flav"]]
    if include_osss:
        common_axes.append(axes["osss"])

    for d in disc_list:
        if d not in jet_fields:
            continue
        disc_axes = {
            "btag": Hist.axis.Regular(50, 0.0, 1, name="discr", label=d),
            "Bprob": Hist.axis.Regular(50, 0, 10, name="discr", label=d),
            "Res": Hist.axis.Regular(40, 0, 1, name="discr", label=d),
            "Corr": Hist.axis.Regular(40, 0, 2, name="discr", label=d),
        }

        if c_wf:
            if "btag" in d or "ProbaN" == d:
                hists[d] = Hist.Hist(
                    *common_axes,
                    disc_axes["btag"],
                    Hist.storage.Weight(),
                )
            elif "Bprob" in d:
                hists[d] = Hist.Hist(
                    *common_axes,
                    disc_axes["Bprob"],
                    Hist.storage.Weight(),
                )
            elif "PNetRegPtRawRes" == d:
                hists[d] = Hist.Hist(
                    *common_axes,
                    disc_axes["PNetRegPtRawRes"],
                    Hist.storage.Weight(),
                )
            elif "PNetRegPtRawCorr" in d:
                hists[d] = Hist.Hist(
                    *common_axes,
                    disc_axes["PNetRegPtRawCorr"],
                    Hist.storage.Weight(),
                )

        for i in range(njet):
            if "btag" in d or "ProbaN" == d:
                hists[d] = Hist.Hist(
                    *common_axes,
                    disc_axes["btag"],
                    Hist.storage.Weight(),
                )
            elif "Bprob" in d:
                hists[d] = Hist.Hist(
                    *common_axes,
                    disc_axes["Bprob"],
                    Hist.storage.Weight(),
                )

            if include_osss:
                if "Res" in d:
                    hists[d] = Hist.Hist(
                        *common_axes,
                        disc_axes["Res"],
                        Hist.storage.Weight(),
                    )
                elif "Corr" in d:
                    hists[d] = Hist.Hist(
                        *common_axes,
                        disc_axes["Corr"],
                        Hist.storage.Weight(),
                    )

            else:
                if "PNetRegPtRawRes" == d:
                    hists[d] = Hist.Hist(
                        *common_axes,
                        disc_axes["Res"],
                        Hist.storage.Weight(),
                    )
                elif "PNetRegPtRawCorr" in d:
                    hists[d] = Hist.Hist(
                        *common_axes,
                        disc_axes["Corr"],
                        Hist.storage.Weight(),
                    )

    return hists
