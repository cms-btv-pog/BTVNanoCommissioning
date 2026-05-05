from BTVNanoCommissioning.helpers.func import update
from BTVNanoCommissioning.utils.selection import btag_wp_dict
import numpy as np
import awkward as ak


def missing_branch(events, campaign=None):
    """
    Add missing branches or rename branches in the `events` object.

    This function adds missing branches or renames existing branches in the `events` object using the `missing_branch` parameter.

    Usage:
    Use the `hasattr` function to check for missing branches.

    Deprecated:
    The `add_jec` function is deprecated. Please use the `JME_shifts` function in the `correction` module instead.

    Example:
    ```python
    events.fixedGridRhoFastjetAll = (
        events.fixedGridRhoFastjetAll
        if hasattr(events, "fixedGridRhoFastjetAll")
        else events.Rho.fixedGridRhoFastjetAll
    )
    ```

    Parameters:
    events (coffea.nanoaodevents): The events object to update.
    missing_branch (str): The name of the missing branch to add or rename.

    Returns:
    events (coffea.nanoaodevents): Events with updated branches.

    """
    # Function implementation here
    if hasattr(events, "fixedGridRhoFastjetAll"):
        pass  # already a top-level branch, nothing to do
    elif hasattr(events, "Rho") and hasattr(events.Rho, "fixedGridRhoFastjetAll"):
        events["fixedGridRhoFastjetAll"] = events.Rho.fixedGridRhoFastjetAll
    else:
        raise AttributeError(
            "fixedGridRhoFastjetAll not found as top-level branch or under Rho collection. "
            "This branch is required for JEC corrections. "
            "Check that the input NanoAOD file contains Rho_fixedGridRhoFastjetAll."
        )
    ## calculate missing nodes

    ## calculate missing nodes
    if not hasattr(events.Jet, "btagDeepFlavB"):
        jets = events.Jet
        jets["btagDeepFlavB"] = (
            events.Jet.btagDeepFlavB_b
            + events.Jet.btagDeepFlavB_bb
            + events.Jet.btagDeepFlavB_lepb
        )
        events.Jet = update(
            events.Jet,
            {"btagDeepFlavB": jets.btagDeepFlavB},
        )
    if (
        hasattr(events.Jet, "btagDeepFlavCvL")
        and hasattr(events.Jet, "btagDeepFlavUDS")
        and not hasattr(events.Jet, "btagDeepFlavC")
    ):
        jets = events.Jet
        jets["btagDeepFlavC"] = (
            events.Jet.btagDeepFlavCvL / (1.0 - events.Jet.btagDeepFlavCvL)
        ) * (events.Jet.btagDeepFlavG + events.Jet.btagDeepFlavUDS)
        events.Jet = update(
            events.Jet,
            {"btagDeepFlavC": jets.btagDeepFlavC},
        )
    if hasattr(events.Jet, "btagDeepFlavCvB") and not hasattr(
        events.Jet, "btagDeepFlavC"
    ):
        jets = events.Jet
        jets["btagDeepFlavC"] = (
            events.Jet.btagDeepFlavCvB / (1.0 - events.Jet.btagDeepFlavCvB)
        ) * (events.Jet.btagDeepFlavB)
        events.Jet = update(
            events.Jet,
            {"btagDeepFlavC": jets.btagDeepFlavC},
        )
    if hasattr(events.Jet, "btagDeepFlavC") and not hasattr(
        events.Jet, "btagDeepFlavCvL"
    ):
        jets = events.Jet
        jets["btagDeepFlavCvL"] = np.maximum(
            np.minimum(
                np.where(
                    ((events.Jet.btagDeepFlavC / (1.0 - events.Jet.btagDeepFlavB)) > 0)
                    & (events.Jet.pt > 15),
                    (events.Jet.btagDeepFlavC / (1.0 - events.Jet.btagDeepFlavB)),
                    -1,
                ),
                0.999999,
            ),
            -1,
        )
        jets["btagDeepFlavCvB"] = np.maximum(
            np.minimum(
                np.where(
                    (
                        (
                            events.Jet.btagDeepFlavC
                            / (events.Jet.btagDeepFlavC + events.Jet.btagDeepFlavB)
                        )
                        > 0
                    )
                    & (events.Jet.pt > 15),
                    (
                        events.Jet.btagDeepFlavC
                        / (events.Jet.btagDeepFlavC + events.Jet.btagDeepFlavB)
                    ),
                    -1,
                ),
                0.999999,
            ),
            -1,
        )
        events.Jet = update(
            events.Jet,
            {
                "btagDeepFlavCvL": jets.btagDeepFlavCvL,
                "btagDeepFlavCvB": jets.btagDeepFlavCvB,
            },
        )
    if not hasattr(events.Jet, "btagPNetCvNotB") and hasattr(events.Jet, "btagPNetB"):
        jets = events.Jet
        jets["btagPNetCvNotB"] = (
            jets.btagPNetCvB * jets.btagPNetB / (1.0 - jets.btagPNetB) ** 2
        )
        events.Jet = update(
            events.Jet,
            {"btagPNetCvNotB": jets.btagPNetCvNotB},
        )
    if not hasattr(events.Jet, "btagRobustParTAK4CvNotB") and hasattr(
        events.Jet, "btagRobustParTAK4B"
    ):
        jets = events.Jet
        jets["btagRobustParTAK4CvNotB"] = (
            jets.btagRobustParTAK4CvB
            * jets.btagRobustParTAK4B
            / (1.0 - jets.btagRobustParTAK4B) ** 2
        )
        events.Jet = update(
            events.Jet,
            {"btagRobustParTAK4CvNotB": jets.btagRobustParTAK4CvNotB},
        )
    for tagger in ("UParTAK4",):
        bvc_name = f"btag{tagger}BvC"
        cvb_name = f"btag{tagger}CvB"
        if not hasattr(events.Jet, bvc_name) and hasattr(events.Jet, cvb_name):
            jets = events.Jet
            jets[bvc_name] = ak.where(
                jets[cvb_name] < 0.0,
                -1,
                1.0 - jets[cvb_name],
            )
            jets[f"{bvc_name}t"] = ak.where(
                jets[bvc_name] < 0.0,
                -1,
                1.0 - np.sqrt(1.0 - jets[bvc_name]),
            )
            for lo, hi in ((25.0, 35.0), (35.0, 50.0), (50.0, 70.0), (70.0, 90.0), (90.0, 120.0)):
                jets[f"{bvc_name}_pt{int(lo)}to{int(hi)}"] = ak.where(
                    (jets.pt < lo) | (jets.pt >= hi),
                    -2,
                    jets[bvc_name],
                )
            jets[f"{bvc_name}_pt120toinf"] = ak.where(
                jets.pt < 120.0,
                -2,
                jets[bvc_name],
            )

            update_map = {bvc_name: jets[bvc_name], f"{bvc_name}t": jets[f"{bvc_name}t"]}
            for lo, hi in ((25, 35), (35, 50), (50, 70), (70, 90), (90, 120)):
                update_map[f"{bvc_name}_pt{lo}to{hi}"] = jets[f"{bvc_name}_pt{lo}to{hi}"]
            update_map[f"{bvc_name}_pt120toinf"] = jets[f"{bvc_name}_pt120toinf"]
            events.Jet = update(events.Jet, update_map)

        hfvlf_name = f"btag{tagger}HFvLF"
        if "UParTAK4" in tagger:
            if (
                not hasattr(events.Jet, hfvlf_name)
                and hasattr(events.Jet, f"btag{tagger}UDG")
                and hasattr(events.Jet, f"btag{tagger}SvUDG")
                and hasattr(events.Jet, f"btag{tagger}CvL")
                and hasattr(events.Jet, f"btag{tagger}CvB")
            ):
                jets = events.Jet
                udg = jets[f"btag{tagger}UDG"]
                svudg = jets[f"btag{tagger}SvUDG"]
                cvl = jets[f"btag{tagger}CvL"]
                cvb = jets[f"btag{tagger}CvB"]

                jets[f"btag{tagger}S"] = ak.where(
                    svudg < 0.0,
                    -1,
                    svudg * udg / (1.0 - svudg),
                )
                jets[f"btag{tagger}C"] = ak.where(
                    svudg < 0.0,
                    -1,
                    cvl * (jets[f"btag{tagger}S"] + udg) / (1.0 - cvl),
                )
                jets[f"btag{tagger}B"] = ak.where(
                    svudg < 0.0,
                    -1,
                    (1.0 - cvb) * jets[f"btag{tagger}C"] / cvb,
                )
                jets[hfvlf_name] = ak.where(
                    svudg < 0.0,
                    -1,
                    (jets[f"btag{tagger}B"] + jets[f"btag{tagger}C"])
                    / (
                        jets[f"btag{tagger}B"]
                        + jets[f"btag{tagger}C"]
                        + jets[f"btag{tagger}S"]
                        + udg
                    ),
                )
                jets[f"{hfvlf_name}t"] = ak.where(
                    jets[hfvlf_name] < 0.0,
                    -1,
                    1.0 - np.sqrt(1.0 - jets[hfvlf_name]),
                )

                for lo, hi in ((25.0, 35.0), (35.0, 50.0), (50.0, 70.0), (70.0, 90.0), (90.0, 120.0)):
                    jets[f"{hfvlf_name}_pt{int(lo)}to{int(hi)}"] = ak.where(
                        (jets.pt < lo) | (jets.pt >= hi),
                        -2,
                        jets[hfvlf_name],
                    )
                jets[f"{hfvlf_name}_pt120toinf"] = ak.where(
                    jets.pt < 120.0,
                    -2,
                    jets[hfvlf_name],
                )

                update_map = {hfvlf_name: jets[hfvlf_name], f"{hfvlf_name}t": jets[f"{hfvlf_name}t"]}
                for lo, hi in ((25, 35), (35, 50), (50, 70), (70, 90), (90, 120)):
                    update_map[f"{hfvlf_name}_pt{lo}to{hi}"] = jets[f"{hfvlf_name}_pt{lo}to{hi}"]
                update_map[f"{hfvlf_name}_pt120toinf"] = jets[f"{hfvlf_name}_pt120toinf"]
                events.Jet = update(events.Jet, update_map)

        elif tagger in ("DeepFlav", "PNet"):
            if (
                not hasattr(events.Jet, hfvlf_name)
                and hasattr(events.Jet, f"btag{tagger}B")
                and hasattr(events.Jet, f"btag{tagger}CvL")
            ):
                jets = events.Jet
                b = jets[f"btag{tagger}B"]
                cvl = jets[f"btag{tagger}CvL"]
                jets[f"btag{tagger}C"] = ak.where(
                    b < 0.0,
                    -1,
                    cvl * (1 - b),
                )
                jets[hfvlf_name] = ak.where(
                    b < 0.0,
                    -1,
                    (b + jets[f"btag{tagger}C"]),
                )
                jets[f"{hfvlf_name}t"] = ak.where(
                    jets[hfvlf_name] < 0.0,
                    -1,
                    1.0 - np.sqrt(1.0 - jets[hfvlf_name]),
                )

                for lo, hi in ((25.0, 35.0), (35.0, 50.0), (50.0, 70.0), (70.0, 90.0), (90.0, 120.0)):
                    jets[f"{hfvlf_name}_pt{int(lo)}to{int(hi)}"] = ak.where(
                        (jets.pt < lo) | (jets.pt >= hi),
                        -2,
                        jets[hfvlf_name],
                    )
                jets[f"{hfvlf_name}_pt120toinf"] = ak.where(
                    jets.pt < 120.0,
                    -2,
                    jets[hfvlf_name],
                )

                update_map = {hfvlf_name: jets[hfvlf_name], f"{hfvlf_name}t": jets[f"{hfvlf_name}t"]}
                for lo, hi in ((25, 35), (35, 50), (50, 70), (70, 90), (90, 120)):
                    update_map[f"{hfvlf_name}_pt{lo}to{hi}"] = jets[f"{hfvlf_name}_pt{lo}to{hi}"]
                update_map[f"{hfvlf_name}_pt120toinf"] = jets[f"{hfvlf_name}_pt120toinf"]
                events.Jet = update(events.Jet, update_map)

        bin2d_name = f"btag{tagger}2Dbin"
        if (
            not hasattr(events.Jet, bin2d_name)
            and hasattr(events.Jet, hfvlf_name)
            and hasattr(events.Jet, bvc_name)
        ):
            jets = events.Jet
            hfvlf = ak.flatten(jets[hfvlf_name])
            bvc = ak.flatten(jets[bvc_name])
            nj = ak.num(jets)
            
            ihfvlf = np.digitize(hfvlf, btag_wp_dict[campaign][tagger]["2D"]["HFvLF"])
            ibvc = np.digitize(bvc, btag_wp_dict[campaign][tagger]["2D"]["BvC"])

            jets[bin2d_name] = ak.unflatten(
                [
                    (
                        -1
                        if hfvlf[i] == -1
                        else btag_wp_dict[campaign][tagger]["2D"]["mapping"][ihfvlf[i]][ibvc[i]]
                    )
                    for i in range(len(ihfvlf))
                ],
                nj,
            )

            for lo, hi in ((25.0, 35.0), (35.0, 50.0), (50.0, 70.0), (70.0, 90.0), (90.0, 120.0)):
                jets[f"{bin2d_name}_pt{int(lo)}to{int(hi)}"] = ak.where(
                    (jets.pt < lo) | (jets.pt >= hi),
                    -2,
                    jets[bin2d_name],
                )
            jets[f"{bin2d_name}_pt120toinf"] = ak.where(
                jets.pt < 120.0,
                -2,
                jets[bin2d_name],
            )

            update_map = {bin2d_name: jets[bin2d_name]}
            for lo, hi in ((25, 35), (35, 50), (50, 70), (70, 90), (90, 120)):
                update_map[f"{bin2d_name}_pt{lo}to{hi}"] = jets[f"{bin2d_name}_pt{lo}to{hi}"]
            update_map[f"{bin2d_name}_pt120toinf"] = jets[f"{bin2d_name}_pt120toinf"]
            events.Jet = update(events.Jet, update_map)

    if hasattr(events, "METFixEE2017"):
        events.MET = events.METFixEE2017
    if hasattr(events.PuppiMET, "ptUnclusteredUp") and not hasattr(
        events.PuppiMET, "MetUnclustEnUpDeltaX"
    ):
        met = events.PuppiMET
        met["MetUnclustEnUpDeltaX"] = (met.ptUnclusteredUp - met.pt) * np.cos(
            met.phiUnclusteredUp
        )
        met["MetUnclustEnUpDeltaY"] = (met.ptUnclusteredUp - met.pt) * np.sin(
            met.phiUnclusteredUp
        )
        events.PuppiMET = update(
            events.PuppiMET,
            {
                "MetUnclustEnUpDeltaX": met.MetUnclustEnUpDeltaX,
                "MetUnclustEnUpDeltaY": met.MetUnclustEnUpDeltaY,
            },
        )
    return events
