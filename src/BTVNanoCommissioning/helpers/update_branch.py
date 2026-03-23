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
    events["fixedGridRhoFastjetAll"] = (
        events.fixedGridRhoFastjetAll
        if hasattr(events, "fixedGridRhoFastjetAll")
        else events.Rho.fixedGridRhoFastjetAll
    )

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
    if not hasattr(events.Jet, "btagUParTAK4BvC") and hasattr(
        events.Jet, "btagUParTAK4CvB"
    ):
        jets = events.Jet
        jets["btagUParTAK4BvC"] = ak.where(
            jets.btagUParTAK4CvB < 0.0,
            -1,
            1.0 - jets.btagUParTAK4CvB,
        )
        jets["btagUParTAK4BvCt"] = ak.where(
            jets.btagUParTAK4BvC < 0.0,
            -1,
            1.0 - np.sqrt(1.0 - jets.btagUParTAK4BvC),
        )
        jets["btagUParTAK4BvC_pt25to35"] = ak.where(
            (jets.pt < 25.0) | (jets.pt >= 35.0),
            -2,
            jets.btagUParTAK4BvC,
        )
        jets["btagUParTAK4BvC_pt35to50"] = ak.where(
            (jets.pt < 35.0) | (jets.pt >= 50.0),
            -2,
            jets.btagUParTAK4BvC,
        )
        jets["btagUParTAK4BvC_pt50to70"] = ak.where(
            (jets.pt < 50.0) | (jets.pt >= 70.0),
            -2,
            jets.btagUParTAK4BvC,
        )
        jets["btagUParTAK4BvC_pt70to90"] = ak.where(
            (jets.pt < 70.0) | (jets.pt >= 90.0),
            -2,
            jets.btagUParTAK4BvC,
        )
        jets["btagUParTAK4BvC_pt90to120"] = ak.where(
            (jets.pt < 90.0) | (jets.pt >= 120.0),
            -2,
            jets.btagUParTAK4BvC,
        )
        jets["btagUParTAK4BvC_pt120toinf"] = ak.where(
            jets.pt < 120.0,
            -2,
            jets.btagUParTAK4BvC,
        )
        events.Jet = update(
            events.Jet,
            {
                "btagUParTAK4BvC": jets.btagUParTAK4BvC,
                "btagUParTAK4BvCt": jets.btagUParTAK4BvCt,
                "btagUParTAK4BvC_pt25to35": jets.btagUParTAK4BvC_pt25to35,
                "btagUParTAK4BvC_pt35to50": jets.btagUParTAK4BvC_pt35to50,
                "btagUParTAK4BvC_pt50to70": jets.btagUParTAK4BvC_pt50to70,
                "btagUParTAK4BvC_pt70to90": jets.btagUParTAK4BvC_pt70to90,
                "btagUParTAK4BvC_pt90to120": jets.btagUParTAK4BvC_pt90to120,
                "btagUParTAK4BvC_pt120toinf": jets.btagUParTAK4BvC_pt120toinf,
            },
        )
    if (
        not hasattr(events.Jet, "btagUParTAK4HFvLF")
        and hasattr(events.Jet, "btagUParTAK4UDG")
        and hasattr(events.Jet, "btagUParTAK4SvUDG")
        and hasattr(events.Jet, "btagUParTAK4CvL")
        and hasattr(events.Jet, "btagUParTAK4CvB")
    ):
        jets = events.Jet
        jets["btagUParTAK4S"] = ak.where(
            jets.btagUParTAK4SvUDG < 0.0,
            -1,
            jets.btagUParTAK4SvUDG
            * jets.btagUParTAK4UDG
            / (1.0 - jets.btagUParTAK4SvUDG),
        )
        jets["btagUParTAK4C"] = ak.where(
            jets.btagUParTAK4SvUDG < 0.0,
            -1,
            jets.btagUParTAK4CvL
            * (jets.btagUParTAK4S + jets.btagUParTAK4UDG)
            / (1.0 - jets.btagUParTAK4CvL),
        )
        jets["btagUParTAK4B"] = ak.where(
            jets.btagUParTAK4SvUDG < 0.0,
            -1,
            (1.0 - jets.btagUParTAK4CvB) * jets.btagUParTAK4C / jets.btagUParTAK4CvB,
        )
        jets["btagUParTAK4HFvLF"] = ak.where(
            jets.btagUParTAK4SvUDG < 0.0,
            -1,
            (jets.btagUParTAK4B + jets.btagUParTAK4C)
            / (
                jets.btagUParTAK4B
                + jets.btagUParTAK4C
                + jets.btagUParTAK4S
                + jets.btagUParTAK4UDG
            ),
        )
        jets["btagUParTAK4HFvLFt"] = ak.where(
            jets.btagUParTAK4HFvLF < 0.0,
            -1,
            1.0 - np.sqrt(1.0 - jets.btagUParTAK4HFvLF),
        )
        jets["btagUParTAK4HFvLF_pt25to35"] = ak.where(
            (jets.pt < 25.0) | (jets.pt >= 35.0),
            -2,
            jets.btagUParTAK4HFvLF,
        )
        jets["btagUParTAK4HFvLF_pt35to50"] = ak.where(
            (jets.pt < 35.0) | (jets.pt >= 50.0),
            -2,
            jets.btagUParTAK4HFvLF,
        )
        jets["btagUParTAK4HFvLF_pt50to70"] = ak.where(
            (jets.pt < 50.0) | (jets.pt >= 70.0),
            -2,
            jets.btagUParTAK4HFvLF,
        )
        jets["btagUParTAK4HFvLF_pt70to90"] = ak.where(
            (jets.pt < 70.0) | (jets.pt >= 90.0),
            -2,
            jets.btagUParTAK4HFvLF,
        )
        jets["btagUParTAK4HFvLF_pt90to120"] = ak.where(
            (jets.pt < 90.0) | (jets.pt >= 120.0),
            -2,
            jets.btagUParTAK4HFvLF,
        )
        jets["btagUParTAK4HFvLF_pt120toinf"] = ak.where(
            jets.pt < 120.0,
            -2,
            jets.btagUParTAK4HFvLF,
        )
        events.Jet = update(
            events.Jet,
            {
                "btagUParTAK4HFvLF": jets.btagUParTAK4HFvLF,
                "btagUParTAK4HFvLFt": jets.btagUParTAK4HFvLFt,
                "btagUParTAK4HFvLF_pt25to35": jets.btagUParTAK4HFvLF_pt25to35,
                "btagUParTAK4HFvLF_pt35to50": jets.btagUParTAK4HFvLF_pt35to50,
                "btagUParTAK4HFvLF_pt50to70": jets.btagUParTAK4HFvLF_pt50to70,
                "btagUParTAK4HFvLF_pt70to90": jets.btagUParTAK4HFvLF_pt70to90,
                "btagUParTAK4HFvLF_pt90to120": jets.btagUParTAK4HFvLF_pt90to120,
                "btagUParTAK4HFvLF_pt120toinf": jets.btagUParTAK4HFvLF_pt120toinf,
            },
        )
    if (
        not hasattr(events.Jet, "btagUParTAK42Dbin")
        and hasattr(events.Jet, "btagUParTAK4HFvLF")
        and hasattr(events.Jet, "btagUParTAK4BvC")
        and campaign is not None
    ):
        jets = events.Jet
        jets_pt = ak.flatten(jets.pt)
        hfvlf = ak.flatten(jets.btagUParTAK4HFvLF)
        bvc = ak.flatten(jets.btagUParTAK4BvC)
        nj = ak.num(jets)
        ihfvlf = np.digitize(hfvlf, btag_wp_dict[campaign]["UParTAK4"]["2D"]["HFvLF"])
        ibvc = np.digitize(bvc, btag_wp_dict[campaign]["UParTAK4"]["2D"]["BvC"])
        jets["btagUParTAK42Dbin"] = ak.unflatten(
            [
                (
                    -1
                    if hfvlf[i] == -1
                    else btag_wp_dict[campaign]["UParTAK4"]["2D"]["mapping"][ihfvlf[i]][
                        ibvc[i]
                    ]
                )
                for i in range(len(ihfvlf))
            ],
            nj,
        )
        jets["btagUParTAK42Dbin_pt25to35"] = ak.where(
            (jets.pt < 25.0) | (jets.pt >= 35.0),
            -2,
            jets.btagUParTAK42Dbin,
        )
        jets["btagUParTAK42Dbin_pt35to50"] = ak.where(
            (jets.pt < 35.0) | (jets.pt >= 50.0),
            -2,
            jets.btagUParTAK42Dbin,
        )
        jets["btagUParTAK42Dbin_pt50to70"] = ak.where(
            (jets.pt < 50.0) | (jets.pt >= 70.0),
            -2,
            jets.btagUParTAK42Dbin,
        )
        jets["btagUParTAK42Dbin_pt70to90"] = ak.where(
            (jets.pt < 70.0) | (jets.pt >= 90.0),
            -2,
            jets.btagUParTAK42Dbin,
        )
        jets["btagUParTAK42Dbin_pt90to120"] = ak.where(
            (jets.pt < 90.0) | (jets.pt >= 120.0),
            -2,
            jets.btagUParTAK42Dbin,
        )
        jets["btagUParTAK42Dbin_pt120toinf"] = ak.where(
            jets.pt < 120.0,
            -2,
            jets.btagUParTAK42Dbin,
        )
        events.Jet = update(
            events.Jet,
            {
                "btagUParTAK42Dbin": jets.btagUParTAK42Dbin,
                "btagUParTAK42Dbin_pt25to35": jets.btagUParTAK42Dbin_pt25to35,
                "btagUParTAK42Dbin_pt35to50": jets.btagUParTAK42Dbin_pt35to50,
                "btagUParTAK42Dbin_pt50to70": jets.btagUParTAK42Dbin_pt50to70,
                "btagUParTAK42Dbin_pt70to90": jets.btagUParTAK42Dbin_pt70to90,
                "btagUParTAK42Dbin_pt90to120": jets.btagUParTAK42Dbin_pt90to120,
                "btagUParTAK42Dbin_pt120toinf": jets.btagUParTAK42Dbin_pt120toinf,
            },
        )
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
