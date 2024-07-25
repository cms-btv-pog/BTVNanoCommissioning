from BTVNanoCommissioning.helpers.func import update
from BTVNanoCommissioning.utils.correction import add_jec_variables
import numpy as np


def missing_branch(events):
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
    if hasattr(events.Jet, "btagDeepFlavCvL") and not hasattr(
        events.Jet, "btagDeepFlavC"
    ):
        jets = events.Jet
        jets["btagDeepFlavC"] = (
            events.Jet.btagDeepFlavCvL / (1.0 - events.Jet.btagDeepFlavCvL)
        ) * events.Jet.btagDeepFlavB
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
    if not hasattr(events.Jet, "btagPNetCvNotB"):
        jets = events.Jet
        jets["btagPNetCvNotB"] = (
            jets.btagPNetCvB * jets.btagPNetB / (1.0 - jets.btagPNetB) ** 2
        )
        events.Jet = update(
            events.Jet,
            {"btagPNetCvNotB": jets.btagPNetCvNotB},
        )
    if not hasattr(events.Jet, "btagRobustParTAK4CvNotB"):
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
    if hasattr(events, "METFixEE2017"):
        events.MET = events.METFixEE2017
    if hasattr(events.PuppiMET, "ptUnclusteredUp") and not hasattr(
        events.PuppiMET, "MetUnclustEnUpDeltaX"
    ):
        met = events.PuppiMET
        met["MetUnclustEnUpDeltaX"] = met.ptUnclusteredUp * np.cos(met.phiUnclusteredUp)
        met["MetUnclustEnUpDeltaY"] = met.ptUnclusteredUp * np.sin(met.phiUnclusteredUp)
        events.PuppiMET = update(
            events.PuppiMET,
            {
                "MetUnclustEnUpDeltaX": met.MetUnclustEnUpDeltaX,
                "MetUnclustEnUpDeltaY": met.MetUnclustEnUpDeltaY,
            },
        )
    return events
