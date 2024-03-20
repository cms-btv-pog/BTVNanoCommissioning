from BTVNanoCommissioning.helpers.func import update
from BTVNanoCommissioning.utils.correction import add_jec_variables
import numpy as np


def missing_branch(events):
    events["fixedGridRhoFastjetAll"] = (
        events.fixedGridRhoFastjetAll
        if hasattr(events, "fixedGridRhoFastjetAll")
        else events.Rho.fixedGridRhoFastjetAll
    )
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
    if hasattr(events, "METFixEE2017"):
        events.MET = events.METFixEE2017
    if not hasattr(events.PuppiMET, "MetUnclustEnUpDeltaX"):
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


### Not used anymore
def add_jec(events, campaign, jmestuff):
    dataset = events.metadata["dataset"]
    jet_factory = jmestuff["jet_factory"]
    met_factory = jmestuff["met_factory"]
    isRealData = not hasattr(events, "genWeight")
    if isRealData:
        if "un" in dataset:
            jecname = dataset[dataset.find("un") + 6]
        else:
            print("No valid jec name")
            raise NameError
        if "2016_" in campaign:
            if "B" == jecname or "C" == jecname or "D" == jecname:
                jecname = "BCD"
            elif "E" == jecname or "F" == jecname:
                jecname = "EF"
            elif "F" == jecname or "G" == jecname or "H" == jecname:
                jecname = "FGH"
            else:
                raise NameError
        elif campaign == "Rereco17_94X":
            jecname = ""
        jets = jet_factory[f"data{jecname}"].build(
            add_jec_variables(events.Jet, events.fixedGridRhoFastjetAll),
            lazy_cache=events.caches[0],
        )
    else:
        jets = jet_factory["mc"].build(
            add_jec_variables(events.Jet, events.fixedGridRhoFastjetAll),
            lazy_cache=events.caches[0],
        )

    if "Run3" not in campaign:
        met = met_factory.build(events.MET, jets, {})
        events = update(events, {"Jet": jets, "MET": met})
    else:
        met = met_factory.build(events.PuppiMET, jets, {})
        events = update(events, {"Jet": jets, "PuppiMET": met})
    return events
