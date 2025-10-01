import hist as Hist
import awkward as ak
import numpy as np


def _flavor_label(flav):
    absflavs = np.abs(flav)
    labels = ak.where(
        (absflavs == 1) | (absflavs == 2),
        0,
        ak.where(
            absflavs == 3,
            1,
            ak.where(
                absflavs == 4,
                2,
                ak.where(absflavs == 5, 3, ak.where(absflavs == 21, 4, 5)),
            ),
        ),
    )
    return labels


def get_histograms(axes, **kwargs):
    hists = {}

    is_dijet = kwargs.get("is_dijet", False)
    jet_fields = kwargs.get("jet_fields", [])

    # Taggers
    taggers = [
        # "btagDeepFlavB",
        # "btagDeepFlavCvB",
        # "btagDeepFlavCvL",
        "btagDeepFlavQG",
        # "btagPNetB",
        # "btagPNetCvB",
        # "btagPNetCvL",
        "btagPNetQvG",
        # "btagRobustParTAK4B",
        # "btagRobustParTAK4CvB",
        # "btagRobustParTAK4CvL",
        "btagRobustParTAK4QG",
        # "btagUParTAK4B",
        # "btagUParTAK4CvB",
        # "btagUParTAK4CvNotB",
        # "btagUParTAK4CvL",
        "btagUParTAK4QG",
        "btagUParTAK4QvG",
        "btagUParTAK4SvCB",
        "btagUParTAK4SvUDG",
    ]

    objs = ["Tag", "SelJet"]
    if is_dijet:
        objs.extend(["FwdJet", "CenJet", "RndJet"])

    for obj in objs:
        if obj == "Tag":
            obj_axes = [
                axes["syst"],
            ]
        else:
            obj_axes = [axes["syst"], axes["flav"]]

        for tagger in taggers:
            if tagger not in jet_fields:
                continue
            # hists[f"Obj{obj}_Var{tagger}"] = Hist.Hist(
            # *obj_axes,
            # Hist.axis.Regular(50, 0, 1, name=tagger, label=tagger),
            # storage=Hist.storage.Weight(),
            # )
            hists[f"Obj{obj}_Var{tagger}_pteta"] = Hist.Hist(
                *obj_axes,
                Hist.axis.Regular(50, 0, 1, name=tagger, label=tagger),
                axes["pt"],
                axes["eta"],
                storage=Hist.storage.Weight(),
            )

        hists[f"Obj{obj}_Varpt"] = Hist.Hist(
            *obj_axes,
            axes["pt"],
            storage=Hist.storage.Weight(),
        )
        hists[f"Obj{obj}_Vareta"] = Hist.Hist(
            *obj_axes,
            axes["eta"],
            storage=Hist.storage.Weight(),
        )
        hists[f"Obj{obj}_Varphi"] = Hist.Hist(
            *obj_axes,
            axes["phi"],
            storage=Hist.storage.Weight(),
        )
        hists[f"Obj{obj}_Varmass"] = Hist.Hist(
            *obj_axes,
            axes["mass"],
            storage=Hist.storage.Weight(),
        )

    return hists


def qg_writer(
    events,
    output,
    weights,
    systematics: list,
    isSyst: bool,
    SF_map: dict,
):
    for syst in systematics:
        if not isSyst and syst != "nominal":
            break
        weight = (
            weights.weight()
            if syst == "nominal" or syst not in list(weights.variations)
            else weights.weight(modifier=syst)
        )
        # weight = weight * weights.partial_weight(include=["psweight"])

        for histname, hist in output.items():
            if histname in ["sumw", "fname", "run", "lumi", "processed", "out"]:
                continue
            hobj = histname.split("_Var")[0].replace("Obj", "")
            var = histname.split("_Var")[1].split("_")[0]
            is_pteta = histname.endswith("_pteta")
            if hobj not in events.fields:
                continue
            if var not in events[hobj].fields:
                continue

            obj_axes = {
                "syst": syst,
                var: ak.flatten(events[hobj][var], axis=None),
                # "weight": weight,
            }
            if is_pteta:
                obj_axes["pt"] = ak.flatten(events[hobj]["pt"], axis=None)
                obj_axes["eta"] = ak.flatten(events[hobj]["eta"], axis=None)

            if hobj != "Tag":
                if "partonFlavour" not in events[hobj].fields:
                    obj_axes["flav"] = ak.zeros_like(
                        ak.flatten(events[hobj].pt, axis=None), dtype=int
                    )
                else:
                    obj_axes["flav"] = ak.flatten(
                        _flavor_label(events[hobj].partonFlavour), axis=None
                    )

            w = ak.flatten(ak.broadcast_arrays(weight, events[hobj][var])[0], axis=None)
            obj_axes["weight"] = w

            output[histname].fill(**obj_axes)

    return output
