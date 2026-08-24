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
    is_trijet = kwargs.get("is_trijet", False)
    is_photondijet = kwargs.get("is_photondijet", False)
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
    if is_trijet:
        # SoftJet (j3) is the gluon-enriched probe of arXiv:1104.1175
        objs.extend(["SoftJet", "LeadJet", "SubleadJet"])
    if is_photondijet:
        # SoftJet (j2) is the quark-enriched probe of arXiv:1104.1175
        objs.extend(["SoftJet", "LeadJet"])

    for obj in objs:
        if obj == "Tag":
            obj_axes = [
                axes["syst"],
            ]
        else:
            obj_axes = [axes["syst"], axes["pflav"]]

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
                Hist.axis.Regular(
                    1,
                    0,
                    1,
                    name=tagger,
                    label=tagger,
                    underflow=False,
                    overflow=False,
                ),
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

    def _event_hists(prefix, spec):
        # Event-level kinematics, split by the parton flavour of the probe jet
        # so the quark/gluon purity can be read off directly.
        for name, ax in spec.items():
            hists[f"{prefix}_{name}"] = Hist.Hist(
                axes["syst"],
                axes["pflav"],
                axes[ax],
                storage=Hist.storage.Weight(),
            )
        hists[f"{prefix}_njet"] = Hist.Hist(
            axes["syst"],
            axes["n"],
            storage=Hist.storage.Weight(),
        )

    if is_trijet:
        _event_hists(
            "trijet",
            {
                "disc": "qgdisc",
                "deta12": "deta",
                "abseta3": "abseta",
                "mass": "trijetmass",
            },
        )
    if is_photondijet:
        _event_hists(
            "photondijet",
            {
                "disc": "gammadisc",
                "etaprod": "etaprod",
                "drgj2": "dr",
                "mass": "trijetmass",
            },
        )

    return hists


def _event_level_writer(
    events,
    output,
    weights,
    systematics: list,
    isSyst: bool,
    prefix: str,
    variables: dict,
    probe: str = "SoftJet",
):
    """Fill the event-level histograms of a quark/gluon workflow.

    The per-object histograms are handled by `qg_writer`; this only covers the
    cross-object quantities the workflow attached to `pruned_ev`. Everything is
    split by the parton flavour of `probe`, the jet the selection is meant to
    purify.

    variables: {histogram suffix: axis name} - the value is read from the
    `f"{prefix}_{suffix}"` field of `events`.
    """
    if f"{prefix}_njet" not in output:
        return output

    if "partonFlavour" in events[probe].fields:
        flav = _flavor_label(events[probe].partonFlavour)
    else:
        flav = ak.zeros_like(events[probe].pt, dtype=int)
    flav = ak.to_numpy(flav)

    filled = {
        f"{prefix}_{name}": (
            axis,
            ak.to_numpy(ak.fill_none(events[f"{prefix}_{name}"], np.nan)),
        )
        for name, axis in variables.items()
        if f"{prefix}_{name}" in output
    }
    njet = ak.to_numpy(events["njet"])

    for syst in systematics:
        if not isSyst and syst != "nominal":
            break
        weight = (
            weights.weight()
            if syst == "nominal" or syst not in list(weights.variations)
            else weights.weight(modifier=syst)
        )
        for name, (axis, value) in filled.items():
            output[name].fill(syst=syst, flav=flav, **{axis: value}, weight=weight)
        output[f"{prefix}_njet"].fill(syst=syst, n=njet, weight=weight)

    return output


def trijet_writer(events, output, weights, systematics: list, isSyst: bool):
    """Event-level histograms of the gluon-enriched trijet selection."""
    return _event_level_writer(
        events,
        output,
        weights,
        systematics,
        isSyst,
        "trijet",
        {"disc": "disc", "deta12": "deta", "abseta3": "abseta", "mass": "mass"},
    )


def photondijet_writer(events, output, weights, systematics: list, isSyst: bool):
    """Event-level histograms of the quark-enriched photon+2jet selection."""
    return _event_level_writer(
        events,
        output,
        weights,
        systematics,
        isSyst,
        "photondijet",
        {"disc": "disc", "etaprod": "etaprod", "drgj2": "dr", "mass": "mass"},
    )


def qg_writer(
    events,
    output,
    weights,
    systematics: list,
    isSyst: bool,
    SF_map: dict,
):
    # The axis values read off `events` do not depend on the systematic, and
    # within one object they are shared by every tagger. Flatten each of them
    # once up front rather than once per (systematic, histogram): that repeated
    # flattening dominated the runtime of systematics runs.
    flat = {}  # (obj, field) -> flattened field
    abs_eta = {}  # obj -> flattened |eta|
    flav = {}  # obj -> flattened flavour label
    template = {}  # obj -> array carrying the object's jagged shape
    hist_specs = []

    def flatten_field(hobj, field):
        if (hobj, field) not in flat:
            flat[(hobj, field)] = ak.flatten(events[hobj][field], axis=None)
        return flat[(hobj, field)]

    for histname in output:
        if "Var" not in histname or "Obj" not in histname:
            continue
        hobj = histname.split("_Var")[0].replace("Obj", "")
        var = histname.split("_Var")[1].split("_")[0]
        is_pteta = histname.endswith("_pteta")
        if hobj not in events.fields:
            continue
        if var not in events[hobj].fields:
            continue

        flatten_field(hobj, var)
        template.setdefault(hobj, events[hobj].pt)
        if is_pteta:
            flatten_field(hobj, "pt")
            if hobj not in abs_eta:
                abs_eta[hobj] = ak.flatten(np.abs(events[hobj]["eta"]), axis=None)
        if hobj != "Tag" and hobj not in flav:
            if "partonFlavour" not in events[hobj].fields:
                flav[hobj] = ak.zeros_like(flatten_field(hobj, "pt"), dtype=int)
            else:
                flav[hobj] = ak.flatten(
                    _flavor_label(events[hobj].partonFlavour), axis=None
                )

        hist_specs.append((histname, hobj, var, is_pteta))

    for syst in systematics:
        if not isSyst and syst != "nominal":
            break
        weight = (
            weights.weight()
            if syst == "nominal" or syst not in list(weights.variations)
            else weights.weight(modifier=syst)
        )
        # weight = weight * weights.partial_weight(include=["psweight"])

        # One broadcast per object instead of one per histogram: every tagger of
        # a given object shares that object's shape.
        flat_weight = {
            hobj: ak.flatten(ak.broadcast_arrays(weight, tmpl)[0], axis=None)
            for hobj, tmpl in template.items()
        }

        for histname, hobj, var, is_pteta in hist_specs:
            obj_axes = {
                "syst": syst,
                var: flat[(hobj, var)],
            }
            if is_pteta:
                obj_axes["pt"] = flat[(hobj, "pt")]
                obj_axes["eta"] = abs_eta[hobj]

            if hobj != "Tag":
                obj_axes["flav"] = flav[hobj]

            obj_axes["weight"] = flat_weight[hobj]

            output[histname].fill(**obj_axes)

    return output
