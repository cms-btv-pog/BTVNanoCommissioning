from functools import lru_cache

from BTVNanoCommissioning.helpers.func import uproot_writeable
from BTVNanoCommissioning.helpers.xsection import xsection
from BTVNanoCommissioning.helpers.xsection_13TeV import xsection_13TeV
import numpy as np
import awkward as ak
import os, uproot

# b/c-tag SF names weight_manager() (correction.py) may add via btagSFs(); only one
# is active at a time there, kept as a list so ftagWeight resolves to whichever is.
_FTAG_SF_NAMES = [
    "UParTAK4BC",
    "PNetBC",
    "DeepJetC",
    "DeepJetB",
    "DeepCSVB",
    "DeepCSVC",
]

# Central weight branches array_writer() always computes (see below). Listed
# explicitly here -- rather than relying on the "weight" substring match used for
# other branches -- since "genWeight" collides with a native NanoAOD field name and
# would otherwise be silently dropped by the setdiff1d-based default branch list.
_WEIGHT_BRANCHES = [
    "genWeight",
    "ttbar_weight",
    "puwei",
    "muWeight",
    "eleWeight",
    "ftagWeight",
]

arraySchema = {
    "CFM": _WEIGHT_BRANCHES
    + [
        "xSecWeight",
        "SelJet_btag",
        "SelJet_pt",
        "SelJet_muEF",
        "SelJet_chEmEF",
        "SelJet_chHEF",
        "SelJet_chMultiplicity",
        "SelJet_electronIdx1",
        "SelJet_hfEmEF",
        "SelJet_hfHEF",
        "SelJet_muonIdx1",
        "SelJet_nConstituents",
        "SelJet_nMuons",
        "SelJet_nElectrons",
        "SelJet_nSVs",
        "SelJet_neEmEF",
        "SelJet_neHEF",
        "SelJet_neMultiplicity",
        "SelJet_puIdDisc",
        "SelJet_rawFactor",
        "SelJet_eta",
        "SelJet_phi",
        "SelJet_mass",
        "SelJet_hadronFlavour",
        "SelJet_partonFlavour",
        "njet",
        "PuppiMET_pt",
        "PuppiMET_phi",
        "PV_npvs",
        "PV_npvsGood",
        "dilep_mass",
        "dilep_pt",
        "dilep_eta",
        "dilep_phi",
        "SoftMuon_dxySig",
        "MuonJet_muneuEF",
        "soft_l_ptratio",
        "osss",
        "W_transmass",
        "W_pt",
        "Pileup_nTrueInt",
        "Pileup_nPU",
    ],
    # Minimal per-jet content shared by DY/Wc/ttdilep: kinematics, flavour,
    # associated SV mass and multiplicity, the raw Jet-collection index (jetIdx,
    # meaningless/absent for Wc where SelJet is the single muon-tagged jet, not a
    # per-jet list), and only the inclusive HFvLF/BvC scores of the taggers built
    # for the 2D method (PNet, UParTAK4, NoPIDParTAK4) -- not their
    # pT-binned/2Dbin sub-variants.
    #
    # NB: the SV branches (SelJet_svmass/svd3dsig), SelJet_jetIdx and the
    # HFvLF/BvC scores are not produced on this branch yet -- they arrive with
    # the SV_link / 2D-tagger work. They are kept here verbatim so the schema
    # lights up unchanged once that lands; until then array_writer() warns
    # about each one it could not find rather than silently thinning the tree.
    "SLIM": _WEIGHT_BRANCHES
    + [
        "xSecWeight",
        "njet",
        "osss",  # OS(+1)/SS(-1) tag; needed for Wc to reproduce the OS-SS
        # subtraction applied to the histogram-based templates
        "SelJet_pt",
        "SelJet_eta",
        "SelJet_phi",
        "SelJet_mass",
        "SelJet_hadronFlavour",
        "SelJet_partonFlavour",
        "SelJet_svmass",
        "SelJet_nSVs",
        "SelJet_jetIdx",
        "SelJet_btagPNetHFvLF",
        "SelJet_btagPNetBvC",
        "SelJet_btagUParTAK4HFvLF",
        "SelJet_btagUParTAK4BvC",
        "SelJet_btagNoPIDParTAK4HFvLF",
        "SelJet_btagNoPIDParTAK4BvC",
    ],
}


@lru_cache(maxsize=1)
def _xs_lookup():
    """process_name -> cross section (times kFactor where one is given).

    Built once per process: array_writer() runs per shift per chunk, and the
    two cross-section tables together hold ~10^3 entries.
    """
    xs_dict = {}
    for obj in xsection + xsection_13TeV:
        xs = float(obj["cross_section"])
        if obj.get("kFactor"):
            xs *= float(obj["kFactor"])
        xs_dict[obj["process_name"]] = xs
    return xs_dict


def _flatten_to_perjet(ev, ref_collection):
    """Flatten a per-event branch dict (as produced by uproot_writeable) to one row
    per jet, using `ref_collection` (e.g. "SelJet") to determine jets/event.
    If `ref_collection` is already one row per event (not jagged), `ev` is returned
    unchanged since it's already effectively per-jet.

    uproot_writeable hands back nested collections as a single record-array entry
    (ev["SelJet"] with fields pt/eta/...), not pre-flattened "SelJet_pt" keys --
    uproot's mktree does that flattening itself on write. Flattening a jagged
    record array yields a flat record array, so the written branch names are
    unchanged either way.
    """
    if ref_collection in ev:
        ref = ev[ref_collection]
    else:
        prefix = f"{ref_collection}_"
        sub_keys = [k for k in ev if k.startswith(prefix)]
        if not sub_keys:
            raise ValueError(
                f"perJet=True requires '{ref_collection}' among the output branches"
            )
        ref = ev[sub_keys[0]]
    try:
        ref_counts = ak.num(ref, axis=1)
    except ValueError:  # numpy AxisError subclasses ValueError
        return ev
    # broadcast template: a bare per-jet field, so scalar branches can be
    # broadcast against the jet multiplicity without dragging the record along
    template = ref[ref.fields[0]] if ref.fields else ref

    out = {}
    for name, arr in ev.items():
        try:
            arr_counts = ak.num(arr, axis=1)
        except ValueError:  # numpy AxisError subclasses ValueError
            # scalar-per-event branch: repeat its value across each jet
            broadcasted, _ = ak.broadcast_arrays(arr, template)
            out[name] = ak.flatten(broadcasted, axis=1)
            continue
        if ak.all(arr_counts == ref_counts):
            out[name] = ak.flatten(arr, axis=1)
        # else: jagged with a per-event length unrelated to the jet
        # collection (e.g. LHEReweightingWeight's own weight-variation
        # count) -- cannot be represented in a one-row-per-jet tree, drop it
    return out


def array_writer(
    processor_class,  # the NanoProcessor class ("self")
    pruned_event,  # the event with specific calculated variables stored
    nano_event,  # entire NanoAOD/PFNano event with many variables
    weights,  # weight for the event
    systname,  # name of systematic shift
    dataset,  # dataset name
    isRealData,  # boolean
    out_dir_base="",  # string
    remove=[
        "SoftMuon",
        "MuonJet",
        "dilep",
        "OtherJets",
        "Jet",
    ],  # remove from variable list
    kinOnly=[
        "Muon",
        "Jet",
        "SoftMuon",
        "dilep",
        "charge",
        "MET",
    ],  # variables for which only kinematic properties are kept
    kins=[
        "pt",
        "eta",
        "phi",
        "mass",
        "pfRelIso04_all",
        "pfRelIso03_all",
        "dxy",
        "dz",
    ],  # kinematic propoerties for the above variables
    othersData=[
        "PFCands_*",
        "MuonJet_*",
        "SV_*",
        "PV_npvs",
        "PV_npvsGood",
        "Rho_*",
        "SoftMuon_dxySig",
        "Muon_sip3d",
    ],  # other fields, for Data and MC
    doOnly=None,
    schema=None,
    othersMC=["Pileup_nTrueInt", "Pileup_nPU"],  # other fields, for MC only
    empty=False,
    perJet=None,  # write one row per jet instead of one row per event
    perJet_collection="SelJet",  # jet collection used to determine jets/event
):
    # How much systematic content the trees carry, from runner.py --array-systs
    # (set on the processor instance; a processor constructed directly has no
    # such attribute and gets the same default as the CLI).
    #   none    nominal shift only, central weight components only
    #   weights adds one weight_syst_<variation> branch per weight variation
    #   shifts  also writes a tree per JES/JER shift, under its own directory
    #   both    both of the above  [default: what array_writer always did]
    # The two are independent axes and are deliberately NOT crossed: the
    # weight_* branches are written on the nominal shift only, so enabling
    # "both" gives shifted trees plus one set of varied weights, not their
    # product.
    syst_mode = getattr(processor_class, "arraySysts", "both")
    if systname[0] != "nominal" and syst_mode not in ("shifts", "both"):
        return

    # runner.py --array-perjet, same instance-attribute route as --array-systs.
    # An explicit perJet= from a caller always wins over the CLI default.
    if perJet is None:
        perJet = getattr(processor_class, "arrayPerJet", False)

    if weights is not None:
        pruned_event["weight"] = weights.weight()
        # One branch per weight component, named <component>_weight as on master.
        for ind_wei in weights.weightStatistics.keys():
            pruned_event[f"{ind_wei}_weight"] = weights.partial_weight(
                include=[ind_wei]
            )

        # Grouped central weights, in addition to the per-component branches
        # above, for the CFM/SLIM schemas.
        added = set(weights.weightStatistics.keys())

        def _partial(names):
            names = [n for n in names if n in added]
            if not names:
                return np.ones(len(pruned_event))
            return weights.partial_weight(include=names)

        pruned_event["genWeight"] = _partial(["genweight"])
        pruned_event["ttbar_weight"] = _partial(["ttbar_weight"])
        pruned_event["puwei"] = _partial(["puweight"])
        pruned_event["muWeight"] = _partial([n for n in added if n.startswith("mu_")])
        pruned_event["eleWeight"] = _partial([n for n in added if n.startswith("ele_")])
        pruned_event["ftagWeight"] = _partial(_FTAG_SF_NAMES)

        # Varied total weights, one weight_syst_<variation> branch each, as on
        # master. Only on the nominal shift: a JES-shifted tree's weight
        # variations are the same set, and writing them there would store the
        # cross product for no gain. The "weight" substring is what gets these
        # picked up by the out_branch selection below.
        if syst_mode in ("weights", "both") and systname[0] == "nominal":
            for variation in sorted(weights.variations):
                pruned_event[f"weight_syst_{variation}"] = weights.weight(
                    modifier=variation
                )

    if not isRealData:
        xs_dict = _xs_lookup()
        if dataset in xs_dict:
            pruned_event["xSecWeight"] = np.full(len(pruned_event), xs_dict[dataset])
        else:
            print(
                f"WARNING: '{dataset}' not found in xsection database, xSecWeight not written"
            )

    schema_requested = []
    if empty:
        print("WARNING: No events selected. Writing blank file.")
        out_branch = []
    elif doOnly is not None:
        if "weight" not in doOnly:
            doOnly.extend([b for b in pruned_event.fields if "weight" in b])
        out_branch = np.array(doOnly)
        # The named weight branches do not all contain the "weight" substring
        # ("genWeight" is capitalised, "puwei" is truncated), so append them
        # explicitly the same way the default branch list below does.
        out_branch = np.append(out_branch, _WEIGHT_BRANCHES)
        if not isRealData:
            out_branch = np.append(out_branch, othersMC)
            out_branch = np.append(out_branch, ["xSecWeight"])
    elif schema is not None:
        schema_requested = arraySchema[schema]
        netout = schema_requested + [b for b in pruned_event.fields if "weight" in b]
        out_branch = np.array(netout)
    else:
        # Get only the variables that were added newly
        out_branch = np.setdiff1d(
            np.array(pruned_event.fields), np.array(nano_event.fields)
        )

        # Handle kinOnly vars
        remove = remove + ["PFCands", "hl", "sl", "posl", "negl"]
        for v in remove:
            out_branch = np.delete(out_branch, np.where((out_branch == v)))

        for kin in kins:
            for obj in kinOnly:
                if "MET" in obj and ("pt" != kin or "phi" != kin):
                    continue
                if (obj != "SelMuon" and obj != "SoftMuon") and (
                    "pfRelIso04_all" == kin or "d" in kin
                ):
                    continue
                out_branch = np.append(out_branch, [f"{obj}_{kin}"])

        # Handle data vars
        out_branch = np.append(out_branch, othersData)
        # Weight branches: added explicitly since "genWeight" collides with a
        # native NanoAOD field name and would be missed by the setdiff1d above
        out_branch = np.append(out_branch, _WEIGHT_BRANCHES)

        if not isRealData:
            out_branch = np.append(out_branch, othersMC)
            out_branch = np.append(out_branch, ["xSecWeight"])

    # Write to root files
    outdir = f"{out_dir_base}{processor_class.name}/{systname[0]}/{dataset}/"
    os.makedirs(outdir, exist_ok=True)

    outfile = f"{outdir}/{nano_event.metadata['filename'].split('/')[-1].replace('.root','')}_{int(nano_event.metadata['entrystop']/processor_class.chunksize)}.root"

    with uproot.recreate(outfile) as fout:
        # uproot>=5.7 writes RNTuples for dict-like assignment, mktree keeps TTrees
        if not empty:
            events_out = uproot_writeable(pruned_event, include=out_branch)
            if perJet:
                events_out = _flatten_to_perjet(events_out, perJet_collection)
            fout.mktree("Events", events_out)
        fout.mktree(
            "TotalEventCount",
            ak.Array(
                [nano_event.metadata["entrystop"] - nano_event.metadata["entrystart"]]
            ),
        )
        if not isRealData:
            fout.mktree("TotalEventWeight", ak.Array([ak.sum(nano_event.genWeight)]))

    # Report what actually ended up in the tree rather than what was requested:
    # uproot_writeable silently drops requested names that do not exist on
    # pruned_ev, and silently adds extras (wildcard entries, bare collection
    # names, cross-reference "Idx"/"Flavor" fields).
    with uproot.open(outfile) as fin:
        written = sorted(fin["Events"].keys()) if "Events" in fin else []
    if schema_requested:
        # Only the schema path lists literal branch names, so it is the only one
        # where "requested but absent" is unambiguous (the other paths carry
        # wildcards and bare collection names).
        missing = [b for b in schema_requested if b not in written]
        if missing:
            print(
                f"WARNING: schema '{schema}' requested {len(missing)} branches absent "
                f"from pruned_ev:",
                missing,
            )
    print(
        f"Writing arrays to {os.path.abspath(outfile)} - branches ({len(written)}):",
        written,
    )
