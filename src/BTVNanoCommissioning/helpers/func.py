import awkward as ak
import numpy as np
from coffea import processor
import psutil, os
import uproot


def memory_usage_psutil():
    # return the memory usage in MB
    process = psutil.Process(os.getpid())
    mem = psutil.Process(os.getpid()).memory_info().rss / 1024**2
    return mem


def flatten(ar):  # flatten awkward into a 1d array to hist
    return ak.flatten(ar, axis=None)


def normalize(val, cut):
    if cut is None:
        ar = ak.to_numpy(ak.fill_none(val, np.nan))
        return ar
    else:
        ar = ak.to_numpy(ak.fill_none(val[cut], np.nan))
        return ar


def update(events, collections):
    """Return a shallow copy of events array with some collections swapped out"""
    out = events
    for name, value in collections.items():
        out = ak.with_field(out, value, name)
    return out


# return run & lumiblock in pairs
def dump_lumi(events, output):
    pairs = np.vstack((events.run.to_numpy(), events.luminosityBlock.to_numpy()))
    # remove replicas
    pairs = np.unique(np.transpose(pairs), axis=0)
    pairs = pairs[
        np.lexsort(([pairs[:, i] for i in range(pairs.shape[1] - 1, -1, -1)]))
    ]
    output["fname"] = processor.set_accumulator([events.metadata["filename"]])
    output["run"] = processor.column_accumulator(pairs[:, 0])
    output["lumi"] = processor.column_accumulator(pairs[:, 1])
    return output


def num(ar):
    return ak.num(ak.fill_none(ar[~ak.is_none(ar)], 0), axis=0)


## Based on https://github.com/CoffeaTeam/coffea/discussions/735
def _is_rootcompat(a):
    """Is it a flat or 1-d jagged array?"""
    t = ak.type(a)
    if isinstance(t, ak._ext.ArrayType):
        if isinstance(t.type, ak._ext.PrimitiveType):
            return True
        if isinstance(t.type, ak._ext.ListType) and isinstance(
            t.type.type, ak._ext.PrimitiveType
        ):
            return True
    return False


def uproot_writeable(events, include=["events", "run", "luminosityBlock"]):
    ev = {}
    include = np.array(include)
    no_filter = False
    if len(include) == 1 and include[0] == "*":
        no_filter = False
    for bname in events.fields:
        if not events[bname].fields:
            if not no_filter and bname not in include:
                continue
            ev[bname] = ak.fill_none(
                ak.packed(ak.without_parameters(events[bname])), -99
            )
        else:
            b_nest = {}
            no_filter_nest = False
            if all(np.char.startswith(include, bname) == False):
                continue
            include_nest = [
                i[i.find(bname) + len(bname) + 1 :]
                for i in include
                if i.startswith(bname)
            ]

            if len(include_nest) == 1 and include_nest[0] == "*":
                no_filter_nest = True

            if not no_filter_nest:
                mask_wildcard = np.char.find(include_nest, "*") != -1
                include_nest = np.char.replace(include_nest, "*", "")
            for n in events[bname].fields:
                ## make selections to the filter case, keep cross-ref ("Idx")
                if (
                    not no_filter_nest
                    and all(np.char.find(n, include_nest) == -1)
                    and "Idx" not in n
                    and "Flavor" not in n
                ):
                    continue
                evnums = ak.num(events[bname][n], axis=0)
                if not isinstance(evnums, int):
                    continue
                if not _is_rootcompat(events[bname][n]) and evnums != len(
                    flatten(events[bname][n])
                ):
                    continue
                # skip IdxG
                if "IdxG" in n:
                    continue
                b_nest[n] = ak.fill_none(
                    ak.packed(ak.without_parameters(events[bname][n])), -99
                )
            ev[bname] = ak.zip(b_nest)
    return ev


def array_writer(
    processor_class,  # the NanoProcessor class ("self")
    pruned_event,  # the event with specific calculated variables stored
    nano_event,  # entire NanoAOD/PFNano event with many variables
    systname,  # name of systematic shift
    dataset,  # dataset name
    isRealData,  # boolean
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
    othersMC=["Pileup_nTrueInt", "Pileup_nPU"],  # other fields, for MC only
    empty=False
):
    if empty:
        print("WARNING: No events selected. Writing blank file.")
        out_branch = []
    else:
        # Get only the variables that were added newly
        out_branch = np.setdiff1d(
            np.array(pruned_event.fields), np.array(nano_event.fields)
        )

        # Handle kinOnly vars
        for v in remove:
            out_branch = np.delete(out_branch, np.where((out_branch == v)))

        for kin in kins:
            for obj in kinOnly:
                if "MET" in obj and ("pt" != kin or "phi" != kin):
                    continue
                if (obj != "Muon" and obj != "SoftMuon") and (
                    "pfRelIso04_all" == kin or "d" in kin
                ):
                    continue
                out_branch = np.append(out_branch, [f"{obj}_{kin}"])

        # Handle data vars
        out_branch = np.append(out_branch, othersData)

        if not isRealData:
            out_branch = np.append(out_branch, othersMC)

    # Write to root files
    print("Branches to write:", out_branch)
    outdir = f"{processor_class.name}/{systname}/{dataset}/"
    os.system(f"mkdir -p {outdir}")

    with uproot.recreate(
        f"{outdir}/{nano_event.metadata['filename'].split('/')[-1].replace('.root','')}_{int(nano_event.metadata['entrystop']/processor_class.chunksize)}.root"
    ) as fout:
        if not empty: fout["Events"] = uproot_writeable(pruned_event, include=out_branch)
        fout["TotalEventCount"] = ak.Array(
            [nano_event.metadata["entrystop"] - nano_event.metadata["entrystart"]]
        )
        if not isRealData:
            fout["TotalEventWeight"] = ak.Array([ak.sum(nano_event.genWeight)])
