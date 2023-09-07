from coffea.util import load
import awkward as ak
import numpy as np
import sys, json, glob


def dump_lumi(output, fname):
    lumi, run = [], []
    for m in output.keys():
        for f in output[m].keys():
            lumi.extend(output[m][f]["lumi"].value)
            run.extend(output[m][f]["run"].value)

    lumi, run = np.array(lumi), np.array(run)
    dicts = {}
    for r in list(set(run)):
        dicts[str(r)] = lumi[r == run]
    for r in dicts.keys():
        ar = ak.singletons(ak.Array(dicts[r]))
        ars = ak.concatenate([ar, ar], axis=-1)
        dicts[r] = ak.values_astype(ars, int).tolist()

    with open(f"{fname}_lumi.json", "w") as outfile:
        json.dump(dicts, outfile, indent=2)


def dump_dataset(output, fname, alljson):
    jsonlist = glob.glob(alljson) if "*" in alljson else alljson.split(",")
    print("Original jsons:", jsonlist)
    oldf, newf = {}, {}
    for j in jsonlist:
        old = json.load(open(j))
        for o in old.keys():
            oldf[o] = old[o]
    for m in output.keys():
        for f in output[m].keys():
            newf[f] = list(output[m][f]["fname"])
    failed = {}
    for t in oldf.keys():
        failed[t] = list(set(oldf[t]) - set(newf[t]))

    with open(f"{fname}_failed_dataset.json", "w") as outfile:
        json.dump(failed, outfile, indent=4)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Dump procssed luminosity & failed coffea"
    )
    parser.add_argument(
        "-t",
        "--type",
        default="all",
        choices=["all", "lumi", "failed"],
        help="Choose the function for dump luminosity(`lumi`)/failed files(`failed`) into json",
    )
    parser.add_argument(
        "-c",
        "--coffea",
        required=True,
        help="Processed coffea files, splitted by ,. Wildcard option * available as well.",
    )
    parser.add_argument(
        "-n", "--fname", required=True, help="Output name of jsons(with _lumi/_dataset)"
    )
    parser.add_argument(
        "-j",
        "--jsons",
        type=str,
        help="Original json files, splitted by ,. Wildcard option * available as well. ",
    )
    args = parser.parse_args()
    if len(args.coffea.split(",")) > 1:
        output = {i: load(i) for i in args.coffea.split(",")}
    elif "*" in args.coffea:
        args.coffea = glob.glob(args.coffea)
        output = {i: load(i) for i in args.coffea}
    else:
        output = {args.coffea: load(args.coffea)}

    if args.type == "all" or args.type == "lumi":
        print("===>Dump Processed Luminosity")
        dump_lumi(output, args.fname)
    if args.type == "all" or args.type == "failed":
        print("===>Dump Failed Files")
        dump_dataset(output, args.fname, args.jsons)
