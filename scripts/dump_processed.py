from coffea.util import load
import awkward as ak
import numpy as np
import sys, json, glob
import os


def dump_lumi(output, fname):
    lumi, run = [], []
    for m in output.keys():
        for f in output[m].keys():
            lumi.extend(output[m][f]["lumi"].value)
            run.extend(output[m][f]["run"].value)

    lumi, run = np.array(sorted(lumi)), np.array(sorted(run))
    dicts = {}
    for r in list(set(run)):
        dicts[str(r)] = lumi[r == run]
    for r in dicts.keys():
        ar = ak.singletons(ak.Array(dicts[r]))
        ars = ak.concatenate([ar, ar], axis=-1)
        dicts[r] = ak.values_astype(ars, int).tolist()

    with open(f"{fname}_lumi.json", "w") as outfile:
        json.dump(dicts, outfile, indent=2)

    lumi_in_pb = os.popen(
        f"export PATH=$HOME/.local/bin:/cvmfs/cms-bril.cern.ch/brilconda3/bin:$PATH; brilcalc lumi -c web -i {fname}_lumi.json -u /pb "
    ).read()
    lumi_in_pb = lumi_in_pb[
        lumi_in_pb.find("#Summary:") : lumi_in_pb.find("#Check JSON:")
    ]
    lumi_in_pb = float(lumi_in_pb.split("\n")[-3].split("|")[-2])

    print(f"Luminosity in pb: {lumi_in_pb}")


def dump_dataset(output, fname, alljson):
    jsonlist = glob.glob(alljson) if "*" in alljson else alljson.split(",")
    print("Original jsons:", jsonlist)
    original_list, list_from_coffea = {}, {}
    for j in jsonlist:
        old = json.load(open(j))
        for o in old.keys():
            if o not in original_list.keys():
                original_list[o] = []
            original_list[o].append(old[o])

    for m in output.keys():
        for f in output[m].keys():
            if f not in list_from_coffea.keys():
                list_from_coffea[f] = []
            else:
                list_from_coffea[f] += list(set(output[m][f]["fname"]))
    failed = {}
    for t in original_list.keys():
        failed[t] = []
        if t not in list_from_coffea.keys():
            failed[t] = original_list[t]
            continue
        for f in original_list[t]:
            if not f in list_from_coffea[t]:
                failed[t].append(f)

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
