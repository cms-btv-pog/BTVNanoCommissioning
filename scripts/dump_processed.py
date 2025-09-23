from coffea.util import load
import awkward as ak
import numpy as np
import sys, json, glob
import os


def dump_lumi(output, fname, year):
    lumi, run = [], []
    for m in output.keys():
        for f in output[m].keys():
            lumi.extend(output[m][f]["lumi"].value)
            run.extend(output[m][f]["run"].value)

    # Sort runs and keep lumisections matched
    run, lumi = np.array(run), np.array(lumi)
    sorted_indices = np.lexsort((lumi, run))  # Sort by run first, then lumi
    run = run[sorted_indices]
    lumi = lumi[sorted_indices]
    # Create dictionary with ls values for each run
    dicts = {}
    for r in np.unique(run):
        dicts[str(r)] = lumi[run == r]
    # Convert to format for brilcalc
    for r in dicts.keys():
        ar = ak.singletons(ak.Array(dicts[r]))
        ars = ak.concatenate([ar, ar], axis=-1)
        dicts[r] = ak.values_astype(ars, int).tolist()

    with open(f"{fname}_lumi.json", "w") as outfile:
        json.dump(dicts, outfile, indent=2)

    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/LumiRecommendationsRun3
    if year in ["2022", "2023"]:
        brilcalc_cmd = f"source /cvmfs/cms-bril.cern.ch/cms-lumi-pog/brilws-docker/brilws-env; eval 'brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -c web -i {fname}_lumi.json -u /pb '"
    elif year in ["2024", "2025"]:
        # Using recommended temporary Run 3 normtag
        brilcalc_cmd = f"source /cvmfs/cms-bril.cern.ch/cms-lumi-pog/brilws-docker/brilws-env; eval 'brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_BRIL.json -c web -i {fname}_lumi.json -u /pb --datatag online '"
    lumi_in_pb = os.popen(brilcalc_cmd).read()
    lumi_in_pb = lumi_in_pb[
        lumi_in_pb.find("#Summary:") : lumi_in_pb.find("#Check JSON:")
    ]
    lumi_in_pb = float(lumi_in_pb.split("\n")[4].split("|")[-2])

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
            original_list[o].extend(old[o])
    for m in output.keys():
        for f in output[m].keys():
            if f not in list_from_coffea.keys():
                list_from_coffea[f] = list(output[m][f]["fname"])
            else:
                list_from_coffea[f] += list(set(output[m][f]["fname"]))

    failed = {}
    for t in original_list.keys():
        failed[t] = []
        if t not in list_from_coffea.keys():
            failed[t] = original_list[t]
            continue
        for f in original_list[t]:
            if f not in list_from_coffea[t]:
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
        "-y",
        "--year",
        default="2024",
        choices=["2022", "2023", "2024", "2025"],
        help="The data-taking year to process the luminosity for",
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
        help="Original json files, splitted by ,. Wildcard option * available as well.",
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
        dump_lumi(output, args.fname, args.year)
    if args.type == "all" or args.type == "failed":
        print("===>Dump Failed Files")
        dump_dataset(output, args.fname, args.jsons)
