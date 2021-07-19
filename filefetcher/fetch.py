import os
import json
import argparse
import pprint

xrootd_pfx = {
    "Americas": "root://cmsxrootd.fnal.gov/",
    "Eurasia": "root://xrootd-cms.infn.it/",
    "Yolo": "root://cms-xrd-global.cern.ch/",
}


def get_parser():
    parser = argparse.ArgumentParser(
        description="Query dasgoclient for dataset file lists"
    )

    parser.add_argument(
        "-i",
        "--input",
        help="What input dataset definition file to process.",
        required=True,
    )
    parser.add_argument(
        "-w",
        "--where",
        help="Where are you running your jobs? (default: %(default)s)",
        default="Americas",
        choices=["Americas", "Eurasia", "Yolo"],
    )
    parser.add_argument(
        "-x",
        "--xrootd",
        help="Override xrootd prefix with the one given.",
        default=None,
    )
    parser.add_argument(
        "--dbs-instance",
        dest="instance",
        help="The DBS instance to use for querying datasets. (default: %(default)s)",
        type=str,
        default="prod/global",
        choices=["prod/global", "prod/phys01", "prod/phys02", "prod/phys03"],
    )

    return parser


if __name__ == "__main__":

    parser = get_parser()
    args = parser.parse_args()

    if ".txt" not in args.input:
        raise Exception("Input file must have '.txt' extension and be a text file!")

    fset = []
    with open(args.input) as fp:
        for i, line in enumerate(fp.readlines()):
            if line.strip().startswith("#"):
                continue
            fset.append(tuple(line.strip().split()))
            if len(fset[-1]) != 2:
                raise Exception(
                    f"Text file format should be '<short name> <dataset path>' and nothing else.\nInvalid spec on line {i+1}: '{line}'"
                )

    fdict = {}

    for name, dataset in fset:
        flist = (
            os.popen(
                ("dasgoclient -query='instance={} file dataset={}'").format(
                    args.instance, dataset
                )
            )
            .read()
            .split("\n")
        )
        if name not in fdict:
            fdict[name] = [args.xrootd + f for f in flist if len(f) > 1]
        else:  # needed to collect all data samples into one common key "Data" (using append() would introduce a new element for the key)
            fdict[name].extend([args.xrootd + f for f in flist if len(f) > 1])

    # pprint.pprint(fdict, depth=1)
    with open(args.input[: args.input.rfind(".txt")] + ".json", "w") as fp:
        json.dump(fdict, fp, indent=4)
