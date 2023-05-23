import sys
from os import popen, listdir, makedirs, path, system
import json
import argparse

parser = argparse.ArgumentParser(
    description="Run analysis on baconbits files using processor coffea files"
)
parser.add_argument(
    "-i",
    "--input",
    default=None,
    type=str,
    required=True,
    help="List of samples in DAS (default: %(default)s)",
)
parser.add_argument(
    "-o",
    "--output",
    default=r"test_my_samples.json",
    help="Site (default: %(default)s)",
)
parser.add_argument(
    "--xrd",
    default="root://xrootd-cms.infn.it//",
    type=str,
    help="xrootd prefix string (default: %(default)s)",
)

parser.add_argument(
    "--from_path",
    action="store_true",
    help="For samples that are not published on DAS. If this option is set then the format of the --inpit file must be adjusted. It should be: \n dataset_name path_to_files.",
    default=False,
)

args = parser.parse_args()


def getFilesFromDas(args):
    fset = []
    with open(args.input) as fp:
        lines = fp.readlines()
        for line in lines:
            fset.append(line)

    fdict = {}

    for dataset in fset:
        if dataset.startswith("#") or dataset.strip() == "":
            # print("we skip this line:", line)
            continue

        dsname = dataset.strip().split("/")[1]  # Dataset first name

        Tier = dataset.strip().split("/")[
            3
        ]  # NANOAODSIM for regular samples, USER for private
        if "SIM" not in Tier:
            dsname = dataset.strip().split("/")[1] + "_" + dataset.split("/")[2]
        instance = "prod/global"
        if Tier == "USER":
            instance = "prod/phys03"
        print("Creating list of files for dataset", dsname, Tier, instance)
        flist = (
            popen(
                (
                    "/cvmfs/cms.cern.ch/common/dasgoclient -query='instance={} file dataset={}'"
                ).format(instance, fset[fset.index(dataset)].rstrip())
            )
            .read()
            .split("\n")
        )

        if dsname not in fdict:
            fdict[dsname] = [args.xrd + f for f in flist if len(f) > 1]
        else:  # needed to collect all data samples into one common key "Data" (using append() would introduce a new element for the key)
            fdict[dsname].extend([args.xrd + f for f in flist if len(f) > 1])

    # pprint.pprint(fdict, depth=1)
    return fdict


def getFilesFromPath(args, lim=None):
    fdict = {}
    fset = []
    with open(args.input) as fp:
        lines = fp.readlines()
        for line in lines:
            if line.startswith("#") or line.strip() == "":
                continue
            if line.startswith("/"):
                print(
                    "You are trying to read files from path, but providing a dataset in DAS:\n",
                    line,
                )
                print("That's not gonna work, so we exit here")
                sys.exit(1)
            ds = line.strip().split()
            print("ds=", ds)
            dataset = ds[0]
            fdict[ds[0]] = getRootFilesFromPath(ds[1])

    return fdict


def getRootFilesFromPath(d, lim=None):
    import subprocess

    if "xrootd" in d:
        sp = d.split("/")
        siteIP = "/".join(sp[0:4])
        pathToFiles = "/".join(sp[3:]) + "/"
        allfiles = str(
            subprocess.check_output(["xrdfs", siteIP, "ls", pathToFiles]), "utf-8"
        ).split("\n")
        # rootfiles = [siteIP+'/'+f for i,f in enumerate(allfiles) if f.endswith(".root") and (lim==None or i<lim)]
    else:
        siteIP = ""
        pathToFiles = d
        allfiles = [
            path.join(d, f) for i, f in enumerate(listdir(d)) if f.endswith(".root")
        ]

    # print(siteIP, pathToFiles)
    rootfiles = []
    for file_or_dir in allfiles:
        # print(file_or_dir)
        if file_or_dir == "" or file_or_dir == pathToFiles:
            continue
        file_or_dir = siteIP + file_or_dir
        if file_or_dir.endswith(".root"):
            if lim == None or len(rootfiles) < lim:
                rootfiles.append(file_or_dir)

        elif not "log" in file_or_dir and not file_or_dir[-1] == "/":
            file_or_dir = file_or_dir + "/"
            print("file or dir:", file_or_dir)
            if lim == None:
                rootfiles.extend(getRootFilesFromPath(file_or_dir))
            elif len(rootfiles) < lim:
                rootfiles.extend(
                    getRootFilesFromPath(file_or_dir, lim - len(rootfiles))
                )

    # print("Input path:", d)
    # print("List of root files to be processed:\n",rootfiles)

    return rootfiles


def main(args):
    if args.from_path:
        print("do it from path: ")

        fdict = getFilesFromPath(args)

    else:
        fdict = getFilesFromDas(args)

    # Check the any file lists empty
    empty = True
    for dsname, flist in fdict.items():
        if len(flist) == 0:
            print(dsname, "is empty!!!!")
            empty = False
    assert empty, "you have empty lists"
    output_file = "./%s" % (args.output)
    with open(output_file, "w") as fp:
        json.dump(fdict, fp, indent=4)
        print("The file is saved at: ", output_file)


if __name__ == "__main__":
    print("This is the __main__ part")
    main(args)
