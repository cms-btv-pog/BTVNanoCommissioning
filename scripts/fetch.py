import sys
import os
import json
import argparse
from collections import defaultdict
import uproot
import numpy as np
from BTVNanoCommissioning.workflows import workflows
from BTVNanoCommissioning.utils.sample import predefined_sample

# Adapt some developments from Andrey Pozdnyakov in CoffeaRunner https://github.com/cms-rwth/CoffeaRunner/blob/master/filefetcher/fetch.py
parser = argparse.ArgumentParser(
    description="Run analysis on baconbits files using processor coffea files"
)
parser.add_argument(
    "-i",
    "--input",
    default=None,
    type=str,
    # required=True,
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
    default=None,
    type=str,
    help="xrootd prefix string otherwise get from available sites",
)

parser.add_argument(
    "--from_path",
    action="store_true",
    help="For samples that are not published on DAS. If this option is set then the format of the --input file must be adjusted. It should be: \n dataset_name path_to_files.",
    default=False,
)
parser.add_argument(
    "--from_dataset",
    help="input dataset only",
    action="store_true",
    default=False,
)
parser.add_argument(
    "-wf",
    "--from_workflow",
    help="Use the predefined workflows",
    choices=list(workflows.keys()),
    default=None,
)
parser.add_argument(
    "--testfile",
    action="store_true",
    help="Construct file list in the test directory. Specify the test directory path, create the json file for individual dataset",
    default=False,
)
parser.add_argument(
    "--whitelist_sites",
    help="White list fot sites",
    default=None,
)
parser.add_argument(
    "--blacklist_sites",
    help="Black list for sites",
    default=None,
)
parser.add_argument(
    "--limit", help="Limit numbers of file to create json", default=None, type=int
)

parser.add_argument(
    "-r",
    "--redirector",
    help="xrootd ridirector in case sites are not found",
    choices=["infn", "fnal", "cern"],
    default="infn",
)
parser.add_argument(
    "-j", "--ncpus", help="Number of CPUs to use for validation", default="4"
)
parser.add_argument(
    "--skipvalidation",
    action="store_true",
    help="If true, the readability of files will not be validated.",
    default=False,
)
parser.add_argument(
    "--overwrite",
    action="store_true",
    help="Overwrite existing file?",
    default=False,
)

parser.add_argument(
    "--DAS_campaign",
    help="campaign info, specifying dataset name in DAS. If you are running with ```from_workflow`` option, please do ```data_camapgin,mc_campaign``` split by ,",
    default=None,
    type=str,
)
parser.add_argument(
    "-c",
    "--campaign",
    help="campaign name (same as the campaign in runner.py)",
    default=None,
    require=True,
    type=str,
)
parser.add_argument("--year", help="year", default=None, type=str)

args = parser.parse_args()


# Based on https://github.com/PocketCoffea/PocketCoffea/blob/main/pocket_coffea/utils/rucio.py
def get_xrootd_sites_map():
    sites_xrootd_access = defaultdict(dict)
    if not os.path.exists(".sites_map.json"):
        print("Loading SITECONF info")
        sites = [
            (s, "/cvmfs/cms.cern.ch/SITECONF/" + s + "/storage.json")
            for s in os.listdir("/cvmfs/cms.cern.ch/SITECONF/")
            if s.startswith("T")
        ]
        for site_name, conf in sites:
            if not os.path.exists(conf):
                continue
            try:
                data = json.load(open(conf))
            except:
                continue
            for site in data:
                if site["type"] != "DISK":
                    continue
                if site["rse"] == None:
                    continue
                for proc in site["protocols"]:
                    if proc["protocol"] == "XRootD":
                        if proc["access"] not in ["global-ro", "global-rw"]:
                            continue
                        if "prefix" not in proc:
                            if "rules" in proc:
                                for rule in proc["rules"]:
                                    sites_xrootd_access[site["rse"]][rule["lfn"]] = (
                                        rule["pfn"]
                                    )
                        else:
                            sites_xrootd_access[site["rse"]] = proc["prefix"]
        json.dump(sites_xrootd_access, open(".sites_map.json", "w"))

    return json.load(open(".sites_map.json"))


def getFilesFromDas(args):
    fset = []
    with open(args.input) as fp:
        lines = fp.readlines()
        for line in lines:
            fset.append(line)

    fdict = {}
    if args.blacklist_sites is not None:
        args.blacklist_sites = args.blacklist_sites.split(",")
        print("blacklist sites:", args.blacklist_sites)
    if args.whitelist_sites is not None:
        args.whitelist_sites = args.whitelist_sites.split(",")
        print("whitelist sites:", args.whitelist_sites)
    for dataset in fset:
        if dataset.startswith("#") or dataset.strip() == "":
            continue

        dsname = dataset.strip().split("/")[1]  # Dataset first name

        Tier = dataset.strip().split("/")[3]
        # NANOAODSIM for regular samples, USER for private
        if "Run" in dataset and "pythia8" not in dataset:
            dsname = dataset.strip().split("/")[1] + dataset.strip().split("/")[2]
            if "BTV" in dataset:
                dsname = (
                    dataset.strip().split("/")[1]
                    + dataset.strip().split("/")[2][
                        dataset.strip().split("/")[2].find("-") + 1 :
                    ]
                )
                dsname = dsname[: dsname.find("_BTV")]
        instance = "prod/global"
        if Tier == "USER":
            instance = "prod/phys03"
        print("Creating list of files for dataset", dsname, Tier, instance)
        flist = (
            os.popen(
                (
                    "/cvmfs/cms.cern.ch/common/dasgoclient -query='instance={} file dataset={}'"
                ).format(instance, fset[fset.index(dataset)].rstrip())
            )
            .read()
            .split("\n")
        )
        print("Number of files: ", len(flist))
        import json

        dataset = dataset[:-1] if "\n" in dataset else dataset
        fetchsite = json.loads(
            (
                os.popen(
                    f"/cvmfs/cms.cern.ch/common/dasgoclient -query='instance={instance}  dataset={dataset} site' -json"
                ).read()
            )
        )

        if instance == "prod/phys03":
            possible_sites = [fetchsite[0]["site"][0]["name"]]
        else:
            possible_sites = [
                s["site"][0]["name"]
                for s in fetchsite
                if (s["site"][0]["replica_fraction"] == "100.00%")
                & (s["site"][0]["block_completion"] == "100.00%")
                & (s["site"][0]["kind"] == "DISK")
            ]

        sites_xrootd_prefix = get_xrootd_sites_map()
        xrd = None

        if args.xrd == None:
            if args.whitelist_sites is not None:
                # get first site in whitelist_site
                for w_site in args.whitelist_sites:
                    if xrd is not None:
                        break
                    for site in possible_sites:
                        if site == w_site and type(sites_xrootd_prefix[site]) == str:
                            xrd = sites_xrootd_prefix[site]
            else:
                for site in possible_sites:
                    # skip sites without xrd prefix
                    if (
                        site == "T2_IT_Pisa"
                        or site == "T2_IT_Bari"
                        or site == "T2_BE_UCL"
                        or site == "T2_IT_Rome"
                        or site == "T2_FR_GRIF"
                    ):
                        continue
                    if (
                        args.blacklist_sites is not None
                        and site in args.blacklist_sites
                    ):
                        continue
                    # get first site in possible_sites
                    elif type(sites_xrootd_prefix[site]) == str:
                        xrd = sites_xrootd_prefix[site]
        else:
            xrd = args.xrd

        if xrd is None:
            print(
                f"No SITE available in the whitelist for file {dsname}, change to global redirector: {args.redirector}"
            )
            redirector = {
                "infn": "root://xrootd-cms.infn.it//",
                "fnal": "root://cmsxrootd.fnal.gov/",
                "cern": "root://cms-xrd-global.cern.ch/",
            }
            xrd = redirector[args.redirector]
        if args.limit is not None:
            flist = flist[: args.limit]
        if dsname not in fdict:
            fdict[dsname] = [xrd + f for f in flist if len(f) > 1]
        else:  # needed to collect all data samples into one common key "Data" (using append() would introduce a new element for the key)
            fdict[dsname].extend([xrd + f for f in flist if len(f) > 1])
    return fdict


def getFilesFromPath(args):
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
            fdict[ds[0]] = getRootFilesFromPath(ds[1], args.limit)

    return fdict


def getTestlist(args):
    fdict = {}
    with open(args.input) as fp:
        lines = fp.readlines()
        for line in lines:
            if line.startswith("#") or line.strip() == "":
                continue
            if not line.endswith("/"):
                line = line + "/"
            if "test" not in line:
                print("You are not getting files in test directory")

            dirs_in_test = os.popen(f"gfal-ls {line}").read().split("\n")
            for s in dirs_in_test:
                if s == "":
                    continue
                print("dataset: ", s)
                fdict[s] = getRootFilesFromPath(line + s, 1)
    return fdict


def getRootFilesFromPath(d, lim=None):
    import subprocess

    if "root://" in d:
        sp = d.split("/")
        siteIP = "/".join(sp[0:4])
        pathToFiles = "/".join(sp[3:]) + "/"
        allfiles = str(
            subprocess.check_output(["xrdfs", siteIP, "ls", pathToFiles]), "utf-8"
        ).split("\n")
    else:
        siteIP = ""
        pathToFiles = d
        allfiles = [
            os.path.join(d, f)
            for i, f in enumerate(os.listdir(d))
            if f.endswith(".root")
        ]

    rootfiles = []
    for file_or_dir in allfiles:
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

    return rootfiles


def validate(file):
    n_tries = 0
    check_path = os.popen(f"gfal-ls {file}").read()
    if check_path == "":
        return f"NotFound: {file}"
    not_able_open = True
    while n_tries <= 5:
        try:
            fin = uproot.open(file)
            not_able_open = False
            return fin["Events"].num_entries
        except:
            print("retries", n_tries, file)
            n_tries += n_tries
    if not_able_open:
        return f"FailedRetries: {file}"


def remove_bad_files(sample_dict, outname, remove_bad=True):
    from p_tqdm import p_map

    all_invalid = []
    bad_sample_dict = {}
    f = open(f"{outname.replace('.json','')}_file_status.txt", "w")
    for sample in sample_dict.keys():
        _rmap = p_map(
            validate,
            sample_dict[sample],
            num_cpus=int(args.ncpus),
            desc=f"Validating {sample[:20]}...",
        )

        _results = list(_rmap)
        counts = np.sum([r for r in _results if np.isreal(r) and isinstance(r, int)])
        all_invalid += [r for r in _results if isinstance(r, str)]
        if len(all_invalid) > 0:
            bad_sample_dict[sample] = all_invalid
        f.write(f"{sample} Events: {np.sum(counts)}\n")

    if len(all_invalid) == 0:
        f.write("No bad files found!")
    else:
        print(f"Found {len(all_invalid)} bad files.")
        f.write(f"Found {len(all_invalid)} bad files.")
        if remove_bad == True:
            f.write(
                "\n==========================BAD FILES==========================\n "
            )

            for sample in bad_sample_dict.keys():
                for bad_file in bad_sample_dict[sample]:
                    f.write(bad_file + "\n")
                    if bad_file[bad_file.find("root://") :] in sample_dict[sample]:
                        sample_dict[sample].remove(bad_file[bad_file.find("root://") :])

    return sample_dict


def main(args):

    if args.from_workflow:
        for sample in predefined_sample[args.from_workflow].keys():
            if (
                os.path.exists(
                    f"metadata/{args.campaign}/{sample}_{args.campaign}_{args.year}_{args.from_workflow}.json"
                )
                and args.overwrite == False
            ):
                raise Exception(
                    f"metadata/{args.campaign}/{sample}_{args.campaign}_{args.year}_{args.from_workflow}.json exists"
                )
    elif os.path.exists(args.output) and args.overwrite == False:
        raise Exception(f"{args.output} exists")

    ## If you only provide dataset from the dataset name(DAS) or do from_workflow
    if args.from_dataset or args.from_workflow is not None:
        if args.from_dataset:
            f = open(args.input)
            lines = f.readlines()
        else:
            lines = []
            for sample in predefined_sample[args.from_workflow].keys():
                lines += predefined_sample[args.from_workflow][sample]
            args.input = args.from_workflow + "_predef"
            data_campaign, mc_campaign = (
                args.DAS_campaign.split(",")[0],
                args.DAS_campaign.split(",")[1],
            )

        if args.DAS_campaign is None:
            raise ("Please provide the campaign info when input dataset")
        args.input = args.input + "_DAS_" + args.campaign
        outf = open(args.input, "w")

        for l in lines:
            l = l.replace("\n", "")
            # read campaigns from two inputs
            if args.from_workflow is not None:
                args.DAS_campaign = (
                    data_campaign
                    if l in predefined_sample[args.from_workflow]["data"]
                    else mc_campaign
                )

            dataset = (
                os.popen(
                    f"/cvmfs/cms.cern.ch/common/dasgoclient -query='instance=prod/global dataset=/{l}/*{args.DAS_campaign}*/NANOAOD*'"
                )
                .read()[:-1]
                .split("\n")
            )
            if dataset[0] == "":
                print(l, "not Found! List all campaigns")
                dataset = (
                    os.popen(
                        f"/cvmfs/cms.cern.ch/common/dasgoclient -query='instance=prod/global dataset=/{l}/**/NANOAOD*'"
                    )
                    .read()[:-1]
                    .split("\n")
                )
                if dataset[0] == "":
                    print(f"{l} is not a valid dataset")
                    continue
                campaigns = [
                    d.split("/")[2]
                    for d in dataset
                    if "CMSSW" not in d and "Tier0" not in d and "ECAL" not in d
                ]
                args.DAS_campaign = input(f"which campaign? \n {campaigns} \n")
                dataset = (
                    os.popen(
                        f"/cvmfs/cms.cern.ch/common/dasgoclient -query='instance=prod/global dataset=/{l}/*{args.DAS_campaign}*/NANOAOD*'"
                    )
                    .read()[:-1]
                    .split("\n")
                )

                outf.write(dataset[0] + "\n")
                # continue

            elif len(dataset) > 1:
                campaigns = [d.split("/")[2] for d in dataset]
                if args.from_workflow is None or dataset[0].endswith("SIM"):
                    args.DAS_campaign = input(
                        f"{l} is which campaign? \n {campaigns} \n"
                    )

                    dataset = (
                        os.popen(
                            f"/cvmfs/cms.cern.ch/common/dasgoclient -query='instance=prod/global dataset=/{l}/*{args.DAS_campaign}*/NANOAOD*'"
                        )
                        .read()[:-1]
                        .split("\n")
                    )

                    outf.write(dataset[0] + "\n")
                else:
                    for d in dataset:
                        outf.write(d + "\n")

            else:
                outf.write(dataset[0] + "\n")
        outf.close()
    ## If put the path
    if args.from_path:
        print("do it from path: ")
        fdict = getFilesFromPath(args)
    elif args.testfile:
        fdict = getTestlist(args)
    else:
        fdict = getFilesFromDas(args)
    # Check the any file lists empty
    empty = True
    for dsname, flist in fdict.items():
        if len(flist) == 0:
            print(dsname, "is empty!!!!")
            empty = False
    assert empty, "you have empty lists"
    ## Remove files if not exist
    if not args.skipvalidation:
        fdict = remove_bad_files(fdict, args.output, True)  # remove bad files

    ## Create JSONs
    # create according to workflow
    if args.from_workflow:
        os.system(f"mkdir -p metadata/{args.campaign}/")
        for sample in predefined_sample[args.from_workflow].keys():
            reduced_fdict = {}
            for dataset in fdict.keys():
                for s in predefined_sample[args.from_workflow][sample]:

                    if s in dataset:
                        reduced_fdict[dataset] = fdict[dataset]

            with open(
                f"metadata/{args.campaign}/{sample}_{args.campaign}_{args.year}_{args.from_workflow}.json",
                "w",
            ) as fp:

                json.dump(reduced_fdict, fp, indent=4)

    else:
        os.system(f"mkdir -p metadata/{args.campaign}/")
        with open(f"metadata/{args.campaign}/{args.output}", "w") as fp:
            json.dump(fdict, fp, indent=4)
            print(
                "The file is saved at: metadata/", {args.campaign}, "/", {args.output}
            )


if __name__ == "__main__":
    print("This is the __main__ part")
    main(args)
