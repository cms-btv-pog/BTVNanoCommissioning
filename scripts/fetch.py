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
    required=True,
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

def run_das_command(cmd):
    """Run a DAS command with proper environment in micromamba"""
    import os
    import subprocess
    
    # Check if we're in GitLab CI
    in_ci = 'CI' in os.environ or 'GITLAB_CI' in os.environ
    
    if in_ci:
        # Set up full environment for GitLab CI with micromamba
        full_cmd = f"""
        source /cvmfs/cms.cern.ch/cmsset_default.sh &&
        which dasgoclient || export PATH=$PATH:/cvmfs/cms.cern.ch/common &&
        {cmd}
        """
        
        # Use subprocess to run command and capture output
        try:
            result = subprocess.run(['bash', '-c', full_cmd], 
                                    capture_output=True, 
                                    text=True)
            if result.returncode != 0:
                print(f"Command failed: {result.stderr}")
                return []
            return [line for line in result.stdout.strip().split('\n') if line]
        except Exception as e:
            print(f"Error executing command: {e}")
            return []
    else:
        # Normal execution for local environments
        return os.popen(cmd).read().splitlines()

def getFilesFromDas(args):
    """Improved getFilesFromDas with multiple fallback strategies"""
    fset = []
    with open(args.input) as fp:
        lines = fp.readlines()
        for line in lines:
            fset.append(line)

    fdict = {}
    # Site handling setup
    if args.blacklist_sites is not None:
        args.blacklist_sites = args.blacklist_sites.split(",")
        print("blacklist sites:", args.blacklist_sites)
    if args.whitelist_sites is not None:
        args.whitelist_sites = args.whitelist_sites.split(",")
        print("whitelist sites:", args.whitelist_sites)
        
    # Track failed datasets for summary report
    failed_datasets = []
    
    for dataset in fset:
        if dataset.startswith("#") or dataset.strip() == "":
            continue

        dataset = dataset.strip()
        print(f"\n===== Processing dataset: '{dataset}' =====")
        
        # Validate dataset format
        if not dataset.startswith('/'):
            print(f"WARNING: Dataset '{dataset}' does not start with '/' - trying anyway")
        
        # Try to get dsname safely
        try:
            dsname = dataset.split('/')[1]  # Primary dataset name
        except IndexError:
            print(f"ERROR: Cannot parse dataset name from '{dataset}'")
            print("Using the entire string as the dataset name")
            dsname = dataset.replace('/', '_').strip()
        
        # Try multiple DAS query approaches
        flist = []
        das_queries = [
            f"dasgoclient -query=\"file dataset={dataset}\"",
            f"dasgoclient -query=\"file dataset={dataset} instance=prod/global\"",
            f"dasgoclient -query=\"file dataset={dataset} instance=prod/phys03\"",
            f"dasgoclient -query=\"file dataset=*{dataset}*\""
        ]
        
        for query in das_queries:
            if flist:  # If we already have files, no need to try other queries
                break
                
            print(f"Trying DAS query: {query}")
            try:
                flist = run_das_command(query)
                if flist:
                    print(f"Found {len(flist)} files with query: {query}")
                    break
            except Exception as e:
                print(f"Error with DAS query: {e}")
        
        if not flist:
            print(f"WARNING: No files found for dataset '{dataset}' after trying all queries")
            failed_datasets.append(dataset)
            continue
        
        # Get sites with multiple fallback strategies
        sites = []
        site_queries = [
            f"dasgoclient -query=\"site dataset={dataset}\"",
            f"dasgoclient -query=\"site dataset={dataset} instance=prod/global\"",
            f"dasgoclient -query=\"site dataset={dataset} instance=prod/phys03\""
        ]
        
        for query in site_queries:
            if sites:  # If we already have sites, no need to try other queries
                break
                
            print(f"Trying site query: {query}")
            try:
                sites = run_das_command(query)
                if sites:
                    print(f"Found {len(sites)} sites with query: {query}")
                    break
            except Exception as e:
                print(f"Error with site query: {e}")
        
        # Fallback to default redirector if no sites found
        if not sites:
            print(f"WARNING: No sites found for dataset '{dataset}', using global redirector")
            redirector = {
                "infn": "root://xrootd-cms.infn.it//",
                "fnal": "root://cmsxrootd.fnal.gov/",
                "cern": "root://cms-xrd-global.cern.ch/"
            }
            xrd = redirector[args.redirector]
            if args.limit is not None:
                flist = flist[:args.limit]
            if dsname not in fdict:
                fdict[dsname] = [xrd + f for f in flist if len(f) > 1]
            else:
                fdict[dsname].extend([xrd + f for f in flist if len(f) > 1])
            continue
        
        # Process sites with careful error handling
        xrd = None
        sites_xrootd_prefix = get_xrootd_sites_map()
        
        for site in sites:
            if not site:
                continue
                
            # Handle site blacklisting/whitelisting
            if args.blacklist_sites is not None and site in args.blacklist_sites:
                print(f"Site {site} is blacklisted, skipping")
                continue
            if args.whitelist_sites is not None and site not in args.whitelist_sites:
                print(f"Site {site} is not in whitelist, skipping")
                continue
            
            # Safely handle site redirector lookup
            try:
                if site in sites_xrootd_prefix:
                    if isinstance(sites_xrootd_prefix[site], list):
                        xrd = sites_xrootd_prefix[site][0]
                    elif isinstance(sites_xrootd_prefix[site], str):
                        xrd = sites_xrootd_prefix[site]
                    
                    if xrd:
                        print(f"Using redirector for site {site}: {xrd}")
                        break
            except Exception as e:
                print(f"Error processing site {site}: {e}")
        
        # Fallback to global redirector if no valid site found
        if xrd is None:
            print(f"No valid site found for {dsname}, using global redirector: {args.redirector}")
            redirector = {
                "infn": "root://xrootd-cms.infn.it//",
                "fnal": "root://cmsxrootd.fnal.gov/",
                "cern": "root://cms-xrd-global.cern.ch/"
            }
            xrd = redirector[args.redirector]
        
        if args.limit is not None:
            flist = flist[:args.limit]
        
        # Add files to dictionary
        if dsname not in fdict:
            fdict[dsname] = [xrd + f for f in flist if len(f) > 1]
        else:
            fdict[dsname].extend([xrd + f for f in flist if len(f) > 1])
        
        print(f"Added {len(fdict[dsname])} files for dataset {dsname}")
    
    # Report on failures
    if failed_datasets:
        print("\n===== SUMMARY OF FAILED DATASETS =====")
        for ds in failed_datasets:
            print(f"- {ds}")
        print(f"Total: {len(failed_datasets)} failed datasets out of {len(fset)}")
    
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
