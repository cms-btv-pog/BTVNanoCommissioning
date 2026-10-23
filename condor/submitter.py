import os, re, sys
import json
import shutil
import subprocess
import tarfile
import argparse


def framework_toplevel(source_dir):
    """Top-level entries that make up the framework, i.e. the ones git tracks.

    Local job output (jobs_*, p3_*, arrays_*, the campaign output directories)
    lives beside the code and runs to hundreds of MB; it must not end up in the
    tarball that is transferred to every worker node. Returns None when this is
    not a git checkout, in which case the caller ships everything as before.
    """
    try:
        entries = subprocess.run(
            ["git", "-C", source_dir, "ls-tree", "--name-only", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
            timeout=60,
        ).stdout.split()
    except Exception:
        return None
    entries = [e for e in entries if os.path.exists(os.path.join(source_dir, e))]
    return entries or None


def make_tarfile(output_filename, source_dir, exclude_dirs=[]):
    toplevel = framework_toplevel(source_dir)
    if toplevel is None:
        print("  Not a git checkout: tarring the whole directory.")
        toplevel = sorted(os.listdir(source_dir))
    with tarfile.open(output_filename, "w:gz") as tar:
        for entry in toplevel:
            if entry in exclude_dirs:
                continue
            entry_path = os.path.join(source_dir, entry)
            if os.path.isfile(entry_path):
                tar.add(entry_path, arcname=entry)
                continue
            for root, dirs, files in os.walk(entry_path):
                dirs[:] = [
                    d for d in dirs if d not in exclude_dirs and d != "__pycache__"
                ]  # Exclude specified directories
                for file in files:
                    file_path = os.path.join(root, file)
                    tar.add(file_path, arcname=os.path.relpath(file_path, source_dir))


def get_condor_submitter_parser(parser):
    parser.add_argument(
        "--jobName",
        help="Condor job name to make the job directory",
        required=True,
    )
    parser.add_argument(
        "-nCPU",
        "--nCPU",
        default=1,
        type=int,
        help="Number of CPUs to request for each condor job (default: %(default)s). Job memory scales as nCPU*3GB, adjust if necessary.",
    )
    parser.add_argument(
        "-n",
        "--condorFileSize",
        type=int,
        default=50,
        help="Number of files proceed per condor job (default: %(default)s)",
    )
    parser.add_argument(
        "--outputDir",
        help="Output directory",
        required=True,
    )
    parser.add_argument(
        "--remoteRepo",
        default=None,
        help="If specified, access BTVNanoCommsioning from a remote tarball (downloaded via https), instead of from a transferred sandbox",
    )
    parser.add_argument(
        "--jobqueue",
        default="tomorrow",
        help="JobFlavour for condor@lxplus. E.g. microcentury, longlunch, workday, tomorrow",
    )
    return parser


def get_main_parser():
    parser = argparse.ArgumentParser(description="Arguments for condor submitter")
    ## Inputs
    parser.add_argument(
        "--wf",
        "--workflow",
        dest="workflow",
        help="Which processor to run",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        default=r"hists.coffea",
        help="Output histogram filename (default: %(default)s)",
    )
    parser.add_argument(
        "--samples",
        "--json",
        dest="samplejson",
        nargs="+",
        default="dummy_samples.json",
        help="JSON file containing dataset and file locations (default: %(default)s)",
    )
    ## Configuations
    parser.add_argument("--year", default="2023", help="Year")
    parser.add_argument(
        "--campaign",
        default="Summer23",
        choices=[
            "Rereco17_94X",
            "Winter22Run3",
            "Summer22",
            "Summer22EE",
            "Summer23",
            "Summer23BPix",
            "Summer24",
            "Prompt25",
            "2018-UL",
            "2017-UL",
            "2016preVFP-UL",
            "2016postVFP-UL",
            "CAMPAIGN_prompt_dataMC",
        ],
        help="Dataset campaign, change the corresponding correction files",
    )
    parser.add_argument(
        "--isSyst",
        default=False,
        type=str,
        choices=[
            "False",
            "all",
            "weight_only",
            "JEC_full",
            "JEC_reduced",
            "JEC_reduced_JER_split",
            "JEC_total",
            "JP_MC",
        ],
        help="Run with systematics (default: %(default)s)",
    )
    parser.add_argument(
        "--ttbar-reweights",
        default="none",
        choices=["none", "hdamp_ml", "full"],
        help=(
            "Enable additional ttbar event reweights in correction.py. "
            "'hdamp_ml' applies hdamp ONNX up/down; 'full' additionally reserves "
            "hooks for frag/decay reweights."
        ),
    )
    parser.add_argument("--isArray", action="store_true", help="Output root files")
    parser.add_argument(
        "--array-systs",
        dest="array_systs",
        default="both",
        choices=["none", "weights", "shifts", "both"],
        help="How much systematic content --isArray trees carry; forwarded to "
        "runner.py. See its --array-systs help. Default: %(default)s",
    )
    parser.add_argument(
        "--array-perjet",
        dest="array_perjet",
        action="store_true",
        help="Write --isArray trees with one row per jet; forwarded to runner.py.",
    )

    parser.add_argument(
        "--noHist", action="store_true", help="Not output coffea histogram"
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Overwrite existing files"
    )
    parser.add_argument(
        "--only",
        type=str,
        default=None,
        help="Only  process/skip part of the dataset. By input list of file",
    )
    parser.add_argument(
        "--voms",
        default=None,
        type=str,
        help="Path to voms proxy, made accessible to worker nodes. By default a copy will be made to $HOME.",
    )

    parser.add_argument(
        "--chunk",
        type=int,
        default=75000,
        metavar="N",
        help="Number of events per process chunk",
    )
    parser.add_argument("--skipbadfiles", action="store_true", help="Skip bad files.")
    parser = get_condor_submitter_parser(parser)
    return parser


if __name__ == "__main__":
    parser = get_main_parser()
    args = parser.parse_args()
    print("Running with the following options:")
    print(args)

    # condor_submit delegates the proxy named in the submit file to the jobs, so
    # take it from X509_USER_PROXY (where voms-proxy-init writes when it is set)
    # rather than assuming a fixed location, and refuse to submit without a valid one.
    proxy = os.getenv("X509_USER_PROXY")
    if not proxy or not os.path.isfile(proxy):
        print(f"X509_USER_PROXY does not point to a proxy file ({proxy}).")
        print("Please run the following commands in your shell:")
        print("export X509_USER_PROXY=$HOME/x509up_u`id -u`")
        print("voms-proxy-init --voms cms --valid 192:00")
        sys.exit(1)
    if shutil.which("voms-proxy-info") is not None:
        valid = subprocess.run(
            ["voms-proxy-info", "--exists", "--file", proxy], capture_output=True
        )
        if valid.returncode != 0:
            print(f"The proxy in X509_USER_PROXY ({proxy}) is expired or invalid.")
            print("Please run the following command in your shell:")
            print("voms-proxy-init --voms cms --valid 192:00")
            sys.exit(1)
    print(f"Using the proxy from X509_USER_PROXY: {proxy}")

    current_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = current_dir.replace("/condor", "")

    # Since Nov 2025 the CERN schedds split by filesystem: the standard ones refuse
    # /eos paths in a submit file and the EosSubmit ones (`module load
    # lxbatch/eossubmit`) refuse AFS paths. Everything below references base_dir, so
    # the checkout's filesystem decides which schedd has to receive the job.
    # https://batchdocs.web.cern.ch/local/eossubmit.html
    on_eos = base_dir.startswith("/eos/")
    eossubmit = os.getenv("_myschedd_POOL") == "eossubmit"
    if on_eos and not eossubmit:
        print(
            f"The checkout is on EOS ({base_dir}), which only the EosSubmit schedds "
            "accept. Please run the following command in your shell and resubmit:"
        )
        print("module load lxbatch/eossubmit")
        sys.exit(1)
    if eossubmit and not on_eos:
        print(
            f"The EosSubmit schedds are selected but the checkout is not on EOS "
            f"({base_dir}). Please run the following command in your shell and "
            "resubmit:"
        )
        print("module unload lxbatch/eossubmit")
        sys.exit(1)
    if on_eos:
        # /eos/home-<x>/ and /eos/user/<x>/ are the same volume; the EosSubmit docs
        # only use the latter.
        base_dir = re.sub(r"^/eos/home-([a-z0-9])/", r"/eos/user/\1/", base_dir)

    if args.remoteRepo is not None:
        print("Will use a remote path to access BTVNanoCommissioning:", args.remoteRepo)
    else:
        print("Tarring BTVNanoCommissioning directory...")

        skip_tar = False
        if os.path.exists("BTVNanoCommissioning.tar.gz"):
            user_input = input(
                "BTVNanoCommissioning.tar.gz already exists, skip the tarring? ([y]/n): "
            )
            if user_input.lower() != "n":
                skip_tar = True
            else:
                skip_tar = False
                os.remove("BTVNanoCommissioning.tar.gz")

        if not skip_tar:
            jobdirs = [d for d in os.listdir(base_dir) if d.startswith("jobs_")]
            make_tarfile(
                "BTVNanoCommissioning.tar.gz",
                base_dir,
                exclude_dirs=["BTVNanoCommissioning.egg-info"] + jobdirs,
            )

    # Create job dir
    job_dir = f"jobs_{args.jobName}"
    if os.path.exists(job_dir):
        user_input = input("Job directory already exists, overwrite? ([y]/n): ")
        if user_input.lower() != "n":
            shutil.rmtree(job_dir)
        else:
            raise Exception("Job exiting...")
    print(f"Job directory created: {job_dir}")
    os.mkdir(job_dir)
    os.mkdir(job_dir + "/log")

    # Store job submission files

    ## store parser arguments
    with open(os.path.join(job_dir, "arguments.json"), "w") as json_file:
        json.dump(vars(args), json_file, indent=4)

    ## split the sample json
    if isinstance(args.samplejson, str):
        samplejson = [args.samplejson]
    else:
        samplejson = args.samplejson
    sample_dict = {}
    for js in samplejson:
        with open(js) as f:
            sample_dict.update(json.load(f))

    split_sample_dict = {}
    counter = 0
    only = []
    if args.only is not None:
        if "*" in args.only:
            only = [
                k
                for k in sample_dict.keys()
                if k.lstrip("/").startswith(args.only.rstrip("*"))
            ]
        else:
            only.append(args.only)

    for sample_name, files in sample_dict.items():
        if len(only) != 0 and sample_name not in only:
            continue
        for ifile in range(
            (len(files) + args.condorFileSize - 1) // args.condorFileSize
        ):
            split_sample_dict[counter] = {
                sample_name: files[
                    ifile * args.condorFileSize : (ifile + 1) * args.condorFileSize
                ]
            }
            counter += 1

    ## store the split sample json file
    with open(os.path.join(job_dir, "split_samples.json"), "w") as json_file:
        json.dump(split_sample_dict, json_file, indent=4)
    ## store the jobnum list (0..jobnum-1)
    with open(os.path.join(job_dir, "jobnum_list.txt"), "w") as f:
        f.write("\n".join([str(i) for i in range(counter)]))

    ## store the jdl file
    submit_dir = f"{base_dir}/{job_dir}"
    sandbox = [
        f"{submit_dir}/arguments.json",
        f"{submit_dir}/split_samples.json",
        f"{submit_dir}/jobnum_list.txt",
    ]
    if not args.remoteRepo:
        sandbox.append(f"{base_dir}/BTVNanoCommissioning.tar.gz")

    jdl_template = """Universe   = vanilla
Executable = {executable}
initialdir = {submit_dir}


Arguments = $(JOBNUM) $(request_cpus)

request_cpus = {nCPU}
x509userproxy = {proxy}

+JobFlavour = "{jobqueue}"

Log        = {log_dir}/job.log_$(Cluster)
Output     = {log_dir}/job.out_$(Cluster)-$(Process)
Error      = {log_dir}/job.err_$(Cluster)-$(Process)

max_retries             = 10
periodic_release        = True
should_transfer_files   = YES
when_to_transfer_output = ON_EXIT_OR_EVICT
transfer_input_files    = {transfer_input_files}
JobBatchName            = {batch_name}
transfer_output_files   = .success

Queue JOBNUM from {jobnum_file}
""".format(
        executable=f"{base_dir}/condor/execute.sh",
        proxy=proxy,
        submit_dir=submit_dir,
        jobqueue=args.jobqueue,
        log_dir=f"{submit_dir}/log",
        transfer_input_files=",".join(sandbox),
        nCPU=args.nCPU,
        batch_name=args.jobName,
        jobnum_file=f"{submit_dir}/jobnum_list.txt",
    )
    with open(os.path.join(job_dir, "submit.jdl"), "w") as f:
        f.write(jdl_template)
    os.system(f"condor_submit {job_dir}/submit.jdl")
    # print(
    #     f"Setup completed. Now submit the condor jobs by:\n  condor_submit {job_dir}/submit.jdl"
    # )
