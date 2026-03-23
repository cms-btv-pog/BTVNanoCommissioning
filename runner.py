import os
import sys
import json
import argparse
import time

import numpy as np

import uproot
from coffea.util import load, save
from coffea import processor
from coffea.nanoevents import PFNanoAODSchema
from BTVNanoCommissioning.workflows import workflows


def validate(file):
    try:
        fin = uproot.open(file)
        return fin["Events"].num_entries
    except:
        print("Corrupted file: {}".format(file))
        return file


def validate_dataset_structure(fileset):
    """Check dataset files and return a filtered fileset with only valid files."""
    import uproot
    import logging
    from copy import deepcopy

    # Critical branches that must be present
    required_branches = [
        "Jet_pt",
        "Muon_pt",
        "Electron_pt",
        "PV_npvsGood",
        "Jet_eta",
        "Jet_phi",
        "Muon_eta",
        "Muon_phi",
        "SV_pt",
    ]
    optional_met_branches = ["MET_pt", "PuppiMET_pt", "PFMET_pt"]

    filtered_fileset = deepcopy(fileset)
    all_files_invalid = True

    # Check each sample in the fileset
    for sample_name, files in fileset.items():
        valid_files = []

        # Check each file in the sample
        for filename in files:
            try:
                # print(f"Validating file: {filename}")
                file = uproot.open(filename)
                events = file["Events"]
                branches = set(events.keys())

                # Check branch count
                if len(branches) < 20:
                    print(
                        f"WARNING: File has too few branches ({len(branches)}): {filename}"
                    )
                    # print(f"Branches present: {branches}")
                    continue

                # Check required branches
                missing = [b for b in required_branches if b not in branches]
                if missing:
                    print(
                        f"WARNING: File missing critical branches {missing}: {filename}"
                    )
                    continue

                if not any(met in branches for met in optional_met_branches):
                    print(
                        "WARNING: File missing MET branch; expected one of "
                        f"{optional_met_branches}: {filename}"
                    )
                    continue

                # Check event count
                if len(events) < 10:
                    print(
                        f"WARNING: File has too few events ({len(events)}): {filename}"
                    )
                    continue

                # File passed all checks
                valid_files.append(filename)
                all_files_invalid = False
                # print(f"File validation successful: {filename}")

            except Exception as e:
                print(f"ERROR validating file: {filename}, {e}")
                continue

        # Update the filtered fileset with valid files
        if valid_files:
            filtered_fileset[sample_name] = valid_files
        else:
            # Remove the sample if no valid files
            del filtered_fileset[sample_name]

    # Summary
    if all_files_invalid:
        print("WARNING: All files in dataset failed validation!")
        return None
    else:
        print(
            f"Dataset validation complete. {len(filtered_fileset)} valid samples remaining."
        )
        return filtered_fileset


def check_port(port):
    import socket

    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    try:
        sock.bind(("0.0.0.0", port))
        available = True
    except:
        available = False
    sock.close()
    return available


def retry_handler(exception, task_record):
    from parsl.executors.high_throughput.interchange import ManagerLost

    if isinstance(exception, ManagerLost):
        return 0.1
    else:
        return 1


## From condor/submitter.py https://github.com/cms-btv-pog/BTVNanoCommissioning/blob/9edb9ed6bb0b28730b8de9e5aa1142ec4fdf74b7/condor/submitter.py
# create t


def make_tarfile(output_filename, source_dir, exclude_dirs=[]):
    import tarfile

    with tarfile.open(output_filename, "w:gz") as tar:
        for root, dirs, files in os.walk(source_dir):
            dirs[:] = [
                d for d in dirs if d not in exclude_dirs
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
        "--condorFileSize",
        type=int,
        default=50,
        help="Number of files proceed per condor job",
    )
    parser.add_argument(
        "--outputXrootdDir",
        help="Output directory for xrootd files. Start with root://",
        required=True,
    )
    parser.add_argument(
        "--remoteRepo",
        default=None,
        help="If specified, access BTVNanoCommsioning from a remote tarball (downloaded via https), instead of from a transferred sandbox",
    )
    return parser


def IO_parser(parser):
    ## Inputs
    parser.add_argument(
        "--wf",
        "--workflow",
        dest="workflow",
        choices=list(workflows.keys()),
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
        "--outputdir",
        default=None,
        help="Output directory name (default: %(default)s)",
    )
    parser.add_argument(
        "--json",
        dest="samplejson",
        default="dummy_samples.json",
        help="JSON file containing dataset and file locations (default: %(default)s)",
    )
    return parser


def config_parser(parser):
    ## Configuations
    parser.add_argument("--year", default="2023", required=True, help="Year")
    parser.add_argument(
        "--campaign",
        required=True,
        default="Summer23",
        help="Dataset campaign, change the corresponding correction files",
    )
    parser.add_argument(
        "--isSyst",
        default="False",
        type=str,
        choices=["False", "all", "weight_only", "JERC_split", "JP_MC"],
        help="Run with systematics, all, weights_only(no JERC uncertainties included),JERC_split, None",
    )
    parser.add_argument("--isArray", action="store_true", help="Output root files")
    parser.add_argument(
        "--noHist", action="store_true", help="Not output coffea histogram"
    )
    parser.add_argument(
        "--overwrite", action="store_true", help="Overwrite existing files"
    )
    return parser


def scaleout_parser(parser):
    # Scale out
    parser.add_argument(
        "--executor",
        choices=[
            "iterative",
            "futures",
            "parsl/slurm",
            "parsl/condor",
            "parsl/condor/naf_lite",
            "dask/condor",
            "dask/condor/brux",
            "dask/slurm",
            "dask/lpc",
            "dask/lxplus",
            "dask/casa",
            "condor_standalone",
        ],
        default="futures",
        help="The type of executor to use (default: %(default)s). Other options can be implemented. "
        "For example see https://parsl.readthedocs.io/en/stable/userguide/configuring.html"
        "- `parsl/slurm` - tested at DESY/Maxwell"
        "- `parsl/condor` - tested at DESY, RWTH"
        "- `parsl/condor/naf_lite` - tested at DESY"
        "- `dask/condor/brux` - tested at BRUX (Brown U)"
        "- `dask/slurm` - tested at DESY/Maxwell"
        "- `dask/condor` - tested at DESY, RWTH"
        "- `dask/lpc` - custom lpc/condor setup (due to write access restrictions)"
        "- `dask/lxplus` - custom lxplus/condor setup (due to port restrictions)",
    )
    parser.add_argument(
        "-j",
        "--workers",
        type=int,
        default=3,
        help="Number of workers (cores/threads) to use for multi-worker executors "
        "(e.g. futures or condor) (default: %(default)s)",
    )
    parser.add_argument(
        "-s",
        "--scaleout",
        type=int,
        default=6,
        help="Number of nodes to scale out to if using slurm/condor. Total number of "
        "concurrent threads is ``workers x scaleout`` (default: %(default)s)",
    )
    parser.add_argument(
        "--memory",
        type=float,
        default=4.0,
        help="Memory used in jobs default ``(default: %(default)s)",
    )
    parser.add_argument(
        "--disk",
        type=float,
        default=4,
        help="Disk used in jobs default ``(default: %(default)s)",
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
        default=200000,
        metavar="N",
        help="Number of events per process chunk",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=25,
        metavar="N",
        help="Number of retries for coffea processor",
    )
    parser.add_argument(
        "--fsize",
        type=int,
        default=50,
        help="(Specific for dask/lxplus file splitting, default: %(default)s)\n Numbers of files processed per dask-worker",
    )
    parser.add_argument(
        "--index",
        type=str,
        default="0,0",
        help=f"(Specific for dask/lxplus file splitting, default: %(default)s)\n   Format: $dict_index_start,$file_index_start,$dict_index_stop,$file_index_stop. Stop indices are optional. $dict_index refers to the index, splitted $dict_index and $file_index with ','"
        "$dict_index refers to the sample dictionary of the samples json file. $file_index refers to the N-th batch of files per dask-worker, with its size being defined by the option --index. The job will start (stop) submission from (with) the corresponding indices.",
    )
    return parser


def debug_parser(parser):
    # Debugging
    parser.add_argument(
        "--validate",
        action="store_true",
        help="Do not process, just check all files are accessible",
    )
    parser.add_argument("--skipbadfiles", action="store_true", help="Skip bad files.")
    parser.add_argument(
        "--only", type=str, default=None, help="Only process specific dataset or file"
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        metavar="N",
        help="Limit to the first N files of each dataset in sample JSON",
    )
    parser.add_argument(
        "--max",
        type=int,
        default=None,
        metavar="N",
        help="Max number of chunks to run in total",
    )
    return parser


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run analysis on baconbits files using processor coffea files"
    )
    parser = IO_parser(parser)
    parser = config_parser(parser)
    parser = scaleout_parser(parser)
    parser = debug_parser(parser)
    args = parser.parse_args()
    if args.isSyst == "False":
        args.isSyst = False
    print("Running with the following options:")
    print(args)
    ogoutput = args.output
    histoutdir = ogoutput.split(".")[0]
    outdir = "arrays_" + histoutdir
    outdirprefix = ""
    if args.outputdir:
        outdirprefix = f"{args.outputdir}/"
        histoutdir = f"{outdirprefix}{histoutdir}"
        outdir = f"{outdirprefix}{outdir}"
    coffeaoutput = f"{histoutdir}/{ogoutput}"
    basename = ogoutput.replace(".coffea", "").replace("hists_", "")
    if args.output == parser.get_default("output"):
        index = args.samplejson.rfind("/") + 1
        sample_json = args.samplejson[index:]
        histoutdir = (
            f"{outdirprefix}hists_{args.workflow}_{sample_json.rstrip('.json')}"
        )
        outdir = f"{outdirprefix}arrays_{args.workflow}_{sample_json.rstrip('.json')}"
        coffeaoutput = (
            f'{histoutdir}/hists_{args.workflow}_{(sample_json).rstrip(".json")}.coffea'
        )
    os.system(f"mkdir -p {histoutdir}")
    # load dataset
    with open(args.samplejson) as f:
        sample_dict = json.load(f)
    for key in sample_dict.keys():
        sample_dict[key] = sample_dict[key][: args.limit]
    if args.executor == "dask/casa":
        for key in sample_dict.keys():
            sample_dict[key] = [
                path.replace("xrootd-cms.infn.it/", "xcache")
                for path in sample_dict[key]
            ]
    # check file dict size - avoid large memory consumption for local machine
    filesize = np.sum(np.array([len(sample_dict[key]) for key in sample_dict.keys()]))
    splitjobs = False
    if filesize > 200 and "lxplus" in args.executor:
        splitjobs = True

    # For debugging
    if args.only is not None:
        if args.only in sample_dict.keys():  # is dataset
            sample_dict = dict([(args.only, sample_dict[args.only])])
            coffeaoutput = coffeaoutput.replace(".coffea", f"_{args.only}.coffea")
        elif args.only.isdigit():
            isamp = int(args.only)
            nsamp = len(sample_dict.keys())
            if isamp >= nsamp:
                print(
                    f"There are {nsamp} datasets, please use --only n with n<{nsamp}."
                )
            key = list(sample_dict.keys())[isamp]
            print(f"Will process only {key} instead of all {nsamp} datasets.")
            sample_dict = dict([(key, sample_dict[key])])
            coffeaoutput = coffeaoutput.replace(".coffea", f"_{key}.coffea")
        elif "*" in args.only:  # wildcard for datasets
            _new_dict = {}
            print("Will only proces the following datasets:")
            for k, v in sample_dict.items():
                if args.only.replace("*", "") in k:
                    print("    ", k)
                    _new_dict[k] = v
            sample_dict = _new_dict
        else:  # is file
            for key in sample_dict.keys():
                if args.only in sample_dict[key]:
                    sample_dict = dict([(key, [args.only])])

    # Scan if files can be opened
    if args.validate:
        start = time.time()
        from p_tqdm import p_map

        all_invalid = []
        for sample in sample_dict.keys():
            _rmap = p_map(
                validate,
                sample_dict[sample],
                num_cpus=args.workers,
                desc=f"Validating {sample[:20]}...",
            )
            _results = list(_rmap)
            counts = np.sum(
                [r for r in _results if np.isreal(r) and isinstance(r, int)]
            )
            all_invalid += [r for r in _results if isinstance(r, str)]
            print("Events:", np.sum(counts))
        print("Bad files:")
        for fi in all_invalid:
            print(f"  {fi}")
        end = time.time()
        print("TIME:", time.strftime("%H:%M:%S", time.gmtime(end - start)))
        if len(all_invalid) == 0:
            print("No bad files found!")
        else:
            if input("Remove bad files? (y/n): ") == "y":
                print("Removing...")
                json = args.samplejson
                jsonnew = json.replace(".json", "") + "_backup.json"
                os.system("mv %s %s" % (json, jsonnew))
                inf = open(jsonnew, "r")
                outf = open(json, "w")
                for line in inf:
                    foundline = False
                    for fi in all_invalid:
                        if fi in line:
                            print(f"Removing: {fi}")
                            foundline = True
                            break
                    if not foundline:
                        outf.write(line)
        sys.exit(0)

    # load workflow
    processor_instance = workflows[args.workflow](
        args.year,
        args.campaign,
        outdir,
        args.isSyst,
        args.isArray,
        args.noHist,
        args.chunk,
    )

    filtered_sample_dict = validate_dataset_structure(sample_dict)

    if filtered_sample_dict is None:
        print(f"Creating empty output file for incompatible dataset: {args.samplejson}")
        # Create minimal empty output
        from coffea.util import save

        empty_output = {}  # Minimal dict output
        outname = os.path.join(args.outputdir, args.output)
        os.makedirs(os.path.dirname(outname), exist_ok=True)
        save(empty_output, outname)
        print(f"Empty output file created successfully at {outname}")
        sys.exit(0)  # Exit with success code
    else:
        # Continue with processing using the filtered_sample_dict
        sample_dict = filtered_sample_dict

    ## create tmp directory and check file exist or not
    from os import path

    if path.exists(f"{coffeaoutput}") and args.overwrite == False:
        raise Exception(f"{coffeaoutput} exists")

    if args.isArray:
        if path.exists(outdir) and args.overwrite == False and args.only is None:
            raise Exception(f"Directory {outdir} exists")
        else:
            os.system(f"mkdir -p {outdir}")

    if args.executor not in ["futures", "iterative", "dask/lpc", "dask/casa"]:
        """
        dask/parsl needs to export x509 to read over xrootd
        dask/lpc uses custom jobqueue provider that handles x509
        """
        if args.voms is not None:
            _x509_path = args.voms
        else:
            try:
                _x509_localpath = (
                    [
                        l
                        for l in os.popen("voms-proxy-info").read().split("\n")
                        if l.startswith("path")
                    ][0]
                    .split(":")[-1]
                    .strip()
                )
            except:
                raise RuntimeError(
                    "x509 proxy could not be parsed, try creating it with 'voms-proxy-init'"
                )
            _x509_path = os.environ["HOME"] + f'/.{_x509_localpath.split("/")[-1]}'
            os.system(f"cp {_x509_localpath} {_x509_path}")

        job_script_prologue = [
            "export XRD_RUNFORKHANDLER=1",
            f"export X509_USER_PROXY={_x509_path}",
            f'export X509_CERT_DIR={os.environ["X509_CERT_DIR"]}',
            f"export PYTHONPATH=$PYTHONPATH:{os.getcwd()}",
        ]
        pathvar = [i for i in os.environ["PATH"].split(":") if "envs/btv_coffea" in i][
            0
        ]
        condor_extra = [
            f"export PATH={pathvar}:$PATH",
            f'source {os.environ["HOME"]}/.bashrc',
        ]
        if "brux" in args.executor:
            job_script_prologue.append(f"cd {os.getcwd()}")
        else:
            condor_extra.append(f"cd {os.getcwd()}")

            # Check if Conda is available
            conda_check_command = "command -v conda"
            conda_available = os.system(conda_check_command) == 0

            # Check if Mamba is available
            mamba_check_command = "command -v micromamba"
            mamba_available = os.system(mamba_check_command) == 0

            # Set up environment based on availability
            if conda_available and mamba_available:
                use_conda = True  # Set to False if you prefer Micromamba
                if use_conda:
                    condor_extra.append(f'conda activate {os.environ["CONDA_PREFIX"]}')
                else:
                    condor_extra.append(
                        f"micromamba activate {os.environ['CONDA_PREFIX']}"
                    )
            elif conda_available:
                condor_extra.append(f'conda activate {os.environ["CONDA_PREFIX"]}')
            elif mamba_available:
                condor_extra.append(f"micromamba activate {os.environ['CONDA_PREFIX']}")
            else:
                # Handle the case when neither Conda nor Micromamba is available
                print(
                    "Neither Conda nor Micromamba is available in the environment. At least install one of them."
                )

    #########
    # Execute
    if args.executor in ["futures", "iterative", "standalone_condor"]:
        if args.executor == "iterative":
            _exec = processor.iterative_executor
        else:
            _exec = processor.futures_executor
        if args.executor != "standalone_condor":
            output = processor.run_uproot_job(
                sample_dict,
                treename="Events",
                processor_instance=processor_instance,
                executor=_exec,
                executor_args={
                    "skipbadfiles": args.skipbadfiles,
                    "schema": PFNanoAODSchema,
                    "workers": args.workers,
                    "xrootdtimeout": 900,
                },
                chunksize=args.chunk,
                maxchunks=args.max,
            )

        ## standalone condor from https://github.com/cms-btv-pog/BTVNanoCommissioning/blob/9edb9ed6bb0b28730b8de9e5aa1142ec4fdf74b7/condor/submitter.py
        else:
            current_dir = os.path.dirname(os.path.abspath(__file__))
            base_dir = current_dir.replace("/condor", "")
            if args.remoteRepo is not None:
                print(
                    "Will use a remote path to access BTVNanoCommissioning:",
                    args.remoteRepo,
                )
            else:
                print("Tarring BTVNanoCommissioning directory...")

                skip_tar = False
                if os.path.exists("BTVNanoCommissioning.tar.gz"):
                    user_input = input(
                        "BTVNanoCommissioning.tar.gz already exists, skip the tarring? (y/n): "
                    )
                    if user_input.lower() == "y":
                        skip_tar = True
                    elif user_input.lower() == "n":
                        skip_tar = False
                        os.remove("BTVNanoCommissioning.tar.gz")
                    else:
                        raise Exception("Invalid input, exiting")

                if not skip_tar:
                    make_tarfile(
                        "BTVNanoCommissioning.tar.gz",
                        base_dir,
                        exclude_dirs=[
                            "jsonpog-integration",
                            "BTVNanoCommissioning.egg-info",
                        ],
                    )
            import shutil

            # Create job dir
            job_dir = f"jobs_{args.jobName}"
            if os.path.exists(job_dir):
                user_input = input("Job directory already exists, overwrite? (y/n): ")
                if user_input.lower() == "y":
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
            with open(args.samplejson) as f:
                sample_dict = json.load(f)
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
                            ifile
                            * args.condorFileSize : (ifile + 1)
                            * args.condorFileSize
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
            jdl_template = """Universe   = vanilla
        Executable = {executable}

        +ProjectName="cms.org.cern"

        Arguments = $(JOBNUM)

        request_cpus = 1
        request_memory = 2000
        use_x509userproxy = true

        +JobFlavour = "tomorrow"

        Log        = {log_dir}/job.log_$(Cluster)
        Output     = {log_dir}/job.out_$(Cluster)-$(Process)
        Error      = {log_dir}/job.err_$(Cluster)-$(Process)

        should_transfer_files   = YES
        when_to_transfer_output = ON_EXIT_OR_EVICT
        transfer_input_files    = {transfer_input_files}
        transfer_output_files   = .success

        Queue JOBNUM from {jobnum_file}
        """.format(
                executable=f"{base_dir}/condor/execute.sh",
                log_dir=f"{base_dir}/{job_dir}/log",
                transfer_input_files=f"{base_dir}/{job_dir}/arguments.json,{base_dir}/{job_dir}/split_samples.json,{base_dir}/{job_dir}/jobnum_list.txt"
                + (
                    ""
                    if args.remoteRepo
                    else f",{base_dir}/BTVNanoCommissioning.tar.gz"
                ),
                jobnum_file=f"{base_dir}/{job_dir}/jobnum_list.txt",
            )
            with open(os.path.join(job_dir, "submit.jdl"), "w") as f:
                f.write(jdl_template)

            print(
                f"Setup completed. Now submit the condor jobs by:\n  condor_submit {job_dir}/submit.jdl"
            )

    elif "parsl" in args.executor:
        import parsl
        from parsl.providers import LocalProvider, CondorProvider, SlurmProvider
        from parsl.channels import LocalChannel
        from parsl.config import Config
        from parsl.executors import HighThroughputExecutor
        from parsl.launchers import SrunLauncher
        from parsl.addresses import address_by_hostname, address_by_query

        if "slurm" in args.executor:
            htex_config = Config(
                executors=[
                    HighThroughputExecutor(
                        label="coffea_parsl_slurm",
                        address=address_by_hostname(),
                        prefetch_capacity=0,
                        provider=SlurmProvider(
                            channel=LocalChannel(script_dir="logs_parsl"),
                            launcher=SrunLauncher(),
                            mem_per_node=args.memory,
                            max_blocks=(args.scaleout) + 10,
                            init_blocks=args.scaleout,
                            partition="all",
                            worker_init="\n".join(job_script_prologue),
                            walltime="00:120:00",
                        ),
                    )
                ],
                retries=args.retries,
            )
            if splitjobs:
                htex_config = Config(
                    executors=[
                        HighThroughputExecutor(
                            label="run",
                            address=address_by_hostname(),
                            prefetch_capacity=0,
                            provider=SlurmProvider(
                                channel=LocalChannel(script_dir="logs_parsl"),
                                launcher=SrunLauncher(),
                                mem_per_node=args.memory,
                                max_blocks=(args.scaleout) + 10,
                                init_blocks=args.scaleout,
                                partition="all",
                                worker_init="\n".join(job_script_prologue),
                                walltime="00:120:00",
                            ),
                        ),
                        HighThroughputExecutor(
                            label="merge",
                            address=address_by_hostname(),
                            prefetch_capacity=0,
                            provider=SlurmProvider(
                                channel=LocalChannel(script_dir="logs_parsl"),
                                launcher=SrunLauncher(),
                                mem_per_node=args.memory,
                                max_blocks=(args.scaleout) + 10,
                                init_blocks=args.scaleout,
                                partition="all",
                                worker_init="\n".join(job_script_prologue),
                                walltime="00:30:00",
                            ),
                        ),
                    ],
                    retries=args.retries,
                    retry_handler=retry_handler,
                )
        elif "condor" in args.executor:
            if "naf_lite" in args.executor:
                ## code source: https://github.com/cms-rwth/CoffeaRunner/commit/d5ef86f76723e75b67bb212c3644c4012cae05be (Annika Stein)
                htex_config = Config(
                    executors=[
                        HighThroughputExecutor(
                            label="coffea_parsl_condor",
                            address=address_by_query(),
                            max_workers=1,
                            worker_debug=True,
                            provider=CondorProvider(
                                nodes_per_block=1,
                                cores_per_slot=args.workers,
                                mem_per_slot=2,  # lite job / opportunistic can only use this much
                                init_blocks=args.scaleout,
                                max_blocks=args.scaleout + 5,
                                worker_init="\n".join(
                                    job_script_prologue + condor_extra
                                ),
                                walltime="03:00:00",  # lite / short queue requirement
                            ),
                        )
                    ],
                    retries=args.retries,
                    retry_handler=retry_handler,
                )
                if splitjobs:
                    htex_config = Config(
                        executors=[
                            HighThroughputExecutor(
                                label="run",
                                address=address_by_query(),
                                max_workers=1,
                                worker_debug=True,
                                provider=CondorProvider(
                                    nodes_per_block=1,
                                    cores_per_slot=args.workers,
                                    mem_per_slot=2,  # lite job / opportunistic can only use this much
                                    init_blocks=args.scaleout,
                                    max_blocks=args.scaleout + 5,
                                    worker_init="\n".join(
                                        job_script_prologue + condor_extra
                                    ),
                                    walltime="03:00:00",  # lite / short queue requirement
                                ),
                            ),
                            HighThroughputExecutor(
                                label="merge",
                                address=address_by_query(),
                                max_workers=1,
                                worker_debug=True,
                                provider=CondorProvider(
                                    nodes_per_block=1,
                                    cores_per_slot=args.workers,
                                    mem_per_slot=2,  # lite job / opportunistic can only use this much
                                    init_blocks=args.scaleout,
                                    max_blocks=args.scaleout + 5,
                                    worker_init="\n".join(
                                        job_script_prologue + condor_extra
                                    ),
                                    walltime="00:30:00",  # lite / short queue requirement
                                ),
                            ),
                        ],
                        retries=args.retries,
                        retry_handler=retry_handler,
                    )

            else:
                htex_config = Config(
                    executors=[
                        HighThroughputExecutor(
                            label="coffea_parsl_condor",
                            address=address_by_query(),
                            max_workers=1,
                            provider=CondorProvider(
                                nodes_per_block=1,
                                cores_per_slot=args.workers,
                                mem_per_slot=args.memory,
                                init_blocks=args.scaleout,
                                max_blocks=(args.scaleout) + 10,
                                worker_init="\n".join(
                                    job_script_prologue + condor_extra
                                ),
                                walltime="03:00:00",
                            ),
                        )
                    ],
                    retries=args.retries,
                )
                if splitjobs:
                    htex_config = Config(
                        executors=[
                            HighThroughputExecutor(
                                label="run",
                                address=address_by_query(),
                                max_workers=1,
                                provider=CondorProvider(
                                    nodes_per_block=1,
                                    cores_per_slot=args.workers,
                                    mem_per_slot=args.memory,
                                    init_blocks=args.scaleout,
                                    max_blocks=(args.scaleout) + 10,
                                    worker_init="\n".join(
                                        job_script_prologue + condor_extra
                                    ),
                                    walltime="03:00:00",
                                ),
                            ),
                            HighThroughputExecutor(
                                label="merge",
                                address=address_by_query(),
                                max_workers=1,
                                provider=CondorProvider(
                                    nodes_per_block=1,
                                    cores_per_slot=args.workers,
                                    mem_per_slot=args.memory,
                                    init_blocks=args.scaleout,
                                    max_blocks=(args.scaleout) + 10,
                                    worker_init="\n".join(
                                        job_script_prologue + condor_extra
                                    ),
                                    walltime="03:00:00",
                                ),
                            ),
                        ],
                        retries=args.retries,
                        retry_handler=retry_handler,
                    )

        else:
            raise NotImplementedError

        dfk = parsl.load(htex_config)
        if not splitjobs:
            output = processor.run_uproot_job(
                sample_dict,
                treename="Events",
                processor_instance=processor_instance,
                executor=processor.parsl_executor,
                executor_args={
                    "skipbadfiles": args.skipbadfiles,
                    "schema": PFNanoAODSchema,
                    "config": None,
                },
                chunksize=args.chunk,
                maxchunks=args.max,
            )
        else:
            output = processor.run_uproot_job(
                sample_dict,
                treename="Events",
                processor_instance=processor_instance,
                executor=processor.parsl_executor,
                executor_args={
                    "skipbadfiles": args.skipbadfiles,
                    "schema": PFNanoAODSchema,
                    "merging": True,
                    "merges_executors": ["merge"],
                    "jobs_executors": ["run"],
                    "config": None,
                },
                chunksize=args.chunk,
                maxchunks=args.max,
            )
    elif "dask" in args.executor:
        from dask_jobqueue import SLURMCluster, HTCondorCluster
        from distributed import Client
        from dask.distributed import performance_report

        if "lpc" in args.executor:
            job_script_prologue = [f"export PYTHONPATH=$PYTHONPATH:{os.getcwd()}"]
            from lpcjobqueue import LPCCondorCluster

            cluster = LPCCondorCluster(
                transfer_input_files="/srv/src/",
                ship_env=True,
                job_script_prologue=job_script_prologue,
            )
        elif "lxplus" in args.executor:
            # details: https://batchdocs.web.cern.ch/specialpayload/dask.html
            n_port = 8786
            if not check_port(n_port):
                raise RuntimeError(
                    f"Port '{n_port}' is not occupied on this node. Try another one."
                )
            import socket

            daskjobnm = (
                "dask_job_output" + os.uname()[1].split("lxplus")[-1].split(".cern")[0]
            )
            cluster = HTCondorCluster(
                cores=1,
                memory="2GB",  # hardcoded
                disk="1GB",
                death_timeout="300",
                nanny=False,
                scheduler_options={"port": n_port, "host": socket.gethostname()},
                job_extra_directives={
                    "log": daskjobnm + ".log",
                    "output": daskjobnm + ".out",
                    "error": daskjobnm + ".err",
                    "should_transfer_files": "Yes",
                    "when_to_transfer_output": "ON_EXIT",
                    "+JobFlavour": '"workday"',
                },
                worker_extra_args=["--worker-port 10000:10100"],
                job_script_prologue=job_script_prologue,
            )
        elif "slurm" in args.executor:
            cluster = SLURMCluster(
                queue="all",
                cores=args.workers,
                processes=args.scaleout,
                memory=f"{args.memory}GB",
                disk=f"{args.disk}GB",
                retries=args.retries,
                walltime="00:30:00",
                job_script_prologue=job_script_prologue,
            )
        elif "condor" in args.executor:
            portopts = {}
            if "brux" in args.executor:
                import socket

                portopts = {"host": socket.gethostname()}
            cluster = HTCondorCluster(
                cores=args.workers,
                memory=f"{args.memory}GB",
                disk=f"{args.disk}GB",
                scheduler_options=portopts,
                job_script_prologue=job_script_prologue,
            )

        if args.executor == "dask/casa":
            client = Client("tls://localhost:8786")
            import shutil

            shutil.make_archive("workflows", "zip", base_dir="workflows")
            client.upload_file("workflows.zip")
        else:
            if "brux" in args.executor:
                cluster.adapt(minimum=args.scaleout, maximum=336)
            else:
                cluster.adapt(minimum=args.scaleout)
            client = Client(cluster)
            print("Waiting for at least one worker...")
            client.wait_for_workers(1)
        with performance_report(filename="dask-report.html"):
            if args.executor != "dask/lxplus":
                output = processor.run_uproot_job(
                    sample_dict,
                    treename="Events",
                    processor_instance=processor_instance,
                    executor=processor.dask_executor,
                    executor_args={
                        "client": client,
                        "skipbadfiles": args.skipbadfiles,
                        "schema": PFNanoAODSchema,
                        "retries": args.retries,
                    },
                    chunksize=args.chunk,
                    maxchunks=args.max,
                )

            else:
                findex = int(args.index.split(",")[1])
                for sindex, sample in enumerate(sample_dict.keys()):
                    if sindex < int(args.index.split(",")[0]):
                        continue
                    if int(args.index.split(",")[1]) == findex:
                        mins = findex * args.fsize
                    else:
                        mins = 0
                        findex = 0
                    while mins < len(sample_dict[sample]):
                        splitted = {}
                        maxs = mins + args.fsize
                        splitted[sample] = sample_dict[sample][mins:maxs]
                        mins = maxs
                        findex = findex + 1
                        if (
                            len(args.index.split(",")) == 4
                            and findex > int(args.index.split(",")[3])
                            and sindex > int(args.index.split(",")[2])
                        ):
                            break
                        output = processor.run_uproot_job(
                            splitted,
                            treename="Events",
                            processor_instance=processor_instance,
                            executor=processor.dask_executor,
                            executor_args={
                                "client": client,
                                "skipbadfiles": args.skipbadfiles,
                                "schema": PFNanoAODSchema,
                                "retries": args.retries,
                            },
                            chunksize=args.chunk,
                            maxchunks=args.max,
                        )
                        if args.noHist == False:
                            save(
                                output,
                                coffeaoutput.replace(
                                    ".coffea", f"_{sindex}_{findex}.coffea"
                                ),
                            )
    if not "lxplus" in args.executor:
        if args.noHist == False:
            save(output, coffeaoutput)
    if args.noHist == False:
        # print(output)
        print(f"Saving output to {coffeaoutput}")
