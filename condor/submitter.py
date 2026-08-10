import os, sys
import json
import shutil
import tarfile
import argparse
import re
import subprocess
import time

CONDOR_SUBMIT_SUCCESS_RE = re.compile(
    r"\b(?P<count>\d+)\s+job\(s\)\s+submitted\s+to\s+cluster\s+" r"(?P<cluster>\d+)\b",
    re.IGNORECASE,
)
CONDOR_SUBMIT_PERMANENT_ERROR_RE = re.compile(
    r"(failed to open .*\.jdl|no such file|invalid submit|submit file .*error|"
    r"syntax error)",
    re.IGNORECASE,
)
DEFAULT_SUBMIT_ATTEMPTS = 5
MAX_SUBMIT_RETRY_DELAY_SECONDS = 300


def make_tarfile(output_filename, source_dir, exclude_dirs=[]):
    temporary_filename = f"{output_filename}.tmp"
    try:
        with tarfile.open(temporary_filename, "w:gz") as tar:
            for item in os.listdir(source_dir):
                if item in exclude_dirs:
                    continue
                item_path = os.path.join(source_dir, item)
                if os.path.isdir(item_path):
                    for root, dirs, files in os.walk(item_path):
                        # Ensure we also skip any nested excluded directories
                        dirs[:] = [d for d in dirs if d not in exclude_dirs]
                        for file in files:
                            file_path = os.path.join(root, file)
                            tar.add(
                                file_path,
                                arcname=os.path.relpath(file_path, source_dir),
                            )
                else:
                    # Add top-level files
                    tar.add(item_path, arcname=item)
        os.replace(temporary_filename, output_filename)
    finally:
        if os.path.exists(temporary_filename):
            os.remove(temporary_filename)


def get_condor_submitter_parser(parser, require_job_args=True):
    parser.add_argument(
        "--jobName",
        help="Condor job name to make the job directory",
        required=require_job_args,
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
        required=require_job_args,
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
    parser.add_argument(
        "--reuseTarball",
        action="store_true",
        help="Reuse an existing BTVNanoCommissioning.tar.gz without prompting to recreate it.",
    )
    parser.add_argument(
        "--rebuildTarball",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--submitRetries",
        type=int,
        default=DEFAULT_SUBMIT_ATTEMPTS,
        help=(
            "Maximum number of condor_submit attempts before failing "
            "(default: %(default)s)."
        ),
    )
    return parser


def validate_x509_proxy():
    uid = os.getuid()
    homedir = os.getenv("HOME")
    expected_value = f"{homedir}/x509up_u{uid}"
    current_value = os.getenv("X509_USER_PROXY")
    if current_value != expected_value:
        print("X509_USER_PROXY is NOT set correctly.")
        print("Please run the following command in your shell:")
        print("export X509_USER_PROXY=$HOME/x509up_u`id -u`")
        sys.exit(1)


def prompt_rebuild_tarball():
    while True:
        user_input = input(
            "BTVNanoCommissioning.tar.gz already exists, skip the tarring? ([y]/n): "
        ).strip()
        if user_input in {"y", "Y"}:
            return True
        if user_input in {"n", "N"}:
            return False
        print("Please answer with y/Y or n/N.")


def submit_condor_with_retry(
    submit_jdl_path,
    max_attempts=DEFAULT_SUBMIT_ATTEMPTS,
    retry_delay_seconds=30,
    expected_jobs=None,
):
    if max_attempts < 1:
        raise ValueError("max_attempts must be at least 1")

    last_output = ""
    last_returncode = None
    for attempt in range(1, max_attempts + 1):
        try:
            result = subprocess.run(
                ["condor_submit", submit_jdl_path],
                capture_output=True,
                text=True,
                check=False,
            )
        except OSError as exc:
            raise RuntimeError(f"Could not execute condor_submit: {exc}") from exc

        output = "\n".join(
            part for part in [result.stdout.strip(), result.stderr.strip()] if part
        )
        last_output = output
        last_returncode = result.returncode
        if output:
            print(output)

        match = CONDOR_SUBMIT_SUCCESS_RE.search(output)
        if result.returncode == 0 and match:
            submitted_jobs = int(match.group("count"))
            if expected_jobs is not None and submitted_jobs != expected_jobs:
                raise RuntimeError(
                    "condor_submit reported "
                    f"{submitted_jobs} submitted jobs; expected {expected_jobs}"
                )
            return int(match.group("cluster"))

        if CONDOR_SUBMIT_PERMANENT_ERROR_RE.search(output):
            raise RuntimeError(
                f"condor_submit failed with a permanent error on attempt {attempt}: "
                f"{output or 'no output'}"
            )

        if attempt == max_attempts:
            break

        delay = min(
            retry_delay_seconds * (2 ** (attempt - 1)),
            MAX_SUBMIT_RETRY_DELAY_SECONDS,
        )
        print(
            f"condor_submit did not report a successful submission on attempt {attempt}; "
            f"retrying in {delay}s."
        )
        time.sleep(delay)

    raise RuntimeError(
        f"condor_submit failed after {max_attempts} attempts "
        f"(last return code: {last_returncode}): {last_output or 'no output'}"
    )


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
            "all_withJESTotal",
            "weight_only",
            "JEC_full",
            "JEC_reduced",
            "JEC_reduced_JER_split",
            "JEC_total",
            "JERC_full",
            "JERC_reduced",
            "JERC_total",
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
    if args.isSyst in {"JERC_full", "JERC_reduced", "JERC_total"}:
        args.isSyst = args.isSyst.replace("JERC", "JEC")
    print("Running with the following options:")
    print(args)
    validate_x509_proxy()

    current_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = current_dir.replace("/condor", "")
    if base_dir not in sys.path:
        sys.path.insert(0, base_dir)
    from condor_lxplus import dashboard as job_dashboard

    # Extract version for setuptools_scm fallback in worker nodes
    try:
        with open(os.path.join(base_dir, "src/BTVNanoCommissioning/version.py")) as f:
            v_content = f.read()
            v_match = re.search(r"version = ['\"]([^'\"]+)['\"]", v_content)
            scm_version = v_match.group(1) if v_match else "0.1"
    except Exception:
        scm_version = "0.1"
    setattr(args, "scm_version", scm_version)

    if args.remoteRepo is not None:
        print("Will use a remote path to access BTVNanoCommissioning:", args.remoteRepo)
    else:
        print("Tarring BTVNanoCommissioning directory...")

        if args.reuseTarball and args.rebuildTarball:
            raise ValueError(
                "--reuseTarball and --rebuildTarball are mutually exclusive"
            )

        skip_tar = False
        if os.path.exists("BTVNanoCommissioning.tar.gz"):
            if args.rebuildTarball:
                print("Rebuilding BTVNanoCommissioning.tar.gz")
            elif args.reuseTarball:
                print("Reusing existing BTVNanoCommissioning.tar.gz")
                skip_tar = True
            else:
                skip_tar = prompt_rebuild_tarball()

        if not skip_tar:
            exclude_list = ["jsonpog-integration", "BTVNanoCommissioning.egg-info"]
            for d in os.listdir(base_dir):
                if (
                    d.startswith("jobs_")
                    or d.startswith("arrays_")
                    or d.startswith("hists_")
                    or d.startswith("condor")
                    or d.startswith(".")
                ):
                    exclude_list.append(d)
            make_tarfile(
                "BTVNanoCommissioning.tar.gz",
                base_dir,
                exclude_dirs=exclude_list,
            )

    # Create job dir
    job_dir = f"jobs_{args.jobName}_{args.campaign}"
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
    jdl_template = """Universe   = vanilla
Executable = {executable}


Arguments = $(JOBNUM) $(request_cpus)

request_cpus = {nCPU}
use_x509userproxy = true

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
        jobqueue=args.jobqueue,
        log_dir=f"{base_dir}/{job_dir}/log",
        transfer_input_files=f"{base_dir}/{job_dir}/arguments.json,{base_dir}/{job_dir}/split_samples.json,{base_dir}/{job_dir}/jobnum_list.txt"
        + ("" if args.remoteRepo else f",{base_dir}/BTVNanoCommissioning.tar.gz"),
        nCPU=args.nCPU,
        batch_name=args.jobName,
        jobnum_file=f"{base_dir}/{job_dir}/jobnum_list.txt",
    )
    with open(os.path.join(job_dir, "submit.jdl"), "w") as f:
        f.write(jdl_template)
    cluster_id = submit_condor_with_retry(
        f"{job_dir}/submit.jdl",
        max_attempts=args.submitRetries,
        expected_jobs=counter,
    )
    job_dashboard.record_submission(
        job_dir,
        counter,
        job_name=args.jobName,
        output_dir=args.outputDir,
        cluster_id=cluster_id,
    )
    # print(
    #     f"Setup completed. Now submit the condor jobs by:\n  condor_submit {job_dir}/submit.jdl"
    # )
