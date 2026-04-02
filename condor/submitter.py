import os, sys
import json
import shutil
import tarfile
import argparse


def make_tarfile(output_filename, source_dir, exclude_dirs=[]):
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
        "-n",
        "--condorFileSize",
        type=int,
        default=50,
        help="Number of files proceed per condor job",
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
            "JERC_full",
            "JERC_reduced",
            "JERC_total",
            "JP_MC",
        ],
        help="Run with systematics (default: %(default)s)",
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
    print("Running with the following options:")
    print(args)

    uid = os.getuid()
    homedir = os.getenv("HOME")
    expected_value = f"{homedir}/x509up_u{uid}"
    current_value = os.getenv("X509_USER_PROXY")
    if current_value != expected_value:
        print("X509_USER_PROXY is NOT set correctly.")
        print(f"Please run the following command in your shell:")
        print(f"export X509_USER_PROXY=$HOME/x509up_u`id -u`")
        sys.exit(1)

    current_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = current_dir.replace("/condor", "")

    if args.remoteRepo is not None:
        print("Will use a remote path to access BTVNanoCommissioning:", args.remoteRepo)
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
            jobdirs = [d for d in os.listdir(base_dir) if d.startswith("jobs_")]
            make_tarfile(
                "BTVNanoCommissioning.tar.gz",
                base_dir,
                exclude_dirs=["BTVNanoCommissioning.egg-info"] + jobdirs,
            )

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


Arguments = $(JOBNUM)

request_cpus = 1
request_memory = 2000
use_x509userproxy = true

+JobFlavour = "{jobqueue}"

Log        = {log_dir}/job.log_$(Cluster)
Output     = {log_dir}/job.out_$(Cluster)-$(Process)
Error      = {log_dir}/job.err_$(Cluster)-$(Process)

should_transfer_files   = YES
when_to_transfer_output = ON_EXIT_OR_EVICT
transfer_input_files    = {transfer_input_files}

Queue JOBNUM from {jobnum_file}
""".format(
        executable=f"{base_dir}/condor/execute.sh",
        jobqueue=args.jobqueue,
        log_dir=f"{base_dir}/{job_dir}/log",
        transfer_input_files=f"{base_dir}/{job_dir}/arguments.json,{base_dir}/{job_dir}/split_samples.json,{base_dir}/{job_dir}/jobnum_list.txt"
        + ("" if args.remoteRepo else f",{base_dir}/BTVNanoCommissioning.tar.gz"),
        jobnum_file=f"{base_dir}/{job_dir}/jobnum_list.txt",
    )
    with open(os.path.join(job_dir, "submit.jdl"), "w") as f:
        f.write(jdl_template)
    os.system(f"condor_submit {job_dir}/submit.jdl")
    # print(
    #     f"Setup completed. Now submit the condor jobs by:\n  condor_submit {job_dir}/submit.jdl"
    # )
