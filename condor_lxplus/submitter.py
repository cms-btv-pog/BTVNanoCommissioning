import os
import json
import shutil
import tarfile
import argparse
import subprocess


def make_tarfile(output_filename, source_dir, exclude_dirs=[]):
    with tarfile.open(output_filename, "w:gz") as tar:
        for root, dirs, files in os.walk(source_dir):
            dirs[:] = [
                d for d in dirs if d not in exclude_dirs
            ]  # Exclude specified directories
            for file in files:
                file_path = os.path.join(root, file)
                tar.add(file_path, arcname=os.path.relpath(file_path, source_dir))


# From https://github.com/PocketCoffea/PocketCoffea/blob/main/pocket_coffea/utils/network.py
def get_proxy_path() -> str:
    """
    Checks if the VOMS proxy exists and if it is valid
    for at least 1 hour.
    If it exists, returns the path of it"""
    try:
        subprocess.run("voms-proxy-info -exists -valid 0:20", shell=True, check=True)
    except subprocess.CalledProcessError:
        raise Exception(
            "VOMS proxy expirend or non-existing: please run `voms-proxy-init -voms cms -rfc --valid 168:0`"
        )

    # Now get the path of the certificate
    proxy = subprocess.check_output(
        "voms-proxy-info -path", shell=True, text=True
    ).strip()
    return proxy


def get_condor_submitter_parser(parser):
    parser.add_argument(
        "--jobName",
        help="Condor job name to make the job directory",
        required=True,
    )
    parser.add_argument(
        "--condorFileSize",
        type=int,
        default=1,
        help="Number of files to process per condor job",
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
            "2018_UL",
            "2017_UL",
            "2016preVFP_UL",
            "2016postVFP_UL",
            "CAMPAIGN_prompt_dataMC",
        ],
        help="Dataset campaign, change the corresponding correction files",
    )
    parser.add_argument(
        "--isSyst",
        default=False,
        type=str,
        choices=[False, "all", "weight_only", "JERC_split", "JP_MC"],
        help="Run with systematics, all, weights_only(no JERC uncertainties included),JERC_split, None",
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
    parser.add_argument("--submit", action="store_true", help="Submit condor jobs.")
    parser = get_condor_submitter_parser(parser)
    return parser


if __name__ == "__main__":
    parser = get_main_parser()
    args = parser.parse_args()
    print("Running with the following options:")
    print(args)

    current_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = current_dir.replace("/condor_lxplus", "")
    """
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
            make_tarfile(
                "BTVNanoCommissioning.tar.gz",
                base_dir,
                exclude_dirs=["jsonpog-integration", "BTVNanoCommissioning.egg-info"],
            )
    """
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

    # Handle voms proxy
    proxy_file = get_proxy_path()
    os.system(f"scp {proxy_file} proxy")
    print(f"Copied proxy file {proxy_file} to local directory.")

    # Find conda/mamba environment
    envpath = "/eos/home-m/milee/miniforge3/envs/btv_coffea/bin"
    pathvarlist = [i for i in os.environ["PATH"].split(":") if "envs/btv_coffea" in i]
    if len(pathvarlist) == 0:
        print(
            f"You did not source the btv_coffea conda/mamba environment. Proceed with the central conda environment:\n{envpath} ?"
        )
        response = input("(y/yes): ").strip().lower()
        if response == "y" or response == "yes":
            pass
        else:
            print(
                "First source the conda environment with which 'pip install -e .' was run. Then proceed again."
            )
            exit()
    else:
        envpath = pathvarlist[0]
        print(
            f"Using conda installation in\n{envpath}.\nThis has to be the conda installation with which 'pip install -e .' was run. If not, please source the correct environment and run again."
        )

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


Arguments = $(JOBNUM) {base_dir} {outputDir} {envpath}

request_cpus = 1
request_memory = 2000

+JobFlavour = "longlunch"

Log        = {log_dir}/job.log_$(Cluster)
Output     = {log_dir}/job.out_$(Cluster)-$(Process)
Error      = {log_dir}/job.err_$(Cluster)-$(Process)

should_transfer_files   = YES
when_to_transfer_output = ON_EXIT_OR_EVICT
transfer_input_files    = {transfer_input_files}
transfer_output_files   = .success

Queue JOBNUM from {jobnum_file}
""".format(
        executable=f"{base_dir}/condor_lxplus/execute.sh",
        base_dir=base_dir,
        outputDir=args.outputDir,
        envpath=envpath,
        log_dir=f"{base_dir}/{job_dir}/log",
        transfer_input_files=f"{base_dir}/{job_dir}/arguments.json,{base_dir}/{job_dir}/split_samples.json,{base_dir}/{job_dir}/jobnum_list.txt",
        jobnum_file=f"{base_dir}/{job_dir}/jobnum_list.txt",
    )
    with open(os.path.join(job_dir, "submit.jdl"), "w") as f:
        f.write(jdl_template)

    if "/eos/" in base_dir:
        print(
            "WARNING: Submitting from /eos. Log files will NOT be created.\nTo debug/retrieve logs, use `condor_transfer_data <job id>` after a job fails/terminates."
        )
        spool = "-spool"
    else:
        spool = ""
    if args.submit:
        os.system(f"condor_submit {spool} {job_dir}/submit.jdl")
    else:
        print(
            f"Setup completed. Now submit the condor jobs by:\n  condor_submit {spool} {job_dir}/submit.jdl"
        )
