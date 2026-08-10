import os, argparse, subprocess, shutil, glob, tarfile
from BTVNanoCommissioning.workflows import workflows
from BTVNanoCommissioning.utils.sample import predefined_sample
from BTVNanoCommissioning.utils.AK4_parameters import correction_config
import os, sys, inspect, shlex

current_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)

from runner import config_parser, scaleout_parser, debug_parser
from condor.submitter import get_condor_submitter_parser, validate_x509_proxy

CONDOR_ONLY_ARGS = {
    "--condorOutputBase",
    "--datafiles",
    "--mcfiles",
    "--jobqueue",
    "--jobName",
    "--outputDir",
    "--remoteRepo",
    "--nCPU",
    "-nCPU",
    "--condorFileSize",
    "-n",
    "--reuseTarball",
    "--rebuildTarball",
    "--submitRetries",
}


def requested_condor_only_options(argv):
    return sorted(
        {
            arg.split("=", 1)[0]
            for arg in argv
            if arg.split("=", 1)[0] in CONDOR_ONLY_ARGS
        }
    )


def get_condor_submission_name(args, workflow, sample_type):
    return f"condor_{workflow}_{sample_type}"


def get_condor_output_dir(args, workflow, sample_type):
    return os.path.join(
        args.condorOutputBase, f"{workflow}_{sample_type}{args.version}"
    )


def get_condor_file_size(args, sample_type, condor_file_size_override):
    if condor_file_size_override:
        return args.condorFileSize
    if sample_type == "data":
        return args.datafiles
    return args.mcfiles


def should_skip_mc_family(sample_type, mc_family):
    tokens = set(sample_type.upper().split("_"))
    if mc_family == "lo":
        return "NLO" in tokens
    if mc_family == "nlo":
        return "LO" in tokens and "NLO" not in tokens
    return False


def get_condor_job_dir(args, workflow, sample_type):
    return f"jobs_{get_condor_submission_name(args, workflow, sample_type)}_{args.campaign}"


RUNNER_SKIP_KEYS = {
    "workflow",
    "json",
    "campaign",
    "year",
    "scheme",
    "DAS_campaign",
    "version",
    "condor",
    "local",
    "debug",
    "limit_MC",
    "limit_MC_Wc",
    "validate_workflow",
    "mc",
    "jobName",
    "outputDir",
    "condorOutputBase",
    "datafiles",
    "mcfiles",
    "condorFileSize",
    "remoteRepo",
    "jobqueue",
    "nCPU",
    "reuseTarball",
    "rebuildTarball",
    "submitRetries",
    "response",
}


CONDOR_ARG_KEYS = {
    "isSyst",
    "isArray",
    "noHist",
    "overwrite",
    "only",
    "voms",
    "chunk",
    "skipbadfiles",
    "submitRetries",
}


def is_running_in_ci():
    """Check if running in GitLab CI environment"""
    return os.environ.get("GITLAB_CI") == "true"


def should_refresh_dataset(json_file, max_age_minutes=10):
    """Check if the dataset JSON file needs refreshing based on its age"""
    import os
    import time

    if not os.path.exists(json_file):
        return True  # File doesn't exist, must fetch

    file_mtime = os.path.getmtime(json_file)
    current_time = time.time()
    age_in_minutes = (current_time - file_mtime) / 60

    if age_in_minutes > max_age_minutes:
        print(
            f"⚠️ Dataset file {json_file} is {age_in_minutes:.1f} minutes old, needs refreshing"
        )
        return True
    return False


def workflow_sample_json_path(args, workflow_tag, sample_type):
    return f"metadata/{args.campaign}/{sample_type}_{args.campaign}_{args.year}_{workflow_tag}.json"


def should_refresh_workflow_datasets(
    args,
    workflow_tag,
    allowed_sample_types,
    response_opt,
    overwrite,
    sample_types=None,
):
    """Check whether any sample JSON needed by a workflow is missing or stale."""
    refresh_required = bool(args.overwrite)
    checked_any = False

    candidate_sample_types = (
        sample_types
        if sample_types is not None
        else predefined_sample[workflow_tag].keys()
    )

    for sample_type in candidate_sample_types:
        if allowed_sample_types is not None and sample_type not in allowed_sample_types:
            continue
        if should_skip_mc_family(sample_type, args.mc):
            continue
        if (sample_type != "data" and sample_type != "MC") and (
            args.scheme == "Validation" or args.validate_workflow
        ):
            continue

        checked_any = True
        json_file = workflow_sample_json_path(args, workflow_tag, sample_type)
        if should_refresh_dataset(json_file):
            refresh_required = True

    if not checked_any:
        return False

    if refresh_required:
        fetch_cmd = (
            f"python scripts/fetch.py -c {args.campaign} --from_workflow {workflow_tag} "
            f"--DAS_campaign {args.DAS_campaign} --year {args.year} {overwrite} "
            f"--skipvalidation --overwrite --executor futures -j {args.workers} {response_opt}"
        )
        subprocess.run(shlex.split(fetch_cmd), check=True)

    return refresh_required


def get_condor_sample_plan(args, wf, workflow_tag, allowed_sample_types):
    plan = []
    skipped_existing = []
    for sample_type in predefined_sample[workflow_tag].keys():
        if allowed_sample_types is not None and sample_type not in allowed_sample_types:
            continue
        if should_skip_mc_family(sample_type, args.mc):
            continue
        if (sample_type != "data" and sample_type != "MC") and (
            args.scheme == "Validation" or args.validate_workflow
        ):
            continue

        job_dir = get_condor_job_dir(args, wf, sample_type)
        if os.path.exists(job_dir):
            skipped_existing.append((sample_type, job_dir))
            continue
        plan.append(sample_type)

    return plan, skipped_existing


def pretty_rule(title):
    width = 78
    title = f" {title} "
    side = max((width - len(title)) // 2, 3)
    return f"{'═' * side}{title}{'═' * side}"


def pretty_status(emoji, label, message):
    print(f"{emoji} {label}: {message}")


def validate_reusable_tarball(path="BTVNanoCommissioning.tar.gz"):
    if not os.path.isfile(path):
        raise RuntimeError(f"Expected reusable tarball does not exist: {path}")
    if not os.access(path, os.R_OK):
        raise RuntimeError(f"Expected reusable tarball is not readable: {path}")
    try:
        with tarfile.open(path, "r:gz") as archive:
            next(iter(archive), None)
    except (OSError, tarfile.TarError) as exc:
        raise RuntimeError(f"Reusable tarball is invalid: {path}") from exc


def run_local_smoke_test(args, wf, json_path, sample_type):
    test_outputdir = os.path.join(".suball_test_outputs", f"{wf}_{sample_type}")
    shutil.rmtree(test_outputdir, ignore_errors=True)
    os.makedirs(os.path.dirname(test_outputdir), exist_ok=True)

    command = [
        sys.executable,
        "runner.py",
        "--wf",
        wf,
        "--json",
        json_path,
        "--campaign",
        args.campaign,
        "--year",
        str(args.year),
        "--executor",
        "iterative",
        "--limit",
        "1",
        "--max",
        "1",
        "--outputdir",
        test_outputdir,
    ]
    if args.isArray:
        command.append("--isArray")
    if args.skipbadfiles:
        command.append("--skipbadfiles")
    if args.only is not None:
        command.extend(["--only", args.only])
    if args.isSyst != "False":
        command.extend(["--isSyst", args.isSyst])

    print()
    print(pretty_rule(f" Smoke test for {wf} / {sample_type} "))
    print(f"🔎 Local command: {' '.join(shlex.quote(part) for part in command)}")
    result = subprocess.run(command)

    coffea_files = glob.glob(
        os.path.join(test_outputdir, "**", "*.coffea"), recursive=True
    )
    root_files = glob.glob(os.path.join(test_outputdir, "**", "*.root"), recursive=True)

    if result.returncode != 0:
        print(pretty_rule(f" Smoke test failed for {wf} / {sample_type} "))
        pretty_status("❌", "Exit code", str(result.returncode))
        print(
            f"❌ Submission check failed for {wf}; aborting before Condor submission."
        )
        shutil.rmtree(test_outputdir, ignore_errors=True)
        sys.exit(1)
    pretty_status("✅", "Exit code", "0")
    if len(coffea_files) == 0:
        print(pretty_rule(f" Smoke test failed for {wf} / {sample_type} "))
        print(f"❌ Output .coffea exists: no")
        print(
            f"❌ Submission check failed for {wf}; aborting before Condor submission."
        )
        shutil.rmtree(test_outputdir, ignore_errors=True)
        sys.exit(1)
    pretty_status("✅", "Output .coffea exists", "yes")
    if args.isArray and len(root_files) == 0:
        print(pretty_rule(f" Smoke test failed for {wf} / {sample_type} "))
        print(f"❌ Output .root exists: no")
        print(
            f"❌ Submission check failed for {wf}; aborting before Condor submission."
        )
        shutil.rmtree(test_outputdir, ignore_errors=True)
        sys.exit(1)
    if args.isArray:
        pretty_status("✅", "Output .root exists", "yes")
    else:
        pretty_status("ℹ️", "Output .root exists", "not requested (--isArray not set)")

    shutil.rmtree(test_outputdir, ignore_errors=True)
    parent = os.path.dirname(test_outputdir)
    if os.path.isdir(parent) and not os.listdir(parent):
        os.rmdir(parent)
    print(pretty_rule(f" Smoke test passed for {wf} / {sample_type} "))
    print()


# Get lumi
def get_lumi_from_web(year):
    import requests
    import re

    year = str(year)
    # Define the URL of the directory
    url = (
        f"https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions{year[2:]}/"
    )

    # Send a request to fetch the HTML content of the webpage
    response = requests.get(url)
    html_content = response.text

    # Use regex to find all href links that contain 'Golden.json' but do not contain 'era'
    # Ensures it only captures the URL part within href="..." and not any other content.
    goldenjson_files = re.findall(r'href="([^"]*Golden\.json[^"]*)"', html_content)

    # Filter out any matches that contain 'era' in the filename
    goldenjson_files = [file for file in goldenjson_files if "era" not in file]

    # If there are any such files, find the latest one (assuming the files are sorted lexicographically)
    if goldenjson_files:
        latest_file = sorted(goldenjson_files)[
            -1
        ]  # Assuming lexicographical sorting works for the dates
        os.system(f"wget {url}/{latest_file}")
        os.system(f"mv {latest_file} src/BTVNanoCommissioning/data/DC/.")
        return latest_file
    else:
        raise (
            f"No files for Year{year} containing 'Golden.json' (excluding 'era') were found."
        )


### Manage workflow in one script
# Local example:
# python scripts/suball.py --scheme default_comissioning --campaign Summer23 --DAS_campaign "*Run2023D*Sep2023*,*Run3Summer23BPixNanoAODv12-130X*" --year 2023
# Condor example:
# python scripts/suball.py --condor --condorOutputBase /eos/path --scheme CFM --campaign Summer23 --DAS_campaign "*Run2023D*,*Run3Summer23*" --year 2023
# prerequest a new campaign should create a entry in AK4_parameters.py
#############     #############      ##########     ########
#  dataset  #     #   Run     #      #  Dump  #     #      #
#           # ==> #   coffea  #  ==> #        # ==> # Plot #
#  creation #     # processor #      #  Lumi  #     #      #
#############     #############      ##########     ########
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Mastering workflow submission")
    parser = config_parser(parser)
    paser = scaleout_parser(parser)
    paser = debug_parser(parser)
    parser = get_condor_submitter_parser(parser, require_job_args=False)
    parser.add_argument(
        "--condor",
        action="store_true",
        help=(
            "Submit workflows through HTCondor. Requires --condorOutputBase. "
            "Without this flag, suball runs locally and rejects Condor-only options."
        ),
    )
    parser.add_argument(
        "-sc",
        "--scheme",
        default="Validation",
        choices=list(workflows.keys())
        + ["Validation", "Validation_tt", "SF", "default_comissioning", "CFM"],
        help="Choose the function for dump luminosity(`lumi`)/failed files(`failed`) into json",
    )

    parser.add_argument(
        "-dc",
        "--DAS_campaign",
        required=True,
        help="Input the campaign name for DAS to search appropriate campaigns, use in dataset construction , please do `data_camapgin,mc_campaign` split by `,`, e.g. `*Run2023D*Sep2023*,*Run3Summer23BPixNanoAODv12-130X*` ",
    )
    parser.add_argument("-v", "--version", default="", help="version postfix")
    parser.add_argument(
        "--local",
        action="store_true",
        help="not transfered to https://btvweb.web.cern.ch/Commissioning/dataMC/",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Run local debug test with small set of dataset with iterative executor",
    )
    parser.add_argument(
        "--limit_MC_Wc",
        action="store_true",
        help="Limit MC samples to 100 files regardless of workflow/scheme",
    )
    parser.add_argument(
        "--limit_MC",
        action="store_true",
        help="Limit MC samples to 50 files regardless of workflow/scheme",
    )
    parser.add_argument(
        "--validate_workflow",
        "-vw",
        action="store_true",
        help="Run only data and MC samples for the workflow, skip minor MC samples",
    )
    parser.add_argument(
        "--mc",
        type=lambda value: value.lower(),
        default="nlo",
        choices=["lo", "nlo"],
        help="Choose which MC family to keep when a workflow splits MC into LO/NLO buckets (default: %(default)s).",
    )
    parser.add_argument(
        "--condorOutputBase",
        default=None,
        help="Required base output directory for --condor submissions.",
    )
    parser.add_argument(
        "--datafiles",
        type=int,
        default=15,
        help="Number of data files to put into each Condor job (default: %(default)s).",
    )
    parser.add_argument(
        "--mcfiles",
        type=int,
        default=5,
        help="Number of non-data files to put into each Condor job (default: %(default)s).",
    )
    parser.add_argument(
        "--response",
        default=None,
        help="Forward a non-interactive answer to fetch.py dataset selection prompts.",
    )

    raw_argv = sys.argv[1:]
    args = parser.parse_args()
    use_condor = args.condor
    requested_condor_options = requested_condor_only_options(raw_argv)
    condor_file_size_override = any(
        arg.split("=", 1)[0] in {"--condorFileSize", "-n"} for arg in raw_argv
    )

    if not use_condor and requested_condor_options:
        parser.error(
            "Condor-only option(s) require --condor: "
            + ", ".join(requested_condor_options)
        )

    if use_condor and not args.condorOutputBase:
        parser.error(
            "--condor requires --condorOutputBase.\n"
            "Example:\n"
            f"  {sys.executable} scripts/suball.py --condor --scheme CFM "
            f"--campaign {args.campaign} --year {args.year} "
            f"--DAS_campaign '{args.DAS_campaign}' --jobqueue workday "
            "--condorOutputBase /eos/path --isArray --skipbadfiles"
        )

    if use_condor:
        validate_x509_proxy()

    # summarize diffeerent group for study
    scheme = {
        # scale factor workflows
        "SF": ["BTA_ttbar", "BTA_addPFMuons"],
        # Use for prompt data MC checks for analysis
        "Validation": ["ttdilep_sf", "ctag_Wc_sf"],
        "Validation_tt": ["ttdilep_sf"],
        "Validation_ctag": ["ctag_Wc_sf"],
        # commissioning workflows
        "default_comissioning": [
            "ttdilep_sf",
            "ttsemilep_sf",
            "ctag_Wc_sf",
            "ctag_DY_sf",
            "QCD_sf",
            "QCD_mu_sf",
        ],
        # CFM submission bundle from submit_DYele.sh
        "CFM": [
            "ctag_Wc_noMuVeto_sf",
            "ectag_Wc_sf",
            "ctag_DY_sf",
            "ectag_DY_sf",
            "ctag_ttsemilep_noMuVeto_sf",
            "ectag_ttsemilep_sf",
        ],
    }
    if args.scheme in workflows.keys():
        scheme[args.scheme] = [args.scheme]
        # scheme["test"] = [args.scheme]
        # args.scheme = "test"

    cfm_workflow_aliases = {
        # These two workflows share the same dataset tags and JSON naming in submit_DYele.sh
        "ctag_Wc_noMuVeto_sf": "ctag_Wc_sf",
        "ctag_ttsemilep_noMuVeto_sf": "ctag_ttsemilep_sf",
    }
    cfm_sample_types = {
        "ctag_Wc_noMuVeto_sf": {"data", "MC", "MC_LO", "minor_MC"},
        "ectag_Wc_sf": {"data", "MC", "MC_LO", "minor_MC"},
        "ctag_DY_sf": {"data", "MC_LO", "minor_MC"},
        "ectag_DY_sf": {"data", "MC_LO", "minor_MC"},
        "ctag_ttsemilep_noMuVeto_sf": {"data", "MC", "MC_LO", "minor_MC"},
        "ectag_ttsemilep_sf": {"data", "MC", "MC_LO", "minor_MC"},
    }

    if args.scheme == "CFM":
        missing_cfm_flags = []
        if not args.isArray:
            missing_cfm_flags.append("--isArray")
        if args.isSyst != "all_withJESTotal":
            missing_cfm_flags.append("--isSyst all_withJESTotal")
        if missing_cfm_flags:
            print(
                "⚠️ CFM usually expects " + ", ".join(missing_cfm_flags) + " to be set."
            )
            response = (
                input("Proceed with the CFM bundle anyway? [y/N]: ").strip().lower()
            )
            if response not in {"y", "yes"}:
                print("Aborting CFM submission at user request.")
                sys.exit(1)

    # Check lumiMask exists and replace the Validation
    campaign_config = correction_config[args.campaign]
    campaign_default_config = campaign_config.get("default", campaign_config)
    input_lumi_json = campaign_default_config["DC"]
    if args.campaign != "prompt_dataMC" and not os.path.exists(
        f"src/BTVNanoCommissioning/data/DC/{input_lumi_json}"
    ):
        raise f"src/BTVNanoCommissioning/data/DC/{input_lumi_json} not exist"

    if (
        args.campaign == "prompt_dataMC"
        and campaign_default_config["DC"] == "$PROMPT_DATAMC"
    ):
        input_lumi_json = get_lumi_from_web(args.year)
        os.system(
            f"sed -i 's/$PROMPT_DATAMC/{input_lumi_json}/g' src/BTVNanoCommissioning/utils/AK4_parameters.py"
        )
        print(f"======>{input_lumi_json} is used for {args.year}")

    condor_submit_count = 0
    reusable_tarball_ready = False
    for wf in scheme[args.scheme]:
        workflow_tag = cfm_workflow_aliases.get(wf, wf)
        allowed_sample_types = (
            cfm_sample_types.get(wf) if args.scheme == "CFM" else None
        )
        if args.validate_workflow:
            print(
                f"ℹ️ Running workflow '{wf}' in validation mode (only data and MC samples)"
            )
        if args.debug:
            print(f"======{wf} in {args.scheme}=====")
        overwrite = "--overwrite" if args.overwrite else ""
        response_opt = (
            f"--response {shlex.quote(str(args.response))}"
            if args.response is not None
            else ""
        )
        sample_types_to_submit = list(predefined_sample[workflow_tag].keys())
        if use_condor:
            sample_types_to_submit, skipped_existing_jobs = get_condor_sample_plan(
                args, wf, workflow_tag, allowed_sample_types
            )
            for sample_type, job_dir in skipped_existing_jobs:
                print(
                    f"⚠️ Condor job directory already exists for {wf}/{sample_type}: {job_dir}. "
                    "Skipping this job and moving on."
                )
            if not sample_types_to_submit:
                print(
                    f"⚠️ All condor job directories already exist for workflow '{wf}'. "
                    "Skipping fetch and submission for this workflow."
                )
                continue
        ## create or refresh datasets once before any submissions for this workflow
        if args.debug:
            print(
                f"Checking workflow datasets for {wf}: python scripts/fetch.py -c {args.campaign} --from_workflow {workflow_tag} --DAS_campaign {args.DAS_campaign} --year {args.year} {overwrite} --skipvalidation --overwrite --executor futures -j {args.workers} {response_opt}"
            )
        should_refresh_workflow_datasets(
            args,
            workflow_tag,
            allowed_sample_types,
            response_opt,
            overwrite,
            sample_types=sample_types_to_submit,
        )
        if args.debug:
            os.system(f"ls metadata/{args.campaign}/*.json")

        ## Run the workflows
        for types in sample_types_to_submit:
            if allowed_sample_types is not None and types not in allowed_sample_types:
                if args.debug:
                    print(f"⚠️ Skipping sample type '{types}' for CFM workflow '{wf}'")
                continue

            if should_skip_mc_family(types, args.mc):
                if args.debug:
                    print(
                        f"⚠️ Skipping sample type '{types}' because --mc {args.mc} was selected"
                    )
                continue

            if (types != "data" and types != "MC") and (
                args.scheme == "Validation" or args.validate_workflow
            ):
                print(f"⚠️ Skipping minor sample type '{types}' due to validation mode")
                continue
            print(
                f"hists_{wf}_{types}_{args.campaign}_{args.year}_{wf}/hists_{wf}_{types}_{args.campaign}_{args.year}_{wf}.coffea"
            )
            if (
                not os.path.exists(
                    f"hists_{wf}_{types}_{args.campaign}_{args.year}_{wf}/hists_{wf}_{types}_{args.campaign}_{args.year}_{wf}.coffea"
                )
                or args.overwrite
            ):
                if not os.path.exists(
                    f"metadata/{args.campaign}/{types}_{args.campaign}_{args.year}_{workflow_tag}.json"
                ):
                    raise Exception(
                        f"metadata/{args.campaign}/{types}_{args.campaign}_{args.year}_{workflow_tag}.json not exist"
                    )
                json_file = f"metadata/{args.campaign}/{types}_{args.campaign}_{args.year}_{workflow_tag}.json"

                if use_condor and condor_submit_count == 0:
                    run_local_smoke_test(args, wf, json_file, types)

                json_path = f"metadata/{args.campaign}/{types}_{args.campaign}_{args.year}_{workflow_tag}.json"
                base_command = (
                    [sys.executable, "condor/submitter.py"]
                    if use_condor
                    else [sys.executable, "runner.py"]
                )
                command = base_command + [
                    "--wf",
                    wf,
                    "--json",
                    json_path,
                    "--campaign",
                    args.campaign,
                    "--year",
                    str(args.year),
                ]
                limit_added = False  # Track if we've already added a limit flag

                for key, value in vars(args).items():
                    if use_condor:
                        if key not in CONDOR_ARG_KEYS:
                            continue
                    elif key in RUNNER_SKIP_KEYS:
                        continue

                    # Handle boolean flags
                    if key in [
                        "isArray",
                        "noHist",
                        "overwrite",
                        "validate",
                        "skipbadfiles",
                    ]:
                        if value == True:
                            command.append(f"--{key}")
                    elif value is not None:
                        if key == "limit":
                            command.extend([f"--{key}", str(value)])
                            limit_added = True
                        else:
                            command.extend([f"--{key}", str(value)])

                # Add limit for MC validation if not already present
                if types == "MC" and not limit_added:
                    # Apply limit if it's Validation or the limit_MC flag is set
                    if (
                        "Validation" == args.scheme
                        or "Validation_tt" == args.scheme
                        or args.limit_MC
                    ):
                        if not use_condor:
                            command.extend(["--limit", "50"])
                        limit_added = True
                        print(f"⚠️ Running with 50 files limit for MC samples")
                    elif args.limit_MC_Wc:
                        if not use_condor:
                            command.extend(["--limit", "100"])
                        limit_added = True
                        print(f"⚠️ Running with 100 files limit for MC samples")

                if use_condor:
                    if condor_submit_count == 0 and args.remoteRepo is None:
                        command.append("--rebuildTarball")
                    elif condor_submit_count > 0:
                        if args.remoteRepo is None:
                            if not reusable_tarball_ready:
                                raise RuntimeError(
                                    "Refusing --reuseTarball before the first successful "
                                    "submission created and validated the tarball"
                                )
                            validate_reusable_tarball()
                            command.append("--reuseTarball")
                    command.extend(
                        [
                            "--jobName",
                            get_condor_submission_name(args, wf, types),
                            "--outputDir",
                            get_condor_output_dir(args, wf, types),
                            "--condorFileSize",
                            str(
                                get_condor_file_size(
                                    args, types, condor_file_size_override
                                )
                            ),
                        ]
                    )
                    if args.remoteRepo is not None:
                        command.extend(["--remoteRepo", args.remoteRepo])
                    if args.jobqueue is not None:
                        command.extend(["--jobqueue", args.jobqueue])
                    if args.nCPU is not None:
                        command.extend(["--nCPU", str(args.nCPU)])

                runner_config = shlex.join(command)
                if args.debug:
                    print(f"run the workflow: {runner_config}")
                subprocess.run(command, check=True)
                if use_condor:
                    condor_submit_count += 1
                    if args.remoteRepo is None and condor_submit_count == 1:
                        validate_reusable_tarball()
                        reusable_tarball_ready = True

                with open(
                    f"config_{args.year}_{args.campaign}_{args.scheme}_{args.version}.txt",
                    "w",
                ) as config_list:
                    config_list.write(runner_config)

        if args.debug:
            print(f"workflow is finished for {wf}!")

        if use_condor:
            print(
                f"Condor submissions finished for {wf}; skipping local lumi and plotting."
            )
            continue

        if is_running_in_ci():
            import numpy as np
            import awkward as ak
            import json

            print(f"Running in CI environment - creating lumi JSON for {wf}")

            # Extract luminosity JSON for use in CI
            coffea_file = f"hists_{wf}_data_{args.campaign}_{args.year}_{wf}/hists_{wf}_data_{args.campaign}_{args.year}_{wf}.coffea"
            if os.path.exists(coffea_file):
                try:
                    from coffea.util import load

                    # Create lumi_bril directory if needed
                    os.makedirs("lumi_bril", exist_ok=True)

                    # Create the same structure as dump_processed.py expects
                    # The key is the coffea file path, value is the loaded coffea file
                    output = {coffea_file: load(coffea_file)}

                    # Now follow the exact same logic as in dump_lumi
                    lumi, run = [], []
                    for m in output.keys():  # m is the coffea file path
                        for f in output[
                            m
                        ].keys():  # f is the dataset key inside the coffea file
                            if "lumi" in output[m][f] and "run" in output[m][f]:
                                try:
                                    lumi.extend(output[m][f]["lumi"].value)
                                    run.extend(output[m][f]["run"].value)
                                except Exception as e:
                                    print(f"  Error extracting run/lumi from {f}: {e}")

                    if len(run) > 0 and len(lumi) > 0:
                        print(f"Found {len(run)} run/lumi pairs, creating JSON")

                        # Sort runs and keep lumisections matched
                        run, lumi = np.array(run), np.array(lumi)
                        sorted_indices = np.lexsort(
                            (lumi, run)
                        )  # Sort by run first, then lumi
                        run = run[sorted_indices]
                        lumi = lumi[sorted_indices]

                        # Create dictionary with ls values for each run
                        dicts = {}
                        for r in np.unique(run):
                            # Make sure to cast to int to avoid string keys
                            dicts[str(int(r))] = lumi[run == r]

                        # Convert to format for brilcalc (exactly as in dump_lumi)
                        for r in dicts.keys():
                            ar = ak.singletons(ak.Array(dicts[r]))
                            ars = ak.concatenate([ar, ar], axis=-1)
                            dicts[r] = ak.values_astype(ars, int).tolist()

                        # Save JSON file for brilcalc
                        json_path = f"lumi_bril/{wf}_bril_lumi.json"
                        with open(json_path, "w") as outfile:
                            json.dump(dicts, outfile, indent=2)

                        print(f"Created luminosity JSON file for {wf} at {json_path}")
                    else:
                        print(f"No run/lumi information found for {wf}")
                        # Create an empty file so we know we tried
                        with open(f"lumi_bril/{wf}_bril_empty.json", "w") as outfile:
                            json.dump({}, outfile)
                except Exception as e:
                    print(f"ERROR creating luminosity JSON for {wf}: {e}")
                    import traceback

                    traceback.print_exc()

            # Skip the rest of luminosity calculation and plotting in CI
            print(f"Skipping luminosity calculation and plotting for {wf}")
            continue
        # Get luminosity
        if (
            os.path.exists(
                f"hists_{wf}_data_{args.campaign}_{args.year}_{wf}/hists_{wf}_data_{args.campaign}_{args.year}_{wf}.coffea"
            )
            or args.overwrite
        ):
            if args.debug:
                print(
                    f"Get the luminosity from hists_{wf}_data_{args.campaign}_{args.year}_{wf}/hists_{wf}_data_{args.campaign}_{args.year}_{wf}.coffea"
                )
            if not os.path.exists(
                f"hists_{wf}_data_{args.campaign}_{args.year}_{wf}/hists_{wf}_data_{args.campaign}_{args.year}_{wf}.coffea"
            ):
                raise Exception(
                    f"hists_{wf}_data_{args.campaign}_{args.year}_{wf}/hists_{wf}_data_{args.campaign}_{args.year}_{wf}.coffea not exist"
                )
            lumi = os.popen(
                f"python scripts/dump_processed.py -t all -c hists_{wf}_data_{args.campaign}_{args.year}_{wf}/hists_{wf}_data_{args.campaign}_{args.year}_{wf}.coffea --json metadata/{args.campaign}/data_{args.campaign}_{args.year}_{wf}.json -n {args.campaign}_{args.year}_{wf}"
            ).read()
            print(lumi)
            lumi = int(
                round(
                    float(
                        lumi[
                            lumi.find("Luminosity in pb:")
                            + 18 : lumi.find("===>Dump Failed Files")
                            - 1
                        ]
                    ),
                    0,
                )
            )
            if os.path.exists(
                f"hists_{wf}_MC_{args.campaign}_{args.year}_{wf}/hists_{wf}_MC_{args.campaign}_{args.year}_{wf}.coffea"
            ) and os.path.exists(
                f"hists_{wf}_data_{args.campaign}_{args.year}_{wf}/hists_{wf}_data_{args.campaign}_{args.year}_{wf}.coffea"
            ):
                if args.debug:
                    print(f"Plot the dataMC for {wf}")
                os.system(
                    f'python scripts/plotdataMC.py -i "hists_{wf}_*_{args.campaign}_{args.year}_{wf}/hists_{wf}_*_{args.campaign}_{args.year}_{wf}.coffea" --lumi {lumi} -p {wf} -v all --ext {args.campaign}_{args.year}{args.version}'
                )
                ## Inspired from Uttiya, create remote directory
                # https://github.com/cms-btv-pog/BTVNanoCommissioning/blob/14e654feeb4b4d738ee43ab913efb343ea65fd1d/scripts/submit/createremotedir.sh
                # create remote direcotry
                if args.debug:
                    print(f"Upload plots&coffea to eos: {wf}")
                if not args.local:
                    os.system(f"mkdir -p {args.campaign}{args.version}/{wf}")
                    os.system(f"cp scripts/index.php {args.campaign}{args.version}/.")
                    os.system(
                        f"xrdcp -r  {args.campaign}{args.version}/ root://eosuser.cern.ch//eos/user/b/btvweb/www/Commissioning/dataMC/{args.scheme}/."
                    )
                    os.system(f"cp scripts/index.php {args.campaign}/{wf}/.")
                    os.system(
                        f"cp hists_{wf}_*_{args.campaign}_{args.year}_{wf}/*.coffea {args.campaign}/{wf}/."
                    )
                    os.system(
                        f"cp plot/{wf}_{args.campaign}_{args.year}{args.version}/* {args.campaign}{args.version}/{wf}/."
                    )
                    overwrite = "-f " if args.overwrite else ""
                    os.system(
                        f"xrdcp -r -p {overwrite} {args.campaign}{args.version}/{wf} root://eosuser.cern.ch//eos/user/b/btvweb/www/Commissioning/dataMC/{args.scheme}/{args.campaign}{args.version}/."
                    )
            else:
                raise Exception(
                    f"No input coffea hists_{wf}_data_{args.campaign}_{args.year}_{wf}/hists_{wf}_data_{args.campaign}_{args.year}_{wf}.coffea or hists_{wf}_MC_{args.campaign}_{args.year}_{wf}/hists_{wf}_MC_{args.campaign}_{args.year}_{wf}.coffea"
                )
    # revert prompt_dataMC lumimask
    if args.campaign == "prompt_dataMC":
        os.system(
            f"sed -i 's/{input_lumi_json}/$PROMPT_DATAMC/g' src/BTVNanoCommissioning/utils/AK4_parameters.py"
        )
