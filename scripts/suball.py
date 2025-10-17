import os, argparse
from BTVNanoCommissioning.workflows import workflows
from BTVNanoCommissioning.utils.sample import predefined_sample
from BTVNanoCommissioning.utils.AK4_parameters import correction_config
import os, sys, inspect

current_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)

from runner import config_parser, scaleout_parser, debug_parser


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
# EXAMPLE: python scripts/suball.py --scheme default_comissioning --campaign Summer23  --DAS_campaign "*Run2023D*Sep2023*,*Run3Summer23BPixNanoAODv12-130X*" --year 2023
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
    parser.add_argument(
        "-sc",
        "--scheme",
        default="Validation",
        choices=list(workflows.keys())
        + ["Validation", "Validation_tt", "SF", "default_comissioning"],
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

    args = parser.parse_args()
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
    }
    if args.scheme in workflows.keys():
        scheme[args.scheme] = [args.scheme]
        # scheme["test"] = [args.scheme]
        # args.scheme = "test"

    # Check lumiMask exists and replace the Validation
    input_lumi_json = correction_config[args.campaign]["DC"]
    if args.campaign != "prompt_dataMC" and not os.path.exists(
        f"src/BTVNanoCommissioning/data/DC/{input_lumi_json}"
    ):
        raise f"src/BTVNanoCommissioning/data/DC/{input_lumi_json} not exist"

    if (
        args.campaign == "prompt_dataMC"
        and correction_config[args.campaign]["DC"] == "$PROMPT_DATAMC"
    ):
        input_lumi_json = get_lumi_from_web(args.year)
        os.system(
            f"sed -i 's/$PROMPT_DATAMC/{input_lumi_json}/g' src/BTVNanoCommissioning/utils/AK4_parameters.py"
        )
        print(f"======>{input_lumi_json} is used for {args.year}")

    for wf in scheme[args.scheme]:
        if args.validate_workflow:
            print(
                f"ℹ️ Running workflow '{wf}' in validation mode (only data and MC samples)"
            )
        if args.debug:
            print(f"======{wf} in {args.scheme}=====")
        overwrite = "--overwrite" if args.overwrite else ""
        ## creating dataset
        if (
            not os.path.exists(
                f"metadata/{args.campaign}/MC_{args.campaign}_{args.year}_{wf}.json"
            )
            or args.overwrite
        ):
            if args.debug:
                print(
                    f"Creating MC dataset: python scripts/fetch.py -c {args.campaign} --from_workflow {wf} --DAS_campaign {args.DAS_campaign} --year {args.year} {overwrite} --skipvalidation --executor futures"
                )

            os.system(
                f"python scripts/fetch.py -c {args.campaign} --from_workflow {wf} --DAS_campaign {args.DAS_campaign} --year {args.year} {overwrite} --skipvalidation --executor futures"
            )
            if args.debug:
                os.system(f"ls metadata/{args.campaign}/*.json")

        ## Run the workflows
        for types in predefined_sample[wf].keys():

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
                    f"metadata/{args.campaign}/{types}_{args.campaign}_{args.year}_{wf}.json"
                ):
                    raise Exception(
                        f"metadata/{args.campaign}/{types}_{args.campaign}_{args.year}_{wf}.json not exist"
                    )
                json_file = f"metadata/{args.campaign}/{types}_{args.campaign}_{args.year}_{wf}.json"

                # Check if dataset needs refreshing because of age
                if should_refresh_dataset(json_file):
                    print(f"Refreshing dataset: {json_file}")
                    fetch_cmd = f"python scripts/fetch.py -c {args.campaign} --from_workflow {wf} --DAS_campaign {args.DAS_campaign} --year {args.year} {overwrite} --skipvalidation --overwrite --executor futures"
                    os.system(fetch_cmd)

                runner_config_required = f"python runner.py --wf {wf} --json metadata/{args.campaign}/{types}_{args.campaign}_{args.year}_{wf}.json {overwrite} --campaign {args.campaign} --year {args.year}"
                runner_config = ""
                limit_added = False  # Track if we've already added a limit flag

                for key, value in vars(args).items():
                    # Required info or things not relevant for runner, skip
                    # Skip keys that don't apply to runner
                    if key in [
                        "workflow",
                        "json",
                        "campaign",
                        "year",
                        "scheme",
                        "DAS_campaign",
                        "version",
                        "local",
                        "debug",
                        "limit_MC",
                        "limit_MC_Wc",
                        "validate_workflow",
                    ]:
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
                            runner_config += f" --{key}"
                    elif value is not None:
                        if key == "limit":
                            runner_config += f" --{key}={value}"
                            limit_added = True
                        else:
                            runner_config += f" --{key}={value}"

                # Add limit for MC validation if not already present
                if types == "MC" and not limit_added:
                    # Apply limit if it's Validation or the limit_MC flag is set
                    if (
                        "Validation" == args.scheme
                        or "Validation_tt" == args.scheme
                        or args.limit_MC
                    ):
                        runner_config += " --limit 50"
                        limit_added = True
                        print(f"⚠️ Running with 50 files limit for MC samples")
                    elif args.limit_MC_Wc:
                        runner_config += " --limit 100"
                        limit_added = True
                        print(f"⚠️ Running with 100 files limit for MC samples")
                runner_config = runner_config_required + runner_config
                if args.debug:
                    print(f"run the workflow: {runner_config}")
                os.system(runner_config)

                with open(
                    f"config_{args.year}_{args.campaign}_{args.scheme}_{args.version}.txt",
                    "w",
                ) as config_list:
                    config_list.write(runner_config)

        if args.debug:
            print(f"workflow is finished for {wf}!")

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
