import os, argparse
from BTVNanoCommissioning.workflows import workflows
from BTVNanoCommissioning.utils.sample import predefined_sample
from BTVNanoCommissioning.utils.AK4_parameters import correction_config
import os, sys, inspect

current_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)

from runner import config_parser, scaleout_parser, debug_parser

# Add this function to your suball.py script
def check_xrootd_site_availability(site_url):
    """Check if a XRootD site is responding properly"""
    import subprocess
    import re
    import os
    
    # Check if xrdfs is available first
    if os.system("which xrdfs >/dev/null 2>&1") != 0:
        # Try to install xrootd-clients if not available
        print("xrdfs command not found, attempting to install xrootd-clients...")
        os.system("apt-get update -qq && apt-get install -y xrootd-client >/dev/null 2>&1")
        
        # Check again
        if os.system("which xrdfs >/dev/null 2>&1") != 0:
            print("Failed to install xrootd-client. Using global redirector as fallback.")
            # Just assume cms-xrd-global.cern.ch is working
            return site_url.startswith("root://cms-xrd-global.cern.ch")
    
    
    # Extract the base URL without protocol and path
    match = re.search(r'root://([^/]+)', site_url)
    if not match:
        print(f"Invalid XRootD URL format: {site_url}")
        return False
        
    base_url = match.group(1)
    check_cmd = f"xrdfs {base_url} ping"
    
    try:
        result = subprocess.run(check_cmd, shell=True, capture_output=True, timeout=10)
        if result.returncode == 0:
            print(f"✅ Site {base_url} is responding")
            return True
        else:
            print(f"❌ Site {base_url} is not responding: {result.stderr.decode()}")
            return False
    except subprocess.TimeoutExpired:
        print(f"❌ Site {base_url} timed out")
        return False
    except Exception as e:
        print(f"❌ Error checking site {base_url}: {str(e)}")
        return False

# Add this function to sanitize and validate your dataset JSON
def validate_dataset_redirectors(json_path, fallback_redirectors=None):
    """Validate and potentially fix redirectors in a dataset JSON file"""
    import json
    import re
    import os
    
    if fallback_redirectors is None:
        fallback_redirectors = [
            "cms-xrd-global.cern.ch",
            "xrootd-cms.infn.it",
            "cmsxrootd.fnal.gov"
        ]
    
    print(f"Validating redirectors in {json_path}...")
    
    if not os.path.exists(json_path):
        print(f"Error: JSON file {json_path} not found")
        return False
        
    with open(json_path, 'r') as f:
        data = json.load(f)
    
    modified = False
    
    # Process each dataset
    for dataset_name, files in data.items():
        print(f"Checking dataset: {dataset_name}")
        valid_files = []
        
        for file_path in files:
            if not file_path.startswith("root://"):
                print(f"  Skipping non-XRootD path: {file_path}")
                valid_files.append(file_path)
                continue
                
            # Extract redirector from path
            match = re.search(r'root://([^/]+)', file_path)
            if not match:
                print(f"  Invalid XRootD URL: {file_path}")
                continue
                
            redirector = match.group(1)
            
            # Check if the original redirector is responsive
            if check_xrootd_site_availability(f"root://{redirector}"):
                valid_files.append(file_path)
                continue
                
            # Try fallback redirectors
            store_path = file_path[file_path.find("/store/"):]
            for fallback in fallback_redirectors:
                fallback_url = f"root://{fallback}{store_path}"
                if check_xrootd_site_availability(fallback_url):
                    print(f"  Replaced {redirector} with {fallback}")
                    valid_files.append(fallback_url)
                    modified = True
                    break
            else:
                print(f"  ⚠️ Warning: Could not find working redirector for {file_path}")
        
        # Update the dataset with valid files only
        if len(valid_files) < len(files):
            print(f"  Removed {len(files) - len(valid_files)} inaccessible files")
            modified = True
        
        if not valid_files:
            print(f"  ⚠️ Warning: No accessible files found for dataset {dataset_name}")
            
        data[dataset_name] = valid_files
    
    # Write back the updated JSON if modified
    if modified:
        with open(json_path, 'w') as f:
            json.dump(data, f, indent=4)
        print(f"Updated {json_path} with working redirectors")
    else:
        print(f"No changes needed for {json_path}")
    
    return True

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
        os.system(f"mv {latest_file} src/BTVNanoCommissioning/data/lumiMasks/.")
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
        choices=list(workflows.keys()) + ["Validation", "SF", "default_comissioning"],
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

    args = parser.parse_args()
    # summarize diffeerent group for study
    scheme = {
        # scale factor workflows
        "SF": ["BTA_ttbar", "BTA_addPFMuons"],
        # Use for prompt data MC checks for analysis
        "Validation": ["ttdilep_sf", "ctag_Wc_sf"],
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
        scheme["test"] = [args.scheme]
        args.scheme = "test"
    # Check lumiMask exists and replace the Validation
    input_lumi_json = correction_config[args.campaign]["lumiMask"]
    if args.campaign != "prompt_dataMC" and not os.path.exists(
        f"src/BTVNanoCommissioning/data/lumiMasks/{input_lumi_json}"
    ):
        raise f"src/BTVNanoCommissioning/data/lumiMasks/{input_lumi_json} not exist"

    if (
        args.campaign == "prompt_dataMC"
        and correction_config[args.campaign]["lumiMask"] == "$PROMPT_DATAMC"
    ):
        input_lumi_json = get_lumi_from_web(args.year)
        os.system(
            f"sed -i 's/$PROMPT_DATAMC/{input_lumi_json}/g' src/BTVNanoCommissioning/utils/AK4_parameters.py"
        )
        print(f"======>{input_lumi_json} is used for {args.year}")

    for wf in scheme[args.scheme]:
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
                    f"Creating MC dataset: python scripts/fetch.py -c {args.campaign} --from_workflow {wf} --DAS_campaign {args.DAS_campaign} --year {args.year} {overwrite} --skipvalidation"
                )

            os.system(
                f"python scripts/fetch.py -c {args.campaign} --from_workflow {wf} --DAS_campaign {args.DAS_campaign} --year {args.year} {overwrite} --skipvalidation"
            )
            if args.debug:
                os.system(f"ls metadata/{args.campaign}/*.json")

        ## Run the workflows
        for types in predefined_sample[wf].keys():

            if (types != "data" and types != "MC") and args.scheme == "Validation":
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
                runner_config_required = f"python runner.py --wf {wf} --json metadata/{args.campaign}/{types}_{args.campaign}_{args.year}_{wf}.json {overwrite} --campaign {args.campaign} --year {args.year}"
                runner_config = ""
                for key, value in vars(args).items():
                    # Required info or things not relevant for runner, skip
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
                    ]:
                        continue
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
                        if (
                            "Validation" == args.scheme
                            and types == "MC"
                            and "limit" not in key
                        ):
                            runner_config += " --limit 50"

                        else:
                            runner_config += f" --{key}={value}"
                if 'CI' in os.environ or 'GITLAB_CI' in os.environ:
                    # Validate all redirectors in the JSON and try fix if necessary
                    json_path = f"metadata/{args.campaign}/{types}_{args.campaign}_{args.year}_{wf}.json"
                    validate_result = validate_dataset_redirectors(json_path)
                    
                    if not validate_result and types == "data":
                        print(f"⚠️ WARNING: Could not validate data redirectors in {json_path}")
                        print(f"Due to data integrity requirements, skipping data workflow")
                        continue  # Skip this workflow if validation fails for data
                    import tempfile
                    
                    # Create a temporary bash script to handle proxy detection and variable expansion
                    with tempfile.NamedTemporaryFile(mode='w', suffix='.sh', delete=False) as script_file:
                        script_path = script_file.name
                        script_file.write('#!/bin/bash\n\n')
                        script_file.write('echo "Starting proxy detection script"\n\n')
                        
                        # Try multiple proxy paths with proper shell expansion
                        script_file.write('PROXY_FOUND=false\n')
                        script_file.write('for PROXY_PATH in "/cms-analysis/btv/software-and-algorithms/autobtv/proxy/x509_proxy" "/btv/software-and-algorithms/autobtv/proxy/x509_proxy" "${CI_PROJECT_DIR}/proxy/x509_proxy"; do\n')
                        script_file.write('    if [ -f "$PROXY_PATH" ]; then\n')
                        script_file.write('        echo "Found proxy at: $PROXY_PATH"\n')
                        script_file.write('        export X509_USER_PROXY="$PROXY_PATH"\n')
                        script_file.write('        echo "Set X509_USER_PROXY=$X509_USER_PROXY"\n')
                        script_file.write('        ls -la "$X509_USER_PROXY" || echo "Warning: Cannot list proxy file"\n')
                        script_file.write('        PROXY_FOUND=true\n')
                        script_file.write('        break\n')
                        script_file.write('    fi\n')
                        script_file.write('done\n\n')
                        
                        # Add warning if no proxy found
                        script_file.write('if [ "$PROXY_FOUND" = false ]; then\n')
                        script_file.write('    echo "WARNING: No proxy found in any standard location!"\n')
                        script_file.write('fi\n\n')
                        
                        # Add the runner command with all environment variables preserved
                        script_file.write(f'echo "Executing: {runner_config_required}{runner_config}"\n')
                        script_file.write(f'{runner_config_required}{runner_config}\n')
                        
                        # Exit with the runner's exit code
                        script_file.write('exit $?\n')
                    
                    # Make script executable
                    os.chmod(script_path, 0o755)
                    
                    # Execute the script
                    if args.debug:
                        print(f"Executing proxy detection script: {script_path}")
                    os.system(script_path)
                    
                    # Clean up
                    os.unlink(script_path)
                else:
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
                f"hists_{wf}_data_{args.campaign}_{args.year}_{wf}/hists_{wf}_data_{args.campaign}_{args.year}_{wf}.coffea "
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
