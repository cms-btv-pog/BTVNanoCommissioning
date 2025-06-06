import os, argparse
from BTVNanoCommissioning.workflows import workflows
from BTVNanoCommissioning.utils.sample import predefined_sample
from BTVNanoCommissioning.utils.AK4_parameters import correction_config
import os, sys, inspect
from datetime import datetime

current_dir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)

from runner import config_parser, scaleout_parser, debug_parser
  

def create_eos_upload_script():
    """Create a bash script for reliable file uploads using xrdcp only"""
    upload_script = os.path.join(os.path.dirname(__file__), "upload_to_eos.sh")
    
    print("Creating upload script using xrdcp...")
    with open(upload_script, "w") as f:
        f.write("""#!/bin/bash
# Simple Upload Script for BTVNanoCommissioning using xrdcp only
set -e

CAMPAIGN=$1
VERSION=$2
WF=$3
SCHEME=$4
FORCE=$5
YEAR="2024"  # You can make this a parameter if needed

# Parse force flag
if [[ "$FORCE" == "--force" ]]; then
    OVERWRITE="-f"
else
    OVERWRITE=""
fi

# Colors for better output
GREEN='\\033[0;32m'
YELLOW='\\033[1;33m'
RED='\\033[0;31m'
NC='\\033[0m' # No Color

echo -e "${GREEN}=== BTV Commissioning Upload ===${NC}"
echo "Campaign: $CAMPAIGN"
echo "Workflow: $WF"
echo "Scheme: $SCHEME"
echo "Target: /eos/user/b/btvweb/www/Commissioning/dataMC/$SCHEME/$CAMPAIGN$VERSION"

# Set up authentication for uploads
setup_auth() {
    echo -e "${YELLOW}Setting up authentication...${NC}"
    
    # Check for existing proxy files in common locations
    if [ -f "/tmp/x509up_u$(id -u)" ]; then
        echo "Found proxy in /tmp/x509up_u$(id -u)"
        export X509_USER_PROXY="/tmp/x509up_u$(id -u)"
    fi
    
    # Check GitLab CI specific location
    if [ -f "/builds/cms-analysis/btv/software-and-algorithms/autobtv/proxy/x509_proxy" ]; then
        echo "Found GitLab CI proxy file"
        export X509_USER_PROXY="/builds/cms-analysis/btv/software-and-algorithms/autobtv/proxy/x509_proxy"
    fi
    
    # If we find a proxy file in the current directory or parent directories
    PROXY_FILE=$(find $HOME -name "x509_proxy" -o -name "x509up_*" 2>/dev/null | head -1)
    if [ -n "$PROXY_FILE" ]; then
        echo "Found proxy file: $PROXY_FILE"
        export X509_USER_PROXY="$PROXY_FILE"
    fi
    
    # Display proxy info if available
    if [ -n "$X509_USER_PROXY" ]; then
        voms-proxy-info -all 2>/dev/null || echo "Cannot read proxy info"
    else
        echo -e "${YELLOW}No proxy file found, will try default authentication${NC}"
    fi
}

# Find all relevant files and prepare upload structure
prepare_files() {
    echo -e "${YELLOW}Preparing files for upload...${NC}"
    
    # Create temporary directories for collecting files
    mkdir -p ./upload_tmp/coffea/$CAMPAIGN/$WF
    mkdir -p ./upload_tmp/plots/$CAMPAIGN$VERSION/$WF
    
    # Copy index files
    if [ -f "scripts/index.php" ]; then
        cp scripts/index.php ./upload_tmp/coffea/$CAMPAIGN/
        cp scripts/index.php ./upload_tmp/coffea/$CAMPAIGN/$WF/
        cp scripts/index.php ./upload_tmp/plots/$CAMPAIGN$VERSION/
        cp scripts/index.php ./upload_tmp/plots/$CAMPAIGN$VERSION/$WF/
    else
        echo "Warning: index.php not found in scripts directory"
    fi
    
    # Find all coffea files - FIXED pattern to match actual file structure
    echo "Finding coffea files..."
    COFFEA_FILES_COUNT=0
    
    # Direct pattern for exact coffea files
    echo "Looking for coffea files matching exact pattern: ./hists_${WF}_*_${CAMPAIGN}_${YEAR}_${WF}/*.coffea"
    for COFFEA_FILE in $(find . -name "hists_${WF}_*_${CAMPAIGN}_${YEAR}_${WF}.coffea" 2>/dev/null); do
        echo "Found coffea file: $COFFEA_FILE"
        cp "$COFFEA_FILE" ./upload_tmp/coffea/$CAMPAIGN/$WF/
        COFFEA_FILES_COUNT=$((COFFEA_FILES_COUNT + 1))
    done
    
    # If no files found, try broader patterns
    if [ $COFFEA_FILES_COUNT -eq 0 ]; then
        echo "No coffea files found with exact pattern, trying broader patterns..."
        # Try to find coffea files in subdirectories
        for COFFEA_FILE in $(find . -path "*/hists_${WF}_*/*.coffea" 2>/dev/null); do
            echo "Found coffea file with broader pattern: $COFFEA_FILE"
            cp "$COFFEA_FILE" ./upload_tmp/coffea/$CAMPAIGN/$WF/
            COFFEA_FILES_COUNT=$((COFFEA_FILES_COUNT + 1))
        done
    fi
    
    # Find all plot files - use multiple patterns
    echo "Finding plot files..."
    PLOT_FILES_COUNT=0
    
    # Look in the plot directory with flexible patterns
    # Try specific pattern first (workflow_campaign_year)
    for PLOT_FILE in $(find ./plot -path "*/${WF}_${CAMPAIGN}_${YEAR}*/*" -name "*.png" -o -name "*.pdf" 2>/dev/null); do
        echo "Found plot file: $PLOT_FILE"
        cp "$PLOT_FILE" ./upload_tmp/plots/$CAMPAIGN$VERSION/$WF/
        PLOT_FILES_COUNT=$((PLOT_FILES_COUNT + 1))
    done
    
    # Try a more general pattern based on workflow name
    if [ $PLOT_FILES_COUNT -eq 0 ]; then
        for PLOT_FILE in $(find ./plot -path "*/${WF}_*/*" -name "*.png" -o -name "*.pdf" 2>/dev/null); do
            echo "Found plot file with broader pattern: $PLOT_FILE"
            cp "$PLOT_FILE" ./upload_tmp/plots/$CAMPAIGN$VERSION/$WF/
            PLOT_FILES_COUNT=$((PLOT_FILES_COUNT + 1))
        done
    fi
    
    # Last resort - check if any files in the plot directory match the workflow name pattern
    if [ $PLOT_FILES_COUNT -eq 0 ]; then
        for PLOT_FILE in $(find ./plot -path "*${WF}*/*" -name "*.png" -o -name "*.pdf" 2>/dev/null); do
            echo "Found plot file with workflow in path: $PLOT_FILE"
            cp "$PLOT_FILE" ./upload_tmp/plots/$CAMPAIGN$VERSION/$WF/
            PLOT_FILES_COUNT=$((PLOT_FILES_COUNT + 1))
        done
    fi
    
    # Check if the workflow name is part of a different format
    if [ $PLOT_FILES_COUNT -eq 0 ]; then
        # Look for files with ttdilep_sf (for example) somewhere in the path
        for PLOT_FILE in $(find ./plot -type f \( -name "*.png" -o -name "*.pdf" \) | grep -i "${WF}" 2>/dev/null); do
            echo "Found plot file with workflow in filename: $PLOT_FILE"
            cp "$PLOT_FILE" ./upload_tmp/plots/$CAMPAIGN$VERSION/$WF/
            PLOT_FILES_COUNT=$((PLOT_FILES_COUNT + 1))
        done
    fi
    
    # Report found files
    echo "Found $COFFEA_FILES_COUNT coffea files and $PLOT_FILES_COUNT plot files"
    
    # List files for verification
    echo "Coffea files to upload:"
    ls -la ./upload_tmp/coffea/$CAMPAIGN/$WF/
    
    echo "Plot files to upload:"
    ls -la ./upload_tmp/plots/$CAMPAIGN$VERSION/$WF/
    
    # Check if we have any files to upload
    if [ $COFFEA_FILES_COUNT -eq 0 ] && [ $PLOT_FILES_COUNT -eq 0 ]; then
        echo -e "${RED}No files found to upload!${NC}"
        return 1
    fi
    
    return 0
}

# Upload using xrdcp - FIXED to not use unsupported --mkdir option
upload_files() {
    echo -e "${YELLOW}Uploading files...${NC}"
    UPLOAD_SUCCESS=0
    
    # Upload coffea files if any exist
    if [ -n "$(ls -A ./upload_tmp/coffea/$CAMPAIGN/$WF/ 2>/dev/null)" ]; then
        echo "Uploading coffea files..."
        TARGET_PATH="root://eosuser.cern.ch//eos/user/b/btvweb/www/Commissioning/dataMC/$SCHEME/$CAMPAIGN/"
        
        echo "DEBUG: xrdcp source: ./upload_tmp/coffea/$CAMPAIGN/$WF/"
        echo "DEBUG: xrdcp target: $TARGET_PATH"
        
        # Create directories first using xrdfs
        echo "Creating target directories using xrdfs..."
        xrdfs root://eosuser.cern.ch/ mkdir -p /eos/user/b/btvweb/www/Commissioning/dataMC/$SCHEME/$CAMPAIGN/$WF
        
        # Then upload using xrdcp without --mkdir
        xrdcp -r -p $OVERWRITE ./upload_tmp/coffea/$CAMPAIGN/$WF/ "$TARGET_PATH"
        XRDCP_STATUS=$?
        
        if [ $XRDCP_STATUS -eq 0 ]; then
            echo -e "${GREEN}Coffea files uploaded successfully${NC}"
            UPLOAD_SUCCESS=1
        else
            echo -e "${RED}Failed to upload coffea files, status: $XRDCP_STATUS${NC}"
        fi
    else
        echo "No coffea files to upload"
    fi
    
    # Upload plot files if any exist
    if [ -n "$(ls -A ./upload_tmp/plots/$CAMPAIGN$VERSION/$WF/ 2>/dev/null)" ]; then
        echo "Uploading plot files..."
        TARGET_PATH="root://eosuser.cern.ch//eos/user/b/btvweb/www/Commissioning/dataMC/$SCHEME/$CAMPAIGN$VERSION/"
        
        echo "DEBUG: xrdcp source: ./upload_tmp/plots/$CAMPAIGN$VERSION/$WF/"
        echo "DEBUG: xrdcp target: $TARGET_PATH"
        
        # Create directories first using xrdfs
        echo "Creating target directories using xrdfs..."
        xrdfs root://eosuser.cern.ch/ mkdir -p /eos/user/b/btvweb/www/Commissioning/dataMC/$SCHEME/$CAMPAIGN$VERSION/$WF
        
        # Then upload using xrdcp without --mkdir
        xrdcp -r -p $OVERWRITE ./upload_tmp/plots/$CAMPAIGN$VERSION/$WF/ "$TARGET_PATH"
        XRDCP_STATUS=$?
        
        if [ $XRDCP_STATUS -eq 0 ]; then
            echo -e "${GREEN}Plot files uploaded successfully${NC}"
            UPLOAD_SUCCESS=1
        else
            echo -e "${RED}Failed to upload plot files, status: $XRDCP_STATUS${NC}"
        fi
    else
        echo "No plot files to upload"
    fi
    
    # Cleanup
    echo "Cleaning up temporary files..."
    rm -rf ./upload_tmp
    
    if [ $UPLOAD_SUCCESS -eq 1 ]; then
        return 0
    else
        return 1
    fi
}

# Main execution flow
main() {
    # Set up authentication
    setup_auth
    
    # Prepare files for upload
    prepare_files
    PREP_STATUS=$?
    
    if [ $PREP_STATUS -ne 0 ]; then
        echo -e "${RED}Failed to prepare files for upload${NC}"
        return 1
    fi
    
    # Upload files
    upload_files
    UPLOAD_STATUS=$?
    
    if [ $UPLOAD_STATUS -eq 0 ]; then
        echo -e "${GREEN}Upload completed successfully${NC}"
        return 0
    else
        echo -e "${RED}Upload failed${NC}"
        return 1
    fi
}

# Run the main function
main
""")
    os.chmod(upload_script, 0o755)
    print(f"Created upload script: {upload_script}")
    return upload_script
    
def add_file_level_fallbacks():
    """Monkey-patch uproot to use fallback redirectors for individual files"""
    import os
    
    try:
        import uproot
        from uproot.source.xrootd import XRootDSource
        import functools
        
        # Store the original open method
        original_open = XRootDSource._open
        
        # Create a list of fallback redirectors
        fallbacks = [
            "root://cms-xrd-global.cern.ch/",
            "root://xrootd-cms.infn.it/",
            "root://cmsxrootd.fnal.gov/"
        ]
        
        # Define our wrapper function with fallbacks
        @functools.wraps(original_open)
        def open_with_fallbacks(self):
            try:
                return original_open(self)
            except Exception as original_error:
                # Extract the file path after /store/
                if "/store/" in self._file_path:
                    file_path = "/store/" + self._file_path.split("/store/")[1]
                    
                    # Try all fallback redirectors
                    for fallback in fallbacks:
                        try:
                            print(f"⚠️ Primary redirector failed, trying fallback: {fallback}")
                            # Temporarily modify the file path
                            original_path = self._file_path
                            self._file_path = f"{fallback}{file_path}"
                            
                            # Try to open with this redirector
                            result = original_open(self)
                            print(f"✅ Successfully accessed file via fallback: {fallback}")
                            return result
                        except Exception as fallback_error:
                            # Restore the original path before trying the next fallback
                            self._file_path = original_path
                            print(f"❌ Fallback {fallback} also failed: {str(fallback_error)}")
                    
                # If all fallbacks fail or file path doesn't contain /store/, raise original error
                raise original_error
        
        # Apply our monkey patch
        XRootDSource._open = open_with_fallbacks
        print("✅ Applied XRootD fallback system")
        
    except (ImportError, AttributeError) as e:
        print(f"⚠️ Could not set up XRootD fallbacks: {e}")

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
        print(f"⚠️ Dataset file {json_file} is {age_in_minutes:.1f} minutes old, needs refreshing")
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
    # Add test_lumi flag
    # TODO: remove after testing
    parser.add_argument(
        "--test_lumi",
        action="store_true",
        help="Test luminosity calculation with limited dataset (20 files)",
    )
    
    args = parser.parse_args()
    # summarize diffeerent group for study
    scheme = {
        # scale factor workflows
        "SF": ["BTA_ttbar", "BTA_addPFMuons"],
        # Use for prompt data MC checks for analysis
        #"Validation": ["ttdilep_sf", "ctag_Wc_sf"],
        "Validation": ["ttdilep_sf"], ### TODO: change to all after testing
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
                    
                json_file = f"metadata/{args.campaign}/{types}_{args.campaign}_{args.year}_{wf}.json"

                # Check if dataset needs refreshing because of age
                if should_refresh_dataset(json_file):
                    print(f"Refreshing dataset: {json_file}")
                    fetch_cmd = f"python scripts/fetch.py -c {args.campaign} --from_workflow {wf} --DAS_campaign {args.DAS_campaign} --year {args.year} {overwrite} --skipvalidation --overwrite"
                    os.system(fetch_cmd)
                    
                runner_config_required = f"python runner.py --wf {wf} --json metadata/{args.campaign}/{types}_{args.campaign}_{args.year}_{wf}.json {overwrite} --campaign {args.campaign} --year {args.year}"
                runner_config = ""
                limit_added = False  # Track if we've already added a limit flag

                for key, value in vars(args).items():
                    # Required info or things not relevant for runner, skip
                    # Skip keys that don't apply to runner
                    if key in ["workflow", "json", "campaign", "year", "scheme", "DAS_campaign", 
                            "version", "local", "debug", "test_lumi"]:
                        continue
                        
                    # Handle boolean flags
                    if key in ["isArray", "noHist", "overwrite", "validate", "skipbadfiles"]:
                        if value == True:
                            runner_config += f" --{key}"
                    elif value is not None:
                        if key == "limit":
                            runner_config += f" --{key}={value}"
                            limit_added = True
                        else:
                            runner_config += f" --{key}={value}"

                # Add limit for MC validation if not already present
                if "Validation" == args.scheme and types == "MC" and not limit_added:
                    runner_config += " --limit 10" ###TODO: change to 50 after testing
                    imit_added = True
                    print(f"⚠️ Running Validation with 50 files limit for MC")
                # Add test_lumi limited processing (20 files) if flag is set
                elif args.test_lumi and not limit_added:
                    runner_config += " --limit 5"
                    limit_added = True
                    print(f"⚠️ Running in test_lumi mode with 20 files limit for {types}")
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
                    
                    if args.debug:
                        print(f"Preparing to upload files for {wf}...")
                        print("Searching for files to upload:")
                        
                        # Look for coffea files
                        print("\nLooking for coffea files:")
                        os.system(f"find . -name 'hists_{wf}_*' -type d | sort")
                        os.system(f"find . -name '*.coffea' | grep {wf} | head -10")
                        
                        # Look for plot files with very flexible patterns
                        print("\nLooking for plot files:")
                        os.system(f"find ./plot -type d | sort")
                        os.system(f"find ./plot -name '*.png' | head -5")
                        os.system(f"find ./plot -path '*{wf}*' -name '*.png' | head -5")
                        
                    # Create and get the upload script path
                    upload_script = create_eos_upload_script()
                    force_flag = "--force" if args.overwrite else ""
                    
                    # Execute the upload script with the correct parameters
                    upload_cmd = f"{upload_script} {args.campaign} {args.version} {wf} {args.scheme} {force_flag}"
                    print(f"Executing: {upload_cmd}")
                    exit_code = os.system(upload_cmd)
                    
                    if exit_code != 0:
                        print(f"Warning: Upload script exited with code {exit_code}")
                        # Continue with local copies
                        print("Local copies of files have been created")
            else:
                raise Exception(
                    f"No input coffea hists_{wf}_data_{args.campaign}_{args.year}_{wf}/hists_{wf}_data_{args.campaign}_{args.year}_{wf}.coffea or hists_{wf}_MC_{args.campaign}_{args.year}_{wf}/hists_{wf}_MC_{args.campaign}_{args.year}_{wf}.coffea"
                )
    # revert prompt_dataMC lumimask
    if args.campaign == "prompt_dataMC":
        os.system(
            f"sed -i 's/{input_lumi_json}/$PROMPT_DATAMC/g' src/BTVNanoCommissioning/utils/AK4_parameters.py"
        )
