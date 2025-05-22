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
def validate_and_fix_redirectors(json_path, fallback_redirectors=None):
    """Validate and fix redirectors in dataset JSON file"""
    import json
    import re
    import os
    import tempfile
    import subprocess
    from datetime import datetime
    
    print(f"🔍 [{datetime.now().strftime('%H:%M:%S')}] Validating and fixing redirectors in {json_path}")
    
    if fallback_redirectors is None:
        fallback_redirectors = ["cms-xrd-global.cern.ch", "xrootd-cms.infn.it", "cmsxrootd.fnal.gov"]
    
    if not os.path.exists(json_path):
        print(f"❌ Error: JSON file {json_path} not found")
        return False
    
    # Create a script to find working redirectors and update JSON file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.sh', delete=False) as script_file:
        script_path = script_file.name
        script_file.write('#!/bin/bash\n\n')
        
        # Try to source potential grid environments
        script_file.write('echo "🔎 Searching for XRootD environment..."\n')
        script_file.write('[ -f /cvmfs/grid.cern.ch/etc/profile.d/setup-cvmfs-ui.sh ] && source /cvmfs/grid.cern.ch/etc/profile.d/setup-cvmfs-ui.sh &>/dev/null && echo "✅ Sourced CVMFS grid environment"\n')
        script_file.write('[ -f /cvmfs/cms.cern.ch/cmsset_default.sh ] && source /cvmfs/cms.cern.ch/cmsset_default.sh &>/dev/null && echo "✅ Sourced CMSSW environment"\n')
        
        # Check for Python XRootD module
        script_file.write('if python3 -c "from XRootD import client; print(\'XRootD Python module available\')" &>/dev/null; then\n')
        script_file.write('    echo "✅ Found Python XRootD module"\n')
        script_file.write('    HAS_PYTHON_XROOTD=true\n')
        script_file.write('else\n')
        script_file.write('    echo "❌ Python XRootD module not available"\n')
        script_file.write('    HAS_PYTHON_XROOTD=false\n')
        script_file.write('fi\n\n')
        
        # Locate xrdfs
        script_file.write('XRDFS_CMD=""\n')
        script_file.write('for xrd_path in $(which -a xrdfs 2>/dev/null) /cvmfs/grid.cern.ch/*/bin/xrdfs /usr/bin/xrdfs; do\n')
        script_file.write('    if [ -x "$xrd_path" ]; then\n')
        script_file.write('        echo "✅ Found xrdfs at: $xrd_path"\n')
        script_file.write('        XRDFS_CMD="$xrd_path"\n')
        script_file.write('        break\n')
        script_file.write('    fi\n')
        script_file.write('done\n\n')
        
        # Function to check if a redirector is working
        script_file.write('function check_redirector() {\n')
        script_file.write('    local redirector=$1\n')
        script_file.write('    echo "🔄 Testing redirector: $redirector..."\n')
        
        # Try xrdfs ping
        script_file.write('    if [ -n "$XRDFS_CMD" ]; then\n')
        script_file.write('        if "$XRDFS_CMD" "$redirector" ping &>/dev/null; then\n')
        script_file.write('            echo "✅ Redirector $redirector is responding via xrdfs"\n')
        script_file.write('            return 0\n')
        script_file.write('        else\n')
        script_file.write('            echo "⚠️ Redirector $redirector is NOT responding via xrdfs"\n')
        script_file.write('        fi\n')
        script_file.write('    fi\n')
        
        # Try Python XRootD
        script_file.write('    if [ "$HAS_PYTHON_XROOTD" = true ]; then\n')
        script_file.write('        if python3 -c "from XRootD import client; exit(0 if client.FileSystem(\'root://$redirector\').ping()[0].ok else 1)" &>/dev/null; then\n')
        script_file.write('            echo "✅ Redirector $redirector is responding via Python XRootD"\n')
        script_file.write('            return 0\n')
        script_file.write('        else\n')
        script_file.write('            echo "⚠️ Redirector $redirector is NOT responding via Python XRootD"\n')
        script_file.write('        fi\n')
        script_file.write('    fi\n')
        
        script_file.write('    echo "❌ Redirector $redirector is not responding via any method"\n')
        script_file.write('    return 1\n')
        script_file.write('}\n\n')
        
        # Function to check if a file exists
        script_file.write('function check_file_exists() {\n')
        script_file.write('    local redirector=$1\n')
        script_file.write('    local filepath=$2\n')
        
        # Try xrdfs stat
        script_file.write('    if [ -n "$XRDFS_CMD" ]; then\n')
        script_file.write('        if "$XRDFS_CMD" "$redirector" stat "$filepath" &>/dev/null; then\n')
        script_file.write('            return 0\n')
        script_file.write('        fi\n')
        script_file.write('    fi\n')
        
        # Try Python XRootD
        script_file.write('    if [ "$HAS_PYTHON_XROOTD" = true ]; then\n')
        script_file.write('        if python3 -c "from XRootD import client; exit(0 if client.FileSystem(\'root://$redirector\').stat(\'$filepath\')[0].ok else 1)" &>/dev/null; then\n')
        script_file.write('            return 0\n')
        script_file.write('        fi\n')
        script_file.write('    fi\n')
        
        script_file.write('    return 1\n')
        script_file.write('}\n\n')
        
        # JSON file path
        script_file.write(f'JSON_FILE="{json_path}"\n')
        script_file.write('JSON_BACKUP="${JSON_FILE}.bak"\n')
        script_file.write('echo "📄 Reading dataset file: $JSON_FILE"\n')
        script_file.write('cp "$JSON_FILE" "$JSON_BACKUP"\n\n')
        
        # Test all redirectors including fallbacks
        script_file.write('declare -A REDIRECTOR_STATUS\n')
        
        # First check existing redirectors in the JSON
        script_file.write('echo "🔎 Extracting redirectors from dataset JSON..."\n')
        script_file.write('REDIRECTORS=$(python3 -c \'import json, re, sys; data = json.load(open(sys.argv[1])); urls = set(); redirector_pattern = re.compile(r"root://([^/]+)"); [urls.add(redirector_pattern.search(url).group(1)) for dataset in data.values() for url in dataset if url.startswith("root://") and redirector_pattern.search(url)]; print("\\n".join(urls))\' "$JSON_FILE" 2>/dev/null)\n\n')
        
        script_file.write('echo "⏳ Testing existing redirectors..."\n')
        script_file.write('while read -r redirector; do\n')
        script_file.write('    if [ -n "$redirector" ]; then\n')
        script_file.write('        if check_redirector "$redirector"; then\n')
        script_file.write('            REDIRECTOR_STATUS["$redirector"]="ok"\n')
        script_file.write('        else\n')
        script_file.write('            REDIRECTOR_STATUS["$redirector"]="fail"\n')
        script_file.write('        fi\n')
        script_file.write('    fi\n')
        script_file.write('done <<< "$REDIRECTORS"\n\n')
        
        # Then check fallback redirectors
        for fallback in fallback_redirectors:
            script_file.write(f'echo "⏳ Testing fallback redirector: {fallback}"\n')
            script_file.write(f'if check_redirector "{fallback}"; then\n')
            script_file.write(f'    REDIRECTOR_STATUS["{fallback}"]="ok"\n')
            script_file.write('else\n')
            script_file.write(f'    REDIRECTOR_STATUS["{fallback}"]="fail"\n')
            script_file.write('fi\n\n')
        
        # Find working fallback redirectors
        script_file.write('WORKING_FALLBACKS=("")\n')
        for fallback in fallback_redirectors:
            script_file.write(f'if [ "${{REDIRECTOR_STATUS["{fallback}"]}}" = "ok" ]; then\n')
            script_file.write(f'    WORKING_FALLBACKS+=("{fallback}")\n')
            script_file.write('fi\n')
        
        # Summary of redirector status
        script_file.write('echo "📊 Redirector Status Summary:"\n')
        script_file.write('for redirector in "${!REDIRECTOR_STATUS[@]}"; do\n')
        script_file.write('    if [ "${REDIRECTOR_STATUS[$redirector]}" = "ok" ]; then\n')
        script_file.write('        echo "  ✅ $redirector: WORKING"\n')
        script_file.write('    else\n')
        script_file.write('        echo "  ❌ $redirector: FAILING"\n')
        script_file.write('    fi\n')
        script_file.write('done\n\n')
        
        # Now process each URL and replace failing redirectors
        script_file.write('echo "🔄 Processing dataset file URLs..."\n')
        script_file.write('python3 -c \'\n')
        script_file.write('import json, re, sys, os\n')
        script_file.write('import subprocess\n')
        script_file.write('\n')
        script_file.write('# Load redirector status from the environment\n')
        script_file.write('redirector_status = {}\n')
        script_file.write('for key, value in os.environ.items():\n')
        script_file.write('    if key.startswith("REDIRECTOR_STATUS_"):\n')
        script_file.write('        redirector = key[18:]\n')
        script_file.write('        redirector_status[redirector] = value\n')
        script_file.write('\n')
        script_file.write('# Load fallbacks from environment\n')
        script_file.write('fallbacks = []\n')
        script_file.write('fallback_env = os.environ.get("WORKING_FALLBACKS", "")\n')
        script_file.write('if fallback_env:\n')
        script_file.write('    fallbacks = fallback_env.split(":")\n')
        script_file.write('\n')
        script_file.write('# Function to check if file exists using subprocess\n')
        script_file.write('def check_file_exists(redirector, filepath):\n')
        script_file.write('    check_cmd = f"bash -c \'source $0; check_file_exists {redirector} {filepath}\'" "{sys.argv[3]}"\n')
        script_file.write('    try:\n')
        script_file.write('        result = subprocess.run(check_cmd, shell=True, capture_output=True)\n')
        script_file.write('        return result.returncode == 0\n')
        script_file.write('    except Exception:\n')
        script_file.write('        return False\n')
        script_file.write('\n')
        script_file.write('# Load the JSON file\n')
        script_file.write('with open(sys.argv[1], "r") as f:\n')
        script_file.write('    data = json.load(f)\n')
        script_file.write('\n')
        script_file.write('redirector_pattern = re.compile(r"root://([^/]+)")\n')
        script_file.write('modified = False\n')
        script_file.write('stats = {"total": 0, "fixed": 0, "failed": 0}\n')
        script_file.write('\n')
        script_file.write('# Process each dataset\n')
        script_file.write('for dataset, urls in data.items():\n')
        script_file.write('    print(f"  Dataset: {dataset} - {len(urls)} files")\n')
        script_file.write('    new_urls = []\n')
        script_file.write('    \n')
        script_file.write('    for url in urls:\n')
        script_file.write('        stats["total"] += 1\n')
        script_file.write('        if not url.startswith("root://"):\n')
        script_file.write('            new_urls.append(url)  # Non-root URL, keep as is\n')
        script_file.write('            continue\n')
        script_file.write('            \n')
        script_file.write('        match = redirector_pattern.search(url)\n')
        script_file.write('        if not match:\n')
        script_file.write('            new_urls.append(url)  # Invalid format, keep as is\n')
        script_file.write('            continue\n')
        script_file.write('            \n')
        script_file.write('        redirector = match.group(1)\n')
        script_file.write('        current_status = os.environ.get(f"REDIRECTOR_STATUS_{redirector}", "unknown")\n')
        script_file.write('        \n')
        script_file.write('        # If redirector is working, keep URL as is\n')
        script_file.write('        if current_status == "ok":\n')
        script_file.write('            new_urls.append(url)\n')
        script_file.write('        else:\n')
        script_file.write('            # Extract file path\n')
        script_file.write('            path_start = url.find("/store/")\n')
        script_file.write('            if path_start == -1:\n')
        script_file.write('                path_start = url.find("/", url.find(redirector) + len(redirector))\n')
        script_file.write('            \n')
        script_file.write('            if path_start == -1:\n')
        script_file.write('                new_urls.append(url)  # Can\'t find path, keep as is\n')
        script_file.write('                stats["failed"] += 1\n')
        script_file.write('                print(f"    ⚠️ Can\'t parse path in: {url}")\n')
        script_file.write('                continue\n')
        script_file.write('                \n')
        script_file.write('            filepath = url[path_start:]\n')
        script_file.write('            fixed = False\n')
        script_file.write('            \n')
        script_file.write('            # Try each working fallback\n')
        script_file.write('            for fallback in fallbacks:\n')
        script_file.write('                if not fallback:  # Skip empty strings\n')
        script_file.write('                    continue\n')
        script_file.write('                    \n')
        script_file.write('                new_url = f"root://{fallback}{filepath}"\n')
        script_file.write('                print(f"    Testing {url} => {new_url}")\n')
        script_file.write('                \n')
        script_file.write('                if check_file_exists(fallback, filepath):\n')
        script_file.write('                    print(f"    ✅ Replaced {redirector} with {fallback}")\n')
        script_file.write('                    new_urls.append(new_url)\n')
        script_file.write('                    fixed = True\n')
        script_file.write('                    stats["fixed"] += 1\n')
        script_file.write('                    modified = True\n')
        script_file.write('                    break\n')
        script_file.write('            \n')
        script_file.write('            if not fixed:\n')
        script_file.write('                print(f"    ⚠️ No working redirector found for {url}")\n')
        script_file.write('                new_urls.append(url)  # Keep original if no working alternative\n')
        script_file.write('                stats["failed"] += 1\n')
        script_file.write('    \n')
        script_file.write('    # Update dataset with new URLs\n')
        script_file.write('    data[dataset] = new_urls\n')
        script_file.write('\n')
        script_file.write('# Write the updated JSON\n')
        script_file.write('if modified:\n')
        script_file.write('    with open(sys.argv[2], "w") as f:\n')
        script_file.write('        json.dump(data, f, indent=4)\n')
        script_file.write('    print(f"\\n✅ Updated JSON saved to {sys.argv[2]}")\n')
        script_file.write('    print(f"📊 Stats: {stats[\\"total\\"]} URLs processed, {stats[\\"fixed\\"]} fixed, {stats[\\"failed\\"]} failed")\n')
        script_file.write('else:\n')
        script_file.write('    print("✅ No changes needed to JSON file")\n')
        script_file.write('\' "$JSON_FILE" "$JSON_FILE" "$(realpath "$0")"\n\n')
        
        # Export redirector status for Python script
        script_file.write('# Export redirector status for Python\n')
        script_file.write('for redirector in "${!REDIRECTOR_STATUS[@]}"; do\n')
        script_file.write('    # Safe variable name by replacing special chars\n')
        script_file.write('    safe_redirector=$(echo "$redirector" | tr \':.:\' \'__\')\n')
        script_file.write('    export "REDIRECTOR_STATUS_$safe_redirector=${REDIRECTOR_STATUS[$redirector]}"\n')
        script_file.write('done\n\n')
        
        # Export working fallbacks
        script_file.write('export WORKING_FALLBACKS="')
        for i, fallback in enumerate(fallback_redirectors):
            if i > 0:
                script_file.write(':')
            script_file.write(fallback)
        script_file.write('"\n\n')
        
        # Verify script completed successfully
        script_file.write('echo "✅ Redirector validation completed"\n')
        script_file.write('exit 0\n')
    
    # Make script executable
    os.chmod(script_path, 0o755)
    
    # Run the script
    print(f"📋 Running redirector validation script: {script_path}")
    exit_code = os.system(script_path)
    
    # Clean up
    os.unlink(script_path)
    
    return exit_code == 0


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
                    print(f"🔄 [{datetime.now().strftime('%H:%M:%S')}] Validating redirectors in {json_path}")
                    
                    validate_success = validate_and_fix_redirectors(json_path)
                    if not validate_success:
                        print(f"⚠️ [{datetime.now().strftime('%H:%M:%S')}] Warning: Redirector validation had issues")
                        if types == "data":
                            print(f"⚠️ WARNING: Skipping data workflow due to redirector validation issues")
                            continue
                    
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
