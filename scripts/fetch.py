import sys
import os
import json
import argparse
from collections import defaultdict
import uproot
import numpy as np
import time
import re
from BTVNanoCommissioning.workflows import workflows
from BTVNanoCommissioning.utils.sample import predefined_sample
from BTVNanoCommissioning.utils.xrootdtools import get_xrootd_sites_map

# Adapt some developments from Andrey Pozdnyakov in CoffeaRunner https://github.com/cms-rwth/CoffeaRunner/blob/master/filefetcher/fetch.py
parser = argparse.ArgumentParser(
    description="Run analysis on baconbits files using processor coffea files"
)
parser.add_argument(
    "-i",
    "--input",
    default=None,
    type=str,
    help="List of samples in DAS (default: %(default)s)",
)
parser.add_argument(
    "-o",
    "--output",
    default=r"test_my_samples.json",
    help="Site (default: %(default)s)",
)
parser.add_argument(
    "--xrd",
    default=None,
    type=str,
    help="xrootd prefix string otherwise get from available sites",
)

parser.add_argument(
    "--from_path",
    action="store_true",
    help="For samples that are not published on DAS. If this option is set then the format of the --input file must be adjusted. It should be: \n dataset_name path_to_files.",
    default=False,
)
parser.add_argument(
    "--from_dataset",
    help="input dataset only",
    action="store_true",
    default=False,
)
parser.add_argument(
    "-wf",
    "--from_workflow",
    help="Use the predefined workflows",
    choices=list(workflows.keys()),
    default=None,
)
parser.add_argument(
    "--testfile",
    action="store_true",
    help="Construct file list in the test directory. Specify the test directory path, create the json file for individual dataset",
    default=False,
)
parser.add_argument(
    "--whitelist_sites",
    help="White list fot sites",
    default=None,
)
parser.add_argument(
    "--blacklist_sites",
    help="Black list for sites",
    default=None,
)
parser.add_argument(
    "--limit", help="Limit numbers of file to create json", default=None, type=int
)

parser.add_argument(
    "-r",
    "--redirector",
    help="xrootd ridirector in case sites are not found",
    choices=["infn", "fnal", "cern"],
    default="infn",
)
parser.add_argument(
    "-j", "--ncpus", help="Number of CPUs to use for validation", default="4"
)
parser.add_argument(
    "--skipvalidation",
    action="store_true",
    help="If true, the readability of files will not be validated.",
    default=False,
)
parser.add_argument(
    "--overwrite",
    action="store_true",
    help="Overwrite existing file?",
    default=False,
)

parser.add_argument(
    "--DAS_campaign",
    help="campaign info, specifying dataset name in DAS. If you are running with ```from_workflow`` option, please do ```campaign1,campaign2,campaign3``` split by `,`. E.g. `--DAS_campaign Run2022C*Sep2023,Run2022D*Sep2023,Run3Summer22NanoAODv12-130X`",
    default=None,
    type=str,
)
parser.add_argument(
    "-c",
    "--campaign",
    help="campaign name (same as the campaign in runner.py)",
    default=None,
    required=True,
    type=str,
)
parser.add_argument("--year", help="year", default=None, type=str)

parser.add_argument(
    "--verbose",
    "-v",
    action="store_true",
    help="Show detailed processing information (default: show only essential info)",
)

parser.add_argument(
    "--executor",
    choices=["iterative", "futures"],
    default="iterative",
    help="The type of executor to use for parallelization",
)
parser.add_argument(
    "--workers",
    "-w",
    type=int,
    default=4,
    help="Number of workers to use for futures executor",
)

args = parser.parse_args()


def conditional_print(message, level="INFO", always_print=False):
    """Print message based on verbosity setting"""
    if always_print or args.verbose or level in ["WARNING", "ERROR"]:
        if level == "WARNING":
            print(f"WARNING: {message}")
        elif level == "ERROR":
            print(f"ERROR: {message}")
        else:
            print(message)


# Dataset status helpers - these are always printed
def print_dataset_start(dataset):
    """Print the start of dataset processing - always shown"""
    print(f"📊 Processing dataset: {dataset}")


def print_dataset_success(dataset, dsname, file_count):
    """Print successful dataset processing - always shown"""
    print(f"✅ Successfully processed: {dataset} → {dsname} ({file_count} files)")


def print_dataset_failure(dataset, reason):
    """Print failed dataset processing - always shown"""
    print(f"❌ Failed to process: {dataset} - {reason}")


def print_dataset_empty(dataset):
    """Print when a dataset has no files - always shown"""
    print(f"⚠️ No files found for: {dataset}")


site_url_formats = {}


# Based on https://github.com/PocketCoffea/PocketCoffea/blob/main/pocket_coffea/utils/rucio.py
def get_xrootd_sites_map():

    in_ci_env = "CI" in os.environ or "GITLAB_CI" in os.environ

    sites_xrootd_access = defaultdict(dict)
    # Check if the cache file has been modified in the last 10 minutes
    cache_valid = False
    if os.path.exists(".sites_map.json"):
        file_time = os.path.getmtime(".sites_map.json")
        current_time = time.time()
        # ten_minutes_ago = current_time - 600
        twenty_minutes_ago = current_time - 1200
        sixty_minutes_ago = current_time - 3600
        if file_time > sixty_minutes_ago:
            cache_valid = True

    if not os.path.exists(".sites_map.json") or not cache_valid or in_ci_env:
        if args.verbose:
            print("Loading SITECONF info")
        sites = [
            (s, "/cvmfs/cms.cern.ch/SITECONF/" + s + "/storage.json")
            for s in os.listdir("/cvmfs/cms.cern.ch/SITECONF/")
            if s.startswith("T")
        ]
        for site_name, conf in sites:
            if not os.path.exists(conf):
                continue
            try:
                data = json.load(open(conf))
            except:
                continue
            for site in data:
                if site["type"] != "DISK":
                    continue
                if site["rse"] == None:
                    continue
                for proc in site["protocols"]:
                    if proc["protocol"] == "XRootD":
                        if proc["access"] not in ["global-ro", "global-rw"]:
                            continue
                        if "prefix" not in proc:
                            if "rules" in proc:
                                for rule in proc["rules"]:
                                    sites_xrootd_access[site["rse"]][rule["lfn"]] = (
                                        rule["pfn"]
                                    )
                        else:
                            sites_xrootd_access[site["rse"]] = proc["prefix"]
        json.dump(sites_xrootd_access, open(".sites_map.json", "w"))

    return json.load(open(".sites_map.json"))


if args.DAS_campaign == "auto":
    DAS_campaign_map = {
        "Summer22_2022": "Run2022C*Sep2023,Run2022D*Sep2023,Run3Summer22NanoAODv12-130X",
        "Summer22EE_2022": "Run2022E*Sep2023,Run2022F*Sep2023,Run2022G*Sep2023,Run3Summer22EENanoAODv12-130X",
        "Summer23_2023": "Run2023C*Sep2023,Run3Summer23NanoAODv12",
        "Summer23BPix_2023": "Run2023D*Sep2023,Run3Summer23BPixNanoAODv12",
        "Summer24_2024": "Run2024*MINIv6,RunIII2024Summer24NanoAODv15",
    }
    key = f"{args.campaign}_{args.year}"
    if key in DAS_campaign_map:
        args.DAS_campaign = DAS_campaign_map[key]
        print(
            f"Automatically selected DAS_campaign {args.DAS_campaign} based on {key}."
        )
    else:
        raise (
            f"Cannot automatically assign DAS_campaign based on {key}. Valid keys:",
            list(DAS_campaign_map.keys()),
        )


def _get_pfn_for_site(path, rules):
    """
    Utility function that converts the file path to a valid pfn matching
    the file path with the site rules (regexes).
    """
    if isinstance(rules, dict):
        for rule, pfn in rules.items():
            if m := re.match(rule, path):
                grs = m.groups()
                for i in range(len(grs)):
                    pfn = pfn.replace(f"${i+1}", grs[i])
                return pfn
    else:
        # Remove leading slash from path if present
        if path.startswith("/"):
            path = path[1:]

        # Check if adding our slash would create a triple slash
        if rules.endswith("//"):
            # Redirector ends with double slash, be careful not to add another
            return rules + path
        else:
            # Normal case - add a slash between redirector and path
            return rules + "/" + path


def normalize_xrootd_url(url):
    """
    Create port and non-port variants of an XRootD URL
    """

    if not url.startswith("root://"):
        print(f"Not an XRootD URL, returning original: {url}")
        return [url]

    # Handle port removal only - nothing else should change
    urls = []
    # Always include original URL as an option
    urls.append(url)

    # Check if there's a port number in the server part
    server_part = url.split("/")[2]
    if ":" in server_part:
        try:
            # Extract the port number pattern
            port_pattern = server_part.split(":")[1]
            # Remove ONLY the port, keeping everything else identical
            no_port_url = url.replace(f":{port_pattern}", "")
            urls.append(no_port_url)
            # print(f"Created no-port variant: {no_port_url}")
        except Exception as e:
            print(f"Error extracting port: {e}")
    # else:
    #    print(f"No port found in URL, skipping variant creation")

    return urls


def access_xrootd_file(url, xrootd_tools, check_all_variants=False):
    """Try to access file with variants of the URL using both stat and open"""
    from XRootD import client
    import threading

    # Define a timeout exception and handler using threading instead of signals
    class TimeoutException(Exception):
        pass

    def with_timeout(func, args=(), kwargs={}, timeout_duration=40):
        """Run a function with a timeout using threading instead of signals"""
        result = [None]
        exception = [None]
        finished = [False]

        def worker():
            try:
                result[0] = func(*args, **kwargs)
                finished[0] = True
            except Exception as e:
                exception[0] = e
                finished[0] = True

        thread = threading.Thread(target=worker)
        thread.daemon = True
        thread.start()
        thread.join(timeout_duration)

        if finished[0]:
            if exception[0]:
                raise exception[0]
            return result[0]
        else:
            raise TimeoutException(
                f"Operation timed out after {timeout_duration} seconds"
            )

    timeout = 30  # 40 seconds timeout
    url_variants = normalize_xrootd_url(url)
    success_results = []

    for variant in url_variants:
        variant_success = False

        # First try direct file open with timeout
        try:
            file_client = client.File()

            # Use threading-based timeout instead of signals
            def open_file():
                return file_client.open(variant)

            open_status, _ = with_timeout(open_file, timeout_duration=timeout)

            if open_status.ok:
                file_client.close()
                variant_success = True
                success_results.append((True, variant, "open"))
                if not check_all_variants:
                    return True, variant
        except TimeoutException:
            print(f"⚠️ Open timed out after {timeout}s with URL variant {variant}")
        except Exception as e:
            print(f"❌ Open failed with URL variant {variant}: {e}")

        # If open fails, try stat method with timeout
        if not variant_success:
            try:
                xrd_client = client.FileSystem(variant)

                # Use threading-based timeout instead of signals
                def stat_file():
                    return xrd_client.stat(variant)

                status, response = with_timeout(stat_file, timeout_duration=timeout)

                if status.ok:
                    variant_success = True
                    success_results.append((True, variant, "stat"))
                    if not check_all_variants:
                        return True, variant
            except TimeoutException:
                print(f"⚠️ Stat timed out after {timeout}s with URL variant {variant}")
            except Exception as e:
                print(f"❌ Stat failed with URL variant {variant}: {e}")

    # If we get here and check_all_variants=True, we need to check if any variant succeeded
    if success_results:
        # Return the first success (or you could choose based on other criteria)
        return success_results[0][0], success_results[0][1]

    return False, None


def determine_execution_mode():
    """
    Determine whether to use CI execution mode or local execution mode.
    Returns True for CI mode, False for local mode.
    """
    import os

    # First check if we're in CI at all
    in_ci_environment = "CI" in os.environ or "GITLAB_CI" in os.environ

    # If not in CI, always use local execution mode
    if not in_ci_environment:
        return False

    # If in CI, check if we're in the btv_coffea environment
    conda_env = os.environ.get("CONDA_DEFAULT_ENV", "")
    conda_prefix = os.environ.get("CONDA_PREFIX", "")

    in_btv_coffea = "btv_coffea" in conda_env or "btv_coffea" in conda_prefix

    # If we're in the btv_coffea environment in CI, use local execution mode
    if in_btv_coffea:
        if args.verbose:
            print(
                "📌 Detected btv_coffea environment in CI - using local execution mode"
            )
        return False

    # Otherwise, use CI execution mode
    return True


def get_dasgoclient_command():
    """Get the full path to dasgoclient executable"""
    import os
    import shutil

    # Check if explicitly set in environment
    das_path = os.environ.get("DASGOCLIENT_PATH", "")
    if das_path and os.path.exists(das_path):
        return das_path

    # Check if it's in PATH
    in_path = shutil.which("dasgoclient")
    if in_path:
        return "dasgoclient"

    # Try standard locations
    for path in [
        "/cvmfs/cms.cern.ch/common/dasgoclient",
        "/cms.cern.ch/common/dasgoclient",
        "/usr/local/bin/dasgoclient",
    ]:
        if os.path.exists(path):
            if args.verbose:
                print(f"Found dasgoclient at: {path}")
            return path

    # Fallback to the name and let the system handle it
    print(f"WARNING: Could not find dasgoclient, using default name")
    return "dasgoclient"


def run_das_command(cmd):
    """Run a DAS command with proper environment in micromamba"""
    import os
    import subprocess
    import re
    import tempfile

    # Add debug info

    # print(f"\n==== DAS Command Debug Information ====")
    # print(f"Original command: {cmd}")

    # Check if we're in GitLab CI
    in_ci = determine_execution_mode()

    if in_ci:
        # Extract the query part
        match = re.search(r'-query="([^"]+)"', cmd)
        if match:
            query = match.group(1)
            escaped_query = query.replace("*", "\\*")
            escaped_cmd = cmd.replace(f'-query="{query}"', f'-query="{escaped_query}"')
        else:
            escaped_cmd = cmd

        # Fix paths for CI
        if (
            "/common/dasgoclient" in escaped_cmd
            and not "/cms.cern.ch/common/dasgoclient" in escaped_cmd
        ):
            escaped_cmd = escaped_cmd.replace(
                "/common/dasgoclient", "/cms.cern.ch/common/dasgoclient"
            )

        # COMPLETELY DIFFERENT APPROACH: Write a script file and execute it
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".sh", delete=False
        ) as script_file:
            script_path = script_file.name
            script_file.write("#!/bin/bash\n\n")
            script_file.write('echo "Starting DAS query from temp script"\n')
            script_file.write('echo "Command: ' + escaped_cmd + '"\n')

            # Add this near the beginning of your script generation, after the shebang
            script_file.write('echo "Checking for proxy..."\n')
            script_file.write(
                'if [ -n "$X509_USER_PROXY" ] && [ -f "$X509_USER_PROXY" ]; then\n'
            )
            script_file.write('    echo "Using proxy from $X509_USER_PROXY"\n')
            script_file.write(
                'elif [ -f "${CI_PROJECT_DIR}/proxy/x509_proxy" ]; then\n'
            )
            script_file.write(
                '    export X509_USER_PROXY="${CI_PROJECT_DIR}/proxy/x509_proxy"\n'
            )
            script_file.write(
                '    echo "Found and using proxy at ${X509_USER_PROXY}"\n'
            )
            script_file.write("else\n")
            script_file.write(
                '    echo "WARNING: No proxy found! DAS queries may fail."\n'
            )
            script_file.write("fi\n\n")

            # Add all possible CMS environment setup paths
            script_file.write("if [ -f /cms.cern.ch/cmsset_default.sh ]; then\n")
            script_file.write("    source /cms.cern.ch/cmsset_default.sh\n")
            script_file.write('    echo "Sourced /cms.cern.ch/cmsset_default.sh"\n')
            script_file.write(
                "elif [ -f /cvmfs/cms.cern.ch/cmsset_default.sh ]; then\n"
            )
            script_file.write("    source /cvmfs/cms.cern.ch/cmsset_default.sh\n")
            script_file.write(
                '    echo "Sourced /cvmfs/cms.cern.ch/cmsset_default.sh"\n'
            )
            script_file.write("else\n")
            script_file.write('    echo "WARNING: Could not find cmsset_default.sh"\n')
            script_file.write("fi\n\n")

            script_file.write('echo "Searching for dasgoclient:"\n')
            script_file.write(
                'DASGOCLIENT=""\n'
            )  # Changed variable name to match usage below
            script_file.write("if [ -f /cms.cern.ch/common/dasgoclient ]; then\n")
            script_file.write('    echo "Found at /cms.cern.ch/common/dasgoclient"\n')
            script_file.write('    DASGOCLIENT="/cms.cern.ch/common/dasgoclient"\n')
            script_file.write(
                "elif [ -f /cvmfs/cms.cern.ch/common/dasgoclient ]; then\n"
            )
            script_file.write(
                '    echo "Found at /cvmfs/cms.cern.ch/common/dasgoclient"\n'
            )
            script_file.write(
                '    DASGOCLIENT="/cvmfs/cms.cern.ch/common/dasgoclient"\n'
            )
            script_file.write("fi\n\n")

            # Add error checking for dasgoclient
            script_file.write('if [ -z "$DASGOCLIENT" ]; then\n')
            script_file.write('    echo "ERROR: dasgoclient not found!"\n')
            script_file.write("    exit 1\n")
            script_file.write("fi\n\n")

            # Extract the query and run it - properly handle the extraction
            query_match = re.search(r'-query="([^"]+)"', cmd)
            if query_match:
                query = query_match.group(1)
                # Remove double asterisks which cause problems
                query = query.replace("**", "*")
                # Add quotes around the query and execute with proper syntax
                script_file.write(
                    f'echo "Executing command: $DASGOCLIENT -query=\\"{query}\\""\n'
                )
                script_file.write(f'$DASGOCLIENT -query="{query}"\n')
            else:
                script_file.write(f"{cmd}\n")

        # Make script executable
        os.chmod(script_path, 0o755)
        # print(f"Created temporary script at: {script_path}")

        # Execute the script
        try:
            print(f"Executing script: {script_path}")
            result = subprocess.run([script_path], capture_output=True, text=True)

            print(f"Script return code: {result.returncode}")

            if result.stdout:
                print(f"Script stdout (first 800 chars): {result.stdout[:800]}")
            if result.stderr:
                print(f"Script stderr: {result.stderr}")

            if result.returncode != 0:
                print(f"Script failed with code {result.returncode}")
                return []

            output = [line for line in result.stdout.strip().split("\n") if line]
            # Remove the script's debug lines from output
            output = [
                line
                for line in output
                if not line.startswith("Starting")
                and not line.startswith("Sourced")
                and not line.startswith("Found")
                and not line.startswith("Using")
                and not line.startswith("Command")
                and not line.startswith("Checking")
                and not line.startswith("Searching")
                and not line.startswith("Trying")
                and not line.startswith("Executing")
            ]
            return output

        except Exception as e:
            print(f"Exception executing script: {e}")
            return []
        finally:
            # Clean up the temp script
            os.unlink(script_path)

    else:
        # Local environment - unchanged

        dasgoclient_path = get_dasgoclient_command()

        # Replace dasgoclient with full path if needed
        if dasgoclient_path != "dasgoclient" and "dasgoclient" in cmd:
            cmd = cmd.replace("dasgoclient", dasgoclient_path)
            # print(f"Using full dasgoclient path: {cmd}")

        # print("Executing local command with os.popen:", cmd)
        result = os.popen(cmd).read().splitlines()
        return result


def run_python_xrootd_ping(server, site, timeout=10):
    """Run Python XRootD ping operation through a bash script to ensure proxy is properly set"""
    import tempfile
    import subprocess
    import os

    # Check if we're in CI
    in_ci = determine_execution_mode()

    if in_ci:
        # Create a temporary Python script that will be executed by bash
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".py", delete=False
        ) as python_file:
            python_path = python_file.name

            # Write a self-contained Python script
            python_file.write("""
import sys
import os
import time
from XRootD import client

# Get server from command line
server = sys.argv[1]
site = sys.argv[2]

# Print environment info
print("XRootD environment settings:")
xrd_vars = [var for var in os.environ if var.startswith('XRD_')]
for var in xrd_vars:
    print(f"  {var}={os.environ[var]}")

if 'X509_USER_PROXY' in os.environ:
    print(f"Using X509_USER_PROXY: {os.environ['X509_USER_PROXY']}")
    if os.path.exists(os.environ['X509_USER_PROXY']):
        print(f"  Proxy file exists, size: {os.path.getsize(os.environ['X509_USER_PROXY'])} bytes")
    else:
        print(f"  WARNING: Proxy file doesn't exist!")

# Try the connection
print("Attempting XRootD connection...")
start_time = time.time()

try:
    # Create a connection with detailed flags
    fs = client.FileSystem(f'root://{server}')
    
    # Execute a ping operation
    status, response = fs.ping()
    response_time = time.time() - start_time
    
    if status.ok:
        print(f"SUCCESS: Site {site} redirector {server} responsive ({response_time:.2f}s)")
        sys.exit(0)
    else:
        print(f"FAILED: {status.message}")
        sys.exit(1)
except Exception as e:
    print(f"ERROR: {str(e)}")
    sys.exit(2)
""")

        # Create the bash wrapper script that sets up the environment
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".sh", delete=False
        ) as script_file:
            script_path = script_file.name
            script_file.write("#!/bin/bash\n\n")
            script_file.write(
                'echo "Starting XRootD Python operation from temp script"\n\n'
            )

            # Use the same proxy detection as run_das_command
            script_file.write('echo "Checking for proxy..."\n')
            script_file.write(
                'if [ -n "$X509_USER_PROXY" ] && [ -f "$X509_USER_PROXY" ]; then\n'
            )
            script_file.write('    echo "Using proxy from $X509_USER_PROXY"\n')
            script_file.write(
                'elif [ -f "${CI_PROJECT_DIR}/proxy/x509_proxy" ]; then\n'
            )
            script_file.write(
                '    export X509_USER_PROXY="${CI_PROJECT_DIR}/proxy/x509_proxy"\n'
            )
            script_file.write(
                '    echo "Found and using proxy at ${X509_USER_PROXY}"\n'
            )
            script_file.write("else\n")
            script_file.write(
                '    echo "WARNING: No proxy found! XRootD operations may fail."\n'
            )
            script_file.write("fi\n\n")

            # Set up XRootD environment variables
            script_file.write("# Set up XRootD environment\n")
            script_file.write("export XRD_CONNECTIONRETRY=3\n")
            script_file.write("export XRD_REQUESTTIMEOUT=60\n")
            script_file.write("export XRD_PLUGIN=gsi\n")
            script_file.write("export XRD_SECPROTOCOLS=gsi\n\n")

            # Source CMS environment
            script_file.write("if [ -f /cvmfs/cms.cern.ch/cmsset_default.sh ]; then\n")
            script_file.write("    source /cvmfs/cms.cern.ch/cmsset_default.sh\n")
            script_file.write('    echo "Sourced CMS environment"\n')
            script_file.write("fi\n\n")

            # Execute the Python script
            script_file.write(f"python3 {python_path} {server} {site}\n")
            script_file.write("XROOTD_EXIT_CODE=$?\n")
            script_file.write("exit $XROOTD_EXIT_CODE\n")

        # Make scripts executable
        os.chmod(script_path, 0o755)
        os.chmod(python_path, 0o644)

        try:
            # Execute the bash script
            result = subprocess.run(
                [script_path], capture_output=True, text=True, timeout=timeout
            )

            success = result.returncode == 0
            return {
                "success": success,
                "returncode": result.returncode,
                "stdout": result.stdout,
                "stderr": result.stderr,
            }
        finally:
            # Clean up temp files
            os.unlink(script_path)
            os.unlink(python_path)
    else:
        # For local execution, use direct Python approach
        from XRootD import client
        import time

        start_time = time.time()
        try:
            fs = client.FileSystem(f"root://{server}")
            status, response = fs.ping()
            response_time = time.time() - start_time

            success = status.ok
            message = "" if success else status.message
            return {
                "success": success,
                "message": message,
                "response_time": response_time,
            }
        except Exception as e:
            return {
                "success": False,
                "message": str(e),
                "response_time": time.time() - start_time,
            }


def run_xrdfs_command(server, timeout=5):
    """Run xrdfs command directly in bash script for CI environment"""
    import tempfile
    import subprocess
    import os
    import time

    # Check if we're in CI
    in_ci = determine_execution_mode()

    if in_ci:
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".sh", delete=False
        ) as script_file:
            script_path = script_file.name
            script_file.write("#!/bin/bash\n\n")
            script_file.write('echo "Starting xrdfs operation from temp script"\n\n')

            # Use the same proxy detection as run_das_command
            script_file.write('echo "Checking for proxy..."\n')
            script_file.write(
                'if [ -n "$X509_USER_PROXY" ] && [ -f "$X509_USER_PROXY" ]; then\n'
            )
            script_file.write('    echo "Using proxy from $X509_USER_PROXY"\n')
            script_file.write(
                'elif [ -f "${CI_PROJECT_DIR}/proxy/x509_proxy" ]; then\n'
            )
            script_file.write(
                '    export X509_USER_PROXY="${CI_PROJECT_DIR}/proxy/x509_proxy"\n'
            )
            script_file.write(
                '    echo "Found and using proxy at ${X509_USER_PROXY}"\n'
            )
            script_file.write("else\n")
            script_file.write(
                '    echo "WARNING: No proxy found! XRootD operations may fail."\n'
            )
            script_file.write("fi\n\n")

            # Set up XRootD environment variables
            script_file.write("# Set up XRootD environment\n")
            script_file.write("export XRD_CONNECTIONRETRY=3\n")
            script_file.write("export XRD_REQUESTTIMEOUT=60\n")
            script_file.write("export XRD_PLUGIN=gsi\n")
            script_file.write("export XRD_SECPROTOCOLS=gsi\n\n")

            # Source CMS environment
            script_file.write("if [ -f /cvmfs/cms.cern.ch/cmsset_default.sh ]; then\n")
            script_file.write("    source /cvmfs/cms.cern.ch/cmsset_default.sh\n")
            script_file.write('    echo "Sourced CMS environment"\n')
            script_file.write("fi\n\n")

            # Find xrdfs executable
            script_file.write("# Find xrdfs executable\n")
            script_file.write("XRDFS_PATH=$(which xrdfs 2>/dev/null)\n")
            script_file.write('if [ -z "$XRDFS_PATH" ]; then\n')
            script_file.write(
                '    for candidate in "/cvmfs/cms.cern.ch/slc7_amd64_gcc820/cms/xrootd/4.8.5-pafccj3/bin/xrdfs" "/usr/bin/xrdfs"; do\n'
            )
            script_file.write('        if [ -f "$candidate" ]; then\n')
            script_file.write('            XRDFS_PATH="$candidate"\n')
            script_file.write('            echo "Found xrdfs at $XRDFS_PATH"\n')
            script_file.write("            break\n")
            script_file.write("        fi\n")
            script_file.write("    done\n")
            script_file.write("else\n")
            script_file.write('    echo "Using xrdfs from PATH: $XRDFS_PATH"\n')
            script_file.write("fi\n\n")

            script_file.write('if [ -z "$XRDFS_PATH" ]; then\n')
            script_file.write('    echo "ERROR: xrdfs not found!"\n')
            script_file.write("    exit 1\n")
            script_file.write("fi\n\n")

            # Execute the xrdfs command
            script_file.write('echo "Executing xrdfs ping command..."\n')
            script_file.write(f"$XRDFS_PATH {server} ping\n")
            script_file.write("XRDFS_EXIT_CODE=$?\n")
            script_file.write('echo "xrdfs exit code: $XRDFS_EXIT_CODE"\n')
            script_file.write("exit $XRDFS_EXIT_CODE\n")

        # Make script executable
        os.chmod(script_path, 0o755)

        try:
            # Execute the bash script
            start_time = time.time()
            result = subprocess.run(
                [script_path], capture_output=True, text=True, timeout=timeout
            )
            response_time = time.time() - start_time

            return {
                "success": result.returncode == 0,
                "returncode": result.returncode,
                "stdout": result.stdout,
                "stderr": result.stderr,
                "response_time": response_time,
            }
        finally:
            # Clean up temp files
            os.unlink(script_path)
    else:
        # For local execution, this function should not be called
        # Return a clear error
        return {
            "success": False,
            "returncode": -1,
            "message": "run_xrdfs_command called in local mode",
            "response_time": 0,
        }


def check_xrootd_availability():
    """Check if XRootD tools are available in the current environment"""
    import subprocess
    import shutil

    # Check for xrdfs command
    xrdfs_available = shutil.which("xrdfs") is not None

    # Check for Python XRootD module
    try:
        import XRootD

        python_xrootd_available = True
    except ImportError:
        python_xrootd_available = False

    # Check for gfal tools as alternative
    gfal_available = shutil.which("gfal-ls") is not None

    if args.verbose:
        print(
            f"Environment check: xrdfs: {'✅' if xrdfs_available else '❌'}, "
            f"Python XRootD: {'✅' if python_xrootd_available else '❌'}, "
            f"gfal: {'✅' if gfal_available else '❌'}"
        )

    return {
        "xrdfs": xrdfs_available,
        "python_xrootd": python_xrootd_available,
        "gfal": gfal_available,
    }


def check_site_completion(site, dataset, xrootd_tools):
    """Check what percentage of the dataset is available at the site by directly counting files"""
    try:
        # Get total file list once for the entire dataset
        in_ci = determine_execution_mode()

        if in_ci:
            total_files_cmd = f'/cms.cern.ch/common/dasgoclient -query="file dataset={dataset} | grep file.name"'
        else:
            total_files_cmd = (
                f'dasgoclient -query="file dataset={dataset} | grep file.name"'
            )

        # print(f"Getting file list for {dataset}...")
        total_files_result = run_das_command(total_files_cmd)
        total_files_count = len(
            [f for f in total_files_result if f and not f.startswith("ERROR:")]
        )

        if total_files_count == 0:
            print_dataset_empty(dataset)
            return 0.0

        if args.verbose:
            print(f"Dataset has {total_files_count} files in total")

        # Now get files at the specific site
        if in_ci:
            site_files_cmd = f'/cms.cern.ch/common/dasgoclient -query="file dataset={dataset} site={site} | grep file.name"'
        else:
            site_files_cmd = f'dasgoclient -query="file dataset={dataset} site={site} | grep file.name"'

        # print(f"Checking files available at {site}...")
        site_files_result = run_das_command(site_files_cmd)
        site_files_count = len(
            [f for f in site_files_result if f and not f.startswith("ERROR:")]
        )

        # Calculate completion percentage
        if total_files_count > 0 and site_files_count > 0:
            completion_pct = (site_files_count / total_files_count) * 100
            if args.verbose:
                print(
                    f"Site {site} has {site_files_count}/{total_files_count} files ({completion_pct:.2f}%)"
                )
            return completion_pct
        else:
            print(f"WARNING: Site {site} has no valid files for this dataset")
            return 0.0

    except Exception as e:
        print(f"ERROR: Error checking completion for site {site}: {e}")
        import traceback

        traceback.print_exc()
        return 0.0


def check_site_responsiveness(
    site, sites_xrootd_prefix, xrootd_tools, dataset_files=None
):
    """Test if the site's redirector is currently responsive using available tools"""
    import time
    import os
    import subprocess  # Add missing import for gfal-ls
    import random

    # Skip if site not in our map
    if site not in sites_xrootd_prefix:
        print(f"Site {site} not found in redirector map")
        return False

    # If we have actual files, test ALL files in the dataset
    if dataset_files and len(dataset_files) > 0:
        total_files = len(dataset_files)
        if args.verbose:
            print(
                f"Testing accessibility of all {total_files} files for site {site}..."
            )

        # Determine print frequency based on dataset size
        if total_files < 20:
            print_frequency = 5
        elif total_files < 120:
            print_frequency = 10
        else:
            print_frequency = 50

        # Test each file
        successful_files = 0
        failed_files = []
        for idx, test_file in enumerate(dataset_files):
            # Construct the full PFN using the site's redirector rules
            full_pfn = _get_pfn_for_site(test_file, sites_xrootd_prefix[site])
            # Use the existing access_xrootd_file function
            success, working_url = access_xrootd_file(
                full_pfn, xrootd_tools, check_all_variants=True
            )

            if success:
                successful_files += 1
                # Only print status for certain files based on frequency
                if (
                    (idx + 1) % print_frequency == 0
                    or idx == 0
                    or idx == total_files - 1
                ):
                    if args.verbose:
                        print(
                            f"✅ File {idx+1}/{total_files} accessible: {working_url}"
                        )

                # Save the first successful redirector format
                if successful_files == 1:
                    parts = working_url.split("store/")
                    redirector = parts[0]
                    site_url_formats[site] = {
                        "url": working_url,
                        "redirector": redirector,
                    }
            else:
                print(f"❌ File {idx+1}/{total_files} NOT accessible")
                # Fail fast - stop checking as soon as one file fails
                print(
                    f"❌ Site {site} fails - file {idx+1}/{total_files} not accessible"
                )
                return False

        # If we've checked all files and they're all accessible, the site passes
        if args.verbose:
            print(f"✅ Site {site} passes with all {total_files} files accessible")
        return True
    else:
        # Get the redirector for this site
        xrd = None
        if isinstance(sites_xrootd_prefix[site], list):
            xrd = sites_xrootd_prefix[site][0]
        elif isinstance(sites_xrootd_prefix[site], dict):
            # For dictionary redirectors, extract a redirector string for testing
            # Pick the first redirector value from the dictionary
            first_pattern = list(sites_xrootd_prefix[site].values())[0]
            # Extract just the server part (up to the first path component)
            # Check how the redirector is formatted in the dictionary
            if "//$1" in first_pattern or "/$1" in first_pattern:
                # Create a complete test URL with a path for testing
                test_path = "/store/mc"  # Common valid path that should exist
                complete_url = _get_pfn_for_site(test_path, sites_xrootd_prefix[site])

                # Extract just the server part for the ping test (xrdfs needs server only)
                parts = complete_url.split("//")
                if len(parts) >= 2:
                    server = parts[1].split("/")[0]
                    # Set xrd to the server with proper root:// prefix
                    xrd = f"root://{server}"

                print(
                    f"Using server with path for dictionary redirector: {server} (from {complete_url})"
                )
            else:
                # For simpler dictionary patterns
                parts = first_pattern.split("//")
                if len(parts) >= 2:
                    protocol = parts[0] + "//"  # root://
                    hostname = parts[1].split("/")[0]  # hostname:port
                    xrd = f"{protocol}{hostname}"
                    print(
                        f"Using server part of dictionary redirector for {site}: {xrd}"
                    )
                else:
                    print(f"Invalid redirector format for {site}: {first_pattern}")
                    return False
        elif isinstance(sites_xrootd_prefix[site], str):
            xrd = sites_xrootd_prefix[site]

        if not xrd:
            print(f"No redirector found for site {site}")
            return False

        # Extract the server part
        server = xrd.replace("root://", "").split("/")[0]
        # print(f"Testing redirector {server} for site {site}...")

        # Try using xrdfs if available
        if xrootd_tools["xrdfs"]:
            try:
                start_time = time.time()

                # Run xrdfs command directly in local mode
                in_ci = determine_execution_mode()
                if not in_ci:
                    # Direct local execution
                    cmd = ["xrdfs", server, "ping"]
                    result = subprocess.run(cmd, capture_output=True, timeout=5)
                    response_time = time.time() - start_time

                    if result.returncode == 0:
                        print(
                            f"✅ Site {site} redirector {server} is responsive via xrdfs (response time: {response_time:.2f}s)"
                        )
                        return True
                    else:
                        print(f"❌ Site {site} redirector {server} failed xrdfs ping")
                else:
                    # Use the bash script wrapper for CI
                    result = run_xrdfs_command(server, timeout=5)
                    response_time = result.get("response_time", 0)

                    if result.get("success", False):
                        print(
                            f"✅ Site {site} redirector {server} is responsive via xrdfs (response time: {response_time:.2f}s)"
                        )
                        return True
                    else:
                        print(f"❌ Site {site} redirector {server} failed xrdfs ping")
                        if "stdout" in result:
                            print(result["stdout"])
            except Exception as e:
                print(f"Error with xrdfs ping for {site}: {e}")

        # Try using Python XRootD
        if xrootd_tools["python_xrootd"]:
            try:
                print("Attempting XRootD connection through bash script...")

                # Use the bash script wrapper for Python XRootD
                result = run_python_xrootd_ping(server, site)

                if result.get("success", False):
                    print(
                        f"✅ Site {site} redirector {server} is responsive via Python XRootD"
                    )
                    # Only try to print stdout if it exists, otherwise print the message
                    if "stdout" in result:
                        print(result["stdout"])
                    elif "message" in result:
                        print(f"Response info: {result.get('message', '')}")
                    return True
                else:
                    print(
                        f"❌ Site {site} redirector {server} failed Python XRootD ping:"
                    )
                    # Safely access all potential keys
                    for key in ["stdout", "stderr", "message"]:
                        if key in result:
                            print(f"{key}: {result[key]}")
            except Exception as e:
                print(f"Error with Python XRootD ping for {site}: {e}")

        # Try using gfal if available
        if xrootd_tools["gfal"]:
            try:
                start_time = time.time()
                cmd = ["gfal-ls", f"root://{server}//"]
                result = subprocess.run(cmd, capture_output=True, timeout=5)
                response_time = time.time() - start_time

                if result.returncode == 0:
                    print(
                        f"✅ Site {site} redirector {server} is responsive via gfal-ls (response time: {response_time:.2f}s)"
                    )
                    return True
                else:
                    print(f"❌ Site {site} redirector {server} failed gfal-ls")
            except Exception as e:
                print(f"Error with gfal-ls for {site}: {e}")

        # If all methods failed
        print(f"❌ Site {site} redirector is not responsive via any available method")
        return False


def validate_xrootd_path(path):
    """Validate an XRootD path and fix common issues"""
    if not path.startswith("root://"):
        return path

    # Split into server and path parts
    parts = path.split("//")

    if len(parts) < 2:
        return path

    server = parts[0] + "//"
    remaining = "//".join(parts[1:])

    # Fix triple slash issue
    if remaining.startswith("/"):
        # Remove the extra leading slash
        remaining = remaining.lstrip("/")
        return server + remaining

    return path


def getFilesFromDas(args):
    """Improved getFilesFromDas with multiple fallback strategies"""
    # Check if we're in GitLab CI
    in_ci = determine_execution_mode()

    fset = []
    with open(args.input) as fp:
        lines = fp.readlines()
        for line in lines:
            fset.append(line)

    fdict = {}
    # Site handling setup
    if args.blacklist_sites is not None:
        args.blacklist_sites = args.blacklist_sites.split(",")
        if args.verbose:
            print(f"Blacklist sites: {args.blacklist_sites}")
    if args.whitelist_sites is not None:
        args.whitelist_sites = args.whitelist_sites.split(",")
        if args.verbose:
            print(f"Blacklist sites: {args.blacklist_sites}")

    # Track failed datasets for summary report
    failed_datasets = []

    for dataset in fset:
        if dataset.startswith("#") or dataset.strip() == "":
            continue

        dataset = dataset.strip()
        print_dataset_start(dataset)

        # Validate dataset format
        if not dataset.startswith("/"):
            print(
                f"WARNING: Dataset '{dataset}' does not start with '/' - trying anyway"
            )
        # Try to get dsname safely
        try:

            # Extract primary dataset name
            primary_name = dataset.split("/")[1]
            # More precise check for data: should contain "/Run20XX" pattern
            # Data pattern: /MuonEG/Run2024G-MINIv6NANOv15-v3/NANOAOD
            data_run_pattern = re.search(r"/Run20\d\d[A-Z]", dataset)

            if data_run_pattern:
                # This is actual data with a real run period
                run_period = dataset.split("/")[2].split("-")[0]  # Gets 'Run2024X'
                dsname = f"{primary_name}{run_period}"
            else:
                # This is MC - use only the primary dataset name
                dsname = primary_name

        except IndexError:
            print(f"ERROR: Cannot parse dataset name from '{dataset}'")
            if args.verbose:
                print("Using the entire string as the dataset name")
            dsname = dataset.replace("/", "_").strip()

        # Try multiple DAS query approaches
        flist = []
        if in_ci:
            das_queries = [
                f'/cms.cern.ch/common/dasgoclient -query="file dataset={dataset}"',
                f'/cms.cern.ch/common/dasgoclient -query="file dataset={dataset} instance=prod/global"',
                f'/cms.cern.ch/common/dasgoclient -query="file dataset={dataset} instance=prod/phys03"',
                f'/cms.cern.ch/common/dasgoclient -query="file dataset=*{dataset}*"',
            ]
        else:
            das_queries = [
                f'dasgoclient -query="file dataset={dataset}"',
                f'dasgoclient -query="file dataset={dataset} instance=prod/global"',
                f'dasgoclient -query="file dataset={dataset} instance=prod/phys03"',
                f'dasgoclient -query="file dataset=*{dataset}*"',
            ]

        for query in das_queries:
            if flist:  # If we already have files, no need to try other queries
                break

            # print(f"Trying DAS query: {query}")

            try:
                flist = run_das_command(query)
                if flist:
                    if args.verbose:
                        print(f"Found {len(flist)} files with query: {query}")
                    break
            except Exception as e:
                print(f"Error with DAS query: {e}")

        if not flist:
            print(
                f"WARNING: No files found for dataset '{dataset}' after trying all queries"
            )
            failed_datasets.append(dataset)
            continue

        # Get sites with multiple fallback strategies
        sites = []
        if in_ci:
            site_queries = [
                f'/cms.cern.ch/common/dasgoclient -query="site dataset={dataset}"',
                f'/cms.cern.ch/common/dasgoclient -query="site dataset={dataset} instance=prod/global"',
                f'/cms.cern.ch/common/dasgoclient -query="site dataset={dataset} instance=prod/phys03"',
            ]
        else:
            site_queries = [
                f'dasgoclient -query="site dataset={dataset}"',
                f'dasgoclient -query="site dataset={dataset} instance=prod/global"',
                f'dasgoclient -query="site dataset={dataset} instance=prod/phys03"',
            ]

        for query in site_queries:
            if sites:  # If we already have sites, no need to try other queries
                break

            # print(f"Trying site query: {query}")

            try:
                sites = run_das_command(query)
                if sites:
                    if args.verbose:
                        print(f"Found {len(sites)} sites with query: {query}")
                    break
            except Exception as e:
                print(f"Error with site query: {e}")

        # Fallback to default redirector if no sites found
        if not sites:
            print(
                f"WARNING: No sites found for dataset '{dataset}', using global redirector"
            )
            redirector = {
                "infn": "root://xrootd-cms.infn.it//",
                "fnal": "root://cmsxrootd.fnal.gov/",
                "cern": "root://cms-xrd-global.cern.ch/",
            }
            xrd = redirector[args.redirector]
            if args.limit is not None:
                flist = flist[: args.limit]
            if dsname not in fdict:
                fdict[dsname] = [xrd + f for f in flist if len(f) > 1]
            else:
                fdict[dsname].extend([xrd + f for f in flist if len(f) > 1])
            continue

        # Process sites with careful error handling
        xrootd_tools = check_xrootd_availability()
        xrd = None
        sites_xrootd_prefix = get_xrootd_sites_map()

        # First gather completion information for all sites without checking responsiveness
        site_completions = []
        if args.verbose:
            print(
                f"Evaluating completion for {len(sites)} sites for dataset {dataset}..."
            )

        for site in sites:
            if not site:
                continue

            # Handle site blacklisting/whitelisting
            if args.blacklist_sites is not None and site in args.blacklist_sites:
                if args.verbose:
                    print(f"Site {site} is blacklisted, skipping")
                continue
            if args.whitelist_sites is not None and site not in args.whitelist_sites:
                if args.verbose:
                    print(f"Site {site} is not in whitelist, skipping")
                continue

            # Skip if site has no XRootD access in our map
            if site not in sites_xrootd_prefix:
                if args.verbose:
                    print(f"Site {site} has no XRootD access in our map, skipping")
                continue

            # Only check completion first
            completion = check_site_completion(site, dataset, xrootd_tools)

            # Save all sites with files for sorting
            if completion > 0:
                site_completions.append((site, completion))

        # Sort sites by completion percentage (highest first)
        site_completions.sort(key=lambda x: x[1], reverse=True)

        # Now check responsiveness only for the top sites in order of completion
        site_candidates = []
        if args.verbose:
            print(
                f"Found {len(site_completions)} sites with files, checking responsiveness in order of completion..."
            )

        # Try sites in order of completion until we find a responsive one
        for site, completion in site_completions:
            if args.verbose:
                print(
                    f"Checking responsiveness for site {site} with {completion:.2f}% completion..."
                )

            if (
                xrootd_tools["xrdfs"]
                or xrootd_tools["python_xrootd"]
                or xrootd_tools["gfal"]
            ):
                is_responsive = check_site_responsiveness(
                    site, sites_xrootd_prefix, xrootd_tools, dataset_files=flist
                )
                if is_responsive:
                    site_candidates.append((site, completion))
                    print(
                        f"Found responsive site: {site} with {completion:.2f}% completion"
                    )
                    break  # Stop after finding the first responsive site
            else:
                # If no XRootD tools available, just assume site is responsive
                print(
                    f"⚠️ No XRootD tools available - assuming site {site} is responsive"
                )
                site_candidates.append((site, completion))
                break

        best_site = None
        # If no sites were responsive, try the global redirector
        if not site_candidates:
            print(f"No responsive sites found for {dsname}, using global redirector")

            redirector = {
                "infn": "root://xrootd-cms.infn.it//",
                "fnal": "root://cmsxrootd.fnal.gov/",
                "cern": "root://cms-xrd-global.cern.ch/",
            }
            xrd = redirector[args.redirector]

        else:
            # Use the first (and only) responsive site
            best_site, completion = site_candidates[0]
            if args.verbose:
                print(
                    f"Selected site {best_site} with {completion:.2f}% dataset completion"
                )

            # Get redirector for the best site
            if isinstance(sites_xrootd_prefix[best_site], list):
                xrd = sites_xrootd_prefix[best_site][0]
            elif isinstance(sites_xrootd_prefix[best_site], dict):
                # For dictionary redirectors, we need to pass the whole dictionary
                # to _get_pfn_for_site which knows how to handle regex rules
                # print(f"Using dictionary redirector for {best_site}: {sites_xrootd_prefix[best_site]}")
                xrd = sites_xrootd_prefix[best_site]
            elif isinstance(sites_xrootd_prefix[best_site], str):
                xrd = sites_xrootd_prefix[best_site]
            else:
                print(f"No redirector found for site {best_site}")
                # Fall back to global redirector
                redirector = {
                    "infn": "root://xrootd-cms.infn.it//",
                    "fnal": "root://cmsxrootd.fnal.gov/",
                    "cern": "root://cms-xrd-global.cern.ch/",
                }
                xrd = redirector[args.redirector]

        if args.limit is not None:
            flist = flist[: args.limit]

        # Add files to dictionary
        if dsname not in fdict:

            if best_site in site_url_formats:
                redirector = site_url_formats[best_site]["redirector"]
                if args.verbose:
                    print(f"Using redirector for site {best_site}: {redirector}")

                # Apply the redirector to all files - clean and simple
                paths = [f"{redirector}{f.lstrip('/')}" for f in flist if len(f) > 1]
                if args.verbose:
                    print(f"Sample path: {paths[0]}")
            else:
                # No cached format, use the standard approach
                paths = [_get_pfn_for_site(f, xrd) for f in flist if len(f) > 1]

            validated_paths = [validate_xrootd_path(p) for p in paths]
            fdict[dsname] = validated_paths
        else:

            # Same logic for the else clause
            if best_site and best_site in site_url_formats:
                redirector = site_url_formats[best_site]["redirector"]
                if args.verbose:
                    print(f"Using redirector for site {best_site}: {redirector}")

                # Apply the redirector to all files - clean and simple
                paths = [f"{redirector}{f.lstrip('/')}" for f in flist if len(f) > 1]
                if args.verbose:
                    print(f"Sample path: {paths[0]}")
            else:
                # No cached format, use the standard approach
                paths = [_get_pfn_for_site(f, xrd) for f in flist if len(f) > 1]

            validated_paths = [validate_xrootd_path(p) for p in paths]
            fdict[dsname].extend(validated_paths)

        print(f"Added {len(fdict[dsname])} files for dataset {dsname}")

    # Report on failures
    if failed_datasets:
        print("\n===== SUMMARY OF FAILED DATASETS =====")
        for ds in failed_datasets:
            print(f"- {ds}")
        print(f"Total: {len(failed_datasets)} failed datasets out of {len(fset)}")

    return fdict


def process_single_dataset(dataset, args_dict):
    """Process a single dataset, returning (dsname, files)"""
    import os
    import re
    import uuid
    import time

    # Convert args_dict back to an object
    class Args:
        def __init__(self, **kwargs):
            for key, value in kwargs.items():
                setattr(self, key, value)

    args = Args(**args_dict)

    # Create a TRULY unique temp filename for this thread
    # Using UUID and thread-specific timestamp to avoid collisions
    temp_input_file = f"temp_{os.getpid()}_{uuid.uuid4()}_{time.time()}.txt"

    # Write just this dataset to the temp file
    with open(temp_input_file, "w") as f:
        f.write(dataset + "\n")

    # Temporarily replace input file
    original_input = args.input
    args.input = temp_input_file

    try:
        # Force proper run period extraction from this dataset
        run_pattern = re.search(r"Run20\d\d([A-Z])", dataset)
        expected_run_full = run_pattern.group(1) if run_pattern else None

        # Also extract version information to distinguish between multiple versions
        version_pattern = re.search(r"-([^-]+)-v(\d+)/", dataset)
        version_info = (
            f"_{version_pattern.group(1)}_v{version_pattern.group(2)}"
            if version_pattern
            else ""
        )

        # Process just this one dataset using existing function
        result = getFilesFromDas(args)

        # If we got results, extract and validate
        if result and len(result) > 0:
            dsname = next(iter(result.keys()))
            file_list = result[dsname]

            # Fix dsname if needed, preserving the actual run year and adding version
            if (
                expected_run_full
                and "MuonEG" in dsname
                and not expected_run_full in dsname
            ):
                correct_dsname = f"MuonEG{expected_run_full}{version_info}"
                print(f"⚠️ Fixing incorrect dsname: {dsname} → {correct_dsname}")
                return correct_dsname, file_list
            else:
                # For existing dsnames, append version info if it's a data sample
                if (
                    expected_run_full
                    and "MuonEG" in dsname
                    and not version_info in dsname
                ):
                    correct_dsname = f"{dsname}{version_info}"
                    return correct_dsname, file_list
                return dsname, file_list
        else:
            # Build a dsname from the dataset name as fallback
            try:
                primary_name = dataset.split("/")[1]
                if expected_run_full:
                    dsname = f"{primary_name}Run2024{expected_run_full}"
                else:
                    dsname = primary_name
            except:
                dsname = dataset.replace("/", "_").strip()

            return dsname, []
    finally:
        # Always clean up
        try:
            os.remove(temp_input_file)
        except:
            pass
        args.input = original_input


def getFilesFromPath(args):
    fdict = {}
    fset = []
    with open(args.input) as fp:
        lines = fp.readlines()
        for line in lines:
            if line.startswith("#") or line.strip() == "":
                continue
            if line.startswith("/"):
                print(
                    "You are trying to read files from path, but providing a dataset in DAS:\n",
                    line,
                )
                print("That's not gonna work, so we exit here")
                sys.exit(1)
            ds = line.strip().split()
            print("ds=", ds)
            dataset = ds[0]
            fdict[ds[0]] = getRootFilesFromPath(ds[1], args.limit)

    return fdict


def getTestlist(args):
    fdict = {}
    with open(args.input) as fp:
        lines = fp.readlines()
        for line in lines:
            if line.startswith("#") or line.strip() == "":
                continue
            if not line.endswith("/"):
                line = line + "/"
            if "test" not in line:
                print("You are not getting files in test directory")

            dirs_in_test = os.popen(f"gfal-ls {line}").read().split("\n")
            for s in dirs_in_test:
                if s == "":
                    continue
                print("dataset: ", s)
                fdict[s] = getRootFilesFromPath(line + s, 1)
    return fdict


def getRootFilesFromPath(d, lim=None):
    import subprocess

    if "root://" in d:
        sp = d.split("/")
        siteIP = "/".join(sp[0:4])
        pathToFiles = "/".join(sp[3:]) + "/"
        allfiles = str(
            subprocess.check_output(["xrdfs", siteIP, "ls", pathToFiles]), "utf-8"
        ).split("\n")
    else:
        siteIP = ""
        pathToFiles = d
        allfiles = [
            os.path.join(d, f)
            for i, f in enumerate(os.listdir(d))
            if f.endswith(".root")
        ]

    rootfiles = []
    for file_or_dir in allfiles:
        if file_or_dir == "" or file_or_dir == pathToFiles:
            continue
        file_or_dir = siteIP + file_or_dir
        if file_or_dir.endswith(".root"):
            if lim == None or len(rootfiles) < lim:
                rootfiles.append(file_or_dir)

        elif not "log" in file_or_dir and not file_or_dir[-1] == "/":
            file_or_dir = file_or_dir + "/"
            print("file or dir:", file_or_dir)
            if lim == None:
                rootfiles.extend(getRootFilesFromPath(file_or_dir))
            elif len(rootfiles) < lim:
                rootfiles.extend(
                    getRootFilesFromPath(file_or_dir, lim - len(rootfiles))
                )

    return rootfiles


def validate(file):
    n_tries = 0
    check_path = os.popen(f"gfal-ls {file}").read()
    if check_path == "":
        return f"NotFound: {file}"
    not_able_open = True
    while n_tries <= 5:
        try:
            fin = uproot.open(file)
            not_able_open = False
            return fin["Events"].num_entries
        except:
            print("retries", n_tries, file)
            n_tries += 1
    if not_able_open:
        return f"FailedRetries: {file}"


def remove_bad_files(sample_dict, outname, remove_bad=True):
    from p_tqdm import p_map

    all_invalid = []
    bad_sample_dict = {}
    f = open(f"{outname.replace('.json','')}_file_status.txt", "w")
    for sample in sample_dict.keys():
        _rmap = p_map(
            validate,
            sample_dict[sample],
            num_cpus=int(args.ncpus),
            desc=f"Validating {sample[:20]}...",
        )

        _results = list(_rmap)
        counts = np.sum([r for r in _results if np.isreal(r) and isinstance(r, int)])
        all_invalid += [r for r in _results if isinstance(r, str)]
        if len(all_invalid) > 0:
            bad_sample_dict[sample] = all_invalid
        f.write(f"{sample} Events: {np.sum(counts)}\n")

    if len(all_invalid) == 0:
        f.write("No bad files found!")
    else:
        print(f"Found {len(all_invalid)} bad files.")
        f.write(f"Found {len(all_invalid)} bad files.")
        if remove_bad == True:
            f.write(
                "\n==========================BAD FILES==========================\n "
            )

            for sample in bad_sample_dict.keys():
                for bad_file in bad_sample_dict[sample]:
                    f.write(bad_file + "\n")
                    if bad_file[bad_file.find("root://") :] in sample_dict[sample]:
                        sample_dict[sample].remove(bad_file[bad_file.find("root://") :])

    return sample_dict


def direct_das_query(dataset_name, campaign_pattern):
    """
    Execute a DAS query using the existing run_das_command function
    Returns list of matching datasets or empty list if none found
    """
    # Remove extra asterisks which cause problems
    clean_pattern = campaign_pattern.replace("**", "*")

    # Build query for dataset discovery
    query = f"dataset=/{dataset_name}/*{clean_pattern}*/NANOAOD*"

    if args.verbose:
        print(f"Executing direct DAS query: {query}")

    try:
        # Check if we're in CI environment

        in_ci = determine_execution_mode()

        if not in_ci:
            # For local environment - use direct dasgoclient call
            cmd = f'dasgoclient -query="instance=prod/global {query}"'

            if args.verbose:
                print(f"Local command: {cmd}")
            # Use the already-working run_das_command function instead of os.popen
            output = run_das_command(cmd)
        else:
            # For CI environment - use the two-step approach
            # First query without wildcards in the campaign pattern
            basic_query = f'/cms.cern.ch/common/dasgoclient -query="dataset=/{dataset_name}/*{clean_pattern}*/NANOAOD* instance=prod/global"'
            print(f"CI basic query: {basic_query}")

            # Use the existing run_das_command function
            all_datasets = run_das_command(basic_query)

            # if not all_datasets:
            #    print(f"No datasets found for {dataset_name}")
            #    return []

            # Now filter datasets matching our pattern criteria
            # pattern_core = clean_pattern.replace("*", "")
            # print(f"Filtering with pattern core: '{pattern_core}'")

            # matching_datasets = [
            #    ds for ds in all_datasets
            #    if pattern_core in ds and "NANOAOD" in ds
            # ]

            output = all_datasets

        # Check results
        if output and output[0]:
            # Filter out validation errors
            valid_outputs = [ds for ds in output if not "Validation error" in ds]
            if valid_outputs:
                print(f"Query found {len(valid_outputs)} datasets:")
                if args.verbose:
                    for i, ds in enumerate(valid_outputs[:3]):
                        print(f"  {i+1}: {ds}")
                    if len(valid_outputs) > 3:
                        print(f"  ... and {len(valid_outputs) - 3} more")
                return valid_outputs
            else:
                print("All results were validation errors")
                return []
        else:
            print(f"WARNING: No datasets found matching: {query}")
            return []

    except Exception as e:
        print(f"ERROR: Error executing DAS query: {e}")
        import traceback

        traceback.print_exc()
        return []


def main(args):

    if args.from_workflow:
        for sample in predefined_sample[args.from_workflow].keys():
            if (
                os.path.exists(
                    f"metadata/{args.campaign}/{sample}_{args.campaign}_{args.year}_{args.from_workflow}.json"
                )
                and args.overwrite == False
            ):
                raise Exception(
                    f"metadata/{args.campaign}/{sample}_{args.campaign}_{args.year}_{args.from_workflow}.json exists"
                )
    elif os.path.exists(args.output) and args.overwrite == False:
        raise Exception(f"{args.output} exists")

    ## If you only provide dataset from the dataset name(DAS) or do from_workflow
    if args.from_dataset or args.from_workflow is not None:
        if args.from_dataset:
            f = open(args.input)
            lines = f.readlines()
        else:
            lines = []
            for sample in predefined_sample[args.from_workflow].keys():
                lines += predefined_sample[args.from_workflow][sample]
            args.input = args.from_workflow + "_predef"
            campaign_list = args.DAS_campaign.split(",")

        if args.DAS_campaign is None:
            raise ("Please provide the campaign info when input dataset")
        args.input = args.input + "_DAS_" + args.campaign
        outf = open(args.input, "w")

        for l in lines:
            l = l.replace("\n", "")
            # read campaigns from two inputs
            if args.from_workflow is not None:
                dataset = []
                for campaign in campaign_list:
                    dataset.extend(direct_das_query(l, campaign))

            if not dataset:
                print(
                    f"WARNING: No datasets found for {l} among campaigns {campaign_list}\n"
                )
                continue  # Skip this dataset if query fails

            if len(dataset) > 1:
                while True:
                    print(f"Found multiple datasets for {l}:")
                    for i, d in enumerate(dataset):
                        print(f"  {i+1}: {d}")
                    campaigns = [d.split("/")[2] for d in dataset]
                    campaign_input = input(
                        f"{l} is which campaign? [Enter integer corresponding to above list. Use ',' for multiple]: "
                    )
                    camp_idxs = []
                    for camp_idx in campaign_input.split(","):
                        try:
                            idx = int(camp_idx) - 1
                            campaigns[idx]
                            camp_idxs.append(idx)
                        except:
                            print(f"{camp_idx} is not a valid input. Try again!\n")
                            continue
                    break

                dataset = []
                for idx in camp_idxs:
                    dataset.extend(direct_das_query(l, campaigns[idx]))

                if dataset and not any("ERROR:" in d for d in dataset):
                    for d in dataset:
                        outf.write(d + "\n")
                else:
                    print(f"WARNING: Skipping invalid dataset result for {l}")
                # continue

            elif len(dataset) > 1:
                if args.verbose:
                    print(f"Found multiple datasets for {l}:")
                    for i, d in enumerate(dataset):
                        print(f"  {i+1}: {d}")
                campaigns = [d.split("/")[2] for d in dataset]
                if args.from_workflow is None or dataset[0].endswith("SIM"):
                    args.DAS_campaign = input(
                        f"{l} is which campaign? \n {campaigns} \n"
                    )

                    dataset = direct_das_query(l, args.DAS_campaign)

                    if dataset and not any("ERROR:" in d for d in dataset):
                        outf.write(dataset[0] + "\n")
                    else:
                        print(f"WARNING: Skipping invalid dataset result for {l}")
                else:
                    for d in dataset:
                        if not "ERROR:" in d:
                            outf.write(d + "\n")
                        else:
                            print(f"WARNING: Skipping invalid dataset result: {d}")

            else:
                if dataset and not any("ERROR:" in d for d in dataset):
                    outf.write(dataset[0] + "\n")
                else:
                    print(f"WARNING: Skipping invalid dataset result for {l}")

            print()
        outf.close()
    ## If put the path
    if args.from_path:
        print("do it from path: ")
        fdict = getFilesFromPath(args)
    elif args.testfile:
        fdict = getTestlist(args)
    else:
        # This is where we'll add parallelization
        print(f"Processing datasets with {args.executor} executor")

        # Read dataset list
        fset = []
        with open(args.input) as fp:
            lines = fp.readlines()
            for line in lines:
                if not line.startswith("#") and line.strip() != "":
                    fset.append(line.strip())

        # Choose executor based on args
        if args.executor == "iterative":
            # Process sequentially
            fdict = getFilesFromDas(args)

        elif args.executor == "futures":
            # Process with concurrent.futures
            import concurrent.futures

            # Convert args to dict for serialization
            args_dict = vars(args)
            fdict = {}

            # Read dataset list and deduplicate
            fset = []
            seen_datasets = set()
            with open(args.input) as fp:
                lines = fp.readlines()
                for line in lines:
                    line = line.strip()
                    if not line.startswith("#") and line and line not in seen_datasets:
                        seen_datasets.add(line)
                        fset.append(line)

            print(
                f"Processing {len(fset)} unique datasets with futures executor (workers={args.workers})"
            )

            with concurrent.futures.ThreadPoolExecutor(
                max_workers=args.workers
            ) as executor:
                # Create a mapping of futures to datasets
                future_to_dataset = {}
                for dataset in fset:
                    future = executor.submit(process_single_dataset, dataset, args_dict)
                    future_to_dataset[future] = dataset
                    if args.verbose:
                        print(f"Submitted dataset: {dataset}")

                # Process results as they complete
                for future in concurrent.futures.as_completed(future_to_dataset):
                    dataset = future_to_dataset[future]
                    try:
                        dsname, file_paths = future.result()
                        if dsname and file_paths:
                            if dsname in fdict:
                                print(
                                    f"⚠️ Accumulating files for {dsname} (adding {len(file_paths)} files to existing {len(fdict[dsname])} files)"
                                )
                                fdict[dsname].extend(file_paths)
                            else:
                                fdict[dsname] = file_paths
                            print_dataset_success(dataset, dsname, len(fdict[dsname]))
                        else:
                            print_dataset_empty(dataset)
                    except Exception as exc:
                        print_dataset_failure(dataset, str(exc))
                        import traceback

                        traceback.print_exc()

    # Check the any file lists empty
    empty = True
    for dsname, flist in fdict.items():
        if len(flist) == 0:
            print(f"WARNING: {dsname} is empty!!!!")
            empty = False
    assert empty, "you have empty lists"
    ## Remove files if not exist
    if not args.skipvalidation:
        fdict = remove_bad_files(fdict, args.output, True)  # remove bad files

    ## Create JSONs
    # create according to workflow
    if args.from_workflow:
        os.system(f"mkdir -p metadata/{args.campaign}/")
        for sample in predefined_sample[args.from_workflow].keys():
            reduced_fdict = {}
            for dataset in fdict.keys():
                for s in predefined_sample[args.from_workflow][sample]:

                    if s in dataset:
                        reduced_fdict[dataset] = fdict[dataset]

            with open(
                f"metadata/{args.campaign}/{sample}_{args.campaign}_{args.year}_{args.from_workflow}.json",
                "w",
            ) as fp:

                json.dump(reduced_fdict, fp, indent=4)

    else:
        os.system(f"mkdir -p metadata/{args.campaign}/")
        with open(f"metadata/{args.campaign}/{args.output}", "w") as fp:
            json.dump(fdict, fp, indent=4)
            print(f"The file is saved at: metadata/{args.campaign}/{args.output}")


if __name__ == "__main__":
    print("Starting fetch.py")
    main(args)
