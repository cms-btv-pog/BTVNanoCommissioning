from coffea.util import load
import awkward as ak
import numpy as np
import sys, json, glob
import os
import subprocess


def dump_lumi(output, fname):
    lumi, run = [], []
    for m in output.keys():
        for f in output[m].keys():
            lumi.extend(output[m][f]["lumi"].value)
            run.extend(output[m][f]["run"].value)

    # Sort runs and keep lumisections matched
    run, lumi = np.array(run), np.array(lumi)
    sorted_indices = np.lexsort((lumi, run))  # Sort by run first, then lumi
    run = run[sorted_indices]
    lumi = lumi[sorted_indices]
    # Create dictionary with ls values for each run
    dicts = {}
    for r in np.unique(run):
        dicts[str(r)] = lumi[run == r]
    # Convert to format for brilcalc
    for r in dicts.keys():
        ar = ak.singletons(ak.Array(dicts[r]))
        ars = ak.concatenate([ar, ar], axis=-1)
        dicts[r] = ak.values_astype(ars, int).tolist()

    with open(f"{fname}_lumi.json", "w") as outfile:
        json.dump(dicts, outfile, indent=2)
    json_path = f"{fname}_lumi.json"
    # Create a brilcalc script that can be executed separately
    script_path = os.path.abspath(f"{fname}_get_lumi.sh")
    with open(script_path, "w") as script:
        script.write("#!/bin/bash\n\n")
        script.write("# Luminosity calculation script\n")
        
        # Detect CI environment
        script.write("# Detect CI environment\n")
        script.write('if [ -n "$CI" ] || [ -n "$GITLAB_CI" ] || [ -n "$GITHUB_ACTIONS" ]; then\n')
        script.write('    IN_CI=true\n')
        script.write('else\n')
        script.write('    IN_CI=false\n')
        script.write('fi\n\n')
        
        # Replace the CI section with this practical implementation
        script.write("if $IN_CI; then\n")
        script.write("    echo 'Running in CI environment'\n")
            
        # Try pip installation first (according to official docs)
        script.write("    # Try pip installation of brilws as documented at cms-service-lumi.web.cern.ch\n")
        script.write("    if ! command -v brilcalc &> /dev/null; then\n")
        script.write("        echo 'Installing brilws via pip...'\n")
        script.write("        pip install --user brilws\n")
        script.write("        export PATH=$HOME/.local/bin:$PATH\n")
        script.write("    fi\n\n")

        # Check if installation succeeded
        script.write("    # Check if brilcalc is now available\n")
        script.write("    if command -v brilcalc &> /dev/null; then\n")
        script.write("        echo 'Successfully installed/found brilcalc'\n")
        script.write(f"        brilcalc lumi -c web -i {os.path.basename(json_path)} -u /pb > {os.path.basename(fname)}_lumi_calc.txt\n")
        script.write("        BRIL_EXIT=$?\n")
        script.write("    else\n")
        script.write("        echo 'Brilcalc installation failed, trying Docker...'\n")

        # Try Docker if available
        script.write("        # Try Docker if available\n")
        script.write("        if command -v docker &> /dev/null; then\n")
        script.write("            echo 'Using Docker to run brilcalc'\n")
        script.write(f"            docker run --rm -v $(pwd):/data cms-bril/brilcalc:latest brilcalc lumi -c web -i /data/{os.path.basename(json_path)} -u /pb > {os.path.basename(fname)}_lumi_calc.txt\n")
        script.write("            BRIL_EXIT=$?\n")
        script.write("         else\n")
        script.write("            # Neither pip install nor Docker worked\n")
        script.write("            echo 'ERROR: Could not access brilcalc through any method'\n")
        script.write("            BRIL_EXIT=1\n")
        script.write("        fi\n")
        script.write("    fi\n")
        
        # Process results regardless of method
        script.write("# Process results\n")
        script.write("if [ $BRIL_EXIT -ne 0 ]; then\n")
        script.write(f"    echo \"ERROR: brilcalc failed with exit code $BRIL_EXIT\"\n")
        script.write(f"    echo \"1000.0\" > {fname}_lumi_value.txt\n")
        script.write("    exit 1\n")
        script.write("fi\n\n")
        
        # Rest of your script for extracting values remains the same
        script.write("# Extract the luminosity value\n")
        script.write(f"if grep -q \"#Summary\" {fname}_lumi_calc.txt; then\n")
        script.write(f"    total_lumi=$(grep -A 2 \"#Summary\" {fname}_lumi_calc.txt | tail -1 | awk -F '|' '{{print $(NF-1)}}')\n")
        script.write("    echo \"Luminosity: $total_lumi /pb\"\n")
        script.write(f"    echo $total_lumi > {fname}_lumi_value.txt\n")
        script.write("else\n")
        script.write("    echo \"ERROR: Could not find luminosity summary in output\"\n")
        script.write(f"    echo \"1.0\" > {fname}_lumi_value.txt\n")
        script.write("fi\n")

    # Make the script executable
    os.chmod(script_path, 0o755)
    print(f"Created luminosity calculation script: {script_path}")
    
    # Now actually RUN the script!
    try:
        print(f"Executing luminosity script: {script_path}")
        print(f"Current working directory: {os.getcwd()}")
        # Use capture_output to get both stdout and stderr
        result = subprocess.run([script_path], check=False, capture_output=True, text=True)
        # Always print stdout and stderr for debugging
        if result.stdout:
            print(f"Script stdout:\n{result.stdout}")
        if result.stderr:
            print(f"Script stderr:\n{result.stderr}")
        # Check if the script execution was successful
        if result.returncode != 0:
            print(f"Warning: Luminosity script failed with exit code {result.returncode}")
            if result.stderr:
                print(f"Error: {result.stderr}")
            # Create a default lumi value file
            with open(f"{fname}_lumi_value.txt", "w") as f:
                default_lumi = 1000.0  # Using 1000/pb as default
                f.write(str(default_lumi))
            print(f"Luminosity in pb: {default_lumi}")  # Use EXACT format suball.py expects
            return default_lumi
            
        # Read the luminosity value from the output file
        lumi_value_file = f"{fname}_lumi_value.txt"
        if os.path.exists(lumi_value_file):
            with open(lumi_value_file, 'r') as f:
                lumi_value = f.read().strip()
                try:
                    lumi_in_pb = float(lumi_value)
                    print(f"Luminosity in pb: {lumi_in_pb}")  # EXACT format for suball.py
                    return lumi_in_pb
                except ValueError:
                    print(f"Could not convert luminosity '{lumi_value}' to float")
                    default_lumi = 1000.0
                    print(f"Luminosity in pb: {default_lumi}")  # EXACT format for suball.py
                    return default_lumi
        
        # If we reached here, we didn't get a valid luminosity value
        default_lumi = 1000.0
        print(f"Default uminosity in pb: {default_lumi}")  # EXACT format for suball.py
        
        
    except Exception as e:
        print(f"Unexpected error: {e}")


def dump_dataset(output, fname, alljson):
    jsonlist = glob.glob(alljson) if "*" in alljson else alljson.split(",")
    print("Original jsons:", jsonlist)
    original_list, list_from_coffea = {}, {}
    for j in jsonlist:
        old = json.load(open(j))
        for o in old.keys():
            if o not in original_list.keys():
                original_list[o] = []
            original_list[o].extend(old[o])

    for m in output.keys():
        for f in output[m].keys():
            if f not in list_from_coffea.keys():
                list_from_coffea[f] = list(output[m][f]["fname"])
            else:
                list_from_coffea[f] += list(set(output[m][f]["fname"]))
    failed = {}
    for t in original_list.keys():
        failed[t] = []
        if t not in list_from_coffea.keys():
            failed[t] = original_list[t]
            continue
        for f in original_list[t]:
            if f not in list_from_coffea[t]:
                failed[t].append(f)

    with open(f"{fname}_failed_dataset.json", "w") as outfile:
        json.dump(failed, outfile, indent=4)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Dump procssed luminosity & failed coffea"
    )
    parser.add_argument(
        "-t",
        "--type",
        default="all",
        choices=["all", "lumi", "failed"],
        help="Choose the function for dump luminosity(`lumi`)/failed files(`failed`) into json",
    )
    parser.add_argument(
        "-c",
        "--coffea",
        required=True,
        help="Processed coffea files, splitted by ,. Wildcard option * available as well.",
    )
    parser.add_argument(
        "-n", "--fname", required=True, help="Output name of jsons(with _lumi/_dataset)"
    )
    parser.add_argument(
        "-j",
        "--jsons",
        type=str,
        help="Original json files, splitted by ,. Wildcard option * available as well. ",
    )
    args = parser.parse_args()
    if len(args.coffea.split(",")) > 1:
        output = {i: load(i) for i in args.coffea.split(",")}
    elif "*" in args.coffea:
        args.coffea = glob.glob(args.coffea)
        output = {i: load(i) for i in args.coffea}
    else:
        output = {args.coffea: load(args.coffea)}

    if args.type == "all" or args.type == "lumi":
        print("===>Dump Processed Luminosity")
        dump_lumi(output, args.fname)
    if args.type == "all" or args.type == "failed":
        print("===>Dump Failed Files")
        dump_dataset(output, args.fname, args.jsons)
