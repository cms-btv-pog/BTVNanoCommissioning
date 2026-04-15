import os, sys, json, shutil
from glob import glob
from alive_progress import alive_bar
from rich import print
from concurrent.futures import ThreadPoolExecutor, as_completed

import argparse
parser = argparse.ArgumentParser(description="Create hadd commands and submission scripts.")
parser.add_argument("outputdir", help="Output directory containing array results")
parser.add_argument("--condor", "-c", action="store_true", help="Automatically run the final condor command")
parser.add_argument("-n", "--ncpu", type=int, default=2, help="Number of CPUs for hadd parallelism")
args = parser.parse_args()

outputdir = os.path.abspath(args.outputdir)

arraydirs = [i for i in os.listdir(outputdir) if i.startswith("arrays_")]
systlist = list(set([i.split("/")[-1] for i in glob(f"{outputdir}/arrays_*/*")]))

if "/hadd" not in outputdir:        
    root_files = glob(f"{outputdir}/hadd/**/*.root", recursive=True)
else:
    root_files = glob(f"{outputdir}/**/*.root", recursive=True)
if not root_files:
    check = False
else:
    check = True
    print(f"[cyan]Found[/] [bold cyan]{len(root_files)}[/] [cyan].root files in[/] [bold white]{outputdir}[/][cyan]. Checking if all files are valid.[/]")

if check:
    import uproot
    cmdlist = []    

    hadd_sh = f"{outputdir}/dohadd.sh"
    if os.path.exists(hadd_sh):
        with open(hadd_sh, "r") as f:
            expected_cmds = [line.strip() for line in f if line.strip()]
    else:
        # Fallback to checking existing files if dohadd.sh is missing
        expected_cmds = []
        for filepath in root_files:
            relpath = os.path.relpath(filepath, f"{outputdir}/hadd")
            syst = os.path.dirname(relpath)
            samp = os.path.basename(relpath).replace(".root", "")
            expected_cmds.append(
                f"hadd -j {args.ncpu} -v 0 {filepath} {outputdir}/arrays_*/{syst}/{samp}/*.root"
            )

    def check_root_file(cmd):
        parts = cmd.split(" ")
        filepath = next((p for p in parts if p.endswith(".root")), None)
        if not filepath:
            return None
        
        is_fine = os.path.exists(filepath)
        if is_fine:
            try:
                with uproot.open(filepath) as f:
                    if "Events" not in f and "Event" not in f:
                        is_fine = False
            except Exception:
                is_fine = False

        if not is_fine:
            # Add -f for retry
            return cmd if "hadd -f" in cmd else cmd.replace("hadd ", "hadd -f ")
        return None

    with alive_bar(len(expected_cmds), title="Checking root files") as bar:
        with ThreadPoolExecutor(max_workers=16) as executor:
            futures = {executor.submit(check_root_file, cmd): cmd for cmd in expected_cmds}
            for future in as_completed(futures):
                res = future.result()
                if res:
                    cmdlist.append(res)
                bar()

    if len(cmdlist) == 0:
        print("[green][b]All files are fine![/][/]")
        exit()
    
    print(f"\n[red][b]Found {len(cmdlist)} incomplete/corrupt files![/][/]\n")
    haddfile = f"{outputdir}/dohadd_retry.sh"
    retry_suffix = "_retry"

else:
    cmdlist = []
    def create_cmds(syst):
        samplist = list(
            set([i.split("/")[-1] for i in glob(f"{outputdir}/arrays_*/{syst}/*")])
        )
        newoutdir = f"{outputdir}/hadd/{syst}"
        os.makedirs(newoutdir, exist_ok=True)
        local_cmds = []
        for samp in samplist:
            local_cmds.append(f"hadd -j {args.ncpu} -v 0 {newoutdir}/{samp}.root {outputdir}/arrays_*/{syst}/{samp}/*.root")
        return local_cmds

    with alive_bar(len(systlist), title="Creating hadd commands") as bar:
        with ThreadPoolExecutor(max_workers=16) as executor:
            futures = {executor.submit(create_cmds, syst): syst for syst in systlist}
            for future in as_completed(futures):
                cmdlist.extend(future.result())
                bar()

    if len(cmdlist) == 0:
        print("[red][b]Found zero files to hadd. Are you pointing to the right output directory?[/][/]")
        exit()

    haddfile = f"{outputdir}/dohadd.sh"
    retry_suffix = ""

with open(haddfile, "w") as outfl:
    for cmd in cmdlist:
        outfl.write(cmd + "\n")
        print(f"[yellow]{cmd}[/]")

if not args.condor:
    print("\n[b]Now switch to a root env and run:[/]")
    print(f"[cyan]bash {haddfile}[/]")
    print("or")
    print(f"[cyan]parallel :::: {haddfile}[/]")

# Create condor submission files
current_dir = os.path.dirname(os.path.abspath(__file__)) + "/haddscripts"
os.system("mkdir -p " + current_dir)
suffix = os.path.basename(outputdir.rstrip("/")) + retry_suffix
hadd_sh = f"{current_dir}/hadd_wrapper_{suffix}.sh"
hadd_sub = f"{current_dir}/hadd_{suffix}.sub"
haddfile_copy = f"{current_dir}/{os.path.basename(haddfile).replace('.sh', '')}_{suffix}.sh"
shutil.copy(haddfile, haddfile_copy)

os.system(f"mkdir -p {current_dir}/hadd_logs")

with open(hadd_sh, "w") as f:
    f.write("#!/bin/bash\n")
    f.write("source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc11-opt/setup.sh\n")
    f.write("LINE_NUM=$(($1 + 1))\n")
    f.write("COMMAND=$(sed -n \"${LINE_NUM}p\" $2)\n")
    f.write('echo "Executing: $COMMAND"\n')
    f.write('eval "$COMMAND"\n')
os.chmod(hadd_sh, 0o755)

with open(hadd_sub, "w") as f:
    f.write(f"executable = {hadd_sh}\n")
    f.write(f"arguments = $(ProcID) {haddfile_copy}\n")
    f.write(f"log = {current_dir}/hadd_logs/hadd.log\n")
    f.write(f"output = {current_dir}/hadd_logs/hadd.out_$(ProcID)\n")
    f.write(f"error = {current_dir}/hadd_logs/hadd.err_$(ProcID)\n")
    f.write(f"request_cpus = {args.ncpu}\n")
    f.write("getenv = True\n")
    if retry_suffix == "_retry":
        f.write('+JobFlavour = "workday"\n')
    else:
        f.write('+JobFlavour = "longlunch"\n')
    f.write(f"queue {len(cmdlist)}\n")

if args.condor:
    print("\n[b]Submitting to condor...[/]")
    os.system(f"condor_submit {hadd_sub}")
else:
    print("or submit to condor:")
    print(f"[yellow]condor_submit {hadd_sub}[/]")
