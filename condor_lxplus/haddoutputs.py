import os, sys, json, shutil, subprocess
from glob import glob
from alive_progress import alive_bar
from rich import print
from concurrent.futures import ThreadPoolExecutor, as_completed

import argparse
import dashboard as job_dashboard

parser = argparse.ArgumentParser(
    description="Create hadd commands and submission scripts."
)
parser.add_argument(
    "input",
    help="Output directory containing array results, or the jobs_condor_* directory",
)
parser.add_argument(
    "--condor",
    "-c",
    action="store_true",
    help="Automatically run the final condor command",
)
parser.add_argument(
    "-n", "--ncpu", type=int, default=2, help="Number of CPUs for hadd parallelism"
)
args = parser.parse_args()

input_path = os.path.abspath(args.input)
jobdir = None
outputdir = input_path

if os.path.isdir(input_path) and os.path.basename(input_path).startswith(
    ("jobs_", "job_")
):
    jobdir = input_path
    status = job_dashboard.load_status(jobdir, create=False)
    outputdir = status.get("output_dir") or ""
    if not outputdir:
        args_path = os.path.join(jobdir, "arguments.json")
        if os.path.isfile(args_path):
            with open(args_path) as f:
                outputdir = json.load(f).get("outputDir", "")
    if not outputdir:
        raise ValueError(
            f"Could not determine output directory from {jobdir}/job_status.json or arguments.json"
        )
    outputdir = os.path.abspath(outputdir)
else:
    jobdir = job_dashboard.find_job_dir_by_output(outputdir)

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
    print(
        f"[cyan]Found[/] [bold cyan]{len(root_files)}[/] [cyan].root files in[/] [bold white]{outputdir}[/][cyan]. Checking if all files are valid.[/]"
    )

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
            futures = {
                executor.submit(check_root_file, cmd): cmd for cmd in expected_cmds
            }
            for future in as_completed(futures):
                res = future.result()
                if res:
                    cmdlist.append(res)
                bar()

    if len(cmdlist) == 0:
        if jobdir:
            job_dashboard.record_hadd_checked(jobdir, len(expected_cmds), 0)
        print("[green][b]All files are fine![/][/]")
        exit()

    print(f"\n[red][b]Found {len(cmdlist)} incomplete/corrupt files![/][/]\n")
    if jobdir:
        job_dashboard.record_hadd_checked(jobdir, len(expected_cmds), len(cmdlist))
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
            local_cmds.append(
                f"hadd -j {args.ncpu} -v 0 {newoutdir}/{samp}.root {outputdir}/arrays_*/{syst}/{samp}/*.root"
            )
        return local_cmds

    with alive_bar(len(systlist), title="Creating hadd commands") as bar:
        with ThreadPoolExecutor(max_workers=16) as executor:
            futures = {executor.submit(create_cmds, syst): syst for syst in systlist}
            for future in as_completed(futures):
                cmdlist.extend(future.result())
                bar()

    if len(cmdlist) == 0:
        print(
            "[red][b]Found zero files to hadd. Are you pointing to the right output directory?[/][/]"
        )
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
    sys.exit(0)

# Create condor submission files
if not jobdir:
    raise ValueError(
        "Could not identify the jobs_condor_* directory for this output directory. "
        "Pass the job directory directly, or make sure its status/arguments JSON points to this output directory."
    )
current_dir = os.path.join(os.path.abspath(jobdir), "haddscripts")
os.makedirs(current_dir, exist_ok=True)
suffix = os.path.basename(outputdir.rstrip("/")) + retry_suffix
hadd_sh = f"{current_dir}/hadd_wrapper_{suffix}.sh"
hadd_sub = f"{current_dir}/hadd_{suffix}.sub"
haddfile_copy = (
    f"{current_dir}/{os.path.basename(haddfile).replace('.sh', '')}_{suffix}.sh"
)
shutil.copy(haddfile, haddfile_copy)

hadd_log_dir = f"{current_dir}/hadd_logs"
os.makedirs(hadd_log_dir, exist_ok=True)

with open(hadd_sh, "w") as f:
    f.write("#!/bin/bash\n")
    f.write(
        "source /cvmfs/sft.cern.ch/lcg/views/LCG_105/x86_64-el9-gcc11-opt/setup.sh\n"
    )
    f.write("LINE_NUM=$(($1 + 1))\n")
    f.write('COMMAND=$(sed -n "${LINE_NUM}p" $2)\n')
    f.write('echo "Executing: $COMMAND"\n')
    f.write('eval "$COMMAND"\n')
os.chmod(hadd_sh, 0o755)

with open(hadd_sub, "w") as f:
    f.write(f"executable = {hadd_sh}\n")
    f.write(f"arguments = $(ProcID) {haddfile_copy}\n")
    f.write(f"log = {hadd_log_dir}/hadd.log_$(Cluster)\n")
    f.write(f"output = {hadd_log_dir}/hadd.out_$(Cluster)-$(ProcID)\n")
    f.write(f"error = {hadd_log_dir}/hadd.err_$(Cluster)-$(ProcID)\n")
    f.write(f"request_cpus = {args.ncpu}\n")
    f.write("getenv = True\n")
    if retry_suffix == "_retry":
        f.write('+JobFlavour = "workday"\n')
    else:
        f.write('+JobFlavour = "longlunch"\n')
    f.write(f"queue {len(cmdlist)}\n")

print("\n[b]Submitting to condor...[/]")
completed = subprocess.run(
    ["condor_submit", hadd_sub],
    capture_output=True,
    text=True,
    check=False,
)
submit_output = "\n".join(
    part for part in [completed.stdout.strip(), completed.stderr.strip()] if part
)
if submit_output:
    print(submit_output)
if completed.returncode != 0:
    raise SystemExit(f"condor_submit failed with return code {completed.returncode}")
submission = job_dashboard.parse_condor_submit_output(submit_output)
if submission is None:
    raise SystemExit(
        "condor_submit returned success but did not confirm a submitted cluster"
    )
if submission["count"] != len(cmdlist):
    raise SystemExit(
        "condor_submit reported "
        f"{submission['count']} submitted jobs; expected {len(cmdlist)}"
    )
job_dashboard.record_hadd_submitted(
    jobdir,
    len(cmdlist),
    cluster_id=submission["cluster_id"],
)
