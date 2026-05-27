import os, sys, json, concurrent.futures
from copy import deepcopy
from collections import defaultdict
import time
import argparse
from rich import print
from BTVNanoCommissioning.utils.xrootdtools import get_xrootd_sites_map, find_other_file
from glob import glob
from alive_progress import alive_bar

parser = argparse.ArgumentParser()
parser.add_argument("job_dir", type=str, help="Input string")
parser.add_argument(
    "-u",
    "--update-xrootd",
    action="store_true",
    help="Update xrootd paths of failed jobs.",
)
parser.add_argument(
    "-c",
    "--condor",
    action="store_true",
    help="Automatically run the final condor command",
)
args = parser.parse_args()

sitemap = get_xrootd_sites_map()
jobdir = args.job_dir
updatexrootd = args.update_xrootd
run_condor = args.condor

uid = os.getuid()
homedir = os.getenv("HOME")
expected_value = f"{homedir}/x509up_u{uid}"
current_value = os.getenv("X509_USER_PROXY")
if current_value != expected_value:
    print("[red][b]X509_USER_PROXY is NOT set correctly.[/][/]")
    print(f"Please run the following command in your shell:")
    print(f"export X509_USER_PROXY=$HOME/x509up_u`id -u`")
    sys.exit(1)

numlist = [i.strip() for i in open(f"{jobdir}/jobnum_list.txt", "r").readlines()]

with open(f"{jobdir}/arguments.json") as f:
    args = json.load(f)
outdir = args["outputDir"]


def check_job(n, outdir):
    outfile = f"{outdir}/hists_{n}/hists_{n}.coffea"
    if not os.path.isfile(outfile):
        return n, f"[red]Job {n}[/] does not have an output histogram file."
    arraylist = glob(f"{outdir}/arrays_hists_{n}/*/*/*.root")
    for rootfile in arraylist:
        try:
            filesize = os.path.getsize(rootfile)
            if filesize < 10 * 1024:  # Files smaller than 10 kB are likely zombies
                return n, f"[red]Job {n}[/] has at least one zombie root output file."
        except OSError:
            return n, f"[red]Job {n}[/] error accessing root output file."
    return None, None


toresubmit = []

with alive_bar(len(numlist), title="Checking jobs") as bar:
    with concurrent.futures.ThreadPoolExecutor(max_workers=20) as executor:
        futures = {executor.submit(check_job, n, outdir): n for n in numlist}
        for future in concurrent.futures.as_completed(futures):
            res, msg = future.result()
            if res:
                print(msg)
                toresubmit.append(res)
            bar()

toresubmit = sorted(list(set(toresubmit)), key=lambda x: int(x))

if len(toresubmit) == 0:
    print("[green][b]All jobs complete. Nothing to resubmit![/][/]\n")
    print("[b]To hadd, run:[/]")
    print(f"[yellow]python condor_lxplus/haddoutputs.py {outdir}[/]")
    exit()

if updatexrootd:
    with open(f"{jobdir}/split_samples.json") as f:
        samples = json.load(f)
    newsamples = deepcopy(samples)
    for r in toresubmit:
        sampdict = samples[r]
        for samp, fllist in sampdict.items():
            newfllist = []
            for fl in fllist:
                newfl = find_other_file(fl, sitemap)
                newfllist.append(newfl)
                if newfl == fl:
                    print(
                        f"[red][b]WARNING:[/] Could not find replacement for {fl}.[/]"
                    )
                else:
                    print(f"[yellow]{fl}[/] [cyan]->[/] [green]{newfl}[/].")
            newsamples[r][samp] = newfllist
    with open(os.path.join(jobdir, "split_samples_resubmit.json"), "w") as json_file:
        json.dump(newsamples, json_file, indent=4)


numlist2 = open(f"{jobdir}/jobnum_list_resubmit.txt", "w")
for r in toresubmit:
    numlist2.write(f"{r}\n")
numlist2.close()

jdl = open(f"{jobdir}/submit.jdl", "r").readlines()

queues = [
    "espresso",
    "microcentury",
    "longlunch",
    "workday",
    "tomorrow",
    "testmatch",
    "nextweek",
]
current_flavour = None
for line in jdl:
    if "+JobFlavour" in line:
        current_flavour = line.split("=")[1].strip().replace('"', "")
        break

next_flavour = current_flavour
if current_flavour in queues:
    idx = queues.index(current_flavour)
    if idx < len(queues) - 1:
        next_flavour = queues[idx + 1]
        print(
            f"[yellow]Bumping JobFlavour: [red]{current_flavour}[/] -> [green]{next_flavour}[/][/]"
        )

with open(f"{jobdir}/resubmit.jdl", "w") as jdlnew:
    for line in jdl:
        towrite = line.replace("jobnum_list.txt", "jobnum_list_resubmit.txt")
        if updatexrootd:
            towrite = towrite.replace(
                "split_samples.json", "split_samples_resubmit.json"
            )
        if "+JobFlavour" in line and current_flavour:
            towrite = towrite.replace(f'"{current_flavour}"', f'"{next_flavour}"')
        jdlnew.write(towrite)

print(f"[yellow]Found {len(toresubmit)} missing outputs: {toresubmit}[/]")
if run_condor:
    print("\n[b]Submitting to condor...[/]")
    os.system(f"condor_submit {jobdir}/resubmit.jdl")
else:
    print("[b]Resubmit with:[/]")
    print(f"condor_submit {jobdir}/resubmit.jdl")
