import os, sys, json
from copy import deepcopy
from collections import defaultdict
import time
import argparse
from rich import print
from BTVNanoCommissioning.utils.xrootdtools import get_xrootd_sites_map, find_other_file

parser = argparse.ArgumentParser()
parser.add_argument("job_dir", type=str, help="Input string")
parser.add_argument(
    "-u",
    "--update-xrootd",
    action="store_true",
    help="Update xrootd paths of failed jobs.",
)
args = parser.parse_args()

sitemap = get_xrootd_sites_map()
jobdir = args.job_dir
updatexrootd = args.update_xrootd

numlist = [i.strip() for i in open(f"{jobdir}/jobnum_list.txt", "r").readlines()]

with open(f"{jobdir}/arguments.json") as f:
    args = json.load(f)
outdir = args["outputDir"]

toresubmit = []
for n in numlist:
    outfile = f"{outdir}/hists_{n}/hists_{n}.coffea"
    if not os.path.isfile(outfile):
        toresubmit.append(n)

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
jdlnew = open(f"{jobdir}/resubmit.jdl", "w")
for line in jdl:
    towrite = line.replace("jobnum_list.txt", "jobnum_list_resubmit.txt")
    if updatexrootd:
        towrite = towrite.replace("split_samples.json", "split_samples_resubmit.json")
    jdlnew.write(towrite)

print(f"[yellow]Found {len(toresubmit)} missing outputs: {toresubmit}[/]")
print("[b]Resubmit with:[/]")
print(f"condor_submit {jobdir}/resubmit.jdl")
