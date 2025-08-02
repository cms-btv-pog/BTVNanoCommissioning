import os, sys, json
from glob import glob

if len(sys.argv) < 2:
    raise ValueError("Syntax: python condor_lxplus/haddoutputs.py output_dir")

outputdir = sys.argv[1]

arraydirs = [i for i in os.listdir(outputdir) if i.startswith("arrays_")]
systlist = list(set([i.split("/")[-1] for i in glob(f"{outputdir}/arrays_*/*")]))

cmdlist = []
for syst in systlist:
    samplist = list(
        set([i.split("/")[-1] for i in glob(f"{outputdir}/arrays_*/{syst}/*")])
    )
    newoutdir = f"{outputdir}/hadd/{syst}"
    os.system("mkdir -p " + newoutdir)
    for samp in samplist:
        cmd = f"hadd {newoutdir}/{samp}.root {outputdir}/arrays_*/{syst}/{samp}/*.root"
        cmdlist.append(cmd)

if len(cmdlist) == 0:
    print("Found zero files to hadd. Are you pointing to the right output directory?")
    exit()

haddfile = f"{outputdir}/dohadd.sh"
print()
with open(haddfile, "w") as outfl:
    for cmd in cmdlist:
        outfl.write(cmd + "\n")
        print(cmd)

print("\nNow switch to a root env and run:")
print("bash", haddfile)
print("or")
print(f"parallel :::: {haddfile}")
