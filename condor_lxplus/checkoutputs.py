import os, sys, json

if len(sys.argv) < 2:
    raise ValueError("Syntax: python condor_lxplus/checkoutputs.py job_condor_dir")

jobdir = sys.argv[1]

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
    print("All jobs complete. Nothing to resubmit!")
    exit()

numlist2 = open(f"{jobdir}/jobnum_list_resubmit.txt", "w")
for r in toresubmit:
    numlist2.write(f"{r}\n")
numlist2.close()

jdl = open(f"{jobdir}/submit.jdl", "r").readlines()
jdlnew = open(f"{jobdir}/resubmit.jdl", "w")
for line in jdl:
    jdlnew.write(line.replace("jobnum_list.txt", "jobnum_list_resubmit.txt"))

print(f"Found {len(toresubmit)} missing outputs:", toresubmit)
print("Resubmit with:")
print(f"condor_submit {jobdir}/resubmit.jdl")
