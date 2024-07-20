import os, sys
from glob import glob

indir = sys.argv[1]

systs = os.listdir(indir)

outfile = open("hadd.sh", "w")

for syst in systs:
    roots = glob(f"{indir}/{syst}/*/*.root")
    if len(roots) == 0:
        print(f"Skipping {indir}/{syst}. Not the right directory structure.")
        continue
    samps = os.listdir(f"{indir}/{syst}")
    for samp in samps:
        if len(glob(f"{indir}/{syst}/{samp}/*.root")) == 0:
            continue
        outfile.write(
            f"hadd -v 0 {indir}/{syst}/{samp}.root {indir}/{syst}/{samp}/*.root\n"
        )

print(
    "Now run `parallel :::: hadd.sh` from an environment with ROOT installed. E.g. \nconda activate rootenv\nparallel :::: hadd.sh\nconda activate btv_coffea"
)
