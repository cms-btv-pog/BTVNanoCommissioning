import os, sys
from glob import glob

indir = sys.argv[1]
subdirs = [i for i in os.listdir(indir) if i.startswith("arrays_")]
systs = os.listdir(indir)

outfile = open("hadd.sh", "w")

for subdir in subdirs:
    systs = os.listdir(f"{indir}/{subdir}")
    for syst in systs:
        roots = glob(f"{indir}/{subdir}/{syst}/*/*.root")
        if len(roots) == 0:
            print(f"Skipping {indir}/{subdir}/{syst}. Not the right directory structure.")
            continue
        samps = os.listdir(f"{indir}/{subdir}/{syst}")
        for samp in samps:
            if len(glob(f"{indir}/{subdir}/{syst}/{samp}/*.root")) == 0:
                continue
            outfile.write(
                f"hadd -v 0 {indir}/{subdir}/{syst}/{samp}.root {indir}/{subdir}/{syst}/{samp}/*.root\n"
            )

print(
    "Now run `parallel :::: hadd.sh` from an environment with ROOT installed. E.g. \nconda activate rootenv\nparallel :::: hadd.sh\nconda activate btv_coffea"
)
