import os, re, sys
from glob import glob
from coffea.processor import accumulate
from coffea.util import load, save

if len(sys.argv) < 2:
    raise ValueError("Syntax: python condor_lxplus/combinecoffeaoutputs.py output_dir")

outputdir = sys.argv[1].rstrip("/")

# One file per job: <outputdir>/hists_<n>/hists_<n>.coffea . Job directories that
# are empty (job still running or failed) are skipped, and so is any previously
# combined output living in the same directory.
coffealist = []
for cdir in sorted(os.listdir(outputdir)):
    if re.fullmatch(r"hists_\d+", cdir) is None:
        continue
    cfile = f"{outputdir}/{cdir}/{cdir}.coffea"
    if not os.path.isfile(cfile):
        print(f"Skipping {outputdir}/{cdir}, no {cdir}.coffea inside.")
        continue
    coffealist.append(cfile)

if len(coffealist) == 0:
    print("Found zero coffea files. Are you pointing to the right output directory?")
    exit()

extras = [f for f in glob(f"{outputdir}/*/*.coffea") if f not in coffealist]
if len(extras) > 0:
    print(f"\nIgnoring {len(extras)} file(s) not matching hists_<n>/hists_<n>.coffea:")
    for f in extras:
        print(f"  {f}")

print(f"\nCombining {len(coffealist)} coffea files from {outputdir}")
# Fold the files in one at a time. Loading them all up front costs
# n_files x (size of one merged output), which at the current histogram
# granularity is hundreds of GB for a full MC campaign; folding keeps the peak
# at roughly twice the merged result.
combined_hist = load(coffealist[0])
for i, cfile in enumerate(coffealist[1:], 2):
    combined_hist = accumulate([combined_hist, load(cfile)])
    if i % 25 == 0 or i == len(coffealist):
        print(f"  merged {i}/{len(coffealist)}")

outfile = f"{outputdir}/combined_hist.coffea"
save(combined_hist, outfile)
print(f"Datasets in output: {len(combined_hist)}")
print(f"Wrote {outfile}")
