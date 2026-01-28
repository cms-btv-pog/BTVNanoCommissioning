import os, sys, json
from glob import glob
from coffea.processor import accumulate
from coffea.util import load, save

if len(sys.argv) < 2:
    raise ValueError("Syntax: python condor_lxplus/combinecoffearoutputs.py output_dir")

outputdir = sys.argv[1]

coffeadirs = os.listdir(outputdir)

dict_hist = {}
for cdir in coffeadirs:
    dict_hist[cdir] = load(f"{outputdir}/{cdir}/{cdir}.coffea")

combined_hist = accumulate([dict_hist[key] for key in dict_hist.keys()])

save(combined_hist, f"{outputdir}/combined_hist.coffea")