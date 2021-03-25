import os
import json
import argparse
import pprint

fset = []
with open("UL16APVcheck.txt") as fp: 
    lines = fp.readlines() 
    for line in lines: 
        fset.append(line)

fdict = {}

instance = 'prod/phys03'

xrd = 'root://xrootd-cms.infn.it//'

for dataset in fset:
    flist = os.popen(("dasgoclient -query='instance={} file dataset={}'").format(instance,fset[fset.index(dataset)].rstrip())).read().split('\n')
    if "Summer19" in dataset:
        dictname = "ttbar19"
    elif "Summer20" in dataset:
        dictname = "ttbar20"
    else:
        dictname = dataset.rstrip()
    if dictname not in fdict:
        fdict[dictname] = [xrd+f for f in flist if len(f) > 1]
    else: #needed to collect all data samples into one common key "Data" (using append() would introduce a new element for the key)
        fdict[dictname].extend([xrd+f for f in flist if len(f) > 1])

#pprint.pprint(fdict, depth=1)
with open('samplesUL16APV.json', 'w') as fp:
    json.dump(fdict, fp, indent=4)
