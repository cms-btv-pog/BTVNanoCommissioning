import os
import json
import argparse
import pprint

fset = []
with open("list.txt") as fp: 
    lines = fp.readlines() 
    for line in lines: 
        fset.append(line)

fdict = {}

instance = 'prod/phys03'

xrd = 'root://xrootd-cms.infn.it//'

for dataset in fset:
    flist = os.popen(("dasgoclient -query='instance={} file dataset={}'").format(instance,fset[fset.index(dataset)].rstrip())).read().split('\n')
    if "TTToSemi" in dataset:
        dictname = "ttbarSL"
    if "TTTo2L2Nu" in dataset:
        dictname = "ttbarDL"
    elif "DY" in dataset:#TODO split in HT bins (xsecs?)
        dictname = "DYJets"
    elif "WJets" in dataset:#TODO split in HT bins (xsecs?)
        dictname = "WJets"
    elif "ST_t-channel_top" in dataset:
        dictname = "st_tch_top"
    elif "ST_t-channel_antitop" in dataset:
        dictname = "st_tch_atop"
    elif "ST_tW_top" in dataset:
        dictname = "st_tw_top"
    elif "ST_tW_antitop" in dataset:
        dictname = "st_tw_atop"
    elif "QCD" in dataset:#TODO split in HT or pT bins (xsecs?)
        dictname = "QCD"
    elif "Run201" in dataset: #all data samples have this snippet in common
        dictname = "Data"
    else:
        dictname = dataset.rstrip()
    if dictname not in fdict:
        fdict[dictname] = [xrd+f for f in flist if len(f) > 1]
    else: #needed to collect all data samples into one common key "Data" (using append() would introduce a new element for the key)
        fdict[dictname].extend([xrd+f for f in flist if len(f) > 1])

#pprint.pprint(fdict, depth=1)
with open('datasets.json', 'w') as fp:
    json.dump(fdict, fp, indent=4)
