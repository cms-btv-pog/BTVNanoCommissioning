import os
import json
import argparse
parser = argparse.ArgumentParser(description='Run analysis on baconbits files using processor coffea files')
parser.add_argument('-i', '--input', default=r'singlemuon', help='List of samples in DAS (default: %(default)s)')
parser.add_argument('-s', '--site', default=r'global', help='Site (default: %(default)s)')
parser.add_argument('-o', '--output', default=r'singlemuon', help='Site (default: %(default)s)')
args = parser.parse_args()
fset = []
<<<<<<< HEAD
with open(args.input) as fp: 
=======
with open("tt5.txt") as fp: 
>>>>>>> up/master
    lines = fp.readlines() 
    for line in lines: 
        fset.append(line)

fdict = {}

<<<<<<< HEAD
instance = 'prod/'+args.site
=======
instance = 'prod/global'
>>>>>>> up/master

xrd = 'root://xrootd-cms.infn.it//'

for dataset in fset:
    print(fset)
    flist = os.popen(("/cvmfs/cms.cern.ch/common/dasgoclient -query='instance={} file dataset={}'").format(instance,fset[fset.index(dataset)].rstrip())).read().split('\n')
<<<<<<< HEAD
    
=======
    #if "Summer19" in dataset:
    #    dictname = "ttbar19"
    #elif "Summer20" in dataset:
    #    dictname = "ttbar20"
    #else:
>>>>>>> up/master
    dictname = dataset.rstrip()
    if dictname not in fdict:
        fdict[dictname] = [xrd+f for f in flist if len(f) > 1]
    else: #needed to collect all data samples into one common key "Data" (using append() would introduce a new element for the key)
        fdict[dictname].extend([xrd+f for f in flist if len(f) > 1])

#pprint.pprint(fdict, depth=1)
<<<<<<< HEAD
with open('../metadata/%s.json'%(args.output), 'w') as fp:
=======
with open('../muEG_nano.json', 'w') as fp:
>>>>>>> up/master
    json.dump(fdict, fp, indent=4)
