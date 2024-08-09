#!/bin/bash -xe

JOBID=$1

export HOME=`pwd`
if [ -d /afs/cern.ch/user/${USER:0:1}/$USER ]; then
  export HOME=/afs/cern.ch/user/${USER:0:1}/$USER  # crucial on lxplus condor but cannot set on cmsconnect
fi
env


WORKDIR=`pwd`


# Set up mamba environment
## Interactive bash script with fallback pointing to $HOME, hence setting $PWD of worker node as $HOME
export HOME=`pwd`
wget -L micro.mamba.pm/install.sh
chmod +x install.sh
## FIXME parsing arguments does not work. will use defaults in install.sh instead, see https://github.com/mamba-org/micromamba-releases/blob/main/install.sh 
## Tried solutions listed in https://stackoverflow.com/questions/14392525/passing-arguments-to-an-interactive-program-non-interactively
./install.sh <<< $'bin\nY\nY\nmicromamba\n' 
source .bashrc

export PATH=$WORKDIR:$PATH


if [ ! -d /afs/cern.ch/user/${USER:0:1}/$USER ]; then
    ## install necessary packages if on cmsconnect
    micromamba install -c conda-forge jq --yes
fi

# Get arguments
declare -A ARGS
for key in workflow output samplejson year campaign isSyst isArray noHist overwrite voms chunk skipbadfiles outputXrootdDir remoteRepo; do
    ARGS[$key]=$(jq -r ".$key" $WORKDIR/arguments.json)
done

# Create base env with python=3.10 and setuptools<=70.1.1
micromamba activate 
micromamba install python=3.10 -c conda-forge xrootd --yes
micromamba activate base
micromamba install setuptools=70.1.1

# Install BTVNanoCommissioning
mkdir BTVNanoCommissioning
cd BTVNanoCommissioning
if [ ! -f $WORKDIR/BTVNanoCommissioning.tar.gz ]; then
    ## clone the BTVNanoCommissioning repo only, no submodule
    git clone ${ARGS[remoteRepo]} .
else
    tar xaf $WORKDIR/BTVNanoCommissioning.tar.gz
fi
rm -rf src/BTVNanoCommissioning/jsonpog-integration
ln -s /cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration src/BTVNanoCommissioning/jsonpog-integration  # link jsonpog-integration

pip install -e .

## other dependencies
pip install psutil

# Build the sample json given the job id
python -c "import json; json.dump(json.load(open('$WORKDIR/split_samples.json'))['$JOBID'], open('sample.json', 'w'), indent=4)"

# Unparse arguments and send to runner.py
OPTS="--wf ${ARGS[workflow]} --year ${ARGS[year]} --campaign ${ARGS[campaign]} --chunk ${ARGS[chunk]}"
if [ "${ARGS[voms]}" != "null" ]; then
    OPTS="$OPTS --voms ${ARGS[voms]}"
fi
if [ "${ARGS[isSyst]}" != "false" ]; then
    OPTS="$OPTS --isSyst ${ARGS[isSyst]}"
fi

for key in  isArray noHist overwrite skipbadfiles; do
    if [ "${ARGS[$key]}" == true ]; then
        OPTS="$OPTS --$key"
    fi
done
OPTS="$OPTS --output ${ARGS[output]//.coffea/_$JOBID.coffea}"  # add a suffix to output file name
OPTS="$OPTS --json sample.json"  # use the sample json for this JOBID
OPTS="$OPTS --worker 1"  # use number of worker = 1
OPTS="$OPTS --executor iterative"

# Launch
echo "Now launching: python runner.py $OPTS"
python runner.py $OPTS

# Transfer output
if [[ ${ARGS[outputXrootdDir]} == root://* ]]; then

    xrdcp --silent -p -f -r hists_* ${ARGS[outputXrootdDir]}/
    if [[ "$OPTS" == *"isArray"* ]]; then
	xrdcp --silent -p -f -r arrays_* ${ARGS[outputXrootdDir]}/
    fi
else
    mkdir -p ${ARGS[outputXrootdDir]}
    cp -p -f -r hists_* ${ARGS[outputXrootdDir]}/
    if [[ "$OPTS" == *"isArray"* ]]; then
	cp -p -f -r arrays_* ${ARGS[outputXrootdDir]}/
    fi
fi

### one can also consider origanizing the root files in the subdirectory structure ###
# for filename in `\ls *.root`; do
#     SAMPLENAME=$(echo "$filename" | sed -E 's/(.*)_f[0-9-]+_[0-9]+\.root/\1/')
#     # SAMPLENAME=$(echo "$filename" | sed -E 's/(.*)_[0-9a-z]{9}-[0-9a-z]{4}-.*\.root/\1/')
#     xrdcp --silent -p -f $filename ${ARGS[outputXrootdDir]}/$SAMPLENAME/
# done

touch $WORKDIR/.success
