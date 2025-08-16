#!/bin/bash -xe

JOBID=$1
COMMDIR=$2

export HOME=`pwd`
if [ -d /afs/cern.ch/user/${USER:0:1}/$USER ]; then
  export HOME=/afs/cern.ch/user/${USER:0:1}/$USER  # crucial on lxplus condor but cannot set on cmsconnect
fi
env

WORKDIR=`pwd`

cd $COMMDIR
export X509_USER_PROXY=$COMMDIR/proxy
voms-proxy-info

export PATH="$4:$PATH" 

# Build the sample json given the job id
python -c "import json, os; flname = 'split_samples.json' if os.path.isfile(f'$WORKDIR/split_samples.json') else 'split_samples_resubmit.json';  json.dump(json.load(open(f'$WORKDIR/{flname}'))['$JOBID'], open('$WORKDIR/sample.json', 'w'), indent=4)"

declare -A ARGS
for key in workflow output samplejson year campaign isSyst isArray noHist overwrite voms chunk skipbadfiles outputDir remoteRepo; do
    ARGS[$key]=$(jq -r ".$key" $WORKDIR/arguments.json)
done

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
OPTS="$OPTS --json $WORKDIR/sample.json"  # use the sample json for this JOBID
OPTS="$OPTS --worker 1"  # use number of worker = 1
OPTS="$OPTS --executor iterative --overwrite --outputdir $3"

# Launch
echo "Now launching: python runner.py $OPTS"
python runner.py $OPTS

touch $WORKDIR/.success

