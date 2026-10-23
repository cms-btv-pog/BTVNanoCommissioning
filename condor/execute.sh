#!/bin/bash -x

JOBID=$1
NCPU=$2

export HOME=`pwd`
if [ -d /afs/cern.ch/user/${USER:0:1}/$USER ]; then
  export HOME=/afs/cern.ch/user/${USER:0:1}/$USER  # crucial on lxplus condor but cannot set on cmsconnect
fi
env


WORKDIR=`pwd`

# Get arguments
declare -A ARGS
for key in workflow output samplejson year campaign isSyst ttbar_reweights isArray array_systs array_perjet noHist overwrite voms chunk skipbadfiles outputDir remoteRepo; do
    echo $(jq -r ".$key" $WORKDIR/arguments.json)
    ARGS[$key]=$(jq -r ".$key" $WORKDIR/arguments.json)
done

# Set up mamba environment
## Interactive bash script with fallback pointing to $HOME, hence setting $PWD of worker node as $HOME
export HOME=`pwd`

echo "Setting up mamba environment"
if [[ ${ARGS[remoteRepo]} != "" ]]; then
    echo "remoteRepo is set to ${ARGS[remoteRepo]}"
    for i in {1..10}; do
        wget -L micro.mamba.pm/install.sh
        if [ $? -eq 0 ]; then
            break
        fi
	echo "Failed attempt #"$i" to download mamba installer. Will retry."
        sleep 30
    done
    chmod +x install.sh
    ## FIXME parsing arguments does not work. will use defaults in install.sh instead, see https://github.com/mamba-org/micromamba-releases/blob/main/install.sh 
    ## Tried solutions listed in https://stackoverflow.com/questions/14392525/passing-arguments-to-an-interactive-program-non-interactively
    ./install.sh <<< $'bin\nY\nY\nmicromamba\n' 
    source .bashrc
fi

export PATH=$WORKDIR:$PATH

if [ ! -d /afs/cern.ch/user/${USER:0:1}/$USER ]; then
    ## install necessary packages if on cmsconnect
    micromamba install -c conda-forge jq --yes
fi

# Create base env with Python 3.12 (the version validated against coffea 2026;
# coffea 2026 requires >=3.10). setuptools is installed unpinned -- the old
# setuptools<=70.1.1 cap was a coffea 0.7 constraint and is no longer needed.
micromamba activate
micromamba install python=3.12 -c conda-forge xrootd setuptools --yes
micromamba activate base

# Install BTVNanoCommissioning
mkdir BTVNanoCommissioning
cd BTVNanoCommissioning
if [ ! -f $WORKDIR/BTVNanoCommissioning.tar.gz ]; then
    ## clone the BTVNanoCommissioning repo only, no submodule
    git clone ${ARGS[remoteRepo]} .
else
    tar xaf $WORKDIR/BTVNanoCommissioning.tar.gz
fi

# setuptools_scm infers the version from git metadata, which the sandbox tarball
# deliberately does not carry (shipping .git would add ~100 MB to every job).
# Hand it the version recorded when the tarball was made; without this the
# editable install fails with "unable to detect version".
SCM_VERSION=$(sed -n "s/^__version__ = version = '\(.*\)'$/\1/p" src/BTVNanoCommissioning/version.py 2>/dev/null | head -1)
export SETUPTOOLS_SCM_PRETEND_VERSION=${SCM_VERSION:-0.0.0}
echo "Using SETUPTOOLS_SCM_PRETEND_VERSION=$SETUPTOOLS_SCM_PRETEND_VERSION"

pip install -e .

## other dependencies
pip install psutil

# Build the sample json given the job id
python -c "import json, os; flname = 'split_samples.json' if os.path.isfile(f'$WORKDIR/split_samples.json') else 'split_samples_resubmit.json';  json.dump(json.load(open(f'$WORKDIR/{flname}'))['$JOBID'], open('$WORKDIR/sample.json', 'w'), indent=4)"
cp $WORKDIR/sample.json $WORKDIR/BTVNanoCommissioning/sample.json

ls -lah $WORKDIR
ls -lah $WORKDIR/BTVNanoCommissioning

# Unparse arguments and send to runner.py
OPTS="--wf ${ARGS[workflow]} --year ${ARGS[year]} --campaign ${ARGS[campaign]} --chunk ${ARGS[chunk]}"
if [ "${ARGS[voms]}" != "null" ]; then
    OPTS="$OPTS --voms ${ARGS[voms]}"
fi
if [ "${ARGS[isSyst]}" != "false" ]; then
    OPTS="$OPTS --isSyst ${ARGS[isSyst]}"
fi
if [ "${ARGS[ttbar_reweights]}" != "none" ]; then
    OPTS="$OPTS --ttbar-reweights ${ARGS[ttbar_reweights]}"
fi
# Forward every explicit value, "none" included: runner.py's default is not
# "none", so dropping it here would silently turn the trees' systematics back on.
if [ "${ARGS[array_systs]}" != "null" ]; then
    OPTS="$OPTS --array-systs ${ARGS[array_systs]}"
fi
if [ "${ARGS[array_perjet]}" == true ]; then
    OPTS="$OPTS --array-perjet"
fi
for key in  isArray noHist overwrite skipbadfiles; do
    if [ "${ARGS[$key]}" == true ]; then
        OPTS="$OPTS --$key"
    fi
done
OPTS="$OPTS --output ${ARGS[output]//.coffea/_$JOBID.coffea}"  # add a suffix to output file name
OPTS="$OPTS --json sample.json"  # use the sample json for this JOBID

# Check the number of CPUs requested and set the worker accordingly.
# If nCPU > 1, use futures executor with nCPU workers. If nCPU = 1, use iterative executor with 1 worker.
if [ $NCPU -gt 1 ]; then
    OPTS="$OPTS --workers $NCPU"  # use number of worker = nCPU
    OPTS="$OPTS --executor futures"
else
    OPTS="$OPTS --workers 1"  # use number of worker = 1
    OPTS="$OPTS --executor iterative"
fi

# Launch
echo "Now launching: python runner.py $OPTS"
python runner.py $OPTS
RUNNER_STATUS=$?

# Transfer output
COPY_STATUS=0
if [[ ${ARGS[outputDir]} == root://* ]]; then

    xrdcp --silent -p -f -r hists_* ${ARGS[outputDir]}/ || COPY_STATUS=$?
    if [[ "$OPTS" == *"isArray"* ]]; then
	xrdcp --silent -p -f -r arrays_* ${ARGS[outputDir]}/ || COPY_STATUS=$?
    fi
else
    mkdir -p ${ARGS[outputDir]}
    cp -p -f -r hists_* ${ARGS[outputDir]}/ || COPY_STATUS=$?
    if [[ "$OPTS" == *"isArray"* ]]; then
	cp -p -f -r arrays_* ${ARGS[outputDir]}/ || COPY_STATUS=$?
    fi
fi

### one can also consider origanizing the root files in the subdirectory structure ###
# for filename in `\ls *.root`; do
#     SAMPLENAME=$(echo "$filename" | sed -E 's/(.*)_f[0-9-]+_[0-9]+\.root/\1/')
#     # SAMPLENAME=$(echo "$filename" | sed -E 's/(.*)_[0-9a-z]{9}-[0-9a-z]{4}-.*\.root/\1/')
#     xrdcp --silent -p -f $filename ${ARGS[outputDir]}/$SAMPLENAME/
# done

# `.success` is what condor transfers back, so it is always created; but the
# job's exit status must still reflect reality. Exiting 0 unconditionally meant a
# failed runner -- or a worker with no EOS mount, where the copy above fails --
# was reported to condor as a success, so `max_retries` never fired and the loss
# only ever showed up as a missing output file.
touch $WORKDIR/.success

if [ $RUNNER_STATUS -ne 0 ] || [ $COPY_STATUS -ne 0 ]; then
    echo "JOB FAILED: runner exit=$RUNNER_STATUS, output copy exit=$COPY_STATUS"
    exit 1
fi
exit 0
