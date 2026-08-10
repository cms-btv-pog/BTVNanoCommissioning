#!/bin/bash -x
set -euo pipefail

JOBID=$1
NCPU=$2

HOME=$(pwd)
export HOME
if [[ -d "/afs/cern.ch/user/${USER:0:1}/$USER" ]]; then
  HOME="/afs/cern.ch/user/${USER:0:1}/$USER"  # crucial on lxplus condor but cannot set on cmsconnect
  export HOME
fi
env


WORKDIR=$(pwd)
rm -f "$WORKDIR/.success"

# Get arguments
declare -A ARGS
for key in workflow output samplejson year campaign isSyst ttbar_reweights isArray noHist overwrite voms chunk skipbadfiles outputDir remoteRepo scm_version; do
    ARGS[$key]=$(jq -r ".$key" "$WORKDIR/arguments.json")
done

# Set up mamba environment
## Interactive bash script with fallback pointing to $HOME, hence setting $PWD of worker node as $HOME
HOME=$(pwd)
export HOME

echo "Setting up mamba environment"
if [[ -n "${ARGS[remoteRepo]}" && "${ARGS[remoteRepo]}" != "null" ]]; then
    echo "remoteRepo is set to ${ARGS[remoteRepo]}"
    installer_downloaded=false
    for i in {1..10}; do
        if wget -L micro.mamba.pm/install.sh; then
            installer_downloaded=true
            break
        fi
        echo "Failed attempt #${i} to download mamba installer. Will retry."
        sleep 30
    done
    if [[ "$installer_downloaded" != true ]]; then
        echo "Failed to download the micromamba installer after 10 attempts." >&2
        exit 1
    fi
    chmod +x install.sh
    ## FIXME parsing arguments does not work. will use defaults in install.sh instead, see https://github.com/mamba-org/micromamba-releases/blob/main/install.sh 
    ## Tried solutions listed in https://stackoverflow.com/questions/14392525/passing-arguments-to-an-interactive-program-non-interactively
    ./install.sh <<< $'bin\nY\nY\nmicromamba\n' 
    # shellcheck source=/dev/null
    source .bashrc
fi

export PATH=$WORKDIR:$PATH

if [[ ! -d "/afs/cern.ch/user/${USER:0:1}/$USER" ]]; then
    ## install necessary packages if on cmsconnect
    micromamba install -c conda-forge jq --yes
fi

# Create base env with python=3.10 and setuptools<=70.1.1
micromamba activate 
micromamba install python=3.10 -c conda-forge xrootd --yes
micromamba activate base
micromamba install setuptools=70.1.1

# Install BTVNanoCommissioning
mkdir BTVNanoCommissioning
cd BTVNanoCommissioning
if [[ ! -f "$WORKDIR/BTVNanoCommissioning.tar.gz" ]]; then
    ## clone the BTVNanoCommissioning repo only, no submodule
    git clone "${ARGS[remoteRepo]}" .
else
    tar xaf "$WORKDIR/BTVNanoCommissioning.tar.gz"
fi

export SETUPTOOLS_SCM_PRETEND_VERSION=${ARGS[scm_version]}
pip install -e .

## other dependencies
pip install psutil

# Build the sample json given the job id
python -c "import json, os; flname = 'split_samples.json' if os.path.isfile(f'$WORKDIR/split_samples.json') else 'split_samples_resubmit.json';  json.dump(json.load(open(f'$WORKDIR/{flname}'))['$JOBID'], open('$WORKDIR/sample.json', 'w'), indent=4)"
cp "$WORKDIR/sample.json" "$WORKDIR/BTVNanoCommissioning/sample.json"

ls -lah "$WORKDIR"
ls -lah "$WORKDIR/BTVNanoCommissioning"

# Unparse arguments and send to runner.py
OPTS=(
    --wf "${ARGS[workflow]}"
    --year "${ARGS[year]}"
    --campaign "${ARGS[campaign]}"
    --chunk "${ARGS[chunk]}"
)
if [ "${ARGS[voms]}" != "null" ]; then
    OPTS+=(--voms "${ARGS[voms]}")
fi
if [ "${ARGS[isSyst]}" != "false" ]; then
    OPTS+=(--isSyst "${ARGS[isSyst]}")
fi
if [ "${ARGS[ttbar_reweights]}" != "none" ]; then
    OPTS+=(--ttbar-reweights "${ARGS[ttbar_reweights]}")
fi
for key in  isArray noHist overwrite skipbadfiles; do
    if [ "${ARGS[$key]}" == true ]; then
        OPTS+=("--$key")
    fi
done
OPTS+=(--output "${ARGS[output]//.coffea/_$JOBID.coffea}")  # add a suffix to output file name
OPTS+=(--json sample.json)  # use the sample json for this JOBID

# Check the number of CPUs requested and set the worker accordingly.
# If nCPU > 1, use futures executor with nCPU workers. If nCPU = 1, use iterative executor with 1 worker.
if [[ "$NCPU" -gt 1 ]]; then
    OPTS+=(--worker "$NCPU")  # use number of worker = nCPU
    OPTS+=(--executor futures)
else
    OPTS+=(--worker 1)  # use number of worker = 1
    OPTS+=(--executor iterative)
fi

# Launch
echo "Now launching: python runner.py ${OPTS[*]}"
python runner.py "${OPTS[@]}"

# Transfer output
if [[ ${ARGS[outputDir]} == root://* ]]; then
    if [[ "${ARGS[noHist]}" != true ]]; then
        xrdcp --silent -p -f -r hists_* "${ARGS[outputDir]}/"
    fi
    if [[ "${ARGS[isArray]}" == true ]]; then
        xrdcp --silent -p -f -r arrays_* "${ARGS[outputDir]}/"
    fi
else
    mkdir -p "${ARGS[outputDir]}"
    if [[ "${ARGS[noHist]}" != true ]]; then
        cp -p -f -r hists_* "${ARGS[outputDir]}/"
    fi
    if [[ "${ARGS[isArray]}" == true ]]; then
        cp -p -f -r arrays_* "${ARGS[outputDir]}/"
    fi
fi

### one can also consider origanizing the root files in the subdirectory structure ###
# for filename in `\ls *.root`; do
#     SAMPLENAME=$(echo "$filename" | sed -E 's/(.*)_f[0-9-]+_[0-9]+\.root/\1/')
#     # SAMPLENAME=$(echo "$filename" | sed -E 's/(.*)_[0-9a-z]{9}-[0-9a-z]{4}-.*\.root/\1/')
#     xrdcp --silent -p -f $filename ${ARGS[outputDir]}/$SAMPLENAME/
# done

touch "$WORKDIR/.success"
