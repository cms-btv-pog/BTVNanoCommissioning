#!/bin/bash

# Check if the arguments file exists
if [ ! -f "./scripts/submit/arguments.txt" ]; then
    echo "Error: arguments.txt not found."
    return 1
fi

# Read arguments from the file
args=$(<./scripts/submit/arguments.txt)
# Split the arguments into an array
IFS=$'\n' read -d '' -r -a arg_array <<< "$args"

# Check if enough arguments are provided
if [ "${#arg_array[@]}" -lt 4 ]; then
    echo "Error: Insufficient arguments provided in arguments.txt."
    exit 1
fi

# Assign arguments to variables
echo "Running with arguments: ${arg_array[@]}"
arg_campaign="${arg_array[0]}"
arg_year="${arg_array[1]}"
arg_executor="${arg_array[2]}"
arg_lumi="${arg_array[3]}"
scaleout="${arg_array[4]}"
overwrite="${arg_array[5]}"

chunksize=200000
limit=10
max=10

# check files
basepath="/net/data_cms3a-1/fzinn/BTV/btag_sf/nobackup"

datafile="hists_btag_ttbar_sf_data_"$arg_campaign"_"$arg_year"_em_BTV_Run3_"$arg_year"_Comm_MINIAODv4_NanoV12/hists_btag_ttbar_sf_data_"$arg_campaign"_"$arg_year"_em_BTV_Run3_"$arg_year"_Comm_MINIAODv4_NanoV12.coffea"
mcfile="hists_btag_ttbar_sf_MC_"$arg_campaign"_"$arg_year"_BTV_Run3_"$arg_year"_Comm_MINIAODv4_NanoV12/hists_btag_ttbar_sf_MC_"$arg_campaign"_"$arg_year"_BTV_Run3_"$arg_year"_Comm_MINIAODv4_NanoV12.coffea"

mcfile_storage="$basepath/$arg_campaign/hists_MC_${arg_campaign}.coffea"
datafile_storage="$basepath/$arg_campaign/hists_data_${arg_campaign}.coffea"

# Main running scripts
echo "
Processing MC..."
mcjson="metadata/$arg_campaign/MC_"$arg_campaign"_"$arg_year"_BTV_Run3_"$arg_year"_Comm_MINIAODv4_NanoV12.json"
# mcjson="test_my_samples.json"
if [[ ! -f $mcfile_storage ]]; then
    python runner.py --json $mcjson --workflow btag_ttbar_sf --campaign $arg_campaign --year $arg_year --executor $arg_executor --retries 1 -s $scaleout --chunk $chunksize --limit $limit --max $max --overwrite
    mv $mcfile $mcfile_storage
elif [[ "$overwrite" == "true" ]]; then
    python runner.py --json $mcjson --workflow btag_ttbar_sf --campaign $arg_campaign --year $arg_year --executor $arg_executor --retries 1 -s $scaleout --chunk $chunksize --limit $limit --max $max --overwrite
    mv $mcfile $mcfile_storage
else
    echo "MC file already exists. Skipping..."
fi

# Check if the python command ran successfully
if [[ $? -ne 0 ]]; then
    echo "Error: MC Processing failed."
    return 1
fi

echo "
Processing data..."
if [[ ! -f $datafile_storage ]]; then
    python runner.py --json metadata/$arg_campaign/data_"$arg_campaign"_"$arg_year"_em_BTV_Run3_"$arg_year"_Comm_MINIAODv4_NanoV12.json --workflow btag_ttbar_sf --campaign $arg_campaign --year $arg_year --executor $arg_executor --retries 1 -s $scaleout --chunk $chunksize --limit $limit --max $max --overwrite
    mv $datafile $datafile_storage
elif [[ "$overwrite" == "true" ]]; then
    python runner.py --json metadata/$arg_campaign/data_"$arg_campaign"_"$arg_year"_em_BTV_Run3_"$arg_year"_Comm_MINIAODv4_NanoV12.json --workflow btag_ttbar_sf --campaign $arg_campaign --year $arg_year --executor $arg_executor --retries 1 -s $scaleout --chunk $chunksize --limit $limit --max $max --overwrite
    mv $datafile $datafile_storage
else
    echo "Data file already exists. Skipping..."
fi

if [[ $? -ne 0 ]]; then
    echo "Error: Data Processing failed."
    return 1
fi

# plots
echo "Plotting..."
python scripts/plotdataMC.py -i $mcfile_storage,$datafile_storage --lumi $arg_lumi -p btag_ttbar_sf --ext $arg_campaign -v *pt,*eta,*phi --log

# Check if the python command ran successfully
if [[ $? -ne 0 ]]; then
    echo "Error: Plotting script failed."
    return 1
fi
plotdir="$basepath/$arg_campaign/plots/"
mkdir --parents $plotdir
mv plot/btag_ttbar_sf_${arg_campaign}/* $plotdir
echo "Plots saved in $plotdir"