#!/bin/bash
DATE=$(date +%y_%m_%d)
echo "running on $DATE"

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

chunksize=75000

# check files
datafile="hists_btag_ttbar_sf_data_"$arg_campaign"_"$arg_year"_em_BTV_Run3_"$arg_year"_Comm_MINIAODv4_NanoV12/hists_btag_ttbar_sf_data_"$arg_campaign"_"$arg_year"_em_BTV_Run3_"$arg_year"_Comm_MINIAODv4_NanoV12.coffea"
mcfile="hists_btag_ttbar_sf_MC_"$arg_campaign"_"$arg_year"_BTV_Run3_"$arg_year"_Comm_MINIAODv4_NanoV12/hists_btag_ttbar_sf_MC_"$arg_campaign"_"$arg_year"_BTV_Run3_"$arg_year"_Comm_MINIAODv4_NanoV12.coffea"

# Main running scripts
# echo "Processing MC..."
# if [[ ! -f $mcfile ]]; then
#     python runner.py --json metadata/MC_"$arg_campaign"_"$arg_year"_BTV_Run3_"$arg_year"_Comm_MINIAODv4_NanoV12.json --workflow btag_ttbar_sf --campaign $arg_campaign --year $arg_year --executor $arg_executor --retries 1 -s $scaleout --chunk $chunksize --limit 10  # --max 1
# elif [[ "$overwrite" == "true" ]]; then
#     python runner.py --json metadata/MC_"$arg_campaign"_"$arg_year"_BTV_Run3_"$arg_year"_Comm_MINIAODv4_NanoV12.json --workflow btag_ttbar_sf --campaign $arg_campaign --year $arg_year --executor $arg_executor --retries 1 -s $scaleout --chunk $chunksize --limit 10 --overwrite  # --max 1
# fi

echo "Processing data..."
if [[ ! -f $datafile ]]; then
    python runner.py --json metadata/data_"$arg_campaign"_"$arg_year"_em_BTV_Run3_"$arg_year"_Comm_MINIAODv4_NanoV12.json --workflow btag_ttbar_sf --campaign $arg_campaign --year $arg_year --executor $arg_executor --retries 1 -s $scaleout --chunk $chunksize --limit 10  # --max 1 
elif [[ "$overwrite" == "true" ]]; then
    python runner.py --json metadata/data_"$arg_campaign"_"$arg_year"_em_BTV_Run3_"$arg_year"_Comm_MINIAODv4_NanoV12.json --workflow btag_ttbar_sf --campaign $arg_campaign --year $arg_year --executor $arg_executor --retries 1 -s $scaleout --chunk $chunksize --limit 20 --overwrite  # --max 1
fi

# plots
echo "Plotting..."
python scripts/plotdataMC.py -i $mcfile,$datafile --lumi $arg_lumi -p btag_ttbar_sf --ext $arg_campaign -v all --log


# move files
echo "Moving files..."
basepath="/net/data_cms3a-1/fzinn/BTV/btag_sf/nobackup"
if [[ -f $mcfile ]]; then
    mv $mcfile $basepath/$arg_campaign/hists_MC_${arg_year}_BTV_em.coffea
fi
if [[ -f $datafile ]]; then
    mv $datafile $basepath/$arg_campaign/hists_data_${arg_year}_em.coffea
fi
plotdir="$basepath/$arg_campaign/plots/"
mkdir --parents $plotdir
mv plot/BTV/btag_ttbar_sf_${arg_campaign}_${DATE}/* $basepath/$arg_campaign/plots/