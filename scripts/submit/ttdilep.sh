#!/bin/bash
DATE=$(date +%y_%m_%d)
echo "running on $DATE"
# Check if the arguments file exists
if [ ! -f "arguments.txt" ]; then
    echo "Error: arguments.txt not found."
    exit 1
fi

# Read arguments from the file
args=$(<arguments.txt)
# Split the arguments into an array
IFS=$'\n' read -d '' -r -a arg_array <<< "$args"

# Check if enough arguments are provided
if [ "${#arg_array[@]}" -lt 4 ]; then
    echo "Error: Insufficient arguments provided in arguments.txt."
    exit 1
fi

# Assign arguments to variables
arg1="${arg_array[0]}"
arg2="${arg_array[1]}"
arg3="${arg_array[2]}"
arg4="${arg_array[3]}"

# Main running scripts
python runner.py --json metadata/MC_"$arg1"_"$arg2"_BTV_Run3_"$arg2"_Comm_MINIAODv4_NanoV12.json --workflow ttdilep_sf --campaign $arg1 --year $arg2 --executor $arg3 --retries 1 -s 5 --limit 1 --max 1
python runner.py --json metadata/data_"$arg1"_"$arg2"_em_BTV_Run3_"$arg2"_Comm_MINIAODv4_NanoV12.json --workflow ttdilep_sf --campaign $arg1 --year $arg2 --executor $arg3 --retries 1 -s 5 --limit 1 --max 1 
python scripts/plotdataMC.py -i hists_ttdilep_sf_MC_"$arg1"_"$arg2"_BTV_Run3_"$arg2"_Comm_MINIAODv4_NanoV12/hists_ttdilep_sf_MC_"$arg1"_"$arg2"_BTV_Run3_"$arg2"_Comm_MINIAODv4_NanoV12.coffea,hists_ttdilep_sf_data_"$arg1"_"$arg2"_em_BTV_Run3_"$arg2"_Comm_MINIAODv4_NanoV12/hists_ttdilep_sf_data_"$arg1"_"$arg2"_em_BTV_Run3_"$arg2"_Comm_MINIAODv4_NanoV12.coffea --lumi $arg4 -p ttdilep_sf --ext $arg1 -v all
mv hists_ttdilep_sf_MC_"$arg1"_"$arg2"_BTV_Run3_"$arg2"_Comm_MINIAODv4_NanoV12/hists_ttdilep_sf_MC_"$arg1"_"$arg2"_BTV_Run3_"$arg2"_Comm_MINIAODv4_NanoV12.coffea plot/BTV/ttdilep_sf_"$arg1"_"$DATE"/hists_ttdilep_sf_MC_"$arg1"_"$arg2"_BTV_Run3_"$arg2"_Comm_MINIAODv4_NanoV12.coffea
mv hists_ttdilep_sf_data_"$arg1"_"$arg2"_em_BTV_Run3_"$arg2"_Comm_MINIAODv4_NanoV12/hists_ttdilep_sf_data_"$arg1"_"$arg2"_em_BTV_Run3_"$arg2"_Comm_MINIAODv4_NanoV12.coffea plot/BTV/ttdilep_sf_"$arg1"_"$DATE"/hists_ttdilep_sf_data_"$arg1"_"$arg2"_em_BTV_Run3_"$arg2"_Comm_MINIAODv4_NanoV12.coffea
mv plot/BTV/ttdilep_sf_"$arg1"_"$DATE"/ BTV_Run3_"$arg2"_Comm_MinoAODv4_"$DATE"/"$arg1"/ 
xrdcp plot/BTV/ttdilep_sf_"$arg1"_"$DATE"/* root://eosuser.cern.ch//eos/user/b/btvweb/www/Commissioning/dataMC/"$DATE"/"$arg1"/ttdilep_sf/
