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
python runner.py --json metadata/MC_"$arg1"_"$arg2"_BTV_Run3_"$arg2"_Comm_MINIAODv4_NanoV12.json --workflow ctag_Wc_sf --campaign $arg1 --year $arg2 --executor $arg3 --retries 1 -s 5 --limit 1 --max 1
python runner.py --json metadata/data_"$arg1"_"$arg2"_mu_BTV_Run3_"$arg2"_Comm_MINIAODv4_NanoV12.json --workflow ctag_Wc_sf --campaign $arg1 --year $arg2 --executor $arg3 --retries 1 -s 5 --limit 1 --max 1 
python scripts/plotdataMC.py -i hists_ctag_Wc_sf_MC_"$arg1"_"$arg2"_BTV_Run3_"$arg2"_Comm_MINIAODv4_NanoV12/hists_ctag_Wc_sf_MC_"$arg1"_"$arg2"_BTV_Run3_"$arg2"_Comm_MINIAODv4_NanoV12.coffea,hists_ctag_Wc_sf_data_"$arg1"_"$arg2"_mu_BTV_Run3_"$arg2"_Comm_MINIAODv4_NanoV12/hists_ctag_Wc_sf_data_"$arg1"_"$arg2"_mu_BTV_Run3_"$arg2"_Comm_MINIAODv4_NanoV12.coffea --lumi $arg4 -p ctag_Wc_sf --ext $arg1 -v njet,mujet_eta,mujet_phi,mujet_pt,soft_l_pt,soft_l_eta,soft_l_phi,hl_pt,hl_eta,hl_phi,MET_pt,MET_phi,DeepCSV_trackSip3dSigAboveCharm,DeepCSV_trackSip3dValAboveCharm,DeepCSV_trackSip2dSigAboveCharm,DeepCSV_trackSip2dValAboveCharm,DeepJet_sv_d3dsig_0,DeepJet_sv_d3d_0,DeepJet_sv_ntracks_0,DeepJet_sv_mass_0,DeepJet_sv_normchi2_0,DeepJet_Cpfcan_BtagPf_trackSip2dSig_0,DeepJet_Cpfcan_BtagPf_trackSip2dVal_0,DeepJet_Cpfcan_BtagPf_trackSip3dSig_0,DeepJet_Cpfcan_BtagPf_trackSip3dVal_0,btagDeepFlavCvL_0,btagDeepFlavCvB_0,btagPNetCvL_0,btagPNetCvB_0,btagRobustParTAK4CvL_0,btagRobustParTAK4CvB_0,btagNegDeepFlavCvL_0,btagNegDeepFlavCvB_0,btagNegPNetCvL_0,btagNegPNetCvB_0,btagNegRobustParTAK4CvL_0,btagNegRobustParTAK4CvB_0
mv hists_ctag_Wc_sf_MC_"$arg1"_"$arg2"_BTV_Run3_"$arg2"_Comm_MINIAODv4_NanoV12/hists_ctag_Wc_sf_MC_"$arg1"_"$arg2"_BTV_Run3_"$arg2"_Comm_MINIAODv4_NanoV12.coffea plot/BTV/ctag_Wc_sf_"$arg1"_"$DATE"/hists_ctag_Wc_sf_MC_"$arg1"_"$arg2"_BTV_Run3_"$arg2"_Comm_MINIAODv4_NanoV12.coffea
mv hists_ctag_Wc_sf_data_"$arg1"_"$arg2"_mu_BTV_Run3_"$arg2"_Comm_MINIAODv4_NanoV12/hists_ctag_Wc_sf_data_"$arg1"_"$arg2"_mu_BTV_Run3_"$arg2"_Comm_MINIAODv4_NanoV12.coffea plot/BTV/ctag_Wc_sf_"$arg1"_"$DATE"/hists_ctag_Wc_sf_data_"$arg1"_"$arg2"_mu_BTV_Run3_"$arg2"_Comm_MINIAODv4_NanoV12.coffea
xrdcp plot/BTV/ctag_Wc_sf_"$arg1"_"$DATE"/* root://eosuser.cern.ch//eos/user/b/btvweb/www/Commissioning/dataMC/"$DATE"/"$arg1"/ctag_Wc_sf/
