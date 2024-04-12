#!/bin/bash
#####Based on Ming-Yan's createdir.py script#######
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
campaign="${arg_array[0]}"

# Define workflows and campaigns
workflows=("ctag_Wc_sf" "ctag_DY_sf" "ttdilep_sf" "ttsemilep_sf")
campaigns=("Winter22" "Summer22" "Summer22EE" "Summer23" "Summer23BPix" "all")

# Check if campaign is provided and valid
if [[ -z "$campaign" ]]; then
    echo "Campaign name is required."
    exit 1
elif [[ ! " ${campaigns[@]} " =~ " ${campaign} " ]]; then
    echo "Invalid campaign name."
    exit 1
fi

# Create directories
for wf in "${workflows[@]}"; do
    if [[ ! -d "dataMC/$DATE/$campaign/$wf" ]]; then
        mkdir -p "dataMC/$DATE/$campaign/$wf"
        xrdcp root://eosuser.cern.ch//eos/user/b/btvweb/www/Commissioning/index.php "dataMC/$DATE/$campaign/$wf/"
    fi
done

# Sending the dir structure to BTV web
xrdcp -r -f dataMC/"$DATE"/ root://eosuser.cern.ch//eos/user/b/btvweb/www/Commissioning/dataMC/