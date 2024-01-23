#!/bin/bash
#Usage: ./run.sh ${executor} ${Lumi}
echo "Starting the BTV Commissioning task........"
echo "Fetching Dataset DAS paths from https://btvweb.web.cern.ch/BTVNanoProduction/csvoutputs/prod_tables.php"
python input_dataset.py
python scripts/fetch.py --input "$1_$2.txt" --output "$1_$2_Run3_BTV_Comm_v1.json"
echo "Running workflows........................................................."
echo "Running over Data"
echo "Running Dileptonic ttbar phase space : check performance for btag SFs, emu channel"
python runner.py --workflow ttdilep_sf --json "metadata/$1_$2_Run3_BTV_Comm_v1.json"  --campaign $3 --year $4 --executor $5 
echo "Running Semileptonic ttbar phase space : check performance for btag SFs, muon channel"
python runner.py --workflow ttsemilep_sf --json metadata/.json --campaign $3 --year $4 --executor $5
echo "Running Dileptonic ttbar phase space : check performance for charm SFs, bjets enriched SFs, muon channel"
python runner.py --workflow ctag_ttdilep_sf --json metadata/data_Summer22_Run3_2022_mu_BTV_Comm_v2_NanoV12_noPF.json  --campaign Summer22Run3 --year 2022(--executor ${scaleout_site})
echo "Running Semileptonic ttbar phase space : check performance for charm SFs, bjets enriched SFs, muon channel"
python runner.py --workflow ctag_ttsemilep_sf --json metadata/data_Summer22_Run3_2022_mu_BTV_Comm_v2_NanoV12_noPF.json  --campaign Summer22Run3 --year 2022(--executor ${scaleout_site})
echo "Running W+c phase space : check performance for charm SFs, cjets enriched SFs, muon channel"
python runner.py --workflow ctag_Wc_sf --json metadata/data_Summer22_Run3_2022_mu_BTV_Comm_v2_NanoV12_noPF.json  --campaign Summer22Run3 --year 2022(--executor ${scaleout_site})
echo "Running DY phase space : check performance for charm SFs, light jets enriched SFs, muon channel"
python runner.py --workflow ctag_DY_sf --json metadata/data_Summer22_Run3_2022_mu_BTV_Comm_v2_NanoV12_noPF.json  --campaign Summer22Run3 --year 2022(--executor ${scaleout_site})



