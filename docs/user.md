# For user : Running the workflows

![structure](figs/comm_wf.png)


All in one script : `scripts/suball.py`

All the steps are summarized in the [`suball.py`](#all-in-one-script--scriptssuballpy) scripts for the existing workflows. You can simply just run

```python
python scripts/suball.py --scheme ${SCHEME_FOR_STUDY} --campaign ${CAMPAIGN_FOR_SF} --year ${YEAR}  --DAS_campaign "$DATA_CAMPAIGN_RGX,$MC_CAMPAIGN_RGX"
#Example with 2023 Summer23 campaign
python scripts/suball.py --scheme default_comissioning --campaign Summer23 --year 2023  --DAS_campaign "*Run2023D*Sep2023*,*Run3Summer23BPixNanoAODv12-130X*" 
```
This wrap up the steps mentioned above as a streamline to obtained the required info

The only missing item need to do manually is to change the updated correction in `AK4_parameter.py` as written [here](#correction-files-configurations)
Each steps are also explained in detailed below, this can be obtain by 

## 0.  Make the dataset json files 

Use `fetch.py` in folder `scripts/` to obtain your samples json files for the predefined workflow with the refined MC. For more flexible usage please find [details](scripts.md#fetchpy--create-input-json)

The fetch script reads the predefine data & MC samples dataset name and output the json file to `metadata/$CAMPAIGN/`, but to find the exact dataset for BTV studies, we usually need to specify the `DAS_campaign`

```
python scripts/fetch.py -c {campaign} --year {args.year}  --from_workflow {wf} --DAS_campaign {DAS_campaign} {overwrite} 
# campaign :  the campaign name like Summer23,Winter22
# year : data taking years 2022/2023...
# wf: workflow name like ttdilep_sf, ctag_Wc_sf
# DAS_campaign: Input the campaign name for DAS to search appropriate campaigns, use in dataset construction , please do `campaign1,campaign2,campaign3`. Also supports "auto" (hard-coded!) if campaign and year are specified.
# overwrite (bool): recreate the exist json file 
```


:::{caution}
Do not make the file list greater than 4k files to avoid scaleout issues in various site (file open limit)
:::

:::{tip}
If `gfal-ls` does not work on your machine, reset the gfal-python with `GFAL_PYTHONBIN=/usr/bin/python3`
:::

## 1. Correction files configurations & add new correction files (Optional)

If the correction files are not supported yet by jsonpog-integration, you can still try with custom input data.

All the `lumiMask`, correction files (SFs, pileup weight), and JEC, JER files are under  `BTVNanoCommissioning/src/data/` following the substructure `${type}/${year}_${campaign}/${files}`(except `lumiMasks` and `Prescales`)

| Type        | File type |  Comments|
| :---:   | :---: | :---: |
| `lumiMasks` |`.json` | Masked good lumi-section used for physics analysis|
| `Prescales` | `.json.` | HLT paths for prescaled triggers|
| `PU`  | `.pkl.gz` or `.histo.root` | Pileup reweight files, matched MC to data| 
| `LSF` | `.histo.root` | Lepton ID/Iso/Reco/Trigger SFs|
| `BTV` | `.csv` or `.root` | b-tagger, c-tagger SFs|
| `JME` | `.txt` | JER, JEC files|
| `JPCalib` | `.root` | Jet probablity calibration, used in LTSV methods|

Create a `dict` entry under `correction_config` with dedicated campaigns in `BTVNanoCommissioning/src/utils/AK4_parameters.py`.


  
  
The official correction files collected in [jsonpog-integration](https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration) is updated by POG, except `lumiMask` and `JME` still updated by by the BTVNanoCommissioning framework user/developer.  For centrally maintained correction files, no input files have to be defined anymore in the `correction_config`. The example to implemented new corrections from POG can be found in [git](https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/blob/master/examples/), and the contents of the correction files are in the [summary](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/)


  ```python
    "2017_UL": {
          # Same with custom config
          "lumiMask": "Cert_294927-306462_13TeV_UL2017_Collisions17_MuonJSON.json",

          "JME": {
             "MC": "Summer19UL17_V5_MC",
            "Run2017F": "Summer19UL17_RunF_V5_DATA",
          },
          ### Alternatively, take the txt files in  https://github.com/cms-jet/JECDatabase/tree/master/textFiles
          "JME": {
                      # specified the name of JEC
                      "name": "V1_AK4PFPuppi",
                      # dictionary of jec text files
                      "MC": [
                          "Summer23Prompt23_V1_MC_L1FastJet_AK4PFPuppi",
                          "Summer23Prompt23_V1_MC_L2Relative_AK4PFPuppi",
                          "Summer23Prompt23_V1_MC_L2Residual_AK4PFPuppi",
                          "Summer23Prompt23_V1_MC_L3Absolute_AK4PFPuppi",
                          "Summer23Prompt23_V1_MC_UncertaintySources_AK4PFPuppi",
                          "Summer23Prompt23_V1_MC_Uncertainty_AK4PFPuppi",
                          "Summer23Prompt23_JRV1_MC_SF_AK4PFPuppi",
                          "Summer23Prompt23_JRV1_MC_PtResolution_AK4PFPuppi",
                      ],
                      "dataCv123": [
                          "Summer23Prompt23_RunCv123_V1_DATA_L1FastJet_AK4PFPuppi",
                          "Summer23Prompt23_RunCv123_V1_DATA_L2Relative_AK4PFPuppi",
                          "Summer23Prompt23_RunCv123_V1_DATA_L3Absolute_AK4PFPuppi",
                          "Summer23Prompt23_RunCv123_V1_DATA_L2L3Residual_AK4PFPuppi",
                      ],
                      "dataCv4": [
                          "Summer23Prompt23_RunCv4_V1_DATA_L1FastJet_AK4PFPuppi",
                          "Summer23Prompt23_RunCv4_V1_DATA_L2Relative_AK4PFPuppi",
                          "Summer23Prompt23_RunCv4_V1_DATA_L3Absolute_AK4PFPuppi",
                          "Summer23Prompt23_RunCv4_V1_DATA_L2L3Residual_AK4PFPuppi",
                      ],
                  },
          ###
          # no config need to be specify for PU weights
          "PU": None,
          # Alternatively, take root file as input
          "PU": "puwei_Summer23.histo.root",
          # Btag SFs - specify $TAGGER : $TYPE-> find [$TAGGER_$TYPE] in json file
          "BTV": {"deepCSV": "shape", "deepJet": "shape"},
          "roccor": None,
          # JMAR, IDs from JME- Following the scheme: "${SF_name}": "${WP}"
          "JMAR": {"PUJetID_eff": "L"},
          "LSF": {
          # Electron SF - Following the scheme: "${SF_name} ${year}": "${WP}"
          # https://github.com/cms-egamma/cms-egamma-docs/blob/master/docs/EgammaSFJSON.md
              "ele_ID 2017": "wp90iso",
              "ele_Reco 2017": "RecoAbove20",

          # Muon SF - Following the scheme: "${SF_name} ${year}": "${WP}"

              "mu_Reco 2017_UL": "NUM_TrackerMuons_DEN_genTracks",
              "mu_HLT 2017_UL": "NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight",
              "mu_ID 2017_UL": "NUM_TightID_DEN_TrackerMuons",
              "mu_Iso 2017_UL": "NUM_TightRelIso_DEN_TightIDandIPCut",
          },
          # use for BTA production, jet probablity
        "JPCalib": {
            "Run2022E": "calibeHistoWrite_Data2022F_NANO130X_v1.root",
            "Run2022F": "calibeHistoWrite_Data2022F_NANO130X_v1.root",
            "Run2022G": "calibeHistoWrite_Data2022G_NANO130X_v1.root",
            "MC": "calibeHistoWrite_MC2022EE_NANO130X_v1.root",
        },
      },
  ```




## 2. Run the workflow to get coffea files

The `runner.py` handles the options to select the workflow with dedicated configuration for each campaign. The miniumum required info is 
```
python runner.py --wf {wf} --json metadata/{args.campaign}/{types}_{args.campaign}_{args.year}_{wf}.json {overwrite} --campaign {args.campaign} --year {args.year} 
```

:::{tip}
- In case just to test your program, you can limit only one file with one chunk using iterative executor to avoid overwriting error message by `--max 1 --limit 1 --executor iterative`
- In case you only want to run particular sample in your json `--only $dataset_name`, i.e. `--only TT*` or `--only MuonEG_Run2023A`
- Change the numbers of scale job by `-s $NJOB`
- Store the arrays by setting the flag `--isArray`
- Modifying chunk size in case the jobs is to big `--chunk $N_EVENTS_PER_CHUNK`
- Sometimes the global redirector is insufficient, you can increase the numbers of retries (only in parsl/dask) `--retries 30`, or skip the files `--skipbadfiles` and later reprocess the missing info by create the json with skipped files. Methods to create the json files discussed in the next part.
:::

Other options detail can be found here 

<details><summary>runner options</summary>
<p>

```python
### ====> REQUIRED <=======
# --wf {validation,ttdilep_sf,ttsemilep_sf,c_ttsemilep_sf,emctag_ttdilep_sf,ctag_ttdilep_sf,ectag_ttdilep_sf,ctag_ttsemilep_sf,ectag_ttsemilep_sf,QCD_sf,QCD_smu_sf,ctag_Wc_sf,ectag_Wc_sf,ctag_DY_sf,ectag_DY_sf,BTA,BTA_addPFMuons,BTA_addAllTracks,BTA_ttbar}, --workflow {validation,ttdilep_sf,ttsemilep_sf,c_ttsemilep_sf,emctag_ttdilep_sf,ctag_ttdilep_sf,ectag_ttdilep_sf,ctag_ttsemilep_sf,ectag_ttsemilep_sf,QCD_sf,QCD_smu_sf,ctag_Wc_sf,ectag_Wc_sf,ctag_DY_sf,ectag_DY_sf,BTA,BTA_addPFMuons,BTA_addAllTracks,BTA_ttbar}
#                         Which processor to run
#   -o OUTPUT, --output OUTPUT
#                         Output histogram filename (default: hists.coffea)
#   --json SAMPLEJSON     JSON file containing dataset and file locations (default: dummy_samples.json)
#   --year YEAR           Year
#   --campaign CAMPAIGN   Dataset campaign, change the corresponding correction files

#=======Optional======
# ==> configurations for storing info
#   --isSyst {False,all,weight_only,JERC_split,JP_MC}
#                         Run with systematics, all, weights_only(no JERC uncertainties included),JERC_split, None
#   --isArray             Output root files
#   --noHist              Not output coffea histogram
#   --overwrite           Overwrite existing files
# ==> scale out options
#   --executor {iterative,futures,parsl/slurm,parsl/condor,parsl/condor/naf_lite,dask/condor,dask/condor/brux,dask/slurm,dask/lpc,dask/lxplus,dask/casa,condor_standalone}
#                         The type of executor to use (default: futures). Other options can be implemented. For example see https://parsl.readthedocs.io/en/stable/userguide/configuring.html-
#                         `parsl/slurm` - tested at DESY/Maxwell- `parsl/condor` - tested at DESY, RWTH- `parsl/condor/naf_lite` - tested at DESY- `dask/condor/brux` - tested at BRUX (Brown U)-
#                         `dask/slurm` - tested at DESY/Maxwell- `dask/condor` - tested at DESY, RWTH- `dask/lpc` - custom lpc/condor setup (due to write access restrictions)- `dask/lxplus` - custom
#                         lxplus/condor setup (due to port restrictions)
#   -j WORKERS, --workers WORKERS
#                         Number of workers (cores/threads) to use for multi-worker executors (e.g. futures or condor) (default: 3)
#   -s SCALEOUT, --scaleout SCALEOUT
#                         Number of nodes to scale out to if using slurm/condor. Total number of concurrent threads is ``workers x scaleout`` (default: 6)
#   --memory MEMORY       Memory used in jobs default ``(default: 4.0)
#   --disk DISK           Disk used in jobs default ``(default: 4)
#   --voms VOMS           Path to voms proxy, made accessible to worker nodes. By default a copy will be made to $HOME.
#   --chunk N             Number of events per process chunk
#   --retries N           Number of retries for coffea processor
#   --fsize FSIZE         (Specific for dask/lxplus file splitting, default: 50) Numbers of files processed per dask-worker
#   --index INDEX         (Specific for dask/lxplus file splitting, default: 0,0) Format: $dict_index_start,$file_index_start,$dict_index_stop,$file_index_stop. Stop indices are optional. $dict_index
#                         refers to the index, splitted $dict_index and $file_index with ','$dict_index refers to the sample dictionary of the samples json file. $file_index refers to the N-th batch
#                         of files per dask-worker, with its size being defined by the option --index. The job will start (stop) submission from (with) the corresponding indices.
# ==> debug option
#   --validate            Do not process, just check all files are accessible
#   --skipbadfiles        Skip bad files.
#   --only ONLY           Only process specific dataset or file
#   --limit N             Limit to the first N files of each dataset in sample JSON
#   --max N               Max number of chunks to run in total
```

</p>
</details>


## 3. Dump processed information to obtain luminoisty and processed files

After obtained `coffea` file, we can check the processed files and obtain the luminoisty.

Get the run & luminosity information for the processed events from the coffea output files. When you use `--skipbadfiles`, the submission will ignore files not accesible(or time out) by xrootd. This script helps you to dump the processed luminosity into a json file which can be calculated by brilcalc tool and provide a list of failed lumi sections by comparing the original json input to the one from the `.coffea` files.
We will see the luminosity info in `/pb` and the skipped files as new json for resubmission.


```bash
python scripts/dump_processed.py -t all -c INPUT_COFFEA --json ORIGINAL_JSON_INPUT -n {args.campaign}_{args.year}_{wf}
#   -t {all,lumi,failed}, --type {all,lumi,failed}
#                         Choose the function for dump luminosity(`lumi`)/failed files(`failed`) into json
#   -c COFFEA, --coffea COFFEA
#                         Processed coffea files, splitted by ,. Wildcard option * available as well.
#   -n FNAME, --fname FNAME
#                         Output name of jsons(with _lumi/_dataset)
#   -j JSONS, --jsons JSONS
#                         Original input json files, splitted by ,. Wildcard option * available as well.
```

## 4. Obtain data/MC plots 

We can obtain data/MC plots from coffea via the plotting scripts: 

You can specify `-v all` to plot all the variables in the `coffea` file, or use wildcard options (e.g. `-v "*DeepJet*"` for the input variables containing `DeepJet`)

:new: non-uniform rebinning is possible, specify the bins with  list of edges `--autorebin 50,80,81,82,83,100.5`

```bash
python scripts/plotdataMC.py -i a.coffea,b.coffea --lumi 41500 -p ttdilep_sf -v z_mass,z_pt  
python scripts/plotdataMC.py -i "test*.coffea" --lumi 41500 -p ttdilep_sf -v z_mass,z_pt # with wildcard option need ""

```

<details><summary>options</summary>
<p>
```
options:
  -h, --help            show this help message and exit
  --lumi LUMI           luminosity in /pb
  --com COM             sqrt(s) in TeV
  -p {ttdilep_sf,ttsemilep_sf,ctag_Wc_sf,ctag_DY_sf,ctag_ttsemilep_sf,ctag_ttdilep_sf}, --phase {dilep_sf,ttsemilep_sf,ctag_Wc_sf,ctag_DY_sf,ctag_ttsemilep_sf,ctag_ttdilep_sf}
                        which phase space
  --log LOG             log on y axis
  --norm NORM           Use for reshape SF, scale to same yield as no SFs case
  -v VARIABLE, --variable VARIABLE
                        variables to plot, splitted by ,. Wildcard option * available as well. Specifying `all` will run through all variables.
  --SF                  make w/, w/o SF comparisons
  --ext EXT             prefix name
  -i INPUT, --input INPUT
                        input coffea files (str), splitted different files with ','. Wildcard option * available as well.
   --autorebin AUTOREBIN
                        Rebin the plotting variables, input `int` or `list`. int: merge N bins. list of number: rebin edges(non-uniform bin is possible)
   --xlabel XLABEL      rename the label for x-axis
   --ylabel YLABEL      rename the label for y-axis
   --splitOSSS SPLITOSSS 
                        Only for W+c phase space, split opposite sign(1) and same sign events(-1), if not specified, the combined OS-SS phase space is used
   --xrange XRANGE      custom x-range, --xrange xmin,xmax
   --flow FLOW 
                        str, optional {None, 'show', 'sum'} Whether plot the under/overflow bin. If 'show', add additional under/overflow bin. If 'sum', add the under/overflow bin content to first/last bin.
   --split {flavor,sample,sample_flav}
                        Decomposition of MC samples. Default is split to jet flavor(udsg, pu, c, b), possible to split by group of MC
                        samples. Combination of jetflavor+ sample split is also possible 
```

</p>
</details>



## Reading coffea `hist`  


Quick tutorial to go through coffea files. Example coffea files can be found in `testfile/` 


### Structure of the file

The coffea contains histograms  wrapped in a dictionary with `$dataset:{$histname:hist}`, the `hist` is the object using
[hist](https://hist.readthedocs.io/en/latest/) which allows multidimensional bins with different types of array
```python
{'WW_TuneCP5_13p6TeV-pythia8':{
'btagDeepFlavB_b_0': Hist(
  IntCategory([0, 1, 4, 5, 6], name='flav', label='Genflavour'),
  IntCategory([1, -1], name='osss', label='OS(+)/SS(-)'),
  StrCategory(['noSF'], growth=True, name='syst'),
  Regular(30, -0.2, 1, name='discr', label='btagDeepFlavB_b'),
  storage=Weight()) # Sum: WeightedSum(value=140, variance=140), 'btagDeepFlavB_bb_0': Hist(
  IntCategory([0, 1, 4, 5, 6], name='flav', label='Genflavour'),
  IntCategory([1, -1], name='osss', label='OS(+)/SS(-)'),
  StrCategory(['noSF'], growth=True, name='syst'),
  Regular(30, -0.2, 1, name='discr', label='btagDeepFlavB_bb'),
  storage=Weight()) # Sum: WeightedSum(value=140, variance=140), 'btagDeepFlavB_lepb_0': Hist(
  IntCategory([0, 1, 4, 5, 6], name='flav', label='Genflavour'),
  IntCategory([1, -1], name='osss', label='OS(+)/SS(-)'),
  StrCategory(['noSF'], growth=True, name='syst'),
  Regular(30, -0.2, 1, name='discr', label='btagDeepFlavB_lepb'),
  storage=Weight()) # Sum: WeightedSum(value=140, variance=140)}}
```
There are also `column_array` stores the processed file and lumi/run info in each dataset for data. The information are used in [dump_processed info](user.md#3-dump-processed-information-to-obtain-luminoisty-and-processed-files) 



The histogram is a multidimentinal histogram, with all the axis listed
```python
Hist(
  IntCategory([0, 1, 4, 5, 6], name='flav', label='Genflavour'),# different genflavor, 0 for light, 1 for PU, 2 for c, 3 for b. Always 0 for data.
  IntCategory([1, -1], name='osss', label='OS(+)/SS(-)'),# opposite sign or same sign, only appears in W+c workflow
  StrCategory(['noSF','PUUp','PUDown'], growth=True, name='syst'),# systematics variations,
  Regular(30, -0.2, 1, name='discr', label='btagDeepFlavB_lepb'),# discriminator distribution, the last axis is always the variable
  storage=Weight()) # Sum: WeightedSum(value=140, variance=140)# Value is sum of the entries, Variances is sum of the variances.
```

### Read coffea files and explore the histogram

```python
from coffea.util import load
# open coffea file
output=load("hists_ctag_Wc_sf_VV.coffea"）
# get the histogram and read the info
hist=output['WW_TuneCP5_13p6TeV-pythia8']['btagDeepFlavB_lepb_0']
# addition for two histogram is possible if the axis is the same
histvv=output['WW_TuneCP5_13p6TeV-pythia8']['btagDeepFlavB_lepb_0']+
       output['WZ_TuneCP5_13p6TeV-pythia8']['btagDeepFlavB_lepb_0']+
       output['ZZ_TuneCP5_13p6TeV-pythia8']['btagDeepFlavB_lepb_0']
# To get 1D histogram, we need to reduce the dimention
# we can specify the axis we want to read, e.g. read charm jet, opposite sign events with noSF
axis={'flav':3,'os':0,'syst':'noSF'}
hist1d=hist[axis] #--> this is the 1D histogram Hist
# you can also sum over the axis, e.g. here shows no jet flavor split and sum os+ss
axis={'flav':sum,'os':sum,'syst':'noSF'}
# rebin the axis is also possible, rebin the discrimnator by merged two bins into one
axis={'flav':sum,'os':sum,'syst':'noSF','discr':hist.rebin(2)}
```

### Plot the histogram
You can simply plot the histogram using [mplhep](https://mplhep.readthedocs.io/en/latest/)
```python
import mplhep as hep
import matplotlib.pyplot as plt
# set the plot style like tdr style
plt.style.use(hep.style.ROOT)
# make 1D histogram plot
hep.histplot(hist1D) 
```
### convert coffea hist to ROOT TH1

 `scripts/make_template.py` does the work to convert the coffea hist into 1D/2D ROOT histogram:

```python
python scripts/make_template.py -i $INPUT_COFFEA --lumi $LUMI_IN_invPB -o $ROOTFILE_NAME -v $VARIABLE -a $HIST_AXIS 
## Example
python scripts/make_template.py -i "testfile/*.coffea" --lumi 7650 -o test.root -v mujet_pt -a '{"flav":0,"osss":"sum"}' --mergemap 
#####
# -h, --help            show this help message and exit
#  -i INPUT, --input INPUT
#                        Input coffea file(s), can be a regular expression contains '*'
#  -v VARIABLE, --variable VARIABLE
#                        Variables to store(histogram name)
#  -a AXIS, --axis AXIS  dict, put the slicing of histogram, specify 'sum' option as string
#  --lumi LUMI           Luminosity in /pb , normalize the MC yields to corresponding luminosity
#  -o OUTPUT, --output OUTPUT
#                        output root file name
#  --mergemap MERGEMAP   Specify mergemap as dict, '{merge1:[dataset1,dataset2]...}' Also works with the json file with dict

#### EXAMPLE MERGEMAP
{
    "WJets": ["WJetsToLNu_TuneCP5_13p6TeV-madgraphMLM-pythia8"],
    "VV": [ "WW_TuneCP5_13p6TeV-pythia8", "WZ_TuneCP5_13p6TeV-pythia8", "ZZ_TuneCP5_13p6TeV-pythia8"],
    "TT": [ "TTTo2J1L1Nu_CP5_13p6TeV_powheg-pythia8", "TTTo2L2Nu_CP5_13p6TeV_powheg-pythia8"],
    "ST":[ "TBbarQ_t-channel_4FS_CP5_13p6TeV_powheg-madspin-pythia8", "TbarWplus_DR_AtLeastOneLepton_CP5_13p6TeV_powheg-pythia8", "TbarBQ_t-channel_4FS_CP5_13p6TeV_powheg-madspin-pythia8", "TWminus_DR_AtLeastOneLepton_CP5_13p6TeV_powheg-pythia8"],
"data":[ "Muon_Run2022C-PromptReco-v1", "SingleMuon_Run2022C-PromptReco-v1", "Muon_Run2022D-PromptReco-v1", "Muon_Run2022D-PromptReco-v2"]
}
```
 





