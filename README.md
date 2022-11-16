
# BTVNanoCommissioning
[![Linting](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/python_linting.yml/badge.svg)](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/python_linting.yml)
[![TTbar](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/ttbar_workflow.yml/badge.svg)](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/ttbar_workflow.yml)
[![TTbar DL+SL](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/ttbar_SL_DL_workflow.yml/badge.svg)](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/ttbar_SL_DL_workflow.yml)
[![ctag DY+jets Workflow](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/ctag_DY_workflow.yml/badge.svg)](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/ctag_DY_workflow.yml)
[![ctag W+c Workflow](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/ctag_Wc_workflow.yml/badge.svg)](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/ctag_Wc_workflow.yml)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

Repository for Commissioning studies in the BTV POG based on (custom) nanoAOD samples

## Requirements
### Setup 

:heavy_exclamation_mark: suggested to install under `bash` environment

```
# only first time 
git clone git@github.com:cms-btv-pog/BTVNanoCommissioning.git 

# activate enviroment once you have coffea framework 
conda activate btv_coffea
```
### Coffea installation with Miniconda
For installing Miniconda, see also https://hackmd.io/GkiNxag0TUmHnnCiqdND1Q#Local-or-remote
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
# Run and follow instructions on screen
bash Miniconda3-latest-Linux-x86_64.sh
```
NOTE: always make sure that conda, python, and pip point to local Miniconda installation (`which conda` etc.).

You can either use the default environment `base` or create a new one:
```
# create new environment with python 3.7, e.g. environment of name `btv_coffea`
conda create --name btv_coffea python=3.7
# activate environment `btv_coffea`
conda activate btv_coffea
```
You could simply create the environment through the existing `test_env.yml` under your conda environment
```
conda env create -f test_env.yml -p ${conda_dir}/envs/btv_coffea
```

Or install manually for the required packages, coffea, xrootd, and more:
```
pip install git+https://github.com/CoffeaTeam/coffea.git #latest published release with `pip install coffea`
conda install -c conda-forge xrootd
conda install -c conda-forge ca-certificates
conda install -c conda-forge ca-policy-lcg
conda install -c conda-forge dask-jobqueue
conda install -c anaconda bokeh 
conda install -c conda-forge 'fsspec>=0.3.3'
conda install dask
conda install -c anaconda 'openssl==1.1.1s'
conda install -c conda-forge parsl
```

Once the environment is set up, compile the python package:
```
pip install -e .
```

### Other installation options for coffea
See https://coffeateam.github.io/coffea/installation.html

## Structure
Example worfkflow for ttbar is included. 

Each workflow can be a separate "processor" file, creating the mapping from NanoAOD to
the histograms we need. Workflow processors can be passed to the `runner.py` script 
along with the fileset these should run over. Multiple executors can be chosen 
(for now iterative - one by one, uproot/futures - multiprocessing and dask-slurm). 

To test a small set of files to see whether the workflows run smoothly, run:
```
python runner.py --workflow ttsemilep_sf --json metadata/test_w_dj_mu.json --campaign Rereco17_94X --year 2017
```

More options for `runner.py` 
<details><summary>more options</summary>
<p>

```
--wf {validation,ttcom,ttdilep_sf,ttsemilep_sf,emctag_ttdilep_sf,ctag_ttdilep_sf,ectag_ttdilep_sf,ctag_ttsemilep_sf,ectag_ttsemilep_sf,ctag_Wc_sf,ectag_Wc_sf,ctag_DY_sf,ectag_DY_sf}, --workflow {validation,ttcom,ttdilep_sf,ttsemilep_sf,emctag_ttdilep_sf,ctag_ttdilep_sf,ectag_ttdilep_sf,ctag_ttsemilep_sf,ectag_ttsemilep_sf,ctag_Wc_sf,ectag_Wc_sf,ctag_DY_sf,ectag_DY_sf}
                        Which processor to run
  -o OUTPUT, --output OUTPUT
                        Output histogram filename (default: hists.coffea)
  --samples SAMPLEJSON, --json SAMPLEJSON
                        JSON file containing dataset and file locations
                        (default: dummy_samples.json)
  --year YEAR           Year
  --campaign CAMPAIGN   Dataset campaign, change the corresponding correction
                        files
  --isCorr              Run with SFs
  --isJERC              JER/JEC implemented to jet
  --executor {iterative,futures,parsl/slurm,parsl/condor,parsl/condor/naf_lite,dask/condor,dask/slurm,dask/lpc,dask/lxplus,dask/casa}
                        The type of executor to use (default: futures). 
  -j WORKERS, --workers WORKERS
                        Number of workers (cores/threads) to use for multi- worker executors (e.g. futures or condor) (default:
                        12)
  -s SCALEOUT, --scaleout SCALEOUT
                        Number of nodes to scale out to if using slurm/condor.
                        Total number of concurrent threads is ``workers x
                        scaleout`` (default: 6)
  --memory MEMORY       Memory used in jobs  ``(default: 4GB)
  --disk DISK           Disk used in jobs  ``(default: 4GB)
  --voms VOMS           Path to voms proxy, made accessible to worker nodes.
                        By default a copy will be made to $HOME.
  --chunk N             Number of events per process chunk
  --retries N           Number of retries for coffea processor
 --index INDEX         (Specific for dask/lxplus file splitting, default:0,0) 
                        Format: $dictindex,$fileindex. $dictindex refers to the index of the file list split to 50 files per dask-worker.
                        The job will start submission from the corresponding indices
  --validate            Do not process, just check all files are accessible
  --skipbadfiles        Skip bad files.
  --only ONLY           Only process specific dataset or file
  --limit N             Limit to the first N files of each dataset in sample
                        JSON
  --max N               Max number of chunks to run in total
```
</p>
</details>

### Roadmap for running the tool

1. Is the `.json` file ready? If not, create it following the instructions in the  [Make the json files](#make-the-json-files) section. Please use the correct naming scheme

2. Put the `lumiMask`, correction files (SFs, pileup weight), and JER, JEC files under the dict entry in `BTVNanoCommissioning/src/utils/AK4_parameters.py`. See details in [Correction files configurations](#correction-files-configurations)

3. If the JERC file `jec_compiled.pkl.gz` is missing in the `data/JME/${campaign}` directory, create it through [Create compiled JERC file](#create-compiled-jerc-filepklgz)

4. Run the workflow with dedicated input and campaign name. Example commands for Run 3 can be found [here](#commands-for-different-phase-space). For first usage, the JERC file needs to be recompiled first, see [Create compiled JERC file](#create-compiled-jerc-filepklgz)

5. Once you obtain the `.coffea` file(s), you can make plots using the [plotting scripts](#plotting-code), if the xsection for your sample is missing, please add to `src/BTVNanoCommissioning/helpers/xsection.py`

### Commands for different phase space

After a small test, you can run the full campaign for a dedicated phase space, separately for data and for MC.

#### b-SFs 

<details><summary>details</summary>
<p>

- Dileptonic ttbar phase space : check performance for btag SFs, emu channel

```
 python runner.py --workflow ttdilep_sf --json metadata/data_Winter22_emu_BTV_Run3_2022_Comm_v1.json  --campaign Winter22Run3 --year 2022 --isJERC --isCorr  (--executor ${scaleout_site}) 
```

- Semileptonic ttbar phase space : check performance for btag SFs, muon channel

```
python runner.py --workflow ttsemilep_sf --json metadata/data_Winter22_mu_BTV_Run3_2022_Comm_v1.json --campaign Winter22Run3 --year 2022 --isJERC --isCorr  (--executor ${scaleout_site})
```

</p>
</details>

#### c-SFs
<details><summary>details</summary>
<p>

- Dileptonic ttbar phase space : check performance for charm SFs, bjets enriched SFs, muon channel

```
python runner.py --workflow ctag_ttdilep_sf --json metadata/data_Winter22_mumu_BTV_Run3_2022_Comm_v1.json --campaign Winter22Run3 --year 2022 --isJERC --isCorr (--executor ${scaleout_site})
```


- Semileptonic ttbar phase space : check performance for charm SFs, bjets enriched SFs, muon channel

```
python runner.py --workflow ctag_ttsemilep_sf --json metadata/data_Winter22_mu_BTV_Run3_2022_Comm_v1.json --campaign Winter22Run3 --year 2022 --isJERC --isCorr (--executor ${scaleout_site})
```

- W+c phase space : check performance for charm SFs, cjets enriched SFs, muon  channel

```
python runner.py --workflow ctag_Wc_sf --json metadata/data_Winter22_mu_BTV_Run3_2022_Comm_v1.json --campaign Winter22Run3 --year 2022 --isJERC --isCorr (--executor ${scaleout_site})
```

- DY phase space : check performance for charm SFs, light jets enriched SFs, muon channel

```
python runner.py --workflow ctag_DY_sf --json metadata/data_Winter22_mumu_BTV_Run3_2022_Comm_v1.json --campaign Winter22Run3 --year 2022 --isJERC --isCorr (--executor ${scaleout_site})
```

</p>
</details>

#### Validation - check different campaigns

<details><summary>details</summary>
<p>

Only basic jet selections(PUID, ID, pT, $\eta$) applied. Put the json files with different campaigns

```
python runner.py --workflow valid --json metadata/$json file
```

</p>
</details>



## Scale-out (Sites)

Scale out can be notoriously tricky between different sites. Coffea's integration of `slurm` and `dask`
makes this quite a bit easier and for some sites the ``native'' implementation is sufficient, e.g Condor@DESY.
However, some sites have certain restrictions for various reasons, in particular Condor @CERN and @FNAL.

### Condor@FNAL (CMSLPC)
Follow setup instructions at https://github.com/CoffeaTeam/lpcjobqueue. After starting 
the singularity container run with 
```bash
python runner.py --wf ttcom --executor dask/lpc
```

### Condor@CERN (lxplus)
Only one port is available per node, so its possible one has to try different nodes until hitting
one with `8786` being open. Other than that, no additional configurations should be necessary.

```bash
python runner.py --wf ttcom --executor dask/lxplus
```

jobs automatically split to 50 files per jobs to avoid job failure due to crowded cluster on lxplus with the naming scheme `hist_$workflow_$json_$dictindex_$fileindex.coffea`. The `.coffea` files can be then combined at plotting level


:exclamation: The optimal scaleout options on lxplus are `-s 50 --chunk 50000`

To deal with unstable condor cluster and dask worker on lxplus, you can resubmit failure jobs via `--index $dictindex,$fileindex` option. `$dictindex` refers to the index in the `.json dict`, `$fileindex` refers to the index of the file list split to 50 files per dask-worker. The total number of files of each dict can be computed by `math.ceil(len($filelist)/50)` The job will start from the corresponding indices.

### Coffea-casa (Nebraska AF)
Coffea-casa is a JupyterHub based analysis-facility hosted at Nebraska. For more information and setup instuctions see
https://coffea-casa.readthedocs.io/en/latest/cc_user.html

After setting up and checking out this repository (either via the online terminal or git widget utility run with
```bash
python runner.py --wf ttcom --executor dask/casa
```
Authentication is handled automatically via login auth token instead of a proxy. File paths need to replace xrootd redirector with "xcache", `runner.py` does this automatically.


### Condor@DESY 
```bash
python runner.py --wf ttcom --executor dask/condor
```

### Maxwell@DESY 
```bash
python runner.py --wf ttcom --executor parsl/slurm
```


## Make the json files

Use the `fetch.py` in `filefetcher`, the `$input_DAS_list` is the info extract from DAS, and output json files in `metadata/`

```
python fetch.py --input ${input_DAS_list} --output ${output_json_name} --site ${site}
```
The `output_json_name` must contain the BTV name tag (e.g. `BTV_Run3_2022_Comm_v1`).

You might need to rename the json key name with following name scheme:

For the data sample please use the naming scheme,
```
$dataset_$Run
#i.e.
SingleMuon_Run2022C-PromptReco-v1
```
For MC, please be consistent with the dataset name in CMS DAS, as it cannot be mapped to the cross section otherwise.
```
$dataset
#i.e.
WW_TuneCP5_13p6TeV-pythia8
```

:exclamation: Do not make the file list greater than 4k files to avoid scaleout issues in various site

## Correction files configurations

All the `lumiMask`, correction files (SFs, pileup weight), and JEC, JER files are under  `BTVNanoCommissioning/src/data/` following the substructure `${type}/${campaign}/${files}`(except `lumiMasks` and `Prescales`)

| Type        | File type |  Comments|
| :---:   | :---: | :---: | 
| `lumiMasks` |`.json` | Masked good lumi-section used for physics analysis|
| `Prescales` | `.txt` | HLT paths for prescaled triggers|
| `PU`  | `.pkl.gz` or `.histo.root` | Pileup reweight files, matched MC to data| 
| `LSF` | `.histo.root` | Lepton ID/Iso/Reco/Trigger SFs|
| `BTV` | `.csv` or `.root` | b-tagger, c-tagger SFs|
| `JME` | `.txt` | JER, JEC files|

Create a `dict` entry under `correction_config` with dedicated campaigns in `BTVNanoCommissioning/src/utils/AK4_parameters.py`.

Take `Rereco17_94X` as an example.

```python
# specify campaign
"Rereco17_94X": 
{
      ##Load files with dedicated type:filename
        "lumiMask": "Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt",
        "PU": "94XPUwei_corrections.pkl.gz",
        "JME": "jec_compiled.pkl.gz",
      ## Btag SFs- create dict specifying SFs for DeepCSV b-tag(DeepCSVB),  DeepJet b-tag(DeepJetB),DeepCSV c-tag(DeepCSVC),  DeepJet c-tag(DeepJetC),
        "BTV": {
            ### b-tag 
            "DeepCSVB": "DeepCSV_94XSF_V5_B_F.csv",
            "DeepJetB": "DeepFlavour_94XSF_V4_B_F.csv",
            ### c-tag
            "DeepCSVC": "DeepCSV_ctagSF_MiniAOD94X_2017_pTincl_v3_2_interp.root",
            "DeepJetC": "DeepJet_ctagSF_MiniAOD94X_2017_pTincl_v3_2_interp.root",
        },
        ## lepton SF - create dict specifying SFs for electron/muon ID/ISO/RECO SFs
        "LSF": {
            ### Following the scheme "${SF_name} ${histo_name_in_root_file}": "${file}"
            "ele_Trig TrigSF": "Ele32_L1DoubleEG_TrigSF_vhcc.histo.root",
            "ele_ID EGamma_SF2D": "ElectronIDSF_94X_MVA80WP.histo.root",
            "ele_Rereco EGamma_SF2D": "ElectronRecoSF_94X.histo.root",
            "mu_ID NUM_TightID_DEN_genTracks_pt_abseta": "RunBCDEF_MuIDSF.histo.root",
            "mu_ID_low NUM_TightID_DEN_genTracks_pt_abseta": "RunBCDEF_MuIDSF_lowpT.histo.root",
            "mu_Iso NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta": "RunBCDEF_MuISOSF.histo.root",
        },
    },
```

## Create compiled JERC file(`pkl.gz`)

:exclamation: In case existing correction file doesn't work for you due to the incompatibility of `cloudpickle` in different python versions. Please recompile the file to get new pickle file.

Under `compile_jec.py` you need to create dedicated jet factory files with different campaigns. Following the name scheme with `mc` for MC and `data${run}` for data.

Compile correction pickle files for a specific JEC campaign by changing the dict of jet_factory, and define the MC campaign and the output file name by passing it as arguments to the python script:

```
python -m BTVNanoCommissioning.utils.compile_jec ${campaign} jec_compiled
e.g. python -m BTVNanoCommissioning.utils.compile_jec Winter22Run3 jec_compiled
```



## Plotting code

- data/MC comparisons
:exclamation_mark: If using wildcard for input, do not forget the quoatation marks! (see 2nd example below)

You can specify `-v all` to plot all the variables in the `coffea` file, or use wildcard options (e.g. `-v "*DeepJet*"` for the input variables containing `DeepJet`)

```
python plotdataMC.py -i a.coffea,b.coffea --lumi 41500 -p dilep_sf -v z_mass,z_pt
python plotdataMC.py -i "test*.coffea" --lumi 41500 -p dilep_sf -v z_mass,z_pt

options:
  -h, --help            show this help message and exit
  --lumi LUMI           luminosity in /pb
  --com COM             sqrt(s) in TeV
  -p {dilep_sf,ttsemilep_sf,ctag_Wc_sf,ctag_DY_sf,ctag_ttsemilep_sf,ctag_ttdilep_sf}, --phase {dilep_sf,ttsemilep_sf,ctag_Wc_sf,ctag_DY_sf,ctag_ttsemilep_sf,ctag_ttdilep_sf}
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
                        Rebin the plotting variables by merging N bins in case the current binning is too fine for you 
   --xlabel XLABEL      rename the label for x-axis
```
- data/data, MC/MC comparisons

You can specify `-v all` to plot all the variables in the `coffea` file, or use wildcard options (e.g. `-v "*DeepJet*"` for the input variables containing `DeepJet`)
:exclamation_mark: If using wildcard for input, do not forget the quoatation marks! (see 2nd example below)

```
python comparison.py -i a.coffea,b.coffea -p ttsemilep_sf -r SingleMuon_Run2017B-106X_PFNanov1 -c DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8 -v DeepJet_Cpfcan_BtagPf_trackJetDistVal_0 --shortref Run2017B --shortcomp DYJets (--sepflav True/False)
python comparison.py -i "test*.coffea" -p ttsemilep_sf -r SingleMuon_Run2017B-106X_PFNanov1 -c DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8 -v DeepJet_Cpfcan_BtagPf_trackJetDistVal_0 --shortref Run2017B --shortcomp DYJets (--sepflav True/False)

options:
  -h, --help            show this help message and exit
  -p {dilep_sf,ttsemilep_sf,ctag_Wc_sf,ctag_DY_sf,ctag_ttsemilep_sf,ctag_ttdilep_sf}, --phase {dilep_sf,ttsemilep_sf,ctag_Wc_sf,ctag_DY_sf,ctag_ttsemilep_sf,ctag_ttdilep_sf}
                        which phase space
  -i INPUT, --input INPUT
                        input coffea files (str), splitted different files with ','. Wildcard option * available as well.
  -r REF, --ref REF     referance dataset
  -c COMPARED, --compared COMPARED
                        compared datasets, splitted by ,
  --sepflav SEPFLAV     seperate flavour(b/c/light)
  --log                 log on y axis
  -v VARIABLE, --variable VARIABLE
                        variables to plot, splitted by ,. Wildcard option * available as well. Specifying `all` will run through all variables.
  --ext EXT             prefix name
  --com COM             sqrt(s) in TeV
  --shortref SHORTREF   short name for reference dataset for legend
  --shortcomp SHORTCOMP
                        short names for compared datasets for legend, split by ','
   --autorebin AUTOREBIN
                        Rebin the plotting variables by merging N bins in case the current binning is too fine for you 
   --xlabel XLABEL      rename the label for x-axis
```


### Running jupyter remotely
See also https://hackmd.io/GkiNxag0TUmHnnCiqdND1Q#Remote-jupyter

1. On your local machine, edit `.ssh/config`:
```
Host lxplus*
  HostName lxplus7.cern.ch
  User <your-user-name>
  ForwardX11 yes
  ForwardAgent yes
  ForwardX11Trusted yes
Host *_f
  LocalForward localhost:8800 localhost:8800
  ExitOnForwardFailure yes
```
2. Connect to remote with `ssh lxplus_f`
3. Start a jupyter notebook:
```
jupyter notebook --ip=127.0.0.1 --port 8800 --no-browser
```
4. URL for notebook will be printed, copy and open in local browser


