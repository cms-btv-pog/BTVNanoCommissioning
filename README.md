
# BTVNanoCommissioning
[![Linting](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/python_linting.yml/badge.svg)](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/python_linting.yml)
[![btag ttbar](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/ttbar_workflow.yml/badge.svg)](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/ttbar_workflow.yml)
[![ctag ttbar](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/ctag_ttbar_workflow.yml/badge.svg)](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/ctag_ttbar_workflow.yml)
[![ctag DY+jets Workflow](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/ctag_DY_workflow.yml/badge.svg)](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/ctag_DY_workflow.yml)
[![ctag W+c Workflow](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/ctag_Wc_workflow.yml/badge.svg)](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/ctag_Wc_workflow.yml)
[![BTA Workflow](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/BTA_workflow.yml/badge.svg)](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/BTA_workflow.yml)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

Repository for Commissioning studies in the BTV POG based on (custom) nanoAOD samples
Detailed documentation in [btv-wiki](https://btv-wiki.docs.cern.ch/SoftwareAlgorithms/BTVNanoCommissioning/)

## Requirements
### Setup 

:heavy_exclamation_mark: suggested to install under `bash` environment

```
# only first time, including submodules
git clone --recursive git@github.com:cms-btv-pog/BTVNanoCommissioning.git 

# activate enviroment once you have coffea framework 
conda activate btv_coffea
```
### Coffea installation with Micromamba
For installing Micromamba, see [[here](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)]
```
wget -L micro.mamba.pm/install.sh
# Run and follow instructions on screen
bash install.sh
```
NOTE: always make sure that conda, python, and pip point to local micromamba installation (`which conda` etc.).

You can simply create the environment through the existing `test_env.yml` under your micromamba environment using micromamba, and activate it
```
micromamba env create -f test_env.yml 
micromamba activate btv_coffea
```

Once the environment is set up, compile the python package:
```
pip install -e .
pip install -e .[dev] # for developer
```

### Other installation options for coffea
See https://coffeateam.github.io/coffea/installation.html

## Structure
 
Each workflow can be a separate "processor" file, creating the mapping from NanoAOD to
the histograms we need. Workflow processors can be passed to the `runner.py` script 
along with the fileset these should run over. Multiple executors can be chosen 
(for now iterative - one by one, uproot/futures - multiprocessing and dask-slurm). 

To test a small set of files to see whether the workflows run smoothly, run:
```
python runner.py --workflow ttsemilep_sf --json metadata/test_bta_run3.json --campaign Summer22EERun3 --year 2022
```

More options for `runner.py` 
<details><summary>more options</summary>
<p>

```
--wf {validation,ttcom,ttdilep_sf,ttsemilep_sf,emctag_ttdilep_sf,ctag_ttdilep_sf,ectag_ttdilep_sf,ctag_ttsemilep_sf,ectag_ttsemilep_sf,ctag_Wc_sf,ectag_Wc_sf,ctag_DY_sf,ectag_DY_sf,BTA,BTA_addPFMuons,BTA_addAllTracks,BTA_ttbar}, --workflow {validation,ttcom,ttdilep_sf,ttsemilep_sf,emctag_ttdilep_sf,ctag_ttdilep_sf,ectag_ttdilep_sf,ctag_ttsemilep_sf,ectag_ttsemilep_sf,ctag_Wc_sf,ectag_Wc_sf,ctag_DY_sf,ectag_DY_sf,BTA,BTA_addPFMuons,BTA_addAllTracks,BTA_ttbar}
                        Which processor to run
  -o OUTPUT, --output OUTPUT
                        Output histogram filename (default: hists.coffea)
  --samples SAMPLEJSON, --json SAMPLEJSON
                        JSON file containing dataset and file locations
                        (default: dummy_samples.json)
  --year YEAR           Year
  --campaign CAMPAIGN   Dataset campaign, change the corresponding correction
                        files{ "Rereco17_94X","Winter22Run3","Summer22Run3","Summer22EERun3","2018_UL","2017_UL","2016preVFP_UL","2016postVFP_UL"}
  --isSyst              Run with systematics, all, weight_only(no JERC uncertainties included),JERC_split, None(not extract)
  --isArray             Output root files
  --noHist              Not save histogram coffea files
  --overwrite           Overwrite existing files
   --executor {iterative,futures,parsl/slurm,parsl/condor,parsl/condor/naf_lite,dask/condor,dask/condor/brux,dask/slurm,dask/lpc,dask/lxplus,dask/casa}
                        The type of executor to use (default: futures). Other options can be implemented. For
                        example see https://parsl.readthedocs.io/en/stable/userguide/configuring.html-
                        `parsl/slurm` - tested at DESY/Maxwell- `parsl/condor` - tested at DESY, RWTH-
                        `parsl/condor/naf_lite` - tested at DESY- `dask/condor/brux` - tested at BRUX (Brown U)-
                        `dask/slurm` - tested at DESY/Maxwell- `dask/condor` - tested at DESY, RWTH- `dask/lpc` -
                        custom lpc/condor setup (due to write access restrictions)- `dask/lxplus` - custom
                        lxplus/condor setup (due to port restrictions)
  -j WORKERS, --workers WORKERS
                        Number of workers (cores/threads) to use for multi- worker executors (e.g. futures or condor) (default:
                        3)
  -s SCALEOUT, --scaleout SCALEOUT
                        Number of nodes to scale out to if using slurm/condor.
                        Total number of concurrent threads is ``workers x
                        scaleout`` (default: 6)
  --memory MEMORY       Memory used in jobs (in GB) ``(default: 4GB)
  --disk DISK           Disk used in jobs  ``(default: 4GB)
  --voms VOMS           Path to voms proxy, made accessible to worker nodes.
                        By default a copy will be made to $HOME.
  --chunk N             Number of events per process chunk
  --retries N           Number of retries for coffea processor
  --fsize FSIZE         (Specific for dask/lxplus file splitting, default: 50) Numbers of files processed per
                        dask-worker
  --index INDEX         (Specific for dask/lxplus file splitting, default: 0,0) Format:
                        $dict_index_start,$file_index_start,$dict_index_stop,$file_index_stop. Stop indices are
                        optional. $dict_index refers to the index, splitted $dict_index and $file_index with ','
                        $dict_index refers to the sample dictionary of the samples json file. $file_index refers to the N-th batch of files per dask-worker, with its size being defined by the option --index. The job will start (stop) submission from (with) the corresponding indices.
  --validate            Do not process, just check all files are accessible
  --skipbadfiles        Skip bad files.
  --only ONLY           Only process specific dataset or file
  --limit N             Limit to the first N files of each dataset in sample
                        JSON
  --max N               Max number of chunks to run in total
```
</p>
</details>

### Roadmap for running the tool (for commissioning tasks)

1. Is the `.json` file ready? If not, create it following the instructions in the  [Make the json files](#make-the-json-files) section. Please use the correct naming scheme

2. Add the `lumiMask`, correction files (SFs, pileup weight), and JER, JEC files under the dict entry in `utils/AK4_parameters.py`. See details in [Correction files configurations](#correction-files-configurations)

3. If the JERC file `jec_compiled.pkl.gz` is missing in the `data/JME/${campaign}` directory, create it through [Create compiled JERC file](#create-compiled-jerc-filepklgz)

4. If selections and output histogram/arrays need to be changed, modify the dedicated `workflows`

5. Run the workflow with dedicated input and campaign name. Example commands for Run 3 can be found [here](#commands-for-different-phase-space). For first usage, the JERC file needs to be recompiled first, see [Create compiled JERC file](#create-compiled-jerc-filepklgz). You can also specify `--isArray` to store the skimmed root files

6. Fetch the failed files to obtain events that have been processed and events that have to be resubmitted using `scripts/dump_processed.py`. Check the luminosity of the processed dataset used for the plotting script and re-run failed jobs if needed (details in [get procssed info](#get-processed-information))

7. Once you obtain the `.coffea` file(s), you can make plots using the [plotting scripts](#plotting-code), if the xsection for your sample is missing, please add it to `src/BTVNanoCommissioning/helpers/xsection.py`

Check out [notes for developer](https://btv-wiki.docs.cern.ch/SoftwareAlgorithms/BTVNanoCommissioning/#notes-for-developers) for more info!

### Commands for different phase space

After a small test, you can run the full campaign for a dedicated phase space, separately for data and for MC.

#### b-SFs 

<details><summary>details</summary>
<p>

- Dileptonic ttbar phase space : check performance for btag SFs, emu channel

```
 python runner.py --workflow ttdilep_sf --json metadata/data_Summer22_Run3_2022_em_BTV_Comm_v2_NanoV12_noPF.json  --campaign Summer22Run3 --year 2022 (--executor ${scaleout_site}) 
```

- Semileptonic ttbar phase space : check performance for btag SFs, muon channel

```
python runner.py --workflow ttsemilep_sf --json metadata/data_Summer22_Run3_2022_mu_BTV_Comm_v2_NanoV12_noPF.json --campaign Summer22Run3 --year 2022 (--executor ${scaleout_site})
```

</p>
</details>

#### c-SFs
<details><summary>details</summary>
<p>

- Dileptonic ttbar phase space : check performance for charm SFs, bjets enriched SFs, muon channel

```
python runner.py --workflow ctag_ttdilep_sf --json metadata/data_Summer22_Run3_2022_mu_BTV_Comm_v2_NanoV12_noPF.json  --campaign Summer22Run3 --year 2022(--executor ${scaleout_site})
```


- Semileptonic ttbar phase space : check performance for charm SFs, bjets enriched SFs, muon channel

```
python runner.py --workflow ctag_ttsemilep_sf --json metadata/data_Summer22_Run3_2022_mu_BTV_Comm_v2_NanoV12_noPF.json  --campaign Summer22Run3 --year 2022(--executor ${scaleout_site})
```

- W+c phase space : check performance for charm SFs, cjets enriched SFs, muon  channel

```
python runner.py --workflow ctag_Wc_sf --json metadata/data_Summer22_Run3_2022_mu_BTV_Comm_v2_NanoV12_noPF.json  --campaign Summer22Run3 --year 2022(--executor ${scaleout_site})
```

- DY phase space : check performance for charm SFs, light jets enriched SFs, muon channel

```
python runner.py --workflow ctag_DY_sf --json metadata/data_Summer22_Run3_2022_mu_BTV_Comm_v2_NanoV12_noPF.json  --campaign Summer22Run3 --year 2022(--executor ${scaleout_site})
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


#### BTA - BTagAnalyzer Ntuple producer

Based on Congqiao's [development](notebooks/BTA_array_producer.ipynb) to produce BTA ntuples based on PFNano.

:exclamation: Only the newest version [BTV_Run3_2022_Comm_v2](https://github.com/cms-jet/PFNano/tree/13_0_7_from124MiniAOD) ntuples work. Example files are given in [this](metadata/test_bta_run3.json) json. Optimize the chunksize(`--chunk`) in terms of the memory usage. This depends on sample, if the sample has huge jet collection/b-c hardons. The more info you store, the more memory you need. I would suggest to test with `iterative` to estimate the size.

<details><summary>details</summary>
<p>

Run with the nominal `BTA` workflow to include the basic event variables, jet observables, and GEN-level quarks, hadrons, leptons, and V0 variables. 
```
python runner.py --wf BTA --json metadata/test_bta_run3.json --campaign Summer22EERun3 --isJERC
```

Run with the `BTA_addPFMuons` workflow to additionally include the `PFMuon` and `TrkInc` collection, used by the b-tag SF derivation with the QCD(μ) methods.
```
python runner.py --wf BTA_addPFMuons --json metadata/test_bta_run3.json --campaign Summer22EERun3 --isJERC
```

Run with the `BTA_addAllTracks` workflow to additionally include the `Tracks` collection, used by the JP variable calibration.
```
python runner.py --wf BTA_addAllTracks --json metadata/test_bta_run3.json --campaign Summer22EERun3 --isJERC
```

</p>
</details>

## Scale-out 

Scale out can be notoriously tricky between different sites. Coffea's integration of `slurm` and `dask`
makes this quite a bit easier and for some sites the ``native'' implementation is sufficient, e.g Condor@DESY.
However, some sites have certain restrictions for various reasons, in particular Condor @CERN and @FNAL. The scaleout scheme is named as follows: `$cluster_schedule_system/scheduler/site`. The existing sites are documented in [sites configuration](#sites-configuration-with-daskparsl-schedular) while [standalone condor submission](#standalone-condor-jobslxpluscmsconnect) is possible and strongly suggested when working on lxplus. 


Memory usage is also useful to adapt to cluster. Check the memory by calling  `memory_usage_psutil()` from `helpers.func.memory_usage_psutil` to optimize job size. Example with `ectag_Wc_sf` summarized below.
 Type        |Array+Hist |  Hist only| Array Only|
| :---:   | :---: | :---: | :---: |
DoubleMuon (BTA,BTV_Comm_v2)| 1243MB |	848MB	|1249MB|
DoubleMuon (PFCands, BTV_Comm_v1)|1650MB	|1274MB	|1632MB|
DoubleMuon (Nano_v11)|1183MB|	630MB	|1180MB|
WJets_inc (BTA,BTV_Comm_v2)| 1243MB	|848MB	|1249MB|
WJets_inc (PFCands, BTV_Comm_v1)|1650MB	|1274MB	|1632MB
WJets_inc (Nano_v11)|1183MB	|630MB	|1180MB|


### Sites configuration with dask/parsl schedular

<details><summary>details</summary>
<p>

#### Condor@FNAL (CMSLPC)
Follow setup instructions at https://github.com/CoffeaTeam/lpcjobqueue. After starting 
the singularity container run with 
```bash
python runner.py --wf ttcom --executor dask/lpc
```

#### Condor@CERN (lxplus)
Only one port is available per node, so its possible one has to try different nodes until hitting
one with `8786` being open. Other than that, no additional configurations should be necessary.

```bash
python runner.py --wf ttcom --executor dask/lxplus
```

jobs automatically split to 50 files per jobs to avoid job failure due to crowded cluster on lxplus with the naming scheme `hist_$workflow_$json_$dictindex_$fileindex.coffea`. The `.coffea` files can be then combined at plotting level


:exclamation: The optimal scaleout options on lxplus are `-s 50 --chunk 50000`

To deal with unstable condor cluster and dask worker on lxplus, you can resubmit failure jobs via `--index $dictindex,$fileindex` option. `$dictindex` refers to the index in the `.json dict`, `$fileindex` refers to the index of the file list split to 50 files per dask-worker. The total number of files of each dict can be computed by `math.ceil(len($filelist)/50)` The job will start from the corresponding indices.

#### Coffea-casa (Nebraska AF)
Coffea-casa is a JupyterHub based analysis-facility hosted at Nebraska. For more information and setup instuctions see
https://coffea-casa.readthedocs.io/en/latest/cc_user.html

After setting up and checking out this repository (either via the online terminal or git widget utility run with
```bash
python runner.py --wf ttcom --executor dask/casa
```
Authentication is handled automatically via login auth token instead of a proxy. File paths need to replace xrootd redirector with "xcache", `runner.py` does this automatically.


#### Condor@DESY 
```bash
python runner.py --wf ttcom --executor dask/condor(parsl/condor)
```

#### Maxwell@DESY 
```bash
python runner.py --wf ttcom --executor parsl/slurm
```

</p>
</details>

### Standalone condor jobs@lxplus/cmsconnect

:heavy_exclamation_mark: :heavy_exclamation_mark: :heavy_exclamation_mark: Strongly suggest to use this in lxplus :heavy_exclamation_mark: :heavy_exclamation_mark: :heavy_exclamation_mark: 
You have the option to run the framework through "standalone condor jobs", bypassing the native coffea-supported job submission system. Within each job you submit, a standalone script will execute the following on the worker node:

 - Set up a minimal required Python environment.
 - Retrieve the BTVNanoCommissioning repository, either from a git link or transferred locally.
 - Launch the `python runner.py ...` command to execute the coffea framework in the iterative executor mode.
 
This utility is currently adapted for the lxplus and cmsconnect condor systems. To generate jobs for launching, replace `python runner.py` with `python condor/submitter.py`, append the existing arguments, and add the following arguments in addition:

 - `--jobName`: Specify the desired condor job name. A dedicated folder will be generated, including all submission-related files.
 - `--outputXrootdDir`: Indicate the XRootD directory's path (starting with `root://`) where the produced .coffea (and .root) files from each worker node will be transferred to.
 - `--condorFileSize`: Define the number of files to process per condor job (default is 50). The input file list will be divided based on this count.
 - `--remoteRepo` (optional, but recommended): Specify the path and branch of the remote repository to download the BTVNanoCommissioning repository. If not specified, the local directory will be packed and transferred as the condor input, potentially leading to higher loads for condor transfers. Use the format e.g. `--remoteRepo 'https://github.com/cms-btv-pog/BTVNanoCommissioning.git -b master'`.

After executing the command, a new folder will be created, preparing the submission. Follow the on-screen instructions and utilize `condor_submit ...` to submit the jdl file. The output will be transferred to the designated XRootD destination.

<details><summary>Frequent issues for standalone condor jobs submission
</summary>
<p>

1. CMS Connect provides a condor interface where one can submit jobs to all resources available in the CMS Global Pool. See [WorkBookCMSConnect Twiki](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCMSConnect#Requesting_different_Operative_S) for the instructions if you use it for the first time.
2. The submitted jobs are of the kind which requires a proper setup of the X509 proxy, to use the XRootD service to access and store data. In the generated `.jdl` file, you may see a line configured for this purpose `use_x509userproxy = true`. If you have not submitted jobs of this kind on lxplus condor, we recommend you to add a line
   ```bash
   export X509_USER_PROXY=$HOME/x509up_u`id -u`
   ```
   to `.bashrc` and run it so the proxy file will be stored in your AFS folder instead of in your `/tmp/USERNAME` folder. For submission on cmsconnect, no specific action is required.

</p>
</details>


## Make the dataset json files

Use `fetch.py` in folder `scripts/` to obtain your samples json files. You can create `$input_list` ,which can be a list of datasets taken from CMS DAS , and create the json contains `dataset_name:[filelist]`. One can specify the local path in that input list for samples not published in CMS DAS.
`$output_json_name$` is the name of your output samples json file.

The `--whitelist_sites, --blacklist_sites` are considered for fetch dataset if multiple sites are available

```
## File publish in DAS
python fetch.py --input ${input_DAS_list} --output ${output_json_name} (--xrd {prefix_forsite})

## Not publish case, specify site by --xrd prefix
python fetch.py --input ${input_list} --output ${output_json_name} --xrd {prefix_forsite}
# where the input list should contains
$DATASET_NAME $PATH_TO_FILE
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


## Get processed information

Get the run & luminosity information for the processed events from the coffea output files. When you use `--skipbadfiles`, the submission will ignore files not accesible(or time out) by `xrootd`. This script helps you to dump the processed luminosity into a json file which can be calculated by `brilcalc` tool and provide a list of failed lumi sections by comparing the original json input to the one from the `.coffea` files.


```bash
# all is default, dump lumi and failed files, if run -t lumi only case. no json file need to be specified
python scripts/dump_processed.py -c $COFFEA_FILES -n $OUTPUT_NAME (-j $ORIGINAL_JSON -t [all,lumi,failed])
```


## Correction files configurations
:heavy_exclamation_mark:  If the correction files are not supported yet by jsonpog-integration, you can still try with custom input data.

### Options with custom input data 

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


<details><summary>Take `Rereco17_94X` as an example.</summary>
<p>

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

</p>
</details>

### Use central maintained jsonpog-integration 
The official correction files collected in [jsonpog-integration](https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration) is updated by POG except `lumiMask` and `JME` still updated by maintainer. No longer to request input files in the `correction_config`.  

<details><summary>See the example with `2017_UL`.</summary>
<p>

```python
  "2017_UL": {
        # Same with custom config
        "lumiMask": "Cert_294927-306462_13TeV_UL2017_Collisions17_MuonJSON.txt",
        "JME": "jec_compiled.pkl.gz",
        # no config need to be specify for PU weights
        "PU": None,
        # Btag SFs - specify $TAGGER : $TYPE-> find [$TAGGER_$TYPE] in json file
        "BTV": {"deepCSV": "shape", "deepJet": "shape"},
        "roccor": None,
         # JMAR, IDs from JME- Following the scheme: "${SF_name}": "${WP}"
        "JMAR": {"PUJetID_eff": "L"},
        "LSF": {
        # Electron SF - Following the scheme: "${SF_name} ${SF_map} ${year}": "${WP}"
        # https://github.com/cms-egamma/cms-egamma-docs/blob/master/docs/EgammaSFJSON.md
            "ele_ID 2017 UL-Electron-ID-SF": "wp90iso",
            "ele_Reco 2017 UL-Electron-ID-SF": "RecoAbove20",

        # Muon SF - Following the scheme: "${SF_name} ${year}": "${WP}"
        # WPs : ['NUM_GlobalMuons_DEN_genTracks', 'NUM_HighPtID_DEN_TrackerMuons', 'NUM_HighPtID_DEN_genTracks', 'NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight', 'NUM_LooseID_DEN_TrackerMuons', 'NUM_LooseID_DEN_genTracks', 'NUM_LooseRelIso_DEN_LooseID', 'NUM_LooseRelIso_DEN_MediumID', 'NUM_LooseRelIso_DEN_MediumPromptID', 'NUM_LooseRelIso_DEN_TightIDandIPCut', 'NUM_LooseRelTkIso_DEN_HighPtIDandIPCut', 'NUM_LooseRelTkIso_DEN_TrkHighPtIDandIPCut', 'NUM_MediumID_DEN_TrackerMuons', 'NUM_MediumID_DEN_genTracks', 'NUM_MediumPromptID_DEN_TrackerMuons', 'NUM_MediumPromptID_DEN_genTracks', 'NUM_Mu50_or_OldMu100_or_TkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose', 'NUM_SoftID_DEN_TrackerMuons', 'NUM_SoftID_DEN_genTracks', 'NUM_TightID_DEN_TrackerMuons', 'NUM_TightID_DEN_genTracks', 'NUM_TightRelIso_DEN_MediumID', 'NUM_TightRelIso_DEN_MediumPromptID', 'NUM_TightRelIso_DEN_TightIDandIPCut', 'NUM_TightRelTkIso_DEN_HighPtIDandIPCut', 'NUM_TightRelTkIso_DEN_TrkHighPtIDandIPCut', 'NUM_TrackerMuons_DEN_genTracks', 'NUM_TrkHighPtID_DEN_TrackerMuons', 'NUM_TrkHighPtID_DEN_genTracks']

            "mu_Reco 2017_UL": "NUM_TrackerMuons_DEN_genTracks",
            "mu_HLT 2017_UL": "NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight",
            "mu_ID 2017_UL": "NUM_TightID_DEN_TrackerMuons",
            "mu_Iso 2017_UL": "NUM_TightRelIso_DEN_TightIDandIPCut",
        },
    },
```

</p>
</details>

## Create compiled JERC file(`pkl.gz`)

:exclamation: In case existing correction file doesn't work for you due to the incompatibility of `cloudpickle` in different python versions. Please recompile the file to get new pickle file.

Under `compile_jec.py` you need to create dedicated jet factory files with different campaigns. Following the name scheme with `mc` for MC and `data${run}` for data.

Compile correction pickle files for a specific JEC campaign by changing the dict of jet_factory, and define the MC campaign and the output file name by passing it as arguments to the python script:

```
python -m BTVNanoCommissioning.utils.compile_jec ${campaign} jec_compiled
e.g. python -m BTVNanoCommissioning.utils.compile_jec Summer22Run3 jec_compiled
```


## Plotting code

- data/MC comparisons
:exclamation_mark: If using wildcard for input, do not forget the quoatation marks! (see 2nd example below)

You can specify `-v all` to plot all the variables in the `coffea` file, or use wildcard options (e.g. `-v "*DeepJet*"` for the input variables containing `DeepJet`)

:new: non-uniform rebinning is possible, specify the bins with  list of edges `--autorebin 50,80,81,82,83,100.5`

```bash
python scripts/plotdataMC.py -i a.coffea,b.coffea --lumi 41500 -p ttdilep_sf -v z_mass,z_pt  
python scripts/plotdataMC.py -i "test*.coffea" --lumi 41500 -p ttdilep_sf -v z_mass,z_pt # with wildcard option need ""
```

<details><summary>more arguments</summary>
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

</details>
</p>

- data/data, MC/MC comparisons

You can specify `-v all` to plot all the variables in the `coffea` file, or use wildcard options (e.g. `-v "*DeepJet*"` for the input variables containing `DeepJet`)
:exclamation_mark: If using wildcard for input, do not forget the quoatation marks! (see 2nd example below)

```bash
# with merge map, compare ttbar with data
python scripts/comparison.py -i "*.coffea" --mergemap '{"ttbar": ["TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8","TTto4Q_TuneCP5_13p6TeV_powheg-pythia8","TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8],"data":["MuonRun2022C-27Jun2023-v1","MuonRun2022D-27Jun2023-v1"]}' -r ttbar -c data -v mu_pt  -p ttdilep_sf
# if no  mergemap, take the key name directly
python scripts/comparison.py -i datac.coffea,datad.coffea -r MuonRun2022C-27Jun2023-v1 -c MuonRun2022D-27Jun2023-v1 -v mu_pt  -p ttdilep_sf

```

<details><summary>more arguments</summary>
<p>

 ```
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
  --mergemap MERGEMAP
                        Group list of sample(keys in coffea) as reference/compare set as dictionary format. Keys would be the new lables of the group
   --autorebin AUTOREBIN
                        Rebin the plotting variables, input `int` or `list`. int: merge N bins. list of number: rebin edges(non-uniform bin is possible)
   --xlabel XLABEL      rename the label for x-axis
   --ylabel YLABEL      rename the label for y-axis
   --norm               compare shape, normalized yield to reference
   --xrange XRANGE       custom x-range, --xrange xmin,xmax
   --flow FLOW 
                        str, optional {None, 'show', 'sum'} Whether plot the under/overflow bin. If 'show', add additional under/overflow bin. If 'sum', add the under/overflow bin content to first/last bin.
```

</details>
</p>

## Store histograms from coffea file

Use `scripts/make_template.py` to dump 1D/2D histogram from `.coffea` to `TH1D/TH2D` with hist. MC histograms can be reweighted to according to luminosity value given via `--lumi`. You can also merge several files 

```python
python scripts/make_template.py -i "testfile/*.coffea" --lumi 7650 -o test.root -v mujet_pt -a '{"flav":0,"osss":"sum"}'
```

<details><summary>more arguments</summary>
<p>

```
  -i INPUT, --input INPUT
                        Input coffea file(s)
  -v VARIABLE, --variable VARIABLE
                        Variables to store (histogram name)
  -a AXIS, --axis AXIS  dict, put the slicing of histogram, specify 'sum' option as string
  --lumi LUMI           Luminosity in /pb
  -o OUTPUT, --output OUTPUT
                        output root file name
  --mergemap MERGEMAP   Specify mergemap as dict, '{merge1:[dataset1,dataset2]...}' Also works with the json file with dict
```

</details>
</p>



<details><summary>mergemap example</summary>
<p>

```json
{
    "WJets": ["WJetsToLNu_TuneCP5_13p6TeV-madgraphMLM-pythia8"],
    "VV": [ "WW_TuneCP5_13p6TeV-pythia8", "WZ_TuneCP5_13p6TeV-pythia8", "ZZ_TuneCP5_13p6TeV-pythia8"],
    "TT": [ "TTTo2J1L1Nu_CP5_13p6TeV_powheg-pythia8", "TTTo2L2Nu_CP5_13p6TeV_powheg-pythia8"],
    "ST":[ "TBbarQ_t-channel_4FS_CP5_13p6TeV_powheg-madspin-pythia8", "TbarWplus_DR_AtLeastOneLepton_CP5_13p6TeV_powheg-pythia8", "TbarBQ_t-channel_4FS_CP5_13p6TeV_powheg-madspin-pythia8", "TWminus_DR_AtLeastOneLepton_CP5_13p6TeV_powheg-pythia8"],
"data":[ "Muon_Run2022C-PromptReco-v1", "SingleMuon_Run2022C-PromptReco-v1", "Muon_Run2022D-PromptReco-v1", "Muon_Run2022D-PromptReco-v2"]
}
```

</p>
</details>



## Notes for developers
The BTV tutorial for coffea part is under `notebooks` and the template to construct new workflow is `src/BTVNanoCommissioning/workflows/example.py`
Here are some tips provided for developers working on their forked version of this repository. Also some useful git commands can be found [here](https://codimd.web.cern.ch/wY3IrOBBT3i3GXIQxLMWPA)
### Setup CI pipeline for fork branch
Since the CI pipelines involve reading files via `xrootd` and access gitlab.cern.ch, you need to save some secrets in your forked directory. 

Yout can find the secret configuration in the direcotry : `Settings>>Secrets>>Actions`, and create the following secrets:

- `GIT_CERN_SSH_PRIVATE`: 
  1. Create a ssh key pair with `ssh-keygen -t rsa -b 4096` (do not overwrite with your local one), add the public key to your CERN gitlab account
  2. Copy the private key to the entry
- `GRID_PASSWORD`: Add your grid password to the entry.
- `GRID_USERCERT` & `GRID_USERKEY`:  Encrypt your grid user certification `base64 -i ~/.globus/userkey.pem | awk NF=NF RS= OFS=` and `base64 -i ~/.globus/usercert.pem | awk NF=NF RS= OFS=` and copy the output to the entry. 

Special commit head messages could run different commands in actions (add the flag in front of your commit)
The default configureation is doing 
```
python runner.py --workflow emctag_ttdilep_sf --json metadata/test_bta_run3.json --limit 1 --executor iterative --campaign Summer22Run3 --isArray --isSyst all
```

- `[skip ci]`: not running ci at all in the commit message
- `ci:skip array` : remove `--isArray` option
- `ci:skip syst` : remove `--isSyst all` option
- `ci:JERC_split` : change systematic option to split JERC uncertainty sources `--isSyst JERC_split`
- `ci:weight_only` : change systematic option to weight only variations `--isSyst weight_only`
 
### Running jupyter remotely
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
