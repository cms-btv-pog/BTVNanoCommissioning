
![hackmd-github-sync-badge](https://hackmd.io/OanqNrsDRLqvP6P5dCh9eQ/badge)
[TOC]

# BTVNanoCommissioning
Repository for Commissioning studies in the BTV POG based on (custom) nanoAOD samples

## Requirements
### Setup 
```
# only first time 
git clone git@github.com:cms-btv-pog/BTVNanoCommissioning.git 

# activate enviroment once you have coffea framework 
conda activate coffea
```
### Coffea installation with Miniconda
For installing Miniconda, see also https://hackmd.io/GkiNxag0TUmHnnCiqdND1Q#Local-or-remote
```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
# Run and follow instructions on screen
bash Miniconda3-latest-Linux-x86_64.sh
```
NOTE: always make sure that conda, python, and pip point to local Miniconda installation (`which conda` etc.).

You can either use the default environment`base` or create a new one:
```
# create new environment with python 3.7, e.g. environment of name `coffea`
conda create --name coffea python=3.7
# activate environment `coffea`
conda activate coffea
```
Install coffea, xrootd, and more:
```
pip install git+https://github.com/CoffeaTeam/coffea.git #latest published release with `pip install coffea`
conda install -c conda-forge xrootd
conda install -c conda-forge ca-certificates
conda install -c conda-forge ca-policy-lcg
conda install -c conda-forge dask-jobqueue
conda install -c anaconda bokeh 
conda install -c conda-forge 'fsspec>=0.3.3'
conda install dask
```
### Other installation options for coffea
See https://coffeateam.github.io/coffea/installation.html
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



## Structure
Example worfkflow for ttbar is included. 

Each workflow can be a separate "processor" file, creating the mapping from NanoAOD to
the histograms we need. Workflow processors can be passed to the `runner.py` script 
along with the fileset these should run over. Multiple executors can be chosen 
(for now iterative - one by one, uproot/futures - multiprocessing and dask-slurm). 

To run the example, run:
```
python runner.py --workflow ttcom
```

Example plots can be found in ` make_some_plots.ipynb` though we might want to make
that more automatic in the end.

---------------------
The current studies are base on the 2017 Rereco dataset with corresponding MC

Purpose | phase space | workflow | json
--- | --- | --- | --- 
btag SF | dilepton ttbar | ttdilep_sf | Rereco17_doublemu.json
btag SF | semilepton ttbar | ttsemilep_sf | Rereco17_singlemu.json
ctag SF | b-enriched, dilepton ttbar, muon channel | ctag_ttdilep_sf | Rereco17_doublemu.json
ctag SF | b-enriched, dilepton ttbar, electron channel | ectag_ttdilep_sf | dilep_AK4_runD_ele.json
ctag SF | b-enriched, dilepton ttbar, emu channel | emctag_ttdilep_sf | dilep_AK4_runD_m.json
ctag SF | b-enriched, semilepton ttbar, muon channel | ctag_semilep_sf | Rereco17_singlemu.json
ctag SF | b-enriched, semilepton ttbar | ctag_semilep_sf | Rereco17_singlemu.json
ctag SF | c-enriched,W+c | ctag_Wc_sf, electron channel | Rereco17_singlemu.json
ctag SF | udsg-enriched, DY+jets | ctag_DY_sf, electron channel | ctag_DY.json

### Validation - check different campaigns


Only basic jet selections(PUID, ID, pT, $\eta$) applied. Put the json files with different campaigns

```
python runner.py --workflow valid --json {}
```

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




## Corrections
Corrections are stored in `correction.py`

- Pileup correction:
```python=
# Replace L44 in corrections
with gzip.open("data/corrections.pkl.gz") as fin
# notified the dict in workflows
weights.add('puweight', compiled['2017_pileupweight'](events.Pileup.nPU))
```
- b/c tag SFs:
```python=
# Replace L53-60 in correction files & methods
deepcsvb_sf = BTagScaleFactor("data/DeepCSV_94XSF_V5_B_F.csv",BTagScaleFactor.RESHAPE,methods='iterativefit,iterativefit,iterativefit')
#reshape can replace as WP
btag_sf = BTagScaleFactor("data/DeepCSV_94XSF_V5_B_F.csv", "medium")
deepcsvc_sf = "data/DeepCSV_ctagSF_MiniAOD94X_2017_pTincl_v3_2_interp.root"
deepjetb_sf = BTagScaleFactor("data/DeepFlavour_94XSF_V4_B_F.csv",BTagScaleFactor.RESHAPE,methods='iterativefit,iterativefit,iterativefit')
deepjetc_sf = "data/DeepJet_ctagSF_MiniAOD94X_2017_pTincl_v3_2_interp.root"
```
- Lepton SFs:
```python=
# Replace L47-50 in correction files  
ext.add_weight_sets(["ele_Trig TrigSF data/Ele32_L1DoubleEG_TrigSF_vhcc.histo.root"])
ext.add_weight_sets(["ele_ID EGamma_SF2D data/ElectronIDSF_94X_MVA80WP.histo.root"])
ext.add_weight_sets(["ele_Rereco EGamma_SF2D data/ElectronRecoSF_94X.histo.root"])           
ext.add_weight_sets(["mu_ID NUM_TightID_DEN_genTracks_pt_abseta data/RunBCDEF_SF_ID.histo.root"])
ext.add_weight_sets(["mu_ID_low NUM_TightID_DEN_genTracks_pt_abseta data/RunBCDEF_SF_MuID_lowpT.histo.root"])
ext.add_weight_sets(["mu_Iso NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta data/RunBCDEF_SF_ISO.histo.root"])
```
## Plotting code

- data/MC comparison code

```python=
python plotdataMC.py --lumi ${lumi} --phase ctag_ttdilep_sf --output ctag_ttdilep_sf (--discr zmass --log True/False --data data_runD)
# lumi in /pb
# phase = workflow 
# output coffea file output = hist_$output$.coffea 
# discr = input variables, the defaults are the discriminators, can add multiple variables with space
# log = logorithm on y-axis
# data = data name
```

- data/data, MC/MC comparison

```python=
python comparison.py --phase ctag_ttdilep_sf --output ctag_ttdilep_sf -ref 2017_runB --compared 2017_runC 2017_runD (--discr zmass --log True/False --sepflav True/False)
# phase = workflow 
# output coffea file output = hist_$output$.coffea 
# ref = reference data/MC sample
# comapred = 
# discr = input variables, the defaults are the discriminators, can add multiple variables with space
# log = logorithm on y-axis
# sepflav = separate the jets into different flavor
```

