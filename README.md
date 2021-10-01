# BTV Nano comissioning framework


[TOC]


## Introdcution

Contruct by PFNano and coffea, a slimmer, user-friendly framework

To run the tool, run:
```
python runner.py --workflow ${workflow} --json ${json}
```

The current studies are base on the 2017 Rereco dataset with corresponding MC

Purpose | phase space | workflow | json
--- | --- | --- | --- 
btag SF | dilepton ttbar | ttdilep_sf | Rereco17_doublemu.json
btag SF | semilepton ttbar | ttdilep_sf | Rereco17_singlemu.json
ctag SF | b-enriched, dilepton ttbar | ctag_ttdilep_sf | Rereco17_doublemu.json
ctag SF | b-enriched, semilepton ttbar | ctag_semilep_sf | Rereco17_singlemu.json
ctag SF | c-enriched,W+c | ctag_Wc_sf | Rereco17_singlemu.json
ctag SF | udsg-enriched, DY+jets | ctag_DY_sf | ctag_DY.json

### General variables



## Selections

### Btag SF - Dilepton ttbar 

### Btag SF - Semilepton ttbar 

### Ctag SF - b-enriched, Dilepton ttbar

### Ctag SF - b-enriched, Semilepton ttbar 

### Ctag SF - c-enriched, W+c

### Ctag SF - l-enriched, DY+jets

## MC samples 


sample | Rereco17_doublemuon | Rereco17_singlemuon | ctag_DY 
--- | --- | --- | --- 
SingleMuon ||✔||
DoubleMuon |✔||✔|
DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8 |✔|✔|✔|
WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8 |✔|✔|✔|
ST_s-channel_4f_hadronicDecays_TuneCP5_13TeV-amcatnlo-pythia8|✔|✔||
ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8|✔|✔||
ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8|✔|✔||
ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8|✔|✔||
ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8|✔|✔||
ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8|✔|✔||
ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8|✔|✔||
TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8|✔|✔||
TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8|✔|✔||
TTToHadronic_TuneCP5_13TeV-powheg-pythia8|✔|✔||
WW_TuneCP5_13TeV-pythia8|✔|✔||
WZ_TuneCP5_13TeV-pythia8|✔|✔||
ZZ_TuneCP5_13TeV-pythia8|✔|✔||

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
# Replace L47-50 in corrections file & methods
deepcsvb_sf = BTagScaleFactor("data/DeepCSV_94XSF_V5_B_F.csv",BTagScaleFactor.RESHAPE,methods='iterativefit,iterativefit,iterativefit')
#reshape can replace as WP
btag_sf = BTagScaleFactor("data/DeepCSV_94XSF_V5_B_F.csv", "medium")
deepcsvc_sf = "data/DeepCSV_ctagSF_MiniAOD94X_2017_pTincl_v3_2_interp.root"
deepjetb_sf = BTagScaleFactor("data/DeepFlavour_94XSF_V4_B_F.csv",BTagScaleFactor.RESHAPE,methods='iterativefit,iterativefit,iterativefit')
deepjetc_sf = "data/DeepJet_ctagSF_MiniAOD94X_2017_pTincl_v3_2_interp.root"
```
- Lepton SFs:
```python=
# Replace L47-50 in corrections file & methods
deepcsvb_sf = BTagScaleFactor("data/DeepCSV_94XSF_V5_B_F.csv",BTagScaleFactor.RESHAPE,methods='iterativefit,iterativefit,iterativefit')
#reshape can replace as WP
btag_sf = BTagScaleFactor("data/DeepCSV_94XSF_V5_B_F.csv", "medium")
deepcsvc_sf = "data/DeepCSV_ctagSF_MiniAOD94X_2017_pTincl_v3_2_interp.root"
deepjetb_sf = BTagScaleFactor("data/DeepFlavour_94XSF_V4_B_F.csv",BTagScaleFactor.RESHAPE,methods='iterativefit,iterativefit,iterativefit')
deepjetc_sf = "data/DeepJet_ctagSF_MiniAOD94X_2017_pTincl_v3_2_interp.root"
```
## Plotting code

- data/MC comparison code
`python plotdataMC.py --lumi ${lumi} --phase ttdilep --output ttdilep_sf_charm`
