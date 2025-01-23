# Index of the framework 

<!-### Need a figure?  -->

The main ingredients of the framework are wrapped in `src/` directories with supported directories in the root path. 

## `src/data` : custom corrections

Stores the customize corrections used in analysis. e.g. jet probality calibration, custom scale factors...etc. It has a structure similar to [`jsonpog-intergration`](https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/) split by POGs, corrections.
| Type        | File type |  Comments|
| :---:   | :---: | :---: |
| `lumiMasks` |`.json` | Masked good lumi-section used for physics analysis|
| `Prescales` | `.json.` | HLT paths for prescaled triggers|
| `PU`  | `.pkl.gz` or `.histo.root` | Pileup reweight files, matched MC to data| 
| `LSF` | `.histo.root` | Lepton ID/Iso/Reco/Trigger SFs|
| `BTV` | `.csv` or `.root` | b-tagger, c-tagger SFs|
| `JME` | `.txt` | JER, JEC files|
| `JPCalib` | `.root` | Jet probablity calibration, used in LTSV methods|

## `src/utils`: configurations of frameworks

### `histogrammer.py`: collections of hisograms & hisogram writter
### `selection.py`: collections of common selections
### `correction.py`: `coffea` corrections used in analysis

**supported corrections in the current FW**

| Type | corrections        | systematics |  Description| 
| -----  | ---  | :---: | ----- |
| scale<br>variation | event<br>based<br>syst | $\mu_R/\mu_F$<br>0.5,2|  Evaluate renoramlization/ factorization<br>up and down independently<br> based on NanoAOD `LHEScaleWeight` |
| PDF | event<br>based<br>syst | PDF sets| Evaluate PDF uncertainties of <br> `NNPDF31_nnlo_hessian_pdfas`<br>, PDF, alphaS, and PDF+alphaS<br>based on NanoAOD `LHEPDFWeight` |
| parton<br>showering | event<br>based<br>syst | weights|  Evaluate parton showering ISR/FSR<br>uncertainties<br> based on NanoAOD `PSWeight` |
| top pT | SF | event<br>based<br>SF+syst |  no SF | ttbar pT reweighting <br> based on [top PAG]( https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting#TOP_PAG_corrections_based_on_the ) |
| pileup | event<br>based<br>SF+syst | up/down | pileup SF from LUM POG <br>from jsonpog-integration<br> or self derived root file| 
| electron<br>ID/ISO/Reco  | ele<br>SF+syst | up/down | electron ID/Iso/Reco SFs  <br>from EGM using jsonpog-integration | 
| muon<br>ID/Iso   | muon<br>SF+syst  | up/down | muon ID/Iso SFs <br>from MUO in jsonpog-integration | 
| jet<br>veto | --  | -- | veto jets in problematic region<br> suggested by JME in jsonpog-integration | 
| jet<br>prob. | jet trk<br>prob. | apply MC<br>to data | veto jets in problematic <br>region suggested by<br> JME in jsonpog-integration | 
| JER/JEC | jet<br>energy<br>correction<br> or smearing | total<br>up/down<br>energy<br>mass<br>variation | JERC uncertainties from JME, support both txt-like <br>or files from jsonpog-integration | 

### `sample.py`: refined MC sample list for each workflow
### `AK4_parameters.py`: correction, lumi configuration for each campaign
### `plot_utils.py`: plotting utilities
### `array_writter.py`: write out events into root file

## `src/helpers`: functionalities of the framework

### `xsection(_13TeV).py`: xsection diectionary
### `particle*.csv`: particle mass info
### `defintions.py`: definitions of histogram name 
### `BTA_helper.py`: special tools for BTA workflow
### `func.py`: useful functionality
### `update_branch.py`: update missing branch (tagger likelihood ratio)
### `cTagSFReader.py`(deprecated): csv reader of cSF

## `src/workflows`: collections of workflow

Collections of different selections used in commissioning and scale factor. Find the detail in [workflow section](./wf.md).


## `runner.py`: running coffea jobs
## `condor/`: standalone submission

standalone condor submission script with futures executor. See the [details](scaleout.md#standalone-condor-jobs@lxplus/cmsconnect)

## `scripts`: plotting, post-processing, all in one scripts

### `fetch.py`: obtain dataset json file 
### `suball.py`: all in one script to run commissioning/ quick data MC check....
### Output post-processing
    ### `dump_prescale.py`: dump prescale info by `brilcalc`
    ### `dump_processed.py`: dump processed info from output coffea file: lumi & processed files 
    ### `make_template.py`: convert coffea to ROOT hist
    ### `do_hadd.py`: hadd processed root file
    ### `missingFiles.py`: **for customiuzed ** check missing files not include and recreate new submission scripts

### Plotting scripts: 
    ### `comparison.py`: data/data, MC/MC comparison 
    ### `plotdataMC.py`: data/MC comparison
    ### `validation_plot.py`: plot  ROC curve & efficiency
    ### `2Dhistograms.py`: plot 2D histogram from root file
    ### `correction_plot.py`: plot  corelation matrix from root file

## `metadata`: collections of dataset json files

Collections of json file for different campaigns. Split directory by campaign name.

## `testfile`: example coffea files






