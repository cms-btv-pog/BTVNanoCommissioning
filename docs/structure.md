## Structure of the framework 

<!-- Need a figure?  -->

The main ingredients of the framework are wrapped in `src/` directories with supported directories in the root path. 

### `src/data` : custom correctiosn

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

### `src/utils`: configurations of frameworks

- `histogrammer.py`: collections of hisograms & hisogram writter
- `selection.py`: collections of common selections
- `correction.py`: `coffea` corrections used in analysis
- `sample.py`: refined MC sample list for each workflow
- `AK4_parameters.py`: correction, lumi configuration for each campaign
- `plot_utils.py`: plotting utilities
- `array_writter.py`: write out events into root file

### `src/helpers`: functionalities of the framework

- `xsection(_13TeV).py`: xsection diectionary
- `particle*.csv`: particle mass info
- `defintions.py`: definitions of histogram name 
- `BTA_helper.py`: special tools for BTA workflow
- `func.py`: useful functionality
- `update_branch.py`: update missing branch (tagger likelihood ratio)
- `cTagSFReader.py`(deprecated): csv reader of cSF

### `src/workflows`: collections of workflow

Collections of different selections used in commissioning and scale factor. Find the detail in [workflow section](./wf.md).


### `runner.py`: running coffea jobs
### `condor/`: standalone submission

standalone condor submission script with futures executor. See the [details](scaleout.md#standalone-condor-jobs@lxplus/cmsconnect)

### `scripts`: plotting, post-processing, all in one scripts

- `fetch.py`: obtain dataset json file 
- `suball.py`: all in one script to run commissioning/ quick data MC check....
- Output post-processing
    - `dump_prescale.py`: dump prescale info by `brilcalc`
    - `dump_processed.py`: dump processed info from output coffea file: lumi & processed files 
    - `make_template.py`: convert coffea to ROOT hist
    - `do_hadd.py`: hadd processed root file
    - `missingFiles.py`: **for customiuzed ** check missing files not include and recreate new submission scripts

- Plotting scripts: 
    - `comparison.py`: data/data, MC/MC comparison 
    - `plotdataMC.py`: data/MC comparison
    - `validation_plot.py`: plot  ROC curve & efficiency
    - `2Dhistograms.py`: plot 2D histogram from root file
    - `correction_plot.py`: plot  corelation matrix from root file

### `metadata`: collections of dataset json files

Collections of json file for different campaigns. Split directory by campaign name.

### `testfile`: example coffea files






