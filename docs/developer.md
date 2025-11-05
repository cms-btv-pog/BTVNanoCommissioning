# For developers: Add new workflow


The BTV tutorial for the coffea part at [`notebooks/BTV_commissioning_tutorial-coffea.ipynb`](https://github.com/cms-btv-pog/BTVNanoCommissioning/tree/master/notebooks/BTV_commissioning_tutorial-coffea.ipynb) and the template to construct new workflow is located at [`src/BTVNanoCommissioning/workflows/example.py`](https://github.com/cms-btv-pog/BTVNanoCommissioning/blob/master/src/BTVNanoCommissioning/workflows/example.py)

## 0. Add new workflow info to `workflows/__init__.py` 


```python
# if the workflow is in a new files, you need to import your workflow.py
from BTVNanoCommissioning.workflows.new_workflow import (
  NanoProcessor as new_WORKFLOWProcessor
)
# And then include the processor into the modules with the name of workflow. The name is used when specifying --workflow when running the jobs
workflows["name_workflow"] = new_WORKFLOWProcessor
# IF the workflow is based on the modifier to existing workflow, put the selectionModifier used in the existing workflow
workflows["ctag_ttsemilep_sf"] = partial(
    CTAGWcTTValidSFProcessor, selectionModifier="semittM"
)
```
Notice that if you are working on a WP SFs, please put **WP** in the name.

## 1. Add histogram collections to `utils/histogramming/histograms`

The framework uses [`hist`](https://hist.readthedocs.io/en/latest/) histograms, which can be easily to converted to ROOT histogram by `uproot` or to `numpy` histograms. For a quickstart to hist see [this page](https://hist.readthedocs.io/en/latest/user-guide/quickstart.html).

There are a common axes located in `utils/histogramming/axes/common.py`, which are used by many of the workflows. You can define your own axes collection in the same directory and then add it to `utils/histogramming/hist_helpers.py` to make it available for the histogrammer. The `name` and dictionary key of the axis should be consistent with axes defined in the histogram.

An example collection of axes:
```python
axes = {
    "flav": hist.axis.IntCategory([0, 1, 4, 5, 6], name="flav", label="Genflavour"),
    "syst": hist.axis.StrCategory([], name="syst", growth=True),
    "pt": hist.axis.Regular(60, 0, 300, name="pt", label=" $p_{T}$ [GeV]"),
}
```

Regardless of the axis or histogram collection, the histograms need to be contained in a `dict`, and each histogram should at least contain as its first axis a **syst_axis** and a **Hist.storage.Weight()** as its last axis. The key for the histogram is suggested to use the format `$OBJECT_$VAR`, with `OBJECT` being an object defined during the workflow and `$VAR` a variable for that object, for example:

```python
_hist_dict["mujet_pt"] = Hist.Hist(
            syst_axis, flav_axis, dr_axis, Hist.storage.Weight()
        )  # create cutstomize histogram
```

After choosing and defining the axis and histogram collections you want to use in your workflow, use the [histogrammer](https://github.com/cms-btv-pog/BTVNanoCommissioning/blob/5c56c0affa6209c7b3f99855ed52c41fe5e225ec/src/BTVNanoCommissioning/utils/histogramming/histogrammer.py#L7) function to create the output histograms. Note that there are also common histogram collections for creating common histograms ([common.py](https://github.com/cms-btv-pog/BTVNanoCommissioning/blob/master/src/BTVNanoCommissioning/utils/histogramming/histograms/common.py)) and 4-vector histograms for a given `obj_list` ([fourvec.py](https://github.com/cms-btv-pog/BTVNanoCommissioning/blob/master/src/BTVNanoCommissioning/utils/histogramming/histograms/fourvec.py)). An example histogrammer call could be:
```python
output = histogrammer(
    jet_fields=events.Jet.fields, # Pass the jet fields to check available jet variables
    obj_list=["posl", "negl", "dilep", "jet0"], # Pass the objects to define 4-vector quantities for
    hist_collections=["common", "fourvec", "DY"], # Choose the histogram collections to use
    include_m=isMu, # Pass a collection specific kwarg
    # This call uses the common axis collection
)
```

##  2. Selections: Implement selections on events (`workflow/`)

Create `boolean` arrays along event axis. Also check whether some common selctions already in `utils/selection.py`

```python
      ## HLT- put trigger paths
      triggers = [
          "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
          "Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
          "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
      ]
      req_trig = HLT_helper(events, triggers)

      ##### Add some selections
      ## Muon cuts
      # muon twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2

      muon_sel = (events.Muon.pt > 15) & (mu_idiso(events, self._campaign)) # applied selection on pT, and `mu_idiso` is predefined selection in `selection.py` which refers to cut-based tight muon ID+Iso
      event_mu = events.Muon[muon_sel] # Pruned the muon collections with the selection
      req_muon = ak.num(event_mu.pt) == 1 # Check each event has exact one muon
      .... # other selections
      ## Apply all selections
      event_level = (
          req_trig & req_lumi & req_jet & req_muon & req_ele & req_leadlep_pt
      )
```

In case you are modifying an existing workflow, you need to add to  `__init__` and in the selections create a `selectionModifier`:

```python
# in init
self.selMod = selectionModifier
# In selection 
if self.selMod=="WcM":
  event_level = req_trig & req_lumi & req_jet & req_muon & req_ele & req_leadlep_pt& req_Wc
```

##  3. Selected objects: Pruned objects with reduced event_level
Store the selected objects to event-based arrays. If using the common [histogram writer](https://github.com/cms-btv-pog/BTVNanoCommissioning/blob/5c56c0affa6209c7b3f99855ed52c41fe5e225ec/src/BTVNanoCommissioning/utils/histogramming/histogrammer.py#L63), the selected object must contain **Sel**, and the muon-enriched jet and soft muon are called **MuonJet** and **SoftMu**, respectively. You can then use the defined objects to derive new variables, ie.

```python
  # Keep the structure of events and pruned the object size
  pruned_ev = events[event_level]
  pruned_ev["SelJet"] = event_jet[event_level][:, 0]
  pruned_ev["SelMuon"] = event_mu[event_level][:, 0]
  pruned_ev["mujet_ptratio"] = event_mu[event_level].pt / pruned_ev.SelJet.pt # notice that the cross-object need to be created specificaly
  pruned_ev["mujet_dr"] =  event_mu[event_level].delta_r(pruned_ev.SelJet) 
```


 The pruned information is then stored into the defined histograms and the output arrays, and used to evaluate the event weights. In case you have a custom object for [corrections](#add-additional-weight-or-uncertainty-information), or a new common [variable](#add-new-common-variables) that you need to add, please go to dedicated section.
 
See details below for the usage of `pruned_ev`

<details><summary>Output section</summary>
<p>

```python
####################
#     Output       #
####################
# Configure SFs - read pruned objects from the pruned_ev and apply SFs and call the systematics
weights = weight_manager(pruned_ev, self.SF_map, self.isSyst)
# Configure systematics shifts 
if shift_name is None:
    systematics = ["nominal"] + list(weights.variations) # nominal + weight variation systematics
else:
    systematics = [shift_name] # JES/JER systematics

# Fill the weight to output arrys
if not isRealData:
    pruned_ev["weight"] = weights.weight()
    for ind_wei in weights.weightStatistics.keys():
        pruned_ev[f"{ind_wei}_weight"] = weights.partial_weight(
            include=[ind_wei]
        )
# Configure histograms- fill the histograms with pruned objects
if not self.noHist:
    output = histo_writter(
        pruned_ev, output, weights, systematics, self.isSyst, self.SF_map
    )
# Output arrays - store the pruned objects in the output arrays
if self.isArray:
    array_writer(self, pruned_ev, events, weights, systematics,dataset, isRealData)
  ```

</p>
</details>


## 4. Setup CI pipeline `.github/workflow`

The actions are checking the changes would break the framework. The actions are collected in `.github/workflow`
You can simply include a workflow by adding the entries with name

```yaml
- name:  semileptonic + c ttbar workflows with correctionlib
      run: |
        string=$(git log -1 --pretty=format:'%s')
        if [[ $string == *"ci:skip array"* ]]; then
        opts=$(echo "$opts" | sed 's/--isArray //g')
        fi
        if [[ $string == *"ci:skip syst"* ]]; then
            opts=$(echo "$opts" | sed 's/--isSyst all//g')
        elif [[ $string == *"ci:JERC_split"* ]]; then
            opts=$(echo "$opts" | sed 's/--isSyst all/--isSyst JERC_split/g')
        elif [[ $string == *"ci:weight_only"* ]]; then
            opts=$(echo "$opts" | sed 's/--isSyst all/--isSyst weight_only/g') 
        fi
        python runner.py --workflow c_ttsemilep_sf --json metadata/test_bta_run3.json --limit 1 --executor iterative --campaign Summer23 --year 2023  $opts
```

Special commit head messages could run different commands in actions (add the flag in front of your commit)
The default configureation is doing 
```python
python runner.py --workflow emctag_ttdilep_sf --json metadata/test_bta_run3.json --limit 1 --executor iterative --campaign Summer23 --isArray --isSyst all
```

- `[skip ci]`: not running ci at all in the commit message
- `ci:skip array` : remove `--isArray` option
- `ci:skip syst` : remove `--isSyst all` option
- `ci:JERC_split` : change systematic option to split JERC uncertainty sources `--isSyst JERC_split`
- `ci:weight_only` : change systematic option to weight only variations `--isSyst weight_only`

<details><summary>Set CI in your github account</summary>
<p>

Since the CI pipelines involve reading files via `xrootd` and access gitlab.cern.ch, you need to save some secrets in your forked directory. 

Yout can find the secret configuration in the direcotry : `Settings>>Secrets>>Actions`, and create the following secrets:

- `GIT_CERN_SSH_PRIVATE`: 
  1. Create a ssh key pair with `ssh-keygen -t rsa -b 4096` (do not overwrite with your local one), add the public key to your CERN gitlab account
  2. Copy the private key to the entry
- `GRID_PASSWORD`: Add your grid password to the entry.
- `GRID_USERCERT` & `GRID_USERKEY`:  Encrypt your grid user certification `base64 -i ~/.globus/userkey.pem | awk NF=NF RS= OFS=` and `base64 -i ~/.globus/usercert.pem | awk NF=NF RS= OFS=` and copy the output to the entry. 

</p>
</details>


## 5. Refine used MC as input `sample.py`
The `sample.py` collects the samples (dataset name) used in the workflow. This collections are use to create the dataset json file.
- `data` : data sample (MuonEG, Muon0....)
- `MC`: main MC used for the workflow
- `minor_MC` : minor MC samples use for background events
- `syst_MC`: systematic MC samples (TTbar sample mass, Hdamp ... variations)

Here's the example for BTA_ttbar
```python
"BTA_ttbar": {
        "data": ["MuonEG"],
        "MC": ["TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8"],
        "minor_MC": [
            "TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",
            "TWminusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
            "TbarWplusto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
            "TWminustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",
            "TbarWplustoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",
            "TbarBQ_t-channel_4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8",
            "TBbarQ_t-channel_4FS_TuneCP5_13p6TeV_powheg-madspin-pythia8",
            "WWto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
            "ZZto2L2Q_TuneCP5_13p6TeV_powheg-pythia8",
            "WZto3LNu_TuneCP5_13p6TeV_powheg-pythia8",
            "WZtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",
        ],
        "syst_MC": [
            "TTto2L2Nu_MT-171p5_TuneCP5_13p6TeV_powheg-pythia8",
            "TTto2L2Nu_MT-175p5_TuneCP5_13p6TeV_powheg-pythia8",
            "TTto2L2Nu_Hdamp-158_TuneCP5_13p6TeV_powheg-pythia8",
            "TTto2L2Nu_Hdamp-418_TuneCP5_13p6TeV_powheg-pythia8",
            "TTto2L2Nu_TuneCP5Down_13p6TeV_powheg-pythia8",
            "TTto2L2Nu_TuneCP5Up_13p6TeV_powheg-pythia8",
        ],
    },
```

## Optional changes
### Add workflow to `scripts/suball.py` 
The `suball.py` summarize the steps to obtain the result.
In case your task requires to run several workflows, you can wrapped them as `dict` of the workflows
```python
scheme = {
        # scale factor workflows
        "SF": ["BTA_ttbar", "BTA_addPFMuons"],
        # Use for prompt data MC checks for analysis
        "Validation": ["ttdilep_sf", "ctag_Wc_sf"],
        # commissioning workflows
        "default_comissioning": [
            "ttdilep_sf",
            "ttsemilep_sf",
            "ctag_Wc_sf",
            "ctag_DY_sf",
            "QCD_sf",
            "QCD_mu_sf"
        ],
    }
```
### Add new common variables in `helper/definition.py`

In the `definition.py` we collect the axis definition, name and label of tagger scores/input variables 
```python
disc_list=[....] # tagger score
definitions_dict = {
    "DeepCSV_jetNSelectedTracks": # name used in tree 
    {
        "displayname": "Jet N Selected Tracks", # axis name
        "manual_ranges": [0.0, 25],
        "ylabel_text": "Jets",
        "format_unit": "2f",
        "format_unit_digits": 2,
        "bins": 25,
        "inputVar_units": None,
    },
    ...
}
```
### Additional corrections and uncertainty variations not in the framework
The corrections are collected in `utils/correction.py`.  There are two types of the variation: weight varations, i.e. SFs, ueps weight, or object energy scale/resolution variations: JES/JER. Here's an example to add new corrections 

1. Add new info `utils/AK4_parameter.py` 
```python
"JPCalib": {
            "Run2023D-22Sep2023_v1": "calibeHistoWrite_Data2023D-22Sep2023_v1.root",
            "Run2023D-22Sep2023_v2": "calibeHistoWrite_Data2023D-22Sep2023_v2.root",
            "MC": "calibeHistoWrite_MC2023_Summer23BPix.root",
        },        
```
2. Add new collections to `load_SF` in `utils/correction.py`
Depends on corrections file type, read differently from its definition. See details in: [correctionlib official](https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/examples), or other custom way used in [coffea](https://coffea-hep.readthedocs.io/en/latest/notebooks/applying_corrections.html). This load all the correction information as `evaluator` can be use to extract weight information later
```python
for SF in config[campaign].keys():
        if SF == "DC":
            continue
        ## pileup weight
        if SF == "LUM":
            ## Check whether files in jsonpog-integration exist
            if os.path.exists(
                f"/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/LUM/{year}_{campaign}"
            ):
                correct_map["LUM"] = correctionlib.CorrectionSet.from_file(
                    f"/cvmfs/cms.cern.ch/rsync/cms-nanoAOD/jsonpog-integration/POG/LUM/{year}_{campaign}/puWeights.json.gz"
                )
            ## Otherwise custom files
            else:
                _pu_path = f"BTVNanoCommissioning.data.PU.{campaign}"
                with importlib.resources.path(
                    _pu_path, config[campaign]["LUM"]
                ) as filename:
                    if str(filename).endswith(".pkl.gz"):
                        with gzip.open(filename) as fin:
                            correct_map["LUM"] = cloudpickle.load(fin)[
                                "2017_pileupweight"
                            ]
                    elif str(filename).endswith(".json.gz"):
                        correct_map["LUM"] = correctionlib.CorrectionSet.from_file(
                            str(filename)
                        )
                    elif str(filename).endswith(".histo.root"):
                        ext = extractor()
                        ext.add_weight_sets([f"* * {filename}"])
                        ext.finalize()
                        correct_map["LUM"] = ext.make_evaluator()

```
3.1 Add weight based correction

In the `utils/correction` create the reader to get the weight 

Create your function to readout the weight from the evaluator stored in the `correction_map`, the idea is to add the weight/systematic information to the event and return to the workflow



```python
def btagSFs(jet, correct_map, weights, SFtype, syst=False):
    .....
    if i == 0 and syst == False:
            weights.add(SFtype, sfs)
       
    if syst == True:
        weights.add_multivariation(
            SFtype,
            sfs,
            systlist,
            np.array(list(sfs_up_all.values())),
            np.array(list(sfs_down_all.values())),
        )
        # in case you only have the up/down variation
        weights.add(
            SFtype,# name of the weight
            sfs,# nominal SFs
            sf_unc_up,#up varition 
            sf_unc_down, #down varition 
        )

```

In case it's a common correction, add to the `weight_manager` in `utils/correction` otherwise directly to the workflow (weight based)

```python
def weight_manager(pruned_ev, SF_map, isSyst):
    weights = Weights(len(pruned_ev), storeIndividual=True)
    ...
    btagSFs(pruned_ev.SelJet, SF_map, weights, "DeepJetC", syst_wei)
    btagSFs(pruned_ev.SelJet, SF_map, weights, "DeepJetB", syst_wei)
    btagSFs(pruned_ev.SelJet, SF_map, weights, "DeepCSVB", syst_wei)
    btagSFs(pruned_ev.SelJet, SF_map, weights, "DeepCSVC", syst_wei)
    return weights
```

3.2 Add object variations

For the object scale / resolution variation we shift object energy/positions as a list of `shifts` to the original object. 
The `shifts` is a list of shifted object after the corrctions are applied to the objects

```python
# JES uncertainty 
if "JES_Total" in jets.fields:
    shifts += [
        (
            {
                "Jet": jets.JES_Total.up, # change the objects to JES up variation
                "MET": met.JES_Total.up,
            },
            "JESUp",  # the uncertainty variation name
        ),
        (
            {
                "Jet": jets.JES_Total.down,
                "MET": met.JES_Total.down,
            },
            "JESDown",
        ),
    ]
```

In case the shifts are in common , put to `common_shifts`:
```python
if "JME" in self.SF_map.keys():
        syst_JERC = self.isSyst
        if self.isSyst == "JERC_split":
            syst_JERC = "split"
        shifts = JME_shifts(
            shifts, self.SF_map, events, self._campaign, isRealData, syst_JERC
        )
    else:
        if int(self._year) < 2020:
            shifts = [
                ({"Jet": events.Jet, "MET": events.MET, "Muon": events.Muon}, None)
            ]
        else:
            shifts = [
                (
                    {
                        "Jet": events.Jet,
                        "MET": events.PuppiMET,
                        "Muon": events.Muon,
                    },
                    None,
                )
            ]
```


otherwise in your workflow `process(self, events)` add new shifts
```python
def process(self, events):
        events = missing_branch(events)
        shifts = common_shifts(self, events)
        shifts+=[({obj:variation},name)]
```