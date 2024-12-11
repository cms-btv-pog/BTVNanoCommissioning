# Example coffea file 

Quick tutorial to go through coffea files from commissioning framework.

## Structure of the file

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

The histogram is a multidimentinal histogram, with all the axis listed
```python
Hist(
  IntCategory([0, 1, 4, 5, 6], name='flav', label='Genflavour'),# different genflavor, 0 for light, 1 for PU, 2 for c, 3 for b. Always 0 for data.
  IntCategory([1, -1], name='osss', label='OS(+)/SS(-)'),# opposite sign or same sign, only appears in W+c workflow
  StrCategory(['noSF','PUUp','PUDown'], growth=True, name='syst'),# systematics variations,
  Regular(30, -0.2, 1, name='discr', label='btagDeepFlavB_lepb'),# discriminator distribution, the last axis is always the variable
  storage=Weight()) # Sum: WeightedSum(value=140, variance=140)# Value is sum of the entries, Variances is sum of the variances.
```

## Read coffea files and explore the histogram

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
hist1d=hist[axis]
# you can also sum over the axis, e.g. here shows no jet flavor split and sum os+ss
axis={'flav':sum,'os':sum,'syst':'noSF'}
# rebin the axis is also possible, rebin the discrimnator by merged two bins into one
axis={'flav':sum,'os':sum,'syst':'noSF','discr':hist.rebin(2)}
```
