
# BTVNanoCommissioning
[![Linting](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/python_linting.yml/badge.svg)](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/python_linting.yml)
[![btag ttbar](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/ttbar_workflow.yml/badge.svg)](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/ttbar_workflow.yml)
[![ctag ttbar](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/ctag_ttbar_workflow.yml/badge.svg)](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/ctag_ttbar_workflow.yml)
[![ctag DY+jets Workflow](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/ctag_DY_workflow.yml/badge.svg)](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/ctag_DY_workflow.yml)
[![ctag W+c Workflow](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/ctag_Wc_workflow.yml/badge.svg)](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/ctag_Wc_workflow.yml)
[![BTA Workflow](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/BTA_workflow.yml/badge.svg)](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/BTA_workflow.yml)
[![QCD Workflow](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/QCD_workflow.yml/badge.svg)](https://github.com/cms-btv-pog/BTVNanoCommissioning/actions/workflows/QCD_workflow.yml)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

Repository for Commissioning studies in the BTV POG based on (custom) nanoAOD samples


This framework is based on [coffea](https://coffeateam.github.io/coffea/) and using [btvnano](https://btv-wiki.docs.cern.ch/SoftwareAlgorithms/PFNano/) as input. The framework is also used as frontend for the btv automation task [autobtv](https://gitlab.cern.ch/cms-analysis/btv/software-and-algorithms/autobtv)

This framework is based on [coffea processor](https://coffeateam.github.io/coffea/concepts.html#coffea-processor). Each workflow can be a separate **processor** file in the `workflows`, creating the mapping from `PFNano` to the histograms as `coffea` file or creating `.root` files by saving awkward arrays. Workflow processors can be passed to the `runner.py` script along with the fileset these should run over. Multiple executors can be chosen
(`iterative` - one by one, `futures` - multiprocessing). Scale out to clusters depend on facilities. Obtain the histograms as plot(`.pdf`) or save to template `.root` file with dedicated scripts

The minimum requirement commands are shown in follow, specified the selections, datataset, campaign and year
```
python runner.py --workflow ttsemilep_sf --json metadata/test_bta_run3.json --campaign Summer22EERun3 --year 2022
```
- Detailed documentation [here](https://btvnanocommissioning.readthedocs.io/en/latest/)
- To running the commissioning task or producing the template: go to [Preparation for commissioning/SFs tasks](https://btvnanocommissioning.readthedocs.io/en/latest/user.html)
- To develop new workflow, the instruction can be found in [Instruction for developers](https://btvnanocommissioning.readthedocs.io/en/latest/user.html)
- Current working in progress [issues](https://gitlab.cern.ch/cms-btv-coordination/tasks/-/issues/?label_name%5B%5D=Software%3A%3A%20BTVnano%20%26CommFW)



## Setup 

You can install your [standalone conda envrionment](#standalone-conda-environment) via `yaml` or on the lxplus you can directly jump to [setup](#setup-the-framework)
### Standalone conda environment
> [!Caution]
> suggested to install under `bash` environment


For installing Micromamba, see [[here](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html)]
```
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
# Run and follow instructions on screen
bash Miniforge3-$(uname)-$(uname -m).sh

micromamba activate
```
NOTE: always make sure that conda, python, and pip point to local micromamba installation (`which conda` etc.).


You can simply create the environment through the existing `test_env.yml` under your micromamba environment using micromamba, and activate it
```
micromamba env create -f test_env.yml 

```
### Setup the framework

```bash
# activate enviroment once you have coffea framework 
conda/micromamba activate btv_coffea

conda/micromamba activate  /eos/home-m/milee/miniforge3/envs/btv_coffea

# only first time, including submodules
git clone git@github.com:cms-btv-pog/BTVNanoCommissioning.git 
# Once the environment is set up, compile the python package:
pip install -e .
pip install -e .[dev, doc] # for developer
```

You can still install additional packages itself by `pip install $PACKAGE`

`conda/micromamba activate btv_coffea` is required to setup

