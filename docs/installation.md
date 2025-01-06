## Setup the environment

You can install your [standalone conda envrionment](#standalone-conda-environment) via `yaml` or on the lxplus you can directly jump to [setup](#setup-the-framework)

### Standalone conda environment
:::{caution}
suggested to install under `bash` environment
:::

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


### Other installation options for coffea
See [https://coffeateam.github.io/coffea/installation.html](https://coffeateam.github.io/coffea/installation.html)
