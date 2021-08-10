# hgg-coffea
Tools for running the CMS Higgs to Two Photons Analysis on NanoAOD
  
[![Actions Status][actions-badge]][actions-link]
[![Documentation Status][rtd-badge]][rtd-link]
[![Code style: black][black-badge]][black-link]

[![PyPI version][pypi-version]][pypi-link]
[![Conda-Forge][conda-badge]][conda-link]
[![PyPI platforms][pypi-platforms]][pypi-link]

[![GitHub Discussion][github-discussions-badge]][github-discussions-link]
[![Gitter][gitter-badge]][gitter-link]

[actions-badge]:            https://github.com/lgray/hgg-coffea/workflows/CI/badge.svg
[actions-link]:             https://github.com/lgray/hgg-coffea/actions
[black-badge]:              https://img.shields.io/badge/code%20style-black-000000.svg
[black-link]:               https://github.com/psf/black
[conda-badge]:              https://img.shields.io/conda/vn/conda-forge/hgg-coffea
[conda-link]:               https://github.com/conda-forge/hgg-coffea-feedstock
[github-discussions-badge]: https://img.shields.io/static/v1?label=Discussions&message=Ask&color=blue&logo=github
[github-discussions-link]:  https://github.com/lgray/hgg-coffea/discussions
[gitter-badge]:             https://badges.gitter.im/https://github.com/lgray/hgg-coffea/community.svg
[gitter-link]:              https://gitter.im/https://github.com/lgray/hgg-coffea/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge
[pypi-link]:                https://pypi.org/project/hgg-coffea/
[pypi-platforms]:           https://img.shields.io/pypi/pyversions/hgg-coffea
[pypi-version]:             https://badge.fury.io/py/hgg-coffea.svg
[rtd-badge]:                https://readthedocs.org/projects/hgg-coffea/badge/?version=latest
[rtd-link]:                 https://hgg-coffea.readthedocs.io/en/latest/?badge=latest
[sk-badge]:                 https://scikit-hep.org/assets/images/Scikit--HEP-Project-blue.svg

## Structure
Example worfkflow for drell-yan studies is included. 

Each workflow can be a separate "processor" file, creating the mapping from NanoAOD to
the histograms we need. Workflow processors can be passed to the `runner.py` script 
along with the fileset these should run over. Multiple executors can be chosen 
(for now iterative - one by one, uproot/futures - multiprocessing and dask-slurm). 

To run the example, run:
```
python runner.py --workflow dystudies
```

Example plots can be found in ` make_some_plots.ipynb` though we might want to make
that more automatic in the end.

## Requirements
### Coffea installation with Miniconda
For installing Miniconda, see also https://hackmd.io/GkiNxag0TUmHnnCiqdND1Q#Local-or-remote
```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
# Run and follow instructions on screen
bash Miniconda3-latest-Linux-x86_64.sh
```
NOTE: always make sure that conda, python, and pip point to local Miniconda installation (`which conda` etc.).

You can either use the default environment`base` or create a new one:
```bash
# create new environment with python 3.7, e.g. environment of name `coffea`
conda create --name coffea python=3.7
# activate environment `coffea`
conda activate coffea
```
Install coffea, xrootd, and more:
```bash
conda install -c conda-forge coffea # pip install git+https://github.com/CoffeaTeam/coffea.git # for bleeding edge
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
```bash
jupyter notebook --ip=127.0.0.1 --port 8800 --no-browser
```
4. URL for notebook will be printed, copy and open in local browser

## Scale-out (Sites)

Scale out can be notoriously tricky between different sites. Coffea's integration of `slurm` and `dask`
makes this quite a bit easier and for some sites the ``native'' implementation is sufficient, e.g Condor@DESY.
However, some sites have certain restrictions for various reasons, in particular Condor @CERN and @FNAL.

### Condor@FNAL (CMSLPC)
Follow setup instructions at https://github.com/CoffeaTeam/lpcjobqueue, run them from within the `hgg-coffea` directory that you have checked out.
After starting the singularity container run a test with 
```bash
python runner.py --meta Era2017_legacy_v1.json --wf dystudies -d root://cmseos.fnal.gov//store/user/$USER/hgg_test/ --executor dask/lpc --samples filefetcher/dystudies.json --chunk=100000 --max=5
```

### Condor@CERN (lxplus)
Only one port is available per node, so its possible one has to try different nodes until hitting
one with `8786` being open. Other than that, no additional configurations should be necessary.

```bash
python runner.py --wf dystudies --executor dask/lxplus
```

### Coffea-casa (Nebraska AF)
Coffea-casa is a JupyterHub based analysis-facility hosted at Nebraska. For more information and setup instuctions see
https://coffea-casa.readthedocs.io/en/latest/cc_user.html

After setting up and checking out this repository (either via the online terminal or git widget utility run with
```bash
python runner.py --wf dystudies --executor dask/casa
```
Authentication is handled automatically via login auth token instead of a proxy. File paths need to replace xrootd redirector with "xcache", `runner.py` does this automatically.

### Canned example on cms-lpc scaleout
```bash
python runner.py --meta Era2017_legacy_v1.json --wf dystudies -d root://cmseos.fnal.gov//store/user/$USER/hgg_test/ --executor dask/lpc --samples filefetcher/dystudies.json --chunk=100000 --scaleout 1 --limit 2 --only DYJets-M50 --ts DummyTagger1 DummyTagger2
```
