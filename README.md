# BTVNanoCommissioning
Repository for Commissioning studies in the BTV POG based on (custom) nanoAOD samples


## Structure
Example worfkflow for ttbar is included. 

Each workflow can be a separete "processor" file, creating the mapping from NanoAOD to
the histograms we need. Workflow processors can be passed to the `runner.py` script 
along with the fileset these should run over. Multiple executors can be chosen 
(for now iterative - one by one, uproot/futures - multiprocessing and dask-slurm). 

To run the example, run:
```
python runner.py --workflow ttcom
```

Example plots can be found in ` make_some_plots.ipynb` though we might want to make
that more automatic in the end.
