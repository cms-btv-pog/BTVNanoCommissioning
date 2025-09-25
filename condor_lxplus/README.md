# Light Condor submission script for lxplus

The jobs submitted by this script rely on eos/afs being mounted on the condor nodes.

- It does not create a new installation of conda/mamba and instead uses a preinstalled conda env (you can replace this with your own conda/mamba path in the PATH variable, if needed).
- It does not create a tarball of the BTVNanoComm code either, simply `cd`s to the working directory in eos/afs.
- Copies the proxy locally and reads it directly from the condor node.

## Example
```
voms-proxy-init --voms cms --valid 192:00
conda activate btv_coffea   # Or conda activate /eos/home-m/milee/miniforge3/envs/btv_coffea
python condor_lxplus/submitter.py --workflow ctag_DY_sf --json fetched_list.json --campaign Summer22 --year 2022 --isArray --skipbadfiles --jobName condor_1 --outputDir /preferably/eos/output/directory --submit
```

## Check outputs
Check if all jobs succeeded with:
```
python condor_lxplus/checkoutputs.py <job_condor_dir>
```
Then manually resubmit failed jobs with the newly created config.

## hadd outputs
hadd all outputs using
```
python condor_lxplus/haddoutputs.py <job_output_dir>
```


## Pros
- This is likely faster and will run out of the box.
- Proxy handling works even if condor's native `user_proxy` method fails.

## Cons
- Relies on eos/afs mount, hence jobs will fail to run on condor nodes where the mount is unstable. **You may need to keep releasing jobs in that case.**
- Will not work where eos is not mounted, e.g. on CMSConnect nodes.
