## Scale-out (Sites & file split scheme) 

Scale out can be notoriously tricky between different sites. Coffea's integration of `slurm` and `dask`
makes this quite a bit easier and for some sites the ``native'' implementation is sufficient, e.g Condor@DESY.
However, some sites have certain restrictions for various reasons, in particular Condor @CERN and @FNAL. The scaleout scheme is named as follows: `$cluster_schedule_system/scheduler/site`. The existing sites are documented in [sites configuration](#sites-configuration-with-daskparsl-schedular) while [standalone condor submission](#standalone-condor-jobslxpluscmsconnect) is possible and strongly suggested when working on lxplus. 


::: {admonition} file splitting scheme
   Here we also provide some file split schemes provided in the framework, the later two schemes can be used together to boost the processing time. **All in one** case is the default any non-lxplus clusters which avoid complications for merging several files in plotting step. **Dataset split** scheme can be used to save intermediate steps while the input json file need to be created by users. **File split** scheme is the default scheme for lxplus cluster(see details in [condor-dask@lxplus](#condorcern-lxplus)) where the file list of each dataset in the json file are split by a certain numbers of file `--fsize` and process sequencially automatically. You can also combine both **dataset split** and **file split** to parallize the job into different processes(notice in lxplus you need to login different machine due to port restriction) with quick job.

    

Memory usage is also useful to adapt to cluster. Check the memory by calling  `memory_usage_psutil()` from `helpers.func.memory_usage_psutil` to optimize job size. Example with `ectag_Wc_sf` summarized below.

| Type        |Array+Hist |  Hist only| Array Only|
| :---:   | :---: | :---: | :---: |
DoubleMuon (BTA,BTV_Comm_v2)| 1243MB |	848MB	|1249MB|
DoubleMuon (PFCands, BTV_Comm_v1)|1650MB	|1274MB	|1632MB|
DoubleMuon (Nano_v11)|1183MB|	630MB	|1180MB|
WJets_inc (BTA,BTV_Comm_v2)| 1243MB	|848MB	|1249MB|
WJets_inc (PFCands, BTV_Comm_v1)|1650MB	|1274MB	|1632MB
WJets_inc (Nano_v11)|1183MB	|630MB	|1180MB|

:::

### Sites configuration with dask/parsl schedular



#### dask: Condor@FNAL (CMSLPC)
Follow setup instructions at https://github.com/CoffeaTeam/lpcjobqueue. After starting 
the singularity container run with 
```bash
python runner.py --wf ttcom --executor dask/lpc
```

#### daskLCondor@CERN (lxplus)
Only one port is available per node, so its possible one has to try different nodes until hitting
one with `8786` being open. Other than that, no additional configurations should be necessary.

```bash
python runner.py --wf ttcom --executor dask/lxplus
```

jobs automatically split to 50 files per jobs to avoid job failure due to crowded cluster on lxplus with the naming scheme `hist_$workflow_$json_$dictindex_$fileindex.coffea`. The `.coffea` files can be then combined at plotting level


:::{caution}
The optimal scaleout options on lxplus are `-s 50 --chunk 50000`
:::

To deal with unstable condor cluster and dask worker on lxplus, you can resubmit failure jobs via `--index $dictindex,$fileindex` option. `$dictindex` refers to the index in the `.json dict`, `$fileindex` refers to the index of the file list split to 50 files per dask-worker. The total number of files of each dict can be computed by `math.ceil(len($filelist)/50)` The job will start from the corresponding indices.

#### Coffea-casa (Nebraska AF)
Coffea-casa is a JupyterHub based analysis-facility hosted at Nebraska. For more information and setup instuctions see
https://coffea-casa.readthedocs.io/en/latest/cc_user.html

After setting up and checking out this repository (either via the online terminal or git widget utility run with
```bash
python runner.py --wf ttcom --executor dask/casa
```
Authentication is handled automatically via login auth token instead of a proxy. File paths need to replace xrootd redirector with "xcache", `runner.py` does this automatically.


#### parsl/dask with Condor 
```bash
python runner.py --wf ttcom --executor dask/condor(parsl/condor)
```

:::{tip}
if the jobs hang for a long time, you can check the logs of the condor jobs `runinfo/$ID/submit_scripts/parsl.parsl.run.blockxxx.err`
:::


### Standalone condor jobs@lxplus/cmsconnect

:::{caution}
Strongly suggest to use this on lxplus.
Check the end of this subsection for an alternative, lighter lxplus condor submitter with some caveats.
:::

You have the option to run the framework through "standalone condor jobs", bypassing the native coffea-supported job submission system. Within each job you submit, a standalone script will execute the following on the worker node:

 - Set up a minimal required Python environment.
 - Retrieve the BTVNanoCommissioning repository, either from a git link or transferred locally.
 - Launch the `python runner.py ...` command to execute the coffea framework in the iterative executor mode.
 
This utility is currently adapted for the lxplus and cmsconnect condor systems. To generate jobs for launching, replace `python runner.py` with `python condor/submitter.py`, append the existing arguments, and add the following arguments in addition:

 - `--jobName`: Specify the desired condor job name. A dedicated folder will be generated, including all submission-related files.
 - `--outputXrootdDir`: Indicate the XRootD directory's path (starting with `root://`) where the produced .coffea (and .root) files from each worker node will be transferred to.
 - `--condorFileSize`: Define the number of files to process per condor job (default is 50). The input file list will be divided based on this count.
 - `--remoteRepo` (optional, but recommended): Specify the path and branch of the remote repository to download the BTVNanoCommissioning repository. If not specified, the local directory will be packed and transferred as the condor input, potentially leading to higher loads for condor transfers. Use the format e.g. `--remoteRepo 'https://github.com/cms-btv-pog/BTVNanoCommissioning.git -b master'`.

After executing the command, a new folder will be created, preparing the submission. Follow the on-screen instructions and utilize `condor_submit ...` to submit the jdl file. The output will be transferred to the designated XRootD destination.

::: {admonition} Frequent issues for standalone condor jobs submission

1. CMS Connect provides a condor interface where one can submit jobs to all resources available in the CMS Global Pool. See [WorkBookCMSConnect Twiki](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCMSConnect#Requesting_different_Operative_S) for the instructions if you use it for the first time.
2. The submitted jobs are of the kind which requires a proper setup of the X509 proxy, to use the XRootD service to access and store data. In the generated `.jdl` file, you may see a line configured for this purpose `use_x509userproxy = true`. If you have not submitted jobs of this kind on lxplus condor, we recommend you to add a line
   ```bash
   export X509_USER_PROXY=$HOME/x509up_u`id -u`
   ```
   to `.bashrc` and run it so the proxy file will be stored in your AFS folder instead of in your `/tmp/USERNAME` folder. For submission on cmsconnect, no specific action is required.
:::

<details>
  <summary>A lighter version that works *only* on lxplus</summary>

   The jobs submitted by this script rely on eos/afs being mounted on the condor nodes.

   - It does not create a new installation of conda/mamba and instead uses a preinstalled conda env (you can replace this with your own conda/mamba path in the PATH variable, if needed).
   - It does not create a tarball of the BTVNanoComm code either, simply `cd`s to the working directory in eos/afs.
   - Copies the proxy locally and reads it directly from the condor node.

   **Example**
   ```
   voms-proxy-init --voms cms --valid 192:00
   conda activate btv_coffea   # Or conda activate /eos/home-m/milee/miniforge3/envs/btv_coffea
   python condor_lxplus/submitter.py --workflow ctag_DY_sf --json fetched_list.json --campaign Summer22 --year 2022 --isArray --skipbadfiles --jobName condor_1 --outputDir /preferably/eos/output/directory --submit
   ```

   **Check outputs**
   Check if all jobs succeeded with:
   ```
   python condor_lxplus/checkoutputs.py <job_condor_dir>
   ```
   Then manually resubmit failed jobs with the newly created config.

   **hadd outputs**
   hadd all outputs using
   ```
   python condor_lxplus/haddoutputs.py <job_output_dir>
   ```

   **Pros**
   - This is likely faster and will run out of the box.
   - Proxy handling works even if condor's native `user_proxy` method fails.

   **Cons**
   - Relies on eos/afs mount, hence jobs will fail to run on condor nodes where the mount is unstable. **You may need to keep releasing jobs in that case.**
   - Will not work where eos is not mounted, e.g. on CMSConnect nodes.

</details>

### FAQ for submission

- All jobs held: might indicate environment setup issue→ check out the condor err/out for parsl jobs the info are in `runinfo/JOBID/submit_scripts/`
- Exit without complain: might be huge memory consumption: 
   - Reduce `--chunk`, especially JERC variation are memory intense
   - check the memory usage by calling `memory_usage_psutil`
- partially failed/held: 	
   - could be temporarily unavailable of the files/site. If the retries not work, considering obtained failure file list and resubmit.
   - error of certain files→ check the failed files and run it locally with `--executor iterative`
