"""Merge a condor output directory one dataset at a time, resumably.

``combinecoffeaoutputs.py`` folds every job file into one accumulator, which peaks
at roughly three times the size of the merged result -- ~20 GB for a full HT-binned
MC campaign at the current histogram granularity, enough to be OOM-killed on a
shared lxplus node.  Each job file holds exactly one dataset (``split_samples.json``
splits within a dataset), so folding per dataset and then taking the union of the
per-dataset results gives the identical output with a peak of one dataset instead
of the whole campaign (~3 GB).

Each finished dataset is written to ``<outdir>/_mergeparts/<dataset>.coffea`` before
the next one starts, so a run killed halfway (lxplus reaps processes when the login
session ends) resumes where it stopped instead of starting over.

Usage::

    python condor_lxplus/combine_by_dataset.py <output_dir> [--job-dir DIR] [--keep-parts]

The job directory defaults to ``jobs_<output_dir>``; every ``jobs_<output_dir>*``
directory is read, so job numbers from a ``*_retry`` submission are picked up too.
"""

import argparse
import json
import os
import re
import time
from glob import glob

PARTS = "_mergeparts"


def job_dataset_map(job_dirs: list[str]) -> dict[str, str]:
    """``{job number: dataset}`` from every split_samples.json."""
    jobmap: dict[str, str] = {}
    for d in job_dirs:
        path = os.path.join(d, "split_samples.json")
        if not os.path.isfile(path):
            continue
        for jid, samp in json.load(open(path)).items():
            names = list(samp)
            if len(names) != 1:
                raise ValueError(f"job {jid} in {path} spans {len(names)} datasets")
            jobmap[jid] = names[0]
    if not jobmap:
        raise FileNotFoundError(f"no split_samples.json in {job_dirs}")
    return jobmap


def files_by_dataset(outdir: str, jobmap: dict[str, str]) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    for d in sorted(os.listdir(outdir)):
        if not re.fullmatch(r"hists_\d+", d):
            continue
        f = os.path.join(outdir, d, f"{d}.coffea")
        if not os.path.isfile(f):
            print(f"Skipping {d}, no output inside.")
            continue
        jid = d.split("_")[1]
        if jid not in jobmap:
            raise KeyError(f"{d} is not in any split_samples.json; wrong --job-dir?")
        out.setdefault(jobmap[jid], []).append(f)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("outdir")
    ap.add_argument("--job-dir", default=None, help="default: jobs_<outdir>*")
    ap.add_argument(
        "--keep-parts", action="store_true", help="keep _mergeparts/ afterwards"
    )
    args = ap.parse_args()

    from coffea.processor import accumulate
    from coffea.util import load, save

    outdir = args.outdir.rstrip("/")
    job_dirs = (
        [args.job_dir]
        if args.job_dir
        else sorted(glob(f"jobs_{os.path.basename(outdir)}*"))
    )
    by_ds = files_by_dataset(outdir, job_dataset_map(job_dirs))
    parts = os.path.join(outdir, PARTS)
    os.makedirs(parts, exist_ok=True)
    print(
        f"{sum(map(len, by_ds.values()))} job files in {len(by_ds)} dataset(s), from {job_dirs}"
    )

    for ds, files in sorted(by_ds.items()):
        part = os.path.join(parts, f"{ds}.coffea")
        if os.path.isfile(part):
            print(f"  {ds}: part exists, skipping")
            continue
        t0 = time.time()
        acc = load(files[0])
        for f in files[1:]:
            acc = accumulate([acc, load(f)])
        if list(acc) != [ds]:
            raise ValueError(f"{ds}: part holds {list(acc)}")
        save(acc, part + ".tmp")
        os.replace(part + ".tmp", part)  # only a complete part ever gets the final name
        print(f"  {ds}: {len(files)} files -> part in {time.time() - t0:.0f}s")

    combined = {}
    for ds in sorted(by_ds):
        combined.update(load(os.path.join(parts, f"{ds}.coffea")))  # keys are disjoint
    out = os.path.join(outdir, "combined_hist.coffea")
    save(combined, out + ".tmp")
    os.replace(out + ".tmp", out)
    print(f"Wrote {out}: {len(combined)} datasets, {os.path.getsize(out) / 1e6:.0f} MB")

    if not args.keep_parts:
        for ds in by_ds:
            os.remove(os.path.join(parts, f"{ds}.coffea"))
        os.rmdir(parts)


if __name__ == "__main__":
    main()
