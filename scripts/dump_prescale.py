import json, os, argparse, pandas, io, subprocess, numpy
import correctionlib.schemav2 as cs
from BTVNanoCommissioning.helpers.BTA_helper import BTA_HLT
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
import numpy as np, pandas as pd

### NOTICE The scripts only works on lxplus...

parser = argparse.ArgumentParser(description="Create prescale weights(lxplus)")

parser.add_argument(
    "-l",
    "--lumimask",
    default="src/BTVNanoCommissioning/data/DC/Cert_Collisions2022_355100_362760_Golden.json",
    help="lumimask to generate prescale weights",
)
parser.add_argument(
    "-H",
    "--HLT",
    default=BTA_HLT,
    type=str,
    help="Which HLT is used; comma separated for multiple.",
)
parser.add_argument("-v", "--verbose", action="store_true", help="debugging")
parser.add_argument("-t", "--test", action="store_true", help="test with only 5 runs")
parser.add_argument("-f", "--force", action="store_true", help="recreate .csv")
parser.add_argument(
    "-i", "--ignore_csv_output", action="store_true", help="Ignore writing the .csv"
)
parser.add_argument(
    "-n",
    "--nthreads",
    default=None,
    type=int,
    help="Number of threads for parallel run processing (default: auto)",
)

### NOTICE The scripts only works on lxplus...


def process_run(run_input):
    run, trg = run_input

    bc_alias = "singularity -s exec  --env PYTHONPATH=/home/bril/.local/lib/python3.10/site-packages /cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-cloud/brilws-docker:latest brilcalc"
    brilcall = [
        f"{bc_alias} trg -r {run} --prescale --hltpath HLT_{trg}_v* --output-style csv"
    ]

    command = subprocess.run(
        brilcall,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        shell=True,
        executable="/bin/bash",
    )

    if command.returncode != 0:
        print(f"Error: {command.stderr}")
        print(
            "Check that you have sourced brilcalc with the command `source /cvmfs/cms-bril.cern.ch/cms-lumi-pog/brilws-docker/brilws-env`"
        )
        return None

    csv_output = command.stdout
    df = pandas.read_csv(
        io.StringIO(csv_output),
        usecols=["cmsls", "totprescval", "# run", "hltpath/prescval"],
        dtype={
            "cmsls": numpy.int32,
            "totprescval": numpy.float64,
            "# run": numpy.int32,
            "hltpath/prescval": str,
        },
    )

    return df


def get_prescale(
    HLT,
    lumimask,
    verbose=False,
    test=False,
    force=False,
    ignore_csv_output=False,
    nthreads=None,
):
    prescales = pandas.DataFrame()
    runs = json.load(open(lumimask))
    runs = list(runs.keys())
    if test:
        runs = runs[: min(len(runs), 10)]

    outcsv = f"src/BTVNanoCommissioning/data/Prescales/HLTinfo_{HLT}_run{min(runs)}_{max(runs)}.csv"
    if force or not os.path.exists(outcsv):
        with ThreadPoolExecutor(max_workers=nthreads) as executor:
            dfs = list(
                tqdm(
                    executor.map(process_run, [(run, HLT) for run in runs]),
                    total=len(runs),
                )
            )

        prescales = pandas.concat(dfs, ignore_index=True)
        if not ignore_csv_output:
            prescales.to_csv(outcsv)

        if verbose:
            print("prescales :", prescales)
    else:
        prescales = pandas.read_csv(outcsv)
    return prescales


### code from Hsin-Wei: https://github.com/cms-btv-pog/BTVNanoCommissioning/blob/f2a5db0e325c9b26d220089a49ddb8f73682f846/prescales.ipynb
## read prescale


def get_ps(ps, verbose=False):
    if len(ps) != 1:
        print(ps)
        print("Length of ps after selection ", len(ps))
        raise ValueError(ps)
    if verbose:
        print("Final prescale weight: ", ps.iloc[0]["totprescval"])
    return float(ps.iloc[0]["totprescval"])


def build_lumibins(ps, verbose=False):
    ##### to sort as bin edges properly, starting lumi sections need to be stored as floats
    if verbose:
        print("Path: ", ps["hltpath/prescval"], ps["totprescval"])
    edges = sorted(set(ps["cmsls"]))
    if len(edges) == 1:
        return get_ps(ps)
    elif len(edges) > 1:
        edges.append("inf")
        if verbose:
            print("Lumi bin edges: ", list(zip(edges[:-1], edges[1:])))
        content = [
            get_ps(
                ps[
                    (ps["cmsls"].astype(float) >= lo)
                    & (ps["cmsls"].astype(float) < float(hi))
                ]
            )
            for lo, hi in zip(edges[:-1], edges[1:])
        ]
        if verbose:
            print("Prescales: ", content)
        return cs.Binning.model_validate(
            {
                "nodetype": "binning",
                "input": "lumi",
                "edges": edges,
                "content": content,
                "flow": "clamp",
            }
        )


def build_runs(ps, all_runs, HLT_paths, verbose=False):
    df = ps
    runs_np = np.array(sorted(df["# run"].unique())).astype(int)
    for run in all_runs:
        if run in runs_np:
            continue
        index = np.abs(runs_np - int(run)).argmin()
        closest_run = runs_np[index]
        new_row = df.loc[df["# run"] == closest_run].iloc[0].copy()
        new_row["# run"] = int(run)
        df = pd.concat([df, new_row.to_frame().T], ignore_index=True)
        print(
            f"[WARNING] Could not find info for Run #{run}! Replacing with info from closest run, #{closest_run}."
        )

    runs = sorted(df["# run"].unique())
    if verbose:
        print("Selected ", len(runs), ": ", runs)
    return cs.Category.model_validate(
        {
            "nodetype": "category",
            "input": "run",
            "content": [
                {
                    "key": int(run),
                    "value": build_paths(df[df["# run"] == run], HLT_paths, verbose),
                }
                for run in runs
            ],
        }
    )


def build_paths(ps, HLT_paths, verbose=False):
    if verbose:
        print("Run: ", ps["# run"].iloc[0], type(ps["# run"].iloc[0]))
    # paths are unique bc of hltpath/lumi --> make array of path name separate
    paths = [HLT_paths]
    if verbose:
        print("Type of path key: ", type(paths[0]), paths)
    return cs.Category.model_validate(
        {
            "nodetype": "category",
            "input": "path",
            "content": [
                {"key": str(path), "value": build_lumibins(ps, verbose)}
                for path in paths
            ],
        }
    )


if __name__ == "__main__":
    args = parser.parse_args()
    if "," in args.HLT:
        args.HLT = args.HLT.split(",")
    else:
        args.HLT = [args.HLT]

    os.system(
        "source /cvmfs/cms-bril.cern.ch/cms-lumi-pog/brilws-docker/brilws-env; which brilcalc"
    )

    # Get all_runs from the lumimask
    with open(args.lumimask) as f:
        lumimask_data = json.load(f)
    runs_in_mask = set([int(run) for run in lumimask_data.keys()])
    print(f"Total runs in lumimask: {len(runs_in_mask)}")

    for HLT in args.HLT:
        if HLT.startswith("HLT_"):
            HLT = HLT[4:]
        print("HLT : ", HLT)
        ps_csvData = get_prescale(
            HLT,
            args.lumimask,
            args.verbose,
            args.test,
            args.force,
            args.ignore_csv_output,
            args.nthreads,
        )
        psCorr = cs.Correction.model_validate(
            {
                "version": 2,
                "name": "prescaleWeight",
                "inputs": [
                    {"name": "run", "type": "int"},
                    {"name": "path", "type": "string"},
                    {"name": "lumi", "type": "real"},
                ],
                "output": {"name": "weight", "type": "real"},
                "data": build_runs(
                    ps_csvData, runs_in_mask, "HLT_" + HLT, args.verbose
                ),
            }
        )
        cset = cs.CorrectionSet(
            schema_version=2,
            corrections=[psCorr],
            description=f"prescales for HLT_{HLT}",
        )
        dumpfile = f"src/BTVNanoCommissioning/data/Prescales/ps_weight_{HLT}_run{min(runs_in_mask)}_{max(runs_in_mask)}.json"
        with open(
            dumpfile,
            "w",
        ) as f:
            f.write(cset.model_dump_json(exclude_unset=True))
