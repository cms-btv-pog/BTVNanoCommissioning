import json, os, argparse, pandas
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
parser.add_argument("-H", "--HLT", default=None, type=str, help="Which HLT is used; comma separated for multiple; do not add the 'HLT_' prefix")
parser.add_argument("-v", "--verbose", action="store_true", help="debugging")
parser.add_argument("-t", "--test", action="store_true", help="test with only 5 runs")
parser.add_argument("-f", "--force", action="store_true", help="recreate .csv")

### NOTICE The scripts only works on lxplus...


def process_run(ir_run):
    ir, run = ir_run
    tmpfile = f".tmp_{ir}.csv"
    os.system(
        f"singularity -s exec --env PYTHONPATH=/home/bril/.local/lib/python3.10/site-packages "
        f"/cvmfs/unpacked.cern.ch/gitlab-registry.cern.ch/cms-cloud/brilws-docker:latest "
        f"brilcalc trg --prescale --hltpath 'HLT_{HLT}_v*' -r {run} --output-style csv &>{tmpfile}"
    )
    return pandas.read_csv(tmpfile)


def get_prescale(HLT, lumimask, verbose=False, test=False, force=False):
    # os.system("source /cvmfs/cms-bril.cern.ch/cms-lumi-pog/brilws-docker/brilws-env")
    prescales = pandas.DataFrame()
    runs = json.load(open(lumimask))
    runs = list(runs.keys())
    if test:
        runs = runs[: min(len(runs), 100)]

    outcsv = f"src/BTVNanoCommissioning/data/Prescales/HLTinfo_{HLT}_run{runs[0]}_{runs[-1]}.csv"
    if force or not os.path.exists(outcsv):
        with ThreadPoolExecutor() as executor:
            dfs = list(
                tqdm(executor.map(process_run, enumerate(runs)), total=len(runs))
            )

        prescales = pandas.concat(dfs, ignore_index=True)
        # prescales= prescales[prescales['totprescval']!=0]
        prescales.to_csv(outcsv)
        os.system(f"rm -rf .tmp_*.csv")

        if verbose:
            print("prescales :", prescales)
    else:
        prescales = pandas.read_csv(outcsv)
    return runs, prescales


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
    edges = sorted(set(ps["cmsls"].astype(float)))
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


def build_runs(ps, allruns, HLT_paths, verbose=False):
    runs = sorted(ps["# run"].unique())
    runs_np = np.array(runs).astype(int)
    for run in allruns:
        if int(run) in runs_np:
            continue
        index = np.abs(runs_np - int(run)).argmin()
        closest_run= runs_np[index]
        new_row = ps.loc[ps['# run'] == closest_run].iloc[0].copy()
        new_row['# run'] = int(run)
        ps = pd.concat([ps, new_row.to_frame().T], ignore_index=True)
        print(f"[WARNING] Could not find info for Run #{run}! Replacing with info from closest run, #{closest_run}.")

    runs = sorted(ps["# run"].unique())
    if verbose:
        print("Selected ", len(runs), ": ", runs)
    return cs.Category.model_validate(
        {
            "nodetype": "category",
            "input": "run",
            "content": [
                {
                    "key": int(run),
                    "value": build_paths(ps[ps["# run"] == run], HLT_paths, verbose),
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
    if args.HLT is None:
        args.HLT = BTA_HLT
    else:
        if "," in args.HLT:
            args.HLT = args.HLT.split(",")
        else:
            args.HLT = [args.HLT]

    os.system(
        "source /cvmfs/cms-bril.cern.ch/cms-lumi-pog/brilws-docker/brilws-env; which brilcalc"
    )

    for HLT in args.HLT:
        print("HLT : ", HLT)
        runs, ps_csvData = get_prescale(
            HLT, args.lumimask, args.verbose, args.test, args.force
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
                "data": build_runs(ps_csvData, runs, "HLT_" + HLT, args.verbose),
            }
        )
        cset = cs.CorrectionSet(
            schema_version=2,
            corrections=[psCorr],
            description=f"prescales for HLT_{HLT}",
        )
        runs = json.load(open(args.lumimask))

        runs = list(runs.keys())
        dumpfile = f"src/BTVNanoCommissioning/data/Prescales/ps_weight_{HLT}_run{runs[0]}_{runs[-1]}.json"
        with open(
            dumpfile,
            "w",
        ) as f:
            f.write(cset.model_dump_json(exclude_unset=True))

        print(f"Dumped prescales in {dumpfile}.")
