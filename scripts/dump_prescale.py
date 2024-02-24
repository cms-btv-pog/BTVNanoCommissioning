import json, os, argparse, pandas
import correctionlib.schemav2 as cs
from BTVNanoCommissioning.helpers.BTA_helper import BTA_HLT

### NOTICE The scripts only works on lxplus...

parser = argparse.ArgumentParser(description="Create prescale weights(lxplus)")

parser.add_argument(
    "--lumimask",
    default="src/BTVNanoCommissioning/data/lumiMasks/Cert_Collisions2022_355100_362760_Golden.json",
    help="lumimask to generate prescale weights",
)
parser.add_argument("--HLT", default=None, type=str, help="Which HLT is used")
parser.add_argument("-v", "--verbose", action="store_true", help="debugging")


### NOTICE The scripts only works on lxplus...
def get_prescale(HLT, lumimask, verbose=False):
    # os.system("source /cvmfs/cms-bril.cern.ch/cms-lumi-pog/brilws-docker/brilws-env")
    prescales = pandas.DataFrame()
    runs = json.load(open(lumimask))
    runs = list(runs.keys())
    os.system("source /cvmfs/cms-bril.cern.ch/cms-lumi-pog/brilws-docker/brilws-env")
    if not os.path.exists(
        f"src/BTVNanoCommissioning/data/Prescales/HLTinfo_{HLT}_run{runs[0]}_{runs[-1]}.csv"
    ):

        for run in runs:

            os.system(
                f"brilcalc trg --prescale --hltpath 'HLT_{HLT}*'  -r {run}  --output-style csv &>tmp.csv"
            )
            ps = pandas.read_csv("tmp.csv")
            prescales = prescales._append(ps)
        # prescales= prescales[prescales['totprescval']!=0]
        prescales.to_csv(
            f"src/BTVNanoCommissioning/data/Prescales/HLTinfo_{HLT}_run{runs[0]}_{runs[-1]}.csv"
        )

        if verbose:
            print("prescaled :", prescales)
    else:
        prescales = pandas.read_csv(
            f"src/BTVNanoCommissioning/data/Prescales/HLTinfo_{HLT}_run{runs[0]}_{runs[-1]}.csv"
        )
    return prescales


### code from Hsin-Wei: https://github.com/cms-btv-pog/BTVNanoCommissioning/blob/f2a5db0e325c9b26d220089a49ddb8f73682f846/prescales.ipynb
## read prescale


def get_ps(ps, verbose=False):
    if len(ps) != 1:
        if verbose:
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
        edges.append(float("inf"))
        if verbose:
            print("Lumi bin edges: ", list(zip(edges[:-1], edges[1:])))
        content = [
            get_ps(
                ps[(ps["cmsls"].astype(float) >= lo) & (ps["cmsls"].astype(float) < hi)]
            )
            for lo, hi in zip(edges[:-1], edges[1:])
        ]
        if verbose:
            print("Prescales: ", content)
        return cs.Binning.parse_obj(
            {
                "nodetype": "binning",
                "input": "lumi",
                "edges": edges,
                "content": content,
                "flow": "clamp",
            }
        )


def build_runs(ps, HLT_paths, verbose=False):
    runs = sorted(ps["# run"].unique())
    if verbose:
        print("Selected ", len(runs), ": ", runs)
    return cs.Category.parse_obj(
        {
            "nodetype": "category",
            "input": "run",
            "content": [
                {
                    "key": int(run),
                    "value": build_paths(ps[ps["# run"] == run], HLT_paths),
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
    return cs.Category.parse_obj(
        {
            "nodetype": "category",
            "input": "path",
            "content": [
                {"key": str(path), "value": build_lumibins(ps)} for path in paths
            ],
        }
    )


if __name__ == "__main__":
    args = parser.parse_args()
    if args.HLT is None:
        args.HLT = BTA_HLT
    else:
        if args.HLT.split(","):
            args.HLT = args.HLT
        args.HLT = [args.HLT]

    for HLT in args.HLT:
        print("HLT : ", HLT)
        ps_csvData = get_prescale(HLT, args.lumimask, args.verbose)
        psCorr = cs.Correction.parse_obj(
            {
                "version": 2,
                "name": "prescaleWeight",
                "inputs": [
                    {"name": "run", "type": "int"},
                    {"name": "path", "type": "string"},
                    {"name": "lumi", "type": "real"},
                ],
                "output": {"name": "weight", "type": "real"},
                "data": build_runs(ps_csvData, "HLT_" + HLT),
            }
        )
        cset = cs.CorrectionSet(
            schema_version=2,
            corrections=[psCorr],
            description=f"prescales for HLT_{HLT}",
        )
        runs = json.load(open(args.lumimask))

        runs = list(runs.keys())
        with open(
            f"src/BTVNanoCommissioning/data/Prescales/ps_weight_{HLT}_run{runs[0]}_{runs[-1]}.json",
            "w",
        ) as f:
            f.write(cset.json(exclude_unset=True))
