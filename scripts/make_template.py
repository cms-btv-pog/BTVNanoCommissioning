from coffea.util import load
import numpy as np
import uproot, sys, os, argparse, hist
from BTVNanoCommissioning.helpers.xs_scaler import collate, scaleSumW
from BTVNanoCommissioning.utils.plot_utils import (
    rebin_hist,
)
import numpy as np

parser = argparse.ArgumentParser(description="Make templates from coffea files")
parser.add_argument(
    "-i",
    "--input",
    type=str,
    required=True,
    help="Input coffea file(s)",
)
parser.add_argument(
    "-v",
    "--variable",
    type=str,
    required=True,
    help="Variables to store(histogram name)",
)
parser.add_argument(
    "-a",
    "--axis",
    type=str,
    required=True,
    help="dict, put the slicing of histogram, specify 'sum' option as string",
)
parser.add_argument(
    "--lumi",
    type=float,
    required=True,
    help="Luminosity in /pb",
)
parser.add_argument(
    "-o",
    "--output",
    type=str,
    required=True,
    help="output root file name",
)
parser.add_argument(
    "--mergemap",
    type=str,
    default=None,
    help="Specify mergemap as dict, '{merge1:[dataset1,dataset2]...}' Also works with the json file with dict",
)
parser.add_argument(
    "--autorebin",
    default=None,
    help="Rebin the plotting variables, input `int` or `list`. int: merge N bins. list of number: rebin edges(non-uniform bin is possible)",
)


def create_template(inputs, variable, mergemap, axis, lumi, output):
    inputs = scaleSumW(inputs, lumi)
    collated = collate(inputs, mergemap)
    mergemap_all = {"all": [s for s in inputs]}
    merge_all = collate(inputs, {"all": [s for s in inputs]})
    hall = merge_all["all"][variable][axis]

    fout = uproot.recreate(f"{variable}_{output}")
    subtract = False
    if "osss" in axis.keys():
        if "sum" in str(axis["osss"]):
            subtract = True
            axis.pop("osss")
    for data in collated.keys():
        if subtract:
            collated[data][variable] = (
                collated[data][variable][{"osss": 0}]
                + collated[data][variable][{"osss": 1}] * -1
            )
        if arg.autorebin is not None:
            if arg.autorebin.isdigit():
                rebin = int(arg.autorebin)
            else:
                rebin = np.array([float(i) for i in arg.autorebin.split(",")])
            fout[data.replace("-", "_")] = rebin_hist(
                collated[data][variable][axis],
                collated[data][variable].axes[-1].name,
                rebin,
            )
        else:
            fout[data.replace("-", "_")] = collated[data][variable][axis]
        if variable not in collated[data].keys():
            raise f"{variable} not in {data}"
        elif np.sum(collated[data][variable].values()) == 0:
            print(f"0 value in {variable} histogram of {data}, set empty histogram")
            fout[data.replace("-", "_")] = hall
        else:
            fout[data.replace("-", "_")] = collated[data][variable][axis]


if __name__ == "__main__":
    import glob

    arg = parser.parse_args()
    if len(arg.input.split(",")) > 1:
        inputs = {i: load(i) for i in arg.input.split(",")}
    elif "*" in arg.input:
        files = glob.glob(arg.input)
        inputs = {i: load(i) for i in files}
    else:
        inputs = {arg.input: load(arg.input)}
    import json

    if arg.mergemap != None:
        mergemap = (
            json.load(open(arg.mergemap))
            if arg.mergemap.endswith(".json")
            else json.loads(arg.mergemap)
        )
    # if not specify the mergemap, stored per dataset
    else:
        flists = []
        for f in inputs.keys():
            for d in inputs[f].keys():
                flists.append(d)
        mergemap = {d: d for d in flists}

    axis = json.loads(arg.axis)
    # modify axis to hist methods
    for a in axis:
        if axis[a] == "sum":
            axis[a] = sum
    create_template(inputs, arg.variable, mergemap, axis, arg.lumi, arg.output)
