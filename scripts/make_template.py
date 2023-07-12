from coffea.util import load
import uproot, sys, os, argparse, hist
from BTVNanoCommissioning.helpers.xs_scaler import collate, scaleSumW

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


def create_template(inputs, variable, mergemap, axis, lumi, output):
    inputs = scaleSumW(inputs, lumi)
    collated = collate(inputs, mergemap)
    fout = uproot.recreate(f"{variable}_{output}")
    for data in collated.keys():
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
