import pathlib
import glob
from coffea.util import load, save
import argparse
import tqdm

arg = argparse.ArgumentParser("Combine coffea output files")

arg.add_argument(
    "-i",
    "--input",
    required=True,
    help="Input coffea files split by a comma. Glob wildcards are also supported.",
)
arg.add_argument(
    "-o", "--output", required=True, help="Output combined coffea file name."
)

args = arg.parse_args()

files = args.input.split(",")
files = [glob.glob(f) for f in files]
files = [item for sublist in files for item in sublist]

output = {}

for f in tqdm.tqdm(files, desc="Combining coffea files"):
    try:
        data = load(f)
    except Exception as e:
        print(f"Error loading {f}: {e}")
        continue
    for key in data:
        if key not in output:
            output[key] = data[key]
        else:
            # Check if a dictionary
            if isinstance(data[key], dict):
                for subkey in data[key]:
                    if subkey not in output[key]:
                        output[key][subkey] = data[key][subkey]
                    else:
                        output[key][subkey] += data[key][subkey]
            else:
                output[key] += data[key]

pathlib.Path(args.output).parent.mkdir(parents=True, exist_ok=True)
save(output, args.output)