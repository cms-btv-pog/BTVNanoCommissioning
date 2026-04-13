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
    help="Input coffea files separated by a comma. Glob wildcards are also supported.",
)
arg.add_argument(
    "-o", "--output", required=True, help="Output combined coffea file name."
)

args = arg.parse_args()

files = args.input.split(",")
files = [glob.glob(f) for f in files]
files = [item for sublist in files for item in sublist]

output = {}

keys = set()
for f in tqdm.tqdm(files, desc="Combining coffea files"):
    try:
        data = load(f)
    except Exception as e:
        print(f"Error loading {f}: {e}")
        continue
    for key in data:
        print(f"Processing key: {key} from file: {f}")
        keys.add(key)
        if key not in output:
            output[key] = data[key]
        else:
            # Check if a dictionary
            if isinstance(data[key], dict):
                for subkey in data[key]:
                    if subkey not in output[key]:
                        output[key][subkey] = data[key][subkey]
                    else:
                        try:
                            output[key][subkey] += data[key][subkey]
                        except ValueError as e:
                            print(f"ValueError combining {key}[{subkey}] from {f}")
                            print(f"Output type: {type(output[key][subkey])}, Data type: {type(data[key][subkey])}")
                            print(f"Output value: {output[key][subkey]},\nData value: {data[key][subkey]}")
                            print(f"Error: {e}")
            else:
                try:
                    output[key] += data[key]
                except ValueError as e:
                    print(f"ValueError combining {key} from {f}")
                    print(f"Output type: {type(output[key])}, Data type: {type(data[key])}")
                    print(f"Output value: {output[key]},\nData value: {data[key]}")
                    print(f"Error: {e}")

pathlib.Path(args.output).parent.mkdir(parents=True, exist_ok=True)
save(output, args.output)

sortkeys = sorted(keys)
print(f"Combined {len(files)} files with keys: {sortkeys} into {args.output}")