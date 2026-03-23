import psutil
import os
import uproot
import warnings
import argparse
import glob
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import awkward as ak
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
from BTVNanoCommissioning.helpers.definitions import get_definitions
from BTVNanoCommissioning.helpers.definitions import get_discriminators

# Suppress the specific FutureWarning from uproot
warnings.filterwarnings("ignore", category=FutureWarning, module="uproot")
warnings.filterwarnings("ignore", category=FutureWarning, module="seaborn")

definitions_dict = get_definitions()
disc_list = get_discriminators()

filtered_names = [
    name for name in disc_list if name.endswith("B") and not name.endswith("CvB")
]


def is_within_range(array, min_value, max_value):
    """Check if all values in the array are within the specified range."""
    return ak.all((array >= min_value) & (array <= max_value))


def load_single_file(
    file_path,
    base_dir,
    chunk_size=100000,
    temp_subdir="temp_data",
    limit_inputs=False,
    limit_outputs=False,
    SMu=False,
    flavour_split=False,
):
    temp_dir = os.path.join(base_dir, temp_subdir)
    if not os.path.exists(temp_dir):
        os.makedirs(temp_dir, exist_ok=True)

    with uproot.open(file_path) as file:
        if "Events" not in file:
            print(f"Skipping file without 'Events' key: {file_path}")
            return

        tree = file["Events"]
        if tree.num_entries == 0:
            print(f"Skipping empty file: {file_path}")
            return

        # Apply limit_inputs filter
        if limit_inputs:
            filtered_definitions_dict = {
                k: v
                for k, v in definitions_dict.items()
                if any(sub in k for sub in ["Npfcan", "Cpfcan", "sv"])
                and k.endswith("_0")
                or not any(sub in k for sub in ["Npfcan", "Cpfcan", "sv"])
            }
        else:
            filtered_definitions_dict = definitions_dict

        # Initialize filtered_names
        filtered_names = [
            name
            for name in disc_list
            if name.endswith("B") and not name.endswith("CvB")
        ]

        # Apply limit_outputs filter
        if limit_outputs:
            filtered_names = [
                name
                for name in filtered_names
                if name.endswith("B")
                and "Neg" not in name
                and "btagPNetProbB" not in name
            ]

        # Filter branches to include only those with 'SelJets' in their name, containing keys from filtered_definitions_dict, and matching filtered_names
        branches = [
            branch
            for branch in tree.keys()
            if "MuonJet" in branch
            and (
                any(key in branch for key in filtered_definitions_dict.keys())
                or any(name in branch for name in filtered_names)
            )
        ]

        if limit_outputs:
            branches = [
                branch
                for branch in branches
                if not any(name in branch for name in filtered_names)
                or any(
                    name in branch for name in filtered_names if branch.endswith("B")
                )
            ]

        # Include branches that have 'SMu' in their name if the SMu flag is True
        if SMu:
            smu_branches = [branch for branch in tree.keys() if "SoftMuon" in branch]
            smu_branches.append("MuonJet_muEF")
            if limit_inputs:
                smu_branches = [
                    "SoftMuon_tunepRelPt",
                    "SoftMuon_pfRelIso03_chg",
                    "SoftMuon_eta",
                    "SoftMuon_phi",
                    "SoftMuon_jetPtRelv2",
                    "SoftMuon_dxy",
                    "SoftMuon_dxyErr",
                    "SoftMuon_jetRelIso",
                    "SoftMuon_sip3d",
                    "SoftMuon_dzErr",
                    "SoftMuon_pfRelIso04_all",
                    "SoftMuon_ip3d",
                    "SoftMuon_pt",
                    "SoftMuon_ptErr",
                    "SoftMuon_tkRelIso",
                    "SoftMuon_dz",
                    "SoftMuon_pfRelIso03_all",
                    "MuonJet_muEF",
                    "SoftMuon_charge",
                ]
            branches.extend(smu_branches)
            branches = list(set(branches))  # Remove duplicates

        # Include the flavour column if flavour_split is enabled
        if flavour_split:
            branches.append("MuonJet_hadronFlavour")

        # Extract manual ranges for x and y columns if they don't start with SelJet_ or SMu_
        x_limits_dict = {}
        for x_col in branches:
            if not x_col.startswith(("SelJet_", "SoftMuon_", "MuonJet_")):
                x_limits = definitions_dict.get(x_col, {}).get("manual_ranges", None)
            else:
                x_col_str = (
                    x_col.lstrip("SelJet_").lstrip("SoftMuon_").lstrip("MuonJet_")
                )
                x_limits = definitions_dict.get(x_col_str, {}).get(
                    "manual_ranges", None
                )

            x_limits_dict[x_col] = x_limits

        for i, data_chunk in enumerate(
            tree.iterate(branches, library="pd", step_size=chunk_size)
        ):
            if isinstance(data_chunk, tuple):
                data_chunk = data_chunk[0]  # Extract the DataFrame from the tuple

            # Filter the data based on MuJet_Cpfcan_ptrel_0
            ###data_chunk = data_chunk[data_chunk['MuJet_DeepJet_Cpfcan_ptrel_0'] >= -20]

            # Filter the data based on x_limits
            for col, limits in x_limits_dict.items():
                if limits is not None:
                    min_value, max_value = limits
                    data_chunk = data_chunk[
                        (data_chunk[col] >= min_value) & (data_chunk[col] <= max_value)
                    ]
            temp_file_path = os.path.join(
                temp_dir, f"{os.path.basename(file_path)}_chunk_{i}.parquet"
            )
            data_chunk.to_parquet(temp_file_path)


def inspect_first_file(
    file_path, limit_inputs=False, limit_outputs=False, SMu=False, flavour_split=False
):
    with uproot.open(file_path) as file:
        print(f"Inspecting file: {file_path}")
        print("Keys:", file.keys())
        if "Events" in file:
            tree = file["Events"]
            print("Number of entries in 'Events':", tree.num_entries)

            # Apply limit_inputs filter
            if limit_inputs:
                filtered_definitions_dict = {
                    k: v
                    for k, v in definitions_dict.items()
                    if any(sub in k for sub in ["Npfcan", "Cpfcan", "sv"])
                    and k.endswith("_0")
                    or not any(sub in k for sub in ["Npfcan", "Cpfcan", "sv"])
                }
            else:
                filtered_definitions_dict = definitions_dict

            # Initialize filtered_names
            filtered_names = [
                name
                for name in disc_list
                if name.endswith("B") and not name.endswith("CvB")
            ]

            # Apply limit_outputs filter
            if limit_outputs:
                filtered_names = [
                    name
                    for name in filtered_names
                    if name.endswith("B")
                    and "Neg" not in name
                    and "btagPNetProbB" not in name
                ]

            print(tree.keys())
            # Filter branches to include only those with 'SelJets' in their name, containing keys from definitions_dict, and matching filtered_names

            branches = [
                branch
                for branch in tree.keys()
                if "MuJet" in branch
                and (
                    any(key in branch for key in filtered_definitions_dict.keys())
                    or any(name in branch for name in filtered_names)
                )
            ]
            if limit_outputs:
                branches = [
                    branch
                    for branch in branches
                    if not any(name in branch for name in filtered_names)
                    or any(
                        name in branch
                        for name in filtered_names
                        if branch.endswith(name)
                    )
                ]
            print("Filtered Branches:", branches)
            print(len(branches))
            if SMu:
                smu_branches = [
                    branch for branch in tree.keys() if "SoftMuon" in branch
                ]
                smu_branches.append("MuonJet_muEF")
                print(smu_branches)
                if limit_inputs:
                    smu_branches = [
                        "SoftMuon_tunepRelPt",
                        "SoftMuon_pfRelIso03_chg",
                        "SoftMuon_eta",
                        "SoftMuon_phi",
                        "SoftMuon_jetPtRelv2",
                        "SoftMuon_dxy",
                        "SoftMuon_dxyErr",
                        "SoftMuon_jetRelIso",
                        "SoftMuon_sip3d",
                        "SoftMuon_dzErr",
                        "SoftMuon_pfRelIso04_all",
                        "SoftMuon_ip3d",
                        "SoftMuon_pt",
                        "SoftMuon_ptErr",
                        "SoftMuon_tkRelIso",
                        "SoftMuon_dz",
                        "SoftMuon_pfRelIso03_all",
                        "MuonJet_muEF",
                        "SoftMuon_charge",
                    ]

                branches.extend(smu_branches)
                branches = list(set(branches))  # Remove duplicates
                print("Branches with SMu:", smu_branches)

            # Include the flavour column if flavour_split is enabled
            if flavour_split:
                branches.append("MuonJet_hadronFlavour")

            print("Branches:", branches)

            # Find branches that contain any of the keys from definitions_dict
            matching_keys = [
                branch
                for branch in branches
                if any(key in branch for key in definitions_dict.keys())
            ]
            print("Branches containing keys from definitions_dict:", matching_keys)
            print(len(matching_keys))

            matching_keys = [
                branch
                for branch in branches
                if any(name in branch for name in filtered_names)
            ]
            print("Branches containing keys from definitions_dict:", matching_keys)
            print(len(matching_keys))

            # Print the filtered names
            print("Filtered Names:", filtered_names)

            # Extract manual ranges for x and y columns if they don't start with SelJet_ or SMu_
            x_limits_dict = {}
            for x_col in branches:
                if not x_col.startswith(("SelJet_", "SoftMuon_", "MuonJet_")):
                    x_limits = definitions_dict.get(x_col, {}).get(
                        "manual_ranges", None
                    )
                else:
                    x_col_str = (
                        x_col.lstrip("SelJet_").lstrip("SoftMuon_").lstrip("MuonJet_")
                    )
                    x_limits = definitions_dict.get(x_col_str, {}).get(
                        "manual_ranges", None
                    )

                x_limits_dict[x_col] = x_limits
                # Debugging prints
                # print(f"x_col: {x_col}")
                # print(f"Manual ranges for {x_col}: {x_limits}")

            # Load the first chunk to inspect memory usage
            data_chunk = next(tree.iterate(branches, library="pd", step_size=100000))
            if isinstance(data_chunk, tuple):
                data_chunk = data_chunk[0]  # Extract the DataFrame from the tuple

            # Filter the data based on x_limits
            for col, limits in x_limits_dict.items():
                if limits is not None:
                    min_value, max_value = limits
                    data_chunk = data_chunk[
                        (data_chunk[col] >= min_value) & (data_chunk[col] <= max_value)
                    ]

            memory_usage = data_chunk.memory_usage(deep=True).sum()
            print(
                f"Memory usage of the first DataFrame chunk: {memory_usage / (1024 ** 2):.2f} MB"
            )
        else:
            print("'Events' key not found in the file.")


def load_data(
    file_paths,
    base_dir,
    batch_size=2,
    temp_subdir="temp_data",
    output_subdir="parquet_data",
    max_files=None,
    limit_inputs=False,
    limit_outputs=False,
    SMu=False,
    flavour_split=False,
):
    # Step 1: Inspect the nominal folder
    nominal_folder = os.path.join(base_dir, "nominal")
    nominal_subfolders = {}
    if os.path.exists(nominal_folder):
        for subfolder in os.listdir(nominal_folder):
            subfolder_path = os.path.join(nominal_folder, subfolder)
            if os.path.isdir(subfolder_path):
                # Count the number of .root files in the subfolder
                root_files = [
                    f for f in os.listdir(subfolder_path) if f.endswith(".root")
                ]
                root_file_count = len(root_files)
                # Print the count of .root files
                print(
                    f"Subfolder: {subfolder}, Number of .root files: {root_file_count}"
                )
                # Store the subfolder path for later use
                nominal_subfolders[subfolder] = subfolder_path

    if max_files is not None:
        file_paths = file_paths[:max_files]

    output_dir = os.path.join(base_dir, output_subdir)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with ThreadPoolExecutor(max_workers=batch_size) as executor:
        futures = [
            executor.submit(
                load_single_file,
                file_path,
                base_dir,
                temp_subdir=temp_subdir,
                limit_inputs=limit_inputs,
                limit_outputs=limit_outputs,
                SMu=SMu,
                flavour_split=flavour_split,
            )
            for file_path in file_paths
        ]
        for future in tqdm(
            as_completed(futures), total=len(futures), desc="Processing files"
        ):
            future.result()  # Ensure any exceptions are raised

    combined_data = {}
    temp_files = glob.glob(os.path.join(base_dir, temp_subdir, "*.parquet"))
    print(f"Temporary files found: {len(temp_files)}")  # Debugging line
    event_counts = []
    for temp_file in temp_files:
        df = pd.read_parquet(temp_file)
        key = os.path.basename(temp_file).replace(".parquet", "")
        combined_data[key] = df  # Store the DataFrame directly
        event_counts.append((key, len(df)))  # Store the number of events and the key

    # Define ranking factors for each subfolder
    ### FIXME: sumw has to change, whenever you run a different set. Pay attention!!!
    ### FIXME: The former factor is the xs, the second - sumw!
    ranking_factors = {
        "QCD_PT-15to20_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8": 295600 / 142083,
        "QCD_PT-20to30_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8": 2689000 / 5926,
        "QCD_PT-30to50_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8": 1442000 / 339153,
        "QCD_PT-50to80_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8": 405800 / 94942,
        "QCD_PT-80to120_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8": 96060 / 488561,
        "QCD_PT-120to170_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8": 23230 / 210233,
        "QCD_PT-1700to300_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8": 7763 / 216617,
        "QCD_PT-300to470_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8": 701.4 / 196409,
        "QCD_PT-470to600_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8": 68.24 / 254794,
        "QCD_PT-600to800_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8": 21.23 / 34096,
        "QCD_PT-800to1000_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8": 3.9 / 155007,
        "QCD_PT-1000_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8": 1.323 / 237819,
    }

    # Sort event counts by the number of events in descending order
    event_counts.sort(key=lambda x: x[1], reverse=True)

    # Display the ranking
    print("Ranking by number of events:")
    for rank, (key, count) in enumerate(event_counts, start=1):
        # Strip the _chunk_0 suffix
        stripped_key = key.replace("_chunk_0", "")
        # Determine the parent subfolder
        parent_subfolder = None
        for subfolder, path in nominal_subfolders.items():
            if stripped_key in os.listdir(path):
                parent_subfolder = subfolder
                break
        print(
            f"{rank}. File: {key}, Number of events: {count}, Parent subfolder: {parent_subfolder}"
        )

    # Calculate weighted event counts
    weighted_event_counts = []
    for key, count in event_counts:
        stripped_key = key.replace("_chunk_0", "")
        parent_subfolder = None
        for subfolder, path in nominal_subfolders.items():
            if stripped_key in os.listdir(path):
                parent_subfolder = subfolder
                break
        factor = ranking_factors.get(
            parent_subfolder, 1.0
        )  # Default factor is 1.0 if not specified
        weighted_count = count * factor
        weighted_event_counts.append((key, weighted_count, parent_subfolder))

    # Sort weighted event counts by the weighted number of events in descending order
    weighted_event_counts.sort(key=lambda x: x[1], reverse=True)

    # Display the ranking
    ### Ranking done by the yields: from the highest to the lowest
    print("Ranking by weighted number of events:")
    for rank, (key, weighted_count, parent_subfolder) in enumerate(
        weighted_event_counts, start=1
    ):
        print(
            f"{rank}. File: {key}, Weighted number of events: {weighted_count}, Parent subfolder: {parent_subfolder}"
        )

    print(f"Combined data keys: {list(combined_data.keys())}")  # Debugging line
    return combined_data


variables_to_zoom = {
    "MuonJet_DeepCSV_flightDistance2dSig": [2, 4],
    "MuonJet_DeepCSV_flightDistance2dVal": [0, 0.4],
    "MuonJet_DeepCSV_flightDistance3dSig": [1, 10],
    "MuonJet_DeepCSV_flightDistance3dVal": [0, 0.5],
    "MuonJet_DeepCSV_jetNSelectedTracks": [5, 12],
    "MuonJet_DeepCSV_trackJetPt": [25, 60],
    "MuonJet_DeepCSV_trackSip2dSigAboveCharm": [0, 2],
    "MuonJet_DeepCSV_trackSip2dValAboveCharm": [0, 0.01],
    "MuonJet_DeepCSV_trackSip3dSigAboveCharm": [0, 2],
    "MuonJet_DeepCSV_trackSip3dValAboveCharm": [0, 0.01],
    "MuonJet_DeepCSV_trackSumJetDeltaR": [0, 0.03],
    "MuonJet_DeepCSV_vertexEnergyRatio": [0, 0.4],
    "MuonJet_DeepCSV_vertexJetDeltaR": [0, 0.07],
    "SoftMuon_ip3d": [0, 0.05],
    "SoftMuon_jetPtRelv2": [0, 1],
    "SoftMuon_jetRelIso": [0, 30],
    "SoftMuon_pfRelIso03_all": [0, 20],
    "SoftMuon_pfRelIso03_chg": [0, 20],
    "SoftMuon_pfRelIso04_all": [0, 20],
    "SoftMuon_sip3d": [0, 15],
}


# Function to print keys of definitions_dict
def print_definitions_keys():
    keys = definitions_dict.keys()
    print("Keys in definitions_dict:", keys)
    print(len(keys))


def plot_2d_histogram(
    data, x_col, y_col, output_dir, x_limits=None, y_limits=None, x_bins=50, y_bins=50
):
    plt.figure(figsize=(10, 8))
    hist = sns.histplot(
        data, x=x_col, y=y_col, bins=(x_bins, y_bins), pthresh=0.1, cmap="viridis"
    )

    if x_limits is not None:
        plt.xlim(x_limits)
    if y_limits is not None:
        plt.ylim(y_limits)

    plt.xlabel(x_col)
    plt.ylabel(y_col)
    plt.title(f"2D Histogram of {x_col} vs {y_col}")
    plt.colorbar(hist.collections[0], label="Counts")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"{x_col}_vs_{y_col}.png"))
    plt.close()


def plot_all_histograms(
    data,
    definitions_dict,
    filtered_names,
    output_dir,
    include_smu=False,
    flavour_split=False,
    split_region_b=False,
    zoom=False,
):
    print("Data columns:", data.columns)
    print("Definitions dict keys:", definitions_dict.keys())
    print("Filtered names:", filtered_names)

    # Extend definitions_dict with SMu columns if include_smu is True
    if include_smu:
        smu_columns = [col for col in data.columns if "SMu" in col]
        for col in smu_columns:
            if col not in definitions_dict:
                definitions_dict[col] = col  # Add SMu columns to definitions_dict

    print("Definitions dict keys:", definitions_dict.keys())

    if flavour_split:
        flavour_column = "MuJet_hadronFlavour"

        flavour_groups = {
            "light": data[data[flavour_column] == 0],
            "charm": data[data[flavour_column] == 4],
            "bottom": data[data[flavour_column] == 5],
        }
        for flavour in flavour_groups.keys():
            flavour_data = flavour_groups[flavour]
            flavour_dir = os.path.join(output_dir, f"flavour_{flavour}")
            os.makedirs(flavour_dir, exist_ok=True)
            plot_histograms_for_data(
                flavour_data, definitions_dict, filtered_names, flavour_dir, zoom
            )

    elif split_region_b:
        deepflavb_column = "MuonJet_btagDeepFlavB"
        high_probB_data = data[data[deepflavb_column] > 0.5]
        low_probB_data = data[data[deepflavb_column] <= 0.5]

        high_probB_dir = os.path.join(output_dir, "high_probB")
        low_probB_dir = os.path.join(output_dir, "low_probB")

        os.makedirs(high_probB_dir, exist_ok=True)
        os.makedirs(low_probB_dir, exist_ok=True)

        plot_histograms_for_data(
            high_probB_data, definitions_dict, filtered_names, high_probB_dir, zoom
        )
        plot_histograms_for_data(
            low_probB_data, definitions_dict, filtered_names, low_probB_dir, zoom
        )

    else:
        plot_histograms_for_data(
            data, definitions_dict, filtered_names, output_dir, zoom
        )


def plot_histograms_for_data(
    data, definitions_dict, filtered_names, output_dir, zoom=False
):
    for input_col in definitions_dict.keys():
        for output_col in filtered_names:
            matching_input_cols = [col for col in data.columns if input_col in col]
            matching_output_cols = [col for col in data.columns if output_col in col]

            for x_col in matching_input_cols:
                for y_col in matching_output_cols:
                    print(f"Plotting {x_col} vs {y_col}")

                    # Extract manual ranges for x and y columns if they don't start with SelJet_ or SMu_
                    if not x_col.startswith(("SelJet_", "SMu_", "MuJet_")):
                        x_limits = definitions_dict.get(x_col, {}).get(
                            "manual_ranges", None
                        )
                        x_bins = definitions_dict.get(x_col, {}).get("bins", 50)
                    else:
                        x_col_str = (
                            x_col.lstrip("SelJet_").lstrip("SMu_").lstrip("MuJet_")
                        )
                        x_limits = definitions_dict.get(x_col_str, {}).get(
                            "manual_ranges", None
                        )
                        x_bins = definitions_dict.get(x_col_str, {}).get("bins", 50)

                    if not y_col.startswith(("SelJet_", "SMu_", "MuJet_")):
                        y_limits = definitions_dict.get(y_col, {}).get(
                            "manual_ranges", None
                        )
                        y_bins = definitions_dict.get(x_col, {}).get("bins", 50)
                    else:
                        y_col_str = (
                            y_col.lstrip("SelJet_").lstrip("SMu_").lstrip("MuJet_")
                        )
                        y_limits = definitions_dict.get(y_col_str, {}).get(
                            "manual_ranges", None
                        )
                        y_bins = definitions_dict.get(y_col_str, {}).get("bins", 50)

                    # Debugging prints
                    print(f"x_col: {x_col}, y_col: {y_col}")
                    print(
                        f"definitions_dict[{x_col}]: {definitions_dict.get(x_col, 'Not found')}"
                    )
                    print(
                        f"definitions_dict[{y_col}]: {definitions_dict.get(y_col, 'Not found')}"
                    )
                    print(f"Manual ranges for {x_col}: {x_limits}")
                    print(f"Manual ranges for {y_col}: {y_limits}")

                    print(f"Bin count for {x_col}: {x_bins}")
                    print(f"Bin count for {y_col}: {y_bins}")
                    # If zoom is active, create zoomed histograms
                    if zoom and x_col in variables_to_zoom:
                        zoomed_dir = os.path.join(output_dir, "zoomed_hists")
                        os.makedirs(zoomed_dir, exist_ok=True)

                        x_zoom = variables_to_zoom[x_col]
                        plot_2d_histogram(
                            data,
                            x_col,
                            y_col,
                            zoomed_dir,
                            x_limits=x_zoom,
                            y_limits=y_limits,
                            x_bins=40,
                            y_bins=y_bins,
                        )

                    plot_2d_histogram(
                        data,
                        x_col,
                        y_col,
                        output_dir,
                        x_limits=x_limits,
                        y_limits=y_limits,
                        x_bins=x_bins,
                        y_bins=y_bins,
                    )


def main():
    parser = argparse.ArgumentParser(description="Generate correlation plots.")
    parser.add_argument(
        "folder", type=str, help="The folder containing the data files."
    )
    parser.add_argument(
        "--max_files",
        type=int,
        default=None,
        help="Maximum number of files to process.",
    )
    parser.add_argument(
        "--limit_inputs",
        action="store_true",
        help="Limit inputs to variables with npf, cpf, sv and _0",
    )
    parser.add_argument(
        "--limit_outputs",
        action="store_true",
        help="Limit outputs to names ending with B and not containing Neg",
    )
    parser.add_argument(
        "--SMu", action="store_true", help="Include branches with SMu in their name"
    )
    parser.add_argument(
        "--flavour_split",
        action="store_true",
        help="Split correlations by hadron flavour",
    )
    parser.add_argument(
        "--split_region_b", action="store_true", help="Split correlations by region B"
    )
    parser.add_argument(
        "--zoom", action="store_true", help="Zoom in on the correlation plots"
    )
    args = parser.parse_args()

    data_file_paths = glob.glob(
        os.path.join(args.folder, "**", "*.root"), recursive=True
    )
    print(f"Data file paths: {data_file_paths}")  # Debugging line
    # print(data_file_paths)
    # print_definitions_keys()

    if data_file_paths:
        inspect_first_file(
            data_file_paths[0],
            limit_inputs=args.limit_inputs,
            limit_outputs=args.limit_outputs,
            SMu=args.SMu,
            flavour_split=args.flavour_split,
        )

    parquet_dir = os.path.join(args.folder, "parquet_data")
    histograms_dir = os.path.join(args.folder, "2D_histograms")
    os.makedirs(histograms_dir, exist_ok=True)

    data_dict = load_data(
        data_file_paths,
        args.folder,
        max_files=args.max_files,
        limit_inputs=args.limit_inputs,
        limit_outputs=args.limit_outputs,
        SMu=args.SMu,
        flavour_split=args.flavour_split,
    )

    # Concatenate all DataFrames if load_data returns a dictionary of DataFrames
    if isinstance(data_dict, dict):
        data = pd.concat(data_dict.values(), ignore_index=True)
    else:
        data = data_dict

    # Print data columns for debugging
    print("Data columns after loading and concatenation:", data.columns)

    # Plot 2D histograms of inputs vs outputs
    plot_all_histograms(
        data,
        definitions_dict,
        filtered_names,
        histograms_dir,
        include_smu=args.SMu,
        flavour_split=args.flavour_split,
        split_region_b=args.split_region_b,
        zoom=args.zoom,
    )


if __name__ == "__main__":
    main()
