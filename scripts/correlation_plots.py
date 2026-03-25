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
from BTVNanoCommissioning.helpers.definitions import get_definitions, get_discriminators

# Suppress the specific FutureWarning from uproot
warnings.filterwarnings("ignore", category=FutureWarning, module="uproot")

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
    specify_MC=False,
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

    print(base_dir)
    if "MC" in base_dir:
        print("MC file detected")

    # Check if the base directory name contains "MC" and the specify_MC flag is active
    if "MC" in base_dir and specify_MC:
        # Filter file_paths to only include files from a specific subfolder
        print("Filtering file paths to include files from a specific subfolder")
        # specific_subfolder = "QCD_PT-50to80_MuEnrichedPt5_TuneCP5_13p6TeV_pythia8"  # Replace with the actual subfolder name
        specific_subfolder = "QCD_PT-80to120_TuneCP5_13p6TeV_pythia8"
        file_paths = [fp for fp in file_paths if specific_subfolder in fp]

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
        "QCD_PT-50to80_TuneCP5_13p6TeV_pythia8": 16760000.0 / 306000,
        "QCD_PT-80to120_TuneCP5_13p6TeV_pythia8": 2514000.0 / 362000,
        "QCD_PT-120to170_TuneCP5_13p6TeV_pythia8": 442300.0 / 352000,
        "QCD_PT-170to300_TuneCP5_13p6TeV_pythia8": 113400.0 / 334000,
        "QCD_PT-300to470_TuneCP5_13p6TeV_pythia8": 7610.0 / 128000,
        "QCD_PT-470to600_TuneCP5_13p6TeV_pythia8": 626.6 / 223000,
        "QCD_PT-600to800_TuneCP5_13p6TeV_pythia8": 179.5 / 230000,
        "QCD_PT-800to1000_TuneCP5_13p6TeV_pythia8": 30.69 / 38000,
        "QCD_PT-1000to1400_TuneCP5_13p6TeV_pythia8": 8.956 / 143000,
        "QCD_PT-1400to1800_TuneCP5_13p6TeV_pythia8": 0.8112 / 219000,
        "QCD_PT-1800to2400_TuneCP5_13p6TeV_pythia8": 0.1152 / 211000,
        "QCD_PT-2400to3200_TuneCP5_13p6TeV_pythia8": 0.007574 / 334000,
        "QCD_PT-3200_TuneCP5_13p6TeV_pythia8": 0.0002313 / 44000,
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


def compute_correlations(data, SMu=False, flavour_split=False, split_region_B=False):
    if not data:
        raise ValueError("No data available to compute correlations.")

    # Combine all data into a single DataFrame
    combined_df = pd.concat(data.values(), axis=0)
    print(f"Combined DataFrame shape: {combined_df.shape}")  # Debugging line
    print(f"Combined DataFrame columns: {combined_df.columns}")  # Debugging line

    combined_df = combined_df[sorted(combined_df.columns)]
    print(f"Sorted DataFrame columns: {combined_df.columns}")  # Debugging line

    # Separate columns into two groups
    filtered_columns = [
        col
        for col in combined_df.columns
        if any(name in col for name in filtered_names)
    ]
    definition_columns = [
        col
        for col in combined_df.columns
        if any(key in col for key in definitions_dict.keys())
    ]

    # Include SMu columns if the SMu flag is True
    if SMu:
        smu_columns = [col for col in combined_df.columns if "SoftMuon" in col]
        smu_columns.append("MuonJet_muEF")  # Ensure MuJet_muEF is included
        definition_columns.extend(smu_columns)
        definition_columns = list(set(definition_columns))  # Remove duplicates

        # Ensure SMu columns are at the end of the definition_columns list
        smu_columns = [col for col in definition_columns if "SoftMuon" in col]
        non_smu_columns = [col for col in definition_columns if "SoftMuon" not in col]
        definition_columns = non_smu_columns + smu_columns

    # Compute the full correlation matrix
    all_columns = sorted(filtered_columns + definition_columns)
    if flavour_split:
        flavour_column = "MuonJet_hadronFlavour"
        if flavour_column not in combined_df.columns:
            raise ValueError(f"Column '{flavour_column}' not found in data.")

        flavour_groups = {
            "light": combined_df[combined_df[flavour_column] == 0],
            "charm": combined_df[combined_df[flavour_column] == 4],
            "bottom": combined_df[combined_df[flavour_column] == 5],
        }

        correlation_matrices = {}
        for flavour, df in flavour_groups.items():
            correlation_matrices[flavour] = df[all_columns].corr()

        return correlation_matrices, filtered_columns, definition_columns
    elif split_region_B:
        deepflavb_column = "MuonJet_btagDeepFlavB"
        if deepflavb_column not in combined_df.columns:
            raise ValueError(f"Column '{deepflavb_column}' not found in data.")

        # Split the DataFrame based on DeepFlavB
        high_probB_df = combined_df[combined_df[deepflavb_column] > 0.5]
        low_probB_df = combined_df[combined_df[deepflavb_column] <= 0.5]

        # Calculate correlation matrices for high and low DeepFlavB
        high_probB_corr_matrix = high_probB_df[all_columns].corr()
        low_probB_corr_matrix = low_probB_df[all_columns].corr()

        # Return the correlation matrices along with the filtered and definition columns
        return (
            {"high_probB": high_probB_corr_matrix, "low_probB": low_probB_corr_matrix},
            filtered_columns,
            definition_columns,
        )
    else:
        correlation_matrix = combined_df[all_columns].corr()
        return correlation_matrix, filtered_columns, definition_columns


# Function to print keys of definitions_dict
def print_definitions_keys():
    keys = definitions_dict.keys()
    print("Keys in definitions_dict:", keys)
    print(len(keys))


def create_label_mapping(columns):
    sorted_columns = sorted(columns)
    return {col: f"Label {i+1}" for i, col in enumerate(sorted_columns)}


def plot_correlations(
    correlation_data,
    filtered_columns,
    definition_columns,
    output_dir,
    flavour_split=False,
    split_region_B=False,
    threshold=0.5,
):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    def plot_heatmap(matrix, title, filename, label_mapping, threshold):
        num_labels = len(matrix.columns)
        font_size = (
            10 if num_labels < 30 else 6
        )  # Adjust font size based on number of labels

        plt.figure(figsize=(24, 16))  # Increase the size of the plot
        sns.heatmap(
            matrix,
            annot=False,
            cmap="coolwarm",
            fmt=".2f",
            cbar_kws={"label": "Correlation Coefficient"},
        )

        # Add annotations manually
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                value = matrix.iloc[i, j]
                if abs(value) >= threshold:
                    plt.text(
                        j + 0.5,
                        i + 0.5,
                        f"{value:.2f}",
                        ha="center",
                        va="center",
                        fontsize=font_size,
                    )

        plt.title(title)

        # Replace x and y labels with numbers or symbols
        ax = plt.gca()
        ax.set_xticklabels(
            [
                label_mapping.get(label.get_text(), label.get_text())
                for label in ax.get_xticklabels()
            ]
        )
        ax.set_yticklabels(
            [
                label_mapping.get(label.get_text(), label.get_text())
                for label in ax.get_yticklabels()
            ]
        )

        # Add legend to the right of the colorbar
        colorbar = ax.collections[0].colorbar
        colorbar.ax.set_position([0.85, 0.1, 0.03, 0.8])  # Adjust position as needed

        # Create legend handles and sort them alphabetically
        sorted_label_mapping = dict(sorted(label_mapping.items()))
        legend_handles = [
            plt.Line2D([0], [0], color="w", label=f"{v}: {k}")
            for k, v in sorted_label_mapping.items()
        ]

        # Determine the number of columns for the legend
        ncol = 2 if len(legend_handles) > 20 else 1

        # Add legend
        plt.legend(
            handles=legend_handles,
            loc="center left",
            bbox_to_anchor=(1.2, 0.5),
            title="Legend",
            ncol=ncol,
        )

        plt.tight_layout()
        plt.savefig(filename, bbox_inches="tight")
        plt.close()

    # Create label mappings
    filtered_label_mapping = create_label_mapping(filtered_columns)
    definition_label_mapping = create_label_mapping(definition_columns)

    # Sort columns and index alphabetically
    sorted_filtered_columns = sorted(filtered_columns)
    sorted_definition_columns = sorted(definition_columns)

    if flavour_split:
        flavour_output_dir = os.path.join(output_dir, "flavour_splitting")
        if not os.path.exists(flavour_output_dir):
            os.makedirs(flavour_output_dir, exist_ok=True)

        for flavour, correlation_matrix in correlation_data.items():
            flavour_dir = os.path.join(flavour_output_dir, flavour)
            if not os.path.exists(flavour_dir):
                os.makedirs(flavour_dir, exist_ok=True)

            sorted_correlation_matrix = correlation_matrix.loc[
                sorted_filtered_columns + sorted_definition_columns,
                sorted_filtered_columns + sorted_definition_columns,
            ]
            plot_heatmap(
                sorted_correlation_matrix,
                f"{flavour.capitalize()} Flavour Correlation Matrix",
                os.path.join(flavour_dir, f"{flavour}_full_correlation_matrix.png"),
                {**filtered_label_mapping, **definition_label_mapping},
                threshold,
            )

            # Plot the filtered_names x filtered_names submatrix
            filtered_corr_matrix = correlation_matrix.loc[
                sorted_filtered_columns, sorted_filtered_columns
            ]
            plot_heatmap(
                filtered_corr_matrix,
                f"{flavour.capitalize()} Filtered Names Correlation Matrix",
                os.path.join(
                    flavour_dir, f"{flavour}_filtered_names_correlation_matrix.png"
                ),
                filtered_label_mapping,
                0.0,
            )

            # Plot the definition_dict x definition_dict submatrix
            definition_corr_matrix = correlation_matrix.loc[
                sorted_definition_columns, sorted_definition_columns
            ]
            plot_heatmap(
                definition_corr_matrix,
                f"{flavour.capitalize()} Definition Dict Correlation Matrix",
                os.path.join(
                    flavour_dir, f"{flavour}_definition_dict_correlation_matrix.png"
                ),
                definition_label_mapping,
                0.75,
            )

            # Plot the definition_dict x filtered_names submatrix
            def_filtered_corr_matrix = correlation_matrix.loc[
                sorted_definition_columns, sorted_filtered_columns
            ]
            plot_heatmap(
                def_filtered_corr_matrix,
                f"{flavour.capitalize()} Definition Dict x Filtered Names Correlation Matrix",
                os.path.join(
                    flavour_dir,
                    f"{flavour}_definition_dict_filtered_names_correlation_matrix.png",
                ),
                {**definition_label_mapping, **filtered_label_mapping},
                0.12,
            )

            # Plot the filtered_names x definition_dict submatrix
            filtered_def_corr_matrix = correlation_matrix.loc[
                sorted_filtered_columns, sorted_definition_columns
            ]
            plot_heatmap(
                filtered_def_corr_matrix,
                f"{flavour.capitalize()} Filtered Names x Definition Dict Correlation Matrix",
                os.path.join(
                    flavour_dir,
                    f"{flavour}_filtered_names_definition_dict_correlation_matrix.png",
                ),
                {**filtered_label_mapping, **definition_label_mapping},
                threshold,
            )
    elif split_region_B:
        for region, matrix in correlation_data.items():
            region_dir = os.path.join(output_dir, region)
            if not os.path.exists(region_dir):
                os.makedirs(region_dir, exist_ok=True)

            sorted_correlation_matrix = matrix.loc[
                sorted_filtered_columns + sorted_definition_columns,
                sorted_filtered_columns + sorted_definition_columns,
            ]
            plot_heatmap(
                sorted_correlation_matrix,
                f"Correlation Matrix for {region} DeepFlavB",
                os.path.join(region_dir, f"correlation_matrix_{region}.png"),
                {**filtered_label_mapping, **definition_label_mapping},
                threshold,
            )

            # Plot the filtered_names x filtered_names submatrix
            filtered_corr_matrix = matrix.loc[
                sorted_filtered_columns, sorted_filtered_columns
            ]
            plot_heatmap(
                filtered_corr_matrix,
                f"Filtered Names Correlation Matrix for {region} DeepFlavB",
                os.path.join(
                    region_dir, f"filtered_names_correlation_matrix_{region}.png"
                ),
                filtered_label_mapping,
                0.0,
            )

            # Plot the definition_dict x definition_dict submatrix
            definition_corr_matrix = matrix.loc[
                sorted_definition_columns, sorted_definition_columns
            ]
            plot_heatmap(
                definition_corr_matrix,
                f"Definition Dict Correlation Matrix for {region} DeepFlavB",
                os.path.join(
                    region_dir, f"definition_dict_correlation_matrix_{region}.png"
                ),
                definition_label_mapping,
                0.75,
            )

            # Plot the definition_dict x filtered_names submatrix
            def_filtered_corr_matrix = matrix.loc[
                sorted_definition_columns, sorted_filtered_columns
            ]
            plot_heatmap(
                def_filtered_corr_matrix,
                f"Definition Dict x Filtered Names Correlation Matrix for {region} DeepFlavB",
                os.path.join(
                    region_dir,
                    f"definition_dict_filtered_names_correlation_matrix_{region}.png",
                ),
                {**definition_label_mapping, **filtered_label_mapping},
                0.12,
            )

            # Plot the filtered_names x definition_dict submatrix
            filtered_def_corr_matrix = matrix.loc[
                sorted_filtered_columns, sorted_definition_columns
            ]
            plot_heatmap(
                filtered_def_corr_matrix,
                f"Filtered Names x Definition Dict Correlation Matrix for {region} DeepFlavB",
                os.path.join(
                    region_dir,
                    f"filtered_names_definition_dict_correlation_matrix_{region}.png",
                ),
                {**filtered_label_mapping, **definition_label_mapping},
                threshold,
            )

    else:
        # Plot the full correlation matrix
        sorted_correlation_matrix = correlation_data.loc[
            sorted_filtered_columns + sorted_definition_columns,
            sorted_filtered_columns + sorted_definition_columns,
        ]
        plot_heatmap(
            sorted_correlation_matrix,
            "Full Correlation Matrix",
            os.path.join(output_dir, "full_correlation_matrix.png"),
            {**filtered_label_mapping, **definition_label_mapping},
            threshold,
        )

        # Plot the filtered_names x filtered_names submatrix
        filtered_corr_matrix = correlation_data.loc[
            sorted_filtered_columns, sorted_filtered_columns
        ]
        plot_heatmap(
            filtered_corr_matrix,
            "Filtered Names Correlation Matrix",
            os.path.join(output_dir, "filtered_names_correlation_matrix.png"),
            filtered_label_mapping,
            0.0,
        )

        # Plot the definition_dict x definition_dict submatrix
        definition_corr_matrix = correlation_data.loc[
            sorted_definition_columns, sorted_definition_columns
        ]
        plot_heatmap(
            definition_corr_matrix,
            "Definition Dict Correlation Matrix",
            os.path.join(output_dir, "definition_dict_correlation_matrix.png"),
            definition_label_mapping,
            0.75,
        )

        # Plot the definition_dict x filtered_names submatrix
        def_filtered_corr_matrix = correlation_data.loc[
            sorted_definition_columns, sorted_filtered_columns
        ]
        plot_heatmap(
            def_filtered_corr_matrix,
            "Definition Dict x Filtered Names Correlation Matrix",
            os.path.join(
                output_dir, "definition_dict_filtered_names_correlation_matrix.png"
            ),
            {**definition_label_mapping, **filtered_label_mapping},
            0.12,
        )

        # Plot the filtered_names x definition_dict submatrix
        filtered_def_corr_matrix = correlation_data.loc[
            sorted_filtered_columns, sorted_definition_columns
        ]
        plot_heatmap(
            filtered_def_corr_matrix,
            "Filtered Names x Definition Dict Correlation Matrix",
            os.path.join(
                output_dir, "filtered_names_definition_dict_correlation_matrix.png"
            ),
            {**filtered_label_mapping, **definition_label_mapping},
            threshold,
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
        "--specify_MC",
        action="store_true",
        help="Specify MC file to make better data MC comparison",
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

    output_plot_path = os.path.join(args.folder, "output_plot.png")
    parquet_dir = os.path.join(args.folder, "parquet_data")

    data = load_data(
        data_file_paths,
        args.folder,
        max_files=args.max_files,
        limit_inputs=args.limit_inputs,
        limit_outputs=args.limit_outputs,
        SMu=args.SMu,
        flavour_split=args.flavour_split,
        specify_MC=args.specify_MC,
    )
    correlation_matrix, filtered_cols, def_cols = compute_correlations(
        data,
        SMu=args.SMu,
        flavour_split=args.flavour_split,
        split_region_B=args.split_region_b,
    )

    # Create a separate subfolder for plots if specify_MC is active
    if args.specify_MC:
        output_dir = os.path.join(args.folder, "correlation_plots_MC_specified")
    else:
        output_dir = os.path.join(args.folder, "correlation_plots")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    plot_correlations(
        correlation_matrix,
        filtered_cols,
        def_cols,
        output_dir,
        flavour_split=args.flavour_split,
        split_region_B=args.split_region_b,
    )


if __name__ == "__main__":
    main()
