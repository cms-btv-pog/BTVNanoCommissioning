import os
import pandas as pd
import sys

# Function to filter files based on criteria and extract information
def process_csv_files(directory, category, era):
    selected_files = []

    for filename in os.listdir(directory):
        if filename.endswith(".csv"):
            if category == "data" and era in filename:
                selected_files.append(filename)
            elif category == "mc" and era in filename.lower():
                selected_files.append(filename)

    return selected_files

# Function to extract information from the 3rd column of CSV files and write to a txt file
def extract_third_column_and_write_txt(files, directory, output_filename):
    with open(output_filename, 'w') as output_file:
        for filename in files:
            file_path = os.path.join(directory, filename)
            df = pd.read_csv(file_path)
            print("Adding following files")
            if len(df.columns) >= 3:  # Ensure at least 3 columns exist
                third_column_data = df.iloc[:, 2].astype(str).str.cat(sep='\n')  # Assuming 0-based indexing
                print(third_column_data)
                output_file.write(f"{third_column_data}\n")

# Check if the correct number of command-line arguments are provided
if len(sys.argv) != 3:
    print("Usage: python script.py <data/mc> <2022/2023/Summer22/Summer22EE>")
    sys.exit(1)

# Set your directory path
directory_path = '/eos/user/b/btvweb/www/BTVNanoProduction/csvoutputs'

# Process files based on command-line arguments
category_arg = sys.argv[1]
era_arg = sys.argv[2]

selected_files = process_csv_files(directory_path, category_arg, era_arg)

if category_arg == "data":
    output_filename = f'{category_arg}_{era_arg}.txt'
else:
    output_filename = f'{category_arg}_{era_arg}.txt'

extract_third_column_and_write_txt(selected_files, directory_path, output_filename)
print(f"New txt files created containing DAS path names: {output_filename}")
