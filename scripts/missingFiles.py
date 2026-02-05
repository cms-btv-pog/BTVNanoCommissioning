import os, sys
import argparse

parser = argparse.ArgumentParser(description="Check for missing hists_N.coffea files.")
parser.add_argument(
    "--jobName",
    "-j",
    type=str,
    required=True,
    help="Path to the folder containing jobnum_list.txt",
    default="jobs_DY_MC",
)
parser.add_argument(
    "--outputXrootdDir",
    "-o",
    type=str,
    required=True,
    help="Path to the folder containing hists_N.coffea files",
    default="DY_MC",
)
parser.add_argument(
    "--missingfilename", "-f", type=str, help="Name outputfile", default=None
)
parser.add_argument(
    "--updateJDL", "-u", action="store_true", help="Update submit.jdl file"
)  # , action=store_true)
parser.add_argument("--test", "-t", action="store_true", help="test behaviour")

args = parser.parse_args()

# Read the jobnum_list.txt and get all the job numbers
jobFolder = "jobs_" + args.jobName + "/"
jobnum_list_file = jobFolder + "jobnum_list.txt"
if not os.path.isfile(jobnum_list_file):
    print(
        f"The jobnum_list.txt file does not exist at the provided path: {jobnum_list_file}"
    )
    exit(1)

if not os.path.isdir(args.outputXrootdDir):
    print(f"The folder path provided does not exist: {args.outputXrootdDir}")
    exit(1)

with open(jobnum_list_file, "r") as file:
    job_numbers = file.read().splitlines()

# List all the files in the folder
files_in_folder = os.listdir(args.outputXrootdDir)

# Check for each number if the corresponding file exists
missing_files = []
for job_number in job_numbers:
    expected_folder = f"hists_{job_number}"
    fol_i = (
        os.listdir(args.outputXrootdDir + "/" + expected_folder)
        if os.path.isdir(args.outputXrootdDir + "/" + expected_folder)
        else [""]
    )
    expected_file_name = expected_folder + ".coffea"
    if expected_file_name not in files_in_folder and expected_file_name not in fol_i:
        missing_files.append(job_number)

missingfilename = (
    args.missingfilename.replace(".txt", "") + ".txt"
    if args.missingfilename
    else "missing_files_" + args.outputXrootdDir.replace("/", "_") + ".txt"
)
missingfileloc = jobFolder + missingfilename

# Save the list of missing files to missing_files.txt
if len(missing_files) < 1:
    print("All histograms in folder, file not being created")
    exit()
else:
    print("Job numbers missing:", missing_files, jobFolder)
if args.test:
    exit()
with open(missingfileloc, "w") as file:
    for missing_file in missing_files:
        file.write(str(missing_file) + "\n")

print(f"Missing files have been saved to ", missingfileloc)

# Update the jdl file if -u option is on
if args.updateJDL:
    print("am i storing", args.updateJDL)
    jdl_file = "submit.jdl"
    jdl_loc = jobFolder + jdl_file
    os.system("cp " + jdl_loc + " " + jdl_loc.replace(".jdl", "_all.jdl"))
    with open(jdl_loc, "r") as file:
        filedata = file.read()

    filedata = filedata.replace("jobnum_list.txt", missingfilename)

    with open(jdl_loc, "w") as file:
        file.write(filedata)

    print(f"The file {jdl_loc} has been updated with the new missing filename.")
