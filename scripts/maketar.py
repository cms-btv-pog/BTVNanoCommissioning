import os
import tarfile


def make_tarfile(output_filename, source_dir, exclude_dirs=[]):
    with tarfile.open(output_filename, "w:gz") as tar:
        for root, dirs, files in os.walk(source_dir):
            dirs[:] = [
                d for d in dirs if d not in exclude_dirs
            ]  # Exclude specified directories
            for file in files:
                file_path = os.path.join(root, file)
                tar.add(file_path, arcname=os.path.relpath(file_path, source_dir))


base_dir = "."
jobdirs = [d for d in os.listdir(base_dir) if d.startswith("jobs_")]
make_tarfile(
    "BTVNanoCommissioning.tar.gz",
    base_dir,
    exclude_dirs=["jsonpog-integration", "BTVNanoCommissioning.egg-info"] + jobdirs,
)
