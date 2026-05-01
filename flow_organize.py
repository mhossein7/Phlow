import os
import shutil
import re
import uuid

def organize_fcs_files(mother_folder, labels:str,num_conds=4):
    """
    mother_folder: path to the directory containing .fcs files
    labels: list of label names (strings)
    """

    # Define ranges for each label
    # 1–4 → label 0, 5–8 → label 1, 9–12 → label 2, 13–16 → label 3
    num_labels = len(labels)
    
    ranges={}
    j = 1
    for i in range(num_labels):
        ranges[i] = range(j,j+num_conds)
        j += num_conds

    # Create subfolders
    subfolders = []
    for label in labels:
        subfolder = os.path.join(mother_folder, label)
        os.makedirs(subfolder, exist_ok=True)
        subfolders.append(subfolder)

    # Regex to capture ending "-#" before .fcs
    pattern = re.compile(r"-([0-9]+)\.fcs$", re.IGNORECASE)

    # Loop over all fcs files
    for filename in os.listdir(mother_folder):
        if filename.lower().endswith(".fcs"):
            match = pattern.search(filename)
            if match:
                num = int(match.group(1))
                
                # Determine which label it belongs to
                for idx, r in ranges.items():
                    if num in r:
                        src = os.path.join(mother_folder, filename)
                        dst = os.path.join(subfolders[idx], filename)
                        shutil.move(src, dst)
                        print(f"Moved {filename} → {labels[idx]}")
                        break
                else:
                    # Number not in 1–16
                    print(f"Skipping {filename}: number {num} outside expected range")


def normalize_filenames(folder):
    """
    Renames all .fcs files in the folder so they end with -1.fcs ... -4.fcs.
    Order is decided by sorting by the original trailing number.
    """

    pattern = re.compile(r"-([0-9]+)\.fcs$", re.IGNORECASE)

    files = []
    for filename in os.listdir(folder):
        if filename.lower().endswith(".fcs"):
            match = pattern.search(filename)
            if match:
                original_num = int(match.group(1))
                files.append((original_num, filename))

    if len(files) == 0:
        print(f"No .fcs files found in: {folder}")
        return

    # Sort by the original number so renaming is consistent
    files.sort(key=lambda x: x[0])

    # Rename sequentially to -1, -2, -3, ...
    for new_num, (_, old_name) in enumerate(files, start=1):
        old_path = os.path.join(folder, old_name)

        # Create new filename by replacing the old "-#.fcs"
        new_name = re.sub(r"-[0-9]+\.fcs$", f"-{new_num}.fcs", old_name, flags=re.IGNORECASE)
        new_path = os.path.join(folder, new_name)

        os.rename(old_path, new_path)
        print(f"{old_name} → {new_name}")

def normalize_in_all_subfolders(mother_folder):
    """
    Walks through all subfolders inside the mother folder and normalizes filenames in each.
    """

    for entry in os.listdir(mother_folder):
        path = os.path.join(mother_folder, entry)
        if os.path.isdir(path):
            print(f"\nNormalizing files in folder: {entry}")
            normalize_filenames(path)


def reverse_fcs_numbering(folder):
    """
    Reverses the numbering of .fcs files in a folder.
    Example: -1 <-> -12, -2 <-> -11, etc.
    """

    pattern = re.compile(r"(.+)-([0-9]+)\.fcs$", re.IGNORECASE)

    files = []
    for fname in os.listdir(folder):
        match = pattern.match(fname)
        if match:
            prefix = match.group(1)
            num = int(match.group(2))
            files.append((fname, prefix, num))

    if len(files) < 2:
        print("Not enough files to reverse.")
        return

    nums = [f[2] for f in files]
    min_n, max_n = min(nums), max(nums)

    # Step 1: rename to temporary unique names to avoid collisions
    temp_map = {}
    for fname, prefix, num in files:
        tmp_name = f"{prefix}-TMP-{uuid.uuid4().hex}.fcs"
        os.rename(
            os.path.join(folder, fname),
            os.path.join(folder, tmp_name)
        )
        temp_map[tmp_name] = (prefix, num)

    # Step 2: rename to final reversed numbers
    for tmp_name, (prefix, old_num) in temp_map.items():
        new_num = min_n + max_n - old_num
        new_name = f"{prefix}-{new_num}.fcs"

        os.rename(
            os.path.join(folder, tmp_name),
            os.path.join(folder, new_name)
        )
        print(f"{prefix}-{old_num}.fcs → {new_name}")
