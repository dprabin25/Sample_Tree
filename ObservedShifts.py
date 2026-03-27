#!/usr/bin/env python3
import os
import re
import itertools
import pandas as pd

print("\n BioShift - Automated Merge + Combination Engine")
print("=================================================================\n")

BASE_DIR = os.getcwd()
GROUP_DIR_OUT = os.path.join(BASE_DIR, "Observed_Shifts_by_group")
COMBO_OUT = os.path.join(BASE_DIR, "Observed_Shifts")

os.makedirs(GROUP_DIR_OUT, exist_ok=True)
os.makedirs(COMBO_OUT, exist_ok=True)

print("Working directory:", BASE_DIR)
print("Merged outputs  ->", GROUP_DIR_OUT)
print("Final combos    ->", COMBO_OUT, "\n")


# ================================================================
# PART 1 - MERGE Input_*.csv INSIDE EACH group_*
# ================================================================
def merge_group_inputs(group_dir):
    files = [f for f in os.listdir(group_dir)
             if f.startswith("Input_") and f.endswith(".csv")]

    if not files:
        print(f"No Input_*.csv in: {group_dir}")
        return None

    dfs = []
    for f in files:
        path = os.path.join(group_dir, f)
        try:
            df = pd.read_csv(path)
        except Exception:
            print("Error reading:", path)
            continue

        if "Element" not in df.columns or "Observed Shift" not in df.columns:
            print("Invalid columns:", f)
            continue

        dfs.append(df[["Element", "Observed Shift"]])

    if not dfs:
        print(f"No valid Input files in {group_dir}")
        return None

    merged = pd.concat(dfs, ignore_index=True)

    # Resolve conflicts: identical -> keep, different -> 0
    merged = merged.groupby("Element")["Observed Shift"].apply(
        lambda x: x.iloc[0] if len(set(x)) == 1 else 0
    ).reset_index()

    base = os.path.basename(group_dir)

    # Try to extract node+type from folder name: e.g. group_1_node201_Cell
    m = re.search(r"node(\d+)_([A-Za-z0-9_.-]+)$", base)
    if m:
        node, tname = m.group(1), m.group(2)
        # Compact, standard-style name: 201_Cell.csv
        outname = f"{node}_{tname}.csv"
    else:
        # Fallback: just use the folder name as a short filename
        outname = f"{base}.csv"

    outpath = os.path.join(GROUP_DIR_OUT, outname)
    merged.to_csv(outpath, index=False)

    print(f"MERGED -> {outname}")
    return outname


# ---- Find and merge all group_* folders ----
print("Searching for group_* folders...\n")

group_dirs = []
for root, dirs, files in os.walk(BASE_DIR):
    for d in dirs:
        if d.lower().startswith("group_"):
            group_dirs.append(os.path.join(root, d))

if not group_dirs:
    print("No group_* folders found.")
    exit()

print(f"Found {len(group_dirs)} groups:\n")
for g in group_dirs:
    print("  -", g)

print("\nMerging groups...\n")

merged_files = []
for g in group_dirs:
    out = merge_group_inputs(g)
    if out:
        merged_files.append(out)

if not merged_files:
    print("No merged files created. Aborting combinations.")
    exit()

print("\nMerging completed.")
print("---------------------------------------------------------------\n")


# ================================================================
# PART 2 - HIERARCHICAL COMBINATION LOGIC
# ================================================================
print("Starting Combination Engine\n")


# ---------- A. Parse merged filenames ----------
def parse_merge_filename(fname):
    """
    Supported patterns:
      1) 201_Cell.csv -> node=201, type=Cell
      2) group_1_node201_Cell.csv -> node=201, type=Cell
      3) group_1.csv -> node=1, type='group'
    Returns: (node, typename) as strings
    """
    # 1) Format: 201_Cell.csv
    m = re.match(r"(\d+)_([A-Za-z0-9_.-]+)\.csv$", fname)
    if m:
        return m.group(1), m.group(2)

    # 2) Format: group_1_node201_Cell.csv
    m = re.search(r"node(\d+)_([A-Za-z0-9_.-]+)\.csv$", fname)
    if m:
        return m.group(1), m.group(2)

    # 3) Format: group_1.csv  -> treat as generic "group"
    m = re.match(r"group_(\d+)\.csv$", fname, re.IGNORECASE)
    if m:
        return m.group(1), "group"

    return None, None


# ---------- B. Organize files by type ----------
type_to_files = {}

for fname in merged_files:
    node, tname = parse_merge_filename(fname)

    if node is None:
        print("Skipping unrecognized filename:", fname)
        continue

    if tname not in type_to_files:
        type_to_files[tname] = []

    type_to_files[tname].append((node, fname))

print("Input types detected:")
for t in type_to_files:
    print(f"  - {t} ({len(type_to_files[t])})")
print()


# ---------- C. Determine highest valid combination size ----------
types = list(type_to_files.keys())
N = len(types)

if N >= 3:
    K = 3
elif N == 2:
    K = 2
else:
    K = 1

print(f"Highest valid combination size = {K}-way\n")


# ---------- D. Load CSV helper ----------
def load_shift_csv(path):
    df = pd.read_csv(path)
    return df[["Element", "Observed Shift"]]


# ---------- E. Combine with conflict rule ----------
def combine_list(filelist):
    dfs = []
    for tname, node, path in filelist:
        df = load_shift_csv(path).copy()
        df.rename(columns={"Observed Shift": f"{tname}_{node}"}, inplace=True)
        dfs.append(df)

    merged = dfs[0]
    for df in dfs[1:]:
        merged = merged.merge(df, on="Element", how="outer")

    def resolve(row):
        vals = [v for v in row if pd.notna(v)]
        if len(vals) == 0:
            return None
        if len(set(vals)) == 1:
            return vals[0]
        return 0

    merged["Observed Shift"] = merged.drop(columns=["Element"]).apply(resolve, axis=1)
    return merged[["Element", "Observed Shift"]]


# ---------- F. Perform combinations ----------
combo_counter = 0

type_groups = list(type_to_files.keys())

# Generate all K-way combinations of types
for tset in itertools.combinations(type_groups, K):

    # Build product of file lists
    filegrid = [type_to_files[t] for t in tset]

    for choice in itertools.product(*filegrid):
        combo_counter += 1

        filelist = []
        for tname, (node, fname) in zip(tset, choice):
            fpath = os.path.join(GROUP_DIR_OUT, fname)
            filelist.append((tname, node, fpath))

        # Build short tags like C201, P10, M305 from tname + node
        short_tags = []
        for tname, node, _ in filelist:
            abbrev = tname[0].upper() if tname else "T"  # first letter of type
            short_tags.append(f"{abbrev}{node}")

        # Final compact name: 001_C201_P10_M305.csv
        k = len(tset)  # number of types (just in case you care later)
        outname = f"{combo_counter:03d}_" + "_".join(short_tags) + ".csv"
        outpath = os.path.join(COMBO_OUT, outname)

        combined = combine_list(filelist)
        combined.to_csv(outpath, index=False)

        print("Saved:", outname)


print("\nALL DONE")
print(f"Saved {combo_counter} combined CSVs into:")
print("   ->", COMBO_OUT)
print("=================================================================")
