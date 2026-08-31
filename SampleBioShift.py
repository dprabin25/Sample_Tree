# -*- coding: utf-8 -*-
"""
@author: prabindawadi
SampleBioShift.py -- the one script that runs everything: `python
SampleBioShift.py FolderName`.
  0. Pre-flight validation of the input dataset(s) and target.txt --
     see run_preflight_validation() -- runs once, before anything else.
  1. tree_pipeline.R, once per block in methods.txt (single block or
     "# FILE N" batch -- see split_methods_blocks()).
  2. build_observed_shifts() -- merges tree_pipeline.R's clade-based TREND
     output across representations into BioShift_Req/Observed_shifts/.
  3. BioShift.py (disease + healthy passes), reading Observed_shifts/.
Steps 1-2 always run and never touch OpenAI. Step 3 needs a working API
key, checked right before it runs -- a bad key still leaves the tree/
cluster/differential-analysis results from steps 1-2 intact.
Note: target.txt needs real "Y" values for the samples you're
targeting -- if every row is "N", no clade will ever be found.
"""
import os
import re
import sys
import csv
import glob
import shutil
import argparse
import subprocess
import itertools
import pandas as pd
from collections import defaultdict
from datetime import datetime
BIOSHIFT_MODE = "full_with_graphviz"
BIOSHIFT_DIR_NAME = "BioShift_Req"      # source-side folder (BioShift.py, config, graphviz/)
BIOSHIFT_ARCHIVE_NAME = "BioShift"      # name used once archived into <FolderName>/RunN/
# Distance metrics requiring a phylogenetic tree -- must match
# tree_pipeline.R's `phylo_metrics` vector. These write to a different
# Outputs/ path shape (no Normalization/SD0 nesting).
PHYLO_METRICS = {"unifrac", "unifracw", "mpd", "mpdw", "mntd", "mntdw"}
# ------------------------------------------------------------------
# PRE-FLIGHT VALIDATION -- runs once, before step 1 (tree_pipeline.R)
# ever executes. Catches bad input data / target.txt early instead of
# burning a full tree_pipeline.R run before failing.
#
# Checks:
#   1. input file exists
#   2. first column contains non-empty sample IDs
#   3. sample IDs are unique
#   4. all feature columns are numeric
#   5. no unexpected NA values were created
#   6. target.txt contains Sample and Target
#   7. target sample IDs are unique
#   8. all target sample IDs exist in the input dataset
#   9. Target values are only accepted values (Y/N)
# methods.txt combination validity is intentionally NOT checked here --
# that stays tree_pipeline.R's job.
# ------------------------------------------------------------------
ACCEPTED_TARGET_VALUES = {"Y", "N"}
def _find_project_file(base_dir, filename):
    """Locates a file referenced by methods.txt/target.txt on disk.
    Checked in order: base_dir itself, then a few common data
    subfolders, then a recursive search under base_dir (skipping
    Outputs/, the BioShift_Req folder, and Run* archive folders)."""
    candidates = [
        os.path.join(base_dir, filename),
        os.path.join(base_dir, "Data", filename),
        os.path.join(base_dir, "Input", filename),
        os.path.join(base_dir, "Inputs", filename),
    ]
    for c in candidates:
        if os.path.isfile(c):
            return c
    skip = {"Outputs", BIOSHIFT_DIR_NAME, BIOSHIFT_ARCHIVE_NAME}
    for root, dirs, files in os.walk(base_dir):
        dirs[:] = [d for d in dirs if d not in skip and not (d.startswith("Run") and d[3:].isdigit())]
        if filename in files:
            return os.path.join(root, filename)
    return None
def _validate_input_csv(base_dir, path, filename, errors):
    """Checks 1-5 against one input file. Returns the loaded DataFrame
    (so target.txt validation can cross-check sample IDs against it),
    or None if it couldn't be loaded / is unusable for later checks."""
    # 1. input file exists
    if path is None:
        errors.append(
            f"[{filename}] input file not found (looked in {base_dir}, "
            f"Data/, Input/, Inputs/, and subfolders)."
        )
        return None
    try:
        df = pd.read_csv(path)
    except Exception as e:
        errors.append(f"[{filename}] could not be read as CSV: {e}")
        return None
    if df.shape[1] < 2:
        errors.append(
            f"[{filename}] has fewer than 2 columns -- need a sample-ID "
            f"column plus at least one feature column."
        )
        return None
    id_col = df.columns[0]
    ids = df[id_col]
    # 2. first column contains non-empty sample IDs
    blank_mask = ids.isna() | (ids.astype(str).str.strip() == "")
    if blank_mask.any():
        bad_rows = df.index[blank_mask].tolist()
        shown = bad_rows[:10]
        errors.append(
            f"[{filename}] {int(blank_mask.sum())} row(s) have an empty sample ID "
            f"in column '{id_col}' (rows: {shown}{'...' if len(bad_rows) > 10 else ''})."
        )
    # 3. sample IDs are unique
    dupe_mask = ids.duplicated(keep=False) & ~blank_mask
    if dupe_mask.any():
        dupe_vals = sorted(set(ids[dupe_mask].astype(str)))
        shown = dupe_vals[:10]
        errors.append(
            f"[{filename}] duplicate sample IDs in column '{id_col}': "
            f"{shown}{'...' if len(dupe_vals) > 10 else ''}."
        )
    # 4. all feature columns are numeric
    feature_cols = list(df.columns[1:])
    coerced = df[feature_cols].apply(lambda c: pd.to_numeric(c, errors="coerce"))
    for col in feature_cols:
        newly_na = df[col].notna() & coerced[col].isna()
        if newly_na.any():
            bad_vals = sorted(set(df.loc[newly_na, col].astype(str)))
            shown = bad_vals[:5]
            errors.append(
                f"[{filename}] feature column '{col}' has non-numeric value(s): "
                f"{shown}{'...' if len(bad_vals) > 5 else ''}."
            )
    # 5. no unexpected NA values were created (by coercing feature
    #    columns to numeric -- distinct from NAs already present in the
    #    original file, which are left alone)
    original_na = int(df[feature_cols].isna().sum().sum())
    coerced_na = int(coerced.isna().sum().sum())
    if coerced_na > original_na:
        errors.append(
            f"[{filename}] converting feature columns to numeric produced "
            f"{coerced_na - original_na} unexpected NA value(s) beyond the "
            f"{original_na} already present in the file -- see the non-numeric "
            f"value(s) reported above."
        )
    return df
def _read_target_txt(target_path):
    """target.txt may be comma- or tab-delimited -- sniff the delimiter
    rather than assuming one."""
    with open(target_path, encoding="utf-8") as f:
        sample = f.read(4096)
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        sep = dialect.delimiter
    except csv.Error:
        sep = "\t" if "\t" in sample else ","
    tdf = pd.read_csv(target_path, sep=sep)
    tdf.columns = [c.strip() for c in tdf.columns]
    return tdf
def _validate_target_txt(base_dir, id_cols_by_file, errors):
    """Checks 6-9 against target.txt. id_cols_by_file is
    {filename: pd.Series of sample IDs} for every successfully-loaded
    input file, used to cross-check target.txt's Sample column."""
    target_path = _find_project_file(base_dir, "target.txt")
    if target_path is None:
        errors.append(f"target.txt not found (looked in {base_dir}, Data/, Input/, Inputs/, and subfolders).")
        return
    try:
        tdf = _read_target_txt(target_path)
    except Exception as e:
        errors.append(f"target.txt could not be read: {e}")
        return
    # 6. target.txt contains Sample and Target
    if "Sample" not in tdf.columns or "Target" not in tdf.columns:
        errors.append(
            f"target.txt must contain 'Sample' and 'Target' columns; "
            f"found {list(tdf.columns)}."
        )
        return
    sample_ids = tdf["Sample"].astype(str).str.strip()
    # 7. target sample IDs are unique
    dupe_mask = sample_ids.duplicated(keep=False)
    if dupe_mask.any():
        dupe_vals = sorted(set(sample_ids[dupe_mask]))
        shown = dupe_vals[:10]
        errors.append(
            f"target.txt has duplicate Sample IDs: {shown}{'...' if len(dupe_vals) > 10 else ''}."
        )
    # 8. all target sample IDs exist in the input dataset
    if id_cols_by_file:
        known_ids = set()
        for ids in id_cols_by_file.values():
            known_ids |= set(ids.astype(str).str.strip())
        missing = sorted(set(sample_ids) - known_ids)
        if missing:
            shown = missing[:10]
            errors.append(
                f"target.txt Sample ID(s) not found in any input dataset "
                f"({list(id_cols_by_file.keys())}): {shown}{'...' if len(missing) > 10 else ''}."
            )
    # 9. Target values are only accepted values (Y/N)
    target_vals = tdf["Target"].astype(str).str.strip().str.upper()
    bad_targets = sorted(set(target_vals) - ACCEPTED_TARGET_VALUES)
    if bad_targets:
        shown = bad_targets[:10]
        errors.append(
            f"target.txt 'Target' column has value(s) outside "
            f"{sorted(ACCEPTED_TARGET_VALUES)}: {shown}{'...' if len(bad_targets) > 10 else ''}."
        )
def run_preflight_validation(base_dir, blocks, log_fh):
    """Runs checks 1-9 once, before tree_pipeline.R touches anything.
    Collects every problem found instead of stopping at the first, so a
    single re-run can fix everything at once. Aborts (sys.exit) if any
    check fails; otherwise prints/logs confirmation and returns.
    Intentionally does NOT validate methods.txt combination validity --
    that stays tree_pipeline.R's job."""
    header = "\n" + "=" * 80 + "\n[RUN] Pre-flight validation\n" + "=" * 80 + "\n"
    print(header, end="")
    log_fh.write(header)
    errors = []
    # One Filename can repeat across blocks (same file, different
    # Normalization/Distance_Metric) -- validate each distinct file once.
    filenames = []
    for block in blocks:
        fields = parse_methods_fields(block)
        fn = fields.get("Filename")
        if fn and fn not in filenames:
            filenames.append(fn)
    if not filenames:
        errors.append("methods.txt has no 'Filename = ...' entries -- nothing to validate.")
    id_cols_by_file = {}
    for fn in filenames:
        path = _find_project_file(base_dir, fn)
        df = _validate_input_csv(base_dir, path, fn, errors)
        if df is not None:
            id_cols_by_file[fn] = df[df.columns[0]]
    _validate_target_txt(base_dir, id_cols_by_file, errors)
    if errors:
        msg = "\n[ABORT] Pre-flight validation failed:\n" + "\n".join(f"  - {e}" for e in errors) + "\n"
        print(msg)
        log_fh.write(msg)
        sys.exit(msg)
    ok = "Pre-flight validation OK -- input file(s) and target.txt look valid.\n\n"
    print(ok, end="")
    log_fh.write(ok)
# ------------------------------------------------------------------
# build_observed_shifts: reads tree_pipeline.R's native clade-based
# TREND output and combines it across representations into
# BioShift_Req/Observed_shifts/.
# ------------------------------------------------------------------
_FIELD_LINE_RE = re.compile(r"^\s*([^=\n#]+?)\s*=\s*(.*)$")
def parse_methods_fields(block_text):
    """Parses one methods.txt block's `Key = value` lines into a dict."""
    fields = {}
    for ln in block_text.splitlines():
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        m = _FIELD_LINE_RE.match(s)
        if not m:
            continue
        key = m.group(1).strip()
        val = re.split(r"\s+#", m.group(2), maxsplit=1)[0].strip()
        fields[key] = val
    return fields
def rep_specs_from_blocks(blocks):
    """One (file_base, Normalization, Distance_Metric) per methods.txt
    block with Differential analysis != No -- these are the jobs that
    can have written trend_outputs/ for build_observed_shifts to read."""
    specs = []
    for block in blocks:
        fields = parse_methods_fields(block)
        filename = fields.get("Filename")
        dist = fields.get("Distance_Metric")
        if not filename or not dist:
            continue
        file_base = os.path.splitext(filename)[0]
        norm = fields.get("Normalization") or "None"
        diff_raw = fields.get("Differential analysis", "Yes").strip().lower()
        if diff_raw in ("no", "n", "false", "0"):
            print(f"{file_base}: Differential analysis = No -- skipping (no trend_outputs to read).")
            continue
        specs.append((file_base, norm, dist))
    return specs
def _build_rep_tags(rep_names):
    """Short, unambiguous tags for combined-file names. Not
    rep[0].upper(): two names can share a first letter and collide.
    First three letters, lengthened only on collision."""
    rep_names = list(rep_names)
    n = 3
    tags = {rep: rep[:n] for rep in rep_names}
    while len(set(tags.values())) < len(tags):
        n += 1
        tags = {rep: rep[:n] for rep in rep_names}
    return tags
def merge_input_csvs(group_dir):
    """Merges every Input_*.csv in one clade's group_ folder into a
    single Element/Observed-Shift table: identical value everywhere ->
    keep it, disagreement -> 0."""
    files = [f for f in glob.glob(os.path.join(group_dir, "Input_*.csv"))
             if not os.path.basename(f).startswith("Input_MERGED_")]
    dfs = []
    for f in files:
        try:
            df = pd.read_csv(f)
        except Exception:
            continue
        if "Element" in df.columns and "Observed Shift" in df.columns:
            dfs.append(df[["Element", "Observed Shift"]])
    if not dfs:
        return None
    merged = pd.concat(dfs, ignore_index=True)
    merged = merged.groupby("Element")["Observed Shift"].apply(
        lambda x: x.iloc[0] if len(set(x)) == 1 else 0
    ).reset_index()
    return merged
def build_observed_shifts(base_dir, bioshift_dir, rep_specs):
    out_dir = os.path.join(bioshift_dir, "Observed_shifts")
    os.makedirs(out_dir, exist_ok=True)
    for f in glob.glob(os.path.join(out_dir, "*.csv")):
        os.remove(f)
    if not rep_specs:
        print("\nNo representations with Differential analysis = Yes in methods.txt -- nothing to combine.")
        return
    rep_tag = _build_rep_tags([r[0] for r in rep_specs])
    # Collect each representation's clade(s) in memory -- only the final
    # combined file(s) get written to Observed_shifts/.
    type_to_dfs = {}   # {rep: [(node, df), ...]}
    for rep, norm, dist in rep_specs:
        file_base = rep
        if dist.strip().lower() in PHYLO_METRICS:
            trend_root = os.path.join(base_dir, "Outputs", file_base, dist, "trend_outputs")
        else:
            trend_root = os.path.join(base_dir, "Outputs", file_base, norm, dist, "SD0", "trend_outputs")
        if not os.path.isdir(trend_root):
            print(f"{rep}: no {trend_root} yet -- has tree_pipeline.R run for this job "
                  f"with a real Y-target in target.txt? Treated as no clade.")
            continue
        group_dirs = sorted(glob.glob(os.path.join(trend_root, "group_*")))
        if not group_dirs:
            print(f"{rep}: no clade found (target.txt Y samples didn't form one -- "
                  f"0 group_* folders in {trend_root}).")
            continue
        found = []
        for g_dir in group_dirs:
            merged = merge_input_csvs(g_dir)
            if merged is None or merged.empty:
                continue
            m = re.search(r"node(\d+)", os.path.basename(g_dir))
            node = m.group(1) if m else os.path.basename(g_dir)
            found.append((node, merged))
            print(f"{rep}: clade node{node} found ({len(merged)} element(s))")
        if found:
            type_to_dfs[rep] = found
    reps_with_clades = [r for r in type_to_dfs if type_to_dfs[r]]
    if not reps_with_clades:
        print("\nNo clades found in any representation -- nothing to combine. "
              "(target.txt needs real Y-flagged samples that actually form a clade "
              "meeting min_targeted/max_others.)")
        return
    # Combine across representation types. If every representation found
    # exactly one clade, this is a single 3-way merge. If a representation
    # found more than one clade, each of ITS clades pairs separately with
    # the other representations' single clade(s), so multiple combos come
    # out as distinct biological groupings, not noise to collapse.
    n_types = len(reps_with_clades)
    k = 3 if n_types >= 3 else (2 if n_types == 2 else 1)
    def resolve(row):
        vals = [v for v in row if pd.notna(v)]
        if not vals:
            return None
        return vals[0] if len(set(vals)) == 1 else 0
    combos = []
    for tset in itertools.combinations(reps_with_clades, k):
        filegrid = [type_to_dfs[t] for t in tset]
        for choice in itertools.product(*filegrid):
            combos.append((tset, choice))
    combo_counter = 0
    for tset, choice in combos:
        combo_counter += 1
        dfs, tags = [], []
        for rep, (node, df) in zip(tset, choice):
            dfs.append(df.rename(columns={"Observed Shift": f"{rep}_{node}"}))
            tags.append(f"{rep_tag[rep]}{node}")
        merged = dfs[0]
        for df in dfs[1:]:
            merged = merged.merge(df, on="Element", how="outer")
        merged["Observed Shift"] = merged.drop(columns=["Element"]).apply(resolve, axis=1)
        # One combo overall -> plain name; multiple combos -> numbered.
        out_name = ("Combined_Observed_Shifts.csv" if len(combos) == 1
                    else f"Combined_{combo_counter:02d}_" + "_".join(tags) + ".csv")
        merged[["Element", "Observed Shift"]].to_csv(os.path.join(out_dir, out_name), index=False)
        print(f"Combined -> {out_name}")
    print(f"\n{combo_counter} combined file(s) -> {out_dir}")
# ------------------------------------------------------------------
# API key validation -- checked right before step 3, never blocks 1-2.
# ------------------------------------------------------------------
def validate_api_key(bioshift_dir):
    config_path = os.path.join(bioshift_dir, "config_bioshift.txt")
    cfg = {}
    if os.path.exists(config_path):
        for raw in open(config_path, encoding="utf-8"):
            line = raw.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            k, v = line.split("=", 1)
            cfg[k.strip().upper()] = v.strip()
    key = cfg.get("KEY", "").strip() or os.getenv("OPENAI_API_KEY", "").strip()
    if not key or key.upper().startswith("XXXXXXXX"):
        sys.exit(
            f"[ABORT] No API key set in {config_path} (KEY=...) or OPENAI_API_KEY env var. "
            f"tree_pipeline.R + build_observed_shifts already ran successfully -- only the "
            f"BioShift.py (LLM) step is blocked. Fix the key and re-run."
        )
    try:
        import openai
    except Exception:
        sys.exit("[ABORT] 'openai' package not installed. Run: pip install openai")
    print("Validating OpenAI API key before starting BioShift.py...")
    try:
        try:
            from openai import OpenAI
            client = OpenAI(api_key=key)
            client.models.list()
        except ImportError:
            openai.api_key = key
            openai.Model.list()
    except Exception as e:
        sys.exit(
            f"[ABORT] API key check failed -- {e}\n"
            f"tree_pipeline.R + build_observed_shifts already ran successfully -- only the "
            f"BioShift.py (LLM) step is blocked. Fix KEY in {config_path} and try again."
        )
    print("API key OK.\n")
# ------------------------------------------------------------------
# Utilities
# ------------------------------------------------------------------
def run_or_die(cmd, label, log_fh, cwd=None):
    header = (
        "\n" + "=" * 80 +
        f"\n[RUN] {label}\nCommand: " + " ".join(cmd) +
        "\n" + "=" * 80 + "\n"
    )
    print(header, end="")
    log_fh.write(header)
    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        cwd=cwd,
    )
    if result.stdout:
        print(result.stdout, end="")
        log_fh.write(result.stdout)
    if result.returncode != 0:
        msg = f"\n[ERROR] {label} failed with exit code {result.returncode}\n"
        print(msg)
        log_fh.write(msg)
        sys.exit(result.returncode)
_FILE_MARKER_RE = re.compile(r"^\s*#\s*FILE\b", re.IGNORECASE)
_FILENAME_FIELD_RE = re.compile(r"^\s*Filename\s*=\s*(.+?)\s*(?:#.*)?$", re.IGNORECASE)
def split_methods_blocks(text):
    """Splits methods.txt into blocks on "# FILE 1", "# FILE 2", ... lines.
    No markers -> returns [text] (single-job run, as before)."""
    lines = text.splitlines(keepends=True)
    marker_idxs = [i for i, ln in enumerate(lines) if _FILE_MARKER_RE.match(ln)]
    if not marker_idxs:
        return [text]
    blocks = []
    for n, start in enumerate(marker_idxs):
        end = marker_idxs[n + 1] if n + 1 < len(marker_idxs) else len(lines)
        blocks.append("".join(lines[start:end]))
    return blocks
def block_label(block_text, index):
    for ln in block_text.splitlines():
        m = _FILENAME_FIELD_RE.match(ln)
        if m:
            return m.group(1).strip()
    return f"block {index + 1}"
def create_run_folder(base_output_dir):
    os.makedirs(base_output_dir, exist_ok=True)
    runs = [
        d for d in os.listdir(base_output_dir)
        if d.startswith("Run") and d[3:].isdigit()
    ]
    run_id = max([int(d[3:]) for d in runs], default=0) + 1
    run_folder = os.path.join(base_output_dir, f"Run{run_id}")
    os.makedirs(run_folder, exist_ok=True)
    return run_folder
# ------------------------------------------------------------------
# ARCHIVING -- always sweep these three output-only locations into
# the run folder (they never hold anything else, so this is safe
# unconditionally, every run).
# ------------------------------------------------------------------
ARCHIVE_LOCATIONS = [
    ("Outputs", "Outputs"),
    (os.path.join(BIOSHIFT_DIR_NAME, "Observed_shifts"), os.path.join(BIOSHIFT_ARCHIVE_NAME, "Observed_shifts")),
    (os.path.join(BIOSHIFT_DIR_NAME, "BioShiftOutputs"), os.path.join(BIOSHIFT_ARCHIVE_NAME, "BioShiftOutputs")),
]
def archive_run_outputs(base_dir, run_folder, log_fh):
    print("\n[INFO] Archiving this run's outputs\n")
    log_fh.write("\n[INFO] Archiving this run's outputs\n")
    for rel_src, rel_dst in ARCHIVE_LOCATIONS:
        src = os.path.join(base_dir, rel_src)
        if not os.path.isdir(src) or not os.listdir(src):
            msg = f"[INFO] Nothing to archive at {rel_src} (empty or missing)\n"
            print(msg, end="")
            log_fh.write(msg)
            continue
        dst = os.path.join(run_folder, rel_dst)
        try:
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.move(src, dst)
            msg = f"MOVED: {rel_src} -> {os.path.relpath(dst, base_dir)}\n"
            print(msg, end="")
            log_fh.write(msg)
        except Exception as e:
            err = f"[WARN] Could not move {rel_src}: {e}\n"
            print(err, end="")
            log_fh.write(err)
# ------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("folder_name", help="Run folder (e.g., Run_2026_08_26)")
    parser.add_argument("--samples", nargs="*", default=None)
    parser.add_argument("--python-exe", default="python")
    parser.add_argument("--rscript-exe", default="Rscript")
    args = parser.parse_args()
    base_dir = os.path.dirname(os.path.abspath(__file__))       # SampleBioShift/ root
    bioshift_dir = os.path.join(base_dir, BIOSHIFT_DIR_NAME)     # BioShift_Req/ subfolder
    base_output_dir = os.path.join(base_dir, args.folder_name)
    print("Working directory:", base_dir)
    run_folder = create_run_folder(base_output_dir)
    log_path = os.path.join(run_folder, "log.txt")
    with open(log_path, "a", encoding="utf-8") as log_fh:
        log_fh.write("#" * 80 + "\n")
        log_fh.write(f"Run started: {datetime.now().isoformat()}\n")
        log_fh.write(f"Run folder: {run_folder}\n")
        log_fh.write("#" * 80 + "\n")
        print("\n=== PIPELINE CONFIGURATION ===")
        print("tree_pipeline.R (whatever methods.txt currently holds) -> build_observed_shifts() -> BioShift.py")
        print("BioShift mode:", BIOSHIFT_MODE)
        try:
            _run_pipeline_steps(base_dir, bioshift_dir, methods_path_args=args, log_fh=log_fh)
        finally:
            # Always archive, even if a step above failed or exited early.
            archive_run_outputs(base_dir, run_folder, log_fh)
    print("\nPipeline complete (see log for whether every step succeeded).")
    print("All generated outputs moved to:", run_folder)
def _run_pipeline_steps(base_dir, bioshift_dir, methods_path_args, log_fh):
        args = methods_path_args
        # ---- STEP 1: tree_pipeline.R, once per methods.txt block ----
        methods_path = os.path.join(base_dir, "methods.txt")
        with open(methods_path, encoding="utf-8") as f:
            original_methods_text = f.read()
        blocks = split_methods_blocks(original_methods_text)
        print(f"methods.txt: {len(blocks)} combo(s) to run\n")
        # ---- STEP 0: pre-flight validation, before analysis begins ----
        run_preflight_validation(base_dir, blocks, log_fh)
        try:
            for i, block in enumerate(blocks):
                label = block_label(block, i)
                with open(methods_path, "w", encoding="utf-8") as f:
                    f.write(block)
                run_or_die([args.rscript_exe, os.path.join(base_dir, "tree_pipeline.R")],
                           f"tree_pipeline.R ({label})", log_fh)
        finally:
            with open(methods_path, "w", encoding="utf-8") as f:
                f.write(original_methods_text)
        # ---- STEP 2: build_observed_shifts ----
        print("\n" + "=" * 80 + "\n[RUN] build_observed_shifts\n" + "=" * 80 + "\n")
        log_fh.write("\n" + "=" * 80 + "\n[RUN] build_observed_shifts\n" + "=" * 80 + "\n")
        rep_specs = rep_specs_from_blocks(blocks)
        build_observed_shifts(base_dir, bioshift_dir, rep_specs)
        # ---- STEP 3: BioShift.py (LLM), gated on the API key ----
        validate_api_key(bioshift_dir)
        for ctx in ["disease", "healthy"]:
            if args.samples:
                for s in args.samples:
                    run_or_die(
                        [args.python_exe, "BioShift.py",
                         "--context", ctx,
                         "--mode", BIOSHIFT_MODE,
                         "--sample", s],
                        f"BioShift {ctx} {s}", log_fh, cwd=bioshift_dir
                    )
            else:
                run_or_die(
                    [args.python_exe, "BioShift.py",
                     "--context", ctx,
                     "--mode", BIOSHIFT_MODE],
                    f"BioShift {ctx}", log_fh, cwd=bioshift_dir
                )
if __name__ == "__main__":
    main()
