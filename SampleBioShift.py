# -*- coding: utf-8 -*-
"""
SampleBioShift.py  -- the ONE script that controls everything in this
folder. Lives at the root (sibling of tree_pipeline.R), not inside
BioShift_Req/, since it's the orchestrator, not part of the LLM-interpretation
package.

Pipeline, one call to `python SampleBioShift.py FolderName`:
  1. tree_pipeline.R, one run per block in methods.txt (single block or
     "# FILE N" batch -- see split_methods_blocks()). methods.txt has no
     CLI config-path argument, it's always read from "methods.txt" next
     to tree_pipeline.R, so this script writes each block into that file
     in turn and restores the original content when done.
  2. build_observed_shifts() (below -- inlined here, no longer a separate
     script). tree_pipeline.R's own run_trend_for_job() writes clade-based
     TREND results natively to
       Outputs/<file_base>/<Norm>/<Dist>/SD0/trend_outputs/group_*/Input_*.csv
     whenever target.txt's Y-targeted samples form a real clade (per
     min_targeted/max_others in sampletree_control.txt; tolerated N
     contamination is folded into that same clade, everyone else is
     "other" -- there's no per-patient grouping anymore, just the one
     global targeted cohort). This function merges each clade's Input_*.csv,
     then combines across representations (whatever Filename
     blocks methods.txt currently has, e.g. Cell/Pro1log10/BacCount --
     derived from methods.txt itself via rep_specs_from_blocks(), never a
     hardcoded mapping) the same way the original package's
     ObservedShifts.py did, and writes the result into
     BioShift_Req/Observed_shifts/.
  3. BioShift.py (disease + healthy passes), which reads whatever's in
     BioShift_Req/Observed_shifts/.

Steps 1-2 run UNCONDITIONALLY every time -- they're the real clustering /
differential-analysis pipeline and don't touch OpenAI at all. Only step 3
needs a working API key, so validate_api_key() is checked right before
step 3, not at the top of main() -- a bad key still lets you get the real
tree/cluster results out of a run, it just skips the LLM interpretation.

Note: target.txt's Target column must have real "Y" values for the
samples you actually want to interrogate -- if every row is "N", no
clade will ever be found and step 2 will report "no clades found in any
representation" every time.
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
# Name of the ROOT-level subfolder holding BioShift.py/config_bioshift.txt/
# graphviz/ (renamed from "BioShift" to "BioShift_Req" to avoid reading as
# a collision with the root-level "Outputs" folder). Single source of
# truth for the source-side name -- change here only.
BIOSHIFT_DIR_NAME = "BioShift_Req"
# Name used for that same folder ONCE ARCHIVED into <FolderName>/RunN/ --
# kept as plain "BioShift" there, since inside a run folder there's no
# "Outputs" to collide with.
BIOSHIFT_ARCHIVE_NAME = "BioShift"

# Distance metrics that require a phylogenetic tree (PhyloTree=<file>.nwk).
# Must match tree_pipeline.R's own `phylo_metrics` vector exactly -- these
# jobs write to a DIFFERENT Outputs/ path shape (no Normalization/SD0
# nesting; see rep_trend_root() below), so build_observed_shifts() has to
# know which ones they are to find their trend_outputs/ at all.
PHYLO_METRICS = {"unifrac", "unifracw", "mpd", "mpdw", "mntd", "mntdw"}


# ------------------------------------------------------------------
# Step 2: build_observed_shifts (inlined -- was its own build_observed_shifts.py;
# folded in here so there's exactly one script that controls everything)
#
# Reads tree_pipeline.R's own NATIVE clade-based TREND output --
# Outputs/<file_base>/<Norm>/<Dist>/SD0/trend_outputs/group_<tag>_<method>_
# node<N>_<label>/Input_<label>_<library>.csv -- written by
# run_trend_for_job() whenever target.txt's Y-targeted samples form a real
# clade (min_targeted/max_others from sampletree_control.txt). There is no
# per-patient concept here anymore: target.txt's Y column defines ONE
# global targeted cohort, and each representation (Cell/Cytokines/
# Microbes) either finds a clade for it or doesn't -- so this collects
# whatever clade(s) each representation found and combines them across
# representations the same way the original package's ObservedShifts.py
# did (merge Input_*.csv within a clade, then every combination across
# representation types, since there's no patient-level alignment to key
# off anymore).
# ------------------------------------------------------------------

_FIELD_LINE_RE = re.compile(r"^\s*([^=\n#]+?)\s*=\s*(.*)$")


def parse_methods_fields(block_text):
    """Parses one methods.txt block's `Key = value` lines into a dict,
    stripping trailing '## comment' / '# comment' text. This is the ONLY
    place representation identity (Filename/Normalization/Distance_Metric)
    comes from -- no hand-maintained Python dict to keep in sync. Rename a
    Filename in methods.txt (Micro.csv -> BacCount.csv, add a 4th block,
    whatever) and this picks it up automatically next run."""
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
    block that actually ran differential analysis. file_base IS the
    representation's name (derived straight from its Filename, e.g.
    "Cell.csv" -> "Cell") -- these are exactly the jobs that can have
    written Outputs/<file_base>/<Norm>/<Dist>/SD0/trend_outputs/. Blocks
    with Differential analysis = No are skipped: tree_pipeline.R stops
    after the plain sample tree for those and never writes trend_outputs/
    at all, so there is nothing here for build_observed_shifts to read."""
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
    """Short, unambiguous tags for combined-file names, derived from the
    representation names themselves. NOT rep[0].upper(): two file_bases
    can share a first letter (e.g. "Cell" and "Cytokines" style names),
    which would silently collide into the same tag ("C201" could mean
    either). First three letters is usually enough; if two reps still
    collide on that, fall back to lengthening just the colliding ones
    until they're distinct.
    """
    rep_names = list(rep_names)
    n = 3
    tags = {rep: rep[:n] for rep in rep_names}
    while len(set(tags.values())) < len(tags):
        n += 1
        tags = {rep: rep[:n] for rep in rep_names}
    return tags


def merge_input_csvs(group_dir):
    """Merges every Input_*.csv (excluding any prior Input_MERGED_*.csv)
    inside one clade's group_ folder into a single Element/Observed-Shift
    table -- same conflict rule throughout this pipeline: identical value
    everywhere -> keep it, disagreement -> 0."""
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

    # ---- collect each representation's clade(s), if any -- kept IN
    # MEMORY only. Observed_shifts/ should end up holding just the final
    # combined file(s) that actually feed BioShift.py, not these
    # per-representation intermediates. ----
    type_to_dfs = {}   # {rep: [(node, df), ...]}
    for rep, norm, dist in rep_specs:
        file_base = rep
        # tree_pipeline.R writes phylo-distance jobs (UniFrac/UniFracW/MPD/
        # MPDw/MNTD/MNTDw) straight to Outputs/<file_base>/<dist>/ -- no
        # Normalization/SD0 nesting, since there's no SD sweep for those
        # (support comes from raw-count bootstrap instead). Non-phylo jobs
        # keep the full Outputs/<file_base>/<norm>/<dist>/SD0/ nesting.
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

    # ---- combine across representation types. When every representation
    # found exactly one clade (the common case), this is one 3-way merge
    # -> exactly one combined file. When a representation found more than
    # one clade (the Y-targeted cohort split into distinct sub-clusters
    # for that data type), each of ITS clades pairs separately with the
    # single clade from the other representations -- e.g. 2 Cell clades x
    # 1 Cytokines clade x 1 Microbes clade -> combo1, combo2, each a
    # genuinely different biological grouping, not noise to collapse. ----
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
        # Single combo overall -> plain name, no pointless "001_" prefix.
        # Multiple combos (a representation split into sub-clusters) ->
        # numbered so each stays distinguishable.
        out_name = ("Combined_Observed_Shifts.csv" if len(combos) == 1
                    else f"Combined_{combo_counter:02d}_" + "_".join(tags) + ".csv")
        merged[["Element", "Observed Shift"]].to_csv(os.path.join(out_dir, out_name), index=False)
        print(f"Combined -> {out_name}")

    print(f"\n{combo_counter} combined file(s) -> {out_dir}")


# ------------------------------------------------------------------
# API key validation -- checked right before step 3 (BioShift.py), not
# at the top of main(), so a missing/wrong key never blocks steps 1-2.
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
    """Splits methods.txt content into one or more blocks on lines like
    "# FILE 1", "# FILE 2", ... Supports BOTH formats:
      - single combo (no "# FILE" markers at all) -> returns [whole text],
        unchanged from the original one-job-per-run behavior.
      - multiple "# FILE N" blocks -> tree_pipeline.R runs once per block,
        automatically, in one SampleBioShift.py call (this is what makes
        "run multiple, then merged for the inputs" work again, matching
        the old package's multi-row methods.txt, just written as stacked
        blocks in the same single methods.txt instead of separate rows)."""
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
# ARCHIVING -- always sweep the three known, exclusively-output
# locations into this run's folder. (Previously this diffed a
# before/after filesystem snapshot to find "new" paths, but that broke
# the moment a prior run's output was left un-archived at the root: the
# next run overwrites the SAME paths in place, so they're never "new"
# relative to before_snapshot and silently never get moved. Since
# Outputs/, BioShift_Req/Observed_shifts/, and BioShift_Req/BioShiftOutputs/
# never hold anything except this pipeline's own output -- no permanent
# assets live in any of them -- it's simplest and most robust to just
# always move all three, every run, full stop.)
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
            # ARCHIVE OUTPUTS -- always, even if a step above failed (e.g. a
            # placeholder/invalid API key aborting step 3) or exited early.
            # Whatever got generated still belongs in this run's folder,
            # not left at the root.
            archive_run_outputs(base_dir, run_folder, log_fh)

    print("\nPipeline complete (see log for whether every step succeeded).")
    print("All generated outputs moved to:", run_folder)


def _run_pipeline_steps(base_dir, bioshift_dir, methods_path_args, log_fh):
        args = methods_path_args
        # ---------------- STEP 1: tree_pipeline.R ----------------
        # methods.txt can hold ONE combo (single Filename/Normalization/...
        # block, no markers -- runs once, as before) OR several stacked
        # "# FILE 1" / "# FILE 2" / ... blocks -- tree_pipeline.R has no
        # native batch support (no CLI config-path arg, always reads
        # "methods.txt"), so each block gets copied over methods.txt in
        # turn and tree_pipeline.R runs once per block, automatically, in
        # this one call. The original multi-block methods.txt is restored
        # afterward so your saved file isn't left holding just the last
        # block.
        methods_path = os.path.join(base_dir, "methods.txt")
        with open(methods_path, encoding="utf-8") as f:
            original_methods_text = f.read()
        blocks = split_methods_blocks(original_methods_text)
        print(f"methods.txt: {len(blocks)} combo(s) to run\n")
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

        # ---------------- STEP 2: build_observed_shifts (inlined, no subprocess) ----------------
        # rep_specs comes straight from the SAME methods.txt blocks Step 1
        # just ran -- not a separate hardcoded mapping -- so it can never
        # drift out of sync with what tree_pipeline.R actually produced.
        print("\n" + "=" * 80 + "\n[RUN] build_observed_shifts\n" + "=" * 80 + "\n")
        log_fh.write("\n" + "=" * 80 + "\n[RUN] build_observed_shifts\n" + "=" * 80 + "\n")
        rep_specs = rep_specs_from_blocks(blocks)
        build_observed_shifts(base_dir, bioshift_dir, rep_specs)

        # ---------------- STEP 3: BioShift.py (LLM) -- gated on the API key ----------------
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
