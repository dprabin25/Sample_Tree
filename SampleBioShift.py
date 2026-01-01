# -*- coding: utf-8 -*-

import os
import subprocess
import sys
import argparse
from datetime import datetime
import shutil

BIOSHIFT_MODE = "full_with_graphviz"

# ------------------------------------------------------------------
# Utilities
# ------------------------------------------------------------------

def run_or_die(cmd, label, log_fh):
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
        text=True
    )

    if result.stdout:
        print(result.stdout, end="")
        log_fh.write(result.stdout)

    if result.returncode != 0:
        msg = f"\n[ERROR] {label} failed with exit code {result.returncode}\n"
        print(msg)
        log_fh.write(msg)
        sys.exit(result.returncode)


def create_run_folder(base_output_dir):
    os.makedirs(base_output_dir, exist_ok=True)
    runs = [
        d for d in os.listdir(base_output_dir)
        if os.path.isdir(os.path.join(base_output_dir, d)) and d.startswith("Run")
    ]
    nums = [int(d[3:]) for d in runs if d[3:].isdigit()]
    run_id = max(nums, default=0) + 1
    run_folder = os.path.join(base_output_dir, f"Run{run_id}")
    os.makedirs(run_folder, exist_ok=True)
    return run_folder


# ------------------------------------------------------------------
# SNAPSHOT LOGIC (FIXED)
# ------------------------------------------------------------------

def snapshot_files(root, exclude_prefixes=()):
    """Track files by size + mtime"""
    exts = {
        ".csv", ".tsv", ".txt", ".log",
        ".dot", ".gv", ".gml",
        ".json", ".yaml", ".yml",
        ".pdf", ".png", ".svg",
        ".xlsx", ".xls", ".html",
        ".nwk"
    }
    snap = {}
    for dp, _, files in os.walk(root):
        if any(os.path.commonpath([dp, ex]) == ex for ex in exclude_prefixes):
            continue
        for f in files:
            if os.path.splitext(f)[1].lower() not in exts:
                continue
            full = os.path.join(dp, f)
            try:
                st = os.stat(full)
                snap[full] = (st.st_size, st.st_mtime)
            except OSError:
                pass
    return snap


def move_and_delete(src, dst, log_fh):
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if os.path.isdir(src):
        shutil.copytree(src, dst)
        shutil.rmtree(src)
    else:
        shutil.copy2(src, dst)
        os.remove(src)
    msg = f"Moved: {src} -> {dst}\n"
    print(msg)
    log_fh.write(msg)


def move_new_outputs(base_dir, run_folder, before, after, log_fh):
    changed = [
        p for p in after
        if p not in before or before[p] != after[p]
    ]

    if not changed:
        msg = "\n[INFO] No new outputs detected — nothing moved.\n"
        print(msg)
        log_fh.write(msg)
        return False

    print("\n[INFO] Moving outputs into run folder\n")
    log_fh.write("\n[INFO] Moving outputs into run folder\n")

    for src in changed:
        rel = os.path.relpath(src, base_dir)
        dst = os.path.join(run_folder, rel)
        move_and_delete(src, dst, log_fh)

    return True


def cleanup_generated_folders(base_dir, protected_dirs):
    for dp, _, _ in os.walk(base_dir, topdown=False):
        if dp in protected_dirs or dp == base_dir:
            continue
        try:
            if not os.listdir(dp):
                os.rmdir(dp)
        except Exception:
            pass


# ------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("folder_name", help="Patient folder (e.g., Pat1)")
    parser.add_argument("--samples", nargs="*", default=None)
    parser.add_argument("--python-exe", default="python")
    parser.add_argument("--rscript-exe", default="Rscript")
    args = parser.parse_args()

    base_dir = os.path.dirname(os.path.abspath(__file__))
    base_output_dir = os.path.join(base_dir, args.folder_name)

    print("Working directory:", base_dir)

    # 🔑 FIX: SNAPSHOT BEFORE RUN FOLDER EXISTS
    before_snapshot = snapshot_files(base_dir, exclude_prefixes=[base_output_dir])

    run_folder = create_run_folder(base_output_dir)
    log_path = os.path.join(run_folder, "log.txt")

    with open(log_path, "a", encoding="utf-8") as log_fh:
        log_fh.write("#" * 80 + "\n")
        log_fh.write(f"Run started: {datetime.now().isoformat()}\n")
        log_fh.write(f"Run folder: {run_folder}\n")
        log_fh.write("#" * 80 + "\n")

        print("\n=== PIPELINE CONFIGURATION ===")
        print("SampleTree → ObservedShifts → BioShift")
        print("BioShift mode:", BIOSHIFT_MODE)

        run_or_die([args.rscript_exe, "sampletree_simple.R"], "SampleTree", log_fh)
        run_or_die([args.python_exe, "ObservedShifts.py"], "ObservedShifts", log_fh)

        for ctx in ["disease", "healthy"]:
            if args.samples:
                for s in args.samples:
                    run_or_die(
                        [args.python_exe, "BioShift.py",
                         "--context", ctx,
                         "--mode", BIOSHIFT_MODE,
                         "--sample", s],
                        f"BioShift {ctx} {s}", log_fh
                    )
            else:
                run_or_die(
                    [args.python_exe, "BioShift.py",
                     "--context", ctx,
                     "--mode", BIOSHIFT_MODE],
                    f"BioShift {ctx}", log_fh
                )

        after_snapshot = snapshot_files(base_dir, exclude_prefixes=[base_output_dir])

        moved = move_new_outputs(
            base_dir, run_folder,
            before_snapshot, after_snapshot, log_fh
        )

        if moved:
            cleanup_generated_folders(
                base_dir,
                protected_dirs={base_dir, run_folder, base_output_dir}
            )

    print("\n✅ Pipeline complete.")
    print("All outputs are safely stored in:", run_folder)


if __name__ == "__main__":
    main()
