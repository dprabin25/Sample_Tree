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
        if d.startswith("Run") and d[3:].isdigit()
    ]
    run_id = max([int(d[3:]) for d in runs], default=0) + 1
    run_folder = os.path.join(base_output_dir, f"Run{run_id}")
    os.makedirs(run_folder, exist_ok=True)
    return run_folder


# ------------------------------------------------------------------
# SNAPSHOT LOGIC
# ------------------------------------------------------------------

def snapshot_state(base_dir):
    """
    Snapshot ALL filesystem paths (files + folders) recursively.
    """
    paths = set()
    for root, dirs, files in os.walk(base_dir):
        for d in dirs:
            paths.add(os.path.join(root, d))
        for f in files:
            paths.add(os.path.join(root, f))
    return paths


def move_only_new_outputs(base_dir, run_folder, log_fh, before_snapshot):
    """
    Move ALL newly generated outputs (files AND folders) created during the run,
    while preserving workspace integrity.
    """

    after_snapshot = snapshot_state(base_dir)
    new_paths = after_snapshot - before_snapshot

    print("\n[INFO] Moving ALL newly generated outputs and folders\n")
    log_fh.write("\n[INFO] Moving ALL newly generated outputs and folders\n")

    # Determine TOP-LEVEL new items only
    top_level_items = set()
    for p in new_paths:
        rel = os.path.relpath(p, base_dir)
        top = rel.split(os.sep)[0]
        top_level_items.add(top)

    patient_root = os.path.dirname(run_folder)

    for item in sorted(top_level_items):
        src = os.path.join(base_dir, item)
        dst = os.path.join(run_folder, item)

        # --- SAFETY RULES (CRITICAL) ---

        # 1. Never move the run folder itself
        if os.path.abspath(src) == os.path.abspath(run_folder):
            continue

        # 2. Never move the patient workspace root
        if os.path.abspath(src) == os.path.abspath(patient_root):
            continue

        # 3. Never move pipeline scripts
        if item.endswith((".py", ".R")):
            continue

        # --------------------------------

        try:
            shutil.move(src, dst)   # CUT + MOVE (folders included)
            msg = f"MOVED: {item}\n"
            print(msg, end="")
            log_fh.write(msg)
        except Exception as e:
            err = f"[WARN] Could not move {item}: {e}\n"
            print(err, end="")
            log_fh.write(err)


# ------------------------------------------------------------------
# MAIN
# ------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("folder_name", help="Patient folder (e.g., DemoPrabin)")
    parser.add_argument("--samples", nargs="*", default=None)
    parser.add_argument("--python-exe", default="python")
    parser.add_argument("--rscript-exe", default="Rscript")
    args = parser.parse_args()

    base_dir = os.path.dirname(os.path.abspath(__file__))
    base_output_dir = os.path.join(base_dir, args.folder_name)

    print("Working directory:", base_dir)

    run_folder = create_run_folder(base_output_dir)
    log_path = os.path.join(run_folder, "log.txt")

    # SNAPSHOT BEFORE RUN
    before_snapshot = snapshot_state(base_dir)

    with open(log_path, "a", encoding="utf-8") as log_fh:
        log_fh.write("#" * 80 + "\n")
        log_fh.write(f"Run started: {datetime.now().isoformat()}\n")
        log_fh.write(f"Run folder: {run_folder}\n")
        log_fh.write("#" * 80 + "\n")

        print("\n=== PIPELINE CONFIGURATION ===")
        print("SampleTree → ObservedShifts → BioShift")
        print("BioShift mode:", BIOSHIFT_MODE)

        # ---------------- PIPELINE ----------------
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
        # ---------------- PIPELINE ----------------

        # ARCHIVE OUTPUTS
        move_only_new_outputs(
            base_dir,
            run_folder,
            log_fh,
            before_snapshot
        )

    print("\nPipeline complete.")
    print("All generated outputs moved to:", run_folder)


if __name__ == "__main__":
    main()
