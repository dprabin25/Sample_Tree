# -*- coding: utf-8 -*-
import os
import subprocess
import sys
import argparse
from datetime import datetime
import shutil
import fnmatch

BIOSHIFT_MODE = "full_with_graphviz"

# Files that should NEVER be moved (including wildcard patterns)
SKIP_FILES = {
    "bioshift.py",
    "observedshifts.py",
    "runallfinal.py",
    "runall.py",
    "sampletree_simple.r",
    "sampletree_control.txt",
    "config.txt",
    "methods.txt",
    "target.txt",
    "*.nwk",  # Wildcard pattern for files ending in .nwk
    "*.csv",  # Wildcard pattern for files ending in .csv
}

def move_and_delete(src, dest, log_fh):
    """Move file or folder and delete the original, logging the move."""
    try:
        if os.path.isdir(src):
            shutil.copytree(src, dest)
        else:
            shutil.copy2(src, dest)
        if os.path.isdir(src):
            shutil.rmtree(src)
        else:
            os.remove(src)
        msg = f"Moved and deleted: {src} -> {dest}\n"
        print(msg)
        log_fh.write(msg)
        log_fh.flush()
    except Exception as e:
        msg = f"Error while moving {src}: {e}\n"
        print(msg)
        log_fh.write(msg)
        log_fh.flush()

def should_skip_file(file_name):
    """Check if a file should be skipped based on wildcard patterns in SKIP_FILES."""
    for pattern in SKIP_FILES:
        if fnmatch.fnmatch(file_name, pattern):
            return True
    return False

def run_or_die(cmd, label, log_fh):
    """Run the command and stop if it fails."""
    header = (
        "\n" + "=" * 80 +
        f"\n[RUN] {label}\nCommand: " + " ".join(cmd) +
        "\n" + "=" * 80 + "\n"
    )
    print(header, end="")
    log_fh.write(header)
    log_fh.flush()

    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True
    )

    if result.stdout:
        print(result.stdout, end="")
        log_fh.write(result.stdout)
        log_fh.flush()

    if result.returncode != 0:
        msg = f"\n[ERROR] {label} failed with exit code {result.returncode}\n"
        print(msg)
        log_fh.write(msg)
        log_fh.flush()
        sys.exit(result.returncode)

def create_run_folder(base_output_dir):
    """Create a numbered Run folder under the patient directory."""
    os.makedirs(base_output_dir, exist_ok=True)
    existing_runs = [
        d for d in os.listdir(base_output_dir)
        if os.path.isdir(os.path.join(base_output_dir, d)) and d.startswith("Run")
    ]
    run_numbers = [int(d[3:]) for d in existing_runs if d[3:].isdigit()]
    next_run = max(run_numbers, default=0) + 1
    run_folder = os.path.join(base_output_dir, f"Run{next_run}")
    os.makedirs(run_folder, exist_ok=True)
    return run_folder

def snapshot_files(root, exclude_prefixes=()):
    """Snapshot of output files to detect what’s new."""
    exts = {".csv", ".tsv", ".txt", ".log", ".dot", ".gv", ".gml",
            ".json", ".yaml", ".yml", ".pdf", ".png", ".svg",
            ".xlsx", ".xls", ".html"}
    files = {}
    for dirpath, dirnames, filenames in os.walk(root):
        skip = any(os.path.commonpath([dirpath, ex]) == ex for ex in exclude_prefixes)
        if skip:
            continue
        for fn in filenames:
            ext = os.path.splitext(fn)[1].lower()
            if ext not in exts or fn.endswith(".py") or fn.endswith(".R"):
                continue
            full = os.path.join(dirpath, fn)
            try:
                stat = os.stat(full)
                files[full] = (stat.st_size, stat.st_mtime)
            except OSError:
                continue
    return files

def move_new_outputs(base_dir, run_folder, before_snapshot, after_snapshot, log_fh):
    """Move new/changed outputs and new folders into run folder, logging the process."""
    new_or_changed = [
        path for path, meta in after_snapshot.items()
        if path not in before_snapshot or before_snapshot[path] != meta
    ]
    before_dirs = {d for d, _, _ in os.walk(base_dir)}
    after_dirs = {d for d, _, _ in os.walk(base_dir)}
    new_dirs = [d for d in after_dirs if d not in before_dirs and d != run_folder]

    if not new_or_changed and not new_dirs:
        msg = "\n[INFO] No new/changed outputs or folders detected.\n"
        print(msg)
        log_fh.write(msg)
        return

    msg = "\n[INFO] Moving new/changed outputs and folders into run folder...\n"
    print(msg)
    log_fh.write(msg)

    for src in new_or_changed:
        rel = os.path.relpath(src, base_dir)
        dst = os.path.join(run_folder, rel)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        move_and_delete(src, dst, log_fh)  # Pass log_fh here to log the move

    for folder in new_dirs:
        name = os.path.basename(folder)
        dst = os.path.join(run_folder, name)
        move_and_delete(folder, dst, log_fh)  # Pass log_fh here to log the move

    for dirpath, dirnames, filenames in os.walk(base_dir, topdown=False):
        if dirpath in [base_dir, run_folder]:
            continue
        if not os.listdir(dirpath):
            try:
                os.rmdir(dirpath)
            except OSError:
                pass

def cleanup_generated_folders(base_dir, exclude_prefixes=()):
    """Delete all folders generated after running the pipeline."""
    for dirpath, dirnames, filenames in os.walk(base_dir, topdown=False):
        if any(os.path.commonpath([dirpath, ex]) == ex for ex in exclude_prefixes):
            continue
        if dirpath != base_dir:  # Don't delete the base folder
            try:
                shutil.rmtree(dirpath)
                print(f"Deleted folder: {dirpath}")
            except Exception as e:
                print(f"Error deleting folder {dirpath}: {e}")

def main():
    parser = argparse.ArgumentParser(
        description=("Collect output folders/files and move them into a run folder."
                     " Afterwards, delete the original files and folders.")
    )
    parser.add_argument("folder_name", help="Patient/output folder name (e.g., Pat1)")
    parser.add_argument("--samples", nargs="*", help="Optional BioShift samples")
    parser.add_argument("--python-exe", default="python", help='Python executable')
    parser.add_argument("--rscript-exe", default="Rscript", help='Rscript executable')
    args = parser.parse_args()

    base_dir = os.path.dirname(os.path.abspath(__file__))
    print("Working directory:", base_dir)

    base_output_dir = os.path.join(base_dir, args.folder_name)
    run_folder = create_run_folder(base_output_dir)
    print(f"\nAll outputs for this run will be saved in:\n  {run_folder}\n")

    before_snapshot = snapshot_files(base_dir, exclude_prefixes=[base_output_dir])
    log_path = os.path.join(run_folder, "log.txt")

    with open(log_path, "a", encoding="utf-8") as log_fh:
        log_fh.write("\n" + "#" * 80 + "\n")
        log_fh.write(f"Run started: {datetime.now().isoformat()}\n")
        log_fh.write(f"Patient folder: {base_output_dir}\n")
        log_fh.write(f"Run folder: {run_folder}\n")
        log_fh.write("#" * 80 + "\n")

        config = [
            "=== PIPELINE CONFIGURATION ===",
            "  Step 1: sampletree_simple.R : WILL RUN",
            "  Step 2: ObservedShifts.py   : WILL RUN",
            "  Step 3: BioShift.py         : WILL RUN [disease, healthy]",
            f"  BioShift mode (fixed)       : {BIOSHIFT_MODE}",
            "==========================================\n"
        ]
        if args.samples:
            config.insert(4, "  BioShift samples            : " + ", ".join(args.samples))
        else:
            config.insert(4, "  BioShift samples            : (none specified; all handled)")
        log_fh.write("\n".join(config) + "\n")
        print("\n".join(config))

        run_or_die([args.rscript_exe, "sampletree_simple.R"], "Sampletree", log_fh)
        run_or_die([args.python_exe, "ObservedShifts.py"], "ObservedShifts", log_fh)

        for context in ["disease", "healthy"]:
            if args.samples:
                for sample in args.samples:
                    run_or_die(
                        [args.python_exe, "BioShift.py", "--context", context,
                         "--mode", BIOSHIFT_MODE, "--sample", sample],
                        f"BioShift ({context}, sample={sample})",
                        log_fh
                    )
            else:
                run_or_die(
                    [args.python_exe, "BioShift.py", "--context", context,
                     "--mode", BIOSHIFT_MODE],
                    f"BioShift ({context}, all samples)",
                    log_fh
                )

        final_msg = "\nAll steps finished successfully for both contexts.\n"
        print(final_msg)
        log_fh.write(final_msg)

        after_snapshot = snapshot_files(base_dir, exclude_prefixes=[base_output_dir])
        move_new_outputs(base_dir, run_folder, before_snapshot, after_snapshot, log_fh)

        # Cleanup generated folders after moving outputs
        cleanup_generated_folders(base_dir, exclude_prefixes=[base_output_dir, run_folder])

    print("\n✅ Pipeline complete. All outputs are in:", run_folder)

if __name__ == "__main__":
    main()
