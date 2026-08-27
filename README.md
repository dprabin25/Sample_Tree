# SampleBioShift

## Description

**SampleBioShift** runs the **SampleTree** workflow and **BioShift** together, from one command.

`SampleBioShift.py` is the orchestrator — the only script you run. It calls `tree_pipeline.R` once per representation in `methods.txt`, builds the combined `Observed_shifts/` input BioShift needs, then calls `BioShift_Req/BioShift.py` for the disease and healthy interpretation passes.

SampleTree groups samples by similarity and looks for **clades of targeted (`Y`) samples** in the resulting tree. There is no per-patient grouping: `target.txt` defines one global targeted cohort (every sample flagged `Y`), and clade detection runs across that whole cohort at once — mirroring the real-world case where you don't know in advance which patient a target sample came from. BioShift then interprets which elements shift within the clades that were found.

### Run modes

`python SampleBioShift.py FolderName` always executes the same three steps, but how far each run actually gets depends on `methods.txt` and `BioShift_Req/config_bioshift.txt` — no separate mode flag is needed:

1. **SampleTree only** — set `Differential analysis = No` on every block in `methods.txt`. `tree_pipeline.R` stops right after the plain `support_tree.nwk/pdf/png`: no clade detection, no highlighted tree, no differential analysis. Use this to just see how samples cluster before deciding anything else.
2. **SampleTree → df** — set `Differential analysis = Yes` (the default), but leave `BioShift_Req/config_bioshift.txt`'s `KEY` unset or invalid. Clade detection, `support_tree_highlighted.jpeg`, and the differential analysis (TREND: limma/MaAsLin2) all run, and `Observed_shifts/` gets built — `SampleBioShift.py` then aborts cleanly at the API-key check, before touching OpenAI. Use this to get the full clustering + differential-analysis result without spending API credits.
3. **SampleTree → df → BioShift** — `Differential analysis = Yes` and a valid `KEY` in `config_bioshift.txt`. The full pipeline runs end to end, including the LLM interpretation passes in `BioShift/BioShiftOutputs/`.

Because steps 1–2 never touch OpenAI, a run always reaches at least mode 2's output even if mode 3 was intended but the API key turns out to be missing or wrong — nothing already computed is lost.

### Important note

The pipeline does not perform input scaling. It is up to the user to decide whether to use scaled or unscaled data (`Normalization` in `methods.txt` handles this per representation).

## Dependencies

### 1. Environment

- Install Anaconda: https://www.anaconda.com/distribution/

```bash
conda create -n Sample_Tree python=3.12.2 r-base=4.5.3 pandas -y -c conda-forge

conda activate Sample_Tree

conda install -c conda-forge r-recommended
```

#### R packages

```
Rscript -e "options(repos=c(CRAN='https://cloud.r-project.org')); install.packages(c('ape','vegan','picante','dplyr','stringr','progress','ggplot2','ggnewscale','gtools','reshape2','ggrepel','uwot','BiocManager'))"
```

Install the Bioconductor packages

```
Rscript -e "BiocManager::install(c('phyloseq','ggtree','limma','Maaslin2'), ask=FALSE, update=FALSE)"
```

Verify load packages

```
Rscript -e "pkgs <- c('optparse','ape','vegan','picante','phangorn','progress','ggplot2','dplyr','readr','stringr','tibble','reshape2','data.table','tidyr','pbapply','matrixStats','Hmisc','quantreg','lme4','lmerTest','Rcpp','phyloseq','ggtree','treeio','limma','Maaslin2'); print(data.frame(Package=pkgs, Loads=sapply(pkgs, function(p) require(p, character.only=TRUE, quietly=TRUE))))"
```

#### Python packages

Installing pandas, numpy and graphviz packages

```
conda install -c conda-forge pandas=2.2.2 numpy=1.26.4
conda install -c conda-forge graphviz python-graphviz
```

Installing openai package

```
conda install conda-forge::openai==2.30.0
```

### 2. API key

1. Be signed up for OpenAI.

   https://platform.openai.com/

(NOTE: For new users: Sign up, create an account, generate an API key by providing an API Key Name and a Project Name when prompted. Copy the generated key and store it in a safe, secure location — you'll need it to access the API.)

2. Once logged in, click your profile icon (top-right corner) → Manage Account → Billing.

3. In the Billing section, set up Prepaid Billing or Auto Recharge.
   Prepaid: manually add credit (e.g., $5, $10).
   Auto Recharge: automatically top up when balance is low.

4. Check your usage: open Usage from the left-hand menu to monitor your monthly spend and remaining balance. Pricing: https://openai.com/api/pricing/

5. Go to OpenAI API keys: https://platform.openai.com/api-keys

6. Click "Create new secret key" → copy the key (it looks like `sk-...`).

Important: treat this key like a password — never share it or commit it to public code repositories.

---

## Files Required

Your working directory (`SampleBioShift/`) should contain the following:

```
SampleBioShift/
├── SampleBioShift.py
├── tree_pipeline.R
├── target.txt
├── sampletree_control.txt
├── methods.txt
├── Inputs/
│   ├── Cell.csv
│   ├── Pro1log10.csv
│   ├── BacCount.csv
│   └── BacTree.nwk
└── BioShift_Req/
    ├── BioShift.py
    ├── config_bioshift.txt
    └── graphviz/
        ├── Graphviz1.txt
        ├── Graphviz2.txt
        └── Graphviz3.txt
```

There is no separate `differential_analysis.R` — differential analysis is built directly into `tree_pipeline.R` (see below).

### Input files (Inputs/)

All CSV files contain a `Sample` column with sample IDs; other columns hold frequency, abundance, or scaled expression values depending on the data type. Any number of representations can be used — `SampleBioShift.py` reads whichever `Filename`s are listed in `methods.txt` at run time, so adding, removing, or renaming a representation only ever requires editing `methods.txt`, never the Python or R code.

1. **Cell.csv** — cell frequency data. Used for sample-tree clustering and differential analysis.
2. **Pro1log10.csv** — protein/cytokine expression data (log10). Used for sample-tree clustering and differential analysis.
3. **BacCount.csv** — bacterial abundance data (count). Used for sample-tree clustering (including phylogenetic methods) and differential analysis.
4. **BacTree.nwk** — reference phylogenetic tree, required only when a representation's `Distance_Metric` is a phylogenetic one (UniFrac, UniFracW, MPD, MPDw, MNTD, MNTDw).

---

## File Descriptions

1. **SampleBioShift.py**

   The orchestrator — run this one. For each block in `methods.txt` it runs `tree_pipeline.R`, then builds `Observed_shifts/` from whatever clades were found (`build_observed_shifts()`, inlined — no separate script), then calls `BioShift_Req/BioShift.py` for the disease and healthy passes. All three steps write into a fresh `RunN/` folder under whatever `FolderName` you pass on the command line; archiving happens even if a later step fails (e.g. a bad API key), so the tree/clustering results are never lost.

2. **tree_pipeline.R**

   Clusters samples (phylogenetic and/or non-phylogenetic methods), detects clades of targeted (`Y`) samples against the global `target.txt` cohort, and — when `Differential analysis = Yes` — runs the differential analysis (limma or MaAsLin2, via the built-in `run_trend_for_job()`) on every clade it finds. There is no separate `differential_analysis.R` file anymore; this is all one script.

3. **BioShift_Req/BioShift.py**

   Data curation and calling large language models for interpretation.
   For details: [BioShift on GitHub](https://github.com/dprabin25/BioShift).

4. **methods.txt**

   The config file `tree_pipeline.R` reads. Supports either a single block (one representation per run) or several stacked `# FILE 1` / `# FILE 2` / ... blocks — `SampleBioShift.py` runs `tree_pipeline.R` once per block automatically and restores your original file afterward. Fields per block:

   - `Filename` — input CSV under `Inputs/`. This name is also the representation's identity everywhere downstream (output folder, combined-file tags) — nothing else needs to be told about it.
   - `Input_Type` — `log10` | `Count` | `Frequency`
   - `Distance_Metric` — non-phylo: `Bray` | `Euclidean`; phylo: `UniFrac` | `UniFracW` | `MPD` | `MPDw` | `MNTD` | `MNTDw`
   - `PhyloTree` — a `.nwk` file, or `None` (required only for phylo distance metrics)
   - `Normalization` — depends on `Input_Type`: `log10` → `MinMax` | `Zscore` (Euclidean only); `Count` → `TSS` | `Hellinger` | `CLR` (Euclidean only) | `Rarefaction`; `Frequency` → `None` (Bray) | `Hellinger` (Bray) | `Zscore` (Euclidean only)
   - `Support/bootstrap` — bootstrap replicate count
   - `SD` — `0` = no injected noise (control run) | `>0` = Gaussian/Dirichlet/Poisson noise magnitude
   - `Differential analysis` — `Yes` | `No`. **This is a real gate, not decorative.** `Yes`: clade detection runs, `support_tree_highlighted.jpeg` is drawn, and TREND (limma/MaAsLin2) runs on every kept clade. `No`: the pipeline stops right after the plain `support_tree.nwk/pdf/png` — no clade detection, no highlighted tree, no TREND output for that block, and `build_observed_shifts()` skips it entirely (there is nothing in `Outputs/` for it to read).
   - `Input df`, `Profile_Library` — kept in the file for readability/documentation but not read by `tree_pipeline.R`: the differential-analysis library is chosen automatically from `Input_Type` (`log10` → limma, `Count`/`Frequency` → MaAsLin2), so `Profile_Library` should just describe what will actually run.

   See the invalid-combination notes at the bottom of `methods.txt` itself for combinations the pipeline will refuse (e.g. `CLR + Bray`, `MinMax + Count`, a phylo metric with `PhyloTree = None`).

5. **sampletree_control.txt**

   Configures clade detection: `min_targeted`, `max_clade_size`, `max_others`, `assign_policy`.
   `assign_policy` options:
   - `best` (default): ranks clades by most targeted → fewer others → smaller total tips.
   - `first`: selects the first qualifying clade encountered.
   - `largest`: prefers clades with more total tips.
   - `smallest`: prefers clades with fewer total tips.

6. **target.txt**

   Two columns: `Sample`, `Target`. Assign `Y` to every sample of interest — together these form the one global targeted cohort that clade detection is run against, with no notion of which patient a sample came from. A `Patient` column is not required and, if present, is purely informational; it plays no role in clade detection.

7. **BioShift_Req/config_bioshift.txt**

   Needs to be updated with your API key and the large language model version you want to use. `SampleBioShift.py` validates this key with a live API call right before the `BioShift.py` step, so a missing or wrong key aborts there without losing the tree/clustering/differential-analysis results already produced in steps 1–2.

8. **BioShift_Req/graphviz/**

   Pre-built DOT-format pathway diagrams. `BioShift.py` colors any node matching a significant element and renders a highlighted image per diagram per combined-clade result.

---

## Running the Script

1. Go to the working directory containing all the required files (`Inputs/`, `methods.txt`, `SampleBioShift.py`, `tree_pipeline.R`, `target.txt`, `sampletree_control.txt`, `BioShift_Req/`) on the Anaconda terminal.

2. Run the script with your chosen folder name:

   ```
   python SampleBioShift.py FolderName
   ```

   You can assign any name to `FolderName`. Its `RunN/` subfolder will contain every output from that run, plus `log.txt` with the full run log.

3. Re-running with the same `FolderName` creates the next sequential run folder — `Run1`, `Run2`, `Run3`, etc. — without touching earlier runs.

Optional flags: `--samples` (restrict BioShift to specific sample IDs), `--python-exe` / `--rscript-exe` (override the Python/R executables if `python`/`Rscript` aren't on your PATH).

## Output Structure

Everything from a run lands under `FolderName/RunN/`:

```
FolderName/RunN/
├── log.txt
├── Outputs/
│   └── <file_base>/<Normalization>/<Distance_Metric>/
│       ├── SD0/            (Control run — no injected noise, always present)
│       │   ├── support_tree.nwk / .pdf / .png
│       │   ├── support_tree_highlighted.jpeg      (only if Differential analysis = Yes and a clade was found)
│       │   ├── clades/                             (which samples make up each detected clade)
│       │   └── trend_outputs/group_<...>_node<N>/   (only if Differential analysis = Yes)
│       │       ├── Input_<file_base>_<library>.csv
│       │       ├── sig_features_<file_base>_<library>.csv
│       │       └── box plots
│       └── SD<value>/      (Test run — only when methods.txt's SD > 0; same structure as SD0,
│                             a noise-robustness check run in parallel)
└── BioShift/
    ├── Observed_shifts/
    │   └── Combined_Observed_Shifts.csv         (or Combined_01_..., Combined_02_..., ...)
    └── BioShiftOutputs/
        ├── Disease/<combo_stem>/
        │   ├── elements/
        │   ├── prompts/
        │   ├── tables/
        │   ├── graphviz/
        │   └── Prompt_Co_Output/
        └── Healthy/<combo_stem>/
            └── (same structure)
```

**`Outputs/`** — one subfolder per representation × normalization × distance metric configured in `methods.txt`, each holding a `SD0/` Control run (always) and, only if that block's `SD > 0`, a parallel `SD<value>/` Test run as a noise-robustness check. `support_tree_highlighted.jpeg` uses a colorblind-safe legend: **blue = Target (Y)**, **vermillion = Contamination** (non-target samples admitted into a clade), **grey = Other**; node labels show support/bootstrap values. `trend_outputs/` holds one `group_*` folder per clade found, with the raw differential-analysis output for that clade.

**`BioShift/Observed_shifts/`** — built only from each representation's **`SD0` (Control)** results; the `SD<value>` Test run is written for reference but is not fed into BioShift. This is the combined input BioShift reads. `build_observed_shifts()` merges each clade's `Input_*.csv` (identical value everywhere → kept; disagreement → 0), then combines across representations. If every representation found exactly one clade, this is a single `Combined_Observed_Shifts.csv`. If a representation's targeted cohort split into more than one clade, each of that representation's clades pairs separately with the other representations' clade(s), producing `Combined_01_...`, `Combined_02_...`, etc. — each numbered file is a genuinely distinct biological grouping, not noise to collapse into one.

**`BioShift/BioShiftOutputs/`** — LLM interpretation of each combined file, split into `Disease/` and `Healthy/` passes, one subfolder per combined-file stem.

---

## Reference

[1] Prabin Dawadi, Ryan M Tobin, Jorge Frias-Lopez, Alpdogan Kantarci, Flavia Teles, Sayaka Miura. Uncovering Periodontal Ecosystem Complexity with Sample Trees. (2025) Under Review

---

## Copyright 2025, Authors and University of Mississippi

BSD 3-Clause "New" or "Revised" License, which is a permissive license similar to the BSD 2-Clause License except that it prohibits others from using the name of the project or its contributors to promote derived products without written consent.
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
