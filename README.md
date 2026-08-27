# SampleTree

## Description

**SampleTree finds where your targeted samples cluster, tests which features differ significantly in those clusters, and gives the biological interpretation — in one command.**

It works on any sample data: each CSV in `Inputs/` is a set of features (columns) measured across your samples (rows) — cell types, proteins, microbial taxa, or anything else you're tracking. It works in two layers:

- Clustering: builds a similarity tree from your samples and detects **clades of targeted (`Y`) samples** within it. Every sample flagged `Y` in `target.txt` belongs to one shared targeted cohort, and clade detection runs across that whole cohort at once — whether those samples cluster together or not.
- Interpreting : interprets the features found significant within each clade and uses an large language model to explain what those shifts mean biologically.

`SampleBioShift.py` is the orchestrator — the only script you run. It calls `tree_pipeline.R` once per representation in `methods.txt`, builds the combined `Observed_shifts/` input BioShift needs, then calls `BioShift_Req/BioShift.py` for the disease and healthy interpretation passes.

### Run modes

`python SampleBioShift.py FolderName` always executes the same three steps, but how far each run actually gets depends on `methods.txt` and `BioShift_Req/config_bioshift.txt` — no separate mode flag is needed:

1. **SampleTree only** — set `Differential analysis = No` on every block in `methods.txt`. `tree_pipeline.R` stops right after the plain `support_tree.nwk/pdf/png`: no clade detection, no highlighted tree, no differential analysis. Use this to just see how samples cluster before deciding anything else.
2. **SampleTree → df** — set `Differential analysis = Yes` (the default), but leave `BioShift_Req/config_bioshift.txt`'s `KEY` unset or invalid. Clade detection, `support_tree_highlighted.jpeg`, and the differential analysis (TREND: limma/MaAsLin2) all run, and `Observed_shifts/` gets built — `SampleBioShift.py` then aborts cleanly at the API-key check, before touching OpenAI. Use this to get the full clustering + differential-analysis result without spending API credits.
3. **SampleTree → df → BioShift** — `Differential analysis = Yes` and a valid `KEY` in `config_bioshift.txt`. The full pipeline runs end to end, including the LLM interpretation passes in `BioShift/BioShiftOutputs/`.

Because steps 1–2 never touch OpenAI, a run always reaches at least mode 2's output even if mode 3 was intended but the API key turns out to be missing or wrong — nothing already computed is lost.

Independently of run mode, `methods.txt` can hold **one file or several**: a single block runs that one representation through whichever mode above is active, while stacking multiple `# FILE 1` / `# FILE 2` / ... blocks runs each representation in turn in the same call and — if `Differential analysis = Yes` and each finds a clade — combines their results into `Observed_shifts/` (see `methods.txt` below).


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

### Example input files (Inputs/)

The method is built around **your own data**, not these specific files — `SampleBioShift.py` reads whichever `Filename`s are listed in `methods.txt` at run time, so adding, removing, or renaming a representation only ever requires editing `methods.txt`, never the Python or R code. The files below, shipped in this repo's `Inputs/`, are examples showing the shape your own CSVs should take: a `Sample` column with sample IDs, and other columns holding frequency, abundance, or scaled expression values depending on the data type. Anyone can point this pipeline at their own files by matching this format.

1. **Cell.csv** — example cell frequency data. Used for sample-tree clustering and differential analysis.

E.g.

 <img width="525" height="345" alt="image" src="https://github.com/user-attachments/assets/7233e7cd-8e96-47e1-8b04-6792f78f1664" />

3. **Pro1log10.csv** — example protein/cytokine expression data (log10). Used for sample-tree clustering and differential analysis.
4. **Micro.csv** — example bacterial relative-abundance data (frequency-type), an alternate representation of the same kind of microbiome data as `BacCount.csv`.
5. **BacCount.csv** — example bacterial abundance data (count-type). Used for sample-tree clustering (including phylogenetic methods) and differential analysis.
6. **BacTree.nwk** — example reference phylogenetic tree, required only when a representation's `Distance_Metric` is a phylogenetic one (UniFrac, UniFracW, MPD, MPDw, MNTD, MNTDw).


---

## File Descriptions

1. **SampleBioShift.py**

   The orchestrator — run this one. For each block in `methods.txt` it runs `tree_pipeline.R`, then builds `Observed_shifts/` from whatever clades were found (`build_observed_shifts()`, inlined — no separate script), then calls `BioShift_Req/BioShift.py` for the disease and healthy passes. All three steps write into a fresh `RunN/` folder under whatever `FolderName` you pass on the command line; archiving happens even if a later step fails (e.g. a bad API key), so the tree/clustering results are never lost.

2. **tree_pipeline.R**

   Clusters samples (phylogenetic and/or non-phylogenetic methods), detects clades of targeted (`Y`) samples against the global `target.txt` cohort, and — when `Differential analysis = Yes` — runs the differential analysis (limma or MaAsLin2, via the built-in `run_trend_for_job()`) on every clade it finds. 

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
   - `Profile_Library` — `limma` | `MaAsLin2`. **This is a real gate.** Set explicitly, it controls which library `run_trend_for_job()` uses. Left blank, it falls back to the `Input_Type`-derived default (`log10` → limma, `Count`/`Frequency` → MaAsLin2). An unrecognized value stops the run with an error.
   - `Input df` — the differential-analysis is for the clade detected runs in this input (Clustered samples vs Other Samples).

   See the invalid-combination notes at the bottom of `methods.txt` itself for combinations the pipeline will refuse (e.g. `CLR + Bray`, `MinMax + Count`, a phylo metric with `PhyloTree = None`).

5. **sampletree_control.txt**

   Configures clade detection: `min_targeted`, `max_clade_size`, `max_others`, `assign_policy`.
   `assign_policy` options:
   - `best` (default): ranks clades by most targeted → fewer others → smaller total tips.
   - `first`: selects the first qualifying clade encountered.
   - `largest`: prefers clades with more total tips.
   - `smallest`: prefers clades with fewer total tips.

E.g. 

<img width="938" height="105" alt="image" src="https://github.com/user-attachments/assets/83281168-d528-496e-8f64-698be552df9f" />



6. **target.txt**

   Two columns: `Sample`, `Target`. Assign `Y` to every sample of interest — together these form the one global targeted cohort that clade detection is run against, with no notion of which patient a sample came from. A `Patient` column is not required and, if present, is purely informational; it plays no role in clade detection.
   
   E.g.
   
   <img width="234" height="419" alt="image" src="https://github.com/user-attachments/assets/03315b4c-9b5f-43f8-bfb5-a7e42ecdae08" />


8. **BioShift_Req/config_bioshift.txt**

   Needs to be updated with your API key and the large language model version you want to use. `SampleBioShift.py` validates this key with a live API call right before the `BioShift.py` step, so a missing or wrong key aborts there without losing the tree/clustering/differential-analysis results already produced in steps 1–2.

   E.g.
   
   <img width="471" height="205" alt="image" src="https://github.com/user-attachments/assets/37448f42-eec3-430e-bbf4-dc0d61a79c61" />


10. **BioShift_Req/graphviz/**

   Pre-built DOT-format pathway diagrams. `BioShift.py` colors any node matching a significant element and renders a highlighted image per diagram per combined-clade result.

---

## Running the Script

1. Go to the working directory containing all the required files (`Inputs/`, `methods.txt`, `SampleBioShift.py`, `tree_pipeline.R`, `target.txt`, `sampletree_control.txt`, `BioShift_Req/`) on the Anaconda terminal.

E.g. 

<img width="892" height="157" alt="image" src="https://github.com/user-attachments/assets/ffb0ada6-d848-4e0e-9bc6-1ee3394ba746" />


3. Run the script with your chosen folder name:

   ```
   python SampleBioShift.py FolderName
   ```

E.g.

<img width="463" height="31" alt="image" src="https://github.com/user-attachments/assets/16db97cf-b2e1-46ad-9ad0-813a668ef35c" />



You can assign any name to `FolderName`. Its `RunN/` subfolder will contain every output from that run, plus `log.txt` with the full run log.

3. Re-running with the same `FolderName` creates the next sequential run folder — `Run1`, `Run2`, `Run3`, etc. — without touching earlier runs.

Optional flags: `--python-exe` / `--rscript-exe` (override the Python/R executables if `python`/`Rscript` aren't on your PATH).

## Output Structure

Everything from a run lands under `FolderName/RunN/`:

```
FolderName/RunN/
├── log.txt
├── Outputs/
│   ├── <file_base>/<Normalization>/<Distance_Metric>/       (non-phylo: Bray | Euclidean)
│   │   ├── SD0/            (Control run — no injected noise, always present)
│   │   │   ├── support_tree.nwk / .pdf / .png
│   │   │   ├── support_tree_highlighted.jpeg      (only if Differential analysis = Yes and a clade was found)
│   │   │   ├── clades/                             (which samples make up each detected clade)
│   │   │   └── trend_outputs/group_<...>_node<N>/   (only if Differential analysis = Yes)
│   │   │       ├── Input_<file_base>_<library>.csv
│   │   │       ├── sig_features_<file_base>_<library>.csv
│   │   │       └── box plots
│   │   └── SD<value>/      (Test run — only when methods.txt's SD > 0; same structure as SD0,
│   │                         a noise-robustness check run in parallel)
│   └── <file_base>/<Distance_Metric>/                        (phylo: UniFrac | UniFracW | MPD | MPDw | MNTD | MNTDw)
│       ├── support_tree.nwk / .pdf / .png
│       ├── support_tree_highlighted.jpeg          (only if Differential analysis = Yes and a clade was found)
│       ├── clades/
│       └── trend_outputs/group_<...>_node<N>/       (only if Differential analysis = Yes)
│           ├── Input_<file_base>_<library>.csv
│           ├── sig_features_<file_base>_<library>.csv
│           └── box plots
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

**`Outputs/`** — one subfolder per representation configured in `methods.txt`, shaped differently depending on `Distance_Metric`:
- **Non-phylo** (`Bray`, `Euclidean`) — nested under `<Normalization>/<Distance_Metric>/`, with a `SD0/` Control run (always) and, only if that block's `SD > 0`, a parallel `SD<value>/` Test run as a noise-robustness check. Support comes from SD-noise replicates.
- **Phylo** (`UniFrac`, `UniFracW`, `MPD`, `MPDw`, `MNTD`, `MNTDw`) — output goes straight to `Outputs/<file_base>/<Distance_Metric>/`, with no `Normalization` or `SD0` nesting: there's no SD sweep for these (`SD` in `methods.txt` is ignored), support instead comes from raw-count bootstrap replicates.

Either way, `support_tree_highlighted.jpeg` uses a colorblind-safe legend: **blue = Target (Y)**, **vermillion = Contamination** (non-target samples admitted into a clade), **grey = Other**; node labels show support/bootstrap values. `trend_outputs/` holds one `group_*` folder per clade found, with the raw differential-analysis output for that clade.

E.g. Cluster 

<img width="1173" height="564" alt="image" src="https://github.com/user-attachments/assets/eefadf46-57d3-4688-a368-1adf835a6b25" />

E.g. differential-analysis

<img width="443" height="251" alt="image" src="https://github.com/user-attachments/assets/489b694d-4f06-4e5f-b655-d9232023ba83" />




**`BioShift/Observed_shifts/`** — built only from each representation's **Control** results (non-phylo: the `SD0/` run; phylo: the one run there is, since phylo jobs don't have a Control/Test split). For non-phylo jobs, a `SD<value>` Test run is written for reference but is not fed into BioShift. This is the combined input BioShift reads. `build_observed_shifts()` merges each clade's `Input_*.csv` (identical value everywhere → kept; disagreement → 0), then combines across representations. If every representation found exactly one clade, this is a single `Combined_Observed_Shifts.csv`. If a representation's targeted cohort split into more than one clade, each of that representation's clades pairs separately with the other representations' clade(s), producing `Combined_01_...`, `Combined_02_...`, etc. — each numbered file is a genuinely distinct biological grouping, not noise to collapse into one.

E.g. 

<img width="199" height="235" alt="image" src="https://github.com/user-attachments/assets/7f38fd10-f38a-49a2-ab12-8d32b90a8064" />

_1= increase, -1 = decrease_


**`BioShift/BioShiftOutputs/`** — LLM interpretation of each combined file, split into `Disease/` and `Healthy/` passes, one subfolder per combined-file stem.

E.g. 

<img width="1287" height="95" alt="image" src="https://github.com/user-attachments/assets/96f22c08-b132-4938-89e7-1cbb11188b35" />


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
