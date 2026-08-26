# Sample_Tree

## Description
The **SampleBioShift** pipeline runs the **SampleTree** workflow and **BioShift** simultaneously.  
`SampleBioShift.py` is a small orchestration script that runs the full analysis pipeline for a targeted sample set. 

SampleTree groups samples into clusters based on their similarity, using similarity matrices to identify which samples are most alike. BioShift, on the other hand, focuses on identifying and defining the changes (or shifts) of elements within these sample clusters. Essentially, SampleTree organizes the samples, while BioShift analyzes the shifts that occur within those groups.

Note: Users can also run only Sample tree approach skipping the BioShift if they don't update config.txt. 

### Important note:
The pipeline does not perform input scaling. It is up to the user to decide whether to use scaled or unscaled data.

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

Installing opanai package
```
conda install conda-forge::openai==2.30.0   
```

### 2. API key
1. Be signed up for OpenAI.

   https://platform.openai.com/

(NOTE: For new users: Sign up, create an account, generate an API key by providing an API Key Name and a Project Name when prompted. Copy the generated key and store it in a safe, secure location — you’ll need it to access the API.)

3. Once logged in, click your profile icon (top-right corner) → Manage Account → Billing.
   
2. In the Billing section, set up Prepaid Billing or Auto Recharge
   Prepaid: Manually add credit (e.g., $5, $10).
   Auto Recharge: Automatically top up when balance is low.
   
 3. Check Your Usage
   Open Usage from the left-hand menu to monitor your monthly spend and remaining balance
   Link for pricing: https://openai.com/api/pricing/

 4. Go to OpenAI API keys: https://platform.openai.com/api-keys

 5. Click “Create new secret key” → Copy the key (it looks like sk-...).
    
Important: Treat this key like a password — never share it or commit it to public code repositories.
---

##  Files Required

Your working directory should contain the following files:

##  Files Required
Your working directory (SampleBioShift/) should contain the following:

├── Inputs/

   └── File1.csv              

   └── File2.csv 
   
   └──File3.csv               

   └──FileX.csv               
.
.

└── FileX.nwk              

├── SampleBioShift.py      

├── tree_pipeline.R    

├── differential_analysis.R    

├── target.txt       

├── sampletree_control.txt    

├── methods.txt     

├── BioShift/

   └── BioShift.py             
    
   └── config_bioshift.txt     
    
   └── graphviz/              



### Input files according to ExampleInputs (that we shared here in the repository)
All csv files contains Sample column that have sample IDs and other columns contain frequency, abundance, scaled expression values depending upon type of input. 

Example input files

1. Pro1log10.csv: Expression data from protein dataset (log10) for sample tree clustering and differential analysis for targeted samples clustered in sample trees.
 
<img width="449" height="191" alt="image" src="https://github.com/user-attachments/assets/2be0cfe1-e04e-40b0-9110-323c8c5fad94" />


2. Cell.csv : Cell frequency dataset  - Required for sample tree clustering and differential analysis

<img width="414" height="191" alt="image" src="https://github.com/user-attachments/assets/ec0aa0bc-5aff-4af6-bc9f-c53a941bcaa1" />


3. BacCount.csv : Bacteria abundance dataset (Count) - Required for sample tree clustering using phylogenetic methods and differential analysis

<img width="567" height="198" alt="image" src="https://github.com/user-attachments/assets/dd0e8b3a-cbb7-4c31-99a0-1d506ed0b3df" />


4. BacTree.nwk : Required for making sample tree using phyolgenetic methods

<img width="592" height="470" alt="image" src="https://github.com/user-attachments/assets/2a900187-a518-42b2-87a7-9c0429ff2b33" />





##  File Descriptions

##  File Descriptions

1. **SampleBioShift.py**
   
   The orchestrator — run this one. Controls the whole pipeline in one call: runs `tree_pipeline.R` against whatever `methods.txt` currently holds, builds `Observed_shifts/` internally, then calls `BioShift/BioShift.py` for the disease and healthy passes.
2. **tree_pipeline.R**
   
   For clustering samples based on phylogenetic and/or non-phylogenetic methods, and detecting clades of targeted (`Y`) samples.
3. **differential_analysis.R**
   
   Sourced automatically by `tree_pipeline.R` when `Differential analysis = Yes` in `methods.txt`. Runs limma or MaAsLin2 on each detected clade and writes the `Element` / `Observed Shift` table `SampleBioShift.py` reads to build `Observed_shifts/`.
4. **BioShift.py**
   
   Data curation and calling large language models for interpretation.  
   For details: [BioShift on GitHub](https://github.com/dprabin25/BioShift).
5. **methods.txt**
   
   Single config file `tree_pipeline.R` reads — edit it by hand before each run to select one dataset/combo at a time. Fields: `Filename`, `Input_Type` (log10 | Count | Frequency), `Distance_Metric` (Bray | Euclidean | UniFrac | UniFracW | MPD | MPDw | MNTD | MNTDw), `PhyloTree`, `Normalization`, `Support/bootstrap`, `SD`, `Differential analysis` (Yes/No), `Input df`, `Profile_Library`.  
   - Note: You need to provide count data to run the Bray-Curtis dissimilarity method.  
   - For other methods (MPD, MPDw, MNTD, MNTDw, UniFrac, UniFracW), you need both count data and a tree file.  
   Users can choose the package to work with: `"limma"` or `"MaAslin2"` based on their data type.
   You may see methods.txt on our repository for an example.
6. **sampletree_control.txt**
   
   Configure the clade assignment with `min_targeted`, `max_clade_size`, `max_others`, and `assign_policy`.  
   `assign_policy` options:
   - `best`: (default) Ranks clades by most targeted → fewer others → smaller total tips.
   - `first`: Selects the first qualifying clade encountered.
   - `largest`: Prefers clades with more total tips.
   - `smallest`: Prefers clades with fewer total tips.
7. **target.txt**
   
   This file should contain `Sample`, `Target`, and `Patient` columns.  
   Assign `Y` for the samples of interest; `Patient` groups a sample's timepoints so clades are detected per patient rather than across one global cohort.
8. **BioShift/config_bioshift.txt**
   
   This needs to be updated with your API key and the version of the large language model you want to use. `SampleBioShift.py` validates this key with a live API call before the `BioShift.py` step runs, so a missing or wrong key aborts there without losing the tree/clustering results already produced.
9. **BioShift/graphviz/**
   
   Pre-built DOT-format pathway diagrams. `BioShift.py` colors any node matching a significant element green (increase) or blue (decrease) and renders a highlighted image per diagram per patient-group.

---

## Running the Script

1. Go to the working directory containing all the required files (All inputs in csv format along with Methods.txt, ObservedShifts.py, SampleBioShift.py, config.txt, sampletree_control.txt, sampletree_simple.R and target.txt) on the command terminal.
<img width="150" height="235" alt="image" src="https://github.com/user-attachments/assets/c2f554f3-4e1d-4a45-ab08-5898cfc439fb" />


3. Run the script with your chosen folder name:

`python SampleBioShift.py FolderName`

<img width="885" height="203" alt="image" src="https://github.com/user-attachments/assets/04372e2a-ada5-46b9-a536-2f781ad171a1" />



You can assign any name to "FolderName". This folder will contain:

- All output files

- target.txt you used

- log.txt with details about the run in Run1 inside FolderName.

3. If you use the same name for "FolderName" for subsequent runs, it will create sequential run folders: Run1, Run2, Run3, etc.

## Output Structure

The output is saved inside the "FolderName" you assigned. 

The folder contains:

### SampleTree Outputs:

Outputs for each input file, with the suffix of the method used.

Example:

- *_UniFrac

- *_Bray

ExampleOutputs have these folders:
1. BacCount_UniFrac
   
3. Cell_Bray
   
5. Pro1log10MinMax_Bray

Inside this folder, you will find a subfolder called "clades", where you can see which elements are shifting significantly in the targeted samples that formed clades, as demonstrated in the boxplots and statistical summary.

### Observed_Shifts_by_group:

Contains elements increasing or decreasing in different clades from input files.

ExampleOutputs shows folder name "Observed_Shifts_by_group" 



### Observed_Shifts:

The output consists of combined results with combinations from CSV files from "Observed_Shifts_by_group" folder which have input for BioShift.  

ExampleOutputs shows folder name "Observed_Shifts" 

#### Process:

One clade is selected from each CSV file individually, as applicable.

These selected clades are merged to create various combinations.

#### Folder Structure:

Each combination of selected clades is stored as a separate CSV file, containing the merged results for:

- Element

- Observed Shift

Note: The number of CSV files used in each combination depends on the user's input. Some combinations may involve just two or one file, based on statistically significant features. 

### BioShiftOutputs

Contains data interpretation for disease and healthy for each specific combinations of the clades. 

ExampleOutputs shows folder name "BioShiftOutputs" 


--------
## Reference

[1] Prabin Dawadi, Ryan M Tobin, Jorge Frias-Lopez, Alpdogan Kantarci, Flavia Teles, Sayaka Miura.  Uncovering Periodontal Ecosystem Complexity with Sample Trees. (2025) Under Review

--------
## Copyright 2025, Authors and University of Mississippi

BSD 3-Clause "New" or "Revised" License, which is a permissive license similar to the BSD 2-Clause License except that that it prohibits others from using the name of the project or its contributors to promote derived products without written consent. 
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
