# Sample_Tree

## Description
The **SampleBioShift** pipeline runs the **SampleTree** workflow and **BioShift** simultaneously.  
`SampleBioShift.py` is a small orchestration script that runs the full analysis pipeline for a targeted sample set. 

SampleTree groups samples into clusters based on their similarity, using similarity matrices to identify which samples are most alike. BioShift, on the other hand, focuses on identifying and defining the changes (or shifts) of elements within these sample clusters. Essentially, SampleTree organizes the samples, while BioShift analyzes the shifts that occur within those groups.

Note: Users can also run only Sample tree approach skipping the BioShift if they don't update config.txt. 

### Important note:
The pipeline does not perform input scaling. It is up to the user to decide whether to use scaled or unscaled data.

## Dependencies

### 1. Anaconda and R

- Open Anaconda terminal and then create conda environment for Bioshift.

Please install Anaconda : https://www.anaconda.com/distribution/

To create the environment, include to install correct Python version and R.

- To create the environment, install the correct Python version and R:

```bash
conda create -n bioshift python=3.12.2 r-base=4.5.0 pandas -y

conda activate bioshift

```

- Make sure you have these packages installed for R:
  

```r
install.packages(c(
  "optparse", "ape", "vegan", "picante", "phangorn", "progress",
  "ggplot2", "dplyr", "readr", "stringr", "tibble", "reshape2",
  "data.table", "tidyr", "pbapply", "matrixStats", "Hmisc",
  "quantreg", "lme4", "lmerTest", "Rcpp", "RcppEigen"
))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "phyloseq", "ggtree", "treeio", "limma", "MaAsLin2"
))
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

├── File1.csv ## input CSV file

├── File2.csv ## input CSV file

├── File3.csv ## input CSV file

├── FileX.csv ## input CSV file (select based on your data types)

├── sampletree_simple.R

├── BioShift.py

├── ObservedShifts.py

├── SampleBioShift.py

├── target.txt

├── config.txt

├── methods.txt

├── graphviz/


### Input files according to ExampleInputs (that we shared here in the repository)
All csv files contains Sample column that have sample IDs and other columns contain frequency, abundance, scaled expression values depending upon type of input. 

Example input files
1. Pro1log10MinMax.csv : MinMax scaled expression data from protein dataset for making sample tree

<img width="706" height="240" alt="image" src="https://github.com/user-attachments/assets/d2bd62a7-cb75-4193-a097-d8e7297b3f57" />
 

2. Pro1log10.csv: Expression data from protein dataset (unscaled) for differential analysis for targeted samples clustered in sample trees.
 
<img width="702" height="243" alt="image" src="https://github.com/user-attachments/assets/8c2e1f47-3822-4d65-aa47-d05b007d6efd" />

3. Cell.csv : Cell frequency dataset (unscaled) - Required for sample tree clustering and differential analysis

<img width="673" height="245" alt="image" src="https://github.com/user-attachments/assets/d8aa499d-cbdb-4096-9149-c449b79eb6f5" />

4. BacCount.csv : Bacteria abundance dataset (unscaled) - Required for sample tree clustering using phylogenetic methods and differential analysis

<img width="809" height="241" alt="image" src="https://github.com/user-attachments/assets/f1527507-bbe1-4587-add1-d2f328d03fc4" />

5. BacFreq.csv : MinMax scaled dataset from bacterial abundance. Required for differential analysis.

<img width="807" height="241" alt="image" src="https://github.com/user-attachments/assets/3ba1bc82-1a30-43b4-9d26-07bdeed7fa62" />

6. BacTree.nwk : Required for making sample tree using phyolgenetic methods

<img width="592" height="470" alt="image" src="https://github.com/user-attachments/assets/2a900187-a518-42b2-87a7-9c0429ff2b33" />





##  File Descriptions

1. **sampletree_simple.R**
   
   For clustering samples based on phylogenetic and/or non-phylogenetic methods.

2. **ObservedShifts.py**
   
   Merges outputs from **BioShift.py**.

3. **BioShift.py**
   
   Data curation and calling large language models for interpretation.  
   For details: [BioShift on GitHub](https://github.com/dprabin25/BioShift).

5. **methods.txt**
   
   Contains four columns: `File`, `Method`, `Boot`, `Library`, `Tree`.  
   Users can enter the inputs file name with the `.csv` extension that they want to run with the phylogenetic method.  
   - Note: You need to provide count data to run the Bray-Curtis dissimilarity method.  
   - For other methods (MPD, MPDw, MNTD, MNTDw, UniFrac, UniFracW), you need both count data and a tree file.  

   Users can choose the package to work with: `"limma"` or `"MaAslin2"` based on their data type.

   You may see methods.txt on our repository for an example. 


6. **sampletree_control.txt**
   
Configure the clade assignment with `min_targeted`, `max_other_samples`, `max_total_samples`, and `assign_policy`.  
`assign_policy` options:
- `best`: (default) Ranks clades by most targeted → fewer others → smaller total tips.
- `first`: Selects the first qualifying clade encountered.
- `largest`: Prefers clades with more total tips.
- `smallest`: Prefers clades with fewer total tips.

7. **target.txt**
   
This file should contain `Sample` and `Target` columns.  
Assign `Y` for the samples of interest.

8. **config.txt**
   
This needs to be updated with your API key and the version of the large language model you want to use.

---

## Running the Script

1. Go to the working directory containing all the required files (All inputs in csv format along with Methods.txt, ObservedShifts.py, SampleBioShift.py, config.txt, sampletree_control.txt, sampletree_simple.R and target.txt) on the command terminal.

2. Run the script with your chosen folder name:

`python SampleBioShift.py FolderName`


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

### Observed_Shifts_by_group:

Contains elements increasing or decreasing in different clades from input files.

### Observed_Shifts:

The output consists of combined results with combinations from CSV files from "Observed_Shifts_by_group" folder. 

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
