# Sample_Tree

## Description
The **SampleBioShift** pipeline runs the **SampleTree** workflow and **BioShift** simultaneously.  
`SampleBioShift.py` is a small orchestration script that runs the full analysis pipeline for a given patient/sample set, keeping the project directory clean and reproducible.

Make sure you have both **R** and **Python** installed, with the necessary packages. Please check `requirements.txt` for installation details.

## Dependencies

### 1a. Anaconda
Please install Anaconda: https://www.anaconda.com/distribution/

Open Anaconda terminal and then create conda environment for Bioshift. 
Please install Anaconda: https://www.anaconda.com/distribution/
Note: We tested in 3.12.1 and 3.12.2 

- Make sure you have these packages installed:
  
`argparse, pandas, shutil, fnmatch, re, subprocess, sys, itertools, datetime`


### 1b. R
Please install R from https://posit.co/download/rstudio-desktop/

Note: We tested in 4.5.0

- Make sure you have these packages installed:
  
`optparse, ape, vegan, picante, phyloseq, phangorn, progress, ggplot2, ggtree, treeio, dplyr, readr, stringr, tibble, tools, reshape2, limma`



---

## Files Required

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


## File Descriptions

1. **sampletree_simple.R**
   
   For clustering samples based on phylogenetic and/or non-phylogenetic methods.

2. **ObservedShifts.py**
   
   Merges outputs from **BioShift.py**.

3. **BioShift.py**
   
   Data curation and calling large language models for interpretation.  
   For details: [BioShift on GitHub](https://github.com/dprabin25/BioShift).

4. **methods.txt**
   
   Contains four columns: `File`, `Method`, `Boot`, `Library`, `Tree`.  
   Users can enter the scaled input file name with the `.csv` extension that they want to run with the phylogenetic method.  
   - Note: You need to provide count data to run the Bray-Curtis dissimilarity method.  
   - For other methods (MPD, MPDw, MNTD, MNTDw, UniFrac, UniFracW), you need both count data and a tree file.  

   Users can choose the package to work with: `"limma"` or `"MaAslin2"` based on their data type.

   Example:
   File Method Boot Library Tree
Pro1log10_with_ID.csv Bray 0 limma NA
Cell_with_ID.csv Bray 0 MaAsLin2 NA
BacLog10Freq_with_ID.csv UniFrac 100 MaAsLin2 eHOMD_Ribosomal_Protein_Tree_1_pruneReName.nwk


5. **sampletree_control.txt**
   
Configure the clade assignment with `min_targeted`, `max_other_samples`, `max_total_samples`, and `assign_policy`.  
`assign_policy` options:
- `best`: (default) Ranks clades by most targeted → fewer others → smaller total tips.
- `first`: Selects the first qualifying clade encountered.
- `largest`: Prefers clades with more total tips.
- `smallest`: Prefers clades with fewer total tips.

6. **target.txt**
   
This file should contain `Sample` and `Target` columns.  
Assign `Y` for the samples of interest.

7. **config.txt**
   
This needs to be updated with your API key and the version of the large language model you want to use.

---

## Running the Script

1. Go to the working directory containing all the required files from the command terminal.

2. Run the script with your chosen folder name:

`python SampleBioShift.py FolderName`


- You can assign any name to FolderName. This folder will contain:

The folder will contain

-All output files

-target.txt you used

-log.txt with details about the run in Run1 inside FolderName.

3. If you use the same FolderName for subsequent runs, it will create sequential run folders: Run1, Run2, Run3, etc.

## Output Structure

The output is saved inside the FolderName you assigned. It contains:

### SampleTree Outputs:

Outputs for each input file, with the suffix of the method used.

Example:

-BacLog10Freq_with_ID_UniFrac

-Pro1log10_with_ID_Bray

### Observed_Shifts_by_group:

Contains elements increasing or decreasing in different clades from input files.

### Observed_Shifts:

The output consists of combined results where one clade is selected at a time from each applicable input CSV file (CSV1, CSV2, CSV3, etc.) and then merged.

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

[1] Sayaka Miura, Prabin Dawadi, Ryan M Tobin, Jorge Frias-Lopez, Alpdogan Kantarci, Flavia Teles.  Uncovering Periodontal Ecosystem Complexity with Sample Trees. (2025) Under Review

--------
## Copyright 2025, Authors and University of Mississippi

BSD 3-Clause "New" or "Revised" License, which is a permissive license similar to the BSD 2-Clause License except that that it prohibits others from using the name of the project or its contributors to promote derived products without written consent. 
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
