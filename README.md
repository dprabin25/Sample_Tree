
# SampleTreeBioShift

The **SampleBioShift** pipeline runs the **SampleTree** workflow and **BioShift** simultaneously.  
`SampleBioShift.py` is a small orchestration script that executes the full analysis pipeline for a given patient/sample set, keeping the project directory clean and reproducible.

### Requirements

Make sure you have both **R** and **Python** installed, with the necessary packages. Please check `requirements.txt` for installation details.


### Files Required

Your working directory should contain the following files:

WorkingDirectory/
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

less
Copy code

### File Descriptions

1. **sampletree_simple.R**  
   For clustering samples based on phylogenetic and/or non-phylogenetic methods.
   
2. **ObservedShifts.py**  
   Merges outputs from **BioShift.py**.

3. **BioShift.py**  
   Data curation and calling large language models for interpretation.  
   For details: [BioShift on GitHub](https://github.com/dprabin25/BioShift).

4. **methods.txt**  
   Contains four columns: `File`, `Method`, `Boot`, `Library`, `Tree`.  
   Users specify input files and methods like **Bray-Curtis**, **UniFrac**, etc., along with the package to use (`limma` or `MaAsLin2`).

   Example:
File Method Boot Library Tree
Pro1log10_with_ID.csv Bray 0 limma NA
Cell_with_ID.csv Bray 0 MaAsLin2 NA
BacLog10Freq_with_ID.csv UniFrac 100 MaAsLin2 eHOMD_Ribosomal_Protein_Tree_1_pruneReName.nwk

r
Copy code

5. **sampletree_control.txt**  
Configure clade selection policies, such as `min_targeted`, `max_other_samples`, and `assign_policy`.  
`assign_policy` options:
- `best`: Ranks by most targeted → fewer others → smaller total tips.
- `first`: Selects the first qualifying clade.
- `largest`: Prefers clades with more total tips.
- `smallest`: Prefers clades with fewer total tips.

6. **target.txt**  
Should contain `Sample` and `Target` columns.  
Set `Y` for samples of interest and `N` for others.

7. **config.txt**  
Contains your API key and the version of the large language model to use.

---

### Running the Script

1. Go to the working directory containing all the required files from the command terminal.

2. Run the script with your chosen folder name:
```bash
python SampleBioShift.py FolderName
You can assign any name to FolderName. This folder will contain:

All output files

target.txt you used

log.txt with details about the run in Run1 inside FolderName.

If you use the same FolderName for subsequent runs, it will create sequential run folders: Run1, Run2, Run3, etc.

Output Structure
The output is saved inside the FolderName you assigned. It contains:

SampleTree Outputs:
Outputs for each input file, with the suffix of the method used.
Example:

BacLog10Freq_with_ID_UniFrac

Pro1log10_with_ID_Bray

Observed_Shifts_by_group:
Contains elements increasing or decreasing in different clades from input files.

Observed_Shifts:
The output consists of combined results where one clade is selected at a time from each applicable input CSV file (CSV1, CSV2, CSV3, etc.) and then merged.

Process:

One clade is selected from each CSV file individually, as applicable.

These selected clades are merged to create various combinations.

Folder Structure:

Each combination of selected clades is stored as a separate CSV file, containing the merged results for:

-Element

-Observed Shift

Note: The number of CSV files used in each combination depends on the user's input. Some combinations may involve just two or one file, based on statistically significant features.
