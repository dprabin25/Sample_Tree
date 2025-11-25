# Sample_Tree
SampleBioShift pipeline run Sampletree workflow and BioShift simulateneously. 
SampleBioShift.py is a small `RunAll.py` is a small orchestration script that runs the full analysis pipeline for a given patient/sample set and keeps the project directory clean and reproducible.

Make sure you have both R and python installed with necessary packages (Please check requirements.txt). 

## How do you run
Your working directory should contain these files
1. sampletree_simple.R
- For clustering samples based on phylogenetic and/or non-phylogenetic methods
2. ObservedShifts.py
- For merging outputs for BioShift.py
3. BioShift.py
- Data curation and calling large language model for interpretation (For details: https://github.com/dprabin25/BioShift)
5. methods.txt
  It contains four columns   File, Method, Boot, Library	Tree

Users can enter scaled input file name with csv extension that they want to run with phylogenetic and run phylogenetic method. Note: users need to provide count data to run Bray-Curtis dissimilarity method while they need count data and tree file for running other methods (MPD, MPDw, MNTD, MNTDw, UniFrac, UniFracW).

Users also can choose the package they want to work with "limma" or "MaAslin2"

### Example, 

File                     Method   Boot  Library	Tree
Pro1log10_with_ID.csv          Bray     0     limma	NA
Cell_with_ID.csv    Bray     0     MaAsLin2	NA
BacLog10Freq_with_ID.csv UniFrac  100   MaAsLin2	eHOMD_Ribosomal_Protein_Tree_1_pruneReName.nwk



6. sampletree_control.txt
7. target.txt
8. config.txt

## What it does

For each run, the script 

1. **Creates a patient-specific folder**  
   - You pass a folder name (e.g. `Pat1`, `Patient_10737`).
   - Inside that folder it creates sequential run folders: `Run1`, `Run2`, `Run3`, …
