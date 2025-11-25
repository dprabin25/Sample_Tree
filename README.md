# Sample_Tree
SampleBioShift pipeline run Sampletree workflow and BioShift simulateneously. 
SampleBioShift.py is a small `RunAll.py` is a small orchestration script that runs the full analysis pipeline for a given patient/sample set and keeps the project directory clean and reproducible.

Make sure you have both R and python installed with necessary packages (Please check requirements.txt). 

## Files required

WorkingDirectory/

├── File1.csv/
├
├── File2.csv/

-------------
├── File2.csv/
├── ProLog10freq_with_ID_Bray/
├── BioShift.py
├── ObservedShifts.py
├── RunAllFinal.py
├── sampletree_simple.R
├── target.txt
├── config.txt
├── methods.txt
├── graphviz/
├── inputs/
└── (other permanent files)

Your working directory should contain these files
1. sampletree_simple.R
- For clustering samples based on phylogenetic and/or non-phylogenetic methods
2. ObservedShifts.py
- For merging outputs for BioShift.py
3. BioShift.py
- Data curation and calling large language model for interpretation
- For details: https://github.com/dprabin25/BioShift
5. methods.txt
  
  It contains four columns   File, Method, Boot, Library	Tree

Users can enter scaled input file name with csv extension that they want to run with phylogenetic and run phylogenetic method. Note: Users need to provide count data to run Bray-Curtis dissimilarity method while they need count data and tree file for running other methods (MPD, MPDw, MNTD, MNTDw, UniFrac, UniFracW).

Users also can choose the package they want to work with "limma" or "MaAslin2" that they want to run aganist their csv file. 

### Example, 

File                    Method   Boot  Library	Tree
Pro1log10_with_ID.csv          Bray     0     limma	NA
Cell_with_ID.csv    Bray     0     MaAsLin2	NA
BacLog10Freq_with_ID.csv UniFrac  100   MaAsLin2	eHOMD_Ribosomal_Protein_Tree_1_pruneReName.nwk

You may have any number of csv files as your input. It doesn't have to be always 3. 

6. sampletree_control.txt

You can assign min_targeted for your clade of sample interst and max other samples that clade can hold, max total samples the clade can hold and also assign_policy (choose best, first, larget, smallest).

Note about setting assign_policy
best : (default) Ranks candidate clades by: more targeted → fewer others → smaller total tips. Use when you want the cleanest, most on-target clades with minimal contamination.

first : Keeps the first qualifying clade encountered and excludes overlapping ones after. Use for maximum speed and fully deterministic runs (order = tree traversal).

largest: Prefers clades with more total tips (still must meet thresholds). Use when you want broader clades capturing more samples (even if slightly noisier).

smallest: Prefers clades with fewer total tips. Use when you want tighter, highly specific clades.
  
8. target.txt
This .txt should contain Sample and Target columns. Sample columns contains your samples while Target contains 'Y' and 'N' alphabets.
You need to assign 'Y' to the samples of your interst.

10. config.txt
This needs to be updated with your API key and the version of large language model you want to use.


# Running the script
Go to working directory with all the required files from command terminal.

Then run the script SampleBioShift.py as
'python SampleBioShift.py FolderName"

You can assign any name to "FolderName". This folder will contain all the outputs from that run, contain target.txt that you used and have log.txt that contains details about the run in "Run1" within "FolderName". 

If you use the same FolderName, inside that folder it creates sequential run folders: `Run1`, `Run2`, `Run3. This make easy to manipulate Target.txt with User's intension. 


## Output Structure




