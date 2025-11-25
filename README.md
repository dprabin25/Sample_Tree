# Sample_Tree

The **SampleBioShift** pipeline runs the **SampleTree** workflow and **BioShift** simultaneously.  
`SampleBioShift.py` is a small orchestration script that runs the full analysis pipeline for a given patient/sample set, keeping the project directory clean and reproducible.

Make sure you have both **R** and **Python** installed, with the necessary packages. Please check `requirements.txt` for installation details.

---

### Files Required

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
   Users can enter the scaled input file name with the `.csv` extension that they want to run with the phylogenetic method.  
   - Note: You need to provide count data to run the Bray-Curtis dissimilarity method.  
   - For other methods (MPD, MPDw, MNTD, MNTDw, UniFrac, UniFracW), you need both count data and a tree file.  

   Users can choose the package to work with: `"limma"` or `"MaAslin2"` based on their data type.

   Example:
