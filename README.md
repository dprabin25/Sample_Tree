# Sample_Tree

## Description
The **SampleBioShift** pipeline runs the **SampleTree** workflow and **BioShift** simultaneously.  
`SampleBioShift.py` is a small orchestration script that runs the full analysis pipeline for a targeted sample set.

SampleTree groups samples into clusters based on their similarity, using similarity matrices to identify which samples are most alike. BioShift, on the other hand, focuses on identifying and defining the changes (or shifts) of elements within these sample clusters. Essentially, SampleTree organizes the samples, while BioShift analyzes the shifts that occur within those groups.

Note: Users can also run only Sample tree approach skipping the BioShift if they don't update `config.txt`.

### Important note:
The pipeline does not perform input scaling. It is up to the user to decide whether to use scaled or unscaled data.

## Dependencies

### 1a. Anaconda

Please install Anaconda: https://www.anaconda.com/distribution/

Open the Anaconda Prompt (terminal), then create and activate a conda environment for Bioshift.

Note: We tested Bioshift with Python 3.12.1 / 3.12.2 and R 4.5.0.


 `conda create -n bioshift python=3.12.2 r-base=4.5.0 pandas -y`

`conda activate bioshift` 

Make sure you have these packages installed:
argparse, pandas, shutil, fnmatch, re, subprocess, sys, itertools, datetime
