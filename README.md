# DIA-SIFT
#### A post-processing tool to improve protein quantification when using Stable Isotope Labeling by Amino acids in Cell culture (SILAC) and LC-IMS-MS acquisition. From the Martin Lab at the University of Michigan.

**[Martin Lab Homepage](https://sites.google.com/umich.edu/brentmartinlab/home)**

<br />  

## DIA-SIFT instructions:
*(requires a python environment with pandas and numpy libraries)*:

After searching 2-plex SILAC .raw files with ProteinLynx Global Server (PLGS), move all final_fragment and final_peptide .csv files to the same directory for filtering.

<br />

DIA-SIFT is performed in two steps:  

(1) fragment ion SILAC ratios (Light / Heavy) are extracted and technical replicates are merged (using the script `merge_and_preprocess_PLGS_outputs.py`)

(2) Light / Heavy SILAC ratios are quartile filtered for each protein (using the script `DIA-SIFT_filtering.py`)

<br /> 

### How to run DIA-SIFT:

**1.** download `merge_and_preprocess_PLGS_outputs.py` and `DIA-SIFT_filtering.py`

**2.** in `merge_and_preprocess_PLGS_outputs.py`, set the current working directory (`cwd`) to the folder containing PLGS peptide & fragment files, and set the file name `prefix` and suffixes (`p_suffix` and `f_suffix`) as indicated in the docstring (copied below).

```
This script pre-processes ProteinLynx Global Server (PLGS) output .csv files for SILAC quantile filtering (DIA-SIFT).

INPUTS (each of these must be defined below):
- [cwd] folder containing peptide and fragment files from 3 replicate injections of each sample
- [prefix] include only the first part of the filenames right up until a unique sample identifier
- [p_suffix] include only the part of the peptide filenames immediately following the unique sample identifier
- [f_suffix] include only the part of the fragment filenames immediately following the unique sample identifier

EXAMPLE INPUTS:
SILAC_HeLa_SampleA_replicate001_IA_final_peptide.csv
SILAC_HeLa_SampleA_replicate001_IA_final_fragment.csv

[SILAC_HeLa_]SampleA_replicate001[_IA_final_peptide.csv]
      ^                                     ^
    prefix                               p_suffix
    
[SILAC_HeLa_]SampleA_replicate001[_IA_final_fragment.csv]
      ^                                     ^
    prefix                               f_suffix
    
prefix = 'SILAC_HeLa_'
p_suffix = '_IA_final_peptide.csv'
f_suffix = '_IA_final_fragment.csv'
(the unique sample identifier would be 'SampleA_replicate001')


OUTPUTS:
- merged and pre-processed .csv files for peptides & fragments:
    "[SampleName]_peptideMergedAndPreprocessed.csv"
    "[SampleName]_fragmentMergedAndPreprocessed.csv"
```

**3.** in `DIA-SIFT_filtering.py`, set the current working directory (`cwd`) to the folder containing the pre-processed PLGS peptide & fragment files (from steps 1 - 2 above). Docstring is copied below.

```
This script performs the main quantile filtering function of DIA-SIFT. 

INPUT:
- [cwd] folder containing processed .csv files (from merge_and_preprocess_PLGS_outputs.py)
    "[SampleName]_peptideMergedAndPreprocessed.csv"
    "[SampleName]_fragmentMergedAndPreprocessed.csv"
    
OUTPUT:
- filtered & unfiltered .csv files in the same directory
    "[SampleName]_unfiltered.csv"
    "[SampleName]_DIA-SIFTed.csv"
    "[SampleName]_unfiltered_protein.csv"
    "[SampleName]_DIA-SIFTed_protein.csv"
```

The first two files contain all peptide and fragment observations, before and after DIA-SIFT filtering. (For this reason, fragment-only information will be missing from rows containing peptide observations, and vice versa.) The `_protein.csv` files contain protein-level quantification results. All values refer to Light / Heavy SILAC ratios.
