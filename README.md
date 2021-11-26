# AlphScore

### Overview
Improving pathogenicity prediction of missense variants by using AlphaFold-derived features

This Repository is structured into two parts: 
- extract features
Extract features from Alphafold2 derived structures ( https://alphafold.ebi.ac.uk/ ) using the tools DSSP, FEATURE, protinter and biop... . Finally these Features are combined with data from dbNSFP 4.2.
![alt text](https://github.com/Ax-Sch/AlphScore/blob/main/Overview.png?raw=true)


- Analysis
Here the extracted and merged features are analysed and used for prediction of pathogenicity. This part is written in R.


### Install
This was tested on Ubuntu / Xubuntu Linux using miniconda3.

Make sure to have conda installed, if not install it (e.g. https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html ).

Clone or download the repository to a directory with plenty of free storage, then change directory into the AlphScore-main folder. Then set up the conda environment.

```
cd [folder where the repository was cloned/downloaded to]/AlphScore
conda env create --file=environment/alphafold_environment.yaml
conda activate alphafold1
chmod +x tools/dssp/mkdssp
chmod +x tools/feature/feature-3.1.0/bin/featurize
```

### Run
```
snakemake --cores 1
```
