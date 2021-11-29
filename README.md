# AlphScore

### Overview
Improving pathogenicity prediction of missense variants by using AlphaFold-derived features

This Repository is structured into two parts: 

- Extract features:

Extract features from Alphafold2 derived structures ( https://alphafold.ebi.ac.uk/ ) using the tools DSSP, FEATURE, protinter ( https://github.com/maxibor/protinter ) and the python modules biopython, biopandas and networkx. Finally the generated features are combined added to the data of dbNSFP 4.2.
![alt text](https://github.com/Ax-Sch/AlphScore/blob/main/Overview.png?raw=true)


- Analysis:

Here the extracted and merged features are analysed and used for prediction of pathogenicity. This part is written in R.

abcdef

### Install dependencies of feature extraction
This workwas tested on Ubuntu / Xubuntu Linux using miniconda3.

Make sure to have conda installed, if not install it (e.g. https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html ).

Clone or download the repository to a directory with plenty of free storage, then change directory into the AlphScore-main folder. Then set up the conda environment:

```
cd [folder where the repository was cloned/downloaded to]/AlphScore
conda env create --file=environment/alphafold_environment.yaml
conda activate alphafold1
chmod +x tools/dssp/mkdssp
chmod +x tools/feature/feature-3.1.0/bin/featurize
```

Note that DSSP, the FEATURE framework and protinter are practically not available via conda. Currently, these tools are uploaded to this repository to facilitate installing.

### Run feature extraction
Feature extraction can be run using the follwoing commands:

```
cd extract_features
snakemake --cores [number of cores you have]

```

### Run analysis
After the feature extraction is finished, the R-scripts in analysis/scripts can be run consecutively, as specified e.g. in the file preprocess_first_prediction.sh.

```
cd analysis
bash preprocess_first_prediction.sh

```
The prediction file has several parameters that can be set. To get an overview of all possible parameters run:

```
Rscript scripts/03_prediction.R --help
```

