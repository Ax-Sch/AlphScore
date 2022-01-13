# AlphScore

This code belongs to the project "Improving pathogenicity prediction of missense variants by using AlphaFold-derived features" and is structured into two folders: 

### extract_features

Extracts features of each amino acid from AlphaFold2 structures ( https://alphafold.ebi.ac.uk/ ) using a snakemake workflow. For feature extraction the tools DSSP, FEATURE, protinter ( https://github.com/maxibor/protinter ) and the python modules biopython, biopandas and networkx are used. Finally, the generated features are added to the data of dbNSFP 4.2 and dbNSFP is used to extract variants relevant for the analysis part (see below).
Workflow:
![alt text](https://github.com/Ax-Sch/AlphScore/blob/main/dag.jpg?raw=true)

### analysis:

Here the extracted variants with the AlphFold derived features are analysed and used for prediction of pathogenicity.

### Install dependencies of feature extraction
This code was tested on Ubuntu / Xubuntu Linux using miniconda3.

Conda should be installed, if you need to install it (e.g. https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html ).

Clone or download the repository to a directory with plenty of free storage, then change the directory into the AlphScore-main folder. Then set up the conda environment and make the tools executable:

```
cd [folder where the repository was cloned/downloaded to]/AlphScore
conda env create --file=environment/alphafold_environment.yaml
conda activate alphafold1
chmod +x tools/dssp/mkdssp
chmod +x tools/feature/feature-3.1.0/bin/featurize
```

Note that DSSP, the FEATURE framework and protinter are practically not available via conda. Currently, these tools in the tools folder of this repository to facilitate installing.

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

