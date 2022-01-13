# AlphScore

This code belongs to the project "Predicting the pathogenicity of missense variants based on AlphaFold-derived features". If you simply would like to obtain pre-calculated scores, please visit XXX.

### Short workflow description

In the first part of the workflow features of each amino acid from AlphaFold structures ( https://alphafold.ebi.ac.uk/ ) are extracted using a snakemake workflow. For feature extraction the tools DSSP, FEATURE, protinter ( https://github.com/maxibor/protinter ) and the python modules biopython, biopandas and networkx are used. The extracted features are then added to the data of dbNSFP 4.2 and dbNSFP is subsequently used to form variant sets (gnomAD, ClinVar). Then these variants are used to train tree-based machine learning classifiers. 

Workflow:
![alt text](https://github.com/Ax-Sch/AlphScore/blob/main/dag.jpg?raw=true)

### Install the requiered dependencies
First conda should be installed (e.g. https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html ).
Clone or download the repository to a directory with plenty of free storage and change the directory into the AlphScore-main folder. Then set up the conda environment and make the tools executable:

```
cd [folder where the repository was cloned/downloaded to]/AlphScore
conda env create --file=environment/alphafold_environment_from_history.yaml
conda activate alphafold1
chmod +x tools/dssp/mkdssp
chmod +x tools/feature/feature-3.1.0/bin/featurize
```

Note that DSSP, the FEATURE framework and protinter are practically not available via conda. Currently, these tools were placed in the tools folder of this repository to facilitate installing.

### Run snakemake pipeline
Feature extraction can be run using the follwoing commands:

```
snakemake --cores [number of cores you have]

```
The snakemake pipeline will execute >100,000 jobs. For testing, you can set the variable testing=True at the beginning of the pipeline. This code was tested on a HPC cluster with CentOS Linux 7, miniconda3 and the job scheduler Slurm. If you have any questions just message me.
