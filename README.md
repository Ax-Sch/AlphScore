# AlphScore

This code belongs to the project "Predicting the pathogenicity of missense variants based on AlphaFold-derived features". If you simply would like to obtain pre-calculated scores, please visit https://uni-bonn.sciebo.de/s/iE2GcXYUPoYgWHl and download the file all_possible_values_concat.csv.gz. 

### Short description of the workflow

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

### Run the snakemake pipeline
Feature extraction can be run using the follwoing commands:

```
snakemake --cores [number of cores you have]

```
The snakemake pipeline will execute >100,000 jobs. To avoid all the computation, you can set the variable testing=True at the beginning of the pipeline. After running the pipeline in testing mode, you could replace the file data/train_testset1/gnomad_extracted_prepro_rec.csv.gz by the one from https://uni-bonn.sciebo.de/s/iE2GcXYUPoYgWHl . Then you would have the full data set for fitting machine learning models yourself. The evaluation in ClinVar variants will be on fewer proteins, however.

This code was tested on a HPC cluster with CentOS Linux 7, miniconda3 and the job scheduler Slurm. If you have any questions you can also contact me.
