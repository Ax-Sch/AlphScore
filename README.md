# Source code of AlphScore

This code belongs to the project "Predicting the pathogenicity of missense variants using features derived from AlphaFold2". This is the source code of the project; precalculated scores are deposited under the DOI 10.5281/zenodo.6288139 . 

### Workflow

In the first part of the workflow features of each amino acid from AlphaFold structures ( https://alphafold.ebi.ac.uk/ ) are extracted. For feature extraction the tools DSSP, FEATURE, protinter and the python modules biopython, biopandas and networkx are used. The extracted features are then added to the data of dbNSFP 4.2a and dbNSFP is subsequently used to form variant sets (gnomAD, ClinVar). Then these variants are used to train tree-based machine learning classifiers. 

Workflow:
![alt text](https://github.com/Ax-Sch/AlphScore/blob/main/dag.jpg?raw=true)

### Install the requiered dependencies
First clone or download the repository and the modified version of the tool protinter to a directory with plenty (~500Gb) of free storage by. Protinter was modified from its original version to increase its speed and to be able to obtain residues that are spatially close. 

```
git clone https://github.com/Ax-Sch/AlphScore
cd AlphScore/tools
git clone https://github.com/Ax-Sch/protinter.git
chmod +x protinter/protinter
cd ../
```

If you do not have conda installed, please now obtain it e.g. from https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html ) and set up the main conda environment and the conda environment in which the FEATURE framework will run:

```
conda env update --file workflow/envs/AlphScore.yaml
conda env update --file workflow/envs/dssp.yaml
```

Additionally the FEATURE framework version 3.1.0 needs to be obtained from the following website: https://simtk.org . To install the feature framework place the file feature-3.1.0-src.tar.gz (obtained from simtk.org) to the folder tools/feature. Then execute the following commands:

```
cd tools/feature
conda activate DSSP
tar -xf feature-3.1.0-src.tar.gz 
cd feature-3.1.0/
make
cd ../../..
```
After this, the executable of FEATURE should be located in the following path: tools/feature/feature-3.1.0/bin/featurize 

### Run the snakemake pipeline
The snakemake pipeline can then be run by executing the following command.

```
conda activate AlphScore
snakemake --cores [number of cores you would like to use] --use-conda --conda-frontend conda
```

At the top of the file workflow/Snakefile the variable testing is set to True (testing=True). This reduces the computation to a set of about 210 proteins, to test if everything works as expected without the need of processing all proteins. You can set the variable testing to False (change testing=True to testing=False) to run all proteins, which will result in >100,000 jobs. This code was tested on a HPC cluster with CentOS Linux 7, miniconda3 and the job scheduler Slurm. 

If you have any questions do not hesitate to contact me.

