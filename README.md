# AlphScore

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
chmod +x tools/dssp/mkdssp
```
Additionally, dssp (version 3) is required. For ease of installation the dssp executable and the required Boost libraries are provided within the folder tools/dssp. Note that the Boost software license applies, which can be found in tools/dssp/LICENSE_1_0.txt . DSSP has additional dependencies (e.g. libc6, libgcc1, libstdc++6) which might not be installed on your system; if this is the case you could try to install dssp/mkdssp via the package manager of your (linux) operating system. To check if the provided dssp executable works, execute the following commands in the AlphScore folder:

```
export LD_LIBRARY_PATH="tools/dssp/"
tools/dssp/mkdssp
```

Additionally the feature framework version 3.1.0 needs to be obtained from the following website: https://simtk.org . Its executable should be located in the following path: tools/feature/feature-3.1.0/bin/featurize 
Ensure that the executable works proberly:

```
chmod +x tools/feature/feature-3.1.0/bin/featurize
tools/feature/feature-3.1.0/bin/featurize
```

If you do not have conda installed, obtain it e.g. from https://docs.conda.io/projects/conda/en/latest/user-guide/install/linux.html ) and set up the conda environment:

```
conda env update --file workflow/envs/AlphScore.yaml
conda activate AlphScore
```


### Run the snakemake pipeline
The snakemake pipeline can then be run by executing the following command.

```
snakemake --cores [number of cores you have]
```

The snakemake pipeline will execute >100,000 jobs. To perform the computation on a reduced set of ~110 proteins, you can set the variable testing=True at the beginning of the pipeline. This code was tested on a HPC cluster with CentOS Linux 7, miniconda3 and the job scheduler Slurm. If you have any questions do not hesitate to contact me.


