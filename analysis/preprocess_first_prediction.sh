#!/bin/bash
# The first two scripts need to be run to be able to run the prediction.
# To see the command line options of the 01_prediction script, add the --help flag.

Rscript 01_create_validation_set.R
Rscript scripts/01_preprocess.R -i /media/axel/Dateien/Arbeit_Gen/alphafold2/data_from_xcat_v2/gnomad_extracted_v2.csv.gz ### REPLACE this path with your relative path
Rscript scripts/01_preprocess.R -i data/validation/validation_set.csv.gz
Rscript scripts/02_recalibrate_variants.R
Rscript scripts/03_prediction.R -w TRUE -v data/preprocess/validation_set.csv.gzpreprocessed.csv.gz
