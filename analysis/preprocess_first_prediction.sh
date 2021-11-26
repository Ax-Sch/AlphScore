#!/bin/bash
# The first two scripts need to be run to be able to run the prediction.
# To see the command line options of the 01_prediction script, add the --help flag.

Rscript scripts/01_preprocess.R
Rscript scripts/02_recalibrate_variants.R
Rscript scripts/01_prediction.R
