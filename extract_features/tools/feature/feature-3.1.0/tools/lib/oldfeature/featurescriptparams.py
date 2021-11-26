
#$Id: featurescriptparams.py,v 1.3 2002/04/10 05:23:02 mliang Exp $
# Copyright 2001 (c) Mike Liang.  All Rights Reserved
#
# Globals and common functions used in FeatureScripts
#

import sys
import string
import re
import os
import os.path
import urllib
import ConfigParser
import fetchproteinfile
from featureexts import *


DEFAULT_PARAMETER_FILE=os.environ.get("FEATURE_PARAMETER_FILE","")

FEATURE_PARAM_SECTION="Feature"
FEATURE_TRAINING_PARAM_SECTION="Feature Training"
FEATURE_SCAN_PARAM_SECTION="Feature Scan"
FEATURE_EVAL_PARAM_SECTION="Feature Evaluation"

DEFAULT_FEATURE_PARAMS="""
[%s]
ProteinDir: Proteins
OutputDir: Model

[%s]
ExcludedResidues:
PLevel:                  0.01
OutputDir:               Scan-%%(scanName)s
""" % (FEATURE_TRAINING_PARAM_SECTION, FEATURE_SCAN_PARAM_SECTION)



featureBinDir=os.environ.get("FEATURE_PATH","")

# Feature program configurations
TRAINING_EXE=os.path.join(featureBinDir, "training")
TRAINING_PARAM_FILENAME="training_params.txt"
TRAINING_CONFIG_FILENAME="training.ini"

SCAN_EXE=os.path.join(featureBinDir, "scan")
SCAN_PARAM_FILENAME="scan_params.txt"
SCAN_CONFIG_FILENAME="scan.ini"


# ----------------------------------------------------------------------
def DisplayParameterFileFormat():
# ----------------------------------------------------------------------
    """ Displays Parameter File Format
    """
    
    print """
    Parameter file format: 
        [%s]

        [%s]
        SiteFile:                    Full path to site file
        NumberOfShells:              Number of shells to build feature
        ShellThickness:              Thickness of each shell (angstroms)

        OutputDir:                   Directory to output model

        [%s]
        ProteinList:                 File with list of PDB ids OR
                                     Protein site file (.site)
        GridSize:                    Size of grid (angstroms)
        ExcludedResidues:            List of residues to exclude

        PLevel:                      p-value for significance
        TrainingModelFile:           Training model
        OutputDir:                   Directory to output scan results

        [%s]
        Cutoff:                      Cutoff threshold
        Radius:                      Radius threshold
        ScanHitsDir:                 Directory containing Scan hit results
    """ % (FEATURE_PARAM_SECTION, FEATURE_TRAINING_PARAM_SECTION,
           FEATURE_SCAN_PARAM_SECTION, FEATURE_EVAL_PARAM_SECTION)



from sitefile import SiteFile
# ----------------------------------------------------------------------
def ReadSites(sitefile):
# ----------------------------------------------------------------------
    """ Returns list of PDBids used in sitefile
    """
    sfile = SiteFile(sitefile)
    pdbids = sfile.sites.keys()
    pdbids.sort()

    return pdbids



# ----------------------------------------------------------------------
def FetchPDB(pdbid, filePath, cache=1):
# ----------------------------------------------------------------------
    """ Gets PDB file from Local Cache/Global Cache/URL
    """
    return fetchproteinfile.fetchPDBFile(pdbid, filePath, cache)



# ----------------------------------------------------------------------
def FetchDSSP(pdbid, filePath, cache=1):
# ----------------------------------------------------------------------
    """ Gets DSSP file from Local Cache/Global Cache/URL
    """
    return fetchproteinfile.fetchDSSPFile(pdbid, filePath, cache)


# ----------------------------------------------------------------------
def VerifyOrCreateDir(dirpath):
# ----------------------------------------------------------------------
    """ Verifies that dirpath is directory, or creates the directory if it does not exist.  Returns 0 if error, or 1 if okay
    """
    return fetchproteinfile.verifyOrCreateDir(dirpath)



# Load Parameter File
def UseConfig(config):
    pass

