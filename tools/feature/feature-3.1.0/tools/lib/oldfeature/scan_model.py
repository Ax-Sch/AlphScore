#! /usr/bin/env python

# $Id: scan_model.py,v 1.10 2004/02/22 21:22:46 mliang Exp $
# Copyright (c) 2001 Mike Liang.  All rights reserved.
#
# Uses model of scan proteins
#

import sys
import os
import os.path
import string
import getopt
import ConfigParser
import StringIO
from featurescriptparams import *

global skipfetch
skipfetch = 0



# ----------------------------------------------------------------------
def Usage(extra=0):
# ----------------------------------------------------------------------
    """ Displays program Usage
    """
    
    print """
    Usage: scan_model [options] [<parameterfile.ini>]

    options:
        -h, -?, --help               Display Help File
        --xhelp                      Display Extra Help Information

        --proteinlist,
        -f <proteinlist>             File with list of PDB ids OR
                                     Protein site file (.site)

        --gridsize,
        -g <gridsize>                Size of grid (angstroms)

        --excludedresidues,
        -x <excluded residues>       List of residues to exclude
    """

    if extra:
        print """
        --plevel,
        -p <plevel>                  p-value for significance

        --trainingmodel,
        -t <training model file>     Training model

        --numshells,
        -n <number of shells>        Number of shells used in training model

        --shellthickness,
        -w <shell thickness>         Shell thickness used in training model

        --outputdir,
        -O <output directory>        Directory to store output files
        """
        DisplayParameterFileFormat()



# ----------------------------------------------------------------------
def ParseCommandLineArguments(args):
# ----------------------------------------------------------------------
    """ Parses command line arguments
    """

    # Parse User Options
    try:
        (optlist, paramFileList) = getopt.getopt(args, "h?f:p:g:x:t:n:w:O:",
            ("help","xhelp","proteinlist=","plevel=","gridsize=","excludedresidues=","trainingmodel=","numshells=","shellthickness=","outputdir=","skipfetch"))
    except getopt.GetoptError:
        Usage()
        sys.exit(2)

    scanArgs = StringIO.StringIO()
    print >> scanArgs, "[%s]" % FEATURE_SCAN_PARAM_SECTION
    trainingArgs = StringIO.StringIO()
    print >> trainingArgs, "[%s]" % FEATURE_TRAINING_PARAM_SECTION
    for (option, value) in optlist:
        if option in ("-h","-?","--help"):
            Usage()
            sys.exit()
        elif option in ("--skipfetch",):
            global skipfetch
            skipfetch = 1
        elif option in ("--xhelp",):
            Usage(1)
            sys.exit()
        elif option in ("-f", "--proteinlist"):
            print >> scanArgs, "ProteinList: %s" % value
        elif option in ("-p", "--plevel"):
            print >> scanArgs, "PLevel:  %s" % value
        elif option in ("-g", "--gridsize"):
            print >> scanArgs, "GridSize: %s" % value
        elif option in ("-x", "--excludedresidues"):
            print >> scanArgs, "ExcludedResidues: %s" % value
        elif option in ("-t", "--trainingmodel"):
            print >> scanArgs, "TrainingModelFile: %s" % value
        elif option in ("-O", "--outputdir"):
            print >> scanArgs, "OutputDir: %s" % value
        elif option in ("-n", "--numshells"):
            print >> trainingArgs, "NumberOfShells: %s" % value
        elif option in ("-w", "--shellthickness"):
            print >> trainingArgs, "ShellThickness: %s" %value
    trainingArgs.seek(0)
    scanArgs.seek(0)

    # Load Parameter Files
    paramFileList.insert(0, DEFAULT_PARAMETER_FILE)
    config = ConfigParser.ConfigParser()
    try:
        config.readfp(StringIO.StringIO(DEFAULT_FEATURE_PARAMS))
        config.read(paramFileList)
        config.readfp(trainingArgs)
        config.readfp(scanArgs)
    except ConfigParser.ParsingError, e:
        print >> sys.stderr, "Invalid Parameter File: ", e
        sys.exit(1)



    return config



# ----------------------------------------------------------------------
# __main__()
# ----------------------------------------------------------------------

# Parse User Arguments
config = ParseCommandLineArguments(sys.argv[1:])
UseConfig(config)


# Optional parameters
# if trainingmodel not specified, then sitefile must be specified
# if both are specified, they must be the same
if config.has_option(FEATURE_TRAINING_PARAM_SECTION, string.lower("SiteFile")):
    siteFile=config.get(FEATURE_TRAINING_PARAM_SECTION, "SiteFile")
    (siteName,siteExt)=os.path.splitext(siteFile)
    if config.has_option(FEATURE_SCAN_PARAM_SECTION, string.lower("TrainingModelFile")):
        # they better be the same
        modelFile=config.get(FEATURE_SCAN_PARAM_SECTION, "TrainingModelFile")
        (modelName,modelExt)=os.path.splitext(modelFile)
        if os.path.basename(siteName) != os.path.basename(modelName):
            print >> sys.stderr, "TrainingModelFile does not match Feature Training parameters"
            sys.exit(1)
    else:
        # set training file to sitefile
        modelDir=config.get(FEATURE_TRAINING_PARAM_SECTION, "OutputDir")
        modelFile="%s%s" % (os.path.basename(siteName), SCOREFILE_EXT)
        modelFile=os.path.join(modelDir,modelFile)
        config.set(FEATURE_SCAN_PARAM_SECTION, string.lower("TrainingModelFile"), modelFile)


# Verify Parameters Specified
for param in ("ProteinList", "PLevel", "GridSize", "ExcludedResidues", "TrainingModelFile", "OutputDir"):
    if not config.has_option(FEATURE_SCAN_PARAM_SECTION, string.lower(param)):
        print "Missing %s Parameter: %s" % (FEATURE_SCAN_PARAM_SECTION, param)
        sys.exit(2)
for param in ("NumberOfShells", "ShellThickness"):
    if not config.has_option(FEATURE_TRAINING_PARAM_SECTION, string.lower(param)):
        print "Missing %s Parameter: %s" % (FEATURE_TRAINING_PARAM_SECTION, param)
        sys.exit(2)


# Extract Parameters from Config
proteinFile = config.get(FEATURE_SCAN_PARAM_SECTION, "ProteinList")
pLevel = config.getfloat(FEATURE_SCAN_PARAM_SECTION, "PLevel")
gridSize = config.getfloat(FEATURE_SCAN_PARAM_SECTION, "GridSize")
excludedResidues = config.get(FEATURE_SCAN_PARAM_SECTION, "ExcludedResidues")
trainingModel = config.get(FEATURE_SCAN_PARAM_SECTION, "TrainingModelFile")
outputDir = config.get(FEATURE_SCAN_PARAM_SECTION, "OutputDir", 1)
numShells = config.getint(FEATURE_TRAINING_PARAM_SECTION, "NumberOfShells")
shellThickness = config.getfloat(FEATURE_TRAINING_PARAM_SECTION, "ShellThickness")


# Check if input files exist
if not os.access(proteinFile, os.R_OK):
    print >> sys.stderr, "Could not read protein file: %s" % proteinFile
    sys.exit(1)
if not os.access(trainingModel, os.R_OK):
    print >> sys.stderr, "Could not read training model: %s" % trainingModel
    sys.exit(1)


# Extract fields from trainingModel
(modelFile, modelExt) = os.path.splitext(os.path.normpath(trainingModel))
if modelExt != SCOREFILE_EXT:
    print >> sys.stderr, "Model filename does not have extension %s" % SCOREFILE_EXT
    sys.exit(1)


# Splice in scanName variable, if necessary
(scanName, ext) = os.path.splitext(os.path.basename(proteinFile))
outputDir = outputDir % {"scanName": scanName}
proteinDir = os.path.join(outputDir, "Proteins")
hitDir = os.path.join(outputDir, "Hits")


# Create Directories if they don't exist
print "Creating directories..."
for dir in [outputDir, proteinDir, hitDir]:
    if not VerifyOrCreateDir(dir):
        print >> sys.stderr, "Error creating directory %s" % dir
        sys.exit(1)


# Fetch Protein Files
if skipfetch:
    print "Skipping fetch protein files..."
else:
    print "Fetching protein files...",
    (base,ext)=os.path.splitext(proteinFile)
    if ext == SITEFILE_EXT:
        pdblist=ReadSites(proteinFile)
    else:
        try:
            pdblist=map(string.strip, open(proteinFile).readlines())
        except:
            print >> sys.stderr, "Error reading protein list %s" % proteinFile
            sys.exit(1)

    for pdbid in pdblist:
        sys.stdout.write(".")
        # Fetch PDB file
        FetchPDB(pdbid, os.path.join(proteinDir, "%s%s" % (pdbid, PDB_EXT)), 0)
        # Fetch DSSP file
        FetchDSSP(pdbid, os.path.join(proteinDir, "%s%s" % (pdbid, DSSP_EXT)), 0)
    print


# Generate Scan Parameter File
print "Generating Scan Parameter file..."
scanParamFile=os.path.join(outputDir, SCAN_PARAM_FILENAME)
spFile = open(scanParamFile, "w")
print >> spFile, """
P-LEVEL:           %g
DELTA-VALUE:       %g
NUM-OF-SHELLS:     %d
SHELL-THICKNESS:   %g
SCORE-FILE:        %s
PROTEINS-PATH:     %s
OUTPUT-PATH:       %s
EXCLUDED-RESIDUES: %s
""" % (pLevel, gridSize, numShells, shellThickness, modelFile, proteinDir, hitDir, excludedResidues)
spFile.close()


# Display Configuration
print "Storing Configuration file..."
cfgOut = open(os.path.join(outputDir,SCAN_CONFIG_FILENAME), "w")
config.write(cfgOut)
cfgOut.close()


# Run Feature-Scan
print "Running Feature scan..."
os.system("%s %s" % (SCAN_EXE, scanParamFile))
