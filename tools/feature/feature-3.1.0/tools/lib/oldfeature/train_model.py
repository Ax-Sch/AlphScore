#! /usr/bin/env python

# $Id: train_model.py,v 1.8 2004/02/22 21:22:46 mliang Exp $
# Copyright (c) 2001 Mike Liang.  All rights reserved.
#
# Creates model of site
#

import sys
import os
import os.path
import string
import getopt
import ConfigParser
import StringIO
from featurescriptparams import *



# ----------------------------------------------------------------------
def Usage(extra=0):
# ----------------------------------------------------------------------
    """ Displays program Usage
    """
    
    print """
    Usage: train_model [options] [<parameterfile.ini>]

    options:
        -h, -?, --help               Display Help File
        --xhelp                      Display Extra Help

        --sitefile,
        -s <sitefile>                Full path to site file

        --numshells,
        -n <number of shells>        Number of shells to build feature

        --shellthickness,
        -w <shell thickness>         Thickness of each shell (angstroms)

        --outputdir,
        -O <output directory>        Directory to store output files
    """

    if extra:
        DisplayParameterFileFormat()



# ----------------------------------------------------------------------
def ParseCommandLineArguments(args):
# ----------------------------------------------------------------------
    """ Parses command line arguments
    """

    # Parse User Options
    try:
        (optlist, paramFileList) = getopt.getopt(args, "h?s:n:w:O:",
            ("help","xhelp","sitefile=","numshells=","shellthickness=","outputdir="))
    except getopt.GetoptError:
        Usage()
        sys.exit(2)

    cmdArgs = StringIO.StringIO()
    print >> cmdArgs, "[%s]" % FEATURE_TRAINING_PARAM_SECTION
    for (option, value) in optlist:
        if option in ("-h","-?","--help"):
            Usage()
            sys.exit()
        elif option in ("--xhelp",):
            Usage(1)
            sys.exit()
        elif option in ("-s", "--sitefile"):
            print >> cmdArgs, "SiteFile: %s" % value
        elif option in ("-n", "--numshells"):
            print >> cmdArgs, "NumberOfShells:  %s" % value
        elif option in ("-w", "--shellthickness"):
            print >> cmdArgs, "ShellThickness: %s" % value
        elif option in ("-O", "--outputdir"):
            print >> cmdArgs, "OutputDir: %s" % value
    cmdArgs.seek(0)

    # Load Parameter Files
    paramFileList.insert(0, DEFAULT_PARAMETER_FILE)
    config = ConfigParser.ConfigParser()
    try:
        config.readfp(StringIO.StringIO(DEFAULT_FEATURE_PARAMS))
        config.read(paramFileList)
        config.readfp(cmdArgs)
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


# Verify Parameters Specified
for param in ("SiteFile", "NumberOfShells", "ShellThickness", "OutputDir"):
    if not config.has_option(FEATURE_TRAINING_PARAM_SECTION, string.lower(param)):
        print "Missing Parameter: %s" % param
        sys.exit(2)


# Extract Parameters from Config
siteFile = config.get(FEATURE_TRAINING_PARAM_SECTION, "SiteFile")
numShells = config.getint(FEATURE_TRAINING_PARAM_SECTION, "NumberOfShells")
shellThickness = config.getfloat(FEATURE_TRAINING_PARAM_SECTION, "ShellThickness")
outputDir = config.get(FEATURE_TRAINING_PARAM_SECTION, "OutputDir")
proteinDir = os.path.join(outputDir, "Proteins")


# Check if sitefile exists
if not os.access(siteFile, os.R_OK):
    print >> sys.stderr, "Could not read sitefile: %s" % siteFile
    sys.exit(1)


# Extract fields from sitefilename
(siteDir, siteName) = os.path.split(os.path.normpath(siteFile))
(siteName, siteExt) = os.path.splitext(siteName)
if siteDir == "":
    siteDir = "."
if siteExt != SITEFILE_EXT:
    print >> sys.stderr, "Site filename does not end in %s" % SITEFILE_EXT
    sys.exit(1)


# Create Directories if they don't exist
print "Creating directories..."
if not VerifyOrCreateDir(proteinDir):
    print >> sys.stderr, "Error creating directory %s" % proteinDir
    sys.exit(1)
if not VerifyOrCreateDir(outputDir):
    print >> sys.stderr, "Error creating directory %s" % outputDir
    sys.exit(1)


# Fetch Protein Files
print "Fetching protein files...",
pdblist=ReadSites(siteFile)
for pdbid in pdblist:
    sys.stdout.write(".")
    # Fetch PDB file
    FetchPDB(pdbid, os.path.join(proteinDir, "%s%s" % (pdbid, PDB_EXT)), 0)
    # Fetch DSSP file
    FetchDSSP(pdbid, os.path.join(proteinDir, "%s%s" % (pdbid, DSSP_EXT)), 0)
print


# Generate Training Parameter File
print "Generating Training Parameter file..."
trainingParamFile=os.path.join(outputDir, TRAINING_PARAM_FILENAME)
tpFile = open(trainingParamFile, "w")
print >> tpFile, """
NUM-OF-SHELLS:     %d
SHELL-THICKNESS:   %g
ANALYSIS-NAME:     %s
ALL-SITES-PATH:    %s
PROTEINS-PATH:     %s
STAT-FILES-PATH:   %s
""" % (numShells, shellThickness, siteName, siteDir, proteinDir, outputDir)
tpFile.close()


# Display Configuration
print "Storing Configuration file..."
cfgOut = open(os.path.join(outputDir, TRAINING_CONFIG_FILENAME), "w")
config.write(cfgOut)
cfgOut.close()


# Run Feature-Training
print "Running Feature training..."
os.system("%s %s" % (TRAINING_EXE, trainingParamFile))
