#! /usr/bin/env python

# $Id: leaveoneout.py,v 1.3 2004/02/22 21:22:46 mliang Exp $
# Copyright (c) 2001 Mike Liang.  All rights reserved.

# Perform cross validation by leaving one protein out

# Requires all parameters necessary for training and scan
#
# [Feature Training]
# SiteFile:
# NumberOfShells:
# ShellThickness:
#
# [Feature Scan]
# ProteinList:
# GridSize:
# ExcludedResidues:


import os
import sys
import thread
from sitefile import SiteFile, Site


# =================================================================
# TRAINING COMMAND
def FeatureTrainCommand(dirname, numshells, shellthickness):
    return 'cd %s; train_model.py -s train.site -n %d -w %g' % (dirname, numshells, shellthickness)

# =================================================================
# SCAN COMMAND
def FeatureScanCommand(dirname, gridsize, excludedresidues):
    return 'cd %s; scan_model.py -f test.site -g %g -x "%s" Model/training.ini' % (dirname, gridsize, excludedresidues)

# =================================================================
# Executes Feature Commands
def RunCommands(dirname, numshells, shellthickness, gridsize, excludedresidues):
    os.system(FeatureTrainCommand(dirname, numshells, shellthickness))
    os.system(FeatureScanCommand(dirname, gridsize, excludedresidues))

# =================================================================
# Makes a copy of object attributes
def copyAttributes(src,dst,attrList):
    for attr in attrList:
        setattr(dst,attr,getattr(src,attr))

# =================================================================
# Makes a copy of sitefile, w/o site list
def copySiteFile(orig):
    aCopy = SiteFile()
    copyAttributes(orig, aCopy, ("name", "radius", "filename"))
    return aCopy

# =================================================================
# Write SiteFile to a file
def writeSiteFile(filename, sFile):
    fh = open(filename, "w")
    fh.write(str(sFile))
    fh.close()

# =================================================================
# Checks if any of the sites are SITES
def containsSite(aList):
    for site in aList:
        if site.isSite():
            return 1
    return 0

# =================================================================
# 
def LeaveOneOut(sfile):
    global NumberOfShells, ShellThickness, GridSize, ExcludedResidues

    # get site lists
    siteLists = sfile.sites.items()

    for idx in range(len(siteLists)):
        (pdbid, testSites) = siteLists[idx]

        if not containsSite(testSites):
            continue

        tmpList = siteLists[:idx]+siteLists[idx+1:]
        trainSites = reduce(lambda a,b:(None, a[1]+b[1]), tmpList)[1]

        # generate new sites
        tstFile = copySiteFile(sfile)
        tstFile.addSites(testSites)

        trnFile = copySiteFile(sfile)
        trnFile.addSites(trainSites)

        # Write out sitefiles
        dirname = "LOO_%s" % pdbid
        if not os.path.exists(dirname):
            os.mkdir(dirname)
        writeSiteFile(os.path.join(dirname, "test.site"), tstFile)
        writeSiteFile(os.path.join(dirname, "train.site"), trnFile)

        # Run FEATURE Commands
        thread.start_new_thread(RunCommands,(dirname, NumberOfShells, ShellThickness, GridSize, ExcludedResidues))


# Train on one block, validate on the others


def main():
    global NumberOfShells, ShellThickness, GridSize, ExcludedResidues
    NumberOfShells = 6
    ShellThickness = 1.25
    GridSize = 1.652
    ExcludedResidues = "CA"

    args = sys.argv[1:]
    if len(args) != 1:
        print "Usage: leaveoneout <sitefile>"
        sys.exit(2)
    sitefilename = args.pop(0)
    LeaveOneOut(SiteFile(sitefilename))



if __name__ == "__main__":
    main()
