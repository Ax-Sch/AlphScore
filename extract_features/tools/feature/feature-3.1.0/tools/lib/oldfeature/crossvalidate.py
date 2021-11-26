#! /usr/bin/env python

# $Id: crossvalidate.py,v 1.5 2004/02/22 21:22:46 mliang Exp $
# Copyright (c) 2001 Mike Liang.  All rights reserved.

# Perform k-cross validation

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
import random
import getopt
import threading
import time
from feature.sitefile import SiteFile, Site

PROGNAME = sys.argv[0]
ARGS = sys.argv[1:]


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
    global threadSem
    threadSem.acquire()
    os.system(FeatureTrainCommand(dirname, numshells, shellthickness))
    os.system(FeatureScanCommand(dirname, gridsize, excludedresidues))
    threadSem.release()

# =================================================================
# Splits list into specified number of blocks
def splitList(aList, num):
    retList = []
    size = len(aList)/num
    for idx in range(num-1):
        retList.append(aList[idx*size:(idx+1)*size])
    retList.append(aList[(idx+1)*size:])
    return retList

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
# Split SiteFile into separate blocks
def SplitSiteFile(sfile, numblocks, prefix="CV_"):
    global NumberOfShells, ShellThickness, GridSize, ExcludedResidues, OutputDir
    
    # get site lists
    sites = sfile.getSites(Site.SITE_TYPE)
    random.shuffle(sites)
    sites = splitList(sites,numblocks)
    nonsites = sfile.getSites(Site.NONSITE_TYPE)
    random.shuffle(nonsites)
    nonsites = splitList(nonsites,numblocks)

    for count in range(numblocks):
        tstSites = sites[count]
        tstNonSites = nonsites[count]
        trnSites = reduce(lambda x,y:x+y, sites[:count]+sites[count+1:])
        trnNonSites = reduce(lambda x,y:x+y, nonsites[:count]+nonsites[count+1:])

        # generate new sites
        tstFile = copySiteFile(sfile)
        tstFile.addSites(tstSites)
        tstFile.addSites(tstNonSites)
        tstFile.name += "_test%02d" % count
        trnFile = copySiteFile(sfile)
        trnFile.addSites(trnSites)
        trnFile.addSites(trnNonSites)
        trnFile.name += "_train%02d" % count

        # Write out sitefiles
        dirname = os.path.join(OutputDir, "%s%02d" % (prefix, count))
        CreateDir(dirname)
        writeSiteFile(os.path.join(dirname, "test.site"), tstFile)
        writeSiteFile(os.path.join(dirname, "train.site"), trnFile)

        threading.Thread(target=RunCommands,args=(dirname, NumberOfShells, ShellThickness, GridSize, ExcludedResidues)).start()


# Train on one block, validate on the others
def Usage(ec=None):
    print """
    Usage: %s [options] <sitefile> <numblocks>
        Options:
        -O <outputdir>
        -P <prefix>
        -t <max num threads>

        Training:
        -n <numshells>
        -w <shellthickness>

        Scan:
        -g <gridsize>
        -x <excludedresidues>
""" % sys.argv[0]
    if ec != None:
        sys.exit(ec)

def CreateDir(dirname):
    if not os.path.exists(dirname):
        os.mkdir(dirname)

def main():
    global ARGS
    global NumberOfShells, ShellThickness, GridSize, ExcludedResidues, OutputDir
    global MaxNumThreads
    global threadSem
    
    NumberOfShells = 6
    ShellThickness = 1.25
    GridSize = 1.652
    ExcludedResidues = "CA"

    OutputDir = "."
    prefix = "CV_"
    MaxNumThreads = 4

    (opts, ARGS) = getopt.getopt(ARGS, "O:P:t:n:w:g:x:", 
        ("numshells=","shellthickness=", "gridsize=", "excludedresidues="))
    for (opt, arg) in opts:
        if opt in ("-O",):
            OutputDir = arg
        elif opt in ("-P",):
            prefix = arg
        elif opt in ("-t",):
            MaxNumThreads = int(arg)
        elif opt in ("-n","--numshells"):
            NumberOfShells = int(arg)
        elif opt in ("-w","--shellthickness"):
            ShellThickness = float(arg)
        elif opt in ("-g","--gridsize"):
            GridSize = float(arg)
        elif opt in ("-x","--excludedresidues"):
            ExcludedResidues = arg
    if len(ARGS) != 2:
        Usage(2)
    sitefilename = ARGS.pop(0)
    numblocks = int(ARGS.pop(0))

    threadSem = threading.Semaphore(MaxNumThreads)

    CreateDir(OutputDir)
    SplitSiteFile(SiteFile(sitefilename), numblocks, prefix)

    # Wait for threads to end
    main = threading.currentThread()
    children = threading.enumerate()
    children.remove(main)
    for child in children:
        child.join()


if __name__ == "__main__":
    main()
