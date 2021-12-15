#! /usr/bin/env python

# $Id: retrain.py,v 1.3 2004/02/22 21:22:46 mliang Exp $
# Copyright (c) 2001 Mike Liang.  All rights reserved.

# Takes the output of evaluation to retrain the set

import os
import sys
import getopt
import cStringIO
import re
from sitefile import SiteFile, Site

PROGNAME = sys.argv[0]
ARGS = sys.argv[1:]



# ---------------------------------------------------------------------
# User Parameters
class USER:
    outputDir = "Retrain"
    minSens = 0.80
    topFalse = 20
    numIter = None
    siteFile = None
    numBlocks = None
    evaluationSiteFile = None

    genSiteFilename = "retrain.site"
    retrainPrefix = "iter_"
    crossvalidatePrefix = "cv_"

# ---------------------------------------------------------------------
class EvaluationResults:
    def __init__(self, data=None):
        self.info = []
        self.nonsites = []
        if data:
            self.parseData(data)

    def parseData(self, data):
        # read from file handle or string
        if type(data) == type(""):
            input = cStringIO.StringIO(data)
        else:
            input = data

        # Skip everything up to Nonsites
        line = input.readline()
        while line and not re.match("NEW NONSITES", line):
            if re.match("RAD", line):
                words = line.split()
                rad = float(words[1])
                cut = float(words[3])
                sens = float(words[5])
                spec = float(words[7])
                ppv = float(words[9])
                self.info.append((rad, cut, sens, spec, ppv))
            line = input.readline()

        # Parse nonsites
        line = input.readline()
        while line:
            (siteline, score) = line.split(";;")
            site = Site(siteline)
            site.score = float(score)
            self.nonsites.append(site)
            line = input.readline()

    def writeInfo(self, param, closefile=0):
        if type(param) == type(""):
            file = open(param, "w")
            closefile = 1
        else:
            file = param
        print >>file, "RAD:\tCUT:\tSENS:\tSPEC:\tPPV:"
        for entry in self.info:
            print >>file, "%g\t%g\t%g\t%g\t%g" % entry
        if closefile:
            file.close()

    def writeNonsites(self, param, closefile=0):
        if type(param) == type(""):
            file = open(param, "w")
            closefile = 1
        else:
            file = param
        for site in self.nonsites:
            print >>file, "%s ;; %g" % (str(site), site.score)
        if closefile:
            file.close()
        

# ---------------------------------------------------------------------
def CrossValidateCommand(outputdir, sitefile, numblocks, prefix):
    return "crossvalidate.py -t 7 -O %s -P %s %s %d" % \
           (outputdir, prefix, sitefile, numblocks)

# ---------------------------------------------------------------------
def EvaluationCommand(outputdir, prefix, minSens, topFalse, siteFile):
    return "evaluateCrossValidation.py -P %s -n %g -t %d -s %s %s" % \
           (prefix, minSens, topFalse, siteFile, outputdir)

# ---------------------------------------------------------------------
def Usage(ec=None):
    print """
    Usage: %s [options] <numiter> <sitefile> <numblocks>
        Options:
        -O <outputdir>        Output directory
        -t <numfp>            Number of top false positives to retrain
        -m <sens>             Minimum Sensitivity (%%) for false positives
        -e <sitefile>         Evaluation Sitefile
""" % PROGNAME
    if ec != None:
        sys.exit(ec)

# ---------------------------------------------------------------------
def ParseArgs(args):
    # Command Line Arguments
    (opts, args) = getopt.getopt(args, "O:t:m:e:")
    for (opt, arg) in opts:
        if opt in ("-O",):
            USER.outputDir = arg
        elif opt in ("-t",):
            USER.topFalse = int(arg)
        elif opt in ("-m",):
            USER.minSens = float(arg)
        elif opt in ("-e",):
            USER.evaluationSiteFile = arg
    if len(args) != 3:
        Usage(2)
    USER.numIter = int(args.pop(0))
    USER.siteFile = args.pop(0)
    USER.numBlocks = int(args.pop(0))
    if not USER.evaluationSiteFile:
        USER.evaluationSiteFile = USER.siteFile

# ---------------------------------------------------------------------
def CreateDir(dirname):
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    

# ---------------------------------------------------------------------
def main():
    global ARGS

    # Parse Command Line Arguments
    ARGS = ParseArgs(ARGS)

    # Create Output Directory
    CreateDir(USER.outputDir)

    # Load Sitefile
    sFile = SiteFile(USER.siteFile)

    for iter in range(USER.numIter):
        # Create Iteration Directory
        iterdir = os.path.join(USER.outputDir,"%s%d" % (USER.retrainPrefix, iter))
        CreateDir(iterdir)

        # Generate Sitefile
        siteFilename = os.path.join(iterdir, USER.genSiteFilename)
        sFile.unduplicate()
        sFile.write(siteFilename)

        # Cross Validate
        os.system(CrossValidateCommand(iterdir, siteFilename, USER.numBlocks, USER.crossvalidatePrefix))

        # Evaluate
        fh = os.popen(EvaluationCommand(iterdir, USER.crossvalidatePrefix, USER.minSens, USER.topFalse, USER.evaluationSiteFile))
        data = fh.read()
        fh.close()
        evalInfo = EvaluationResults(data)
        sFile.addSites(evalInfo.nonsites)
        evalInfo.writeInfo(os.path.join(iterdir, "eval.out"))
        evalInfo.writeNonsites(os.path.join(iterdir, "new.site"))


if __name__ == "__main__":
    main()
