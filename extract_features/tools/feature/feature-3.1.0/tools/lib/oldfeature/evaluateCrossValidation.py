#! /usr/bin/env python

# $Id: evaluateCrossValidation.py,v 1.5 2004/05/22 01:37:42 mliang Exp $
# Copyright (c) 2001 Mike Liang.  All rights reserved.

# Evaluate Cross Validation

import os
import sys
import glob
import getopt
import cStringIO
import re
from Numeric import *
from feature.sitefile import Site
from feature.hitsfile import Hit

VERBOSITY = 0
PRINTFORMAT = 'tabular'


# =================================================================
def fraction(a,b):
    if not a: return 0.0
    return float(a)/(a+b)


# =================================================================
class EvalResults:
    def __init__(self, data=None):
        self.radius = None
        self.cutoff = None
        self.siteTP = 0
        self.siteFN = 0
        self.siteTN = 0
        self.siteFP = 0
        self.pointTP = 0
        self.pointFP = 0
        self.info = {}
        self.hits = []
        self.sites = []

        if data:
            self.parseData(data)

    def parseData(self, data):
        input = cStringIO.StringIO(data)
        line = input.readline()

        while line and not re.match("Table Content", line):
            self.parseVerboseInfo(line[:-1])
            line = input.readline()

        try:
            words = input.readline().split()
            self.cutoff = float(words[-1])
            words = input.readline().split()
            self.radius = float(words[-1])
            words = input.readline().split()
            self.siteTP = float(words[1])
            self.siteFN = float(words[3])
            words = input.readline().split()
            self.siteTN = float(words[1])
            self.siteFP = float(words[3])
            words = input.readline().split()
            self.pointTP = float(words[1])
            self.pointFP = float(words[3])
        except IndexError:
            print "Bad data:"
            sys.stdout.write(data)

    def parseVerboseInfo(self, line):
        if re.match("PDB:", line):
            self.parseSummaryLine(line)
        elif re.match("HIT", line):
            self.parseHitLine(line)
        elif re.match("SITE|NONSITE", line):
            self.parseSiteLine(line)

    def parseSummaryLine(self, line):
        words = line.split()
        (key, pdbid) = (words.pop(0)[:-1], words.pop(0))
        self.info[pdbid] = {}
        while words:
            (key, value) = (words.pop(0)[:-1], words.pop(0))
            try: value = int(value)
            except ValueError: pass
            self.info[pdbid][key] = value

    def parseHitLine(self, line):
        words = line.split()
        recid = words.pop(0)

        hit = Hit()
        hit.pdbid = words.pop(0)
        hit.numSites = int(words.pop(0))
        hit.x = float(words.pop(0))
        hit.y = float(words.pop(0))
        hit.z = float(words.pop(0))
        hit.score = float(words.pop(0))

        self.hits.append(hit)

    def parseSiteLine(self, line):
        words = line.split()
        recid = words.pop(0)

        site = Site()
        site.pdbid = words.pop(0)
        site.numHits = int(words.pop(0))
        site.x = float(words.pop(0))
        site.y = float(words.pop(0))
        site.z = float(words.pop(0))
        site.type = words.pop(0)

        self.sites.append(site)

    def sensitivity(self):
        try: return float(self.siteTP)/(self.siteTP+self.siteFN)
        except ZeroDivisionError: return None
            
        
    def specificity(self):
        try: return float(self.siteTN)/(self.siteTN+self.siteFP)
        except ZeroDivisionError: return None

    def precision(self):
        try: return float(self.pointTP)/(self.pointTP+self.pointFP)
        except ZeroDivisionError: return None

    def tuple(self):
        return (self.siteTP, self.siteFN, self.siteTN, self.siteFP,
                self.pointTP, self.pointFP)


# =================================================================
# EVALUATION COMMAND
def FeatureEvaluationCommand(radius, cutoff, hitsdir, sitefile):
    return 'evaluation -v -r %g -c %g -p "%s" "%s"' % (radius, cutoff, hitsdir, sitefile)

# =================================================================
# Accumulates a list of tuples
def accumulate(aList):
    if not aList:
        return aList
    transpose = apply(zip, aList)
    sum = map(lambda l:reduce(lambda a,b:a+b,l), transpose)
    return sum
    
# =================================================================
# Averages a list of tuples
def average(aList):
    if not aList:
        return aList
    count = len(aList)
    sum = accumulate(aList)
    average = map(lambda v,n=float(count):v/n, accumulate)
    return average
    
# =================================================================
# Performs one evaluation on a particular directory
def evaluateRun(radius, cutoff, hitsDir, siteFile):
    cmd = FeatureEvaluationCommand(radius, cutoff, hitsDir, siteFile)
    if VERBOSITY > 2:
        print >>sys.stderr, cmd
    fh = os.popen(cmd)
    data = fh.read()
    fh.close()
    return EvalResults(data)

def getFalsePositives(resultList):
    falsePositives = {}

    for result in resultList:
        for hit in result.hits:
            if hit.numSites == 0 and hit.score >= result.cutoff:
                key = (hit.pdbid,hit.x,hit.y,hit.z)
                falsePositives.setdefault(key,[]).append((hit.score,hit))

    retval = []
    for aList in falsePositives.values():
        retval.append(max(aList)[1])

    return retval


def analyzeInfo(resultList):
    cumulative = {"gsT":{}, "gsF":{}}

    # Accumulate information
    # For each run (cross validation)
    for result in resultList:
        infoMap = result.info
        # For each PDB in the Run
        for (pdb, data) in infoMap.items():
            if data["gsT"]:
                for (key, value) in data.items():
                    cumulative["gsT"][key] = cumulative["gsT"].get(key,0) + value
            if data["gsF"]:
                for (key, value) in data.items():
                    cumulative["gsF"][key] = cumulative["gsF"].get(key,0) + value
    # Combined info
    combined = {}
    for data in cumulative.values():
        for (key, value) in data.items():
            combined[key] = combined.get(key,0) + value
    for key in ("sTP", "sFN", "sTN", "sFP", "hTP", "hFP"):
        exec "%s = float(combined.get('%s', 0))" % (key, key)
    sens = 100*fraction(sTP,sFN)
    spec = 100*fraction(sTN,sFP)
    ppv = 100*fraction(hTP,hFP)
    
    # Calculate PPV
    for (key,data) in cumulative.items():
        tp = data.get("hTP", 0)
        fp = data.get("hFP", 0)


# =================================================================
# Performs one evaluation on each CV* directories
def evaluate(radius, cutoff, dirname, sitefile=None, crossvalidate="CV_*",
             scandir="Scan-test"):
    if crossvalidate:
        dirlist = glob.glob(os.path.join(dirname, crossvalidate))
    else:
        dirlist = [dirname]

    genSiteFile = 1
    if sitefile:
        genSiteFile = 0

    # get data for all runs
    data = []
    for runDir in dirlist:
        hitsDir = os.path.join(runDir, "%s/Hits" % scandir)
        if genSiteFile:
            sitefile = os.path.join(runDir, "test.site")
        results = evaluateRun(radius, cutoff, hitsDir, sitefile)
        data.append(results)

    # calculate performance
    (sTP, sFN, sTN, sFP, pTP, pFP) = accumulate(map(lambda d:d.tuple(), data))
    sens = fraction(sTP,sFN)
    spec = fraction(sTN,sFP)
    ppv = fraction(pTP,pFP)
    if PRINTFORMAT=='tabular':
        sys.stdout.write("%g\t%g\t%.4f\t%.4f\t%.4f" % (radius, cutoff, sens, spec, ppv))
        if VERBOSITY > 0:
            sys.stdout.write("\t%d\t%d\t%d" % (sTP+sFN, sTN+sFP, pTP+pFP))
        if VERBOSITY > 1:
            sys.stdout.write("\t%d\t%d\t%d\t%d\t%d\t%d" % (sTP,sFN,sTN,sFP,pTP,pFP))
        print
    else:
        print "RAD: %4g CUT: %4g SENS: %5.4f SPEC: %5.4f PPV: %5.4f" % (radius, cutoff, sens, spec, ppv),
        if VERBOSITY > 0:
            print "SITES: %4d nSITES: %4d HITS: %6d" % (sTP+sFN, sTN+sFP, pTP+pFP),
        if VERBOSITY > 1:
            print "sTP: %3d sFN: %3d sTN: %3d sFP: %3d pTP: %3d pFP: %3d" % (sTP,sFN,sTN,sFP,pTP,pFP),
        print
    
    return (sens, getFalsePositives(data))


def Usage(ec=None):
    print """
    Usage: %s [options] <dirname>
        Options:
        -r <radius>                site radius
        -l <lower cutoff>          cutoff lower limit
        -u <upper cutoff>          cutoff upper limit
        -i <increment>             cutoff increment
        -s <sitefile>              reference sitefile instead of relative
        -d                         don't perform cross validation
                                   <dirname> is scan directory
        -P <prefix>                cross validation directory prefix
        -n <sensitivity>           minimum sensitivity (0..1)
        -t <numfp>                 number of top false positives to display
        -v                         increase verbosity
        -b                         decrease verbosity
""" % sys.argv[0]
    if ec != None:
        sys.exit(ec)

def main():
    global VERBOSITY

    # default parameters
    radius = 7
    lower = 0
    upper = 120
    increment = 5
    dirname = None
    minsensitivity = 0.80
    topfalse = None
    getscandir = 0
    printHeader = 1

    evalKwArgs = {"crossvalidate":"CV*"}

    # process command line arguments
    args = sys.argv[1:]
    (opts, args) = getopt.getopt(args, "r:l:u:i:s:dP:n:t:vbH")
    for (opt, arg) in opts:
        if opt in ("-r",):
            radius = float(arg)
        elif opt in ("-l",):
            lower = float(arg)
        elif opt in ("-u",):
            upper = float(arg)
        elif opt in ("-i",):
            increment = float(arg)
        elif opt in ("-s",):
            evalKwArgs["sitefile"] = os.path.abspath(arg)
        elif opt in ("-d",):
            evalKwArgs["crossvalidate"] = None
        elif opt in ("-P",):
            evalKwArgs["crossvalidate"] = arg + "*"
        elif opt in ("-n",):
            minsensitivity = float(arg)
            # assume that it is percentage if greater than 1.0
            if minsensitivity > 1.0:
                minsensitivity = float(arg) / 100.0
        elif opt in ("-t",):
            topfalse = int(arg)
        elif opt in ("-v",):
            VERBOSITY += 1
        elif opt in ("-b",):
            VERBOSITY -= 1
        elif opt in ("-H",):
            printHeader = 0
    if len(args) != 1:
        Usage(2)
    dirname = args.pop(0)

    # Strip trailing '/'
    if dirname[-1] == '/':
        dirname = dirname[:-1]

    # dirname specifies scan directory if not performing cross validation
    if evalKwArgs["crossvalidate"] == None:
        (dirname, scandir) = os.path.split(dirname)
        if not dirname:
            dirname = "."
        evalKwArgs["scandir"] = scandir

    if printHeader and PRINTFORMAT=='tabular':
        sys.stdout.write("RAD\tCUT\tSENS\tSPEC\tPPV")
        if VERBOSITY > 0:
            sys.stdout.write("\tSITES\tnSITES\tHITS")
        if VERBOSITY > 1:
            sys.stdout.write("\tsTP\tsFN\tsTN\tsFP\tpTP\tpFP")
        sys.stdout.write("\n")
        

    # perform evaluation at each setting
    data = []
    for cutoff in arange(lower,upper+increment/2.0,increment):
        data.append(apply(evaluate,(radius, cutoff, dirname),evalKwArgs))

    if topfalse:
        MINSENS = minsensitivity
        # consider only false positives with at least minimum sensitivity
        fpList = filter(lambda (s,f),ms=MINSENS:s >= ms, data)
        try:
            (sens,fp) = min(fpList)
        except ValueError:
            (sens,fp) = (1.0, [])

        print "NEW NONSITES (Sensitivity=%5.4g)" % sens
        fp.sort(lambda a,b:cmp(b.score, a.score))
        for hit in fp[:topfalse]:
            print '("%s" X %g Y %g Z %g %s)  ;; %g' % \
                  (hit.pdbid, hit.x, hit.y, hit.z, "NIL", hit.score)


if __name__ == "__main__":
    main()
