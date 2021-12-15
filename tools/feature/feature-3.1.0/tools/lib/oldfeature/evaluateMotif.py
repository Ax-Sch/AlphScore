#! /usr/bin/env python

# $Id: evaluateMotif.py,v 1.4 2004/02/22 21:22:46 mliang Exp $
# Copyright (c) 2001 Mike Liang.  All rights reserved.
#
# Evaluate EFHand sequence motifs and FEATURE hits
#

# Description:
#
# Gold Standard: Presence of Calcium in proximity of fragment
#
# Sequences:
#   True positive: motif match has calcium within 5 Angstroms of alpha-carbon
#   False positive: motif match does not have calcium within 5 Angstroms
#   False negative: no motif match near calcium
#   True negative: ???
#
# Feature:
#   True positive: hiscore within radius of calcium
#   False positive: hiscore not within radius of calcium
#   False negative: no hiscore within radius of calcium
#   True negative: ???
#   Ignore hiscore hits not within radius of calcium but within radius
#   of motif match


# -----------------------------------------------------------------------
# Imports
import sys
import os.path
import glob
import time
import getopt
import string

import Motif
from fsutil import *
from featureexts import *
from nngrid import Grid
from sitefile import SiteFile
from hitfile import HitFile
import proteindatabankdbh


# -----------------------------------------------------------------------
# Global Constants
# FEATURE Filename Extensions

STDOUT = sys.stdout



# -----------------------------------------------------------------------
# Functions
def HiScores(aList, cutoff):
    return filter(lambda x:x.score >= cutoff, aList)

def Usage(ec=None):
    print """
    Usage: %s [options] <radius> [<cutoff>] <scandir> <sitefile>
        Options:
        -m motif only
        -M <patternfile>
        -O <outputdir>
""" % (sys.argv[0],)
    if ec != None:
        sys.exit(ec)

def siteOnly(aList):
    return filter(lambda s:s.isSite(), aList)
    
def notsiteOnly(aList):
    return filter(lambda s:not s.isSite(), aList)


# -----------------------------------------------------------------------
# Input:
hitsPath = None
patternfile = None
radius = None
cutoff = None
outputdir = None
motifonly = 0

args = sys.argv[1:]
optlist, args = getopt.getopt(args, "M:O:m")
for (option, value) in optlist:
    if option in ("-m",):
        motifonly = 1
    elif option in ("-M",):
        if value != "-":
            patternfile = value
        else:
            patternfile = None
    elif option in ("-O",):
        if value != "-":
            outputdir = value
        else:
            outputdir = None

if not motifonly and len(args) == 4:
    radius = float(args.pop(0))
    cutoff = float(args.pop(0))
    scandir = args.pop(0)
    sitefilename = args.pop(0)
elif motifonly and len(args) == 3:
    radius = float(args.pop(0))
    cutoff = 0.0
    scandir = args.pop(0)
    sitefilename = args.pop(0)
else:
    Usage(2)


if outputdir:
    outputfile = os.path.join(outputdir, "eval.r%g.c%g.txt" % (radius, cutoff))
    sys.stdout = open(outputfile,"w")

hitsPath = os.path.join(scandir, "Hits")


# Get list of PDBs
hitfilelist = glob.glob(os.path.join(hitsPath,"*"+HITSFILE_EXT))
pdblist = [os.path.splitext(os.path.basename(x))[0] for x in hitfilelist]
pdblist.sort()

# Load motifs
if patternfile:
    patternlist = [x.strip() for x in open(patternfile).readlines()]
else:
    patternlist = []




class SummaryAnalysis:
    attributes = [
        "num.Sites.GP",
        "num.Motif.Sites.TP",
        "num.Feature.Sites.TP",
        "num.Both.Sites.TP",
        "num.Neither.Sites.TP",

        "num.Sites.GN",
        "num.Motif.Sites.FP",
        "num.Feature.Sites.FP",
        "num.Both.Sites.FP",
        "num.Neither.Sites.FP",

        "num.Motif.Hits",
        "num.Motif.Hits.TP",
        "num.Motif.Hits.FP",

        "num.Feature.Hits",
        "num.Feature.Hits.TP",
        "num.Feature.Hits.TP.Motif",
        "num.Feature.Hits.TP.NotMotif",

        "num.Feature.Hits.FP",
        "num.Feature.Hits.FP.Motif",
        "num.Feature.Hits.FP.NotMotif"
        ]
    
    def __init__(self, pdbid = None, cutoff = 0.0, siteList = [], matchlist = [], hitlist = []):
        self.pdbid = pdbid
        self.cutoff = cutoff
        self.allsitesList = siteList
        self.matchlist = matchlist
        self.hitlist = hitlist
        self.properties = {}

        self.siteList = siteOnly(self.allsitesList)
        self.notsiteList = notsiteOnly(self.allsitesList)


        if self.siteList:
            self.properties['num.Sites.GP'] = len(self.siteList)
            list1 = filter(lambda x:x.motifHits, self.siteList)
            self.properties['num.Motif.Sites.TP'] = len(list1)
            list2 = filter(lambda x,c=self.cutoff:HiScores(x.featureHits,c), self.siteList)
            self.properties['num.Feature.Sites.TP'] = len(list2)
            list = filter(lambda x,l1=list1,l2=list2:x in l1 and x in l2, self.siteList)
            self.properties['num.Both.Sites.TP'] = len(list)
            list = filter(lambda x,l1=list1,l2=list2:x not in l1 and x not in l2, self.siteList)
            self.properties['num.Neither.Sites.TP'] = len(list)

        if self.notsiteList:
            self.properties['num.Sites.GN'] = len(self.notsiteList)
            list1 = filter(lambda x:x.motifHits, self.notsiteList)
            self.properties['num.Motif.Sites.FP'] = len(list1)
            list2 = filter(lambda x,c=self.cutoff:HiScores(x.featureHits,c), self.notsiteList)
            self.properties['num.Feature.Sites.FP'] = len(list2)
            list = filter(lambda x,l1=list1,l2=list2:x in l1 and x in l2, self.notsiteList)
            self.properties['num.Both.Sites.FP'] = len(list)
            list = filter(lambda x,l1=list1,l2=list2:x not in l1 and x not in l2, self.notsiteList)
            self.properties['num.Neither.Sites.FP'] = len(list)
            

        if self.matchlist:
            self.properties['num.Motif.Hits'] = len(self.matchlist)
            list = filter(lambda x:siteOnly(x.siteHits), self.matchlist)
            self.properties['num.Motif.Hits.TP'] = len(list)
            list = filter(lambda x:not siteOnly(x.siteHits), self.matchlist)
            self.properties['num.Motif.Hits.FP'] = len(list)

        if self.hitlist:
            hiScoreHits = HiScores(self.hitlist,self.cutoff)
            self.properties['num.Feature.Hits'] = len(hiScoreHits)
            list1 = filter(lambda x:siteOnly(x.siteHits), hiScoreHits)
            self.properties['num.Feature.Hits.TP'] = len(list1)
            list2 = filter(lambda x:x.motifHits, list1)
            self.properties['num.Feature.Hits.TP.Motif'] = len(list2)
            list2 = filter(lambda x:not x.motifHits, list1)
            self.properties['num.Feature.Hits.TP.NotMotif'] = len(list2)
            list1 = filter(lambda x:not siteOnly(x.siteHits), hiScoreHits)
            self.properties['num.Feature.Hits.FP'] = len(list1)
            list2 = filter(lambda x:x.motifHits, list1)
            self.properties['num.Feature.Hits.FP.Motif'] = len(list2)
            list2 = filter(lambda x:not x.motifHits, list1)
            self.properties['num.Feature.Hits.FP.NotMotif'] = len(list2)

    def getProp(self, attribute):
        return self.properties.get(attribute,0)

    def accumulate(self, other):
        for attribute in SummaryAnalysis.attributes:
            self.properties[attribute] = self.getProp(attribute) + other.getProp(attribute)

    def displayHeadings(self):
        print "pdbid",
        for attribute in SummaryAnalysis.attributes:
            print attribute,
        print

    def display(self):
        print self.pdbid,
        for attribute in SummaryAnalysis.attributes:
            print self.getProp(attribute),
        print

    def displaySummary(self):
        print "Motif:",
        print "Sens:",
        try: sens = float(self.getProp('num.Motif.Sites.TP'))/self.getProp('num.Sites.GP')
        except ZeroDivisionError: sens = 0.0
        print "%5.4f" % sens,
        print "Spec:",
        try: spec = 1.0-float(self.getProp('num.Motif.Sites.FP'))/self.getProp('num.Sites.GN')
        except ZeroDivisionError: spec = 0.0
        print "%5.4f" % spec,
        print "PPV:",
        try: ppv = float(self.getProp('num.Motif.Hits.TP'))/self.getProp('num.Motif.Hits')
        except ZeroDivisionError: ppv = 0.0
        print "%5.4f" % ppv

        print "Feature:",
        print "Sens:",
        try: sens = float(self.getProp('num.Feature.Sites.TP'))/self.getProp('num.Sites.GP')
        except ZeroDivisionError: sens = 0.0
        print "%5.4f" % sens,
        print "Spec:",
        try: spec = 1.0-float(self.getProp('num.Feature.Sites.FP'))/self.getProp('num.Sites.GN')
        except ZeroDivisionError: spec = 0.0
        print "%5.4f" % spec,
        print "PPV:",
        try: ppv = float(self.getProp('num.Feature.Hits.TP'))/self.getProp('num.Feature.Hits')
        except ZeroDivisionError: ppv = 0.0
        print "%5.4f" % ppv

    
        
# -----------------------------------------------------------------------
# main()

# Generate Motif object from patterns
motiflist = map(lambda x:Motif.Motif(x), patternlist)

globalAnalysis = SummaryAnalysis('TOTAL')
globalAnalysis.displayHeadings()

# Load Sitefile
siteFile = SiteFile(sitefilename)

# Open PDB dbh
pdbdbh = proteindatabankdbh.dbopen()


# For each PDB
for pdbid in pdblist:
    print >>sys.stderr, pdbid

    if motiflist:
        clock = time.clock()
        # Load protein
        protein = pdbdbh.getProtein(pdbid)
        print >>sys.stderr, "Loading protein:", time.clock()-clock

    # Get Sites
    siteList = siteFile.sites.get(pdbid,None)
    if siteList == None:
        print >>sys.stderr, "Warning:",pdbid,"Does not have any Sites!"
        siteList = []
    for site in siteList:
        site.motifHits = []
        site.featureHits = []
    
    clock = time.clock()
    # Locate Motifs <-> Sites
    matchlist = []
    for motif in motiflist:
        # find motif match
        matches = motif.scanProtein(protein).matches
        matchlist.extend(matches)
        for match in matches:
            # get alpha-carbons
            match.alphas = [x.getAtom("CA") for x in match.residues]
            match.alphas = filter(lambda x:x, match.alphas)
            alphaPoints = [atom.getLocation() for atom in match.alphas]
            match.centroid = CalculateCentroid(alphaPoints)
            # identify nearby calcium
            siteHits = GetNearbyObjects(match.alphas, siteList, radius)
            # link them up
            match.featureHits = []
            match.siteHits = siteHits
            for site in match.siteHits:
                site.motifHits.append(match)
    print >>sys.stderr, "Locating motifs<->sites:", time.clock()-clock


    if not motifonly:
        clock = time.clock()
        # Locate Hits <-> Sites
        hitsfilename = os.path.join(hitsPath, pdbid+HITSFILE_EXT)
        hitfile = HitFile(hitsfilename)
        hitlist = hitfile.hits
        grid = Grid(radius/2)
        for hit in hitlist:
            grid.add(hit.getLocation(),hit)
            hit.siteHits = []
            hit.motifHits = []
        for site in siteList:
            # find all hits near calcium
            datalist = grid.query(site.getLocation(), radius)
            # link them up
            for (point,hit) in datalist:
                hit.siteHits.append(site)
                site.featureHits.append(hit)
        print >>sys.stderr, "Locating hits<->sites:", time.clock()-clock
    else:
        hitlist = []
    

    if not motifonly:
        clock = time.clock()
        # Locate Hits <-> Motifs
        for match in matchlist:
            datalist = grid.query(match.centroid,radius)
            for (point,hit) in datalist:
                hit.motifHits.append(match)
                match.featureHits.append(hit)
        print >>sys.stderr, "Locating hits<->motif:", time.clock()-clock


    clock = time.clock()
    # Summary
    localAnalysis = SummaryAnalysis(pdbid, cutoff, siteList, matchlist, hitlist)
    print >>sys.stderr, "Summary:", time.clock()-clock

    localAnalysis.display()

    globalAnalysis.accumulate(localAnalysis)
    
            
# Summaries:
print "================================================================="

print "Patterns: %d" % len(patternlist)
for pattern in patternlist:
    print pattern
print "PDBs: %d" % len(pdblist)
print string.join(pdblist," ")

globalAnalysis.displayHeadings()
globalAnalysis.display()
globalAnalysis.displaySummary()
print "Radius:",radius,"Cutoff:",cutoff

