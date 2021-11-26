#! /usr/bin/env python

# $Id: findCalciums.py,v 1.5 2004/02/22 21:22:46 mliang Exp $
# Copyright (c) 2001 Mike Liang.  All rights reserved.

# Finds all the calciums and creates sitefile

import sys
import os
import glob
import getopt
import proteindatabankdbh
from sitefile import SiteFile, Site


pdbdbh = proteindatabankdbh.dbopen()

def isSiteFile(filename):
    (base,ext) = os.path.splitext(filename)
    if ext == ".site":
        return 1
    return 0

def isPdbIdFile(filename):
    (base,ext) = os.path.splitext(filename)
    if ext == ".ids":
        return 1
    return 0

def getpdbid(filename):
    (head, tail) = os.path.split(filename)
    (base, ext) = os.path.splitext(tail)
    return base

def atomToSite(pdbAtom):
    site = Site()
    (site.x, site.y, site.z) = pdbAtom.getLocation()
    site.type = Site.SITE_TYPE
    return site

def extractCalciums(pdbid):
    print >>sys.stderr, "Extracting: '%s'" % pdbid
    protein = pdbdbh.getProtein(pdbid)
    if not protein:
        return []
    
    caResidues = filter(lambda r: r.resName == "CA", protein.getResidues())
    caAtoms = map(lambda r:r.getAtom("CA"), caResidues)
    sites = map(atomToSite, caAtoms)
    for site in sites:
        site.pdbid = pdbid
    return sites

def main():
    args = sys.argv[1:]

    if len(args) != 1:
        print "Usage: %s <scandirname|filename>" % sys.argv[0]
        sys.exit(2)

    name = args[0]
    if isSiteFile(name):
        siteFile = SiteFile(name)
        pdbids = siteFile.sites.keys()
        siteFile.sites = {}
    else:
        if isPdbIdFile(name):
            pdbids = [line.strip() for line in open(name).readlines()]
        else:
            if os.path.isdir(name):
                filenames = glob.glob(os.path.join(name,"*.hits"))
            else:
                filenames = [name]
            pdbids = map(getpdbid, filenames)

        siteFile = SiteFile()
        siteFile.name = "extractedCalciums"
        siteFile.radius = 0
        siteFile.filename = "<%s>" % sys.argv[0]


    pdbids.sort()
    for pdbid in pdbids:
        siteFile.addSites(extractCalciums(pdbid))
    sys.stdout.write(str(siteFile))

if __name__ == "__main__":
    main()
