#! /usr/bin/env python

# Find distribution of score for sites and nonsites

import sys
import os
from sitefile import SiteFile
from hitfile import HitFile
from nngrid import Grid, distance
from mannwhitney import ranksumtestAppx


# ================================================================
def Usage(error=None):
    print "Usage: %s <sitefile> <hitdir> <radius>" % sys.argv[0]
    if error != None:
        sys.exit(error)


# ================================================================
def main(sitefilename, hitdirname, radius):
    # load sitefile
    sfile = SiteFile(sitefilename)
    pdbids = sfile.sites.keys()
    pdbids.sort()

    # For each PDB
    for pdbid in pdbids:
        hitfilename = os.path.join(hitdirname, "%s.hits" % pdbid.upper())
        # load hitfile
        hitfile = HitFile(hitfilename)

        # create query data structure
        grid = Grid(radius/2.0)
        for hit in hitfile.hits:
            grid.add(hit.getLocation(),hit)

        # for each site of PDB,
        for site in sfile.sites[pdbid]:
            # find all hits around site
            neighbors = grid.query(site.getLocation(), radius, dataOnly=1)
            # display pdbid, location, type, score, distance from site
            for hit in neighbors:
                print "%s %g %g %g %s %g %g" % \
                      (pdbid, hit.x, hit.y, hit.z, site.type, hit.score,
                       distance(hit.getLocation(), site.getLocation()))


# ================================================================
def main2(sitefilename, hitdirname, radius):
    # load sitefile
    sfile = SiteFile(sitefilename)
    pdbids = sfile.sites.keys()
    pdbids.sort()


    # For each PDB
    allHits = {}
    for pdbid in pdbids:
        hitfilename = os.path.join(hitdirname, "%s.hits" % pdbid.upper())
        # load hitfile
        hitfile = HitFile(hitfilename)

        # for each hit,
        neighbors = {}
        for hit in hitfile.hits:
            # find if near site
            type = "OUT"
            for site in sfile.sites[pdbid]:
                dist = distance(hit.getLocation(), site.getLocation())
                if dist < radius:
                    type = site.type
                    break
            if type == "OUT":
                dist = -1
            neighbors.setdefault(type,[]).append((hit, dist))
                

        # display pdbid, location, type, score, distance from site
        for (type,hits) in neighbors.items():
            for (hit, dist) in hits:
                allHits.setdefault(type,[]).append(hit.score)
                # print "%s %g %g %g %s %g %g" % \
                # (pdbid, hit.x, hit.y, hit.z, type, hit.score, dist)

    # Analyze score distribution
    def compare(type1, type2, alpha=0.001):
        opts = []
        if ranksumtestAppx(allHits.get(type1),allHits.get(type2),alpha,opts):
            result = "DIFFERENT"
        else:
            result = "SIMILAR"
        if opts:
            print "%s and %s %s (z=%g p=%g alpha=%g)" % \
                  (type1, type2, result, opts[0], opts[1], alpha)
        else:
            print "%s and %s Not Comparable" % (type1, type2)

    compare('T','NIL')
    compare('T','OUT')
    compare('NIL','OUT')


# ================================================================
if __name__ == "__main__":
    # Process command line arguments
    args = sys.argv[1:]
    if len(args) != 3:
        Usage(2)
    sitefilename = args[0]
    hitdirname = args[1]
    radius = float(args[2])

    main2(sitefilename, hitdirname, radius)
