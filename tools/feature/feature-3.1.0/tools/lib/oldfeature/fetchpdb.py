#! /bin/env python

# Fetches PDB file from RCSB

import sys
import urllib
import htmllib
import formatter


URL="http://www.rcsb.org/pdb/cgi/export.cgi?pdbId=%s;format=PDB"



# ===================================================================
def fetchPDB(pdbid, filename):
    # Get Data from Web
    try:
        data = urllib.urlopen(URL % pdbid).read()
    except IOError:
        return 0

    if not data:
        return 0

    # Save Data to file
    if filename and filename != "-":
        try:
            file = open(filename, "w")
            file.write(data)
            file.close()
        except:
            return 0
    else:
        sys.stdout.write(data)

    return 1



# ===================================================================
if __name__ == "__main__":
    try:
        apply(fetchPDB, sys.argv[1:])
    except TypeError:
        print "Usage: fetchpdb <pdbid> <filename>"
