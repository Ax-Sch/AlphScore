#! /bin/env python

# Fetches DSSP file from SRS

import sys
import urllib
import re



URL="http://srs.ebi.ac.uk/srs6bin/cgi-bin/wgetz?-e+[DSSP:'%s']"



# ===================================================================
def fetchDSSP(pdbid, filename):
    # Get Data from Web
    try:
        data = urllib.urlopen(URL % pdbid).read()
    except IOError:
        return 0

    # If Error
    if data.find("SRS error") != -1:
        return 0

    # Strip HTML tags/entities
    data = re.sub("</?pre>","",data)[:-1]

    # Save data to File
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
        apply(fetchDSSP, sys.argv[1:])
    except TypeError:
        print "Usage: fetchdssp <pdbid> <filename>"

