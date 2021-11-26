# $Id: webload.py,v 1.4 2004/11/11 11:13:10 mliang Exp $
# Copyright (c) 2004 Mike Liang. All rights reserved.
#
# PyMol plugin (requires PyMol version 0.96 or above)
# Provides command line extension for fetching PDB structures from the web
# Use [ Plugin | Install Plugin... ] to install the webload plugin
#
# For versions 0.95 and earlier, include the following in .pymolrc
#     run /path/to/webload/webload.py

from pymol import cmd

_PDB_URL='http://www.rcsb.org/pdb/cgi/export.cgi?format=PDB&pdbId=%(pdbid)s&compression=gz'

def __init__(self=None):
    cmd.extend('webload',webload)

def webload(pdbid):
    import sys,os
    import StringIO
    import gzip
    import urllib

    label = pdbid
    url = _PDB_URL % vars()
    try:
        datafh = StringIO.StringIO(urllib.urlopen(url).read())
        data = gzip.GzipFile(fileobj=datafh).read()
    except Exception, e:
        print "Error: pdbid=%(pdbid)s" % vars(), e
        return
        
    if not data:
        print "Invalid PDB id %(pdbid)s" % vars()
        return
 
    # load it into pymol
    cmd.read_pdbstr(data,label)

if __name__ == 'pymol':
    __init__()
