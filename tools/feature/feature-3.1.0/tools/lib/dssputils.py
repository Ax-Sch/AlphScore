#! /usr/bin/env python

# DSSP utilities

import os
import logutils

_env_vars = [ 'DSSP_DIR', 'LOCAL_DSSP_DIR' ]
for _env_var in _env_vars:
    _DSSP_DIR = os.environ.get(_env_var,'')
    if _DSSP_DIR:
        break
else:
    logutils.warning("None of",_env_vars,"environment variables set")

def gooddsspfilename(pdbid, dsspdir=None):
    # first check if pdbid is really a file
    if os.path.exists(pdbid):
        return pdbid

    # extract pdbid from pdbid
    if dsspdir is None:
        dsspdir = _DSSP_DIR
    pdbid = pdbid[:4]
    pdbidl = pdbid.lower()
    branch = pdbidl[1:3]

    # generate filename variants
    basenames = [x%vars() for x in ( "%(pdbid)s", "pdb%(pdbidl)s" )]
    extensions = ( "", ".dssp", ".DSSP" )
    compressions = ( "", ".gz", ".Z" )

    # generate subdirectory locations
    subdirs = [x%vars() for x in (
        "",
        "%(dsspdir)s",
        os.path.join("%(dsspdir)s","%(branch)s"),
        os.path.join("%(dsspdir)s","divided","%(branch)s"),
    )]

    # search tree
    for subdir in subdirs:
        for cmp in compressions:
            for base in basenames:
                for ext in extensions:
                    filename = os.path.join(subdir,"%(base)s%(ext)s%(cmp)s" % vars())
                    if os.path.exists(filename):
                        return filename

    return None
