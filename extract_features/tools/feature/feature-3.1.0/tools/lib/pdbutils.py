#! /usr/bin/env python

# PDB utilities

import os
import re
import logutils
from utils import uniq

_PDB_DIR = os.environ.get('PDB_DIR','').split(os.pathsep)
_DSSP_DIR = os.environ.get('DSSP_DIR','').split(os.pathsep)

if not _PDB_DIR:
    logutils.warning("PDB_DIR environment not set")
_DSSP_DIR = os.environ.get('DSSP_DIR',_PDB_DIR)

def pdbfilename(pdbid):
    pdbid = pdbid.lower()
    return os.path.join(_PDB_DIR[0],pdbid[1:3],'pdb%s.ent.gz' % pdbid)

def pdbSeqresFilename():
    return os.path.join(_PDB_DIR[0],"../derived_data/pdb_seqres.txt")

def pdbidFromFilename(filename):
    root,ext = os.path.splitext(os.path.basename(filename))
    m = re.search(r'\d\w{3}',root)
    if m:
        return m.group()
    return root

def goodpdbfilename(filename, pdbdirList=None):
    # first check if pdbid is really a file
    if os.path.exists(filename):
        return filename

    # extract pdbid from pdbid
    if pdbdirList is None:
        pdbdirList = _PDB_DIR
    if type(pdbdirList) == str:
        pdbdirList = pdbdirList.split(os.pathsep)
    pdbid = pdbidFromFilename(filename)
    pdbidl = pdbid.lower()
    branch = pdbidl[1:3]

    # generate filename variants
    basenames = [x%vars() for x in ( "%(filename)s", "%(pdbid)s", "%(pdbidl)s", "pdb%(pdbidl)s" )]
    basenames = uniq(basenames)
    extensions = ( "", ".pdb", ".ent", ".FULL" )
    compressions = ( "", ".gz", ".Z" )

    # generate subdirectory locations
    subdirs = []
    for pdbdir in pdbdirList:
        innerSubdirs = [x%vars() for x in (
            "",
            "%(branch)s",
            "%(pdbdir)s",
            os.path.join("%(pdbdir)s","%(branch)s"),
            os.path.join("%(pdbdir)s","divided","%(branch)s"),
            os.path.join("%(pdbdir)s","data","structures","divided","pdb","%(branch)s")
        )]
        subdirs.extend(innerSubdirs)
    subdirs = uniq(subdirs)

    # search tree
    for subdir in subdirs:
        for cmp in compressions:
            for base in basenames:
                for ext in extensions:
                    filename = os.path.join(subdir,"%(base)s%(ext)s%(cmp)s" % vars())
                    if os.path.exists(filename):
                        return filename

    return None

def gooddsspfilename(filename, dsspdirList=None):
    # first check if pdbid is really a file
    if os.path.exists(filename):
        return filename

    # extract pdbid from pdbid
    if dsspdirList is None:
        dsspdirList = _DSSP_DIR
    if type(dsspdirList) == str:
        dsspdirList = _DSSP_DIR.split(os.pathsep)
    pdbid = pdbidFromFilename(filename)
    pdbidl = pdbid.lower()
    branch = pdbidl[1:3]

    # generate filename variants
    basenames = [x%vars() for x in ( "%(filename)s", "%(pdbid)s", "%(pdbidl)s", "pdb%(pdbidl)s" )]
    basenames = uniq(basenames)
    extensions = ( "", ".dssp", ".DSSP" )
    compressions = ( "", ".gz", ".Z" )

    # generate subdirectory locations
    subdirs = []
    for dsspdir in dsspdirList:
        innerSubdirs = [x%vars() for x in (
            "",
            "%(branch)s",
            "%(dsspdir)s",
            os.path.join("%(dsspdir)s","%(branch)s"),
            os.path.join("%(dsspdir)s","divided","%(branch)s"),
            os.path.join("%(dsspdir)s","data","structures","divided","pdb","%(branch)s")
        )]
        subdirs.extend(innerSubdirs)
    subdirs = uniq(subdirs)

    # search tree
    for subdir in subdirs:
        for cmp in compressions:
            for base in basenames:
                for ext in extensions:
                    filename = os.path.join(subdir,"%(base)s%(ext)s%(cmp)s" % vars())
                    if os.path.exists(filename):
                        return filename

    return None
