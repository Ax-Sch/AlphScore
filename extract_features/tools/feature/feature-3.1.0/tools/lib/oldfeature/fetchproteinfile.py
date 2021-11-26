
# $Id: fetchproteinfile.py,v 1.5 2004/02/22 21:23:26 mliang Exp $
# Copyright 2001 (c) Mike Liang.  All Rights Reserved
#
# Fetches Protein Files
#

import re
import os
import shutil
import sys
import urllib
import callable
import fetchpdb
import fetchdssp


# ----------------------------------------------------------------------
def _fetchPDB(pdbid, filePath):
    return fetchpdb.fetchPDB(pdbid, filePath)



# ----------------------------------------------------------------------
def _fetchDSSP(pdbid, filePath):
    return fetchdssp.fetchDSSP(pdbid, filePath)



# ----------------------------------------------------------------------
class _FileType:
    def canonicalFilename(self, id):
        return id

    def localPath(self, id):
        dir = self.localDir
        filename = self.canonicalFilename(id)
        return os.path.join(dir, filename)

    def localPathGz(self, id):
        dir = self.localDir
        filename = self.canonicalFilename(id) + ".gz"
        return os.path.join(dir, filename)

    def localPathDiv(self, id):
        dir = self.localDir
        filename = self.canonicalFilename(id)
        id = id.lower()
        return os.path.join(dir, id[1:3], filename)

    def localPathDivGz(self, id):
        dir = self.localDir
        filename = self.canonicalFilename(id) + ".gz"
        id = id.lower()
        return os.path.join(dir, id[1:3], filename)

    def globalPath(self, id):
        dir = self.globalDir
        filename = self.canonicalFilename(id)
        return os.path.join(dir, filename)



# ----------------------------------------------------------------------
class _PDBFileType(_FileType):
    localDir = os.environ.get("LOCAL_PDB_DIR")
    if not localDir:
        localDir = os.environ.get("PDB_DIR","")
    globalDir = os.environ.get("PROTEIN_PDB_DIR","")
    fetchFile = callable.Callable(_fetchPDB)

    def canonicalFilename(self, id):
        return "pdb%s.ent" % id.lower()
PDBFileType = _PDBFileType()



# ----------------------------------------------------------------------
class _DSSPFileType(_FileType):
    localDir = os.environ.get("LOCAL_DSSP_DIR")
    if not localDir:
        localDir = os.environ.get("DSSP_DIR","")
    globalDir = os.environ.get("PROTEIN_DSSP_DIR","")
    fetchFile = callable.Callable(_fetchDSSP)

    def canonicalFilename(self, id):
        return "pdb%s.dssp" % id.lower()
DSSPFileType = _DSSPFileType()



# ----------------------------------------------------------------------
def _cachedFetch(id, filePath, fileType):
    """ Retrieves a file from cache or fetch it """

    # Check if OS supports symlink, if not
    if not hasattr(os, "symlink"):
        # Fetch file directly
        if not os.access(filePath, os.R_OK):
            fileType.fetchFile(id, filePath)
        return os.access(filePath, os.R_OK)

    paths = [filePath, fileType.localPath(id), fileType.globalPath(id)]
    # Filter for nonempty paths
    paths = filter(lambda p:p, paths)

    # Scan forward for existence of files
    for i in range(len(paths)):
        # if file exists,
        if os.access(paths[i], os.R_OK):
            # backward symlink
            for j in range(i,0,-1):
                if not os.path.exists(paths[j-1]):
                    os.symlink(paths[j], paths[j-1])
            # return
            return os.access(paths[0], os.R_OK)

    # Walk backward to fetch files
    for i in range(len(paths)-1,-1,-1):
        fileType.fetchFile(id,paths[i])
        # if fetch file successfully,
        if os.access(paths[i], os.R_OK):
            # backward symlink
            for j in range(i,0,-1):
                if not os.path.exists(paths[j-1]):
                    os.symlink(paths[j],paths[j-1])
            # return
            return os.access(paths[0], os.R_OK)

    return os.access(paths[0], os.R_OK)
 

# ----------------------------------------------------------------------
def _searchFetch(id, filePath, fileType):
    """ Retrieves a file from a search path or fetch it """

    # if file exists, return
    if os.access(filePath, os.R_OK):
        return os.access(filePath, os.R_OK)

    filenames = [fileType.localPath(id), fileType.localPathGz(id),
                 fileType.localPathDiv(id), fileType.localPathDivGz(id),
                 fileType.globalPath(id)]

    # Search filename list for existence
    for filename in filenames:
        # if file exists
        if os.access(filename, os.R_OK):
            # make a symlink or copy
            if hasattr(os, "symlink"):
                os.symlink(filename, filePath)
            else:
                shutil.copyfile(filename, filePath)
            return os.access(filePath, os.R_OK)

    # fetch file from URL
    fileType.fetchFile(id,filePath)
    return os.access(filePath, os.R_OK)


# ----------------------------------------------------------------------
def verifyOrCreateDir(dirpath):
    """ Verifies that dirpath is directory, or creates the directory if it does not exist.  Returns 0 if error, or 1 if okay
    """

    if os.path.exists(dirpath):
        if not os.path.isdir(dirpath):
            return 0
    else:
        try:
            os.makedirs(dirpath)
        except OSError:
            return 0
    return 1



# ----------------------------------------------------------------------
def fetchPDBFile(id, filepath, cache=1):
    if cache:
        return _cachedFetch(id, filepath, PDBFileType)
    else:
        return _searchFetch(id, filepath, PDBFileType)



# ----------------------------------------------------------------------
def fetchDSSPFile(id, filepath, cache=1):
    if cache:
        return _cachedFetch(id, filepath, DSSPFileType)
    else:
        return _searchFetch(id, filepath, DSSPFileType)



# ===================================================================
if __name__ == "__main__":
    if len(sys.argv[1:]) != 3:
        print "Usage: fetchproteinfile <type> <pdbid> <filename>"
        sys.exit(2)
    else:
        (filetype, pdbid, filename) = sys.argv[1:]
        if filetype == "pdb":
            _fetchPDB(pdbid, filename)
        elif filetype == "dssp":
            _fetchDSSP(pdbid, filename)
        else:
            print "Usage: fetchproteinfile <type> <pdbid> <filename>"
            sys.exit(2)
