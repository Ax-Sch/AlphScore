# $Id: proteindatabankdbh.py,v 1.4 2004/02/22 21:23:41 mliang Exp $
# Copyright 2001 (c) Mike Liang.  All rights reserved.

# Interface to PDB database

# ====================================================================
# IMPORTS
# ====================================================================
import sys,os
import PDB



# ====================================================================
# CONFIGURATION
# ====================================================================
DO_REMOTE = 0
PROTEIN_DATA_BANK_DB = os.environ.get("LOCAL_PDB_DIR")
if not PROTEIN_DATA_BANK_DB:
    PROTEIN_DATA_BANK_DB = os.environ.get("PDB_DIR",
                                      "/home/mliang/local/db/pdb/")


# ====================================================================
# FUNCTIONS
# ====================================================================
def findpdbfile(pdbid):
    pdbid = pdbid.lower()
    divid = pdbid[1:3]
    
    possiblenames = [ "pdb%(pdbid)s.ent",
                      "pdb%(pdbid)s.ent.gz",
                      "%(divid)s/pdb%(pdbid)s.ent",
                      "%(divid)s/pdb%(pdbid)s.ent.gz" ]
    
    for name in possiblenames:
        filename = os.path.join(PROTEIN_DATA_BANK_DB,name % vars())
        if os.path.exists(filename):
            return filename
    return None


# ====================================================================
# CLASSES
# ====================================================================
class ProteinDataBankDatabase:
    def __init__(self, database):
        self.database = database

    def getProtein(self, pdbid, filename=None):
        if filename == None:
            filename = findpdbfile(pdbid)
        if DO_REMOTE and filename == None:
            import fetchproteinfile
            import tempfile
            filename = tempfile.mktemp(".pdb")
            fetchproteinfile.fetchPDBFile(pdbid, filename)
        return PDB.readFile(filename)


# ====================================================================
def dbopen(database=None):
    if database == None:
        database = PROTEIN_DATA_BANK_DB
    if not os.path.isdir(database):
        return None
    return ProteinDataBankDatabase(database)


# ====================================================================
def setdb(database):
    global PROTEIN_DATA_BANK_DB
    PROTEIN_DATA_BANK_DB = database


# ====================================================================
def setremote(flag):
    global DO_REMOTE
    DO_REMOTE = flag
