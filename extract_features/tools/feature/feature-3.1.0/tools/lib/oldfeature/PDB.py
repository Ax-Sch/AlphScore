#!/usr/pubsw/bin/python
# $Id: PDB.py,v 1.7 2003/04/18 13:12:40 mliang Exp $
# PDB file reader
# Copyright (c) 2000 Mike Liang. All rights reserved.

# Status:
# - handles multiple models
# - handles only ATOM/HETATM records
# - ignores alternate atom locations
# - residues assumed to be amino acids

# ====================================================================
# IMPORTS
# ====================================================================
import sys, os
import string
from fsutil import gzopen



# ====================================================================
# CONFIGURATION
# ====================================================================
DEFAULT_MAX_MODELS = 1



# ====================================================================
# GLOBALS
# ====================================================================
READALL_MAX_MODELS = -1

# Map from Residue Abbreviation to Alphabet
resLetters = { "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "ASX": "B",
               "CYS": "C", "GLN": "Q", "GLU": "E", "GLX": "Z", "GLY": "G", 
               "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M",
               "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
               "TYR": "Y", "VAL": "V" }



# ====================================================================
# FUNCTIONS
# ====================================================================
# Residue abbreviation to letter map
def residueLetter(abbreviation):
    return resLetters.get(abbreviation,"?")


# ===========================================================================
def residueListSequence(aList):
    """ Converts a list of 3 letter AA to string """
    return string.join(map(residueLetter,aList),"")


# ====================================================================
# PDB-derived model of protein
# Protein Class
class Protein:
    # contains a list of models
    def __init__(self, filename = None):
        self.models = []
        self.filename = filename
        self.expmethod = "X-ray Diffraction"
        self.pdbID = ""
        if self.filename:
            self._setpdbid()

    def setExperimentalMethod(self, method):
        self.expmethod = method

    def isXrayExp(self):
        return self.expmethod.upper().find("X-RAY") != -1

    def addModel(self, model):
        self.models.append(model)

    def getAtoms(self, modelNo=0):
        return self.models[modelNo].getAtoms()

    def getResidues(self, modelNo=0):
        return self.models[modelNo].getResidues()

    def getChains(self, modelNo=0):
        return self.models[modelNo].getChains()

    def getChain(self, chainid, modelNo=0):
        return self.models[modelNo].getChain(chainid)
        
    def displayInfo(self):
        print "Protein Name: %s" % self.filename
        print "Number of models: %d" % len(self.models)
        map(Model.displayInfo, self.models)

    def _setpdbid(self):
        base = os.path.basename(self.filename)
        if base[:3] == "pdb" and base[-4:] == ".ent":
            self.pdbID = base[3:7].upper()
        else:
            self.pdbID = base[:4].upper()


# ====================================================================
# Model Class
class Model:
    # contains a list of chains
    def __init__(self):
        self.chains = []
        self.atoms = []

    def addChain(self, chain):
        self.chains.append(chain)

    def addAtom(self, atom):
        self.atoms.append(atom)

    def getChains(self):
        return self.chains

    def getChain(self,chainid):
        for chain in self.chains:
            if chain.chainID == chainid:
                return chain
        return None

    def getResidues(self):
        residues = []
        for chain in self.chains:
            residues.extend(chain.getResidues())
        return residues

    def getAtoms(self):
        return self.atoms

    def displayInfo(self):
        print "Number of chains: %d" % len(self.chains)
        map(Chain.displayInfo, self.chains)


# ====================================================================
# Chain Class
class Chain:
    # contains a list of residues
    def __init__(self, chainID = None):
        self.residues = []
        self.chainID = chainID
        self.sequence = ""

    def addResidue(self, residue):
        residue.chain = self
        residue.index = len(self.sequence)

        self.residues.append(residue)
        self.sequence = self.sequence + residueLetter(residue.resName)

    def getResidues(self):
        return self.residues

    def getResidue(self, resSeq, iCode=""):
        for residue in self.residues:
            if (residue.resSeq, residue.iCode) == (resSeq, iCode):
                return residue
        return None

    def getAtoms(self):
        return reduce(lambda a,b:a+b,[r.getAtoms() for r in self.residues])

    def getSequence(self):
        return residueListSequence(map(lambda r:r.resName,self.residues))


    def displayInfo(self):
        print "ChainID: %s" % self.chainID
        print "Number of residues: %d" % len(self.residues)
        print "Sequence: %s" % self.sequence


# ====================================================================
# Residue Class
class Residue:
    # contains a hash of atoms
    def __init__(self, resName = None, resSeq = None, iCode = None, chainID = None):
        self.atomHash = {}
        self.resName = resName
        self.resSeq = resSeq
        self.iCode = iCode
        self.chainID = chainID

        self.chain = None
        self.index = None

    def addAtom(self, atom):
        if not self.atomHash.has_key(atom.name):
            self.atomHash[atom.name]=atom
            atom.residue = self

    def displayInfo(self):
        print self

    def __str__(self):
        return "%s(%d%s:%s) atoms: %d" % (self.resName, self.resSeq,
                                          self.iCode, self.chainID,
                                          len(self.atomHash))

    def getAtoms(self):
        return self.atomHash.values()

    def getAtom(self, atomName):
        return self.atomHash.get(atomName)


# ====================================================================
# Atom Class
class Atom:
    def __init__(self, line):
        self.residue = None
        
        self.recName = string.strip(line[0:6])
        self.serialNo = string.atoi(line[6:11])
        self.name = string.strip(line[12:16])
        self.altLoc = string.strip(line[16:17])

        self.resName = string.strip(line[17:20])
        self.chainID = string.strip(line[21:22])
        self.resSeq = string.atoi(line[22:26])
        self.iCode = string.strip(line[26:27])

        self.x = string.atof(line[30:38])
        self.y = string.atof(line[38:46])
        self.z = string.atof(line[46:54])
        self.occupancy = string.atof(line[54:60])
        self.tempFactor = string.atof(line[60:66])
        self.segID = string.strip(line[72:76])
        self.element = string.strip(line[76:78])
        self.charge = string.strip(line[78:80])

    def getLocation(self):
        return (self.x, self.y, self.z)

    def getFullResidueName(self):
        return "%s%d%s:%s" % (self.resName, self.resSeq,
                              self.iCode, self.chainID)
    
    def getFullAtomName(self):
        return "%s.%s" % (self.getFullResidueName(), self.name)

    def __str__(self):
        return "%-6s%5d %-4s%1s%-3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s" \
               % (self.recName, self.serialNo, self.name, self.altLoc,
                  self.resName, self.chainID, self.resSeq, self.iCode,
                  self.x, self.y, self.z, self.occupancy, self.tempFactor,
                  self.segID, self.element, self.charge)



# ====================================================================
# PDB File parser
def readFile(filename=None, max_models = DEFAULT_MAX_MODELS, fobj=None):
    do_close = 0
    if fobj == None and filename == None:
        raise TypeError("Must specify either file or filename")
    if not fobj and filename:
        try:
            fobj = gzopen(filename)
        except IOError,e:
            print >>sys.stderr, e
            return None
        do_close = 1
    if fobj and not filename:
        filename = fobj.name


    curProtein = Protein(filename)
    curModel = None
    curChain = None
    curResidue = None

    for line in fobj.readlines():
        recType = string.strip(line[0:6])

        # Handle Coordinate Records
        if recType == "END":
            break
        elif recType == "EXPDTA":
            method = line[10:70].strip()
            curProtein.setExperimentalMethod(method)
        elif recType == "MODEL":
            curModel = None
        elif recType == "ENDMDL":
            curModel = None
            if max_models != READALL_MAX_MODELS \
               and len(curProtein.models) >= max_models:
                break
        elif recType == "TER":
            curChain = None
        elif recType == "ATOM" or recType == "HETATM":
            atom = Atom(line)

            # create model if necessary
            if curModel == None:
                curModel = Model()
                curProtein.addModel(curModel)
                curChain = None
            # create chain if necessary
            if curChain == None or atom.chainID != curChain.chainID:
                curChain = Chain(atom.chainID)
                curModel.addChain(curChain)
                curResidue = None
            # create residue if necessary
            if curResidue == None or (atom.resSeq != curResidue.resSeq or \
                                      atom.iCode != curResidue.iCode):
                curResidue = Residue(atom.resName, atom.resSeq, atom.iCode, \
                                     atom.chainID)
                curChain.addResidue(curResidue)

            # add atom
            curResidue.addAtom(atom)
            curModel.addAtom(atom)
        else:
            pass

    if do_close:
        fobj.close()
    return curProtein



polarResidues = (
    )
def isPolar(atom):
    pass


chargedResidues = {
    }
def isCharged(atom):
    pass
