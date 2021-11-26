#! /usr/bin/env python

# Simple PDB parser
# extracts atom coordinates
# uses only first model
# ignores hetatms
# ignores alternate locations
# ignores hydrogens
import sys,os
import gzip

LETTER_CODE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
      'A': 'a',   'G': 'g',   'T': 't',   'C': 'c',   'U': 'u',
}
def letterCode(resName):
    return LETTER_CODE.get(resName,'?')

def openFile(filename):
    basename,ext = os.path.splitext(filename)
    if ext in ('.ent','.pdb','.FULL'):
        infh = file(filename)
    elif ext in ('.gz',):
        infh = gzip.open(filename)
    elif ext in ('.Z',):
        infh = os.popen('gunzip -c %s' % filename)
    else:
        raise IOError('Unknown file type')
    return infh
    

class PDBResidueId:
    def __init__(self,field=None):
        if isinstance(field,str):
            self.parseLine(field)
        elif isinstance(field,PDBAtom):
            self.resName = field.resName
            self.chainId = field.chainId
            self.resSeq = field.resSeq
            self.iCode = field.iCode

    def parseLine(self,line):
        self.resName = line[17:20].strip()
        self.chainId = line[21:22].strip()
        self.resSeq = int(line[22:26])
        self.iCode = line[26:27].strip()

    def __iter__(self):
        return iter((self.resName,self.chainId,self.resSeq,self.iCode))

    def __getitem__(self,index):
        return tuple(self)[index]

    def __cmp__(self,other):
        if other == None:
            return 1
        return cmp(tuple(self),tuple(other))

    def __str__(self):
        return "%(resName)s%(resSeq)s%(iCode)s:%(chainId)s" % vars(self)

class PDBResidue:
    def __init__(self, residueId=None):
        self.residueId = residueId
        self.atoms = []
        self._atomsByName = {}

    def __iter__(self):
        return iter(self.atoms)

    def append(self,atom):
        self.atoms.append(atom)

    def __getattr__(self, name):
        return getattr(self.residueId,name)

    def __getitem__(self, index):
        if type(index) == type(""):
            if not self._atomsByName:
                self.buildAtomDict()
            return self._atomsByName[index]
        return self.atoms[index]

    def get(self,index,default=None):
        try:
            return self[index]
        except KeyError:
            return default
    getAtom = get

    def __setitem__(self, index, value):
        self.atoms[index] = value

    def __repr__(self):
        return "%(resName)s" % vars(self.residueId)

    def __len__(self):
        return len(self.atoms)

    def buildAtomDict(self):
        self._atomsByName = dict([(a.name,a) for a in self.atoms])

    def __cmp__(self,other):
        return cmp(self.residueId,other.residueId)

    def getAtoms(self):
        return self.atoms

class PDBChain:
    def __init__(self, chainId=None):
        self.chainId = chainId
        self.residues = []
        self._sequence = ""

    def __iter__(self):
        return iter(self.residues)

    def append(self,residue):
        return self.residues.append(residue)

    def __getitem__(self, index):
        return self.residues[index]

    def __setitem__(self, index, value):
        self.residues[index] = value

    def __repr__(self):
        if not self.chainId: return '_'
        return "%(chainId)s" % vars(self)

    def __len__(self):
        return len(self.residues)

    def getSequence(self):
        if not self._sequence:
            self.buildSequence()
        return self._sequence
    sequence = getSequence

    def getResidues(self):
        return self.residues

    def buildSequence(self):
        resNames = [residue.resName for residue in self.residues]
        resLetters = map(letterCode,resNames)
        self._sequence = "".join(resLetters)


class PDBAtom:
    def __init__(self,line=None):
        if line:
            self.parse(line)

    def copy(self,other):
        self.__dict__.update(other.__dict__)
        return self

    def parse(self,line):
        self.line = line
        self.rectype = line[:6].strip()
        self.parseAtomId(line)
        self.parseResId(line)
        self.parseCoord(line)

    def parseAtomId(self,line):
        self.serial = int(line[6:11])
        self.name = line[12:16].strip()
        self.altLoc = line[16:17].strip()

    def parseResId(self,line):
        self.resName = line[17:20].strip()
        self.chainId = line[21:22].strip()
        self.resSeq = int(line[22:26])
        self.iCode = line[26:27].strip()

    def parseCoord(self,line):
        self.x = float(line[30:38])
        self.y = float(line[38:46])
        self.z = float(line[46:54])

    def parseMisc(self,line):
        self.occupancy = float(line[54:60])
        self.tempFactor = float(line[60:66])
        self.segId = line[72:76].strip()
        self.element = line[76:78].strip()
        self.charge = line[78:80].strip()

    def getCoord(self):
        return (self.x,self.y,self.z)
    coord = getCoord

    def getResidueId(self):
        return PDBResidueId(self)
#        return (self.resName,self.chainId,self.resSeq,self.iCode)
    resId = getResidueId

    def resIdString(self):
        return "%(resName)s%(resSeq)s%(iCode)s:%(chainId)s" % vars(self)

    def getAtomId(self):
        residstr = self.resIdString()
        atomname = self.name
        altloc = ""
        if self.altLoc:
            altloc = "'%s" % self.altLoc
        return "%(residstr)s@%(atomname)s%(altloc)s" % locals()
    fullIdString = getAtomId

    def __iter__(self):
        return iter((self.resName,self.chainId,self.resSeq,self.iCode))

    def __repr__(self):
        return "%(name)s" % vars(self)

    def __getattr__(self,name):
        if name == 'residueId':
            return self.resId()
        raise AttributeError(name)

    def appxEqu(self,other):
        if self.name == other.name and \
           self.resSeq == other.resSeq and \
           self.iCode == other.iCode and \
           self.resName == other.resName and \
           self.chainId == other.chainId:
            return 1
        return 0


class PDBFile:
    def __init__(self,filename=None):
        self._sequence = ""
        self._atomCoords = []
        self._residues = []
        self._chains = []
        self._chainsByName = {}
        self.atoms = []
        if filename:
            self.load(filename)

    def load(self,filename=None,fh=None):
        self.atoms = []

        if not fh:
            self.filename = filename
            infh = openFile(filename)
        else:
            try:
                self.filename = fh.name
            except AttributeError:
                self.filename = fh.filename
            infh = fh

        lastatom = None
        started = 0
        for line in infh.readlines():
            rectype = line[:6].strip()
            # get only first model
            if started and rectype in ('ENDMDL','END','MODEL'):
                break
            # ignore everything but ATOM
            if rectype not in ('ATOM',):
                continue
            started = 1
            atom = PDBAtom(line)
            # if atom same as previous, ignore
            if lastatom and atom.appxEqu(lastatom):
                continue
            # ignore hydrogen/deuterium atoms
            if atom.name[0] in ('H','D'):
                continue
            self.atoms.append(atom)
            lastatom = atom
        if not fh:
            infh.close()

    def atomCoords(self):
        if not self._atomCoords:
            self.buildAtomCoords()
        return self._atomCoords

    def buildAtomCoords(self):
        self._atomCoords = [atom.coord() for atom in self.atoms]

    def getResidues(self):
        if not self._residues:
            self.buildResidues()
        return self._residues
    residues = getResidues

    def buildResidues(self):
        residues = []
        lastResId = None
        for atom in self.atoms:
            resId = atom.resId()
            if resId != lastResId:
                lastResId = resId
                residues.append(PDBResidue(resId))
            residues[-1].append(atom)
        self._residues = residues

    def getChains(self):
        if not self._chains:
            self.buildChains()
        return self._chains
    chains = getChains

    def getChain(self,chainId):
        if not self._chainsByName:
            self.buildChains()
        return self._chainsByName.get(chainId,[None])[0]

    def getChainIds(self):
        if not self._chainsByName:
            self.buildChains()
        ids = self._chainsByName.keys()
        ids.sort()
        return ids

    def buildChains(self):
        chains = []
        chainsByName = {}
        lastChainId = None
        for residue in self.residues():
            chainId = residue.chainId
            if lastChainId != chainId:
                lastChainId = chainId
                chain = PDBChain(chainId)
                chains.append(chain)
                chainsByName.setdefault(chainId,[]).append(chain)
            chains[-1].append(residue)
        self._chains = chains
        self._chainsByName = chainsByName
            
    def sequence(self):
        if not self._sequence:
            self.buildSequence()
        return self._sequence

    def buildSequence(self):
        resAtoms = self.residues()
        resNames = [atomList[0].resName for atomList in resAtoms]
        resLetters = map(letterCode,resNames)
        self._sequence = "".join(resLetters)
