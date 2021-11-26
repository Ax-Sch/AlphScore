#! /usr/bin/env python

import re
from pdbutils import goodpdbfilename
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Structure import Structure


LETTER_CODE = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
      'A': 'a',   'G': 'g',   'T': 't',   'C': 'c',   'U': 'u',
}
def get_letter_code(resName):
    return LETTER_CODE.get(resName,'?')

def get_models(structure,firstmodel=1):
    models = []
    for model in structure.get_iterator():
        models.append(model)
    if firstmodel:
        return models[:1]
    return models

def get_chains(structure):
    chains = []
    for model in get_models(structure):
        for chain in model.get_iterator():
            chains.append(chain)
    return chains

def get_residues(structure):
    residues = []
    for chain in get_chains(structure):
        for residue in chain.get_iterator():
            residues.append(residue)
    return residues

def get_atoms(structure):
    atoms = []
    for residue in get_residues(structure):
        for atom in residue.get_iterator():
            atoms.append(atom)
    return atoms

def is_het(residue):
    if residue.level == 'A':
        residue = residue.parent
    return residue.id[0] != ' '

def is_wat(residue):
    if residue.level == 'A':
        residue = residue.parent
    return residue.id[0] == 'W'

def get_element(atom):
    element = re.sub(r'\W',' ',atom.name).strip()[:1]
    return element

def is_hydrogen(atom):
    return get_element(atom) in ('H','D','Q')

def is_backbone(atom):
    return not is_het(atom) and atom.name in ('N','CA','C','O','OXT')

def get_sequence(chain,usehetatm=0):
    letters = [get_letter_code(residue.resname) for residue in chain.get_iterator() if usehetatm or not is_het(residue)]
    return ''.join(letters)

def get_structure(pdbid):
    pdbfilename = goodpdbfilename(pdbid)
    if not pdbfilename:
        return
    return PDBParser().get_structure(pdbid,pdbfilename)
    
