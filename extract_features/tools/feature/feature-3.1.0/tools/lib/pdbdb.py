#! /usr/bin/env python
# $Id: pdbdb.py,v 1.2 2004/09/06 22:18:13 mliang Exp $

import sys,os
import re
from Bio import Fasta
from indexdb import IndexedDatabase
from simplepdb import PDBFile
from pdbutils import goodpdbfilename
import utils

PDB_DIR = os.environ.get('PDB_DIR','/Users/mliang/db/pdb')
STRUCTURE_DATABASE_DIR = PDB_DIR
SEQRES_DATABASE_FILE = utils.pathsearch(PDB_DIR,'derived_data/pdb_seqres.txt')

PDBID_RE = re.compile(r'\d\w{3}')

class Callable:
    def __init__(self,callable):
        self.__call__ = callable


def SeqresEntry_getattr(self,name):
    if name == 'name':
        self.name = self.title.split()[0]
        return self.name
    elif name == 'pdbid':
        m = PDBID_RE.search(self.title)
        if m:
            self.pdbid = m.group()
            return self.pdbid
    raise AttributeError(name)


class SeqresDatabase(IndexedDatabase):
    KEY_TYPES = ['name','pdbid']
    DEFAULT_KEY_TYPE = 'pdbid'
    PARSER_CLASS = Fasta.RecordParser
    ITERATOR_CLASS = Fasta.Iterator
    GET_ATTR = Callable(SeqresEntry_getattr)

    def __init__(self,filename=None):
        IndexedDatabase.__init__(self,filename)

class StructureDatabase:
    def __init__(self,dirname=None):
        self.dirname = dirname

    def getEntry(self,key):
        filename = goodpdbfilename(key,self.dirname)
        if not filename:
            return
        entry = PDBFile(filename)
        return entry
    
seqresDB = SeqresDatabase(SEQRES_DATABASE_FILE)
structureDB = StructureDatabase(STRUCTURE_DATABASE_DIR)
