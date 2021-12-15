#! /usr/bin/env python

import sys,os
import time
from logutils import *
from Bio import Fasta

DEFAULT_DB_DIR = '/project1/structure/mliang/pdb/derived_data'
DEFAULT_INDEX_FILENAME = 'pdb_seqres.idx'
DEFAULT_SEQRES_FILENAME = 'pdb_seqres.txt'

class PdbSeqresLookup:
    def __init__(self,**kwargs):
        self.db_dir = DEFAULT_DB_DIR
        self.index_filename = DEFAULT_INDEX_FILENAME
        self.seqres_filename = DEFAULT_SEQRES_FILENAME
        for key,value in kwargs:
            if key in ['db_dir']:
                self.db_path = arg
            elif key in ['index_filename']:
                self.index_filename = arg
            elif key in ['seqres_filename']:
                self.seqres_filename = arg
        self.full_index_filename = os.path.join(self.db_dir,self.index_filename)
        self.full_seqres_filename = os.path.join(self.db_dir,self.seqres_filename)
        self.offsets = {}
        self.namesByPdbid = {}
        self.seqres_file = file(self.full_seqres_filename)
        self.fasta_parser = Fasta.RecordParser()
        self.load_index_file()

    def init_index_file(self):
        raise IndexError("Not implemented yet")

    def save_index_file(self):
        raise IndexError("Not implemented yet")

    def load_index_file(self):
        index_file = file(self.full_index_filename)
        for line in index_file:
            name,offset = line.split(None,1)
            name = name.lower()
            pdbid = name[:4]
            self.offsets[name] = int(offset)
            self.namesByPdbid.setdefault(pdbid,[]).append(name)
        index_file.close()

    def query_gen(self,ids):
        if type(ids) == str:
            ids = [ids]
        idlist = []
        for name in ids:
            if len(name) == 4:
                idlist.extend(self.namesByPdbid.get(name,[]))
            else:
                idlist.append(name)
        for name in idlist:
            offset = self.offsets.get(name)
            if offset is None:
                continue
            self.seqres_file.seek(offset)
            record = self.fasta_parser.parse(self.seqres_file)
            yield record

def main():
    args = sys.argv[1:]
    if not args:
        idlist = [l.strip() for l in sys.stdin]
    else:
        idlist = args
    idlist = [x.lower() for x in idlist]

    lookup = PdbSeqresLookup()
    for id in idlist:
        for e in lookup.query_gen(id):
            print e

if __name__ == '__main__':
    main()
