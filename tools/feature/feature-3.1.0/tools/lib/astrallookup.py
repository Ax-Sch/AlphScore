#! /usr/bin/env python
# $Id: astrallookup.py,v 1.1 2004/11/17 09:29:08 mliang Exp $
# Copyright (c) 2004 Mike Liang. All rights reserved.
# Looks up astral domain id for pdbid
import sys,os
import fileinput

DEFAULT_ASTRAL_ID_FILENAME = os.path.join(os.environ.get('ASTRAL_DIR',''),'scopseq-1.65','astral-scopdom-seqres-all-1.65.id')

class AstralIdDatabase:
    def __init__(self,filename):
        self.entries = {}
        self.Load(filename)

    def Load(self,filename):
        self.filename = filename
        infh = file(filename)
        for line in infh:
            id = line.strip()
            self.Add(id)

    def Add(self,id):
        pdbid = id[1:5].lower()
        self.entries.setdefault(pdbid,[]).append(id)

    def get(self,key,default=None):
        key = key.lower()
        if key in self.entries:
            return self.entries[key]
        return default

def DefaultAstralIdDatabase():
    aidFilename = DEFAULT_ASTRAL_ID_FILENAME 
    aid = AstralIdDatabase(aidFilename)
    return aid

def main():
    aid = DefaultAstralIdDatabase()
    args = sys.argv[1:]
    infh = fileinput.input(args)
    for line in infh:
        pdbid = line.strip()
        for astralid in aid.get(pdbid,[]):
            print astralid

if __name__ == '__main__':
    main()
