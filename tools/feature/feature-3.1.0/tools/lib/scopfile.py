#! /usr/bin/env python

# Scop file formats

import sys,os
from utils import chomp,sort,uniq

DEFAULT_SCOP_DIR = '/project1/structure/mliang/scop'
DEFAULT_CLASS_FILENAME = os.path.join(DEFAULT_SCOP_DIR,'dir.cla.scop.txt_1.65')

class ScopClassFileEntry:
    def __init__(self,line=None):
        self.isinit = False
        if line:
            self.parse(line)

    def parse(self,line):
        line = chomp(line)
        self.line = line
        fields = line.split('\t')
        try:
            self.domainid = fields[0]
            self.pdbid = fields[1]
            self.chainidField = fields[2]
            self.scopfamily = fields[3]
            self.scopid = fields[4]
            self.scopinfoField = fields[5]

            self.chainids = []
            self.chainrange = {}
            for chainid in self.chainidField.split(','):
                chainfields = chainid.split(':')
                if len(chainfields) == 1:
                    if chainfields[0] == '-':
                        chainfields = ['','']
                    else:
                        chainfields = ['',chainfields[0]]
                self.chainids.append(chainfields[0])
                self.chainrange[chainfields[0]] = chainfields[1]
        except:
            return
        self.isinit = True
        
    def __nonzero__(self):
        return self.isinit

    def __str__(self):
        return self.line

class ScopClassFile:
    def __init__(self,filename=None):
        if filename:
            self.read(filename)

    def read(self,filename):
        infh = file(filename)
        self.entries = []
        self.entriesByPdbid = {}
        self.entriesByDomainid = {}
        for line in infh:
            if line[0] == '#':
                continue
            entry = ScopClassFileEntry(line)
            if entry:
                self.entries.append(entry)
                self.entriesByPdbid.setdefault(entry.pdbid,[]).append(entry)
                self.entriesByDomainid.setdefault(entry.domainid,[]).append(entry)

    def getScopFamily(self,id,chainid=None,resnum=None):
        if len(id) == 4:
            families = []
            for entry in self.entriesByPdbid.get(id,[]):
                if chainid is not None:
                    # check chainid
                    pass
                if resnum is not None:
                    # check resnum
                    pass
                families.append(entry.scopfamily)
            families = sort(uniq(families))
            return ';'.join(families)
        return ''
