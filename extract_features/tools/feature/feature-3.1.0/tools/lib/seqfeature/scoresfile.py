#! /usr/bin/env python
# $Id: scoresfile.py,v 1.3 2004/10/27 00:08:39 mliang Exp $
# Copyright (c) 2004 by Mike Liang. All rights reserved.

# Represents seqfeature .scores output file

import sys,os
from utils import chomp

class ScoresEntry:
    def __init__(self,line=None):
        if line:
            self.Parse(line)

    def Parse(self,line):
        self.line = line
        fields = chomp(line).split('\t')
        self.label = fields[0]
        subfields = fields[0].split('_')
        self.pdbid = '_'.join(subfields[1:-1])
        self.index = int(subfields[-1])
        self.score = float(fields[1])
        self.x = float(fields[2])
        self.y = float(fields[3])
        self.z = float(fields[4])
        self.description = '\t'.join(fields[5:])

    def __cmp__(self,rhs):
        return cmp(self.score,rhs.score)

    def __str__(self):
        return self.line

class ScoresFile:
    DEFAULT_TYPE = 'label'
    def __init__(self,filename=None):
        self.entries = []
        self.indexMap = {'pdbid':{},'label':{}}
        if filename:
            self.Load(filename)

    def Load(self,filename):
        infh = file(filename)
        for line in infh:
            entry = ScoresEntry(line)
            self.Add(entry)

    def Add(self,entry):
        self.entries.append(entry)
        self.indexMap['label'].setdefault(entry.label,[]).append(entry)
        self.indexMap['pdbid'].setdefault(entry.pdbid,[]).append(entry)

    def SortAscending(self):
        self.entries.sort()

    def SortDescending(self):
        self.entries.sort()
        self.entries.reverse()

    def __len__(self):
        return len(self.entries)

    def __iter__(self):
        return iter(self.entries)

    def Get(self,key,type='label'):
        return self.indexMap[type][key][0]

    def __getitem__(self,key):
        return self.indexMap[self.DEFAULT_TYPE][key][0]
