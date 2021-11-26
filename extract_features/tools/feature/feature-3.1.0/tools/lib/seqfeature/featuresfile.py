#! /usr/bin/env python
# $Id: featuresfile.py,v 1.2 2004/10/27 00:08:39 mliang Exp $

import sys,os
from utils import chomp

class FeaturesEntry:
    def __init__(self,line):
        self.label = None
        self.pdbid = None
        self.index = None
        self.values = []
        self.comments = None
        self.Parse(line)

    def Parse(self,line):
        fields = chomp(line).split('#',1)
        data = fields[0].strip()
        if len(fields) > 1:
            self.comments = fields[1].strip()
        fields = data.split('\t')
        self.label = fields[0]
        self.values = map(float,fields[1:])
        subfields = self.label.split('_')
        self.pdbid = '_'.join(subfields[1:-1])
        self.index = int(subfields[-1])

class FeaturesFile:
    DEFAULT_TYPE = 'label'

    def __init__(self,filename=None):
        self.entries = []
        self.indexMap = {'pdbid':{},'label':{}}
        if filename:
            self.Load(filename)

    def Load(self,filename):
        infh = file(filename)
        for line in infh:
            entry = FeaturesEntry(line)
            self.Add(entry)

    def Add(self,entry):
        self.entries.append(entry)
        self.indexMap['label'].setdefault(entry.label,[]).append(entry)
        self.indexMap['pdbid'].setdefault(entry.pdbid,[]).append(entry)

    def Get(self,key,type='label'):
        return self.indexMap[type][key][0]

    def __getitem__(self,key):
        return self.indexMap[self.DEFAULT_TYPE][key][0]
