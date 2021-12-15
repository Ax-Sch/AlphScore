#! /usr/bin/env python
# $Id: modelfile.py,v 1.2 2004/10/14 00:34:26 mliang Exp $

# SeqFeature .model parser

import sys,os
from utils import chomp

class Parameter:
    SPECIAL_PARAMETERS = {
        'NUM_BINS': ('numBins',int),
        'P_SITE': ('pSite',float),
        'NUM_PROPERTIES': ('numProperties',int),
        'NUM_SHELLS': ('numShells',int)
    }

    def __init__(self,line):
        self.Parse(line)

    def Parse(self,line):
        fields = chomp(line).split('\t')
        name = fields[1]
        value = fields[2]
        if name not in self.SPECIAL_PARAMETERS:
            return None

        name,func = self.SPECIAL_PARAMETERS[name]
        self.name = name
        self.value = func(value)
        return self

class Property:
    NUM_BINS = 5
    EXTENDED_FIELDNAMES = [
        ('tvalue',float)
    ]

    def __init__(self,line):
        self.Parse(line)

    def Parse(self,line):
        self.line = line
        fields = chomp(line).split('\t')
        self.name = fields[0]
        subfields = self.name.split('-')
        self.property = '-'.join(subfields[:-1])
        self.volume = int(subfields[-1])
        self.pvalue = float(fields[1])
        self.minval = float(fields[2])
        self.binsize = float(fields[3])
        self.scores = map(float,fields[4:4+self.NUM_BINS])

        extendedFields = fields[4+self.NUM_BINS:]
        if extendedFields and extendedFields[0] == '#':
            for idx,value in enumerate(extendedFields[1:]):
                if not idx < len(self.EXTENDED_FIELDNAMES):
                    warning('Too many extended fields')
                    break
                name,func = self.EXTENDED_FIELDNAMES[idx]
                setattr(self,name,func(value))

    def __str__(self):
        return self.line

class ModelFile:
    def __init__(self,filename=None):
        self.propertyList = []
        self.indexMap = {'property':{}, 'volume':{}}

        self.properties = []
        self.volumes = []

        if filename:
            self.Load(filename)

    def Load(self,filename):
        self.filename = filename
        infh = file(filename)
        for line in infh:
            if line[0] == '#':
                param = Parameter(line)
                setattr(self,param.name,param.value)
                if param.name == 'numBins':
                    Property.NUM_BINS = param.value
                continue
            property = Property(line)
            self.Add(property)

    def Add(self,prop):
        self.propertyList.append(prop)

        # keep track of unique properties and unique volumes
        if prop.property not in self.indexMap['property']:
            self.properties.append(prop.property)
        if prop.volume not in self.indexMap['volume']:
            self.volumes.append(prop.volume)

        # add property into appropriate index
        self.indexMap['property'].setdefault(prop.property,[]).append(prop)
        self.indexMap['volume'].setdefault(prop.volume,[]).append(prop)

    def Get(self,property=None,volume=None):
        if property is None and volume is None:
            return self.propertyList
        if volume is None:
            return self.indexMap['property'].get(property,[])
        if property is None:
            return self.indexMap['volume'].get(volume,[])
        return [ p for p in self.indexMap['property'].get(property,[]) if p.volume == volume ]

    def GetPropertyNameList(self):
        return self.properties

    def GetShellList(self):
        return self.volumes

    def GetNumProperties(self):
        return len(self.properties)

    def GetNumShells(self):
        return len(self.volumes)
        

if __name__ == '__main__':
    outfh = sys.stdout
    args = sys.argv[1:]
    m = ModelFile(args[0])
    outfh.write(''.join(map(str,m.Get(property='ATOM-TYPE-IS-Ca'))))
