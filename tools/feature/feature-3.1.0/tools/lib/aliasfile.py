#! /usr/bin/env python

# Represents a mapping of alias to pdbid
# File format:
# <canonical> <aliases>...
import re

class AliasFile:
    def __init__(self,filename=None):
        self.map = {}
        if filename:
            self.load(filename)

    def load(self,filename):
        self.filename = filename
        infh = file(self.filename)
        for line in infh.readlines():
            if not line.strip():
                continue
            if line[:1] == "#":
                continue
            fields = line.split()
            for field in fields:
                self.map.setdefault(field,[]).append(fields[0])

    def get(self,index,default=None):
        return self.map.get(index,[default])[0]

    def getall(self,index):
        return self.map.get(index,[])

    def getmatch(self,index):
        return [value[0] for key,value in self.map.items() if re.match(index,key)]

    def getallmatch(self,index):
        return [value for key,value in self.map.items() if re.match(index,key)]

    def __getitem__(self,index):
        return self.map[index]
