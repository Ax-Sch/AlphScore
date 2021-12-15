#! /usr/env/bin python

# Generates index file for file
import cPickle as pickle

class IndexEntry:
    def __init__(self):
        self.name = None
        self.offset = None

class IndexFile:
    def __init__(self):
        self.entries = []
        self.entriesByName = {}

    def add(self,name,offset):
        entry = IndexEntry(name,offset)
        self.entries.append(entry)
        self.entriesByName[entry.name] = entry

    def __getitem__(self,key):
        return self.entriesByName[key]

def write(indexFile,outfh):
    pickle.dump(indexFile,outfh,-1)

def read(infh):
    return pickle.load(infh)
