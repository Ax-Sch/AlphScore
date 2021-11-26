#! /usr/bin/env python

# Represents the Cleaver input file formats
# this includes files that contain all sites all nonsites or mixed

class InputFile:
    def Entry:
        def __init__(self,label,values):
            self.label = label
            self.values = values

    def __init__(self, filename=None):
        if filename:
            self.load(filename)

    def load(self,filename):
        self.filename = filename
        self.fieldnames = []
        self.entries = []

        infh = file(filename)
        self.fieldnames = infh.readline()[:-1].split('\t')
        numfields = len(self.fieldnames)
        for line in infile.readlines():
            entry = line[:-1].split('\t')[:numfields]
            label = entry[0]
            values = entry[1:]
            self.entries.append(InputFile.Entry(label,values))
        infh.close()

    def import(self,values,labels,fieldnames):
        """ Imports tab delimited file """
        assert(values and labels and fieldnames)
        assert(len(values[0]) == len(fieldnames))
        assert(len(values) == len(labels))
        self.fieldnames = ["LABEL"] + fieldnames
        for idx in range(len(labels)):
            self.entries.append(InputFile.Entry(labels[idx],values[idx]))
        
