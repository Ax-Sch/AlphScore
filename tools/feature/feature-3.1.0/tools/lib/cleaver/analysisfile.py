#! /usr/bin/env python

# Represents Cleaver's analysisdata.txt file

import sys,os
import re

class AnalysisFile:
    """ Represents Cleaver's analysisdata.txt raw output file """
    
    class Score:
        """ Represents a Score entry in analysisdata.txt """
        def __init__(self, name, score):
            self.name = name
            self.score = float(score)

        def fieldnames(self):
            return ("name","score")

        def __iter__(self):
            return iter([self.__dict__[field] for field in self.fieldnames()])

        def __repr__(self):
            return str(tuple(self))

    class Feature:
        """ Represents a Feature entry in analysisdata.txt """
        def __init__(self, feature, power, weight):
            self.feature = feature
            self.power = float(power)
            self.weight = float(weight)

        def fieldnames(self):
            return ("feature","power","weight")

        def __iter__(self):
            return iter([self.__dict__[field] for field in self.fieldnames()])

        def __repr__(self):
            return str(tuple(self))

    def __init__(self, filename=None):
        if filename:
            self.load(filename)

    def load(self, filename):
        self.filename = filename
        self.accuracy = None
        self.scores = []
        self.features = []
        
        infh = file(filename)
        # Get accuracy
        while 1:
            line = infh.readline()
            if not line: break
            if line.find('accuracy of prediction') != -1:
                match = re.search(r'(\d*\.?\d*) %',line)
                if match:
                    self.accuracy = float(match.group(1))/100
                    break
        
        # Get scores
        while 1:
            line = infh.readline()
            if not line: break
            if line.find('name\t\tscore') != -1:
                break
        while 1:
            line = infh.readline()
            if not line: break
            try:
                name,score = line.strip().split()
                self.scores.append(AnalysisFile.Score(name,score))
            except ValueError:
                break
            
        # Get features
        while 1:
            line = infh.readline()
            if not line: break
            if line.find('feature\t\tpredictive power\t\tweight') != -1:
                break
        while 1:
            line = infh.readline()
            if not line: break
            try:
                feature,power,weight = line.strip().split()
                self.features.append(AnalysisFile.Feature(feature,power,weight))
            except ValueError:
                break

        infh.close()

    def exportFeatures(self,outfh=sys.stdout,sep="\t"):
        for feature in self.features:
            print >>outfh, sep.join(map(str,feature))

    def exportScores(self,outfh=sys.stdout,sep="\t"):
        for score in self.scores:
            print >>outfh, sep.join(map(str,score))

    def export(self,outfh=sys.stdout,sep="\t"):
        print >>outfh, "# Accuracy"
        print >>outfh, self.accuracy
        print >>outfh, "# Scores"
        print >>outfh, "#",sep.join(self.scores.fieldnames())
        self.exportScores(outfh,sep)
        print >>outfh, "# Features"
        print >>outfh, "#",sep.join(self.features.fieldnames())
        self.exportFeatures(outfh,sep)
