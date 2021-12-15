#!/usr/pubsw/bin/python
# $Id: Motif.py,v 1.4 2002/04/05 04:08:39 mliang Exp $
# Represents Sequence Motifs
# Copyright (c) 2000 Mike Liang.  All rights reserved.

# Motif class
import re
import string

class Motif:
    def __init__(self, pattern):
        self.motif = pattern
        self.pattern = re.compile(self.motif.upper())

    # Returns MotifScanResults
    def scanProtein(self, protein, modelNo=0):
        results = MotifScanResults(protein, self.motif)
        model = protein.models[modelNo]
        for chain in model.chains:
            pos = 0
            match = self.pattern.search(chain.sequence, pos)
            while match != None:
                pos = match.start() + 1
                results.addMatch(modelNo, chain, match.start(), match.end())
                match = self.pattern.search(chain.sequence, pos)
        return results

class MotifMatch:
    def __init__(self, protein, modelNo, chainID, start, end,
                 sequence, residues, pattern):
        self.protein = protein
        self.modelNo = modelNo
        self.chainID = chainID
        self.start = start
        self.end = end
        self.sequence = sequence
        self.residues = residues
        self.pattern = pattern

    def displayInfo(self):
        print "Protein=%s modelNo=%d chainID=%s span=[%d,%d] seq=%s" % \
              (self.protein.filename, self.modelNo, self.chainID,
               self.start, self.end, self.sequence)
        
    def __str__(self):
        return "%s %s %d %d %s %s" % (self.protein.pdbID, self.chainID,
                                      self.start, self.end,
                                      self.pattern,
                                      self.sequence)
        
        
class MotifScanResults:
    def __init__(self, protein, pattern):
        self.pattern = pattern
        self.protein = protein
        self.matches = []
        
    def addMatch(self, modelNo, chain, start, end):
        match = MotifMatch(self.protein, modelNo, chain.chainID,
                           start, end,
                           chain.sequence[start:end],
                           chain.residues[start:end], self.pattern)
        self.matches.append(match)

    def displayInfo(self):
        map(MotifMatch.displayInfo, self.matches)


class MotifHit:
    def __init__(self, line=None):
        self.pdbid = ""
        self.chainid = ""
        self.start = 0
        self.end = 0
        self.pattern = ""
        self.sequence = ""
        
        if line:
            self.parseLine(line)

    def __str__(self):
        return "%s %s %d %d %s %s" % \
               (self.pdbid, self.chainid, self.start, self.end, self.pattern,
                self.sequence)

    def parseLine(self, line):
        words = line.split(" ")
        self.pdbid = words[0]
        self.chainid = words[1]
        self.start = int(words[2])
        self.end = int(words[3])
        self.pattern = words[4]
        self.sequence = words[5]
