
# $Id: threemotifrecord.py,v 1.2 2002/04/05 04:10:53 mliang Exp $
# Copyright (c) 2001 Mike Liang.  All rights reserved.

# 3Motif Record

class ThreeMotifRecord:
    class Match:
        def __init__(self, line=None):
            self.pdbid = ""
            self.chainid = ""
            self.offset = None
            self.sequence = ""
            self.startRes = None
            self.stopRes = None
            if line:
                self.loadFromString(line)

        def loadFromString(self, line):
            (pdbchain,dummy,remains) = line.split()
            (self.pdbid, self.chainid) = pdbchain.split("_")

            fields = remains.split("|")
            (dummy, self.offset) = fields.pop(0).split(":")
            self.offset = int(self.offset)
            self.sequence = fields.pop(0)
            self.startRes = int(fields.pop(0))
            self.stopRes = int(fields.pop(0))


    def __init__(self, file=None):
        self.expectation = None
        self.specificity = None
        self.percent = None
        self.blockac = ""
        self.description = ""
        self.pattern = ""
        self.matches = []
        if file:
            self.loadFromFile(file)

    def loadFromFile(self, file):
        line = file.readline()
        if not line or line[0] != ">":
            return None

        fields = line[1:-1].split("|")
        self.expectation = float(fields.pop(0))
        self.specificity = float(fields.pop(0))
        self.percent = float(fields.pop(0)[:-1])/100
        self.blockac = fields.pop(0)
        self.description = fields.pop(0)
        self.pattern = fields.pop(0)

        line = file.readline()
        while line:
            if line[0] == ">":
                file.seek(-len(line),1)
                break
            match = ThreeMotifRecord.Match()
            match.loadFromString(line[:-1])
            self.matches.append(match)
            line = file.readline()
    
