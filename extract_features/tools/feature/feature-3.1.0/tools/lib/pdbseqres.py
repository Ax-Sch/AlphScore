#! /usr/bin/env python
# $Id: pdbseqres.py,v 1.2 2004/05/22 21:50:15 mliang Exp $
# Copyright (c) 2003 Mike Liang. All rights reserved.

# Find the specified pattern in the PDB
# Uses the pdb_seqres.txt file

import sys,os
import re


class SeqResRecord:
    def __init__(self,line=None):
        self.header = None
        self.lines = []
        if line:
            self.parseHeader(line)

    def parseHeader(self,line):
        self.header = line.strip()
        fields = self.header.split(None,3)
        self.pdbId = fields[0][1:5]
        self.chainId = fields[0][6:7]
        self.type = fields[1].split(':')[1]
        self.length = fields[2].split(':')[1]
        self.description = fields[3]

    def append(self,line):
        self.lines.append(line.strip())

    def sequence(self):
        return "".join(self.lines)

class SeqResFile:
    def __init__(self,filename=None):
        if filename:
            self.parse(filename)

    def parse(self,filename):
        self.filename = filename
        self.records = []
        fh = file(filename)
        for line in fh.readlines():
            if line[0] == '>':
                self.records.append(SeqResRecord(line))
            else:
                self.records[-1].append(line)

def eusage():
    print "Usage: %s PATTERN" % os.path.basename(sys.argv[0])
    sys.exit(1)

def main():
    import pdbutils
    args = sys.argv[1:]
    try:
        pattern, = args
    except ValueError:
        eusage()

    regex = re.compile(pattern.upper())
    srf = SeqResFile(pdbutils.pdbSeqresFilename())
    for record in srf.records:
        m = regex.search(record.sequence())
        if m:
            print m.span(), record.pdbId, record.chainId, record.description

if __name__ == "__main__":
    main()
