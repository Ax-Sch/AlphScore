#! /usr/bin/env python

# Parses DSSP File

import sys,os
import re

class Record:
    def __init__(self,line=None):
        if line:
            self.parse(line)

    def parse(self,line):
        # handle chain breaks
        self.chainbreak = 0
        self.chainborder = 0
        if line[13:14] == '!':
            self.chainbreak = 1
            if line[14:15] == '*':
                self.chainborder = 1
            return

        self.resnum = int(line[0:5]) # dssp resnum

        self.resseq = int(line[5:10])
        self.icode = line[10:11].strip()
        self.chainid = line[11:12].strip()
        self.resname = line[13:14].strip() # single letter AA

        self.structure = line[16:25] # structure

        self.bp1 = int(line[25:29]) # bridge partner 1 resnum
        self.bp2 = int(line[29:33]) # bridge partner 2 resnum
        self.bplabel = line[33:34].strip() # sheet label

        self.acc = int(line[34:38]) # sq Angstroms

        def parsehbond(field):
            vals = field.split(',')
            return (int(vals[0]),float(vals[1]))

        self.nho1 = parsehbond(line[39:50]) # offset, energy
        self.ohn1 = parsehbond(line[50:61])
        self.nho2 = parsehbond(line[61:72])
        self.ohn2 = parsehbond(line[72:83])

        self.tco = float(line[83:91])
        self.kappa = float(line[91:97])
        self.alpha = float(line[97:103])

        self.phi = float(line[103:109])
        self.psi = float(line[109:115])

        self.xca = float(line[115:122])
        self.yca = float(line[122:129])
        self.zca = float(line[129:136])

    def getResId(self):
        return (self.resseq,self.icode,self.chainid)
    resId = getResId

    def getCoord(self):
        return (self.xca, self.yca, self.zca)
    coord = getCoord
    

class DSSPFile:
    def __init__(self,filename=None):
        self.headers = []
        self.records = []
        if filename:
            self.load(filename)

    def load(self,filename=None):
        self.filename = filename
        infh = file(filename)

        # ignore most of the header information
        self.headers = []
        while 1:
            line = infh.readline()
            if not line:
                break
            self.headers.append(line)
            if line.find('  #  RESIDUE AA') != -1:
                break

        self.records = []
        while 1:
            line = infh.readline()
            if not line:
                break
            record = Record(line)
            self.records.append(record)


if __name__ == '__main__':
    import sys,os
    def eusage():
        print "Usage: %s DSSPFILE" % os.path.basename(sys.argv[0])
        sys.exit(1)

    args = sys.argv[1:]
    if len(args) != 1:
        eusage()

    dsspfilename = args[0]
    dsspfile = DSSPFile(dsspfilename)
    print "Num records:",len(dsspfile.records)
    for record in dsspfile.records:
        if not record.chainbreak:
            print ".".join(map(str,record.getResId())),record.structure[0],record.acc
