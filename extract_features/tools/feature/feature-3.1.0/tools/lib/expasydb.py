#! /usr/bin/env python
# $Id: expasydb.py,v 1.2 2004/05/22 21:50:15 mliang Exp $
# Copyright (c) 2003 Mike Liang. All rights reserved.

# expasy database format

def getline(fh):
    line = fh.readline()
    if not line: raise IOError('no more record')
    return line

class ExpasyRecord:
    def __init__(self, fh=None):
        if fh:
            self.loadFromFile(fh)

    def loadFromFile(self, fh):
        self.lines = []
        self.linesByType = {}

        # Scan for ID line, get offset
        offset = fh.tell()
        self.ignored = 0
        line = getline(fh)
        while line[:2] != 'ID':
            self.ignored += 1
            offset = fh.tell()
            line = getline(fh)
        self.offset = offset
        self.addLine(line)

        # Load rest of the lines
        while 1:
            line = getline(fh)
            self.addLine(line)
            if line[:2] == '//':
                break

    def addLine(self, line):
        type = line[:2]
        self.lines.append(line)
        self.linesByType.setdefault(type,[]).append(line)

class SwissProtRecord(ExpasyRecord):
    class CrossReferenceEntry:
        def __init__(self,field=None):
            if field:
                self.parse(field)

        def parse(self,field):
            fields = [f.strip() for f in field.split(';')]
            self.db_id = fields[0]
            self.ids = fields[1:]

    def __init__(self, fh=None):
        ExpasyRecord.__init__(self, fh)
        self.parseID()
        self.parseAC()
        self.parseDR()

    def parseID(self):
        for line in self.linesByType.get('ID',[]):
            fields = line[2:].split(';')
            self.entryName = fields[0].strip().split()[0]
            break

    def parseAC(self):
        self.accessionNumbers = []
        for line in self.linesByType.get('AC',[]):
            fields = line[2:].split(';')
            ids = [id.strip() for id in fields]
            ids = filter(lambda a:a, ids)
            self.accessionNumbers.extend(ids)

    def parseDR(self):
        self.databaseReferences = []
        self.databaseReferencesByType = {}
        for line in self.linesByType.get('DR',[]):
            entry = self.CrossReferenceEntry(line[2:])
            self.databaseReferences.append(entry)
            self.databaseReferencesByType.setdefault(entry.db_id,[]).append(entry)

class PrositeRecord(ExpasyRecord):
    class SwissProtCrossReference:
        def __init__(self,field=None):
            if field:
                self.parse(field)

        def parse(self,field):
            self.ac_nb,self.entry_name,self.type = [f.strip() for f in field.split(',')]

    def __init__(self, fh=None):
        ExpasyRecord.__init__(self, fh)
        self.parseID()
        self.parseAC()
        self.parseDR()

    def parseID(self):
        for line in self.linesByType.get('ID',[]):
            fields = line[2:].split(';')
            self.entryName = fields[0].strip()

    def parseAC(self):
        self.accessionNumbers = []
        for line in self.linesByType.get('AC',[]):
            fields = line[2:].split(';')
            ids = [id.strip() for id in fields]
            ids = filter(lambda a:a, ids)
            self.accessionNumbers.extend(ids)

    def parseDR(self):
        self.swissProtEntries = []
        self.swissProtEntriesByType = {}
        for line in self.linesByType.get('DR',[]):
            fields = line[2:].split(';')
            for field in fields:
                if not field.strip():
                    continue
                entry = self.SwissProtCrossReference(field)
                self.swissProtEntries.append(entry)
                self.swissProtEntriesByType.setdefault(entry.type,[]).append(entry)

class IndexFile:
    def __init__(self,filename=None):
        self.indices = {}
        if filename:
            self.load(filename)

    def load(self,filename):
        self.filename = filename
        fh = file(filename)
        for line in fh.readlines():
            key,offset = line[:-1].split('\t')
            self.indices[key] = int(offset)
        fh.close()

    def __getitem__(self,key):
        return self.indices[key]

class DatabaseFile:
    def __init__(self,dataFilename=None,indexFilename=None):
        if dataFilename:
            self.dbOpen(dataFilename)
        if indexFilename:
            self.idxOpen(indexFilename)

    def dbOpen(self,dataFilename):
        self.dataFilename = dataFilename
        self.datafh = file(dataFilename)

    def idxOpen(self,idxFilename):
        self.idxFilename = idxFilename
        self.index = IndexFile(idxFilename)

class PrositeDatabaseFile(DatabaseFile):
    def __init__(self,dataFilename=None,indexFilename=None):
        DatabaseFile.__init__(self,dataFilename,indexFilename)

    def getRecordByOffset(self,offset):
        self.datafh.seek(offset)
        return PrositeRecord(self.datafh)

    def getNextRecord(self):
        return PrositeRecord(self.datafh)

    def getRecordByIndex(self,index):
        self.datafh.seek(self.index[index])
        return PrositeRecord(self.datafh)
