import os
import cPickle as pickle

class IndexedDatabase:
    KEY_TYPES = ['name']
    DEFAULT_KEY_TYPE = KEY_TYPES[0]
    PARSER_CLASS = None
    ITERATOR_CLASS = None
    FILENAME_KEY = '__FILENAME__'
    GET_ATTR=getattr
    
    def __init__(self,filename):
        self.entries = []
        self.entryMap = {}
        self.indexMap = {}
        for keyType in self.KEY_TYPES:
            self.entryMap[keyType] = {}
            self.indexMap[keyType] = {}
        self.parser = self.PARSER_CLASS()
        
        indexfilename = self.getIndexFilename(filename)
        if os.path.exists(indexfilename):
            self.loadIndex(indexfilename)
        else:
            self.createIndex(filename)
            self.saveIndex(indexfilename)

    def addEntry(self,entry):
        self.entries.append(entry)
        for keyType in self.KEY_TYPES:
            key = self.GET_ATTR(entry,keyType)
            self.entryMap[keyType].setdefault(key,[]).append(entry)

    def getKeys(self,keyType=None):
        if keyType is None:
            keyType = self.DEFAULT_KEY_TYPE
        keys = self.indexMap[keyType].keys()
        keys.sort()
        return keys

    def __iter__(self):
        for key in self.getKeys():
            entries = self.getEntries(key)
            for entry in entries:
                yield entry

    def getEntry(self,*args):
        entries = self.getEntries(*args)
        if entries:
            return entries[0]
        return None

    def getEntries(self,key,keyType=None):
        if keyType is None:
            keyType = self.DEFAULT_KEY_TYPE
        try:
            entries = self.entryMap[keyType][key]
        except KeyError:
            self.getEntriesFromIndex(key,keyType)
            entries = self.entryMap[keyType][key]
        return entries

    def getEntriesFromIndex(self,key,keyType):
        offsets = self.indexMap[keyType][key]
        for offset,size in offsets:
            self.infh.seek(offset)
            recordString = self.infh.read(size)
            record = self.parser.parse_str(recordString)
            self.addEntry(record)

    def getIndexFilename(self,filename):
        return filename + '.idx'

    def createIndex(self,filename):
        self.infh = file(filename)
        self.indexMap[self.FILENAME_KEY] = filename
        iterator = self.ITERATOR_CLASS(self.infh,self.parser)
        while True:
            offset = iterator._uhandle.tell()
            record = iterator.next()
            size = iterator._uhandle.tell() - offset
            
            if record is None:
                break
            for keyType in self.KEY_TYPES:
                self.indexMap[keyType].setdefault(self.GET_ATTR(record,keyType),[]).append((offset,size))
        
    def loadIndex(self,indexfilename):
        self.indexMap = pickle.load(file(indexfilename))
        self.infh = file(self.indexMap[self.FILENAME_KEY])

    def saveIndex(self,indexfilename):
        pickle.dump(self.indexMap,file(indexfilename,'w'),pickle.HIGHEST_PROTOCOL)
