#! /usr/bin/env python

# Process sprot42.dat file
# Generate table of ID-OFFSET
# Generate FASTA file

import sys,os
import getopt
import time
from logutils import *


DEFAULT_SPROT_DIR = os.environ.get('SPROT_DIR',"/home/mliang/project1/swissprot")
DEFAULT_SPROT_FILENAME = "sprot42.dat"
DEFAULT_SPROT_INDEXNAME = "sprot42.idx"

class FeatureEntry:
    def __init__(self,line=None):
        if line:
            self.parse(line)

    def parse(self,line):
        line = line[5:]    # get rid of junk in front
        self.name = line[0:8].strip()
        try:
            self.from_res = int(line[9:15])
        except ValueError:
            self.from_res = line[9:15].strip()
        try:
            self.to_res = int(line[16:22])
        except ValueError:
            self.to_res = line[16:22].strip()
        self.description = line[29:70].strip()
        #if there is a feature_id (FTId), store it away
        if line[29:35]==r"/FTId=":
            self.ft_id = line[35:70].strip()[:-1]
        else:
            self.ft_id =""

    def __iadd__(self,other):
        self.description += " "+other.description
        self.description.strip()
        return self

    def __str__(self):
        return "%s\t%s\t%s" % (self.name,self.from_res,self.to_res)

class Record:
    def __init__(self):
        self.entry_name = None
        self.features = []

class RecordParser:
    def __init__(self):
        pass

    def parse(self,infh):
        record = Record()
        sequence = []
        for line in infh:
            if line[:2] == 'ID':
                record.entry_name = line.split()[1]
            elif line[:2] == 'FT':
                ftentry = FeatureEntry(line)
                if not ftentry.name:
                    assert record.features
                    record.features[-1] += ftentry
                else:
                    record.features.append(ftentry)
            elif line[:2] == '  ':
                sequence.append(line[:-1])
            elif line[:2] == '//':
                break
        record.sequence = ''.join(sequence).replace(' ','')
        return record

class LineParser:
    def parse(self,infh):
        lines = []
        for line in infh:
            lines.append(line)
            if line[:2] == '//':
                break
        return lines

class SprotLookup:
    def __init__(self,**kwargs):
        self.sprot_dir = DEFAULT_SPROT_DIR
        self.sprot_filename = DEFAULT_SPROT_FILENAME
        self.sprot_indexname = DEFAULT_SPROT_INDEXNAME
        for key,arg in kwargs.items():
            if key in ['sprot_dir']:
                self.sprot_dir = arg
            elif key in ['sprot_filename']:
                self.sprot_filename = arg
            elif key in ['sprot_indexname']:
                self.sprot_indexname = arg
        self.full_sprot_filename = os.path.join(self.sprot_dir,self.sprot_filename)
        self.full_sprot_indexname = os.path.join(self.sprot_dir,self.sprot_indexname)
        self.sprot_file = file(self.full_sprot_filename)
        self.sprot_parser = RecordParser()
        self.offsets = {}
        self.load_index_file()

    def load_index_file(self):
        debug(0,"Loading index file")
        stime = time.time()
        index_file = file(self.full_sprot_indexname)
        for line in index_file:
            id,offset = line.split('\t',2)
            self.offsets[id] = int(offset)
        etime = time.time() - stime
        debug(0,"Finished loading index file",etime)
        

    def query(self,id):
        if self.setoffset(id) < 0:
            return
        return self.sprot_parser.parse(self.sprot_file)

    def view(self,id):
        if self.setoffset(id) < 0:
            return
        return ''.join(LineParser().parse(self.sprot_file))

    def setoffset(self,id):
        offset = self.offsets.get(id)
        if offset is None:
            return -1
        self.sprot_file.seek(offset)
        return offset

def GetActiveSiteSerine(record):
    class Entry:
        def __init__(self,*args):
            self.sprotName = args[0]
            self.featureName = args[1]
            self.featureStart = args[2]
            self.featureStop = args[3]
            self.featureResidue = args[4]
            self.featureFragment = args[5]

        def __iter__(self):
            return iter([self.sprotName,self.featureName,self.featureStart,self.featureStop,self.featureResidue,self.featureFragment])

        def __str__(self):
            return '\t'.join(map(str,list(self)))

    for item in record.features:
        if item.name == 'ACT_SITE':
            residue = record.sequence[item.from_res-1:item.to_res]
            if residue == 'S':
                frag_start = max(0,item.from_res-1-5)
                frag_end = min(len(record.sequence),item.to_res+5)
                fragment = record.sequence[frag_start:frag_end]
                yield Entry(record.entry_name,item.name,item.from_res,item.to_res,residue,fragment)
            else:
                frag_start = max(0,item.from_res-1-5)
                frag_end = min(len(record.sequence),item.to_res+5)
                fragment = record.sequence[frag_start:frag_end]
                debug(1,'not serine',item,item.from_res,fragment,record)
        else:
            debug(1,'not active',item)

def GetActiveSite(record):
    class Entry:
        def __init__(self,*args):
            self.sprotName = args[0]
            self.featureName = args[1]
            self.featureStart = args[2]
            self.featureStop = args[3]
            self.featureResidue = args[4]
            self.featureFragment = args[5]

        def __iter__(self):
            return iter([self.sprotName,self.featureName,self.featureStart,self.featureStop,self.featureResidue,self.featureFragment])

        def __str__(self):
            return '\t'.join(map(str,list(self)))

    for item in record.features:
        if item.name == 'ACT_SITE':
            residue = record.sequence[item.from_res-1:item.to_res]
            frag_start = max(0,item.from_res-1-5)
            frag_end = min(len(record.sequence),item.to_res+5)
            fragment = record.sequence[frag_start:frag_end]
            yield Entry(record.entry_name,item.name,item.from_res,item.to_res,residue,fragment)

def GetInterestingFeatures(record,goodNames=['ACT_SITE']):
    class Entry:
        def __init__(self,*args):
            self.sprotName = args[0]
            self.featureName = args[1]
            self.featureStart = args[2]
            self.featureStop = args[3]
            self.featureResidue = args[4]
            self.featureFragment = args[5]

        def __iter__(self):
            return iter([self.sprotName,self.featureName,self.featureStart,self.featureStop,self.featureResidue,self.featureFragment])

        def __str__(self):
            return '\t'.join(map(str,list(self)))

    for item in record.features:
        if item.name in goodNames:
            residue = record.sequence[item.from_res-1:item.to_res]
            frag_start = max(0,item.from_res-1-5)
            frag_end = min(len(record.sequence),item.to_res+5)
            fragment = record.sequence[frag_start:frag_end]
            yield Entry(record.entry_name,item.name,item.from_res,item.to_res,residue,fragment)


def main():
    args = sys.argv[1:]
    if not args:
        idlist = [x.strip() for x in sys.stdin]
    else:
        idlist = args
    idlist = [x.upper() for x in idlist]

    lookup = SprotLookup()
    for id in idlist:
        record = lookup.query(id)
        if not record:
            continue
        for x in GetActiveSite(record):
            print '\t'.join(map(str,x))


if __name__ == '__main__':
    main()
