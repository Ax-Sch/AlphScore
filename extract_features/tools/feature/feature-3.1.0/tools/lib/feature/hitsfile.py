#! /usr/bin/env python

# Represents the FEATURE .hits file
# Provides method to export to other formats

import sys,os
import re
from cStringIO import StringIO
from utils import fhfilename


HIT_RECORD = re.compile(r'\(\((?P<coord>[^)]+)\) (?P<score>[^)]+)\)')

class Hit:
    def __init__(self,line=None):
        if line is not None:
            self.parse(line)

    def parse(self,line):
        m = HIT_RECORD.search(line)
        if not m:
            return
        self.x,self.y,self.z = map(float,m.group('coord').split())
        self.score = float(m.group('score'))

    def coord(self):
        return (self.x,self.y,self.z)

class HitsFile:
    def __init__(self,filename=None,fh=None):
        self.hits = []
        if filename:
            self.loadFile(filename)
        elif fh:
            self.filename = fhfilename(fh)
            if not self.filename:
                self.filename='<FH>'
            self.load(fh)

    def loadFile(self,filename):
        self.filename = filename
        fh = file(filename)
        self.load(fh)
        fh.close()

    def load(self,fh):
        data = fh.read()
        self.hits = []
        for m in HIT_RECORD.finditer(data):
            hit = Hit()
            self.hits.append(hit)

            hit.x,hit.y,hit.z = map(float,m.group('coord').split())
            hit.score = float(m.group('score'))

    def writeDelimited(self,outfh=sys.stdout,sep='\t',header=0):
        if header:
            print >>outfh, sep.join(["X","Y","Z","SCORE"])
        for hit in self.hits:
            print >>outfh, sep.join(["%g" % v for v in (hit.x,hit.y,hit.z,hit.score)])

    def write(self,outfh=sys.stdout):
        for hit in self.hits:
            print >>outfh, "((%(x)s %(y)s %(z)s) %(score)s)" % vars(hit)



from Ft.Xml import MarkupWriter

def u(s):
    if s is None:
        return s
    return unicode(s)

class HitsXmlDoc:
    def __init__(self,hfile):
        sio = StringIO()
        writer = MarkupWriter(sio,indent=u'yes',omitXmlDeclaration=u'yes')
        writer.startDocument()
        writer.startElement(u'hitsfile')
        writer.simpleElement(u'filename',content=u(hfile.filename))
        writer.startElement(u'hits')
        for hit in hfile.hits:
            writer.startElement(u'hit')
            writer.simpleElement(u'score',content=u(hit.score))
            writer.startElement(u'coord')
            writer.simpleElement(u'x',content=u(hit.x))
            writer.simpleElement(u'y',content=u(hit.y))
            writer.simpleElement(u'z',content=u(hit.z))
            writer.endElement(u'coord')
            writer.endElement(u'hit')
        writer.endElement(u'hits')
        writer.endElement(u'hitsfile')
        writer.endDocument()
        self.data = sio.getvalue()

    def __str__(self):
        return self.data

