#! /usr/bin/env python

# Represents the FEATURE .site class
# Provides methods to export to other formats
import sys,os
import re
from cStringIO import StringIO
from utils import fhfilename

class Site:
    SITE_TYPE = 'site'
    NONSITE_TYPE = 'nonsite'

    def __init__(self):
        pass

    def fieldnames(self):
        return ("pdbid","x","y","z","label","comment")

    def __iter__(self):
        return iter([self.__dict__[field] for field in self.fieldnames()])

    def __repr__(self):
        return str(tuple(self))

    def coord(self):
        return (self.x,self.y,self.z)

RELABEL_MAP = {'T':Site.SITE_TYPE,'NIL':Site.NONSITE_TYPE}
def relabel(label):
    return RELABEL_MAP.get(label,label)

SITEFILE_HEADER = re.compile(r'''
\(:SITE-NAME\s+"(?P<name>[^"]*)"\s*\)\s*
\(:SITE-RADIUS\s+(?P<radius>[0-9.+-]+)\s*\)
''', re.VERBOSE)

SITEFILE_RECORD = re.compile(r'''
\(\s*"(?P<pdbid>[^"]+)"\s+
X\s+(?P<x>[0-9.+-]+)\s+
Y\s+(?P<y>[0-9.+-]+)\s+
Z\s+(?P<z>[0-9.+-]+)\s+
(?P<label>\w+)\s*
\)
(\s*\#\s*(?P<comment>.*))?
''', re.VERBOSE)

class SiteFile:
    def __init__(self,filename=None,fh=None):
        self.filename = "<None>"
        self.name = "<None>"
        self.radius = "0.0"
        self.sites = []

        self.bypdb = {}
        self.bylabel = {}

        if fh:
            self.filename = fhfilename(fh)
            if not self.filename:
                self.filename='<FH>'
            self.load(fh)
        elif filename:
            self.loadfile(filename)

    def loadfile(self,filename):
        self.filename = filename
        fh = file(filename)
        self.load(fh)
        fh.close()

    def load(self,infh):
        data = infh.read()
        m = SITEFILE_HEADER.search(data)
        if m:
            self.name = m.group('name')
            self.radius = float(m.group('radius'))
        for m in SITEFILE_RECORD.finditer(data):
            site = Site()
            site.pdbid = m.group('pdbid')
            site.x = float(m.group('x'))
            site.y = float(m.group('y'))
            site.z = float(m.group('z'))
            site.originalLabel = m.group('label')
            site.label = relabel(m.group('label'))
            site.comment = m.group('comment')
            
            self.sites.append(site)
            self.bypdb.setdefault(site.pdbid,[]).append(site)
            self.bylabel.setdefault(site.label,[]).append(site)

    def getPdbids(self):
        keys = self.bypdb.keys()
        keys.sort()
        return keys

    def getByPdbid(self,pdbid):
        return self.bypdb[pdbid]

    def getByLabel(self,label):
        return self.bylabel[label]

    def writeDelimited(self,outfh=sys.stdout,sep="\t",**options):
        for opt,arg in options.items():
            if opt == 'header' and arg:
                print >>outfh, "# Sites"
                print >>outfh, "#",sep.join(self.sites[0].fieldnames())

        for site in self.sites:
            print >>outfh,sep.join(map(str,site))

    def write(self,outfh=sys.stdout):
        print >>outfh, '''((:SITE-NAME "%(name)s") (:SITE-RADIUS %(radius)g) (:SITES (''' % vars(self)
        for site in self.sites:
            print >>outfh,'''("%(pdbid)s" X %(x)g Y %(y)g Z %(z)g %(originalLabel)s)''' % vars(site)
        print >>outfh, ''')))'''

    def __str__(self):
        sio = StringIO()
        self.write(sio)
        return sio.getvalue()

    def addSites(self,siteList):
        for site in siteList:
            self.sites.append(site)
            self.bypdb.setdefault(site.pdbid,[]).append(site)
            self.bylabel.setdefault(site.label,[]).append(site)

    def getSites(self,key,index='bylabel'):
        return getattr(self,index).get(key,[])


from Ft.Xml import MarkupWriter

def u(s):
    if s is None:
        return s
    return unicode(s)

class SiteXmlDoc:
    def __init__(self,sfile):
        sio = StringIO()
        writer = MarkupWriter(sio,indent=u'yes',omitXmlDeclaration=u'yes')
        writer.startDocument()
        writer.startElement(u'sitefile')
        writer.simpleElement(u'filename',content=u(sfile.filename))
        writer.simpleElement(u'name',content=u(sfile.name))
        writer.simpleElement(u'radius',content=u(sfile.radius))
        writer.startElement(u'sites')
        for site in sfile.sites:
            writer.startElement(u'site')
            writer.simpleElement(u'pdbid',content=u(site.pdbid))
            writer.startElement(u'coord')
            writer.simpleElement(u'x',content=u(site.x))
            writer.simpleElement(u'y',content=u(site.y))
            writer.simpleElement(u'z',content=u(site.z))
            writer.endElement(u'coord')
            writer.simpleElement(u'label',content=u(site.label))
            writer.simpleElement(u'originalLabel',content=u(site.originalLabel))
            writer.simpleElement(u'comment',content=u(site.comment))
            writer.endElement(u'site')
        writer.endElement(u'sites')
        writer.endElement(u'sitefile')
        writer.endDocument()
        self.data = sio.getvalue()

    def __str__(self):
        return self.data
