#! /usr/bin/env python

# Represents the FEATURE .sitedataaf file
# Provides methods to export to other formats

import sys
import re
from cStringIO import StringIO
from utils import fhfilename

STATFILE_RECORD = re.compile(r'''
\(:PROTEIN\ NIL\)\s+
\(:PROPERTY\ (?P<property>[^)]+)\)\s+
\(:COLLECTOR\ SHELL\)\s+
\(:VOLUME\ (?P<volume>[^)]+)\)\s+
\(:SITE-VALUES\ \((?P<sites>[^)]*)\)\)\s+
\(:NONSITE-VALUES\ \((?P<nonsites>[^)]*)\)\)\s+
\(:P-LEVEL\ (?P<pvalue>[^)]+)\)\s+
\(:T-RESULT\s+
\(\(:T-VALUE\ (?P<tvalue>[^)]+)\)\s+
\(:DEGREES-OF-FREEDOM\ (?P<dof>[^)]+)\)\s+
\(:MEAN\ (?P<mean>[^)]+)\)\s+
\(:SD\ (?P<stdev>[^)]+)\)\s+
\(:SUM-OF-RANKS\ (?P<sumOfRanks>[^)]+)\)\s+
\(:NS\ (?P<numSmall>[^)]+)\)\s+
\(:NB\ (?P<numBig>[^)]+)\)\)\)\)
''', re.VERBOSE)

def sort(alist):
    blist = alist[:]
    blist.sort()
    return blist

class StatFileRecord:
    def __init__(self, matchobj):
        attr = matchobj.groupdict()
        self.property = attr['property'].strip()
        self.volume = int(attr['volume'])
        self.sites = map(float,attr['sites'].split())
        self.nonsites = map(float,attr['nonsites'].split())
        self.pvalue = float(attr['pvalue'])
        self.tvalue = float(attr['tvalue'])
        self.dof = float(attr['dof'])
        self.mean = float(attr['mean'])
        self.stdev = float(attr['stdev'])
        self.sumOfRanks = float(attr['sumOfRanks'])
        self.numSmall = float(attr['numSmall'])
        self.numBig = float(attr['numBig'])

class StatFile:
    def __init__(self,filename=None,fh=None):
        self.records = []
        self.recordsOrdered = []
        self.recordsByName = {}
        self.recordsByVolume = {}
        self.numProperties = 0
        self.numVolumes = 0
        self.numSites = 0
        self.numNonSites = 0
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

        lastPropertyName = None
        for m in STATFILE_RECORD.finditer(data):
            record = StatFileRecord(m)
            if record.property != lastPropertyName:
                lastPropertyName = record.property
                self.numProperties += 1
                self.recordsOrdered.append([])
            self.records.append(record)
            self.recordsOrdered[-1].append(record)
            self.recordsByName.setdefault(record.property,[]).append(record)
            self.recordsByVolume.setdefault(record.volume,[]).append(record)

        if self.recordsOrdered:
            self.numVolumes = len(self.recordsOrdered[0])
            if self.numVolumes:
                self.numSites = len(self.records[0].sites)
                self.numNonSites = len(self.records[0].nonsites)

        self.properties = []
        for i,record in enumerate(self.records):
            if record.property not in self.properties:
                self.properties.append(record.property)

    def GetProperties(self):
        return self.properties

    def writeDelimited(self,outfh,delim="\t"):
        # Write Header
        outfh.write("PROPERTY-VOLUME")
        if self.numSites:
            outfh.write(delim)
            outfh.write(delim.join(["SITE_%d" % i for i in range(self.numSites)]))
        if self.numNonSites:
            outfh.write(delim)
            outfh.write(delim.join(["NONSITE_%d" % i for i in range(self.numNonSites)]))
        outfh.write("\n")
        # Write values
        for volume in sort(self.recordsByVolume.keys()):
            for record in self.recordsByVolume[volume]:
                outfh.write("%s-%d" % (record.property,record.volume))
                if record.sites:
                    outfh.write(delim)
                    outfh.write(delim.join(["%g" % v for v in record.sites]))
                if record.nonsites:
                    outfh.write(delim)
                    outfh.write(delim.join(["%g" % v for v in record.nonsites]))
                outfh.write("\n")


MARKER_BLANK=" "
MARKER_LOW="."
MARKER_HIGH="O"
MARKER_MAP = { -1: MARKER_LOW, 0: MARKER_BLANK, 1: MARKER_HIGH }

PROPERTY_FMT = "%32s "
VOLUME_FMT = "%2s"
def ShowModel(sfile,pthresh=0.01,outfh=sys.stdout,ignoreNonSig=1):
    def sigtype(record):
        if record.pvalue <= pthresh:
            if record.tvalue < 0:
                return -1
            elif record.tvalue > 0:
                return 1
        return 0
        
    # Write Header
    outfh.write(PROPERTY_FMT % " ")
    for volNum in range(sfile.numVolumes):
        outfh.write(VOLUME_FMT % volNum)
    outfh.write("\n")

    # Display Markers
    for recordList in sfile.recordsOrdered:
        propName = recordList[0].property
        significance = [sigtype(record) for record in recordList]
        if ignoreNonSig and 1 not in significance and -1 not in significance:
            continue
        outfh.write(PROPERTY_FMT % propName)
        for sig in significance:
            outfh.write(VOLUME_FMT % MARKER_MAP[sig])
        outfh.write("\n")




from Ft.Xml import MarkupWriter

def u(s):
    if s is None:
        return s
    return unicode(s)

class StatXmlDoc:
    def __init__(self,stfile):
        sio = StringIO()
        writer = MarkupWriter(sio,indent=u'yes',omitXmlDeclaration=u'yes')
        writer.startDocument()
        writer.startElement(u'statfile')
        writer.simpleElement(u'filename',content=u(stfile.filename))

        writer.startElement(u'properties')
        for propertyIndex,record in enumerate(stfile.records):
            writer.startElement(u'property')
            writer.attribute(u'index',u(propertyIndex))
            writer.simpleElement(u'name',content=u(record.property))
            writer.simpleElement(u'volume',content=u(record.volume))
            writer.simpleElement(u'pvalue',content=u(record.pvalue))
            writer.simpleElement(u'tvalue',content=u(record.tvalue))
            writer.simpleElement(u'dof',content=u(record.dof))
            writer.simpleElement(u'mean',content=u(record.mean))
            writer.simpleElement(u'stdev',content=u(record.stdev))
            writer.simpleElement(u'sumOfRanks',content=u(record.sumOfRanks))
            writer.simpleElement(u'numSmall',content=u(record.numSmall))
            writer.simpleElement(u'numBig',content=u(record.numBig))
            writer.endElement(u'property')
        writer.endElement(u'properties')

        writer.startElement(u'sites')
        for siteIndex in range(stfile.numSites):
            writer.startElement(u'site')
            writer.attribute(u'index',u(siteIndex))
            writer.attribute(u'label',u'site')
            for propertyIndex,record in enumerate(stfile.records):
                writer.startElement(u'value')
                writer.attribute(u'index',u(propertyIndex))
                writer.text(u(record.sites[siteIndex]))
                writer.endElement(u'value')
            writer.endElement(u'site')
        for nonSiteIndex in range(stfile.numNonSites):
            writer.startElement(u'site')
            writer.attribute(u'index',u(nonSiteIndex))
            writer.attribute(u'label',u'nonsite')
            for propertyIndex,record in enumerate(stfile.records):
                writer.startElement(u'value')
                writer.attribute(u'index',u(propertyIndex))
                writer.text(u(record.nonsites[nonSiteIndex]))
                writer.endElement(u'value')
            writer.endElement(u'site')
        writer.endElement(u'sites')

        writer.endElement(u'statfile')
        writer.endDocument()
        self.data = sio.getvalue()

    def __str__(self):
        return self.data






StatFile.writeModel = ShowModel


if __name__=="__main__":
    import sys,os
    def eusage():
        print "Usage: %s FILENAME" % os.path.basename(sys.argv[0])
        sys.exit(1)

    args = sys.argv[1:]
    if len(args) != 1:
        eusage()
    filename = args[0]

    sf = StatFile(filename)
    ShowModel(sf)
