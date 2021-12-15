# $Id: sitedata.py,v 1.1 2002/04/05 04:11:50 mliang Exp $
# Copyright (c) 2002 Mike Liang.  All rights reserved.

# Feature .sitedataaf file parser

import cStringIO
import sitedata_format
from xml.sax import handler


class SiteDataFileParser(handler.ContentHandler):
    def __init__(self):
        handler.ContentHandler.__init__(self)
        self.data = None

        self.obj_stack = []
        self.curtext = ""

        self.parser = sitedata_format.format.make_parser()
        self.parser.setContentHandler(self)
        self.parser.setErrorHandler(handler.ErrorHandler())

    def parseStream(self, stream, data):
        self.data = data
        self.parser.parseFile(stream)

    def startElement(self, name, attrs):
        self.curtext = ""
        if name == "record":
            self.obj_stack.append(SiteDataRecord())
        elif name == "sitevalues":
            self.obj_stack.append([])
        elif name == "nonsitevalues":
            self.obj_stack.append([])
        elif name == "tresults":
            self.obj_stack.append(TResultRecord())

    def characters(self, content):
        self.curtext += content

    def endElement(self, name):
        try: top_obj = self.obj_stack[-1]
        except: top_obj = None
        if name == "protein": top_obj.protein = self.curtext
        elif name == "property": top_obj.property = self.curtext
        elif name == "collector": top_obj.collector = self.curtext
        elif name == "volume": top_obj.volume = int(self.curtext)
        elif name == "siteval": top_obj.append(float(self.curtext))
        elif name == "nonsiteval": top_obj.append(float(self.curtext))
        elif name == "plevel": top_obj.plevel = float(self.curtext)
        elif name == "tvalue": top_obj.tvalue = float(self.curtext)
        elif name == "dof": top_obj.dof = float(self.curtext)
        elif name == "mean": top_obj.mean = float(self.curtext)
        elif name == "stdev": top_obj.stdev = float(self.curtext)
        elif name == "ranksum": top_obj.ranksum = float(self.curtext)
        elif name == "numsmall": top_obj.numsmall = int(self.curtext)
        elif name == "numbig": top_obj.numbig = int(self.curtext)
        elif name == "sitevalues":
            values = self.obj_stack.pop()
            self.obj_stack[-1].sitevalues = values
        elif name == "nonsitevalues":
            values = self.obj_stack.pop()
            self.obj_stack[-1].nonsitevalues = values
        elif name == "tresults":
            value = self.obj_stack.pop()
            self.obj_stack[-1].tresult = value
        elif name == "record":
            value = self.obj_stack.pop()
            self.data.records.append(value)
            
        self.curtext = ""


class TResultRecord:
    def __init__(self):
        pass

    def __str__(self):
        sio = cStringIO.StringIO()
        print >>sio, '(:T-RESULT (',
        print >>sio, '(:T-VALUE %g)' % self.tvalue,
        print >>sio, '(:DEGREES-OF-FREEDOM %g)' % self.dof,
        print >>sio, '(:MEAN %g)' % self.mean,
        print >>sio, '(:SD %g)' % self.stdev,
        print >>sio, '(:SUM-OF-RANKS %g)' % self.ranksum,
        print >>sio, '(:NS %d)' % self.numsmall,
        print >>sio, '(:NB %d)' % self.numbig,
        print >>sio, '))',
        return sio.getvalue()

class SiteDataRecord:
    def __init__(self):
        pass

    def __str__(self):
        sio = cStringIO.StringIO()
        print >>sio, '(',
        print >>sio, '(:PROTEIN %s)' % self.protein,
        print >>sio, '(:PROPERTY %s)' % self.property,
        print >>sio, '(:COLLECTOR %s)' % self.collector,
        print >>sio, '(:VOLUME %d)' % self.volume,
        print >>sio, '(:SITE-VALUES (',
        for val in self.sitevalues:
            print >>sio, '%g' % val,
        print >>sio, '))',
        print >>sio, '(:NONSITE-VALUES (',
        for val in self.nonsitevalues:
            print >>sio, '%g' % val,
        print >>sio, '))',
        print >>sio, '(:P-LEVEL %g)' % self.plevel,
        print >>sio, self.tresult,
        print >>sio, ')',
        return sio.getvalue()

class SiteDataFile:
    def __init__(self, filename = None):
        self.records = []
        
        if filename:
            self.loadFile(filename)

    def loadFile(self, filename):
        SiteDataFileParser().parseStream(open(filename), self)

    def __str__(self):
        sio = cStringIO.StringIO()
        print >>sio, '('
        for record in self.records:
            print >>sio, record
        print >>sio, ')'
        return sio.getvalue()
