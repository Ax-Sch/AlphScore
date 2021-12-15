
# $Id: hitfile.py,v 1.1 2002/04/05 04:11:50 mliang Exp $
# Copyright (c) 2002 Mike Liang.  All rights reserved.

# Feature .hits file parser

# =================================================================
# IMPORTS
# =================================================================
import cStringIO
import hitfile_format
from xml.sax import handler


# =================================================================
# CLASSES
# =================================================================
class HitFileParser(handler.ContentHandler):
    # =========================================================
    def __init__(self):
        handler.ContentHandler.__init__(self)
        self.data = None

        self.elem_stack = []
        self.obj_stack = []
        self.curtext = ""

        self.parser = hitfile_format.format.make_parser()
        self.parser.setContentHandler(self)
        self.parser.setErrorHandler(handler.ErrorHandler())

    # =========================================================
    def parseStream(self, stream, data):
        self.data = data
        self.parser.parseFile(stream)

    # =========================================================
    def startElement(self, name, attrs):
        self.curtext = ""
        if name == "hit":
            self.obj_stack.append(Hit())

    # =========================================================
    def characters(self, content):
        self.curtext += content

    # =========================================================
    def endElement(self, name):
        try: top_obj = self.obj_stack[-1]
        except: top_obj = None
        if name == "x": top_obj.x = float(self.curtext)
        elif name == "y": top_obj.y = float(self.curtext)
        elif name == "z": top_obj.z = float(self.curtext)
        elif name == "score": top_obj.score = float(self.curtext)
        elif name == "hit": self.data.hits.append(self.obj_stack.pop())


# =================================================================
class Hit:
    # =========================================================
    def __init__(self):
        self.x = None
        self.y = None
        self.z = None
        self.score = None

    # =========================================================
    def __str__(self):
        return '((%g %g %g) %g)' % (self.x, self.y, self.z, self.score)
    
    # =========================================================
    def __cmp__(self, other):
        if isinstance(other, Hit):
            return cmp(self.score, other.score)
        return cmp(self.score, other)

    # =========================================================
    def getScore(self):
        return self.score

    # =========================================================
    def getLocation(self):
        return (self.x, self.y, self.z)


# =================================================================
class HitFile:
    # =========================================================
    def __init__(self, filename = None):
        self.hits = []

        if filename:
            self.loadFile(filename)

    # =========================================================
    def loadFile(self, filename):
        HitFileParser().parseStream(open(filename), self)

    # =========================================================
    def __str__(self):
        sio = cStringIO.StringIO()
        for hit in self.hits:
            print >>sio, hit
        return sio.getvalue()
