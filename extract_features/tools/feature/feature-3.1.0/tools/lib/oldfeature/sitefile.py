
# $Id: sitefile.py,v 1.4 2003/04/18 09:25:26 mliang Exp $
# Copyright (c) 2002 Mike Liang.  All rights reserved.

# Feature .site file parser

# =================================================================
# IMPORTS
# =================================================================
import cStringIO
import re



# =================================================================
# GLOBALS
# =================================================================
INFO_RE = re.compile(r'\(:SITE-NAME\s+"(?P<name>.+?)"\s*\)\s*\(:SITE-RADIUS\s+(?P<radius>.+?)\)')
SITE_RE = re.compile(r'\(\s*"(?P<pdbid>\S+?)"\s+X\s+(?P<x>\S+?)\s+Y\s+(?P<y>\S+?)\s+Z\s+(?P<z>\S+?)\s+(?P<type>T|NIL)\s*\)')



# =================================================================
# CLASSES
# =================================================================
class SiteFileParser:
    # =========================================================
    def __init__(self):
        pass

    # =========================================================
    def parseStream(self, file, obj):
        self.obj = obj

        # get info from first line
        line = file.readline()
        m = INFO_RE.search(line)
        if m:
            self.obj.name = m.group('name')
            self.obj.radius = float(m.group('radius'))

        # read one entry per line
        for line in file.xreadlines():
            site = Site()
            if site.parse(line):
                self.obj.sites.setdefault(site.pdbid, []).append(site)


# =================================================================
class Site:
    SITE_TYPE = "T"
    NONSITE_TYPE = "NIL"
    TYPEMAP = {
            SITE_TYPE: "SITE",
            NONSITE_TYPE: "NONSITE"
        }

    # =========================================================
    def __init__(self, pdbid="", x=None, y=None, z=None, type="", weight=1.0):
        self.pdbid = pdbid
        self.x = x
        self.y = y
        self.z = z
        self.type = type
        self.weight = weight

    # =========================================================
    def __str__(self):
        return '("%s" X %g Y %g Z %g %s)' % (self.pdbid, self.x, self.y, self.z, self.type)

    # =========================================================
    def altStr(self):
        return "%s (%g, %g, %g)" % (Site.TYPEMAP[self.type], self.x, self.y, self.z)

    # =========================================================
    def isSite(self):
        return (self.type == Site.SITE_TYPE)

    # =========================================================
    def getLocation(self):
        return (self.x, self.y, self.z)

    # =========================================================
    def parse(self,line):
        m = SITE_RE.search(line)
        if not m:
            return None

        self.pdbid = m.group('pdbid')
        self.x = float(m.group('x'))
        self.y = float(m.group('y'))
        self.z = float(m.group('z'))
        self.type = m.group('type')

        return self
        


# =================================================================
class SiteFile:
    # =========================================================
    def __init__(self, filename = None):
        self.name = ""
        self.radius = 0.0
        self.sites = {}
        self.filename = filename

        if filename:
            self.loadFile(filename)

    # =========================================================
    def addSite(self, site):
        self.sites.setdefault(site.pdbid.upper(), []).append(site)

    # =========================================================
    def getSites(self, pdbid):
        return self.sites.get(pdbid.upper(), [])
    
    # =========================================================
    def loadFile(self, filename):
        SiteFileParser().parseStream(open(filename), self)

    # =========================================================
    def __str__(self):
        sio = cStringIO.StringIO()
        print >>sio, '((:SITE-NAME "%s") (:SITE-RADIUS %g) (:SITES (' % (self.name, self.radius)
        pdbids = self.sites.keys()
        pdbids.sort()
        for pdbid in pdbids:
            for site in self.sites[pdbid]:
                print >>sio, site
        print >>sio, ')))'
        return sio.getvalue()
