#! /usr/bin/env python
# $Id: multifile.py,v 1.2 2004/09/22 04:31:52 mliang Exp $

from iterutils import xzip

class MultiFile:
    def __init__(self,*filenames,**kwargs):
        self.filenames = filenames
        self.mode = kwargs.get('mode','r')

        self.fhList = [file(filename,self.mode) for filename in self.filenames]

    def __iter__(self):
        return iter(xzip(*self.fhList))

    def write(self,*lines):
        for line,fh in xzip(lines,self.fhList):
            fh.write(line)

