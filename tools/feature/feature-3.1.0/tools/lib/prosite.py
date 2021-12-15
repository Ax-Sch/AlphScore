#! /usr/bin/env python
# $Id: prosite.py,v 1.3 2004/09/06 22:08:42 mliang Exp $

import sys,os
import re
from Bio import Prosite
from Bio.Prosite.Pattern import prosite_to_grouped_re
from indexdb import IndexedDatabase

PROSITE_DIR = os.environ.get('PROSITE_DIR','/Users/mliang/db/prosite')
PROSITE_DATABASE_FILE = os.path.join(PROSITE_DIR,'prosite.dat')

PrositeEntry = Prosite.Record

class SearchablePrositeEntry:
    def __init__(self,prositeEntry):
        if prositeEntry.type != 'PATTERN':
            raise ValueError('Not PATTERN type')
        self.prositeEntry = prositeEntry
        self.repattern = re.compile(prosite_to_grouped_re(prositeEntry.pattern),re.I)

    def __getattr__(self,name):
        if name in ['search','findall','finditer']:
            return getattr(self.repattern,name)
        return getattr(self.prositeEntry,name)


class PrositeDatabase(IndexedDatabase):
    # Customize the IndexedDatabase
    KEY_TYPES = ['name']
    PARSER_CLASS = Prosite.RecordParser
    ITERATOR_CLASS = Prosite.Iterator
    
    def __init__(self,filename=None):
        IndexedDatabase.__init__(self,filename)
        

class SearchablePrositeDatabase:
    def __init__(self,prositeDB):
        self.db = prositeDB

    def getEntry(self,*args):
        entries = self.getEntries(*args)
        if entries:
            return entries[0]
        return None

    def getEntries(self,*args):
        entries = self.db.getEntries(*args)
        if type(entries) == list:
            return [SearchablePrositeEntry(e) for e in entries]
        return entries

    def __iter__(self):
        for e in self.db:
            yield SearchablePrositeEntry(e)


prositeDB = PrositeDatabase(PROSITE_DATABASE_FILE)
