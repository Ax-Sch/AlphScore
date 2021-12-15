#! /usr/bin/env python

# Common utilities
import os,sys
import re
import gzip
import imp

from logutils import *

# =========================================================================
# List Functions
# =========================================================================
def uniq_ordered(alist):
    amap = {}
    idx = 0
    for item in alist:
        alist.setdefault(item,idx)
        idx += 1
    return [a[0] for a in amap.items().sort(lambda a,b:cmp(a[1],b[1]))]

def uniq(alist):
    amap = {}
    for item in alist:
        amap[item] = 1
    return amap.keys()

def uniq_count(alist):
    amap = {}
    for item in alist:
        amap[item] = amap.get(item,0) + 1
    return amap.items()

def sort(alist):
    alist.sort()
    return alist

def sortuniq(alist):
    return sort(uniq(alist))

def sort_copy(alist):
    nlist = alist[:]
    nlist.sort()
    return nlist

def collapse(alist):
    return reduce(lambda a,b:a+b,alist)

def take(alist,indices):
    return [alist[i] for i in indices]

# =========================================================================
# File Functions
# =========================================================================
def openfile(filename):
    if type(filename) != str:
        return filename
    if filename == '-':
        return sys.stdin
    basename,ext = os.path.splitext(filename)
    if ext in ('.gz',):
        infh = gzip.open(filename)
    elif ext in ('.Z',):
        infh = os.popen('gunzip -c %s' % filename)
    elif ext in ('.bz2',):
        infh = os.popen('bunzip2 -c %s' % filename)
    else:
        infh = file(filename)
    return infh

def isfh(obj):
    if not hasattr(obj,'read'):
        return 0
    if not hasattr(obj,'tell'):
        return 0
    if not hasattr(obj,'readline'):
        return 0
    return 1

def fhfilename(fh,default=None):
    for attr in ['filename','name']:
        if hasattr(fh,attr):
            return getattr(fh,attr)
    return default

def pathsearch(pathdirs,filenames):
    if type(pathdirs) == str:
        pathdirs = pathdirs.split(os.pathsep)
    if type(filenames) == str:
        filenames = filenames.split(os.pathsep)

    for filename in filenames:
        for pathdir in pathdirs:
            fullname = os.path.join(pathdir,filename)
            if os.path.isfile(fullname):
                return fullname

# =========================================================================
# String Functions
# =========================================================================
def translate(str,map1,map2):
    assert len(map1) == len(map2)
    table = range(256)
    for l1,l2 in zip(map1,map2):
        table[ord(l1)] = ord(l2)
    return str.translate("".join(map(chr,table)))

def importfile(filename):
    root,base = os.path.split(filename)
    name,ext = os.path.splitext(base)
    args=imp.find_module(name,[root])
    return imp.load_module(name,*args)

CHOMP_RE = re.compile(r'\r?\n$')
def chomp(line):
    return CHOMP_RE.sub('',line)
