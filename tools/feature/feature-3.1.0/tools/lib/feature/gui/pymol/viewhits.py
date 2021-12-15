'''
$Id: viewhits.py,v 1.2 2005/05/12 05:47:33 mliang Exp $
Copyright (c) 2003 Mike Liang.  All rights reserved.

DESCRIPTION

    "viewhits" will visualize hitfiles from Feature.  By default, it
    will display the top 100 scores.  Color can be specified by name
    or parenthesized, comma separted RGB tuple.  There is no way to
    specify "top X" scores yet, only cutoff scores.  

USAGE

    viewhits <hitfile> [, <cutoff>] [, <color>] [, <radius>]

'''

from pymol import cmd
from pymol.cgo import *
from pymol import viewing
from pymol import querying
import os
import re


# =================================================================
class _Hit:
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
        if isinstance(other, _Hit):
            return cmp(self.score, other.score)
        return cmp(self.score, other)

    # =========================================================
    def getScore(self):
        return self.score

    # =========================================================
    def getLocation(self):
        return (self.x, self.y, self.z)


# approximate floating point expression
FLOAT_RE = re.compile('[0-9.-]+')

# =================================================================
class _HitFile:
    # =========================================================
    def __init__(self, filename = None):
        self.hits = []

        if filename:
            self.loadFile(filename)

    # =========================================================
    def loadFile(self, filename):
        for line in open(filename).xreadlines():
            hit = _Hit()
            hit.x, hit.y, hit.z, hit.score = map(float,FLOAT_RE.findall(line))
            self.hits.append(hit)
        
    # =========================================================
    def __str__(self):
        sio = cStringIO.StringIO()
        for hit in self.hits:
            print >>sio, hit
        return sio.getvalue()


# =================================================================
def _viewhits(hitfile, cutoff=None, color=(1.0,0.0,0.0), radius=0.3 ):
    hit_objs = []
    try:
        hf = _HitFile(hitfile)
    except IOError:
        print "Error opening file: %s" % hitfile
        return
    if not hf: return

    print 'hitfile',hitfile,'cutoff',cutoff,'color',color,'radius',radius

    if cutoff == None:    
        hf.hits.sort()
        hits = hf.hits[-100:]
        print 'Cutoff Score:',hits[0].score
    else:
        hits = filter(lambda h,c=cutoff:h.score>=c, hf.hits)

    print "Number of Hits:", len(hits)
    hit_objs += [ COLOR, color[0], color[1], color[2] ]    
    for hit in hits:
        hit_objs += [ SPHERE, hit.x, hit.y, hit.z, radius ]

    pdbid = os.path.splitext(os.path.basename(hitfile))[0].lower()
    obj_name = '%s hits' % pdbid

    # supress view reset on load
    cur_view = viewing.get_view(0)
    cmd.delete(obj_name)
    cmd.load_cgo(hit_objs, obj_name)
    viewing.set_view(cur_view)
    
    if hf.hits:
        histogram([h.score for h in hf.hits], 10)


# =================================================================
def histogram(alist, numbins):
    minval = min(alist)
    maxval = max(alist)
    binsize = float(maxval-minval)/numbins
    bins = [0]*numbins

    for item in alist:
        idx = int((item-minval)/binsize)
        if idx == numbins: idx = numbins-1
        bins[idx] += 1

    print 'Min:',minval, 'Max:',maxval
    for idx in range(numbins):
        print '%6.2f-%6.2f: %d' % (minval+idx*binsize, minval+(idx+1)*binsize, bins[idx])

    return bins, minval, binsize


def _color_lookup(color):
    color_tuple = cmd._cmd.get_color(color,0)
    if color_tuple is None:
        try:
            color_tuple = eval(color)
        except NameError:
            print "Unknown color: %s" % color
            return None
    return color_tuple


# =================================================================
def viewhits(hitfile, cutoff='', color='(1.0,0.0,0.0)', radius='0.3'):
    '''
    DESCRIPTION

    "viewhits" will visualize hitfiles from Feature.  By default, it
    will display the top 100 scores.  Color can be specified by name
    or parenthesized, comma separted RGB tuple.  There is no way to
    specify "top X" scores yet, only cutoff scores.  

    USAGE

    viewhits <hitfile> [, <cutoff>] [, <color>] [, <radius>]
    '''
    
    if not cutoff: cutoff = None
    else: cutoff = eval(cutoff)

    # Handle color lookup
    color = _color_lookup(color)
    if color is None:
        return

    radius = eval(radius)
    _viewhits(hitfile, cutoff, color, radius)

def __init__(self=None):
    cmd.extend('viewhits',viewhits)
    

# =================================================================
if __name__ == 'pymol':
    __init__()
