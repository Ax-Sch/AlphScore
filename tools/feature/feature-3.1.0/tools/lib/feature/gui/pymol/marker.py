# $Id: marker.py,v 1.1 2005/05/12 05:47:44 mliang Exp $
# Copyright (c) 2004 Mike Liang. All rights reserved.
#
# PyMol plugin (requires PyMol version 0.96 or above)
# Provides command line extension for adding marker
# Use [ Plugin | Install Plugin... ] to install the webload plugin
#
# For versions 0.95 and earlier, include the following in .pymolrc
#     run /path/to/marker/marker.py

from pymol import cmd
from pymol.cgo import *

def __init__(self=None):
    cmd.extend('marker',marker)
    
def _color_lookup(color):
    color_tuple = cmd._cmd.get_color(color,0)
    if color_tuple is None:
        try:
            color_tuple = eval(color)
        except NameError:
            print "Unknown color: %s" % color
            return None
    return color_tuple

def marker(x,y,z,c='red',r='1.0',p='marker'):
    (x,y,z,radius) = map(eval,(x,y,z,r))
    c = _color_lookup(c)
    if c is None:
        return
    (cr,cg,cb) = c

    color = [ COLOR, cr, cg, cb ]
    marker = [ SPHERE, x, y, z, radius ]

    nameList = cmd.get_names()
    name = p
    index = 1

    while name in nameList:
        name = '%s%d' % (p,index)
        index = index + 1

    v = cmd.get_view()
    cmd.load_cgo(color+marker,name)
    cmd.set_view(v)
    
if __name__ == 'pymol':
    __init__()
