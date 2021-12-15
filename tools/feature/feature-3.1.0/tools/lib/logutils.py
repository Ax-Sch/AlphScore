#! /usr/bin/env python
# $Id: logutils.py,v 1.2 2004/05/22 21:50:15 mliang Exp $
# Copyright (c) 2003 Mike Liang. All rights reserved.

# Logging utilities

import sys

_debuglevel = 0
_error_fh = sys.stderr
_warning_fh = sys.stderr
_debug_fh = sys.stderr

# =========================================================================
# Logging Functions
# =========================================================================
def setdebuglevel(level):
    global _debuglevel
    _debuglevel = level

def printfh(fh,*args):
    print >>fh," ".join(map(str,args))

def error(*args):
    print >>_error_fh, "Error:"," ".join(map(str,args))

def fatal_error(*args):
    print >>_error_fh, "Fatal Error:"," ".join(map(str,args))
    sys.exit(2)

def warning(*args):
    print >>_warning_fh, "Warning:"," ".join(map(str,args))

def debug(level,*args):
    if _debuglevel >= level:
        print >>_debug_fh, " ".join(map(str,args))

