# $Id: format_utils.py,v 1.1 2002/04/05 04:11:50 mliang Exp $
# Copyright (c) 2002 Mike Liang.  All rights reserved.

# Useful variables for Martel formats

from Martel import *


def UnquotedString(name = None):
    """(name = None) -> match an alphanumeric string

    If 'name' is not None, the matching text will be put inside of a
    group of the given name.
    """
    exp = Re(r"[a-zA-Z0-9_-]+")
    if name is None:
        return exp
    else:
        return Group(name, exp)


def QuotedString(name = None):
    """(name = None) -> match quoted string, can be single or double quoted

    If 'name' is not None, the matching text will be put inside of a
    group of the given name.
    """
    if name is None:
        single_exp = Re(r"'[^']*'")
        double_exp = Re(r'"[^"]*"')
    else:
        single_exp = Re(r"'(?P<%s>[^']*)'" % name)
        double_exp = Re(r'"(?P<%s>[^"]*)"' % name)
    return Alt(single_exp, double_exp)


def AnyString(name=None):
    return Alt(QuotedString(name), UnquotedString(name))


WhiteSpace = Re(r"\s+")
