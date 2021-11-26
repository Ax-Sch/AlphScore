#! /usr/bin/env python
# $Id: iterutils.py,v 1.1 2004/09/22 04:31:52 mliang Exp $

# =========================================================================
# Iter Functions
# =========================================================================
def xzip(*iterators):
    """Iterative version of builtin 'zip'."""
    iterators = map(iter, iterators)
    while 1:
        yield tuple([x.next() for x in iterators])
