
# $Id: callable.py,v 1.2 2001/10/28 23:00:06 mliang Exp $
#
# Callable interface (wrapper)
#

# ----------------------------------------------------------------------
class Callable:
# ----------------------------------------------------------------------
    def __init__(self, anycallable):
        self.__call__ = anycallable



