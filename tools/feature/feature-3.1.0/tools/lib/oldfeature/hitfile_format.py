# $Id: hitfile_format.py,v 1.1 2002/04/05 04:11:50 mliang Exp $
# Copyright (c) 2002 Mike Liang.  All rights reserved.

# .hits format for Martel

from format_utils import *
from Martel import *

coord_expr = Str('(') + Opt(WhiteSpace) + \
             Float('x') + WhiteSpace + \
             Float('y') + WhiteSpace + \
             Float('z') + Opt(WhiteSpace) + Str(')')

hit_expr = Group('hit', Str('(') + Opt(WhiteSpace) + coord_expr + \
           Opt(WhiteSpace) + Float('score') + Opt(WhiteSpace) + Str(')'))

format = Opt(WhiteSpace) + Rep(hit_expr + Opt(WhiteSpace))
