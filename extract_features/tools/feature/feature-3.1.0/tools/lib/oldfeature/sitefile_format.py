# $Id: sitefile_format.py,v 1.1 2002/04/05 04:11:50 mliang Exp $
# Copyright (c) 2002 Mike Liang.  All rights reserved.

# .site file format for Martel

from Martel import *
from format_utils import *

sitename_expr = Str('(:SITE-NAME') + WhiteSpace + QuotedString('name') + Opt(WhiteSpace) + Str(')')
siteradius_expr = Str('(:SITE-RADIUS') + WhiteSpace + Float('radius') + Opt(WhiteSpace) + Str(')')
site_expr = Group('site',
                  Str('(') + Opt(WhiteSpace) + QuotedString('pdbid') + WhiteSpace + \
                  Str('X') + WhiteSpace + Float('x') + WhiteSpace + \
                  Str('Y') + WhiteSpace + Float('y') + WhiteSpace + \
                  Str('Z') + WhiteSpace + Float('z') + WhiteSpace + \
                  UnquotedString('type') + Opt(WhiteSpace) + Str(')')
                  )
sitelist_expr = Str('(:SITES') + WhiteSpace + Str('(') + Opt(WhiteSpace) + \
                Rep1(site_expr + Opt(WhiteSpace)) + Str(')') + Opt(WhiteSpace) + \
                Str(')')

format = Opt(WhiteSpace) + Str('(') + Opt(WhiteSpace) + sitename_expr + Opt(WhiteSpace) + \
         siteradius_expr + Opt(WhiteSpace) + \
         sitelist_expr + Opt(WhiteSpace) + Str(')') + Opt(WhiteSpace)

