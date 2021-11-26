# $Id: sitedata_format.py,v 1.1 2002/04/05 04:11:50 mliang Exp $
# Copyright (c) 2002 Mike Liang.  All rights reserved.

# .sitedataaf format for Martel

from format_utils import *
from Martel import *


def StringProperty(property, groupname):
    return Str('(:%s' % property) + WhiteSpace + \
           AnyString(groupname) + Opt(WhiteSpace) + Str(')')


def FloatProperty(property, groupname):
    return Str('(:%s' % property) + WhiteSpace + \
           Float(groupname) + Opt(WhiteSpace) + Str(')')


def IntegerProperty(property, groupname):
    return Str('(:%s' % property) + WhiteSpace + \
           Integer(groupname) + Opt(WhiteSpace) + Str(')')


def FloatList(elementname=None):
    return Str('(') + Opt(WhiteSpace) + \
           Rep(Float(elementname) + Opt(WhiteSpace)) + Str(')')


protein_expr = StringProperty('PROTEIN', 'protein')
property_expr = StringProperty('PROPERTY', 'property')
collector_expr = StringProperty('COLLECTOR', 'collector')
volume_expr = IntegerProperty('VOLUME', 'volume')
sitevalue_expr = Group('sitevalues', Str('(:SITE-VALUES') + WhiteSpace + \
                 FloatList('siteval') + Opt(WhiteSpace) + Str(')'))
nonsitevalue_expr = Group('nonsitevalues', Str('(:NONSITE-VALUES') + WhiteSpace + \
                    FloatList('nonsiteval') + Opt(WhiteSpace) + Str(')'))
plevel_expr = FloatProperty('P-LEVEL', 'plevel')

tvalue_expr = FloatProperty('T-VALUE', 'tvalue')
dof_expr = FloatProperty('DEGREES-OF-FREEDOM', 'dof')
mean_expr = FloatProperty('MEAN', 'mean')
stdev_expr = FloatProperty('SD', 'stdev')
ranksum_expr = FloatProperty('SUM-OF-RANKS', 'ranksum')
numsmall_expr = IntegerProperty('NS', 'numsmall')
numbig_expr = IntegerProperty('NB', 'numbig')

tresult_fields = (
    tvalue_expr + Opt(WhiteSpace),
    dof_expr + Opt(WhiteSpace),
    mean_expr + Opt(WhiteSpace),
    stdev_expr + Opt(WhiteSpace),
    ranksum_expr + Opt(WhiteSpace),
    numsmall_expr + Opt(WhiteSpace),
    numbig_expr + Opt(WhiteSpace)
    )
tresult_expr = Group('tresults', Str('(:T-RESULT') + WhiteSpace + \
               Str('(') + Opt(WhiteSpace) + Seq(*tresult_fields) + Str(')') + \
               Opt(WhiteSpace) + Str(')'))


record_fields = (
    protein_expr + Opt(WhiteSpace),
    property_expr + Opt(WhiteSpace),
    collector_expr + Opt(WhiteSpace),
    volume_expr + Opt(WhiteSpace),
    sitevalue_expr + Opt(WhiteSpace),
    nonsitevalue_expr + Opt(WhiteSpace),
    plevel_expr + Opt(WhiteSpace),
    tresult_expr + Opt(WhiteSpace)
    )
record_expr = Group('record', Str('(') + Opt(WhiteSpace) + Seq(*record_fields) + Str(')'))

format = Opt(WhiteSpace) + Str('(') + Opt(WhiteSpace) + \
         Rep(record_expr + Opt(WhiteSpace)) + Str(')') + Opt(WhiteSpace)
