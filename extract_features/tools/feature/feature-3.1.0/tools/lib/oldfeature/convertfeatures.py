#! /usr/bin/env python

# $Id: convertfeatures.py,v 1.3 2004/02/22 21:22:46 mliang Exp $
# Copyright (c) 2002.  Mike Liang.  All rights reserved
#
# Converts FEATURE .sitedataaf files to tab delimited files for Cluster
#

import re
import sys

def LoadWeights(filename):
    weights = {}
    FIELDS_RE = re.compile(r'{.*"(?P<name>.*?)",.*?,(?P<scale>.*)}')
    
    file = open(filename)

    # find start of property list
    while 1:
        line = file.readline()
        if not line:
            break
        if line.find("[]") != -1:
            break

    # process property list
    while 1:
        line = file.readline()
        if not line:
            break
        if line.find("};") != -1:
            break
        
        m = FIELDS_RE.search(line)
        if not m:
            continue

        # add weight to dictionary
        name = m.group('name').strip()
        scale = float(m.group('scale'))
        weights[name] = scale

    file.close()
    
    return weights


# Gets: property name, volume number, site values, nonsite values
FIELDS_RE = re.compile(r'\(:PROPERTY (?P<name>.+?)\).*\(:VOLUME (?P<volume>.+?)\).*\(:SITE-VALUES \((?P<sites>.*?)\)\).*\(:NONSITE-VALUES \((?P<nonsites>.*?)\)\)')


# Process command line arguments
args = sys.argv[1:]
if len(args) not in (1,2):
    print "Usage: %s <sitedataaf> [<propertylist.cc>]" % sys.argv[0]
    sys.exit(2)
    
filename = args[0]
weights = None
if len(args) == 2:
    weightfile = args[1]
    weights = LoadWeights(weightfile)
file = open(filename)

output = sys.stdout

header = 0
for line in file.readlines():
    # match fields
    m = FIELDS_RE.search(line)
    if not m:
        continue

    # Parse fields
    name = m.group('name')
    volume = m.group('volume')
    sites = map(float, m.group('sites').split())
    nonsites = map(float, m.group('nonsites').split())

    # print header first time
    if not header:
        # Print header
        output.write('UNIQID')
        output.write('\tNAME')
        if weights:
            output.write('\tGWEIGHT')
        for idx in range(len(sites)):
            output.write('\tSITE%d' % idx)
        for idx in range(len(nonsites)):
            output.write('\tNONSITE%d' % idx)
        output.write('\n')
        header = 1

    # print out values
    output.write('%s-%s' % (name, volume))
    output.write('\t%s-%s' % (name, volume))
    if weights:
        output.write('\t%g' % (1.0/weights.get(name,1)))
    for entry in sites:
        output.write('\t%g' % entry)
    for entry in nonsites:
        output.write('\t%g' % entry)
    output.write('\n')
