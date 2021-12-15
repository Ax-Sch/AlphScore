#! /usr/bin/env python

# Views motifs from motif file of format:
# 0 1      2      3         4       5
# # Match: 1kwm:A (147,160) [56:68] GQNKPAIFMDCGF

# from pymol.cgo import *
# from pymol import viewing
# from pymol import querying
import os
import re


def viewmotif(motiffile, pdbid=''):
    '''
    DESCRIPTION

    "viewmotif" displays motifs from a motif file.

    USAGE

    viewmotif <motiffile> [<pdbid>]
    '''
    class Entry:
        pass

    fh = file(motiffile)
    entryMap = {}
    for line in fh.readlines():
        entry = Entry()
        fields = line.split()
        molid = fields[2]
        span = fields[3]
        resSpan = fields[4]

        entry.pdbid,entry.chainid = molid.split(':')
        entry.start,entry.end = map(int,re.search(r'\((\d+),(\d+)\)',span).groups())
        entry.resStart,entry.resStop = map(int,re.search(r'\[(\d+):(\d+)\]',resSpan).groups())
        entry.seq = fields[5]

        entryMap.setdefault(entry.pdbid,[]).append(entry)

    if not pdbid:
        pdbids = entryMap.keys()
        pdbids.sort()
        for key in pdbids:
            print key,len(entryMap[key])
    else:
        relevantEntries = entryMap[pdbid]
        idx = 0
        for entry in relevantEntries:
            from pymol import cmd
            selname = 'motif%d' % idx
            if not entry.chainid:
                cmd.do('select %s,(resi %d-%d)' % (selname,entry.resStart,entry.resStop))
            else:
                cmd.do('select %s,(resi %d-%d and chain %s)' % (selname,entry.resStart,entry.resStop,entry.chainid))
            cmd.do('show cartoon,%s' % selname)
            cmd.do('color red,%s' % selname)
            print entry.pdbid, entry.chainid, entry.seq
            idx += 1

if __name__ == '__pymol__':
    from pymol import cmd
    cmd.extend('viewmotif',viewmotif)
