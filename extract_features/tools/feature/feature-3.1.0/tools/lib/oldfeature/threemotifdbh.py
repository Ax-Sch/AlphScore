# $Id: threemotifdbh.py,v 1.3 2002/04/05 04:09:32 mliang Exp $
# Copyright 2001 (c) Mike Liang.  All rights reserved.

# Interface to 3Motif database

import os.path
from threemotifrecord import ThreeMotifRecord


THREE_MOTIF_DB = "/home/mliang/local/db/3motif/"
DBFILENAME = "pdb.blocksplus.new"


class ThreeMotifDatabase:
    def __init__(self, database, dbfilename):
        self.fieldnames = [ "blockac", "motif" ]
        # Open Handler
        self.database = database
        self.dbfilename = dbfilename
        self.datah = open(os.path.join(self.database, self.dbfilename))
        # Load Indices
        self.index = {}
        for fieldname in self.fieldnames:
            self._loadIndex(fieldname)

    def __del__(self):
        # Close Handler
        self.datah.close()

    def _loadIndex(self, fieldname):
        self.index[fieldname] = {}
        index = self.index[fieldname]

        idxh = open(os.path.join(self.database, "%s.idx" % fieldname))
        for line in idxh.readlines():
            fields = line[:-1].split("\t")
            key = fields[0]
            value = map(int, fields[1:])
            index[key] = value
        idxh.close()
            
    def getThreeMotifRecords(self, key, fieldname="blockac"):
        offsets = self.index[fieldname].get(key,[])
        records = []
        for offset in offsets:
            self.datah.seek(offset)
            tmrec = ThreeMotifRecord(self.datah)
            records.append(tmrec)
        return records


def dbopen(database=THREE_MOTIF_DB, dbfilename=DBFILENAME):
    return ThreeMotifDatabase(database, dbfilename)

