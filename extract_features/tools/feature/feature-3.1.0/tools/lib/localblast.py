#! /usr/bin/env python

# Run blast locally

import tempfile
from Bio.Blast import NCBIStandalone
from Bio import File

DEFAULT_BLAST_EXE = '/home/mliang/opt/blast/blastall'
DEFAULT_BLAST_DB = '/home/mliang/opt/blast/db/sprot41'
DEFAULT_BLAST_MODE = 'blastp'

class LocalBlastP:
    def __init__(self,dbname=None,blastexe=None,mode=None,parser=None):
        if dbname is None:
            dbname = DEFAULT_BLAST_DB
        if blastexe is None:
            blastexe = DEFAULT_BLAST_EXE
        if mode is None:
            mode = DEFAULT_BLAST_MODE
        if parser is None:
            parser = NCBIStandalone.BlastParser()
        self.dbname = dbname
        self.blastexe = blastexe
        self.parser = parser
        self.mode = mode

    def blastfile(self,filename):
        # run blast
        b_out,e_info = NCBIStandalone.blastall(self.blastexe,self.mode,self.dbname,filename)
        data = b_out.read()
        if not data:
            raise ValueError, 'BLAST error: %s' % e_info.read()

        return data

    def blast(self,sequence):
        # generate temporary file to contain sequence
        tmpfile = tempfile.NamedTemporaryFile()
        print >>tmpfile, "> tmpfile"
        print >>tmpfile, sequence
        tmpfile.flush()

        # run blast
        data = self.blastfile(tmpfile.name)
        tmpfile.close()

        # parse output
        if self.parser:
            b_record = self.parser.parse(File.StringHandle(data))

        # return record
        return b_record

