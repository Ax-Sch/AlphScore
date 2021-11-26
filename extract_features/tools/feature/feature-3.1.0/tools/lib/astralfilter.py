#! /usr/bin/env python

# Filters pdbids by astral homolog

import sys,os
import getopt
import utils
from logutils import *

MISSING_MARKER="***"
DEFAULT_IDENTITY = 30
DEFAULT_OUTFH = sys.stdout
DEFAULT_ID_TYPE = 'pdb'
DEFAULT_DISPLAY_FORMAT = 'single'
DEFAULT_ONLY_UNIQUE = False
DEFAULT_LOOKUP_TYPE = 'pdb'
DISPLAY_FORMATS = ('single','all')
ID_TYPES = ('scop','pdb')
LOOKUP_TYPES = ('pdb','chain','domain')

def getidtype(key):
    idtype = { 4: 'pdb', 5: 'chain', 6: 'domain' }
    return idtype.get(len(key))


class AstralHomologFile:
    HOMOLOG_DIR = os.environ.get('ASTRAL_DIR','.')
    HOMOLOG_FILE_FMT = 'homolog%s.txt'
    SCORE_FILE = 'scores.txt'

    def __init__(self,identity):
        self.loadscores()
        self.identity = identity
        self.filename = os.path.join(self.HOMOLOG_DIR,self.HOMOLOG_FILE_FMT % identity)
        self.loadmap(self.filename)
        self.only_unique=DEFAULT_ONLY_UNIQUE
        self.display_format=DEFAULT_DISPLAY_FORMAT
        self.id_type=DEFAULT_ID_TYPE


    def loadmap(self,filename):
        self.domainidToCanonicalid = {}
        self.pdbidToDomainid = {}
        self.chainidToDomainid = {}
        for line in file(filename):
            domainids = line[:-1].lower().split('\t')
            canonicalid = domainids[0]
            for domainid in domainids:
                pdbid = domainid[:4]
                chainid = domainid[:5]
                self.pdbidToDomainid.setdefault(pdbid,[]).append(domainid)
                self.chainidToDomainid.setdefault(chainid,[]).append(domainid)
                self.domainidToCanonicalid.setdefault(domainid,[]).append(canonicalid)


    def loadscores(self):
        filename = os.path.join(self.HOMOLOG_DIR,self.SCORE_FILE)
        self.scores = {}
        for line in file(filename):
            domainid,score = line.split()
            self.scores[domainid] = float(score)


    def __getitem__(self,key):
        """ returns the canonical domainid of a given domainid
            if multiple canonical id, return the first - assumes 
            already ordered by best group """
        return self.domainidToCanonicalid[key.lower()][0]


    def get(self,key,default=None):
        """ returns the canonical domainid of a given domainid
            if multiple canonical id, return the first - assumes 
            already ordered by best group """
        ids = self.domainidToCanonicalid.get(key.lower())
        if ids is None:
            return default
        return ids[0]



    def getbytype(self,id,idtype=None):
        """ expands the id to all domainid, then returns canonical 
            ids of those domainids """
        if idtype is None:
            idtype = getidtype(id)

        id = id.lower()

        # perform domainid expansion
        ids = self.expand(id,idtype)

        # for each of the expanded domainids
        retval = []
        for id in ids:
            # if exists canonical, add it to return list
            try:
                retval.append(self[id])
            except KeyError:
                pass
        return retval


    def expand(self,id,idtype=None):
        if idtype is None:
            idtype = getidtype(id)

        id = id.lower()
        if idtype == 'pdb':
            domainids = self.pdbidToDomainid.get(id,[])
        elif idtype == 'chain':
            domainids = self.chainidToDomainid.get(id,[])
            # if not found, maybe chainid is '.'
            if not domainids:
                id = id[:4]+'.'
                domainids = self.chainidToDomainid.get(id,[])
        else:
            domainids = [id]
        return domainids


    def clusterbycanonical(self,ids,idtype=None):
        if idtype is None:
            idtype = getidtype(ids[0])

        def entrysorter(a,b):
            return cmp(self.scores.get(a[0]),self.scores.get(b[0]))
        retdict = {}
        for id in ids:
            domainids = self.expand(id,idtype)
            for domainid in domainids:
                # get first canonicalid referred by domainid
                canonicalid = self.get(domainid)
                # add domainid (and referring id)
                retdict.setdefault(canonicalid,[]).append((domainid,id))
        # sort entries by score
        for entries in retdict.values():
            entries.sort(entrysorter)
        return retdict
        

    def getcanonicalpdb(self,pdbids):
        uniq_map = {}

        for pdbid in pdbids:
            homolog = self.getbytype(pdbid,self.lookup_type)
            if homolog:
                homolog = utils.sort(utils.uniq(homolog))
                if self.id_type == 'pdb':
                    homolog = [id[:4] for id in homolog]
                if self.display_format == 'all':
                    homolog = "\t".join(homolog)
                else:
                    homolog = homolog[0]
            else:
                homolog = pdbid+MISSING_MARKER
            if not self.only_unique or homolog not in uniq_map:
                line = "%s\t%s" % (pdbid,homolog)
                yield line
            uniq_map[homolog] = 1


def eusage():
    print "Usage: %s [OPTIONS] PDBID..." % os.path.basename(sys.argv[0])
    print """
Options:
    -f FILENAME
        Load PDBIDs from FILENAME
    -n IDENTITY
        Use identity cutoff IDENTITY (Default: %(DEFAULT_IDENTITY)s)
    -d DISPLAY_FORMAT
        Use DISPLAY_FORMAT (Default: %(DEFAULT_DISPLAY_FORMAT)s)
        formats: %(DISPLAY_FORMATS)s
    -u  Display uniq (Default: %(DEFAULT_ONLY_UNIQUE)s)
    -t ID_TYPE
        Use ID_TYPE  (Default: %(DEFAULT_ID_TYPE)s)
        types: %(ID_TYPES)s
    -l LOOKUP_TYPE
        Use LOOKUP_TYPE (Default: %(DEFAULT_LOOKUP_TYPE)s)
        types: %(LOOKUP_TYPES)s
""" % globals()
    sys.exit(1)

def main():
    identity = DEFAULT_IDENTITY
    outfh = DEFAULT_OUTFH
    display_format = DEFAULT_DISPLAY_FORMAT
    only_unique = DEFAULT_ONLY_UNIQUE
    id_type = DEFAULT_ID_TYPE
    lookup_type = DEFAULT_LOOKUP_TYPE
    pdbids = []
    args = sys.argv[1:]
    try:
        optlist,args = getopt.getopt(args,"f:n:d:ut:al:")
    except getopt.GetoptError,e:
        eusage()
    for opt,arg in optlist:
        if opt in ['-f']:
            if not os.path.exists(arg):
                fatal_error("Could not open file",arg)
            pdbids.extend([x.strip() for x in file(arg)])
        elif opt in ['-n']:
            identity = int(arg)
        elif opt in ['-d']:
            if arg in DISPLAY_FORMATS:
                display_format = arg
            else:
                warning('Invalid display format',arg)
        elif opt in ['-u']:
            only_unique = True
        elif opt in ['-a']:
            display_format = 'all'
        elif opt in ['-t']:
            if arg in ID_TYPES:
                id_type = arg
            else:
                warning('Invalid id type',arg)
        elif opt in ['-l']:
            if arg in LOOKUP_TYPES:
                lookup_type = arg
            else:
                warning('Invalid lookup type',arg)
    pdbids.extend(args)
    if not pdbids:
        pdbids.extend([x.strip() for x in sys.stdin])

    try:
        homologfile = AstralHomologFile(identity)
    except IOError:
        fatal_error("Could not find homolog file for identity",identity)

    homologfile.only_unique = only_unique
    homologfile.display_format = display_format
    homologfile.id_type = id_type
    homologfile.lookup_type = lookup_type
    for line in homologfile.getcanonicalpdb(pdbids):
        print >>outfh,line


if __name__ == '__main__':
    main()
