#! /usr/bin/env python

# Converts statfile to scorefile

import sys,os
from math import log
import re
from cStringIO import StringIO
from utils import fhfilename

# Globals/Defaults:
# numBins
NUM_BINS = 5
# lowerBoundFudgeFactor
LOWER_BOUND_FUDGE_FACTOR = -1E-5
# P(Site)
SITE_PRIOR = 0.01
# P(NonSite)
NONSITE_PRIOR = 1-SITE_PRIOR
# minScore
MIN_SCORE = log(SITE_PRIOR)

SCORE_HEADER = re.compile(r'''
\(:NUM-BINS\ (?P<numbins>[^)]+)\)\s+
\(:PROB-SITE\ (?P<psite>[^)]+)\)\s+
\(:PROB-NONSITE\ (?P<pnonsite>[^)]+)\)
''',re.VERBOSE)

SCORE_RECORD = re.compile(r'''
\(\(:PROPERTY\ (?P<property>[^)]+)\)\s+
\(:VOLUME\ (?P<volume>[^)]+)\)\s+
\(:P-VALUE\ (?P<pvalue>[^)]+)\)\s+
\(:BIN-SIZE\ (?P<binsize>[^)]+)\)\s+
\(:LOW-VALUE\ (?P<lowvalue>[^)]+)\)\s+
\(:SCORES\ \((?P<scores>[^)]*)\)\)\)
''',re.VERBOSE)

class ScoreInfo:
    def __init__(self,*args):
        if len(args) == 2:
            propvol,numBins = args
            self.property = propvol.property
            self.volume = propvol.volume
            self.pValue = propvol.pvalue
            self.numBins = numBins
            self.scores = [0.0] * numBins

            self._calculateBinSize(propvol)

    def _calculateBinSize(self,propvol):
        maxSiteVal = max(propvol.sites)
        maxNonSiteVal = max(propvol.nonsites)
        minSiteVal = min(propvol.sites)
        minNonSiteVal = min(propvol.nonsites)

        upperBound = max(maxSiteVal,maxNonSiteVal)
        lowerBound = min(minSiteVal,minNonSiteVal)
        lowerBound += LOWER_BOUND_FUDGE_FACTOR
        binSize = (upperBound-lowerBound)/self.numBins

        self.lowValue = lowerBound
        self.binSize = binSize

    def getBin(self,value):
        bin = int((value-self.lowValue)/self.binSize)
        if bin < 0:
            bin = 0
        elif bin >= self.numBins:
            bin = self.numBins-1
        return bin


class ScoreFile:
    def __init__(self,filename=None,fh=None):
        self.numBins = 0
        self.pSite = 0
        self.pNonSite = 1-self.pSite
        self.scores = []

        if fh:
            self.filename = fhfilename(fh)
            if not self.filename:
                self.filename='<FH>'
            self.load(fh)
        elif filename:
            self.loadfile(filename)

    def loadfile(self,filename):
        self.filename = filename
        fh = file(filename)
        self.load(fh)
        fh.close()

    def load(self,fh):
        data = fh.read()
        m = SCORE_HEADER.search(data)
        if m:
            self.numBins = int(m.group('numbins'))
            self.pSite = float(m.group('psite'))
            self.pNonSite = float(m.group('pnonsite'))
        for m in SCORE_RECORD.finditer(data):
            scoreInfo = ScoreInfo()
            scoreInfo.numBins = self.numBins
            scoreInfo.property = m.group('property').strip()
            scoreInfo.volume = int(m.group('volume'))
            scoreInfo.pValue = float(m.group('pvalue'))
            scoreInfo.binSize = float(m.group('binsize'))
            scoreInfo.lowValue = float(m.group('lowvalue'))
            scoreInfo.scores = [float(v) for v in m.group('scores').split()]
            self.scores.append(scoreInfo)
            
            

    def convert(self,statFile,numBins=NUM_BINS,pSite=SITE_PRIOR,pNonSite=NONSITE_PRIOR):
        self.numBins = numBins
        self.pSite = pSite
        self.pNonSite = pNonSite
        self.scores = []

        # For each record (prop-vol)
        for recordList in statFile.recordsOrdered:
            for record in recordList:
                # Determine bins for the values
                scoreInfo = ScoreInfo(record,numBins)
                self.scores.append(scoreInfo)

                # Determine frequencies for the bins
                # p(Bin i|Site)
                fBinSite = [0]*numBins
                tBinSite = 0
                for val in record.sites:
                    fBinSite[scoreInfo.getBin(val)] += 1
                    tBinSite += 1
                # p(Bin i|NonSite)
                fBinNonSite = [0]*numBins
                tBinNonSite = 0
                for val in record.nonsites:
                    fBinNonSite[scoreInfo.getBin(val)] += 1
                    tBinNonSite += 1

                # For each bin
                for bin in range(numBins):
                    # Calculate probabilty of bin from frequency counts
                    pBinSite = fBinSite[bin]/float(tBinSite)
                    pBinNonSite = fBinNonSite[bin]/float(tBinNonSite)
                    # p(Site|Bin i)
                    if pBinSite == 0 and pBinNonSite == 0:
                        pSiteBin = pSite
                    else:
                        pSiteBin = pBinSite*pSite/(pBinSite*pSite+pBinNonSite*pNonSite)
                    # Score = log( p(Site|Bin i) / p(Site) )
                    if pSiteBin == 0:
                        score = MIN_SCORE
                    else:
                        score = log(pSiteBin/pSite)
                        if score < MIN_SCORE:
                            score = MIN_SCORE
                    scoreInfo.scores[bin] = score

        return self

    def write(self,outfh=sys.stdout):
        # write header
        outfh.write("""(
(:NUM-BINS %(numBins)d)
(:PROB-SITE %(pSite)g)
(:PROB-NONSITE %(pNonSite)g)
""" % vars(self))
        # write scores
        for scoreInfo in self.scores:
            # write score header
            outfh.write("""((:PROPERTY %(property)s) (:VOLUME %(volume)d) (:P-VALUE %(pValue)g) (:BIN-SIZE %(binSize)g) (:LOW-VALUE %(lowValue)g) (:SCORES (""" % vars(scoreInfo))
            outfh.write(" ".join(["%g" % score for score in scoreInfo.scores]))
            # write score footer
            outfh.write(""")))
""")
        # write footer
        outfh.write(""")
""")

    def writeDelimited(self,outfh=sys.stdout,sep='\t',**options):
        for opt,arg in options.items():
            if opt == 'header' and arg:
                outfh.write('PROPERTY-VOLUME')
                outfh.write(sep)
                outfh.write(sep.join(['LOW_VALUE','BIN_SIZE','P_VALUE']))
                outfh.write(sep)
                outfh.write(sep.join(["BIN_%d" % i for i in range(self.numBins)]))
                outfh.write('\n')

        fields = {}
        for scoreInfo in self.scores:
            sio = StringIO()
            sio.write('%(property)s-%(volume)s' % vars(scoreInfo))
            sio.write(sep)
            sio.write(sep.join(["%g" % v for v in [scoreInfo.pValue,scoreInfo.lowValue,scoreInfo.binSize]]))
            sio.write(sep)
            sio.write(sep.join(["%g" % score for score in scoreInfo.scores]))
            sio.write('\n')
            fields.setdefault(scoreInfo.volume,[]).append(sio.getvalue())
        volumes = fields.keys()
        volumes.sort()
        for vol in volumes:
            for line in fields[vol]:
                outfh.write(line)


if __name__=="__main__":
    import sys,os
    from feature.statfile import StatFile
    def eusage():
        print "Usage: %s FILENAME" % os.path.basename(sys.argv[0])
        sys.exit(1)

    args = sys.argv[1:]
    if len(args) != 1:
        eusage()
    filename = args[0]

    statFile = StatFile(filename)
    scoreFile = ScoreFile().convert(statFile)
    scoreFile.write(sys.stdout)
