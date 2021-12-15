#! /usr/bin/env python
# $Id: modelinfofile.py,v 1.1 2004/11/17 09:29:30 mliang Exp $

# Prepares models for webfeature
import sys,os
import glob
import shutil
import re
from sets import Set
from utils import chomp
from seqfeature.perffile import PerformanceFile
from seqfeature.scoresfile import ScoresEntry

MODEL_INFO_FILENAME = 'model.info'
PERFORMANCE_FILENAME = 'assesscv.40.out'
CUTOFF_FILENAME = 'cutoff_0.99_txt'

DEFAULT_NUM_SHELLS = 6
DEFAULT_SHELL_WIDTH = 1.25

MODEL_EXT = '.model'
XMODEL_EXT = '.xmodel'

def FloatFmt(val):
    if type(val) == float:
        return '%.3f' % val
    return str(val)

class CutoffFile:
    def __init__(self,filename):
        self.cutoff = None
        self.sensitivity = None
        self.specificity = None
        self.falseNegatives = []
        self.Load(filename)

    def Load(self,filename):
        infh = file(filename)
        for line in infh:
            if line[0] == '#':
                fields = chomp(line).split('\t')
                self.cutoff = float(fields[1])
                self.specificity = float(fields[2])
                self.sensitivity = float(fields[3])
            else:
                self.falseNegatives.append(ScoresEntry(line))

    def GetCutoff(self):
        return self.cutoff

    def GetSpecificity(self):
        return self.specificity

    def GetSensitivity(self):
        return self.sensitivity

class WebFeatureModel:
    def __init__(self,modelFile):
        self.modelFile = modelFile
        dirname = os.path.dirname(modelFile)

        self.performanceFilename = os.path.join(dirname,PERFORMANCE_FILENAME)
        self.performanceFile = PerformanceFile(self.performanceFilename)
        self.auc = self.performanceFile.GetAuc()

        self.cutoffFilename = os.path.join(dirname,CUTOFF_FILENAME)
        self.cutoffFile = CutoffFile(self.cutoffFilename)
        self.cutoff = FloatFmt(self.cutoffFile.GetCutoff())
        self.cutoffSens = FloatFmt(self.cutoffFile.GetSensitivity())
        self.cutoffSpec = FloatFmt(self.cutoffFile.GetSpecificity())

        self.modelName = os.path.basename(dirname)
        fields = self.modelName.split('.')
        self.motifName = fields[0]
        self.ccSitePos = int(fields[1])
        self.resName = fields[2]
        self.atom = fields[3]

        self.numShells = DEFAULT_NUM_SHELLS
        self.shellWidth = DEFAULT_SHELL_WIDTH
        
        self.baseModelFilename = self.modelName + MODEL_EXT
        self.baseExtModelFilename = self.modelName + XMODEL_EXT

        self.extModelFile = re.sub(r'\%s$' % MODEL_EXT,XMODEL_EXT,self.modelFile)

    def WriteModel(self,modelpath):
        outputDir = os.path.join(modelpath,self.modelName)
        if not os.path.exists(outputDir):
            os.makedirs(outputDir)

        modelInfoFilename = os.path.join(outputDir,MODEL_INFO_FILENAME)
        modelInfoFile = file(modelInfoFilename,'w')
        # create model info
        print >>modelInfoFile, '''[Model Info]
name: %(modelName)s
description: Built from %(motifName)s ProSite pattern hits on Astral40 centered on %(resName)s@%(atom)s at pattern position %(ccSitePos)s. Performance=%(auc)s (AUC). At %(cutoff)s cutoff, specificity=%(cutoffSpec)s and sensitivity=%(cutoffSens)s.
type: protein-seqfeature
docfile: model.html

[SeqFeature Model]
residue=%(resName)s
atom=%(atom)s
modelfilename=%(baseModelFilename)s
numshells=%(numShells)s
shellwidth=%(shellWidth)s
motifName=%(motifName)s
ccSitePos=%(ccSitePos)s
auc=%(auc)s
cutoff=%(cutoff)s
''' % vars(self)

        # copy model file
        newFilename = os.path.join(outputDir,self.baseModelFilename)
        shutil.copy(self.modelFile,newFilename)
        newFilename = os.path.join(outputDir,self.baseExtModelFilename)
        shutil.copy(self.extModelFile,newFilename)
        shutil.copy(self.performanceFilename,outputDir)
        shutil.copy(self.cutoffFilename,outputDir)

class ModelFilter:
    def __init__(self,filename):
        self.Load(filename)

    def Load(self,filename):
        infh = file(filename)
        self.disallowedEntries = Set()
        for line in infh:
            fields = line.split()
            entry = tuple(fields[0].split('/')[:3])
            self.disallowedEntries.add(entry)

    def Allow(self,filename):
        entry = tuple(os.path.dirname(filename).split('/')[-3:])
        return entry in self.disallowedEntries



inputPath = 'raw_models'
outputPath = 'models'
#filterFilename = 'assesscv.40.maxauc.out'
#modelFilter = ModelFilter(filterFilename)
modelFileList = glob.glob(os.path.join(inputPath,'*/*.model'))
for modelFile in modelFileList:
#    if not modelFilter.Allow(modelFile):
#        continue
    model = WebFeatureModel(modelFile)
    print >>sys.stderr,model.modelName,model.auc
    model.WriteModel(outputPath)
