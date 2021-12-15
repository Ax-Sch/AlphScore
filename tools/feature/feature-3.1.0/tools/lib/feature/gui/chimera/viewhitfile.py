
# $Id: viewhitfile.py,v 1.1.1.1 2004/05/22 01:18:15 mliang Exp $
# Copyright 2001 (c) Mike Liang.  All rights reserved.

# Visualize FEATURE .hits file in Chimera

# =================================================================
# IMPORTS
# =================================================================
import os
import bisect
import Tkinter
import chimera
import CGLutil.vrml
import Midas
import tkFileDialog
import myvrml
#import fetchproteinfile
from hitfile import HitFile, Hit



# =================================================================
# CLASSES
# =================================================================
class HitFileDialog(Tkinter.Toplevel):
    INTERVAL = 0.5
    
    # =========================================================
    def __init__(self, hitfile, master=None):
        Tkinter.Toplevel.__init__(self, master)

        self.hitfile = hitfile
        self.hitNodeList = HitNodeList(self.hitfile.hits)
        self.maxScore = 0
        self.models = []

        basename = os.path.basename(hitfile.filename)
        (pdbid, ext) = os.path.splitext(basename)
        self.pdbid = pdbid


        # Find MaxScore for slider
        if self.hitfile.hits:
            self.maxScore = max(self.hitfile.hits).score

        # Open Molecule if it isn't already open
        self.mol = None
        if myvrml.getModelNo(pdbid.lower()) == None:
            self.mol = myvrml.openMolecule(pdbid)
            # highlight hetatms
            Midas.represent("sphere", "#%d:.het" % self.mol.id)
            Midas.represent("sphere", "#%d:ca" % self.mol.id)

        self._buildGui()

    # =========================================================
    def _buildGui(self):
        # Setup window
        top = self
        top.title('FEATURE hits - %s' % self.pdbid)
        top.protocol('WM_DELETE_WINDOW', self.destroy)

        label = Tkinter.Label(top, text="%s Score cutoff" % self.pdbid)
        label.grid(row=0,column=0)
        
        self.scale = Tkinter.Scale(top, orient=Tkinter.HORIZONTAL, command=self.setCutoff, to=self.maxScore, resolution=HitFileDialog.INTERVAL)
        self.scale.grid(row=1,column=0)
        self.scale.set(self.maxScore)
        
        label = Tkinter.Label(top, text="Number Hits:")
        label.grid(row=2,column=0)
        
        self.numHiHits = Tkinter.StringVar(top)
        label = Tkinter.Entry(top, textvariable=self.numHiHits)
        label.grid(row=3,column=0)

    # =========================================================
    def setCutoff(self, cutoff):
        cutoff = float(cutoff)
        self.closeModels()
        numHits = self.hitNodeList.setCutoff2(cutoff)
        self.numHiHits.set(numHits)
        self.addModels(myvrml.useVRML(self.hitNodeList, "Sites - %s" % self.pdbid))

    # =========================================================
    def destroy(self):
        self.closeModels()
        if self.mol:
            if self.mol in chimera.openModels.list():
                chimera.openModels.close(self.mol)
        Tkinter.Toplevel.destroy(self)

    # =========================================================
    def closeModels(self):
        Midas.freeze()
        for model in self.models:
            if model in chimera.openModels.list():
                chimera.openModels.close([model])
        self.models = []
        
    # =========================================================
    def addModels(self, modelList):
        self.models.extend(modelList)


# =================================================================
class HitNode(CGLutil.vrml.Transform):
    # =========================================================
    def __init__(self, hit, radius=0.1, color="red", transparency=None):
        CGLutil.vrml.Transform.__init__(self, translation = hit.getLocation())
        self.hit = hit
        sphere = myvrml.newSphere(radius, color)
        if transparency:
            myvrml.setTransparency(sphere, transparency)
        self.addChild(sphere)

    # =========================================================
    def __cmp__(self, other):
        if isinstance(other, HitNode):
            return cmp(self.hit, other.hit)
        return cmp(self.hit, other)

    # =========================================================
    def getScore(self):
        return self.hit.getScore()

    # =========================================================
    def getLocation(self):
        return self.hit.getLocation()


# =================================================================
class HitNodeList(myvrml.GroupNode):
    # =========================================================
    def __init__(self, hits, cutoff=None, minR=0.1, maxR=0.5):
        myvrml.GroupNode.__init__(self)

        self.hits = hits
        self.hits.sort()
        self.hitNodes = map(HitNode, self.hits)
        self.cutoff = cutoff
        self.children = []

        self.propRange = {
            'radius': (minR, maxR),
            'colorR': (0, 1),
            'colorG': (0, 0),
            'colorB': (0.5, 0),
            'transparency': (1, 0)
        }

        self.setCutoff(self.cutoff)

    # =========================================================
    def setCutoff(self, cutoff):
        self.cutoff = cutoff
        if cutoff == None:
            self.children = self.hitNodes
        else:
            idx = bisect.bisect(self.hitNodes, cutoff)
            self.children = self.hitNodes[idx:]
        return len(self.children)

    # =========================================================
    def setCutoff2(self, cutoff):
        self.cutoff = cutoff
        return self.updateNodes()

    # =========================================================
    def updateNodes(self):
        # Initialize property scales
        propScale = {}
        if self.cutoff:
            for prop in self.propRange.keys():
                propScale[prop] = float(self.propRange[prop][1] - self.propRange[prop][0])/self.cutoff
        else:
            for prop in self.propRange.keys():
                propScale[prop] = 0
            
        hitNodes = []
        for hit in self.hits:
            # Calculate property values
            value = min(self.cutoff, hit.score)
            propValue = {}
            for prop in self.propRange.keys():
                propValue[prop] = self.propRange[prop][0] + propScale[prop] * value

            # Set values
            radius = propValue['radius']
            color = (propValue['colorR'], propValue['colorG'], propValue['colorB'])
            transparency = propValue['transparency']

            # Create node
            node = HitNode(hit, radius=radius, color=color, transparency=transparency)

            # Add to list
            hitNodes.append(node)

        self.children = hitNodes
        return len(self.children)


# =================================================================
# FUNCTIONS
# =================================================================
def getHitFilename():
    filetypes = [("FEATURE results", "*.hits"), ("All Files", "*")]
    d = tkFileDialog.Open(master=None, title="FEATURE Scan Results",
                            filetypes=filetypes)
    return d.show()


# =================================================================
def loadFile(hitfilename=None, master=None):
    if hitfilename == None:
        hitfilename = getHitFilename()

    print "Loading file '%s'" % hitfilename
    hitfile = HitFile(hitfilename)
    return HitFileDialog(hitfile, master)
