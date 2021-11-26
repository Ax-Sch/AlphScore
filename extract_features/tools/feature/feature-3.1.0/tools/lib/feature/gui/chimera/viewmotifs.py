#! /usr/bin/env python

# $Id: viewmotifs.py,v 1.2 2004/05/22 21:50:15 mliang Exp $
# Copyright (c) 2001 Mike Liang.  All rights reserved.

# Visualizes motif information for a specific PDB structure

# =================================================================
# IMPORTS
# =================================================================
import os
import string
import chimera
import Midas
import myvrml
import Pmw
import Tkinter
import tkFileDialog
import sys



# =================================================================
# GLOBALS
# =================================================================
gHaveCore = 0
gVerbose = 0



# =================================================================
# CLASSES
# =================================================================
class Motif:
    if gHaveCore:
        maxsplits = 6
    else:
        maxsplits = 4
        
    # =========================================================
    def __init__(self, line):
        self.parseLine(line)

    # =========================================================
    def parseLine(self, line):
        fields = line.rstrip().split(" ",self.maxsplits)

        try:
            self.pdbid = fields.pop(0)
            self.chainid = fields.pop(0)
            self.start = int(fields.pop(0))
            self.end = int(fields.pop(0))
            if gHaveCore:
                self.coreStart = int(fields.pop(0))
                self.coreEnd = int(fields.pop(0))
            self.motif = fields.pop(0)
        except IndexError, e:
            print >>sys.stderr, "Error on line:", line
            raise e

    # =========================================================
    def __str__(self):
        if gHaveCore:
            return string.join(map(str,(self.pdbid, self.chainid, self.start, self.end, self.coreStart, self.coreEnd, self.motif)))
        return string.join(map(str,(self.pdbid, self.chainid, self.start, self.end, self.motif)))

    # =========================================================
    def __cmp__(self, other):
        return cmp(self.tuple(), other.tuple())

    # =========================================================
    def tuple(self):
        return (self.pdbid, self.chainid, self.start, self.end, self.motif)

    # =========================================================
    def prettyDisplay(self):
        return "%4s %1s %3d %3d %s" % self.tuple()


# =================================================================
class MotifFile:
    # =========================================================
    def __init__(self, filename = None):
        self.motifs = {}
        self.filename = filename

        if filename:
            self.loadFile(filename)

    # =========================================================
    def loadFile(self, filename):
        self.filename = filename
        file = open(filename)
        for line in file.readlines():
            motif = Motif(line)
            self.motifs.setdefault(motif.pdbid, []).append(motif)
        file.close()


# =================================================================
class PDBMotifDialog(Tkinter.Toplevel):
    # =========================================================
    def __init__(self, pdbid, motifs, master=None):
        Tkinter.Toplevel.__init__(self, master)

        self.pdbid = pdbid
        self.motifs = motifs
        self.models = []
        self.motifs.sort()

        # Open Molecule if it isn't already open
        self.mol = None
        if myvrml.getModelNo(pdbid.lower()) == None:
            self.mol = myvrml.openMolecule(pdbid)
            # highlight hetatms
            Midas.represent("sphere", "#%d:.het" % self.mol.id)
            Midas.represent("sphere", "#%d:ca" % self.mol.id)

        self._displayGui()

    # =========================================================
    def _displayGui(self):
        # Setup window
        top = self
        top.title('Motif Browser - %s' % self.pdbid)
        top.protocol('WM_DELETE_WINDOW', self.destroy)

        # List of Motifs
        frame = Tkinter.Frame(top)
        items = map(Motif.prettyDisplay, self.motifs)
        self.listbox = Pmw.ScrolledListBox(frame,
                items=items,
                labelpos='nw',
                label_text='Motifs',
                selectioncommand = self.motifSelect,
                usehullsize = 1,
                hull_width = 300,
                hull_height = 200,
        )
        self.listbox._listbox.config(font="courier 8")
        self.listbox._listbox.configure(selectmode='extended')
        self.listbox.pack(fill = 'both', expand = 1, padx = 5, pady = 5)
        frame.pack(fill='both',expand=1)
      
    # =========================================================
    def motifSelect(self):
        self.closeModels()
        self._clearModelColor()
        sel_idxs = map(int, self.listbox.curselection())
        if not sel_idxs:
            return
        node = myvrml.GroupNode()
        for idx in sel_idxs:
            motif = self.motifs[idx]
            m = self._displayMotif(motif)
            node.addChild(m)

        # self.addModels(myvrml.useVRML(node, "%s - Motifs" % motif.pdbid))

    # =========================================================
    def _clearModelColor(self):
        # obtain modelNo
        modelNo = myvrml.getModelNo(self.pdbid.lower())
        if modelNo == None:
            return

        # Reset model color
        Midas.uncolor(None, "#%d" % modelNo)
        Midas.unrescolor(None, "#%d" % modelNo)
        
    # =========================================================
    def _displayMotif(self, motif):
        # obtain modelNo
        modelNo = myvrml.getModelNo(motif.pdbid.lower())
        if modelNo == None:
            print "No Model"
            return

        group = myvrml.GroupNode()

        # obtain chain residue numbers
        chainRes = []
        mol = chimera.openModels.list(id=modelNo)[0]
        residues = mol.residues
        residues.sort(lambda a,b:cmp(a.id.position,b.id.position))
        for res in residues:
            if res.id.chainId.strip() == motif.chainid:
                chainRes.append(res)
        
        # Hilight motif
        motifSel = "#%d:%d-%d.%s" % (modelNo, chainRes[motif.start].id.position,
                                     chainRes[motif.end-1].id.position, motif.chainid)
        if gVerbose:
            print "MOTIFSEL:",motifSel
            print "START RES:",chainRes[motif.start].type,
            print "STOP RES:",chainRes[motif.end-1].type
        Midas.color("red", motifSel)
        Midas.rescolor("red", motifSel)

        if gHaveCore:
            # Hilight coreMotif
            coreMotifSel = "#%d:%d-%d.%s" % (modelNo, chainRes[motif.coreStart].id.position,
                                             chainRes[motif.coreEnd-1].id.position, motif.chainid)
            print "CORESEL:",coreMotifSel
            Midas.color("yellow", coreMotifSel)

        print motif
        
        return group

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
class MotifFileDialog(Tkinter.Toplevel):
    # =========================================================
    def __init__(self, motiffile, master=None):
        Tkinter.Toplevel.__init__(self, master)

        self.motiffile = motiffile
        self.prevWin = None
        
        self._buildGui()

    # =========================================================
    def _buildGui(self):
        # Setup window
        top = self
        top.title('Protein Visualizer')
        top.protocol('WM_DELETE_WINDOW', self.destroy)

        # List of PDBs
        frame = Tkinter.Frame(top)
        items = self.motiffile.motifs.keys()
        items.sort()
        self.listbox = Pmw.ScrolledListBox(frame,
                items=items,
                labelpos='nw',
                label_text='PDB ids',
                selectioncommand = self.pdbSelect,
                usehullsize = 1,
                hull_width = 100,
                hull_height = 200,
        )
        self.listbox.pack(fill = 'both', expand = 1, padx = 5, pady = 5)
        frame.pack(fill='both',expand=1)

    # =========================================================
    def pdbSelect(self):
        sels = self.listbox.getcurselection()
        if not sels:
            return
        if self.prevWin != None and self.prevWin.winfo_exists():
            self.prevWin.destroy()
        pdbid = sels[0]
        self.prevWin = PDBMotifDialog(pdbid, self.motiffile.motifs[pdbid], master=self)
        
    # =========================================================
    def destroy(self):
        Tkinter.Toplevel.destroy(self)



# =================================================================
# FUNCTIONS
# =================================================================
def getMotifFileName():
    filetypes = [("FEATURE motifs", "*.motifs"), ("All Files", "*")]
    d = tkFileDialog.Open(master=None, title="Motif Files", filetypes=filetypes)
    return d.show()

    
# =================================================================
def loadFile(motiffilename=None, master=None, verbose=None):
    global gVerbose
    
    if motiffilename == None:
        motiffilename = getMotifFileName()

    if verbose != None:
        gVerbose = verbose

    print "Loading file '%s'" % motiffilename
    motiffile = MotifFile(motiffilename)
    return MotifFileDialog(motiffile, master)
