
# $Id: viewhitfiledir.py,v 1.1.1.1 2004/05/22 01:18:15 mliang Exp $
# Copyright 2001 (c) Mike Liang.  All rights reserved.

# Visualize FEATURE .hits directory in Chimera

# =================================================================
# IMPORTS
# =================================================================
import os
import glob
import Tkinter
import Pmw
import tkFileDialog
import viewhitfile



# =================================================================
# CLASSES
# =================================================================
class HitDirDialog(Tkinter.Toplevel):
    # =========================================================
    def __init__(self, dirname, master=None):
        Tkinter.Toplevel.__init__(self, master)
        
        self.dirname = dirname
        self.models = []
        self.prevWin = None

        # Load hitfilenames
        self.hitfilenames = {}
        for hitfilename in glob.glob(os.path.join(self.dirname, "*.hits")):
            (pdbid, ext) = os.path.splitext(os.path.basename(hitfilename))
            self.hitfilenames[pdbid] = hitfilename

        self._buildGui()

    # =========================================================
    def _buildGui(self):
        # Setup window
        top = self
        top.title('FEATURE Hits Visualizer - %s' % self.dirname)
        top.protocol('WM_DELETE_WINDOW', self.destroy)

        # List of PDBs
        frame = Tkinter.Frame(top)
        items = self.hitfilenames.keys()
        items.sort()
        self.listbox = Pmw.ScrolledListBox(frame,
                items=items,
                labelpos='nw',
                label_text='PDB ids',
                selectioncommand = self.pdbSelect,
                usehullsize = 1,
                hull_width = 200,
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
        hfile = viewhitfile.HitFile(self.hitfilenames[pdbid])
        self.prevWin = viewhitfile.HitFileDialog(hfile, master=self)
        
    # =========================================================
    def destroy(self):
        Tkinter.Toplevel.destroy(self)



# =================================================================
# FUNCTIONS
# =================================================================
def getHitDirname():
    d = tkFileDialog.Open(master=None, title="FEATURE Scan Directory")
    dirname, file = os.path.split(d.show())
    return dirname

    
# =================================================================
def loadDir(hitdirname=None):
    if hitdirname == None:
        hitdirname = getHitDirname()

    print "Loading directory",hitdirname
    hitdir = HitDirDialog(hitdirname, master=None)
