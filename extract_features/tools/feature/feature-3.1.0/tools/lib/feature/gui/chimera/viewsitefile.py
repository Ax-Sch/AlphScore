
# $Id: viewsitefile.py,v 1.1.1.1 2004/05/22 01:18:15 mliang Exp $
# Copyright 2001 (c) Mike Liang.  All Rights Reserved.

# Visualize FEATURE .site file in Chimera

# =================================================================
# IMPORTS
# =================================================================
import sys
import os
import Tkinter
import Pmw
import chimera
import tkFileDialog
import myvrml
import Midas
import feature.sitefile
SiteFile = feature.sitefile.SiteFile


def PrettySiteString(site):
    return "%(pdbid)4s %(x)8.3f %(y)8.3f %(z)8.3f %(label)8s" % vars(site)


# =================================================================
# CLASSES
# =================================================================
class PDBSiteDialog(Tkinter.Toplevel):
    # =========================================================
    def __init__(self, pdbid, sites, master=None):
        Tkinter.Toplevel.__init__(self, master)
        

        self._buildGui()

    # =========================================================
    def _buildGui(self):
        # Setup window
        top = self
        top.title('FEATURE sites - %s' % self.pdbid)
        top.protocol('WM_DELETE_WINDOW', self.destroy)

        # List of Sites
        frame = Tkinter.Frame(top)
        frame.pack(fill='both',expand=1)





# =================================================================
class SiteFileDialog(Tkinter.Toplevel):
    # =========================================================
    def __init__(self, sitefile, master=None):
        Tkinter.Toplevel.__init__(self, master)
        
        self.sitefile = sitefile # SiteFile()
        self.pdbid = None        # current pdbid
        self.sites = []          # list of Site() for pdbid
        self.models = []         # list of models opened
        self.mol = None          # current open molecule

        self._buildGui()

    # =========================================================
    def _buildGui(self):
        # Setup window
        top = self
        top.title('FEATURE sites - %s' % os.path.basename(self.sitefile.filename))
        top.protocol('WM_DELETE_WINDOW', self.destroy)

        # List of PDBs
        frame = Tkinter.Frame(top)
        pdbids = self.sitefile.getPdbids()
        self.pdbListBox = Pmw.ScrolledListBox(frame,
                items=pdbids,
                labelpos='nw',
                label_text='PDB ids',
                selectioncommand = self.pdbSelect,
                usehullsize = 1,
                hull_width = 100,
                hull_height = 200,
        )
        self.pdbListBox.grid(row=0,column=0,padx=5,pady=5)

        # List of Sites
        self.siteStrings = []
        self.siteListBox = Pmw.ScrolledListBox(frame,
                items=self.siteStrings,
                labelpos='nw',
                label_text='Sites',
                selectioncommand = self.siteSelect,
                usehullsize = 1,
                hull_width = 300,
                hull_height = 200,
        )
        self.siteListBox._listbox.configure(selectmode='multiple')
        self.siteListBox.grid(row=0,column=1,sticky='e',padx=5,pady=5)
        
        frame.pack(fill='both',expand=1)


    # =========================================================
    def destroy(self):
        self.cleanupPdbSelect()
        Tkinter.Toplevel.destroy(self)

    def cleanupPdbSelect(self):
        # close open Site vrml objects
        self.closeModels()
        # close any molecule we opened
        if self.mol:
            if self.mol in chimera.openModels.list():
                chimera.openModels.close(self.mol)
        self.mol = None

    # =========================================================
    def pdbSelect(self):
        # Identify selected pbid
        sels = self.pdbListBox.getcurselection()
        if not sels:
            return
        pdbid = sels[0]

        # Set pdbid
        self.pdbid = pdbid

        # clean up previous pdbid selection
        self.cleanupPdbSelect()

        # Open Molecule if it isn't already open
        if myvrml.getModelNo(pdbid.lower()) == None:
            self.mol = myvrml.openMolecule(pdbid)
            # highlight hetatms
            Midas.represent("sphere", "#%d:.het" % self.mol.id)
            Midas.represent("sphere", "#%d:ca" % self.mol.id)

        # Get sites for the pdbid
        self.sites = self.sitefile.getByPdbid(pdbid)

        # Update listbox
        siteStrings = [PrettySiteString(s) for s in self.sites]
        self.siteListBox.setlist(siteStrings)
        
        # Select all sites
        self.siteListBox._listbox.select_set(0,'end')
        self.siteSelect()

        # Update siteListBox label_text
        

    # =========================================================
    def siteSelect(self):
        self.closeModels()
        # identify selected sites
        sel_idxs = map(int, self.siteListBox.curselection())
        if not sel_idxs:
            return

        # create vrml node for sites
        node = myvrml.GroupNode()
        for idx in sel_idxs:
            site = self.sites[idx]
            t = myvrml.newTranslation(site.x, site.y, site.z)
            m = myvrml.MarkerNode(4, self.typeColor(site.label), thickness=.5)
            t.addChild(m)
            node.addChild(t)
        self.addModels(myvrml.useVRML(node, "Sites - %s" % self.pdbid))

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

    # =========================================================
    COLORMAP = {
            'site': "green",
            'nonsite': "blue"
        }
    def typeColor(self, typeStr):
        return self.COLORMAP[typeStr]



# =================================================================
# FUNCTIONS
# =================================================================
def getSiteFileName():
    filetypes = [("FEATURE sitefiles", "*.site"), ("All Files", "*")]
    d = tkFileDialog.Open(master=None, title="FEATURE .site Files",
                          filetypes=filetypes)
    return d.show()

    
# =================================================================
def loadFile(sitefilename=None, master=None):
    if sitefilename == None:
        sitefilename = getSiteFileName()

    print "Loading file '%s'" % sitefilename
    sitefile = SiteFile(sitefilename)
    return SiteFileDialog(sitefile, master)
