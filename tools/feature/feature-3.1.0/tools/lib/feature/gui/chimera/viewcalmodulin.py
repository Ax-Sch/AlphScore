
# $Id: viewcalmodulin.py,v 1.1.1.1 2004/05/22 01:18:15 mliang Exp $
# Copyright 2001 (c) Mike Liang.  All rights reserved.

import os
import chimera
import Midas
from myvrml import *

pdbid = "1CLL"
pdbdir = os.environ.get("PROTEIN_PDB_DIR","C:/Databases/PDB")
pdbfilename = os.path.join(pdbdir, "pdb%s.ent" % pdbid.lower())
hitfilename = "G:/CurrentResearch/PyFindSite/centroid2/Scan-centroid2/Hits/1CLL.hits"



def MethodOne():
    Midas.close("all")
    openMolecule(pdbid)

    efhand = "#0:117-147"
    calcium = "#0:151.het"
    backboneAtoms = "@n,c,ca,o"
    negRes = "::asp,glu"

    # Display EF-Hand
    Midas.represent("bs", efhand)
    Midas.show(efhand + "&" + negRes)
    Midas.color("green", efhand + "&" + negRes)
    Midas.color("cyan", "#0:135")
    Midas.color("byatom", efhand + backboneAtoms)
    # Midas.resrepr("edged", efhand)
    Midas.resrepr("smooth", "#0")

    # Display Calcium
    Midas.display(calcium)
    Midas.color("cyan", calcium)
    Midas.vdwdefine("ca", "ca", 0.5, calcium)
    Midas.represent("sphere", calcium)

    # Adjust View
    Midas.center(calcium)
    Midas.cofr(calcium)
    chimera.viewer.scaleFactor = 7


    # Add Feature Shells
    s = FeatureShellNode(6, 1, sphereColor=(.2,0,.4),sphereTransparency=.1)
    t = apply(newTranslation, getLocation(calcium))
    t.addChild(s.spheres)
    useVRML(t,name="FEATURE Spheres")
    t = apply(newTranslation, getLocation(calcium))
    t.addChild(s.circles)
    useVRML(t,name="FEATURE Rings")


    # Add Markers
    markerSize = 1
    markerColor = "yellow"

    alphaLoc = (26.768,15.445,3.308)
    betaLoc = (24.006,18.672,1.650)
    centroidLoc = (27.1862,17.9809,1.31815)

    g = GroupNode()
    g.addChild(newMarker(alphaLoc, markerSize, markerColor))
    g.addChild(newMarker(betaLoc, markerSize, markerColor))
    g.addChild(newMarker(centroidLoc, markerSize, markerColor))
    #useVRML(g,name="FEATURE Markers")

    # Show Water
    water = "#0:231.water"
    Midas.display(water)
    Midas.vdwdefine("hoh", "o", 0.2, water)
    Midas.represent("sphere", water)
    

    # Add Hits
    #import viewhitfile
    #viewhitfile.loadFile(hitfilename)


def MethodTwo():
    Midas.close("all")

    # Parameters
    featureNumShells = 6
    featureShellThickness = 1.0
    
    # Load Molecule
    modelSel = "#0"
    [mol] = chimera.openModels.open(pdbfilename)
    chimera.viewer.scaleFactor = 4.5

    # Draw as Ribbon
    Midas.resrepr("smooth", modelSel)
    for res in mol.residues:
	res.drawMode = chimera.Residue.Ribbon_Edged

    # Display Calciums
    calciumsSel = modelSel + "::ca"
    Midas.vdwdefine("ca", "ca", 0.5, calciumsSel)
    Midas.show(calciumsSel)
    Midas.represent("sphere", calciumsSel)
    
    # Display Shells
    calciums = Midas._selectedAtoms(calciumsSel)
    featureModels = []
    for calcium in calciums:
        location = calcium.coord().xyz.data()
        t = apply(newTranslation, location)
        f = FeatureShellNode(featureNumShells, featureShellThickness)
        t.addChild(f)
        featureModels += useVRML(t,name="FEATURE Shell")

    # Load hitsfile
    import viewhitfile
    viewhitfile.loadFile(hitfilename)
