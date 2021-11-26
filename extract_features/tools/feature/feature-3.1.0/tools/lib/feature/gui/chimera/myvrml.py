
# $Id: myvrml.py,v 1.1.1.1 2004/05/22 01:18:15 mliang Exp $
# Copyright 2001 (c) Mike Liang.  All rights reserved.

import chimera
import CGLutil.vrml
from StringIO import StringIO
import Midas
import math
import types
import re
import os


def getMolecules():
    return chimera.openModels.list(modelTypes=[chimera.Molecule])

def getResidues(mol):
    return mol.residues

def getResidueName(res):
    return res.type

def selectChain(residues,chain):
    return filter(lambda x,chain=chain:x.id.chainId==chain,residues)

def selectResidues(residues,resName):
    return filter(lambda x,resName=resName:x.type==resName,residues)

def getAtoms(mol):
    return mol.atoms

def getLocation(specifier):
    [atom] = Midas._selectedAtoms(specifier)
    return atom.coord().xyz.data()

def lookupColor(color):
    if type(color) != types.StringType:
        return tuple(color)
    return chimera.Color_lookup(color).rgba()[:3]

def newNode():
    return CGLutil.vrml.Transform()

def newTranslation(x,y,z):
    return CGLutil.vrml.Transform(translation=(x,y,z))

def newRotation(x,y,z,angle):
    return CGLutil.vrml.Transform(rotation=(x,y,z,angle))

def newSphere(radius,color):
    return CGLutil.vrml.Sphere(radius=radius,color=lookupColor(color))

def newBox(length,width,thickness,color):
    return CGLutil.vrml.Box(size=(length, width, thickness), color=lookupColor(color))

def newCylinder(radius,height,color):
    return CGLutil.vrml.Cylinder(radius=radius, height=height, color=lookupColor(color))

def newCone(bottomRadius,height,color):
    return CGLutil.vrml.Cone(bottomRadius=bottomRadius, height=height, color=lookupColor(color))

def setTransparency(node, transparency):
    if transparency != None:
        setattr(node,'transparency',transparency)

def setColor(node, color):
    setattr(node,'color',lookupColor(color))

def useVRML(wrl, name=None, **kwargs):
    if type(wrl) == types.StringType:
        vrmlStr = wrl
    else:
        vrmlStr = CGLutil.vrml.vrml(wrl)

    kwargs['type'] = 'VRML'

    models = apply(chimera.openModels.open,(vrmlStr,), kwargs)

    if name:
        map(lambda x,name=name:setattr(x,'name',name), models)

    return models


class FontStyle(CGLutil.vrml._Node):
    def writeNode(self, f, prefix=''):
        Indent=CGLutil.vrml.Indent
        f.write('%sFontStyle {\n' % prefix)
        p = prefix + Indent
        self.writeAttribute(f, p, 'family', 'family ["%s"]')
        self.writeBooleanAttribute(f, p, 'horizontal', 'horizontal %s')
        self.writeAttribute(f, p, 'justify', 'justify ["%s"]')
        self.writeAttribute(f, p, 'language', 'language %s')
        self.writeBooleanAttribute(f, p, 'leftToRight', 'leftToRight %s')
        self.writeAttribute(f, p, 'size', 'size %g')
        self.writeAttribute(f, p, 'spacing', 'spacing %g')
        self.writeAttribute(f, p, 'style', 'style %s')
        self.writeBooleanAttribute(f, p, 'topToBottom', 'topToBottom %s')
        f.write('%s}\n' % prefix)

    def __str__(self):
        out = StringIO()
        self.writeNode(out)
        return out.getvalue()

class Text(CGLutil.vrml._Node):
    def writeNode(self, f, prefix=''):
        Indent=CGLutil.vrml.Indent
        f.write('%sShape {\n' % prefix)
        self.writeAppearance(f, prefix)
        f.write('%s%sgeometry Text {\n' % (prefix, Indent))
        p = prefix + Indent + Indent
        self.writeAttribute(f, p, 'string', 'string ["%s"]')
        self.writeAttribute(f, p, 'fontStyle', 'fontStyle %s')
        self.writeAttribute(f, p, 'length', 'length [%g]')
        self.writeAttribute(f, p, 'maxExtent', 'maxExtent %g')
        f.write('%s%s}\n' % (prefix, Indent))
        f.write('%s}\n' % prefix)


class GroupNode(CGLutil.vrml._Node):
    def writeNode(self, f, prefix=''):
        f.write('%sGroup {\n' % prefix)
        p = prefix + CGLutil.vrml.Indent
        self.writeChildren(f,p)
        f.write('%s}\n' % prefix)


class FeatureShellNode(GroupNode):
    def __init__(self, numShells, shellThickness,
                 sphereColor=(0,0,.5), sphereTransparency=.5,
                 circleColor="yellow"):
        GroupNode.__init__(self)
        self.numShells = numShells
        self.shellThickness = shellThickness

        # sphere settings
        self.sphereColor = sphereColor
        self.sphereTransparency = sphereTransparency
        self.spheres = self._createSpheres()
        
        # circle settings
        self.circleColor = circleColor
        self.circleNumSegments = 40
        self.circles = self._createCircles()

        # add to node        
        self.addChild(self.spheres)
        self.addChild(self.circles)

    def _createSpheres(self):
        g = GroupNode()
        for shell in range(self.numShells):
            s = newSphere((1+shell)*self.shellThickness, self.sphereColor)
            setTransparency(s, self.sphereTransparency)
            g.addChild(s)
        return g

    def _createCircles(self):
        g = GroupNode()
        for shell in range(self.numShells):
            c = PolygonNode((1+shell)*self.shellThickness, self.circleNumSegments, self.circleColor)
            g.addChild(c)
        return g


class PolygonNode(GroupNode):
    def __init__(self, radius, numSides, color, thickness = 0.05):
        if numSides < 3: raise ValueError, numSides
        
        GroupNode.__init__(self)

        # save circle settings
        self.radius = radius
        self.color = color
        self.thickness = thickness
        self.numSides = numSides

        # Create circle        
        numSegments = numSides
        angle = 2.0*math.pi/numSides
        length = 2*radius * math.tan(angle/2.0)
        for count in range(numSides):
            r=newRotation(0, 0, 1, angle*count)
            t=newTranslation(radius, 0, 0)
            c=newCylinder(thickness, length, color)
            self.addChild(r)
            r.addChild(t)
            t.addChild(c)
        

class MarkerNode(GroupNode):
    def __init__(self, size, color, thickness = 0.05):
        GroupNode.__init__(self)

        # save Marker settings
        self.size = size
        self.color = color
        self.thickness = thickness
        
        # Create Marker 
        angle = math.pi/2.0
        for axis in [(0, 0, 1), (0, 1, 0), (1, 0, 0)]:
            r = apply(newRotation, axis, {'angle':angle})
            c = newCylinder(thickness, size, color)
            self.addChild(r)
            r.addChild(c)

def newMarker(location, size, color, thickness = 0.1):
    t = apply(newTranslation, location)
    m = MarkerNode(size, color, thickness)
    t.addChild(m)
    return t

def getModelNo(modelName):
    for m in chimera.openModels.list(modelTypes=[chimera.Molecule]):
        if m.name.find(modelName) != -1:
            return m.id
    return None

def openMolecule(pdbid, filename=None):
    import tempfile

    if not filename:
        import pdbutils
        filename = pdbutils.pdbfilename(pdbid)

    # Handle gzipped files
    tmpfile = None
    (root,ext) = os.path.splitext(filename)
    if ext in ('.gz',):
        import sys
        sys.path.insert(0,r'C:\Python22\dlls')
        import gzip
        sys.path = sys.path[1:]
        (base,suffix) = os.path.splitext(root)
        tmpfile = tempfile.TemporaryFile(suffix=suffix)
        tmpfile.write(gzip.open(filename).read())
        filename = tmpfile.name

    m = chimera.openModels.open(filename)

    # Close any temporary files
    if tmpfile:
        tmpfile.close()
    # Close all but the first model
    if len(m) > 1:
        chimera.openModels.close(m[1:])
    m = m[0]
    m.name = pdbid
    return m
