
# ====================================================================
# IMPORTS
# ====================================================================
import sys, os
import types
import math
import re
import gzip



# ====================================================================
# FUNCTIONS
# ====================================================================
def uniqList(aList):
    retList = [ aList[0] ]
    for item in aList[1:]:
        if item != retList[-1]:
            retList.append(item)
    return retList


# ====================================================================
def sum(aList):
    if not aList:
        return None
    value = reduce(lambda a,b:a+b, aList)
    return value


# ====================================================================
def mean(aList):
    if not aList:
        return None
    total = float(sum(aList))
    return total/len(aList)


# ====================================================================
def stdev(aList, meanVal=None):
    if len(aList) < 2:
        return None
    if meanVal == None:
        meanVal = mean(aList)
    sumsqr = float(sum([(val-meanVal)**2 for val in aList]))
    return math.sqrt(sumsqr/(len(aList)-1))


# ====================================================================
def point(obj):
    return (obj.x,obj.y,obj.z)


# ====================================================================
def efile(filename,*args,**keywords):
    try:
        return file(filename,*args,**keywords)
    except IOError,e:
        print e
        sys.exit(2)


# ====================================================================
def gzopen(filename,*args,**keywords):
    fileext = os.path.splitext(filename)[1]
    if fileext == ".gz":
        return gzip.open(filename,*args,**keywords)
    elif fileext == '.Z':
        return os.popen("gunzip -c %s" % filename)
    else:
        return file(filename,*args,**keywords)


# ====================================================================
def CalculateCentroid(points):
    # get rid of Empty points
    points = filter(lambda x:x,points)

    # if no more points, return
    if not points:
        return None
    
    centroid=[0,0,0]
    count=0
    # sum points
    for point in points:
        centroid[0] += point[0]
        centroid[1] += point[1]
        centroid[2] += point[2]
        count+=1

    # divide by total number of points
    centroid[0] /= 1.0*count
    centroid[1] /= 1.0*count
    centroid[2] /= 1.0*count
    return centroid


# ====================================================================
def GetAtomLocation(residue, atomNames):
    if type(atomNames) == types.StringType:
        atomNames = [atomNames]
    for atomName in atomNames:
        atom = residue.getAtom(atomName)
        if atom:
            return atom.getLocation()
    return None


# ====================================================================
def PatternToMask(pattern):
    def expandRepeat(m):
        r = int(m.group(3))
        return m.group(1)*r
    
    mask = re.sub("\[.+?\]","2",pattern)
    mask = re.sub("\.","0",mask)
    mask = re.sub("(.){\s*(\d+\s*,\s*)*(\d+)\s*}",expandRepeat,mask)
    mask = re.sub("[^0]","1",mask)
    return mask


# ====================================================================
def MaskResidues(residues, mask):
    retlist = []
    for idx in range(len(residues)):
        if mask[idx] != "0":
            retlist.append(residues[idx])
    return retlist


# ====================================================================
def distance(pointA, pointB):
    sum = 0
    for idx in range(len(pointA)):
        sum += (pointA[idx]-pointB[idx])*(pointA[idx]-pointB[idx])
    return math.sqrt(sum)


# ====================================================================
def GetNearestAtom(atoms, point):
    if len(atoms) < 1:
        return None

    nearest = (atoms[0], distance(point, atoms[0].getLocation()))
    
    for atom in atoms:
        dist = distance(point, atom.getLocation())
        if dist < nearest[1]:
            nearest = (atom, dist)

    return nearest


# ====================================================================
def GetNearestResidue(residues, point, atomname):
    atoms = []
    for residue in residues:
        if atomname == "*":
            atoms.extend(residue.getAtoms())
        else:
            atom = residue.getAtom(atomname)
            if atom:
                atoms.append(atom)
    
    nearest = GetNearestAtom(atoms, point)
    if nearest:
        return nearest[0].residue
    return None


# ====================================================================
def ResidueListSequence(residues):
    import PDB
    seq = ""
    for r in residues:
        seq += PDB.residueLetter(r.resName)
    return seq


# ====================================================================
def GetNearbyObjects(list1, list2, radius):
    retlist = []
    for obj1 in list1:
        for obj2 in list2:
            if distance(point(obj1),point(obj2)) <= radius:
                retlist.append(obj2)
    return retlist


