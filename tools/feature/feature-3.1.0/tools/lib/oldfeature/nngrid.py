
# $Id: nngrid.py,v 1.4 2004/02/22 21:24:04 mliang Exp $
#
# Data structure to store 3D point data and perform quick queries on
# points within radius of query
#

import math


def distance(pointA, pointB):
    retval = 0
    for idx in range(len(pointA)):
        diff = (pointA[idx]-pointB[idx])
        retval += diff*diff
    return math.sqrt(retval)

def addValue(point, val):
    retval = range(len(point))
    for idx in range(len(point)):
        retval[idx] = point[idx] + val
    return retval

def mulValue(point, val):
    retval = range(len(point))
    for idx in range(len(point)):
        retval[idx] = point[idx] * val
    return retval

def addPoint(pointA, pointB):
    retval = range(len(pointA))
    for idx in range(len(pointA)):
        retval[idx] = pointA[idx]+pointB[idx]
    return retval


class Grid:
    def __init__(self, gridSize):
        self.gridSize = float(gridSize)
        self.table = {}

    def add(self, point, data):
        center=self.gridPoint(point)
        cell=self.table.setdefault(center,[])
        cell.append((point,data))

    def remove(self,point,data):
        cell=self.table.get(point,[])
        try:
            cell.remove(data)
        except ValueError:
            return False
        return True

    def query(self, point, radius, dataOnly=0):
        # returns list of (point, data) pairs
        retlist = []

        # define bounding box
        upper=self.gridPoint(addValue(point,radius))
        lower=self.gridPoint(addValue(point,-radius))

        # traverse grids in box
        for x in range(lower[0], upper[0]+1):
            for y in range(lower[1], upper[1]+1):
                for z in range(lower[2], upper[2]+1):
                    gridpoint=(x,y,z)

                    # add points within radius
                    retlist.extend(self.checkCell(gridpoint, point, radius))

        if dataOnly:
            return [entry[1] for entry in retlist]
        return retlist

    def checkCell(self, gridpoint, point, radius):
        retlist=[]
        cell=self.table.get(gridpoint,[])
        for data in cell:
            if distance(data[0], point) < radius:
                retlist.append(data)
        return retlist

    def gridPoint(self, point):
        gridval = range(len(point))
        for idx in range(len(point)):
            gridval[idx] = self.gridValue(point[idx])
        return tuple(gridval)

    def gridValue(self, val):
        return int(math.floor(val / self.gridSize))


    def quadPoint(self, gridpoint, location):
        # retval: gridPoint*gridSize + (location+1)*0.5*gridSize
        lower = mulValue(gridpoint, self.gridSize)
        offset = addValue(location, 1)
        offset = mulValue(offset, 0.5*self.gridSize)
        return addPoint(lower, offset)

