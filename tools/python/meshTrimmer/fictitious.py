#! /usr/bin/pvpython
#------------------------------------------------------------------------------#
# Load system modules
import sys
import os
import paraview
from paraview import servermanager as sm
from math import *
#------------------------------------------------------------------------------#
# load the additional modules
fsiroot = os.environ.get( "OPENFTLROOT" )
modules = fsiroot + '/tools/python/modules'
sys.path.append(modules)
from XMLReader import *
from XMLWriter import *
from utils import *
import numpy

#------------------------------------------------------------------------------#
#Get the inputfile
InputfileName = sys.argv[1]
OutputfileName = 'trimmed_' + InputfileName

#------------------------------------------------------------------------------#
# Access data from input file
data = Input(InputfileName)

#------------------------------------------------------------------------------#
# Read data from input 
# - geometry
points,nPoints = GetPoints(data)
cells,nCells   = GetCells(data)
# - data
pointdata = GetPointData(data)
celldata  = GetCellData(data)

#------------------------------------------------------------------------------#
# Assume that the Distance function is named 'Distance'
distance_name = 'Distance'
distance = GetPointDataArray(data,pointdata,distance_name)

#------------------------------------------------------------------------------#
# Identify the different types of cells
cutCells, inActiveCells, activeCells = CellIdentifier(cells,distance)

#------------------------------------------------------------------------------#
# Copy the active cells to a new array
offset = 0
outputCells = []
for cell in activeCells:
	offset = offset + len(cell[0])
	outputCells.append( [cell[0], offset, cell[2] ] )	

# Array of raw data of point data
pointDataArray = []
for array in pointdata :
        pointDataArray.append(GetPointDataArray(data,pointdata,array[0]))

# Array of raw data of cell data
cellDataArray = []
for array in celldata:
	cellDataArray.append(GetCellDataArray(data,celldata,array[0]))


#------------------------------------------------------------------------------#
offset = outputCells[-1][1]
newPoints    = []
newCells     = []
# go through all cut cells
for cell in cutCells:
        # order cell vertices (quadratic elements)
        oCell = orderCell(cell)
        # new cell vertices
        new_vtx = []
        # number of vertices
        nv = len(oCell[0])
        # go through all edges
        for e in range(0,nv):
                # first vertex with distance
                v1 = oCell[0][e]
                phi1 = distance[v1][0]
                # second vertex with distance
                v2 = oCell[0][0]
                if e+1 < nv: v2 = oCell[0][e+1]
                phi2 = distance[v2][0]
                # store positive vertices
                if phi1 >= 0: new_vtx.append(v1)
                # check for cut
                if (phi1 * phi2) <0:
                        # dimensionless paramter
                        s = - phi1/(phi2 - phi1)
                        # geometry of edge
                        X1 = numpy.array(points[v1])
                        X2 = numpy.array(points[v2])
                        dX = X2 - X1
                        # compute new point
                        Xnew = X1 + s * dX
                        # store new node
                        points.append(Xnew)
                        new_vtx.append( nPoints)
                        nPoints = nPoints +1
                        # compute the point data of the new point
                        for pointDatum in pointDataArray:
                                # vertex values
                                pd1 = numpy.array(pointDatum[v1])
                                pd2 = numpy.array(pointDatum[v2])
                                pdnew = pd1 + s * (pd2-pd1)
                                # append to array
                                pointDatum.append( pdnew )
                        # store the new cell
        offset = offset + len( new_vtx )
        outputCells.append( [new_vtx,offset,7] )

#------------------------------------------------------------------------------#
# Gather the indices of all cells for the output
activeCellIds = [cells.index(cell) for cell in activeCells ]
cutCellIds    = [cells.index(cell) for cell in cutCells ]
outputCellIds = activeCellIds + cutCellIds

#------------------------------------------------------------------------------#
# Sample the cell data
outputCellData = []
for array in cellDataArray:
        outputCellData.append( [array[cellId] for cellId in outputCellIds] )

#------------------------------------------------------------------------------#
# Write output fill
outFile = Output(OutputfileName)

WriteHeader( outFile, len(points), len(outputCells) )
WritePoints( outFile, points )
WriteCells(  outFile, outputCells )

WritePointDataHeader(outFile)
for array,name in zip(pointDataArray,pointdata):
	WritePointData(outFile,array,name[0])
WritePointDataFooter(outFile)

WriteCellDataHeader(outFile)
for array,name in zip(outputCellData,celldata):
	WriteCellData(outFile,array,name[0])
WriteCellDataFooter(outFile)

WriteFooter(outFile)  

outFile.close()
#------------------------------------------------------------------------------#
