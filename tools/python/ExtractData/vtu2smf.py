###! /usr/bin/pvpython
#------------------------------------------------------------------#
#------------------------------------------------------------------#
# vtu2smf: generates an smf file given a vtu file
#
# Usage:  ./vtu2smf.py input.vtu
# 
# Important: The pvpython interpreter has to be set in the first
#            line of this script. You can use your local paraview
#            installation or a system installation for this purpose.
#------------------------------------------------------------------#
#------------------------------------------------------------------#

#----------------------------------------------------------------#
import os
import sys
from math import *
# get the path of the extra modules needed
openFTLroot = os.environ.get( "OPENFTLROOT" )
modules = openFTLroot + '/tools/python/modules'
sys.path.append(modules)
# import paraview stuff
import paraview
from paraview import servermanager as sm
from XMLReader import *

#----------------------------------------------------------------#
# input file is first argument, output is constructed
Inputfile  = sys.argv[1]
Outputfile = Inputfile[:-4] + '.smf'
# connect to server manager
connection = sm.Connect()
# read unstructured grid
reader = sm.sources.XMLUnstructuredGridReader()
reader.FileName = Inputfile
# get points and cells
data = sm.Fetch( reader )
Points,nPoints = GetPoints(data)
Cells,nCells = GetCells(data)
#----------------------------------------------------------------#
# output
f = open( Outputfile, 'w' )
f.write( '%i' % nPoints )
f.write( '  ' )
f.write( '%i' % nCells )
f.write( '\n' )
# write points
for i in range(nPoints):
        for j in range(len(Points[i])):
                f.write('%f' % Points[i][j])
                f.write('  ')
        f.write( '\n' )
# write cells
for i in range(nCells):
        for j in range(len(Cells[i][0])):
                f.write('%i' % Cells[i][0][j])
                f.write('  ')
        f.write( '\n' )

f.close()	
