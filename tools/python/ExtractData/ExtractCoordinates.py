###! /usr/bin/pvpython
#------------------------------------------------------------------#
#------------------------------------------------------------------#
# ExtractCoordinates: extracts coordinates and writes the to a file
#
# Usage:  ./ExtractCoordinates input.vtu
# 
# Important: The pvpython interpreter has to be set in the first
#            line of this script. You can use your local paraview
#            installation or a system installation for this purpose.
#------------------------------------------------------------------#
#------------------------------------------------------------------#

#----------------------------------------------------------------#
# system and math import
import os
from math import *
import sys
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
Outputfile = Inputfile[:-4] + '.coords'

# connect to servermanager
connection=sm.Connect()

# read unstructured grid
reader = sm.sources.XMLUnstructuredGridReader()
reader.FileName = Inputfile

# get points
data = sm.Fetch( reader )
PointData = GetPointData(data) 


Points,nPoints = GetPoints(data)

#----------------------------------------------------------------#
# write point coordinates 
f = open( Outputfile, 'w' )
f.write( repr( nPoints ) )
f.write( '\n' )
for i in range(nPoints):
        for j in range( len( Points[i] ) ):
            f.write( repr( Points[i][j] ) )
            f.write('  ')
        f.write( '\n' )

f.close()	
