###! /usr/bin/pvpython
#------------------------------------------------------------------#
#------------------------------------------------------------------#
# ExtractCoordinatesAndData: extracts coordinates and specified
#                            data fields in adjacent columns
#
# Usage:  ./ExtractCoordinatesAndData input.vtu
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
Outputfile = Inputfile[:-4] + '.data'

# connect to servermanager
connection=sm.Connect()

# read unstructured grid
reader = sm.sources.XMLUnstructuredGridReader()
reader.FileName = Inputfile

# get data
data = sm.Fetch( reader )
PointData = GetPointData(data) 

# get points
Points,nPoints = GetPoints(data)

#----------------------------------------------------------------#
# Extract specific array
print 'Which array do you want to extract ?'
for info in PointData: 
        print  PointData.index(info), info[0]
        
ArrayId   = int (raw_input('Enter number of array: '))
ArrayName = PointData[ArrayId][0]
Array     = GetPointDataArray( data, PointData, ArrayName )

#----------------------------------------------------------------#
# write point coordinates and array values next to each other
f=open(Outputfile,'w')
##f.write( repr( nPoints ) )
##f.write( '\n' )
for i in range( nPoints ):
        for j in range( len( Points[i] ) ):
            f.write( repr( Points[i][j] ) )
            f.write('  ')
        for j in range( len( Array[i] ) ):
            f.write( repr( Array[i][j] ) )
            f.write('  ')

        f.write( '\n' )

f.close()	
