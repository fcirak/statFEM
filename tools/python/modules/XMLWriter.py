def Output(OUTPUTFILE):
	f_out =  open (OUTPUTFILE,'w')
	return f_out
#-------------------------------------------------------------------------------------------#
def WriteHeader(f,NumberOfPoints=0,NumberOfcells=0):
	f.write('\
<?xml version="1.0"?>\n\
<VTKFile type="UnstructuredGrid">\n\
\t<UnstructuredGrid>\n\
\t<Piece NumberOfPoints="'+repr(NumberOfPoints)+'"  NumberOfCells="'+repr(NumberOfcells)+'">\n\
')
#-------------------------------------------------------------------------------------------#
def WritePoints(f,PointsCoordinates):
	f.write('\
\t\t<Points>\n\
\t\t\t<DataArray type="Float64" NumberOfComponents="3" format="ascii">\n\
')
	for i in range(0,len(PointsCoordinates)):
		for j in range(0,3):
			f.write(repr(PointsCoordinates[i][j])+' ')
		f.write('\n')
	f.write('\
\t\t\t</DataArray>\n\
\t\t</Points>\n\
')
#-------------------------------------------------------------------------------------------#
def WriteCells(f,cells):
	f.write('\
\t\t<Cells>\n\
\t\t\t<DataArray type="Int32" Name="connectivity" format="ascii">\n\
')
#The connectivity
	for i in range(len(cells)):
		for j in range(len(cells[i][0])):
			f.write(repr(cells[i][0][j])+' ')
		f.write('\n')
	f.write('\
\t\t\t</DataArray>\n\
')  
	f.write('\
\t\t\t<DataArray type="Int32" Name="offsets" format="ascii">\n\
')
# The offsets
	for i in range(0,len(cells)):
			f.write(repr(cells[i][1])+' ')
			f.write('\n')
	f.write('\
\t\t\t</DataArray>\n\
')    
	f.write('\
\t\t\t<DataArray type="UInt8" Name="types" format="ascii">\n\
')
#The types 
	for i in range(0,len(cells)):
			f.write(repr(cells[i][2])+' ')
			f.write('\n')
	f.write('\
\t\t\t</DataArray>\n\
\t\t</Cells>\n')
#-------------------------------------------------------------------------------------------#
def WritePointDataHeader(f):
	f.write('\
\t\t<PointData>\n\
')
def WritePointDataFooter(f):
	f.write('\
\t\t</PointData>\n\
')
def WriteCellDataHeader(f):
	f.write('\
\t\t<CellData>\n\
')
def WriteCellDataFooter(f):
	f.write('\
\t\t</CellData>\n\
')
#-------------------------------------------------------------------------------------------#
def WritePointData(f,PointsDatas,Name):
	NumberOfComponents = len(PointsDatas[0])
	f.write('\
\t\t\t<DataArray type="Float64" NumberOfComponents="'+repr(NumberOfComponents)+'" Name="'+Name+'" format="ascii">\n\
')
	for i in range(0,len(PointsDatas)):
		for j in range(0,NumberOfComponents):
			f.write(repr(PointsDatas[i][j])+' ')
		f.write('\n')
	f.write('\
\t\t\t</DataArray>\n\
')
#-------------------------------------------------------------------------------------------#
def WriteCellData(f,CellDatas,Name):
	NumberOfComponents = len(CellDatas[0])
	f.write('\
\t\t\t<DataArray type="Float64" NumberOfComponents="'+repr(NumberOfComponents)+'" Name="'+Name+'" format="ascii">\n\
')
	for i in range(0,len(CellDatas)):
		for j in range(0,len(CellDatas[0])):
			f.write(repr(CellDatas[i][j])+' ')
		f.write('\n')
	f.write('\
\t\t\t</DataArray>\n\
')
#-------------------------------------------------------------------------------------------#
def WriteFooter(f):
	f.write('\
\t</Piece>\n\
\t</UnstructuredGrid>\n\
</VTKFile>\n\
')
