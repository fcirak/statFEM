def Input(INPUTFILE):
	import paraview
	from paraview import servermanager as sm

	# connect to right files
	connection=sm.Connect()

	# create a reader
	reader = sm.sources.XMLUnstructuredGridReader( FileName=INPUTFILE )
	reader.UpdatePipeline()

	# get the data i.e moves data from the server to the client
	data = sm.Fetch( reader )
	return data
#-------------------------------------------------------------------------------------------#
def GetPoints(data):
	points = data.GetPoints()
	npoints = data.GetNumberOfPoints()
	Points=[]
	for i in range(npoints):
		a,b,c=points.GetPoint(i)
		Points.append([a,b,c])
	return Points,npoints
#-------------------------------------------------------------------------------------------#
def GetCells(data):
    ncells = data.GetNumberOfCells()
    Cells = []
    offsets =0
    for i in range(0,ncells):
        cell=data.GetCell(i)
        types = cell.GetCellType()
        offsets = offsets + cell.GetNumberOfPoints()
        pointid = cell.GetPointIds()
        pointpercell = cell.GetNumberOfPoints()
        connectivity=[]
        for j in range(0,pointpercell):
            connectivity.append(pointid.GetId(j))
        Cells.append([connectivity,offsets,types])
    return Cells,ncells 
#-------------------------------------------------------------------------------------------#
def GetPointData(data):
	pdata = data.GetPointData()
	nparrays = pdata.GetNumberOfArrays()
	PointData=[]
	for i in range(0,nparrays):
		array = pdata.GetArray(i)
		name = array.GetName()
		numberofcomponents = array.GetNumberOfComponents()
		PointData.append([name,numberofcomponents])
	return PointData
#-------------------------------------------------------------------------------------------#
def GetPointDataArray(data,PointData,name):
	pdata = data.GetPointData()
	nparrays = pdata.GetNumberOfArrays()
	npoints = data.GetNumberOfPoints()
	place = 0 #default value
	for i in range(0,nparrays):
		if PointData[i][0] == name:
			place = i
			break

	dataarray = pdata.GetArray(place)
	DataArray =[ ]
	if PointData[place][1] == 1:	
		for i in range(0,npoints):
			a = dataarray.GetTuple1(i)
			DataArray.append([a])
	elif PointData[place][1] == 2:	
		for i in range(0,npoints):
			a,b = dataarray.GetTuple2(i)
			DataArray.append([a,b])
	elif PointData[place][1] == 3:	
		for i in range(0,npoints):
			a,b,c = dataarray.GetTuple3(i)
			DataArray.append([a,b,c])
	elif PointData[place][1] == 4:	
		for i in range(0,npoints):
			a,b,c,d = dataarray.GetTuple4(i)
			DataArray.append([a,b,c,d])
	elif PointData[place][1] == 9:	
		for i in range(0,npoints):
			a,b,c,d,e,f,g,h,i = dataarray.GetTuple9(i)
			DataArray.append([a,b,c,d,e,f,g,h,i])
	else:
		print 'error'

	return DataArray
#-------------------------------------------------------------------------------------------#
def GetCellData(data):
	cdata = data.GetCellData()
	ncarrays = cdata.GetNumberOfArrays()
	CellData=[]
	for i in range(0,ncarrays):
		array = cdata.GetArray(i)
		name = array.GetName()
		numberofcomponents = array.GetNumberOfComponents()
		CellData.append([name,numberofcomponents])
	return CellData
#CellData has the name and the number of components of each field
#-------------------------------------------------------------------------------------------#
def GetCellDataArray(data,CellData,name):
	cells = data.GetCells()
	ncells = cells.GetNumberOfCells()
	cdata = data.GetCellData()
	ncarrays = cdata.GetNumberOfArrays()
	place = 0 #default value
	for i in range(0,ncarrays):
		if CellData[i][0] == name:
			place = i
			break

	dataarray = cdata.GetArray(place)
	CellDataArray =[ ]
	if CellData[place][1] == 1:	
		for i in range(0,ncells):
			a = dataarray.GetTuple1(i)
			CellDataArray.append([a])
	elif CellData[place][1] == 2:	
		for i in range(0,ncells):
			a,b = dataarray.GetTuple2(i)
			CellDataArray.append([a,b])
	elif CellData[place][1] == 3:	
		for i in range(0,ncells):
			a,b,c = dataarray.GetTuple3(i)
			CellDataArray.append([a,b,c])
	elif CellData[place][1] == 4:	
		for i in range(0,ncells):
			a,b,c,d = dataarray.GetTuple4(i)
			CellDataArray.append([a,b,c,d])
	elif CellData[place][1] == 9:	
		for i in range(0,ncells):
			a,b,c,d,e,f,g,h,i = dataarray.GetTuple9(i)
			CellDataArray.append([a,b,c,d,e,f,g,h,i])
	else:
		print 'error'

	return CellDataArray
