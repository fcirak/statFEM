#------------------------------------------------------------------------------#
# Divide the Cells array into three categories: 
#   -  plusCells which are entirely on the positive side of the interface
#   - minusCells which are entirely on the negative side of the interface
#   -   cutCells which are traversed by the interface
# argin: 
#   - Cells:    The list of cells
#   - Distance: Array of the nodal signed distances to the interface
# returns: three lists of cells
def CellIdentifier(Cells,Distance):
	plusCells    = []
	minusCells   = []
	cutCells     = []
	for cell in Cells:
                # retrieve the distance for each cell vertex
		dcell=[Distance[j][0] for j in cell[0]]
		if max(dcell)*min(dcell)<=0 :
			cutCells.append(cell)
		elif min(dcell)>=0 :
			plusCells.append(cell)
		elif max(dcell)<=0 :
			minusCells.append(cell)		
	return cutCells, minusCells, plusCells

#------------------------------------------------------------------------------#
# Change the order of the cell vertices for quadratic elements. The input cell
# is ordered according to VTK specifications. The output cell will be ordered
# such that two adjacent vertices are adjacent in the list of vertices.
# argin:
#  - cell:  A cell
# returns: reordered cell
def orderCell(cell):
        type = cell[2]
        new_cell = cell;
        if type == 22:
                vertices = [cell[0][0], cell[0][3], cell[0][1], \
                                    cell[0][4], cell[0][2], cell[0][5]]
                new_cell = [vertices, cell[1], cell[2]]
        elif type == 23:
                vertices = [cell[0][0], cell[0][4], cell[0][1], cell[0][5],\
                                    cell[0][2], cell[0][6], cell[0][3], cell[0][7]]
                new_cell = [vertices, cell[1], cell[2]]
        return new_cell

#------------------------------------------------------------------------------#










