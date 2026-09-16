# check environment variable
if [ ! -n "$STATFEMROOT" ] ; then
    echo "\$STATFEMROOT=$STATFEMROOT not set"
    exit
fi

# Gmsh input file name
geo=input/bridge.geo

# make mesh
msh=${geo%.*}.msh
#/Applications/Gmsh.app/Contents/MacOS/Gmsh -v 0 -2 -order 1 $geo -o $msh
gmsh -v 0 -2 -order 1 -format msh2 $geo -o $msh

echo "---> created MSH file"

# determine tri xor quad mesh
if [ `grep '^[0-9]*[ ][2][ ]' -c $msh` -ne "0" ] ; then
    simplex=TRIANGLE
    numNodesPerElem=3
    numNodesPerFace=2
elif [ `grep '^[0-9]*[ ][3][ ]' -c $msh` -ne "0" ] ; then
    simplex=QUADRILATERAL
    numNodesPerElem=4
    numNodesPerFace=2
else
    echo "Cannot handle Gmsh file"
    exit
fi

# make SMF file
smf=${geo%.*}.smf
$STATFEMROOT/tools/input/gmsh/gmsh2smf -e $numNodesPerElem -f $numNodesPerFace -i $msh
echo "---> created SMF file"

# make VTU file
vtu=${geo%.*}.vtu
$STATFEMROOT/tools/input/smf2vtu/smf2vtu $smf $vtu
echo "---> created VTU file"
