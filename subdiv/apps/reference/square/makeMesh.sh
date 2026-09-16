# environment
if [ ! -n "$OPENFTLROOT" ] ; then
    if [ ! -n "$STATFEMROOT" ] ; then
        echo "Neither \$OPENFTLROOT nor \$STATFEMROOT is set"
        exit
    fi
    OPENFTLROOT=$STATFEMROOT
fi

# Gmsh input file name
for name in input/meshQuad input/meshTria
do
    
    echo "---> processing $name"

    # Generate mesh
    #/Applications/Gmsh.app/Contents/MacOS/Gmsh -v 0 -2 -order 1 $name.geo
    gmsh -v 0 -2 -order 1 $name.geo

    echo "---> mesh generated in gmsh"

    # determine tri xor quad mesh
    msh=$name.msh
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

    # convert to SMF
    $OPENFTLROOT/tools/input/gmsh/gmsh2smf -e $numNodesPerElem -f $numNodesPerFace -i $name.msh 
    echo "---> smf mesh generated"
    # convert to VTU
    $OPENFTLROOT/tools/input/smf2vtu/smf2vtu $name.smf $name.vtu
    echo "---> vtu visualisation file generated" 
done


 
