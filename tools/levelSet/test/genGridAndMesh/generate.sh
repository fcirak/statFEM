#! /bin/bash

## make sure user calls correctly
if test $# -ne 1; then
    echo "Usage:   generate.sh <file-name>.<ext>"
    echo ""
    echo "Example: generate.sh box2D.dat          to generate a unit grid in 2D"
    echo "Example: generate.sh bridge.geo         to generate a mesh for the bridge geometry"
    exit
fi

## extract basename and extension
FILENAME="$1"
BASENAME=${FILENAME%.*}
EXTENAME=${FILENAME#*.}

if [ "${EXTENAME}" == "dat" ]; then
    ## generate grids
    ##--------------------------------------------------------------------------
    GENGRID="../../../input/meshing/genGrid/gridFile2Smf"

    ${GENGRID} ${FILENAME}
    mv ${BASENAME}.smf grid.smf

    echo "Generated grid for ${BASENAME}, wrote to grid.smf"
elif [ "${EXTENAME}" == "geo" ]; then
    ## generate meshes
    ##--------------------------------------------------------------------------
    GMSHTOSMF="../../../input/gmsh/gmsh2smf"

    ## 2D surfaces
    if [ "${BASENAME}" == "bridge" ] || [ "${BASENAME}" == "sphere" ] || [ "${BASENAME}" == "plate" ]; then
	gmsh -2 -order 1 ${FILENAME}
	${GMSHTOSMF} -e 3 -f 2 -i ${BASENAME}.msh
	mv ${BASENAME}.smf mesh.smf
    ## 1D lines
    elif [ "${BASENAME}" == "circle" ] || [ "${BASENAME}" == "line" ] || [ "${BASENAME}" == "zigzag" ]; then
	gmsh -2 -order 1 ${FILENAME}
	${GMSHTOSMF} -e 2 -f 2 -i ${BASENAME}.msh
	mv ${BASENAME}.smf mesh.smf
    else
	echo "Error: Mesh not recognised"
	exit
    fi

    echo "Generated mesh for ${BASENAME}, wrote to mesh.smf"
    rm *cnstr *msh
else
    echo "Incorrect file extension, should be .dat or .geo"
fi