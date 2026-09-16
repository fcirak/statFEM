#! /bin/bash

if [ -e *.cnstr ]; then
    rm *cnstr
fi

if [ -e *.msh ]; then
    rm *msh
fi

if [ -e *.smf ]; then
    rm *smf
fi

if [ -e core ]; then
    rm core
fi
