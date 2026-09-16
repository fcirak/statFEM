#######################################
# MAKEFILE for this directory         #
#######################################

include $(STATFEMROOT)/statFEM_definitions.mk

# compile/link executables and flags
SHELL = /bin/bash
CXX   = $(CPLUSPLUS)

INCLUDES  = -I$(STATFEMROOT)
INCLUDES += -I$(BOOST_DIR)
INCLUDES += -I$(EIGEN_DIR)

# common flags
CPPFLAGS = $(INCLUDES)

all : smf2off

rebuild : clean all

clean :
	rm -f smf2off

smf2off : smf2off.cpp
	$(CXX) $(CPPFLAGS) smf2off.cpp -o smf2off $(LIBS)
