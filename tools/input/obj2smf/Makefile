#######################################
# MAKEFILE for this directory         #
#######################################

include $(STATFEMROOT)/statFEM_definitions.mk

INCLUDES  = -I$(STATFEMROOT) -I$(EIGEN_DIR)

# common flags
CPPFLAGS = $(INCLUDES)

all : obj2smf

rebuild : clean all

clean :
	rm -f obj2smf

obj2smf : obj2smf.cpp
	g++ -g3 $(CPPFLAGS) obj2smf.cpp -o obj2smf $(LIBS)
