# 
include $(STATFEMROOT)/statFEM_definitions.mk

# includes
INCLUDES = -I$(STATFEMROOT) -I$(EIGEN_DIR)


# common flags
CPPFLAGS = $(INCLUDES)

all : abaqus2smf

rebuild : clean all

clean :
	rm -f abaqus2smf

abaqus2smf : main.cpp
	g++ $(CPPFLAGS) main.cpp -o abaqus2smf 
