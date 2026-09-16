# 
include $(STATFEMROOT)/statFEM_definitions.mk

# compile/link executables and flags
SHELL = /bin/bash
CXX   = $(CPLUSPLUS)

# includes
INCLUDES = -I$(STATFEMROOT)/ 
INCLUDES += -I$(BOOST_DIR)

# common flags
CPPFLAGS = $(INCLUDES)

all : off2smf

rebuild : clean all

clean :
	rm -f off2smf

off2smf : off2smf.cpp
	$(CXX) $(CPPFLAGS) off2smf.cpp -o off2smf $(LIBS)
