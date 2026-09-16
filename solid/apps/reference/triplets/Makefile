#######################################
# MAKEFILE for this directory         #
#######################################

# determine modus of compilation
DEBUG=

# pull this in _after_ modi (above) are set
include $(STATFEMROOT)/statFEM_definitions.mk

# compile/link executables and flags
SHELL = /bin/bash
CXX   = $(CPLUSPLUS)
LD    = $(CPLUSPLUS)

# includes
INCLUDES  = -I$(STATFEMROOT) 
INCLUDES += -I$(BOOST_DIR)

# libraries
LIBS = 

# common flags
CPPFLAGS = $(INCLUDES) 

# set flags depending on mode
ifeq ($(DEBUG),YES) 
	CPPFLAGS += $(DEBDEFINES) $(DEBCPPFLAGS) 
	LDFLAGS  += $(DEBLDFLAGS)
else
	CPPFLAGS += $(RELDEFINES) $(RELCPPFLAGS) 
	LDFLAGS  += $(RELLDFLAGS)
endif

# for cleaning up
RM = rm -f

#######################################
# the rules
#
all: testTriplet

# link the targets
testTriplet: testTriplet.cpp
	$(CXX) $(CPPFLAGS) $(LDFLAGS) -o $@ $< $(LIBS) 

# clean
clean:
	$(RM) testTriplet

