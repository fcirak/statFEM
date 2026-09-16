#######################################
# MAKEFILE for this directory         #
#######################################

include $(STATFEMROOT)/statFEM_definitions.mk

# determine modus of compilation
DEBUG   =

# choose if you want extra output
OUTPUT  =
VERBOSE = YES

# generate info for profiler
PROFILE = 

# compile/link executables and flags
SHELL = /bin/bash
CXX   = $(CPLUSPLUS)
LD    = $(CPLUSPLUS)

# includes
INCLUDES = -I. -I$(STATFEMROOT) -I$(BOOST_DIR) -I$(EIGEN_DIR)

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

# add extra file output
ifeq ($(strip $(OUTPUT)),YES)
	CPPFLAGS += -DOUTPUT
endif

# add extra statistics
ifeq ($(strip $(VERBOSE)),YES)
	CPPFLAGS += -DVERBOSE
endif

# add extra flags
ifeq ($(strip $(PROFILE)),YES)
	CPPFLAGS += -p -pg
endif

LIBS = -lboost_timer -lboost_chrono


# detect local headers to set up rules properly
LOCHEADERS = $(STATFEMROOT)/tools/input/smf2vtu/readwrite.hpp
HEADERS = $(LOCHEADERS)

#######################################
# the rules
#
all: smfMidsides

smfMidsides: smfMidsides.cpp $(HEADERS)
	$(CXX) $(CPPFLAGS) -o $@ $< $(LIBS)

# test executables
test:
	(cd tests ; $(MAKE))

# clean
clean:
	-rm -f $(OBJS) smfMidsides

