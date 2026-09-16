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
CXX   = $(MPICPLUSPLUS)
LD    = $(MPICPLUSPLUS)

# includes
INCLUDES = -I. -I$(STATFEMROOT) -I$(BOOST_DIR) -I$(EIGEN_DIR) -I$(METIS_INC)

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

LIBS = -Wl,-rpath,$(METIS_LIB_PATH) -L$(METIS_LIB_PATH) $(METIS_LIB) -lboost_timer -lboost_chrono


# detect local headers to set up rules properly
LOCHEADERS = PartitionerMetis.hpp  \
             $(STATFEMROOT)/tools/input/smf2vtu/readwrite.hpp
HEADERS = $(LOCHEADERS)

#######################################
# the rules
#
all: partitioner psmf2smf

partitioner: partitioner.cpp $(HEADERS)
	$(CXX) $(CPPFLAGS) -o $@ $< $(LIBS)

psmf2smf: ../utils/psmf2smf.cpp $(HEADERS)
	$(CXX) $(CPPFLAGS) -o $@ $< $(LIBS)

# test executables
test:
	(cd tests ; $(MAKE))

# clean
clean:
	-rm -f $(OBJS) partitioner psmf2smf

