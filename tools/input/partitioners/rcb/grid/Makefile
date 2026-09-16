#######################################
# MAKEFILE for this directory         #
#######################################

include $(STATFEMROOT)/statFEM_definitions.mk

# determine modus of compilation
DEBUG=NO

# compile/link executables and flags
SHELL = /bin/bash
CXX   = $(CPLUSPLUS)
CC    = $(CPLUSPLUS)

# includes
INCLUDES  = -I$(STATFEMROOT)
INCLUDES += -I$(BOOST_DIR)
INCLUDES += -I$(EIGEN_DIR)

# libraries
LDLIBS  = -lboost_timer -lboost_chrono

# common flags
CPPFLAGS = $(INCLUDES) 

# set flags depending on mode
ifeq ($(strip $(DEBUG)),YES) 
	CPPFLAGS += $(DEBDEFINES) $(DEBCPPFLAGS) 
	LDFLAGS  += $(DEBLDFLAGS)
else
	CPPFLAGS += $(RELDEFINES) $(RELCPPFLAGS) 
	LDFLAGS  += $(RELLDFLAGS)
endif

# for cleanup
RM = -rm -f

TARGETS = partitioner

#######################################
# the rules
#
all: $(TARGETS)

# clean: 
clean:
	$(RM) *.o $(TARGETS) core

# implict rules:
# $(CC) $(LDFLAGS) n.o $(LOADLIBES) $(LDLIBS)
# $(CXX) -c $(CPPFLAGS) $(CXXFLAGS)