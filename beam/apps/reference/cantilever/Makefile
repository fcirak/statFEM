#######################################
# MAKEFILE for this directory         #
#######################################

include $(STATFEMROOT)/statFEM_definitions.mk

# determine modus of compilation
DEBUG = YES
HIRES = 0

# flags for the debug mode
DEBCPPFLAGS = -g  -Wall
DEBLDFLAGS  = -g  -Wall
DEBDEFINES  = -DDEBUG

# flags for release mode
RELCPPFLAGS = -O3
RELLDFLAGS  = -O3
RELDEFINES  = -DNDEBUG

# compile/link executables and flags
SHELL = /bin/bash
CXX = $(CPLUSPLUS)
LD = $(CPLUSPLUS)

# includes
INCLUDES = -I$(STATFEMROOT)  \
           -I$(BOOST_DIR)  \
           -I$(EIGEN_DIR)

# libraries
LIBS = 

# common flags
CPPFLAGS = $(INCLUDES) -fpermissive

ifeq ($(HIRES),1)
	CPPFLAGS += -DD_HIRES
endif

# set flags depending on mode
ifeq ($(DEBUG),YES) 
	CPPFLAGS += $(DEBDEFINES) $(DEBCPPFLAGS) 
	LDFLAGS  += $(DEBLDFLAGS)
else
	CPPFLAGS += $(RELDEFINES) $(RELCPPFLAGS) 
	LDFLAGS  += $(RELLDFLAGS)
endif

# std library flag
LDFLAGS  += -lstdc++

# dependency generation
DEPEND = makedepend $(CPPFLAGS)
DEPFLAGS = -Y

# the final target
TARGET = $(shell basename $(PWD))

# find the sources and object files
CPPFILES = $(wildcard *.cpp)
CPPOBJS = $(CPPFILES:.cpp=.o)

#define sources and objects
SRCS = $(CPPFILES)
OBJS = $(CPPOBJS)

# for cleaning up
RM = -rm -f

#######################################
# the rules
#

# if nothing is given, make the target
all: $(TARGET)

# link the target
$(TARGET): $(OBJS) 
	$(CXX) $(LDFLAGS) -o $@ $(OBJS) $(LIBS) 

# the source dependencies
depend: 
	$(DEPEND) $(DEPFLAGS) $(SRCS)

# implicit rule for makin object files
.cpp.o:	
	$(CXX) $(CPPFLAGS) -c $<

# clean
clean:
	$(RM) $(OBJS) $(TARGET)

#############################################################################
# DO NOT DELETE
