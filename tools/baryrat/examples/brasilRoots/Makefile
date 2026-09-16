#######################################
# MAKEFILE for this directory         #
#######################################
include $(STATFEMROOT)/statFEM_definitions.mk

# determine modus of compilation
DEBUG=YES

# compile/link executables and flags
SHELL = /bin/bash
CXX   = $(CPLUSPLUS)
LD    = $(CPLUSPLUS)

# includes
INCLUDES = -I$(STATFEMROOT)

# libraries
LIBS  =

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

# the final target
TARGET = brasilRoots

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

# implicit rule for makin object files
.cpp.o:	
	$(CXX) $(CPPFLAGS) -c $<

clean:
	$(RM) $(OBJS) $(TARGET)

