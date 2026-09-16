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
INCLUDES = -I$(STATFEMROOT) -I$(BOOST_DIR) -I$(EIGEN_DIR)

# LU solver
INCLUDES += -I$(SUPERLU_INC)
LDLIBS += -L$(SUPERLU_LIB_PATH) $(SUPERLU_LIB)

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
TARGET = potentialFlow

# find the sources and object files
CPPFILES = main.cpp
CPPOBJS = main.o

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
	$(CXX) $(LDFLAGS) -o $@ $(OBJS) $(LIBS) $(LDLIBS)

# implicit rule for makin object files
.cpp.o:	
	$(CXX) $(CPPFLAGS) -c $<

clean:
	$(RM) $(OBJS) $(TARGET)

