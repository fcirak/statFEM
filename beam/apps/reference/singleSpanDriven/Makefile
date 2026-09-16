#######################################
# MAKEFILE for this directory         #
#######################################

include $(STATFEMROOT)/statFEM_definitions.mk

DEBUG  = YES	
SOLVER = 
HIRES  =

# compile/link executables and flags
SHELL = /bin/bash
CXX   = $(CPLUSPLUS)
CC    = $(CPLUSPLUS)

# includes
INCLUDES = -I$(STATFEMROOT)  \
           -I$(BOOST_DIR)  \
           -I$(EIGEN_DIR)
#           -I$(SUPERLU_INC)

# libraries
LIBS = 
#LDLIBS = -L$(SUPERLU_LIB_PATH) $(SUPERLU_LIB) 

# common flags
CPPFLAGS = $(INCLUDES)

# high resolution output
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

# for cleaning up
RM = -rm -f

TARGETS = singleSpanDriven

#######################################
# the rules
#

# if nothing is given, make the target
all: $(TARGETS)

# link the target
##$(TARGET): $(OBJS) 
##	$(CXX) $(LDFLAGS) -o $@ $(OBJS) $(LIBS) 
##
### implicit rule for makin object files
##.cpp.o:	
##	$(CXX) $(CPPFLAGS) -c $<

# clean
clean:
	$(RM) *.o $(TARGETS)

