#######################################
# MAKEFILE for this directory         #
#######################################

include $(STATFEMROOT)/statFEM_definitions.mk


# set compiler: INTEL,[unset implies GNU]
ifneq ($(origin STATFEMCOMPENV),environment)
  STATFEMCOMPENV = 
endif
# determine modus of compilation: YES,[unset implies release]
ifneq ($(origin DEBUG),environment)
  DEBUG = 	
endif
# set special solver: SUPERLU,PETSC,PARDISO,[unset implies uBLAS Cholesky]
ifneq ($(origin SOLVER),environment)
  SOLVER = 
endif


# compile/link executables and flags
SHELL = /bin/bash
CXX   = $(CPLUSPLUS)
LD    = $(CPLUSPLUS)

# includes
INCLUDES = -I$(STATFEMROOT)  \
           -I$(BOOST_DIR)  \
           -I$(EIGEN_DIR)  

# libraries
LIBS = 

# common flags
CPPFLAGS = $(INCLUDES) -fpermissive

# upsampled output
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

# the final target
TARGET = pendulum

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
