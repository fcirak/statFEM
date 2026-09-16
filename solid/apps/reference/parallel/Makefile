#######################################
# MAKEFILE for this directory         #
#######################################

# env
# set compiler: INTEL,[unset implies GNU]
ifeq ($(origin OPENFTLCOMPENV),undefined)
  OPENFTLCOMPENV =
endif

# determine modus of compilation
DEBUG =
# linear solver: PETSC,BDDCML,[unset implies errors]
SOLVER = PETSC

# pull this in _after_ modi (above) are set
include $(STATFEMROOT)/statFEM_definitions.mk

# compile/link executables and flags
SHELL = /bin/bash

# includes
INCLUDES  = -I$(STATFEMROOT)
INCLUDES += -I$(BOOST_DIR)
INCLUDES += -I$(EIGEN_DIR)

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

# set solver environment
ifeq ($(strip $(SOLVER)),PETSC)
        # set PETSC environment
        PETSC_CC_INCLUDES = $(shell pkg-config --cflags PETSc)
        PETSC_LIB = $(shell pkg-config --libs PETSc)
        CPPFLAGS += $(PETSC_CC_INCLUDES) -DPETSC
        LIBS += $(PETSC_LIB)
endif
ifeq ($(strip $(SOLVER)),BDDCML)
        # set BDDCML library and required packages ( MUMPS, ParMETIS, ScaLAPACK, BLOPEX, LAPACK )
        include $(BDDCML_DIR)/make.inc
        CPPFLAGS += $(BDDCML_INC) -DBDDCML
        LIBS += $(BDDCML_LINK_SEQUENCE)

        # for cross-compilation working well
        LIBS += -L$(FORTRAN_LIB_PATH) $(FORTRAN_LIB)
        # add preprocessor definition to produce correct names of Fortran symbols
        CPPFLAGS += $(C2FDEF)
endif

# profiling
ifeq ($(PROFILE),YES)
        CPPFLAGS += -DSYSTEMSOLVEPETSC_PROFILE
endif

# executable
TARGET = parallel

# find the sources and object files
CPPFILES = $(wildcard *.cpp)
CPPOBJS = $(CPPFILES:.cpp=.o)

#define sources and objects
SRCS = $(CPPFILES)
OBJS = $(CPPOBJS)

# for cleaning up
RM = rm -f

#######################################
# the rules
#
all: $(TARGET)

# link the target
$(TARGET): $(OBJS) 
	$(MPICPLUSPLUS) $(LDFLAGS) -o $@ $(OBJS) $(LIBS) 

# implicit rule for makin object files
.cpp.o:	
	$(MPICPLUSPLUS) $(CPPFLAGS) -c $<

# clean: PETSC has its own clean
clean:
	$(RM) $(OBJS) $(TARGET)

