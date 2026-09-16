#######################################
# MAKEFILE for this directory         #
#######################################


# include system-specific configuration
include $(STATFEMROOT)/statFEM_definitions.mk

# set compiler: INTEL,[unset implies GNU]
ifeq ($(origin OPENFTLCOMPENV),undefined)
  OPENFTLCOMPENV =
endif
#OPENFTLCOMPENV = 

# determine modus of compilation: YES,[unset implies release]
DEBUG =	

# set special solver: SUPERLU,PETSC,PARDISO,[unset implies uBLAS Cholesky]
SOLVER = 

# compile/link executables and flags
SHELL = /bin/bash
CXX   = $(CPLUSPLUS)
LD    = $(CPLUSPLUS)

# includes
INCLUDES  = -I$(STATFEMROOT)
INCLUDES += -I$(BOOST_DIR)
INCLUDES += -I$(EIGEN_DIR)

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

# set PETSC stuff
ifeq ($(SOLVER),PETSC)
  include $(PETSC_DIR)/$(PETSC_CONF_DIR)/variables
  CPPFLAGS += $(PETSC_INCLUDE) -DPETSC
  LIBS += $(PETSC_LIB)
  CXX   = $(MPICPLUSPLUS)
  LD    = $(MPICPLUSPLUS)
endif
# set SuperLU stuff
ifeq ($(SOLVER),SUPERLU)
  CPPFLAGS += -I$(SUPERLU_INC) -DSUPERLU
  LIBS += -L$(SUPERLU_LIB_PATH) $(SUPERLU_LIB)
endif
# set Pardiso stuff
ifeq ($(SOLVER),PARDISO)
  ifeq (,$(findstring INTEL,$(OPENFTLCOMPENV)))
    $(error Pardiso solver only available with Intel compiler: set OPENFTLCOMPENV = INTEL??)  # no tab before $(error ...)
  endif
  CPPFLAGS += -DPARDISO_MKL
  LDFLAGS += -mkl=parallel
endif

# the final target
TARGET = static

# find the sources and object files
CPPFILES = $(wildcard *.cpp)
CPPOBJS = $(CPPFILES:.cpp=.o)

#define sources and objects
SRCS = $(CPPFILES)
OBJS = $(CPPOBJS)

# for cleanup
RM = -rm -f

#######################################
# the rules
#
all: $(TARGET)

# link the target
$(TARGET): $(OBJS) 
	$(LD) $(LDFLAGS) -o $@ $(OBJS) $(LIBS) 

# implicit rule for making the object file
.cpp.o:	
	$(CXX) $(CPPFLAGS) -c $<

# clean: 
clean:
	$(RM) $(OBJS) $(TARGET)

