#######################################
# MAKEFILE for this directory         #
#######################################


# include system-specific configuration
include $(STATFEMROOT)/statFEM_definitions.mk

# set compiler: INTEL,[unset implies GNU]
ifneq ($(origin OPENFTLCOMPENV),environment)
  OPENFTLCOMPENV = 
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



# the final targets
TARGETS = staticDrivenFine staticDrivenMedium staticDrivenCoarse

# for cleanup
RM = -rm -f

#######################################
# the rules
#
all: $(TARGETS)

# link the target
staticDrivenFine: staticDrivenFine.o
	$(LD) $(LDFLAGS) -o $@ $< $(LIBS) 

staticDrivenMedium: staticDrivenMedium.o
	$(LD) $(LDFLAGS) -o $@ $< $(LIBS) 

staticDrivenCoarse: staticDrivenCoarse.o
	$(LD) $(LDFLAGS) -o $@ $< $(LIBS) 

# implicit rule for making the object file
.cpp.o:	
	$(CXX) $(CPPFLAGS) -c $<

# clean: 
clean:
	$(RM) $(OBJS) $(TARGETS)

