#######################################
# MAKEFILE for this directory         #
#######################################

include $(STATFEMROOT)/statFEM_definitions.mk

# determine modus of compilation {YES, NO}
DEBUG = NO

# determine the type of shape functions 
# {CATMULL_CLARK, LOOP}
SHAPEFUN = LOOP

# dynamic mode
# {BE,NEWMARK}
DYN = BE

# set compiler: INTEL,[unset implies GNU]
#OPENFTLCOMPENV = 

# set special solver: SUPERLU,PETSC,PARDISO,[unset implies Eigen Sparse Cholesky]
#SOLVER = SUPERLU

# compile/link executables and flags
SHELL = /bin/bash
CXX = $(CPLUSPLUS)
LD  = $(CPLUSPLUS)

# header locations
INCLUDES  = -I$(STATFEMROOT) -I$(EIGEN_DIR)
INCLUDES += -I$(BOOST_DIR)

# common flags
CPPFLAGS = $(INCLUDES) -fpermissive

# libraries
LIBS += -lboost_timer -lboost_chrono

# Catmull Clark subdivision
ifeq ($(strip $(SHAPEFUN)),CATMULL_CLARK)
  CPPFLAGS += -DCATMULL_CLARK_
else 
  # Loop subdivision 
  ifeq ($(strip $(SHAPEFUN)),LOOP)
    CPPFLAGS += -DLOOP_
  endif
endif

# dynamics
ifeq ($(strip $(DYN)),BE)
  CPPFLAGS += -DDYNAMICS_BACKWARDEULER
endif
ifeq ($(strip $(DYN)),NEWMARK)
  CPPFLAGS += -DDYNAMICS_NEWMARK
endif

# set flags depending on mode
ifeq ($(strip $(DEBUG)),YES) 
  CPPFLAGS += $(DEBDEFINES) $(DEBCPPFLAGS) -pg 
  LDFLAGS  += $(DEBLDFLAGS) -pg -g
else
  CPPFLAGS += $(RELDEFINES) $(RELCPPFLAGS) 
  LDFLAGS  += $(RELLDFLAGS)
endif

# set PETSC stuff
#ifeq ($(SOLVER),PETSC)
#  include $(PETSC_DIR)/$(PETSC_CONF_DIR)/variables
#  CPPFLAGS += $(PETSC_INCLUDE) -DPETSC
#  LIBS += $(PETSC_LIB)
#  CXX   = mpicxx
#  LD    = mpicxx
#endif
## set SuperLU stuff
#ifeq ($(SOLVER),SUPERLU)
#  CPPFLAGS += -I$(SUPERLU_INC) -DSUPERLU
#  LIBS += -L$(SUPERLU_LIB_PATH) $(SUPERLU_LIB)
#endif
## set Pardiso stuff
#ifeq ($(SOLVER),PARDISO)
#  ifeq (,$(findstring INTEL,$(OPENFTLCOMPENV)))
#    $(error Pardiso solver only available with Intel compiler: set OPENFTLCOMPENV=INTEL)  # _no_ tab before $(error ...)
#  endif
#  CPPFLAGS += -DPARDISO_MKL -openmp
#  LDFLAGS += -mkl=parallel
#endif

# find the sources and object files
CPPFILES = $(wildcard *.cpp)
CPPOBJS = $(CPPFILES:.cpp=.o)

# the final target
TARGET = $(CPPFILES:.cpp=)

# define sources and objects
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

# clean: PETSC has its own clean
clean:
	$(RM) $(OBJS) $(TARGET)

# temp
bla:
	echo $(LOCALDEFSFILE)

#######################################
# EOF
