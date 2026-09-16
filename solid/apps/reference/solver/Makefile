#######################################
# MAKEFILE for this directory         #
#######################################

include $(STATFEMROOT)/statFEM_definitions.mk

# determine modus of compilation
DEBUG=

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
CPPFLAGS = $(INCLUDES) -fpermissive

# set flags depending on mode
ifeq ($(DEBUG),YES) 
	CPPFLAGS += $(DEBDEFINES) $(DEBCPPFLAGS) 
	LDFLAGS  += $(DEBLDFLAGS)
else
	CPPFLAGS += $(RELDEFINES) $(RELCPPFLAGS) 
	LDFLAGS  += $(RELLDFLAGS)
endif

# set PETSC stuff
PETSC_CC_INCLUDES = $(shell pkg-config --cflags PETSc)
PETSC_LIB = $(shell pkg-config --libs PETSc)
PETSCINCLUDES = $(PETSC_CC_INCLUDES) #defined in file above

TARGETS = solver_snes solver_ublas solver_petsc

# all: make all the targets
all: $(TARGETS)

#------------------------------------------------------------
# solver with petsc's snes wrapper			
solver_snes: solver_snes.o
	$(MPICPLUSPLUS) $(LDFLAGS) -o $@ solver_snes.o $(LIBS) $(PETSC_LIB)

solver_snes.o:
	$(MPICPLUSPLUS) $(CPPFLAGS) $(PETSCINCLUDES) -c solver_snes.cpp

#------------------------------------------------------------
# solver based on a ublas cholesky decomposition
solver_ublas: solver_ublas.o
	$(LD) $(LDFLAGS) -o $@ solver_ublas.o $(LIBS)

solver_ublas.o:
	$(CXX) $(CPPFLAGS) -c solver_ublas.cpp

#------------------------------------------------------------
# solver with petsc for the linear solve
solver_petsc: solver_petsc.o
	$(MPICPLUSPLUS) $(LDFLAGS) -o $@ solver_petsc.o $(LIBS) $(PETSC_LIB)

solver_petsc.o:
	$(MPICPLUSPLUS) $(CPPFLAGS) $(PETSCINCLUDES) -c solver_petsc.cpp

#------------------------------------------------------------

# for cleaning up
RM = -rm -f

clean:
	$(RM) *.o $(TARGETS)

