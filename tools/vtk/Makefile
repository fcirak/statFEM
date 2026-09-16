# ------------------------------------------------------------------------------
include $(STATFEMROOT)/statFEM_definitions.mk
# ------------------------------------------------------------------------------
DEBUG=

# compile/link executables and flags
SHELL = /bin/bash
CXX   = $(CPLUSPLUS)
CC    = $(CPLUSPLUS)

# initialise includes
INCLUDES  = -I$(STATFEMROOT) -I$(BOOST_DIR)

# VTK headers and libraries
ifneq ($(origin VTK_PATH),undefined)
  VTK_INC = $(VTK_PATH)/include/vtk-5.6
  VTK_LIB_PATH = $(VTK_PATH)/lib/vtk-5.6
  VTK_LIB = -Wl,-rpath,$(VTK_LIB_PATH) -L$(VTK_LIB_PATH)
else
  VTK_PATH = /usr/include
  VTK_INC = $(VTK_PATH)/$(shell basename $(VTK_PATH)/vtk-5.*)
endif

INCLUDES += -I$(VTK_INC)
VTK_LIB += -lvtkHybrid -lvtkIO -lvtkFiltering -lvtkCommon -lvtkGraphics

# compilation flags
CPPFLAGS = $(INCLUDES) -g3  #-03

# libraries
LIBS = $(VTK_LIB)

# location of your TECIO header
TECINCL = $(TEC_360_2010)/include

# location of your tecio.a library
TECLIB = $(TEC_360_2010)/lib/tecio.a

# set flags depending on mode
ifeq ($(DEBUG),YES) 
	CPPFLAGS += $(DEBDEFINES) $(DEBCPPFLAGS) 
	LDFLAGS  += $(DEBLDFLAGS)
else
	CPPFLAGS += $(RELDEFINES) $(RELCPPFLAGS) 
	LDFLAGS  += $(RELLDFLAGS)
endif

# libraries
LDLIBS  = $(LIBS)

# ------------------------------------------------------------------------------
# standard targets, compilable with implicit rules
TARGETS = boxClip deformedShape extractPointDatumOverTime plotOverLine \
					vorticity vtkAddVector vtpClean vtpProbe vtu2vtp \
					vtuComp vtuDiff vtuMerge vtuPointLocator vtx2vtxbin 	

# target which requires mpi
MPITARGETS = processFilesParallel 

# tecplot targets
TECTARGETS = vtu2tec vtp2tec vtu2tecTransient vts2tec

# ------------------------------------------------------------------------------
all: $(TARGETS) $(MPITARGETS) $(TECTARGETS)

clean:
	rm -f  *.o $(TARGETS) $(MPITARGETS) $(TECTARGETS) core

# ------------------------------------------------------------------------------
# special rules with tecplot header and library
vtu2tec: 
	$(CPLUSPLUS) vtu2tec.cpp -o $@ -I$(TECINCL) $(CPPFLAGS) $(LIBS) $(TECLIB)

vtp2tec: 
	$(CPLUSPLUS) vtp2tec.cpp -o $@ -I$(TECINCL) $(CPPFLAGS) $(LIBS) $(TECLIB)

vtu2tecTransient: 
	$(CPLUSPLUS) vtu2tecTransient.cpp -o $@ -I$(TECINCL) $(CPPFLAGS) $(LIBS) $(TECLIB)

vts2tec: 
	$(CPLUSPLUS) vts2tec.cpp -o $@ -I$(TECINCL) $(CPPFLAGS) $(LIBS) $(TECLIB)

# special rule with mpi
processFilesParallel: 
	$(MPICPLUSPLUS) processFilesParallel.cpp -o $@ $(CPPFLAGS) $(LIBS) 

