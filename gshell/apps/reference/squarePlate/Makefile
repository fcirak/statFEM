#######################################
# MAKEFILE for this directory         #
#######################################

include $(STATFEMROOT)/statFEM_definitions.mk

# determine modus of compilation {YES, NO}
DEBUG = YES

# determine the type of shape functions 
# {CATMULL_CLARK, LOOP}
SHAPEFUN = LOOP

# mesh file { , A, B, D, CC}
MESH = 

# shell kinematics
# {SHEARFLEXIBLE, KIRCHHOFFLOVE}
KINEM = SHEARFLEXIBLE

# dense linear algebra
# {UBLAS,C}
LINALG = C

# compile/link executables and flags
SHELL = /bin/bash
CXX = $(CPLUSPLUS)
LD  = $(CPLUSPLUS)

# header locations
INCLUDES  = -I$(STATFEMROOT)
INCLUDES += -I$(EIGEN_DIR)
INCLUDES += -I$(BOOST_DIR)
#INCLUDES += -I$(SUPERLU_INC)

# libraries
#LIBS += -L$(SUPERLU_LIB_PATH)
#LIBS += $(SUPERLU_LIB)

# common flags
CPPFLAGS = $(INCLUDES) -fpermissive

# Catmull Clark subdivision
ifeq ($(strip $(SHAPEFUN)),CATMULL_CLARK)
  CPPFLAGS += -DCATMULL_CLARK_
else 
  # Loop subdivision 
  ifeq ($(strip $(SHAPEFUN)),LOOP)
    CPPFLAGS += -DLOOP_
  endif
endif

# Mesh type
ifeq ($(strip $(MESH)),A)
  CPPFLAGS += -DMESH_A
else
  ifeq ($(strip $(MESH)),B)
    CPPFLAGS += -DMESH_B
  else
    ifeq ($(strip $(MESH)),D)
      CPPFLAGS += -DMESH_D
    else
      ifeq ($(strip $(MESH)),CC)
        CPPFLAGS += -DMESH_CC
      endif
    endif
  endif
endif

# kinematics
ifeq ($(strip $(KINEM)),SHEARFLEXIBLE)
  CPPFLAGS += -DKINEMATICS_SHEARFLEXIBLE
endif
ifeq ($(strip $(KINEM)),KIRCHHOFFLOVE)
  CPPFLAGS += -DKINEMATICS_KIRCHHOFFLOVE
endif

# linalg
ifeq ($(strip $(LINALG)),C)
  CPPFLAGS += -DLINALG_C
else
  CPPFLAGS += -DLINALG_UBLAS
endif

# profiling
ifeq ($(PROFILING),YES)
  CPPFLAGS += -pg 
  LDFLAGS  += -pg
endif

# set flags depending on mode
ifeq ($(DEBUG),YES) 
  CPPFLAGS += $(DEBDEFINES) $(DEBCPPFLAGS)
  LDFLAGS  += $(DEBLDFLAGS) -g
else
  CPPFLAGS += $(RELDEFINES) $(RELCPPFLAGS) 
  LDFLAGS  += $(RELLDFLAGS)
endif

# std library flag
LDFLAGS  += -lstdc++ 

# find the sources and object files
CPPFILES = $(wildcard *.cpp)
CPPOBJS = $(CPPFILES:.cpp=.o)

# the final target
TARGET = squarePlate

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

# implicit rule for makin object files
.cpp.o:	
	$(CXX) $(CPPFLAGS) -c $<

# clean: PETSC has its own clean
clean:
	$(RM) $(OBJS) $(TARGET) core

# temp
bla:
	echo $(LOCALDEFSFILE)

#######################################
# EOF
