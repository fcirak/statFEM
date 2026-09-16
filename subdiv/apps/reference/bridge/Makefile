#######################################
# MAKEFILE for this directory         #
#######################################

include $(STATFEMROOT)/statFEM_definitions.mk

# determine modus of compilation {YES, NO}
DEBUG = NO

# determine the type of shape functions 
# {CATMULL_CLARK, LOOP}
SHAPEFUN = LOOP 

# compile/link executables and flags
SHELL = /bin/bash
CXX = $(CPLUSPLUS)
LD  = $(CPLUSPLUS)

# header locations
INCLUDES  = -I$(STATFEMROOT) -I$(EIGEN_DIR)
INCLUDES += -I$(BOOST_DIR)

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

# set flags depending on mode
ifeq ($(DEBUG),YES) 
  CPPFLAGS += $(DEBDEFINES) $(DEBCPPFLAGS) -pg 
  LDFLAGS  += $(DEBLDFLAGS) -pg -g
else
  CPPFLAGS += $(RELDEFINES) $(RELCPPFLAGS) 
  LDFLAGS  += $(RELLDFLAGS)
endif

# the final target
TARGET = bridge

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

# clean: PETSC has its own clean
clean:
	$(RM) $(OBJS) $(TARGET)

# temp
bla:
	echo $(LOCALDEFSFILE)

#######################################
# EOF
