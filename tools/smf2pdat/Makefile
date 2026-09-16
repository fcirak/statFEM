# 
include $(STATFEMROOT)/statFEM_definitions.mk

# determine modus of compilation
DEBUG=NO

# compile/link executables and flags
SHELL = /bin/bash
CXX   = $(CPLUSPLUS)
CC    = $(CPLUSPLUS)
LD    = $(CPLUSPLUS)

# includes
INCLUDES = -I$(STATFEMROOT)/
INCLUDES += -I$(BOOST_DIR)
INCLUDES += -I$(EIGEN_DIR)

# libraries
LIBS = -Wl,-rpath,$(BOOST_LIB_PATH) -L$(BOOST_LIB_PATH) \
       -lboost_program_options

# common flags
CPPFLAGS = $(INCLUDES) 

# set flags depending on mode
ifeq ($(strip $(DEBUG)),YES) 
	CPPFLAGS += $(DEBDEFINES) $(DEBCPPFLAGS) 
	LDFLAGS  += $(DEBLDFLAGS)
else
	CPPFLAGS += $(RELDEFINES) $(RELCPPFLAGS) 
	LDFLAGS  += $(RELLDFLAGS)
endif


# the final target
TARGET = smf2pdat

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
