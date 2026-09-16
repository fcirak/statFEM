# 
include $(STATFEMROOT)/statFEM_definitions.mk

# determine modus of compilation
DEBUG=YES

# verbosity
VERBOSE = 
#VERBOSE = -DVERBOSE

# includes
INCLUDES  = -I$(STATFEMROOT) -I$(BOOST_DIR)
INCLUDES += -I$(EIGEN_DIR)

ifeq ($(DEBUG),YES)
CPPFLAGS += $(DEBDEFINES) $(DEBCPPFLAGS)
LDFLAGS += $(DEBLDFLAGS)
else
CPPFLAGS += $(RELDEFINES) $(RELCPPFLAGS)
LDFLAGS += $(RELLDFLAGS)
endif

# common flags
CPPFLAGS = $(INCLUDES) $(VERBOSE)

LIBS = 

TARGETS = smf2vtu psmf2vtu

all : $(TARGETS)

rebuild : clean all

smf2vtu : smf2vtu.cpp readwrite.hpp
	$(CPLUSPLUS) $(CPPFLAGS) $< -o $@ $(LIBS)

psmf2vtu : psmf2vtu.cpp readwrite.hpp
	$(CPLUSPLUS) $(CPPFLAGS) $< -o $@ 

test:
	(cd tests ; $(MAKE))

clean :
	rm -f $(TARGETS)

