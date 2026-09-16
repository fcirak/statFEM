# my configuration
include $(STATFEMROOT)/statFEM_definitions.mk

# determine modus of compilation {YES, NO}
DEBUG = 

# includes
INCLUDES  = -I$(STATFEMROOT) 
INCLUDES += -I$(BOOST_DIR)

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

# libraries
LIBS = 

# makefile build targets
TARGETS = tg2tgl tgl2ptg

all : $(TARGETS)

rebuild : clean all

tg2tgl : tg2tgl.cpp readwrite.hpp
	$(CPLUSPLUS) $(CPPFLAGS) $< -o $@ $(LIBS)

tgl2ptg : tgl2ptg.cpp readwrite.hpp
	$(CPLUSPLUS) $(CPPFLAGS) $< -o $@ $(LIBS)

#test:
#	(cd tests ; $(MAKE))

clean :
	rm -f $(TARGETS)

