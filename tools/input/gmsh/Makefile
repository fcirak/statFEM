include $(STATFEMROOT)/statFEM_definitions.mk

# determine modus of compilation
DEBUG=YES

# includes
INCLUDES  = -I$(STATFEMROOT) 
INCLUDES += -I$(BOOST_DIR)
INCLUDES += -I$(EIGEN_DIR)

ifeq ($(DEBUG),YES)
CPPFLAGS += $(DEBDEFINES) $(DEBCPPFLAGS)
LDFLAGS += $(DEBLDFLAGS)
else
CPPFLAGS += $(RELDEFINES) $(RELCPPFLAGS)
LDFLAGS += $(RELLDFLAGS)
endif

# common flags
CPPFLAGS += $(INCLUDES)

all : gmsh2smf

clean :
	rm -f gmsh2smf

rebuild : clean all

gmsh2smf : gmsh2smf.cpp
	$(CPLUSPLUS)  $(CPPFLAGS) gmsh2smf.cpp -o gmsh2smf $(LIBS)

msh :
	gmsh example.geo -2 -order 1
