# 
include $(STATFEMROOT)/statFEM_definitions.mk

# includes
INCLUDES = -I$(UTILSPATH)

# libraries
LIBS = -L$(UTILSPATH) -lutils

# common flags
CPPFLAGS = $(INCLUDES)

all : genesis2smf

rebuild : clean all

clean :
	rm -f genesis2smf

genesis2smf : main.cpp
	g++ $(CPPFLAGS) main.cpp -o genesis2smf $(LIBS)
