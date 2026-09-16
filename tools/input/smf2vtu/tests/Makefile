#######################################
# MAKEFILE for this directory         #
#######################################

include $(STATFEMROOT)/statFEM_definitions.mk

TESTS = ring32 square

# the rules
#
# if nothing is given, make the target
all: $(TESTS)

ring32:
	../smf2vtu ring32.smf

square:
	../smf2vtu square.smf

# clean
clean:
	$(RM) *.vtu 

