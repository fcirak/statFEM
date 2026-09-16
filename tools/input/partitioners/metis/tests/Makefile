#######################################
# MAKEFILE for this directory         #
#######################################

include $(STATFEMROOT)/statFEM_definitions.mk

NSUB = 8 
OVERLAP = 1
PROBLEM = pinchedHemisphere

# the rules
#
# if nothing is given, make the target
all: partition conversion join join_conversion

# run partitioner
partition: 
	../partitioner $(NSUB) $(OVERLAP) $(PROBLEM)

# run converter
conversion: 
	$(STATFEMROOT)/tools/input/smf2vtu/psmf2vtu $(NSUB) $(PROBLEM)

# run merger
join: 
	../psmf2smf $(NSUB) $(PROBLEM)

# run convertor
join_conversion: 
	$(STATFEMROOT)/tools/input/smf2vtu/smf2vtu $(PROBLEM)_joined.smf

# clean
clean:
	$(RM) $(PROBLEM).*.psmf $(PROBLEM).*.vtu $(PROBLEM).pvd $(TARGET) $(PROBLEM)_joined.smf

