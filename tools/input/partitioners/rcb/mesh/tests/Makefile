#######################################
# MAKEFILE for this directory         #
#######################################

PARTITIONER = $(STATFEMROOT)/tools/input/partitioners/rcb/mesh/partitioner
PSMF2VTU    = $(STATFEMROOT)/tools/input/smf2vtu/psmf2vtu
SMF2VTU     = $(STATFEMROOT)/tools/input/smf2vtu/smf2vtu

PARTITIONS  = 3
OVERLAP     = 1

SMFFILES    = $(wildcard *.smf)
BASENAMES   = $(SMFFILES:.smf=)
VTUFILES    = $(SMFFILES:.smf=.vtu)

#######################################
# the rules
#
all: test vtu

test: $(BASENAMES)

vtu: $(VTUFILES)

$(BASENAMES):
	$(PARTITIONER) $(PARTITIONS) $(OVERLAP) $@

$(VTUFILES):
	$(PSMF2VTU) $(PARTITIONS) $(@:.vtu=)
	$(SMF2VTU) $(@:.vtu=).smf

# clean
clean:
	-rm -f *psmf *vtu *pvd core

