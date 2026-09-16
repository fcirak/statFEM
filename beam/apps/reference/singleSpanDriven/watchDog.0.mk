# makefile to perform watchdog testing
#
# Make sure the target names remain what they are, i.e.
#    'clean' to remove objects, binaries, etc
#    'build' to actually create an execetable
#    'run'   to execute the latter
#    'test'  to test a result versus a reference
# Also make sure to provide the '@echo "WatchDog ..."' statements
# appear at the bottom of every make target. They are
# scanned by the tools/watchdog scripts to automate
# testing.
#
# Note that the STATFEMROOT variable must not be set within
# this Makefile, otherwise it looses its portability.

# diff ASCII VTU subject to tolerance
VTUDIFF = $(STATFEMROOT)/tools/watchdog/vtuDiff.py
DATDIFF = $(STATFEMROOT)/tools/watchdog/datDiff.py

# convenience target to run complete test cycle on the fly
all: clean build run test
.PHONY: all

# main targets as used by watchdog automated testing
clean:
	make -f Makefile clean
	rm -rf animation*
	@echo "WatchDog clean successful"

build:
	make -f Makefile
	@echo "WatchDog build successful"

run:
	./singleSpanDriven
	@echo "WatchDog run successful"

test:
	$(VTUDIFF) --new=animation/output.010.vtu --ref=watchDog.0.output.010.vtu --tol=1.e-6
	$(DATDIFF) --new=animation/nodeMonitor.dat --ref=watchDog.0.nodeMonitor.dat --tol=1.e-6
	@echo "WatchDog test successful"


