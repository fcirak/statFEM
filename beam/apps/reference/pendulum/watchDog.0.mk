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
	./pendulum
	@echo "WatchDog run successful"

test:
	$(VTUDIFF) --new=animation/line8.0005.vtu --ref=watchDog.0.line8.0005.vtu --tol=1.e-6
	@echo "WatchDog test successful"


