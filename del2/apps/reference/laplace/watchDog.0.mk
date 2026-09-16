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

# convenience
VTUDIFF = $(STATFEMROOT)/tools/watchdog/vtuDiff.py
L2ERR   = 0.000996162
ERRTOL  = 1.e-8
ERRFILE = error.dat

# convenience target to run complete test cycle on the fly
all: clean build run test
.PHONY: all

# main targets as used by watchdog automated testing
clean:
	make -f Makefile clean
	rm -f laplace.vtu
	@echo "WatchDog clean successful"

build:
	make -f Makefile DEBUG=NO
	@echo "WatchDog build successful"

run:
	./laplace
	@echo "WatchDog run successful"

# First check:  use vtuDiff to compare output vtu-file with reference
# Second check: read L2-error from file and compare with expected value
#
test:
	$(VTUDIFF) --new=laplace.vtu --ref=watchDog.0.mk.laplace.vtu --tol=${ERRTOL}
	awk '/L2-error/ {exit !(sqrt(($$3-${L2ERR})*($$3-${L2ERR})) < ${ERRTOL}) }' ${ERRFILE}
	@echo "WatchDog test successful"

