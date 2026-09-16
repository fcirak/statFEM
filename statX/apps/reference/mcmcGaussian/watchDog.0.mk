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
DATDIFF         = $(STATFEMROOT)/tools/watchdog/datDiff.py
ERRTOL          = 1.e-8
ACCEPTRATIO     = 0.245920000
ACCEPTRATIOFILE = acceptRatio.dat
SAMPLEFILE      = sample.dat
SAMPLEREF       = watchDog.0.mk.sample.dat
LIKELIHOODFILE  = likelihood.dat
LIKELIHOODREF   = watchDog.0.mk.likelihood.dat

# convenience target to run complete test cycle on the fly
all: clean build run test
.PHONY: all

# main targets as used by watchdog automated testing
clean:
	make -f Makefile clean
	rm -f ${ACCEPTRATIOFILE} ${SAMPLEFILE} ${LIKELIHOODFILE}
	@echo "WatchDog clean successful"

build:
	make -f Makefile DEBUG=NO
	@echo "WatchDog build successful"

run:
	./mcmcGaussian
	@echo "WatchDog run successful"

# Check 1: use datDiff.py to compare samples with reference
# Check 2: use datDiff.py to compare samples' likelihood with reference
# Check 3: read accept ratio from file and compare with reference value
test:
	$(DATDIFF) --new=${SAMPLEFILE} --ref=${SAMPLEREF} --tol=${ERRTOL}
	$(DATDIFF) --new=${LIKELIHOODFILE} --ref=${LIKELIHOODREF} --tol=${ERRTOL}
	awk '/likelihood/ {exit !(sqrt(($$3-${ACCEPTRATIO})*($$3-${ACCEPTRATIO})) < ${ERRTOL}) }' ${ACCEPTRATIOFILE}
	@echo "WatchDog test successful"
	
