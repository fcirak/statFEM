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
DATDIFF        = $(STATFEMROOT)/tools/watchdog/datDiff.py
ERRTOL         = 1.e-8
LIKELIHOOD     = -4.837171259
LIKELIHOODFILE = likelihood.dat
PREDMEANFILE   = predictionMean.dat
PREDMEANREF    = watchDog.0.mk.predictionMean.dat
PREDVARFILE    = predictionVar.dat
PREDVARREF     = watchDog.0.mk.predictionVar.dat

# convenience target to run complete test cycle on the fly
all: clean build run test
.PHONY: all

# main targets as used by watchdog automated testing
clean:
	make -f Makefile clean
	rm -f ${LIKELIHOODFILE} ${PREDMEANFILE} ${PREDVARFILE}
	@echo "WatchDog clean successful"

build:
	make -f Makefile DEBUG=NO
	@echo "WatchDog build successful"

run:
	./gaussianProcess
	@echo "WatchDog run successful"

# Check 1: use datDiff.py to compare output prediction mean with reference
# Check 2: use datDiff.py to compare output prediction variance with reference
# Check 3: read likelihood from file and compare with reference value
test:
	$(DATDIFF) --new=${PREDMEANFILE} --ref=${PREDMEANREF} --tol=${ERRTOL}
	$(DATDIFF) --new=${PREDVARFILE} --ref=${PREDVARREF} --tol=${ERRTOL}
	awk '/likelihood/ {exit !(sqrt(($$3-${LIKELIHOOD})*($$3-${LIKELIHOOD})) < ${ERRTOL}) }' ${LIKELIHOODFILE}
	@echo "WatchDog test successful"
