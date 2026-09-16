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

# tolerance 
TOL          = 1.e-8

# hyperparameter reference values
RHO          = 0.804412
SIGMAD       = 0.0191724
ELLD         = 0.0740285
# hyperparameter inferred values
HPARAMFILE   = inferredHyperParam.dat

# testing the inferred FE against its reference values
DATDIFF      = $(STATFEMROOT)/tools/watchdog/datDiff.py
FEINFERREF   = watchDog.0.mk.inferredFEcondData.dat
FEINFERFILE  = inferredFEcondData.dat

# convenience target to run complete test cycle on the fly
all: clean build run test
.PHONY: all

# main targets as used by watchdog automated testing
clean:
	make -f Makefile clean
	@echo "WatchDog clean successful"

build:
	make -f Makefile DEBUG=YES
	@echo "WatchDog build successful"

run:
	./laplaceRHSinfer
	@echo "WatchDog run successful"

# Check:  read log-likelihood and hyperparamaters from file and compare with the expected values
test:
	awk '/rho/ {exit !(sqrt(($$2-${RHO})*($$2-${RHO})) < ${TOL}) }' ${HPARAMFILE}
	awk '/missmatch.sigma/ {exit !(sqrt(($$2-${SIGMAD})*($$2-${SIGMAD})) < ${TOL}) }' ${HPARAMFILE}
	awk '/missmatch.length/ {exit !(sqrt(($$2-${ELLD})*($$2-${ELLD})) < ${TOL}) }' ${HPARAMFILE}
	$(DATDIFF) --new=${FEINFERFILE} --ref=${FEINFERREF} --tol=${TOL}
	
	@echo "WatchDog test successful"

