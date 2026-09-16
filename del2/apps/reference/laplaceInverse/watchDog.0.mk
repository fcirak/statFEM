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
TOL          = 1.e-6

# hyperparameter inferred values
HPARAMFILE   = hParMean.dat

# hyperparameter reference values
W0       = 0.770050
W1       = 1.060490
W2       = 0.847863
W3       = 0.319244
W4       = 0.705646
EPS      = 0.010782

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
	./laplaceInverse
	@echo "WatchDog run successful"

# Check:  read log-likelihood and hyperparamaters from file and compare with the expected values
#
test:
	awk '/gaussianProcess.w0/ {exit !(sqrt(($$2-${W0})*($$2-${W0})) < ${TOL}) }' ${HPARAMFILE}
	awk '/gaussianProcess.w1/ {exit !(sqrt(($$2-${W1})*($$2-${W1})) < ${TOL}) }' ${HPARAMFILE}
	awk '/gaussianProcess.w2/ {exit !(sqrt(($$2-${W2})*($$2-${W2})) < ${TOL}) }' ${HPARAMFILE}
	awk '/gaussianProcess.w3/ {exit !(sqrt(($$2-${W3})*($$2-${W3})) < ${TOL}) }' ${HPARAMFILE}
	awk '/gaussianProcess.w4/ {exit !(sqrt(($$2-${W4})*($$2-${W4})) < ${TOL}) }' ${HPARAMFILE}
	awk '/mismatch.eps/ {exit !(sqrt(($$2-${EPS})*($$2-${EPS})) < ${TOL}) }' ${HPARAMFILE}
	@echo "WatchDog test successful"

