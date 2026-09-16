# makefile to perform watchdog testing
#
# Make sure the target names remain what they are, i.e.
#    'clean' to remove objects, binaries, etc
#    'build' to actually create an execetable
#    'run'   to execute the latter
#    'test'  to test a result versus a reference
# Also make sure to provide the '@echo "WatchDog ..."' statements
# appear at the bottom of every make target. They are
# looked for by the tools/watchdog scripts to automate
# testing.
#
# Note that the OPENFTLROOT variable must not be set within
# this Makefile, otherwise it looses its portability.

# setup variables to run savely on iwulf
ifeq ($(shell hostname -s),master)
  SOURCE_MY_ENV = source /opt/bin/profile-clean.sh &&  \
	          export OPENFTLCOMPENV=INTEL12 &&  \
                  source /opt/bin/openftlcompenv-setup.sh
else
  SOURCE_MY_ENV = source /dev/null
endif

# convenience
VTUDIFF = $(OPENFTLROOT)/tools/watchdog/vtuDiff.py

# convenience target to run complete test cycle on the fly
all: clean build run test
.PHONY: all

# main targets as used by watchdog automated testing
clean:
	bash -c "${SOURCE_MY_ENV} && make -f Makefile clean"
	rm -f planeSheet/planeSheet.000-0020.vtu
	rm -f planeSheet/planeSheet.001-0020.vtu
	rm -f planeSheet/planeSheet.002-0020.vtu
	rm -f planeSheet/planeSheet.003-0020.vtu
	@echo "WatchDog clean successful"

build:
	bash -c "${SOURCE_MY_ENV} && make -f Makefile"
	@echo "WatchDog build successful"

run:
	bash -c "${SOURCE_MY_ENV} && cd planeSheet && mpirun -np 4 ../parallelDriven && cd .."
	@echo "WatchDog run successful"

test:
	$(VTUDIFF) --new=planeSheet/planeSheet.000-0020.vtu --ref=planeSheet/watchDog.0.planeSheet.000-0020.vtu --comp=rel --tol=1.e-5
	$(VTUDIFF) --new=planeSheet/planeSheet.001-0020.vtu --ref=planeSheet/watchDog.0.planeSheet.001-0020.vtu --comp=rel --tol=1.e-5
	$(VTUDIFF) --new=planeSheet/planeSheet.002-0020.vtu --ref=planeSheet/watchDog.0.planeSheet.002-0020.vtu --comp=rel --tol=1.e-5
	$(VTUDIFF) --new=planeSheet/planeSheet.003-0020.vtu --ref=planeSheet/watchDog.0.planeSheet.003-0020.vtu --comp=rel --tol=1.e-5
	@echo "WatchDog test successful"


