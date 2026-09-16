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
# Note that the OPENFTLROOT variable must not be set within
# this Makefile, otherwise it looses its portability.

# setup variables to run savely on iwulf
ifeq ($(shell hostname -s),master)
    MY_CLEAN = source /opt/bin/profile-clean.sh
    MY_ICC = export OPENFTLCOMPENV=INTEL && source /opt/bin/openftlcompenv-setup.sh
else
    $(error Only available on iwulf)
    MY_CLEAN = /dev/null
    MY_ICC = /dev/null
endif
SOURCE_MY_ENV = "${MY_CLEAN} && ${MY_ICC}"

# convenience
VTUDIFF = $(OPENFTLROOT)/tools/watchdog/vtuDiff.py

# convenience target to run complete test cycle on the fly
all: clean build run test
.PHONY: all

# main targets as used by watchdog automated testing
clean:
	bash -c ${SOURCE_MY_ENV}" && make -f Makefile clean SOLVER=PARDISO"
	rm -rf animation*
	@echo "WatchDog clean successful"

build:
	bash -c ${SOURCE_MY_ENV}" && make -f Makefile SOLVER=PARDISO"
	@echo "WatchDog build successful"

run:
	bash -c ${SOURCE_MY_ENV}" && ./static"
	@echo "WatchDog run successful"

test:
	$(VTUDIFF) --new=animation/planeSheetQQ.smf.0003.vtu --ref=watchDog.0.planeSheetQQ.smf.0003.vtu --tol=1.e-6
	@echo "WatchDog test successful"


