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

# diff ASCII VTU subject to tolerance
VTUDIFF = $(OPENFTLROOT)/tools/watchdog/vtuDiff.py

# convenience target to run complete test cycle on the fly
all: clean build run test
.PHONY: all

# main targets as used by watchdog automated testing
clean:
	make -f Makefile clean
	rm -f 2DGrid_circleMesh_BF.vts
	rm -f 2DGrid_circleMesh_EBB.vts
	rm -f 2DGrid_lineMesh_BF.vts
	rm -f 2DGrid_lineMesh_EBB.vts
	rm -f 2DGrid_zigzagMesh_BF.vts
	rm -f 2DGrid_zigzagMesh_EBB.vts
	rm -f sphereGrid_sphereMesh1_BF.vts
	rm -f sphereGrid_sphereMesh1_EBB.vts
	rm -f sphereGrid_sphereMesh2_BF.vts
	rm -f sphereGrid_sphereMesh2_EBB.vts
	rm -f heartGrid1_heartMesh_BF.vts
	rm -f heartGrid1_heartMesh_EBB.vts
	rm -f heartGrid2_heartMesh_BF.vts
	rm -f heartGrid2_heartMesh_EBB.vts
	rm -f heartGrid3_heartMesh_BF.vts
	rm -f heartGrid3_heartMesh_EBB.vts
	@echo "WatchDog clean successful"

build:
	make -f Makefile levelSetTest DEBUG=NO
	@echo "WatchDog build successful"

run:
	./levelSetTest 2DGrid.sgf circleMesh.smf 1 0 0.2
	./levelSetTest 2DGrid.sgf circleMesh.smf 0 0 0.2
	./levelSetTest 2DGrid.sgf lineMesh.smf   1 0 0.2
	./levelSetTest 2DGrid.sgf lineMesh.smf   0 0 0.2
	./levelSetTest 2DGrid.sgf zigzagMesh.smf 1 1 0.2
	./levelSetTest 2DGrid.sgf zigzagMesh.smf 0 1 0.2
	./levelSetTest sphereGrid.sgf sphereMesh1.smf 1 0 0.15
	./levelSetTest sphereGrid.sgf sphereMesh1.smf 0 0 0.15
	./levelSetTest sphereGrid.sgf sphereMesh2.smf 1 0 0.15
	./levelSetTest sphereGrid.sgf sphereMesh2.smf 0 0 0.15
	./levelSetTest heartGrid1.sgf heartMesh.smf 1 0 5.0
	./levelSetTest heartGrid1.sgf heartMesh.smf 0 0 5.0
	./levelSetTest heartGrid2.sgf heartMesh.smf 1 0 5.0
	./levelSetTest heartGrid2.sgf heartMesh.smf 0 0 5.0
	./levelSetTest heartGrid3.sgf heartMesh.smf 1 0 5.0
	./levelSetTest heartGrid3.sgf heartMesh.smf 0 0 5.0
	@echo "WatchDog run successful"

test:
	$(VTUDIFF) --new=2DGrid_circleMesh_BF.vts  --ref=watchDog.2DGrid_circleMesh_BF.vts  --tol=1.e-6
	$(VTUDIFF) --new=2DGrid_circleMesh_EBB.vts --ref=watchDog.2DGrid_circleMesh_EBB.vts --tol=1.e-6
	$(VTUDIFF) --new=2DGrid_lineMesh_BF.vts    --ref=watchDog.2DGrid_lineMesh_BF.vts    --tol=1.e-6
	$(VTUDIFF) --new=2DGrid_lineMesh_EBB.vts   --ref=watchDog.2DGrid_lineMesh_EBB.vts   --tol=1.e-6
	$(VTUDIFF) --new=2DGrid_zigzagMesh_BF.vts  --ref=watchDog.2DGrid_zigzagMesh_BF.vts  --tol=1.e-6
	$(VTUDIFF) --new=2DGrid_zigzagMesh_EBB.vts --ref=watchDog.2DGrid_zigzagMesh_EBB.vts --tol=1.e-6
	$(VTUDIFF) --new=sphereGrid_sphereMesh1_BF.vts  --ref=watchDog.sphereGrid_sphereMesh1_BF.vts  --tol=1.e-6
	$(VTUDIFF) --new=sphereGrid_sphereMesh1_EBB.vts --ref=watchDog.sphereGrid_sphereMesh1_EBB.vts --tol=1.e-6
	$(VTUDIFF) --new=sphereGrid_sphereMesh2_BF.vts  --ref=watchDog.sphereGrid_sphereMesh2_BF.vts  --tol=1.e-6
	$(VTUDIFF) --new=sphereGrid_sphereMesh2_EBB.vts --ref=watchDog.sphereGrid_sphereMesh2_EBB.vts --tol=1.e-6
	$(VTUDIFF) --new=heartGrid1_heartMesh_BF.vts  --ref=watchDog.heartGrid1_heartMesh_BF.vts  --tol=1.e-6
	$(VTUDIFF) --new=heartGrid1_heartMesh_EBB.vts --ref=watchDog.heartGrid1_heartMesh_EBB.vts --tol=1.e-6
	$(VTUDIFF) --new=heartGrid2_heartMesh_BF.vts  --ref=watchDog.heartGrid2_heartMesh_BF.vts  --tol=1.e-6
	$(VTUDIFF) --new=heartGrid2_heartMesh_EBB.vts --ref=watchDog.heartGrid2_heartMesh_EBB.vts --tol=1.e-6
	$(VTUDIFF) --new=heartGrid3_heartMesh_BF.vts  --ref=watchDog.heartGrid3_heartMesh_BF.vts  --tol=1.e-6
	$(VTUDIFF) --new=heartGrid3_heartMesh_EBB.vts --ref=watchDog.heartGrid3_heartMesh_EBB.vts --tol=1.e-6
	@echo "WatchDog test successful"
