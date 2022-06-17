PROGRAM_NAME	= cpp_IIM_PCG
VERSION			= v2
CC				= g++
CFLAGS			= -g -lm -Wall -Wextra
CFLAGS_fast		= -O1 -lm -Wall -Wextra

compile:
	cd ./src && $(CC) -o $(PROGRAM_NAME)_$(VERSION).x $(PROGRAM_NAME)_main_$(VERSION).cpp $(PROGRAM_NAME)_sub_$(VERSION).cpp $(CFLAGS)

compile_fast:
	cd ./src && $(CC) -o $(PROGRAM_NAME)_$(VERSION).x $(PROGRAM_NAME)_main_$(VERSION).cpp $(PROGRAM_NAME)_sub_$(VERSION).cpp $(CFLAGS_fast)

# for debug
run0:	mkdir_output
	cd ./analysis/output && ../../src/$(PROGRAM_NAME)_$(VERSION).x ../input_glo_3_10x10.txt ../input_loc_3_10x10.txt > result.dat

run1:	mkdir_output
	cd ./analysis/output && ../../src/$(PROGRAM_NAME)_$(VERSION).x ../input_glo_3_20x20.txt ../input_loc_3_20x20.txt

run2:	mkdir_output
	cd ./analysis/output && ../../src/$(PROGRAM_NAME)_$(VERSION).x ../input_glo_3_30x30.txt ../input_loc_3_30x30.txt

run_2D:	mkdir_output
	cd ./analysis/output && ../../src/$(PROGRAM_NAME)_$(VERSION).x ../input_2D/input.txt

run_3D:	mkdir_output
	cd ./analysis/output && ../../src/$(PROGRAM_NAME)_$(VERSION).x ../input_3D/input.txt

# NURBS viewer, 0: local, 1: global
view0:
	./src/NURBS_viewer/NURBS_view.exe ./analysis/output/input_local.txt ./analysis/output/overlay_view.dat ./analysis/output/Displacement_loc.dat > temp.dat

view1:
	./src/NURBS_viewer/NURBS_view.exe ./analysis/output/input_for_NURBS.txt ./analysis/output/view.dat ./analysis/output/Displacement.dat > temp.dat

mkdir_output:
	cd ./analysis && mkdir -p output

# test run
testrun0:	compile_fast run0 view0

testrun1:	compile run1 view0

testrun_2D: compile_fast run_2D view1

testrun_3D: compile run_3D view0