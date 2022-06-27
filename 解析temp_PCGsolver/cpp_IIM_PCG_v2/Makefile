PROGRAM_NAME	= cpp_IIM_PCG
VERSION			= v2
CC				= g++
CFLAGS			= -g -lm -Wall -Wextra
CFLAGS_fast		= -O1 -lm -Wall -Wextra
CFLAGS_so_fast	= -O3 -lm -Wall -Wextra

compile:
	cd ./src && $(CC) -o $(PROGRAM_NAME)_$(VERSION).x $(PROGRAM_NAME)_main_$(VERSION).cpp $(PROGRAM_NAME)_sub_$(VERSION).cpp $(CFLAGS)

compile_fast:
	cd ./src && $(CC) -o $(PROGRAM_NAME)_$(VERSION).x $(PROGRAM_NAME)_main_$(VERSION).cpp $(PROGRAM_NAME)_sub_$(VERSION).cpp $(CFLAGS_fast)

compile_so_fast:
	cd ./src && $(CC) -o $(PROGRAM_NAME)_$(VERSION).x $(PROGRAM_NAME)_main_$(VERSION).cpp $(PROGRAM_NAME)_sub_$(VERSION).cpp $(CFLAGS_so_fast)

# for debug
run0:	mkdir_output
	cd ./analysis/output && ../../src/$(PROGRAM_NAME)_$(VERSION).x ../input_glo_3_10x10.txt ../input_loc_3_10x10.txt

run1:	mkdir_output
	cd ./analysis/output && ../../src/$(PROGRAM_NAME)_$(VERSION).x ../input_glo_3_20x20.txt ../input_loc_3_20x20.txt

run2:	mkdir_output
	cd ./analysis/output && ../../src/$(PROGRAM_NAME)_$(VERSION).x ../input_glo_3_30x30.txt ../input_loc_3_30x30.txt

run_2D:	mkdir_output
	cd ./analysis/output && ../../src/$(PROGRAM_NAME)_$(VERSION).x ../input_2D/input.txt

run_3D:	mkdir_output
	cd ./analysis/output && ../../src/$(PROGRAM_NAME)_$(VERSION).x ../input_3D/input.txt

run_crack:	mkdir_output
	cd ./analysis/output && ../../src/$(PROGRAM_NAME)_$(VERSION).x ../input_crack/input_glo.txt ../input_crack/input_loc.txt

run_2nd_order:	mkdir_output
	cd ./analysis/output && ../../src/$(PROGRAM_NAME)_$(VERSION).x ../input_glo_2_10x10.txt

run_gomimesh:	mkdir_output
	cd ./analysis/output && ../../src/$(PROGRAM_NAME)_$(VERSION).x ../input_glo_3_30x30.txt ../input_loc_3_5x5.txt

# NURBS viewer, 0: local, 1: global
view0:
	cd ./analysis/output && ../../src/NURBS_viewer/NURBS_view.exe ./input_local.txt ./overlay_view.dat ./Displacement_loc.dat > temp.dat

view1:
	cd ./analysis/output && ../../src/NURBS_viewer/NURBS_view.exe ./input_for_NURBS.txt ./view.dat ./Displacement.dat > temp.dat

mkdir_output:
	cd ./analysis && mkdir -p output

# test run
testrun0:	compile_fast run0

testrun1:	compile_so_fast run1 view0

testrun2:	compile_so_fast run2 view0

testrun_2D: compile_so_fast run_2D view1

testrun_3D: compile_so_fast run_3D

testrun_crack: compile_so_fast run_crack view0

testrun_2nd_order: compile_fast run_2nd_order

testrun_gomimesh:	compile_so_fast run_gomimesh view1