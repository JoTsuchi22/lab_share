PROGRAM_NAME	= cpp_IIM_PCG
VERSION			= v1
CC				= g++
CFLAGS			= -O1 -lm -Wall -Wextra

compile:
	cd ./src && $(CC) -o $(PROGRAM_NAME)_$(VERSION).x $(PROGRAM_NAME)_main_$(VERSION).cpp $(PROGRAM_NAME)_sub_$(VERSION).cpp $(CFLAGS)

# for debug
run0:	mkdir_output
	cd ./analysis/output && ../../src/$(PROGRAM_NAME)_$(VERSION).x ../input_glo_3_10x10.txt ../input_loc_3_10x10.txt

run1:	mkdir_output
	cd ./analysis/output && ../../src/$(PROGRAM_NAME)_$(VERSION).x ../input_glo_3_20x20.txt ../input_loc_3_20x20.txt

run2:	mkdir_output
	cd ./analysis/output && ../../src/$(PROGRAM_NAME)_$(VERSION).x ../input_glo_3_30x30.txt ../input_loc_3_30x30.txt

# NURBS viewer, 0: local, 1: global
view0:
	./src/NURBS_viewer/NURBS_view.exe ./analysis/output/input_local.txt ./analysis/output/overlay_view.dat ./analysis/output/Displacement_loc.dat

view1:
	./src/NURBS_viewer/NURBS_view.exe ./analysis/output/input_for_NURBS.txt ./analysis/output/view.dat ./analysis/output/Displacement.dat

mkdir_output:
	cd ./analysis && mkdir -p output

# test run
testrun0:	compile run0 view0

testrun1:	compile run1 view0