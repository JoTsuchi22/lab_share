PROGRAM_NAME	= cpp_IIM_PCG
VERSION			= v3
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
	cd ./analysis/output && ../../src/$(PROGRAM_NAME)_$(VERSION).x ../input_2D_IGA/input.txt

run_3D:	mkdir_output
	cd ./analysis/output && ../../src/$(PROGRAM_NAME)_$(VERSION).x ../input_3D_IGA/input.txt

run_3D_S-IGA:	mkdir_output
	cd ./analysis/output && ../../src/$(PROGRAM_NAME)_$(VERSION).x ../input_3D_S-IGA/input_glo.txt ../input_3D_S-IGA/input_loc.txt

run_3D_S-IGA_2:	mkdir_output
	cd ./analysis/output && ../../src/$(PROGRAM_NAME)_$(VERSION).x ../input_3D_S-IGA_2/input_glo.txt ../input_3D_S-IGA_2/input_loc.txt

run_3D_S-IGA_3:	mkdir_output
	cd ./analysis/output && ../../src/$(PROGRAM_NAME)_$(VERSION).x ../input_3D_S-IGA_3/input.txt

run_crack:	mkdir_output
	cd ./analysis/output && ../../src/$(PROGRAM_NAME)_$(VERSION).x ../input_crack/input_glo.txt ../input_crack/input_loc.txt

run_2nd_order:	mkdir_output
	cd ./analysis/output && ../../src/$(PROGRAM_NAME)_$(VERSION).x ../input_glo_2_10x10.txt

run_gomimesh:	mkdir_output
	cd ./analysis/output && ../../src/$(PROGRAM_NAME)_$(VERSION).x ../input_glo_3_30x30.txt ../input_loc_3_5x5.txt

# mkdir
mkdir_output:
	cd ./analysis && mkdir -p output

# test run
testrun0:	compile_fast run0

testrun1:	compile_so_fast run1

testrun2:	compile_so_fast run2

testrun_2D: compile_so_fast run_2D

testrun_3D: compile_so_fast run_3D

testrun_3D_S-IGA: compile_so_fast run_3D_S-IGA

testrun_3D_S-IGA_2: compile_so_fast run_3D_S-IGA_2

testrun_3D_S-IGA_3: compile_so_fast run_3D_S-IGA_3

testrun_crack: compile_so_fast run_crack

testrun_2nd_order: compile_fast run_2nd_order

testrun_gomimesh:	compile_so_fast run_gomimesh