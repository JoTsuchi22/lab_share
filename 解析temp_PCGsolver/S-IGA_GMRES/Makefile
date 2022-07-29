PROGRAM_NAME	= S-IGA_PCG
VERSION			= v1
CC				= g++
CFLAGS			= -g -lm -Wall -Wextra
CFLAGS_fast		= -O1 -lm -Wall -Wextra
CFLAGS_so_fast	= -O3 -lm -Wall -Wextra

#compile
compile_so_fast:
	cd ./src && $(CC) -o $(PROGRAM_NAME)_$(VERSION).x $(PROGRAM_NAME)_main_$(VERSION).cpp $(PROGRAM_NAME)_sub_$(VERSION).cpp $(CFLAGS_so_fast)

compile_debug:
	cd ./src && $(CC) -o $(PROGRAM_NAME)_$(VERSION).x $(PROGRAM_NAME)_main_$(VERSION).cpp $(PROGRAM_NAME)_sub_$(VERSION).cpp $(CFLAGS)

compile_fast:
	cd ./src && $(CC) -o $(PROGRAM_NAME)_$(VERSION).x $(PROGRAM_NAME)_main_$(VERSION).cpp $(PROGRAM_NAME)_sub_$(VERSION).cpp $(CFLAGS_fast)

# for debug
run:	rm	mkdir_output
	cd ./analysis/output && ../../src/$(PROGRAM_NAME)_$(VERSION).x ../input/input_glo_3_20x20.txt ../input/input_loc_3_20x20.txt

run_ss:	rm	mkdir_output
	cd ./analysis/output && ../../src/$(PROGRAM_NAME)_$(VERSION).x ../input/input_ss.txt

run_VS:	rm	mkdir_output
	cd ./analysis/output && ../../src/VS_GMRES.exe ../input/input_glo_3_20x20.txt ../input/input_loc_3_20x20.txt

test:	rm	mkdir_output
	cd ./analysis/output && ../../src/$(PROGRAM_NAME)_$(VERSION).x ../input/input.txt

# mkdir
mkdir_output:
	cd ./analysis && mkdir -p output

# rm
rm:
	cd ./analysis && rm -r output -f