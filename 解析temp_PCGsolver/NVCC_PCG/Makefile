PROGRAM_NAME	= S-IGA_PCG
VERSION			= v1
CC				= nvcc
CFLAGS			= -Xcompiler "/wd 4819"
# -lcufft

#compile
compile_so_fast:
	cd ./src && $(CC) -O3 $(PROGRAM_NAME)_main_$(VERSION).cu $(PROGRAM_NAME)_sub_$(VERSION).cu $(PROGRAM_NAME)_kernel_$(VERSION).cu -o $(PROGRAM_NAME)_$(VERSION).x $(CFLAGS)

# for debug
run:	rm	mkdir_output
	cd ./analysis/output && ../../src/$(PROGRAM_NAME)_$(VERSION).x ../input/input_glo_3_20x20.txt ../input/input_loc_3_20x20.txt

run_VS:	rm	mkdir_output
	cd ./analysis/output && ../../src/VS_PCG.exe ../input/input_glo_3_20x20.txt ../input/input_loc_3_20x20.txt

test:	rm	mkdir_output
	cd ./analysis/output && ../../src/$(PROGRAM_NAME)_$(VERSION).x ../input/input.txt

# mkdir
mkdir_output:
	cd ./analysis && mkdir -p output

# rm
rm:
	cd ./analysis && rm -r output -f