CC				= g++
CFLAGS			= -O3 -lm -Wall -Wextra

MC3D:	./src/s_IGA_MC3D.cpp
	cd ./src && $(CC) -o s_IGA_MC3D.x s_IGA_MC3D.cpp $(CFLAGS)

# for debug
run_a:	MC3D
	cd ./textfile && ../src/s_IGA_MC3D.x ./test_boundary.txt ./test_MC3D_00_edited.txt ./test_MC3D_01_edited.txt

run_b:	MC3D
	cd ./textfile && ../src/s_IGA_MC3D.x ./test_boundary.txt ./test_MC3D_00_edited.txt ./test_MC3D_01_edited.txt