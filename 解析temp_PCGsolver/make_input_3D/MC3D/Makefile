CC				= g++
CFLAGS			= -O3 -lm -Wall -Wextra

MC3D:	./src/s_IGA_MC3D.cpp
	cd ./src && $(CC) -o s_IGA_MI3D.x s_IGA_MI3D.cpp $(CFLAGS)

# for debug
run0:	MC3D
	cd ./textfile && ../src/s_IGA_MI3D.x ./boundary_test.txt ./test_MI3D_00.txt ./test_MI3D_01.txt