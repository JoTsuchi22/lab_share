CC				= g++
CFLAGS			= -O3 -lm -Wall -Wextra

MI3D:	./src/S_IGA_MI3D.cpp	./src/header_MI3D.h
	cd ./src && $(CC) -o S_IGA_MI3D.x S_IGA_MI3D.cpp $(CFLAGS)

# for debug
run0:	MI3D
	cd ./textfile && ../src/S_IGA_MI3D.x ../textfile/test_MI3D_00.txt ../textfile/test_MI3D_01.txt ../textfile/test_MI3D_02.txt

run1:	MI3D
	cd ./textfile && ../src/S_IGA_MI3D.x ../textfile/test_MC3D_00.txt ../textfile/test_MC3D_01.txt ../textfile/test_MC3D_02.txt ../textfile/test_MC3D_03.txt