CC				= g++
CFLAGS			= -O1 -lm -Wall -Wextra

MI3D:	src/s_IGA_MI3D.cpp	src/header_MI3D.h
	cd ./src && $(CC) -o s_IGA_MI3D.x s_IGA_MI3D.cpp $(CFLAGS)

MC3D:	src/s_IGA_MC3D.cpp
	cd ./src && $(CC) -o s_IGA_MI3D.x s_IGA_MI3D.cpp $(CFLAGS)

# for debug
run0:	MI3D

run1:	MC3D