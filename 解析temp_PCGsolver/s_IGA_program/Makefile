all_program: 02

01:	IIM_PCG_v2.c
	gcc -o IIM_PCG_v2.x IIM_PCG_v2.c -O1 -lm -Wall -Wextra

02:	IIM_PCG_v2_cu.cu
	nvcc -O1 IIM_PCG_v2_cu.cu -o IIM_PCG_v2_cu.x -Xcompiler "/wd 4819"

03:	IIM_PCG_v2_cu.cpp
	g++ -o IIM_PCG_v2_cu.x IIM_PCG_v2_cu.cpp -O0 -lm -Wall -Wextra