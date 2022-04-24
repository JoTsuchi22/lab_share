compile:	mkdir_temp cuda_IIM_PCG_v1.x
	mv ./cuda_IIM_PCG_v1.x ./temp/cuda_IIM_PCG_v1.x
	mv ./cuda_IIM_PCG_v1.exp ./temp/uda_IIM_PCG_v1.exp
	mv ./cuda_IIM_PCG_v1.lib ./temp/cuda_IIM_PCG_v1.lib

cuda_IIM_PCG_v1.x:	cuda_IIM_PCG_main_v1 cuda_IIM_PCG_sub_v1 myheader.h
	nvcc -O1 -Xcompiler "/wd 4819" -o cuda_IIM_PCG_v1.x cuda_IIM_PCG_main_v1.o cuda_IIM_PCG_sub_v1.obj
	# nvcc -O1 -Xcompiler "/wd 4819" -o cuda_IIM_PCG_v1.x cuda_IIM_PCG_main_v1.cpp cuda_IIM_PCG_sub_v1.cu

cuda_IIM_PCG_main_v1:	cuda_IIM_PCG_main_v1.cpp
	g++ -c cuda_IIM_PCG_main_v1.cpp

cuda_IIM_PCG_sub_v1:	cuda_IIM_PCG_sub_v1.cu
	nvcc -c cuda_IIM_PCG_sub_v1.cu

mkdir_temp:
	mkdir -p temp

run:
	./temp/cuda_IIM_PCG_v1.x