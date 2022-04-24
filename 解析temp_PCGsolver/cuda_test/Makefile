compile:	mkdir_temp cuda_IIM_PCG_v1.x
	mv ./cuda_IIM_PCG_v1.x ./temp/cuda_IIM_PCG_v1.x
	mv ./cuda_IIM_PCG_v1.exp ./temp/uda_IIM_PCG_v1.exp
	mv ./cuda_IIM_PCG_v1.lib ./temp/cuda_IIM_PCG_v1.lib

cuda_IIM_PCG_v1.x:	cuda_IIM_PCG_main_v1.cu cuda_IIM_PCG_sub_v1.cpp myheader.h
	nvcc -O2 -Xcompiler "/wd 4819" -o cuda_IIM_PCG_v1.x cuda_IIM_PCG_main_v1.cu cuda_IIM_PCG_sub_v1.cpp

mkdir_temp:
	mkdir -p temp

run:
	./temp/cuda_IIM_PCG_v1.x