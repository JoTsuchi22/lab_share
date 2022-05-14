00:
	echo which make ?

01:	mkdir_temp cuda_IIM_PCG_v1.x
	mv ./cuda_IIM_PCG_v1.x ./temp/cuda_IIM_PCG_v1.x
	mv ./cuda_IIM_PCG_v1.exp ./temp/uda_IIM_PCG_v1.exp
	mv ./cuda_IIM_PCG_v1.lib ./temp/cuda_IIM_PCG_v1.lib

cuda_IIM_PCG_v1.x:	cuda_IIM_PCG_main_v1.cu cuda_IIM_PCG_sub_v1.cpp s_IGA_header.h s_IGA_main.h s_IGA_sub.h
	nvcc -O1 -Xcompiler "/wd 4819" -o cuda_IIM_PCG_v1.x cuda_IIM_PCG_main_v1.cu cuda_IIM_PCG_sub_v1.cpp

02:	mkdir_temp cpp_IIM_PCG_v1.x

cpp_IIM_PCG_v1.x:	cuda_IIM_PCG_main_v1.cpp cuda_IIM_PCG_sub_v1.cpp s_IGA_header.h s_IGA_main.h s_IGA_sub.h
	g++ -o cpp_IIM_PCG_v1.x cuda_IIM_PCG_main_v1.cpp cuda_IIM_PCG_sub_v1.cpp -O1 -g

mkdir_temp:
	mkdir -p temp

run:
	./temp/cuda_IIM_PCG_v1.x