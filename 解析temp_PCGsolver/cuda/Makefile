s_IGA_cuda_cpp.x:	s_IGA_cuda_cpp.cpp
					g++ -o s_IGA_cuda_cpp.x s_IGA_cuda_cpp.cpp -O1

01:	s_IGA_cuda.cu
	nvcc -Xcompiler -O0 s_IGA_cuda.cu -o s_IGA_cuda.x -Xcompiler "/wd 4819"

a:	a.cu
	nvcc -O0 a.cu -o a.x -Xcompiler "/wd 4819" && ./a.x

b:	a.cu
	nvcc -O2 cuda_test_IMM.cu -o cuda_test_IMM.x -Xcompiler "/wd 4819"

c:	cuda_test_IMM.cu
	nvcc -O0 cuda_test_IMM.cu -o cuda_test_IMM.x -Xcompiler "/wd 4819"

d:	cuda_test_IMM.cpp
	nvcc -O0 cuda_test_IMM.cpp -o cuda_test_IMM.x -Xcompiler "/wd 4819"