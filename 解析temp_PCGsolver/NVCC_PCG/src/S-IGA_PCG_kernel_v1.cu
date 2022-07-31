#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// cuda
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

// header
#include "S_IGA_header.h"
#include "S_IGA_sub.h"

using namespace std;

// kernel
int RoundUpFunc(int a, int b)
{
	return ((a + b - 1) / b);
}


void Preprocessing_IGA_pararel(information *info)
{
	// GPU側でメモリを割り当てる
	cudaMalloc((void**)&dev_a, N * N * sizeof(double));
	cudaMalloc((void**)&dev_b, N * N * sizeof(double));
	cudaMalloc((void**)&dev_c, N * N * sizeof(double));

    // 配列aと配列bをGPUにコピーする
	cudaMemcpy(dev_a, a, N * N * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_b, b, N * N * sizeof(double), cudaMemcpyHostToDevice);

	int size_of_parallel_process = N * N;
	int block_dim = 128;
	int grid_dim = RoundUpFunc(size_of_parallel_process, block_dim);					// グリッド内のブロック数 (1次元，各ブロックで同じ値，切り上げ)
	pararel <<<grid_dim, block_dim>>>(dev_a, dev_b, dev_c, N, N, N * N);	// カーネル スレッドあたりの処理・演算が小さくなるように書く (できるだけforループ使わない)

	// 配列cをGPUからCPUにコピーする
	cudaMemcpy(c, dev_c, N * N * sizeof(double), cudaMemcpyDeviceToHost);

    // GPU側で割り当てたメモリを開放する
	cudaFree(dev_a);
	cudaFree(dev_b);
	cudaFree(dev_c);
}

__global__ void pararel()
{

}