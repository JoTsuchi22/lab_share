#include <iostream>
#include <stdio.h>
#include <time.h>
#include <math.h>

using namespace std;

#define rep(i, n) for (i = 0; i < (n); i++)
#define SIZE_OF_MATRIX 1400

void block_calc2(double A[SIZE_OF_MATRIX][SIZE_OF_MATRIX], double B[SIZE_OF_MATRIX][SIZE_OF_MATRIX], double C[SIZE_OF_MATRIX][SIZE_OF_MATRIX])
{
	int i, j, k, ii, jj, kk;
	static const int size = 20;

	for (ii = 0; ii < SIZE_OF_MATRIX; ii += size)
		for (jj = 0; jj < SIZE_OF_MATRIX; jj += size)
			for (kk = 0; kk < SIZE_OF_MATRIX; kk += size)
				for (i = ii; i < ii + size; i += 2)
					for (j = jj; j < jj + size; j += 2)
						for (k = kk; k < kk + size; k++)
						{
							C[i][j]         += A[i][k]     * B[k][j];
							C[i + 1][j]     += A[i + 1][k] * B[k][j];
							C[i][j + 1]     += A[i][k]     * B[k][j + 1];
							C[i + 1][j + 1] += A[i + 1][k] * B[k][j + 1];
						}
	
}

void block_calc(double A[SIZE_OF_MATRIX][SIZE_OF_MATRIX], double B[SIZE_OF_MATRIX][SIZE_OF_MATRIX], double C[SIZE_OF_MATRIX][SIZE_OF_MATRIX])
{
	int i, j, k, ii, jj, kk;
	static const int size = 25;

	for (ii = 0; ii < SIZE_OF_MATRIX; ii += size)
		for (jj = 0; jj < SIZE_OF_MATRIX; jj += size)
			for (kk = 0; kk < SIZE_OF_MATRIX; kk += size)
				for (i = ii; i < ii + size; i++)
					for (k = kk; k < kk + size; k++)
						for (j = jj; j < jj + size; j++)
							C[i][j] += A[i][k] * B[k][j];
}

void calc(double A[SIZE_OF_MATRIX][SIZE_OF_MATRIX], double B[SIZE_OF_MATRIX][SIZE_OF_MATRIX], double C[SIZE_OF_MATRIX][SIZE_OF_MATRIX])
{
	int i, j, k;

	for (i = 0; i < SIZE_OF_MATRIX; i++)
		for (k = 0; k < SIZE_OF_MATRIX; k++)
			for (j = 0; j < SIZE_OF_MATRIX; j++)
				C[i][j] += A[i][k] * B[k][j];
}

int main(int argc, char **argv)
{
	clock_t start, end;

	static double A[SIZE_OF_MATRIX][SIZE_OF_MATRIX], B[SIZE_OF_MATRIX][SIZE_OF_MATRIX];
	static double X1[SIZE_OF_MATRIX][SIZE_OF_MATRIX], X2[SIZE_OF_MATRIX][SIZE_OF_MATRIX];

	int i, j;
	rep(i, SIZE_OF_MATRIX) rep(j, SIZE_OF_MATRIX)
	{
		A[i][j] = i * SIZE_OF_MATRIX + j;
		B[i][j] = i * SIZE_OF_MATRIX + j;
		X1[i][j] = 0.0;
		X2[i][j] = 0.0;
	}

	// calc
	start = clock();
	calc(A, B, X1);
	end = clock();
	cout << "nomal time = " << (double)(end - start) / CLOCKS_PER_SEC << endl;

	// block calc
	start = clock();
	block_calc(A, B, X2);
	end = clock();
	cout << "block time = " << (double)(end - start) / CLOCKS_PER_SEC << endl;
	
	// block calc 2
	start = clock();
	block_calc2(A, B, X2);
	end = clock();
	cout << "block2 time = " << (double)(end - start) / CLOCKS_PER_SEC << endl;
}
