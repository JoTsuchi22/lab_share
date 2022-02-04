#include <stdio.h>
#include <stdlib.h>

void func(double *B, double **D) // 引数注意
{
    int n = 5, m = 5;
    // A, B 関数より上でmallocした場合は引数にポインタの先頭アドレス入れるだけでよい
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            B[i * m + j] = i * m + j;
        }
    }

    // C, D 関数内でmallocする場合は以下のように書いて先頭のポインタアドレスを渡す
    // *D = (double *)malloc(sizeof(double) * n * m);
    *D = malloc(sizeof(double) * n * m);

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            D[0][i * m + j] = (i * m + j) * 10.0; // この関数内で代入する場合はD[0][??]となり多次元配列となる(?)
        }
    }
}

int main()
{
    int n = 5, m = 5;
    double *A, *C; //Cは普通に宣言
    A = (double *)malloc(sizeof(double) * n * m);
    // C = (double *)malloc(sizeof(double) * 1);

    func(A, &C); //&に注意
    
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            // C[i * m + j] = (i * m + j) * 10.0;
            printf("%le\n", A[i * m + j]);
            printf("%le\n", C[i * m + j]);
        }
    }

    free(A), free(C);
    return 0;
}