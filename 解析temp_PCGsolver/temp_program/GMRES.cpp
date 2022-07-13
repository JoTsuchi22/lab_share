#include <iostream>
#include <stdlib.h>

using namespace std;

double dot(int length, double *a, double *b)
{
    int i;
    double result = 0.0;
    
    for (i = 0; i < length; i++)
        result += a[i] * b[i];
    
    return result;
}

void GMRES(int length, double *result_vec, double *right_vec, double *matrix)
{
    int i, j, k, l;
    const int max_length = 100;
    const int max_itr = 100;
    const double tol = 1.0e-14;
    static double H[max_itr + 1][max_itr + 1];
    static double V[max_itr + 1][length];

    double norm;

    double *Ax = (double *)calloc(length, sizeof(double));
    double *r = (double *)malloc(sizeof(double) * length);
    double *v = (double *)malloc(sizeof(double) * length);
    double *v_var = (double *)malloc(sizeof(double) * length);
    double *Av = (double *)malloc(sizeof(double) * length);
    double *h = (double *)malloc(sizeof(double) * length);
    
    // x0
    for (i = 0; i < length; i++)
        x[i] = 0.0;

    // Ax
    for (i = 0; i < length; i++)
        for (j = 0; j < length; j++)
            Ax[i] += matrix[i * length + j] * x[j];

    // r initialization
    for (i = 0; i < length; i++)
        r[i] = right_vec[i] - Ax[i];
    
    // l2 norm
    norm = 0;
    for (i = 0; i < length; i++)
        norm += pow(r[i], 2);
    norm = sqrt(norm);

    // v
    for (i = 0; i < length; i++)
        v[i] = r[i] / norm;
    
    // loop
    for (i = 0; i < max_itr; i++)
    {
        // Av
        for (j = 0; j < length; j++)
            for (k = 0; k < length; k++)
                Av[j] += matrix[j * length + k] * v[k];
        
        // V
        for (j = 0; j < length; j++)
            V[i][j] = v[j];

        // h
        for (j = 0; j < i; j++)
            h[j] = dot(length, Av, V[j]);

        // v_var
        for (j = 0; j < length; j++)
        {
            v_var[j] = Av[j];
            for (k = 0; k < i; k++)
                for (l = 0; l < length; l++)
                    v_var[j] -= h[k] * V[k][l];
        }

        // v_var l2 norm
        double v_var_norm = 0.0;
        for (j = 0; j < length; j++)
            v_var_norm += pow(v_var[j], 2);
        v_var_norm = sqrt(v_var_norm);

        // h
        h[i + 1] = v_var_norm;

        // H
        for (j = 0; j < i + 1)
            H[j][i] = h[j];

        // 
    }
}


int main(int argc, char **argv) {

    int i, j;

    // input
    const int size_matrix = 4 * 4, size_vector = 4;
    double *A = (double *)malloc(sizeof(double) * size_matrix);
    double *x = (double *)calloc(size_vector, sizeof(double));
    double *b = (double *)malloc(sizeof(double) * size_vector);

    A[0] = 3; A[1] = 1; A[2] = 1; A[3] = 2;
    A[4] = 5; A[5] = 1; A[6] = 3; A[7] = 4;
    A[8] = 2; A[9] = 0; A[10] = 1; A[11] = 0;
    A[12] = 1; A[13] = 3; A[14] = 2; A[15] = 1;

    double *A_inv = (double *)malloc(sizeof(double) * size_matrix);
    A_inv[0] = 0.5; A_inv[1] = -0.22727272727273; A_inv[2] = 0.36363636363636; A_inv[3] = -0.090909090909091;
    A_inv[4] = 0.5; A_inv[5] = -0.31818181818182; A_inv[6] = -0.090909090909091; A_inv[7] = 0.27272727272727;
    A_inv[8] = -1; A_inv[9] = 0.45454545454545; A_inv[10] = 0.27272727272727; A_inv[11] = 0.18181818181818;
    A_inv[12] = 0; A_inv[13] = 0.27272727272727; A_inv[14] = -0.63636363636364; A_inv[15] = -0.090909090909091;

    double *x_exact = (double *)calloc(size_vector, sizeof(double));

    b[0] = 1, b[1] = 2, b[2] = 3, b[3] = 4;

    // exact
    for (i = 0; i < size_vector; i++)
        for (j = 0; j < size_vector; j++)
            x_exact[i] += A_inv[i * size_vector + j] * b[j];
    
    for (i = 0; i < size_vector; i++)
        cout << x_exact[i] << endl;
    cout << endl;

    // calc
    GMRES(size_matrix, x, b, A);

    // output
    for (i = 0; i < size_vector; i++)
        cout << x[i] << endl;
}