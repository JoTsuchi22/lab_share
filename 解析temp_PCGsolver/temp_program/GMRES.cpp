#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <vector>

using namespace std;

template<class T1, class T2>
double dot(int length, T1 &a, T2 &b)
{
    int i;
    double result = 0.0;

    for (i = 0; i < length; i++)
        result += a[i] * b[i];

    return result;
}

void GMRES(int length, double *result_vec, double *right_vec, double *matrix)
{
    int i, j, k, l, m;
    static const int max_itr = 60;
    static const double tol = 1.0e-14;

    vector<vector<double>> H(max_itr + 1, vector<double>(max_itr + 1));
    vector<vector<double>> V(max_itr + 1, vector<double>(max_itr + 1));
    vector<vector<double>> Omega(max_itr + 1, vector<double>(max_itr + 1));
    vector<vector<double>> R(max_itr + 1, vector<double>(max_itr + 1));
    vector<vector<double>> R_temp(max_itr + 1, vector<double>(max_itr + 1));
    vector<vector<double>> R_H(max_itr + 1, vector<double>(max_itr + 1));

    vector<double> e1(max_itr + 1);
    vector<double> g(max_itr + 1);
    vector<double> g_temp(max_itr + 1);
    vector<double> y(max_itr + 1);

    vector<double> Ax(length);
    vector<double> r(length);
    vector<double> v(max_itr + 1);
    vector<double> v_hat(max_itr + 1);
    vector<double> Av(length);
    vector<double> h(max_itr + 1);

    double *Ax = (double *)calloc(length, sizeof(double));
    double *r = (double *)malloc(sizeof(double) * length);
    double *v = (double *)malloc(sizeof(double) * (max_itr + 1));
    double *v_hat = (double *)malloc(sizeof(double) * (max_itr + 1));
    double *Av = (double *)malloc(sizeof(double) * length);
    double *h = (double *)malloc(sizeof(double) * (max_itr + 1));

    int itr;
    double norm;


    for (i = 0; i < max_itr + 1; i++)
        for (j = 0; j < max_itr + 1; j++)
        {
            H[i][j] = 0.0;
            V[i][j] = 0.0;
        }

    // x0
    for (i = 0; i < length; i++)
        result_vec[i] = 0.0;

    // Ax0
    for (i = 0; i < length; i++)
        Ax[i] = 0.0;
    for (i = 0; i < length; i++)
        for (j = 0; j < length; j++)
            Ax[i] += matrix[i * length + j] * result_vec[j];

    // r initialization
    for (i = 0; i < length; i++)
        r[i] = right_vec[i] - Ax[i];

    // l2 norm
    norm = 0.0;
    for (i = 0; i < length; i++)
        norm += pow(r[i], 2);
    norm = sqrt(norm);

    // e1
    e1[0] = norm;
    for (i = 1; i < max_itr + 1; i++)
        e1[i] = 0.0;

    // v
    for (i = 0; i < length; i++)
        v[i] = r[i] / norm;

    // loop
    for (i = 0; i < max_itr; i++)
    {
        // Av
        for (j = 0; j < length; j++)
            Av[j] = 0.0;
        for (j = 0; j < length; j++)
            for (k = 0; k < length; k++)
                Av[j] += matrix[j * length + k] * v[k];

        // V
        for (j = 0; j < length; j++)
            V[i][j] = v[j];

        // h
        for (j = 0; j <= i; j++)
            h[j] = dot<double *, vector<double>>(length, Av, V[j]);

        // v_hat
        for (j = 0; j < length; j++)
        {
            v_hat[j] = Av[j];
            for (k = 0; k <= i; k++)
                v_hat[j] -= h[k] * V[k][j];
        }

        // v_hat l2 norm
        norm = 0.0;
        for (j = 0; j < length; j++)
            norm += pow(v_hat[j], 2);
        norm = sqrt(norm);

        // h
        h[i + 1] = norm;

        // H
        for (j = 0; j <= i + 1; j++)
            H[j][i] = h[j];
        
        // v
        for (j = 0; j < length; j++)
            v[j] = v_hat[j] / h[i + 1];

        // givens rotation
        for (j = 0; j <= i; j++)
        {
            double nu = sqrt(pow(H[j][j], 2) + pow(H[j + 1][j], 2));
            double c_i = H[j][j] / nu;
            double s_i = H[j + 1][j] / nu;

            // Omega
            for (k = 0; k <= i + 1; k++)
                for (l = 0; l <= i + 1; l++)
                {
                    if (k == l)
                        Omega[k][l] = 1.0;
                    else
                        Omega[k][l] = 0.0;
                }
            Omega[j][j] = c_i;
            Omega[j][j + 1] = s_i;
            Omega[j + 1][j] = - s_i;
            Omega[j + 1][j + 1] = c_i;

            // R
            if (j == 0)
                for (k = 0; k <= i + 1; k++)
                    for (l = 0; l <= i + 1; l++)
                        R[k][l] = Omega[k][l];
            else
            {
                for (k = 0; k <= i + 1; k++)
                    for (l = 0; l <= i + 1; l++)
                        R_temp[k][l] = 0.0;
                for (k = 0; k <= i + 1; k++)
                    for (m = 0; m <= i + 1; m++)
                        for (l = 0; l <= i + 1; l++)
                            R_temp[k][l] += Omega[k][m] * R[m][l];
                            // R_temp[k][l] += R[k][m] * Omega[m][l];
                for (k = 0; k <= i + 1; k++)
                    for (l = 0; l <= i + 1; l++)
                        R[k][l] = R_temp[k][l];
            }
        }

        // R_H
        for (j = 0; j <= i + 1; j++)
            for (k = 0; k <= i + 1; k++)
                R_H[j][k] = 0.0;
        for (j = 0; j <= i + 1; j++)
            for (l = 0; l <= i + 1; l++)
                for (k = 0; k <= i; k++)
                    R_H[j][k] += R[j][l] * H[l][k];

        // g
        for (j = 0; j <= i + 1; j++)
        {
            g[j] = 0.0;
            g_temp[j] = 0.0;
            for (k = 0; k <= i + 1; k++)
                g[j] += R[j][k] * e1[k];
            for (k = 0; k <= i + 1; k++)
                g_temp[j] = g[j];
        }

        // y
        for (j = i; j >= 0; j--)
        {
            y[j] = g[j] / R_H[j][j];
            for (k = 0; k <= j; k++)
                g[k] -= R_H[k][j] * y[j];
        }

        // r_i norm
        norm = 0.0;
        for (j = 0; j <= i; j++)
        {
            double temp = g_temp[j];
            for (k = 0; k <= i; k++)
                temp -= R_H[j][k] * y[k];
            norm += pow(temp, 2);
        }
        norm += pow(g_temp[i + 1], 2);
        norm = sqrt(norm);

        cout << "itr = " << i <<  "  ";
        cout << "tol = " << norm << endl;
        itr = i;
        if (norm <= tol)
            break;
    }

    // result
    for (i = 0; i < length; i++)
        result_vec[i] = 0.0;
    for (i = 0; i < length; i++)
        for (j = 0; j <= itr; j++)
            result_vec[i] += V[j][i] * y[j];
}

int main()
{
    int i, j;

    // input
    const int size_matrix = 4 * 4, size_vector = 4;
    double *A = (double *)malloc(sizeof(double) * size_matrix);
    double *x = (double *)calloc(size_vector, sizeof(double));
    double *b = (double *)malloc(sizeof(double) * size_vector);

    A[0] = 3;
    A[1] = 1;
    A[2] = 1;
    A[3] = 2;
    A[4] = 5;
    A[5] = 1;
    A[6] = 3;
    A[7] = 4;
    A[8] = 2;
    A[9] = 0;
    A[10] = 1;
    A[11] = 0;
    A[12] = 1;
    A[13] = 3;
    A[14] = 2;
    A[15] = 1;

    double *A_inv = (double *)malloc(sizeof(double) * size_matrix);
    A_inv[0] = 0.5;
    A_inv[1] = -0.22727272727273;
    A_inv[2] = 0.36363636363636;
    A_inv[3] = -0.090909090909091;
    A_inv[4] = 0.5;
    A_inv[5] = -0.31818181818182;
    A_inv[6] = -0.090909090909091;
    A_inv[7] = 0.27272727272727;
    A_inv[8] = -1;
    A_inv[9] = 0.45454545454545;
    A_inv[10] = 0.27272727272727;
    A_inv[11] = 0.18181818181818;
    A_inv[12] = 0;
    A_inv[13] = 0.27272727272727;
    A_inv[14] = -0.63636363636364;
    A_inv[15] = -0.090909090909091;

    double *x_exact = (double *)calloc(size_vector, sizeof(double));

    b[0] = 1, b[1] = 2, b[2] = 3, b[3] = 4;

    // calc
    GMRES(size_vector, x, b, A);

    // output
    cout << "GMRES" << endl;
    for (i = 0; i < size_vector; i++)
        cout << x[i] << endl;

    // exact
    for (i = 0; i < size_vector; i++)
        for (j = 0; j < size_vector; j++)
            x_exact[i] += A_inv[i * size_vector + j] * b[j];

    cout << "exact" << endl;
    for (i = 0; i < size_vector; i++)
        cout << x_exact[i] << endl;
    cout << endl;
}