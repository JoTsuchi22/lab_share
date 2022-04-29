/*****************************************************************************
 * integration_2D.c :
 *****************************************************************************
 * 2次元でメッシュの面積を計算
 * 
 *****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define DIMENSION 2
#define GP_1D 3         // 1方向のガウス点数
#define GP_2D GP_1D * GP_1D
#define ELEMENT_N 100   // 要素数最大値
#define PI 3.141592653589793238462643

static double GP_temp[GP_1D];
static double w_temp[GP_1D];
static double GP[GP_2D * DIMENSION];
static double w[GP_2D * DIMENSION];
double element_GP[GP_1D];

FILE *fp;

void Make_gauss_array();
double det_2D(double a[4]);
double J_2D();
void Integration(int element_num);

int main()
{
    int i, j, k, l;

    Make_gauss_array();
    Integration();

    
}


void Make_gauss_array()
{
    int i, j;

    // ガウス点3x3の場合
    GP_temp[0] =  - sqrt(3.0 / 5.0);
    GP_temp[1] = 0.0;
    GP_temp[2] = sqrt(3.0 / 5.0);
    w_temp[0] = 5.0 / 9.0;
    w_temp[1] = 8.0 / 9.0;
    w_temp[2] = 5.0 / 9.0;

    for (i = 0; i < GP_1D; i++)
    {
        for (j = 0; j < GP_1D; j++)
        {
            GP[(i * GP_1D + j) * DIMENSION + 0] = GP_temp[j];
            GP[(i * GP_1D + j) * DIMENSION + 1] = GP_temp[i];
            w[i * GP_1D + j] = w_temp[i] * w_tmep[j];
        }
    }
}


double J_2D(double temp_delta_el)
{
    int i, j, k;
    double a[2][2];
    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < 2; j++)
        {
            a[i][j] = 0.0;
        }
    }

    for (i = 0; i < 2; i++)
    {
        for (j = 0; j < 2; j++)
        {
            for (k = 0; k < 4; k++)
            {
                a[i][j] += X[k][i] * deriv;
            }
        }
    }
    return det_2D(a); // partial x / partial xi
}

double det_2D(double a[4])
{
    return a[0] * a[3] - a[1] * a[2];
}

double deriv(int axis, double xi, double eta)
{
    double result;
    if (axis == 0) // partial N (xi, eta) / partial xi
    {
        if (num == 0)
        {
            result = - (1.0 - eta) / 4.0;
        }
        else if (num == 1)
        {
            result = (1.0 - eta) / 4.0;
        }
        else if (num == 2)
        {
            result = (1.0 + eta) / 4.0;
        }
        else if (num == 3)
        {
            result = - (1.0 + eta) / 4.0;
        }
    }
    else if (axis == 1) // partial N (xi, eta) / partial eta
    {
        if (num == 0)
        {
            result = - (1.0 - xi) / 4.0;
        }
        else if (num == 1)
        {
            result = - (1.0 + xi) / 4.0;
        }
        else if (num == 2)
        {
            result = (1.0 + xi) / 4.0;
        }
        else if (num == 3)
        {
            result = (1.0 - xi) / 4.0;
        }
    }

    return result;
}

void Integration(double *coordinete)
{
    int i, j;
    double result = 0.0;

    char filename = "result.txt";
    fp = fopen(filename, "w");
    fprintf(fp, "elementnum\tintegration_result\n");

    for (i = 0; i < Element; i++)
    {
        double temp_result = 0.0;
        for (j = 0; j < GP_2D; j++)
        {
            temp_result += J_2D() * w[j];
        }
        result += temp_result;
        fprintf(fp, "%d\t%.15e", i, result);
    }
    fprintf(fp, "Total_integration_esult\n");
    fprintf(fp, "%.15e", result);
    fclose(fp);
}