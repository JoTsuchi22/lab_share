/*****************************************************************************
 * integration_1D.c :
 *****************************************************************************
 * 1次元のガウス積分で刻み幅と誤差の関係を出力
 * 
 *****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define MAX_GP_1D 3     // 1方向のガウス点数
#define MAX_DELTA_EXP 6 // 分割数の指数の上限 n ^ MAX_DELTA_EXP
#define PI 3.141592653589793238462643

static double GP[MAX_GP_1D];
static double w[MAX_GP_1D];
double element_GP[MAX_GP_1D];

FILE *fp;

void Make_gauss_array(int temp_n_int);
double J_1D(double temp_delta_el);

int main()
{
    int i, j, k, l;  
    double scale = PI / 2.0;    // 定積分範囲 ex) 0 <= x <= 2 -> 2

    for (i = 0; i < MAX_GP_1D; i++)
    {
        double n_int = i + 1;
        Make_gauss_array(n_int);

        char filename[256];
        filename[0] = i + 1 + '0';
        filename[1] = '.';
        filename[2] = 't';
        filename[3] = 'x';
        filename[4] = 't';
        filename[5] = '\0';
        fp = fopen(filename, "w");

        fprintf(fp, "積分点:    %d\n\n", i + 1);
        fprintf(fp, "分割数  Error\n");

        for (j = 0; j < MAX_DELTA_EXP; j++)
        {
            double result_int = 0.0;
            int Total_element = (int)pow(10.0, j);
            double delta_el = 1.0 / pow(10.0, j) * scale;   // 刻み幅，1要素の大きさ
            double J = J_1D(delta_el);

            for (k = 0; k < Total_element; k++)
            {
                double element_center = (k + 0.5) * delta_el;
                for (l = 0; l < n_int; l++)
                {
                    element_GP[l] = element_center + (GP[l] * delta_el / 2.0);
                    result_int += sin(element_GP[l]) * J * w[l];
                }
            }

            double Error = fabs(result_int - 1.0); // 解 = 1.0
            fprintf(fp, "%*d    %.20e\n", - (MAX_DELTA_EXP + 2), (int)pow(10.0, j), Error);
        }
        fclose(fp);
    }
}


void Make_gauss_array(int temp_n_int)
{
    if (temp_n_int == 1)
    {
        GP[0] = 0.0;
        w[0] = 2.0;
    }
    else if (temp_n_int == 2)
    {
        GP[0] = - sqrt(1.0 / 3.0);
        GP[1] = sqrt(1.0 / 3.0);
        w[0] = 1.0;
        w[1] = 1.0;
    }
    else if (temp_n_int == 3)
    {
        GP[0] =  - sqrt(3.0 / 5.0);
        GP[1] = 0.0;
        GP[2] = sqrt(3.0 / 5.0);
        w[0] = 5.0 / 9.0;
        w[1] = 8.0 / 9.0;
        w[2] = 5.0 / 9.0;
    }
}


double J_1D(double temp_delta_el)
{
    return temp_delta_el / 2.0; // partial x / partial xi
}