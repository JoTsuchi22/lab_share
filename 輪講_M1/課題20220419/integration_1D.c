/*****************************************************************************
 * integration_1D.c :
 *****************************************************************************
 * 1次元のガウス積分で刻み幅と誤差の関係を出力
 * 
 *****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define MAX_GP_1D 3
#define PI 3.141592653589793238462643

static double GP[MAX_GP_1D];
static double w[MAX_GP_1D];

FILE *fp;


int main()
{
    int i, j, k, l;  
    double scale = PI / 2.0;    // 定積分範囲 ex) 0 <= x <= 2 -> 2

    for (i = 0; i < MAX_GP_1D; i++)
    {
        double n_int = i + 1;
        double element_GP[n_int];
        Make_gauss_array(n_int);

        for (j = 0; j < 5; j++)
        {
            double result_int = 0.0;
            int Total_element = (int)(pow(10.0, j));
            double delta_el = 1.0 / pow(10.0, j) * scale;   // 刻み幅，1要素の大きさ
            printf("delta : %le\n", delta_el);

            // double *te = (double *)malloc(sizeof(double) * Total_element);
            // if (te == NULL)
            // {
            //     printf("Memory allocation error\n");
            //     exit(1);
            // }

            double J = J_1D(delta_el);

            for (k = 0; k < Total_element; k++)
            {
                double element_center = ;
                for (l = 0; l < n_int; l++)
                {
                    element_GP[l] = element_center + (GP[l] * delta_el)
                    result_int += * J * w[l];
                }
            }

            // free(te);
        }
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