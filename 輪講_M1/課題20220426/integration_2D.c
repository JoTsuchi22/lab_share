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
void Integration(int element_num);
double J_2D();

int main()
{
    int i, j, k, l;  
    double scale = PI / 2.0;    // 定積分範囲 ex) 0 <= x <= 2 -> 2

    Make_gauss_array();

    for (i = 0; i < GP_1D; i++)
    {
        double n_int = i + 1;

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
        fclose(fp);
    }
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


double J_1D(double temp_delta_el)
{
    return temp_delta_el / 2.0; // partial x / partial xi
}

void Integration(double *coordinete)
{

}