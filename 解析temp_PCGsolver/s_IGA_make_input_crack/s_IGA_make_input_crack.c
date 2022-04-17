/*****************************************************************************
 * s_IGA_make_input_crack.c :
 *****************************************************************************
 * s_IGAにおける，き裂の特異パッチのインプットデータを作成するためのプログラムです
 * このプログラムのインプットデータについてはREADMEを読んでください :)
 * 
 * 2022.4.14
 *****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define DIMENSION 2                 // 2次元
#define MERGE_DISTANCE 1.0e-13      // コントロールポイントが同じ点と判定する距離
// #define MAX_DISP_CONSTRAINT 10      // 変位指定する変位量の最大個数
// #define MAX_DISP_CONSTRAINT_EDGE 10 // 変位指定する辺の最大個数
// #define MAX_DISTRIBUTED_LOAD 5      // 分布荷重の最大個数

FILE *fp;

static double E_and_nu[2];
static int mode[2];                 // 0: analysis model mode, 1: making mode, 2: sigular patch mode
static double crack_tip[DIMENSION];
static double singular_width;
static double length_each_side[3];  // (outer side, inner side, upper side)
static int cp_after_each_side[3];   // (outer side, inner side, upper side)


void Get_input_a(char *filename);
void Get_input_b(char *filename, int num, double *temp_CP, double *temp_KV, int *temp_CP_info, int *temp_KV_info);
void Make_patch_file(int num);
void Make_model(char **temp_argv);
void Make_CP_KV_quarter_model(double *temp_CP, double *temp_KV, int *temp_CP_info, int *temp_KV_info, char **temp_argv);
void Make_CP_KV_full_model(double *temp_CP, double *temp_KV, int *temp_CP_info, int *temp_KV_info);


int main(int argc, char **argv)
{
    int i;

    // ファイル読み込み
    Get_input_a(argv[1]);

    if (argc == 2)
    {
        // エラー出力
        if (mode[0] != 0 && mode[0] != 1)
        {
            printf("Error analysis model mode must be 0 or 1\n");
            exit(1);
        }
        if (mode[0] == 1 && fabs(crack_tip[0] - length_each_side[2]) >= MERGE_DISTANCE)
        {
            printf("Error\nif auto making mode == 1, upper side must match tip crack coordinate x\n");
            exit(1);
        }

        // 各パッチ作成
        int singular_patch_n, Total_patch_n;
        if (mode[0] == 0)
        {
            singular_patch_n = 4;
            Total_patch_n = 10;
        }
        else if (mode[0] == 1)
        {
            singular_patch_n = 16;
            Total_patch_n = 40;
        }

        for (i = singular_patch_n; i < Total_patch_n; i++)
        {
            Make_patch_file(i);
        }

        return 0;
    }

    // 本解析プログラムのインプットデータ作成
    Make_model(argv);

    return 0;
}


void Get_input_a(char *filename)
{
    int i;
    char s[256];

    int temp_i;
    double temp_d;

    fp = fopen(filename, "r");

    // DIMENSION
    fscanf(fp, "%d", &temp_i);
    if (temp_i != DIMENSION)
    {
        printf("Error, at input file\nDIMENSION must be 2 in this program\n");
        exit(1);
    }

    fgets(s, 256, fp);

    // ヤング率, ポアソン比
    for (i = 0; i < 2; i++)
    {
        fscanf(fp, "%lf", &temp_d);
        E_and_nu[i] = temp_d;
        printf("E_and_nu[%d] = %le\n", i, E_and_nu[i]);
    }

    fgets(s, 256, fp);

    // mode[0]
    fscanf(fp, "%d", &temp_i);
    mode[0] = temp_i;
    printf("mode[0] = %d\n", mode[0]);

    fgets(s, 256, fp);

    // mode[1]
    fscanf(fp, "%d", &temp_i);
    mode[1] = temp_i;
    printf("mode[1] = %d\n", mode[1]);

    fgets(s, 256, fp);

    // crack_tip
    for (i = 0; i < DIMENSION; i++)
    {
        fscanf(fp, "%lf", &temp_d);
        crack_tip[i] = temp_d;
        printf("crack_tip[%d] = %le\n", i, crack_tip[i]);

        if (i == 1 && fabs(crack_tip[1]) >= MERGE_DISTANCE)
        {
            printf("tip crack coordinate y must be '0.0'.\n");
            exit(1);
        }
    }

    fgets(s, 256, fp);

    // singular_width
    fscanf(fp, "%lf", &temp_d);
    singular_width = temp_d;
    printf("singular_width = %le\n", singular_width);

    fgets(s, 256, fp);

    // length_each_side
    for (i = 0; i < 3; i++)
    {
        fscanf(fp, "%lf", &temp_d);
        length_each_side[i] = temp_d;
        printf("length_each_side[%d] = %le\n", i, length_each_side[i]);
    }

    fgets(s, 256, fp);

    // cp_after_each_side
    for (i = 0; i < 3; i++)
    {
        fscanf(fp, "%d", &temp_i);
        cp_after_each_side[i] = temp_i;
        printf("cp_after_each_side[%d] = %d\n", i, cp_after_each_side[i]);
    }

    fclose(fp);
}


void Get_input_b(char *filename, int num, double *temp_CP, double *temp_KV, int *temp_CP_info, int *temp_KV_info)
{
    int i, j;
    char s[256];

    int temp_i;
    double temp_d;

    int singular_patch_n;
    if (mode[0] == 0)
    {
        singular_patch_n = 4;
    }
    else if (mode[0] == 1)
    {
        singular_patch_n = 16;
    }

    int patch_num = num + singular_patch_n;
    int CP_to_here = 0;
    int KV_to_here = 0;
    for (i = 0; i < patch_num; i++)
    {
        CP_to_here += temp_CP_info[i * DIMENSION] * temp_CP_info[i * DIMENSION + 1] * (DIMENSION + 1);
        for (j = 0; j < DIMENSION; j++)
        {
            KV_to_here += temp_KV_info[i * DIMENSION + j];
        }
    }

    fp = fopen(filename, "r");

    // Total CP
    fscanf(fp, "%d", &temp_i);

    fgets(s, 256, fp);

    // Order
    for (i = 0; i < DIMENSION; i++)
    {
        fscanf(fp, "%d", &temp_i);
    }

    fgets(s, 256, fp);

    // Knot Vector length
    for (i = 0; i < DIMENSION; i++)
    {
        fscanf(fp, "%d", &temp_i);
    }

    fgets(s, 256, fp);

    // CP xi eta
    for (i = 0; i < DIMENSION; i++)
    {
        fscanf(fp, "%d", &temp_i);
    }

    fgets(s, 256, fp);

    // Knot Vector
    for (i = 0; i < DIMENSION; i++)
    {
        for (j = 0; j < temp_KV_info[patch_num * DIMENSION + i]; j++)
        {
            fscanf(fp, "%lf", &temp_d);
            temp_KV[KV_to_here + j] = temp_d;
        }
        KV_to_here += temp_KV_info[patch_num * DIMENSION + i];
    }

    fgets(s, 256, fp);

    // CP
    for (i = 0; i < temp_CP_info[patch_num * DIMENSION] * temp_CP_info[patch_num * DIMENSION + 1]; i++)
    {
        fscanf(fp, "%d", &temp_i);
        for (j = 0; j < DIMENSION + 1; j++)
        {
            fscanf(fp, "%lf", &temp_d);
            temp_CP[CP_to_here + j] = temp_d;
        }
        CP_to_here += DIMENSION + 1;
    }

    fclose(fp);
}


void Make_patch_file(int num)
{
    int counter = 0;
    int array[2];
    char str[256];

    if (num < 10)
    {
        str[0] = num + '0';
        counter++;
    }
    else if (num >= 10)
    {
        array[0] = num / 10;
        array[1] = num % 10;
        str[0] = array[0] + '0';
        str[1] = array[1] + '0';
        counter += 2;
    }

    str[counter] = '.';
    str[counter + 1] = 't';
    str[counter + 2] = 'x';
    str[counter + 3] = 't';
    str[counter + 4] = '\0';
    
    fp = fopen(str, "w");

    // DIMENSION
    fprintf(fp, "%d\n\n", DIMENSION);

    // Total CP
    fprintf(fp, "%d\n\n", 4);

    // Order
    fprintf(fp, "%d   %d\n\n", 1, 1);

    // Knot Vector length
    fprintf(fp, "%d   %d\n\n", 4, 4);

    // CP xi eta
    fprintf(fp, "%d   %d\n\n", 2, 2);

    // Knot Vector
    fprintf(fp, "%.17e   %.17e   %.17e   %.17e\n", 0.0, 0.0, 1.0, 1.0);
    fprintf(fp, "%.17e   %.17e   %.17e   %.17e\n\n", 0.0, 0.0, 1.0, 1.0);

    // CP coordinate
    if (mode[0] == 0)
    {
        if (num == 4)
        {
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] + singular_width, crack_tip[1], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] + length_each_side[0], crack_tip[1], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] + singular_width, crack_tip[1] + singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] + length_each_side[0], crack_tip[1] + singular_width, 1.0);
        }
        else if (num == 5)
        {
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] + singular_width, crack_tip[1] + singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] + length_each_side[0], crack_tip[1] + singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] + singular_width, crack_tip[1] + length_each_side[2], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] + length_each_side[0], crack_tip[1] + length_each_side[2], 1.0);
        }
        else if (num == 6)
        {
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0], crack_tip[1] + singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] + singular_width, crack_tip[1] + singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0], crack_tip[1] + length_each_side[2], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] + singular_width, crack_tip[1] + length_each_side[2], 1.0);
        }
        else if (num == 7)
        {
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] - singular_width, crack_tip[1] + singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0], crack_tip[1] + singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] - singular_width, crack_tip[1] + length_each_side[2], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0], crack_tip[1] + length_each_side[2], 1.0);
        }
        else if (num == 8)
        {
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] - length_each_side[1], crack_tip[1] + singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] - singular_width, crack_tip[1] + singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] - length_each_side[1], crack_tip[1] + length_each_side[2], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] - singular_width, crack_tip[1] + length_each_side[2], 1.0);
        }
        else if (num == 9)
        {
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] - length_each_side[1], crack_tip[1], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] - singular_width, crack_tip[1], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] - length_each_side[1], crack_tip[1] + singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] - singular_width, crack_tip[1] + singular_width, 1.0);
        }
    }
    else if (mode[1] == 1)
    {
        if (num == 16)
        {
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] + singular_width, crack_tip[1], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] + length_each_side[0], crack_tip[1], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] + singular_width, crack_tip[1] + singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] + length_each_side[0], crack_tip[1] + singular_width, 1.0);
        }
        else if (num == 17)
        {
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] + singular_width, crack_tip[1] + singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] + length_each_side[0], crack_tip[1] + singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] + singular_width, crack_tip[1] + length_each_side[2], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] + length_each_side[0], crack_tip[1] + length_each_side[2], 1.0);
        }
        else if (num == 18)
        {
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0], crack_tip[1] + singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] + singular_width, crack_tip[1] + singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0], crack_tip[1] + length_each_side[2], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] + singular_width, crack_tip[1] + length_each_side[2], 1.0);
        }
        else if (num == 19)
        {
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] - singular_width, crack_tip[1] + singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0], crack_tip[1] + singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] - singular_width, crack_tip[1] + length_each_side[2], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0], crack_tip[1] + length_each_side[2], 1.0);
        }
        else if (num == 20)
        {
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] - length_each_side[1], crack_tip[1] + singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] - singular_width, crack_tip[1] + singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] - length_each_side[1], crack_tip[1] + length_each_side[2], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] - singular_width, crack_tip[1] + length_each_side[2], 1.0);
        }
        else if (num == 21)
        {
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] - length_each_side[1], crack_tip[1], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] - singular_width, crack_tip[1], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] - length_each_side[1], crack_tip[1] + singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] - singular_width, crack_tip[1] + singular_width, 1.0);
        }
        else if (num == 22)
        {
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] - length_each_side[1], crack_tip[1] - singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] - singular_width, crack_tip[1] - singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] - length_each_side[1], crack_tip[1], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] - singular_width, crack_tip[1], 1.0);
        }
        else if (num == 23)
        {
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] - length_each_side[1], crack_tip[1] - length_each_side[2], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] - singular_width, crack_tip[1] - length_each_side[2], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] - length_each_side[1], crack_tip[1] - singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] - singular_width, crack_tip[1] - singular_width, 1.0);
        }
        else if (num == 24)
        {
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] - singular_width, crack_tip[1] - length_each_side[2], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0], crack_tip[1] - length_each_side[2], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] - singular_width, crack_tip[1] - singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0], crack_tip[1] - singular_width, 1.0);
        }
        else if (num == 25)
        {
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0], crack_tip[1] - length_each_side[2], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] + singular_width, crack_tip[1] - length_each_side[2], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0], crack_tip[1] - singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] + singular_width, crack_tip[1] - singular_width, 1.0);
        }
        else if (num == 26)
        {
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] + singular_width, crack_tip[1] - length_each_side[2], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] + length_each_side[0], crack_tip[1] - length_each_side[2], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] + singular_width, crack_tip[1] - singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] + length_each_side[0], crack_tip[1] - singular_width, 1.0);
        }
        else if (num == 27)
        {
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] + singular_width, crack_tip[1] - singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] + length_each_side[0], crack_tip[1] - singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] + singular_width, crack_tip[1], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] + length_each_side[0], crack_tip[1], 1.0);
        }
        else if (num == 28)
        {
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] + singular_width - (2.0 * singular_width), crack_tip[1], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] + length_each_side[0] - (2.0 * singular_width), crack_tip[1], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] + singular_width - (2.0 * singular_width), crack_tip[1] + singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] + length_each_side[0] - (2.0 * singular_width), crack_tip[1] + singular_width, 1.0);
        }
        else if (num == 29)
        {
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] + singular_width - (2.0 * singular_width), crack_tip[1] + singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] + length_each_side[0] - (2.0 * singular_width), crack_tip[1] + singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] + singular_width - (2.0 * singular_width), crack_tip[1] + length_each_side[2], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] + length_each_side[0] - (2.0 * singular_width), crack_tip[1] + length_each_side[2], 1.0);
        }
        else if (num == 30)
        {
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] - (2.0 * singular_width), crack_tip[1] + singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] + singular_width - (2.0 * singular_width), crack_tip[1] + singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] - (2.0 * singular_width), crack_tip[1] + length_each_side[2], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] + singular_width - (2.0 * singular_width), crack_tip[1] + length_each_side[2], 1.0);
        }
        else if (num == 31)
        {
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] - singular_width - (2.0 * singular_width), crack_tip[1] + singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] - (2.0 * singular_width), crack_tip[1] + singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] - singular_width - (2.0 * singular_width), crack_tip[1] + length_each_side[2], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] - (2.0 * singular_width), crack_tip[1] + length_each_side[2], 1.0);
        }
        else if (num == 32)
        {
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] - length_each_side[1] - (2.0 * singular_width), crack_tip[1] + singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] - singular_width - (2.0 * singular_width), crack_tip[1] + singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] - length_each_side[1] - (2.0 * singular_width), crack_tip[1] + length_each_side[2], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] - singular_width - (2.0 * singular_width), crack_tip[1] + length_each_side[2], 1.0);
        }
        else if (num == 33)
        {
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] - length_each_side[1] - (2.0 * singular_width), crack_tip[1], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] - singular_width - (2.0 * singular_width), crack_tip[1], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] - length_each_side[1] - (2.0 * singular_width), crack_tip[1] + singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] - singular_width - (2.0 * singular_width), crack_tip[1] + singular_width, 1.0);
        }
        else if (num == 34)
        {
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] - length_each_side[1] - (2.0 * singular_width), crack_tip[1] - singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] - singular_width - (2.0 * singular_width), crack_tip[1] - singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] - length_each_side[1] - (2.0 * singular_width), crack_tip[1], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] - singular_width - (2.0 * singular_width), crack_tip[1], 1.0);
        }
        else if (num == 35)
        {
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] - length_each_side[1] - (2.0 * singular_width), crack_tip[1] - length_each_side[2], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] - singular_width - (2.0 * singular_width), crack_tip[1] - length_each_side[2], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] - length_each_side[1] - (2.0 * singular_width), crack_tip[1] - singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] - singular_width - (2.0 * singular_width), crack_tip[1] - singular_width, 1.0);
        }
        else if (num == 36)
        {
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] - singular_width - (2.0 * singular_width), crack_tip[1] - length_each_side[2], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] - (2.0 * singular_width), crack_tip[1] - length_each_side[2], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] - singular_width - (2.0 * singular_width), crack_tip[1] - singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] - (2.0 * singular_width), crack_tip[1] - singular_width, 1.0);
        }
        else if (num == 37)
        {
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] - (2.0 * singular_width), crack_tip[1] - length_each_side[2], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] + singular_width - (2.0 * singular_width), crack_tip[1] - length_each_side[2], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] - (2.0 * singular_width), crack_tip[1] - singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] + singular_width - (2.0 * singular_width), crack_tip[1] - singular_width, 1.0);
        }
        else if (num == 38)
        {
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] + singular_width - (2.0 * singular_width), crack_tip[1] - length_each_side[2], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] + length_each_side[0] - (2.0 * singular_width), crack_tip[1] - length_each_side[2], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] + singular_width - (2.0 * singular_width), crack_tip[1] - singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] + length_each_side[0] - (2.0 * singular_width), crack_tip[1] - singular_width, 1.0);
        }
        else if (num == 39)
        {
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] + singular_width - (2.0 * singular_width), crack_tip[1] - singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] + length_each_side[0] - (2.0 * singular_width), crack_tip[1] - singular_width, 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] + singular_width - (2.0 * singular_width), crack_tip[1], 1.0);
            fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] + length_each_side[0] - (2.0 * singular_width), crack_tip[1], 1.0);
        }
    }

    // o.e.
    fprintf(fp, "%d   %d\n\n", 1, 1);

    // k.i.
    fprintf(fp, "%d   %d\n\n", 0, 0);
    if (mode[0] == 0)
    {
        // パッチごとのxi方向のコントロールポイント数
        if (num == 4 || num == 5)
        {
            fprintf(fp, "%d   ", cp_after_each_side[0]);
        }
        else if (num == 6 || num == 7)
        {
            fprintf(fp, "%d   ", 3);
        }
        else if (num == 8 || num == 9)
        {
            fprintf(fp, "%d   ", cp_after_each_side[1]);
        }

        // パッチごとのeta方向のコントロールポイント数
        if (num == 4 || num == 9)
        {
            fprintf(fp, "%d\n\n", 3);
        }
        else if (num == 5 || num == 6 || num == 7 || num == 8)
        {
            fprintf(fp, "%d\n\n", cp_after_each_side[2]);
        }
    }
    else if (mode[0] == 1)
    {
        // パッチごとのxi方向のコントロールポイント数
        if (num == 16 || num == 17 || num == 26 || num == 27 || num == 32 || num == 33 || num == 34 || num == 35)
        {
            fprintf(fp, "%d   ", cp_after_each_side[0]);
        }
        else if (num == 18 || num == 19 || num == 24 || num == 25 || num == 30 || num == 31 || num == 36 || num == 37)
        {
            fprintf(fp, "%d   ", 3);
        }
        else if (num == 20 || num == 21 || num == 22 || num == 23 || num == 28 || num == 29 || num == 38 || num == 29)
        {
            fprintf(fp, "%d   ", cp_after_each_side[1]);
        }

        // パッチごとのeta方向のコントロールポイント数
        if (num == 16 || num == 21 || num == 22 || num == 27 || num == 28 || num == 33 || num == 34 || num == 39)
        {
            fprintf(fp, "%d\n\n", 3);
        }
        else if (num == 17  || num == 18 || num == 19 || num == 20 || num == 23 || num == 24 || num == 25 || num == 26 || num == 29 || num == 30 || num == 31 || num == 32 || num == 35 || num == 36 || num == 37 || num == 38)
        {
            fprintf(fp, "%d\n\n", cp_after_each_side[2]);
        }
    }
    fprintf(fp, "%d   %d\n\n", 0, 0);

    fclose(fp);
}


void Make_model(char **temp_argv)
{
    int i, j;
    int Total_patch_loc = 0;
    int singular_patch_n = 0;

    // モデル分岐(1/4 or full)
    if (mode[0] == 0)
    {
        Total_patch_loc = 10;
        singular_patch_n = 4;
    }
    else if (mode[0] == 1)
    {
        Total_patch_loc = 40;
        singular_patch_n = 16;
    }

    // 動的メモリ確保
    int *Order = (int *)malloc(sizeof(int) * Total_patch_loc * DIMENSION);      // int Order[パッチ番号][DIMENSION]
    int *KV_info = (int *)malloc(sizeof(int) * Total_patch_loc * DIMENSION);    // int KV_info[パッチ番号][DIMENSION]
    int *CP_info = (int *)malloc(sizeof(int) * Total_patch_loc * DIMENSION);    // int CP_info[パッチ番号][DIMENSION]

    if (Order == NULL || KV_info == NULL || CP_info == NULL)
    {
        printf("Memory cannot be allocated\n");
        exit(1);
    }

    // 各パッチでの情報を取得
    for (i = 0; i < singular_patch_n; i++)
    {
        for (j = 0; j < DIMENSION; j++)
        {
            Order[i * DIMENSION + j] = 2;   // 基底関数は2次で固定
            if (j == 0)         // xi方向
            {
                // 特異パッチの要素の種類
                if (mode[1] == 1)
                {
                    CP_info[i * DIMENSION + j] = 3;
                }
                else if (mode[1] == 2)
                {
                    CP_info[i * DIMENSION + j] = 4;
                }
            }
            else if (j == 1)    // eta方向
            {
                CP_info[i * DIMENSION + j] = 3;
            }
            KV_info[i * DIMENSION + j] = Order[i * DIMENSION + j] + CP_info[i * DIMENSION + j] + 1;
        }
    }
    for (i = singular_patch_n; i < Total_patch_loc; i++)
    {
        // 全パッチで同一
        for (j = 0; j < DIMENSION; j++)
        {
            Order[i * DIMENSION + j] = 2;   // 基底関数は2次で固定
        }

        // モデル分岐(1/4 or full)
        if (mode[0] == 0)
        {
            // パッチごとのxi方向のコントロールポイント数
            if (i == 4 || i == 5)
            {
                CP_info[i * DIMENSION + 0] = cp_after_each_side[0];
            }
            else if (i == 6 || i == 7)
            {
                CP_info[i * DIMENSION + 0] = 3;
            }
            else if (i == 8 || i == 9)
            {
                CP_info[i * DIMENSION + 0] = cp_after_each_side[1];
            }

            // パッチごとのeta方向のコントロールポイント数
            if (i == 4 || i == 9)
            {
                CP_info[i * DIMENSION + 1] = 3;
            }
            else if (i == 5 || i == 6 || i == 7 || i == 8)
            {
                CP_info[i * DIMENSION + 1] = cp_after_each_side[2];
            }
        }
        else if (mode[0] == 1)
        {
            // パッチごとのxi方向のコントロールポイント数
            if (i == 16 || i == 17 || i == 26 || i == 27 || i == 32 || i == 33 || i == 34 || i == 35)
            {
                CP_info[i * DIMENSION + 0] = cp_after_each_side[0];
            }
            else if (i == 18 || i == 19 || i == 24 || i == 25 || i == 30 || i == 31 || i == 36 || i == 37)
            {
                CP_info[i * DIMENSION + 0] = 3;
            }
            else if (i == 20 || i == 21 || i == 22 || i == 23 || i == 28 || i == 29 || i == 38 || i == 29)
            {
                CP_info[i * DIMENSION + 0] = cp_after_each_side[1];
            }

            // パッチごとのeta方向のコントロールポイント数
            if (i == 16 || i == 21 || i == 22 || i == 27 || i == 28 || i == 33 || i == 34 || i == 39)
            {
                CP_info[i * DIMENSION + 1] = 3;
            }
            else if (i == 17  || i == 18 || i == 19 || i == 20 || i == 23 || i == 24 || i == 25 || i == 26 || i == 29 || i == 30 || i == 31 || i == 32 || i == 35 || i == 36 || i == 37 || i == 38)
            {
                CP_info[i * DIMENSION + 1] = cp_after_each_side[2];
            }
        }

        // Order と CP_info から算出
        for (j = 0; j < DIMENSION; j++)
        {
            KV_info[i * DIMENSION + j] = Order[i * DIMENSION + j] + CP_info[i * DIMENSION + j] + 1;
        }
    }

    int temp1 = 0; // temp1 : 全パッチ含めた総コントロールポイント数
    int temp2 = 0; // temp2 : 4辺のコントロールポイントの和を全パッチ分足した値 * 2
    int temp3 = 0; // temp3 : 全パッチ含めた総ノットベクトル数

    for (i = 0; i < Total_patch_loc; i++)
    {
        temp1 += CP_info[i * DIMENSION] * CP_info[i * DIMENSION + 1];
        temp2 += 2 * (CP_info[i * DIMENSION] + CP_info[i * DIMENSION + 1]);
        temp3 += KV_info[i * DIMENSION] + KV_info[i * DIMENSION + 1];
    }
    printf("Total Control Point = %d\n", temp1);

    // 動的メモリ確保
    double *CP = (double *)malloc(sizeof(double) * temp1 * (DIMENSION + 1));            // double CP[パッチ番号][パッチ内CP番号][xyw -> 3]
    double *CP_result = (double *)malloc(sizeof(double) * temp1 * (DIMENSION + 1));     // double CP_result[通しのコントロールポイント番号(連番)][xyw -> 3]
    int *A = (int *)malloc(sizeof(int) * temp2);                                        // int    A[パッチ番号][辺番号(0~3)][辺内のコネクティビティ]
    double *B = (double *)malloc(sizeof(double) * Total_patch_loc * 16 * (DIMENSION + 1));  // double B[パッチ番号][辺番号(0~7)(正負方向)][各辺の端の2頂点][座標xyw -> 3]
    int *Connectivity = (int *)malloc(sizeof(int) * temp1);                             // int    Connectivity[パッチ番号][パッチ内CP番号]
    double *KV = (double *)malloc(sizeof(double) * Total_patch_loc * temp3);                // double KV[パッチ番号][DIMENSION][ノットベクトル番号]

    if (CP == NULL || CP_result == NULL || A == NULL || B == NULL || Connectivity == NULL || KV == NULL)
    {
        printf("Memory cannot be allocated\n");
        exit(1);
    }

    // 特異パッチのコントロールポイントとノットベクトルを作成
    Make_CP_KV_quarter_model(CP, KV, CP_info, KV_info, temp_argv);

    // 動的メモリ確保
    int *Edge_info = (int *)malloc(sizeof(int) * Total_patch_loc * 32);         // int Edge_info[パッチ番号][own 辺番号(正固定0~3)][opp 辺番号(0~7)]
    int *Opponent_patch_num = (int *)malloc(sizeof(int) * Total_patch_loc * 4); // int Opponent_patch_num[パッチ番号][own 辺番号(正固定0~3]

    if (Edge_info == NULL || Opponent_patch_num == NULL)
    {
        printf("Memory cannot be allocated\n");
        exit(1);
    }

    for (i = 0; i < Total_patch_loc * 32; i++)
    {
        Edge_info[i] = 0;
    }

    // メモリ解放
    free(Order), free(KV_info), free(CP_info);
    free(CP), free(CP_result), free(A), free(B), free(Connectivity), free(KV);
    free(Edge_info), free(Opponent_patch_num);
    // free(length_before), free(length_after), free(Boundary), free(Boundary_result);
}


void Make_CP_KV_quarter_model(double *temp_CP, double *temp_KV, int *temp_CP_info, int *temp_KV_info, char **temp_argv)
{
    int i;
    int CP_to_here = 0;
    int KV_to_here = 0;

    // 特異パッチの要素の種類
    if (mode[1] == 1)
    {
        // パッチ0
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] + singular_width; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0] + singular_width; temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0); temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0] + singular_width; temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 8

        temp_KV[KV_to_here] = 0.0; temp_KV[KV_to_here + 1] = 0.0; temp_KV[KV_to_here + 2] = 0.0; temp_KV[KV_to_here + 3] = 1.0; temp_KV[KV_to_here + 4] = 1.0; temp_KV[KV_to_here + 5] = 1.0; KV_to_here += 6;
        temp_KV[KV_to_here] = 0.0; temp_KV[KV_to_here + 1] = 0.0; temp_KV[KV_to_here + 2] = 0.0; temp_KV[KV_to_here + 3] = 1.0; temp_KV[KV_to_here + 4] = 1.0; temp_KV[KV_to_here + 5] = 1.0; KV_to_here += 6;

        // パッチ1
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] + singular_width; temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 2.0); temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 8

        temp_KV[KV_to_here] = 0.0; temp_KV[KV_to_here + 1] = 0.0; temp_KV[KV_to_here + 2] = 0.0; temp_KV[KV_to_here + 3] = 1.0; temp_KV[KV_to_here + 4] = 1.0; temp_KV[KV_to_here + 5] = 1.0; KV_to_here += 6;
        temp_KV[KV_to_here] = 0.0; temp_KV[KV_to_here + 1] = 0.0; temp_KV[KV_to_here + 2] = 0.0; temp_KV[KV_to_here + 3] = 1.0; temp_KV[KV_to_here + 4] = 1.0; temp_KV[KV_to_here + 5] = 1.0; KV_to_here += 6;

        // パッチ2
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 2.0); temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0] - singular_width; temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 8

        temp_KV[KV_to_here] = 0.0; temp_KV[KV_to_here + 1] = 0.0; temp_KV[KV_to_here + 2] = 0.0; temp_KV[KV_to_here + 3] = 1.0; temp_KV[KV_to_here + 4] = 1.0; temp_KV[KV_to_here + 5] = 1.0; KV_to_here += 6;
        temp_KV[KV_to_here] = 0.0; temp_KV[KV_to_here + 1] = 0.0; temp_KV[KV_to_here + 2] = 0.0; temp_KV[KV_to_here + 3] = 1.0; temp_KV[KV_to_here + 4] = 1.0; temp_KV[KV_to_here + 5] = 1.0; KV_to_here += 6;

        // パッチ3
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] - singular_width; temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0] - singular_width; temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0); temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0] - singular_width; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 8

        temp_KV[KV_to_here] = 0.0; temp_KV[KV_to_here + 1] = 0.0; temp_KV[KV_to_here + 2] = 0.0; temp_KV[KV_to_here + 3] = 1.0; temp_KV[KV_to_here + 4] = 1.0; temp_KV[KV_to_here + 5] = 1.0; KV_to_here += 6;
        temp_KV[KV_to_here] = 0.0; temp_KV[KV_to_here + 1] = 0.0; temp_KV[KV_to_here + 2] = 0.0; temp_KV[KV_to_here + 3] = 1.0; temp_KV[KV_to_here + 4] = 1.0; temp_KV[KV_to_here + 5] = 1.0; KV_to_here += 6;
    }
    else if (mode[1] == 2)
    {
        // パッチ0
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 2.0); temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0] + singular_width; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 2.0); temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 4.0); temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0] + singular_width; temp_CP[CP_to_here + 1] = crack_tip[1] + + (singular_width / 2.0); temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 8
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 9
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 2.0); temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0); temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 10
        temp_CP[CP_to_here] = crack_tip[0] + singular_width; temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 11

        temp_KV[KV_to_here] = 0.0; temp_KV[KV_to_here + 1] = 0.0; temp_KV[KV_to_here + 2] = 0.0; temp_KV[KV_to_here + 3] = 0.5; temp_KV[KV_to_here + 4] = 1.0; temp_KV[KV_to_here + 5] = 1.0; temp_KV[KV_to_here + 6] = 1.0; KV_to_here += 7;
        temp_KV[KV_to_here] = 0.0; temp_KV[KV_to_here + 1] = 0.0; temp_KV[KV_to_here + 2] = 0.0; temp_KV[KV_to_here + 3] = 1.0; temp_KV[KV_to_here + 4] = 1.0; temp_KV[KV_to_here + 5] = 1.0; KV_to_here += 6;

        // パッチ1
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 2.0); temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0); temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0] + singular_width; temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 4.0); temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0); temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 2.0); temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 8
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 9
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0); temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 10
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 11

        temp_KV[KV_to_here] = 0.0; temp_KV[KV_to_here + 1] = 0.0; temp_KV[KV_to_here + 2] = 0.0; temp_KV[KV_to_here + 3] = 0.5; temp_KV[KV_to_here + 4] = 1.0; temp_KV[KV_to_here + 5] = 1.0; temp_KV[KV_to_here + 6] = 1.0; KV_to_here += 7;
        temp_KV[KV_to_here] = 0.0; temp_KV[KV_to_here + 1] = 0.0; temp_KV[KV_to_here + 2] = 0.0; temp_KV[KV_to_here + 3] = 1.0; temp_KV[KV_to_here + 4] = 1.0; temp_KV[KV_to_here + 5] = 1.0; KV_to_here += 6;

        // パッチ2
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0); temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 4.0); temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0); temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 2.0); temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 8
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 9
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 2.0); temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0); temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 10
        temp_CP[CP_to_here] = crack_tip[0] - singular_width; temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 11

        temp_KV[KV_to_here] = 0.0; temp_KV[KV_to_here + 1] = 0.0; temp_KV[KV_to_here + 2] = 0.0; temp_KV[KV_to_here + 3] = 0.5; temp_KV[KV_to_here + 4] = 1.0; temp_KV[KV_to_here + 5] = 1.0; temp_KV[KV_to_here + 6] = 1.0; KV_to_here += 7;
        temp_KV[KV_to_here] = 0.0; temp_KV[KV_to_here + 1] = 0.0; temp_KV[KV_to_here + 2] = 0.0; temp_KV[KV_to_here + 3] = 1.0; temp_KV[KV_to_here + 4] = 1.0; temp_KV[KV_to_here + 5] = 1.0; KV_to_here += 6;

        // パッチ3
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 2.0); temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0); temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0] - singular_width; temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 2.0); temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 4.0); temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0] - singular_width; temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0); temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 8
        temp_CP[CP_to_here] = crack_tip[0]; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 9
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 2.0); temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 10
        temp_CP[CP_to_here] = crack_tip[0] - singular_width; temp_CP[CP_to_here + 1] = crack_tip[1]; temp_CP[CP_to_here + 2] = 1.0; CP_to_here += 3; // point 11

        temp_KV[KV_to_here] = 0.0; temp_KV[KV_to_here + 1] = 0.0; temp_KV[KV_to_here + 2] = 0.0; temp_KV[KV_to_here + 3] = 0.5; temp_KV[KV_to_here + 4] = 1.0; temp_KV[KV_to_here + 5] = 1.0; temp_KV[KV_to_here + 6] = 1.0; KV_to_here += 7;
        temp_KV[KV_to_here] = 0.0; temp_KV[KV_to_here + 1] = 0.0; temp_KV[KV_to_here + 2] = 0.0; temp_KV[KV_to_here + 3] = 1.0; temp_KV[KV_to_here + 4] = 1.0; temp_KV[KV_to_here + 5] = 1.0; KV_to_here += 6;
    }

    int singular_patch_n, Total_patch_n;
    if (mode[0] == 0)
    {
        singular_patch_n = 4;
        Total_patch_n = 10;
    }
    else if (mode[0] == 1)
    {
        singular_patch_n = 16;
        Total_patch_n = 40;
    }

    for (i = 0; i < Total_patch_n - singular_patch_n; i++)
    {
        Get_input_b(temp_argv[i + 2], i, temp_CP, temp_KV, temp_CP_info, temp_KV_info);
    }
}


void Make_CP_KV_full_model(double *temp_CP, double *temp_KV, int *temp_CP_info, int *temp_KV_info)
{
    
}