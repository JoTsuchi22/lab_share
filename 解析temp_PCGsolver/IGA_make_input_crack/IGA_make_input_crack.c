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

#define DIMENSION 2            // 2次元
#define MERGE_DISTANCE 1.0e-13 // コントロールポイントが同じ点と判定する距離

FILE *fp;

static double E_and_nu[2];
static int mode[2]; // 0: analysis model mode, 1: sigular patch mode
static double crack_tip[DIMENSION];
static double singular_width;
static double length_each_side[3]; // (outer side, inner side, upper side)
static double global_length[2];    // (global width, global hight)
static int cp_after_each_side[4];  // (outer side, inner side, upper side, global side)
static double affine[3];
static int counter;
static int CP_result_to_here;
static int CP_to_here_counter;
static double Edge_L[4];
static double Edge_delta[4];

void Get_input_a(char *filename);
void Get_input_b(char *filename, int num, double *temp_CP, double *temp_KV, int *temp_CP_info, int *temp_KV_info);
void Make_patch_file(int num);
void Make_model(char **temp_argv);
void Make_CP_KV_full_model(double *temp_CP, double *temp_KV, int *temp_CP_info, int *temp_KV_info, char **temp_argv);
void Make_B(int num, double *temp_B, int *temp_CP_info, double *temp_CP);
void Check_B(int num_own, int num_opponent, double *temp_B, int *temp_Edge_info, int *temp_Opponent_patch_num);
void Make_connectivity(int num, int *temp_CP_info, int *temp_Edge_info, int *temp_Opponent_patch_num, int *temp_Connectivity, int *temp_A, double *temp_CP, double *temp_CP_result);
void Output_inputdata(int *temp_Order, int *temp_KV_info, int *temp_CP_info, int *temp_Connectivity, double *temp_KV, double *temp_CP_result);
void Output_J(int *temp_CP_info, int *temp_A);
void Output_SVG(double *temp_B, double *temp_CP_result);
void swap(int *a, int *b);
int getLeft(int parent);
int getRight(int parent);
int getParent(int child);
void addHeap(int *a, int size);
void removeHeap(int *a, int size);
void makeHeap(int *a, int num);
void heapSort(int *a, int num);
void Dedupe(int *a, int *num, int *a_new, int *num_new, int n);
double rot_rev(double x, double y, double bool);
double rot(double x, double y, double bool);

int main(int argc, char **argv)
{
    int i;

    // ファイル読み込み
    Get_input_a(argv[1]);

    if (argc == 2)
    {
        // エラー出力
        if (fabs(crack_tip[0] - length_each_side[1]) >= MERGE_DISTANCE)
        {
            printf("Error\n inner side must match tip crack coordinate x\n");
            exit(1);
        }

        // 各パッチ作成
        int singular_patch_n, Total_patch_n;
        singular_patch_n = 16;
        Total_patch_n = 40 + 24;

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
    mode[0] = 1;
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

    // global_length
    for (i = 0; i < 2; i++)
    {
        fscanf(fp, "%lf", &temp_d);
        global_length[i] = temp_d;
        printf("global_length[%d] = %le\n", i, length_each_side[i]);
    }

    fgets(s, 256, fp);

    // cp_after_each_side
    for (i = 0; i < 4; i++)
    {
        fscanf(fp, "%d", &temp_i);
        cp_after_each_side[i] = temp_i;
        printf("cp_after_each_side[%d] = %d\n", i, cp_after_each_side[i]);
    }

    fgets(s, 256, fp);

    // affine infomation
    affine[0] = 0.0;
    affine[1] = 0.0;
    fscanf(fp, "%lf", &temp_d);
    affine[2] = temp_d;
    printf("affine_rot = %le\n", affine[2]);

    fclose(fp);
}

void Get_input_b(char *filename, int num, double *temp_CP, double *temp_KV, int *temp_CP_info, int *temp_KV_info)
{
    int i, j;
    char s[256];

    int temp_i;
    double temp_d;

    int singular_patch_n = 16;

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
    double a, b;
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
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] + singular_width - (2.0 * crack_tip[0]), crack_tip[1], 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] + length_each_side[1] - (2.0 * crack_tip[0]), crack_tip[1], 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] + singular_width - (2.0 * crack_tip[0]), crack_tip[1] + singular_width, 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] + length_each_side[1] - (2.0 * crack_tip[0]), crack_tip[1] + singular_width, 1.0);
    }
    else if (num == 29)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] + singular_width - (2.0 * crack_tip[0]), crack_tip[1] + singular_width, 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] + length_each_side[1] - (2.0 * crack_tip[0]), crack_tip[1] + singular_width, 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] + singular_width - (2.0 * crack_tip[0]), crack_tip[1] + length_each_side[2], 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] + length_each_side[1] - (2.0 * crack_tip[0]), crack_tip[1] + length_each_side[2], 1.0);
    }
    else if (num == 30)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] - (2.0 * crack_tip[0]), crack_tip[1] + singular_width, 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] + singular_width - (2.0 * crack_tip[0]), crack_tip[1] + singular_width, 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] - (2.0 * crack_tip[0]), crack_tip[1] + length_each_side[2], 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] + singular_width - (2.0 * crack_tip[0]), crack_tip[1] + length_each_side[2], 1.0);
    }
    else if (num == 31)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] - singular_width - (2.0 * crack_tip[0]), crack_tip[1] + singular_width, 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] - (2.0 * crack_tip[0]), crack_tip[1] + singular_width, 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] - singular_width - (2.0 * crack_tip[0]), crack_tip[1] + length_each_side[2], 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] - (2.0 * crack_tip[0]), crack_tip[1] + length_each_side[2], 1.0);
    }
    else if (num == 32)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] - length_each_side[0] - (2.0 * crack_tip[0]), crack_tip[1] + singular_width, 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] - singular_width - (2.0 * crack_tip[0]), crack_tip[1] + singular_width, 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] - length_each_side[0] - (2.0 * crack_tip[0]), crack_tip[1] + length_each_side[2], 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] - singular_width - (2.0 * crack_tip[0]), crack_tip[1] + length_each_side[2], 1.0);
    }
    else if (num == 33)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] - length_each_side[0] - (2.0 * crack_tip[0]), crack_tip[1], 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] - singular_width - (2.0 * crack_tip[0]), crack_tip[1], 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] - length_each_side[0] - (2.0 * crack_tip[0]), crack_tip[1] + singular_width, 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] - singular_width - (2.0 * crack_tip[0]), crack_tip[1] + singular_width, 1.0);
    }
    else if (num == 34)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] - length_each_side[0] - (2.0 * crack_tip[0]), crack_tip[1] - singular_width, 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] - singular_width - (2.0 * crack_tip[0]), crack_tip[1] - singular_width, 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] - length_each_side[0] - (2.0 * crack_tip[0]), crack_tip[1], 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] - singular_width - (2.0 * crack_tip[0]), crack_tip[1], 1.0);
    }
    else if (num == 35)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] - length_each_side[0] - (2.0 * crack_tip[0]), crack_tip[1] - length_each_side[2], 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] - singular_width - (2.0 * crack_tip[0]), crack_tip[1] - length_each_side[2], 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] - length_each_side[0] - (2.0 * crack_tip[0]), crack_tip[1] - singular_width, 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] - singular_width - (2.0 * crack_tip[0]), crack_tip[1] - singular_width, 1.0);
    }
    else if (num == 36)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] - singular_width - (2.0 * crack_tip[0]), crack_tip[1] - length_each_side[2], 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] - (2.0 * crack_tip[0]), crack_tip[1] - length_each_side[2], 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] - singular_width - (2.0 * crack_tip[0]), crack_tip[1] - singular_width, 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] - (2.0 * crack_tip[0]), crack_tip[1] - singular_width, 1.0);
    }
    else if (num == 37)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] - (2.0 * crack_tip[0]), crack_tip[1] - length_each_side[2], 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] + singular_width - (2.0 * crack_tip[0]), crack_tip[1] - length_each_side[2], 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] - (2.0 * crack_tip[0]), crack_tip[1] - singular_width, 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] + singular_width - (2.0 * crack_tip[0]), crack_tip[1] - singular_width, 1.0);
    }
    else if (num == 38)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] + singular_width - (2.0 * crack_tip[0]), crack_tip[1] - length_each_side[2], 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] + length_each_side[1] - (2.0 * crack_tip[0]), crack_tip[1] - length_each_side[2], 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] + singular_width - (2.0 * crack_tip[0]), crack_tip[1] - singular_width, 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] + length_each_side[1] - (2.0 * crack_tip[0]), crack_tip[1] - singular_width, 1.0);
    }
    else if (num == 39)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] + singular_width - (2.0 * crack_tip[0]), crack_tip[1] - singular_width, 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, crack_tip[0] + length_each_side[1] - (2.0 * crack_tip[0]), crack_tip[1] - singular_width, 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] + singular_width - (2.0 * crack_tip[0]), crack_tip[1], 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n\n", 3, crack_tip[0] + length_each_side[1] - (2.0 * crack_tip[0]), crack_tip[1], 1.0);
    }
    else if (num == 40)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] + length_each_side[0], crack_tip[1] - length_each_side[2], 1.0);
        a = global_length[0] / 2.0;
        b = rot(crack_tip[0] + length_each_side[0], crack_tip[1] - length_each_side[2], 1);
        Edge_L[0] = (global_length[1] / 2.0) - b;
        Edge_delta[0] = b;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] + length_each_side[0], crack_tip[1] - singular_width, 1.0);
        a = global_length[0] / 2.0;
        b = Edge_delta[0] + Edge_L[0] * (cp_after_each_side[2]) / (cp_after_each_side[2] + singular_width);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 3, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
    }
    else if (num == 41)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] + length_each_side[0], crack_tip[1] - singular_width, 1.0);
        a = global_length[0] / 2.0;
        b = Edge_delta[0] + Edge_L[0] * (cp_after_each_side[2]) / (cp_after_each_side[2] + singular_width);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] + length_each_side[0], crack_tip[1], 1.0);
        a = global_length[0] / 2.0;
        b = global_length[1] / 2.0;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 3, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
    }
    else if (num == 42)
    {
        Edge_L[1] = (global_length[0] / 2.0) - rot(crack_tip[0] + length_each_side[0], crack_tip[1] + length_each_side[2], 0);
        Edge_delta[1] = rot(crack_tip[0] + length_each_side[0], crack_tip[1] + length_each_side[2], 0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] + length_each_side[0], crack_tip[1], 1.0);
        a = global_length[0] / 2.0;
        b = global_length[1] / 2.0;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] + length_each_side[0], crack_tip[1] + singular_width, 1.0);
        a = Edge_delta[1] + Edge_L[1] * (cp_after_each_side[2]) / (cp_after_each_side[2] + singular_width);
        b = global_length[1] / 2.0;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 3, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
    }
    else if (num == 43)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] + length_each_side[0], crack_tip[1] + singular_width, 1.0);
        a = Edge_delta[1] + Edge_L[1] * (cp_after_each_side[2]) / (cp_after_each_side[2] + singular_width);
        b = global_length[1] / 2.0;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] + length_each_side[0], crack_tip[1] + length_each_side[2], 1.0);
        a = rot(crack_tip[0] + length_each_side[0], crack_tip[1] + length_each_side[2], 0);
        b = global_length[1] / 2.0;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 3, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
    }
    else if (num == 44)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] + length_each_side[0], crack_tip[1] + length_each_side[2], 1.0);
        a = rot(crack_tip[0] + length_each_side[0], crack_tip[1] + length_each_side[2], 0);
        b = global_length[1] / 2.0;
        Edge_L[2] = (global_length[0] / 2.0) + a; 
        Edge_delta[2] = a;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] + singular_width, crack_tip[1] + length_each_side[2], 1.0);
        a = Edge_delta[2] - Edge_L[2] * (cp_after_each_side[0]) / (cp_after_each_side[0] + cp_after_each_side[1] + 2 * singular_width);
        b = global_length[1] / 2.0;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 3, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
    }
    else if (num == 45)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] + singular_width, crack_tip[1] + length_each_side[2], 1.0);
        a = Edge_delta[2] - Edge_L[2] * (cp_after_each_side[0]) / (cp_after_each_side[0] + cp_after_each_side[1] + 2.0 * singular_width);
        b = global_length[1] / 2.0;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0], crack_tip[1] + length_each_side[2], 1.0);
        a = Edge_delta[2] - Edge_L[2] * (cp_after_each_side[0] + singular_width) / (cp_after_each_side[0] + cp_after_each_side[1] + 2.0 * singular_width);
        b = global_length[1] / 2.0;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 3, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
    }
    else if (num == 46)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0], crack_tip[1] + length_each_side[2], 1.0);
        a = Edge_delta[2] - Edge_L[2] * (cp_after_each_side[0] + singular_width) / (cp_after_each_side[0] + cp_after_each_side[1] + 2.0 * singular_width);
        b = global_length[1] / 2.0;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0] - singular_width, crack_tip[1] + length_each_side[2], 1.0);
        a = Edge_delta[2] - Edge_L[2] * (cp_after_each_side[0] + 2.0 * singular_width) / (cp_after_each_side[0] + cp_after_each_side[1] + 2.0 * singular_width);
        b = global_length[1] / 2.0;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 3, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
    }
    else if (num == 47)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0] - singular_width, crack_tip[1] + length_each_side[2], 1.0);
        a = Edge_delta[2] - Edge_L[2] * (cp_after_each_side[0] + 2.0 * singular_width) / (cp_after_each_side[0] + cp_after_each_side[1] + 2.0 * singular_width);
        b = global_length[1] / 2.0;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, 0.0, crack_tip[1] + length_each_side[2], 1.0);
        a = - global_length[0] / 2.0;
        b = global_length[1] / 2.0;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 3, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
    }
    else if (num == 48)
    {
        Edge_L[3] = (global_length[1] / 2.0) - rot(- (crack_tip[0] + length_each_side[0]), - (crack_tip[1] - length_each_side[2]), 1);
        Edge_delta[3] = - rot(- (crack_tip[0] + length_each_side[0]), - (crack_tip[1] - length_each_side[2]), 1);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, 0.0, crack_tip[1] + length_each_side[2], 1.0);
        a = - global_length[0] / 2.0;
        b = global_length[1] / 2.0;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, - crack_tip[0] + singular_width, crack_tip[1] + length_each_side[2], 1.0);
        a = - global_length[0] / 2.0;
        b = - Edge_delta[3] + Edge_L[3] * (cp_after_each_side[0] + 2.0 * singular_width) / (cp_after_each_side[0] + cp_after_each_side[1] + 2.0 * singular_width);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 3, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
    }
    else if (num == 49)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, - crack_tip[0] + singular_width, crack_tip[1] + length_each_side[2], 1.0);
        a = - global_length[0] / 2.0;
        b = - Edge_delta[3] + Edge_L[3] * (cp_after_each_side[0] + 2.0 * singular_width) / (cp_after_each_side[0] + cp_after_each_side[1] + 2.0 * singular_width);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, - crack_tip[0], crack_tip[1] + length_each_side[2], 1.0);
        a = - global_length[0] / 2.0;
        b = - Edge_delta[3] + Edge_L[3] * (cp_after_each_side[0] + singular_width) / (cp_after_each_side[0] + cp_after_each_side[1] + 2.0 * singular_width);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 3, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
    }
    else if (num == 50)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, - crack_tip[0], crack_tip[1] + length_each_side[2], 1.0);
        a = - global_length[0] / 2.0;
        b = - Edge_delta[3] + Edge_L[3] * (cp_after_each_side[0] + singular_width) / (cp_after_each_side[0] + cp_after_each_side[1] + 2.0 * singular_width);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, - crack_tip[0] - singular_width, crack_tip[1] + length_each_side[2], 1.0);
        a = - global_length[0] / 2.0;
        b = - Edge_delta[3] + Edge_L[3] * (cp_after_each_side[0]) / (cp_after_each_side[0] + cp_after_each_side[1] + 2.0 * singular_width);;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 3, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
    }
    else if (num == 51)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, - crack_tip[0] - singular_width, crack_tip[1] + length_each_side[2], 1.0);
        a = - global_length[0] / 2.0;
        b = - Edge_delta[3] + Edge_L[3] * (cp_after_each_side[0]) / (cp_after_each_side[0] + cp_after_each_side[1] + 2.0 * singular_width);;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, - crack_tip[0] - length_each_side[0], crack_tip[1] + length_each_side[2], 1.0);
        a = - global_length[0] / 2.0;
        b = rot(- (crack_tip[0] + length_each_side[0]), - (crack_tip[1] - length_each_side[2]), 1);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 3, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
    }
    else if (num == 52)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, - (crack_tip[0] + length_each_side[0]), - (crack_tip[1] - length_each_side[2]), 1.0);
        a = - (global_length[0] / 2.0);
        b = rot(- (crack_tip[0] + length_each_side[0]), - (crack_tip[1] - length_each_side[2]), 1);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, - (crack_tip[0] + length_each_side[0]), - (crack_tip[1] - singular_width), 1.0);
        a = - (global_length[0] / 2.0);
        b = - Edge_delta[0] - Edge_L[0] * (cp_after_each_side[2]) / (cp_after_each_side[2] + singular_width);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 3, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
    }
    else if (num == 53)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, - (crack_tip[0] + length_each_side[0]), - (crack_tip[1] - singular_width), 1.0);
        a = - global_length[0] / 2.0;
        b = - Edge_delta[0] - Edge_L[0] * (cp_after_each_side[2]) / (cp_after_each_side[2] + singular_width);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, - (crack_tip[0] + length_each_side[0]), - crack_tip[1], 1.0);
        a = - global_length[0] / 2.0;
        b = - global_length[1] / 2.0;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 3, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
    }
    else if (num == 54)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, - (crack_tip[0] + length_each_side[0]), - crack_tip[1], 1.0);
        a = - global_length[0] / 2.0;
        b = - global_length[1] / 2.0;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, - (crack_tip[0] + length_each_side[0]), - (crack_tip[1] + singular_width), 1.0);
        a = - Edge_delta[1] - Edge_L[1] * (cp_after_each_side[2]) / (cp_after_each_side[2] + singular_width);
        b = - global_length[1] / 2.0;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 3, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
    }
    else if (num == 55)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, - (crack_tip[0] + length_each_side[0]), - (crack_tip[1] + singular_width), 1.0);
        a = - Edge_delta[1] - Edge_L[1] * (cp_after_each_side[2]) / (cp_after_each_side[2] + singular_width);
        b = - global_length[1] / 2.0;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, - (crack_tip[0] + length_each_side[0]), - (crack_tip[1] + length_each_side[2]), 1.0);
        a = rot(- (crack_tip[0] + length_each_side[0]), - (crack_tip[1] + length_each_side[2]), 0);
        b = - global_length[1] / 2.0;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 3, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
    }
    else if (num == 56)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, - (crack_tip[0] + length_each_side[0]), - (crack_tip[1] + length_each_side[2]), 1.0);
        a = rot(- (crack_tip[0] + length_each_side[0]), - (crack_tip[1] + length_each_side[2]), 0);
        b = - global_length[1] / 2.0;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, - (crack_tip[0] + singular_width), - (crack_tip[1] + length_each_side[2]), 1.0);
        a = - Edge_delta[2] + Edge_L[2] * (cp_after_each_side[0]) / (cp_after_each_side[0] + cp_after_each_side[1] + 2 * singular_width);
        b = - global_length[1] / 2.0;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 3, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
    }
    else if (num == 57)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, - (crack_tip[0] + singular_width), - (crack_tip[1] + length_each_side[2]), 1.0);
        a = - Edge_delta[2] + Edge_L[2] * (cp_after_each_side[0]) / (cp_after_each_side[0] + cp_after_each_side[1] + 2.0 * singular_width);
        b = - global_length[1] / 2.0;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, - crack_tip[0], - (crack_tip[1] + length_each_side[2]), 1.0);
        a = - Edge_delta[2] + Edge_L[2] * (cp_after_each_side[0] + singular_width) / (cp_after_each_side[0] + cp_after_each_side[1] + 2.0 * singular_width);
        b = - global_length[1] / 2.0;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 3, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
    }
    else if (num == 58)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, - crack_tip[0], - (crack_tip[1] + length_each_side[2]), 1.0);
        a = - Edge_delta[2] + Edge_L[2] * (cp_after_each_side[0] + singular_width) / (cp_after_each_side[0] + cp_after_each_side[1] + 2.0 * singular_width);
        b = - global_length[1] / 2.0;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, - (crack_tip[0] - singular_width), - (crack_tip[1] + length_each_side[2]), 1.0);
        a = - Edge_delta[2] + Edge_L[2] * (cp_after_each_side[0] + 2.0 * singular_width) / (cp_after_each_side[0] + cp_after_each_side[1] + 2.0 * singular_width);
        b = - global_length[1] / 2.0;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 3, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
    }
    else if (num == 59)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, - (crack_tip[0] - singular_width), - (crack_tip[1] + length_each_side[2]), 1.0);
        a = - Edge_delta[2] + Edge_L[2] * (cp_after_each_side[0] + 2.0 * singular_width) / (cp_after_each_side[0] + cp_after_each_side[1] + 2.0 * singular_width);
        b = - global_length[1] / 2.0;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, 0.0, - (crack_tip[1] + length_each_side[2]), 1.0);
        a = global_length[0] / 2.0;
        b = - global_length[1] / 2.0;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 3, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
    }
    else if (num == 60)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, 0.0, - (crack_tip[1] + length_each_side[2]), 1.0);
        a = global_length[0] / 2.0;
        b = - global_length[1] / 2.0;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, - (- crack_tip[0] + singular_width), - (crack_tip[1] + length_each_side[2]), 1.0);
        a = global_length[0] / 2.0;
        b = Edge_delta[3] - Edge_L[3] * (cp_after_each_side[0] + 2.0 * singular_width) / (cp_after_each_side[0] + cp_after_each_side[1] + 2.0 * singular_width);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 3, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
    }
    else if (num == 61)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, - (- crack_tip[0] + singular_width), - (crack_tip[1] + length_each_side[2]), 1.0);
        a = global_length[0] / 2.0;
        b = Edge_delta[3] - Edge_L[3] * (cp_after_each_side[0] + 2.0 * singular_width) / (cp_after_each_side[0] + cp_after_each_side[1] + 2.0 * singular_width);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, crack_tip[0], - (crack_tip[1] + length_each_side[2]), 1.0);
        a = global_length[0] / 2.0;
        b = Edge_delta[3] - Edge_L[3] * (cp_after_each_side[0] + singular_width) / (cp_after_each_side[0] + cp_after_each_side[1] + 2.0 * singular_width);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 3, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
    }
    else if (num == 62)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, crack_tip[0], - (crack_tip[1] + length_each_side[2]), 1.0);
        a = global_length[0] / 2.0;
        b = Edge_delta[3] - Edge_L[3] * (cp_after_each_side[0] + singular_width) / (cp_after_each_side[0] + cp_after_each_side[1] + 2.0 * singular_width);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, - (- crack_tip[0] - singular_width), - (crack_tip[1] + length_each_side[2]), 1.0);
        a = global_length[0] / 2.0;
        b = Edge_delta[3] - Edge_L[3] * (cp_after_each_side[0]) / (cp_after_each_side[0] + cp_after_each_side[1] + 2.0 * singular_width);;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 3, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
    }
    else if (num == 63)
    {
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 0, - (- crack_tip[0] - singular_width), - (crack_tip[1] + length_each_side[2]), 1.0);
        a = global_length[0] / 2.0;
        b = Edge_delta[3] - Edge_L[3] * (cp_after_each_side[0]) / (cp_after_each_side[0] + cp_after_each_side[1] + 2.0 * singular_width);;
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 1, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 2, - (- crack_tip[0] - length_each_side[0]), - (crack_tip[1] + length_each_side[2]), 1.0);
        a = global_length[0] / 2.0;
        b = rot(crack_tip[0] + length_each_side[0], crack_tip[1] - length_each_side[2], 1);
        fprintf(fp, "%d   %.17e   %.17e   %.17e\n", 3, rot_rev(a, b, 0), rot_rev(a, b, 1), 1.0);
    }
    fprintf(fp, "\n");

    // o.e.
    fprintf(fp, "%d   %d\n\n", 1, 1);

    // k.i.
    fprintf(fp, "%d   %d\n\n", 0, 0);

    // パッチごとのxi方向のコントロールポイント数
    if (num == 16 || num == 17 || num == 26 || num == 27 || num == 32 || num == 33 || num == 34 || num == 35)
    {
        fprintf(fp, "%d   ", cp_after_each_side[0]);
    }
    else if (num == 18 || num == 19 || num == 24 || num == 25 || num == 30 || num == 31 || num == 36 || num == 37)
    {
        fprintf(fp, "%d   ", 3);
    }
    else if (num == 20 || num == 21 || num == 22 || num == 23 || num == 28 || num == 29 || num == 38 || num == 39)
    {
        fprintf(fp, "%d   ", cp_after_each_side[1]);
    }
    else if (num >= 40 && num <= 63)
    {
        fprintf(fp, "%d   ", cp_after_each_side[3]);
    }

    // パッチごとのeta方向のコントロールポイント数
    if (num == 16 || num == 21 || num == 22 || num == 27 || num == 28 || num == 33 || num == 34 || num == 39)
    {
        fprintf(fp, "%d\n\n", 3);
    }
    else if (num == 17 || num == 18 || num == 19 || num == 20 || num == 23 || num == 24 || num == 25 || num == 26 || num == 29 || num == 30 || num == 31 || num == 32 || num == 35 || num == 36 || num == 37 || num == 38)
    {
        fprintf(fp, "%d\n\n", cp_after_each_side[2]);
    }
    else if (num == 41 || num == 42 || num == 45 || num == 46 || num == 49 || num == 50 || num == 53 || num == 54 || num == 57 || num == 58 || num == 61 || num == 62)
    {
        fprintf(fp, "%d\n\n", 3);
    }
    else if (num == 44 || num == 51 || num == 56 || num == 63)
    {
        fprintf(fp, "%d\n\n", cp_after_each_side[0]); // outer
    }
    else if (num == 47 || num == 48 || num == 59 || num == 60)
    {
        fprintf(fp, "%d\n\n", cp_after_each_side[1]); // inner
    }
    else if (num == 40 || num == 43 || num == 52 || num == 55)
    {
        fprintf(fp, "%d\n\n", cp_after_each_side[2]); // upper
    }
    fprintf(fp, "%d   %d\n\n", 0, 0);

    fclose(fp);
}

void Make_model(char **temp_argv)
{
    int i, j;
    int Total_patch_loc = 0;
    int singular_patch_n = 0;

    Total_patch_loc = 40 + 24;
    singular_patch_n = 16;

    // 動的メモリ確保
    int *Order = (int *)malloc(sizeof(int) * Total_patch_loc * DIMENSION);   // int Order[パッチ番号][DIMENSION]
    int *KV_info = (int *)malloc(sizeof(int) * Total_patch_loc * DIMENSION); // int KV_info[パッチ番号][DIMENSION]
    int *CP_info = (int *)calloc(Total_patch_loc * DIMENSION, sizeof(int));  // int CP_info[パッチ番号][DIMENSION]

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
            Order[i * DIMENSION + j] = 2; // 基底関数は2次で固定
            if (j == 0)                   // xi方向
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
            else if (j == 1) // eta方向
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
            Order[i * DIMENSION + j] = 2; // 基底関数は2次で固定
        }

        // パッチごとのxi方向のコントロールポイント数
        if (i == 16 || i == 17 || i == 26 || i == 27 || i == 32 || i == 33 || i == 34 || i == 35)
        {
            CP_info[i * DIMENSION + 0] = cp_after_each_side[0];
        }
        else if (i == 18 || i == 19 || i == 24 || i == 25 || i == 30 || i == 31 || i == 36 || i == 37)
        {
            CP_info[i * DIMENSION + 0] = 3;
        }
        else if (i == 20 || i == 21 || i == 22 || i == 23 || i == 28 || i == 29 || i == 38 || i == 39)
        {
            CP_info[i * DIMENSION + 0] = cp_after_each_side[1];
        }
        else if (i >= 40 && i <= 63)
        {
            CP_info[i * DIMENSION + 0] = cp_after_each_side[3];
        }

        // パッチごとのeta方向のコントロールポイント数
        if (i == 16 || i == 21 || i == 22 || i == 27 || i == 28 || i == 33 || i == 34 || i == 39)
        {
            CP_info[i * DIMENSION + 1] = 3;
        }
        else if (i == 17 || i == 18 || i == 19 || i == 20 || i == 23 || i == 24 || i == 25 || i == 26 || i == 29 || i == 30 || i == 31 || i == 32 || i == 35 || i == 36 || i == 37 || i == 38)
        {
            CP_info[i * DIMENSION + 1] = cp_after_each_side[2];
        }
        else if (i == 41 || i == 42 || i == 45 || i == 46 || i == 49 || i == 50 || i == 53 || i == 54 || i == 57 || i == 58 || i == 61 || i == 62)
        {
            CP_info[i * DIMENSION + 1] = 3;
        }
        else if (i == 44 || i == 51 || i == 56 || i == 63)
        {
            CP_info[i * DIMENSION + 1] = cp_after_each_side[0]; // outer
        }
        else if (i == 47 || i == 48 || i == 59 || i == 60)
        {
            CP_info[i * DIMENSION + 1] = cp_after_each_side[1]; // inner
        }
        else if (i == 40 || i == 43 || i == 52 || i == 55)
        {
            CP_info[i * DIMENSION + 1] = cp_after_each_side[2]; // upper
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
    double *CP = (double *)malloc(sizeof(double) * temp1 * (DIMENSION + 1));               // double CP[パッチ番号][パッチ内CP番号][xyw -> 3]
    double *CP_result = (double *)malloc(sizeof(double) * temp1 * (DIMENSION + 1));        // double CP_result[通しのコントロールポイント番号(連番)][xyw -> 3]
    int *A = (int *)malloc(sizeof(int) * temp2);                                           // int    A[パッチ番号][辺番号(0~3)][辺内のコネクティビティ]
    double *B = (double *)malloc(sizeof(double) * Total_patch_loc * 16 * (DIMENSION + 1)); // double B[パッチ番号][辺番号(0~7)(正負方向)][各辺の端の2頂点][座標xyw -> 3]
    int *Connectivity = (int *)malloc(sizeof(int) * temp1);                                // int    Connectivity[パッチ番号][パッチ内CP番号]
    double *KV = (double *)malloc(sizeof(double) * Total_patch_loc * temp3);               // double KV[パッチ番号][DIMENSION][ノットベクトル番号]

    if (CP == NULL || CP_result == NULL || A == NULL || B == NULL || Connectivity == NULL || KV == NULL)
    {
        printf("Memory cannot be allocated\n");
        exit(1);
    }

    // 特異パッチのコントロールポイントとノットベクトルを作成
    Make_CP_KV_full_model(CP, KV, CP_info, KV_info, temp_argv);

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

    // パッチコネクティビティの作成
    for (i = 0; i < Total_patch_loc; i++)
    {
        Make_B(i, B, CP_info, CP);
    }

    counter = 0, CP_to_here_counter = 0;
    for (i = 0; i < Total_patch_loc; i++)
    {
        for (j = 0; j < i; j++)
        {
            Check_B(i, j, B, Edge_info, Opponent_patch_num);
        }
        Make_connectivity(i, CP_info, Edge_info, Opponent_patch_num, Connectivity, A, CP, CP_result);
    }

    // アフィン変換
    if (affine[2] >= MERGE_DISTANCE)
    {
        double *Affine_CP = (double *)malloc(sizeof(double) * CP_result_to_here); // アフィン変換後のコントロールポイント

        if (Affine_CP == NULL)
        {
            printf("Memory cannot be allocated\n");
            exit(1);
        }

        double PI = 3.14159265358979323846264338327950288;
        double theta = affine[2] * PI / 180.0;
        double rot[4] = {cos(theta), -sin(theta), sin(theta), cos(theta)};

        int temp_counter = 0;
        for (i = 0; i < (CP_result_to_here + 1) / 3; i++)
        {
            Affine_CP[temp_counter] = rot[0] * CP_result[temp_counter] + rot[1] * CP_result[temp_counter + 1] + affine[0];
            Affine_CP[temp_counter + 1] = rot[2] * CP_result[temp_counter] + rot[3] * CP_result[temp_counter + 1] + affine[1];
            Affine_CP[temp_counter + 2] = CP_result[temp_counter + 2];
            temp_counter += 3;
        }

        for (i = 0; i < CP_result_to_here; i++)
        {
            CP_result[i] = Affine_CP[i];
        }

        free(Affine_CP);
    }

    // 出力
    Output_inputdata(Order, KV_info, CP_info, Connectivity, KV, CP_result);
    Output_J(CP_info, A);

    // 図の出力
    Output_SVG(B, CP_result); // SVG出力

    // メモリ解放
    free(Order), free(KV_info), free(CP_info);
    free(CP), free(CP_result), free(A), free(B), free(Connectivity), free(KV);
    free(Edge_info), free(Opponent_patch_num);
}

void Make_CP_KV_full_model(double *temp_CP, double *temp_KV, int *temp_CP_info, int *temp_KV_info, char **temp_argv)
{
    int i;
    int CP_to_here = 0;
    int KV_to_here = 0;

    // 特異パッチの要素の種類
    if (mode[1] == 1)
    {
        // パッチ0
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] + singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0] + singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0] + singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;

        // パッチ1
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] + singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 2.0);
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;

        // パッチ2
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 2.0);
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0] - singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;

        // パッチ3
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] - singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0] - singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0] - singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;

        // パッチ4
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] - singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0] - singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1] - (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0] - singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;

        // パッチ5
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] - singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 2.0);
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;

        // パッチ6
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 2.0);
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0] + singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;

        // パッチ7
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] + singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0] + singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1] - (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0] + singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;

        // パッチ8
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] + singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0] + singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0] + singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;

        // パッチ9
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] + singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 2.0) - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;

        // パッチ10
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 2.0) - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0] - singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;

        // パッチ11
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] - singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0] - singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0] - singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;

        // パッチ12
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] - singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0] - singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] - (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0] - singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;

        // パッチ13
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] - singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 2.0) - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;

        // パッチ14
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 2.0) - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0] + singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;

        // パッチ15
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] + singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0] + singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] - (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0] + singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;
    }
    else if (mode[1] == 2)
    {
        // パッチ0
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 2.0);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0] + singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 2.0);
        temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 4.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0] + singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1] + +(singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 9
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 2.0);
        temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 10
        temp_CP[CP_to_here] = crack_tip[0] + singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 11

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 0.5;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        temp_KV[KV_to_here + 6] = 1.0;
        KV_to_here += 7;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;

        // パッチ1
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 2.0);
        temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0] + singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 4.0);
        temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 2.0);
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 9
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 10
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 11

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 0.5;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        temp_KV[KV_to_here + 6] = 1.0;
        KV_to_here += 7;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;

        // パッチ2
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 4.0);
        temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 2.0);
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 9
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 2.0);
        temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 10
        temp_CP[CP_to_here] = crack_tip[0] - singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 11

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 0.5;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        temp_KV[KV_to_here + 6] = 1.0;
        KV_to_here += 7;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;

        // パッチ3
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 2.0);
        temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0] - singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 2.0);
        temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 4.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0] - singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 9
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 2.0);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 10
        temp_CP[CP_to_here] = crack_tip[0] - singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 11

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 0.5;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        temp_KV[KV_to_here + 6] = 1.0;
        KV_to_here += 7;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;

        // パッチ4
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 2.0);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0] - singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 2.0);
        temp_CP[CP_to_here + 1] = crack_tip[1] - (singular_width / 4.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0] - singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1] - (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 9
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 2.0);
        temp_CP[CP_to_here + 1] = crack_tip[1] - (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 10
        temp_CP[CP_to_here] = crack_tip[0] - singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 11

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 0.5;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        temp_KV[KV_to_here + 6] = 1.0;
        KV_to_here += 7;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;

        // パッチ5
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 2.0);
        temp_CP[CP_to_here + 1] = crack_tip[1] - (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0] - singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 4.0);
        temp_CP[CP_to_here + 1] = crack_tip[1] - (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 2.0);
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 9
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1] - (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 10
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 11

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 0.5;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        temp_KV[KV_to_here + 6] = 1.0;
        KV_to_here += 7;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;

        // パッチ6
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1] - (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 4.0);
        temp_CP[CP_to_here + 1] = crack_tip[1] - (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 2.0);
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 9
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 2.0);
        temp_CP[CP_to_here + 1] = crack_tip[1] - (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 10
        temp_CP[CP_to_here] = crack_tip[0] + singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 11

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 0.5;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        temp_KV[KV_to_here + 6] = 1.0;
        KV_to_here += 7;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;

        // パッチ7
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 2.0);
        temp_CP[CP_to_here + 1] = crack_tip[1] - (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0] + singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 2.0);
        temp_CP[CP_to_here + 1] = crack_tip[1] - (singular_width / 4.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0] + singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1] - (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8
        temp_CP[CP_to_here] = crack_tip[0];
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 9
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 2.0);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 10
        temp_CP[CP_to_here] = crack_tip[0] + singular_width;
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 11

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 0.5;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        temp_KV[KV_to_here + 6] = 1.0;
        KV_to_here += 7;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;

        // パッチ8
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 2.0) - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0] + singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 2.0) - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 4.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0] + singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] + +(singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 9
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 2.0) - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 10
        temp_CP[CP_to_here] = crack_tip[0] + singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 11

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 0.5;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        temp_KV[KV_to_here + 6] = 1.0;
        KV_to_here += 7;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;

        // パッチ9
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 2.0) - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0] + singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 4.0) - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 2.0) - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 9
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 10
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 11

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 0.5;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        temp_KV[KV_to_here + 6] = 1.0;
        KV_to_here += 7;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;

        // パッチ10
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 4.0) - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 2.0) - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 9
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 2.0) - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 10
        temp_CP[CP_to_here] = crack_tip[0] - singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 11

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 0.5;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        temp_KV[KV_to_here + 6] = 1.0;
        KV_to_here += 7;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;

        // パッチ11
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 2.0) - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0] - singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] + singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 2.0) - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 4.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0] - singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] + (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 9
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 2.0) - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 10
        temp_CP[CP_to_here] = crack_tip[0] - singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 11

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 0.5;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        temp_KV[KV_to_here + 6] = 1.0;
        KV_to_here += 7;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;

        // パッチ12
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 2.0) - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0] - singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 2.0) - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] - (singular_width / 4.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0] - singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] - (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 9
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 2.0) - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] - (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 10
        temp_CP[CP_to_here] = crack_tip[0] - singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 11

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 0.5;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        temp_KV[KV_to_here + 6] = 1.0;
        KV_to_here += 7;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;

        // パッチ13
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 2.0) - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] - (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0] - singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 4.0) - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] - (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0] - (singular_width / 2.0) - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 9
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] - (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 10
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 11

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 0.5;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        temp_KV[KV_to_here + 6] = 1.0;
        KV_to_here += 7;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;

        // パッチ14
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] - (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 4.0) - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] - (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 2.0) - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 9
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 2.0) - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] - (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 10
        temp_CP[CP_to_here] = crack_tip[0] + singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 11

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 0.5;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        temp_KV[KV_to_here + 6] = 1.0;
        KV_to_here += 7;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;

        // パッチ15
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 0
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 1
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 2.0) - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] - (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 2
        temp_CP[CP_to_here] = crack_tip[0] + singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] - singular_width;
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 3
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 4
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 5
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 2.0) - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] - (singular_width / 4.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 6
        temp_CP[CP_to_here] = crack_tip[0] + singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1] - (singular_width / 2.0);
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 7
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 8
        temp_CP[CP_to_here] = crack_tip[0] - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 9
        temp_CP[CP_to_here] = crack_tip[0] + (singular_width / 2.0) - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 10
        temp_CP[CP_to_here] = crack_tip[0] + singular_width - (2.0 * crack_tip[0]);
        temp_CP[CP_to_here + 1] = crack_tip[1];
        temp_CP[CP_to_here + 2] = 1.0;
        CP_to_here += 3; // point 11

        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 0.5;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        temp_KV[KV_to_here + 6] = 1.0;
        KV_to_here += 7;
        temp_KV[KV_to_here] = 0.0;
        temp_KV[KV_to_here + 1] = 0.0;
        temp_KV[KV_to_here + 2] = 0.0;
        temp_KV[KV_to_here + 3] = 1.0;
        temp_KV[KV_to_here + 4] = 1.0;
        temp_KV[KV_to_here + 5] = 1.0;
        KV_to_here += 6;
    }

    int singular_patch_n = 16, Total_patch_n = 40 + 24;

    for (i = 0; i < Total_patch_n - singular_patch_n; i++)
    {
        Get_input_b(temp_argv[i + 2], i, temp_CP, temp_KV, temp_CP_info, temp_KV_info);
    }
}

void Make_B(int num, double *temp_B, int *temp_CP_info, double *temp_CP)
{
    int i;
    int CP_to_here = 0;
    int B_to_here = num * 16 * (DIMENSION + 1);

    for (i = 0; i < num; i++)
    {
        CP_to_here += temp_CP_info[i * DIMENSION] * temp_CP_info[i * DIMENSION + 1] * (DIMENSION + 1);
    }

    // B 配列を作成
    // 辺0 点0
    temp_B[B_to_here] = temp_CP[CP_to_here];
    temp_B[B_to_here + 1] = temp_CP[CP_to_here + 1];
    temp_B[B_to_here + 2] = temp_CP[CP_to_here + 2];
    B_to_here += DIMENSION + 1;

    // 辺0 点1
    temp_B[B_to_here] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] - 1) * (DIMENSION + 1)];
    temp_B[B_to_here + 1] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] - 1) * (DIMENSION + 1) + 1];
    temp_B[B_to_here + 2] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] - 1) * (DIMENSION + 1) + 2];
    B_to_here += DIMENSION + 1;

    // 辺1 点0
    temp_B[B_to_here] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] - 1) * (DIMENSION + 1)];
    temp_B[B_to_here + 1] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] - 1) * (DIMENSION + 1) + 1];
    temp_B[B_to_here + 2] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] - 1) * (DIMENSION + 1) + 2];
    B_to_here += DIMENSION + 1;

    // 辺1 点1
    temp_B[B_to_here] = temp_CP[CP_to_here];
    temp_B[B_to_here + 1] = temp_CP[CP_to_here + 1];
    temp_B[B_to_here + 2] = temp_CP[CP_to_here + 2];
    B_to_here += DIMENSION + 1;

    // 辺2 点0
    temp_B[B_to_here] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] - 1) * (DIMENSION + 1)];
    temp_B[B_to_here + 1] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] - 1) * (DIMENSION + 1) + 1];
    temp_B[B_to_here + 2] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] - 1) * (DIMENSION + 1) + 2];
    B_to_here += DIMENSION + 1;

    // 辺2 点1
    temp_B[B_to_here] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] * temp_CP_info[num * DIMENSION + 1] - 1) * (DIMENSION + 1)];
    temp_B[B_to_here + 1] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] * temp_CP_info[num * DIMENSION + 1] - 1) * (DIMENSION + 1) + 1];
    temp_B[B_to_here + 2] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] * temp_CP_info[num * DIMENSION + 1] - 1) * (DIMENSION + 1) + 2];
    B_to_here += DIMENSION + 1;

    // 辺3 点0
    temp_B[B_to_here] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] * temp_CP_info[num * DIMENSION + 1] - 1) * (DIMENSION + 1)];
    temp_B[B_to_here + 1] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] * temp_CP_info[num * DIMENSION + 1] - 1) * (DIMENSION + 1) + 1];
    temp_B[B_to_here + 2] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] * temp_CP_info[num * DIMENSION + 1] - 1) * (DIMENSION + 1) + 2];
    B_to_here += DIMENSION + 1;

    // 辺3 点1
    temp_B[B_to_here] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] - 1) * (DIMENSION + 1)];
    temp_B[B_to_here + 1] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] - 1) * (DIMENSION + 1) + 1];
    temp_B[B_to_here + 2] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] - 1) * (DIMENSION + 1) + 2];
    B_to_here += DIMENSION + 1;

    // 辺4 点0
    temp_B[B_to_here] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] * (temp_CP_info[num * DIMENSION + 1] - 1)) * (DIMENSION + 1)];
    temp_B[B_to_here + 1] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] * (temp_CP_info[num * DIMENSION + 1] - 1)) * (DIMENSION + 1) + 1];
    temp_B[B_to_here + 2] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] * (temp_CP_info[num * DIMENSION + 1] - 1)) * (DIMENSION + 1) + 2];
    B_to_here += DIMENSION + 1;

    // 辺4 点1
    temp_B[B_to_here] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] * temp_CP_info[num * DIMENSION + 1] - 1) * (DIMENSION + 1)];
    temp_B[B_to_here + 1] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] * temp_CP_info[num * DIMENSION + 1] - 1) * (DIMENSION + 1) + 1];
    temp_B[B_to_here + 2] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] * temp_CP_info[num * DIMENSION + 1] - 1) * (DIMENSION + 1) + 2];
    B_to_here += DIMENSION + 1;

    // 辺5 点0
    temp_B[B_to_here] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] * temp_CP_info[num * DIMENSION + 1] - 1) * (DIMENSION + 1)];
    temp_B[B_to_here + 1] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] * temp_CP_info[num * DIMENSION + 1] - 1) * (DIMENSION + 1) + 1];
    temp_B[B_to_here + 2] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] * temp_CP_info[num * DIMENSION + 1] - 1) * (DIMENSION + 1) + 2];
    B_to_here += DIMENSION + 1;

    // 辺5 点1
    temp_B[B_to_here] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] * (temp_CP_info[num * DIMENSION + 1] - 1)) * (DIMENSION + 1)];
    temp_B[B_to_here + 1] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] * (temp_CP_info[num * DIMENSION + 1] - 1)) * (DIMENSION + 1) + 1];
    temp_B[B_to_here + 2] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] * (temp_CP_info[num * DIMENSION + 1] - 1)) * (DIMENSION + 1) + 2];
    B_to_here += DIMENSION + 1;

    // 辺6 点0
    temp_B[B_to_here] = temp_CP[CP_to_here];
    temp_B[B_to_here + 1] = temp_CP[CP_to_here + 1];
    temp_B[B_to_here + 2] = temp_CP[CP_to_here + 2];
    B_to_here += DIMENSION + 1;

    // 辺6 点1
    temp_B[B_to_here] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] * (temp_CP_info[num * DIMENSION + 1] - 1)) * (DIMENSION + 1)];
    temp_B[B_to_here + 1] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] * (temp_CP_info[num * DIMENSION + 1] - 1)) * (DIMENSION + 1) + 1];
    temp_B[B_to_here + 2] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] * (temp_CP_info[num * DIMENSION + 1] - 1)) * (DIMENSION + 1) + 2];
    B_to_here += DIMENSION + 1;

    // 辺7 点0
    temp_B[B_to_here] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] * (temp_CP_info[num * DIMENSION + 1] - 1)) * (DIMENSION + 1)];
    temp_B[B_to_here + 1] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] * (temp_CP_info[num * DIMENSION + 1] - 1)) * (DIMENSION + 1) + 1];
    temp_B[B_to_here + 2] = temp_CP[CP_to_here + (temp_CP_info[num * DIMENSION] * (temp_CP_info[num * DIMENSION + 1] - 1)) * (DIMENSION + 1) + 2];
    B_to_here += DIMENSION + 1;

    // 辺7 点1
    temp_B[B_to_here] = temp_CP[CP_to_here];
    temp_B[B_to_here + 1] = temp_CP[CP_to_here + 1];
    temp_B[B_to_here + 2] = temp_CP[CP_to_here + 2];
    B_to_here += DIMENSION + 1;

    printf("B on patch %d\n", num);
    printf("point0[x y w] point1[x y w]\n");
    int temp_counter = num * 16 * (DIMENSION + 1);

    for (i = 0; i < 8; i++)
    {
        printf("[%le %le %le]\t[%le %le %le]\n", temp_B[temp_counter], temp_B[temp_counter + 1], temp_B[temp_counter + 2], temp_B[temp_counter + 3], temp_B[temp_counter + 4], temp_B[temp_counter + 5]);
        temp_counter += 6;
    }
}

void Check_B(int num_own, int num_opponent, double *temp_B, int *temp_Edge_info, int *temp_Opponent_patch_num)
{
    int i, j;
    int ii;
    double x_diff[4], y_diff[4], w_diff[4];

    int Check_B_own_to_here = num_own * 16 * (DIMENSION + 1);
    int Check_B_opponent_to_here = num_opponent * 16 * (DIMENSION + 1);

    // printf("自分パッチの先頭x y %le %le\n", temp_B[Check_B_own_to_here], temp_B[Check_B_own_to_here + 1]);
    // printf("相手パッチの先頭x y %le %le\n", temp_B[Check_B_opponent_to_here], temp_B[Check_B_opponent_to_here + 1]);

    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 8; j++)
        {
            ii = 2 * i;
            // 点0
            x_diff[0] = temp_B[Check_B_own_to_here + ii * 2 * (DIMENSION + 1)] - temp_B[Check_B_opponent_to_here + j * 2 * (DIMENSION + 1)];
            y_diff[0] = temp_B[Check_B_own_to_here + ii * 2 * (DIMENSION + 1) + 1] - temp_B[Check_B_opponent_to_here + j * 2 * (DIMENSION + 1) + 1];
            w_diff[0] = temp_B[Check_B_own_to_here + ii * 2 * (DIMENSION + 1) + 2] - temp_B[Check_B_opponent_to_here + j * 2 * (DIMENSION + 1) + 2];

            // 点1
            x_diff[1] = temp_B[Check_B_own_to_here + ii * 2 * (DIMENSION + 1) + (DIMENSION + 1)] - temp_B[Check_B_opponent_to_here + j * 2 * (DIMENSION + 1) + (DIMENSION + 1)];
            y_diff[1] = temp_B[Check_B_own_to_here + ii * 2 * (DIMENSION + 1) + (DIMENSION + 1) + 1] - temp_B[Check_B_opponent_to_here + j * 2 * (DIMENSION + 1) + (DIMENSION + 1) + 1];
            w_diff[1] = temp_B[Check_B_own_to_here + ii * 2 * (DIMENSION + 1) + (DIMENSION + 1) + 2] - temp_B[Check_B_opponent_to_here + j * 2 * (DIMENSION + 1) + (DIMENSION + 1) + 2];

            // 辺0
            x_diff[2] = temp_B[Check_B_own_to_here + ii * 2 * (DIMENSION + 1) + (DIMENSION + 1)] - temp_B[Check_B_own_to_here + ii * 2 * (DIMENSION + 1)];
            y_diff[2] = temp_B[Check_B_own_to_here + ii * 2 * (DIMENSION + 1) + (DIMENSION + 1) + 1] - temp_B[Check_B_own_to_here + ii * 2 * (DIMENSION + 1) + 1];
            w_diff[2] = temp_B[Check_B_own_to_here + ii * 2 * (DIMENSION + 1) + (DIMENSION + 1) + 2] - temp_B[Check_B_own_to_here + ii * 2 * (DIMENSION + 1) + 2];

            // 辺1
            x_diff[3] = temp_B[Check_B_opponent_to_here + j * 2 * (DIMENSION + 1) + (DIMENSION + 1)] - temp_B[Check_B_opponent_to_here + j * 2 * (DIMENSION + 1)];
            y_diff[3] = temp_B[Check_B_opponent_to_here + j * 2 * (DIMENSION + 1) + (DIMENSION + 1) + 1] - temp_B[Check_B_opponent_to_here + j * 2 * (DIMENSION + 1) + 1];
            w_diff[3] = temp_B[Check_B_opponent_to_here + j * 2 * (DIMENSION + 1) + (DIMENSION + 1) + 2] - temp_B[Check_B_opponent_to_here + j * 2 * (DIMENSION + 1) + 2];

            if (sqrt(pow(x_diff[2], 2) + pow(y_diff[2], 2) + pow(w_diff[2], 2)) <= MERGE_DISTANCE || sqrt(pow(x_diff[3], 2) + pow(y_diff[3], 2) + pow(w_diff[3], 2)) <= MERGE_DISTANCE)
            {
                // 辺の両端が重なっている場合はスキップ
            }
            else if ((num_own == 22 && i == 2) || (num_own == 39 && i == 2))
            {
                // full モデルでパッチ番号 (22, 21), (39, 28) のコネクティビティを切断
                // (4, 3), (15, 8)では特殊なコネクティビティを作成する
            }
            else if (sqrt(pow(x_diff[0], 2) + pow(y_diff[0], 2) + pow(w_diff[0], 2)) <= MERGE_DISTANCE && sqrt(pow(x_diff[1], 2) + pow(y_diff[1], 2) + pow(w_diff[1], 2)) <= MERGE_DISTANCE)
            {
                // 辺が一致している場合 Edge_info を True
                temp_Edge_info[num_own * 32 + i * 8 + j] = 1;
                temp_Opponent_patch_num[num_own * 4 + i] = num_opponent;
                printf("own_patch:%d opp_patch:%d own_edge:%d opp_edge:%d\n", num_own, num_opponent, i, j);
                return;
            }
        }
    }
}

void Make_connectivity(int num, int *temp_CP_info, int *temp_Edge_info, int *temp_Opponent_patch_num, int *temp_Connectivity, int *temp_A, double *temp_CP, double *temp_CP_result)
{
    int i, j, k;
    int p, q;
    int temp_CP_n = 0;
    int Edge[4];
    int Edge_to_here = num * 32;

    printf("make connectivity on patch %d\n", num);

    int A_to_own = 0, A_to_opponent = 0;
    for (i = 0; i < num; i++)
    {
        A_to_own += 2 * (temp_CP_info[i * DIMENSION] + temp_CP_info[i * DIMENSION + 1]);
    }

    for (i = 0; i < 4; i++)
    {
        Edge[i] = 0;
    }

    // 重なっている辺の A 配列を作成
    for (i = 0; i < 4; i++)
    {
        // i == 0 ではなにもしない
        if (i == 1)
        {
            A_to_own += temp_CP_info[num * DIMENSION];
        }
        else if (i == 2)
        {
            A_to_own += temp_CP_info[num * DIMENSION + 1];
        }
        else if (i == 3)
        {
            A_to_own += temp_CP_info[num * DIMENSION];
        }

        for (j = 0; j < 8; j++)
        {
            if (temp_Edge_info[Edge_to_here + i * 8 + j] == 1)
            {
                A_to_opponent = 0;
                for (k = 0; k < temp_Opponent_patch_num[num * 4 + i]; k++)
                {
                    A_to_opponent += 2 * (temp_CP_info[k * DIMENSION] + temp_CP_info[k * DIMENSION + 1]);
                }

                printf("Patch num = %d\n", num);
                printf("Edge num = %d\n", j);
                printf("Opponent patch num = %d\n", temp_Opponent_patch_num[num * 4 + i]);

                Edge[i] = 1;

                p = j / 2;
                q = j % 2;
                printf("p = %d\n", p);
                printf("q = %d\n", q);

                if (p == 0)
                {
                    temp_CP_n = temp_CP_info[temp_Opponent_patch_num[num * 4 + i] * DIMENSION];
                }
                else if (p == 1)
                {
                    temp_CP_n = temp_CP_info[temp_Opponent_patch_num[num * 4 + i] * DIMENSION + 1];
                    A_to_opponent += temp_CP_info[temp_Opponent_patch_num[num * 4 + i] * DIMENSION];
                }
                else if (p == 2)
                {
                    temp_CP_n = temp_CP_info[temp_Opponent_patch_num[num * 4 + i] * DIMENSION];
                    A_to_opponent += temp_CP_info[temp_Opponent_patch_num[num * 4 + i] * DIMENSION] + temp_CP_info[temp_Opponent_patch_num[num * 4 + i] * DIMENSION + 1];
                }
                else if (p == 3)
                {
                    temp_CP_n = temp_CP_info[temp_Opponent_patch_num[num * 4 + i] * DIMENSION + 1];
                    A_to_opponent += 2 * temp_CP_info[temp_Opponent_patch_num[num * 4 + i] * DIMENSION] + temp_CP_info[temp_Opponent_patch_num[num * 4 + i] * DIMENSION + 1];
                }

                // printf("A to own = %d\n", A_to_own);
                // printf("temp_CP_n = %d\n", temp_CP_n);
                // printf("A to opponent = %d\n", A_to_opponent);
                // printf("temp_A[A_to_opponent] = %d\n", temp_A[A_to_opponent]);

                if (q == 0)
                {
                    for (k = 0; k < temp_CP_n; k++)
                    {
                        temp_A[A_to_own + k] = temp_A[A_to_opponent + k];
                        // printf("temp_A[A_to_own + k] = %d\n", temp_A[A_to_own + k]);
                    }
                    break;
                }
                else if (q == 1)
                {
                    for (k = 0; k < temp_CP_n; k++)
                    {
                        temp_A[A_to_own + k] = temp_A[A_to_opponent + (temp_CP_n - 1) - k];
                    }
                    break;
                }
            }
        }
    }

    printf("patch %d array A before make connectivity\n", num);
    A_to_own = 0;
    for (i = 0; i < num; i++)
    {
        A_to_own += 2 * (temp_CP_info[i * DIMENSION] + temp_CP_info[i * DIMENSION + 1]);
    }
    for (i = 0; i < 2 * (temp_CP_info[num * DIMENSION] + temp_CP_info[num * DIMENSION + 1]); i++)
    {
        printf("%d\t", temp_A[A_to_own + i]);
    }
    printf("\n");

    // コネクティビティを作成
    int xi, eta;

    A_to_own = 0;
    for (i = 0; i < num; i++)
    {
        A_to_own += 2 * (temp_CP_info[i * DIMENSION] + temp_CP_info[i * DIMENSION + 1]);
    }

    // パッチ番号 (4, 3), (15, 8)では特殊なコネクティビティを作成する
    if ((num == 4) || (num == 15))
    {
        if (num == 4)
        {
            for (eta = 0; eta < temp_CP_info[num * DIMENSION + 1]; eta++)
            {
                for (xi = 0; xi < temp_CP_info[num * DIMENSION]; xi++)
                {
                    if (xi == 0 && eta == 0)
                    {
                        temp_Connectivity[CP_to_here_counter + eta * temp_CP_info[num * DIMENSION] + xi] = temp_A[A_to_own + xi];
                    }
                    else
                    {
                        temp_Connectivity[CP_to_here_counter + eta * temp_CP_info[num * DIMENSION] + xi] = counter;
                        temp_CP_result[CP_result_to_here] = temp_CP[(CP_to_here_counter + eta * temp_CP_info[num * DIMENSION] + xi) * (DIMENSION + 1)];
                        temp_CP_result[CP_result_to_here + 1] = temp_CP[(CP_to_here_counter + eta * temp_CP_info[num * DIMENSION] + xi) * (DIMENSION + 1) + 1];
                        temp_CP_result[CP_result_to_here + 2] = temp_CP[(CP_to_here_counter + eta * temp_CP_info[num * DIMENSION] + xi) * (DIMENSION + 1) + 2];
                        counter++;
                        CP_result_to_here += (DIMENSION + 1);
                    }
                }
            }
        }
        else if (num == 15)
        {
            for (eta = 0; eta < temp_CP_info[num * DIMENSION + 1]; eta++)
            {
                for (xi = 0; xi < temp_CP_info[num * DIMENSION]; xi++)
                {
                    if (eta == 0 && Edge[0] == 1)
                    {
                        temp_Connectivity[CP_to_here_counter + eta * temp_CP_info[num * DIMENSION] + xi] = temp_A[A_to_own + xi];
                    }
                    else if (xi == 0 && eta == temp_CP_info[num * DIMENSION + 1] - 1)
                    {
                        temp_Connectivity[CP_to_here_counter + eta * temp_CP_info[num * DIMENSION] + xi] = temp_A[A_to_own + temp_CP_info[num * DIMENSION] + temp_CP_info[num * DIMENSION + 1] + xi];
                    }
                    else
                    {
                        temp_Connectivity[CP_to_here_counter + eta * temp_CP_info[num * DIMENSION] + xi] = counter;
                        temp_CP_result[CP_result_to_here] = temp_CP[(CP_to_here_counter + eta * temp_CP_info[num * DIMENSION] + xi) * (DIMENSION + 1)];
                        temp_CP_result[CP_result_to_here + 1] = temp_CP[(CP_to_here_counter + eta * temp_CP_info[num * DIMENSION] + xi) * (DIMENSION + 1) + 1];
                        temp_CP_result[CP_result_to_here + 2] = temp_CP[(CP_to_here_counter + eta * temp_CP_info[num * DIMENSION] + xi) * (DIMENSION + 1) + 2];
                        counter++;
                        CP_result_to_here += (DIMENSION + 1);
                    }
                }
            }
        }
    }
    else
    {
        for (eta = 0; eta < temp_CP_info[num * DIMENSION + 1]; eta++)
        {
            for (xi = 0; xi < temp_CP_info[num * DIMENSION]; xi++)
            {
                if (eta == 0 && Edge[0] == 1)
                {
                    temp_Connectivity[CP_to_here_counter + eta * temp_CP_info[num * DIMENSION] + xi] = temp_A[A_to_own + xi];
                }
                else if (eta == temp_CP_info[num * DIMENSION + 1] - 1 && Edge[2] == 1)
                {
                    temp_Connectivity[CP_to_here_counter + eta * temp_CP_info[num * DIMENSION] + xi] = temp_A[A_to_own + temp_CP_info[num * DIMENSION] + temp_CP_info[num * DIMENSION + 1] + xi];
                }
                else if (xi == 0 && Edge[3] == 1)
                {
                    temp_Connectivity[CP_to_here_counter + eta * temp_CP_info[num * DIMENSION] + xi] = temp_A[A_to_own + 2 * temp_CP_info[num * DIMENSION] + temp_CP_info[num * DIMENSION + 1] + eta];
                }
                else if (xi == temp_CP_info[num * DIMENSION] - 1 && Edge[1] == 1)
                {
                    temp_Connectivity[CP_to_here_counter + eta * temp_CP_info[num * DIMENSION] + xi] = temp_A[A_to_own + temp_CP_info[num * DIMENSION] + eta];
                }
                else
                {
                    temp_Connectivity[CP_to_here_counter + eta * temp_CP_info[num * DIMENSION] + xi] = counter;
                    temp_CP_result[CP_result_to_here] = temp_CP[(CP_to_here_counter + eta * temp_CP_info[num * DIMENSION] + xi) * (DIMENSION + 1)];
                    temp_CP_result[CP_result_to_here + 1] = temp_CP[(CP_to_here_counter + eta * temp_CP_info[num * DIMENSION] + xi) * (DIMENSION + 1) + 1];
                    temp_CP_result[CP_result_to_here + 2] = temp_CP[(CP_to_here_counter + eta * temp_CP_info[num * DIMENSION] + xi) * (DIMENSION + 1) + 2];
                    counter++;
                    CP_result_to_here += (DIMENSION + 1);
                }
            }
        }
    }

    // A 配列の作ってない分を作成
    for (eta = 0; eta < temp_CP_info[num * DIMENSION + 1]; eta++)
    {
        for (xi = 0; xi < temp_CP_info[num * DIMENSION]; xi++)
        {
            if (eta == 0 && Edge[0] == 0)
            {
                temp_A[A_to_own + xi] = temp_Connectivity[CP_to_here_counter + eta * temp_CP_info[num * DIMENSION] + xi];
            }
            if (eta == temp_CP_info[num * DIMENSION + 1] - 1 && Edge[2] == 0)
            {
                temp_A[A_to_own + temp_CP_info[num * DIMENSION] + temp_CP_info[num * DIMENSION + 1] + xi] = temp_Connectivity[CP_to_here_counter + eta * temp_CP_info[num * DIMENSION] + xi];
            }
            if (xi == 0 && Edge[3] == 0)
            {
                temp_A[A_to_own + 2 * temp_CP_info[num * DIMENSION] + temp_CP_info[num * DIMENSION + 1] + eta] = temp_Connectivity[CP_to_here_counter + eta * temp_CP_info[num * DIMENSION] + xi];
            }
            if (xi == temp_CP_info[num * DIMENSION] - 1 && Edge[1] == 0)
            {
                temp_A[A_to_own + temp_CP_info[num * DIMENSION] + eta] = temp_Connectivity[CP_to_here_counter + eta * temp_CP_info[num * DIMENSION] + xi];
            }
        }
    }
    CP_to_here_counter += temp_CP_info[num * DIMENSION] * temp_CP_info[num * DIMENSION + 1];

    printf("patch %d array A after make connectivity\n", num);
    A_to_own = 0;
    for (i = 0; i < num; i++)
    {
        A_to_own += 2 * (temp_CP_info[i * DIMENSION] + temp_CP_info[i * DIMENSION + 1]);
    }
    for (i = 0; i < 2 * (temp_CP_info[num * DIMENSION] + temp_CP_info[num * DIMENSION + 1]); i++)
    {
        printf("%d\t", temp_A[A_to_own + i]);
    }
    printf("\n");

    printf("\n");
}

void Output_inputdata(int *temp_Order, int *temp_KV_info, int *temp_CP_info, int *temp_Connectivity, double *temp_KV, double *temp_CP_result)
{
    int i, j, k;
    int Total_patch = 0;
    char str[256] = "input_loc.txt";

    fp = fopen(str, "w");

    // ヤング率
    fprintf(fp, "%d  ", (int)E_and_nu[0]);

    // ポアソン比
    fprintf(fp, "%le\n\n", E_and_nu[1]);

    // パッチ数
    if (mode[0] == 0)
    {
        Total_patch = 10;
    }
    else if (mode[0] == 1)
    {
        Total_patch = 40 + 24;
    }
    fprintf(fp, "%d\n\n", Total_patch);

    // コントロールポイント数
    fprintf(fp, "%d\n\n", (CP_result_to_here + 1) / 3);
    int temp_num = (CP_result_to_here + 1) / 3, temp_counter = 0;
    while (temp_num != 0)
    {
        temp_num = temp_num / 10;
        temp_counter++;
    }
    temp_num = -(temp_counter + 2);

    // 各パッチ内での各方向の次数
    for (i = 0; i < Total_patch; i++)
    {
        for (j = 0; j < DIMENSION; j++)
        {
            if (j == DIMENSION - 1)
            {
                fprintf(fp, "%d", temp_Order[i * DIMENSION + j]);
            }
            else
            {
                fprintf(fp, "%*d", -6, temp_Order[i * DIMENSION + j]);
            }
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // 各パッチ内での各方向のノットベクトルの数
    for (i = 0; i < Total_patch; i++)
    {
        for (j = 0; j < DIMENSION; j++)
        {
            if (j == DIMENSION - 1)
            {
                fprintf(fp, "%d", temp_KV_info[i * DIMENSION + j]);
            }
            else
            {
                fprintf(fp, "%*d", -6, temp_KV_info[i * DIMENSION + j]);
            }
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // 各パッチ内での各方向のコントロールポイントの数
    for (i = 0; i < Total_patch; i++)
    {
        for (j = 0; j < DIMENSION; j++)
        {
            if (j == DIMENSION - 1)
            {
                fprintf(fp, "%d", temp_CP_info[i * DIMENSION + j]);
            }
            else
            {
                fprintf(fp, "%*d", -6, temp_CP_info[i * DIMENSION + j]);
            }
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // パッチコネクティビティ
    int CP_to_here = 0;
    for (i = 0; i < Total_patch; i++)
    {
        for (j = 0; j < temp_CP_info[i * DIMENSION] * temp_CP_info[i * DIMENSION + 1]; j++)
        {
            if (j == temp_CP_info[i * DIMENSION] * temp_CP_info[i * DIMENSION + 1] - 1)
            {
                fprintf(fp, "%d", temp_Connectivity[CP_to_here + j]);
            }
            else
            {
                fprintf(fp, "%*d", temp_num, temp_Connectivity[CP_to_here + j]);
            }
        }
        CP_to_here += temp_CP_info[i * DIMENSION] * temp_CP_info[i * DIMENSION + 1];
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // 変位拘束するコントロールポイントの数
    fprintf(fp, "????   ");

    // 荷重条件を与えるコントロールポイントの数
    fprintf(fp, "0  ");

    // 分布荷重の数
    fprintf(fp, "????\n\n");

    // 各パッチでの各方向のノットベクトル
    int KV_to_here = 0;
    for (i = 0; i < Total_patch; i++)
    {
        for (j = 0; j < DIMENSION; j++)
        {
            for (k = 0; k < temp_KV_info[i * DIMENSION + j]; k++)
            {
                if (k == temp_KV_info[i * DIMENSION + j] - 1)
                {
                    fprintf(fp, "%.16e", temp_KV[KV_to_here + k]);
                }
                else
                {
                    fprintf(fp, "%.16e  ", temp_KV[KV_to_here + k]);
                }
            }
            KV_to_here += temp_KV_info[i * DIMENSION + j];
            fprintf(fp, "\n");
        }
    }
    fprintf(fp, "\n");

    // コントロールポイント
    for (i = 0; i < (CP_result_to_here + 1) / 3; i++)
    {
        fprintf(fp, "%*d", temp_num, i);
        fprintf(fp, "% .16e  ", temp_CP_result[i * 3]);
        fprintf(fp, "% .16e  ", temp_CP_result[i * 3 + 1]);
        fprintf(fp, "% .16e\n", temp_CP_result[i * 3 + 2]);
    }
    fprintf(fp, "\n");

    fclose(fp);
}

void Output_J(int *temp_CP_info, int *temp_A)
{
    int i, j;
    char str[256] = "J_int.txt";

    double PI = 3.14159265358979323846264338327950288;
    double theta = affine[2] * PI / 180.0;
    double rot[4] = {cos(theta), -sin(theta), sin(theta), cos(theta)};

    fp = fopen(str, "w");

    // き裂先端座標
    double x = rot[0] * crack_tip[0] + rot[1] * crack_tip[1];
    double y = rot[2] * crack_tip[0] + rot[3] * crack_tip[1];

    fprintf(fp, "1 %.17e %.17e %.17e\n", x, y, 1.0);

    // q function
    int q_func_edge_n = 0;
    int q_info[8][3]; // パッチ番号, xi or eta, 始点0 or 終点1
    if (mode[0] == 0)
    {
        q_func_edge_n = 4;
        q_info[0][0] = 0;
        q_info[0][1] = 0;
        q_info[0][2] = 0;
        q_info[1][0] = 1;
        q_info[1][1] = 0;
        q_info[1][2] = 0;
        q_info[2][0] = 2;
        q_info[2][1] = 0;
        q_info[2][2] = 0;
        q_info[3][0] = 3;
        q_info[3][1] = 0;
        q_info[3][2] = 0;
    }
    else if (mode[0] == 1)
    {
        q_func_edge_n = 8;
        q_info[0][0] = 0;
        q_info[0][1] = 0;
        q_info[0][2] = 0;
        q_info[1][0] = 1;
        q_info[1][1] = 0;
        q_info[1][2] = 0;
        q_info[2][0] = 2;
        q_info[2][1] = 0;
        q_info[2][2] = 0;
        q_info[3][0] = 3;
        q_info[3][1] = 0;
        q_info[3][2] = 0;
        q_info[4][0] = 4;
        q_info[4][1] = 0;
        q_info[4][2] = 0;
        q_info[5][0] = 5;
        q_info[5][1] = 0;
        q_info[5][2] = 0;
        q_info[6][0] = 6;
        q_info[6][1] = 0;
        q_info[6][2] = 0;
        q_info[7][0] = 7;
        q_info[7][1] = 0;
        q_info[7][2] = 0;
    }

    int before_l = q_func_edge_n * 3;
    int *temp_array = (int *)malloc(sizeof(int) * before_l);
    int *new_array = (int *)malloc(sizeof(int) * before_l);
    if (temp_array == NULL || new_array == NULL)
    {
        printf("Memory cannot be allocated\n");
        exit(1);
    }

    int temp = 0;
    for (i = 0; i < q_func_edge_n; i++)
    {
        int A_to_here = 0;
        for (j = 0; j < q_info[i][0]; j++)
        {
            A_to_here += 2 * (temp_CP_info[j * DIMENSION] + temp_CP_info[j * DIMENSION + 1]);
        }

        // if (q_info[i][1] == 1 && q_info[i][2] == 0) は何もしない
        if (q_info[i][1] == 0 && q_info[i][2] == 1)
        {
            A_to_here += temp_CP_info[q_info[i][0] * DIMENSION];
        }
        else if (q_info[i][1] == 1 && q_info[i][2] == 1)
        {
            A_to_here += temp_CP_info[q_info[i][0] * DIMENSION] + temp_CP_info[q_info[i][0] * DIMENSION + 1];
        }
        else if (q_info[i][1] == 0 && q_info[i][2] == 0)
        {
            A_to_here += 2 * temp_CP_info[q_info[i][0] * DIMENSION] + temp_CP_info[q_info[i][0] * DIMENSION + 1];
        }

        if (q_info[i][1] == 0)
        {
            for (j = 0; j < temp_CP_info[q_info[i][0] * DIMENSION + 1]; j++)
            {
                temp_array[temp] = temp_A[A_to_here + j];
                temp++;

                // printf("patch num %d\n", q_info[i][0]);
                // printf("xi or eta %d\n", q_info[i][1]);
                // printf("start or end %d\n", q_info[i][2]);
                // printf("temp = %d\n", temp);
                // printf("A_to_here + l = %d\n", A_to_here + l);
                // printf("temp_A[A_to_here + l] = %d\n", temp_A[A_to_here + l]);
            }
        }
        else if (q_info[i][1] == 1)
        {
            for (j = 0; j < temp_CP_info[q_info[i][0] * DIMENSION]; j++)
            {
                temp_array[temp] = temp_A[A_to_here + j];
                temp++;

                // printf("patch num %d\n", q_info[i][0]);
                // printf("xi or eta %d\n", q_info[i][1]);
                // printf("start or end %d\n", q_info[i][2]);
                // printf("temp = %d\n", temp);
                // printf("A_to_here + l = %d\n", A_to_here + l);
                // printf("temp_A[A_to_here + l] = %d\n", temp_A[A_to_here + l]);
            }
        }
    }

    heapSort(temp_array, before_l);

    int after_l = 0;
    for (i = 0; i < before_l - 1; i++)
    {
        if (temp_array[i] < temp_array[i + 1])
        {
            new_array[after_l] = temp_array[i];
            after_l++;
        }
    }
    new_array[after_l] = temp_array[before_l - 1];
    after_l++;

    fprintf(fp, "%d\n", after_l);
    for (i = 0; i < after_l; i++)
    {
        fprintf(fp, "%d 1.0 0.0\n", new_array[i]);
    }

    free(temp_array), free(new_array);

    fclose(fp);
}

void Output_SVG(double *temp_B, double *temp_CP_result)
{
    int i;

    char color_vec[10][10] = {"#696969", "#a9a9a9", "#00bfff", "#00fa9a", "#ffff00", "#ff8c00", "#cd5c5c", "#ff7f50", "#ee82ee", "#8a2be2"};
    //  0   darkgray
    //  1   deepskyblue
    //  2   mediumspringgreen
    //  3   yellow
    //  4   darkorange
    //  5   indianred
    //  6   coral
    //  7   violet
    //  8   blueviolet
    //  https://www.colordic.org/

    char num_color[10] = "#dc143c";

    double x_min = 0, x_max = 1;
    double y_min = 0, y_max = 1;

    double position_x, position_y;

    for (i = 0; i < (CP_result_to_here + 1) / 3; i++)
    {
        if (i == 0)
        {
            x_min = temp_CP_result[i * 3];
            x_max = temp_CP_result[i * 3];
            y_min = temp_CP_result[i * 3 + 1];
            y_max = temp_CP_result[i * 3 + 1];
        }
        else
        {
            if (x_min > temp_CP_result[i * 3])
            {
                x_min = temp_CP_result[i * 3];
            }
            else if (x_max < temp_CP_result[i * 3])
            {
                x_max = temp_CP_result[i * 3];
            }

            if (y_min > temp_CP_result[i * 3 + 1])
            {
                y_min = temp_CP_result[i * 3 + 1];
            }
            else if (y_max < temp_CP_result[i * 3 + 1])
            {
                y_max = temp_CP_result[i * 3 + 1];
            }
        }
    }

    printf("x y : %le %le\n", x_max, y_max);
    printf("x y : %le %le\n", x_min, y_min);

    double space = 3.0;
    double scale = 1000.0 / (x_max - x_min + 2.0 * space);

    double width = 1.5 * (x_max - x_min + 2.0 * space) * scale;
    double height = (y_max - x_min + 2.0 * space) * scale;

    printf("width = %le\n", width);
    printf("height = %le\n", height);

    char str[256] = "input.svg";

    fp = fopen(str, "w");

    fprintf(fp, "<?xml version='1.0'?>\n");
    // fprintf(fp, "<svg width='%lept' height='%lept' viewBox='0 0 %le %le' style = 'background: #eee' xmlns='http://www.w3.org/2000/svg' xmlns:xlink='http://www.w3.org/1999/xlink'>\n", width, height, width, height);
    fprintf(fp, "<svg width='%le' height='%le' version='1.1' style='background: #eee' xmlns='http://www.w3.org/2000/svg' xmlns:xlink='http://www.w3.org/1999/xlink'>\n", width, height);

    // パッチ境界を描画
    int temp_color_num = 0;
    int B_to_here = 0;
    int Total_patch = 0;

    if (mode[0] == 0)
    {
        Total_patch = 10;
    }
    else if (mode[0] == 1)
    {
        Total_patch = 40 + 24;
    }

    double PI = 3.14159265358979323846264338327950288;
    double theta = affine[2] * PI / 180.0;
    double rot[4] = {cos(theta), -sin(theta), sin(theta), cos(theta)};
    double temp_x, temp_y;

    for (i = 0; i < Total_patch; i++)
    {
        temp_x = rot[0] * temp_B[B_to_here] + rot[1] * temp_B[B_to_here + 1];
        temp_y = rot[2] * temp_B[B_to_here] + rot[3] * temp_B[B_to_here + 1];
        position_x = ((temp_x + space) * scale) + width * (1.0 / 4.0);
        position_y = height / 2.0 - ((temp_y + space) * scale);
        fprintf(fp, "<path d='M %le %le ", position_x, position_y);
        B_to_here += 4 * (DIMENSION + 1);

        temp_x = rot[0] * temp_B[B_to_here] + rot[1] * temp_B[B_to_here + 1];
        temp_y = rot[2] * temp_B[B_to_here] + rot[3] * temp_B[B_to_here + 1];
        position_x = ((temp_x + space) * scale) + width * (1.0 / 4.0);
        position_y = height / 2.0 - ((temp_y + space) * scale);
        fprintf(fp, "L %le %le ", position_x, position_y);
        B_to_here += 2 * (DIMENSION + 1);

        temp_x = rot[0] * temp_B[B_to_here] + rot[1] * temp_B[B_to_here + 1];
        temp_y = rot[2] * temp_B[B_to_here] + rot[3] * temp_B[B_to_here + 1];
        position_x = ((temp_x + space) * scale) + width * (1.0 / 4.0);
        position_y = height / 2.0 - ((temp_y + space) * scale);
        fprintf(fp, "L %le %le ", position_x, position_y);
        B_to_here += 2 * (DIMENSION + 1);

        temp_x = rot[0] * temp_B[B_to_here] + rot[1] * temp_B[B_to_here + 1];
        temp_y = rot[2] * temp_B[B_to_here] + rot[3] * temp_B[B_to_here + 1];
        position_x = ((temp_x + space) * scale) + width * (1.0 / 4.0);
        position_y = height / 2.0 - ((temp_y + space) * scale);
        fprintf(fp, "L %le %le ", position_x, position_y);
        B_to_here += 4 * (DIMENSION + 1);

        fprintf(fp, "Z' fill='%s'/>\n", color_vec[temp_color_num % 10]);
        B_to_here += 4 * (DIMENSION + 1);

        temp_color_num++;
    }

    // 点と番号を描画
    for (i = 0; i < (CP_result_to_here + 1) / 3; i++)
    {
        position_x = (temp_CP_result[i * 3] + space) * scale + width * (1.0 / 4.0);
        position_y = height / 2.0 - ((temp_CP_result[i * 3 + 1] + space) * scale);
        fprintf(fp, "<circle cx='%le' cy='%le' r='2' fill='%s'/>\n", position_x, position_y, num_color);
        fprintf(fp, "<text x='%le' y='%le' font-family='Verdana' font-size='6' fill='%s' font-weight='700'>\n", position_x + 2, position_y, num_color);
        fprintf(fp, "%d\n", i);
        fprintf(fp, "</text>\n");
    }

    fprintf(fp, "</svg>");
    fclose(fp);
}

/* データを入れ替える関数 */
void swap(int *a, int *b)
{
    int temp;
    temp = *b;
    *b = *a;
    *a = temp;
}

/* 左の子ノードの位置を取得 */
int getLeft(int parent)
{
    return parent * 2 + 1;
}

/* 右の子ノードの位置を取得 */
int getRight(int parent)
{
    return parent * 2 + 2;
}

/* 親ノードの位置を取得 */
int getParent(int child)
{
    return (child - 1) / 2;
}

/* a[size]を二分ヒープに追加し、二分ヒープを再構成する */
void addHeap(int *a, int size)
{
    int add;    /* 追加ノードの位置 */
    int parent; /* 追加ノードの親の位置 */

    /* まだ二分ヒープに追加していないデータの先頭を二分ヒープに追加 */
    add = size;
    if (add == 0)
    {
        /* 追加したノードが根ノードなら二分ヒープへの追加完了 */
        return;
    }

    /* 二分ヒープを満たすまで、追加したノードを根の方向に移動する */
    while (1)
    {
        /* 親ノードの位置を取得 */
        parent = getParent(add);

        if (a[parent] < a[add])
        {
            /* 親と子で大小関係が逆ならデータを交換する */
            swap(&a[parent], &a[add]);

            /* 追加ノードは親ノードの位置に移動する */
            add = parent;
            if (add == 0)
            {
                /* 追加ノードが根ノードまで移動したら二分ヒープへの追加完了 */
                break;
            }
        }
        else
        {
            /* 大小関係が満たされているなら二分ヒープへの追加完了 */
            break;
        }
    }
}

/* 根ノードを二分ヒープから取り出し、二分ヒープを再構成する */
void removeHeap(int *a, int size)
{
    int left;   /* 左の子ノードの位置 */
    int right;  /* 右の子ノードの位置 */
    int large;  /* データが大きい方の子ノードの位置 */
    int parent; /* 親ノードの位置 */

    /* 根ノードをヒープ外に追い出す */
    /* 一時的に木の末端のノードを根ノードに設定する */
    swap(&a[0], &a[size - 1]);

    /* 二分ヒープのサイズを1減らす
        これにより元々の根ノードが「ソート済みのデータ」の先頭に移動することになる */
    size--;

    /* 根ノードから子ノードとの大小関係を確認していく */
    parent = 0;

    /* 二分ヒープを満たすまで、根ノードを葉の方向に移動する */
    while (1)
    {
        /* 子ノードの位置を取得 */
        left = getLeft(parent);
        right = getRight(parent);

        /* 子ノードの大きい値を持つ方の位置を取得 */
        if (left < size && right < size)
        {
            /* 左右両方の子ノードが存在する場合は比較して確認 */
            if (a[left] < a[right])
            {
                large = right;
            }
            else
            {
                large = left;
            }
        }
        else if (left < size)
        {
            /* 左の子ノードしか存在しない場合は左の子ノードを大きい値を持つとみなす */
            large = left;
        }
        else
        {
            /* 両ノードがヒープ内に存在しない場合は終了 */
            /* (右の子ノードしか存在しない場合はあり得ない) */
            break;
        }

        if (a[large] <= a[parent])
        {
            /* すでに親子の大小関係が満たされているので交換不要 */
            break;
        }

        /* 親と子で大小関係が逆ならデータを交換する */
        swap(&a[large], &a[parent]);

        /* 根ノードはデータを交換した子ノードの位置に移動する */
        parent = large;
    }
}

/* 二分ヒープを作成する関数 */
void makeHeap(int *a, int num)
{
    int size; /* 二分ヒープに追加済みのデータの個数 */

    /* 二分ヒープのデータ個数を0にする */
    size = 0;

    /* sizeがソートするデータの個数になるまで二分ヒープにデータ追加 */
    while (size < num)
    {
        /* a[size]を二分ヒープに追加 */
        addHeap(a, size);

        /* 二分ヒープのデータ数が増えたのでsizeも1増やす */
        size++;
    }
}

/* ヒープソートを行う関数 */
void heapSort(int *a, int num)
{
    int size; /* 二分ヒープのノード個数 */

    /* サイズnumの二分ヒープを作成 */
    makeHeap(a, num);

    /* 二分ヒープの根ノードを１つずつ取り出す */
    for (size = num; size > 0; size--)
    {
        /* サイズsizeの二分ヒープからデータを１つ取り出す */
        removeHeap(a, size);
    }
}

void Dedupe(int *a, int *num, int *a_new, int *num_new, int n)
{
    int i;
    int temp_to_here = 0;
    int temp_length = 0;

    for (i = 0; i < n; i++)
    {
        temp_to_here += num[n];
    }

    for (i = 0; i < num[n] - 1; i++)
    {
        if (a[i] < a[i + 1])
        {
            a_new[temp_to_here] = a[i];
            temp_to_here++;
            temp_length++;
        }
    }
    a_new[temp_to_here] = a[num[n] - 1];
    temp_length++;

    num_new[n] = temp_length;
}

double rot_rev(double x, double y, double bool)
{
    static double PI = 3.14159265358979323846264338327950288;
    double result = 0.0;
    double theta = - affine[2] * PI / 180.0;
    double rot[4] = {cos(theta), -sin(theta), sin(theta), cos(theta)};

    if (bool == 0)      // return x
        result = rot[0] * x + rot[1] * y;
    else if (bool == 1) // return y
        result = rot[2] * x + rot[3] * y;

    return result;
}

double rot(double x, double y, double bool)
{
    static double PI = 3.14159265358979323846264338327950288;
    double result = 0.0;
    double theta = affine[2] * PI / 180.0;
    double rot[4] = {cos(theta), -sin(theta), sin(theta), cos(theta)};

    if (bool == 0)      // return x
        result = rot[0] * x + rot[1] * y;
    else if (bool == 1) // return y
        result = rot[2] * x + rot[3] * y;

    return result;
}