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
// #define MERGE_DISTANCE 1.0e-13      // コントロールポイントが同じ点と判定する距離
// #define MAX_DISP_CONSTRAINT 10      // 変位指定する変位量の最大個数
// #define MAX_DISP_CONSTRAINT_EDGE 10 // 変位指定する辺の最大個数
// #define MAX_DISTRIBUTED_LOAD 5      // 分布荷重の最大個数

FILE *fp;

static double E_and_nu[2];
static int mode[3];                 // 0: analysis model mode, 1: making mode, 2: sigular patch mode
static double crack_tip[DIMENSION];
static double singular_width;
static double length_each_side[3];  // (outer side, inner side, upper side)
static int cp_after_each_side[3];      // (outer side, inner side, upper side)
static int glo_order[DIMENSION];
static double glo_size[DIMENSION];
static int glo_cp_after[DIMENSION];


void Get_input(char *filename);
void Make_quarter_model();
void Make_full_model();
void Make_glo_patch();
void Make_quarter_glo_patch();
void Make_full_glo_patch();


int main(int argc, char **argv)
{
    // int i;

    printf("\nargc = %d\n", argc);
    printf("argv[1] = %s\n\n", argv[1]);

    // ファイル読み込み
    Get_input(argv[1]);

    // モデル分岐(1/4 or full)
    if (mode[0] == 0)
    {
        Make_quarter_model();
    }
    else if (mode[0] == 1)
    {
        Make_full_model();
    }

    // global patch 作成分岐
    if (mode[1] == 1 && mode[0] == 0)
    {
        Make_quarter_glo_patch();
    }
    else if (mode[1] == 1 && mode[0] == 1)
    {
        Make_full_glo_patch();
    }




    return 0;
}


void Get_input(char *filename)
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

    // mode[2]
    fscanf(fp, "%d", &temp_i);
    mode[2] = temp_i;
    printf("mode[2] = %d\n", mode[2]);

    fgets(s, 256, fp);

    // crack_tip
    for (i = 0; i < DIMENSION; i++)
    {
        fscanf(fp, "%lf", &temp_d);
        crack_tip[i] = temp_d;
        printf("crack_tip[%d] = %le\n", i, crack_tip[i]);
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

    fgets(s, 256, fp);

    // global patch
    if (mode[1] == 1)
    {
        // order
        for (i = 0; i < DIMENSION; i++)
        {
            fscanf(fp, "%d", &temp_i);
            glo_order[i] = temp_i;
            printf("glo_order[%d] = %d\n", i, glo_order[i]);
        }

        fgets(s, 256, fp);

        // glo_size
        for (i = 0; i < DIMENSION; i++)
        {
            fscanf(fp, "%lf", &temp_d);
            glo_size[i] = temp_d;
            printf("glo_size[%d] = %le\n", i, glo_size[i]);
        }

        fgets(s, 256, fp);

        // glo_cp_after
        for (i = 0; i < DIMENSION; i++)
        {
            fscanf(fp, "%d", &temp_i);
            glo_cp_after[i] = temp_i;
            printf("glo_cp_after[%d] = %d\n", i, glo_cp_after[i]);
        }
        
    }
    fclose(fp);
}


void Make_quarter_model()
{
    double Total_patch_loc = 10;

    // 動的メモリ確保
    int *Order = (int *)malloc(sizeof(int) * Total_patch_loc * DIMENSION);      // int Order[パッチ番号][DIMENSION]
    int *KV_info = (int *)malloc(sizeof(int) * Total_patch_loc * DIMENSION);    // int KV_info[パッチ番号][DIMENSION]
    int *CP_info = (int *)malloc(sizeof(int) * Total_patch_loc * DIMENSION);    // int CP_info[パッチ番号][DIMENSION]

    if (Order == NULL || KV_info == NULL || CP_info == NULL)
    {
        printf("Memory cannot be allocated\n");
        exit(1);
    }

    // 特異パッチの種類
    if (mode[2] == 1)
    {

    }
    else if (mode[2] == 2)
    {

    }
}


void Make_full_model()
{
    double Total_patch_loc = 40;

    // 特異パッチの種類
    if (mode[2] == 1)
    {

    }
    else if (mode[2] == 2)
    {

    }
}


void Make_quarter_glo_patch()
{

}


void Make_full_glo_patch()
{

}


void 

