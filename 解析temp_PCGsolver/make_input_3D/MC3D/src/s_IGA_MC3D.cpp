#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "header_MC3D.h"

int main(int argc, char **argv)
{
    int i, j, k;
    int Total_input = argc - 1;
    information info, *info_ptr;
    info_ptr = &info;

    // ファイル読み込み(1回目)
    Get_inputdata_boundary_0(argv[1], &info); // boundaryのインプットデータ処理

    // 動的メモリ確保
    int *disp_constraint_n = (double *)malloc(sizeof(double) * info.DIMENSION);                                         // disp_constraint_n[DIMENSION]
    int *disp_constraint_edge_n = (double *)malloc(sizeof(double) * info.DIMENSION * info.MAX_DISP_CONSTRAINT);     // disp_constraint_edge_n[DIMENSION][MAX_DISP_CONSTRAINT]
    double *disp_constraint_amount = (double *)malloc(sizeof(double) * info.DIMENSION * info.MAX_DISP_CONSTRAINT);  // disp_constraint_amount[DIMENSION][MAX_DISP_CONSTRAINT]
    double *disp_constraint = (double *)malloc(sizeof(double) * info.DIMENSION * info.MAX_DISP_CONSTRAINT * info.MAX_DISP_CONSTRAINT_EDGE * 3);     // disp_constraint[DIMENSION][MAX_DISP_CONSTRAINT][MAX_DISP_CONSTRAINT_EDGE][3]
    double *distributed_load_info = (double *)malloc(sizeof(double) * info.distributed_load_n * 9);                     // distributed_load_info[MAX_DISTRIBUTED_LOAD][9]
    info_ptr->disp_constraint_n = disp_constraint_n;
    info_ptr->disp_constraint_edge_n = disp_constraint_edge_n;
    info_ptr->disp_constraint_amount = disp_constraint_amount;
    info_ptr->disp_constraint = disp_constraint;
    info_ptr->distributed_load_info = distributed_load_info;

    // ファイル読み込み(2回目)
    Get_inputdata_boundary_1(argv[1], &info); // boundaryのインプットデータ処理

    // 動的メモリ確保
    int *Order = (int *)malloc(sizeof(int) * info.Total_patch * info.DIMENSION);   // int Order[パッチ番号][DIMENSION]
    int *KV_info = (int *)malloc(sizeof(int) * info.Total_patch * info.DIMENSION); // int KV_info[パッチ番号][DIMENSION]
    int *CP_info = (int *)malloc(sizeof(int) * info.Total_patch * info.DIMENSION); // int CP_info[パッチ番号][DIMENSION]
    info_ptr->Order = Order;
    info_ptr->KV_info = KV_info;
    info_ptr->CP_info = CP_info;

    if (Order == NULL || KV_info == NULL || CP_info == NULL)
    {
        printf("Memory cannot be allocated\n");
        exit(1);
    }

    // patchのインプットデータ処理(1回目)
    for (i = 1; i < Total_input; i++)
    {
        Get_inputdata_patch_0(argv[i + 1], &info);
    }

    int temp1 = 0; // temp1 : 全パッチ含めた総コントロールポイント数
    int temp2 = 0; // temp2 : 4辺のコントロールポイントの和を全パッチ分足した値 * 2
    int temp3 = 0; // temp3 : 全パッチ含めた総ノットベクトル数

    for (i = 0; i < info.Total_patch; i++)
    {
        int Total_CP = 1;
        for (j = 0; j < info.DIMENSION; j++)
        {
            Total_CP *= info.CP_info[i * info.DIMENSION + j];
        }
        temp1 += Total_CP;
        if (info.DIMENSION == 2)
        {
            // 4辺のCP数 * パッチ数
            temp2 += 2 * (info.CP_info[i * info.DIMENSION] + info.CP_info[i * info.DIMENSION + 1]);
        }
        else if (info.DIMENSION == 3)
        {
            // 6面のCP数 * パッチ数
            temp2 += 2 * (info.CP_info[i * info.DIMENSION] * info.CP_info[i * info.DIMENSION + 1] + info.CP_info[i * info.DIMENSION + 1] * info.CP_info[i * info.DIMENSION + 2] + info.CP_info[i * info.DIMENSION] * info.CP_info[i * info.DIMENSION + 2]);
        }
        for (j = 0; j < info.DIMENSION; j++)
        {
            temp3 += KV_info[i * info.DIMENSION + j];
        }
    }
    printf("Total Control Point = %d\n", temp1);

    // 動的メモリ確保
    double *CP = (double *)malloc(sizeof(double) * temp1 * (info.DIMENSION + 1));           // double CP[CP番号][DIMENSION + 1]
    double *CP_result = (double *)malloc(sizeof(double) * temp1 * (info.DIMENSION + 1));    // double CP_result[通しのコントロールポイント番号(連番)][DIMENSION + 1]
    int *A = (int *)malloc(sizeof(int) * temp2);                                            // int    A[パッチ番号][面番号(0~6) or 辺番号(0~3)][辺内のコネクティビティ]
    double *B;
    if (info.DIMENSION == 2)
    {
        B = (double *)malloc(sizeof(double) * info.Total_patch * 4 * 2 * 2 * (info.DIMENSION + 1)); // double B[パッチ番号][辺番号(0~4)][正負方向 2][各辺の端の2頂点][座標xyw  -> 3]
    }
    else if (info.DIMENSION == 3)
    {
        B = (double *)malloc(sizeof(double) * info.Total_patch * 6 * 2 * 4 * (info.DIMENSION + 1)); // double B[パッチ番号][面番号(0~6)][正負方向 2][各面の端の4頂点][座標xyzw -> 4]
    }
    int *Connectivity = (int *)malloc(sizeof(int) * temp1);                                 // int    Connectivity[パッチ番号][パッチ内CP番号]
    double *KV = (double *)malloc(sizeof(double) * temp3);                                  // double KV[パッチ番号][DIMENSION][ノットベクトル番号]
    info_ptr->CP = CP;
    info_ptr->CP_result = CP_result;
    info_ptr->A = A;
    info_ptr->B = B;
    info_ptr->Connectivity = Connectivity;
    info_ptr->KV = KV;

    if (CP == NULL || CP_result == NULL || A == NULL || B == NULL || Connectivity == NULL || KV == NULL)
    {
        printf("Memory cannot be allocated\n");
        exit(1);
    }

    // patchのインプットデータ処理(2回目)
    counter = 0;
    KV_to_here = 0, CP_to_here = 0, B_to_here = 0;
    for (i = 1; i < Total_input; i++)
    {
        Get_inputdata_patch_1(argv[i + 1], &info, counter);
        counter++;
    }

    printf("Done get input\n");

    // 動的メモリ確保
    int *Face_Edge_info;
    int *Opponent_patch_num;
    if (info.DIMENSION == 2)
    {
        Face_Edge_info = (int *)malloc(sizeof(int) * info.Total_patch * 32);         // int Face_Edge_info[パッチ番号][own 辺番号(正固定0~3)][opp 辺番号(0~7)]
        Opponent_patch_num = (int *)malloc(sizeof(int) * info.Total_patch * 4);      // int Opponent_patch_num[パッチ番号][own 辺番号(正固定0~3]
    }
    else if (info.DIMENSION == 3)
    {
        Face_Edge_info = (int *)malloc(sizeof(int) * info.Total_patch * 32);         // int Face_Edge_info[パッチ番号][own 面番号(0~5)][opp 面番号(0~5)]
        Opponent_patch_num = (int *)malloc(sizeof(int) * info.Total_patch * 6);      // int Opponent_patch_num[パッチ番号][own 辺番号(正固定0~5]
    }
    info_ptr->Face_Edge_info = Face_Edge_info;
    info_ptr->Opponent_patch_num = Opponent_patch_num;

    if (Face_Edge_info == NULL || Opponent_patch_num == NULL)
    {
        printf("Memory cannot be allocated\n");
        exit(1);
    }

    for (i = 0; i < info.Total_patch * 32; i++)
    {
        Face_Edge_info[i] = 0;
    }

    // パッチコネクティビティの作成
    printf("state: patch connectivity\n");
    counter = 0;
    CP_to_here = 0, CP_result_to_here = 0, B_to_here = 0;
    for (i = 0; i < info.Total_patch; i++)
    {
        for (j = 0; j < i; j++)
        {
            Check_B(i, j, B, Face_Edge_info, Opponent_patch_num);
        }
        Make_connectivity(i, CP_info, Face_Edge_info, Opponent_patch_num, Connectivity, A, CP, CP_result);
    }

    // 動的メモリ確保
    int temp4 = 0, temp5 = 0;
    for (i = 0; i < info->DIMENSION; i++)
    {
        for (j = 0; j < disp_constraint_n[i]; j++)
        {
            for (k = 0; k < disp_constraint_edge_n[i][j]; k++)
            {
                if (disp_constraint[i][j][k][1] == 0)
                {
                    temp4 += CP_info[disp_constraint[i][j][k][0] * DIMENSION + 1];
                }
                else if (disp_constraint[i][j][k][1] == 1)
                {
                    temp4 += CP_info[disp_constraint[i][j][k][0] * DIMENSION];
                }
            }
        }
        temp5 += disp_constraint_n[i];
    }

    // printf("temp5 = %d\n", temp5);
    // printf("temp4 = %d\n", temp4);

    int *length_before = (int *)malloc(sizeof(int) * temp5);   // 各変位量でのマージ前の長さ
    int *length_after = (int *)malloc(sizeof(int) * temp5);    // 各変位量でのマージ後の長さ
    int *Boundary = (int *)malloc(sizeof(int) * temp4);        // 境界条件のコネクティビティ
    int *Boundary_result = (int *)malloc(sizeof(int) * temp4); // ソート・マージ後境界条件のコネクティビティ

    if (length_before == NULL || length_after == NULL || Boundary == NULL || Boundary_result == NULL)
    {
        printf("Memory cannot be allocated\n");
        exit(1);
    }

    // 強制変位・変位固定の境界条件を作成
    Sort(temp5, CP_info, A, Boundary, Boundary_result, length_before, length_after);

    printf("state: output\n");
    Output_inputdata(Order, KV_info, CP_info, Connectivity, KV, CP_result, Boundary_result, length_before, length_after, temp5);

    // 図の出力
    if (info.DIMENSION == 2)
    {
        Output_SVG(B, CP_result);   // SVG出力
    }

    // メモリ解放
    free(Order), free(KV_info), free(CP_info);
    free(CP), free(CP_result), free(A), free(B), free(Connectivity), free(KV);
    free(Face_Edge_info), free(Opponent_patch_num);
    free(length_before), free(length_after), free(Boundary), free(Boundary_result);

    return 0;
}


// get input data
void Get_inputdata_boundary_0(char *filename, information *info)
{
    int i, j, k;
    char s[256];

    int temp_i;
    double temp_d;

    fp = fopen(filename, "r");

    // DIMENSION
    fscanf(fp, "%d", &temp_i);
    info->DIMENSION = temp_i;
    printf("DIMENSION = %d\n", temp_i);

    fgets(s, 256, fp);

    // パッチ数
    fscanf(fp, "%d", &temp_i);
    info->Total_patch = temp_i;
    printf("Total patch = %d\n", temp_i);

    fgets(s, 256, fp);

    // ヤング率, ポアソン比
    for (i = 0; i < 2; i++)
    {
        fscanf(fp, "%lf", &temp_d);
        info->E_and_nu[i] = temp_d;
        printf("E_and_nu[%d] = %le\n", i, E_and_nu[i]);
    }

    fgets(s, 256, fp);

    // 各方向への変位指定する個数
    int temp_MAX_DISP_CONSTRAINT = 0, temp_MAX_DISP_CONSTRAINT_EDGE = 0;
    for (i = 0; i < info->DIMENSION; i++)
    {
        fscanf(fp, "%d", &temp_i);
        if (temp_MAX_DISP_CONSTRAINT < temp_i)
        {
            temp_MAX_DISP_CONSTRAINT = temp_i;
        }

        for (j = 0; j < disp_constraint_n[i]; j++)
        {
            fscanf(fp, "%d", &temp_i);
            if (temp_MAX_DISP_CONSTRAINT_EDGE < temp_i)
            {
                temp_MAX_DISP_CONSTRAINT_EDGE = temp_i;
            }

            fscanf(fp, "%lf", &temp_d);

            for (k = 0; k < disp_constraint_edge_n[i][j]; k++)
            {
                fscanf(fp, "%d", &temp_i);
                fscanf(fp, "%d", &temp_i);
                fscanf(fp, "%d", &temp_i);
            }

            fgets(s, 256, fp);
        }
    }
    info.MAX_DISP_CONSTRAINT = temp_MAX_DISP_CONSTRAINT;
    info.MAX_DISP_CONSTRAINT_EDGE = temp_MAX_DISP_CONSTRAINT_EDGE;

    // 分布荷重の荷重の個数
    fscanf(fp, "%d", &temp_i);
    infoglo->distributed_load_n = temp_i;

    fclose(fp);
}


void Get_inputdata_boundary_1(char *filename, information *info)
{
    int i, j, k;
    char s[256];

    int temp_i;
    double temp_d;

    fp = fopen(filename, "r");

    // DIMENSION
    fscanf(fp, "%d", &temp_i);

    fgets(s, 256, fp);

    // パッチ数
    fscanf(fp, "%d", &temp_i);

    fgets(s, 256, fp);

    // ヤング率, ポアソン比
    for (i = 0; i < info->DIMENSION; i++)
    {
        fscanf(fp, "%lf", &temp_d);
    }

    fgets(s, 256, fp);

    // x, y方向への変位指定する個数
    int a = info->MAX_DISP_CONSTRAINT * info->MAX_DISP_CONSTRAINT_EDGE * 3;
    int b = info->MAX_DISP_CONSTRAINT_EDGE * 3;
    int c = 3;
    for (i = 0; i < info->DIMENSION; i++)
    {
        fscanf(fp, "%d", &temp_i);
        info->disp_constraint_n[i] = temp_i;

        for (j = 0; j < info->disp_constraint_n[i]; j++)
        {
            fscanf(fp, "%d", &temp_i);
            info->disp_constraint_edge_n[i * info.MAX_DISP_CONSTRAINT + j] = temp_i;

            fscanf(fp, "%lf", &temp_d);
            info->disp_constraint_amount[i * info.MAX_DISP_CONSTRAINT + j] = temp_d;

            for (k = 0; k < disp_constraint_edge_n[i * info.MAX_DISP_CONSTRAINT + j]; k++)
            {
                fscanf(fp, "%d", &temp_i);
                info->disp_constraint[i * a + j * b + k * c + 0] = temp_i;

                fscanf(fp, "%d", &temp_i);
                info->disp_constraint[i * a + j * b + k * c + 1] = temp_i;

                fscanf(fp, "%d", &temp_i);
                info->disp_constraint[i * a + j * b + k * c + 2] = temp_i;
            }

            fgets(s, 256, fp);
        }
    }

    // 分布荷重の荷重の個数
    fscanf(fp, "%d", &temp_i);

    for (i = 0; i < info->distributed_load_n; i++)
    {
        for (j = 0; j < 9; j++)
        {
            fscanf(fp, "%lf", &temp_d);
            info->distributed_load_info[i * 9 + j] = temp_d;
        }
    }

    fclose(fp);
}


void Get_inputdata_patch_0(char *filename, information *info)
{
    int i;
    char s[256];

    int temp_counter;
    int temp_i;

    fp = fopen(filename, "r");

    // コントロールポイント数
    fscanf(fp, "%d", &temp_i);
    printf("%d\n", temp_i);

    fgets(s, 256, fp);

    // 次数
    temp_counter = counter;
    for (i = 0; i < info->DIMENSION; i++)
    {
        fscanf(fp, "%d", &temp_i);
        info->Order[temp_counter] = temp_i;
        printf("info->Order[%d] = %d\n", temp_counter, info->Order[temp_counter]);
        temp_counter++;
    }

    fgets(s, 256, fp);

    // ノットベクトルの個数
    temp_counter = counter;
    for (i = 0; i < info->DIMENSION; i++)
    {
        fscanf(fp, "%d", &temp_i);
        info->KV_info[temp_counter] = temp_i;
        printf("info->KV_info[%d] = %d\n", temp_counter, info->KV_info[temp_counter]);
        temp_counter++;
    }

    fgets(s, 256, fp);

    // 各方向のコントロールポイント数
    for (i = 0; i < info->DIMENSION; i++)
    {
        fscanf(fp, "%d", &temp_i);
        info->CP_info[counter] = temp_i;
        printf("info->CP_info[%d] = %d\n", counter, info->CP_info[counter]);
        counter++;
    }

    fclose(fp);
}


void Get_inputdata_patch_1(char *filename, information *info, int num)
{
    int i, j;
    char s[256];

    int temp_i;
    double temp_d;

    fp = fopen(filename, "r");

    // コントロールポイント数
    fscanf(fp, "%d", &temp_i);

    fgets(s, 256, fp);

    // 次数
    for (i = 0; i < info->DIMENSION; i++)
    {
        fscanf(fp, "%d", &temp_i);
    }

    fgets(s, 256, fp);

    // ノットベクトルの個数
    for (i = 0; i < info->DIMENSION; i++)
    {
        fscanf(fp, "%d", &temp_i);
    }

    fgets(s, 256, fp);

    // 各方向のコントロールポイント数
    for (i = 0; i < info->DIMENSION; i++)
    {
        fscanf(fp, "%d", &temp_i);
    }

    fgets(s, 256, fp);

    // ノットベクトル
    for (i = 0; i < info->DIMENSION; i++)
    {
        for (j = 0; j < info->KV_info[num * info->DIMENSION + i]; j++)
        {
            fscanf(fp, "%lf", &temp_d);
            info->KV[KV_to_here + j] = temp_d;
            printf("%le ", info->KV[KV_to_here + j]);
        }
        KV_to_here += info->KV_info[num * info->DIMENSION + i];
        printf("\n");
    }
    printf("\n");

    fgets(s, 256, fp);

    // コントロールポイント
    int temp_CP_to_here = CP_to_here;
    int Total_CP = 1;
    for (i = 0; i < info->DIMENSION; i++)
    {
        Total_CP *= CP_info[num * info->DIMENSION + i];
    }
    for (i = 0; i < Total_CP; i++)
    {
        fscanf(fp, "%d", &temp_i);
        for (j = 0; j < info->DIMENSION + 1; j++)
        {
            fscanf(fp, "%lf", &temp_d);
            info->CP[CP_to_here + j] = temp_d;
            printf("%le ", info->CP[CP_to_here + j]);
        }
        printf("\n");
        CP_to_here += info->DIMENSION + 1;
    }
    printf("\n");

    fclose(fp);

    // B 配列を作成
    if (info->DIMENSION == 2)
    {
        // 辺0 点0
        info->B[B_to_here]     = info->CP[temp_CP_to_here];
        info->B[B_to_here + 1] = info->CP[temp_CP_to_here + 1];
        info->B[B_to_here + 2] = info->CP[temp_CP_to_here + 2];
        B_to_here += info->DIMENSION + 1;

        // 辺0 点1
        info->B[B_to_here]     = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] - 1) * (info->DIMENSION + 1)];
        info->B[B_to_here + 1] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] - 1) * (info->DIMENSION + 1) + 1];
        info->B[B_to_here + 2] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] - 1) * (info->DIMENSION + 1) + 2];
        B_to_here += info->DIMENSION + 1;

        // 辺1 点0
        info->B[B_to_here]     = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] - 1) * (info->DIMENSION + 1)];
        info->B[B_to_here + 1] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] - 1) * (info->DIMENSION + 1) + 1];
        info->B[B_to_here + 2] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] - 1) * (info->DIMENSION + 1) + 2];
        B_to_here += info->DIMENSION + 1;

        // 辺1 点1
        info->B[B_to_here]     = info->CP[temp_CP_to_here];
        info->B[B_to_here + 1] = info->CP[temp_CP_to_here + 1];
        info->B[B_to_here + 2] = info->CP[temp_CP_to_here + 2];
        B_to_here += info->DIMENSION + 1;

        // 辺2 点0
        info->B[B_to_here]     = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] - 1) * (info->DIMENSION + 1)];
        info->B[B_to_here + 1] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] - 1) * (info->DIMENSION + 1) + 1];
        info->B[B_to_here + 2] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] - 1) * (info->DIMENSION + 1) + 2];
        B_to_here += info->DIMENSION + 1;

        // 辺2 点1
        info->B[B_to_here]     = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * info->CP_info[num * info->DIMENSION + 1] - 1) * (info->DIMENSION + 1)];
        info->B[B_to_here + 1] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * info->CP_info[num * info->DIMENSION + 1] - 1) * (info->DIMENSION + 1) + 1];
        info->B[B_to_here + 2] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * info->CP_info[num * info->DIMENSION + 1] - 1) * (info->DIMENSION + 1) + 2];
        B_to_here += info->DIMENSION + 1;

        // 辺3 点0
        info->B[B_to_here]     = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * info->CP_info[num * info->DIMENSION + 1] - 1) * (info->DIMENSION + 1)];
        info->B[B_to_here + 1] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * info->CP_info[num * info->DIMENSION + 1] - 1) * (info->DIMENSION + 1) + 1];
        info->B[B_to_here + 2] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * info->CP_info[num * info->DIMENSION + 1] - 1) * (info->DIMENSION + 1) + 2];
        B_to_here += info->DIMENSION + 1;

        // 辺3 点1
        info->B[B_to_here]     = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] - 1) * (info->DIMENSION + 1)];
        info->B[B_to_here + 1] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] - 1) * (info->DIMENSION + 1) + 1];
        info->B[B_to_here + 2] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] - 1) * (info->DIMENSION + 1) + 2];
        B_to_here += info->DIMENSION + 1;

        // 辺4 点0
        info->B[B_to_here]     = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * (info->CP_info[num * info->DIMENSION + 1] - 1)) * (info->DIMENSION + 1)];
        info->B[B_to_here + 1] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * (info->CP_info[num * info->DIMENSION + 1] - 1)) * (info->DIMENSION + 1) + 1];
        info->B[B_to_here + 2] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * (info->CP_info[num * info->DIMENSION + 1] - 1)) * (info->DIMENSION + 1) + 2];
        B_to_here += info->DIMENSION + 1;

        // 辺4 点1
        info->B[B_to_here]     = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * info->CP_info[num * info->DIMENSION + 1] - 1) * (info->DIMENSION + 1)];
        info->B[B_to_here + 1] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * info->CP_info[num * info->DIMENSION + 1] - 1) * (info->DIMENSION + 1) + 1];
        info->B[B_to_here + 2] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * info->CP_info[num * info->DIMENSION + 1] - 1) * (info->DIMENSION + 1) + 2];
        B_to_here += info->DIMENSION + 1;

        // 辺5 点0
        info->B[B_to_here]     = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * info->CP_info[num * info->DIMENSION + 1] - 1) * (info->DIMENSION + 1)];
        info->B[B_to_here + 1] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * info->CP_info[num * info->DIMENSION + 1] - 1) * (info->DIMENSION + 1) + 1];
        info->B[B_to_here + 2] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * info->CP_info[num * info->DIMENSION + 1] - 1) * (info->DIMENSION + 1) + 2];
        B_to_here += info->DIMENSION + 1;

        // 辺5 点1
        info->B[B_to_here]     = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * (info->CP_info[num * info->DIMENSION + 1] - 1)) * (info->DIMENSION + 1)];
        info->B[B_to_here + 1] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * (info->CP_info[num * info->DIMENSION + 1] - 1)) * (info->DIMENSION + 1) + 1];
        info->B[B_to_here + 2] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * (info->CP_info[num * info->DIMENSION + 1] - 1)) * (info->DIMENSION + 1) + 2];
        B_to_here += info->DIMENSION + 1;

        // 辺6 点0
        info->B[B_to_here]     = info->CP[temp_CP_to_here];
        info->B[B_to_here + 1] = info->CP[temp_CP_to_here + 1];
        info->B[B_to_here + 2] = info->CP[temp_CP_to_here + 2];
        B_to_here += info->DIMENSION + 1;

        // 辺6 点1
        info->B[B_to_here]     = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * (info->CP_info[num * info->DIMENSION + 1] - 1)) * (info->DIMENSION + 1)];
        info->B[B_to_here + 1] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * (info->CP_info[num * info->DIMENSION + 1] - 1)) * (info->DIMENSION + 1) + 1];
        info->B[B_to_here + 2] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * (info->CP_info[num * info->DIMENSION + 1] - 1)) * (info->DIMENSION + 1) + 2];
        B_to_here += info->DIMENSION + 1;

        // 辺7 点0
        info->B[B_to_here]     = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * (info->CP_info[num * info->DIMENSION + 1] - 1)) * (info->DIMENSION + 1)];
        info->B[B_to_here + 1] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * (info->CP_info[num * info->DIMENSION + 1] - 1)) * (info->DIMENSION + 1) + 1];
        info->B[B_to_here + 2] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * (info->CP_info[num * info->DIMENSION + 1] - 1)) * (info->DIMENSION + 1) + 2];
        B_to_here += info->DIMENSION + 1;

        // 辺7 点1
        info->B[B_to_here]     = info->CP[temp_CP_to_here];
        info->B[B_to_here + 1] = info->CP[temp_CP_to_here + 1];
        info->B[B_to_here + 2] = info->CP[temp_CP_to_here + 2];
        B_to_here += info->DIMENSION + 1;
    }
    else if (info->DIMENSION == 3)
    {
        int temp = info->CP_info[num * info->DIMENSION] * info->CP_info[num * info->DIMENSION + 1] * (info->CP_info[num * info->DIMENSION + 2] - 1);
        int temp_B_to_here;
        int temp_point_B[8];

        // 面0 点0
        temp_point_B[0] = B_to_here;
        info->B[B_to_here]     = info->CP[temp_CP_to_here];
        info->B[B_to_here + 1] = info->CP[temp_CP_to_here + 1];
        info->B[B_to_here + 2] = info->CP[temp_CP_to_here + 2];
        info->B[B_to_here + 3] = info->CP[temp_CP_to_here + 3];
        B_to_here += info->DIMENSION + 1;

        // 面0 点1
        temp_point_B[1] = B_to_here;
        info->B[B_to_here]     = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] - 1) * (info->DIMENSION + 1)];
        info->B[B_to_here + 1] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] - 1) * (info->DIMENSION + 1) + 1];
        info->B[B_to_here + 2] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] - 1) * (info->DIMENSION + 1) + 2];
        info->B[B_to_here + 3] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] - 1) * (info->DIMENSION + 1) + 3];
        B_to_here += info->DIMENSION + 1;

        // 面0 点2
        temp_point_B[2] = B_to_here;
        info->B[B_to_here]     = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * info->CP_info[num * info->DIMENSION + 1] - 1) * (info->DIMENSION + 1)];
        info->B[B_to_here + 1] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * info->CP_info[num * info->DIMENSION + 1] - 1) * (info->DIMENSION + 1) + 1];
        info->B[B_to_here + 2] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * info->CP_info[num * info->DIMENSION + 1] - 1) * (info->DIMENSION + 1) + 2];
        info->B[B_to_here + 3] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * info->CP_info[num * info->DIMENSION + 1] - 1) * (info->DIMENSION + 1) + 3];
        B_to_here += info->DIMENSION + 1;

        // 面0 点3
        temp_point_B[3] = B_to_here;
        info->B[B_to_here]     = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * (info->CP_info[num * info->DIMENSION + 1] - 1)) * (info->DIMENSION + 1)];
        info->B[B_to_here + 1] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * (info->CP_info[num * info->DIMENSION + 1] - 1)) * (info->DIMENSION + 1) + 1];
        info->B[B_to_here + 2] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * (info->CP_info[num * info->DIMENSION + 1] - 1)) * (info->DIMENSION + 1) + 2];
        info->B[B_to_here + 3] = info->CP[temp_CP_to_here + (info->CP_info[num * info->DIMENSION] * (info->CP_info[num * info->DIMENSION + 1] - 1)) * (info->DIMENSION + 1) + 3];
        B_to_here += info->DIMENSION + 1;

        // 面1 点0
        temp_point_B[4] = B_to_here;
        info->B[B_to_here]     = info->CP[temp_CP_to_here + temp];
        info->B[B_to_here + 1] = info->CP[temp_CP_to_here + temp + 1];
        info->B[B_to_here + 2] = info->CP[temp_CP_to_here + temp + 2];
        info->B[B_to_here + 3] = info->CP[temp_CP_to_here + temp + 3];
        B_to_here += info->DIMENSION + 1;

        // 面1 点1
        temp_point_B[5] = B_to_here;
        info->B[B_to_here]     = info->CP[temp_CP_to_here + temp + (info->CP_info[num * info->DIMENSION] - 1) * (info->DIMENSION + 1)];
        info->B[B_to_here + 1] = info->CP[temp_CP_to_here + temp + (info->CP_info[num * info->DIMENSION] - 1) * (info->DIMENSION + 1) + 1];
        info->B[B_to_here + 2] = info->CP[temp_CP_to_here + temp + (info->CP_info[num * info->DIMENSION] - 1) * (info->DIMENSION + 1) + 2];
        info->B[B_to_here + 3] = info->CP[temp_CP_to_here + temp + (info->CP_info[num * info->DIMENSION] - 1) * (info->DIMENSION + 1) + 3];
        B_to_here += info->DIMENSION + 1;

        // 面1 点2
        temp_point_B[6] = B_to_here;
        info->B[B_to_here]     = info->CP[temp_CP_to_here + temp + (info->CP_info[num * info->DIMENSION] * info->CP_info[num * info->DIMENSION + 1] - 1) * (info->DIMENSION + 1)];
        info->B[B_to_here + 1] = info->CP[temp_CP_to_here + temp + (info->CP_info[num * info->DIMENSION] * info->CP_info[num * info->DIMENSION + 1] - 1) * (info->DIMENSION + 1) + 1];
        info->B[B_to_here + 2] = info->CP[temp_CP_to_here + temp + (info->CP_info[num * info->DIMENSION] * info->CP_info[num * info->DIMENSION + 1] - 1) * (info->DIMENSION + 1) + 2];
        info->B[B_to_here + 3] = info->CP[temp_CP_to_here + temp + (info->CP_info[num * info->DIMENSION] * info->CP_info[num * info->DIMENSION + 1] - 1) * (info->DIMENSION + 1) + 3];
        B_to_here += info->DIMENSION + 1;

        // 面1 点3
        temp_point_B[7] = B_to_here;
        info->B[B_to_here]     = info->CP[temp_CP_to_here + temp + (info->CP_info[num * info->DIMENSION] * (info->CP_info[num * info->DIMENSION + 1] - 1)) * (info->DIMENSION + 1)];
        info->B[B_to_here + 1] = info->CP[temp_CP_to_here + temp + (info->CP_info[num * info->DIMENSION] * (info->CP_info[num * info->DIMENSION + 1] - 1)) * (info->DIMENSION + 1) + 1];
        info->B[B_to_here + 2] = info->CP[temp_CP_to_here + temp + (info->CP_info[num * info->DIMENSION] * (info->CP_info[num * info->DIMENSION + 1] - 1)) * (info->DIMENSION + 1) + 2];
        info->B[B_to_here + 3] = info->CP[temp_CP_to_here + temp + (info->CP_info[num * info->DIMENSION] * (info->CP_info[num * info->DIMENSION + 1] - 1)) * (info->DIMENSION + 1) + 3];
        B_to_here += info->DIMENSION + 1;

        // 面2 点0
        temp_B_to_here = temp_point_B[0];
        info->B[B_to_here]     = info->B[temp_B_to_here];
        info->B[B_to_here + 1] = info->B[temp_B_to_here + 1];
        info->B[B_to_here + 2] = info->B[temp_B_to_here + 2];
        info->B[B_to_here + 3] = info->B[temp_B_to_here + 3];
        B_to_here += info->DIMENSION + 1;

        // 面2 点1
        temp_B_to_here = temp_point_B[3];
        info->B[B_to_here]     = info->B[temp_B_to_here];
        info->B[B_to_here + 1] = info->B[temp_B_to_here + 1];
        info->B[B_to_here + 2] = info->B[temp_B_to_here + 2];
        info->B[B_to_here + 3] = info->B[temp_B_to_here + 3];
        B_to_here += info->DIMENSION + 1;

        // 面2 点2
        temp_B_to_here = temp_point_B[7];
        info->B[B_to_here]     = info->B[temp_B_to_here];
        info->B[B_to_here + 1] = info->B[temp_B_to_here + 1];
        info->B[B_to_here + 2] = info->B[temp_B_to_here + 2];
        info->B[B_to_here + 3] = info->B[temp_B_to_here + 3];
        B_to_here += info->DIMENSION + 1;
        
        // 面2 点3
        temp_B_to_here = temp_point_B[4];
        info->B[B_to_here]     = info->B[temp_B_to_here];
        info->B[B_to_here + 1] = info->B[temp_B_to_here + 1];
        info->B[B_to_here + 2] = info->B[temp_B_to_here + 2];
        info->B[B_to_here + 3] = info->B[temp_B_to_here + 3];
        B_to_here += info->DIMENSION + 1;

        // 面3 点0
        temp_B_to_here = temp_point_B[0];
        info->B[B_to_here]     = info->B[temp_B_to_here];
        info->B[B_to_here + 1] = info->B[temp_B_to_here + 1];
        info->B[B_to_here + 2] = info->B[temp_B_to_here + 2];
        info->B[B_to_here + 3] = info->B[temp_B_to_here + 3];
        B_to_here += info->DIMENSION + 1;

        // 面3 点1
        temp_B_to_here = temp_point_B[1];
        info->B[B_to_here]     = info->B[temp_B_to_here];
        info->B[B_to_here + 1] = info->B[temp_B_to_here + 1];
        info->B[B_to_here + 2] = info->B[temp_B_to_here + 2];
        info->B[B_to_here + 3] = info->B[temp_B_to_here + 3];
        B_to_here += info->DIMENSION + 1;

        // 面3 点2
        temp_B_to_here = temp_point_B[5];
        info->B[B_to_here]     = info->B[temp_B_to_here];
        info->B[B_to_here + 1] = info->B[temp_B_to_here + 1];
        info->B[B_to_here + 2] = info->B[temp_B_to_here + 2];
        info->B[B_to_here + 3] = info->B[temp_B_to_here + 3];
        B_to_here += info->DIMENSION + 1;

        // 面3 点3
        temp_B_to_here = temp_point_B[4];
        info->B[B_to_here]     = info->B[temp_B_to_here];
        info->B[B_to_here + 1] = info->B[temp_B_to_here + 1];
        info->B[B_to_here + 2] = info->B[temp_B_to_here + 2];
        info->B[B_to_here + 3] = info->B[temp_B_to_here + 3];
        B_to_here += info->DIMENSION + 1;

        // 面4 点0
        temp_B_to_here = temp_point_B[1];
        info->B[B_to_here]     = info->B[temp_B_to_here];
        info->B[B_to_here + 1] = info->B[temp_B_to_here + 1];
        info->B[B_to_here + 2] = info->B[temp_B_to_here + 2];
        info->B[B_to_here + 3] = info->B[temp_B_to_here + 3];
        B_to_here += info->DIMENSION + 1;

        // 面4 点1
        temp_B_to_here = temp_point_B[2];
        info->B[B_to_here]     = info->B[temp_B_to_here];
        info->B[B_to_here + 1] = info->B[temp_B_to_here + 1];
        info->B[B_to_here + 2] = info->B[temp_B_to_here + 2];
        info->B[B_to_here + 3] = info->B[temp_B_to_here + 3];
        B_to_here += info->DIMENSION + 1;

        // 面4 点2
        temp_B_to_here = temp_point_B[6];
        info->B[B_to_here]     = info->B[temp_B_to_here];
        info->B[B_to_here + 1] = info->B[temp_B_to_here + 1];
        info->B[B_to_here + 2] = info->B[temp_B_to_here + 2];
        info->B[B_to_here + 3] = info->B[temp_B_to_here + 3];
        B_to_here += info->DIMENSION + 1;

        // 面4 点3
        temp_B_to_here = temp_point_B[5];
        info->B[B_to_here]     = info->B[temp_B_to_here];
        info->B[B_to_here + 1] = info->B[temp_B_to_here + 1];
        info->B[B_to_here + 2] = info->B[temp_B_to_here + 2];
        info->B[B_to_here + 3] = info->B[temp_B_to_here + 3];
        B_to_here += info->DIMENSION + 1;

        // 面5 点0
        temp_B_to_here = temp_point_B[3];
        info->B[B_to_here]     = info->B[temp_B_to_here];
        info->B[B_to_here + 1] = info->B[temp_B_to_here + 1];
        info->B[B_to_here + 2] = info->B[temp_B_to_here + 2];
        info->B[B_to_here + 3] = info->B[temp_B_to_here + 3];
        B_to_here += info->DIMENSION + 1;

        // 面5 点1
        temp_B_to_here = temp_point_B[2];
        info->B[B_to_here]     = info->B[temp_B_to_here];
        info->B[B_to_here + 1] = info->B[temp_B_to_here + 1];
        info->B[B_to_here + 2] = info->B[temp_B_to_here + 2];
        info->B[B_to_here + 3] = info->B[temp_B_to_here + 3];
        B_to_here += info->DIMENSION + 1;

        // 面5 点2
        temp_B_to_here = temp_point_B[6];
        info->B[B_to_here]     = info->B[temp_B_to_here];
        info->B[B_to_here + 1] = info->B[temp_B_to_here + 1];
        info->B[B_to_here + 2] = info->B[temp_B_to_here + 2];
        info->B[B_to_here + 3] = info->B[temp_B_to_here + 3];
        B_to_here += info->DIMENSION + 1;

        // 面5 点3
        temp_B_to_here = temp_point_B[7];
        info->B[B_to_here]     = info->B[temp_B_to_here];
        info->B[B_to_here + 1] = info->B[temp_B_to_here + 1];
        info->B[B_to_here + 2] = info->B[temp_B_to_here + 2];
        info->B[B_to_here + 3] = info->B[temp_B_to_here + 3];
        B_to_here += info->DIMENSION + 1;
    }
}


// make connectivity
void Check_B(int num_own, int num_opponent, double *temp_B, int *temp_Edge_info, int *temp_Opponent_patch_num)
{
    int i, j;
    int ii;
    int x_diff[2], y_diff[2], w_diff[2];

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

            // 辺が一致している場合 Face_Edge_info を True
            if (sqrt(pow(x_diff[0], 2) + pow(y_diff[0], 2) + pow(w_diff[0], 2)) <= MERGE_DISTANCE && sqrt(pow(x_diff[1], 2) + pow(y_diff[1], 2) + pow(w_diff[1], 2)) <= MERGE_DISTANCE)
            {
                temp_Edge_info[num_own * 32 + i * 8 + j] = 1;
                temp_Opponent_patch_num[num_own * 4 + i] = num_opponent;
                printf("own_patch:%d opp_patch:%d own_edge:%d opp_edge:%d\n", num_own, num_opponent, i, j);
                return;
            }
        }
    }
}


void Make_connectivity(int num, int *info->CP_info, int *temp_Edge_info, int *temp_Opponent_patch_num, int *temp_Connectivity, int *temp_A, double *temp_CP, double *temp_CP_result)
{
    int i, j, k;
    int p, q;
    int temp_CP_n = 0;
    int Edge[4];
    int Edge_to_here = num * 32;

    printf("make connectivity on patch %d\n", num);

    A_to_own = 0;
    for (i = 0; i < num; i++)
    {
        A_to_own += 2 * (info->CP_info[i * DIMENSION] + info->CP_info[i * DIMENSION + 1]);
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
            A_to_own += info->CP_info[num * DIMENSION];
        }
        else if (i == 2)
        {
            A_to_own += info->CP_info[num * DIMENSION + 1];
        }
        else if (i == 3)
        {
            A_to_own += info->CP_info[num * DIMENSION];
        }

        for (j = 0; j < 8; j++)
        {
            if (temp_Edge_info[Edge_to_here + i * 8 + j] == 1)
            {
                A_to_opponent = 0;
                for (k = 0; k < temp_Opponent_patch_num[num * 4 + i]; k++)
                {
                    A_to_opponent += 2 * (info->CP_info[k * DIMENSION] + info->CP_info[k * DIMENSION + 1]);
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
                    temp_CP_n = info->CP_info[temp_Opponent_patch_num[num * 4 + i] * DIMENSION];
                }
                else if (p == 1)
                {
                    temp_CP_n = info->CP_info[temp_Opponent_patch_num[num * 4 + i] * DIMENSION + 1];
                    A_to_opponent += info->CP_info[temp_Opponent_patch_num[num * 4 + i] * DIMENSION];
                }
                else if (p == 2)
                {
                    temp_CP_n = info->CP_info[temp_Opponent_patch_num[num * 4 + i] * DIMENSION];
                    A_to_opponent += info->CP_info[temp_Opponent_patch_num[num * 4 + i] * DIMENSION] + info->CP_info[temp_Opponent_patch_num[num * 4 + i] * DIMENSION + 1];
                }
                else if (p == 3)
                {
                    temp_CP_n = info->CP_info[temp_Opponent_patch_num[num * 4 + i] * DIMENSION + 1];
                    A_to_opponent += 2 * info->CP_info[temp_Opponent_patch_num[num * 4 + i] * DIMENSION] + info->CP_info[temp_Opponent_patch_num[num * 4 + i] * DIMENSION + 1];
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
        A_to_own += 2 * (info->CP_info[i * DIMENSION] + info->CP_info[i * DIMENSION + 1]);
    }
    for (i = 0; i < 2 * (info->CP_info[num * DIMENSION] + info->CP_info[num * DIMENSION + 1]); i++)
    {
        printf("%d\t", temp_A[A_to_own + i]);
    }
    printf("\n");

    // コネクティビティを作成
    int xi, eta;

    A_to_own = 0;
    for (i = 0; i < num; i++)
    {
        A_to_own += 2 * (info->CP_info[i * DIMENSION] + info->CP_info[i * DIMENSION + 1]);
    }

    for (eta = 0; eta < info->CP_info[num * DIMENSION + 1]; eta++)
    {
        for (xi = 0; xi < info->CP_info[num * DIMENSION]; xi++)
        {
            if (eta == 0 && Edge[0] == 1)
            {
                temp_Connectivity[CP_to_here + eta * info->CP_info[num * DIMENSION] + xi] = temp_A[A_to_own + xi];
            }
            else if (eta == info->CP_info[num * DIMENSION + 1] - 1 && Edge[2] == 1)
            {
                temp_Connectivity[CP_to_here + eta * info->CP_info[num * DIMENSION] + xi] = temp_A[A_to_own + info->CP_info[num * DIMENSION] + info->CP_info[num * DIMENSION + 1] + xi];
            }
            else if (xi == 0 && Edge[3] == 1)
            {
                temp_Connectivity[CP_to_here + eta * info->CP_info[num * DIMENSION] + xi] = temp_A[A_to_own + 2 * info->CP_info[num * DIMENSION] + info->CP_info[num * DIMENSION + 1] + eta];
            }
            else if (xi == info->CP_info[num * DIMENSION] - 1 && Edge[1] == 1)
            {
                temp_Connectivity[CP_to_here + eta * info->CP_info[num * DIMENSION] + xi] = temp_A[A_to_own + info->CP_info[num * DIMENSION] + eta];
            }
            else
            {
                temp_Connectivity[CP_to_here + eta * info->CP_info[num * DIMENSION] + xi] = counter;
                temp_CP_result[CP_result_to_here] = temp_CP[(CP_to_here + eta * info->CP_info[num * DIMENSION] + xi) * (DIMENSION + 1)];
                temp_CP_result[CP_result_to_here + 1] = temp_CP[(CP_to_here + eta * info->CP_info[num * DIMENSION] + xi) * (DIMENSION + 1) + 1];
                temp_CP_result[CP_result_to_here + 2] = temp_CP[(CP_to_here + eta * info->CP_info[num * DIMENSION] + xi) * (DIMENSION + 1) + 2];
                counter++;
                CP_result_to_here += (DIMENSION + 1);
            }
        }
    }

    // A 配列の作ってない分を作成
    for (eta = 0; eta < info->CP_info[num * DIMENSION + 1]; eta++)
    {
        for (xi = 0; xi < info->CP_info[num * DIMENSION]; xi++)
        {
            if (eta == 0 && Edge[0] == 0)
            {
                temp_A[A_to_own + xi] = temp_Connectivity[CP_to_here + eta * info->CP_info[num * DIMENSION] + xi];
            }
            if (eta == info->CP_info[num * DIMENSION + 1] - 1 && Edge[2] == 0)
            {
                temp_A[A_to_own + info->CP_info[num * DIMENSION] + info->CP_info[num * DIMENSION + 1] + xi] = temp_Connectivity[CP_to_here + eta * info->CP_info[num * DIMENSION] + xi];
            }
            if (xi == 0 && Edge[3] == 0)
            {
                temp_A[A_to_own + 2 * info->CP_info[num * DIMENSION] + info->CP_info[num * DIMENSION + 1] + eta] = temp_Connectivity[CP_to_here + eta * info->CP_info[num * DIMENSION] + xi];
            }
            if (xi == info->CP_info[num * DIMENSION] - 1 && Edge[1] == 0)
            {
                temp_A[A_to_own + info->CP_info[num * DIMENSION] + eta] = temp_Connectivity[CP_to_here + eta * info->CP_info[num * DIMENSION] + xi];
            }
        }
    }
    CP_to_here += info->CP_info[num * DIMENSION] * info->CP_info[num * DIMENSION + 1];

    printf("patch %d array A after make connectivity\n", num);
    A_to_own = 0;
    for (i = 0; i < num; i++)
    {
        A_to_own += 2 * (info->CP_info[i * DIMENSION] + info->CP_info[i * DIMENSION + 1]);
    }
    for (i = 0; i < 2 * (info->CP_info[num * DIMENSION] + info->CP_info[num * DIMENSION + 1]); i++)
    {
        printf("%d\t", temp_A[A_to_own + i]);
    }
    printf("\n");

    printf("\n");
}


// output
void Output_inputdata(int *temp_Order, int *temp_KV_info, int *info->CP_info, int *temp_Connectivity, double *temp_KV, double *temp_CP_result,
                      int *temp_Boundary_result, int *temp_length_before, int *temp_length_after, int total_disp_constraint_n)
{
    int i, j, k;
    char str[256] = "input.txt";

    fp = fopen(str, "w");

    // ヤング率
    fprintf(fp, "%d  ", (int)E_and_nu[0]);

    // ポアソン比
    fprintf(fp, "%le\n\n", E_and_nu[1]);

    // パッチ数
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
                fprintf(fp, "%d", info->Order[i * DIMENSION + j]);
            }
            else
            {
                fprintf(fp, "%*d", -6, info->Order[i * DIMENSION + j]);
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
                fprintf(fp, "%d", info->KV_info[i * DIMENSION + j]);
            }
            else
            {
                fprintf(fp, "%*d", -6, info->KV_info[i * DIMENSION + j]);
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
                fprintf(fp, "%d", info->CP_info[i * DIMENSION + j]);
            }
            else
            {
                fprintf(fp, "%*d", -6, info->CP_info[i * DIMENSION + j]);
            }
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // パッチコネクティビティ
    CP_to_here = 0;
    for (i = 0; i < Total_patch; i++)
    {
        for (j = 0; j < info->CP_info[i * DIMENSION] * info->CP_info[i * DIMENSION + 1]; j++)
        {
            if (j == info->CP_info[i * DIMENSION] * info->CP_info[i * DIMENSION + 1] - 1)
            {
                fprintf(fp, "%d", temp_Connectivity[CP_to_here + j]);
            }
            else
            {
                fprintf(fp, "%*d", temp_num, temp_Connectivity[CP_to_here + j]);
            }
        }
        CP_to_here += info->CP_info[i * DIMENSION] * info->CP_info[i * DIMENSION + 1];
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // 変位拘束するコントロールポイントの数
    int temp = 0;
    for (i = 0; i < total_disp_constraint_n; i++)
    {
        temp += temp_length_after[i];
    }
    fprintf(fp, "%*d", -6, temp);

    // 荷重条件を与えるコントロールポイントの数
    fprintf(fp, "%*d", -6, 0);

    // 分布荷重の数
    fprintf(fp, "%d\n\n", distributed_load_n);

    // 各パッチでの各方向のノットベクトル
    KV_to_here = 0;
    for (i = 0; i < Total_patch; i++)
    {
        for (j = 0; j < DIMENSION; j++)
        {
            for (k = 0; k < info->KV_info[i * DIMENSION + j]; k++)
            {
                if (k == info->KV_info[i * DIMENSION + j] - 1)
                {
                    fprintf(fp, "%.16e", temp_KV[KV_to_here + k]);
                }
                else
                {
                    fprintf(fp, "%.16e  ", temp_KV[KV_to_here + k]);
                }
            }
            KV_to_here += info->KV_info[i * DIMENSION + j];
            fprintf(fp, "\n");
        }
    }
    fprintf(fp, "\n");

    // コントロールポイント
    for (i = 0; i < (CP_result_to_here + 1) / 3; i++)
    {
        fprintf(fp, "%*d", temp_num, i);
        fprintf(fp, "%.16e  ", temp_CP_result[i * 3]);
        fprintf(fp, "%.16e  ", temp_CP_result[i * 3 + 1]);
        fprintf(fp, "%.16e\n", temp_CP_result[i * 3 + 2]);
    }
    fprintf(fp, "\n");

    // 拘束するコントロールポイント
    temp_counter = 0;
    temp = 0;
    for (i = 0; i < DIMENSION; i++)
    {
        for (j = 0; j < disp_constraint_n[i]; j++)
        {
            for (k = 0; k < temp_length_after[temp_counter]; k++)
            {
                fprintf(fp, "%*d", temp_num, temp_Boundary_result[temp + k]);
                fprintf(fp, "%*d", temp_num, i);
                fprintf(fp, "%le\n", disp_constraint_amount[i][j]);
            }
            temp += temp_length_before[j];
            temp_counter++;
        }
    }
    fprintf(fp, "\n");

    // 分布荷重
    for (i = 0; i < distributed_load_n; i++)
    {
        if (i != 0)
        {
            fprintf(fp, "\n");
        }
        fprintf(fp, "%*d", temp_num, (int)distributed_load_info[i][0]);
        fprintf(fp, "%*d", temp_num, (int)distributed_load_info[i][1]);
        fprintf(fp, "%*d", temp_num, (int)distributed_load_info[i][2]);
        fprintf(fp, "%le  ", distributed_load_info[i][3]);
        fprintf(fp, "%le  ", distributed_load_info[i][4]);
        fprintf(fp, "%le  ", distributed_load_info[i][5]);
        fprintf(fp, "%le  ", distributed_load_info[i][6]);
        fprintf(fp, "%le  ", distributed_load_info[i][7]);
        fprintf(fp, "%le", distributed_load_info[i][8]);
    }

    fclose(fp);
}


void Output_SVG(double *temp_B, double *temp_CP_result)
{
    int i;

    // char color_vec[11][10] = {"#a9a9a9", "#00bfff", "#00fa9a", "#bdb76b", "#ffff00", "#ff8c00", "#cd5c5c", "#ff7f50", "#dc143c", "#ee82ee", "#8a2be2"};
    // //  0   darkgray
    // //  1   deepskyblue
    // //  2   mediumspringgreen
    // //  3   darkkhaki
    // //  4   yellow
    // //  5   darkorange
    // //  6   indianred
    // //  7   coral
    // //  8   crimson
    // //  9   violet
    // //  10  blueviolet
    // //  https://www.colordic.org/

    char color_vec[10][10] = {"#a9a9a9", "#00bfff", "#00fa9a", "#ffff00", "#ff8c00", "#cd5c5c", "#ff7f50", "#dc143c", "#ee82ee", "#8a2be2"};
    //  0   darkgray
    //  1   deepskyblue
    //  2   mediumspringgreen
    //  3   yellow
    //  4   darkorange
    //  5   indianred
    //  6   coral
    //  7   crimson
    //  8   violet
    //  9   blueviolet
    //  https://www.colordic.org/

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

    double space = 3.0;

    double scale = 2000.0 / (x_max - x_min + 2.0 * space);
    // widthが2000になるよう拡大，縮小する
    // 数字や文字が小さくてつぶれる場合はこの値を大きくするとイイ！！
    // それか line 1287 の font-size を小さくするとか

    double width = (x_max - x_min + 2.0 * space) * scale;
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
    B_to_here = 0;
    for (i = 0; i < Total_patch; i++)
    {
        position_x = (temp_B[B_to_here] + space) * scale;
        position_y = height - ((temp_B[B_to_here + 1] + space) * scale);
        fprintf(fp, "<path d='M %le %le ", position_x, position_y);
        B_to_here += 4 * (DIMENSION + 1);

        position_x = (temp_B[B_to_here] + space) * scale;
        position_y = height - ((temp_B[B_to_here + 1] + space) * scale);
        fprintf(fp, "L %le %le ", position_x, position_y);
        B_to_here += 2 * (DIMENSION + 1);

        position_x = (temp_B[B_to_here] + space) * scale;
        position_y = height - ((temp_B[B_to_here + 1] + space) * scale);
        fprintf(fp, "L %le %le ", position_x, position_y);
        B_to_here += 2 * (DIMENSION + 1);

        position_x = (temp_B[B_to_here] + space) * scale;
        position_y = height - ((temp_B[B_to_here + 1] + space) * scale);
        fprintf(fp, "L %le %le ", position_x, position_y);
        B_to_here += 4 * (DIMENSION + 1);

        fprintf(fp, "Z' fill='%s'/>\n", color_vec[temp_color_num % 10]);
        B_to_here += 4 * (DIMENSION + 1);

        if (temp_color_num % 10 == 5)
        {
            temp_color_num += 2;
        }
        else
        {
            temp_color_num++;
        }
    }

    // 点と番号を描画
    for (i = 0; i < (CP_result_to_here + 1) / 3; i++)
    {
        position_x = (temp_CP_result[i * 3] + space) * scale;
        position_y = height - ((temp_CP_result[i * 3 + 1] + space) * scale);
        fprintf(fp, "<circle cx='%le' cy='%le' r='2' fill='%s'/>\n", position_x, position_y, color_vec[7]);
        fprintf(fp, "<text x='%le' y='%le' font-family='Verdana' font-size='6' fill='%s' font-weight='700'>\n", position_x + 2, position_y, color_vec[7]);
        fprintf(fp, "%d\n", i);
        fprintf(fp, "</text>\n");
    }

    fprintf(fp, "</svg>");
    fclose(fp);
}


// heap sort
void Sort(int n, int *info->CP_info, int *temp_A, int *temp_Boundary, int *temp_Boundary_result, int *temp_length_before, int *temp_length_after)
{
    int i, j, k, l;
    int temp = 0;

    for (i = 0; i < DIMENSION; i++)
    {
        for (j = 0; j < disp_constraint_n[i]; j++)
        {
            for (k = 0; k < disp_constraint_edge_n[i][j]; k++)
            {
                if (disp_constraint[i][j][k][1] == 0)
                {
                    temp_length_before[temp] += info->CP_info[disp_constraint[i][j][k][0] * DIMENSION + 1];
                }
                else if (disp_constraint[i][j][k][1] == 1)
                {
                    temp_length_before[temp] += info->CP_info[disp_constraint[i][j][k][0] * DIMENSION];
                }
            }
            temp++;
        }
    }

    temp = 0;

    for (i = 0; i < DIMENSION; i++)
    {
        for (j = 0; j < disp_constraint_n[i]; j++)
        {
            for (k = 0; k < disp_constraint_edge_n[i][j]; k++)
            {
                int A_to_here = 0;
                for (l = 0; l < disp_constraint[i][j][k][0]; l++)
                {
                    A_to_here += 2 * (info->CP_info[l * DIMENSION] + info->CP_info[l * DIMENSION + 1]);
                }

                // if (disp_constraint[i][j][k][1] == 1 && disp_constraint[i][j][k][2] == 0) は何もしない
                if (disp_constraint[i][j][k][1] == 0 && disp_constraint[i][j][k][2] == 1)
                {
                    A_to_here += info->CP_info[disp_constraint[i][j][k][0] * DIMENSION];
                }
                else if (disp_constraint[i][j][k][1] == 1 && disp_constraint[i][j][k][2] == 1)
                {
                    A_to_here += info->CP_info[disp_constraint[i][j][k][0] * DIMENSION] + info->CP_info[disp_constraint[i][j][k][0] * DIMENSION + 1];
                }
                else if (disp_constraint[i][j][k][1] == 0 && disp_constraint[i][j][k][2] == 0)
                {
                    A_to_here += 2 * info->CP_info[disp_constraint[i][j][k][0] * DIMENSION] + info->CP_info[disp_constraint[i][j][k][0] * DIMENSION + 1];
                }

                if (disp_constraint[i][j][k][1] == 0)
                {
                    for (l = 0; l < info->CP_info[disp_constraint[i][j][k][0] * DIMENSION + 1]; l++)
                    {
                        temp_Boundary[temp] = temp_A[A_to_here + l];
                        temp++;

                        // printf("patch num %d\n", disp_constraint[i][j][k][0]);
                        // printf("xi or eta %d\n", disp_constraint[i][j][k][1]);
                        // printf("start or end %d\n", disp_constraint[i][j][k][2]);
                        // printf("temp = %d\n", temp);
                        // printf("A_to_here + l = %d\n", A_to_here + l);
                        // printf("temp_A[A_to_here + l] = %d\n", temp_A[A_to_here + l]);
                    }
                }
                else if (disp_constraint[i][j][k][1] == 1)
                {
                    for (l = 0; l < info->CP_info[disp_constraint[i][j][k][0] * DIMENSION]; l++)
                    {
                        temp_Boundary[temp] = temp_A[A_to_here + l];
                        temp++;

                        // printf("patch num %d\n", disp_constraint[i][j][k][0]);
                        // printf("xi or eta %d\n", disp_constraint[i][j][k][1]);
                        // printf("start or end %d\n", disp_constraint[i][j][k][2]);
                        // printf("temp = %d\n", temp);
                        // printf("A_to_here + l = %d\n", A_to_here + l);
                        // printf("temp_A[A_to_here + l] = %d\n", temp_A[A_to_here + l]);
                    }
                }
            }
        }
    }

    temp = 0;
    for (i = 0; i < n; i++)
    {
        printf("length_before[%d] = %d\n", i, temp_length_before[i]);
        for (j = 0; j < temp_length_before[i]; j++)
        {
            printf("%d\t", temp_Boundary[temp]);
            temp++;
        }
        printf("\n");
    }

    temp = 0;
    for (i = 0; i < n; i++)
    {
        int *temp_array = (int *)malloc(sizeof(int) * temp_length_before[i]);
        if (temp_array == NULL)
        {
            printf("Memory cannot be allocated\n");
            exit(1);
        }

        for (j = 0; j < temp_length_before[i]; j++)
        {
            temp_array[j] = temp_Boundary[temp + j];
            // printf("%d\t", temp_array[j]);
        }
        printf("\n");

        heapSort(temp_array, temp_length_before[i]);
        Dedupe(temp_array, temp_length_before, temp_Boundary_result, temp_length_after, i);

        free(temp_array);
        temp += temp_length_before[i];
    }

    temp = 0;
    for (i = 0; i < n; i++)
    {
        printf("length_after[%d] = %d\n", i, temp_length_after[i]);
        for (j = 0; j < temp_length_after[i]; j++)
        {
            printf("%d\t", temp_Boundary_result[i * temp + j]);
        }
        printf("\n");
        temp += temp_length_before[i];
    }
    printf("\n");
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