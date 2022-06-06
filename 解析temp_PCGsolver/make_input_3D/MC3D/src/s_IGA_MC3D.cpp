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

    // memory allocation
    info_ptr->disp_constraint_n = (int *)malloc(sizeof(int) * info.DIMENSION);                                          // disp_constraint_n[DIMENSION]
    info_ptr->disp_constraint_face_edge_n = (int *)malloc(sizeof(int) * info.DIMENSION * info.MAX_DISP_CONSTRAINT);     // disp_constraint_face_edge_n[DIMENSION][MAX_DISP_CONSTRAINT]
    info_ptr->disp_constraint_amount = (double *)malloc(sizeof(double) * info.DIMENSION * info.MAX_DISP_CONSTRAINT);    // disp_constraint_amount[DIMENSION][MAX_DISP_CONSTRAINT]
    info_ptr->disp_constraint = (int *)malloc(sizeof(int) * info.DIMENSION * info.MAX_DISP_CONSTRAINT * info.MAX_DISP_CONSTRAINT_FACE_EDGE * 3);     // disp_constraint[DIMENSION][MAX_DISP_CONSTRAINT][MAX_DISP_CONSTRAINT_FACE_EDGE][3]
    info_ptr->distributed_load_info = (double *)malloc(sizeof(double) * info.distributed_load_n * 9);                   // distributed_load_info[MAX_DISTRIBUTED_LOAD][9]
    if (info.disp_constraint_n == NULL || info.disp_constraint_face_edge_n == NULL || info.disp_constraint_amount == NULL || info.disp_constraint == NULL || info.distributed_load_info == NULL)
    {
        printf("Cannot allocate memory\n"); exit(1);
    }

    // ファイル読み込み(2回目)
    Get_inputdata_boundary_1(argv[1], &info); // boundaryのインプットデータ処理

    // memory allocation
    info_ptr->Order = (int *)malloc(sizeof(int) * info.Total_patch * info.DIMENSION);   // int Order[パッチ番号][DIMENSION]
    info_ptr->KV_info = (int *)malloc(sizeof(int) * info.Total_patch * info.DIMENSION); // int KV_info[パッチ番号][DIMENSION]
    info_ptr->CP_info = (int *)malloc(sizeof(int) * info.Total_patch * info.DIMENSION); // int CP_info[パッチ番号][DIMENSION]
    if (info.Order == NULL || info.KV_info == NULL || info.CP_info == NULL)
    {
        printf("Cannot allocate memory\n"); exit(1);
    }

    // patchのインプットデータ処理(1回目)
    for (i = 1; i < Total_input; i++)
    {
        Get_inputdata_patch_0(argv[i + 1], &info);
    }

    int temp1 = 0; // temp1 : 全パッチ含めた総コントロールポイント数
    int temp2 = 0; // temp2 : 4辺 or 6面 のコントロールポイントの和を全パッチ分足した値 * 2
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
            temp3 += info.KV_info[i * info.DIMENSION + j];
        }
    }
    printf("Total Control Point = %d\n", temp1);

    // memory allocation
    info_ptr->CP = (double *)malloc(sizeof(double) * temp1 * (info.DIMENSION + 1));             // CP[CP番号][DIMENSION + 1]
    info_ptr->CP_result = (double *)malloc(sizeof(double) * temp1 * (info.DIMENSION + 1));      // CP_result[通しのコントロールポイント番号(連番)][DIMENSION + 1]
    info_ptr->A = (int *)malloc(sizeof(int) * temp2);                                           // A[パッチ番号][面番号(0~5) or 辺番号(0~3)][辺内のコネクティビティ]
    if (info.DIMENSION == 2)
    {
        info_ptr->B = (double *)malloc(sizeof(double) * info.Total_patch * 16 * (info.DIMENSION + 1)); // B[パッチ番号][辺番号(0~4)][正負方向 2][各辺の端の2頂点][座標xyw  -> 3]
    }
    else if (info.DIMENSION == 3)
    {
        info_ptr->B = (double *)malloc(sizeof(double) * info.Total_patch * 24 * (info.DIMENSION + 1)); // B[パッチ番号][面番号(0~5)][各面の端の4頂点][座標xyzw -> 4]
    }
    info_ptr->Connectivity = (int *)malloc(sizeof(int) * temp1);                                // Connectivity[パッチ番号][パッチ内CP番号]
    info_ptr->KV = (double *)malloc(sizeof(double) * temp3);                                    // KV[パッチ番号][DIMENSION][ノットベクトル番号]

    if (info.CP == NULL || info.CP_result == NULL || info.A == NULL || info.B == NULL || info.Connectivity == NULL || info.KV == NULL)
    {
        printf("Cannot allocate memory\n"); exit(1);
    }

    // patch のインプットデータ処理(2回目)
    counter = 0, KV_to_here = 0, CP_to_here = 0, B_to_here = 0;
    for (i = 1; i < Total_input; i++)
    {
        Get_inputdata_patch_1(argv[i + 1], &info, counter);
        counter++;
    }
    printf("Done get input\n");

    // memory allocation
    if (info.DIMENSION == 2)
    {
        info_ptr->Face_Edge_info = (int *)calloc(info.Total_patch * 32, sizeof(int));       // Face_Edge_info[パッチ番号][own 辺番号(正固定0~3)][opp 辺番号(0~7)]
        info_ptr->Opponent_patch_num = (int *)malloc(sizeof(int) * info.Total_patch * 4);   // Opponent_patch_num[パッチ番号][own 辺番号(正固定0~3]
    }
    else if (info.DIMENSION == 3)
    {
        info_ptr->Face_Edge_info = (int *)malloc(info.Total_patch * 36 * sizeof(int));      // Face_Edge_info[パッチ番号][own 面番号(0~5)][opp 面番号(0~5)]
        info_ptr->Opponent_patch_num = (int *)malloc(sizeof(int) * info.Total_patch * 6);   // Opponent_patch_num[パッチ番号][own 辺番号(正固定0~5]
        for (i = 0; i < info.Total_patch * 36; i++)
        {
            info.Face_Edge_info[i] = -1;
        }
    }

    if (info.Face_Edge_info == NULL || info.Opponent_patch_num == NULL)
    {
        printf("Cannot allocate memory\n"); exit(1);
    }

    // パッチコネクティビティの作成
    printf("state: patch connectivity\n");
    counter = 0, CP_to_here = 0, CP_result_to_here = 0, B_to_here = 0;
    for (i = 0; i < info.Total_patch; i++)
    {
        for (j = 0; j < i; j++)
        {
            if (info.DIMENSION == 2)
            {
                Check_B_2D(i, j, &info);
            }
            else if (info.DIMENSION == 3)
            {
                Check_B_3D(i, j, &info);
            }
        }
        if (info.DIMENSION == 2)
        {
            Make_connectivity_2D(i, &info);
        }
        else if (info.DIMENSION == 3)
        {
            Make_connectivity_3D(i, &info);
        }
    }

    // memory allocation
    int temp4 = 0, temp5 = 0;
    int a = info.MAX_DISP_CONSTRAINT * info.MAX_DISP_CONSTRAINT_FACE_EDGE * 3;
    int b = info.MAX_DISP_CONSTRAINT_FACE_EDGE * 3;
    int c = 3;
    if (info.DIMENSION == 2)
    {
        for (i = 0; i < info.DIMENSION; i++)
        {
            for (j = 0; j < info.disp_constraint_n[i]; j++)
            {
                for (k = 0; k < info.disp_constraint_face_edge_n[i * info.MAX_DISP_CONSTRAINT + j]; k++)
                {
                    int patch_num = info.disp_constraint[i * a + j * b + k * c + 0];
                    if (info.disp_constraint[i * a + j * b + k * c + 1] == 0)
                    {
                        temp4 += info.CP_info[patch_num * info.DIMENSION + 1];
                    }
                    else if (info.disp_constraint[i * a + j * b + k * c + 1] == 1)
                    {
                        temp4 += info.CP_info[patch_num * info.DIMENSION];
                    }
                }
            }
            temp5 += info.disp_constraint_n[i];
        }
    }
    else if (info.DIMENSION == 3)
    {
        for (i = 0; i < info.DIMENSION; i++)
        {
            for (j = 0; j < info.disp_constraint_n[i]; j++)
            {
                for (k = 0; k < info.disp_constraint_face_edge_n[i * info.MAX_DISP_CONSTRAINT + j]; k++)
                {
                    int patch_num = info.disp_constraint[i * a + j * b + k * c + 0];
                    if (info.disp_constraint[i * a + j * b + k * c + 1] == 0)
                    {
                        temp4 += info.CP_info[patch_num * info.DIMENSION + 1] * info.CP_info[patch_num * info.DIMENSION + 2];
                    }
                    else if (info.disp_constraint[i * a + j * b + k * c + 1] == 1)
                    {
                        temp4 += info.CP_info[patch_num * info.DIMENSION] * info.CP_info[patch_num * info.DIMENSION + 2];
                    }
                    else if (info.disp_constraint[i * a + j * b + k * c + 1] == 2)
                    {
                        temp4 += info.CP_info[patch_num * info.DIMENSION] * info.CP_info[patch_num * info.DIMENSION + 1];
                    }
                }
            }
            temp5 += info.disp_constraint_n[i];
        }
    }

    info_ptr->length_before = (int *)calloc(temp5, sizeof(int));    // 各変位量でのマージ前の長さ
    info_ptr->length_after = (int *)malloc(sizeof(int) * temp5);    // 各変位量でのマージ後の長さ
    info_ptr->Boundary = (int *)malloc(sizeof(int) * temp4);        // 境界条件のコネクティビティ
    info_ptr->Boundary_result = (int *)malloc(sizeof(int) * temp4); // ソート・マージ後境界条件のコネクティビティ
    if (info.length_before == NULL || info.length_after == NULL || info.Boundary == NULL || info.Boundary_result == NULL)
    {
        printf("Cannot allocate memory\n"); exit(1);
    }

    // 強制変位・変位固定の境界条件を作成
    Sort(temp5, &info);

    printf("state: output\n");
    Output_inputdata(temp5 ,&info);

    // 図の出力
    if (info.DIMENSION == 2)
    {
        Output_SVG(&info);   // SVG出力
    }

    // memory free
    free(info.disp_constraint_n), free(info.disp_constraint_face_edge_n), free(info.disp_constraint_amount), free(info.disp_constraint), free(info.distributed_load_info);
    free(info.Order), free(info.KV_info), free(info.CP_info);
    free(info.CP), free(info.CP_result), free(info.A), free(info.B), free(info.Connectivity), free(info.KV);
    free(info.Face_Edge_info), free(info.Opponent_patch_num);
    free(info.length_before), free(info.length_after), free(info.Boundary), free(info.Boundary_result);

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
        printf("E_and_nu[%d] = %le\n", i, info->E_and_nu[i]);
    }

    fgets(s, 256, fp);

    // 各方向への変位指定する個数
    int temp_MAX_DISP_CONSTRAINT = 0, temp_MAX_DISP_CONSTRAINT_FACE_EDGE = 0, temp_n_0, temp_n_1;
    for (i = 0; i < info->DIMENSION; i++)
    {
        fscanf(fp, "%d", &temp_i);
        temp_n_0 = temp_i;
        if (temp_MAX_DISP_CONSTRAINT < temp_i)
        {
            temp_MAX_DISP_CONSTRAINT = temp_i;
        }

        for (j = 0; j < temp_n_0; j++)
        {
            fscanf(fp, "%d", &temp_i);
            temp_n_1 = temp_i;
            if (temp_MAX_DISP_CONSTRAINT_FACE_EDGE < temp_i)
            {
                temp_MAX_DISP_CONSTRAINT_FACE_EDGE = temp_i;
            }

            fscanf(fp, "%lf", &temp_d);

            for (k = 0; k < temp_n_1; k++)
            {
                fscanf(fp, "%d", &temp_i);
                fscanf(fp, "%d", &temp_i);
                fscanf(fp, "%d", &temp_i);
            }

            fgets(s, 256, fp);
        }
    }
    info->MAX_DISP_CONSTRAINT = temp_MAX_DISP_CONSTRAINT;
    info->MAX_DISP_CONSTRAINT_FACE_EDGE = temp_MAX_DISP_CONSTRAINT_FACE_EDGE;

    // 分布荷重の荷重の個数
    fscanf(fp, "%d", &temp_i);
    info->distributed_load_n = temp_i;

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
    for (i = 0; i < 2; i++)
    {
        fscanf(fp, "%lf", &temp_d);
    }

    fgets(s, 256, fp);

    // 各方向への変位指定する個数
    int a = info->MAX_DISP_CONSTRAINT * info->MAX_DISP_CONSTRAINT_FACE_EDGE * 3;
    int b = info->MAX_DISP_CONSTRAINT_FACE_EDGE * 3;
    int c = 3;
    for (i = 0; i < info->DIMENSION; i++)
    {
        fscanf(fp, "%d", &temp_i);
        info->disp_constraint_n[i] = temp_i;

        for (j = 0; j < info->disp_constraint_n[i]; j++)
        {
            fscanf(fp, "%d", &temp_i);
            info->disp_constraint_face_edge_n[i * info->MAX_DISP_CONSTRAINT + j] = temp_i;

            fscanf(fp, "%lf", &temp_d);
            info->disp_constraint_amount[i * info->MAX_DISP_CONSTRAINT + j] = temp_d;

            for (k = 0; k < info->disp_constraint_face_edge_n[i * info->MAX_DISP_CONSTRAINT + j]; k++)
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
        Total_CP *= info->CP_info[num * info->DIMENSION + i];
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
        int temp = info->CP_info[num * info->DIMENSION] * info->CP_info[num * info->DIMENSION + 1] * (info->CP_info[num * info->DIMENSION + 2] - 1) * (info->DIMENSION + 1);
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
void Check_B_2D(int num_own, int num_opponent, information *info)
{
    int i, j;
    int ii;

    double x_diff[2], y_diff[2], w_diff[2];

    int Check_B_own_to_here = num_own * 16 * (info->DIMENSION + 1);
    int Check_B_opponent_to_here = num_opponent * 16 * (info->DIMENSION + 1);

    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 8; j++)
        {
            ii = 2 * i;
            // 点0
            x_diff[0] = info->B[Check_B_own_to_here + ii * 2 * (info->DIMENSION + 1)] - info->B[Check_B_opponent_to_here + j * 2 * (info->DIMENSION + 1)];
            y_diff[0] = info->B[Check_B_own_to_here + ii * 2 * (info->DIMENSION + 1) + 1] - info->B[Check_B_opponent_to_here + j * 2 * (info->DIMENSION + 1) + 1];
            w_diff[0] = info->B[Check_B_own_to_here + ii * 2 * (info->DIMENSION + 1) + 2] - info->B[Check_B_opponent_to_here + j * 2 * (info->DIMENSION + 1) + 2];

            // 点1
            x_diff[1] = info->B[Check_B_own_to_here + ii * 2 * (info->DIMENSION + 1) + (info->DIMENSION + 1)] - info->B[Check_B_opponent_to_here + j * 2 * (info->DIMENSION + 1) + (info->DIMENSION + 1)];
            y_diff[1] = info->B[Check_B_own_to_here + ii * 2 * (info->DIMENSION + 1) + (info->DIMENSION + 1) + 1] - info->B[Check_B_opponent_to_here + j * 2 * (info->DIMENSION + 1) + (info->DIMENSION + 1) + 1];
            w_diff[1] = info->B[Check_B_own_to_here + ii * 2 * (info->DIMENSION + 1) + (info->DIMENSION + 1) + 2] - info->B[Check_B_opponent_to_here + j * 2 * (info->DIMENSION + 1) + (info->DIMENSION + 1) + 2];

            // 辺が一致している場合 Face_Edge_info を True
            if (sqrt(pow(x_diff[0], 2) + pow(y_diff[0], 2) + pow(w_diff[0], 2)) <= MERGE_DISTANCE && sqrt(pow(x_diff[1], 2) + pow(y_diff[1], 2) + pow(w_diff[1], 2)) <= MERGE_DISTANCE)
            {
                info->Face_Edge_info[num_own * 32 + i * 8 + j] = 1;
                info->Opponent_patch_num[num_own * 4 + i] = num_opponent;
                printf("own_patch:%d opp_patch:%d own_edge:%d opp_edge:%d\n", num_own, num_opponent, i, j);
                return;
            }
        }
    }
}


void Check_B_3D(int num_own, int num_opponent, information *info)
{
    int i, j, k, ii, jj;

    double face_center_x[2], face_center_y[2], face_center_z[2], face_center_w[2];
    double x_diff, y_diff, z_diff, w_diff;

    int Check_B_own_to_here = num_own * 24 * (info->DIMENSION + 1);
    int Check_B_opponent_to_here = num_opponent * 24 * (info->DIMENSION + 1);

    for (i = 0; i < 6; i++)
    {
        for (j = 0; j < 6; j++)
        {
            ii = i * 4 * (info->DIMENSION + 1);
            jj = j * 4 * (info->DIMENSION + 1);

            // 面own の centerpoint
            face_center_x[0] = (info->B[Check_B_own_to_here + ii + 0] + info->B[Check_B_own_to_here + ii + 4 + 0] + info->B[Check_B_own_to_here + ii + 2 * 4 + 0] + info->B[Check_B_own_to_here + ii + 3 * 4 + 0]) / 4.0;
            face_center_y[0] = (info->B[Check_B_own_to_here + ii + 1] + info->B[Check_B_own_to_here + ii + 4 + 1] + info->B[Check_B_own_to_here + ii + 2 * 4 + 1] + info->B[Check_B_own_to_here + ii + 3 * 4 + 1]) / 4.0;
            face_center_z[0] = (info->B[Check_B_own_to_here + ii + 2] + info->B[Check_B_own_to_here + ii + 4 + 2] + info->B[Check_B_own_to_here + ii + 2 * 4 + 2] + info->B[Check_B_own_to_here + ii + 3 * 4 + 2]) / 4.0;
            face_center_w[0] = (info->B[Check_B_own_to_here + ii + 3] + info->B[Check_B_own_to_here + ii + 4 + 3] + info->B[Check_B_own_to_here + ii + 2 * 4 + 3] + info->B[Check_B_own_to_here + ii + 3 * 4 + 3]) / 4.0;

            // 面opp の centerpoint
            face_center_x[1] = (info->B[Check_B_opponent_to_here + jj + 0] + info->B[Check_B_opponent_to_here + jj + 4 + 0] + info->B[Check_B_opponent_to_here + jj + 2 * 4 + 0] + info->B[Check_B_opponent_to_here + jj + 3 * 4 + 0]) / 4.0;
            face_center_y[1] = (info->B[Check_B_opponent_to_here + jj + 1] + info->B[Check_B_opponent_to_here + jj + 4 + 1] + info->B[Check_B_opponent_to_here + jj + 2 * 4 + 1] + info->B[Check_B_opponent_to_here + jj + 3 * 4 + 1]) / 4.0;
            face_center_z[1] = (info->B[Check_B_opponent_to_here + jj + 2] + info->B[Check_B_opponent_to_here + jj + 4 + 2] + info->B[Check_B_opponent_to_here + jj + 2 * 4 + 2] + info->B[Check_B_opponent_to_here + jj + 3 * 4 + 2]) / 4.0;
            face_center_w[1] = (info->B[Check_B_opponent_to_here + jj + 3] + info->B[Check_B_opponent_to_here + jj + 4 + 3] + info->B[Check_B_opponent_to_here + jj + 2 * 4 + 3] + info->B[Check_B_opponent_to_here + jj + 3 * 4 + 3]) / 4.0;
            
            // centerpoint の diff
            x_diff = face_center_x[1] - face_center_x[0];
            y_diff = face_center_y[1] - face_center_y[0];
            z_diff = face_center_z[1] - face_center_z[0];
            w_diff = face_center_w[1] - face_center_w[0];

            // 面の中心点が一致している場合 Face_Edge_info を Mode 番号 (k) に
            if (sqrt(pow(x_diff, 2) + pow(y_diff, 2) + pow(z_diff, 2) + pow(w_diff, 2)) <= MERGE_DISTANCE)
            {
                info->Opponent_patch_num[num_own * 6 + i] = num_opponent;
                for (k = 0; k < 4; k++)
                {
                    if (k == 0)
                    {
                        x_diff = info->B[Check_B_own_to_here + ii] - info->B[Check_B_opponent_to_here + jj];
                        y_diff = info->B[Check_B_own_to_here + ii + 1] - info->B[Check_B_opponent_to_here + jj + 1];
                        z_diff = info->B[Check_B_own_to_here + ii + 2] - info->B[Check_B_opponent_to_here + jj + 2];
                        w_diff = info->B[Check_B_own_to_here + ii + 3] - info->B[Check_B_opponent_to_here + jj + 3];
                    }
                    else if (k == 1)
                    {
                        x_diff = info->B[Check_B_own_to_here + ii] - info->B[Check_B_opponent_to_here + jj + (info->DIMENSION + 1)];
                        y_diff = info->B[Check_B_own_to_here + ii + 1] - info->B[Check_B_opponent_to_here + jj + (info->DIMENSION + 1) + 1];
                        z_diff = info->B[Check_B_own_to_here + ii + 2] - info->B[Check_B_opponent_to_here + jj + (info->DIMENSION + 1) + 2];
                        w_diff = info->B[Check_B_own_to_here + ii + 3] - info->B[Check_B_opponent_to_here + jj + (info->DIMENSION + 1) + 3];
                    }
                    else if (k == 2)
                    {
                        x_diff = info->B[Check_B_own_to_here + ii] - info->B[Check_B_opponent_to_here + jj + 2 * (info->DIMENSION + 1)];
                        y_diff = info->B[Check_B_own_to_here + ii + 1] - info->B[Check_B_opponent_to_here + jj + 2 * (info->DIMENSION + 1) + 1];
                        z_diff = info->B[Check_B_own_to_here + ii + 2] - info->B[Check_B_opponent_to_here + jj + 2 * (info->DIMENSION + 1) + 2];
                        w_diff = info->B[Check_B_own_to_here + ii + 3] - info->B[Check_B_opponent_to_here + jj + 2 * (info->DIMENSION + 1) + 3];
                    }
                    else if (k == 3)
                    {
                        x_diff = info->B[Check_B_own_to_here + ii] - info->B[Check_B_opponent_to_here + jj + 3 * (info->DIMENSION + 1)];
                        y_diff = info->B[Check_B_own_to_here + ii + 1] - info->B[Check_B_opponent_to_here + jj + 3 * (info->DIMENSION + 1) + 1];
                        z_diff = info->B[Check_B_own_to_here + ii + 2] - info->B[Check_B_opponent_to_here + jj + 3 * (info->DIMENSION + 1) + 2];
                        w_diff = info->B[Check_B_own_to_here + ii + 3] - info->B[Check_B_opponent_to_here + jj + 3 * (info->DIMENSION + 1) + 3];
                    }
                    if (sqrt(pow(x_diff, 2) + pow(y_diff, 2) + pow(z_diff, 2) + pow(w_diff, 2)) <= MERGE_DISTANCE)
                    {
                        info->Face_Edge_info[num_own * 36 + i * 6 + j] = k;
                        printf("own_patch:%d opp_patch:%d own_face:%d opp_face:%d mode:%d\n", num_own, num_opponent, i, j, k);
                        return;
                    }
                }
            }
        }
    }
}


void Make_connectivity_2D(int num, information *info)
{
    int i, j, k;
    int p, q;
    int temp_CP_n = 0;
    int Edge[4];
    int Edge_to_here = num * 32;

    printf("make connectivity on patch %d\n", num);

    int A_to_own = 0;
    for (i = 0; i < num; i++)
    {
        A_to_own += 2 * (info->CP_info[i * info->DIMENSION] + info->CP_info[i * info->DIMENSION + 1]);
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
            A_to_own += info->CP_info[num * info->DIMENSION];
        }
        else if (i == 2)
        {
            A_to_own += info->CP_info[num * info->DIMENSION + 1];
        }
        else if (i == 3)
        {
            A_to_own += info->CP_info[num * info->DIMENSION];
        }

        for (j = 0; j < 8; j++)
        {
            if (info->Face_Edge_info[Edge_to_here + i * 8 + j] == 1)
            {
                int A_to_opponent = 0;
                for (k = 0; k < info->Opponent_patch_num[num * 4 + i]; k++)
                {
                    A_to_opponent += 2 * (info->CP_info[k * info->DIMENSION] + info->CP_info[k * info->DIMENSION + 1]);
                }

                printf("Patch num = %d\n", num);
                printf("Edge num = %d\n", j);
                printf("Opponent patch num = %d\n", info->Opponent_patch_num[num * 4 + i]);

                Edge[i] = 1;

                p = j / 2;
                q = j % 2;
                printf("p = %d\n", p);
                printf("q = %d\n", q);

                if (p == 0)
                {
                    temp_CP_n = info->CP_info[info->Opponent_patch_num[num * 4 + i] * info->DIMENSION];
                }
                else if (p == 1)
                {
                    temp_CP_n = info->CP_info[info->Opponent_patch_num[num * 4 + i] * info->DIMENSION + 1];
                    A_to_opponent += info->CP_info[info->Opponent_patch_num[num * 4 + i] * info->DIMENSION];
                }
                else if (p == 2)
                {
                    temp_CP_n = info->CP_info[info->Opponent_patch_num[num * 4 + i] * info->DIMENSION];
                    A_to_opponent += info->CP_info[info->Opponent_patch_num[num * 4 + i] * info->DIMENSION] + info->CP_info[info->Opponent_patch_num[num * 4 + i] * info->DIMENSION + 1];
                }
                else if (p == 3)
                {
                    temp_CP_n = info->CP_info[info->Opponent_patch_num[num * 4 + i] * info->DIMENSION + 1];
                    A_to_opponent += 2 * info->CP_info[info->Opponent_patch_num[num * 4 + i] * info->DIMENSION] + info->CP_info[info->Opponent_patch_num[num * 4 + i] * info->DIMENSION + 1];
                }

                if (q == 0)
                {
                    for (k = 0; k < temp_CP_n; k++)
                    {
                        info->A[A_to_own + k] = info->A[A_to_opponent + k];
                    }
                    break;
                }
                else if (q == 1)
                {
                    for (k = 0; k < temp_CP_n; k++)
                    {
                        info->A[A_to_own + k] = info->A[A_to_opponent + (temp_CP_n - 1) - k];
                    }
                    break;
                }
            }
        }
    }

    // コネクティビティを作成
    int xi, eta;

    A_to_own = 0;
    for (i = 0; i < num; i++)
    {
        A_to_own += 2 * (info->CP_info[i * info->DIMENSION] + info->CP_info[i * info->DIMENSION + 1]);
    }

    for (eta = 0; eta < info->CP_info[num * info->DIMENSION + 1]; eta++)
    {
        for (xi = 0; xi < info->CP_info[num * info->DIMENSION]; xi++)
        {
            if (eta == 0 && Edge[0] == 1)
            {
                info->Connectivity[CP_to_here + eta * info->CP_info[num * info->DIMENSION] + xi] = info->A[A_to_own + xi];
            }
            else if (eta == info->CP_info[num * info->DIMENSION + 1] - 1 && Edge[2] == 1)
            {
                info->Connectivity[CP_to_here + eta * info->CP_info[num * info->DIMENSION] + xi] = info->A[A_to_own + info->CP_info[num * info->DIMENSION] + info->CP_info[num * info->DIMENSION + 1] + xi];
            }
            else if (xi == 0 && Edge[3] == 1)
            {
                info->Connectivity[CP_to_here + eta * info->CP_info[num * info->DIMENSION] + xi] = info->A[A_to_own + 2 * info->CP_info[num * info->DIMENSION] + info->CP_info[num * info->DIMENSION + 1] + eta];
            }
            else if (xi == info->CP_info[num * info->DIMENSION] - 1 && Edge[1] == 1)
            {
                info->Connectivity[CP_to_here + eta * info->CP_info[num * info->DIMENSION] + xi] = info->A[A_to_own + info->CP_info[num * info->DIMENSION] + eta];
            }
            else
            {
                info->Connectivity[CP_to_here + eta * info->CP_info[num * info->DIMENSION] + xi] = counter;
                info->CP_result[CP_result_to_here] = info->CP[(CP_to_here + eta * info->CP_info[num * info->DIMENSION] + xi) * (info->DIMENSION + 1)];
                info->CP_result[CP_result_to_here + 1] = info->CP[(CP_to_here + eta * info->CP_info[num * info->DIMENSION] + xi) * (info->DIMENSION + 1) + 1];
                info->CP_result[CP_result_to_here + 2] = info->CP[(CP_to_here + eta * info->CP_info[num * info->DIMENSION] + xi) * (info->DIMENSION + 1) + 2];
                counter++;
                CP_result_to_here += (info->DIMENSION + 1);
            }
        }
    }

    // A 配列の作ってない分を作成
    for (eta = 0; eta < info->CP_info[num * info->DIMENSION + 1]; eta++)
    {
        for (xi = 0; xi < info->CP_info[num * info->DIMENSION]; xi++)
        {
            if (eta == 0 && Edge[0] == 0)
            {
                info->A[A_to_own + xi] = info->Connectivity[CP_to_here + eta * info->CP_info[num * info->DIMENSION] + xi];
            }
            if (eta == info->CP_info[num * info->DIMENSION + 1] - 1 && Edge[2] == 0)
            {
                info->A[A_to_own + info->CP_info[num * info->DIMENSION] + info->CP_info[num * info->DIMENSION + 1] + xi] = info->Connectivity[CP_to_here + eta * info->CP_info[num * info->DIMENSION] + xi];
            }
            if (xi == 0 && Edge[3] == 0)
            {
                info->A[A_to_own + 2 * info->CP_info[num * info->DIMENSION] + info->CP_info[num * info->DIMENSION + 1] + eta] = info->Connectivity[CP_to_here + eta * info->CP_info[num * info->DIMENSION] + xi];
            }
            if (xi == info->CP_info[num * info->DIMENSION] - 1 && Edge[1] == 0)
            {
                info->A[A_to_own + info->CP_info[num * info->DIMENSION] + eta] = info->Connectivity[CP_to_here + eta * info->CP_info[num * info->DIMENSION] + xi];
            }
        }
    }
    CP_to_here += info->CP_info[num * info->DIMENSION] * info->CP_info[num * info->DIMENSION + 1];

    printf("patch %d array A after make connectivity\n", num);
    A_to_own = 0;
    for (i = 0; i < num; i++)
    {
        A_to_own += 2 * (info->CP_info[i * info->DIMENSION] + info->CP_info[i * info->DIMENSION + 1]);
    }
    for (i = 0; i < 2 * (info->CP_info[num * info->DIMENSION] + info->CP_info[num * info->DIMENSION + 1]); i++)
    {
        printf("%d\t", info->A[A_to_own + i]);
    }
    printf("\n");

    printf("\n");
}


void Make_connectivity_3D(int num, information *info)
{
    int i, j, k, l;
    int Face[6] = {0};
    int Face_to_here = num * 36;
    int own_CP_a = 0, own_CP_b = 0;

    printf("make connectivity on patch %d\n", num);

    int A_to_own = 0;
    for (i = 0; i < num; i++)
    {
        A_to_own += 2 * (info->CP_info[i * info->DIMENSION] * info->CP_info[i * info->DIMENSION + 1] + info->CP_info[i * info->DIMENSION + 1] * info->CP_info[i * info->DIMENSION + 2] + info->CP_info[i * info->DIMENSION] * info->CP_info[i * info->DIMENSION + 2]);
    }

    // 重なっている面の A 配列を作成
    for (i = 0; i < 6; i++)
    {
        if (i == 0)
        {
            own_CP_a = info->CP_info[num * info->DIMENSION];
            own_CP_b = info->CP_info[num * info->DIMENSION + 1];
        }
        else if (i == 1)
        {
            A_to_own += info->CP_info[num * info->DIMENSION] * info->CP_info[num * info->DIMENSION + 1];
            own_CP_a = info->CP_info[num * info->DIMENSION];
            own_CP_b = info->CP_info[num * info->DIMENSION + 1];
        }
        else if (i == 2)
        {
            A_to_own += info->CP_info[num * info->DIMENSION] * info->CP_info[num * info->DIMENSION + 1];
            own_CP_a = info->CP_info[num * info->DIMENSION + 1];
            own_CP_b = info->CP_info[num * info->DIMENSION + 2];
        }
        else if (i == 3)
        {
            A_to_own += info->CP_info[num * info->DIMENSION + 1] * info->CP_info[num * info->DIMENSION + 2];
            own_CP_a = info->CP_info[num * info->DIMENSION];
            own_CP_b = info->CP_info[num * info->DIMENSION + 2];
        }
        else if (i == 4)
        {
            A_to_own += info->CP_info[num * info->DIMENSION] * info->CP_info[num * info->DIMENSION + 2];
            own_CP_a = info->CP_info[num * info->DIMENSION + 1];
            own_CP_b = info->CP_info[num * info->DIMENSION + 2];
        }
        else if (i == 5)
        {
            A_to_own += info->CP_info[num * info->DIMENSION + 1] * info->CP_info[num * info->DIMENSION + 2];
            own_CP_a = info->CP_info[num * info->DIMENSION];
            own_CP_b = info->CP_info[num * info->DIMENSION + 2];
        }

        for (j = 0; j < 6; j++)
        {
            if (info->Face_Edge_info[Face_to_here + i * 6 + j] >= 0)
            {
                int opp_num = info->Opponent_patch_num[num * 6 + i];
                int A_to_opponent = 0;
                for (k = 0; k < opp_num; k++)
                {
                    A_to_opponent += 2 * (info->CP_info[k * info->DIMENSION] * info->CP_info[k * info->DIMENSION + 1] + info->CP_info[k * info->DIMENSION + 1] * info->CP_info[k * info->DIMENSION + 2] + info->CP_info[k * info->DIMENSION] * info->CP_info[k * info->DIMENSION + 2]);
                }

                printf("Patch num = %d\n", num);
                printf("Face num = %d\n", j);
                printf("Opponent patch num = %d\n", opp_num);

                Face[i] = 1;

                for (k = 1; k <= j; k++)
                {
                    if (k == 1)
                    {
                        A_to_opponent += info->CP_info[opp_num * info->DIMENSION] * info->CP_info[opp_num * info->DIMENSION + 1];
                    }
                    else if (k == 2)
                    {
                        A_to_opponent += info->CP_info[opp_num * info->DIMENSION] * info->CP_info[opp_num * info->DIMENSION + 1];
                    }
                    else if (k == 3)
                    {
                        A_to_opponent += info->CP_info[opp_num * info->DIMENSION + 1] * info->CP_info[opp_num * info->DIMENSION + 2];
                    }
                    else if (k == 4)
                    {
                        A_to_opponent += info->CP_info[opp_num * info->DIMENSION] * info->CP_info[opp_num * info->DIMENSION + 2];
                    }
                    else if (k == 5)
                    {
                        A_to_opponent += info->CP_info[opp_num * info->DIMENSION + 1] * info->CP_info[opp_num * info->DIMENSION + 2];
                    }
                }
                
                if (info->Face_Edge_info[Face_to_here + i * 6 + j] == 0)
                {
                    for (k = 0; k < own_CP_b; k++)
                    {
                        for (l = 0; l < own_CP_a; l++)
                        {
                            info->A[A_to_own + k * own_CP_a + l] = info->A[A_to_opponent + k * own_CP_a + l];
                        }
                    }
                    break;
                }
                else if (info->Face_Edge_info[Face_to_here + i * 6 + j]  == 1)
                {
                    for (k = 0; k < own_CP_b; k++)
                    {
                        for (l = 0; l < own_CP_a; l++)
                        {
                            info->A[A_to_own + k * own_CP_a + l] = info->A[A_to_opponent + k * own_CP_a + ((own_CP_a - 1) - l)];
                        }
                    }
                    break;
                }
                else if (info->Face_Edge_info[Face_to_here + i * 6 + j]  == 2)
                {
                    for (k = 0; k < own_CP_b; k++)
                    {
                        for (l = 0; l < own_CP_a; l++)
                        {
                            info->A[A_to_own + k * own_CP_a + l] = info->A[A_to_opponent + ((own_CP_b - 1) - k) * own_CP_a + ((own_CP_a - 1) - l)];
                        }
                    }
                    break;
                }
                else if (info->Face_Edge_info[Face_to_here + i * 6 + j]  == 3)
                {
                    for (k = 0; k < own_CP_b; k++)
                    {
                        for (l = 0; l < own_CP_a; l++)
                        {
                            info->A[A_to_own + k * own_CP_a + l] = info->A[A_to_opponent + ((own_CP_b - 1) - k) * own_CP_a + l];
                        }
                    }
                    break;
                }
            }
        }
    }

    // コネクティビティを作成
    int xi, eta, zeta;

    A_to_own = 0;
    for (i = 0; i < num; i++)
    {
        A_to_own += 2 * (info->CP_info[i * info->DIMENSION] * info->CP_info[i * info->DIMENSION + 1] + info->CP_info[i * info->DIMENSION + 1] * info->CP_info[i * info->DIMENSION + 2] + info->CP_info[i * info->DIMENSION] * info->CP_info[i * info->DIMENSION + 2]);
    }

    int temp1, temp2, temp3, temp4, temp5;
    temp1 =         info->CP_info[num * info->DIMENSION] * info->CP_info[num * info->DIMENSION + 1];
    temp2 = temp1 + info->CP_info[num * info->DIMENSION] * info->CP_info[num * info->DIMENSION + 1];
    temp3 = temp2 + info->CP_info[num * info->DIMENSION + 1] * info->CP_info[num * info->DIMENSION + 2];
    temp4 = temp3 + info->CP_info[num * info->DIMENSION] * info->CP_info[num * info->DIMENSION + 2];
    temp5 = temp4 + info->CP_info[num * info->DIMENSION + 1] * info->CP_info[num * info->DIMENSION + 2];
    int a, b;
    a = info->CP_info[num * info->DIMENSION] * info->CP_info[num * info->DIMENSION + 1];
    b = info->CP_info[num * info->DIMENSION];

    for (zeta = 0; zeta < info->CP_info[num * info->DIMENSION + 2]; zeta++)
    {
        for (eta = 0; eta < info->CP_info[num * info->DIMENSION + 1]; eta++)
        {
            for (xi = 0; xi < info->CP_info[num * info->DIMENSION]; xi++)
            {
                if (zeta == 0 && Face[0] == 1)
                {
                    info->Connectivity[CP_to_here + zeta * a + eta * b + xi] = info->A[A_to_own + eta * info->CP_info[num * info->DIMENSION] + xi];
                }
                else if (zeta == info->CP_info[num * info->DIMENSION + 2] - 1 && Face[1] == 1)
                {
                    info->Connectivity[CP_to_here + zeta * a + eta * b + xi] = info->A[A_to_own + temp1 + eta * info->CP_info[num * info->DIMENSION] + xi];
                }
                else if (eta == 0 && Face[3] == 1)
                {
                    info->Connectivity[CP_to_here + zeta * a + eta * b + xi] = info->A[A_to_own + temp3 + zeta * info->CP_info[num * info->DIMENSION] + xi];
                }
                else if (eta == info->CP_info[num * info->DIMENSION + 1] - 1 && Face[5] == 1)
                {
                    info->Connectivity[CP_to_here + zeta * a + eta * b + xi] = info->A[A_to_own + temp5 + zeta * info->CP_info[num * info->DIMENSION] + xi];
                }
                else if (xi == 0 && Face[2] == 1)
                {
                    info->Connectivity[CP_to_here + zeta * a + eta * b + xi] = info->A[A_to_own + temp2 + zeta * info->CP_info[num * info->DIMENSION + 1] + eta];
                }
                else if (xi == info->CP_info[num * info->DIMENSION] - 1 && Face[4] == 1)
                {
                    info->Connectivity[CP_to_here + zeta * a + eta * b + xi] = info->A[A_to_own + temp4 + zeta * info->CP_info[num * info->DIMENSION + 1] + eta];
                }
                else
                {
                    info->Connectivity[CP_to_here + zeta * a + eta * b + xi] = counter;
                    info->CP_result[CP_result_to_here] = info->CP[(CP_to_here + zeta * a + eta * b + xi) * (info->DIMENSION + 1)];
                    info->CP_result[CP_result_to_here + 1] = info->CP[(CP_to_here + zeta * a + eta * b + xi) * (info->DIMENSION + 1) + 1];
                    info->CP_result[CP_result_to_here + 2] = info->CP[(CP_to_here + zeta * a + eta * b + xi) * (info->DIMENSION + 1) + 2];
                    info->CP_result[CP_result_to_here + 3] = info->CP[(CP_to_here + zeta * a + eta * b + xi) * (info->DIMENSION + 1) + 3];
                    counter++;
                    CP_result_to_here += (info->DIMENSION + 1);
                }
            }
        }
    }

    A_to_own = 0;
    for (i = 0; i < num; i++)
    {
        A_to_own += 2 * (info->CP_info[i * info->DIMENSION] * info->CP_info[i * info->DIMENSION + 1] + info->CP_info[i * info->DIMENSION + 1] * info->CP_info[i * info->DIMENSION + 2] + info->CP_info[i * info->DIMENSION] * info->CP_info[i * info->DIMENSION + 2]);
    }

    // A 配列の作ってない分を作成
    for (zeta = 0; zeta < info->CP_info[num * info->DIMENSION + 2]; zeta++)
    {
        for (eta = 0; eta < info->CP_info[num * info->DIMENSION + 1]; eta++)
        {
            for (xi = 0; xi < info->CP_info[num * info->DIMENSION]; xi++)
            {
                if (zeta == 0 && Face[0] == 0)
                {
                    info->A[A_to_own + eta * info->CP_info[num * info->DIMENSION] + xi] = info->Connectivity[CP_to_here + zeta * a + eta * b + xi];
                }
                if (zeta == info->CP_info[num * info->DIMENSION + 2] - 1 && Face[1] == 0)
                {
                    info->A[A_to_own + temp1 + eta * info->CP_info[num * info->DIMENSION] + xi] = info->Connectivity[CP_to_here + zeta * a + eta * b + xi];
                }
                if (eta == 0 && Face[3] == 0)
                {
                    info->A[A_to_own + temp3 + zeta * info->CP_info[num * info->DIMENSION] + xi] = info->Connectivity[CP_to_here + zeta * a + eta * b + xi];
                }
                if (eta == info->CP_info[num * info->DIMENSION + 1] - 1 && Face[5] == 0)
                {
                    info->A[A_to_own + temp5 + zeta * info->CP_info[num * info->DIMENSION] + xi] = info->Connectivity[CP_to_here + zeta * a + eta * b + xi];
                }
                if (xi == 0 && Face[2] == 0)
                {
                    info->A[A_to_own + temp2 + zeta * info->CP_info[num * info->DIMENSION + 1] + eta] = info->Connectivity[CP_to_here + zeta * a + eta * b + xi];
                }
                if (xi == info->CP_info[num * info->DIMENSION] - 1 && Face[4] == 0)
                {
                    info->A[A_to_own + temp4 + zeta * info->CP_info[num * info->DIMENSION + 1] + eta] = info->Connectivity[CP_to_here + zeta * a + eta * b + xi];
                }
            }
        }
    }
    CP_to_here += info->CP_info[num * info->DIMENSION] * info->CP_info[num * info->DIMENSION + 1] * info->CP_info[num * info->DIMENSION + 2];
}


// output
void Output_inputdata(int total_disp_constraint_n, const information *info)
{
    int i, j, k;
    char str[256] = "input.txt";

    fp = fopen(str, "w");

    // ヤング率
    fprintf(fp, "%d  ", (int)info->E_and_nu[0]);

    // ポアソン比
    fprintf(fp, "%le\n\n", info->E_and_nu[1]);

    // パッチ数
    fprintf(fp, "%d\n\n", info->Total_patch);

    // コントロールポイント数
    fprintf(fp, "%d\n\n", (CP_result_to_here + 1) / (info->DIMENSION + 1));
    int temp_num = (CP_result_to_here + 1) / (info->DIMENSION + 1), temp_counter = 0;
    while (temp_num != 0)
    {
        temp_num = temp_num / 10;
        temp_counter++;
    }
    temp_num = -(temp_counter + 2);

    // 各パッチ内での各方向の次数
    for (i = 0; i < info->Total_patch; i++)
    {
        for (j = 0; j < info->DIMENSION; j++)
        {
            if (j == info->DIMENSION - 1)
            {
                fprintf(fp, "%d", info->Order[i * info->DIMENSION + j]);
            }
            else
            {
                fprintf(fp, "%*d", -6, info->Order[i * info->DIMENSION + j]);
            }
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // 各パッチ内での各方向のノットベクトルの数
    for (i = 0; i < info->Total_patch; i++)
    {
        for (j = 0; j < info->DIMENSION; j++)
        {
            if (j == info->DIMENSION - 1)
            {
                fprintf(fp, "%d", info->KV_info[i * info->DIMENSION + j]);
            }
            else
            {
                fprintf(fp, "%*d", -6, info->KV_info[i * info->DIMENSION + j]);
            }
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // 各パッチ内での各方向のコントロールポイントの数
    for (i = 0; i < info->Total_patch; i++)
    {
        for (j = 0; j < info->DIMENSION; j++)
        {
            if (j == info->DIMENSION - 1)
            {
                fprintf(fp, "%d", info->CP_info[i * info->DIMENSION + j]);
            }
            else
            {
                fprintf(fp, "%*d", -6, info->CP_info[i * info->DIMENSION + j]);
            }
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // パッチコネクティビティ
    CP_to_here = 0;
    for (i = 0; i < info->Total_patch; i++)
    {
        int Total_CP = 1;
        for (j = 0; j < info->DIMENSION; j++)
        {
            Total_CP *= info->CP_info[i * info->DIMENSION + j];
        }
        for (j = 0; j < Total_CP; j++)
        {
            if (j == Total_CP - 1)
            {
                fprintf(fp, "%d", info->Connectivity[CP_to_here + j]);
            }
            else
            {
                fprintf(fp, "%*d", temp_num, info->Connectivity[CP_to_here + j]);
            }
        }
        CP_to_here += Total_CP;
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // 変位拘束するコントロールポイントの数
    int temp = 0;
    for (i = 0; i < total_disp_constraint_n; i++)
    {
        temp += info->length_after[i];
    }
    fprintf(fp, "%*d", -6, temp);

    // 荷重条件を与えるコントロールポイントの数
    fprintf(fp, "%*d", -6, 0);

    // 分布荷重の数
    fprintf(fp, "%d\n\n", info->distributed_load_n);

    // 各パッチでの各方向のノットベクトル
    KV_to_here = 0;
    for (i = 0; i < info->Total_patch; i++)
    {
        for (j = 0; j < info->DIMENSION; j++)
        {
            for (k = 0; k < info->KV_info[i * info->DIMENSION + j]; k++)
            {
                if (k == info->KV_info[i * info->DIMENSION + j] - 1)
                {
                    fprintf(fp, "%.16e", info->KV[KV_to_here + k]);
                }
                else
                {
                    fprintf(fp, "%.16e  ", info->KV[KV_to_here + k]);
                }
            }
            KV_to_here += info->KV_info[i * info->DIMENSION + j];
            fprintf(fp, "\n");
        }
    }
    fprintf(fp, "\n");

    // コントロールポイント
    if (info->DIMENSION == 2)
    {
        for (i = 0; i < (CP_result_to_here + 1) / (info->DIMENSION + 1); i++)
        {
            fprintf(fp, "%*d", temp_num, i);
            fprintf(fp, "% .16e ", info->CP_result[i * (info->DIMENSION + 1)]);
            fprintf(fp, "% .16e ", info->CP_result[i * (info->DIMENSION + 1) + 1]);
            fprintf(fp, "% .16e\n", info->CP_result[i * (info->DIMENSION + 1) + 2]);
        }
    }
    else if (info->DIMENSION == 3)
    {
        for (i = 0; i < (CP_result_to_here + 1) / (info->DIMENSION + 1); i++)
        {
            fprintf(fp, "%*d", temp_num, i);
            fprintf(fp, "% .16e ", info->CP_result[i * (info->DIMENSION + 1)]);
            fprintf(fp, "% .16e ", info->CP_result[i * (info->DIMENSION + 1) + 1]);
            fprintf(fp, "% .16e ", info->CP_result[i * (info->DIMENSION + 1) + 2]);
            fprintf(fp, "% .16e\n", info->CP_result[i * (info->DIMENSION + 1) + 3]);
        }
    }
    fprintf(fp, "\n");

    // 拘束するコントロールポイント
    temp_counter = 0;
    temp = 0;
    for (i = 0; i < info->DIMENSION; i++)
    {
        for (j = 0; j < info->disp_constraint_n[i]; j++)
        {
            for (k = 0; k < info->length_after[temp_counter]; k++)
            {
                fprintf(fp, "%*d", temp_num, info->Boundary_result[temp + k]);
                fprintf(fp, "%*d", temp_num, i);
                fprintf(fp, "%le\n", info->disp_constraint_amount[i * info->MAX_DISP_CONSTRAINT + j]);
            }
            temp += info->length_before[j];
            temp_counter++;
        }
    }
    fprintf(fp, "\n");

    // 分布荷重
    for (i = 0; i < info->distributed_load_n; i++)
    {
        if (i != 0)
        {
            fprintf(fp, "\n");
        }
        fprintf(fp, "%*d", temp_num, (int)info->distributed_load_info[i * 9 + 0]);
        fprintf(fp, "%*d", temp_num, (int)info->distributed_load_info[i * 9 + 1]);
        fprintf(fp, "%*d", temp_num, (int)info->distributed_load_info[i * 9 + 2]);
        fprintf(fp, "%le  ", info->distributed_load_info[i * 9 + 3]);
        fprintf(fp, "%le  ", info->distributed_load_info[i * 9 + 4]);
        fprintf(fp, "%le  ", info->distributed_load_info[i * 9 + 5]);
        fprintf(fp, "%le  ", info->distributed_load_info[i * 9 + 6]);
        fprintf(fp, "%le  ", info->distributed_load_info[i * 9 + 7]);
        fprintf(fp, "%le", info->distributed_load_info[i * 9 + 8]);
    }

    fclose(fp);
}


void Output_SVG(const information *info)
{
    int i;

    // 点サイズ, デフォルト: 2pt
    int point_size = 2;
    // 文字サイズ, 5 ~ 20 程度が適切, デフォルト: 6pt
    int font_size = 6;
    // 画像サイズ, 500 ~ 3000 程度が適切 scale 大 ⇒ 文字 小
    double size = 1000.0;

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

    for (i = 0; i < (CP_result_to_here + 1) / (info->DIMENSION + 1); i++)
    {
        if (i == 0)
        {
            x_min = info->CP_result[i * 3];
            x_max = info->CP_result[i * 3];
            y_min = info->CP_result[i * 3 + 1];
            y_max = info->CP_result[i * 3 + 1];
        }
        else
        {
            if (x_min > info->CP_result[i * 3])
            {
                x_min = info->CP_result[i * 3];
            }
            else if (x_max < info->CP_result[i * 3])
            {
                x_max = info->CP_result[i * 3];
            }

            if (y_min > info->CP_result[i * 3 + 1])
            {
                y_min = info->CP_result[i * 3 + 1];
            }
            else if (y_max < info->CP_result[i * 3 + 1])
            {
                y_max = info->CP_result[i * 3 + 1];
            }
        }
    }

    double x_gap = - x_min;
    double y_gap = - y_min;

    double origin_width = (x_max - x_min) * (22.0 / 20.0);
    double origin_height = (y_max - y_min) * (22.0 / 20.0);

    // 横幅固定，アスペクト比維持
    double width = size, height = size * (origin_height / origin_width);

    double x_scale = width / origin_width;
    double y_scale = x_scale * (origin_height / origin_width);

    double origin_space_x = origin_width * (1.0 / 20.0);
    double origin_space_y = origin_height * (1.0 / 20.0);

    char str[256] = "input.svg";

    fp = fopen(str, "w");

    fprintf(fp, "<?xml version='1.0'?>\n");
    fprintf(fp, "<svg width='%le' height='%le' version='1.1' style='background: #eee' xmlns='http://www.w3.org/2000/svg' xmlns:xlink='http://www.w3.org/1999/xlink'>\n", width, height);

    // パッチ境界を描画
    int temp_color_num = 0;
    B_to_here = 0;
    for (i = 0; i < info->Total_patch; i++)
    {
        position_x = (x_gap + info->B[B_to_here] + origin_space_x) * x_scale;
        position_y = height - ((y_gap + info->B[B_to_here + 1] + origin_space_y) * y_scale);
        fprintf(fp, "<path d='M %le %le ", position_x, position_y);
        B_to_here += 4 * (info->DIMENSION + 1);

        position_x = (x_gap + info->B[B_to_here] + origin_space_x) * x_scale;
        position_y = height - ((y_gap + info->B[B_to_here + 1] + origin_space_y) * y_scale);
        fprintf(fp, "L %le %le ", position_x, position_y);
        B_to_here += 2 * (info->DIMENSION + 1);

        position_x = (x_gap + info->B[B_to_here] + origin_space_x) * x_scale;
        position_y = height - ((y_gap + info->B[B_to_here + 1] + origin_space_y) * y_scale);
        fprintf(fp, "L %le %le ", position_x, position_y);
        B_to_here += 2 * (info->DIMENSION + 1);

        position_x = (x_gap + info->B[B_to_here] + origin_space_x) * x_scale;
        position_y = height - ((y_gap + info->B[B_to_here + 1] + origin_space_y) * y_scale);
        fprintf(fp, "L %le %le ", position_x, position_y);
        B_to_here += 4 * (info->DIMENSION + 1);

        fprintf(fp, "Z' fill='%s'/>\n", color_vec[temp_color_num % 10]);
        B_to_here += 4 * (info->DIMENSION + 1);

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
    for (i = 0; i < (CP_result_to_here + 1) / (info->DIMENSION + 1); i++)
    {
        position_x = (x_gap + info->CP_result[i * 3] + origin_space_x) * x_scale;
        position_y = height - ((y_gap + info->CP_result[i * 3 + 1] + origin_space_y) * y_scale);
        fprintf(fp, "<circle cx='%le' cy='%le' r='%d' fill='%s'/>\n", position_x, position_y, point_size, color_vec[7]);
        fprintf(fp, "<text x='%le' y='%le' font-family='Verdana' font-size='%d' fill='%s' font-weight='700'>\n", position_x + 2, position_y, font_size, color_vec[7]);
        fprintf(fp, "%d\n", i);
        fprintf(fp, "</text>\n");
    }

    fprintf(fp, "</svg>");
    fclose(fp);
}


// heap sort
void Sort(int n, information *info)
{
    int i, j, k, l;
    int a = info->MAX_DISP_CONSTRAINT * info->MAX_DISP_CONSTRAINT_FACE_EDGE * 3;
    int b = info->MAX_DISP_CONSTRAINT_FACE_EDGE * 3;
    int c = 3;
    int temp = 0;

    if (info->DIMENSION == 2)
    {
        for (i = 0; i < info->DIMENSION; i++)
        {
            for (j = 0; j < info->disp_constraint_n[i]; j++)
            {
                for (k = 0; k < info->disp_constraint_face_edge_n[i * info->MAX_DISP_CONSTRAINT + j]; k++)
                {
                    int patch_num = info->disp_constraint[i * a + j * b + k * c + 0];
                    if (info->disp_constraint[i * a + j * b + k * c + 1] == 0)
                    {
                        info->length_before[temp] += info->CP_info[patch_num * info->DIMENSION + 1];
                    }
                    else if (info->disp_constraint[i * a + j * b + k * c + 1] == 1)
                    {
                        info->length_before[temp] += info->CP_info[patch_num * info->DIMENSION];
                    }
                }
                temp++;
            }
        }
    }
    else if (info->DIMENSION == 3)
    {
        for (i = 0; i < info->DIMENSION; i++)
        {
            for (j = 0; j < info->disp_constraint_n[i]; j++)
            {
                for (k = 0; k < info->disp_constraint_face_edge_n[i * info->MAX_DISP_CONSTRAINT + j]; k++)
                {
                    int patch_num = info->disp_constraint[i * a + j * b + k * c + 0];
                    if (info->disp_constraint[i * a + j * b + k * c + 1] == 0)
                    {
                        info->length_before[temp] += info->CP_info[patch_num * info->DIMENSION + 1] * info->CP_info[patch_num * info->DIMENSION + 2];
                    }
                    else if (info->disp_constraint[i * a + j * b + k * c + 1] == 1)
                    {
                        info->length_before[temp] += info->CP_info[patch_num * info->DIMENSION] * info->CP_info[patch_num * info->DIMENSION + 2];
                    }
                    else if (info->disp_constraint[i * a + j * b + k * c + 1] == 2)
                    {
                        info->length_before[temp] += info->CP_info[patch_num * info->DIMENSION] * info->CP_info[patch_num * info->DIMENSION + 1];
                    }
                }
                temp++;
            }
        }
    }

    temp = 0;

    if (info->DIMENSION == 2)
    {
        for (i = 0; i < info->DIMENSION; i++)
        {
            for (j = 0; j < info->disp_constraint_n[i]; j++)
            {
                for (k = 0; k < info->disp_constraint_face_edge_n[i * info->MAX_DISP_CONSTRAINT + j]; k++)
                {
                    int A_to_here = 0;
                    for (l = 0; l < info->disp_constraint[i * a + j * b + k * c + 0]; l++)
                    {
                        A_to_here += 2 * (info->CP_info[l * info->DIMENSION] + info->CP_info[l * info->DIMENSION + 1]);
                    }

                    // if (info->disp_constraint[i * a + j * b + k * c + 1] == 1 && info->disp_constraint[i * a + j * b + k * c + 2] == 0) は何もしない
                    if (info->disp_constraint[i * a + j * b + k * c + 1] == 0 && info->disp_constraint[i * a + j * b + k * c + 2] == 1)
                    {
                        A_to_here += info->CP_info[info->disp_constraint[i * a + j * b + k * c + 0] * info->DIMENSION];
                    }
                    else if (info->disp_constraint[i * a + j * b + k * c + 1] == 1 && info->disp_constraint[i * a + j * b + k * c + 2] == 1)
                    {
                        A_to_here += info->CP_info[info->disp_constraint[i * a + j * b + k * c + 0] * info->DIMENSION] + info->CP_info[info->disp_constraint[i * a + j * b + k * c + 0] * info->DIMENSION + 1];
                    }
                    else if (info->disp_constraint[i * a + j * b + k * c + 1] == 0 && info->disp_constraint[i * a + j * b + k * c + 2] == 0)
                    {
                        A_to_here += 2 * info->CP_info[info->disp_constraint[i * a + j * b + k * c + 0] * info->DIMENSION] + info->CP_info[info->disp_constraint[i * a + j * b + k * c + 0] * info->DIMENSION + 1];
                    }

                    if (info->disp_constraint[i * a + j * b + k * c + 1] == 0)
                    {
                        for (l = 0; l < info->CP_info[info->disp_constraint[i * a + j * b + k * c + 0] * info->DIMENSION + 1]; l++)
                        {
                            info->Boundary[temp] = info->A[A_to_here + l];
                            temp++;
                        }
                    }
                    else if (info->disp_constraint[i * a + j * b + k * c + 1] == 1)
                    {
                        for (l = 0; l < info->CP_info[info->disp_constraint[i * a + j * b + k * c + 0] * info->DIMENSION]; l++)
                        {
                            info->Boundary[temp] = info->A[A_to_here + l];
                            temp++;
                        }
                    }
                }
            }
        }
    }
    else if (info->DIMENSION == 3)
    {
        for (i = 0; i < info->DIMENSION; i++)
        {
            for (j = 0; j < info->disp_constraint_n[i]; j++)
            {
                for (k = 0; k < info->disp_constraint_face_edge_n[i * info->MAX_DISP_CONSTRAINT + j]; k++)
                {
                    int A_to_here = 0;
                    for (l = 0; l < info->disp_constraint[i * a + j * b + k * c + 0]; l++)
                    {
                        A_to_here += 2 * (info->CP_info[l * info->DIMENSION] * info->CP_info[l * info->DIMENSION + 1] + info->CP_info[l * info->DIMENSION + 1] * info->CP_info[l * info->DIMENSION + 2] + info->CP_info[l * info->DIMENSION] * info->CP_info[l * info->DIMENSION + 2]);
                    }

                    int patch_num = info->disp_constraint[i * a + j * b + k * c + 0];
                    int temp1, temp2, temp3, temp4, temp5;
                    temp1 =         info->CP_info[patch_num * info->DIMENSION] * info->CP_info[patch_num * info->DIMENSION + 1];
                    temp2 = temp1 + info->CP_info[patch_num * info->DIMENSION] * info->CP_info[patch_num * info->DIMENSION + 1];
                    temp3 = temp2 + info->CP_info[patch_num * info->DIMENSION] * info->CP_info[patch_num * info->DIMENSION + 2];
                    temp4 = temp3 + info->CP_info[patch_num * info->DIMENSION + 1] * info->CP_info[patch_num * info->DIMENSION + 2];
                    temp5 = temp4 + info->CP_info[patch_num * info->DIMENSION] * info->CP_info[patch_num * info->DIMENSION + 2];

                    // if (info->disp_constraint[i * a + j * b + k * c + 1] == 2 && info->disp_constraint[i * a + j * b + k * c + 2] == 0) は何もしない
                    if (info->disp_constraint[i * a + j * b + k * c + 1] == 2 && info->disp_constraint[i * a + j * b + k * c + 2] == 1)
                    {
                        A_to_here += temp1;
                    }
                    else if (info->disp_constraint[i * a + j * b + k * c + 1] == 1 && info->disp_constraint[i * a + j * b + k * c + 2] == 0)
                    {
                        A_to_here += temp3;
                    }
                    else if (info->disp_constraint[i * a + j * b + k * c + 1] == 1 && info->disp_constraint[i * a + j * b + k * c + 2] == 1)
                    {
                        A_to_here += temp5;
                    }
                    else if (info->disp_constraint[i * a + j * b + k * c + 1] == 0 && info->disp_constraint[i * a + j * b + k * c + 2] == 0)
                    {
                        A_to_here += temp2;
                    }
                    else if (info->disp_constraint[i * a + j * b + k * c + 1] == 0 && info->disp_constraint[i * a + j * b + k * c + 2] == 1)
                    {
                        A_to_here += temp4;
                    }

                    if (info->disp_constraint[i * a + j * b + k * c + 1] == 0)
                    {
                        for (l = 0; l < info->CP_info[patch_num * info->DIMENSION + 1] * info->CP_info[patch_num * info->DIMENSION + 2]; l++)
                        {
                            info->Boundary[temp] = info->A[A_to_here + l];
                            temp++;
                        }
                    }
                    else if (info->disp_constraint[i * a + j * b + k * c + 1] == 1)
                    {
                        for (l = 0; l < info->CP_info[patch_num * info->DIMENSION + 0] * info->CP_info[patch_num * info->DIMENSION + 2]; l++)
                        {
                            info->Boundary[temp] = info->A[A_to_here + l];
                            temp++;
                        }
                    }
                    else if (info->disp_constraint[i * a + j * b + k * c + 1] == 2)
                    {
                        for (l = 0; l < info->CP_info[patch_num * info->DIMENSION] * info->CP_info[patch_num * info->DIMENSION + 1]; l++)
                        {
                            info->Boundary[temp] = info->A[A_to_here + l];
                            temp++;
                        }
                    }
                }
            }
        }
        
    }

    temp = 0;
    for (i = 0; i < n; i++)
    {
        printf("length_before[%d] = %d\n", i, info->length_before[i]);
        for (j = 0; j < info->length_before[i]; j++)
        {
            printf("%d\t", info->Boundary[temp]);
            temp++;
        }
        printf("\n");
    }

    temp = 0;
    for (i = 0; i < n; i++)
    {
        int *Array = (int *)malloc(sizeof(int) * info->length_before[i]);
        if (Array == NULL)
        {
            printf("Cannot allocate memory\n");
            exit(1);
        }

        for (j = 0; j < info->length_before[i]; j++)
        {
            Array[j] = info->Boundary[temp + j];
        }
        printf("\n");

        heapSort(Array, info->length_before[i]);
        Dedupe(Array, info->length_before, info->Boundary_result, info->length_after, i);

        free(Array);
        temp += info->length_before[i];
    }

    temp = 0;
    for (i = 0; i < n; i++)
    {
        printf("length_after[%d] = %d\n", i, info->length_after[i]);
        for (j = 0; j < info->length_after[i]; j++)
        {
            printf("%d\t", info->Boundary_result[i * temp + j]);
        }
        printf("\n");
        temp += info->length_before[i];
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