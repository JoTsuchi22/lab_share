#include <iostream>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "header_MI3D.h"

using namespace std;

int main(int argc, char **argv)
{
    int i, j;

    for (i = 0; i < argc - 1; i++)
    {
        // declare struct
        info_global info_glo, *info_glo_ptr;
        info_glo_ptr = &info_glo;

        // get dim
        Get_DIM(argv[i + 1], &info_glo);

        // declare struct
        info_each_DIMENSION *info = (info_each_DIMENSION *)malloc(sizeof(info_each_DIMENSION) * info_glo.DIMENSION);
        info_each_DIMENSION **info_ptr = (info_each_DIMENSION **)malloc(sizeof(info_each_DIMENSION *) * info_glo.DIMENSION);
        if (info == NULL || info_ptr == NULL)
        {
            printf("Cannot allocate memory\n"); exit(1);
        }

        for (j = 0; j < info_glo.DIMENSION; j++)
        {
            info_ptr[j] = &info[j];
        }

        // ファイル読み込み1回目
        Get_InputData_1(argv[i + 1], &info_glo, info);

        // memory allocation
        int MAX_Knot_Vector = 0;
        int MAX_Coordinate = 1;
        int temp = 0;
        if (info_glo.mode == 1)
        {
            for (j = 0; j < info_glo.DIMENSION; j++)
            {
                MAX_Knot_Vector += 2 * info[j].KI_cp_n;
                MAX_Coordinate *= 2 * info[j].KI_cp_n;
            }
        }
        else if (info_glo.mode == 0 || info_glo.mode == 2)
        {
            for (j = 0; j < info_glo.DIMENSION; j++)
            {
                MAX_Knot_Vector += info[j].knot_n * (info[j].OE_n + 1) + info[j].KI_non_uniform_n;
                MAX_Coordinate *= info[j].knot_n * (info[j].OE_n + 1) + info[j].KI_non_uniform_n;
            }
        }
        for (j = 0; j < info_glo.DIMENSION; j++)
        {
            temp += info[j].KI_non_uniform_n;
        }
        double *Knot_Vector = (double *)malloc(sizeof(double) * MAX_Knot_Vector);
        double *Coordinate = (double *)malloc(sizeof(double) * info_glo.DIMENSION * MAX_Coordinate);
        double *Weight = (double *)malloc(sizeof(double) * MAX_Coordinate);
        double *temp_Knot_Vector = (double *)malloc(sizeof(double) * MAX_Knot_Vector);
        double *temp_Coordinate = (double *)malloc(sizeof(double) * info_glo.DIMENSION * MAX_Coordinate);
        double *temp_Weight = (double *)malloc(sizeof(double) * MAX_Coordinate);
        double *Insert_Knot_Vector = (double *)malloc(sizeof(double) * temp);
        if (Knot_Vector == NULL || Coordinate == NULL || Weight == NULL || temp_Knot_Vector == NULL || temp_Coordinate == NULL || temp_Weight == NULL || Insert_Knot_Vector == NULL)
        {
            printf("Cannot allocate memory\n"); exit(1);
        }
    
        int temp1 = 0, temp2 = 0, temp3 = 0;
        for (j = 0; j < info_glo.DIMENSION; j++)
        {
            info_ptr[j]->KV = &Knot_Vector[temp1];
            info_ptr[j]->CP = &Coordinate[temp2];
            info_ptr[j]->temp_KV = &temp_Knot_Vector[temp1];
            info_ptr[j]->temp_CP = &temp_Coordinate[temp2];
            info_ptr[j]->insert_knot = &Insert_Knot_Vector[temp3];
            if (info_glo.mode == 1)
            {
                temp1 += 2 * info[j].KI_cp_n;
            }
            else if (info_glo.mode == 0 || info_glo.mode == 2)
            {
                temp1 += info[j].knot_n * (info[j].OE_n + 1) + info[j].KI_non_uniform_n;
            }
            temp2 += MAX_Coordinate;
            temp3 += info[j].KI_non_uniform_n;
        }
        info_glo_ptr->Weight = Weight;
        info_glo_ptr->temp_Weight = temp_Weight;

        // ファイル読み込み2回目
        Get_InputData_2(argv[i + 1], &info_glo, info);

        // オーダーエレベーション
        for (j = 0; j < info_glo.DIMENSION; j++)
        {
            if (info[j].OE_n != 0)
            {
                OE(j, &info_glo, info);
            }
        }

        // ノットインサーション
        for (j = 0; j < info_glo.DIMENSION; j++)
        {
            if (info_glo.mode == 0 && info[j].KI_non_uniform_n != 0)
            {
                KI_non_uniform(j, info[j].KI_non_uniform_n, info[j].insert_knot, &info_glo, info);
            }
            else if (info_glo.mode == 1)
            {
                // オープンノットベクトル
                if (info[j].knot_n == 2 * (info[j].Order + 1))
                {
                    KI_cp(j, &info_glo, info);
                }
                else
                {
                    KI_cp_not_open_knot_vec(j, &info_glo, info);
                }
            }
        }

        // 結果を出力
        OutputData(argv[i + 1], &info_glo, info);

        free(Knot_Vector), free(Coordinate), free(Weight), free(temp_Knot_Vector), free(temp_Coordinate), free(temp_Weight), free(Insert_Knot_Vector);
        free(info), free(info_ptr);
    }

    return 0;
}


// DIMENSION 読み込み
void Get_DIM(char *filename, info_global *info_glo)
{
    int temp_i;

    fp = fopen(filename, "r");

    // DIMENSION
    fscanf(fp, "%d", &temp_i);
    info_glo->DIMENSION = temp_i;

    fclose(fp);
}


// ファイル読み込み1回目
void Get_InputData_1(char *filename, info_global *info_glo, info_each_DIMENSION *info)
{
    int i, j, k;
    int temp_i;
    double temp_d;
    char s[256];

    fp = fopen(filename, "r");

    // DIMENSION(スキップ)
    fscanf(fp, "%d", &temp_i);
    fgets(s, 256, fp);

    // 各方向の次数
    for (j = 0; j < info_glo->DIMENSION; j++)
    {
        fscanf(fp, "%d", &temp_i);
        info[j].Order = temp_i;
    }
    fgets(s, 256, fp);

    // 各方向のノットベクトルの個数
    for (j = 0; j < info_glo->DIMENSION; j++)
    {
        fscanf(fp, "%d", &temp_i);
        info[j].knot_n = temp_i;
    }
    fgets(s, 256, fp);

    // 各方向のコントロールポイントの個数
    info_glo->Total_Control_Point = 1;
    for (j = 0; j < info_glo->DIMENSION; j++)
    {
        fscanf(fp, "%d", &temp_i);
        info[j].CP_n = temp_i;
        info_glo->Total_Control_Point *= temp_i;
    }
    fgets(s, 256, fp);

    // 各方向のノットベクトル(スキップ)
    for (j = 0; j < info_glo->DIMENSION; j++)
    {
        for (k = 0; k < info[j].knot_n; k++)
        {
            fscanf(fp, "%lf", &temp_d);
        }
    }
    fgets(s, 256, fp);

    // コントロールポイントの座標(スキップ)
    for (i = 0; i < info_glo->Total_Control_Point; i++)
    {
        fscanf(fp, "%d", &temp_i);
        for (j = 0; j < info_glo->DIMENSION + 1; j++)
        {
            fscanf(fp, "%lf", &temp_d);
        }
    }
    fgets(s, 256, fp);

    // 各方向のo.e.の回数
    for (j = 0; j < info_glo->DIMENSION; j++)
    {
        fscanf(fp, "%d", &temp_i);
        info[j].OE_n = temp_i;
    }
    fgets(s, 256, fp);

    // 各方向のk.i.(コントロールポイントに合わせて等分割)のコントロールポイント数
    for (j = 0; j < info_glo->DIMENSION; j++)
    {
        fscanf(fp, "%d", &temp_i);
        info[j].KI_cp_n = temp_i;
    }
    fgets(s, 256, fp);

    // 各方向のk.i.(任意のノット)の個数
    for (j = 0; j < info_glo->DIMENSION; j++)
    {
        fscanf(fp, "%d", &temp_i);
        info[j].KI_non_uniform_n = temp_i;
    }
    fgets(s, 256, fp);

    // コントロールポイントに合わせて分割と任意のノット挿入を同時に行っていないか確認
    int temp_check_a = 0;
    int temp_check_b = 0;
    int temp_check_total = 0;
    for (j = 0; j < info_glo->DIMENSION; j++)
    {
        if (info[j].KI_non_uniform_n != 0)
        {
            temp_check_a += 1;
        }
        if (info[j].KI_cp_n != 0)
        {
            temp_check_b += 1;
        }
    }
    if (temp_check_a > 0)
    {
        info_glo->mode = 0;
        temp_check_total++;
    }
    if (temp_check_b > 0)
    {
        info_glo->mode = 1;
        temp_check_total++;
    }
    if (temp_check_a == 0 && temp_check_b == 0)
    {
        info_glo->mode = 2;
    }

    if (temp_check_total > 1)
    {
        printf("インプットデータが間違っています．\nノットインサーションの各操作(2種類)は同時に行えません．\n操作を行わない該当するパラメータをすべて0としてください．\n");
        exit(1);
    }

    fclose(fp);
}


// ファイル読み込み2回目
void Get_InputData_2(char *filename, info_global *info_glo, info_each_DIMENSION *info)
{
    int i, j, k;
    int temp_i;
    double temp_d;
    char s[256];

    fp = fopen(filename, "r");

    // DIMENSION(スキップ)
    fscanf(fp, "%d", &temp_i);
    fgets(s, 256, fp);

    // 各方向の次数(スキップ)
    for (j = 0; j < info_glo->DIMENSION; j++)
    {
        fscanf(fp, "%d", &temp_i);
    }
    fgets(s, 256, fp);

    // 各方向のノットベクトルの個数(スキップ)
    for (j = 0; j < info_glo->DIMENSION; j++)
    {
        fscanf(fp, "%d", &temp_i);
    }
    fgets(s, 256, fp);

    // 各方向のコントロールポイントの個数(スキップ)
    for (j = 0; j < info_glo->DIMENSION; j++)
    {
        fscanf(fp, "%d", &temp_i);
    }
    fgets(s, 256, fp);

    // 各方向のノットベクトル
    for (j = 0; j < info_glo->DIMENSION; j++)
    {
        for (k = 0; k < info[j].knot_n; k++)
        {
            fscanf(fp, "%lf", &temp_d);
            info[j].KV[k] = temp_d;
        }
    }
    fgets(s, 256, fp);

    // コントロールポイントの座標
    for (i = 0; i < info_glo->Total_Control_Point; i++)
    {
        fscanf(fp, "%d", &temp_i);
        for (j = 0; j < info_glo->DIMENSION + 1; j++)
        {
            fscanf(fp, "%lf", &temp_d);
            if (j < info_glo->DIMENSION)
            {
                info[j].CP[i] = temp_d;
            }
            else if (j == info_glo->DIMENSION)
            {
                info_glo->Weight[i] = temp_d;
            }
        }
    }
    fgets(s, 256, fp);

    // 各方向のo.e.の回数(スキップ)
    for (j = 0; j < info_glo->DIMENSION; j++)
    {
        fscanf(fp, "%d", &temp_i);
    }
    fgets(s, 256, fp);

    // 各方向のk.i.(コントロールポイントに合わせて等分割)のコントロールポイント数(スキップ)
    for (j = 0; j < info_glo->DIMENSION; j++)
    {
        fscanf(fp, "%d", &temp_i);
    }
    fgets(s, 256, fp);

    // 各方向のk.i.(任意のノット)の個数(スキップ)
    for (j = 0; j < info_glo->DIMENSION; j++)
    {
        fscanf(fp, "%d", &temp_i);
    }
    fgets(s, 256, fp);

    // 各方向の挿入するノット
    for (j = 0; j < info_glo->DIMENSION; j++)
    {
        if (info[j].KI_non_uniform_n != 0)
        {
            for (k = 0; k < info[j].KI_non_uniform_n; k++)
            {
                fscanf(fp, "%lf", &temp_d);
                info[j].insert_knot[k] = temp_d;
            }
        }
    }

    fclose(fp);
}


// Knot Insertion
void KI_non_uniform(int insert_axis, int insert_knot_n, double *insert_knot_in_KI, info_global *info_glo, info_each_DIMENSION *info)
{
    int i, j, k, l;
    int other_axis = 0;
    int other_axis_1 = 0, other_axis_2 = 0;

    // calc knot
    KI_calc_knot_1D(insert_axis, insert_knot_n, insert_knot_in_KI, info);

    if (info_glo->DIMENSION == 2)
    {
        if (insert_axis == 0)
        {
            other_axis = 1;
        }
        else if (insert_axis == 1)
        {
            other_axis = 0;
        }

        for (i = 0; i < info[other_axis].CP_n; i++)
        {
            // memory allocation
            line_coordinate w, *w_ptr;
            line_coordinate *DIM = (line_coordinate *)malloc(sizeof(line_coordinate) * info_glo->DIMENSION);
            line_coordinate **DIM_ptr = (line_coordinate **)malloc(sizeof(line_coordinate *) * info_glo->DIMENSION);
            w_ptr = &w;
            for (j = 0; j < info_glo->DIMENSION; j++)
            {
                DIM_ptr[j] = &DIM[j];
            }
            double *line_x = (double *)malloc(sizeof(double) * info[insert_axis].CP_n);
            double *line_y = (double *)malloc(sizeof(double) * info[insert_axis].CP_n);
            double *line_w = (double *)malloc(sizeof(double) * info[insert_axis].CP_n);
            double *new_line_x = (double *)malloc(sizeof(double) * (info[insert_axis].CP_n + insert_knot_n));
            double *new_line_y = (double *)malloc(sizeof(double) * (info[insert_axis].CP_n + insert_knot_n));
            double *new_line_w = (double *)malloc(sizeof(double) * (info[insert_axis].CP_n + insert_knot_n));
            DIM_ptr[0]->line = line_x;
            DIM_ptr[1]->line = line_y;
            w_ptr->line = line_w;
            DIM_ptr[0]->new_line = new_line_x;
            DIM_ptr[1]->new_line = new_line_y;
            w_ptr->new_line = new_line_w;

            // get line info
            for (j = 0; j < info[insert_axis].CP_n; j++)
            {
                int temp = 0;
                if (insert_axis == 0)
                {
                    temp = i * info[0].CP_n + j;
                }
                else if (insert_axis == 1)
                {
                    temp = j * info[0].CP_n + i;
                }

                for (k = 0; k < info_glo->DIMENSION; k++)
                {
                    DIM[k].line[j] = info[k].CP[temp];
                }
                w.line[j] = info_glo->Weight[temp];
            }

            // calc line
            KI_calc_T_1D(insert_axis, insert_knot_n, info_glo, info, &w, DIM);

            // update
            for (j = 0; j < info[insert_axis].CP_n + insert_knot_n; j++)
            {
                int temp = 0;
                if (insert_axis == 0)
                {
                    temp = i * (info[0].CP_n + insert_knot_n) + j;
                }
                else if (insert_axis == 1)
                {
                    temp = j * info[0].CP_n + i;
                }

                for (k = 0; k < info_glo->DIMENSION; k++)
                {
                    info[k].temp_CP[temp] = DIM[k].new_line[j];
                }
                info_glo->temp_Weight[temp] = w.new_line[j];
            }

            // memory free
            free(line_x), free(line_y), free(line_w);
            free(new_line_x), free(new_line_y), free(new_line_w);
            free(DIM), free(DIM_ptr);
        }
    }
    else if(info_glo->DIMENSION == 3)
    {
        if (insert_axis == 0)
        {
            other_axis_1 = 1, other_axis_2 = 2;
        }
        else if (insert_axis == 1)
        {
            other_axis_1 = 0, other_axis_2 = 2;
        }
        else if (insert_axis == 2)
        {
            other_axis_1 = 0, other_axis_2 = 1;
        }

        for (i = 0; i < info[other_axis_1].CP_n; i++)
        {
            for (j = 0; j < info[other_axis_2].CP_n; j++)
            {
                // memory allocation
                line_coordinate w, *w_ptr;
                line_coordinate *DIM = (line_coordinate *)malloc(sizeof(line_coordinate) * info_glo->DIMENSION);
                line_coordinate **DIM_ptr = (line_coordinate **)malloc(sizeof(line_coordinate *) * info_glo->DIMENSION);
                w_ptr = &w;
                for (k = 0; k < info_glo->DIMENSION; k++)
                {
                    DIM_ptr[k] = &DIM[k];
                }
                double *line_x = (double *)malloc(sizeof(double) * info[insert_axis].CP_n);
                double *line_y = (double *)malloc(sizeof(double) * info[insert_axis].CP_n);
                double *line_z = (double *)malloc(sizeof(double) * info[insert_axis].CP_n);
                double *line_w = (double *)malloc(sizeof(double) * info[insert_axis].CP_n);
                double *new_line_x = (double *)malloc(sizeof(double) * (info[insert_axis].CP_n + insert_knot_n));
                double *new_line_y = (double *)malloc(sizeof(double) * (info[insert_axis].CP_n + insert_knot_n));
                double *new_line_z = (double *)malloc(sizeof(double) * (info[insert_axis].CP_n + insert_knot_n));
                double *new_line_w = (double *)malloc(sizeof(double) * (info[insert_axis].CP_n + insert_knot_n));
                DIM_ptr[0]->line = line_x;
                DIM_ptr[1]->line = line_y;
                DIM_ptr[2]->line = line_z;
                w_ptr->line = line_w;
                DIM_ptr[0]->new_line = new_line_x;
                DIM_ptr[1]->new_line = new_line_y;
                DIM_ptr[2]->new_line = new_line_z;
                w_ptr->new_line = new_line_w;

                // get line info
                for (k = 0; k < info[insert_axis].CP_n; k++)
                {
                    int temp = 0;
                    if (insert_axis == 0)
                    {
                        temp = j * (info[1].CP_n * info[0].CP_n) + i * info[0].CP_n + k;
                    }
                    else if (insert_axis == 1)
                    {
                        temp = j * (info[1].CP_n * info[0].CP_n) + k * info[0].CP_n + i;
                    }
                    else if (insert_axis == 2)
                    {
                        temp = k * (info[1].CP_n * info[0].CP_n) + j * info[0].CP_n + i;
                    }

                    for (l = 0; l < info_glo->DIMENSION; l++)
                    {
                        DIM[l].line[k] = info[l].CP[temp];
                    }
                    w.line[k] = info_glo->Weight[temp];
                }

                // calc line
                KI_calc_T_1D(insert_axis, insert_knot_n, info_glo, info, &w, DIM);

                // update
                for (k = 0; k < info[insert_axis].CP_n + insert_knot_n; k++)
                {
                    int temp = 0;
                    if (insert_axis == 0)
                    {
                        temp = j * (info[1].CP_n * (info[0].CP_n + insert_knot_n)) + i * (info[0].CP_n + insert_knot_n) + k;
                    }
                    else if (insert_axis == 1)
                    {
                        temp = j * ((info[1].CP_n + insert_knot_n) * info[0].CP_n) + k * info[0].CP_n + i;
                    }
                    else if (insert_axis == 2)
                    {
                        temp = k * (info[1].CP_n * info[0].CP_n) + j * info[0].CP_n + i;
                    }

                    for (l = 0; l < info_glo->DIMENSION; l++)
                    {
                        info[l].temp_CP[temp] = DIM[l].new_line[k];
                    }
                    info_glo->temp_Weight[temp] = w.new_line[k];
                }

                // memory free
                free(line_x), free(line_y), free(line_z), free(line_w);
                free(new_line_x), free(new_line_y), free(new_line_z), free(new_line_w);
                free(DIM), free(DIM_ptr);
            }
        }
    }

    // update
    for (i = 0; i < info[insert_axis].knot_n + insert_knot_n; i++)
    {
        info[insert_axis].KV[i] = info[insert_axis].temp_KV[i];
    }
    info[insert_axis].CP_n += insert_knot_n;
    info[insert_axis].knot_n += insert_knot_n;
    int temp_CP_n = 1;
    for (i = 0; i < info_glo->DIMENSION; i++)
    {
        temp_CP_n *= info[i].CP_n;
    }
    info_glo->Total_Control_Point = temp_CP_n;
    for (i = 0; i < info_glo->Total_Control_Point; i++)
    {
        for (j = 0; j < info_glo->DIMENSION; j++)
        {
            info[j].CP[i] = info[j].temp_CP[i];
        }
        info_glo->Weight[i] = info_glo->temp_Weight[i];
    }
}


void KI_calc_knot_1D(int insert_axis, int insert_knot_n, double *insert_knot_in_KI, info_each_DIMENSION *info)
{
    int i;
    int temp1 = 0, temp2 = 0, temp3 = 0;

    for (i = 0; i < info[insert_axis].knot_n + insert_knot_n; i++)
    {
        if (temp2 < insert_knot_n)
        {
            if (info[insert_axis].KV[temp1] < insert_knot_in_KI[temp2])
            {
                info[insert_axis].temp_KV[temp3] = info[insert_axis].KV[temp1];
                temp1++;
                temp3++;
            }
            else if (info[insert_axis].KV[temp1] > insert_knot_in_KI[temp2])
            {
                info[insert_axis].temp_KV[temp3] = insert_knot_in_KI[temp2];
                temp2++;
                temp3++;
            }
            else if (info[insert_axis].KV[temp1] == insert_knot_in_KI[temp2])
            {
                info[insert_axis].temp_KV[temp3] = info[insert_axis].KV[temp1];
                info[insert_axis].temp_KV[temp3 + 1] = insert_knot_in_KI[temp2];
                temp1++;
                temp2++;
                temp3++;
                temp3++;
            }
        }
        else
        {
            info[insert_axis].temp_KV[temp3] = info[insert_axis].KV[temp1];
            temp1++;
            temp3++;
        }
    }
}


void KI_calc_T_1D(int insert_axis, int insert_knot_n, info_global *info_glo, info_each_DIMENSION *info, line_coordinate *w, line_coordinate *DIM)
{
    int i, j, k;

    int l = info[insert_axis].CP_n;
    int n = info[insert_axis].Order;

    // pre processing
    for (i = 0; i < l; i++)
    {
        for (j = 0; j < info_glo->DIMENSION; j++)
        {
            DIM[j].line[i] = DIM[j].line[i] * w->line[i];
        }
    }

    double *T = (double *)calloc((n + 1) * (l + insert_knot_n) * l, sizeof(double));

    double a = 0.0, b = 0.0;
    for (i = 0; i < l + insert_knot_n; i++)
    {
        for (k = 0; k < n + 1; k++)
        {
            for (j = 0; j < l; j++)
            {
                if (k == 0)
                {
                    if (info[insert_axis].KV[j] <= info[insert_axis].temp_KV[i] && info[insert_axis].temp_KV[i] < info[insert_axis].KV[j + 1])
                    {
                        T[k * (l + insert_knot_n) * l + i * l + j] = 1.0;
                    }
                    else
                    {
                        T[k * (l + insert_knot_n) * l + i * l + j] = 0.0;
                    }
                }
                else
                {
                    if (info[insert_axis].KV[j + k] - info[insert_axis].KV[j] == 0)
                    {
                        a = 0.0;
                    }
                    else if (info[insert_axis].KV[j + k] - info[insert_axis].KV[j] != 0)
                    {
                        a = ((info[insert_axis].temp_KV[i + k] - info[insert_axis].KV[j]) / (info[insert_axis].KV[j + k] - info[insert_axis].KV[j])) * T[(k - 1) * (l + insert_knot_n) * l + i * l + j];
                    }
                    if (info[insert_axis].KV[j + k + 1] - info[insert_axis].KV[j + 1] == 0)
                    {
                        b = 0.0;
                    }
                    else if (info[insert_axis].KV[j + k + 1] - info[insert_axis].KV[j + 1] != 0)
                    {
                        b = ((info[insert_axis].KV[j + k + 1] - info[insert_axis].temp_KV[i + k]) / (info[insert_axis].KV[j + k + 1] - info[insert_axis].KV[j + 1])) * T[(k - 1) * (l + insert_knot_n) * l + i * l + j + 1];
                    }
                    T[k * (l + insert_knot_n) * l + i * l + j] = a + b;
                }
            }
        }
    }

    // inner product
    for (i = 0; i < l + insert_knot_n; i++)
    {
        for (j = 0; j < l; j++)
        {
            if (j == 0)
            {
                for (k = 0; k < info_glo->DIMENSION; k++)
                {
                    DIM[k].new_line[i] = 0.0;
                }
                w->new_line[i] = 0.0;
            }

            for (k = 0; k < info_glo->DIMENSION; k++)
            {
                DIM[k].new_line[i] += T[n * (l + insert_knot_n) * l + i * l + j] * DIM[k].line[j];
            }
            w->new_line[i] += T[n * (l + insert_knot_n) * l + i * l + j] * w->line[j];
        }
    }

    // post processing
    for (i = 0; i < l + insert_knot_n; i++)
    {
        for (j = 0; j < info_glo->DIMENSION; j++)
        {
            DIM[j].new_line[i] = DIM[j].new_line[i] / w->new_line[i];
        }
    }

    free(T);
}


void KI_cp(int insert_axis, info_global *info_glo, info_each_DIMENSION *info)
{
    int i;
    int insert_knot_n = info[insert_axis].KI_cp_n - info[insert_axis].CP_n;
    double *insert_knot_in_KI = (double *)malloc(sizeof(double) * insert_knot_n);

    for (i = 0; i < insert_knot_n; i++)
    {
        insert_knot_in_KI[i] = (i + 1.0) / (insert_knot_n + 1.0);
    }
    KI_non_uniform(insert_axis, insert_knot_n, insert_knot_in_KI, info_glo, info);

    free(insert_knot_in_KI);
}


void KI_cp_not_open_knot_vec(int insert_axis, info_global *info_glo, info_each_DIMENSION *info)
{
    int i;
    int insert_knot_n = info[insert_axis].KI_cp_n - (info[insert_axis].CP_n - (info[insert_axis].knot_n - 2 * (info[insert_axis].Order + 1)));
    int removal_knot_n = info[insert_axis].knot_n - 2 * (info[insert_axis].Order + 1);
    double *insert_knot_in_KI = (double *)malloc(sizeof(double) * insert_knot_n);
    double *removal_knot_in_KI = (double *)malloc(sizeof(double) * removal_knot_n);

    for (i = 0; i < insert_knot_n; i++)
    {
        insert_knot_in_KI[i] = (i + 1.0) / (insert_knot_n + 1.0);
    }
    for (i = 0; i < removal_knot_n; i++)
    {
        removal_knot_in_KI[i] = info[insert_axis].KV[i + info[insert_axis].Order + 1];
    }

    KI_non_uniform(insert_axis, insert_knot_n, insert_knot_in_KI, info_glo, info);
    KR_non_uniform(insert_axis, removal_knot_n, removal_knot_in_KI, info_glo, info);

    free(insert_knot_in_KI);
}


// Order Elevation
void OE(int elevation_axis, info_global *info_glo, info_each_DIMENSION *info)
{
    int i, j, k, l;
    int other_axis = 0;
    int other_axis_1 = 0, other_axis_2 = 0;
    int insert_knot_n = 0, removal_knot_n = 0;
    double *insert_knot = (double *)malloc(sizeof(double) * info[elevation_axis].knot_n * (info[elevation_axis].OE_n + 2));
    double *removal_knot = (double *)malloc(sizeof(double) * info[elevation_axis].knot_n * (info[elevation_axis].OE_n + 2));

    // calc knot
    Calc_insert_knot_in_OE(elevation_axis, &insert_knot_n, insert_knot, &removal_knot_n, removal_knot, info);

    // knot insertion
    KI_non_uniform(elevation_axis, insert_knot_n, insert_knot, info_glo, info);

    if (info_glo->DIMENSION == 2)
    {
        if (elevation_axis == 0)
        {
            other_axis = 1;
        }
        else if (elevation_axis == 1)
        {
            other_axis = 0;
        }

        for (i = 0; i < info[other_axis].CP_n; i++)
        {
            // memory allocation
            line_coordinate w, *w_ptr;
            line_coordinate *DIM = (line_coordinate *)malloc(sizeof(line_coordinate) * info_glo->DIMENSION);
            line_coordinate **DIM_ptr = (line_coordinate **)malloc(sizeof(line_coordinate *) * info_glo->DIMENSION);
            w_ptr = &w;
            for (j = 0; j < info_glo->DIMENSION; j++)
            {
                DIM_ptr[j] = &DIM[j];
            }
            double *line_x = (double *)malloc(sizeof(double) * info[elevation_axis].CP_n * (info[elevation_axis].OE_n + 2));
            double *line_y = (double *)malloc(sizeof(double) * info[elevation_axis].CP_n * (info[elevation_axis].OE_n + 2));
            double *line_w = (double *)malloc(sizeof(double) * info[elevation_axis].CP_n * (info[elevation_axis].OE_n + 2));
            double *new_line_x = (double *)malloc(sizeof(double) * info[elevation_axis].CP_n * (info[elevation_axis].OE_n + 2));
            double *new_line_y = (double *)malloc(sizeof(double) * info[elevation_axis].CP_n * (info[elevation_axis].OE_n + 2));
            double *new_line_w = (double *)malloc(sizeof(double) * info[elevation_axis].CP_n * (info[elevation_axis].OE_n + 2));
            DIM_ptr[0]->line = line_x;
            DIM_ptr[1]->line = line_y;
            w_ptr->line = line_w;
            DIM_ptr[0]->new_line = new_line_x;
            DIM_ptr[1]->new_line = new_line_y;
            w_ptr->new_line = new_line_w;

            // get line info
            for (j = 0; j < info[elevation_axis].CP_n; j++)
            {
                int temp = 0;
                if (elevation_axis == 0)
                {
                    temp = i * info[0].CP_n + j;
                }
                else if (elevation_axis == 1)
                {
                    temp = j * info[0].CP_n + i;
                }

                for (k = 0; k < info_glo->DIMENSION; k++)
                {
                    DIM[k].line[j] = info[k].CP[temp];
                }
                w.line[j] = info_glo->Weight[temp];
            }

            // calc line
            Calc_Bezier(elevation_axis, info_glo, info, &w, DIM);

            // update
            int temp_CP_n = info[elevation_axis].OE_n * ((info[elevation_axis].CP_n - 1) / info[elevation_axis].Order);
            for (j = 0; j < info[elevation_axis].CP_n + temp_CP_n; j++)
            {
                int temp = 0;
                if (elevation_axis == 0)
                {
                    temp = i * (info[0].CP_n + temp_CP_n) + j;
                }
                else if (elevation_axis == 1)
                {
                    temp = j * info[0].CP_n + i;
                }

                for (k = 0; k < info_glo->DIMENSION; k++)
                {
                    info[k].temp_CP[temp] = DIM[k].new_line[j];
                }
                info_glo->temp_Weight[temp] = w.new_line[j];
            }

            // memory free
            free(line_x), free(line_y), free(line_w);
            free(new_line_x), free(new_line_y), free(new_line_w);
            free(DIM), free(DIM_ptr);
        }
    }
    else if(info_glo->DIMENSION == 3)
    {
        if (elevation_axis == 0)
        {
            other_axis_1 = 1, other_axis_2 = 2;
        }
        else if (elevation_axis == 1)
        {
            other_axis_1 = 0, other_axis_2 = 2;
        }
        else if (elevation_axis == 2)
        {
            other_axis_1 = 0, other_axis_2 = 1;
        }

        for (i = 0; i < info[other_axis_1].CP_n; i++)
        {
            for (j = 0; j < info[other_axis_2].CP_n; j++)
            {
                // memory allocation
                line_coordinate w, *w_ptr;
                line_coordinate *DIM = (line_coordinate *)malloc(sizeof(line_coordinate) * info_glo->DIMENSION);
                line_coordinate **DIM_ptr = (line_coordinate **)malloc(sizeof(line_coordinate *) * info_glo->DIMENSION);
                w_ptr = &w;
                for (k = 0; k < info_glo->DIMENSION; k++)
                {
                    DIM_ptr[k] = &DIM[k];
                }
                double *line_x = (double *)malloc(sizeof(double) * info[elevation_axis].CP_n * (info[elevation_axis].OE_n + 2));
                double *line_y = (double *)malloc(sizeof(double) * info[elevation_axis].CP_n * (info[elevation_axis].OE_n + 2));
                double *line_z = (double *)malloc(sizeof(double) * info[elevation_axis].CP_n * (info[elevation_axis].OE_n + 2));
                double *line_w = (double *)malloc(sizeof(double) * info[elevation_axis].CP_n * (info[elevation_axis].OE_n + 2));
                double *new_line_x = (double *)malloc(sizeof(double) * info[elevation_axis].CP_n * (info[elevation_axis].OE_n + 2));
                double *new_line_y = (double *)malloc(sizeof(double) * info[elevation_axis].CP_n * (info[elevation_axis].OE_n + 2));
                double *new_line_z = (double *)malloc(sizeof(double) * info[elevation_axis].CP_n * (info[elevation_axis].OE_n + 2));
                double *new_line_w = (double *)malloc(sizeof(double) * info[elevation_axis].CP_n * (info[elevation_axis].OE_n + 2));
                DIM_ptr[0]->line = line_x;
                DIM_ptr[1]->line = line_y;
                DIM_ptr[2]->line = line_z;
                w_ptr->line = line_w;
                DIM_ptr[0]->new_line = new_line_x;
                DIM_ptr[1]->new_line = new_line_y;
                DIM_ptr[2]->new_line = new_line_z;
                w_ptr->new_line = new_line_w;

                // get line info
                for (k = 0; k < info[elevation_axis].CP_n; k++)
                {
                    int temp = 0;
                    if (elevation_axis == 0)
                    {
                        temp = j * (info[1].CP_n * info[0].CP_n) + i * info[0].CP_n + k;
                    }
                    else if (elevation_axis == 1)
                    {
                        temp = j * (info[1].CP_n * info[0].CP_n) + k * info[0].CP_n + i;
                    }
                    else if (elevation_axis == 2)
                    {
                        temp = k * (info[1].CP_n * info[0].CP_n) + j * info[0].CP_n + i;
                    }

                    for (l = 0; l < info_glo->DIMENSION; l++)
                    {
                        DIM[l].line[k] = info[l].CP[temp];
                    }
                    w.line[k] = info_glo->Weight[temp];
                }

                // calc line
                Calc_Bezier(elevation_axis, info_glo, info, &w, DIM);

                // update
                int temp_CP_n = info[elevation_axis].OE_n * ((info[elevation_axis].CP_n - 1) / info[elevation_axis].Order);
                for (k = 0; k < info[elevation_axis].CP_n + temp_CP_n; k++)
                {
                    int temp = 0;
                    if (elevation_axis == 0)
                    {
                        temp = j * (info[1].CP_n * (info[0].CP_n + temp_CP_n)) + i * (info[0].CP_n + temp_CP_n) + k;
                    }
                    else if (elevation_axis == 1)
                    {
                        temp = j * ((info[1].CP_n + temp_CP_n) * info[0].CP_n) + k * info[0].CP_n + i;
                    }
                    else if (elevation_axis == 2)
                    {
                        temp = k * (info[1].CP_n * info[0].CP_n) + j * info[0].CP_n + i;
                    }

                    for (l = 0; l < info_glo->DIMENSION; l++)
                    {
                        info[l].temp_CP[temp] = DIM[l].new_line[k];
                    }
                    info_glo->temp_Weight[temp] = w.new_line[k];
                }

                // memory free
                free(line_x), free(line_y), free(line_z), free(line_w);
                free(new_line_x), free(new_line_y), free(new_line_z), free(new_line_w);
                free(DIM), free(DIM_ptr);
            }
        }
    }

    // update
    info[elevation_axis].CP_n += info[elevation_axis].OE_n * ((info[elevation_axis].CP_n - 1) / info[elevation_axis].Order);
    info[elevation_axis].Order += info[elevation_axis].OE_n;
    int temp_CP_n = 1;
    for (i = 0; i < info_glo->DIMENSION; i++)
    {
        temp_CP_n *= info[i].CP_n;
    }
    info_glo->Total_Control_Point = temp_CP_n;
    for (i = 0; i < info_glo->Total_Control_Point; i++)
    {
        for (j = 0; j < info_glo->DIMENSION; j++)
        {
            info[j].CP[i] = info[j].temp_CP[i];
        }
        info_glo->Weight[i] = info_glo->temp_Weight[i];
    }
    int temp1 = 0;
    double *temp_knot1 = (double *)calloc(info[elevation_axis].knot_n * (info[elevation_axis].OE_n + 2), sizeof(double));
    double *temp_knot2 = (double *)calloc(info[elevation_axis].knot_n * (info[elevation_axis].OE_n + 2), sizeof(double));
    for (i = 0; i < info[elevation_axis].knot_n; i++)
    {
        if (i == 0)
        {
            temp_knot1[temp1] = info[elevation_axis].KV[i];
            temp1++;
        }
        if (temp_knot1[temp1 - 1] < info[elevation_axis].KV[i])
        {
            temp_knot1[temp1] = info[elevation_axis].KV[i];
            temp1++;
        }
        if (temp_knot1[temp1] == 1.0)
        {
            break;
        }
    }
    for (i = 0; i < info[elevation_axis].OE_n; i++)
    {
        int temp2 = 0, temp3 = 0;
        for (j = 0; j < info[elevation_axis].knot_n + temp1; j++)
        {
            if (temp2 < temp1 && temp3 < info[elevation_axis].knot_n)
            {
                if (temp_knot1[temp2] <= info[elevation_axis].KV[temp3])
                {
                    temp_knot2[temp2 + temp3] = temp_knot1[temp2];
                    temp2++;
                }
                else if (temp_knot1[temp2] > info[elevation_axis].KV[temp3])
                {
                    temp_knot2[temp2 + temp3] = info[elevation_axis].KV[temp3];
                    temp3++;
                }
            }
            else if (temp2 == temp1)
            {
                temp_knot2[temp2 + temp3] = info[elevation_axis].KV[temp3];
                temp3++;
            }
            else if (temp3 == info[elevation_axis].knot_n)
            {
                temp_knot2[temp2 + temp3] = temp_knot1[temp2];
                temp2++;
            }

            if (temp2 + temp3 == info[elevation_axis].knot_n + temp1)
            {
                break;
            }
        }

        info[elevation_axis].knot_n += temp1;

        for (j = 0; j < info[elevation_axis].knot_n; j++)
        {
            info[elevation_axis].KV[j] = temp_knot2[j];
        }
    }
    free(temp_knot1), free(temp_knot2);

    // knot removal
    KR_non_uniform(elevation_axis, removal_knot_n, removal_knot, info_glo, info);

    free(insert_knot), free(removal_knot);
}


void Calc_insert_knot_in_OE(int elevation_axis, int *insert_knot_n, double *insert_knot, int *removal_knot_n, double *removal_knot, info_each_DIMENSION *info)
{
    int i, j;
    int m = info[elevation_axis].knot_n;
    int n = info[elevation_axis].Order;

    double *vec_temp1 = (double *)malloc(sizeof(double) * info[elevation_axis].knot_n * (info[elevation_axis].OE_n + 2));
    double *vec_temp2 = (double *)malloc(sizeof(double) * info[elevation_axis].knot_n * (info[elevation_axis].OE_n + 2));
    double *vec_unique = (double *)malloc(sizeof(double) * info[elevation_axis].knot_n * (info[elevation_axis].OE_n + 2));

    int temp1 = 0;
    for (i = 1; i < m - 1; i++)
    {
        if (info[elevation_axis].KV[i - 1] != info[elevation_axis].KV[i] && info[elevation_axis].KV[i] != info[elevation_axis].KV[i + 1])
        {
            vec_unique[temp1] = info[elevation_axis].KV[i];
            temp1++;
        }
    }

    int temp3 = 0, temp4 = 0, temp_count = 0;
    for (i = 0; i < temp1; i++)
    {
        int temp2 = 0;
        for (j = 0; j < m; j++)
        {
            if (vec_unique[i] == info[elevation_axis].KV[j])
            {
                temp2++;
            }
        }

        temp_count = temp2;
        for (;;)
        {
            if (temp_count < n)
            {
                vec_temp1[temp3] = vec_unique[i];
                temp3++;
                temp_count++;
            }
            else
            {
                break;
            }
        }
        temp_count = temp2;
        for (;;)
        {
            if (temp_count < n + info[elevation_axis].OE_n)
            {
                vec_temp2[temp4] = vec_unique[i];
                temp4++;
                temp_count++;
            }
            else
            {
                break;
            }
        }
    }

    for (i = 0; i < temp3; i++)
    {
        insert_knot[i] = vec_temp1[i];
    }
    *insert_knot_n = temp3;

    for (i = 0; i < temp4; i++)
    {
        removal_knot[i] = vec_temp2[i];
    }
    *removal_knot_n = temp4;

    free(vec_temp1), free(vec_temp2), free(vec_unique);
}


void Bezier_Order_Elevation(int elevation_axis, int Bezier_line, int *counter, info_global *info_glo, info_each_DIMENSION *info, line_coordinate *w, line_coordinate *DIM)
{
    int i, j, k;
    int n = info[elevation_axis].Order;

    if (info_glo->DIMENSION == 2)
    {
        // memory allocation
        line_coordinate B_w, *B_w_ptr;
        line_coordinate *B_DIM = (line_coordinate *)malloc(sizeof(line_coordinate) * info_glo->DIMENSION);
        line_coordinate **B_DIM_ptr = (line_coordinate **)malloc(sizeof(line_coordinate *) * info_glo->DIMENSION);
        B_w_ptr = &B_w;
        for (k = 0; k < info_glo->DIMENSION; k++)
        {
            B_DIM_ptr[k] = &B_DIM[k];
        }
        double *line_x = (double *)malloc(sizeof(double) * info[elevation_axis].CP_n * (info[elevation_axis].OE_n + 2));
        double *line_y = (double *)malloc(sizeof(double) * info[elevation_axis].CP_n * (info[elevation_axis].OE_n + 2));
        double *line_w = (double *)malloc(sizeof(double) * info[elevation_axis].CP_n * (info[elevation_axis].OE_n + 2));
        double *new_line_x = (double *)malloc(sizeof(double) * info[elevation_axis].CP_n * (info[elevation_axis].OE_n + 2));
        double *new_line_y = (double *)malloc(sizeof(double) * info[elevation_axis].CP_n * (info[elevation_axis].OE_n + 2));
        double *new_line_w = (double *)malloc(sizeof(double) * info[elevation_axis].CP_n * (info[elevation_axis].OE_n + 2));
        B_DIM_ptr[0]->line = line_x;
        B_DIM_ptr[1]->line = line_y;
        B_w_ptr->line = line_w;
        B_DIM_ptr[0]->new_line = new_line_x;
        B_DIM_ptr[1]->new_line = new_line_y;
        B_w_ptr->new_line = new_line_w;

        for (i = 0; i < n + 1; i++)
        {
            for (k = 0; k < info_glo->DIMENSION; k++)
            {
                B_DIM[k].line[i] = DIM[k].line[Bezier_line + i];
            }
            B_w.line[i] = w->line[Bezier_line + i];
        }

        for (i = 0; i < info[elevation_axis].OE_n; i++)
        {
            double a, b, a_w, b_w;

            for (j = 0; j < n + 1; j++)
            {
                if (j != n)
                {
                    double alpha = (j + 1.0) / (n + 1.0);

                    a_w = (1.0 - alpha) * B_w.line[j + 1];
                    b_w = alpha * B_w.line[j];
                    B_w.new_line[j] = a_w + b_w;
                    for (k = 0; k < info_glo->DIMENSION; k++)
                    {
                        a = (1.0 - alpha) * (B_DIM[k].line[j + 1] * B_w.line[j + 1]);
                        b = alpha * (B_DIM[k].line[j] * B_w.line[j]);
                        B_DIM[k].new_line[j] = (a + b) / B_w.new_line[j];
                    }
                }
                else if (j == n)
                {
                    for (k = 0; k < info_glo->DIMENSION; k++)
                    {
                        B_DIM[k].new_line[j] = B_DIM[k].line[j];
                    }
                    B_w.new_line[j] = B_w.line[j];
                }
            }
            if (i != info[elevation_axis].OE_n - 1 && info[elevation_axis].OE_n != 1)
            {
                for (j = 0; j < n + 1; j++)
                {
                    for (k = 0; k < info_glo->DIMENSION; k++)
                    {
                        B_DIM[k].line[j + 1] = B_DIM[k].new_line[j];
                    }
                    B_w.line[j + 1] = B_w.new_line[j];
                }
                n++;
            }
        }

        for (i = 0; i < n + 1; i++)
        {
            for (k = 0; k < info_glo->DIMENSION; k++)
            {
                DIM[k].new_line[*counter] = B_DIM[k].new_line[i];
            }
            w->new_line[*counter] = B_w.new_line[i];
            *counter += 1;
        }

        free(line_x), free(line_y), free(line_w);
        free(new_line_x), free(new_line_y), free(new_line_w);
        free(B_DIM), free(B_DIM_ptr);
    }
    else if (info_glo->DIMENSION == 3)
    {
        // memory allocation
        line_coordinate B_w, *B_w_ptr;
        line_coordinate *B_DIM = (line_coordinate *)malloc(sizeof(line_coordinate) * info_glo->DIMENSION);
        line_coordinate **B_DIM_ptr = (line_coordinate **)malloc(sizeof(line_coordinate *) * info_glo->DIMENSION);
        B_w_ptr = &B_w;
        for (k = 0; k < info_glo->DIMENSION; k++)
        {
            B_DIM_ptr[k] = &B_DIM[k];
        }
        double *line_x = (double *)malloc(sizeof(double) * info[elevation_axis].CP_n * (info[elevation_axis].OE_n + 2));
        double *line_y = (double *)malloc(sizeof(double) * info[elevation_axis].CP_n * (info[elevation_axis].OE_n + 2));
        double *line_z = (double *)malloc(sizeof(double) * info[elevation_axis].CP_n * (info[elevation_axis].OE_n + 2));
        double *line_w = (double *)malloc(sizeof(double) * info[elevation_axis].CP_n * (info[elevation_axis].OE_n + 2));
        double *new_line_x = (double *)malloc(sizeof(double) * info[elevation_axis].CP_n * (info[elevation_axis].OE_n + 2));
        double *new_line_y = (double *)malloc(sizeof(double) * info[elevation_axis].CP_n * (info[elevation_axis].OE_n + 2));
        double *new_line_z = (double *)malloc(sizeof(double) * info[elevation_axis].CP_n * (info[elevation_axis].OE_n + 2));
        double *new_line_w = (double *)malloc(sizeof(double) * info[elevation_axis].CP_n * (info[elevation_axis].OE_n + 2));
        B_DIM_ptr[0]->line = line_x;
        B_DIM_ptr[1]->line = line_y;
        B_DIM_ptr[2]->line = line_z;
        B_w_ptr->line = line_w;
        B_DIM_ptr[0]->new_line = new_line_x;
        B_DIM_ptr[1]->new_line = new_line_y;
        B_DIM_ptr[2]->new_line = new_line_z;
        B_w_ptr->new_line = new_line_w;

        for (i = 0; i < n + 1; i++)
        {
            for (k = 0; k < info_glo->DIMENSION; k++)
            {
                B_DIM[k].line[i] = DIM[k].line[Bezier_line + i];
            }
            B_w.line[i] = w->line[Bezier_line + i];
        }

        for (i = 0; i < info[elevation_axis].OE_n; i++)
        {
            double a, b, a_w, b_w;

            for (j = 0; j < n + 1; j++)
            {
                if (j != n)
                {
                    double alpha = (j + 1.0) / (n + 1.0);

                    a_w = (1.0 - alpha) * B_w.line[j + 1];
                    b_w = alpha * B_w.line[j];
                    B_w.new_line[j] = a_w + b_w;
                    for (k = 0; k < info_glo->DIMENSION; k++)
                    {
                        a = (1.0 - alpha) * (B_DIM[k].line[j + 1] * B_w.line[j + 1]);
                        b = alpha * (B_DIM[k].line[j] * B_w.line[j]);
                        B_DIM[k].new_line[j] = (a + b) / B_w.new_line[j];
                    }
                }
                else if (j == n)
                {
                    for (k = 0; k < info_glo->DIMENSION; k++)
                    {
                        B_DIM[k].new_line[j] = B_DIM[k].line[j];
                    }
                    B_w.new_line[j] = B_w.line[j];
                }
            }
            if (i != info[elevation_axis].OE_n - 1 && info[elevation_axis].OE_n != 1)
            {
                for (j = 0; j < n + 1; j++)
                {
                    for (k = 0; k < info_glo->DIMENSION; k++)
                    {
                        B_DIM[k].line[j + 1] = B_DIM[k].new_line[j];
                    }
                    B_w.line[j + 1] = B_w.new_line[j];
                }
                n++;
            }
        }

        for (i = 0; i < n + 1; i++)
        {
            for (k = 0; k < info_glo->DIMENSION; k++)
            {
                DIM[k].new_line[*counter] = B_DIM[k].new_line[i];
            }
            w->new_line[*counter] = B_w.new_line[i];
            *counter += 1;
        }

        free(line_x), free(line_y), free(line_z), free(line_w);
        free(new_line_x), free(new_line_y), free(new_line_z), free(new_line_w);
        free(B_DIM), free(B_DIM_ptr);
    }
}


void Calc_Bezier(int elevation_axis, info_global *info_glo, info_each_DIMENSION *info, line_coordinate *w, line_coordinate *DIM)
{
    int i, j;
    int counter = 1;
    int n = info[elevation_axis].Order;
    int l = info[elevation_axis].CP_n;
    int number_of_Bezier_line = (l - 1) / n;
    int *temp_Bezier_array = (int *)malloc(sizeof(int) * number_of_Bezier_line);

    for (i = 0; i < number_of_Bezier_line; i++)
    {
        temp_Bezier_array[i] = i * n;
    }

    for (i = 0; i < info_glo->DIMENSION; i++)
    {
        DIM[i].new_line[0] = DIM[i].line[0];
    }
    w->new_line[0] = w->line[0];

    for (j = 0; j < number_of_Bezier_line; j++)
    {
        Bezier_Order_Elevation(elevation_axis, temp_Bezier_array[j], &counter, info_glo, info, w, DIM);
    }

    free(temp_Bezier_array);
}


// Knot Removal
void KR_non_uniform(int removal_axis, int removal_knot_n, double *removal_knot, info_global *info_glo, info_each_DIMENSION *info)
{
    int i, j, k, l;
    int other_axis = 0;
    int other_axis_1 = 0, other_axis_2 = 0;

    // calc knot
    KR_calc_knot_1D(removal_axis, removal_knot_n, removal_knot, info);

    if (info_glo->DIMENSION == 2)
    {
        if (removal_axis == 0)
        {
            other_axis = 1;
        }
        else if (removal_axis == 1)
        {
            other_axis = 0;
        }

        for (i = 0; i < info[other_axis].CP_n; i++)
        {
            // memory allocation
            line_coordinate w, *w_ptr;
            line_coordinate *DIM = (line_coordinate *)malloc(sizeof(line_coordinate) * info_glo->DIMENSION);
            line_coordinate **DIM_ptr = (line_coordinate **)malloc(sizeof(line_coordinate *) * info_glo->DIMENSION);
            w_ptr = &w;
            for (j = 0; j < info_glo->DIMENSION; j++)
            {
                DIM_ptr[j] = &DIM[j];
            }
            double *line_x = (double *)malloc(sizeof(double) * info[removal_axis].CP_n);
            double *line_y = (double *)malloc(sizeof(double) * info[removal_axis].CP_n);
            double *line_w = (double *)malloc(sizeof(double) * info[removal_axis].CP_n);
            double *new_line_x = (double *)malloc(sizeof(double) * info[removal_axis].CP_n);
            double *new_line_y = (double *)malloc(sizeof(double) * info[removal_axis].CP_n);
            double *new_line_w = (double *)malloc(sizeof(double) * info[removal_axis].CP_n);
            DIM_ptr[0]->line = line_x;
            DIM_ptr[1]->line = line_y;
            w_ptr->line = line_w;
            DIM_ptr[0]->new_line = new_line_x;
            DIM_ptr[1]->new_line = new_line_y;
            w_ptr->new_line = new_line_w;

            // get line info
            for (j = 0; j < info[removal_axis].CP_n; j++)
            {
                int temp = 0;
                if (removal_axis == 0)
                {
                    temp = i * info[0].CP_n + j;
                }
                else if (removal_axis == 1)
                {
                    temp = j * info[0].CP_n + i;
                }

                for (k = 0; k < info_glo->DIMENSION; k++)
                {
                    DIM[k].line[j] = info[k].CP[temp];
                }
                w.line[j] = info_glo->Weight[temp];
            }

            // calc line
            KR_calc_Tinv_1D(removal_axis, removal_knot_n, info_glo, info, &w, DIM);

            // update
            for (j = 0; j < info[removal_axis].CP_n - removal_knot_n; j++)
            {
                int temp = 0;
                if (removal_axis == 0)
                {
                    temp = i * (info[0].CP_n - removal_knot_n) + j;
                }
                else if (removal_axis == 1)
                {
                    temp = j * info[0].CP_n + i;
                }

                for (k = 0; k < info_glo->DIMENSION; k++)
                {
                    info[k].temp_CP[temp] = DIM[k].new_line[j];
                }
                info_glo->temp_Weight[temp] = w.new_line[j];
            }

            // memory free
            free(line_x), free(line_y), free(line_w);
            free(new_line_x), free(new_line_y), free(new_line_w);
            free(DIM), free(DIM_ptr);
        }
    }
    else if(info_glo->DIMENSION == 3)
    {
        if (removal_axis == 0)
        {
            other_axis_1 = 1, other_axis_2 = 2;
        }
        else if (removal_axis == 1)
        {
            other_axis_1 = 0, other_axis_2 = 2;
        }
        else if (removal_axis == 2)
        {
            other_axis_1 = 0, other_axis_2 = 1;
        }

        for (i = 0; i < info[other_axis_1].CP_n; i++)
        {
            for (j = 0; j < info[other_axis_2].CP_n; j++)
            {
                // memory allocation
                line_coordinate w, *w_ptr;
                line_coordinate *DIM = (line_coordinate *)malloc(sizeof(line_coordinate) * info_glo->DIMENSION);
                line_coordinate **DIM_ptr = (line_coordinate **)malloc(sizeof(line_coordinate *) * info_glo->DIMENSION);
                w_ptr = &w;
                for (k = 0; k < info_glo->DIMENSION; k++)
                {
                    DIM_ptr[k] = &DIM[k];
                }
                double *line_x = (double *)malloc(sizeof(double) * info[removal_axis].CP_n);
                double *line_y = (double *)malloc(sizeof(double) * info[removal_axis].CP_n);
                double *line_z = (double *)malloc(sizeof(double) * info[removal_axis].CP_n);
                double *line_w = (double *)malloc(sizeof(double) * info[removal_axis].CP_n);
                double *new_line_x = (double *)malloc(sizeof(double) * info[removal_axis].CP_n);
                double *new_line_y = (double *)malloc(sizeof(double) * info[removal_axis].CP_n);
                double *new_line_z = (double *)malloc(sizeof(double) * info[removal_axis].CP_n);
                double *new_line_w = (double *)malloc(sizeof(double) * info[removal_axis].CP_n);
                DIM_ptr[0]->line = line_x;
                DIM_ptr[1]->line = line_y;
                DIM_ptr[2]->line = line_z;
                w_ptr->line = line_w;
                DIM_ptr[0]->new_line = new_line_x;
                DIM_ptr[1]->new_line = new_line_y;
                DIM_ptr[2]->new_line = new_line_z;
                w_ptr->new_line = new_line_w;

                // get line info
                for (k = 0; k < info[removal_axis].CP_n; k++)
                {
                    int temp = 0;
                    if (removal_axis == 0)
                    {
                        temp = j * (info[1].CP_n * info[0].CP_n) + i * info[0].CP_n + k;
                    }
                    else if (removal_axis == 1)
                    {
                        temp = j * (info[1].CP_n * info[0].CP_n) + k * info[0].CP_n + i;
                    }
                    else if (removal_axis == 2)
                    {
                        temp = k * (info[1].CP_n * info[0].CP_n) + j * info[0].CP_n + i;
                    }

                    for (l = 0; l < info_glo->DIMENSION; l++)
                    {
                        DIM[l].line[k] = info[l].CP[temp];
                    }
                    w.line[k] = info_glo->Weight[temp];
                }

                // calc line
                KR_calc_Tinv_1D(removal_axis, removal_knot_n, info_glo, info, &w, DIM);

                // update
                for (k = 0; k < info[removal_axis].CP_n - removal_knot_n; k++)
                {
                    int temp = 0;
                    if (removal_axis == 0)
                    {
                        temp = j * (info[1].CP_n * (info[0].CP_n - removal_knot_n)) + i * (info[0].CP_n - removal_knot_n) + k;
                    }
                    else if (removal_axis == 1)
                    {
                        temp = j * ((info[1].CP_n - removal_knot_n) * info[0].CP_n) + k * info[0].CP_n + i;
                    }
                    else if (removal_axis == 2)
                    {
                        temp = k * (info[1].CP_n * info[0].CP_n) + j * info[0].CP_n + i;
                    }

                    for (l = 0; l < info_glo->DIMENSION; l++)
                    {
                        info[l].temp_CP[temp] = DIM[l].new_line[k];
                    }
                    info_glo->temp_Weight[temp] = w.new_line[k];
                }

                // memory free
                free(line_x), free(line_y), free(line_z), free(line_w);
                free(new_line_x), free(new_line_y), free(new_line_z), free(new_line_w);
                free(DIM), free(DIM_ptr);
            }
        }
    }

    // update
    for (i = 0; i < info[removal_axis].knot_n - removal_knot_n; i++)
    {
        info[removal_axis].KV[i] = info[removal_axis].temp_KV[i];
    }
    info[removal_axis].CP_n -= removal_knot_n;
    info[removal_axis].knot_n -= removal_knot_n;
    int temp_CP_n = 1;
    for (i = 0; i < info_glo->DIMENSION; i++)
    {
        temp_CP_n *= info[i].CP_n;
    }
    info_glo->Total_Control_Point = temp_CP_n;
    for (i = 0; i < info_glo->Total_Control_Point; i++)
    {
        for (j = 0; j < info_glo->DIMENSION; j++)
        {
            info[j].CP[i] = info[j].temp_CP[i];
        }
        info_glo->Weight[i] = info_glo->temp_Weight[i];
    }
}


void KR_calc_knot_1D(int removal_axis, int removal_knot_n, double *removal_knot, info_each_DIMENSION *info)
{
    int i, j;

    int temp1 = 0, temp2 = 0, temp3 = 0;
    for (i = 0; i < info[removal_axis].knot_n; i++)
    {
        if (temp2 < removal_knot_n)
        {
            if (info[removal_axis].KV[temp1] != removal_knot[temp2])
            {
                info[removal_axis].temp_KV[temp3] = info[removal_axis].KV[temp1];
                temp1++;
                temp3++;
            }
            else if (info[removal_axis].KV[temp1] == removal_knot[temp2])
            {
                temp1++;
                temp2++;
            }
        }
        else if (temp2 == removal_knot_n)
        {
            for (j = 0; j < info[removal_axis].knot_n - temp1; j++)
            {
                info[removal_axis].temp_KV[temp3] = info[removal_axis].KV[temp1];
                temp1++;
                temp3++;
            }
        }
    }
}


void KR_calc_Tinv_1D(int removal_axis, int removal_knot_n, info_global *info_glo, info_each_DIMENSION *info, line_coordinate *w, line_coordinate *DIM)
{
    int i, j, k;

    int l = info[removal_axis].CP_n - removal_knot_n;
    int n = info[removal_axis].Order;

    // pre processing
    for (i = 0; i < l + removal_knot_n; i++)
    {
        for (j = 0; j < info_glo->DIMENSION; j++)
        {
            DIM[j].line[i] = DIM[j].line[i] * w->line[i];
        }
    }

    double *T = (double *)calloc((n + 1) * (l + removal_knot_n) * l, sizeof(double));

    double a = 0.0, b = 0.0;
    for (i = 0; i < l + removal_knot_n; i++)
    {
        for (k = 0; k < n + 1; k++)
        {
            for (j = 0; j < l; j++)
            {
                if (k == 0)
                {
                    if (info[removal_axis].temp_KV[j] <= info[removal_axis].KV[i] && info[removal_axis].KV[i] < info[removal_axis].temp_KV[j + 1])
                    {
                        T[k * (l + removal_knot_n) * l + i * l + j] = 1.0;
                    }
                    else
                    {
                        T[k * (l + removal_knot_n) * l + i * l + j] = 0.0;
                    }
                }
                else
                {
                    if (info[removal_axis].temp_KV[j + k] - info[removal_axis].temp_KV[j] == 0)
                    {
                        a = 0.0;
                    }
                    else if (info[removal_axis].temp_KV[j + k] - info[removal_axis].temp_KV[j] != 0)
                    {
                        a = ((info[removal_axis].KV[i + k] - info[removal_axis].temp_KV[j]) / (info[removal_axis].temp_KV[j + k] - info[removal_axis].temp_KV[j])) * T[(k - 1) * (l + removal_knot_n) * l + i * l + j];
                    }
                    if (info[removal_axis].temp_KV[j + k + 1] - info[removal_axis].temp_KV[j + 1] == 0)
                    {
                        b = 0.0;
                    }
                    else if (info[removal_axis].temp_KV[j + k + 1] - info[removal_axis].temp_KV[j + 1] != 0)
                    {
                        b = ((info[removal_axis].temp_KV[j + k + 1] - info[removal_axis].KV[i + k]) / (info[removal_axis].temp_KV[j + k + 1] - info[removal_axis].temp_KV[j + 1])) * T[(k - 1) * (l + removal_knot_n) * l + i * l + j + 1];
                    }
                    T[k * (l + removal_knot_n) * l + i * l + j] = a + b;
                }
            }
        }
    }

    // 疑似逆行列の計算

    double *A = (double *)malloc(sizeof(double) * (l + removal_knot_n) * l);
    double *A_T = (double *)malloc(sizeof(double) * l * (l + removal_knot_n));
    double *B = (double *)malloc(sizeof(double) * l * l);
    double *B_inv = (double *)malloc(sizeof(double) * l * l);
    double *A_pinv = (double *)malloc(sizeof(double) * (l + removal_knot_n) * l);
    double *A_pinv_T = (double *)malloc(sizeof(double) * l * (l + removal_knot_n));

    for (i = 0; i < l + removal_knot_n; i++)
    {
        for (j = 0; j < l; j++)
        {
            A[i * l + j] = T[n * (l + removal_knot_n) * l + i * l + j];
        }
    }

    for (i = 0; i < l + removal_knot_n; i++)
    {
        for (j = 0; j < l; j++)
        {
            A_T[j * (l + removal_knot_n) + i] = A[i * l + j];
        }
    }

    for (i = 0; i < l; i++)
    {
        for (j = 0; j < l; j++)
        {
            B[i * l + j] = 0.0;
            for (k = 0; k < l + removal_knot_n; k++)
            {
                B[i * l + j] += A_T[i * (l + removal_knot_n) + k] * A[k * l + j];
            }
        }
    }

    double temp;

    for (i = 0; i < l; i++)
    {
        for (j = 0; j < l; j++)
        {
            B_inv[i * l + j] = 0.0;
        }
    }

    for (i = 0; i < l; i++)
    {
        B_inv[i * l + i] = 1.0;
    }

    for (k = 0; k < l; k++)
    {
        temp = B[k * l + k];
        for (i = 0; i < l; i++)
        {
            B[k * l + i] /= temp;
            B_inv[k * l + i] /= temp;
        }

        for (i = 0; i < l; i++)
        {
            if (i != k)
            {
                temp = B[i * l + k];
                for (j = 0; j < l; j++)
                {
                    B[i * l + j] -= B[k * l + j] * temp;
                    B_inv[i * l + j] -= B_inv[k * l + j] * temp;
                }
            }
        }
    }

    for (i = 0; i < l + removal_knot_n; i++)
    {
        for (j = 0; j < l; j++)
        {
            A_pinv[i * l + j] = 0.0;
            for (k = 0; k < l; k++)
            {
                A_pinv[i * l + j] += A[i * l + k] * B_inv[k * l + j];
            }
        }
    }

    for (i = 0; i < l + removal_knot_n; i++)
    {
        for (j = 0; j < l; j++)
        {
            A_pinv_T[j * (l + removal_knot_n) + i] = A_pinv[i * l + j];
        }
    }

    // inner product
    for (i = 0; i < l; i++)
    {
        for (j = 0; j < l + removal_knot_n; j++)
        {
            if (j == 0)
            {
                for (k = 0; k < info_glo->DIMENSION; k++)
                {
                    DIM[k].new_line[i] = 0.0;
                }
                w->new_line[i] = 0.0;
            }

            for (k = 0; k < info_glo->DIMENSION; k++)
            {
                DIM[k].new_line[i] += A_pinv_T[i * (l + removal_knot_n) + j] * DIM[k].line[j];
            }
            w->new_line[i] += A_pinv_T[i * (l + removal_knot_n) + j] * w->line[j];
        }
    }

    // post processing
    for (i = 0; i < l; i++)
    {
        for (j = 0; j < info_glo->DIMENSION; j++)
        {
            DIM[j].new_line[i] = DIM[j].new_line[i] / w->new_line[i];
        }
    }

    free(A), free(A_T), free(B), free(B_inv), free(A_pinv), free(A_pinv_T);
    free(T);
}


// Output
void OutputData(char *filename, info_global *info_glo, info_each_DIMENSION *info)
{
    int a = strlen(filename);
    char str[256];

    int i, j;
    int temp1 = 0;
    for (i = 0; i < a - 4; i++)
    {
        str[i] = filename[i];
        temp1++;
    }
    str[temp1] = '_';
    temp1++;
    str[temp1] = 'e';
    temp1++;
    str[temp1] = 'd';
    temp1++;
    str[temp1] = 'i';
    temp1++;
    str[temp1] = 't';
    temp1++;
    str[temp1] = 'e';
    temp1++;
    str[temp1] = 'd';
    temp1++;
    str[temp1] = '.';
    temp1++;
    str[temp1] = 't';
    temp1++;
    str[temp1] = 'x';
    temp1++;
    str[temp1] = 't';
    temp1++;
    str[temp1] = '\0';
    temp1++;

    fp = fopen(str, "w");

    fprintf(fp, "%d\n\n", info_glo->Total_Control_Point);

    for (j = 0; j < info_glo->DIMENSION; j++)
    {
        if (j == info_glo->DIMENSION - 1)
        {
            fprintf(fp, "%d", info[j].Order);
        }
        else
        {
            fprintf(fp, "%-5d", info[j].Order);
        }
    }
    fprintf(fp, "\n\n");

    for (j = 0; j < info_glo->DIMENSION; j++)
    {
        if (j == info_glo->DIMENSION - 1)
        {
            fprintf(fp, "%d", info[j].knot_n);
        }
        else
        {
            fprintf(fp, "%-5d", info[j].knot_n);
        }
    }
    fprintf(fp, "\n\n");

    for (j = 0; j < info_glo->DIMENSION; j++)
    {
        if (j == info_glo->DIMENSION - 1)
        {
            fprintf(fp, "%d", info[j].CP_n);
        }
        else
        {
            fprintf(fp, "%-5d", info[j].CP_n);
        }
    }
    fprintf(fp, "\n\n");

    for (j = 0; j < info_glo->DIMENSION; j++)
    {
        for (i = 0; i < info[j].knot_n; i++)
        {
            if (i == info[j].knot_n - 1)
            {
                fprintf(fp, "%.16e", info[j].KV[i]);
            }
            else
            {
                fprintf(fp, "%.16e  ", info[j].KV[i]);
            }
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    temp1 = info_glo->Total_Control_Point;
    int digit = 0;
    while(temp1 != 0)
    {
        temp1 = temp1 / 10;
        ++digit;
    }

    for (i = 0; i < info_glo->Total_Control_Point; i++)
    {
        fprintf(fp, "%*d", (int)(-1 * (digit + 2)), i);
        for (j = 0; j < info_glo->DIMENSION; j++)
        {
            fprintf(fp, "  %.16e", info[j].CP[i]);
        }
        fprintf(fp, "  %.16e\n", info_glo->Weight[i]);
    }

    fclose(fp);
}