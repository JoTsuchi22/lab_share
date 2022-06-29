#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAX_N_INPUTFILE 5
#define MAX_DIMENSION 3
#define MAX_ORDER 3
#define MAX_N_Controlpoint_each_parameter 1000
#define MAX_N_Controlpoint_in_Patch MAX_N_Controlpoint_each_parameter * MAX_N_Controlpoint_each_parameter
#define MAX_N_KNOT MAX_N_Controlpoint_each_parameter + MAX_ORDER + 1

static int Dimension[MAX_N_INPUTFILE];
static int Total_Control_Point[MAX_N_INPUTFILE];
static int Order[MAX_N_INPUTFILE][MAX_DIMENSION];
static int Order_before[MAX_N_INPUTFILE][MAX_DIMENSION];
static int knot_n[MAX_N_INPUTFILE][MAX_DIMENSION];
static int Control_point_n[MAX_N_INPUTFILE][MAX_DIMENSION];
static int Control_point_n_before[MAX_N_INPUTFILE][MAX_DIMENSION];
static double knot[MAX_N_INPUTFILE][MAX_DIMENSION][MAX_N_KNOT];
static double x[MAX_N_INPUTFILE][MAX_N_Controlpoint_in_Patch];
static double y[MAX_N_INPUTFILE][MAX_N_Controlpoint_in_Patch];
static double w[MAX_N_INPUTFILE][MAX_N_Controlpoint_in_Patch];
static int OE_n[MAX_N_INPUTFILE][MAX_DIMENSION];
static int KI_uniform_n[MAX_N_INPUTFILE][MAX_DIMENSION];
static int KI_cp_n[MAX_N_INPUTFILE][MAX_DIMENSION];
static int KI_non_uniform_n[MAX_N_INPUTFILE][MAX_DIMENSION];
static double insert_knot[MAX_N_INPUTFILE][MAX_DIMENSION][MAX_N_KNOT];
static double insert_knot_in_KI[MAX_N_INPUTFILE][MAX_DIMENSION][MAX_N_KNOT];
static double removal_knot[MAX_N_INPUTFILE][MAX_DIMENSION][MAX_N_KNOT];
static int vec_length1[MAX_N_INPUTFILE][MAX_DIMENSION];
static int vec_length2[MAX_N_INPUTFILE][MAX_DIMENSION];
static double temp_knot1[MAX_N_KNOT];
static double temp_knot2[MAX_N_KNOT];
static double x_array[MAX_N_Controlpoint_each_parameter][MAX_N_Controlpoint_each_parameter];
static double y_array[MAX_N_Controlpoint_each_parameter][MAX_N_Controlpoint_each_parameter];
static double w_array[MAX_N_Controlpoint_each_parameter][MAX_N_Controlpoint_each_parameter];
static double temp_x_array[MAX_N_Controlpoint_each_parameter];
static double temp_y_array[MAX_N_Controlpoint_each_parameter];
static double temp_w_array[MAX_N_Controlpoint_each_parameter];
static double temp_x1[MAX_N_Controlpoint_in_Patch];
static double temp_y1[MAX_N_Controlpoint_in_Patch];
static double temp_w1[MAX_N_Controlpoint_in_Patch];
static double temp_x2[MAX_N_Controlpoint_in_Patch];
static double temp_y2[MAX_N_Controlpoint_in_Patch];
static double temp_w2[MAX_N_Controlpoint_in_Patch];
static int temp_Bezier_array[MAX_N_Controlpoint_each_parameter];
static int counter;
static double Bezier_x[MAX_N_Controlpoint_each_parameter][MAX_ORDER];
static double Bezier_y[MAX_N_Controlpoint_each_parameter][MAX_ORDER];
static double Bezier_w[MAX_N_Controlpoint_each_parameter][MAX_ORDER];
static int number_of_Bezier_line;
static int error_check = 0;
static int mode;

void Get_InputData(int tm, char *filename);
void Print_error();
void KI_calc_point_array(int tm);
void KI_calc_knot_1D(int tm, int insert_parameter_axis);
void KI_define_temp_point_array(int tm, int line_number, int insert_parameter_axis);
void KI_calc_T_1D(int tm, int insert_parameter_axis);
void KI_update_point_array(int tm, int line_number, int insert_parameter_axis);
void KI_update(int tm, int insert_parameter_axis);
void KI_reset_array();
void KI_non_uniform_2D(int tm, int insert_parameter_axis, int KI_non_uniform);
void Calc_cp_insert_knot(int tm, int insert_parameter_axis);
void KI_cp_2D(int tm, int insert_parameter_axis);
void Calc_uniform_insert_knot(int tm, int insert_parameter_axis);
void KI_uniform_2D(int tm, int insert_parameter_axis);
void OE_calc_point_array(int tm);
void Calc_insert_knot_in_OE(int tm, int elevation_parameter_axis);
void OE_define_temp_point_array(int tm, int line_number, int elevation_parameter_axis);
void Bezier_Order_Elevation_2D(int tm, int elevation_parameter_axis, int Bezier_line_number);
void Bezier_update_point_array(int tm, int line_number, int elevation_parameter_axis);
void Bezier_update(int tm, int elevation_parameter_axis);
void Calc_Bezier_2D(int tm, int elevation_parameter_axis, int other_axis);
void KR_calc_point_array(int tm);
void KR_calc_knot_1D(int tm, int removal_parameter_axis);
void KR_define_temp_point_array(int tm, int line_number, int removal_parameter_axis);
void KR_calc_Tinv_1D(int tm, int removal_parameter_axis);
void KR_update_point_array(int tm, int line_number, int removal_parameter_axis);
void KR_update(int tm, int removal_parameter_axis);
void KR_reset_array();
void KR_non_uniform_2D(int tm, int removal_parameter_axis);
void OE_2D(int tm, int elevation_parameter_axis);
void Debug_printf(int tm, char *section);
void Fix_numerical_error(int tm);
void OutputData(int tm, char *filename);

FILE *fp;

int  main(int argc, char *argv[])
{
    int tm;
    int Total_file;
    int j;

    Total_file = argc - 1;

    for (tm = 0; tm < Total_file; tm++)
    {
        //ファイル読み込み
	    Get_InputData(tm, argv[tm+1]);
        if (error_check == 1)
        {
            Print_error();
            break;
        }
	    Debug_printf(tm, "Get_InputData");

        //オーダーエレベーション
        for (j = 0; j < Dimension[tm]; j++)
        {
            OE_2D(tm, j);
            Debug_printf(tm, "Order Elevation");
        }
        
        //ノットインサーション
        for (j = 0; j < Dimension[tm]; j++)
        {
            if (mode == 0)
            {
                KI_non_uniform_2D(tm, j, 1);
            }
            else if (mode == 1)
            {
                KI_cp_2D(tm, j);
            }
            else if (mode == 2)
            {
                KI_uniform_2D(tm, j);
            }
            Debug_printf(tm, "Knot Insertion");
        }

        //逆行列での数値誤差を修正
        //Fix_numerical_error(tm);

        //結果を出力
        OutputData(tm, argv[tm+1]);
    }
    return 0;
}


void Get_InputData(int tm, char *filename)
{
	int i, j, k;
	char s[256];

    int temp_dim;
    int temp_Total_Control_Point;
    int temp_Order;
    int temp_knot_n;
    int temp_Control_point_n;
    double temp_knot;
    int temp_num;
    double temp_x, temp_y, temp_w;
    int temp_OE_n;
    int temp_KI_uniform_n;
    int temp_KI_cp_n;
    int temp_KI_non_uniform_n;
    double temp_insert_knot;

	fp = fopen(filename, "r");

    //DIMENSION
    fscanf(fp, "%d", &temp_dim);
	fgets(s, 256, fp);
    printf("DIMENSION:%d\n", temp_dim);
    Dimension[tm] = temp_dim;
    if (Dimension[tm] != 2)
    {
        printf("Error, Dimension must be 2 in this program.\n");
    }

	//コントロールポイントの総数
	fscanf(fp, "%d", &temp_Total_Control_Point);
	fgets(s, 256, fp);
	printf("Total_Control_Point:%d\n", temp_Total_Control_Point);
    Total_Control_Point[tm] = temp_Total_Control_Point;
    

	//ξη方向の次数
    for (j = 0; j < Dimension[tm]; j++)
    {
        fscanf(fp, "%d", &temp_Order);
        printf("Order[%d][%d]=%d\n", tm, j, temp_Order);
        Order[tm][j] = temp_Order;
        Order_before[tm][j] = temp_Order;
    }
	fgets(s, 256, fp);

	//ξη方向のノットベクトルの個数
    for (j = 0; j < Dimension[tm]; j++)
    {
        fscanf(fp, "%d", &temp_knot_n);
        printf("knot_n[%d][%d]=%d\n", tm, j, temp_knot_n);
        knot_n[tm][j] = temp_knot_n;
    }
	fgets(s, 256, fp);

	//ξη方向のコントロールポイントの個数
    for (j = 0; j < Dimension[tm]; j++)
    {
        fscanf(fp, "%d", &temp_Control_point_n);
        printf("Control_point_n[%d][%d]:%d\n", tm, j, temp_Control_point_n);
        Control_point_n[tm][j] = temp_Control_point_n;
        Control_point_n_before[tm][j] = temp_Control_point_n;
    }
	fgets(s, 256, fp);

	//ξη方向のノットベクトル
    for (j = 0; j < Dimension[tm]; j++)
    {
        for (k = 0; k < knot_n[tm][j]; k++)
        {
            fscanf(fp, "%lf", &temp_knot);
            printf("%.16e\t", temp_knot);
            knot[tm][j][k] = temp_knot;
        }
        printf("\n");
    }
    fgets(s, 256, fp);

    //コントロールポイントの座標
    for (i = 0; i < Total_Control_Point[tm]; i++)
    {
        fscanf(fp, "%d", &temp_num);
        fscanf(fp, "%lf", &temp_x);
        fscanf(fp, "%lf", &temp_y);
        fscanf(fp, "%lf", &temp_w);
        printf("%d\t", temp_num);
        printf("%.16e\t", temp_x);
        printf("%.16e\t", temp_y);
        printf("%.16e\n", temp_w);
        x[tm][i] = temp_x;
        y[tm][i] = temp_y;
        w[tm][i] = temp_w;
    }
    fgets(s, 256, fp);

    //ξη方向のo.e.の回数
    for (j = 0; j < Dimension[tm]; j++)
    {
        fscanf(fp, "%d", &temp_OE_n);
        printf("OE_n[%d][%d]:%d\n", tm, j, temp_OE_n);
        OE_n[tm][j] = temp_OE_n;
    }
	fgets(s, 256, fp);

    //ξη方向のk.i.(等分割)の回数
    for (j = 0; j < Dimension[tm]; j++)
    {
        fscanf(fp, "%d", &temp_KI_uniform_n);
        printf("KI_uniform_n[%d][%d]:%d\n", tm, j, temp_KI_uniform_n);
        KI_uniform_n[tm][j] = temp_KI_uniform_n;
    }
	fgets(s, 256, fp);

    //ξη方向のk.i.(コントロールポイントに合わせて等分割)の回数
    for (j = 0; j < Dimension[tm]; j++)
    {
        fscanf(fp, "%d", &temp_KI_cp_n);
        printf("KI_cp_n[%d][%d]:%d\n", tm, j, temp_KI_cp_n);
        KI_cp_n[tm][j] = temp_KI_cp_n;
    }
    fgets(s, 256, fp);

    //ξη方向のk.i.(任意のノット)の回数
    for (j = 0; j < Dimension[tm]; j++)
    {
        fscanf(fp, "%d", &temp_KI_non_uniform_n);
        printf("KI_non_uniform_n[%d][%d]:%d\n", tm, j, temp_KI_non_uniform_n);
        KI_non_uniform_n[tm][j] = temp_KI_non_uniform_n;
    }
    fgets(s, 256, fp);

    //ξη方向の挿入するノット
    for (j = 0; j < Dimension[tm]; j++)
    {
        if (KI_non_uniform_n[tm][j] != 0)
        {
            for (k = 0; k < KI_non_uniform_n[tm][j]; k++)
            {
                fscanf(fp, "%lf", &temp_insert_knot);
                printf("%.16e\t", temp_insert_knot);
                insert_knot[tm][j][k] = temp_insert_knot;
            }
            printf("\n");
        }
    }

    //コントロールポイントに合わせて分割と任意のノット挿入を同時に行っていないか確認
    int temp_check_a = 0;
    int temp_check_b = 0;
    int temp_check_c = 0;
    int temp_check_total = 0;
    for (j = 0; j < Dimension[tm]; j++)
    {
        if(KI_non_uniform_n[tm][j] != 0)
        {
            temp_check_a += 1;
        }
        if(KI_cp_n[tm][j] != 0)
        {
            temp_check_b += 1;
        }
        if(KI_uniform_n[tm][j] != 0)
        {
            temp_check_c += 1;
        }
    }
    if (temp_check_a > 0)
    {
        mode = 0;
        temp_check_total++;
    }
    if (temp_check_b > 0)
    {
        mode = 1;
        temp_check_total++;
    }
    if (temp_check_c > 0)
    {
        mode = 2;
        temp_check_total++;
    }

    if (temp_check_total > 1)
    {
        error_check = 1;
    }
    else
    {
        error_check = 0;
    }

    fclose(fp);
}


void Print_error()
{
    printf("インプットデータが間違っています．\nノットインサーションの各操作(3種類)は同時に行えません．\n操作を行わない該当するパラメータを(0  0)としてください．\n");
}


void KI_calc_point_array(int tm)
{
    int i, j;
    int temp1 = 0, temp2 = 0;
    for (i = 0; i < Control_point_n[tm][1]; i++)
    {
        for (j = 0; j < Control_point_n[tm][0]; j++)
        {
            x_array[temp1][temp2] = x[tm][temp1 + (Control_point_n[tm][0] * temp2)];
            y_array[temp1][temp2] = y[tm][temp1 + (Control_point_n[tm][0] * temp2)];
            w_array[temp1][temp2] = w[tm][temp1 + (Control_point_n[tm][0] * temp2)];
            temp1++;
        }
        temp1 = 0;
        temp2++;
    }
}


void KI_calc_knot_1D(int tm, int insert_parameter_axis)
{
    int i;
    for (i = 0; i < MAX_N_KNOT; i++)
    {
        temp_knot1[i] = 0.0;
        temp_knot2[i] = 0.0;
    }

    for (i = 0; i < knot_n[tm][insert_parameter_axis]; i++)
    {
        temp_knot1[i] = knot[tm][insert_parameter_axis][i];
    }

    int temp1 = 0, temp2 = 0, temp3 = 0;
    for (i = 0; i < knot_n[tm][insert_parameter_axis] + vec_length1[tm][insert_parameter_axis]; i++)
    {
        if (temp2 < vec_length1[tm][insert_parameter_axis])
        {
            if (temp_knot1[temp1] < insert_knot_in_KI[tm][insert_parameter_axis][temp2]) 
            {
                temp_knot2[temp3] = temp_knot1[temp1];
                temp1++;
                temp3++;
            }
            else if (temp_knot1[temp1] > insert_knot_in_KI[tm][insert_parameter_axis][temp2]) 
            {
                temp_knot2[temp3] = insert_knot_in_KI[tm][insert_parameter_axis][temp2];
                temp2++;
                temp3++;
            }
            else if (temp_knot1[temp1] == insert_knot_in_KI[tm][insert_parameter_axis][temp2]) 
            {
                temp_knot2[temp3] = temp_knot1[temp1];
                temp_knot2[temp3 + 1] = insert_knot_in_KI[tm][insert_parameter_axis][temp2];
                temp1++;
                temp2++;
                temp3++;
                temp3++;
            }
        }
        else
        {
            temp_knot2[temp3] = temp_knot1[temp1];
            temp1++;
            temp3++;
        }
    }
}


void KI_define_temp_point_array(int tm, int line_number, int insert_parameter_axis)
{
    int i;
    if (insert_parameter_axis == 0)
    {
        for (i = 0; i < Control_point_n[tm][insert_parameter_axis]; i++)
        {
            temp_x_array[i] = x_array[i][line_number];
            temp_y_array[i] = y_array[i][line_number];
            temp_w_array[i] = w_array[i][line_number];
        }
    }
    else if (insert_parameter_axis == 1)
    {
        for (i = 0; i < Control_point_n[tm][insert_parameter_axis]; i++)
        {
            temp_x_array[i] = x_array[line_number][i];
            temp_y_array[i] = y_array[line_number][i];
            temp_w_array[i] = w_array[line_number][i];
        }
    }
}


void KI_calc_T_1D(int tm, int insert_parameter_axis)
{
    int i, j, k;
    int l, n;
    
    l = Control_point_n[tm][insert_parameter_axis];
    n = Order[tm][insert_parameter_axis];

    for (i = 0; i < l; i++)
    {
        temp_x1[i] = temp_x_array[i] * temp_w_array[i];
        temp_y1[i] = temp_y_array[i] * temp_w_array[i];
        temp_w1[i] = temp_w_array[i];
    }

    double T[n + 1][l + vec_length1[tm][insert_parameter_axis]][l];

    for (int i = 0; i < n + 1; i++)
    {
        for (int j = 0; j < l + vec_length1[tm][insert_parameter_axis]; j++)
        {
            for (int k = 0; k < l; k++)
            {
                T[i][j][k] = 0.0;
            }
        }
    }

    double a = 0.0, b = 0.0;
    for (i = 0; i < l + vec_length1[tm][insert_parameter_axis]; i++)
    {
        for (k = 0; k < n + 1; k++)
        {
            for (j = 0; j < l; j++)
            {
                if (k == 0)
                {
                    if (temp_knot1[j] <= temp_knot2[i] && temp_knot2[i] < temp_knot1[j + 1])
                    {
                        T[k][i][j] = 1.0;
                    }
                    else
                    {
                        T[k][i][j] = 0.0;
                    }
                }
                else
                {
                    if (temp_knot1[j + k] - temp_knot1[j] == 0)
                    {
                        a = 0.0;
                    }
                    else if (temp_knot1[j + k] - temp_knot1[j] != 0)
                    {
                        a = ((temp_knot2[i + k] - temp_knot1[j]) / (temp_knot1[j + k] - temp_knot1[j])) * T[k - 1][i][j];
                    }
                    if (temp_knot1[j + k + 1] - temp_knot1[j + 1] == 0)
                    {
                        b = 0.0;
                    }
                    else if (temp_knot1[j + k + 1] - temp_knot1[j + 1] != 0)
                    {
                        b = ((temp_knot1[j + k + 1] - temp_knot2[i + k]) / (temp_knot1[j + k + 1] - temp_knot1[j + 1])) * T[k - 1][i][j + 1];
                    }
                    T[k][i][j] = a + b;
                }
            }
        }
    }

    //内積の計算
    for (i = 0; i < l + vec_length1[tm][insert_parameter_axis]; i++)
    {
        for (j = 0; j < l; j++)
        {
            if (j == 0)
            {
                temp_x2[i] = 0.0;
                temp_y2[i] = 0.0;
                temp_w2[i] = 0.0;
            }
            temp_x2[i] += T[n][i][j] * temp_x1[j];
            temp_y2[i] += T[n][i][j] * temp_y1[j];
            temp_w2[i] += T[n][i][j] * temp_w1[j];
        }
    }

    for (i = 0; i < l + vec_length1[tm][insert_parameter_axis]; i++)
    {
        temp_x_array[i] = temp_x2[i] / temp_w2[i];
        temp_y_array[i] = temp_y2[i] / temp_w2[i];
        temp_w_array[i] = temp_w2[i];
    }
}


void KI_update_point_array(int tm, int line_number, int insert_parameter_axis)
{
    int i;
    if (insert_parameter_axis == 0)
    {
        for (i = 0; i < Control_point_n[tm][insert_parameter_axis] + vec_length1[tm][insert_parameter_axis]; i++)
        {
            x_array[i][line_number] = temp_x_array[i];
            y_array[i][line_number] = temp_y_array[i];
            w_array[i][line_number] = temp_w_array[i];
        }
    }
    else if (insert_parameter_axis == 1)
    {
        for (i = 0; i < Control_point_n[tm][insert_parameter_axis] + vec_length1[tm][insert_parameter_axis]; i++)
        {
            x_array[line_number][i] = temp_x_array[i];
            y_array[line_number][i] = temp_y_array[i];
            w_array[line_number][i] = temp_w_array[i];
        }
    }
}


void KI_update(int tm, int insert_parameter_axis)
{
    int i, j;
    for (i = 0; i < knot_n[tm][insert_parameter_axis] + vec_length1[tm][insert_parameter_axis]; i++)
    {
        knot[tm][insert_parameter_axis][i] = temp_knot2[i];
    }
    Control_point_n[tm][insert_parameter_axis] = Control_point_n[tm][insert_parameter_axis] + vec_length1[tm][insert_parameter_axis];
    Total_Control_Point[tm] = Control_point_n[tm][0] * Control_point_n[tm][1];
    knot_n[tm][insert_parameter_axis] = Control_point_n[tm][insert_parameter_axis] + Order[tm][insert_parameter_axis] + 1;

    int temp1 = 0, temp2 = 0;
    for (i = 0; i < Control_point_n[tm][1]; i++)
    {
        for (j = 0; j < Control_point_n[tm][0]; j++)
        {
            x[tm][temp1 + (Control_point_n[tm][0] * temp2)] = x_array[temp1][temp2];
            y[tm][temp1 + (Control_point_n[tm][0] * temp2)] = y_array[temp1][temp2];
            w[tm][temp1 + (Control_point_n[tm][0] * temp2)] = w_array[temp1][temp2];
            temp1++;
        }
        temp1 = 0;
        temp2++;
    }
}


void KI_reset_array()
{
    int i, j;
    for (i = 0; i < MAX_N_Controlpoint_each_parameter; i++)
    {
        for (j = 0; j < MAX_N_Controlpoint_each_parameter; j++)
        {
            x_array[i][j] = 0.0;
            y_array[i][j] = 0.0;
            w_array[i][j] = 0.0;
        }
    }
    
    for (i = 0; i < MAX_N_Controlpoint_each_parameter; i++)
    {
        temp_x_array[i] = 0.0;
        temp_y_array[i] = 0.0;
        temp_w_array[i] = 0.0;
    }

    for (i = 0; i < MAX_N_Controlpoint_in_Patch; i++)
    {
        temp_x1[i] = 0.0;
        temp_y1[i] = 0.0;
        temp_w1[i] = 0.0;
        temp_x2[i] = 0.0;
        temp_y2[i] = 0.0;
        temp_w2[i] = 0.0;
    }
}


void KI_non_uniform_2D(int tm, int insert_parameter_axis, int KI_non_uniform)
{
    int i;
    int other_axis = 0;
    if (insert_parameter_axis == 0)
    {
        other_axis = 1;
    }
    else if (insert_parameter_axis == 1)
    {
        other_axis = 0;
    }

    if (KI_non_uniform == 1)
    {
        vec_length1[tm][insert_parameter_axis] = KI_non_uniform_n[tm][insert_parameter_axis];
        for (i = 0; i < vec_length1[tm][insert_parameter_axis]; i++)
        {
            insert_knot_in_KI[tm][insert_parameter_axis][i] = insert_knot[tm][insert_parameter_axis][i];
        }
    }

    KI_calc_point_array(tm);

    for (i = 0; i < Control_point_n[tm][other_axis]; i++)
    {
        if (i == 0)
        {
            KI_calc_knot_1D(tm, insert_parameter_axis);
        }

        KI_define_temp_point_array(tm, i, insert_parameter_axis);
        KI_calc_T_1D(tm, insert_parameter_axis);
        KI_update_point_array(tm, i, insert_parameter_axis);
    }

    KI_update(tm, insert_parameter_axis);
    KI_reset_array();
}


void Calc_cp_insert_knot(int tm, int insert_parameter_axis)
{
    int i;
    int cp_n_after_OE = Control_point_n_before[tm][insert_parameter_axis] + OE_n[tm][insert_parameter_axis] * (Control_point_n_before[tm][insert_parameter_axis] - Order_before[tm][insert_parameter_axis]);
    vec_length1[tm][insert_parameter_axis] = KI_cp_n[tm][insert_parameter_axis] - cp_n_after_OE;

    for (i = 0; i < vec_length1[tm][insert_parameter_axis]; i++)
    {
        insert_knot_in_KI[tm][insert_parameter_axis][i] = (i + 1.0) / (vec_length1[tm][insert_parameter_axis] + 1.0);
    }
}


void KI_cp_2D(int tm, int insert_parameter_axis)
{
    Calc_cp_insert_knot(tm, insert_parameter_axis);
    KI_non_uniform_2D(tm, insert_parameter_axis, 0);
}


void Calc_uniform_insert_knot(int tm, int insert_parameter_axis)
{
    int i;
    int temp1 = 0;
    for (i = 0; i < MAX_N_KNOT; i++)
    {
        temp_knot1[i] = 0.0;
        temp_knot2[i] = 0.0;
    }
    for (i = 0; i < knot_n[tm][insert_parameter_axis]; i++)
    {
        if (i == 0)
        {
            temp_knot1[temp1] = knot[tm][insert_parameter_axis][i];
            temp1++;
        }
        if (temp_knot1[temp1 - 1] < knot[tm][insert_parameter_axis][i])
        {
            temp_knot1[temp1] = knot[tm][insert_parameter_axis][i];
            temp1++;
        }
        if (temp_knot1[temp1] == 1.0)
        {
            break;
        }
    }

    int temp2 = 0;
    for (i = 0; i < temp1 - 1; i++)
    {
        temp_knot2[temp2] = temp_knot1[i] + ((temp_knot1[i + 1] - temp_knot1[i]) / 2.0);
        temp2++;
    }

    vec_length1[tm][insert_parameter_axis] = temp2;

    for (i = 0; i < temp2; i++)
    {
        insert_knot_in_KI[tm][insert_parameter_axis][i] = temp_knot2[i];
    }
}


void KI_uniform_2D(int tm, int insert_parameter_axis)
{
    int i;
    for (i = 0; i < KI_uniform_n[tm][insert_parameter_axis]; i++)
    {
        Calc_uniform_insert_knot(tm, insert_parameter_axis);
        KI_non_uniform_2D(tm, insert_parameter_axis, 0);
    }
}


void Calc_insert_knot_in_OE(int tm, int elevation_parameter_axis)
{
    int m, n;

    m = knot_n[tm][elevation_parameter_axis];
    n = Order[tm][elevation_parameter_axis];

    double vec_temp1[MAX_N_KNOT];
    double vec_temp2[MAX_N_KNOT];
    double vec_unique[MAX_N_KNOT];

    int i, j;
    int temp1 = 0;
    for (i = 1; i < m - 1; i++)
    {
        if (knot[tm][elevation_parameter_axis][i - 1] != knot[tm][elevation_parameter_axis][i] && knot[tm][elevation_parameter_axis][i] != knot[tm][elevation_parameter_axis][i + 1])
        {
            vec_unique[temp1] = knot[tm][elevation_parameter_axis][i];
            temp1++;
        }
    }

    int temp3 = 0;
    for (i = 0; i < temp1; i++)
    {
        int temp2 = 0;
        for (j = 0; j < m; j++)
        {
            if (vec_unique[i] == knot[tm][elevation_parameter_axis][j])
            {
                temp2++;
            }
        }
        if (temp2 < n)
        {
            vec_temp1[temp3] = vec_unique[i];
            temp3++;
        }
    }

    int temp4 = 0;
    if (n - 2 > 0)
    { 
        for (i = 0; i < temp3; i++)
        {
            for (j = 0; j < n - 1; j++)
            {
                vec_temp2[temp4] = vec_temp1[j];
                temp4++;
            }
        }
    }

    if (n - 2 > 0)
    {
        for (i = 0; i < temp4; i++)
        {
            insert_knot_in_KI[tm][elevation_parameter_axis][i] = vec_temp2[i];
        }
        vec_length1[tm][elevation_parameter_axis] = temp4;
    }
    else
    {
        for (i = 0; i < temp3; i++)
        {
            insert_knot_in_KI[tm][elevation_parameter_axis][i] = vec_temp1[i];
        }
        vec_length1[tm][elevation_parameter_axis] = temp3;
    }

    for (i = 0; i < temp3; i++)
    {
        removal_knot[tm][elevation_parameter_axis][i] = vec_temp1[i];
    }
    vec_length2[tm][elevation_parameter_axis] = temp3;
}


void OE_calc_point_array(int tm)
{
    int i, j;
    int temp1 = 0, temp2 = 0;
    for (i = 0; i < Control_point_n[tm][1]; i++)
    {
        for (j = 0; j < Control_point_n[tm][0]; j++)
        {
            x_array[temp1][temp2] = x[tm][temp1 + (Control_point_n[tm][0] * temp2)];
            y_array[temp1][temp2] = y[tm][temp1 + (Control_point_n[tm][0] * temp2)];
            w_array[temp1][temp2] = w[tm][temp1 + (Control_point_n[tm][0] * temp2)];
            temp1++;
        }
        temp1 = 0;
        temp2++;
    }
}


void OE_define_temp_point_array(int tm, int line_number, int elevation_parameter_axis)
{
    int i;
    if (elevation_parameter_axis == 0)
    {
        for (i = 0; i < Control_point_n[tm][elevation_parameter_axis]; i++)
        {
            temp_x_array[i] = x_array[i][line_number];
            temp_y_array[i] = y_array[i][line_number];
            temp_w_array[i] = w_array[i][line_number];
        }
    }
    else if (elevation_parameter_axis == 1)
    {
        for (i = 0; i < Control_point_n[tm][elevation_parameter_axis]; i++)
        {
            temp_x_array[i] = x_array[line_number][i];
            temp_y_array[i] = y_array[line_number][i];
            temp_w_array[i] = w_array[line_number][i];
        }
    }
}


void Bezier_Order_Elevation_2D(int tm, int elevation_parameter_axis, int Bezier_line_number)
{
    int i, j;
    int n = Order[tm][elevation_parameter_axis];

    for (i = 0; i < OE_n[tm][elevation_parameter_axis]; i++)
    {
        double alpha[MAX_ORDER];
        double a_x, a_y, a_w, b_x, b_y, b_w;

        for (j = 0; j < n + 1; j++)
        {
            if (j != n)
            {
                alpha[j] = (j + 1.0) / (n + 1.0);
                a_w = (1.0 - alpha[j]) * Bezier_w[Bezier_line_number][j + 1];
                a_x = (1.0 - alpha[j]) * (Bezier_x[Bezier_line_number][j + 1] * Bezier_w[Bezier_line_number][j + 1]);
                a_y = (1.0 - alpha[j]) * (Bezier_y[Bezier_line_number][j + 1] * Bezier_w[Bezier_line_number][j + 1]);
                b_w = alpha[j] * Bezier_w[Bezier_line_number][j];
                b_x = alpha[j] * (Bezier_x[Bezier_line_number][j] * Bezier_w[Bezier_line_number][j]);
                b_y = alpha[j] * (Bezier_y[Bezier_line_number][j] * Bezier_w[Bezier_line_number][j]);
                temp_w2[j] = a_w + b_w;
                temp_x2[j] = (a_x + b_x) / temp_w2[j];
                temp_y2[j] = (a_y + b_y) / temp_w2[j];
            }
            else if (j == n)
            {
                temp_w2[j] = Bezier_w[Bezier_line_number][j];
                temp_x2[j] = Bezier_x[Bezier_line_number][j];
                temp_y2[j] = Bezier_y[Bezier_line_number][j];
            }
        }
        if (i != OE_n[tm][elevation_parameter_axis] - 1 && OE_n[tm][elevation_parameter_axis] != 1)
        {
            double temp_Bezeier_x[MAX_ORDER], temp_Bezeier_y[MAX_ORDER], temp_Bezeier_w[MAX_ORDER];

            temp_Bezeier_x[0] = Bezier_x[Bezier_line_number][0];
            temp_Bezeier_y[0] = Bezier_y[Bezier_line_number][0];
            temp_Bezeier_w[0] = Bezier_w[Bezier_line_number][0];
            for (j = 0; j < n + 1; j++)
            {
                temp_Bezeier_x[j + 1] = temp_x2[j];
                temp_Bezeier_y[j + 1] = temp_y2[j];
                temp_Bezeier_w[j + 1] = temp_w2[j];
            }
            
            for (j = 0; j < MAX_ORDER; j++)
            {
                Bezier_x[Bezier_line_number][j] = temp_Bezeier_x[j];
                Bezier_y[Bezier_line_number][j] = temp_Bezeier_y[j];
                Bezier_w[Bezier_line_number][j] = temp_Bezeier_w[j];
            }

            n++;
        }
    }

    for (i = 0; i < n + 1;i++)
    {
        temp_x1[counter] = temp_x2[i];
        temp_y1[counter] = temp_y2[i];
        temp_w1[counter] = temp_w2[i];
        counter++;
    }
}


void Bezier_update_point_array(int tm, int line_number, int elevation_parameter_axis)
{
    int i;
    if (elevation_parameter_axis == 0)
    {
        for (i = 0; i < Control_point_n[tm][elevation_parameter_axis] + (OE_n[tm][elevation_parameter_axis] * number_of_Bezier_line); i++)
        {
            x_array[i][line_number] = temp_x1[i];
            y_array[i][line_number] = temp_y1[i];
            w_array[i][line_number] = temp_w1[i];
        }
    }
    else if (elevation_parameter_axis == 1)
    {
        for (i = 0; i < Control_point_n[tm][elevation_parameter_axis] + (OE_n[tm][elevation_parameter_axis] * number_of_Bezier_line); i++)
        {
            x_array[line_number][i] = temp_x1[i];
            y_array[line_number][i] = temp_y1[i];
            w_array[line_number][i] = temp_w1[i];
        }
    }
}


void Bezier_update(int tm, int elevation_parameter_axis)
{
    int i, j;
    Control_point_n[tm][elevation_parameter_axis] = Control_point_n[tm][elevation_parameter_axis] + (OE_n[tm][elevation_parameter_axis] * number_of_Bezier_line);
    Total_Control_Point[tm] = Control_point_n[tm][0] * Control_point_n[tm][1];
    Order[tm][elevation_parameter_axis] = Order[tm][elevation_parameter_axis] + OE_n[tm][elevation_parameter_axis];

    int temp1 = 0;
    for (i = 0; i < MAX_N_KNOT; i++)
    {
        temp_knot1[i] = 0.0;
    }
    for (i = 0; i < knot_n[tm][elevation_parameter_axis]; i++)
    {
        if (i == 0)
        {
            temp_knot1[temp1] = knot[tm][elevation_parameter_axis][i];
            temp1++;
        }
        if (temp_knot1[temp1 - 1] < knot[tm][elevation_parameter_axis][i])
        {
            temp_knot1[temp1] = knot[tm][elevation_parameter_axis][i];
            temp1++;
        }
        if (temp_knot1[temp1] == 1.0)
        {
            break;
        }
    }

    for (i = 0; i < OE_n[tm][elevation_parameter_axis]; i++)
    {
        for (j = 0; j < MAX_N_KNOT; j++)
        {
            temp_knot2[i] = 0.0;
        }

        int temp2 = 0, temp3 = 0;
        for (j = 0; j < knot_n[tm][elevation_parameter_axis] + temp1; j++)
        {
            if (temp2 < temp1 && temp3 < knot_n[tm][elevation_parameter_axis])
            {
                if (temp_knot1[temp2] <= knot[tm][elevation_parameter_axis][temp3])
                {
                    temp_knot2[temp2 + temp3] = temp_knot1[temp2];
                    temp2++;
                }
                else if (temp_knot1[temp2] > knot[tm][elevation_parameter_axis][temp3])
                {
                    temp_knot2[temp2 + temp3] = knot[tm][elevation_parameter_axis][temp3];
                    temp3++;
                }
            }
            else if (temp2 == temp1)
            {
                temp_knot2[temp2 + temp3] = knot[tm][elevation_parameter_axis][temp3];
                temp3++;
            }
            else if (temp3 == knot_n[tm][elevation_parameter_axis])
            {
                temp_knot2[temp2 + temp3] = temp_knot1[temp2];
                temp2++;
            }

            if (temp2 + temp3 == knot_n[tm][elevation_parameter_axis] + temp1)
            {
                break;
            }
        }

        knot_n[tm][elevation_parameter_axis] = knot_n[tm][elevation_parameter_axis] + temp1;

        for (j = 0; j < knot_n[tm][elevation_parameter_axis]; j++)
        {
            knot[tm][elevation_parameter_axis][j] = temp_knot2[j];
        }
    }

    int temp4 = 0, temp5 = 0;
    for (i = 0; i < Control_point_n[tm][1]; i++)
    {
        for (j = 0; j < Control_point_n[tm][0]; j++)
        {
            x[tm][temp4 + (Control_point_n[tm][0] * temp5)] = x_array[temp4][temp5];
            y[tm][temp4 + (Control_point_n[tm][0] * temp5)] = y_array[temp4][temp5];
            w[tm][temp4 + (Control_point_n[tm][0] * temp5)] = w_array[temp4][temp5];
            temp4++;
        }
        temp4 = 0;
        temp5++;
    }
}


void Calc_Bezier_2D(int tm, int elevation_parameter_axis, int other_axis)
{
    int i, j, k;
    int n, l, l_other;
    l = Control_point_n[tm][elevation_parameter_axis];
    l_other = Control_point_n[tm][other_axis];
    n = Order[tm][elevation_parameter_axis];

    number_of_Bezier_line = l - n;

    for (i = 0; i < number_of_Bezier_line; i++)
    {
        temp_Bezier_array[i] = i * n;
    }

    for (i = 0; i < l_other; i++)
    {
        OE_define_temp_point_array(tm, i, elevation_parameter_axis);

        temp_x1[0] = temp_x_array[0];
        temp_y1[0] = temp_y_array[0];
        temp_w1[0] = temp_w_array[0];
        counter = 1;

        for (j = 0; j < number_of_Bezier_line; j++)
        {
            for (k = 0; k < n + 1; k++)
            {
                Bezier_x[j][k] = temp_x_array[temp_Bezier_array[j] + k];
                Bezier_y[j][k] = temp_y_array[temp_Bezier_array[j] + k];
                Bezier_w[j][k] = temp_w_array[temp_Bezier_array[j] + k];
            }
        }
        for (j = 0; j < number_of_Bezier_line; j++)
        {
            Bezier_Order_Elevation_2D(tm, elevation_parameter_axis, j);
        }
        Bezier_update_point_array(tm, i, elevation_parameter_axis);
    }
}


void KR_calc_point_array(int tm)
{
    int i, j;
    int temp1 = 0, temp2 = 0;
    for (i = 0; i < Control_point_n[tm][1]; i++)
    {
        for (j = 0; j < Control_point_n[tm][0]; j++)
        {
            x_array[temp1][temp2] = x[tm][temp1 + (Control_point_n[tm][0] * temp2)];
            y_array[temp1][temp2] = y[tm][temp1 + (Control_point_n[tm][0] * temp2)];
            w_array[temp1][temp2] = w[tm][temp1 + (Control_point_n[tm][0] * temp2)];
            temp1++;
        }
        temp1 = 0;
        temp2++;
    }
}


void KR_calc_knot_1D(int tm, int removal_parameter_axis)
{
    int i, j;
    for (i = 0; i < MAX_N_KNOT; i++)
    {
        temp_knot1[i] = 0.0;
        temp_knot2[i] = 0.0;
    }

    for (i = 0; i < knot_n[tm][removal_parameter_axis]; i++)
    {
        temp_knot1[i] = knot[tm][removal_parameter_axis][i];
    }

    int temp1 = 0, temp2 = 0, temp3 = 0;
    for (i = 0; i < knot_n[tm][removal_parameter_axis]; i++)
    {
        if (temp2 < vec_length2[tm][removal_parameter_axis])
        {
            if (temp_knot1[temp1] != removal_knot[tm][removal_parameter_axis][temp2])
            {
                temp_knot2[temp3] = temp_knot1[temp1];
                temp1++;
                temp3++;
            }
            else if (temp_knot1[temp1] == removal_knot[tm][removal_parameter_axis][temp2])
            {
                temp1++;
                temp2++;
            }
        }
        else if (temp2 == vec_length2[tm][removal_parameter_axis])
        {
            for (j = 0; j < knot_n[tm][removal_parameter_axis] - temp1; j++)
            {
                temp_knot2[temp3] = temp_knot1[temp1];
                temp1++;
                temp3++;
            }
        }
    }
}


void KR_define_temp_point_array(int tm, int line_number, int removal_parameter_axis)
{
    int i;
    if (removal_parameter_axis == 0)
    {
        for (i = 0; i < Control_point_n[tm][removal_parameter_axis]; i++)
        {
            temp_x_array[i] = x_array[i][line_number];
            temp_y_array[i] = y_array[i][line_number];
            temp_w_array[i] = w_array[i][line_number];
        }
    }
    else if (removal_parameter_axis == 1)
    {
        for (i = 0; i < Control_point_n[tm][removal_parameter_axis]; i++)
        {
            temp_x_array[i] = x_array[line_number][i];
            temp_y_array[i] = y_array[line_number][i];
            temp_w_array[i] = w_array[line_number][i];
        }
    }
}


void KR_calc_Tinv_1D(int tm, int removal_parameter_axis)
{
    int i, j, k;
    int l, n;
    
    l = Control_point_n[tm][removal_parameter_axis] - vec_length2[tm][removal_parameter_axis];
    n = Order[tm][removal_parameter_axis];

    for (i = 0; i < l + vec_length2[tm][removal_parameter_axis]; i++)
    {
        temp_x1[i] = temp_x_array[i] * temp_w_array[i];
        temp_y1[i] = temp_y_array[i] * temp_w_array[i];
        temp_w1[i] = temp_w_array[i];
    }

    double T[n + 1][l + vec_length2[tm][removal_parameter_axis]][l];

    for(int i = 0; i < n + 1; i++)
    {
        for(int j = 0; j < l + vec_length2[tm][removal_parameter_axis]; j++)
        {
            for(int k = 0; k < l; k++)
            {
                T[i][j][k] = 0.0;
            }
        }
    }

    double a = 0.0, b = 0.0;
    for (i = 0; i < l + vec_length2[tm][removal_parameter_axis]; i++)
    {
        for (k = 0; k < n + 1; k++)
        {
            for (j = 0; j < l; j++)
            {
                if (k == 0)
                {
                    if (temp_knot2[j] <= temp_knot1[i] && temp_knot1[i] < temp_knot2[j + 1])
                    {
                        T[k][i][j] = 1.0;
                    }
                    else
                    {
                        T[k][i][j] = 0.0;
                    }
                }
                else
                {
                    if (temp_knot2[j + k] - temp_knot2[j] == 0)
                    {
                        a = 0.0;
                    }
                    else if (temp_knot2[j + k] - temp_knot2[j] != 0)
                    {
                        a = ((temp_knot1[i + k] - temp_knot2[j]) / (temp_knot2[j + k] - temp_knot2[j])) * T[k - 1][i][j];
                    }
                    if (temp_knot2[j + k + 1] - temp_knot2[j + 1] == 0)
                    {
                        b = 0.0;
                    }
                    else if (temp_knot2[j + k + 1] - temp_knot2[j + 1] != 0)
                    {
                        b = ((temp_knot2[j + k + 1] - temp_knot1[i + k]) / (temp_knot2[j + k + 1] - temp_knot2[j + 1])) * T[k - 1][i][j + 1];
                    }
                    T[k][i][j] = a + b;
                }
            }
        }
    }

    //疑似逆行列の計算
    // if m <= n, A_pinv = A^T (A A^T)^(-1)
    // (A A^T)^(-1)の解がない, A_pinv = (A (A^T A)^(-1))^T

    double A[l + vec_length2[tm][removal_parameter_axis]][l];
    double A_T[l][l + vec_length2[tm][removal_parameter_axis]];
    double B[l][l];
    double B_inv[l][l];
    double A_pinv[l + vec_length2[tm][removal_parameter_axis]][l];
    double A_pinv_T[l][l + vec_length2[tm][removal_parameter_axis]];

    for (i = 0; i < l + vec_length2[tm][removal_parameter_axis]; i++)
    {
        for (j = 0; j < l; j++)
        {
            A[i][j] = T[n][i][j];
        }
    }

    for (i = 0; i < l + vec_length2[tm][removal_parameter_axis]; i++)
    {
        for (j = 0; j < l; j++)
        {
            A_T[j][i] = A[i][j];
        }
    }

    for (i = 0; i < l; i++)
    {
        for (j = 0; j < l; j++)
        {
            B[i][j] = 0.0;
            for (k = 0; k < l + vec_length2[tm][removal_parameter_axis]; k++)
            {
                B[i][j] += A_T[i][k] * A[k][j];
            }
        }
    }

    double temp;

    for (i = 0; i < l; i++)
    {
        for (j = 0; j < l; j++)
        {
            B_inv[i][j] = 0.0;
        }
    }

    for (i = 0; i < l; i++)
    {
        B_inv[i][i] = 1.0;
    }

    for (k = 0; k < l; k++)
    {
        temp = B[k][k];
        for (i = 0; i < l; i++)
        {
            B[k][i] /= temp;
            B_inv[k][i] /= temp;
        }
        
        for (i = 0; i < l; i++)
        {
            if (i != k)
            {
                temp = B[i][k];
                for (j = 0; j < l; j++)
                {
                    B[i][j] -= B[k][j] * temp;
                    B_inv[i][j] -= B_inv[k][j] * temp;
                }
            }
        }
    }

    for (i = 0; i < l + vec_length2[tm][removal_parameter_axis]; i++)
    {
        for (j = 0; j < l; j++)
        {
            A_pinv[i][j] = 0.0;
            for (k = 0; k < l; k++)
            {
                A_pinv[i][j] += A[i][k] * B_inv[k][j];
            }
        }
    }

    for (i = 0; i < l + vec_length2[tm][removal_parameter_axis]; i++)
    {
        for (j = 0; j < l; j++)
        {
            A_pinv_T[j][i] = A_pinv[i][j];
        }
    }

    for (i = 0; i < l; i++)
    {
        for (j = 0; j < l + vec_length2[tm][removal_parameter_axis]; j++)
        {
            if (j == 0)
            {
                temp_x2[i] = 0.0;
                temp_y2[i] = 0.0;
                temp_w2[i] = 0.0;
            }
            temp_x2[i] += A_pinv_T[i][j] * temp_x1[j];
            temp_y2[i] += A_pinv_T[i][j] * temp_y1[j];
            temp_w2[i] += A_pinv_T[i][j] * temp_w1[j];
        }
    }

    for (i = 0; i < l; i++)
    {
        temp_x_array[i] = temp_x2[i] / temp_w2[i];
        temp_y_array[i] = temp_y2[i] / temp_w2[i];
        temp_w_array[i] = temp_w2[i];
    }
}


void KR_update_point_array(int tm, int line_number, int removal_parameter_axis)
{
    int i;
    if (removal_parameter_axis == 0)
    {
        for (i = 0; i < Control_point_n[tm][removal_parameter_axis] - vec_length2[tm][removal_parameter_axis]; i++)
        {
            x_array[i][line_number] = temp_x_array[i];
            y_array[i][line_number] = temp_y_array[i];
            w_array[i][line_number] = temp_w_array[i];
        }
    }
    else if (removal_parameter_axis == 1)
    {
        for (i = 0; i < Control_point_n[tm][removal_parameter_axis] - vec_length2[tm][removal_parameter_axis]; i++)
        {
            x_array[line_number][i] = temp_x_array[i];
            y_array[line_number][i] = temp_y_array[i];
            w_array[line_number][i] = temp_w_array[i];
        }
    }
}


void KR_update(int tm, int removal_parameter_axis)
{
    int i, j;
    for (i = 0; i < knot_n[tm][removal_parameter_axis] - vec_length2[tm][removal_parameter_axis]; i++)
    {
        knot[tm][removal_parameter_axis][i] = temp_knot2[i];
    }
    Control_point_n[tm][removal_parameter_axis] = Control_point_n[tm][removal_parameter_axis] - vec_length2[tm][removal_parameter_axis];
    Total_Control_Point[tm] = Control_point_n[tm][0] * Control_point_n[tm][1];
    knot_n[tm][removal_parameter_axis] = Control_point_n[tm][removal_parameter_axis] + Order[tm][removal_parameter_axis] + 1;

    int temp1 = 0, temp2 = 0;
    for (i = 0; i < Control_point_n[tm][1]; i++)
    {
        for (j = 0; j < Control_point_n[tm][0]; j++)
        {
            x[tm][temp1 + (Control_point_n[tm][0] * temp2)] = x_array[temp1][temp2];
            y[tm][temp1 + (Control_point_n[tm][0] * temp2)] = y_array[temp1][temp2];
            w[tm][temp1 + (Control_point_n[tm][0] * temp2)] = w_array[temp1][temp2];
            temp1++;
        }
        temp1 = 0;
        temp2++;
    }
}


void KR_reset_array()
{
    int i, j;
    for (i = 0; i < MAX_N_Controlpoint_each_parameter; i++)
    {
        for (j = 0; j < MAX_N_Controlpoint_each_parameter; j++)
        {
            x_array[i][j] = 0.0;
            y_array[i][j] = 0.0;
            w_array[i][j] = 0.0;
        }
    }
    
    for (i = 0; i < MAX_N_Controlpoint_each_parameter; i++)
    {
        temp_x_array[i] = 0.0;
        temp_y_array[i] = 0.0;
        temp_w_array[i] = 0.0;
    }

    for (i = 0; i < MAX_N_Controlpoint_in_Patch; i++)
    {
        temp_x1[i] = 0.0;
        temp_y1[i] = 0.0;
        temp_w1[i] = 0.0;
        temp_x2[i] = 0.0;
        temp_y2[i] = 0.0;
        temp_w2[i] = 0.0;
    }
}


void KR_non_uniform_2D(int tm, int removal_parameter_axis)
{
    int i;
    int other_axis = 0;
    if (removal_parameter_axis == 0)
    {
        other_axis = 1;
    }
    else if (removal_parameter_axis == 1)
    {
        other_axis = 0;
    }

    KR_calc_point_array(tm);

    for (i = 0; i < Control_point_n[tm][other_axis]; i++)
    {
        if (i == 0)
        {
            KR_calc_knot_1D(tm, removal_parameter_axis);
        }

        KR_define_temp_point_array(tm, i, removal_parameter_axis);
        KR_calc_Tinv_1D(tm, removal_parameter_axis);
        KR_update_point_array(tm, i, removal_parameter_axis);
    }

    KR_update(tm, removal_parameter_axis);
    KR_reset_array();
}


void OE_2D(int tm, int elevation_parameter_axis)
{
    int other_axis = 0;
    if (elevation_parameter_axis == 0)
    {
        other_axis = 1;
    }
    else if (elevation_parameter_axis == 1)
    {
        other_axis = 0;
    }

    if (OE_n[tm][elevation_parameter_axis] > 0)
    {
        Calc_insert_knot_in_OE(tm, elevation_parameter_axis);
        KI_non_uniform_2D(tm, elevation_parameter_axis, 0);

        OE_calc_point_array(tm);
        Calc_Bezier_2D(tm, elevation_parameter_axis, other_axis);
        Bezier_update(tm, elevation_parameter_axis);

        KR_non_uniform_2D(tm, elevation_parameter_axis);
    }
}


void Debug_printf(int tm, char *section)
{
    int i, j;

    printf("----------------------------------------------------------\n");
    printf("input file number:\t%d\n", tm);
    printf("status:\t\t\t\t%s\n\n", section);

    printf("Dimension\n");
    printf("%d\n\n", Dimension[tm]);

    printf("Total Control Point\n");
    printf("%d\n\n", Total_Control_Point[tm]);

    printf("Order\n");
    for (j = 0; j < Dimension[tm] ; j++)
    {
        if (j == 0)
        {
            printf("%d", Order[tm][j]);
        }
        else
        {
            printf("\t%d", Order[tm][j]);
        }
    }
    printf("\n\n");

    printf("number of knot vector\n");
    for (j = 0; j < Dimension[tm]; j++)
    {
        if (j == 0)
        {
            printf("%d", knot_n[tm][j]);
        }
        else
        {
            printf("\t%d", knot_n[tm][j]);
        }
    }
    printf("\n\n");

    printf("number of Control Point\n");
    for (j = 0; j < Dimension[tm]; j++)
    {
        if (j == 0)
        {
            printf("%d", Control_point_n[tm][j]);
        }
        else
        {
            printf("\t%d", Control_point_n[tm][j]);
        }
    }
    printf("\n\n");

    printf("knot vector\n");
    for (j = 0; j < Dimension[tm]; j++)
    {
        for (i = 0; i < knot_n[tm][j]; i++)
        {
            if (i == 0)
            {
                printf("%.16e", knot[tm][j][i]);
            }
            else
            {
                printf("\t%.16e", knot[tm][j][i]);
            }
        }
        printf("\n");
    }
    printf("\n");

    printf("coordinate of Control Point\n");
    for (i = 0; i < Total_Control_Point[tm]; i++)
    {
        printf("%d\t", i);
        printf("%.16e\t", x[tm][i]);
        printf("%.16e\t", y[tm][i]);
        printf("%.16e\n", w[tm][i]);
    }
    printf("\n");

    printf("Order of Elevation\n");
    for (j = 0; j < Dimension[tm]; j++)
    {
        if (j == 0)
        {
            printf("%d", OE_n[tm][j]);
        }
        else
        {
            printf("\t%d", OE_n[tm][j]);
        }
    }
    printf("\n\n");

    printf("number of uniform Knot Insertion\n");
    for (j = 0; j < Dimension[tm]; j++)
    {
        if (j == 0)
        {
            printf("%d", KI_uniform_n[tm][j]);
        }
        else
        {
            printf("\t%d", KI_uniform_n[tm][j]);
        }
    }
    printf("\n\n");

    printf("number of non-uniform Knot Insertion\n");
    for (j = 0; j < Dimension[tm]; j++)
    {
        if (j == 0)
        {
            printf("%d", KI_non_uniform_n[tm][j]);
        }
        else
        {
            printf("\t%d", KI_non_uniform_n[tm][j]);
        }
    }
    printf("\n\n");

    printf("insert knot\n");
    for (j = 0; j < Dimension[tm]; j++)
    {
        if (KI_non_uniform_n[tm][j] != 0)
        {
            for (i = 0; i < KI_non_uniform_n[tm][j]; i++)
            {
                if (i == 0)
                {
                    printf("%.16e", insert_knot[tm][j][i]);
                }
                else
                {
                    printf("\t%.16e", insert_knot[tm][j][i]);
                }
            }
            printf("\n");
        }
    }
}


void Fix_numerical_error(int tm)
{
    int i;
    for (i = 0; i < Total_Control_Point[tm]; i++)
    {
        if (fabs(w[tm][i] - 1.0) < 1.0e-15)
        w[tm][i] = 1.0;
    }
}


void OutputData(int tm, char *filename)
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

    // fprintf(fp, "input file name: ");
    // fprintf(fp, "%s\n\n", filename);

    // fprintf(fp, "Total Control Point\n");
    fprintf(fp, "%d\n\n", Total_Control_Point[tm]);

    // fprintf(fp, "Order\n");
    for (j = 0; j < Dimension[tm]; j++)
    {
        if (j == Dimension[tm] - 1)
        {
            fprintf(fp, "%d", Order[tm][j]);
        }
        else
        {
            fprintf(fp, "%-5d", Order[tm][j]);
        }
    }
    fprintf(fp, "\n\n");

    // fprintf(fp, "number of Knot Vector\n");
    for (j = 0; j < Dimension[tm]; j++)
    {
        if (j == Dimension[tm] - 1)
        {
            fprintf(fp, "%d", knot_n[tm][j]);
        }
        else
        {
            fprintf(fp, "%-5d", knot_n[tm][j]);
        }
    }
    fprintf(fp, "\n\n");

    // fprintf(fp, "number of Control Point\n");
    for (j = 0; j < Dimension[tm]; j++)
    {
        if (j == Dimension[tm] - 1)
        {
            fprintf(fp, "%d", Control_point_n[tm][j]);
        }
        else
        {
            fprintf(fp, "%-5d", Control_point_n[tm][j]);
        }
    }
    fprintf(fp, "\n\n");

    // fprintf(fp, "knot vector\n");
    for (j = 0; j < Dimension[tm]; j++)
    {
        for (i = 0; i < knot_n[tm][j]; i++)
        {
            if (i == knot_n[tm][j] - 1)
            {
                fprintf(fp, "%.16e", knot[tm][j][i]);
            }
            else
            {
                fprintf(fp, "%.16e  ", knot[tm][j][i]);
            }
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // fprintf(fp, "coordinate of Control Point\n");
    for (i = 0; i < Total_Control_Point[tm]; i++)
    {
        fprintf(fp, "%-5d", i);
        fprintf(fp, "  %.16e", x[tm][i]);
        fprintf(fp, "  %.16e", y[tm][i]);
        fprintf(fp, "  %.16e\n", w[tm][i]);
    }

    fclose(fp);
}