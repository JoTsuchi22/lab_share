/****************************************************************
 * s-IGA :)
 *
 * 仮定条件
 * 	ローカルメッシュ同士は原則被りなしと仮定
 * 	グローバルパッチはシングルパッチ
 * 
 * 要素の重なりの判定はガウス点で行っている
 *
 * ガウス積分の積分点数：4(一部10)
 *
 * inputファイルを複数読み込む, 2つ目以降がローカルメッシュ
 *
 ****************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>

// header
#include "s_IGA_header.h"
#include "s_IGA_main.h"

int main(int argc, char **argv)
{
	clock_t start, end;

	int i, j, k;
	
	int K_Whole_Size = 0;
	double element_loc[DIMENSION];

	Total_mesh = argc - 1;

	// for s-IGA
	int tm;

	// 引数の個数確認
	if (argc <= 1)
	{
		printf("Argument is missing\n");
	}
	if (argc == 2) // 通常IGA: input file 1つ
	{
		printf("IGA carried out.(No local mesh)\n");
	}
	if (argc >= 3) // s-IGA: input file 複数
	{
		printf("s-IGA carried out.(%d local meshes)\n", argc - 2);
	}

	start = clock();

	// memory allocation
	int *Total_Knot_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));
	int *Total_Patch_on_mesh = (int *)malloc(sizeof(int) * (Total_mesh));				// 各メッシュ上のパッチ数
	int *Total_Patch_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));			// メッシュ[]までのパッチ数（メッシュ[]内のパッチ数は含まない）
	int *Total_Control_Point_on_mesh = (int *)malloc(sizeof(int) * (Total_mesh));		// 各メッシュ上のコントロールポイント数
	int *Total_Control_Point_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));	// メッシュ[]までのコントロールポイント数(メッシュ[]内のコントロールポイント数は含まない)
	int *Total_Element_on_mesh = (int *)malloc(sizeof(int) * (Total_mesh));				// 各メッシュ上の要素数
	int *Total_Element_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));			// メッシュ[]までの要素数(メッシュ[]内の要素数は含まない)
	int *real_Total_Element_on_mesh = (int *)malloc(sizeof(int) * (Total_mesh));		// real_Total_Element_on_mesh[MAX_N_MESH]
	int *real_Total_Element_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));		// real_Total_Element_to_mesh[MAX_N_MESH + 1]
	int *Total_Load_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));
	int *Total_Constraint_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));
	int *Total_DistributeForce_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));

	// ファイル読み込み1回目
    for (tm = 0; tm < Total_mesh; tm++)
    {
		Get_Input_1(tm, Total_Knot_to_mesh,
					Total_Patch_on_mesh, Total_Patch_to_mesh,
					Total_Control_Point_on_mesh, Total_Control_Point_to_mesh,
					Total_Load_to_mesh, Total_Constraint_to_mesh, Total_DistributeForce_to_mesh,
					argv);
	}
	
	// memory allocation
	int *Order = (int *)malloc(sizeof(int) * (Total_Patch_to_mesh[Total_mesh] * DIMENSION));				// Order[MAX_N_PATCH][DIMENSION]
	int *No_knot = (int *)malloc(sizeof(int) * (Total_Patch_to_mesh[Total_mesh] * DIMENSION));				// No_knot[MAX_N_PATCH][DIMENSION]
	int *Total_Control_Point_to_patch = (int *)calloc((Total_Patch_to_mesh[Total_mesh] + 1), sizeof(int));	// Total_Control_Point_to_patch[MAX_N_PATCH]
	int *Total_Knot_to_patch_dim = (int *)calloc((Total_Patch_to_mesh[Total_mesh] * DIMENSION + 1), sizeof(int));	// Total_Knot_to_patch_dim[MAX_N_PATCH][DIMENSION]
	double *Position_Knots = (double *)malloc(sizeof(double) * Total_Knot_to_mesh[Total_mesh]);				// Position_Knots[MAX_N_PATCH][DIMENSION][MAX_N_KNOT];
	int *No_Control_point = (int *)malloc(sizeof(int) * (Total_Patch_to_mesh[Total_mesh] * DIMENSION));		// Order[MAX_N_PATCH][DIMENSION]
	int *No_Control_point_in_patch = (int *)malloc(sizeof(int) * (Total_Patch_to_mesh[Total_mesh]));		// No_Control_point_in_patch[MAX_N_PATCH]
	int *Patch_Control_point = (int *)malloc(sizeof(int) * (Total_Control_Point_to_mesh[Total_mesh]));		// Patch_Control_point[MAX_N_PATCH][MAX_N_Controlpoint_in_Patch]
	int *No_Control_point_ON_ELEMENT = (int *)malloc(sizeof(int) * (Total_Patch_to_mesh[Total_mesh]));		// No_Control_point_ON_ELEMENT[MAX_N_PATCH]
	double *Node_Coordinate = (double *)malloc(sizeof(double) * (Total_Control_Point_to_mesh[Total_mesh] * (DIMENSION + 1)));	// Node_Coordinate[MAX_N_NODE][DIMENSION + 1];
	double *Control_Coord_x = (double *)malloc(sizeof(double) * (Total_Control_Point_to_mesh[Total_mesh])); // Control_Coord[DIMENSION][MAX_N_NODE];
	double *Control_Coord_y = (double *)malloc(sizeof(double) * (Total_Control_Point_to_mesh[Total_mesh])); // Control_Coord[DIMENSION][MAX_N_NODE];
	double *Control_Weight = (double *)malloc(sizeof(double) * (Total_Control_Point_to_mesh[Total_mesh]));	// Control_Weight[MAX_N_NODE];
	int *Constraint_Node_Dir = (int *)malloc(sizeof(int) * (Total_Constraint_to_mesh[Total_mesh] * 2));		// Constraint_Node_Dir[MAX_N_CONSTRAINT][2];
	double *Value_of_Constraint = (double *)calloc(Total_Constraint_to_mesh[Total_mesh], sizeof(double));	// Value_of_Constraint[MAX_N_CONSTRAINT];
	int *Load_Node_Dir = (int *)malloc(sizeof(int) * (Total_Load_to_mesh[Total_mesh] * 2));					// Load_Node_Dir[MAX_N_LOAD][2];
	double *Value_of_Load = (double *)calloc(Total_Load_to_mesh[Total_mesh], sizeof(double));				// Value_of_Load[MAX_N_LOAD];
	int *iPatch_array = (int *)malloc(sizeof(int) * Total_DistributeForce_to_mesh[Total_mesh]);				// iPatch_array[MAX_N_DISTRIBUTE_FORCE]
	int *iCoord_array = (int *)malloc(sizeof(int) * Total_DistributeForce_to_mesh[Total_mesh]);				// iCoord_array[MAX_N_DISTRIBUTE_FORCE]
	int *type_load_array = (int *)malloc(sizeof(int) * Total_DistributeForce_to_mesh[Total_mesh]);			// type_load_array[MAX_N_DISTRIBUTE_FORCE]
	double *val_Coord_array = (double *)calloc(Total_DistributeForce_to_mesh[Total_mesh], sizeof(double));	// val_Coord_array[MAX_N_DISTRIBUTE_FORCE]
	double *Range_Coord_array = (double *)calloc((Total_DistributeForce_to_mesh[Total_mesh] * 2), sizeof(double));	// Range_Coord_array[MAX_N_DISTRIBUTE_FORCE][2]
	double *Coeff_Dist_Load_array = (double *)calloc((Total_DistributeForce_to_mesh[Total_mesh] * 3), sizeof(double));	// Coeff_Dist_Load_array[MAX_N_DISTRIBUTE_FORCE][3] 

	// ファイル読み込み2回目
	for (tm = 0; tm < Total_mesh; tm++)
    {
		Get_Input_2(tm, Total_Knot_to_mesh,
					Total_Patch_to_mesh, Total_Control_Point_to_mesh,
					Total_Element_on_mesh, Total_Element_to_mesh,
					Total_Constraint_to_mesh, Total_Load_to_mesh, Total_DistributeForce_to_mesh,
					Order, No_knot, No_Control_point, No_Control_point_in_patch,
					Patch_Control_point, Position_Knots, No_Control_point_ON_ELEMENT,
					Node_Coordinate, Control_Coord_x, Control_Coord_y, Control_Weight,
					Constraint_Node_Dir, Value_of_Constraint,
					Load_Node_Dir, Value_of_Load,
					type_load_array, iPatch_array, iCoord_array,
					val_Coord_array, Range_Coord_array, Coeff_Dist_Load_array,
					Total_Control_Point_to_patch, Total_Knot_to_patch_dim,
					argv);
	}

	// memory allocation
	int *INC = (int *)malloc(sizeof(int) * (Total_Control_Point_to_mesh[Total_mesh] * DIMENSION));		// INC[MAX_N_NODE][DIMENSION]
	int *Controlpoint_of_Element = (int *)malloc(sizeof(int) * (Total_Element_to_mesh[Total_mesh] * MAX_NO_CCpoint_ON_ELEMENT)); // Controlpoint_of_Element[MAX_N_ELEMENT][MAX_NO_CCpoint_ON_ELEMENT]
	int *Element_patch = (int *)malloc(sizeof(int) * Total_Element_to_mesh[Total_mesh]);				// Element_patch[MAX_N_ELEMENT]
	int *Element_mesh = (int *)malloc(sizeof(int) * Total_Element_to_mesh[Total_mesh]);					// Element_mesh[MAX_N_ELEMENT] 要素がどのメッシュ内にあるかを示す配列
	int *line_No_real_element = (int *)malloc(sizeof(int) * (Total_Patch_to_mesh[Total_mesh] * DIMENSION));		// line_No_real_element[MAX_N_PATCH][DIMENSION] ゼロエレメントではない要素列の数
	int *line_No_Total_element = (int *)malloc(sizeof(int) * (Total_Patch_to_mesh[Total_mesh] * DIMENSION));	// line_No_Total_element[MAX_N_PATCH][DIMENSION] ゼロエレメントを含むすべての要素列の数
	double *difference = (double *)calloc(Total_Knot_to_mesh[Total_mesh], sizeof(double));				// difference[MAX_N_PATCH][DIMENSION][MAX_N_KNOT] 隣り合うノットベクトルの差
	int *Total_element_all_ID = (int *)calloc(Total_Element_to_mesh[Total_mesh], sizeof(int));			// Total_element_all_ID[MAX_N_ELEMENT] ゼロエレメントではない要素 = 1, ゼロエレメント = 0
	int *ENC = (int *)malloc(sizeof(int) * (Total_Element_to_mesh[Total_mesh] * DIMENSION));			// ENC[MAX_N_ELEMENT][DIMENSION] ENC[パッチ][全ての要素][0, 1] = x, y方向の何番目の要素か
	int *real_element_line = (int *)malloc(sizeof(int) * (Total_Patch_to_mesh[Total_mesh] * Total_Element_to_mesh[Total_mesh] * DIMENSION)); // real_element_line[MAX_N_PATCH][MAX_N_ELEMENT][DIMENSION] ゼロエレメントではない要素列
	int *real_element = (int *)malloc(sizeof(int) * Total_Element_to_mesh[Total_mesh]);					// real_element[MAX_N_ELEMENT] ゼロエレメントではない要素の番号
	int *real_El_No_on_mesh = (int *)malloc(sizeof(int) * (Total_Patch_to_mesh[Total_mesh] * Total_Element_to_mesh[Total_mesh]));	// real_El_No_on_mesh[MAX_N_MESH][MAX_N_ELEMENT]
	double *Equivalent_Nodal_Force = (double *)calloc(Total_Control_Point_to_mesh[Total_mesh] * DIMENSION, sizeof(double));	// Equivalent_Nodal_Force[MAX_N_NODE][DIMENSION] Equivalent nodal forces arising from the distributed load

	// INC 等の作成
	Make_INC(tm, Total_Patch_on_mesh, Total_Patch_to_mesh,
			 Total_Element_on_mesh, Total_Element_to_mesh,
			 Total_Control_Point_on_mesh, Total_Control_Point_to_mesh,
			 Total_DistributeForce_to_mesh,
			 INC, Patch_Control_point, Total_Control_Point_to_patch,
			 No_Control_point, Controlpoint_of_Element,
			 Element_patch, Element_mesh,
			 difference, Total_Knot_to_patch_dim,
			 Total_element_all_ID, ENC,
			 real_element_line, line_No_real_element,
			 line_No_Total_element, real_element,
			 real_El_No_on_mesh, real_Total_Element_on_mesh,
			 real_Total_Element_to_mesh, Equivalent_Nodal_Force,
			 type_load_array, iPatch_array, iCoord_array,
			 val_Coord_array, Range_Coord_array, Coeff_Dist_Load_array,
			 Order, No_knot);

	// memory free
	free(Patch_Control_point), free(real_element_line), free(Total_element_all_ID), free(difference);

	// memory allocation
	double *Gauss_Coordinate = (double *)calloc(real_Total_Element_to_mesh[Total_mesh] * POW_Ng * DIMENSION, sizeof(double));
	double *Gauss_Coordinate_ex = (double *)calloc(real_Total_Element_to_mesh[Total_mesh] * POW_Ng_extended * DIMENSION, sizeof(double));
	double *Jac = (double *)malloc(sizeof(double) * real_Total_Element_to_mesh[Total_mesh] * POW_Ng);
	double *Jac_ex = (double *)malloc(sizeof(double) * real_Total_Element_to_mesh[Total_mesh] * POW_Ng_extended);
	double *B_Matrix = (double *)malloc(sizeof(double) * real_Total_Element_to_mesh[Total_mesh] * POW_Ng * D_MATRIX_SIZE);
	double *B_Matrix_ex = (double *)malloc(sizeof(double) * real_Total_Element_to_mesh[Total_mesh] * POW_Ng_extended * D_MATRIX_SIZE);
	
	// check_over_parameter
	for (i = 1; i < Total_mesh; i++)
	{
		int mesh_n_org = 0;
		printf("mesh_n_org: 0\tmesh_n_over: %d\n", i);
		// NNLOVER の算出
		Check_coupled_Glo_Loc_element_for_Gauss(i, mesh_n_org);
		Make_Loc_Glo();
	}



    // 全体剛性マトリックスの制作
	K_Whole_Size = Make_Index_Dof(Total_Control_Point_to_mesh[Total_mesh],
								  Total_Constraint_to_mesh[Total_mesh],
								  Constraint_Node_Dir);

    printf("K_Whole_Size=%d\n",K_Whole_Size);
    Make_K_Whole_Ptr_Col(Total_Element_to_mesh[Total_mesh], Total_Control_Point_to_mesh[Total_mesh], K_Whole_Size);
    Make_K_Whole_Val(E, nu, real_Total_Element_to_mesh[Total_mesh], DM);
    printf("Finish Make_K_Whole\n");

	// memory free
	free(ENC);

	// for s-IGA　複数メッシュのループ内に移動
	// 荷重ベクトルの算出部分
	for (i = 0; i < Total_Load_to_mesh[Total_mesh]; i++)
	{
		printf("Value_of_Load;%.20e\n", Value_of_Load[i]);
	}
	printf("pp\n");
	Make_F_Vec(Total_Load_to_mesh[Total_mesh], Load_Node_Dir, Value_of_Load, K_Whole_Size);
	Make_F_Vec_disp_const(tm, Total_Constraint_to_mesh[Total_mesh], Constraint_Node_Dir, Value_of_Constraint, E, nu, DM);
	Add_Equivalent_Nodal_Forec_to_F_Vec(Total_Control_Point_to_mesh[Total_mesh]);
    printf("Finish Make_K_Whole\n");

	// memory free
	free(Equivalent_Nodal_Force);

	// Kマトリックスのsvg出力
	// K_output_svg(K_Whole_Size);

	// 連立一次方程式
    // 反復回数の設定
    int max_itr = K_Whole_Size;

	// CG法
	// Diag_Scaling_CG_pre(K_Whole_Size, 0);
    // printf("Finish 1st Diag_Scaling_CG_Pre\n");
	// CG_Solver(K_Whole_Size, max_itr, EPS, 0);
	// Diag_Scaling_CG_pre(K_Whole_Size, 1);
	// printf("Finish CG_Solver\n");

	// PCG法
    printf("\nStart PCG solver\n");
	PCG_Solver(K_Whole_Size, max_itr, EPS);
	printf("Finish PCG solver\n\n");

	// 変位と歪と応力
    for(i = 0; i < Total_Constraint_to_mesh[Total_mesh]; i++)
    {
		Constraint_ID[Constraint_Node_Dir[i][0] * DIMENSION + Constraint_Node_Dir[i][1]] = 1;
		Displacement[Constraint_Node_Dir[i][0] * DIMENSION + Constraint_Node_Dir[i][1]] = Value_of_Constraint[i];
    }

	for(i = 0; i < Total_Control_Point_to_mesh[Total_mesh]; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			int index = Index_Dof[i * DIMENSION + j];
			if (index >= 0)
				Displacement[i * DIMENSION + j] = sol_vec[index];
		}
	}
	printf("Finish Make_Displacement\n");
	end = clock();
	printf("Analysis time:%.2f[s]\n",(double)(end - start) / CLOCKS_PER_SEC);
	Make_Strain(real_Total_Element_to_mesh[Total_mesh]);
	printf("Finish Make_Strain\n");
	Make_Stress_2D(E, nu, real_Total_Element_to_mesh[Total_mesh], DM);
	printf("Finish Make_Stress\n");
	Make_ReactionForce(Total_Control_Point_to_mesh[Total_mesh]);
	printf("Finish Make_ReactionForce\n");
	puts("sol_vec");
	Make_Parameter_z(real_Total_Element_to_mesh[Total_mesh], E, nu, DM);
	printf("Finish Make_Parameter_z\n");

	// 全てのローカルパッチについて1つのinput.txt作成
	// (重ね合わせ結果出力のため)
	int l_patch, g_patch;
	int l_cntl_p, g_cntl_p;
	int l_cnst, g_cnst;
	int l_load, g_load;
	int l_dist, g_dist;
	int l_temp;

	g_patch = Total_Patch_to_mesh[1];
	g_cntl_p = Total_Control_Point_to_mesh[1];
	g_cnst = Total_Constraint_to_mesh[1];
	g_load = Total_Load_to_mesh[1];
	g_dist = Total_DistributeForce_to_mesh[1];

	fp = fopen("input_local.txt", "w");
	fprintf(fp, "%.3f\t%.6f\n\n", E, nu);

	l_patch = Total_Patch_to_mesh[Total_mesh] - g_patch;
	fprintf(fp, "%d\n\n", l_patch);

	l_cntl_p = Total_Control_Point_to_mesh[Total_mesh] - g_cntl_p;
	fprintf(fp, "%d\n\n", l_cntl_p);

	for (i = 0; i < l_patch; i++)
	{
		fprintf(fp, "%d\t%d\n",
				Order[i + g_patch][0], Order[i + g_patch][1]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < l_patch; i++)
	{
		fprintf(fp, "%d\t%d\n",
				No_knot[i + g_patch][0], No_knot[i + g_patch][1]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < l_patch; i++)
	{
		fprintf(fp, "%d\t%d\n",
				No_Control_point[i + g_patch][0],
				No_Control_point[i + g_patch][1]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < l_patch; i++)
	{
		for (j = 0; j < No_Control_point_in_patch[i + g_patch]; j++)
		{
			l_temp = Patch_Control_point[i + g_patch][j] - g_cntl_p;
			fprintf(fp, "%d\t", l_temp);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp,"\n");

	l_cnst = Total_Constraint_to_mesh[Total_mesh] - g_cnst;
	l_load = Total_Load_to_mesh[Total_mesh] - g_load;
	l_dist = Total_DistributeForce_to_mesh[Total_mesh] - g_dist;
	fprintf(fp,"%d\t%d\t%d\n\n",
				l_cnst, l_load,	l_dist);

	for (i = 0; i < l_patch; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			for (k = 0; k < No_knot[i + g_patch][j]; k++)
			{
				fprintf(fp, "%.5f\t", Position_Knots[i + g_patch][j][k]);
			}
			fprintf(fp, "\n");
		}
	}
	fprintf(fp,"\n");

	for (i = 0; i < l_cntl_p; i++)
	{
		fprintf(fp, "%d\t", i);
		for (j = 0; j < DIMENSION + 1; j++)
		{
			fprintf(fp, "%.10f\t", Node_Coordinate[i + g_cntl_p][j]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp,"\n");

	for (i = 0; i < l_cnst; i++)
	{
		l_temp = Constraint_Node_Dir[i + g_cnst][0] - g_cntl_p;
		fprintf(fp, "%d\t%d\t%.10f\n",
				l_temp,
				Constraint_Node_Dir[i + g_cnst][1],
				Value_of_Constraint[i + g_cnst]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < l_load; i++)
	{
		l_temp = Load_Node_Dir[i + g_load][0] - g_cntl_p;
		fprintf(fp, "%d\t%d\t%.10f\n",
				l_temp,
				Load_Node_Dir[i + g_load][1],
				Value_of_Load[i + g_load]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < l_dist; i++)
	{
		fprintf(fp, "%d\t%d\t%d\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n",
				type_load_array[i + g_dist],
				iPatch_array[i + g_dist],
				iCoord_array[i + g_dist],
				val_Coord_array[i + g_dist],
				Range_Coord_array[i + g_dist][0],
				Range_Coord_array[i + g_dist][1],
				Coeff_Dist_Load_array[i + g_dist][0],
				Coeff_Dist_Load_array[i + g_dist][1],
				Coeff_Dist_Load_array[i + g_dist][2]);
	}
	fclose(fp);

	fp = fopen("input_for_NURBS.txt", "w");
	fprintf(fp, "%.3f\t%.6f\n\n", E, nu);

	fprintf(fp, "%d\n\n", Total_Patch_to_mesh[Total_mesh]);

	fprintf(fp, "%d\n\n", Total_Control_Point_to_mesh[Total_mesh]);

	for (i = 0; i < Total_Patch_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp, "%d\t%d\n", Order[i][0], Order[i][1]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < Total_Patch_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp, "%d\t%d\n", No_knot[i][0], No_knot[i][1]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < Total_Patch_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp, "%d\t%d\n", No_Control_point[i][0], No_Control_point[i][1]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < Total_Patch_to_mesh[Total_mesh]; i++)
	{
		for (j = 0; j < No_Control_point_in_patch[i]; j++)
		{
			fprintf(fp, "%d\t", Patch_Control_point[i][j]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp,"\n");

	fprintf(fp,"%d\t%d\t%d\n\n",
				Total_Constraint_to_mesh[Total_mesh],
				Total_Load_to_mesh[Total_mesh],
				Total_DistributeForce_to_mesh[Total_mesh]);

	for (i = 0; i < Total_Patch_to_mesh[Total_mesh]; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			for (k = 0; k < No_knot[i][j]; k++)
			{
				fprintf(fp, "%.5f\t", Position_Knots[i][j][k]);
			}
			fprintf(fp, "\n");
		}
	}
	fprintf(fp,"\n");

	for (i = 0; i < Total_Control_Point_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp, "%d\t", i);
		for (j = 0; j < DIMENSION + 1; j++)
		{
			fprintf(fp, "%.10f\t", Node_Coordinate[i][j]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp,"\n");

	for (i = 0; i < Total_Constraint_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp, "%d\t%d\t%.10f\n",
				Constraint_Node_Dir[i][0],
				Constraint_Node_Dir[i][1],
				Value_of_Constraint[i]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < Total_Load_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp, "%d\t%d\t%.10f\n",
				Load_Node_Dir[i][0],
				Load_Node_Dir[i][1],
				Value_of_Load[i]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < Total_DistributeForce_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp, "%d\t%d\t%d\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n",
				type_load_array[i],
				iPatch_array[i],
				iCoord_array[i],
				val_Coord_array[i],
				Range_Coord_array[i][0],
				Range_Coord_array[i][1],
				Coeff_Dist_Load_array[i][0],
				Coeff_Dist_Load_array[i][1],
				Coeff_Dist_Load_array[i][2]);
	}

	fclose(fp);

	fp = fopen("coord_data.txt","w");
	for (i = 0; i < Total_Control_Point_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp, "%d\t", i);
		for (j = 0; j < DIMENSION; j++)
		{
			fprintf(fp, "%.10e\t", Node_Coordinate[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	fp = fopen("Displacement.dat", "w");
	fprintf(fp, "label=Displacement\n");
	fprintf(fp, "num_items=%d\n", Total_Control_Point_to_mesh[Total_mesh]);
	fprintf(fp, "\n");
	for (j = 0; j < Total_Control_Point_to_mesh[Total_mesh]; j++)
	{
		fprintf(fp, "%d:	%.16e %.16e ", j, Displacement[j * DIMENSION + 0], Displacement[j * DIMENSION + 1]);
		fprintf(fp, "\n");
	}
	fclose(fp);

	// 重ね合わせをスキップしてJ積分を行う
	if (SKIP_S_IGA != 1)
	{
		// 重ね合わせの結果
		division_ele_xi = DIVISION_ELE_XI;
		division_ele_eta = DIVISION_ELE_ETA;

		if (division_ele_xi > MAX_DIVISION)
		{
			printf("Error!!\n");
			printf("Too many Divsion at xi!\n"
				   "Maximum of division is %d (Now %d)\n"
				   "\n",
				   MAX_DIVISION, division_ele_xi);
			exit(1);
		}
		if (division_ele_eta > MAX_DIVISION)
		{
			printf("Error!!\n");
			printf("Too many Divsion at eta!\n"
				   "Maximum of division is %d (Now %d)\n"
				   "\n",
				   MAX_DIVISION, division_ele_eta);
			exit(1);
		}

		// ローカルメッシュの情報取得
		GetLocData();

		int patch_n_loc = 0, patch_n_glo = 0; // パッチ番号

		ReadFile();
		fp = fopen("view.dat", "w");
		fprintf(fp, "%d\t%d\t%d\n",
				fields_flag, division_ele_xi, division_ele_eta);
		fclose(fp);
		// machino
		fp = fopen("view_r_theta.dat", "w");
		fprintf(fp, "%d\t%d\t%d\n",
				fields_flag, division_ele_xi, division_ele_eta);
		fclose(fp);

		n_patch_glo = patch_n - n_patch_loc;

		// for s-IGA
		// 重ね合わせ結果出力のためのoverlay_view.dat作成
		fp = fopen("overlay_view.dat", "w");
		fprintf(fp, "%d\t%d\t%d\n",
				fields_flag, division_ele_xi, division_ele_eta);
		fclose(fp);
		// machino
		fp = fopen("overlay_view_r_theta.dat", "w");
		fprintf(fp, "%d\t%d\t%d\n",
				fields_flag, division_ele_xi, division_ele_eta);
		fclose(fp);

		// グラフ作成のための出力
		fp = fopen("disp_graph.txt", "w");
		fprintf(fp, "patch_n\tx\ty\tdisp_x\tdisp_y\n");
		fclose(fp);

		fp = fopen("stress_y_graph.txt", "w");
		fprintf(fp, "patch_n\tx\ty\tstress_yy\n");
		fclose(fp);

		fp = fopen("stress_y_graph_0.txt", "w");
		fprintf(fp, "patch_n\tx\ty\tstress_yy\n");
		fclose(fp);

		fp = fopen("stress_vm_graph.txt", "w");
		fprintf(fp, "xi\teta\tx\ty\tstress_vm\n");
		fclose(fp);

		fp = fopen("over_disp_graph.txt", "w");
		fprintf(fp, "patch_n\tx\ty\tdisp_x\tdisp_y\n");
		fclose(fp);

		fp = fopen("over_stress_x_graph.txt", "w");
		fprintf(fp, "x\ty\tstress_xx\n");
		fclose(fp);

		fp = fopen("over_stress_y_graph.txt", "w");
		fprintf(fp, "x\ty\tstress_yy\n");
		fclose(fp);

		fp = fopen("over_stress_y_graph_0.txt", "w");
		fprintf(fp, "x\ty\tstress_yy\n");
		fclose(fp);

		fp = fopen("over_stress_r_theta_graph.txt", "w");
		fprintf(fp, "xi\teta\tx\ty\tstress_r\tstress_theta\n");
		fclose(fp);

		fp = fopen("over_stress_vm_graph.txt", "w");
		fprintf(fp, "xi\teta\tx\ty\tstress_vm\n");
		fclose(fp);

		// 重ね合わせ
		for (i = 0; i < patch_n; i++)
		{

			fp = fopen("stress_vm_graph.txt", "a");
			fprintf(fp, "\npatch_n;%d\n\n", i);
			fclose(fp);

			graph_patch_n = i;

			printf("----------Start calculation at patch %d----------\n\n", i);
			Calculation(order_xi[i], order_eta[i],
						knot_n_xi[i], knot_n_eta[i],
						cntl_p_n_xi[i], cntl_p_n_eta[i],
						knot_vec_xi[i], knot_vec_eta[i],
						cntl_px[i], cntl_py[i],
						disp_cntl_px[i], disp_cntl_py[i],
						weight[i]);
			printf("-----------End calculation at patch %d-----------\n\n", i);

			if (i >= n_patch_glo) // ローカル上のパッチに対しては重合計算行う
			{
				patch_n_loc = i;
				printf("----------Start overlay calculation at patch %d in LOCAL patch----------\n\n", i);
				for (j = 0; j < n_patch_glo; j++)
				{
					patch_n_glo = j;
					Calculation_overlay(order_xi[patch_n_loc], order_eta[patch_n_loc],
										knot_n_xi[patch_n_loc], knot_n_eta[patch_n_loc],
										cntl_p_n_xi[patch_n_loc], cntl_p_n_eta[patch_n_loc],
										knot_vec_xi[patch_n_loc], knot_vec_eta[patch_n_loc],
										cntl_px[patch_n_loc], cntl_py[patch_n_loc],
										weight[patch_n_loc],
										order_xi[patch_n_glo], order_eta[patch_n_glo],
										cntl_p_n_xi[patch_n_glo], cntl_p_n_eta[patch_n_glo],
										knot_vec_xi[patch_n_glo], knot_vec_eta[patch_n_glo],
										cntl_px[patch_n_glo], cntl_py[patch_n_glo],
										disp_cntl_px[patch_n_glo], disp_cntl_py[patch_n_glo],
										weight[patch_n_glo]);
				}
			}
		}
		printf("End S-IGA\n\n");
	}
	else if (SKIP_S_IGA == 1)
	{
		printf("SKIP_S_IGA = 1, skip S-IGA\n");
	}

	if (SKIP_S_IGA == 2) {
		printf("SKIP_S_IGA = 2, exit before J integration\n");
		exit(0);
	}

	// 重ね合わせをスキップした場合ここから
	printf("Start J Integration Mixed Mode\n\n");

	// For the J-integral Evaluation

	free(INC), free(Position_Knots), free(Total_Knot_to_patch_dim);
	return 0;
}