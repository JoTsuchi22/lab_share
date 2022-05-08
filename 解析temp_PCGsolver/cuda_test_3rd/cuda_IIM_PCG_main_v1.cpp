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
			 Order, No_knot, Total_Knot_to_mesh, Node_Coordinate,
			 Position_Knots, No_Control_point_ON_ELEMENT);

	// memory free
	free(Patch_Control_point), free(real_element_line), free(Total_element_all_ID), free(difference);

	// memory allocation
	if (DIMENSION == 2)
	{
		D_MATRIX_SIZE = 3;
	}
	else if (DIMENSION == 3)
	{
		D_MATRIX_SIZE = 6;
	}
	int *NNLOVER = (int *)malloc(sizeof(int) * real_Total_Element_to_mesh[Total_mesh]);
	int *NELOVER = (int *)malloc(sizeof(int) * real_Total_Element_to_mesh[Total_mesh] * MAX_N_ELEMENT_OVER);
	double *Gauss_Coordinate = (double *)calloc(real_Total_Element_to_mesh[Total_mesh] * POW_Ng * DIMENSION, sizeof(double));
	double *Gauss_Coordinate_ex = (double *)calloc(real_Total_Element_to_mesh[Total_mesh] * POW_Ng_extended * DIMENSION, sizeof(double));
	double *Jac = (double *)malloc(sizeof(double) * real_Total_Element_to_mesh[Total_mesh] * POW_Ng);
	double *Jac_ex = (double *)malloc(sizeof(double) * real_Total_Element_to_mesh[Total_mesh] * POW_Ng_extended);
	double *B_Matrix = (double *)malloc(sizeof(double) * real_Total_Element_to_mesh[Total_mesh] * POW_Ng * D_MATRIX_SIZE * MAX_KIEL_SIZE);
	double *B_Matrix_ex = (double *)malloc(sizeof(double) * real_Total_Element_to_mesh[Total_mesh] * POW_Ng_extended * D_MATRIX_SIZE * MAX_KIEL_SIZE);
	double *Loc_parameter_on_Glo = (double *)malloc(real_Total_Element_to_mesh[Total_mesh] * POW_Ng * DIMENSION * sizeof(double));
	double *Loc_parameter_on_Glo_ex = (double *)malloc(real_Total_Element_to_mesh[Total_mesh] * POW_Ng_extended * DIMENSION * sizeof(double));
	
	// check_over_parameter, NNLOVER の算出
	for (i = 1; i < Total_mesh; i++)
	{
		int mesh_n_org = 0;
		printf("mesh_n_org: 0\tmesh_n_over: %d\n", i);
		Check_coupled_Glo_Loc_element_for_Gauss(i, mesh_n_org, NNLOVER, NELOVER,
												Gauss_Coordinate, Gauss_Coordinate_ex, Jac, Jac_ex, B_Matrix, B_Matrix_ex, Loc_parameter_on_Glo, Loc_parameter_on_Glo_ex,
												real_Total_Element_to_mesh, Node_Coordinate, Total_Control_Point_to_mesh, Controlpoint_of_Element,
												INC, Element_patch, Order, No_Control_point_ON_ELEMENT, Position_Knots,
												Total_Knot_to_patch_dim, No_knot, No_Control_point,
												Control_Coord_x, Control_Coord_y, Control_Weight,
												real_Total_Element_on_mesh, real_element, Total_Patch_on_mesh, line_No_Total_element);
		Make_Loc_Glo(real_Total_Element_on_mesh, real_Total_Element_to_mesh, real_element, NNLOVER, NELOVER);
	}


	return 0;
}