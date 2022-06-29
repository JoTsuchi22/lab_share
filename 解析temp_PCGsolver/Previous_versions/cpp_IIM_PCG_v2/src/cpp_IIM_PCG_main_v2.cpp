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
	int i, tm;
	clock_t start, end;
	information info, *info_ptr;
	info_ptr = &info;

	Total_mesh = argc - 1;

	// 引数の個数確認
	if (argc <= 1)
	{
		printf("Argument is missing\n");
	}
	else if (argc == 2) // 通常IGA: input file 1つ
	{
		printf("IGA carried out.(No local mesh)\n");
	}
	else if (argc >= 3) // s-IGA: input file 複数
	{
		printf("s-IGA carried out.(%d local meshes)\n", argc - 2);
	}

	start = clock();

	// memory allocation
	info_ptr->Total_Knot_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));
	info_ptr->Total_Patch_on_mesh = (int *)malloc(sizeof(int) * (Total_mesh));				// 各メッシュ上のパッチ数
	info_ptr->Total_Patch_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));			// メッシュまでのパッチ数（メッシュ内のパッチ数は含まない）
	info_ptr->Total_Control_Point_on_mesh = (int *)malloc(sizeof(int) * (Total_mesh));		// 各メッシュ上のコントロールポイント数
	info_ptr->Total_Control_Point_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));	// メッシュまでのコントロールポイント数(メッシュ内のコントロールポイント数は含まない)
	info_ptr->Total_Element_on_mesh = (int *)malloc(sizeof(int) * (Total_mesh));			// 各メッシュ上の要素数
	info_ptr->Total_Element_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));			// メッシュまでの要素数(メッシュ内の要素数は含まない)
	info_ptr->real_Total_Element_on_mesh = (int *)malloc(sizeof(int) * (Total_mesh));		// real_Total_Element_on_mesh[MAX_N_MESH]
	info_ptr->real_Total_Element_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));	// real_Total_Element_to_mesh[MAX_N_MESH + 1]
	info_ptr->Total_Load_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));
	info_ptr->Total_Constraint_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));
	info_ptr->Total_DistributeForce_to_mesh = (int *)calloc((Total_mesh + 1), sizeof(int));
	if (info.Total_Knot_to_mesh == NULL || info.Total_Patch_on_mesh == NULL || info.Total_Patch_to_mesh == NULL || info.Total_Control_Point_on_mesh == NULL || info.Total_Control_Point_to_mesh == NULL || info.Total_Element_on_mesh == NULL || info.Total_Element_to_mesh == NULL || info.real_Total_Element_on_mesh == NULL || info.real_Total_Element_to_mesh == NULL || info.Total_Load_to_mesh == NULL || info.Total_Constraint_to_mesh == NULL || info.Total_DistributeForce_to_mesh == NULL)
    {
        printf("Cannot allocate memory\n"); exit(1);
    }

	// ファイル読み込み1回目
    for (tm = 0; tm < Total_mesh; tm++)
    {
		Get_Input_1(tm, argv, &info);
	}
	
	// memory allocation
	info_ptr->Order = (int *)malloc(sizeof(int) * (info.Total_Patch_to_mesh[Total_mesh] * DIMENSION));					// Order[MAX_N_PATCH][DIMENSION]
	info_ptr->No_knot = (int *)malloc(sizeof(int) * (info.Total_Patch_to_mesh[Total_mesh] * DIMENSION));				// No_knot[MAX_N_PATCH][DIMENSION]
	info_ptr->Total_Control_Point_to_patch = (int *)calloc((info.Total_Patch_to_mesh[Total_mesh] + 1), sizeof(int));	// Total_Control_Point_to_patch[MAX_N_PATCH]
	info_ptr->Total_Knot_to_patch_dim = (int *)calloc((info.Total_Patch_to_mesh[Total_mesh] * DIMENSION + 1), sizeof(int));	// Total_Knot_to_patch_dim[MAX_N_PATCH][DIMENSION]
	info_ptr->Position_Knots = (double *)malloc(sizeof(double) * info.Total_Knot_to_mesh[Total_mesh]);					// Position_Knots[MAX_N_PATCH][DIMENSION][MAX_N_KNOT];
	info_ptr->No_Control_point = (int *)malloc(sizeof(int) * (info.Total_Patch_to_mesh[Total_mesh] * DIMENSION));		// Order[MAX_N_PATCH][DIMENSION]
	info_ptr->No_Control_point_in_patch = (int *)malloc(sizeof(int) * (info.Total_Patch_to_mesh[Total_mesh]));			// No_Control_point_in_patch[MAX_N_PATCH]
	info_ptr->Patch_Control_point = (int *)malloc(sizeof(int) * (info.Total_Control_Point_to_mesh[Total_mesh]));		// Patch_Control_point[MAX_N_PATCH][MAX_N_Controlpoint_in_Patch]
	info_ptr->No_Control_point_ON_ELEMENT = (int *)malloc(sizeof(int) * (info.Total_Patch_to_mesh[Total_mesh]));		// No_Control_point_ON_ELEMENT[MAX_N_PATCH]
	info_ptr->Node_Coordinate = (double *)malloc(sizeof(double) * (info.Total_Control_Point_to_mesh[Total_mesh] * (DIMENSION + 1)));	// Node_Coordinate[MAX_N_NODE][DIMENSION + 1];
	info_ptr->Control_Coord_x = (double *)malloc(sizeof(double) * (info.Total_Control_Point_to_mesh[Total_mesh]));		// Control_Coord[DIMENSION][MAX_N_NODE];
	info_ptr->Control_Coord_y = (double *)malloc(sizeof(double) * (info.Total_Control_Point_to_mesh[Total_mesh]));		// Control_Coord[DIMENSION][MAX_N_NODE];
	info_ptr->Control_Weight = (double *)malloc(sizeof(double) * (info.Total_Control_Point_to_mesh[Total_mesh]));		// Control_Weight[MAX_N_NODE];
	info_ptr->Constraint_Node_Dir = (int *)malloc(sizeof(int) * (info.Total_Constraint_to_mesh[Total_mesh] * 2));		// Constraint_Node_Dir[MAX_N_CONSTRAINT][2];
	info_ptr->Value_of_Constraint = (double *)calloc(info.Total_Constraint_to_mesh[Total_mesh], sizeof(double));		// Value_of_Constraint[MAX_N_CONSTRAINT];
	info_ptr->Load_Node_Dir = (int *)malloc(sizeof(int) * (info.Total_Load_to_mesh[Total_mesh] * 2));					// Load_Node_Dir[MAX_N_LOAD][2];
	info_ptr->Value_of_Load = (double *)calloc(info.Total_Load_to_mesh[Total_mesh], sizeof(double));					// Value_of_Load[MAX_N_LOAD];
	info_ptr->iPatch_array = (int *)malloc(sizeof(int) * info.Total_DistributeForce_to_mesh[Total_mesh]);				// iPatch_array[MAX_N_DISTRIBUTE_FORCE]
	info_ptr->iCoord_array = (int *)malloc(sizeof(int) * info.Total_DistributeForce_to_mesh[Total_mesh]);				// iCoord_array[MAX_N_DISTRIBUTE_FORCE]
	info_ptr->type_load_array = (int *)malloc(sizeof(int) * info.Total_DistributeForce_to_mesh[Total_mesh]);			// type_load_array[MAX_N_DISTRIBUTE_FORCE]
	info_ptr->val_Coord_array = (double *)calloc(info.Total_DistributeForce_to_mesh[Total_mesh], sizeof(double));		// val_Coord_array[MAX_N_DISTRIBUTE_FORCE]
	info_ptr->Range_Coord_array = (double *)calloc((info.Total_DistributeForce_to_mesh[Total_mesh] * 2), sizeof(double));	// Range_Coord_array[MAX_N_DISTRIBUTE_FORCE][2]
	info_ptr->Coeff_Dist_Load_array = (double *)calloc((info.Total_DistributeForce_to_mesh[Total_mesh] * 3), sizeof(double));	// Coeff_Dist_Load_array[MAX_N_DISTRIBUTE_FORCE][3] 
	if (info.Order == NULL || info.No_knot == NULL || info.Total_Control_Point_to_patch == NULL || info.Total_Knot_to_patch_dim == NULL || info.Position_Knots == NULL || info.No_Control_point == NULL || info.No_Control_point_in_patch == NULL || info.Patch_Control_point == NULL || info.No_Control_point_ON_ELEMENT == NULL || info.Node_Coordinate == NULL || info.Control_Coord_x == NULL || info.Control_Coord_y == NULL || info.Control_Weight == NULL || info.Constraint_Node_Dir == NULL || info.Value_of_Constraint == NULL || info.Load_Node_Dir == NULL || info.Value_of_Load == NULL || info.iPatch_array == NULL || info.iCoord_array == NULL || info.type_load_array == NULL || info.val_Coord_array == NULL || info.Range_Coord_array == NULL || info.Coeff_Dist_Load_array == NULL)
	{
        printf("Cannot allocate memory\n"); exit(1);
    }

	// ファイル読み込み2回目
	for (tm = 0; tm < Total_mesh; tm++)
    {
		Get_Input_2(tm, argv, &info);
	}

	printf("\nFinish Get Input data\n\n");

	// memory allocation
	info_ptr->INC = (int *)malloc(sizeof(int) * (info.Total_Control_Point_to_mesh[Total_mesh] * DIMENSION));			// INC[MAX_N_NODE][DIMENSION]
	info_ptr->Controlpoint_of_Element = (int *)malloc(sizeof(int) * (info.Total_Element_to_mesh[Total_mesh] * MAX_NO_CCpoint_ON_ELEMENT)); // Controlpoint_of_Element[MAX_N_ELEMENT][MAX_NO_CCpoint_ON_ELEMENT]
	info_ptr->Element_patch = (int *)malloc(sizeof(int) * info.Total_Element_to_mesh[Total_mesh]);						// Element_patch[MAX_N_ELEMENT]
	info_ptr->Element_mesh = (int *)malloc(sizeof(int) * info.Total_Element_to_mesh[Total_mesh]);						// Element_mesh[MAX_N_ELEMENT] 要素がどのメッシュ内にあるかを示す配列
	info_ptr->line_No_real_element = (int *)malloc(sizeof(int) * (info.Total_Patch_to_mesh[Total_mesh] * DIMENSION));	// line_No_real_element[MAX_N_PATCH][DIMENSION] ゼロエレメントではない要素列の数
	info_ptr->line_No_Total_element = (int *)malloc(sizeof(int) * (info.Total_Patch_to_mesh[Total_mesh] * DIMENSION));	// line_No_Total_element[MAX_N_PATCH][DIMENSION] ゼロエレメントを含むすべての要素列の数
	info_ptr->difference = (double *)calloc(info.Total_Knot_to_mesh[Total_mesh], sizeof(double));						// difference[MAX_N_PATCH][DIMENSION][MAX_N_KNOT] 隣り合うノットベクトルの差
	info_ptr->Total_element_all_ID = (int *)calloc(info.Total_Element_to_mesh[Total_mesh], sizeof(int));				// Total_element_all_ID[MAX_N_ELEMENT] ゼロエレメントではない要素 = 1, ゼロエレメント = 0
	info_ptr->ENC = (int *)malloc(sizeof(int) * (info.Total_Element_to_mesh[Total_mesh] * DIMENSION));					// ENC[MAX_N_ELEMENT][DIMENSION] ENC[パッチ][全ての要素][0, 1] = x, y方向の何番目の要素か
	info_ptr->real_element_line = (int *)malloc(sizeof(int) * (info.Total_Patch_to_mesh[Total_mesh] * info.Total_Element_to_mesh[Total_mesh] * DIMENSION)); // real_element_line[MAX_N_PATCH][MAX_N_ELEMENT][DIMENSION] ゼロエレメントではない要素列
	info_ptr->real_element = (int *)malloc(sizeof(int) * info.Total_Element_to_mesh[Total_mesh]);							// real_element[MAX_N_ELEMENT] ゼロエレメントではない要素の番号
	info_ptr->real_El_No_on_mesh = (int *)malloc(sizeof(int) * (info.Total_Patch_to_mesh[Total_mesh] * info.Total_Element_to_mesh[Total_mesh]));	// real_El_No_on_mesh[MAX_N_MESH][MAX_N_ELEMENT]
	info_ptr->Equivalent_Nodal_Force = (double *)calloc(info.Total_Control_Point_to_mesh[Total_mesh] * DIMENSION, sizeof(double));	// Equivalent_Nodal_Force[MAX_N_NODE][DIMENSION] Equivalent nodal forces arising from the distributed load
	if (info.INC == NULL || info.Controlpoint_of_Element == NULL || info.Element_patch == NULL || info.Element_mesh == NULL || info.line_No_real_element == NULL || info.line_No_Total_element == NULL || info.difference == NULL || info.Total_element_all_ID == NULL || info.ENC == NULL || info.real_element_line == NULL || info.real_element == NULL || info.real_El_No_on_mesh == NULL || info.Equivalent_Nodal_Force == NULL)
	{
		printf("Cannot allocate memory\n"); exit(1);
	}

	// INC 等の作成
	Make_INC(tm, &info);

	// memory free
	free(info.real_element_line), free(info.Total_element_all_ID), free(info.difference);

	// memory allocation
	if (DIMENSION == 2)
	{
		D_MATRIX_SIZE = 3;
	}
	else if (DIMENSION == 3)
	{
		D_MATRIX_SIZE = 6;
	}
	info_ptr->NNLOVER = (int *)malloc(sizeof(int) * info.real_Total_Element_to_mesh[Total_mesh]);
	info_ptr->NELOVER = (int *)malloc(sizeof(int) * info.real_Total_Element_to_mesh[Total_mesh] * MAX_N_ELEMENT_OVER);
	info_ptr->Gauss_Coordinate = (double *)calloc(info.real_Total_Element_to_mesh[Total_mesh] * POW_NG * DIMENSION, sizeof(double));
	info_ptr->Gauss_Coordinate_ex = (double *)calloc(info.real_Total_Element_to_mesh[Total_mesh] * POW_NG_EXTEND * DIMENSION, sizeof(double));
	info_ptr->Jac = (double *)malloc(sizeof(double) * info.real_Total_Element_to_mesh[Total_mesh] * POW_NG);
	info_ptr->Jac_ex = (double *)malloc(sizeof(double) * info.real_Total_Element_to_mesh[Total_mesh] * POW_NG_EXTEND);
	info_ptr->B_Matrix = (double *)malloc(sizeof(double) * info.real_Total_Element_to_mesh[Total_mesh] * POW_NG * D_MATRIX_SIZE * MAX_KIEL_SIZE);
	info_ptr->B_Matrix_ex = (double *)malloc(sizeof(double) * info.real_Total_Element_to_mesh[Total_mesh] * POW_NG_EXTEND * D_MATRIX_SIZE * MAX_KIEL_SIZE);
	info_ptr->Loc_parameter_on_Glo = (double *)malloc(info.real_Total_Element_to_mesh[Total_mesh] * POW_NG * DIMENSION * sizeof(double));
	info_ptr->Loc_parameter_on_Glo_ex = (double *)malloc(info.real_Total_Element_to_mesh[Total_mesh] * POW_NG_EXTEND * DIMENSION * sizeof(double));
	if (info.NNLOVER == NULL || info.NELOVER == NULL || info.Gauss_Coordinate == NULL || info.Gauss_Coordinate_ex == NULL || info.Jac == NULL || info.Jac_ex == NULL || info.B_Matrix == NULL || info.B_Matrix_ex == NULL || info.Loc_parameter_on_Glo == NULL || info.Loc_parameter_on_Glo_ex == NULL)
	{
		printf("Cannot allocate memory\n"); exit(1);
	}

	// check_over_parameter, NNLOVER の算出
	for (i = 1; i < Total_mesh; i++)
	{
		int mesh_n_org = 0;
		Check_coupled_Glo_Loc_element_for_Gauss(i, mesh_n_org, &info);
		Make_Loc_Glo(&info);
	}
	printf("\nFinish check_over_parameter\n\n");

	// memory allocation
	MAX_K_WHOLE_SIZE = info.Total_Control_Point_to_mesh[Total_mesh] * DIMENSION;
	MAX_NON_ZERO = MAX_K_WHOLE_SIZE * (MAX_K_WHOLE_SIZE / 2) + MAX_K_WHOLE_SIZE;
	info_ptr->D = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * D_MATRIX_SIZE);
	info_ptr->Node_To_Node = (int *)malloc(sizeof(int) * K_DIVISION_LENGE * info.Total_Control_Point_to_mesh[Total_mesh]);	// Node_To_Node[K_DIVISION_LENGE][10000]
	info_ptr->Total_Control_Point_To_Node = (int *)malloc(sizeof(int) * K_DIVISION_LENGE);									// Total_Control_Point_To_Node[K_DIVISION_LENGE] ある節点に関係する節点番号
	info_ptr->Index_Dof = (int *)calloc(MAX_K_WHOLE_SIZE, sizeof(int));			// Index_Dof[MAX_K_WHOLE_SIZE];
	info_ptr->K_Whole_Ptr = (int *)calloc(MAX_K_WHOLE_SIZE + 1, sizeof(int));	// K_Whole_Ptr[MAX_K_WHOLE_SIZE + 1]
	info_ptr->K_Whole_Col = (int *)malloc(sizeof(int) * MAX_NON_ZERO);			// K_Whole_Col[MAX_NON_ZERO]
	info_ptr->K_Whole_Val = (double *)calloc(MAX_NON_ZERO, sizeof(double));		// K_Whole_Val[MAX_NON_ZERO]
	if (info.D == NULL || info.Node_To_Node == NULL || info.Total_Control_Point_To_Node == NULL || info.Index_Dof == NULL || info.K_Whole_Ptr == NULL || info.K_Whole_Col == NULL || info.K_Whole_Val == NULL)
	{
        printf("Cannot allocate memory\n"); exit(1);
    }

    // 全体剛性マトリックスの制作
	Make_D_Matrix(&info);
	Make_Index_Dof(&info);
	Make_K_Whole_Ptr_Col(&info);
    Make_K_Whole_Val(&info);
    printf("\nFinish Make_K_Whole\n\n");

	// memory free
	free(info.ENC), free(info.Node_To_Node), free(info.Total_Control_Point_To_Node);

	// memory allocation
	info_ptr->sol_vec = (double *)calloc(MAX_K_WHOLE_SIZE, sizeof(int));	// sol_vec[MAX_K_WHOLE_SIZE]
	info_ptr->rhs_vec = (double *)calloc(MAX_K_WHOLE_SIZE, sizeof(int));	// rhs_vec[MAX_K_WHOLE_SIZE]
	if (info.sol_vec == NULL ||	info.rhs_vec == NULL)
	{
		printf("Cannot allocate memory\n"); exit(1);
	}

	// 荷重ベクトルの制作
	Make_F_Vec(&info);
	Make_F_Vec_disp_const(&info);
	Add_Equivalent_Nodal_Force_to_F_Vec(&info);
    printf("\nFinish Make_F_Vec\n\n");

	// memory free
	free(info.Equivalent_Nodal_Force);

	// PCG法
    printf("\nStart PCG solver\n\n");
    int max_itr = K_Whole_Size;
	PCG_Solver(max_itr, EPS, &info);
	printf("\nFinish PCG solver\n\n");

	// memory free

	// memory allocation
    info_ptr->Displacement = (double *)malloc(sizeof(double) * MAX_K_WHOLE_SIZE);	// Displacement[MAX_K_WHOLE_SIZE]
	if (info.Displacement == NULL)
	{
		printf("Cannot allocate memory\n"); exit(1);
	}

	// Postprocessing
	Make_Displacement(&info);
	printf("\nFinish Make_Displacement\n\n");
	end = clock();
	printf("\nAnalysis time:%.2f[s]\n\n", (double)(end - start) / CLOCKS_PER_SEC);
	// Make_Strain(real_Total_Element_to_mesh[Total_mesh]);
	printf("\nFinish Make_Strain\n\n");
	// Make_Stress_2D(E, nu, real_Total_Element_to_mesh[Total_mesh], DM);
	printf("\nFinish Make_Stress\n\n");
	// Make_ReactionForce(Total_Control_Point_to_mesh[Total_mesh]);
	printf("\nFinish Make_ReactionForce\n\n");
	// Make_Parameter_z(real_Total_Element_to_mesh[Total_mesh], E, nu, DM);
	printf("\nFinish Make_Parameter_z\n\n");

	// Kマトリックスのsvg出力
	if (OUTPUT_SVG == 1)
	{
		K_output_svg(&info);
		printf("\nFinish K_output_svg\n\n");
	}

	// viewer のための出力
	output_for_viewer(&info);
	printf("\nFinish output_for_viewer\n\n");

	// 重ね合わせをスキップしてJ積分を行う
	if (SKIP_S_IGA != 1)
	{
		printf("\nStart s_IGA overlay\n\n");
		start = clock();

		// memory allocation
		DIVISION_ELEMENT = (DIVISION_ELE + 1) * (DIVISION_ELE + 1);
		info_ptr->coord_x = (double *)malloc(sizeof(double) * 2 * info.Total_Element_to_mesh[Total_mesh] * DIVISION_ELEMENT);		// coord_x[MAX_POINTS][MAX_POINTS]
		info_ptr->coord_y = (double *)malloc(sizeof(double) * 2 * info.Total_Element_to_mesh[Total_mesh] * DIVISION_ELEMENT);		// coord_y[MAX_POINTS][MAX_POINTS]
		info_ptr->disp_x = (double *)malloc(sizeof(double) * 2 * info.Total_Element_to_mesh[Total_mesh] * DIVISION_ELEMENT);		// disp_x[MAX_POINTS][MAX_POINTS]
		info_ptr->disp_y = (double *)malloc(sizeof(double) * 2 * info.Total_Element_to_mesh[Total_mesh] * DIVISION_ELEMENT);		// disp_y[MAX_POINTS][MAX_POINTS]
		info_ptr->strain = (double *)malloc(sizeof(double) * 2 * info.Total_Element_to_mesh[Total_mesh] * DIVISION_ELEMENT * 3);	// strain[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1][3] // 0: xx, 1: yy, 2: xy
		info_ptr->stress = (double *)malloc(sizeof(double) * 2 * info.Total_Element_to_mesh[Total_mesh] * DIVISION_ELEMENT * 3);	// stress[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1][3]
		if (info.coord_x == NULL || info.coord_y == NULL || info.disp_x == NULL || info.disp_y == NULL || info.strain == NULL || info.stress == NULL) 
		{
			printf("Cannot allocate memory\n"); exit(1);
		}

		s_IGA_overlay(&info);
		end = clock();
		printf("\noverlay time:%.2f[s]\n\n", (double)(end - start) / CLOCKS_PER_SEC);
		printf("\nFinish s_IGA overlay\n\n");
	}
	else if (SKIP_S_IGA == 1)
	{
		printf("SKIP_S_IGA = 1, skip S-IGA\n");
	}

	if (SKIP_S_IGA == 2)
	{
		printf("SKIP_S_IGA = 2, exit before J integration\n");
		exit(0);
	}

	// 重ね合わせをスキップした場合ここから
	printf("Start J Integration Mixed Mode\n\n");

	// For the J-integral Evaluation

	// memory free
	free(info.INC), free(info.Position_Knots), free(info.Total_Knot_to_patch_dim);

	return 0;
}