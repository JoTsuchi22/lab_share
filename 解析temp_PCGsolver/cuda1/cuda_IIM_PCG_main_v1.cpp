/****************************************************************
 * s-IGA :)
 * 
 * 仮定条件
 * 	ローカルメッシュ同士は原則被りなしと仮定
 * 	グローバルパッチはシングルパッチ
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
	clock_t start, end, t1;

	int i, j, k;
	int re;
	int Load_Node_Dir[MAX_N_LOAD][2];
	double Value_of_Load[MAX_N_LOAD];
	int Total_Constraint = 0, Constraint_Node_Dir[MAX_N_CONSTRAINT][2];
	double Value_of_Constraint[MAX_N_CONSTRAINT];
	int Total_DistributeForce = 0;
	int K_Whole_Size = 0;
	int El_No = 0;
    double element_loc[DIMENSION];

	int Total_mesh = argc - 1;

    // for s-IGA
    int tm;

    // 引数の個数確認
	if (argc <= 1)
	{
		printf("Argument is missing\n");
	}
    if (argc == 2)  // 通常IGA: input file 1つ
    {
        printf("IGA carried out.(No local mesh)\n");
    }
    if (argc >= 3)  // s-IGA: input file 複数
    {
        printf("s-IGA carried out.(%d local meshes)\n", argc - 2);
    }

	start = clock();

	// memory allocation
	int *Total_Patch_on_mesh = (int *)malloc(sizeof(int) * (Total_mesh));				// 各メッシュ上のパッチ数
	int *Total_Patch_to_mesh = (int *)malloc(sizeof(int) * (Total_mesh + 1));			// メッシュ[]までのパッチ数（メッシュ[]内のパッチ数は含まない）
	int *Total_Control_Point_on_mesh = (int *)malloc(sizeof(int) * (Total_mesh));		// 各メッシュ上のコントロールポイント数
	int *Total_Control_Point_to_mesh = (int *)malloc(sizeof(int) * (Total_mesh + 1));	// メッシュ[]までのコントロールポイント数(メッシュ[]内のコントロールポイント数は含まない)
	int *Total_Element_on_mesh = (int *)malloc(sizeof(int) * (Total_mesh));				// 各メッシュ上の要素数
	int *Total_Element_to_mesh = (int *)malloc(sizeof(int) * (Total_mesh + 1));			// メッシュ[]までの要素数(メッシュ[]内の要素数は含まない)
	int *Total_Constraint_to_mesh = (int *)malloc(sizeof(int) * (Total_mesh + 1));
	int *Total_Load_to_mesh = (int *)malloc(sizeof(int) * (Total_mesh + 1));
	int *Total_DistributeForce_to_mesh = (int *)malloc(sizeof(int) * (Total_mesh + 1));

	// ファイル読み込み1回目
    for (tm = 0; tm < Total_mesh; tm++)
    {
		Get_input_1(tm, E, nu, Total_Patch_on_mesh, Total_Patch_to_mesh, Total_Control_Point_on_mesh, Total_Control_Point_to_mesh,
					Total_Element_on_mesh, Total_Element_to_mesh);
	}
	
	// memory allocation
	int *Order = (int *)malloc(sizeof(int) * (Total_Patch_to_mesh[Total_mesh] * DIMENSION));				// Order[MAX_N_PATCH][DIMENSION]
	int *No_knot = (int *)malloc(sizeof(int) * (Total_Patch_to_mesh[Total_mesh] * DIMENSION));				// No_knot[MAX_N_PATCH][DIMENSION]
	int *No_Control_point = (int *)malloc(sizeof(int) * (Total_Patch_to_mesh[Total_mesh] * DIMENSION));		// Order[MAX_N_PATCH][DIMENSION]
	int *No_Controlpoint_in_patch = (int *)malloc(sizeof(int) * (Total_Patch_to_mesh[Total_mesh]));			// No_Controlpoint_in_patch[MAX_N_PATCH]
	int *Patch_controlpoint = (int *)malloc(sizeof(int) * (Total_Control_Point_to_mesh[Total_mesh]));		// Patch_controlpoint[MAX_N_PATCH][MAX_N_Controlpoint_in_Patch]
	int *No_Control_point_ON_ELEMENT = (int *)malloc(sizeof(int) * (Total_Patch_to_mesh[Total_mesh]));		// No_Control_point_ON_ELEMENT[MAX_N_PATCH]
	double *Node_Coordinate = (double *)malloc(sizeof(double) * (Total_Control_Point_to_mesh[Total_mesh] * (DIMENSION + 1));	// Node_Coordinate[MAX_N_NODE][DIMENSION + 1];
	double *Control_Coord_x = (double *)malloc(sizeof(double) * (Total_Control_Point_to_mesh[Total_mesh])); //Control_Coord[DIMENSION][MAX_N_NODE];
	double *Control_Coord_y = (double *)malloc(sizeof(double) * (Total_Control_Point_to_mesh[Total_mesh])); //Control_Coord[DIMENSION][MAX_N_NODE];
	double *Control_Weight = (double *)malloc(sizeof(double) * (Total_Control_Point_to_mesh[Total_mesh]));	//Control_Weight[MAX_N_NODE];

	// ファイル読み込み2回目
	for (tm = 0; tm < Total_mesh; tm++)
    {
	    Get_input_2(tm, Load_Node_Dir,
		Value_of_Load, &Total_Constraint, Constraint_Node_Dir, Value_of_Constraint, &Total_DistributeForce, argv);
	}

	// memory allocation

	// ファイル読み込み3回目
	for (tm = 0; tm < Total_mesh; tm++)
    {
	    Get_input_3(tm, Load_Node_Dir,
		Value_of_Load, &Total_Constraint, Constraint_Node_Dir, Value_of_Constraint, &Total_DistributeForce, argv);
	}

	// memory free
	free(Patch_controlpoint);

	// ncheck_over_parameter
	printf("\ncheck_over_parameter;%d\n\n", check_over_parameter);
	for (i = 1; i < Total_mesh; i++)
	{
		printf("mesh_n_org;0\tmesh_n_over;%d\n", i);
		// NNLOVER[over_ele][]=org_eleの算出
		if (check_over_parameter == 0)
		{
			Check_coupled_Glo_Loc_element_for_end(element_loc, i, 0);
		}
		if (check_over_parameter == 1)
		{
			Check_coupled_Glo_Loc_element_for_Gauss(element_loc, i, 0);
		}
		// NNLOVER[org_ele][]=over_eleの算出
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
		for(j = 0; j < DIMENSION; j++)
		{
            int index = Index_Dof[i * DIMENSION + j];
			if (index >= 0)
				Displacement[i * DIMENSION + j] = sol_vec[index];
		}
	}
	printf("Finish Make_Displacement\n");
	end = clock();
	printf("Analysis time:%.2f[s]\n",(double)(end-start)/CLOCKS_PER_SEC);
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
		for (j = 0; j < No_Controlpoint_in_patch[i + g_patch]; j++)
		{
			l_temp = Patch_controlpoint[i + g_patch][j] - g_cntl_p;
			fprintf(fp,"%d\t", l_temp);
		}
		fprintf(fp,"\n");
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
				fprintf(fp,"%.5f\t",Position_Knots[i + g_patch][j][k]);
			}
			fprintf(fp,"\n");
		}
	}
	fprintf(fp,"\n");

	for (i = 0; i < l_cntl_p; i++)
	{
		fprintf(fp,"%d\t",i);
		for (j = 0; j < DIMENSION+1; j++)
		{
			fprintf(fp,"%.10f\t",Node_Coordinate[i + g_cntl_p][j]);
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n");

	for (i = 0; i < l_cnst; i++)
	{
		l_temp = Constraint_Node_Dir[i + g_cnst][0] - g_cntl_p;
		fprintf(fp,"%d\t%d\t%.10f\n",
					l_temp,
					Constraint_Node_Dir[i + g_cnst][1],
					Value_of_Constraint[i + g_cnst]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < l_load; i++)
	{
		l_temp = Load_Node_Dir[i + g_load][0] - g_cntl_p;
		fprintf(fp,"%d\t%d\t%.10f\n",
					l_temp,
					Load_Node_Dir[i + g_load][1],
					Value_of_Load[i + g_load]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < l_dist; i++)
	{
		fprintf(fp,"%d\t%d\t%d\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n",
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
		for (j = 0; j < No_Controlpoint_in_patch[i]; j++)
		{
			fprintf(fp,"%d\t",Patch_controlpoint[i][j]);
		}
		fprintf(fp,"\n");
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
				fprintf(fp,"%.5f\t",Position_Knots[i][j][k]);
			}
			fprintf(fp,"\n");
		}
	}
	fprintf(fp,"\n");

	for (i = 0; i < Total_Control_Point_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp,"%d\t",i);
		for (j = 0; j < DIMENSION+1; j++)
		{
			fprintf(fp,"%.10f\t",Node_Coordinate[i][j]);
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n");

	for (i = 0; i < Total_Constraint_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp,"%d\t%d\t%.10f\n",
					Constraint_Node_Dir[i][0],
					Constraint_Node_Dir[i][1],
					Value_of_Constraint[i]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < Total_Load_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp,"%d\t%d\t%.10f\n",
					Load_Node_Dir[i][0],
					Load_Node_Dir[i][1],
					Value_of_Load[i]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < Total_DistributeForce_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp,"%d\t%d\t%d\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n",
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
		fprintf(fp,"%d\t",i);
		for (j = 0; j < DIMENSION; j++)
		{
			fprintf(fp,"%.10e\t",Node_Coordinate[i][j]);
		}
		fprintf(fp,"\n");
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


	if (Total_mesh == 1)
	{
		printf("start GP calc\n");
		// Calculation_at_GP(E, nu);
		printf("end GP calc\n");
	}


	// 重ね合わせをスキップしてJ積分を行う
	if (SKIP_S_IGA != 1)
	{
		// 重ね合わせの結果
		division_ele_xi  = DIVISION_ELE_XI;
		division_ele_eta = DIVISION_ELE_ETA;

		if (division_ele_xi > MAX_DIVISION) {
			printf("Error!!\n");
			printf("Too many Divsion at xi!\n"
				"Maximum of division is %d (Now %d)\n"
				"\n", MAX_DIVISION, division_ele_xi);
			exit(1);
		}
		if (division_ele_eta > MAX_DIVISION) {
			printf("Error!!\n");
			printf("Too many Divsion at eta!\n"
				"Maximum of division is %d (Now %d)\n"
				"\n", MAX_DIVISION, division_ele_eta);
			exit(1);
		}

		// ローカルメッシュの情報取得
		GetLocData();

		int patch_n_loc = 0, patch_n_glo = 0;	// パッチ番号

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

		// S-IGAでローカル上のガウス点で重ね合わせた値をデータ整理する場合に使う
		if (Total_mesh >= 2)
		{
			printf("start GP calc\n");

			// patch_n_glo = 0;
			// Calculation_overlay_at_GP(E, nu,
			// 						  order_xi[patch_n_glo],order_eta[patch_n_glo],
			// 						  knot_n_xi[patch_n_glo], knot_n_eta[patch_n_glo],
			// 						  cntl_p_n_xi[patch_n_glo], cntl_p_n_eta[patch_n_glo],
			// 						  knot_vec_xi[patch_n_glo], knot_vec_eta[patch_n_glo],
			// 						  cntl_px[patch_n_glo], cntl_py[patch_n_glo],
			// 						  disp_cntl_px[patch_n_glo], disp_cntl_py[patch_n_glo],
			// 						  weight[patch_n_glo]);

			printf("end GP calc\n");
		}

		// 重ね合わせ
		for (i = 0; i < patch_n; i++) {

			fp = fopen("stress_vm_graph.txt", "a");
			fprintf(fp, "\npatch_n;%d\n\n",i);
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

			if (i >= n_patch_glo)	// ローカル上のパッチに対しては重合計算行う
			{
				patch_n_loc = i;
				printf("----------Start overlay calculation at patch %d in LOCAL patch----------\n\n", i);
				for (j = 0; j < n_patch_glo; j++)
				{
					patch_n_glo = j;
					Calculation_overlay(order_xi[patch_n_loc],order_eta[patch_n_loc],
										knot_n_xi[patch_n_loc], knot_n_eta[patch_n_loc],
										cntl_p_n_xi[patch_n_loc], cntl_p_n_eta[patch_n_loc],
										knot_vec_xi[patch_n_loc], knot_vec_eta[patch_n_loc],
										cntl_px[patch_n_loc], cntl_py[patch_n_loc],
										weight[patch_n_loc],
										order_xi[patch_n_glo],order_eta[patch_n_glo],
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
	int Location_Crack_Tip_Patch;
	double Location_Local_Coordinates[DIMENSION];
	double Virtual_Crack_Extension_Ct_Pt[MAX_N_NODE][DIMENSION];
	double DeltaA;
	J_Integral_Input_Data(Total_Control_Point_to_mesh[Total_mesh],&Location_Crack_Tip_Patch,Location_Local_Coordinates,Virtual_Crack_Extension_Ct_Pt, &DeltaA);

	Make_Displacement_grad(El_No);
	printf("Finish Make_Displacement_grad\n");
	Make_StrainEnergyDensity_2D();
	printf("Finish Make_StrainEnergyDensity\n");

	// globalでのDisp_gradを求める
	int e, N;
	Make_Displacement_grad_glo_check(real_Total_Element_to_mesh[Total_mesh]);
	
	Make_gauss_array(0);

	for(re = 0; re < real_Total_Element_to_mesh[Total_mesh]; re++){
        e = real_element[re];
		for(N = 0; N < GP_2D; N++){
			// Disp_gradをlocalとglobalで重ね合わせる
			for(i = 0; i <  DIMENSION * DIMENSION; i++){
				Disp_grad_overlay[e][N][i] = Disp_grad[e][N][i] + Disp_grad_glo[e][N][i];
				// printf("Disp_grad[%d][%d][%d] = %.10e\tDisp_grad_glo[%d][%d][%d] = %.10e\n", e, N, i, Disp_grad[e][N][i], e, N, i, Disp_grad_glo[e][N][i]);
				printf("Disp_grad_overlay[%d][%d][%d] = %.10e\n", e, N, i, Disp_grad_overlay[e][N][i]);
			}
			// Strainをlocalとglobalで重ね合わせる
			for (i = 0; i < D_MATRIX_SIZE; i++){
				Strain_overlay[e][N][i] = Strain[e][N][i] + Strain_glo[e][N][i];
				printf("Strain_overlay[%d][%d][%d] = %.10e\n", e, N, i, Strain_overlay[e][N][i]);
			}
		}
	}

	double unit_basis_local[DIMENSION] = {0.0};
	double r_tip = sqrt(Location_Local_Coordinates[0] * Location_Local_Coordinates[0] + Location_Local_Coordinates[1] * Location_Local_Coordinates[1]);

	// x'-y'（き裂先端）座標における単位基底ベクトル
	unit_basis_local[0] = Location_Local_Coordinates[0] / r_tip;
	unit_basis_local[1] = Location_Local_Coordinates[1] / r_tip;
	printf("unit_basis[0] : % 1.8e\n", unit_basis_local[0]);
	printf("unit_basis[1] : % 1.8e\n", unit_basis_local[1]);

	// 2 × 2 の行列として局所座標系を登録
	T[0][0] = unit_basis_local[0];   T[0][1] = unit_basis_local[1];   // T[0][2] = 0.0;
	T[1][0] = -unit_basis_local[1];  T[1][1] = unit_basis_local[0];   // T[1][2] = 0.0;
	// T[2][0] = 0.0;  T[2][1] = 0.0;   T[2][2] = 0.0;
	printf("T[0][0] = %1.8e\tT[0][1] = %1.8e\tT[1][0] = %1.8e\tT[1][1] = %1.8e\n", T[0][0], T[0][1], T[1][0], T[1][1]);

	double TT[DIMENSION][DIMENSION];
	for (i = 0; i < DIMENSION; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			TT[i][j] = T[j][i];
		}
	}

	for(re = 0; re < real_Total_Element_to_mesh[Total_mesh]; re++){
        e = real_element[re];
		for(N = 0; N < GP_2D; N++){
			double TEMP_Strain_overlay[DIMENSION][DIMENSION], TEMP_Strain_overlay_loc[DIMENSION][DIMENSION];
			TEMP_Strain_overlay[0][0] = Strain_overlay[e][N][0];
			TEMP_Strain_overlay[0][1] = Strain_overlay[e][N][2];
			TEMP_Strain_overlay[1][0] = Strain_overlay[e][N][2];
			TEMP_Strain_overlay[1][1] = Strain_overlay[e][N][1];
			Coordinate_Transform(TT, TEMP_Strain_overlay, TEMP_Strain_overlay_loc);
			Strain_overlay_loc[e][N][0] = TEMP_Strain_overlay_loc[0][0];
			Strain_overlay_loc[e][N][2] = TEMP_Strain_overlay_loc[0][1];
			Strain_overlay_loc[e][N][2] = TEMP_Strain_overlay_loc[1][0];
			Strain_overlay_loc[e][N][1] = TEMP_Strain_overlay_loc[1][1];
		}
	}

	Make_Stress_2D_glo(E, nu, Total_Element_to_mesh[Total_mesh], DM);

	Make_gauss_array(0);

	double D[D_MATRIX_SIZE][D_MATRIX_SIZE];
	Make_D_Matrix_2D(D, E, nu, DM);

	for(re = 0; re < real_Total_Element_to_mesh[Total_mesh]; re++){
		e = real_element[re];
		for(N = 0; N < GP_2D; N++){
			for (i = 0; i < D_MATRIX_SIZE; i++){
				Stress_overlay[e][N][i] = Stress[e][N][i] + Stress_glo[e][N][i];
				for (j = 0; j < D_MATRIX_SIZE; j++){
					Stress_overlay_loc[e][N][i] += D[i][j] * Strain_overlay_loc[e][N][j];
				}
			}
		}
	}

	Make_StrainEnergyDensity_2D_overlay();
	printf("Finish Make_StrainEnergyDensity\n");
	Make_Parameter_z_overlay(Total_Element_to_mesh[Total_mesh], E, nu, DM);
	printf("Finish Make_Parameter_z_overlay\n");

	printf("Start J_Integral_Computation\n");

	double Kfact1_ref, Kfact2_ref, ref_Sigma = 1.0, ref_a = 2.5, ref_W = 50.0, Pi = 4.0 * atan(1.0), rad = Pi / 4.0;
	printf("Pi = %.15e E = %.15e  nu = %.15e\n", Pi, E, nu);
	printf("o = %1.8e\ta = %1.8e\tW = %1.8e\n", ref_Sigma, ref_a, ref_W);
	Kfact1_ref = ref_Sigma * sin(rad) * sin(rad) * sqrt(Pi * ref_a);
	Kfact2_ref = ref_Sigma * sin(rad) * cos(rad) * sqrt(Pi * ref_a);

	J_Integral_Computation_Interaction(Total_Control_Point_to_mesh[Total_mesh], Location_Local_Coordinates, Virtual_Crack_Extension_Ct_Pt, DeltaA, E, nu, DM);

	double K_aux_mode1, K_aux_mode2;
	K_aux_mode1 = sqrt(J_integral_value_aux_mode1 * E / (1.0 - nu * nu));
	K_aux_mode2 = sqrt(J_integral_value_aux_mode2 * E / (1.0 - nu * nu));

	// for interaction integral（モデルの対称性によりかける倍率を変える。 1/4モデル → * 2.0、　フルモデル → * 1.0）
	printf("J_value_aux_mode1 = %1.10e\n", J_integral_value_aux_mode1);
	printf("J_value_aux_mode2 = %1.10e\n", J_integral_value_aux_mode2);

	printf("K_aux_mode1 = %.15e\n", K_aux_mode1);
	printf("K_aux_mode2 = %.15e\n", K_aux_mode2);

	printf("**** J-Integral Value 3 by 3 Gauss Integration ***\n");
	printf("K_mode1 = %.15e\n", K_mode1);
	printf("K_mode2 = %.15e\n", K_mode2);

	// 誤差評価
	printf("Kfact1_ref = %.15e\n", Kfact1_ref);
	printf("Kfact2_ref = %.15e\n", Kfact2_ref);
	printf("K_error(I_integral_mode1) = %.15e\n", (K_mode1 - Kfact1_ref) / Kfact1_ref * 100.0);
	printf("K_error(I_integral_mode2) = %.15e\n", (K_mode2 - Kfact2_ref) / Kfact2_ref * 100.0);
	printf("****                                           ***\n");

	t1 = clock();
	printf("Total calculation time:%.2f[s]\n", (double)(t1-start)/CLOCKS_PER_SEC);

	return 0;
}