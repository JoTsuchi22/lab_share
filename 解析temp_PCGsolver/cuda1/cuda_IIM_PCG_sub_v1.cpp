#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>


// header
#include "s_IGA_header.h"
#include "s_IGA_sub.h"


// ファイル読み込み1回目
void Get_Input_1(int tm, int *Total_Knot_to_mesh,
				 int *Total_Patch_on_mesh, int *Total_Patch_to_mesh,
				 int *Total_Control_Point_on_mesh, int *Total_Control_Point_to_mesh,
				 int *Total_Constraint_to_mesh, int *Total_Load_to_mesh, int *Total_DistributeForce_to_mesh,
				 char **argv)
{
	char s[256];
	int temp_i;
	double temp_d;

	int i, j;

	fp = fopen(argv[tm + 1], "r");

	// 材料定数
	fscanf(fp, "%lf %lf", &E, &nu);
	fgets(s, 256, fp);
	printf("E: %le, nu: %le\n", E, nu);

	// パッチ数
	fscanf(fp, "%d", &temp_i);
	fgets(s, 256, fp);
	int No_Patch = temp_i;
	int CP[temp_i * DIMENSION];
	printf("No_Patch: %d\n", temp_i);
	Total_Patch_on_mesh[tm] = temp_i;
	Total_Patch_to_mesh[tm + 1] = Total_Patch_to_mesh[tm] + temp_i;
	printf("Total_Patch_to_mesh[%d] = %d\n", tm, Total_Patch_to_mesh[tm]);

	// コントロールポイント数
	fscanf(fp, "%d", &temp_i);
	fgets(s, 256, fp);
	int Total_Control_Point = temp_i;
	printf("Total_Control_Point: %d\n", temp_i);
	Total_Control_Point_on_mesh[tm] = temp_i;
	Total_Control_Point_to_mesh[tm + 1] = Total_Control_Point_to_mesh[tm] + temp_i;
	printf("Total_Control_Point_to_mesh[%d] = %d\n", tm, Total_Control_Point_to_mesh[tm]);

	// 各方向の次数(スキップ)
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			fscanf(fp, "%d", &temp_i);
		}
	}
	fgets(s, 256, fp);

	// ノット数
	Total_Knot_to_mesh[tm + 1] = Total_Knot_to_mesh[tm];
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			fscanf(fp, "%d", &temp_i);
			Total_Knot_to_mesh[tm + 1] += temp_i;
		}
	}

	// 各パッチ各方向のコントロールポイント数(スキップ)
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			fscanf(fp, "%d", &temp_i);
			CP[i * DIMENSION + j] = temp_i;
		}
	}
	fgets(s, 256, fp);

	// パッチコネクティビティ(スキップ)
	for (i = 0; i < Total_Patch_on_mesh[tm]; i++)
	{
		for (j = 0; j < CP[i * DIMENSION + 0] * CP[i * DIMENSION + 1]; j++)
		{
		fscanf(fp, "%d", &temp_i);
		}
	}
	fgets(s, 256, fp);

	// 境界条件
	int Total_Constraint, Total_Load, Total_DistributeForce;
	fscanf(fp, "%d %d %d", &Total_Constraint, &Total_Load, &Total_DistributeForce);
	Total_Constraint_to_mesh[tm + 1] = Total_Constraint_to_mesh[tm] + Total_Constraint;
	Total_Load_to_mesh[tm + 1] = Total_Load_to_mesh[tm] + Total_Load;
	Total_DistributeForce_to_mesh[tm + 1] = Total_DistributeForce_to_mesh[tm] + Total_DistributeForce;

	printf("Total_Constraint = %d\n", Total_Constraint);
	printf("Total_Load = %d\n", Total_Load);
	printf("Total_DistributedForce = %d\n", Total_DistributeForce);

	fclose(fp);
}


// ファイル読み込み2回目
void Get_Input_2(int tm, int *Total_Knot_to_mesh,
				 int *Total_Patch_to_mesh, int *Total_Control_Point_to_mesh,
				 int *Total_Element_on_mesh, int *Total_Element_to_mesh,
				 int *Total_Constraint_to_mesh, int *Total_Load_to_mesh, int *Total_DistributeForce_to_mesh,
				 int *Order, int *No_knot, int *No_Control_point, int *No_Control_point_in_patch,
				 int *Patch_Control_point, double *Position_Knots, int *No_Control_point_ON_ELEMENT,
				 double *Node_Coordinate, double *Control_Coord_x, double *Control_Coord_y, double *Control_Weight,
				 int *Constraint_Node_Dir, double *Value_of_Constraint,
				 int *Load_Node_Dir, double *Value_of_Load,
				 int *type_load_array, int *iPatch_array, int *iCoord_array,
				 double *val_Coord_array, double *Range_Coord_array, double *Coeff_Dist_Load_array,
				 int *Total_Control_Point_to_patch, int *Total_Knot_to_patch_dim,
				 char **argv)
{
	char s[256];
	int temp_i;
	double temp_d;

	int i, j, k;

	fp = fopen(argv[tm + 1], "r");

	// 材料定数(スキップ)
	fscanf(fp, "%lf %lf", &temp_d, &temp_d);
	fgets(s, 256, fp);

	// パッチ数(スキップ)
	fscanf(fp, "%d", &temp_i);
	fgets(s, 256, fp);
	int No_Patch = temp_i;

	// コントロールポイント数(スキップ)
	fscanf(fp, "%d", &temp_i);
	fgets(s, 256, fp);
	int Total_Control_Point = temp_i;

	// 各方向の次数
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			fscanf(fp, "%d", &temp_i);
			Order[(i + Total_Patch_to_mesh[tm]) * DIMENSION + j] = temp_i;
			printf("Order[%d] = %d\n", (i + Total_Patch_to_mesh[tm]) * DIMENSION + j, Order[(i + Total_Patch_to_mesh[tm]) * DIMENSION + j]);
		}
	}
	fgets(s, 256, fp);

	// ノット数
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			fscanf(fp, "%d", &temp_i);
			No_knot[(i + Total_Patch_to_mesh[tm]) * DIMENSION + j] = temp_i;
			Total_Knot_to_patch_dim[(i + Total_Patch_to_mesh[tm]) * DIMENSION + j + 1] = Total_Knot_to_patch_dim[(i + Total_Patch_to_mesh[tm]) * DIMENSION + j] + temp_i;
			printf("No_knot[%d] = %d\n", (i + Total_Patch_to_mesh[tm]) * DIMENSION + j, No_knot[(i + Total_Patch_to_mesh[tm]) * DIMENSION + j]);
		}
	}
	fgets(s, 256, fp);

	// 各パッチ各方向のコントロールポイント数
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			fscanf(fp, "%d", &temp_i);
			No_Control_point[(i + Total_Patch_to_mesh[tm]) * DIMENSION + j] = temp_i;
			printf("No_Control_point[%d] = %d\n", (i + Total_Patch_to_mesh[tm]) * DIMENSION + j, No_Control_point[(i + Total_Patch_to_mesh[tm]) * DIMENSION + j]);
		}
	}
	fgets(s, 256, fp);

	for (i = 0; i < No_Patch; i++)
	{
		No_Control_point_in_patch[i + Total_Patch_to_mesh[tm]] = 1;
	}

	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			No_Control_point_in_patch[i + Total_Patch_to_mesh[tm]] *= No_Control_point[(i + Total_Patch_to_mesh[tm]) * DIMENSION + j];
		}
	}

	for (i = 0; i < No_Patch; i++)
	{
		Total_Control_Point_to_patch[i + Total_Patch_to_mesh[tm] + 1] = Total_Control_Point_to_patch[i + Total_Patch_to_mesh[tm]] + No_Control_point_in_patch[i + Total_Patch_to_mesh[tm]];
	}

	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			if (No_knot[(i + Total_Patch_to_mesh[tm]) * DIMENSION + j] != No_Control_point[(i + Total_Patch_to_mesh[tm]) * DIMENSION + j] + Order[(i + Total_Patch_to_mesh[tm]) * DIMENSION + j] + 1)
			{
				printf("wrong relationship between the number of knot vector and the number of control_point \n");
				printf("in mesh_No.%d in patch_No.%d direction:%d\n", tm, i, j);
			}
		}
	}

	for (i = 0; i < No_Patch; i++)
	{
		printf("No_Control_point_in_patch[%d] = %d\t", i + Total_Patch_to_mesh[tm], No_Control_point_in_patch[i + Total_Patch_to_mesh[tm]]);
	}
	printf("\n");

	// パッチコネクティビティ
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < No_Control_point_in_patch[i + Total_Patch_to_mesh[tm]]; j++)
		{
			fscanf(fp, "%d", &temp_i);
			Patch_Control_point[Total_Control_Point_to_patch[i + Total_Patch_to_mesh[tm]] + j] = temp_i;
			if (tm > 0)
			{
				Patch_Control_point[Total_Control_Point_to_patch[i + Total_Patch_to_mesh[tm]] + j] += Total_Control_Point_to_mesh[tm];
			}
		}
	}
	fgets(s, 256, fp);

	// 境界条件(スキップ)
	int Total_Constraint, Total_Load, Total_DistributeForce;
	fscanf(fp, "%d %d %d", &Total_Constraint, &Total_Load, &Total_DistributeForce);
	fgets(s, 256, fp);

	// ノットベクトルの読み込み
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			for (k = 0; k < No_knot[(i + Total_Patch_to_mesh[tm]) * DIMENSION + j]; k++)
			{
				fscanf(fp, "%lf", &temp_d);
				Position_Knots[Total_Knot_to_patch_dim[(i + Total_Patch_to_mesh[tm]) * DIMENSION + j] + k] = temp_d;
				printf("%le\t", Position_Knots[Total_Knot_to_patch_dim[(i + Total_Patch_to_mesh[tm]) * DIMENSION + j] + k]);
			}
			printf("\n");
		}
	}
	fgets(s, 256, fp);

	int Total_Element = 0;
	for (i = 0; i < No_Patch; i++)
	{
		if (DIMENSION == 2)
		{
			Total_Element += (No_Control_point[(i + Total_Patch_to_mesh[tm]) * DIMENSION + 0] - Order[(i + Total_Patch_to_mesh[tm]) * DIMENSION + 0])
						   * (No_Control_point[(i + Total_Patch_to_mesh[tm]) * DIMENSION + 1] - Order[(i + Total_Patch_to_mesh[tm]) * DIMENSION + 1]);
			No_Control_point_ON_ELEMENT[i + Total_Patch_to_mesh[tm]] = (Order[(i + Total_Patch_to_mesh[tm]) * DIMENSION + 0] + 1) * (Order[(i + Total_Patch_to_mesh[tm]) * DIMENSION + 1] + 1);
		}
		else if (DIMENSION == 3)
		{
			Total_Element += (No_Control_point[(i + Total_Patch_to_mesh[tm]) * DIMENSION + 0] - Order[(i + Total_Patch_to_mesh[tm]) * DIMENSION + 0])
						   * (No_Control_point[(i + Total_Patch_to_mesh[tm]) * DIMENSION + 1] - Order[(i + Total_Patch_to_mesh[tm]) * DIMENSION + 1])
						   * (No_Control_point[(i + Total_Patch_to_mesh[tm]) * DIMENSION + 2] - Order[(i + Total_Patch_to_mesh[tm]) * DIMENSION + 2]);
			No_Control_point_ON_ELEMENT[i + Total_Patch_to_mesh[tm]] = (Order[(i + Total_Patch_to_mesh[tm]) * DIMENSION + 0] + 1) * (Order[(i + Total_Patch_to_mesh[tm]) * DIMENSION + 1] + 1) * (Order[(i + Total_Patch_to_mesh[tm]) * DIMENSION + 2] + 1);
		}
	}
	printf("Total_Element = %d\n", Total_Element);
	Total_Element_on_mesh[tm] = Total_Element;
	Total_Element_to_mesh[tm + 1] = Total_Element_to_mesh[tm] + Total_Element;
	printf("Total_Element_on_mesh[%d] = %d\n", tm, Total_Element_on_mesh[tm]);

	for (i = 0; i < No_Patch; i++)
	{
		printf("No_Control_point_ON_ELEMENT[%d] = %d\n", i + Total_Patch_to_mesh[tm], No_Control_point_ON_ELEMENT[i + Total_Patch_to_mesh[tm]]);
	}

	// 節点座標
	for (i = 0; i < Total_Control_Point; i++)
	{
		fscanf(fp, "%d", &temp_i);
		for (j = 0; j < DIMENSION + 1; j++)
		{
			fscanf(fp, "%lf", &temp_d);
			Node_Coordinate[(temp_i + Total_Control_Point_to_mesh[tm]) * (DIMENSION + 1) + j] = temp_d;
		}
	}
	for (i = 0; i < Total_Control_Point; i++)
	{
		for (j = 0; j < DIMENSION + 1; j++)
		{
			// コントロールポイント座標・重みの新たな配列(for s-IGA/NewtonLaphson) DIMENSION == 2 の場合のみ記述
			if (j == 0)
			{
				Control_Coord_x[i + Total_Control_Point_to_mesh[tm]] = Node_Coordinate[(i + Total_Control_Point_to_mesh[tm]) * (DIMENSION + 1) + j];
			}
			else if (j == 1)
			{
				Control_Coord_y[i + Total_Control_Point_to_mesh[tm]] = Node_Coordinate[(i + Total_Control_Point_to_mesh[tm]) * (DIMENSION + 1) + j];
			}
			else if (j == DIMENSION)
			{
				Control_Weight[i + Total_Control_Point_to_mesh[tm]] = Node_Coordinate[(i + Total_Control_Point_to_mesh[tm]) * (DIMENSION + 1) + j];
			}
		}
	}
	fgets(s, 256, fp);

	// 拘束
	for (i = 0; i < Total_Constraint; i++)
	{
		fscanf(fp, "%d %d %lf",
			   &Constraint_Node_Dir[(i + Total_Constraint_to_mesh[tm]) * 2 + 0],
			   &Constraint_Node_Dir[(i + Total_Constraint_to_mesh[tm]) * 2 + 1],
			   &Value_of_Constraint[i + Total_Constraint_to_mesh[tm]]);
		Constraint_Node_Dir[(i + Total_Constraint_to_mesh[tm]) * 2 + 0] = Constraint_Node_Dir[(i + Total_Constraint_to_mesh[tm]) * 2 + 0] + Total_Control_Point_to_mesh[tm];

		printf("Constraint_Node_Dir[%d] = %d Constraint_Node_Dir[%d] = %d Value_of_Constraint[%d] = %le \n",
			   (i + Total_Constraint_to_mesh[tm]) * 2 + 0, Constraint_Node_Dir[(i + Total_Constraint_to_mesh[tm]) * 2 + 0],
			   (i + Total_Constraint_to_mesh[tm]) * 2 + 1, Constraint_Node_Dir[(i + Total_Constraint_to_mesh[tm]) * 2 + 1],
			   i + Total_Constraint_to_mesh[tm], Value_of_Constraint[i + Total_Constraint_to_mesh[tm]]);
	}
	fgets(s, 256, fp);

	// 荷重
	for (i = 0; i < Total_Load; i++)
	{
		fscanf(fp, "%d %d %lf",
			   &Load_Node_Dir[(i + Total_Load_to_mesh[tm]) * 2 + 0],
			   &Load_Node_Dir[(i + Total_Load_to_mesh[tm]) * 2 + 1],
			   &Value_of_Load[i + Total_Load_to_mesh[tm]]);
		Load_Node_Dir[(i + Total_Load_to_mesh[tm]) * 2 + 0] = Load_Node_Dir[(i + Total_Load_to_mesh[tm]) * 2 + 0] + Total_Control_Point_to_mesh[tm];

		printf("Load_Node_Dir[%d]= %d Load_Node_Dir[%d]= %d Value_of_Load[%d] = %le\n",
			   (i + Total_Load_to_mesh[tm]) * 2 + 0, Load_Node_Dir[(i + Total_Load_to_mesh[tm]) * 2 + 0],
			   (i + Total_Load_to_mesh[tm]) * 2 + 1, Load_Node_Dir[(i + Total_Load_to_mesh[tm]) * 2 + 1],
			   i + Total_Load_to_mesh[tm], Value_of_Load[i + Total_Load_to_mesh[tm]]);
	}
	fgets(s, 256, fp);

	int type_load, iPatch, iCoord;
	double val_Coord, Range_Coord[2], Coeff_Dist_Load[3];

	for (i = 0; i < Total_DistributeForce; i++)
	{
		fscanf(fp, "%d %d %d %lf %lf %lf %lf %lf %lf", &type_load, &iPatch, &iCoord, &val_Coord, &Range_Coord[0], &Range_Coord[1], &Coeff_Dist_Load[0], &Coeff_Dist_Load[1], &Coeff_Dist_Load[2]);
		printf("Distibuted load nober: %d\n", i);
		printf("type_load: %d  iPatch: %d iCoord: %d  val_Coord: %.15e  Range_Coord: %.15e  %.15e\n Coef_Dist_Load: %.15e %.15e %.15e\n",
			   type_load, iPatch, iCoord, val_Coord, Range_Coord[0], Range_Coord[1], Coeff_Dist_Load[0], Coeff_Dist_Load[1], Coeff_Dist_Load[2]);

		// for s-IGA
		type_load_array[i + Total_DistributeForce_to_mesh[tm]] = type_load;
		iPatch_array[i + Total_DistributeForce_to_mesh[tm]] = iPatch;
		iCoord_array[i + Total_DistributeForce_to_mesh[tm]] = iCoord;
		val_Coord_array[i + Total_DistributeForce_to_mesh[tm]] = val_Coord;
		Range_Coord_array[(i + Total_DistributeForce_to_mesh[tm]) * 2 + 0] = Range_Coord[0];
		Range_Coord_array[(i + Total_DistributeForce_to_mesh[tm]) * 2 + 1] = Range_Coord[1];
		Coeff_Dist_Load_array[(i + Total_DistributeForce_to_mesh[tm]) * 3 + 0] = Coeff_Dist_Load[0];
		Coeff_Dist_Load_array[(i + Total_DistributeForce_to_mesh[tm]) * 3 + 1] = Coeff_Dist_Load[1];
		Coeff_Dist_Load_array[(i + Total_DistributeForce_to_mesh[tm]) * 3 + 2] = Coeff_Dist_Load[2];
	}
	fclose(fp);
}


// INC 等の作成
void Make_INC(int tm, int *Total_Patch_on_mesh, int *Total_Patch_to_mesh,
			  int *Total_Element_on_mesh, int *Total_Element_to_mesh,
			  int *Total_Control_Point_on_mesh, int *Total_Control_Point_to_mesh,
			  int *Total_DistributeForce_to_mesh,
			  int *INC, int *Patch_Control_point, int *Total_Control_Point_to_patch,
			  int *No_Control_point, int *Controlpoint_of_Element,
			  int *Element_patch, int *Element_mesh,
			  double *difference, int *Total_Knot_to_patch_dim,
			  int *Total_element_all_ID, int *ENC,
			  int *real_element_line, int *line_No_real_element,
			  int *line_No_Total_element, int *real_element,
			  int *real_El_No_on_mesh, int *real_Total_Element_on_mesh,
			  int *real_Total_Element_to_mesh, double *Equivalent_Nodal_Force,
			  int *type_load_array, int *iPatch_array, int *iCoord_array,
			  double *val_Coord_array, double *Range_Coord_array, double *Coeff_Dist_Load_array,
			  int *Order, int *No_knot)
{
	// INC の計算(節点番号をξ, ηの番号で表す為の配列)
	for (tm = 0; tm < Total_mesh; tm++)
	{
		int b, B, e, h, i, j, k, l, n, p, q, x, y, ii, jj, kk, iii, kkk, iiloc, jjloc, kkloc, r = 0;
		int type_load, iPatch, iCoord;
		double val_Coord, Range_Coord[2], Coeff_Dist_Load[3];
		int No_Patch = Total_Patch_on_mesh[tm];
		int Total_Patch_to_Now = Total_Patch_to_mesh[tm];
		int Total_Element = Total_Element_on_mesh[tm];
		int Total_Element_to_Now = Total_Element_to_mesh[tm];
		int Total_Control_Point = Total_Control_Point_on_mesh[tm];
		int Total_DistributeForce = Total_DistributeForce_to_mesh[tm + 1] - Total_DistributeForce_to_mesh[tm];
		if (DIMENSION == 2) // for s-IGA
		{
			e = 0;
			for (l = 0; l < No_Patch; l++)
			{
				i = 0;
				for (jj = 0; jj < No_Control_point[(l + Total_Patch_to_Now) * DIMENSION + 1]; jj++)
				{
					for (ii = 0; ii < No_Control_point[(l + Total_Patch_to_Now) * DIMENSION + 0]; ii++)
					{
						INC[Patch_Control_point[Total_Control_Point_to_patch[l + Total_Patch_to_Now] + i] * DIMENSION + 0] = ii;
						INC[Patch_Control_point[Total_Control_Point_to_patch[l + Total_Patch_to_Now] + i] * DIMENSION + 1] = jj;

						if (ii >= Order[(l + Total_Patch_to_Now) * DIMENSION + 0] && jj >= Order[(l + Total_Patch_to_Now) * DIMENSION + 1])
						{

							for (jjloc = 0; jjloc <= Order[(l + Total_Patch_to_Now) * DIMENSION + 1]; jjloc++)
							{
								for (iiloc = 0; iiloc <= Order[(l + Total_Patch_to_Now) * DIMENSION + 0]; iiloc++)
								{
									B = Patch_Control_point[Total_Control_Point_to_patch[l + Total_Patch_to_Now] + i - jjloc * No_Control_point[(l + Total_Patch_to_Now) * DIMENSION + 0] - iiloc];
									b = jjloc * (Order[(l + Total_Patch_to_Now) * DIMENSION + 0] + 1) + iiloc;
									Controlpoint_of_Element[(e + Total_Element_to_Now) * MAX_NO_CCpoint_ON_ELEMENT + b] = B;
								}
							}
							Element_patch[e + Total_Element_to_Now] = l + Total_Patch_to_Now;
							Element_mesh[e + Total_Element_to_Now] = tm;
							e++;
						}
						i++;
					}
				}
			}
		}

		// 2次元のデバッグのための一時的なコメントアウト，3次元のプログラム作成時にはコメントアウトをはずす	
		// if (DIMENSION == 3)
		// {
		// 	e = 0;
		// 	for (l = 0; l < No_Patch; l++)
		// 	{
		// 		i = 0;
		// 		for (kk = 0; kk < No_Control_point[l][2]; kk++)
		// 		{
		// 			for (jj = 0; jj < No_Control_point[l][1]; jj++)
		// 			{
		// 				for (ii = 0; ii < No_Control_point[l][0]; ii++)
		// 				{
		// 					INC[l][Patch_Control_point[l][i]][0] = ii;
		// 					INC[l][Patch_Control_point[l][i]][1] = jj;
		// 					INC[l][Patch_Control_point[l][i]][2] = kk;
		// 					if (ii >= Order[l][0] && jj >= Order[l][1] && kk >= Order[l][2])
		// 					{
		// 						for (kkloc = 0; kkloc < Order[l][2]; kkloc++)
		// 						{
		// 							for (jjloc = 0; jjloc <= Order[l][1]; jjloc++)
		// 							{
		// 								for (iiloc = 0; iiloc <= Order[l][0]; iiloc++)
		// 								{
		// 									B = Patch_Control_point[l][i - jjloc * No_Control_point[l][0] - iiloc];
		// 									b = jjloc * (Order[l][0] + 1) + iiloc;
		// 									Controlpoint_of_Element[e * MAX_NO_CCpoint_ON_ELEMENT + b] = B;
		// 								}
		// 							}
		// 						}
		// 						Element_patch[e] = l + Total_Patch_to_Now;
		// 						e++;
		// 					}
		// 					i++;
		// 				}
		// 			}
		// 		}
		// 	}
		// }

		// for s-IGA line_No_real_elementの初期化
		for (l = 0; l < No_Patch; l++)
		{
			for (j = 0; j < DIMENSION; j++)
			{
				line_No_real_element[(l + Total_Patch_to_Now) * DIMENSION + j] = 0;
			}
		}

		for (l = 0; l < No_Patch; l++)
		{
			for (j = 0; j < DIMENSION; j++)
			{
				line_No_Total_element[(l + Total_Patch_to_Now) * DIMENSION + j] = No_knot[(l + Total_Patch_to_Now) * DIMENSION + j] - 2 * Order[(l + Total_Patch_to_Now) * DIMENSION + j] - 1;

				for (kkk = Order[(l + Total_Patch_to_Now) * DIMENSION + j]; kkk < No_knot[(l + Total_Patch_to_Now) * DIMENSION + j] - Order[(l + Total_Patch_to_Now) * DIMENSION + j] - 1; kkk++)
				{
					difference[Total_Knot_to_patch_dim[(l + Total_Patch_to_Now) * DIMENSION + j] + kkk - Order[(l + Total_Patch_to_Now) * DIMENSION + j]];
					if (difference[Total_Knot_to_patch_dim[(l + Total_Patch_to_Now) * DIMENSION + j] + kkk - Order[(l + Total_Patch_to_Now) * DIMENSION + j]] != 0)
					{
						line_No_real_element[(l + Total_Patch_to_Now) * DIMENSION + j]++;
					}
				}
			}
		}

		// 要素に行番号, 列番号をつける
		if (DIMENSION == 2)
		{
			for (h = 0; h < Total_Element; h++)
			{
				Total_element_all_ID[h] = 0;
			}

			i = 0;
			for (l = 0; l < No_Patch; l++)
			{
				for (y = 0; y < line_No_Total_element[(l + Total_Patch_to_Now) * DIMENSION + 1]; y++)
				{
					for (x = 0; x < line_No_Total_element[(l + Total_Patch_to_Now) * DIMENSION + 0]; x++)
					{
						ENC[(i + Total_Element_to_mesh[tm]) * DIMENSION + 0] = x;
						ENC[(i + Total_Element_to_mesh[tm]) * DIMENSION + 1] = y;
						i++;
					}
				}
			}
		}

		// 必要な要素の行と列の番号を求める
		for (j = 0; j < DIMENSION; j++)
		{
			for (l = 0; l < No_Patch; l++)
			{
				e = 0;

				for (k = 0; k < line_No_Total_element[(l + Total_Patch_to_Now) * DIMENSION + j]; k++)
				{
					if (difference[Total_Knot_to_patch_dim[(l + Total_Patch_to_Now) * DIMENSION + j] + k] != 0)
					{
						real_element_line[(l + Total_Patch_to_Now) * (Total_Element_to_mesh[Total_mesh] * DIMENSION) + e * DIMENSION + j] = k;
						e++;
					}
				}
			}
		}

		// 必要な要素列上の要素のIDを1にする
		if (DIMENSION == 2)
		{
			for (n = 0; n < Total_Element; n++)
			{
				for (p = 0; p < line_No_real_element[(Element_patch[n + Total_Element_to_Now]) * DIMENSION + 0]; p++)
				{
					if (ENC[(n + Total_Element_to_mesh[tm]) * DIMENSION + 0] == real_element_line[(Element_patch[n + Total_Element_to_Now]) * (Total_Element_to_mesh[Total_mesh] * DIMENSION) + p * DIMENSION + 0])
					{
						for (q = 0; q < line_No_real_element[(Element_patch[n + Total_Element_to_Now]) * DIMENSION + 1]; q++)
						{
							if (ENC[(n + Total_Element_to_mesh[tm]) * DIMENSION + 1] == real_element_line[(Element_patch[n + Total_Element_to_Now]) * (Total_Element_to_mesh[Total_mesh] * DIMENSION) + q * DIMENSION + 1])
							{
								Total_element_all_ID[n]++;
							}
						}
					}
				}

				// IDが1の要素に番号を振る
				if (Total_element_all_ID[n] == 1)
				{
					real_element[r + real_Total_Element_to_mesh[tm]] = n + Total_Element_to_Now;
					real_El_No_on_mesh[tm * Total_Element_to_mesh[Total_mesh] + r] = n + Total_Element_to_Now;
					r++;
				}
			}

			// for s-IGA real_Total_Elementの初期化
			int real_Total_Element = 0;

			for (l = 0; l < No_Patch; l++)
			{
				real_Total_Element += line_No_real_element[(l + Total_Patch_to_Now) * DIMENSION + 0] * line_No_real_element[(l + Total_Patch_to_Now) * DIMENSION + 1];
			}
			real_Total_Element_on_mesh[tm] = real_Total_Element;
			real_Total_Element_to_mesh[tm + 1] = real_Total_Element_to_mesh[tm] + real_Total_Element;
		}

		for (i = 0; i < Total_DistributeForce; i++)
		{
			type_load = type_load_array[i + Total_DistributeForce_to_mesh[tm]];
			iPatch = iPatch_array[i + Total_DistributeForce_to_mesh[tm]];
			iCoord = iCoord_array[i + Total_DistributeForce_to_mesh[tm]];
			val_Coord = val_Coord_array[i + Total_DistributeForce_to_mesh[tm]];
			Range_Coord[0] = Range_Coord_array[(i + Total_DistributeForce_to_mesh[tm]) * 2 + 0];
			Range_Coord[1] = Range_Coord_array[(i + Total_DistributeForce_to_mesh[tm]) * 2 + 1];
			Coeff_Dist_Load[0] = Coeff_Dist_Load_array[(i + Total_DistributeForce_to_mesh[tm]) * 3 + 0];
			Coeff_Dist_Load[1] = Coeff_Dist_Load_array[(i + Total_DistributeForce_to_mesh[tm]) * 3 + 1];
			Coeff_Dist_Load[2] = Coeff_Dist_Load_array[(i + Total_DistributeForce_to_mesh[tm]) * 3 + 2];

			printf("type_load:%d\tiPatch:%d\tiCoord:%d\tval_coord:%lf\t", type_load, iPatch, iCoord, val_Coord);
			printf("Range0:%lf\tRange1:%lf\t", Range_Coord[0], Range_Coord[1]);
			printf("Coeff0:%lf\n", Coeff_Dist_Load[0]);
			Setting_Dist_Load_2D(tm, iPatch, Total_Element_to_mesh[tm + 1], iCoord, val_Coord, Range_Coord, type_load, Coeff_Dist_Load);
		}
	}
}


// Distributed Load
void Setting_Dist_Load_2D(int mesh_n, int iPatch, int Total_Element, int iCoord, double val_Coord,
						  double Range_Coord[2], int type_load, double Coeff_Dist_Load[3], int *Total_Knot_to_mesh,
						  int *Controlpoint_of_Element, int *Order, int *No_knot, int *Total_element_all_ID,
						  int *Total_Knot_to_patch_dim, double *Position_Knots, double *Equivalent_Nodal_Force,
						  int *Total_Element_on_mesh, int *Total_Element_to_mesh, int *Element_patch, int *ENC,
						  double *Node_Coordinate, int *INC, int *No_Control_point_ON_ELEMENT, int *Total_Control_Point_to_mesh,
						  int *No_Control_point)
{
	int iii, jjj;
	int N_Seg_Load_Element_iDir = 0, jCoord;
	int iRange_ele[2];
	int iPos[2] = {-10000, -10000}, jPos[2] = {-10000, -10000};
	int No_Element_For_Dist_Load;
	int iX, iY;
	int iControlpoint[MAX_NO_CCpoint_ON_ELEMENT], ic, ig, NNG = 3;
	double val_jCoord_Local = 0.0;
	double GaussPt[3], Weight[3];
	double Gg = pow(3.0 / 5.0, 0.5);

	double Position_Data_param[DIMENSION];
	int *No_Element_for_Integration = (int *)malloc(sizeof(int) * Total_Knot_to_mesh[Total_mesh]); // No_Element_for_Integration[MAX_N_KNOT]
 
	GaussPt[0] = -Gg;
	GaussPt[1] = 0.0;
	GaussPt[2] = Gg;
	Weight[0] = 5.0 / 9.0;
	Weight[1] = 8.0 / 9.0;
	Weight[2] = 5.0 / 9.0;

	// iCoord=0: Load on Eta=Constant
	// iCoord=1: Load on Xi=Constant
	if (iCoord == 0)
		jCoord = 1;
	if (iCoord == 1)
		jCoord = 0;

	// val_Coord: Value of Eta or Xi of the line or surface to give the distributed load
	// Setting elements needed to computed the distributed load

	for (iii = Order[iPatch * DIMENSION + iCoord]; iii < No_knot[iPatch * DIMENSION + iCoord] - Order[iPatch * DIMENSION + iCoord] - 1; iii++)
	{
		double epsi = 0.00000000001;

		if (Position_Knots[Total_Knot_to_patch_dim[iPatch * DIMENSION + iCoord] + iii] - epsi <= Range_Coord[0])
			iPos[0] = iii;
		if (Position_Knots[Total_Knot_to_patch_dim[iPatch * DIMENSION + iCoord] + iii + 1] - epsi <= Range_Coord[1])
			iPos[1] = iii + 1;
	}
	iRange_ele[0] = iPos[0] - Order[iPatch * DIMENSION + iCoord];
	iRange_ele[1] = iPos[1] - Order[iPatch * DIMENSION + iCoord] - 1;
	printf("iPos[0] = %d iPos[1] = %d\n", iPos[0], iPos[1]);
	printf("iRange_ele[0] = %d  iRange_ele[1] = %d\n", iRange_ele[0], iRange_ele[1]);

	if (iPos[0] < 0 || iPos[1] < 0)
	{
		printf("Error (Stop) iPos[0] = %d iPos[1] = %d\n", iPos[0], iPos[1]);
		exit(0);
	}

	for (jjj = Order[iPatch * DIMENSION + jCoord]; jjj < No_knot[iPatch * DIMENSION + jCoord] - Order[iPatch * DIMENSION + jCoord] - 1; jjj++)
	{
		double epsi = 0.00000000001;

		if (Position_Knots[Total_Knot_to_patch_dim[iPatch * DIMENSION + jCoord] + jjj] - epsi <= val_Coord
			&& Position_Knots[Total_Knot_to_patch_dim[iPatch * DIMENSION + jCoord] + jjj + 1] + epsi > val_Coord)
		{
			jPos[0] = jjj;
			jPos[1] = jjj + 1;
			val_jCoord_Local = -1.0 + 2.0 * (val_Coord - Position_Knots[Total_Knot_to_patch_dim[iPatch * DIMENSION + jCoord] + jjj])
							 / (Position_Knots[Total_Knot_to_patch_dim[iPatch * DIMENSION + jCoord] + jjj + 1] - Position_Knots[Total_Knot_to_patch_dim[iPatch * DIMENSION + jCoord] + jjj]);
		}
	}

	if (jPos[0] < 0 || jPos[1] < 0)
	{
		printf("Error (Stop) jPos[0] = %d jPos[1] = %d\n", jPos[0], jPos[1]);
		exit(0);
	}

	for (iii = iPos[0]; iii < iPos[1]; iii++)
	{
		N_Seg_Load_Element_iDir++;
	}

	iii = 0;
	if (iCoord == 1)
	{
		iX = jPos[0] - Order[iPatch * DIMENSION + 0];
		for (iY = iPos[0] - Order[iPatch * DIMENSION + 1]; iY < iPos[1] - Order[iPatch * DIMENSION + 1]; iY++)
		{
			No_Element_for_Integration[iii] = SearchForElement(mesh_n, iPatch, iX, iY, Total_Element_on_mesh, Total_Element_to_mesh, Element_patch, ENC);
			printf("Check No_Element_for_Integration[%d] = %d\n", iii, No_Element_for_Integration[iii]);
			iii++;
		}
	}

	if (iCoord == 0)
	{
		iY = jPos[0] - Order[iPatch * DIMENSION + 1];
		for (iX = iPos[0] - Order[iPatch * DIMENSION + 0]; iX < iPos[1] - Order[iPatch * DIMENSION + 0]; iX++)
		{
			No_Element_for_Integration[iii] = SearchForElement(mesh_n, iPatch, iX, iY, Total_Element_on_mesh, Total_Element_to_mesh, Element_patch, ENC);
			printf("Check No_Element_for_Integration[%d] = %d\n", iii, No_Element_for_Integration[iii]);
			iii++;
		}
	}
	No_Element_For_Dist_Load = iii;

	// Book keeping finished

	for (iii = 0; iii < No_Element_For_Dist_Load; iii++)
	{
		if (Total_element_all_ID[No_Element_for_Integration[iii]] == 1)
		{
			iX = ENC[No_Element_for_Integration[iii] * DIMENSION + 0];
			iY = ENC[No_Element_for_Integration[iii] * DIMENSION + 1];

			for (ic = 0; ic < (Order[iPatch * DIMENSION + 0] + 1) * (Order[iPatch * DIMENSION + 1] + 1); ic++)
				iControlpoint[ic] = Controlpoint_of_Element[No_Element_for_Integration[iii] * MAX_NO_CCpoint_ON_ELEMENT + ic];

			for (ig = 0; ig < NNG; ig++)
			{
				double Local_Coord[2], sfc, dxyzdge[3], detJ, XiEtaCoordParen, valDistLoad;
				int icc;
				Local_Coord[jCoord] = val_jCoord_Local;
				Local_Coord[iCoord] = GaussPt[ig];
				printf("ig = %d   Local_Coord[jCoord] = %f Local_Coord[iCoord] = %f\n", ig, Local_Coord[jCoord], Local_Coord[iCoord]);

				ShapeFunc_from_paren(Position_Data_param, Local_Coord, iCoord, No_Element_for_Integration[iii], INC, Position_Knots, Total_Knot_to_patch_dim, Controlpoint_of_Element, Element_patch);
				XiEtaCoordParen = Position_Data_param[iCoord];
				printf("Check  Coeff_Dist_Load[0] = %f Coeff_Dist_Load[1] = %f  Coeff_Dist_Load[2] = %f  Position_Data_param[iCoord] = %f\n", Coeff_Dist_Load[0], Coeff_Dist_Load[1], Coeff_Dist_Load[2], Position_Data_param[iCoord]);
				valDistLoad = Coeff_Dist_Load[0] + Coeff_Dist_Load[1] * XiEtaCoordParen + Coeff_Dist_Load[2] * XiEtaCoordParen * XiEtaCoordParen;

				dxyzdge[0] = 0.0;
				dxyzdge[1] = 0.0;
				dxyzdge[2] = 0.0;
				for (icc = 0; icc < (Order[iPatch * DIMENSION + 0] + 1) * (Order[iPatch * DIMENSION + 1] + 1); icc++)
				{
					dxyzdge[0] += dShape_func(icc, iCoord, Local_Coord, No_Element_for_Integration[iii], Controlpoint_of_Element, No_Control_point_ON_ELEMENT, Total_Control_Point_to_mesh, Element_patch, Node_Coordinate, INC, Order, Position_Knots, Total_Knot_to_patch_dim, No_knot, No_Control_point)
								* Node_Coordinate[iControlpoint[icc] * (DIMENSION + 1) + 0];
					dxyzdge[1] += dShape_func(icc, iCoord, Local_Coord, No_Element_for_Integration[iii], Controlpoint_of_Element, No_Control_point_ON_ELEMENT, Total_Control_Point_to_mesh, Element_patch, Node_Coordinate, INC, Order, Position_Knots, Total_Knot_to_patch_dim, No_knot, No_Control_point)
								* Node_Coordinate[iControlpoint[icc] * (DIMENSION + 1) + 1];
				}

				detJ = sqrt(dxyzdge[0] * dxyzdge[0] + dxyzdge[1] * dxyzdge[1]);
				printf("Check the value of detJ etc: detJ = %f dxyzdge[0] = %f dxyzdge[1] = %f\n", detJ, dxyzdge[0], dxyzdge[1]);
				if (type_load < 2)
				{
					for (ic = 0; ic < (Order[iPatch * DIMENSION + 0] + 1) * (Order[iPatch * DIMENSION + 1] + 1); ic++)
					{
						sfc = Shape_func(ic, Local_Coord, No_Element_for_Integration[iii], Node_Coordinate, Total_Control_Point_to_mesh, Controlpoint_of_Element, INC, Element_patch, Order, No_Control_point_ON_ELEMENT, Position_Knots, Total_Knot_to_patch_dim, No_knot, No_Control_point);
						Equivalent_Nodal_Force[iControlpoint[ic] * DIMENSION + type_load] += valDistLoad * sfc * detJ * Weight[ig];
					}
				}

				if (type_load == 2)
				{
					double LoadDir[2];
					LoadDir[0] = dxyzdge[1] / detJ;
					LoadDir[1] = -dxyzdge[0] / detJ;
					for (ic = 0; ic < (Order[iPatch * DIMENSION + 0] + 1) * (Order[iPatch * DIMENSION + 1] + 1); ic++)
					{
						sfc = Shape_func(ic, Local_Coord, No_Element_for_Integration[iii], Node_Coordinate, Total_Control_Point_to_mesh, Controlpoint_of_Element, INC, Element_patch, Order, No_Control_point_ON_ELEMENT, Position_Knots, Total_Knot_to_patch_dim, No_knot, No_Control_point);
						Equivalent_Nodal_Force[iControlpoint[ic] * DIMENSION + 0] += LoadDir[0] * valDistLoad * sfc * detJ * Weight[ig];
						Equivalent_Nodal_Force[iControlpoint[ic] * DIMENSION + 1] += LoadDir[1] * valDistLoad * sfc * detJ * Weight[ig];
						printf("Equivalent_Nodal_Force[%d][0]=%le\nEquivalent_Nodal_Force[%d][1]=%le\n",
							   iControlpoint[ic], Equivalent_Nodal_Force[iControlpoint[ic] * DIMENSION + 0],
							   iControlpoint[ic], Equivalent_Nodal_Force[iControlpoint[ic] * DIMENSION + 1]);
					}
				}
				if (type_load == 3)
				{
					double LoadDir[2];
					LoadDir[0] = dxyzdge[0] / detJ;
					LoadDir[1] = dxyzdge[1] / detJ;
					for (ic = 0; ic < (Order[iPatch * DIMENSION + 0] + 1) * (Order[iPatch * DIMENSION + 1] + 1); ic++)
					{
						sfc = Shape_func(ic, Local_Coord, No_Element_for_Integration[iii], Node_Coordinate, Total_Control_Point_to_mesh, Controlpoint_of_Element, INC, Element_patch, Order, No_Control_point_ON_ELEMENT, Position_Knots, Total_Knot_to_patch_dim, No_knot, No_Control_point);
						Equivalent_Nodal_Force[iControlpoint[ic] * DIMENSION + 0] += LoadDir[0] * valDistLoad * sfc * detJ * Weight[ig];
						Equivalent_Nodal_Force[iControlpoint[ic] * DIMENSION + 1] += LoadDir[1] * valDistLoad * sfc * detJ * Weight[ig];
					}
				}
			}
		}
	}
	free(No_Element_for_Integration);
}


int SearchForElement(int mesh_n, int iPatch, int iX, int iY, int *Total_Element_on_mesh, int *Total_Element_to_mesh, int *Element_patch, int *ENC)
{
	int iii;

	for (iii = 0; iii < Total_Element_on_mesh[mesh_n]; iii++)
	{
		if (Element_patch[iii + Total_Element_to_mesh[mesh_n]] == iPatch)
		{
			if (iX == ENC[(iii + Total_Element_to_mesh[mesh_n]) * DIMENSION + 0] && iY == ENC[(iii + Total_Element_to_mesh[mesh_n]) * DIMENSION + 1])
				goto loopend;
		}
	}
	loopend:

	return (iii);
}


// for s_IGA, coupled matrix を求める, 要素の重なりを要素のガウス点から求める
void Check_coupled_Glo_Loc_element_for_Gauss(int mesh_n_over, int mesh_n_org)
{
	int re;
	int i, j, k, m;
	int b, l, ll;
	int n_elements_over_point[POW_Ng_extended];
	int patch_n = 0, itr_n = 0;

	int Check_coupled_No[MAX_N_ELEMENT_OVER];
	double Percent_Check_coupled_No;
	int MAX_NNLOVER = 0;

	double element_loc[DIMENSION];

	// int gauss_1dir = 3;	// 重なり判定のための一方向ガウス点数
	// int no_gauss_pt = gauss_1dir * gauss_1dir;	// 重なり判定のためのガウス点総数

	for (i = 0; i < MAX_N_ELEMENT_OVER; i++)
	{
		Check_coupled_No[i] = 0;
	}

	for (m = 0; m < 2; m++) // 最初 Ng 個のガウス点で重なりを求め, NNLOVER[e] >= 2 の e に対して, 再度 Ng_extended 個のガウス点で重なりを求める
	{
		Make_gauss_array(m);

		// ローカルパッチ(mesh_n_over)各要素の頂点の物理座標算出
		for (re = 0; re < real_Total_Element_on_mesh[mesh_n_over]; re++)
		{
			int e = real_element[re + real_Total_Element_to_mesh[mesh_n_over]];

			int i_gg, i_ee;
			int g_n;

			double output_para[DIMENSION];
			int Total_n_elements;

			if (m == 0 || (m == 1 && NNLOVER[e] >= 2))
			{
				if (m == 1)
				{
					NNLOVER[e] = 0;
					for (i = 0; i < NNLOVER[e]; i++)
					{
						NELOVER[e][i] = 0;
					}
				}

				Total_n_elements = 0;
				k = 0;
				ll = 0;

				for (i_ee = 0; i_ee < GP_1dir; i_ee++)
				{
					for (i_gg = 0; i_gg < GP_1dir; i_gg++)
					{
						double data_result_shape[3] = {0.0};

						g_n = i_ee * GP_1dir + i_gg;
						element_loc[0] = Gxi[g_n][0];
						element_loc[1] = Gxi[g_n][1];

						for (b = 0; b < No_Control_point_ON_ELEMENT[Element_patch[e]]; b++)
						{
							double R_shape_func = Shape_func(b, element_loc, e);
							for (j = 0; j < DIMENSION; j++)
							{
								data_result_shape[j] += R_shape_func * Node_Coordinate[Controlpoint_of_Element[e * MAX_NO_CCpoint_ON_ELEMENT + b] * (DIMENSION + 1) + j];
							}
						}

						// 算出したローカルパッチ各要素の頂点の物理座標のグローバルパッチでの(xi, eta)算出
						// from NURBSviewer/NURBS_view/clickcalc.c/func.:calcXiEtaByNR
						for (i = 0; i < Total_Patch_on_mesh[mesh_n_org]; i++)
						{
							int ii = Calc_xi_eta(data_result_shape[0], data_result_shape[1],
												 Position_Knots[i][0], Position_Knots[i][1],
												 No_Control_point[i][0], No_Control_point[i][1], Order[i][0], Order[i][1],
												 &output_para[0], &output_para[1]);
							patch_n = i;
							itr_n = ii;
						}
						// Newton Laphsonによって出力されたxi,etaから重なる要素を求める
						n_elements_over_point[k] = ele_check(patch_n, output_para);
						if (itr_n == 0) // data_result_shapeがグローバルメッシュ上にないとき
						{
							n_elements_over_point[k] = 0;
						}
						Total_n_elements += n_elements_over_point[k];
						for (l = 0; l < n_elements_over_point[k]; l++)
						{
							element_n_point[ll] = temp_element_n[l];
							ll++;
						}
						k++;
					}
				}

				// 昇順ソート
				sort(Total_n_elements);
				// 重複削除
				NNLOVER[e] = duplicate_delete(Total_n_elements, e); // NNLOVER: 要素 e に重なる要素の総数
			
				free(element_n_point), free(temp_element_n)
			}
		}
	}

	for (re = 0; re < real_Total_Element_on_mesh[mesh_n_over]; re++)
	{
		e = real_element[re + real_Total_Element_to_mesh[mesh_n_over]];
		printf("-----------------------------------------Element_No:%d-----------------------------------------\n", e);

		Check_coupled_No[NNLOVER[e]]++;

		if (MAX_NNLOVER < NNLOVER[e])
		{
			MAX_NNLOVER = NNLOVER[e];
		}

		printf("NNLOVER[%d] = %d\n", e, NNLOVER[e]);

		for (i = 0; i < NNLOVER[e]; i++)
		{
			printf("\tNELOVER[%d][%d] = %d\n", e, i, NELOVER[e][i]); //要素eに重なるi番目の要素番号
		}
	}

	printf("MAX_NNLOVER = %d\n", MAX_NNLOVER);
	for (i = 0; i <= MAX_NNLOVER; i++)
	{
		Percent_Check_coupled_No = (double)Check_coupled_No[i] * 100.0 / (double)real_Total_Element_on_mesh[mesh_n_over];
		printf("Check_coupled_No[%d] = %d\t\t%3.1lf %%\n", i, Check_coupled_No[i], Percent_Check_coupled_No);
	}
	printf("---------------------------------------------------------------------------------------------------------------------------\n");
}


void Make_Loc_Glo()
{
	int i, j, k;
	int jj;
	int e;
	int j_n;
	int count;

	j_n = real_Total_Element_to_mesh[Total_mesh] - real_Total_Element_on_mesh[0];
	// printf("j_n=%d\n",j_n);
	// printf("%d\t%d\n",)
	for (i = 0; i < real_Total_Element_on_mesh[0]; i++)
	{
		e = real_element[i];
		count = 0;

		for (j = 0; j < j_n; j++)
		{
			jj = real_element[real_Total_Element_to_mesh[1] + j]; //ローカルメッシュ上のreal element番号
			// printf("jj=%d\n",jj);
			if (NNLOVER[jj] > 0)
			{
				// printf("jj=%d\n",jj);
				for (k = 0; k < NNLOVER[jj]; k++)
				{
					if (NELOVER[jj][k] == e)
					{
						NELOVER[e][count] = jj;
						// printf("NELOVER[%d][%d]=%d\n",e,count,
						//							  NELOVER[e][count]);
						count++;
					}
				}
			}
		}
		NNLOVER[e] = count;
		// printf("NNLOVER[%d]=%d\n",e,NNLOVER[e]);
	}
}


// Newton Raphsonによって出力されたxi,etaから重なる要素を求める
int ele_check(int patch_n, double para_coord[DIMENSION])
{
	int i;
	int j;
	int k, kk;
	int l, ll;
	int RangeCheck_flag;					 //要素を求め終えたら立てるフラグ
	int temp_ad[DIMENSION][MAX_ORDER + 1]; //要素の位置を求めるための値
	int No_line[DIMENSION];					 // xi,etaが含まれている要素の列数
	int n = 1;

	for (j = 0; j < DIMENSION; j++)
	{
		//初期化
		l = 0;
		No_line[j] = 0;
		for (i = 0; i < MAX_ORDER + 1; i++)
		{
			temp_ad[j][i] = 0;
		}
		RangeCheck_flag = 0;

		for (k = 0; k < No_Control_point[patch_n][j] - 1; k++)
		{
			if (RangeCheck_flag == 1)
				break;
			// Local要素の頂点がGlobalパッチ内にない場合
			if (para_coord[j] < Position_Knots[patch_n][j][0] || para_coord[j] > Position_Knots[patch_n][j][No_knot[patch_n][j] - 1])
			{
				// printf("no over element\n");
				RangeCheck_flag++;
			}
			// Local要素の頂点がGlobal要素内部にある場合
			if (para_coord[j] < Position_Knots[patch_n][j][Order[patch_n][j] + k])
			{
				int kk = 0;
				for (kk = 0; kk < k + 1; kk++)
				{
					if (para_coord[j] > Position_Knots[patch_n][j][Order[patch_n][j] + k - kk])
					{
						temp_ad[j][l] = k - kk;
						// printf("temp_ad[%d][%d]=%d\n",j,l,temp_ad[j][l]);
						l++;
						RangeCheck_flag++;
						break;
					}
				}
			}
			// Local要素の頂点がGlobal要素境界上にある場合
			if (para_coord[j] == Position_Knots[patch_n][j][Order[patch_n][j] + k])
			{
				//頂点の座標がGlobalパッチの始点上にある場合
				if (para_coord[j] == Position_Knots[patch_n][j][0])
				{
					temp_ad[j][l] = k;
					l++;
					break;
				}
				//頂点の座標がGlobalパッチの終点上にある場合
				if (para_coord[j] == Position_Knots[patch_n][j][No_knot[patch_n][j] - 1])
				{
					temp_ad[j][l] = k - 1;
					l++;
					break;
				}
				//頂点の座標がGlobal要素境界上にある場合
				else
				{
					temp_ad[j][l] = k - 1;
					l++;
					temp_ad[j][l] = k;
					l++;
				}
				for (kk = 0; kk < Order[patch_n][j]; kk++)
				{
					if (para_coord[j] == Position_Knots[patch_n][j][Order[patch_n][j] + k + kk + 1])
					//多重ノット(次数分ループ)
					{
						printf("C0 continuity\n");
						temp_ad[j][l] = k + kk;
						l++;
					}
					if (para_coord[j] != Position_Knots[patch_n][j][Order[patch_n][j] + k + kk + 1])
						break;
				}
				RangeCheck_flag++;
			}
		}
		No_line[j] = l;
		
		//各方向のNo_lineを掛け合わせる
		n *= l;
	}

	if (DIMENSION == 2)
	{
		for (l = 0; l < No_line[1]; l++)
		{
			for (ll = 0; ll < No_line[0]; ll++)
			{
				temp_element_n[l * No_line[0] + ll] = temp_ad[0][ll] + temp_ad[1][l] * line_No_Total_element[patch_n][0];
			}
		}
	}
	return n;
}


//昇順ソート
void sort(int total)
{
	int i, j;
	int tmp;

	for (i = 0; i < total; i++)
	{
		for (j = i + 1; j < total; j++)
		{
			if (element_n_point[i] > element_n_point[j])
			{
				tmp = element_n_point[i];
				element_n_point[i] = element_n_point[j];
				element_n_point[j] = tmp;
			}
		}
	}
}


//重複削除
int duplicate_delete(int total, int element_n)
{
	int i, j;

	j = 0;
	NELOVER[element_n][j] = element_n_point[0];
	j++;
	for (i = 1; i < total; i++)
	{
		if (element_n_point[i] != element_n_point[i - 1])
		{
			NELOVER[element_n][j] = element_n_point[i];
			j++;
		}
	}
	// j = 要素element_nに重なる要素の総数
	return j;
}




// Preprocessing








// Shape Function
double Shape_func(int I_No, double Local_coord[DIMENSION], int El_No, double *Node_Coordinate, int *Total_Control_Point_to_mesh,
				  int *Controlpoint_of_Element, int *INC, int *Element_patch, int *Order, int *No_Control_point_ON_ELEMENT,
				  double *Position_Knots, int *Total_Knot_to_patch_dim, int *No_knot, int *No_Control_point)
{
	int i, j;
	double R, weight_func = 0.0;

	double Position_Data_param[DIMENSION];
	double shape_func[MAX_NO_CCpoint_ON_ELEMENT * DIMENSION];

	double *Shape = (double *)malloc(sizeof(double) * DIMENSION * Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER); // Shape[DIMENSION][MAX_N_NODE][10]
	double *dShape = (double  *)malloc(sizeof(double) * DIMENSION * Total_Control_Point_to_mesh[Total_mesh]); // dShape[DIMENSION][MAX_N_NODE]

	for (i = 0; i < MAX_NO_CCpoint_ON_ELEMENT; i++)
	{
		shape_func[i] = 1.0;
	}

	for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			ShapeFunc_from_paren(Position_Data_param, Local_coord, j, El_No, INC, Position_Knots, Total_Knot_to_patch_dim, Controlpoint_of_Element, Element_patch);
			ShapeFunction1D(Position_Data_param, j, El_No, Shape, dShape, No_knot, Total_Control_Point_to_mesh, Position_Knots, Total_Knot_to_patch_dim, Element_patch, Order, No_Control_point);
			shape_func[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i]] *= Shape[j * (Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + INC[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + j] * MAX_ORDER + Order[Element_patch[El_No] * DIMENSION + j]];
		}
		weight_func += shape_func[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i]] * Node_Coordinate[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * (DIMENSION + 1) + DIMENSION];
	}

	free(Shape), free(dShape);

	if (I_No < No_Control_point_ON_ELEMENT[Element_patch[El_No]])
		R = shape_func[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + I_No]] * Node_Coordinate[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + I_No] * (DIMENSION + 1) + DIMENSION] / weight_func;

	else
		R = ERROR;

	return R;
}


void ShapeFunc_from_paren(double Position_Data_param[DIMENSION], double Local_coord[DIMENSION], int j, int e, int *INC,
						  double *Position_Knots, int *Total_Knot_to_patch_dim, int *Controlpoint_of_Element, int *Element_patch)
{
	int i = 0;

	i = INC[Controlpoint_of_Element[e * MAX_NO_CCpoint_ON_ELEMENT + 0] * DIMENSION + j];
	Position_Data_param[j] = ((Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + i + 1] - Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + i]) * Local_coord[j] + (Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + i + 1] + Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + i])) / 2.0;
}


void ShapeFunction1D(double Position_Data_param[DIMENSION], int j, int e, double *Shape, double *dShape, int *No_knot, int *Total_Control_Point_to_mesh,
					 double *Position_Knots, int *Total_Knot_to_patch_dim,int *Element_patch, int *Order, int *No_Control_point)
{
	int ii;
	int p;

	for (ii = 0; ii < No_knot[Element_patch[e] * DIMENSION + j]; ii++)
	{
		if (Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + ii] == Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + ii + 1])
		{
			Shape[j * (Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + 0] = 0.0;
		}
		else if (Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + ii] != Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + ii + 1] && Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + ii] <= Position_Data_param[j] && Position_Data_param[j] < Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + ii + 1])
		{
			Shape[j * (Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + 0] = 1.0;
		}
		else if (Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + ii] != Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + ii + 1] && Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + ii + 1] == Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + (No_knot[Element_patch[e] * DIMENSION + j] - 1)] && Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + ii] <= Position_Data_param[j] && Position_Data_param[j] <= Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + ii + 1])
		{
			Shape[j * (Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + 0] = 1.0;
		}
		else
		{
			Shape[j * (Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + 0] = 0.0;
		}
	}

	for (ii = 0; ii < No_knot[Element_patch[e] * DIMENSION + j]; ii++)
	{
		for (p = 1; p <= Order[Element_patch[e] * DIMENSION + j]; p++)
		{
			Shape[j * (Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + p] = 0.0;
		}
	}

	double left_term, right_term;
	for (p = 1; p <= Order[Element_patch[e] * DIMENSION + j]; p++)
	{
		for (ii = 0; ii < No_knot[Element_patch[e] * DIMENSION + j]; ii++)
		{
			left_term = 0.0;
			right_term = 0.0;

			if ((Position_Data_param[j] - Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + ii]) * Shape[j * (Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + p - 1] == 0 && Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + ii + p] - Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + ii] == 0)
				left_term = 0.0;
			else
			{
				left_term = (Position_Data_param[j] - Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + ii]) / (Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + ii + p] - Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + ii]) * Shape[j * (Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + p - 1];
			}
			if ((Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + ii + p + 1] - Position_Data_param[j]) * Shape[j * (Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + (ii + 1) * MAX_ORDER + p - 1] == 0 && Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + ii + p + 1] - Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + ii + 1] == 0)
				right_term = 0.0;
			else
			{
				right_term = (Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + ii + p + 1] - Position_Data_param[j]) / (Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + ii + p + 1] - Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + ii + 1]) * Shape[j * (Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + (ii + 1) * MAX_ORDER + p - 1];
			}
			Shape[j * (Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + p] = left_term + right_term;
		}
	}

	double dleft_term, dright_term;
	for (ii = 0; ii < No_Control_point[Element_patch[e] * DIMENSION + j] + 1; ii++)
	{
		dleft_term = 0.0;
		dright_term = 0.0;

		if (Order[Element_patch[e] * DIMENSION + j] * Shape[j * (Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + Order[Element_patch[e] * DIMENSION + j] - 1] == 0 && Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + ii + Order[Element_patch[e] * DIMENSION + j]] - Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + ii] == 0)
			dleft_term = 0.0;
		else
			dleft_term = Order[Element_patch[e] * DIMENSION + j] / (Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + ii + Order[Element_patch[e] * DIMENSION + j]] - Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + ii]) * Shape[j * (Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + Order[Element_patch[e] * DIMENSION + j] - 1];

		if (Order[Element_patch[e] * DIMENSION + j] * Shape[j * (Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + (ii + 1) * MAX_ORDER + Order[Element_patch[e] * DIMENSION + j] - 1] == 0 && Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + ii + Order[Element_patch[e] * DIMENSION + j] + 1] - Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + ii + 1] == 0)
			dright_term = 0.0;
		else
			dright_term = Order[Element_patch[e] * DIMENSION + j] / (Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + ii + Order[Element_patch[e] * DIMENSION + j] + 1] - Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + ii + 1]) * Shape[j * (Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + (ii + 1) * MAX_ORDER + Order[Element_patch[e] * DIMENSION + j] - 1];

		dShape[j * Total_Control_Point_to_mesh[Total_mesh] + ii] = dleft_term - dright_term;
	}
}


double dShape_func(int I_No, int xez, double Local_coord[DIMENSION], int El_No, int *Controlpoint_of_Element, int *No_Control_point_ON_ELEMENT,
				   int *Total_Control_Point_to_mesh, int *Element_patch, int *Node_Coordinate, int *INC, int *Order, double *Position_Knots,
				   int *Total_Knot_to_patch_dim, int *No_knot, int *No_Control_point)
{
	double dR;

	double *dShape_func1 = (double *)malloc(sizeof(double) * Total_Control_Point_to_mesh[Total_mesh]); // dShape_func1[MAX_N_NODE];
	double *dShape_func2 = (double *)malloc(sizeof(double) * Total_Control_Point_to_mesh[Total_mesh]); // dShape_func2[MAX_N_NODE];

	NURBS_deriv(Local_coord, El_No, Node_Coordinate, Total_Control_Point_to_mesh, Controlpoint_of_Element, INC, Element_patch, Order, No_Control_point_ON_ELEMENT,
				Position_Knots, Total_Knot_to_patch_dim, No_knot, No_Control_point, dShape_func1, dShape_func2);

	if (xez != 0 && xez != 1)
		dR = ERROR;

	else if (I_No < No_Control_point_ON_ELEMENT[Element_patch[El_No]])
	{
		if (xez == 0)
		{
			dR = dShape_func1[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + I_No]]
			   * dShapeFunc_from_paren(xez, El_No, INC, Controlpoint_of_Element, Position_Knots, Position_Knots, Total_Knot_to_patch_dim, Element_patch);
		}
		else if (xez == 1)
		{
			dR = dShape_func2[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + I_No]]
			   * dShapeFunc_from_paren(xez, El_No, INC, Controlpoint_of_Element, Position_Knots, Position_Knots, Total_Knot_to_patch_dim, Element_patch);
		}
	}
	else
		dR = ERROR;

	free(dShape_func1), free(dShape_func2);
	return dR;
}


void NURBS_deriv(double Local_coord[DIMENSION], int El_No, double *Node_Coordinate, int *Total_Control_Point_to_mesh,
				 int *Controlpoint_of_Element, int *INC, int *Element_patch, int *Order, int *No_Control_point_ON_ELEMENT,
				 double *Position_Knots, int *Total_Knot_to_patch_dim, int *No_knot, int *No_Control_point,
				 double *dShape_func1, double *dShape_func2)
{
	int i, j;
	double weight_func = 0.0;
	double dWeight_func1 = 0.0;
	double dWeight_func2 = 0.0;

	double Position_Data_param[DIMENSION];
	double shape_func[MAX_NO_CCpoint_ON_ELEMENT * DIMENSION];

	double *Shape = (double *)malloc(sizeof(double) * DIMENSION * Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER); // Shape[DIMENSION][MAX_N_NODE][10]
	double *dShape = (double  *)malloc(sizeof(double) * DIMENSION * Total_Control_Point_to_mesh[Total_mesh]); // dShape[DIMENSION][MAX_N_NODE]

	for (i = 0; i < MAX_NO_CCpoint_ON_ELEMENT; i++)
	{
		shape_func[i] = 1.0;
	}

	for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			ShapeFunc_from_paren(Position_Data_param, Local_coord, j, El_No, INC, Position_Knots, Total_Knot_to_patch_dim, Controlpoint_of_Element, Element_patch);
			ShapeFunction1D(Position_Data_param, j, El_No, Shape, dShape, No_knot, Total_Control_Point_to_mesh, Position_Knots, Total_Knot_to_patch_dim, Element_patch, Order, No_Control_point);
			shape_func[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i]] *= Shape[j * (Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + INC[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + j] * MAX_ORDER + Order[Element_patch[El_No] * DIMENSION + j]];
		}
		weight_func += shape_func[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i]] * Node_Coordinate[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * (DIMENSION + 1) + DIMENSION];
	}
	for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++)
	{
		dWeight_func1 += dShape[0 * Total_Control_Point_to_mesh[Total_mesh] + INC[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + 0]] * Shape[1 * (Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + INC[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + 1] * MAX_ORDER + Order[Element_patch[El_No] * DIMENSION + 1]] * Node_Coordinate[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i]][DIMENSION];
		dWeight_func2 += Shape[0 * (Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + INC[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + 0] * MAX_ORDER + Order[Element_patch[El_No] * DIMENSION + 0]] * dShape[1 * Total_Control_Point_to_mesh[Total_mesh] + INC[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + 1]] * Node_Coordinate[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i]][DIMENSION];
	}
	for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++)
	{
		dShape_func1[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i]] = Node_Coordinate[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * (DIMENSION + 1) + DIMENSION] * (weight_func * dShape[0 * Total_Control_Point_to_mesh[Total_mesh] + INC[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + 0]] * Shape[1 * (Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + INC[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + 1] * MAX_ORDER + Order[Element_patch[El_No] * DIMENSION + 1]] - dWeight_func1 * shape_func[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i]]) / (weight_func * weight_func);
		dShape_func2[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i]] = Node_Coordinate[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * (DIMENSION + 1) + DIMENSION] * (weight_func * Shape[0 * (Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + INC[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + 0] * MAX_ORDER + Order[Element_patch[El_No] * DIMENSION + 0]] * dShape[1 * Total_Control_Point_to_mesh[Total_mesh] + INC[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + 1]] - dWeight_func2 * shape_func[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i]]) / (weight_func * weight_func);
	}

	free(Shape), free(dShape);
}


double dShapeFunc_from_paren(int j, int e, int *INC, int *Controlpoint_of_Element, double *Position_Knots,
							 double *Position_Knots, int *Total_Knot_to_patch_dim, int *Element_patch)
{
	int i;
	double dPosition_Data_param;

	i = INC[Controlpoint_of_Element[e * MAX_NO_CCpoint_ON_ELEMENT + 0] * DIMENSION + j];

	dPosition_Data_param = (Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + i + 1] - Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + i]) / 2.0;
	return dPosition_Data_param;
}









// 拘束されている行数を省いた行列の番号の制作
int Make_Index_Dof(int Total_Control_Point, int Total_Constraint, int Constraint_Node_Dir[MAX_N_CONSTRAINT][2])
{
	int i, k = 0;

	// Index_Dofの初期化(複数メッシュ読み込みのため)
	for (i = 0; i < Total_Control_Point * 2; i++)
	{
		Index_Dof[i] = 0;
	}
	// 拘束されている自由度(Degree Of free)をERRORにする
	for (i = 0; i < Total_Constraint; i++)
	{
		Index_Dof[Constraint_Node_Dir[i][0] * DIMENSION + Constraint_Node_Dir[i][1]] = ERROR;
	}
	// ERROR以外に番号を付ける
	for (i = 0; i < Total_Control_Point * DIMENSION; i++)
	{
		if (Index_Dof[i] != ERROR)
		{
			Index_Dof[i] = k;
			k++;
		}
	}
	printf("Max_Index_Dof = %d\n", k);
	return k;
}

void Make_K_Whole_Ptr_Col(int Total_Element,
						  int Total_Control_Point,
						  int K_Whole_Size)
{
	int i, ii, j, jj, k;
	int NE;
	int N, i_index, j_index;

	//初期化
	// for (i = 0; i < Total_Control_Point * DIMENSION; i++)
	// Total_Control_Point_To_Node[i] = 0;
	for (i = 0; i < K_Whole_Size + 1; i++)
		K_Whole_Ptr[i] = 0;

	for (N = 0; N < Total_Control_Point; N += K_DIVISION_LENGE)
	{ //大きく分割するためのループ
		//各節点に接する節点を取得
		for (i = 0; i < K_DIVISION_LENGE; i++)
		{
			Total_Control_Point_To_Node[i] = 0;
		}
		for (i = 0; i < Total_Element; i++)
		{
			for (ii = 0; ii < No_Control_point_ON_ELEMENT[Element_patch[i]]; ii++)
			{
				NE = Controlpoint_of_Element[i * MAX_NO_CCpoint_ON_ELEMENT + ii] - N;
				// printf("NE=%d\n",NE);
				// printf("K_DIVISION_LENGE=%d,N=%d,NE=%d\n",K_DIVISION_LENGE,N,NE);    //K_DIVISION_LENGE=0,N=0,NE=コネクティビティ的な
				if (0 <= NE && NE < K_DIVISION_LENGE)
				{
					for (j = 0; j < No_Control_point_ON_ELEMENT[Element_patch[i]]; j++) //ローカル要素
					{
						// printf("j=%d\n",j);
						//数字がない時
						if (Total_Control_Point_To_Node[NE] == 0)
						{
							//節点番号を取得
							Node_To_Node[NE][0] = Controlpoint_of_Element[i * MAX_NO_CCpoint_ON_ELEMENT + j];
							Total_Control_Point_To_Node[NE]++;
							// printf("Node_To_Node[%d][0]=%d\n",NE,Node_To_Node[NE][0]);
						}
						// printf("②Total_Control_Point_To_Node[%d]=%d\n",NE,Total_Control_Point_To_Node[NE]);
						//同じものがあったら
						// k > 0 以降の取得
						// kのカウント
						for (k = 0; k < Total_Control_Point_To_Node[NE]; k++)
						{
							// printf("Node_To_Node[%d][%d]=%d\n",NE,k,Node_To_Node[NE][k]);
							//
							// printf("k_1=%d\t",k);
							if (Node_To_Node[NE][k] == Controlpoint_of_Element[i * MAX_NO_CCpoint_ON_ELEMENT + j])
							{
								// printf("break\t");
								break;
							}
						}
						// printf("\nk_2=%d\n",k);
						//未設定のNode_To_Node取得
						if (k == Total_Control_Point_To_Node[NE])
						{
							Node_To_Node[NE][k] = Controlpoint_of_Element[i * MAX_NO_CCpoint_ON_ELEMENT + j];
							// printf("Node_To_Node[%d][%d]=%d\n",NE,k,Node_To_Node[NE][k]);
							Total_Control_Point_To_Node[NE]++;
							// printf("③Total_Control_Point_To_Node[%d]=%d\n",NE,Total_Control_Point_To_Node[NE]);
						}
					}
					//別メッシュとの重なりを考慮
					if (NNLOVER[i] > 0)
					{
						for (jj = 0; jj < NNLOVER[i]; jj++)
						{
							for (j = 0; j < No_Control_point_ON_ELEMENT[Element_patch[NELOVER[i][jj]]]; j++) //ローカル要素
							{
								// printf("j=%d\n",j);
								//数字がない時
								if (Total_Control_Point_To_Node[NE] == 0)
								{
									//節点番号を取得
									Node_To_Node[NE][0] = Controlpoint_of_Element[NELOVER[i][jj] * MAX_NO_CCpoint_ON_ELEMENT + j];
									Total_Control_Point_To_Node[NE]++;
									// printf("Node_To_Node[%d][0]=%d\n",NE,Node_To_Node[NE][0]);
								}
								// printf("②Total_Control_Point_To_Node[%d]=%d\n",NE,Total_Control_Point_To_Node[NE]);
								//同じものがあったら
								// k > 0 以降の取得
								// kのカウント
								for (k = 0; k < Total_Control_Point_To_Node[NE]; k++)
								{
									// printf("Node_To_Node[%d][%d]=%d\n",NE,k,Node_To_Node[NE][k]);
									//
									// printf("k_1=%d\t",k);
									if (Node_To_Node[NE][k] == Controlpoint_of_Element[NELOVER[i][jj] * MAX_NO_CCpoint_ON_ELEMENT + j])
									{
										// printf("break\t");
										break;
									}
								}
								// printf("\nk_2=%d\n",k);
								//未設定のNode_To_Node取得
								if (k == Total_Control_Point_To_Node[NE])
								{
									Node_To_Node[NE][k] = Controlpoint_of_Element[NELOVER[i][jj] * MAX_NO_CCpoint_ON_ELEMENT + j];
									// printf("Node_To_Node[%d][%d]=%d\n",NE,k,Node_To_Node[NE][k]);
									Total_Control_Point_To_Node[NE]++;
									// printf("③Total_Control_Point_To_Node[%d]=%d\n",NE,Total_Control_Point_To_Node[NE]);
								}
							}
						}
					}
				}
				// printf("\n");
			}
			// printf("\n");
		}
		//順番に並び替える
		for (i = 0; i < K_DIVISION_LENGE; i++)
		{
			if (N + i < Total_Control_Point)
			{
				// printf("Node[%d] T=%d; \n",N+i, Total_Control_Point_To_Node[ i ]);
				for (j = 0; j < Total_Control_Point_To_Node[i]; j++)
				{
					int Min = Node_To_Node[i][j], No = j;
					for (k = j; k < Total_Control_Point_To_Node[i]; k++)
					{
						if (Min > Node_To_Node[i][k])
						{
							Min = Node_To_Node[i][k];
							No = k;
						}
					}
					for (k = No; k > j; k--)
					{
						Node_To_Node[i][k] = Node_To_Node[i][k - 1];
					}
					Node_To_Node[i][j] = Min;
					//				printf("%d ",Node_To_Node[i][j]);
				}
				//			printf("\n");
			}

			//並べ替えたNode_To_Node確認
			// for (j = 0; j < Total_Control_Point_To_Node[i]; j++)
			// {
			// 	printf("sort_Node_To_Node[%d][%d]=%d\n",i,j,Node_To_Node[i][j]);
			// }
		}

		//節点からcol ptrを求める
		ii = 0;
		k = 0;
		for (i = 0; i < K_DIVISION_LENGE; i++)
		{
			for (ii = 0; ii < DIMENSION; ii++)
			{
				if (N + i < Total_Control_Point)
				{
					i_index = Index_Dof[(N + i) * DIMENSION + ii];
					// printf("i = %d\n", i);
					// printf("N = %d\n", N);
					// printf("i_index = %d\n", i_index);
					k = 0;
					if (i_index >= 0)
					{
						// K_Whole_Ptr[i_index + 1] = K_Whole_Ptr[i_index];
						K_Whole_Ptr[i_index + 1] = K_Whole_Ptr[i_index];
						// printf("K_Whole_Ptr[%d][%d]=%d\n",tm,i_index,K_Whole_Ptr[tm][i_index+1]);
						for (j = 0; j < Total_Control_Point_To_Node[i]; j++)
						{
							// printf("Total_Control_Point_To_Node[%d] = %d\n", i, Total_Control_Point_To_Node[i]);
							for (jj = 0; jj < DIMENSION; jj++)
							{
								j_index = Index_Dof[Node_To_Node[i][j] * DIMENSION + jj];
								if (j_index >= 0 && j_index >= i_index)
								{
									K_Whole_Ptr[i_index + 1]++;
									// col_N[N/K_DIVISION_LENGE][k] = j_index;
									K_Whole_Col[K_Whole_Ptr[i_index] + k] = j_index;
									// printf("K_Whole_Col[%d]=%d\n"
									//        ,K_Whole_Ptr[i_index]+k
									//        ,K_Whole_Col[K_Whole_Ptr[i_index]+k]);
									k++;
									// printf("ptr[%d]=%d,col[%d]=%d\n",i_index+1,K_Whole_Ptr[i_index+1],K_Whole_Ptr[i_index]+k,K_Whole_Col[K_Whole_Ptr[i_index]+k]);
								}
							}
						}
					}
				}
			}
		}
		// col_N[N/K_DIVISION_LENGE][ k ] = -1;
	}
	// for(i=0;i<K_Whole_Size_array[tm]+1;i++)
	/*
	for(i=0;i<K_Whole_Size+1;i++)
	{
		printf("K_Whole_Ptr[%d]=%d\n",
				i,K_Whole_Ptr[i]);
	}*/
	/*
	for( i = 0; i < K_Whole_Size+1; i++ )//printf("K_Whole_Ptr[%d]= %d\n",i,K_Whole_Ptr[i]);
	//col合成
	k = 0;
	for( N = 0; N < Total_Control_Point ; N +=K_DIVISION_LENGE ){
		for(i = 0; col_N[ N/K_DIVISION_LENGE ][i] != -1; i++ ){
			K_Whole_Col[k] = col_N[ N/K_DIVISION_LENGE ][i];
			k++;
		}
	}
	*/
}

// valを求める
void Make_K_Whole_Val(double E, double nu, int Total_Element, int DM)
{
	int i, j, j1, j2, k1, k2, l;
	int a, b, re;

	for (i = 0; i < MAX_NON_ZERO; i++)
	{
		K_Whole_Val[i] = 0.0;
	}

	for (re = 0; re < Total_Element; re++)
	{
		i = real_element[re];
		Check_BDBJ_flag[i] = 0;
	}

	for (re = 0; re < Total_Element; re++)
	{
		// if(i==real_element[re]){
		// printf("re=%d\n",re);
		i = real_element[re];
		// printf("El_No;i=%d\n", real_element[re]);

		if (Element_mesh[i] == 0 && re == 0) // 2つめの条件は効率化のため
		{
			Make_gauss_array(0);
		}
		else if (Element_mesh[i] > 0)
		{
			printf("NNLOVER[%d]:%d\tNNLOVER[%d]:%d\tElement_mesh[%d]:%d\n", i, NNLOVER[i], real_element[re - 1], NNLOVER[real_element[re - 1]], real_element[re - 1], Element_mesh[real_element[re - 1]]);
			if (NNLOVER[i] == 1 && (NNLOVER[real_element[re - 1]] != 1 || Element_mesh[real_element[re - 1]] == 0)) /*2つめ以降の条件は効率化のため*/
			{
				Make_gauss_array(0);
			}
			else if (NNLOVER[i] >= 2 && (NNLOVER[real_element[re - 1]] == 1 || Element_mesh[real_element[re - 1]] == 0)) /*2つめ以降の条件は効率化のため*/
			{
				Make_gauss_array(1);
			}
		}
		// printf("i= %d\tGaussPt_3D=%d\n",i ,GaussPt_3D);

		for (j = 0; j < GP_2D; j++)
		{
			Same_BDBJ_flag[j] = 0;
		}

		KIEL_SIZE = No_Control_point_ON_ELEMENT[Element_patch[i]] * DIMENSION;
		double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], K_EL[MAX_KIEL_SIZE][MAX_KIEL_SIZE];
		// printf("Total_Element=%d\tre=%d\tEl_No=%d\n", Total_Element, re, i);
		//各要素のKelを求める
		for (j1 = 0; j1 < No_Control_point_ON_ELEMENT[Element_patch[i]]; j1++)
		{
			for (j2 = 0; j2 < DIMENSION; j2++)
			{
				X[j1][j2] = Node_Coordinate[Controlpoint_of_Element[i * MAX_NO_CCpoint_ON_ELEMENT + j1] * (DIMENSION + 1) + j2];
			}
		}

		Make_K_EL(i, X, K_EL, E, nu, DM);

		// Valを求める
		// for (j1 = 0; j1 < No_Control_point_ON_ELEMENT[Element_patch[El_No_on_mesh[tm][i]]]; j1++)
		for (j1 = 0; j1 < No_Control_point_ON_ELEMENT[Element_patch[i]]; j1++)
		{
			for (j2 = 0; j2 < DIMENSION; j2++)
			{
				a = Index_Dof[Controlpoint_of_Element[i * MAX_NO_CCpoint_ON_ELEMENT + j1] * DIMENSION + j2];
				if (a >= 0)
				{
					// for (k1 = 0; k1 < No_Control_point_ON_ELEMENT[Element_patch[El_No_on_mesh[tm][i]]]; k1++)
					for (k1 = 0; k1 < No_Control_point_ON_ELEMENT[Element_patch[i]]; k1++)
					{
						for (k2 = 0; k2 < DIMENSION; k2++)
						{
							b = Index_Dof[Controlpoint_of_Element[i * MAX_NO_CCpoint_ON_ELEMENT + k1] * DIMENSION + k2];
							if (b >= 0 && b >= a)
							{
								for (l = K_Whole_Ptr[a]; l < K_Whole_Ptr[a + 1]; l++)
								{
									if (K_Whole_Col[l] == b)
									{
										K_Whole_Val[l] += K_EL[j1 * DIMENSION + j2][k1 * DIMENSION + k2];
										break;
									}
								}
							}
						}
					}
				}
			}
		}

		if (Element_mesh[i] > 0) //ローカルメッシュ上の要素について
		{
			if (NNLOVER[i] > 0) //重なっている要素が存在するとき
			{
				for (j = 0; j < NNLOVER[i]; j++)
				{
					double XG[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION];
					KIEL_SIZE = No_Control_point_ON_ELEMENT[Element_patch[NELOVER[i][j]]] * DIMENSION;
					double coupled_K_EL[MAX_KIEL_SIZE][MAX_KIEL_SIZE];
					for (j1 = 0; j1 < No_Control_point_ON_ELEMENT[Element_patch[NELOVER[i][j]]]; j1++)
					{
						for (j2 = 0; j2 < DIMENSION; j2++)
						{
							XG[j1][j2] = Node_Coordinate[Controlpoint_of_Element[NELOVER[i][j] * MAX_NO_CCpoint_ON_ELEMENT + j1] * (DIMENSION + 1) + j2];
							//重なっている要素の物理座標取得
						}
					}
					Make_coupled_K_EL(i, NELOVER[i][j],
									  X,
									  XG,
									  coupled_K_EL,
									  E, nu, DM);

					Check_BDBJ_flag[i] += Total_BDBJ_flag;
					if (j == NNLOVER[i] - 1)
					{
						for (j1 = 0; j1 < GP_2D; j1++)
						{
							// printf("Same_BDBJ_flag[%d]=%d\n",j1,Same_BDBJ_flag[j1]);
							if (Same_BDBJ_flag[j1] != 1)
							{
								printf("ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR\n");
							}
						}
						printf("-------------------------Check_BDBJ_flag[%d]=%d-------------------------\n", i, Check_BDBJ_flag[i]);
						if (Check_BDBJ_flag[i] != GP_2D)
						{
							printf("ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR\n");
						}
					}

					// Valを求める
					for (j1 = 0; j1 < No_Control_point_ON_ELEMENT[Element_patch[NELOVER[i][j]]]; j1++)
					{
						for (j2 = 0; j2 < DIMENSION; j2++)
						{
							a = Index_Dof[Controlpoint_of_Element[NELOVER[i][j] * MAX_NO_CCpoint_ON_ELEMENT + j1] * DIMENSION + j2];
							if (a >= 0)
							{
								for (k1 = 0; k1 < No_Control_point_ON_ELEMENT[Element_patch[i]]; k1++)
								{
									for (k2 = 0; k2 < DIMENSION; k2++)
									{
										b = Index_Dof[Controlpoint_of_Element[i * MAX_NO_CCpoint_ON_ELEMENT + k1] * DIMENSION + k2];
										if (b >= 0 && b >= a)
										{
											for (l = K_Whole_Ptr[a]; l < K_Whole_Ptr[a + 1]; l++)
											{
												if (K_Whole_Col[l] == b)
												{
													K_Whole_Val[l] += coupled_K_EL[j1 * DIMENSION + j2][k1 * DIMENSION + k2];
													break;
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}


//分布荷重の等価節点力を足す
void Add_Equivalent_Nodal_Forec_to_F_Vec(int *Total_Control_Point, double *Equivalent_Nodal_Force)
{
	int i, j, index;
	for (j = 0; j < DIMENSION; j++)
	{
		for (i = 0; i < Total_Control_Point[Total_mesh]; i++)
		{
			index = Index_Dof[i * DIMENSION + j];
			if (index >= 0)
			{
				rhs_vec[index] += Equivalent_Nodal_Force[i * DIMENSION + j];
			}
		}
	}
}

//荷重の行列を作る
void Make_F_Vec(int Total_Load, int Load_Node_Dir[MAX_N_LOAD][2], double Value_of_Load[MAX_N_LOAD], int K_Whole_Size)
{
	int i, index;
	for (i = 0; i < K_Whole_Size; i++)
		rhs_vec[i] = 0.0;
	for (i = 0; i < Total_Load; i++)
	{
		index = Index_Dof[Load_Node_Dir[i][0] * DIMENSION + Load_Node_Dir[i][1]];
		if (index >= 0)
			rhs_vec[index] += Value_of_Load[i];
	}
}

//強制変位対策
void Make_F_Vec_disp_const(int Mesh_No, int Total_Constraint,
						   int Constraint_Node_Dir[MAX_N_CONSTRAINT][2],
						   double Value_of_Constraint[MAX_N_CONSTRAINT],
						   double E, double nu, int DM)
{
	int ie, idir, inode, jdir, jnode, kk_const;
	int ii, iii, b, bb, jj, j1, j2, ii_local, jj_local;
	int iee;

	int i;

	Make_gauss_array(0);

	for (ie = 0; ie < real_Total_Element_to_mesh[Total_mesh]; ie++)
	{
		i = real_element[ie];

		KIEL_SIZE = No_Control_point_ON_ELEMENT[Element_patch[i]] * DIMENSION;
		double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], K_EL[MAX_KIEL_SIZE][MAX_KIEL_SIZE];

		iii = 0;
		for (idir = 0; idir < DIMENSION; idir++)
		{
			for (inode = 0; inode < No_Control_point_ON_ELEMENT[Element_patch[i]]; inode++)
			{
				b = Index_Dof[Controlpoint_of_Element[i * MAX_NO_CCpoint_ON_ELEMENT + inode] * DIMENSION + idir];
				if (b < 0)
					iii++;
			}
		}
		// printf("iii;%d\n",iii);
		if (iii > 0)
		{
			for (j1 = 0; j1 < No_Control_point_ON_ELEMENT[Element_patch[i]]; j1++)
			{
				for (j2 = 0; j2 < DIMENSION; j2++)
				{
					X[j1][j2] = Node_Coordinate[Controlpoint_of_Element[i * MAX_NO_CCpoint_ON_ELEMENT + j1] * (DIMENSION + 1) + j2];
				} // end for j2
			}	  // end for j1
			iee = i;
			Make_K_EL(iee, X, K_EL, E, nu, DM);
			for (idir = 0; idir < DIMENSION; idir++)
			{
				for (inode = 0; inode < No_Control_point_ON_ELEMENT[Element_patch[i]]; inode++)
				{
					ii = Controlpoint_of_Element[real_El_No_on_mesh[Mesh_No][ie] * MAX_NO_CCpoint_ON_ELEMENT + inode] * DIMENSION + idir;
					b = Index_Dof[ii];
					if (b >= 0)
					{
						ii_local = inode * DIMENSION + idir;
						for (jdir = 0; jdir < DIMENSION; jdir++)
						{
							for (jnode = 0; jnode < No_Control_point_ON_ELEMENT[Element_patch[i]]; jnode++)
							{
								jj = Controlpoint_of_Element[i * MAX_NO_CCpoint_ON_ELEMENT + jnode] * DIMENSION + jdir;
								bb = Index_Dof[jj];
								if (bb < 0)
								{
									jj_local = jnode * DIMENSION + jdir; // printf("%d,%d\n",ie,jnode);
									for (kk_const = 0; kk_const < Total_Constraint; kk_const++)
									{
										if (Controlpoint_of_Element[i * MAX_NO_CCpoint_ON_ELEMENT + jnode] == Constraint_Node_Dir[kk_const][0] && jdir == Constraint_Node_Dir[kk_const][1])
										{
											rhs_vec[b] -= K_EL[ii_local][jj_local] * Value_of_Constraint[kk_const]; // if(kk_const >= 28){printf("%d , %d ,%16.15e\n",ii_local, jj_local ,  K_EL[ii_local][jj_local]);}
										}																			// end if Controlpoint_of_Element[ie][jnode]
									}																				// end for kk_const
								}																					// end if bb
							}																						// end for jnode
						}																							// end for jdir
					}																								// end if b>=0
				}																									// end for inode
			}																										// end for idir
		}																											// end if iii>0
	}																												// end for ie
} // end

void mat_vec_crs(double *vec_result, double *vec, const int ndof)
{
	int i, j, icount = 0;

	// zero clear
	for (i = 0; i < ndof; i++)
		vec_result[i] = 0;
	
	for (i = 0; i < ndof; i++)
	{
		for (j = K_Whole_Ptr[i]; j < K_Whole_Ptr[i + 1]; j++)
		{
			vec_result[i] += K_Whole_Val[icount] * vec[K_Whole_Col[j]];
			if (i != K_Whole_Col[j])
				vec_result[K_Whole_Col[j]] += K_Whole_Val[icount] * vec[i];
			icount++;
		}
	}
}

double inner_product(int ndof, double *vec1, double *vec2)
{
	double rrr = 0.0;
	int i;
	for (i = 0; i < ndof; i++)
	{
		rrr += vec1[i] * vec2[i];
	}
	return (rrr);
}

int check_conv_CG(int ndof, double alphak, double *pp, double eps, int itr)
{
	double rrr1 = 0.0, rrr2 = 0.0, rrr3;
	int i, istop = 0;

	printf("ndof=%d alphak= %15e\t", ndof, alphak);
	for (i = 0; i < ndof; i++)
	{
		rrr1 += pp[i] * pp[i];
		rrr2 += sol_vec[i] * sol_vec[i];
	}
	rrr3 = fabs(alphak) * sqrt(rrr1 / rrr2);
	printf("Iteration# = %d  residual = %15e (%15e)\n", itr, rrr3, eps);
	if (rrr3 < eps)
		istop = 1;

	return (istop);
}

void Diag_Scaling_CG_pre(int ndof, int flag_operation)
{
	int i, j;
	int icount = 0;

	printf("ndof=%d\n", ndof);
	if (flag_operation == 0)
	{
		diag_scaling[0] = 1.0 / sqrt(K_Whole_Val[0]);

		for (i = 1; i < ndof; i++)
		{
			diag_scaling[i] = 1.0 / sqrt(K_Whole_Val[K_Whole_Ptr[i]]);
			printf("diag=%le\n", diag_scaling[i]);
			printf("K_Whole_Val[%d] = %.16e\n", K_Whole_Ptr[i], K_Whole_Val[K_Whole_Ptr[i]]);
			printf("sqrt=%le\n", sqrt(K_Whole_Val[K_Whole_Ptr[i]]));
		}
		for (i = 0; i < ndof; i++)
		{
			for (j = K_Whole_Ptr[i]; j < K_Whole_Ptr[i + 1]; j++)
			{
				K_Whole_Val[icount] = K_Whole_Val[icount] * diag_scaling[i] * diag_scaling[K_Whole_Col[j]];
				icount++;
			}
			printf("rhs_vec_before[%d]:%le diag_scaling[%d]:%le\n", i, rhs_vec[i], i, diag_scaling[i]);
			rhs_vec[i] = rhs_vec[i] * diag_scaling[i];
		}
	}
	if (flag_operation == 1)
		for (i = 0; i < ndof; i++)
		{
			sol_vec[i] = sol_vec[i] * diag_scaling[i];
		}

	printf("\nqq\n");
}

void CG_Solver(int ndof, int max_itr, double eps, int flag_ini_val)
{
	static double gg[MAX_K_WHOLE_SIZE], dd[MAX_K_WHOLE_SIZE], pp[MAX_K_WHOLE_SIZE];
	static double qqq, ppp, rrr;
	static double alphak, betak;
	int i;
	int itr;
	int ii, istop;

	// Program to solve linear equations by using the CG method
	if (flag_ini_val == 0)
		for (i = 0; i < ndof; i++)
			sol_vec[i] = 0.0;
	
	// Initializing the solution vector if it were not given
	mat_vec_crs(dd, sol_vec, ndof);
	for (i = 0; i < ndof; i++)
	{
		gg[i] = rhs_vec[i] - dd[i];
		printf("rhs_vec[%d]=%f dd[%d]=%f\n", i, rhs_vec[i], i, dd[i]);
		pp[i] = gg[i];
	}

	printf("\nrr");

	for (itr = 0; itr < max_itr; itr++)
	{
		ppp = inner_product(ndof, gg, gg);
		mat_vec_crs(dd, pp, ndof);
		rrr = inner_product(ndof, dd, pp);
		alphak = ppp / rrr;
		for (ii = 0; ii < ndof; ii++)
		{
			sol_vec[ii] += alphak * pp[ii];
			gg[ii] -= alphak * dd[ii];
		}
		qqq = inner_product(ndof, gg, dd);
		betak = qqq / rrr;
		for (ii = 0; ii < ndof; ii++)
			pp[ii] = gg[ii] - betak * pp[ii];
		istop = check_conv_CG(ndof, alphak, pp, eps, itr);
		if (istop == 1)
			break;
	}

	printf("\nss");
}

void Make_M(double *M, int *M_Ptr, int *M_Col, int ndof)
{
	int i, j;
	int ndof_glo = 0;

	// グローバルパッチのdofを求める
	for (i = 0; i < Total_Control_Point_on_mesh[0] * DIMENSION; i++)
	{
		if (Index_Dof[i] != ERROR)
		{
			ndof_glo++;
		}
	}
	printf("ndof		%d\n", ndof);
	printf("ndof_glo	%d\n", ndof_glo);

	int counter = 0;

	// M = [[K^G, 0], [0, K^L]] を作成
	M_Ptr[0] = 0;
	for (i = 0; i < ndof; i++)
	{
		M_Ptr[i + 1] = M_Ptr[i];

		for (j = K_Whole_Ptr[i]; j < K_Whole_Ptr[i + 1]; j++)
		{
			if (i < ndof_glo && K_Whole_Col[j] < ndof_glo)
			{
				M[counter] = K_Whole_Val[j];
				M_Col[counter] = K_Whole_Col[j];
				counter++;
				M_Ptr[i + 1]++;
			}
			else if (i >= ndof_glo)
			{
				M[counter] = K_Whole_Val[j];
				M_Col[counter] = K_Whole_Col[j];
				counter++;
				M_Ptr[i + 1]++;
			}
		}
	}
}

void M_mat_vec_crs(double *M, int *M_Ptr, int *M_Col, double *vec_result, double *vec, const int ndof)
{
	int i, j, icount = 0;

	for (i = 0; i < ndof; i++)
		vec_result[i] = 0;
	for (i = 0; i < ndof; i++)
	{
		for (j = M_Ptr[i]; j < M_Ptr[i + 1]; j++)
		{
			vec_result[i] += M[icount] * vec[M_Col[j]];
			if (i != M_Col[j])
				vec_result[M_Col[j]] += M[icount] * vec[i];
			icount++;
		}
	}
}

int M_check_conv_CG(int ndof, double alphak, double *pp, double eps, double *solution_vec)
{
	double rrr1 = 0.0, rrr2 = 0.0, rrr3;
	int i, istop = 0;
	for (i = 0; i < ndof; i++)
	{
		rrr1 += pp[i] * pp[i];
		rrr2 += solution_vec[i] * solution_vec[i];
	}
	rrr3 = fabs(alphak) * sqrt(rrr1 / rrr2);
	if (rrr3 < eps)
		istop = 1;
	return (istop);
}

void CG(int ndof, double *solution_vec, double *M, int *M_Ptr, int *M_Col, double *right_vec)
{
	int i;

	// CG solver
	static double gg[MAX_K_WHOLE_SIZE], dd[MAX_K_WHOLE_SIZE], pp[MAX_K_WHOLE_SIZE];
	static double qqq, ppp, rrr;
	static double alphak, betak;
	int itr;
	int ii, istop;
	int max_itr = ndof;
	double eps = 1.0e-13;

	for (i = 0; i < ndof; i++)
	{
		solution_vec[i] = 0.0;
	}
	M_mat_vec_crs(M, M_Ptr, M_Col, dd, solution_vec, ndof);
	for (i = 0; i < ndof; i++)
	{
		gg[i] = right_vec[i] - dd[i];
		pp[i] = gg[i];
	}
	for (itr = 0; itr < max_itr; itr++)
	{
		ppp = inner_product(ndof, gg, gg);
		M_mat_vec_crs(M, M_Ptr, M_Col, dd, pp, ndof);
		rrr = inner_product(ndof, dd, pp);
		alphak = ppp / rrr;
		for (ii = 0; ii < ndof; ii++)
		{
			solution_vec[ii] += alphak * pp[ii];
			gg[ii] -= alphak * dd[ii];
		}
		qqq = inner_product(ndof, gg, dd);
		betak = qqq / rrr;
		for (ii = 0; ii < ndof; ii++)
			pp[ii] = gg[ii] - betak * pp[ii];
		istop = M_check_conv_CG(ndof, alphak, pp, eps, solution_vec);
		if (istop == 1)
			break;
	}
	printf("\titr %d\n", itr);
}

int RowCol_to_icount(int row, int col)
{
	for (int j = K_Whole_Ptr[row]; j < K_Whole_Ptr[row + 1]; j++)
	{
		if (K_Whole_Col[j] == col)
		{
			return j;
		}
		else if (K_Whole_Col[j] > col)
		{
			return -1;
		}
	}
	return -1;
}

// 前処理付共役勾配法により[K]{d}={f}を解く
void PCG_Solver(int ndof, int max_itetarion, double eps)
{
	int i, j, k;

	double *r = (double *)malloc(sizeof(double) * ndof);
	double *p = (double *)calloc(ndof, sizeof(double));
	double *y = (double *)malloc(sizeof(double) * ndof);
	double *r2 = (double *)calloc(ndof, sizeof(double));

	// 初期化
	for (i = 0; i < ndof; i++)
		sol_vec[i] = 0.0;

	// 前処理行列作成
	double *M = (double *)malloc(sizeof(double) * MAX_NON_ZERO);
	int *M_Ptr = (int *)malloc(sizeof(int) * MAX_K_WHOLE_SIZE + 1);
	int *M_Col = (int *)malloc(sizeof(int) * MAX_NON_ZERO);
	Make_M(M, M_Ptr, M_Col, ndof);

	// 第0近似解に対する残差の計算
	double *ax = (double *)calloc(ndof, sizeof(double));
	for (i = 0; i < ndof; i++)
	{
		for (j = K_Whole_Ptr[i]; j < K_Whole_Ptr[i + 1]; j++)
		{
			ax[i] += K_Whole_Val[j] * sol_vec[K_Whole_Col[j]];
			if (i != K_Whole_Col[j])
			{
				ax[K_Whole_Col[j]] += K_Whole_Val[j] * sol_vec[i];
			}
		}
	}
	for (i = 0; i < ndof; i++)
	{
		r[i] = rhs_vec[i] - ax[i];
	}
	free(ax);

	// 第0近似解に対する残差の計算
	// for (i = 0; i < ndof; i++)
	// {
	// 	r[i] = rhs_vec[i];
	// }

	// p_0 = (LDL^T)^-1 r_0 の計算 <- CG法で M = [[K^G, 0], [0, K^L]] とし,p_0 = (LDL^T)^-1 r_0 = M^-1 r_0
	CG(ndof, p, M, M_Ptr, M_Col, r);

	// double rr0 = inner_product(ndof, r, p), rr1;
	double rr0;
	double alpha, beta;

	double e = 0.0;
	for (k = 0; k < max_itetarion; k++)
	{
		// rr0 の計算
		rr0 = inner_product(ndof, r, p);

		// y = AP の計算
		for (i = 0; i < ndof; i++)
		{
			double *temp_array_K = (double *)calloc(ndof, sizeof(double));
			for (j = 0; j < ndof; j++)
			{
				int temp1;
				if (i <= j)
				{
					temp1 = RowCol_to_icount(i, j); // temp_array_K[i][j]
				}
				else if (i > j)
				{
					temp1 = RowCol_to_icount(j, i); // temp_array_K[i][j] = temp_array_K[j][i]
				}

				if (temp1 != -1)
				{
					temp_array_K[j] = K_Whole_Val[temp1];
				}
			}
			y[i] = inner_product(ndof, temp_array_K, p);
			free(temp_array_K);
		}

		// alpha = r*r/(P*AP)の計算
		double temp_scaler = inner_product(ndof, p, y);
		alpha = rr0 / temp_scaler;
		// printf("alpha %le\n", alpha);

		// 解x,残差rの更新
		for (i = 0; i < ndof; i++)
		{
			sol_vec[i] += alpha * p[i];
			r[i] -= alpha * y[i];
		}

		// (r*r)_(k+1)の計算
		CG(ndof, r2, M, M_Ptr, M_Col, r);

		// rr1 = inner_product(ndof, r, r2); // 旧
		// rr1 = inner_product(ndof, y, r2); // 新
		// printf("rr1 %le\n", rr1);

		// 収束判定 (||r||<=eps)
		// double rr1 = inner_product(ndof, y, r2);
		// e = sqrt(fabs(rr1));
		// if(e < eps)
		// {
		//     k++;
		//     break;
		// }

		// 収束判定 (CG法と同じ)
		double e1 = 0.0, e2 = 0.0;
		for (i = 0; i < ndof; i++)
		{
			e1 += p[i] * p[i];
			e2 += sol_vec[i] * sol_vec[i];
		}
		e = fabs(alpha) * sqrt(e1 / e2);
		if (e < eps)
		{
			k++;
			break;
		}

		// βの計算とPの更新
		// beta = rr1 / rr0; //旧
		// beta = - rr1 / temp_scaler; // 新
		beta = -inner_product(ndof, y, r2) / temp_scaler;

		for (i = 0; i < ndof; i++)
		{
			// p[i] = r2[i] - beta * p[i];
			p[i] = r2[i] + beta * p[i];
		}
		// printf("beta %le\n", beta);

		// (r*r)_(k+1)を次のステップのために確保しておく
		// rr0 = rr1;

		printf("itr %d\t", k);
		printf("eps %.15e", e);
		// if (rr1 < 0)
		// {
		// 	printf("\t rr1 < 0");
		// }
		printf("\n");
	}

	int max_itr_result = k;
	double eps_result = e;

	printf("\nndof = %d\n", ndof);
	printf("itr_result = %d\n", max_itr_result);
	printf("eps_result = %.15e\n", eps_result);

	free(r), free(p), free(y), free(r2);
	free(M), free(M_Ptr), free(M_Col);
}




// 逆行列を元の行列に代入
double InverseMatrix_2D(double M[2][2])
{
	int i, j;
	double a[2][2];
	double det = M[0][0] * M[1][1] - M[0][1] * M[1][0];

	if (det == 0)
		return ERROR;

	for (i = 0; i < 2; i++)
	{
		for (j = 0; j < 2; j++)
			a[i][j] = M[i][j];
	}
	M[0][0] = a[1][1] / det;
	M[0][1] = a[0][1] * (-1) / det;
	M[1][0] = a[1][0] * (-1) / det;
	M[1][1] = a[0][0] / det;

	return det;
}

double InverseMatrix_3X3(double M[3][3])
{
	int i, j;
	double a[3][3];
	double det = M[0][0] * M[1][1] * M[2][2] + M[1][0] * M[2][1] * M[0][2] + M[2][0] * M[0][1] * M[1][2] - M[0][0] * M[2][1] * M[1][2] - M[2][0] * M[1][1] * M[0][2] - M[1][0] * M[0][1] * M[2][2];

	if (det == 0)
		return ERROR;

	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
			a[i][j] = M[i][j];
	}
	M[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) / det;
	M[0][1] = (a[0][2] * a[2][1] - a[0][1] * a[2][2]) / det;
	M[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) / det;
	M[1][0] = (a[1][2] * a[2][0] - a[1][0] * a[2][2]) / det;
	M[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]) / det;
	M[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) / det;
	M[2][0] = (a[1][0] * a[2][1] - a[1][1] * a[2][0]) / det;
	M[2][1] = (a[0][1] * a[2][0] - a[0][0] * a[2][1]) / det;
	M[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / det;
	// printf("det;%le\n", det);
	return det;
}

// Newton-Raphson法
// from NURBSviewer
double BasisFunc(double *knot_vec, int knot_index, int order, double xi,
				 double *output, double *d_output)
{
	int p, j;
	double sum1 = 0.0;
	double sum2 = 0.0;
	double temp_basis[MAX_ORDER][MAX_ORDER];
	(*output) = 0.0;
	(*d_output) = 0.0;

	if (knot_vec[knot_index] <= xi && knot_vec[knot_index + order + 1] >= xi)
	{
		for (j = 0; j <= order; j++)
		{
			if ((knot_vec[knot_index + j] <= xi) && (xi <= knot_vec[knot_index + j + 1]))
			{
				temp_basis[j][0] = 1.0;
			}
			else
			{
				temp_basis[j][0] = 0.0;
			}
		}

		if (order > 0)
		{
			for (p = 1; p <= order; p++)
			{
				for (j = 0; j <= order - p; j++)
				{
					sum1 = 0.0;
					sum2 = 0.0;
					if ((knot_vec[knot_index + j + p] - knot_vec[knot_index + j]) != 0.0)
					{
						sum1 = (xi - knot_vec[knot_index + j]) / (knot_vec[knot_index + j + p] - knot_vec[knot_index + j]) * temp_basis[j][p - 1];
					}
					if ((knot_vec[knot_index + j + p + 1] - knot_vec[knot_index + j + 1]) != 0.0)
					{
						sum2 = (knot_vec[knot_index + j + p + 1] - xi) / (knot_vec[knot_index + j + p + 1] - knot_vec[knot_index + j + 1]) * temp_basis[j + 1][p - 1];
					}
					temp_basis[j][p] = sum1 + sum2;
				}
			}
			sum1 = 0.0;
			sum2 = 0.0;
			if ((knot_vec[knot_index + order] - knot_vec[knot_index]) != 0.0)
			{
				sum1 = order / (knot_vec[knot_index + order] - knot_vec[knot_index]) * temp_basis[0][order - 1];
			}
			if ((knot_vec[knot_index + order + 1] - knot_vec[knot_index + 1]) != 0.0)
			{
				sum2 = order / (knot_vec[knot_index + order + 1] - knot_vec[knot_index + 1]) * temp_basis[1][order - 1];
			}
		}
		(*output) = temp_basis[0][order];
		(*d_output) = sum1 - sum2;
	}
	return (*output);
}

double rBasisFunc(double *knot_vec, int knot_index,
				  int order, double xi,
				  double *output, double *d_output)
{
	int p, j;
	double sum1 = 0.0;
	double sum2 = 0.0;
	double temp_basis[MAX_ORDER][MAX_ORDER];
	(*output) = 0.0;
	(*d_output) = 0.0;

	if (knot_vec[knot_index] <= xi && xi <= knot_vec[knot_index + order + 1])
	{
		if (knot_index == 0)
		{
			for (j = 0; j <= order; j++)
			{
				if ((knot_vec[j] <= xi) && (xi <= knot_vec[j + 1]))
				{
					temp_basis[j][0] = 1.0;
				}
				else
				{
					temp_basis[j][0] = 0.0;
				}
			}
		}
		else
		{
			for (j = 0; j <= order; j++)
			{
				if ((knot_vec[knot_index + j] < xi) && (xi <= knot_vec[knot_index + j + 1]))
				{
					temp_basis[j][0] = 1.0;
				}
				else
				{
					temp_basis[j][0] = 0.0;
				}
			}
		}

		if (order > 0)
		{
			for (p = 1; p <= order; p++)
			{
				for (j = 0; j <= order - p; j++)
				{
					sum1 = 0.0;
					sum2 = 0.0;
					if ((knot_vec[knot_index + j + p] - knot_vec[knot_index + j]) != 0.0)
					{
						sum1 = (xi - knot_vec[knot_index + j]) / (knot_vec[knot_index + j + p] - knot_vec[knot_index + j]) * temp_basis[j][p - 1];
					}
					if ((knot_vec[knot_index + j + p + 1] - knot_vec[knot_index + j + 1]) != 0.0)
					{
						sum2 = (knot_vec[knot_index + j + p + 1] - xi) / (knot_vec[knot_index + j + p + 1] - knot_vec[knot_index + j + 1]) * temp_basis[j + 1][p - 1];
					}
					temp_basis[j][p] = sum1 + sum2;
				}
			}
			sum1 = 0.0;
			sum2 = 0.0;

			// for (int temp_i = 0; temp_i < No_knot[0][0]; temp_i++)
			// {
			// 	printf("knot_vec[%d] = %f\n", temp_i, knot_vec[temp_i]);
			// }

			if ((knot_vec[knot_index + order] - knot_vec[knot_index]) != 0.0)
			{
				sum1 = order / (knot_vec[knot_index + order] - knot_vec[knot_index]) * temp_basis[0][order - 1];
			}
			if ((knot_vec[knot_index + order + 1] - knot_vec[knot_index + 1]) != 0.0)
			{
				sum2 = order / (knot_vec[knot_index + order + 1] - knot_vec[knot_index + 1]) * temp_basis[1][order - 1];
			}
		}
		(*output) = temp_basis[0][order];
		(*d_output) = sum1 - sum2;
	}
	return (*output);
}

double lBasisFunc(double *knot_vec, int knot_index,
				  int cntl_p_n, int order, double xi,
				  double *output, double *d_output)
{
	int p, j;
	double sum1 = 0.0;
	double sum2 = 0.0;
	double temp_basis[MAX_ORDER][MAX_ORDER];
	(*output) = 0.0;
	(*d_output) = 0.0;

	if (knot_vec[knot_index] <= xi && xi <= knot_vec[knot_index + order + 1])
	{
		if (knot_index == cntl_p_n - 1)
		{
			for (j = 0; j <= order; j++)
			{
				if ((knot_vec[cntl_p_n - 1 + j] <= xi) && (xi <= knot_vec[cntl_p_n + j]))
				{
					temp_basis[j][0] = 1.0;
				}
				else
				{
					temp_basis[j][0] = 0.0;
				}
			}
		}
		else
		{
			for (j = 0; j <= order; j++)
			{
				if ((knot_vec[knot_index + j] <= xi) && (xi < knot_vec[knot_index + j + 1]))
				{
					temp_basis[j][0] = 1.0;
				}
				else
				{
					temp_basis[j][0] = 0.0;
				}
			}
		}

		if (order > 0)
		{
			for (p = 1; p <= order; p++)
			{
				for (j = 0; j <= order - p; j++)
				{
					sum1 = 0.0;
					sum2 = 0.0;
					if ((knot_vec[knot_index + j + p] - knot_vec[knot_index + j]) != 0.0)
					{
						sum1 = (xi - knot_vec[knot_index + j]) / (knot_vec[knot_index + j + p] - knot_vec[knot_index + j]) * temp_basis[j][p - 1];
					}
					if ((knot_vec[knot_index + j + p + 1] - knot_vec[knot_index + j + 1]) != 0.0)
					{
						sum2 = (knot_vec[knot_index + j + p + 1] - xi) / (knot_vec[knot_index + j + p + 1] - knot_vec[knot_index + j + 1]) * temp_basis[j + 1][p - 1];
					}
					temp_basis[j][p] = sum1 + sum2;
				}
			}
			sum1 = 0.0;
			sum2 = 0.0;
			if ((knot_vec[knot_index + order] - knot_vec[knot_index]) != 0.0)
			{
				sum1 = order / (knot_vec[knot_index + order] - knot_vec[knot_index]) * temp_basis[0][order - 1];
			}
			if ((knot_vec[knot_index + order + 1] - knot_vec[knot_index + 1]) != 0.0)
			{
				sum2 = order / (knot_vec[knot_index + order + 1] - knot_vec[knot_index + 1]) * temp_basis[1][order - 1];
			}
		}
		(*output) = temp_basis[0][order];
		(*d_output) = sum1 - sum2;
	}
	return (*output);
}

double NURBS_surface(double *input_knot_vec_xi, double *input_knot_vec_eta,
					 double *cntl_px, double *cntl_py,
					 int cntl_p_n_xi, int cntl_p_n_eta,
					 double *weight, int order_xi, int order_eta,
					 double xi, double eta,
					 double *output_x, double *output_y,
					 double *output_dxi_x, double *output_deta_x,
					 double *output_dxi_y, double *output_deta_y)
{
	int i, j, temp_index;
	double temp1, temp2, temp3;
	double molecule_x, molecule_y;
	double dxi_molecule_x, dxi_molecule_y;
	double deta_molecule_x, deta_molecule_y;
	double denominator, dxi_denominator, deta_denominator;
	double temp_output_xi, temp_output_eta;
	double temp_d_output_xi, temp_d_output_eta;
	molecule_x = 0.0;
	molecule_y = 0.0;
	denominator = 0.0;
	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_denominator = 0.0;
	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_denominator = 0.0;

	int index_min_xi = 0;
	// int index_max_xi = cntl_p_n_xi; //2020_09_12
	int index_max_xi = cntl_p_n_xi - 1; // 2020_09_12
	int index_min_eta = 0;
	// int index_max_eta = cntl_p_n_eta; //2020_09_12
	int index_max_eta = cntl_p_n_eta - 1; // 2020_09_12

	for (i = 0; i < cntl_p_n_xi; i++)
	{
		if (input_knot_vec_xi[i + 1] > xi)
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}
	if (index_min_xi < 0)
		index_min_xi = 0;
	if (index_max_xi > cntl_p_n_xi)
		index_max_xi = cntl_p_n_xi; // 2020_09_12

	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if (input_knot_vec_eta[i + 1] > eta)
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0)
		index_min_eta = 0;
	if (index_max_eta > cntl_p_n_eta)
		index_max_eta = cntl_p_n_eta; // 2020_09_12

	for (i = index_min_xi; i <= index_max_xi; i++)
	{
		BasisFunc(input_knot_vec_xi, i, order_xi, xi,
				  &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j <= index_max_eta; j++)
		{
			BasisFunc(input_knot_vec_eta, j, order_eta, eta,
					  &temp_output_eta, &temp_d_output_eta);
			temp_index = i + j * cntl_p_n_xi;
			temp1 = temp_output_xi * temp_output_eta * weight[temp_index];
			temp2 = temp_d_output_xi * temp_output_eta * weight[temp_index];
			temp3 = temp_output_xi * temp_d_output_eta * weight[temp_index];
			molecule_x += temp1 * cntl_px[temp_index];
			molecule_y += temp1 * cntl_py[temp_index];
			denominator += temp1;
			dxi_molecule_x += temp2 * cntl_px[temp_index];
			dxi_molecule_y += temp2 * cntl_py[temp_index];
			dxi_denominator += temp2;
			deta_molecule_x += temp3 * cntl_px[temp_index];
			deta_molecule_y += temp3 * cntl_py[temp_index];
			deta_denominator += temp3;
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;

	temp1 = denominator * denominator;
	(*output_dxi_x) = (dxi_molecule_x * denominator - molecule_x * dxi_denominator) / temp1;
	(*output_dxi_y) = (dxi_molecule_y * denominator - molecule_y * dxi_denominator) / temp1;
	(*output_deta_x) = (deta_molecule_x * denominator - molecule_x * deta_denominator) / temp1;
	(*output_deta_y) = (deta_molecule_y * denominator - molecule_y * deta_denominator) / temp1;
	return denominator;
}

double rNURBS_surface(double *input_knot_vec_xi, double *input_knot_vec_eta,
					  double *cntl_px, double *cntl_py,
					  int cntl_p_n_xi, int cntl_p_n_eta,
					  double *weight, int order_xi, int order_eta,
					  double xi, double eta,
					  double *output_x, double *output_y,
					  double *output_dxi_x, double *output_deta_x,
					  double *output_dxi_y, double *output_deta_y)
{
	int i, j, temp_index;
	double temp1, temp2, temp3;
	double molecule_x, molecule_y;
	double dxi_molecule_x, dxi_molecule_y;
	double deta_molecule_x, deta_molecule_y;
	double denominator, dxi_denominator, deta_denominator;
	double temp_output_xi, temp_output_eta;
	double temp_d_output_xi, temp_d_output_eta;
	molecule_x = 0.0;
	molecule_y = 0.0;
	denominator = 0.0;
	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_denominator = 0.0;
	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_denominator = 0.0;

	int index_min_xi = 0;
	int index_max_xi = cntl_p_n_xi - 1;
	int index_min_eta = 0;
	int index_max_eta = cntl_p_n_eta - 1;

	for (i = 0; i < cntl_p_n_xi; i++)
	{
		if (input_knot_vec_xi[i + 1] >= xi)
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}
	if (index_min_xi < 0)
		index_min_xi = 0;
	if (index_max_xi > cntl_p_n_xi)
		index_max_xi = cntl_p_n_xi;

	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if (input_knot_vec_eta[i + 1] >= eta)
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0)
		index_min_eta = 0;
	if (index_max_eta > cntl_p_n_eta)
		index_max_eta = cntl_p_n_eta;

	for (i = index_min_xi; i <= index_max_xi; i++)
	{
		rBasisFunc(input_knot_vec_xi, i, order_xi, xi,
				   &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j <= index_max_eta; j++)
		{
			rBasisFunc(input_knot_vec_eta, j, order_eta, eta,
					   &temp_output_eta, &temp_d_output_eta);
			temp_index = i + j * cntl_p_n_xi;
			temp1 = temp_output_xi * temp_output_eta * weight[temp_index];
			temp2 = temp_d_output_xi * temp_output_eta * weight[temp_index];
			temp3 = temp_output_xi * temp_d_output_eta * weight[temp_index];
			molecule_x += temp1 * cntl_px[temp_index];
			molecule_y += temp1 * cntl_py[temp_index];
			denominator += temp1;
			dxi_molecule_x += temp2 * cntl_px[temp_index];
			dxi_molecule_y += temp2 * cntl_py[temp_index];
			dxi_denominator += temp2;
			deta_molecule_x += temp3 * cntl_px[temp_index];
			deta_molecule_y += temp3 * cntl_py[temp_index];
			deta_denominator += temp3;
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;

	temp1 = denominator * denominator;
	(*output_dxi_x) = (dxi_molecule_x * denominator - molecule_x * dxi_denominator) / temp1;
	(*output_dxi_y) = (dxi_molecule_y * denominator - molecule_y * dxi_denominator) / temp1;
	(*output_deta_x) = (deta_molecule_x * denominator - molecule_x * deta_denominator) / temp1;
	(*output_deta_y) = (deta_molecule_y * denominator - molecule_y * deta_denominator) / temp1;
	return denominator;
}

double lNURBS_surface(double *input_knot_vec_xi, double *input_knot_vec_eta,
					  double *cntl_px, double *cntl_py,
					  int cntl_p_n_xi, int cntl_p_n_eta,
					  double *weight, int order_xi, int order_eta,
					  double xi, double eta,
					  double *output_x, double *output_y,
					  double *output_dxi_x, double *output_deta_x,
					  double *output_dxi_y, double *output_deta_y)
{
	int i, j, temp_index;
	double temp1, temp2, temp3;
	double molecule_x, molecule_y;
	double dxi_molecule_x, dxi_molecule_y;
	double deta_molecule_x, deta_molecule_y;
	double denominator, dxi_denominator, deta_denominator;
	double temp_output_xi, temp_output_eta;
	double temp_d_output_xi, temp_d_output_eta;
	molecule_x = 0.0;
	molecule_y = 0.0;
	denominator = 0.0;
	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_denominator = 0.0;
	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_denominator = 0.0;

	int index_min_xi = 0;
	int index_max_xi = cntl_p_n_xi - 1;
	int index_min_eta = 0;
	int index_max_eta = cntl_p_n_eta - 1;

	for (i = 0; i < cntl_p_n_xi; i++)
	{
		if (input_knot_vec_xi[i + 1] >= xi)
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}
	if (index_min_xi < 0)
		index_min_xi = 0;
	if (index_max_xi > cntl_p_n_xi)
		index_max_xi = cntl_p_n_xi;

	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if (input_knot_vec_eta[i + 1] >= eta)
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0)
		index_min_eta = 0;
	if (index_max_eta > cntl_p_n_eta)
		index_max_eta = cntl_p_n_eta;

	for (i = index_min_xi; i <= index_max_xi; i++)
	{
		lBasisFunc(input_knot_vec_xi, i,
				   cntl_p_n_xi, order_xi, xi,
				   &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j <= index_max_eta; j++)
		{
			lBasisFunc(input_knot_vec_eta, j,
					   cntl_p_n_eta, order_eta, eta,
					   &temp_output_eta, &temp_d_output_eta);
			temp_index = i + j * cntl_p_n_xi;
			temp1 = temp_output_xi * temp_output_eta * weight[temp_index];
			temp2 = temp_d_output_xi * temp_output_eta * weight[temp_index];
			temp3 = temp_output_xi * temp_d_output_eta * weight[temp_index];
			molecule_x += temp1 * cntl_px[temp_index];
			molecule_y += temp1 * cntl_py[temp_index];
			denominator += temp1;
			dxi_molecule_x += temp2 * cntl_px[temp_index];
			dxi_molecule_y += temp2 * cntl_py[temp_index];
			dxi_denominator += temp2;
			deta_molecule_x += temp3 * cntl_px[temp_index];
			deta_molecule_y += temp3 * cntl_py[temp_index];
			deta_denominator += temp3;
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;

	temp1 = denominator * denominator;
	(*output_dxi_x) = (dxi_molecule_x * denominator - molecule_x * dxi_denominator) / temp1;
	(*output_dxi_y) = (dxi_molecule_y * denominator - molecule_y * dxi_denominator) / temp1;
	(*output_deta_x) = (deta_molecule_x * denominator - molecule_x * deta_denominator) / temp1;
	(*output_deta_y) = (deta_molecule_y * denominator - molecule_y * deta_denominator) / temp1;
	return denominator;
}

double rlNURBS_surface(double *input_knot_vec_xi, double *input_knot_vec_eta,
					   double *cntl_px, double *cntl_py,
					   int cntl_p_n_xi, int cntl_p_n_eta,
					   double *weight, int order_xi, int order_eta,
					   double xi, double eta,
					   double *output_x, double *output_y,
					   double *output_dxi_x, double *output_deta_x,
					   double *output_dxi_y, double *output_deta_y)
{
	int i, j, temp_index;
	double temp1, temp2, temp3;
	double molecule_x, molecule_y;
	double dxi_molecule_x, dxi_molecule_y;
	double deta_molecule_x, deta_molecule_y;
	double denominator, dxi_denominator, deta_denominator;
	double temp_output_xi, temp_output_eta;
	double temp_d_output_xi, temp_d_output_eta;
	molecule_x = 0.0;
	molecule_y = 0.0;
	denominator = 0.0;
	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_denominator = 0.0;
	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_denominator = 0.0;

	int index_min_xi = 0;
	int index_max_xi = cntl_p_n_xi - 1;
	int index_min_eta = 0;
	int index_max_eta = cntl_p_n_eta - 1;

	for (i = 0; i < cntl_p_n_xi; i++)
	{
		if (input_knot_vec_xi[i + 1] >= xi)
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}
	if (index_min_xi < 0)
		index_min_xi = 0;
	if (index_max_xi > cntl_p_n_xi)
		index_max_xi = cntl_p_n_xi;

	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if (input_knot_vec_eta[i + 1] >= eta)
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0)
		index_min_eta = 0;
	if (index_max_eta > cntl_p_n_eta)
		index_max_eta = cntl_p_n_eta;

	for (i = index_min_xi; i <= index_max_xi; i++)
	{
		rBasisFunc(input_knot_vec_xi, i, order_xi, xi,
				   &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j <= index_max_eta; j++)
		{
			lBasisFunc(input_knot_vec_eta, j,
					   cntl_p_n_eta, order_eta, eta,
					   &temp_output_eta, &temp_d_output_eta);
			temp_index = i + j * cntl_p_n_xi;
			temp1 = temp_output_xi * temp_output_eta * weight[temp_index];
			temp2 = temp_d_output_xi * temp_output_eta * weight[temp_index];
			temp3 = temp_output_xi * temp_d_output_eta * weight[temp_index];
			molecule_x += temp1 * cntl_px[temp_index];
			molecule_y += temp1 * cntl_py[temp_index];
			denominator += temp1;
			dxi_molecule_x += temp2 * cntl_px[temp_index];
			dxi_molecule_y += temp2 * cntl_py[temp_index];
			dxi_denominator += temp2;
			deta_molecule_x += temp3 * cntl_px[temp_index];
			deta_molecule_y += temp3 * cntl_py[temp_index];
			deta_denominator += temp3;
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;

	temp1 = denominator * denominator;
	(*output_dxi_x) = (dxi_molecule_x * denominator - molecule_x * dxi_denominator) / temp1;
	(*output_dxi_y) = (dxi_molecule_y * denominator - molecule_y * dxi_denominator) / temp1;
	(*output_deta_x) = (deta_molecule_x * denominator - molecule_x * deta_denominator) / temp1;
	(*output_deta_y) = (deta_molecule_y * denominator - molecule_y * deta_denominator) / temp1;
	return denominator;
}

double lrNURBS_surface(double *input_knot_vec_xi, double *input_knot_vec_eta,
					   double *cntl_px, double *cntl_py,
					   int cntl_p_n_xi, int cntl_p_n_eta,
					   double *weight, int order_xi, int order_eta,
					   double xi, double eta,
					   double *output_x, double *output_y,
					   double *output_dxi_x, double *output_deta_x,
					   double *output_dxi_y, double *output_deta_y)
{
	int i, j, temp_index;
	double temp1, temp2, temp3;
	double molecule_x, molecule_y;
	double dxi_molecule_x, dxi_molecule_y;
	double deta_molecule_x, deta_molecule_y;
	double denominator, dxi_denominator, deta_denominator;
	double temp_output_xi, temp_output_eta;
	double temp_d_output_xi, temp_d_output_eta;
	molecule_x = 0.0;
	molecule_y = 0.0;
	denominator = 0.0;
	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_denominator = 0.0;
	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_denominator = 0.0;

	int index_min_xi = 0;
	int index_max_xi = cntl_p_n_xi - 1;
	int index_min_eta = 0;
	int index_max_eta = cntl_p_n_eta - 1;

	for (i = 0; i < cntl_p_n_xi; i++)
	{
		if (input_knot_vec_xi[i + 1] >= xi)
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}
	if (index_min_xi < 0)
		index_min_xi = 0;
	if (index_max_xi > cntl_p_n_xi)
		index_max_xi = cntl_p_n_xi;

	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if (input_knot_vec_eta[i + 1] >= eta)
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0)
		index_min_eta = 0;
	if (index_max_eta > cntl_p_n_eta)
		index_max_eta = cntl_p_n_eta;

	for (i = index_min_xi; i <= index_max_xi; i++)
	{
		lBasisFunc(input_knot_vec_xi, i,
				   cntl_p_n_xi, order_xi, xi,
				   &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j <= index_max_eta; j++)
		{
			rBasisFunc(input_knot_vec_eta, j, order_eta, eta,
					   &temp_output_eta, &temp_d_output_eta);
			temp_index = i + j * cntl_p_n_xi;
			temp1 = temp_output_xi * temp_output_eta * weight[temp_index];
			temp2 = temp_d_output_xi * temp_output_eta * weight[temp_index];
			temp3 = temp_output_xi * temp_d_output_eta * weight[temp_index];
			molecule_x += temp1 * cntl_px[temp_index];
			molecule_y += temp1 * cntl_py[temp_index];
			denominator += temp1;
			dxi_molecule_x += temp2 * cntl_px[temp_index];
			dxi_molecule_y += temp2 * cntl_py[temp_index];
			dxi_denominator += temp2;
			deta_molecule_x += temp3 * cntl_px[temp_index];
			deta_molecule_y += temp3 * cntl_py[temp_index];
			deta_denominator += temp3;
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;

	temp1 = denominator * denominator;
	(*output_dxi_x) = (dxi_molecule_x * denominator - molecule_x * dxi_denominator) / temp1;
	(*output_dxi_y) = (dxi_molecule_y * denominator - molecule_y * dxi_denominator) / temp1;
	(*output_deta_x) = (deta_molecule_x * denominator - molecule_x * deta_denominator) / temp1;
	(*output_deta_y) = (deta_molecule_y * denominator - molecule_y * deta_denominator) / temp1;
	return denominator;
}

//算出したローカルパッチ各要素の頂点の物理座標のグローバルパッチでの(xi,eta)算出
int Calc_xi_eta(double px, double py,
				double *input_knot_vec_xi, double *input_knot_vec_eta,
				int cntl_p_n_xi, int cntl_p_n_eta, int order_xi, int order_eta,
				double *output_xi, double *output_eta)
{
	double temp_xi, temp_eta;
	double temp_x, temp_y;
	double temp_matrix[2][2];
	double temp_dxi, temp_deta;
	// double temp_tol_x = DBL_MAX;
	// double temp_tol_y = DBL_MAX;
	double temp_tol_x, temp_tol_y;

	(*output_xi) = 0;
	(*output_eta) = 0;

	int i;
	// int repeat = 1000;
	// double tol = 10e-8;
	// int repeat = 10000;
	int repeat = 100;
	double tol = 10e-14;

	//初期値の設定
	temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;
	// printf("r_temp_xi_eta % 1.8e % 1.8e\n", temp_xi, temp_eta);

	for (i = 0; i < repeat; i++)
	{
		rNURBS_surface(Position_Knots[0][0], Position_Knots[0][1],
					   Control_Coord[0], Control_Coord[1],
					   No_Control_point[0][0], No_Control_point[0][1],
					   Control_Weight, Order[0][0], Order[0][1],
					   temp_xi, temp_eta,
					   &temp_x, &temp_y,
					   &temp_matrix[0][0], &temp_matrix[0][1],
					   &temp_matrix[1][0], &temp_matrix[1][1]);

		temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;

		// printf("xi_0:  % 1.8e\n", temp_xi);
		// printf("eta_0: % 1.8e\n", temp_eta);
		// printf("px: % 1.8e\n",px);
		// printf("temp_x: % 1.8e\n",temp_x);
		// printf("py: % 1.8e\n",py);
		// printf("temp_y: % 1.8e\n",temp_y);
		// printf("temp_tol_x:  % 1.8e\n", temp_tol_x);
		// printf("temp_tol_y:  % 1.8e\n", temp_tol_y);

		//収束した場合////////////////////////////////////////////////////////////////
		// if (temp_tol_x < tol && temp_tol_y < tol)
		if (temp_tol_x + temp_tol_y < tol)
		{
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;

			// printf("r_xi:  % 1.8e\n", temp_xi);
			// printf("r_eta: % 1.8e\n", temp_eta);

			return i;
		}

		InverseMatrix_2D(temp_matrix);

		temp_dxi = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y);
		temp_deta = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_deta;
		if (temp_xi < input_knot_vec_xi[0])
			temp_xi = input_knot_vec_xi[0];
		if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
			temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
		if (temp_eta < input_knot_vec_eta[0])
			temp_eta = input_knot_vec_eta[0];
		if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
			temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];

		// printf("r_xi:  % 1.8e\n", temp_xi);
		// printf("r_eta: % 1.8e\n", temp_eta);
		// printf("i=%d\n",i);

		// double temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
		// printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
	}

	//初期値の設定
	temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;
	// printf("l_temp_xi_eta % 1.8e % 1.8e\n", temp_xi, temp_eta);

	for (i = 0; i < repeat; i++)
	{
		lNURBS_surface(Position_Knots[0][0], Position_Knots[0][1],
					   Control_Coord[0], Control_Coord[1],
					   No_Control_point[0][0], No_Control_point[0][1],
					   Control_Weight, Order[0][0], Order[0][1],
					   temp_xi, temp_eta,
					   &temp_x, &temp_y,
					   &temp_matrix[0][0], &temp_matrix[0][1],
					   &temp_matrix[1][0], &temp_matrix[1][1]);

		temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;

		//収束した場合////////////////////////////////////////////////////////////////
		// if (temp_tol_x < tol && temp_tol_y < tol)
		if (temp_tol_x + temp_tol_y < tol)
		{
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;

			// printf("l_xi:  % 1.8e\n", temp_xi);
			// printf("l_eta: % 1.8e\n", temp_eta);
			// printf("i=%d\n",i);

			return i;
		}

		InverseMatrix_2D(temp_matrix);

		temp_dxi = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y);
		temp_deta = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_deta;
		if (temp_xi < input_knot_vec_xi[0])
			temp_xi = input_knot_vec_xi[0];
		if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
			temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
		if (temp_eta < input_knot_vec_eta[0])
			temp_eta = input_knot_vec_eta[0];
		if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
			temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];

		// double temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
		// printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
	}

	//初期値の設定
	temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;
	// printf("rl_temp_xi_eta % 1.8e % 1.8e\n", temp_xi, temp_eta);

	for (i = 0; i < repeat; i++)
	{
		rlNURBS_surface(Position_Knots[0][0], Position_Knots[0][1],
						Control_Coord[0], Control_Coord[1],
						No_Control_point[0][0], No_Control_point[0][1],
						Control_Weight, Order[0][0], Order[0][1],
						temp_xi, temp_eta,
						&temp_x, &temp_y,
						&temp_matrix[0][0], &temp_matrix[0][1],
						&temp_matrix[1][0], &temp_matrix[1][1]);

		temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;

		//収束した場合////////////////////////////////////////////////////////////////
		// if (temp_tol_x < tol && temp_tol_y < tol)
		if (temp_tol_x + temp_tol_y < tol)
		{
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;

			// printf("rl_xi:  % 1.8e\n", temp_xi);
			// printf("rl_eta: % 1.8e\n", temp_eta);
			// printf("i=%d\n",i);

			return i;
		}

		InverseMatrix_2D(temp_matrix);

		temp_dxi = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y);
		temp_deta = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_deta;
		if (temp_xi < input_knot_vec_xi[0])
			temp_xi = input_knot_vec_xi[0];
		if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
			temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
		if (temp_eta < input_knot_vec_eta[0])
			temp_eta = input_knot_vec_eta[0];
		if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
			temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];

		// double temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
		// printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
	}

	//初期値の設定
	temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;
	// printf("lr_temp_xi_eta % 1.8e % 1.8e\n", temp_xi, temp_eta);

	for (i = 0; i < repeat; i++)
	{
		lrNURBS_surface(Position_Knots[0][0], Position_Knots[0][1],
						Control_Coord[0], Control_Coord[1],
						No_Control_point[0][0], No_Control_point[0][1],
						Control_Weight, Order[0][0], Order[0][1],
						temp_xi, temp_eta,
						&temp_x, &temp_y,
						&temp_matrix[0][0], &temp_matrix[0][1],
						&temp_matrix[1][0], &temp_matrix[1][1]);

		temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;

		//収束した場合////////////////////////////////////////////////////////////////
		// if (temp_tol_x < tol && temp_tol_y < tol)
		if (temp_tol_x + temp_tol_y < tol)
		{
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;

			// printf("lr_xi:  % 1.8e\n", temp_xi);
			// printf("lr_eta: % 1.8e\n", temp_eta);
			// printf("i=%d\n",i);

			return i;
		}

		InverseMatrix_2D(temp_matrix);

		temp_dxi = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y);
		temp_deta = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_deta;
		if (temp_xi < input_knot_vec_xi[0])
			temp_xi = input_knot_vec_xi[0];
		if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
			temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
		if (temp_eta < input_knot_vec_eta[0])
			temp_eta = input_knot_vec_eta[0];
		if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
			temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];

		// double temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
		// printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
	}

	// printf("% 1.8e % 1.8e\n", temp_x, temp_y);
	return 0;
}


// 要素剛性マトリックス
int Jacobian(int El_No, double a[DIMENSION][DIMENSION], double Local_coord[DIMENSION], double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION])
{
	int i, j, k;
	for (i = 0; i < DIMENSION; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			a[i][j] = 0.0;
			for (k = 0; k < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; k++)
			{
				a[i][j] += dShape_func(k, j, Local_coord, El_No) * X[k][i];
			}
		}
	}

	return 0;
}

// Bマトリックスを求める関数
int Make_B_Matrix(int El_No,
				  double B[D_MATRIX_SIZE][MAX_KIEL_SIZE],
				  double Local_coord[DIMENSION],
				  double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION],
				  double *J)
{
	double a[DIMENSION][DIMENSION], b[DIMENSION][MAX_NO_CCpoint_ON_ELEMENT];

	int i, j, k;

	Jacobian(El_No, a, Local_coord, X);

	*J = InverseMatrix_2D(a);
	// printf("B_Matri_J:%le\n",*J);
	if (*J <= 0)
		return -999;

	for (i = 0; i < DIMENSION; i++)
	{
		for (j = 0; j < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; j++)
		{
			b[i][j] = 0.0;
			for (k = 0; k < DIMENSION; k++)
			{
				b[i][j] += a[k][i] * dShape_func(j, k, Local_coord, El_No);
			}
		}
	}

	for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++)
	{
		B[0][2 * i] = b[0][i];
		B[0][2 * i + 1] = 0.0;
		B[1][2 * i] = 0.0;
		B[1][2 * i + 1] = b[1][i];
		B[2][2 * i] = b[1][i];
		B[2][2 * i + 1] = b[0][i];
	}

	return 0;
}


//応力歪マトリックス
int Make_D_Matrix_2D(double D[D_MATRIX_SIZE][D_MATRIX_SIZE], double E, double nu, int DM)
{
	int i, j;

	if (DM == 0) //平面応力状態
	{
		// printf("E:%le nu:%le\n",E,nu);
		double Eone = E / (1.0 - nu * nu);
		double D1[D_MATRIX_SIZE][D_MATRIX_SIZE] = {{Eone, nu * Eone, 0}, {nu * Eone, Eone, 0}, {0, 0, (1 - nu) / 2 * Eone}};

		for (i = 0; i < D_MATRIX_SIZE; i++)
			for (j = 0; j < D_MATRIX_SIZE; j++)
				D[i][j] = D1[i][j];
	}

	else if (DM == 1) //平面ひずみ状態(2Dの場合はこっち)
	{
		// printf("E:%le nu:%le\n",E,nu);
		double Eone = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
		double D1[D_MATRIX_SIZE][D_MATRIX_SIZE] = {{Eone, nu / (1.0 - nu) * Eone, 0}, {nu / (1.0 - nu) * Eone, Eone, 0}, {0, 0, (1 - 2 * nu) / 2 / (1.0 - nu) * Eone}};

		for (i = 0; i < D_MATRIX_SIZE; i++)
			for (j = 0; j < D_MATRIX_SIZE; j++)
				D[i][j] = D1[i][j];
	}

	else
		return ERROR;

	return 0;
}

//ガウスの数値積分法の中身
int BDBJ(double B[D_MATRIX_SIZE][MAX_KIEL_SIZE], double D[D_MATRIX_SIZE][D_MATRIX_SIZE], double J, double K_EL[MAX_KIEL_SIZE][MAX_KIEL_SIZE])
{
	int i, j, k;
	double BD[MAX_KIEL_SIZE][D_MATRIX_SIZE];

	//[B]T[D][B]Jの計算
	for (i = 0; i < KIEL_SIZE; i++)
	{
		for (j = 0; j < D_MATRIX_SIZE; j++)
		{
			BD[i][j] = 0;
			for (k = 0; k < D_MATRIX_SIZE; k++)
			{
				// printf("B[%d][%d]=%le D[%d][%d]=%le\n", k,i,B[k][i],k,j,D[k][j] );
				// printf("B[%d][%d]=%le\n", k,i,B[k][i]);
				BD[i][j] += B[k][i] * D[k][j];
				// printf("BD[%d][%d]=%e\n",i,j,BD[i][j] );
			}
		}
	}
	// for( j = 0; j < D_MATRIX_SIZE; j++ )for( i = 0; i < KIEL_SIZE; i++ )printf("B[%d][%d]=%le\n",j,i,B[j][i] );
	for (i = 0; i < KIEL_SIZE; i++)
	{
		for (j = 0; j < KIEL_SIZE; j++)
		{
			K_EL[i][j] = 0;
			for (k = 0; k < D_MATRIX_SIZE; k++)
			{
				K_EL[i][j] += BD[i][k] * B[k][j];
			}
			K_EL[i][j] *= J;
		}
	}
	return 0;
}

//結合ガウスの数値積分法の中身
int coupled_BDBJ(double B[D_MATRIX_SIZE][MAX_KIEL_SIZE],
				 double D[D_MATRIX_SIZE][D_MATRIX_SIZE],
				 double BG[D_MATRIX_SIZE][MAX_KIEL_SIZE],
				 double J, double K_EL[MAX_KIEL_SIZE][MAX_KIEL_SIZE])
{
	int i, j, k;
	double BD[MAX_KIEL_SIZE][D_MATRIX_SIZE];

	//[B]T[D][B]Jの計算
	for (i = 0; i < KIEL_SIZE; i++)
	{
		for (j = 0; j < D_MATRIX_SIZE; j++)
		{
			BD[i][j] = 0;
			for (k = 0; k < D_MATRIX_SIZE; k++)
			{
				// printf("B[%d][%d]=%le D[%d][%d]=%le\n", k,i,B[k][i],k,j,D[k][j] );
				// printf("B[%d][%d]=%le\n", k,i,B[k][i]);
				BD[i][j] += BG[k][i] * D[k][j];
				// printf("BD[%d][%d]=%e\n",i,j,BD[i][j] );
			}
		}
	}
	// for( j = 0; j < D_MATRIX_SIZE; j++ )for( i = 0; i < KIEL_SIZE; i++ )printf("B[%d][%d]=%le\n",j,i,B[j][i] );
	for (i = 0; i < KIEL_SIZE; i++)
	{
		for (j = 0; j < KIEL_SIZE; j++)
		{
			K_EL[i][j] = 0;
			for (k = 0; k < D_MATRIX_SIZE; k++)
			{
				K_EL[i][j] += BD[i][k] * B[k][j];
			}
			K_EL[i][j] *= J;
		}
	}
	return 0;
}

//要素合成マトリックス
int Make_K_EL(int El_No, double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], double K_EL[MAX_KIEL_SIZE][MAX_KIEL_SIZE], double E, double nu, int DM)
{
	int i, j, k, l;

	double K1[MAX_KIEL_SIZE][MAX_KIEL_SIZE], B[D_MATRIX_SIZE][MAX_KIEL_SIZE];
	double D[D_MATRIX_SIZE][D_MATRIX_SIZE];
	double J = 0.0;
	double J_test = 0.0;

	for (i = 0; i < KIEL_SIZE; i++)
	{
		for (j = 0; j < KIEL_SIZE; j++)
		{
			K_EL[i][j] = 0.0;
		}
	}

	Make_D_Matrix_2D(D, E, nu, DM);

	for (i = 0; i < GP_2D; i++)
	{
		Make_B_Matrix(El_No, B, Gxi[i], X, &J);

		BDBJ(B, D, J, K1);
		J_test += J;
		for (k = 0; k < KIEL_SIZE; k++)
		{
			for (l = 0; l < KIEL_SIZE; l++)
			{
				K_EL[k][l] += w[i] * K1[k][l];
			}
		}
	}
	return 0;
}

//結合要素剛性マトリックス
int Make_coupled_K_EL(int El_No_loc, int El_No_glo,
					  double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION],
					  double XG[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION],
					  double K_EL[MAX_KIEL_SIZE][MAX_KIEL_SIZE],
					  double E, double nu, int DM)
{
	int i, j, jj, k, l;
	int BDBJ_flag;

	double K1[MAX_KIEL_SIZE][MAX_KIEL_SIZE];
	double B[D_MATRIX_SIZE][MAX_KIEL_SIZE], BG[D_MATRIX_SIZE][MAX_KIEL_SIZE];
	double D[D_MATRIX_SIZE][D_MATRIX_SIZE];
	double J = 0.0;
	double J_test = 0.0;
	double G_Gxi[POW_Ng_extended][DIMENSION]; //グローバルパッチ上での親要素内座標xi_bar, eta_bar

	Total_BDBJ_flag = 0;

	for (i = 0; i < KIEL_SIZE; i++)
	{
		for (j = 0; j < KIEL_SIZE; j++)
		{
			K_EL[i][j] = 0.0;
		}
	}

	Make_D_Matrix_2D(D, E, nu, DM);

	for (i = 0; i < GP_2D; i++) //ガウス点のループ(local)
	{
		// printf("gauss point number:%d\n", i);

		////ローカルガウス点がグローバル要素に含まれているかの判定
		//ローカル要素ガウス点の物理座標算出
		double data_result_shape[2] = {0.0};
		double output_xi, output_eta;
		int patch_n = 0;

		for (j = 0; j < No_Control_point_ON_ELEMENT[Element_patch[El_No_loc]]; j++)
		{
			double R_shape_func = Shape_func(j, Gxi[i], El_No_loc);

			for (jj = 0; jj < DIMENSION; jj++)
			{
				data_result_shape[jj] += R_shape_func * X[j][jj];
			}
		}

		//ローカル要素ガウス点のグローバルパッチ上のパラメータ空間座標算出
		for (j = 0; j < Total_Patch_on_mesh[0]; j++) //グローバルメッシュ[0]上
		{
			Calc_xi_eta(data_result_shape[0], data_result_shape[1],
						Position_Knots[j][0], Position_Knots[j][1],
						No_Control_point[j][0], No_Control_point[j][1], Order[j][0], Order[j][1],
						&output_xi, &output_eta);
			// printf("  x: % 1.8e\n", data_result_shape[0]);
			// printf("  y: % 1.8e\n", data_result_shape[1]);
			// printf(" xi: % 1.8e\n", output_xi);
			// printf("eta: % 1.8e\n", output_eta);
			// printf("patch_n: %d\n", j);
			patch_n = j;
		}
		//要素内外判定

		if (output_xi >= Position_Knots[patch_n][0][Order[patch_n][0] + ENC[patch_n][El_No_glo][0]] &&
			output_xi < Position_Knots[patch_n][0][Order[patch_n][0] + ENC[patch_n][El_No_glo][0] + 1] &&
			output_eta >= Position_Knots[patch_n][1][Order[patch_n][1] + ENC[patch_n][El_No_glo][1]] &&
			output_eta < Position_Knots[patch_n][1][Order[patch_n][1] + ENC[patch_n][El_No_glo][1] + 1]) //要素内であるとき
		{
			BDBJ_flag = 1;
			// printf("BDBJ_flag\n");

			//親要素座標の算出
			G_Gxi[i][0] = -1.0 + 2.0 * (output_xi - Position_Knots[patch_n][0][Order[patch_n][0] + ENC[patch_n][El_No_glo][0]]) /
									 (Position_Knots[patch_n][0][Order[patch_n][0] + ENC[patch_n][El_No_glo][0] + 1] - Position_Knots[patch_n][0][Order[patch_n][0] + ENC[patch_n][El_No_glo][0]]);
			G_Gxi[i][1] = -1.0 + 2.0 * (output_eta - Position_Knots[patch_n][1][Order[patch_n][1] + ENC[patch_n][El_No_glo][1]]) /
									 (Position_Knots[patch_n][1][Order[patch_n][1] + ENC[patch_n][El_No_glo][1] + 1] - Position_Knots[patch_n][1][Order[patch_n][1] + ENC[patch_n][El_No_glo][1]]);
			// printf("G_Gxi[][]=\n");
		}
		else //要素外であるとき
		{
			BDBJ_flag = 0;
		}

		// printf("i=%d\n",i );

		////結合要素剛性マトリックス計算
		//要素内であるとき,次を計算
		if (BDBJ_flag)
		{
			Total_BDBJ_flag++;
			Same_BDBJ_flag[i]++;
			// printf("BDBJ_flag\ti=%d\n",i );
			//重なるグローバル要素のBマトリックス
			Make_B_Matrix(El_No_glo, BG, G_Gxi[i], XG, &J);
			//ローカル要素のBマトリックス
			Make_B_Matrix(El_No_loc, B, Gxi[i], X, &J);
			// BGTDBLの計算
			coupled_BDBJ(B, D, BG, J, K1);
			J_test += J;
			for (k = 0; k < KIEL_SIZE; k++)
			{
				for (l = 0; l < KIEL_SIZE; l++)
				{
					K_EL[k][l] += w[i] * K1[k][l];
				}
			} // printf("w[%d]=%f\n",i,w[i]);
		}
	}
	// printf("El_No:%d J_test=%e\n", El_No_loc, J_test);
	// printf("G=%f\n",G );
	/*for ( k = 0; k < KIEL_SIZE; k++) {
		for ( l = 0; l < KIEL_SIZE; l++) {
			printf("K_EL[%d][%d]:%le\n",k,l,K_EL[k][l]);
		}
	}*/

	// if (i == GP_2D - 1)
	// {
	// 	printf("-------------------Total_BDBJ_flag=%d-------------------\n", Total_BDBJ_flag);
	// }

	return 0;
}


// 歪と応力, ひずみエネルギ密度, 変位勾配
void Make_Strain(int Total_Element)
{
	static double U[MAX_KIEL_SIZE];
	double B[D_MATRIX_SIZE][MAX_KIEL_SIZE], X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], J;

	int N, e, i, j;
	// printf("Strain\n");

	Make_gauss_array(0);

	for (e = 0; e < Total_Element; e++)
	{
		// printf("\nElementNo.:%d\n",e);
		for (N = 0; N < GP_2D; N++)
			for (i = 0; i < N_STRAIN; i++)
				Strain[e][N][i] = 0.0;
		// Bマトリックスと各要素の変位を取得
		// printf("El_No:%d,No_Control_point_ON_ELEMENT[%d]:%d\n", El_No,Element_patch[El_No],No_Control_point_ON_ELEMENT[Element_patch[El_No]]);
		for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[e]]; i++)
		{
			for (j = 0; j < DIMENSION; j++)
			{
				U[i * DIMENSION + j] = Displacement[Controlpoint_of_Element[e * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + j];
				X[i][j] = Node_Coordinate[Controlpoint_of_Element[e * MAX_NO_CCpoint_ON_ELEMENT + i] * (DIMENSION + 1) + j];
			}
		}
		//歪
		for (N = 0; N < GP_2D; N++)
		{
			// printf("N:%d\n",N);
			Make_B_Matrix(e, B, Gxi[N], X, &J);
			for (i = 0; i < D_MATRIX_SIZE; i++)
				for (j = 0; j < KIEL_SIZE; j++)
				{
					Strain[e][N][i] += B[i][j] * U[j];
					// printf("B[%d][%d]_in_strain:%le * ",i,j,B[i][j]);
					// if(e==1){
					// printf("U[%d]=%le = %le\n",j,U[j],B[i][j]*U[j]);
					// }
				}
		}
	}
}

//応力
void Make_Stress_2D(double E, double nu, int Total_Element, int DM)
{

	static double D[D_MATRIX_SIZE][D_MATRIX_SIZE];
	int e, i, j, k;
	Make_gauss_array(0);
	Make_D_Matrix_2D(D, E, nu, DM);

	for (e = 0; e < Total_Element; e++)
	{
		for (k = 0; k < GP_2D; k++)
			for (i = 0; i < N_STRESS; i++)
				Stress[e][k][i] = 0.0;
		for (k = 0; k < GP_2D; k++)
			for (i = 0; i < D_MATRIX_SIZE; i++)
				for (j = 0; j < D_MATRIX_SIZE; j++)
					Stress[e][k][i] += D[i][j] * Strain[e][k][j];
	}
}


void Make_ReactionForce(int Total_Control_Point)
{
	int e, i, j, k, l, re;
	double B[D_MATRIX_SIZE][MAX_KIEL_SIZE], X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], J;

	Make_gauss_array(0);

	for (i = 0; i < Total_Control_Point * DIMENSION; i++)
		ReactionForce[i] = 0.0;
	// printf("ReactionForce\n");
	for (re = 0; re < real_Total_Element_to_mesh[Total_mesh]; re++)
	{
		e = real_element[re];
		// printf("%d,No_Control_point_ON_ELEMENT[%d]:%d\n",e,Element_patch[e], No_Control_point_ON_ELEMENT[Element_patch[e]]);
		// Bマトリックスを取得
		for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[e]]; i++)
		{
			for (j = 0; j < DIMENSION; j++)
			{
				X[i][j] = Node_Coordinate[Controlpoint_of_Element[e * MAX_NO_CCpoint_ON_ELEMENT + i] * (DIMENSION + 1) + j];
			}
		}
		for (k = 0; k < GP_2D; k++)
		{
			Make_B_Matrix(e, B, Gxi[k], X, &J);
			for (j = 0; j < D_MATRIX_SIZE; j++)
				for (l = 0; l < DIMENSION; l++)
					for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[e]]; i++)
						ReactionForce[Controlpoint_of_Element[e * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + l] += B[j][i * DIMENSION + l] * Stress[e][k][j] * w[k] * J;
			// printf("J:%le\n", J);
		}
	}
}


void Make_Parameter_z(int Total_Element, double E, double nu, int DM)
{
	int e, k;
	Make_gauss_array(0);

	if (DM == 0)
	{
		// Make_strain_z
		for (e = 0; e < Total_Element; e++)
			for (k = 0; k < GP_2D; k++)
				Strain[e][k][3] = 0.0;

		for (e = 0; e < Total_Element; e++)
			for (k = 0; k < GP_2D; k++)
				Strain[e][k][3] = -1.0 * nu / E * (Stress[e][k][0] + Stress[e][k][1]);
	}

	if (DM == 1)
	{
		// Make_stree_z
		for (e = 0; e < Total_Element; e++)
			for (k = 0; k < GP_2D; k++)
				Stress[e][k][3] = 0.0;

		for (e = 0; e < Total_Element; e++)
			for (k = 0; k < GP_2D; k++)
				Stress[e][k][3] = E * nu / (1.0 + nu) / (1 - 2.0 * nu) * (Strain[e][k][0] + Strain[e][k][1]);
	}
}


void Make_Parameter_z_overlay(int Total_Element, double E, double nu, int DM)
{
	int e, k;

	Make_gauss_array(0);

	if (DM == 0)
	{
		// Make_strain_z
		for (e = 0; e < Total_Element; e++)
			for (k = 0; k < GP_2D; k++)
				Strain_overlay[e][k][3] = 0.0;

		for (e = 0; e < Total_Element; e++)
			for (k = 0; k < GP_2D; k++)
				Strain_overlay[e][k][3] = -1.0 * nu / E * (Stress_overlay[e][k][0] + Stress_overlay[e][k][1]);
	}

	if (DM == 1)
	{
		// Make_stree_z
		for (e = 0; e < Total_Element; e++)
			for (k = 0; k < GP_2D; k++)
				Stress_overlay[e][k][3] = 0.0;

		for (e = 0; e < Total_Element; e++)
			for (k = 0; k < GP_2D; k++)
				Stress_overlay[e][k][3] = E * nu / (1.0 + nu) / (1 - 2.0 * nu) * (Strain_overlay[e][k][0] + Strain_overlay[e][k][1]);
	}
}



//重ね合わせた結果の出力(NURBS_input_for_s-IGA)
void GetLocData()
{
	double temp;

	//必要なのはローカルのパッチ数とコントロールポイント数
	printf("Start Get Local Data\n\n");
	fp = fopen("input_local.txt", "r");

	fscanf(fp, "%lf%lf", &temp, &temp);
	fscanf(fp, "\n");
	fscanf(fp, "%d", &n_patch_loc);
	fscanf(fp, "\n");
	fscanf(fp, "%d", &loc_cntl_p_n);
	fscanf(fp, "\n");
	fclose(fp);
	printf("patches(in local):% d\n", n_patch_loc);
	printf("control points(in local):% d\n", loc_cntl_p_n);
	printf("\nFinish Get Local Data\n\n");
}

void ReadFile()
{
	int i, j;
	double temp1, temp2, temp3;
	int temp_int;
	printf("Start Reading input\n\n");
	fp = fopen("input_for_NURBS.txt", "r");

	fscanf(fp, "%lf%lf", &E, &nu);
	fscanf(fp, "\n");
	printf("E nu: % 1.4e % 1.4e\n", E, nu);

	fscanf(fp, "%d", &patch_n);
	fscanf(fp, "\n");
	printf("patches: %d \n", patch_n);
	if (patch_n > MAX_PATCHES)
	{
		printf("Error!!\n");
		printf("Too many patches!\n"
			   "Maximum of patches is %d (Now %d)\n"
			   "\n",
			   MAX_PATCHES, patch_n);
		exit(1);
	}

	fscanf(fp, "%d", &cntl_p_n);
	fscanf(fp, "\n");
	printf("total control points:%d \n", cntl_p_n);
	if (cntl_p_n > MAX_CNRL_P)
	{
		printf("Error!!\n");
		printf("Too many control points!\n"
			   "Maximum of control points is %d (Now %d)\n"
			   "\n",
			   MAX_CNRL_P, cntl_p_n);
		exit(1);
	}

	for (i = 0; i < patch_n; i++)
	{
		fscanf(fp, "%d%d", &order_xi[i], &order_eta[i]);
		fscanf(fp, "\n");
		printf("order %d: %d %d\n", i, order_xi[i], order_eta[i]);
		if (order_xi[i] > MAX_ORDER)
		{
			printf("Error!!\n");
			printf("Order too big at xi!\n"
				   "Maximum of order is %d (Now %d at patch %d)\n"
				   "\n",
				   MAX_ORDER, order_xi[i], i);
			exit(1);
		}
		if (order_eta[i] > MAX_ORDER)
		{
			printf("Error!!\n");
			printf("Order too big at eta!\n"
				   "Maximum of order is %d (Now %d at patch %d)\n"
				   "\n",
				   MAX_ORDER, order_eta[i], i);
			exit(1);
		}
	}

	for (i = 0; i < patch_n; i++)
	{
		fscanf(fp, "%d%d", &knot_n_xi[i], &knot_n_eta[i]);
		fscanf(fp, "\n");
		printf("knots %d: %d %d\n", i, knot_n_xi[i], knot_n_eta[i]);
		if (knot_n_xi[i] > MAX_KNOTS)
		{
			printf("Error!!\n");
			printf("Knot vector too long at xi!\n"
				   "Maximum of knot vector is %d (Now %d at patch %d)\n"
				   "\n",
				   MAX_KNOTS, knot_n_xi[i], i);
			exit(1);
		}
		if (knot_n_eta[i] > MAX_KNOTS)
		{
			printf("Error!!\n");
			printf("Knot vector too long at eta!\n"
				   "Maximum of knot vector is %d (Now %d at patch %d)\n"
				   "\n",
				   MAX_KNOTS, knot_n_eta[i], i);
			exit(1);
		}
	}

	for (i = 0; i < patch_n; i++)
	{
		fscanf(fp, "%d%d", &cntl_p_n_xi[i], &cntl_p_n_eta[i]);
		printf("control points %d: %d %d\n",
			   i, cntl_p_n_xi[i], cntl_p_n_eta[i]);
		fscanf(fp, "\n");
	}
	printf("\n");

	for (i = 0; i < patch_n; i++)
	{
		for (j = 0; j < cntl_p_n_xi[i] * cntl_p_n_eta[i]; j++)
		{
			fscanf(fp, "%d", &temp_index[i][j]);
			printf("%d ", temp_index[i][j]);
		}
		fscanf(fp, "\n");
		printf("\n");
	}
	printf("\n");

	fscanf(fp, "%lf%lf%lf", &temp1, &temp2, &temp3);
	fscanf(fp, "\n");

	for (i = 0; i < patch_n; i++)
	{
		for (j = 0; j < knot_n_xi[i]; j++)
		{
			fscanf(fp, "%le", &knot_vec_xi[i][j]);
			printf("%f\t", knot_vec_xi[i][j]);
		}
		printf("\n");
		for (j = 0; j < knot_n_eta[i]; j++)
		{
			fscanf(fp, "%le", &knot_vec_eta[i][j]);
			printf("%f\t", knot_vec_eta[i][j]);
		}
		printf("\n");
	}
	printf("\n");

	for (i = 0; i < cntl_p_n; i++)
	{
		fscanf(fp, "%d%le%le%le",
			   &temp_int,
			   &temp_cntl_px[i], &temp_cntl_py[i], &temp_weight[i]);
		printf("%d\t%f\t%f\t%f\n",
			   temp_int,
			   temp_cntl_px[i], temp_cntl_py[i], temp_weight[i]);
	}
	printf("\n");

	for (i = 0; i < patch_n; i++)
	{
		for (j = 0; j < cntl_p_n_xi[i] * cntl_p_n_eta[i]; j++)
		{
			cntl_px[i][j] = temp_cntl_px[temp_index[i][j]];
			cntl_py[i][j] = temp_cntl_py[temp_index[i][j]];
			weight[i][j] = temp_weight[temp_index[i][j]];
			printf("%d\t%f\t%f\t%f\n",
				   temp_index[i][j],
				   cntl_px[i][j], cntl_py[i][j], weight[i][j]);
		}
		printf("\n");
	}
	fclose(fp);
	printf("End Reading input\n\n");

	if (fields_flag)
	{
		printf("Start Reading displacement\n\n");
		fp = fopen("Displacement.dat", "r");
		char buff[256];

		fscanf(fp, "%s", buff);
		fscanf(fp, "%s", buff);

		for (i = 0; i < cntl_p_n; i++)
		{
			fscanf(fp, "%d:%le%le",
				   &temp_int, &temp_disp_x[i], &temp_disp_y[i]);
			printf("%d\t%1.6e\t%1.6e\n",
				   temp_int, temp_disp_x[i], temp_disp_y[i]);
		}
		printf("\n");

		for (i = 0; i < patch_n; i++)
		{
			for (j = 0; j < cntl_p_n_xi[i] * cntl_p_n_eta[i]; j++)
			{
				disp_cntl_px[i][j] = temp_disp_x[temp_index[i][j]];
				disp_cntl_py[i][j] = temp_disp_y[temp_index[i][j]];
				printf("%d\t%f\t%f\t%f\n",
					   temp_index[i][j], cntl_px[i][j], cntl_py[i][j], weight[i][j]);
			}
			printf("\n");
		}
		fclose(fp); // ファイルを閉じる
		printf("End Reading displpacement\n\n");
	}

	fp = fopen("Displacement_loc.dat", "w");
	glo_cntl_p_n = cntl_p_n - loc_cntl_p_n;
	fprintf(fp, "label=Displacement\n"
				"num_items=%d\n\n",
			loc_cntl_p_n);
	for (i = 0; i < loc_cntl_p_n; i++)
	{
		// fprintf(fp, "%d:	%le %le \n", i, temp_disp_x[i + glo_cntl_p_n], temp_disp_y[i + glo_cntl_p_n]);
		fprintf(fp, "%d:	%.16e %.16e \n", i, temp_disp_x[i + glo_cntl_p_n], temp_disp_y[i + glo_cntl_p_n]);
	}
	fclose(fp);
}

int CalcXiEtaByNR(double px, double py,
				  double *input_knot_vec_xi, double *input_knot_vec_eta,
				  double *cntl_px, double *cntl_py,
				  double *disp_cntl_px, double *disp_cntl_py,
				  int cntl_p_n_xi, int cntl_p_n_eta,
				  double *weight, int order_xi, int order_eta,
				  double *output_xi, double *output_eta,
				  double *disp_x_glo, double *disp_y_glo,
				  double *strain_xx_glo, double *strain_yy_glo, double *strain_xy_glo)
{
	double temp_xi, temp_eta;
	double temp_x, temp_y;
	double temp_matrix[2][2];
	double temp_dxi, temp_deta;
	double temp_tol_x = DBL_MAX;
	double temp_tol_y = DBL_MAX;

	(*output_xi) = 0;
	(*output_eta) = 0;

	int i;
	// double tol = 10e-8;
	// int repeat = 10000;
	int repeat = 100;
	double tol = 10e-14;

	temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;
	// printf("% 1.8e % 1.8e\n", temp_xi, temp_eta);

	for (i = 0; i < repeat; i++)
	{
		rNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
					   cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
					   weight, order_xi, order_eta,
					   temp_xi, temp_eta,
					   &temp_x, &temp_y,
					   &temp_matrix[0][0], &temp_matrix[0][1],
					   &temp_matrix[1][0], &temp_matrix[1][1]);

		temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;

		//収束した場合////////////////////////////////////////////////////////////////
		// if (temp_tol_x < tol && temp_tol_y < tol) {
		if (temp_tol_x + temp_tol_y < tol)
		{
			// printf("rNURBS\n");
			// printf("repeat = %d\n", i);
			if (temp_xi == input_knot_vec_xi[0] || temp_eta == input_knot_vec_eta[0])
			{
				break;
			}
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;

			int knot_n_xi = cntl_p_n_xi + order_xi + 1;
			int knot_n_eta = cntl_p_n_eta + order_eta + 1;
			double dtilda_xi = 0.0;
			double dtilda_eta = 0.0;
			double disp_x, disp_y;
			double dxi_x, deta_x, dxi_y, deta_y;
			double dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y;
			double temp_matrix2[2][2];
			double temp;
			double strain_xx, strain_yy, strain_xy;
			// double stress_xx, stress_yy, stress_xy;

			for (i = 0; i < knot_n_xi; i++)
			{
				if (input_knot_vec_xi[i] < temp_xi && temp_xi <= input_knot_vec_xi[i + 1])
				{
					dtilda_xi = (input_knot_vec_xi[i + 1] - input_knot_vec_xi[i]) / 2.0;
					// printf("xi%f\n", dtilda_xi);
					break;
				}
			}
			for (i = 0; i < knot_n_eta; i++)
			{
				if (input_knot_vec_eta[i] < temp_eta && temp_eta <= input_knot_vec_eta[i + 1])
				{
					dtilda_eta = (input_knot_vec_eta[i + 1] - input_knot_vec_eta[i]) / 2.0;
					// printf("eta%f\n", dtilda_eta);
					break;
				}
			}

			rNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
						   cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
						   weight, order_xi, order_eta,
						   temp_xi, temp_eta,
						   &temp, &temp,
						   &dxi_x, &deta_x,
						   &dxi_y, &deta_y);
			// printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
			// 	   dxi_x, deta_x, dxi_y, deta_y);

			rNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
						   disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
						   weight, order_xi, order_eta,
						   temp_xi, temp_eta,
						   &disp_x, &disp_y,
						   &dxi_disp_x, &deta_disp_x,
						   &dxi_disp_y, &deta_disp_y);
			// printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
			// 	   dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y);

			temp_matrix2[0][0] = dxi_x * dtilda_xi;
			temp_matrix2[0][1] = dxi_y * dtilda_xi;
			temp_matrix2[1][0] = deta_x * dtilda_eta;
			temp_matrix2[1][1] = deta_y * dtilda_eta;

			InverseMatrix_2D(temp_matrix2);

			strain_xx = temp_matrix2[0][0] * dtilda_xi * dxi_disp_x + temp_matrix2[0][1] * dtilda_eta * deta_disp_x;
			strain_yy = temp_matrix2[1][0] * dtilda_xi * dxi_disp_y + temp_matrix2[1][1] * dtilda_eta * deta_disp_y;
			strain_xy = temp_matrix2[1][0] * dtilda_xi * dxi_disp_x + temp_matrix2[1][1] * dtilda_eta * deta_disp_x + temp_matrix2[0][0] * dtilda_xi * dxi_disp_y + temp_matrix2[0][1] * dtilda_eta * deta_disp_y;

			// 応力の計算を行わないのでコメントアウト
			// double D_matrix[3][3] = {{0.0}};
			// int DM = 1;
			// if (DM == 0) { //平面応力状態
			// 	temp = E * (1.0 - nu * nu);
			// 	D_matrix[0][0] = temp;
			// 	D_matrix[0][1] = nu * temp;
			// 	D_matrix[1][0] = nu * temp;
			// 	D_matrix[1][1] = temp;
			// 	D_matrix[2][2] = (1.0 - nu) / 2.0 * temp;
			// } else if (DM == 1) { //平面ひずみ状態(2Dの場合はこっち)
			// 	temp = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
			// 	D_matrix[0][0] = temp;
			// 	D_matrix[0][1] = nu / (1.0 - nu) * temp;
			// 	D_matrix[1][0] = nu / (1.0 - nu) * temp;
			// 	D_matrix[1][1] = temp;
			// 	D_matrix[2][2] = (1.0 - 2.0 * nu) / 2.0 / (1.0 - nu) * temp;
			// }

			// strainしか使わないのでコメントアウト
			// stress_xx = D_matrix[0][0] * strain_xx + D_matrix[0][1] * strain_yy;
			// stress_yy = D_matrix[1][0] * strain_xx + D_matrix[1][1] * strain_yy;
			// stress_xy = D_matrix[2][2] * strain_xy;

			temp = sqrt(disp_x * disp_x + disp_y * disp_y);
			// printf("x:   % 1.8e\n", px);
			// printf("y:   % 1.8e\n", py);
			// printf("xi:  % 1.8e\n", temp_xi);
			// printf("eta: % 1.8e\n", temp_eta);
			// printf("Displacement x: % 1.8e\n", disp_x);
			// printf("Displacement y: % 1.8e\n", disp_y);
			// printf("Displacement  : % 1.8e\n", temp);
			// printf("Strain xx: % 1.8e\n", strain_xx);
			// printf("Strain yy: % 1.8e\n", strain_yy);
			// printf("Strain xy: % 1.8e\n", strain_xy);
			// printf("Stress xx: % 1.8e\n", stress_xx);
			// printf("Stress yy: % 1.8e\n", stress_yy);
			// printf("Stress xy: % 1.8e\n", stress_xy);

			(*disp_x_glo) = disp_x;
			(*disp_y_glo) = disp_y;

			(*strain_xx_glo) = strain_xx;
			(*strain_yy_glo) = strain_yy;
			(*strain_xy_glo) = strain_xy;

			return i;
		}

		InverseMatrix_2D(temp_matrix);

		temp_dxi = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y);
		temp_deta = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_deta;
		if (temp_xi < input_knot_vec_xi[0])
			temp_xi = input_knot_vec_xi[0];
		if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
			temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
		if (temp_eta < input_knot_vec_eta[0])
			temp_eta = input_knot_vec_eta[0];
		if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
			temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];

		// temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
		// printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
	}
	temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;
	for (i = 0; i < repeat; i++)
	{
		lNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
					   cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
					   weight, order_xi, order_eta,
					   temp_xi, temp_eta,
					   &temp_x, &temp_y,
					   &temp_matrix[0][0], &temp_matrix[0][1],
					   &temp_matrix[1][0], &temp_matrix[1][1]);

		temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;

		//収束した場合////////////////////////////////////////////////////////////////
		// if (temp_tol_x < tol && temp_tol_y < tol) {
		if (temp_tol_x + temp_tol_y < tol)
		{
			// printf("lNURBS\n");
			// printf("repeat = %d\n", i);
			if (temp_xi == input_knot_vec_xi[cntl_p_n_xi + order_xi] || temp_eta == input_knot_vec_eta[cntl_p_n_eta + order_eta])
			{
				break;
			}
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;

			int knot_n_xi = cntl_p_n_xi + order_xi + 1;
			int knot_n_eta = cntl_p_n_eta + order_eta + 1;
			double dtilda_xi = 0.0;
			double dtilda_eta = 0.0;
			double disp_x, disp_y;
			double dxi_x, deta_x, dxi_y, deta_y;
			double dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y;
			double temp_matrix2[2][2];
			double temp;
			double strain_xx, strain_yy, strain_xy;
			// double stress_xx, stress_yy, stress_xy;

			for (i = 0; i < knot_n_xi; i++)
			{
				if (input_knot_vec_xi[i] <= temp_xi && temp_xi < input_knot_vec_xi[i + 1])
				{
					dtilda_xi = (input_knot_vec_xi[i + 1] - input_knot_vec_xi[i]) / 2.0;
					// printf("%f\n", dtilda_xi);
					break;
				}
			}
			for (i = 0; i < knot_n_eta; i++)
			{
				if (input_knot_vec_eta[i] <= temp_eta && temp_eta < input_knot_vec_eta[i + 1])
				{
					dtilda_eta = (input_knot_vec_eta[i + 1] - input_knot_vec_eta[i]) / 2.0;
					// printf("%f\n", dtilda_eta);
					break;
				}
			}

			lNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
						   cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
						   weight, order_xi, order_eta,
						   temp_xi, temp_eta,
						   &temp, &temp,
						   &dxi_x, &deta_x,
						   &dxi_y, &deta_y);
			// printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
			// 	   dxi_x, deta_x, dxi_y, deta_y);

			lNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
						   disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
						   weight, order_xi, order_eta,
						   temp_xi, temp_eta,
						   &disp_x, &disp_y,
						   &dxi_disp_x, &deta_disp_x,
						   &dxi_disp_y, &deta_disp_y);
			// printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
			// 	   dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y);

			temp_matrix2[0][0] = dxi_x * dtilda_xi;
			temp_matrix2[0][1] = dxi_y * dtilda_xi;
			temp_matrix2[1][0] = deta_x * dtilda_eta;
			temp_matrix2[1][1] = deta_y * dtilda_eta;

			InverseMatrix_2D(temp_matrix2);

			strain_xx = temp_matrix2[0][0] * dtilda_xi * dxi_disp_x + temp_matrix2[0][1] * dtilda_eta * deta_disp_x;
			strain_yy = temp_matrix2[1][0] * dtilda_xi * dxi_disp_y + temp_matrix2[1][1] * dtilda_eta * deta_disp_y;
			strain_xy = temp_matrix2[1][0] * dtilda_xi * dxi_disp_x + temp_matrix2[1][1] * dtilda_eta * deta_disp_x + temp_matrix2[0][0] * dtilda_xi * dxi_disp_y + temp_matrix2[0][1] * dtilda_eta * deta_disp_y;

			// 応力の計算を行わないのでコメントアウト
			// double D_matrix[3][3] = {{0.0}};
			// int DM = 1;
			// if (DM == 0) { //平面応力状態
			// 	temp = E * (1.0 - nu * nu);
			// 	D_matrix[0][0] = temp;
			// 	D_matrix[0][1] = nu * temp;
			// 	D_matrix[1][0] = nu * temp;
			// 	D_matrix[1][1] = temp;
			// 	D_matrix[2][2] = (1.0 - nu) / 2.0 * temp;
			// } else if (DM == 1) { //平面ひずみ状態(2Dの場合はこっち)
			// 	temp = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
			// 	D_matrix[0][0] = temp;
			// 	D_matrix[0][1] = nu / (1.0 - nu) * temp;
			// 	D_matrix[1][0] = nu / (1.0 - nu) * temp;
			// 	D_matrix[1][1] = temp;
			// 	D_matrix[2][2] = (1.0 - 2.0 * nu) / 2.0 / (1.0 - nu) * temp;
			// }

			// strainしか使わないのでコメントアウト
			// stress_xx = D_matrix[0][0] * strain_xx + D_matrix[0][1] * strain_yy;
			// stress_yy = D_matrix[1][0] * strain_xx + D_matrix[1][1] * strain_yy;
			// stress_xy = D_matrix[2][2] * strain_xy;

			temp = sqrt(disp_x * disp_x + disp_y * disp_y);
			// printf("x:   % 1.8e\n", px);
			// printf("y:   % 1.8e\n", py);
			// printf("xi:  % 1.8e\n", temp_xi);
			// printf("eta: % 1.8e\n", temp_eta);
			// printf("Displacement x: % 1.8e\n", disp_x);
			// printf("Displacement y: % 1.8e\n", disp_y);
			// printf("Displacement  : % 1.8e\n", temp);
			// printf("Strain xx: % 1.8e\n", strain_xx);
			// printf("Strain yy: % 1.8e\n", strain_yy);
			// printf("Strain xy: % 1.8e\n", strain_xy);
			// printf("Stress xx: % 1.8e\n", stress_xx);
			// printf("Stress yy: % 1.8e\n", stress_yy);
			// printf("Stress xy: % 1.8e\n", stress_xy);

			(*disp_x_glo) = disp_x;
			(*disp_y_glo) = disp_y;

			(*strain_xx_glo) = strain_xx;
			(*strain_yy_glo) = strain_yy;
			(*strain_xy_glo) = strain_xy;

			return i;
		}

		InverseMatrix_2D(temp_matrix);

		temp_dxi = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y);
		temp_deta = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_deta;
		if (temp_xi < input_knot_vec_xi[0])
			temp_xi = input_knot_vec_xi[0];
		if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
			temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
		if (temp_eta < input_knot_vec_eta[0])
			temp_eta = input_knot_vec_eta[0];
		if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
			temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];

		// temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
		// printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
	}
	temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;
	for (i = 0; i < repeat; i++)
	{
		rlNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
						cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
						weight, order_xi, order_eta,
						temp_xi, temp_eta,
						&temp_x, &temp_y,
						&temp_matrix[0][0], &temp_matrix[0][1],
						&temp_matrix[1][0], &temp_matrix[1][1]);

		temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;

		//収束した場合////////////////////////////////////////////////////////////////
		// if (temp_tol_x < tol && temp_tol_y < tol) {
		if (temp_tol_x + temp_tol_y < tol)
		{
			// printf("rlNURBS\n");
			if (temp_xi == input_knot_vec_xi[0])
			{
				break;
			}
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;

			int knot_n_xi = cntl_p_n_xi + order_xi + 1;
			int knot_n_eta = cntl_p_n_eta + order_eta + 1;
			double dtilda_xi = 0.0;
			double dtilda_eta = 0.0;
			double disp_x, disp_y;
			double dxi_x, deta_x, dxi_y, deta_y;
			double dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y;
			double temp_matrix2[2][2];
			double temp;
			double strain_xx, strain_yy, strain_xy;
			// double stress_xx, stress_yy, stress_xy;

			for (i = 0; i < knot_n_xi; i++)
			{
				if (input_knot_vec_xi[i] < temp_xi && temp_xi <= input_knot_vec_xi[i + 1])
				{
					dtilda_xi = (input_knot_vec_xi[i + 1] - input_knot_vec_xi[i]) / 2.0;
					// printf("%f\n", dtilda_xi);
					break;
				}
			}
			for (i = 0; i < knot_n_eta; i++)
			{
				if (input_knot_vec_eta[i] <= temp_eta && temp_eta < input_knot_vec_eta[i + 1])
				{
					dtilda_eta = (input_knot_vec_eta[i + 1] - input_knot_vec_eta[i]) / 2.0;
					// printf("%f\n", dtilda_eta);
					break;
				}
			}

			rlNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
							cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
							weight, order_xi, order_eta,
							temp_xi, temp_eta,
							&temp, &temp,
							&dxi_x, &deta_x,
							&dxi_y, &deta_y);
			// printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
			//	   dxi_x, deta_x, dxi_y, deta_y);

			rlNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
							disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
							weight, order_xi, order_eta,
							temp_xi, temp_eta,
							&disp_x, &disp_y,
							&dxi_disp_x, &deta_disp_x,
							&dxi_disp_y, &deta_disp_y);
			// printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
			//	   dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y);

			temp_matrix2[0][0] = dxi_x * dtilda_xi;
			temp_matrix2[0][1] = dxi_y * dtilda_xi;
			temp_matrix2[1][0] = deta_x * dtilda_eta;
			temp_matrix2[1][1] = deta_y * dtilda_eta;

			InverseMatrix_2D(temp_matrix2);

			strain_xx = temp_matrix2[0][0] * dtilda_xi * dxi_disp_x + temp_matrix2[0][1] * dtilda_eta * deta_disp_x;
			strain_yy = temp_matrix2[1][0] * dtilda_xi * dxi_disp_y + temp_matrix2[1][1] * dtilda_eta * deta_disp_y;
			strain_xy = temp_matrix2[1][0] * dtilda_xi * dxi_disp_x + temp_matrix2[1][1] * dtilda_eta * deta_disp_x + temp_matrix2[0][0] * dtilda_xi * dxi_disp_y + temp_matrix2[0][1] * dtilda_eta * deta_disp_y;

			// 応力の計算を行わないのでコメントアウト
			// double D_matrix[3][3] = {{0.0}};
			// int DM = 1;
			// if (DM == 0) { //平面応力状態
			// 	temp = E * (1.0 - nu * nu);
			// 	D_matrix[0][0] = temp;
			// 	D_matrix[0][1] = nu * temp;
			// 	D_matrix[1][0] = nu * temp;
			// 	D_matrix[1][1] = temp;
			// 	D_matrix[2][2] = (1.0 - nu) / 2.0 * temp;
			// } else if (DM == 1) { //平面ひずみ状態(2Dの場合はこっち)
			// 	temp = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
			// 	D_matrix[0][0] = temp;
			// 	D_matrix[0][1] = nu / (1.0 - nu) * temp;
			// 	D_matrix[1][0] = nu / (1.0 - nu) * temp;
			// 	D_matrix[1][1] = temp;
			// 	D_matrix[2][2] = (1.0 - 2.0 * nu) / 2.0 / (1.0 - nu) * temp;
			// }

			// strainしか使わないのでコメントアウト
			// stress_xx = D_matrix[0][0] * strain_xx + D_matrix[0][1] * strain_yy;
			// stress_yy = D_matrix[1][0] * strain_xx + D_matrix[1][1] * strain_yy;
			// stress_xy = D_matrix[2][2] * strain_xy;

			temp = sqrt(disp_x * disp_x + disp_y * disp_y);
			// printf("x:   % 1.8e\n", px);
			// printf("y:   % 1.8e\n", py);
			// printf("xi:  % 1.8e\n", temp_xi);
			// printf("eta: % 1.8e\n", temp_eta);
			// printf("Displacement x: % 1.8e\n", disp_x);
			// printf("Displacement y: % 1.8e\n", disp_y);
			// printf("Displacement  : % 1.8e\n", temp);
			// printf("Strain xx: % 1.8e\n", strain_xx);
			// printf("Strain yy: % 1.8e\n", strain_yy);
			// printf("Strain xy: % 1.8e\n", strain_xy);
			// printf("Stress xx: % 1.8e\n", stress_xx);
			// printf("Stress yy: % 1.8e\n", stress_yy);
			// printf("Stress xy: % 1.8e\n", stress_xy);

			(*disp_x_glo) = disp_x;
			(*disp_y_glo) = disp_y;

			(*strain_xx_glo) = strain_xx;
			(*strain_yy_glo) = strain_yy;
			(*strain_xy_glo) = strain_xy;

			return i;
		}

		InverseMatrix_2D(temp_matrix);

		temp_dxi = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y);
		temp_deta = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_deta;
		if (temp_xi < input_knot_vec_xi[0])
			temp_xi = input_knot_vec_xi[0];
		if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
			temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
		if (temp_eta < input_knot_vec_eta[0])
			temp_eta = input_knot_vec_eta[0];
		if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
			temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];

		// temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
		// printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
	}
	temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;
	for (i = 0; i < repeat; i++)
	{
		lrNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
						cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
						weight, order_xi, order_eta,
						temp_xi, temp_eta,
						&temp_x, &temp_y,
						&temp_matrix[0][0], &temp_matrix[0][1],
						&temp_matrix[1][0], &temp_matrix[1][1]);

		temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;

		//収束した場合////////////////////////////////////////////////////////////////
		// if (temp_tol_x < tol && temp_tol_y < tol) {
		if (temp_tol_x + temp_tol_y < tol)
		{
			// printf("lrNURBS\n");
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;

			int knot_n_xi = cntl_p_n_xi + order_xi + 1;
			int knot_n_eta = cntl_p_n_eta + order_eta + 1;
			double dtilda_xi = 0.0;
			double dtilda_eta = 0.0;
			double disp_x, disp_y;
			double dxi_x, deta_x, dxi_y, deta_y;
			double dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y;
			double temp_matrix2[2][2];
			double temp;
			double strain_xx, strain_yy, strain_xy;
			// double stress_xx, stress_yy, stress_xy;

			for (i = 0; i < knot_n_xi; i++)
			{
				if (input_knot_vec_xi[i] <= temp_xi && temp_xi < input_knot_vec_xi[i + 1])
				{
					dtilda_xi = (input_knot_vec_xi[i + 1] - input_knot_vec_xi[i]) / 2.0;
					// printf("%f\n", dtilda_xi);
					break;
				}
			}
			for (i = 0; i < knot_n_eta; i++)
			{
				if (input_knot_vec_eta[i] < temp_eta && temp_eta <= input_knot_vec_eta[i + 1])
				{
					dtilda_eta = (input_knot_vec_eta[i + 1] - input_knot_vec_eta[i]) / 2.0;
					// printf("%f\n", dtilda_eta);
					break;
				}
			}

			lrNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
							cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
							weight, order_xi, order_eta,
							temp_xi, temp_eta,
							&temp, &temp,
							&dxi_x, &deta_x,
							&dxi_y, &deta_y);
			// printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
			//	   dxi_x, deta_x, dxi_y, deta_y);

			lrNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
							disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
							weight, order_xi, order_eta,
							temp_xi, temp_eta,
							&disp_x, &disp_y,
							&dxi_disp_x, &deta_disp_x,
							&dxi_disp_y, &deta_disp_y);
			// printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
			//	   dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y);

			temp_matrix2[0][0] = dxi_x * dtilda_xi;
			temp_matrix2[0][1] = dxi_y * dtilda_xi;
			temp_matrix2[1][0] = deta_x * dtilda_eta;
			temp_matrix2[1][1] = deta_y * dtilda_eta;

			InverseMatrix_2D(temp_matrix2);

			strain_xx = temp_matrix2[0][0] * dtilda_xi * dxi_disp_x + temp_matrix2[0][1] * dtilda_eta * deta_disp_x;
			strain_yy = temp_matrix2[1][0] * dtilda_xi * dxi_disp_y + temp_matrix2[1][1] * dtilda_eta * deta_disp_y;
			strain_xy = temp_matrix2[1][0] * dtilda_xi * dxi_disp_x + temp_matrix2[1][1] * dtilda_eta * deta_disp_x + temp_matrix2[0][0] * dtilda_xi * dxi_disp_y + temp_matrix2[0][1] * dtilda_eta * deta_disp_y;

			// 応力の計算を行わないのでコメントアウト
			// double D_matrix[3][3] = {{0.0}};
			// int DM = 1;
			// if (DM == 0) { //平面応力状態
			// 	temp = E * (1.0 - nu * nu);
			// 	D_matrix[0][0] = temp;
			// 	D_matrix[0][1] = nu * temp;
			// 	D_matrix[1][0] = nu * temp;
			// 	D_matrix[1][1] = temp;
			// 	D_matrix[2][2] = (1.0 - nu) / 2.0 * temp;
			// } else if (DM == 1) { //平面ひずみ状態(2Dの場合はこっち)
			// 	temp = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
			// 	D_matrix[0][0] = temp;
			// 	D_matrix[0][1] = nu / (1.0 - nu) * temp;
			// 	D_matrix[1][0] = nu / (1.0 - nu) * temp;
			// 	D_matrix[1][1] = temp;
			// 	D_matrix[2][2] = (1.0 - 2.0 * nu) / 2.0 / (1.0 - nu) * temp;
			// }

			// strainしか使わないのでコメントアウト
			// stress_xx = D_matrix[0][0] * strain_xx + D_matrix[0][1] * strain_yy;
			// stress_yy = D_matrix[1][0] * strain_xx + D_matrix[1][1] * strain_yy;
			// stress_xy = D_matrix[2][2] * strain_xy;

			temp = sqrt(disp_x * disp_x + disp_y * disp_y);
			// printf("x:   % 1.8e\n", px);
			// printf("y:   % 1.8e\n", py);
			// printf("xi:  % 1.8e\n", temp_xi);
			// printf("eta: % 1.8e\n", temp_eta);
			// printf("Displacement x: % 1.8e\n", disp_x);
			// printf("Displacement y: % 1.8e\n", disp_y);
			// printf("Displacement  : % 1.8e\n", temp);
			// printf("Strain xx: % 1.8e\n", strain_xx);
			// printf("Strain yy: % 1.8e\n", strain_yy);
			// printf("Strain xy: % 1.8e\n", strain_xy);
			// printf("Stress xx: % 1.8e\n", stress_xx);
			// printf("Stress yy: % 1.8e\n", stress_yy);
			// printf("Stress xy: % 1.8e\n", stress_xy);

			(*disp_x_glo) = disp_x;
			(*disp_y_glo) = disp_y;

			(*strain_xx_glo) = strain_xx;
			(*strain_yy_glo) = strain_yy;
			(*strain_xy_glo) = strain_xy;

			return i;
		}

		InverseMatrix_2D(temp_matrix);

		temp_dxi = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y);
		temp_deta = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_deta;
		if (temp_xi < input_knot_vec_xi[0])
			temp_xi = input_knot_vec_xi[0];
		if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
			temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
		if (temp_eta < input_knot_vec_eta[0])
			temp_eta = input_knot_vec_eta[0];
		if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
			temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];

		// temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
		// printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
	}
	// printf("% 1.8e % 1.8e\n", temp_x, temp_y);
	return 0;
}

void Calculation(int order_xi, int order_eta,
				 int knot_n_xi, int knot_n_eta,
				 int cntl_p_n_xi, int cntl_p_n_eta,
				 double *input_knot_vec_xi, double *input_knot_vec_eta,
				 double *cntl_px, double *cntl_py,
				 double *disp_cntl_px, double *disp_cntl_py,
				 double *weight)
{
	int i, j, k, l;
	double temp1, temp2, temp3;
	double temp_matrix[2][2];

	//計算するξ,ηの値決定と ∂ξ/∂チルダξ, ∂η/∂チルダη の計算
	double calc_xi[MAX_POINTS];	  //計算するξの値
	double calc_eta[MAX_POINTS];  //計算するηの値
	double dtilda_xi[MAX_KNOTS];  // ∂ξ/∂チルダξ
	double dtilda_eta[MAX_KNOTS]; // ∂η/∂チルダη
	k = 0;
	l = 0;
	for (i = 0; i < knot_n_xi - 1; i++)
	{
		if (input_knot_vec_xi[i] != input_knot_vec_xi[i + 1])
		{
			calc_xi[k] = input_knot_vec_xi[i];
			printf("%d\t%f\n", k, calc_xi[k]);
			dtilda_xi[l] = (input_knot_vec_xi[i + 1] - input_knot_vec_xi[i]) / 2.0;
			printf("%d\t%f\n", k, dtilda_xi[k]);
			k++;
			l++;
			if (division_ele_xi > 1)
			{
				temp1 = (input_knot_vec_xi[i + 1] - input_knot_vec_xi[i]) / (double)division_ele_xi;
				for (j = 1; j < division_ele_xi; j++)
				{
					calc_xi[k] = calc_xi[k - 1] + temp1;
					printf("%d\t%f\n", k, calc_xi[k]);
					k++;
				}
			}
		}
	}
	calc_xi[k] = input_knot_vec_xi[knot_n_xi - 1];
	printf("%d\t%f\n", k, calc_xi[k]);
	// printf("\n");
	division_n_xi = k + 1;
	element_n_xi = l;

	k = 0;
	l = 0;
	for (i = 0; i < knot_n_eta - 1; i++)
	{
		if (input_knot_vec_eta[i] != input_knot_vec_eta[i + 1])
		{
			calc_eta[k] = input_knot_vec_eta[i];
			// printf("%d\t%f\n", k, calc_eta[k]);
			dtilda_eta[l] = (input_knot_vec_eta[i + 1] - input_knot_vec_eta[i]) / 2.0;
			// printf("%d\t%f\n", k, dtilda_eta[k]);
			k++;
			l++;
			if (division_ele_eta > 1)
			{
				temp1 = (input_knot_vec_eta[i + 1] - input_knot_vec_eta[i]) / (double)division_ele_eta;
				for (j = 1; j < division_ele_eta; j++)
				{
					calc_eta[k] = calc_eta[k - 1] + temp1;
					// printf("%d\t%f\n", k, calc_eta[k]);
					k++;
				}
			}
		}
	}
	calc_eta[k] = input_knot_vec_eta[knot_n_eta - 1];
	// printf("%d\t%f\n", k, calc_eta[k]);
	// printf("\n");
	division_n_eta = k + 1;
	element_n_eta = l;

	if (element_n_xi > MAX_ELEMENTS)
	{
		printf("Error!!\n");
		printf("Too many elements at xi!\n"
			   "Maximum of elements is %d (Now %d)\n"
			   "\n",
			   MAX_ELEMENTS, element_n_xi);
		exit(1);
	}
	if (element_n_eta > MAX_ELEMENTS)
	{
		printf("Error!!\n");
		printf("Too many elements at eta!\n"
			   "Maximum of elements is %d (Now %d)\n"
			   "\n",
			   MAX_ELEMENTS, element_n_eta);
		exit(1);
	}

	int ii, jj, kk, ll;

	//メッシュ座標計算
	printf("Start Calculation mesh\n\n");
	for (i = 0; i < division_n_xi; i++)
	{
		ii = i / division_ele_xi;
		kk = i % division_ele_xi;
		for (j = 0; j < division_n_eta; j++)
		{
			jj = j / division_ele_eta;
			ll = j % division_ele_eta;
			lNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
						   cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
						   weight, order_xi, order_eta,
						   calc_xi[i], calc_eta[j],
						   &coord_x[i][j], &coord_y[i][j],
						   &dxi_x[ii][jj][kk][ll], &deta_x[ii][jj][kk][ll],
						   &dxi_y[ii][jj][kk][ll], &deta_y[ii][jj][kk][ll]);
			// printf("[%d][%d] [%d][%d][%d][%d]"
			// 	   "% 1.4e % 1.4e "
			// 	   "% 1.4e % 1.4e\n",
			// 	   i, j, ii, jj, kk, ll,
			// 	   calc_xi[i], calc_eta[j],
			// 	   coord_x[i][j], coord_y[i][j]);
		}
		// printf("\n");
	}
	printf("\n");
	printf("End Calculation mesh\n\n");

	if (fields_flag)
	{
		//変位計算
		printf("Start Calculation displpacement\n\n");
		for (i = 0; i < division_n_xi; i++)
		{
			ii = i / division_ele_xi;
			kk = i % division_ele_xi;
			for (j = 0; j < division_n_eta; j++)
			{
				jj = j / division_ele_eta;
				ll = j % division_ele_eta;
				lNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
							   disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
							   weight, order_xi, order_eta,
							   calc_xi[i], calc_eta[j],
							   &disp_x[i][j], &disp_y[i][j],
							   &dxi_disp_x[ii][jj][kk][ll], &deta_disp_x[ii][jj][kk][ll],
							   &dxi_disp_y[ii][jj][kk][ll], &deta_disp_y[ii][jj][kk][ll]);
				// printf("[%d][%d] [%d][%d][%d][%d]"
				// 	   "% 1.4e % 1.4e "
				// 	   "% 1.4e % 1.4e\n",
				// 	   i, j, ii, jj, kk, ll,
				// 	   calc_xi[i], calc_eta[j],
				// 	   disp_x[i][j], disp_y[i][j]);
			}
			// printf("\n");
		}
		printf("\n");
		printf("End Calculation displpacement\n\n");

		//足りない微分値計算
		for (ii = 0; ii < element_n_xi; ii++)
		{
			for (jj = 0; jj < element_n_eta; jj++)
			{
				kk = division_ele_xi;
				i = (ii + 1) * division_ele_xi;
				j = jj * division_ele_eta;
				for (ll = 1; ll < division_ele_eta; ll++)
				{
					j++;
					rNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
								   cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
								   weight, order_xi, order_eta,
								   calc_xi[i], calc_eta[j],
								   &coord_x[i][j], &coord_y[i][j],
								   &dxi_x[ii][jj][kk][ll], &deta_x[ii][jj][kk][ll],
								   &dxi_y[ii][jj][kk][ll], &deta_y[ii][jj][kk][ll]);
					rNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
								   disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
								   weight, order_xi, order_eta,
								   calc_xi[i], calc_eta[j],
								   &disp_x[i][j], &disp_y[i][j],
								   &dxi_disp_x[ii][jj][kk][ll], &deta_disp_x[ii][jj][kk][ll],
								   &dxi_disp_y[ii][jj][kk][ll], &deta_disp_y[ii][jj][kk][ll]);
					/*
					printf("[%d][%d] [%d][%d][%d][%d]"
						   "% 1.4e % 1.4e % 1.4e % 1.4e"
						   "% 1.4e % 1.4e % 1.4e % 1.4e\n",
						   i, j, ii, jj, kk, ll,
						   dxi_x[ii][jj][kk][ll], deta_x[ii][jj][kk][ll],
						   dxi_y[ii][jj][kk][ll], deta_y[ii][jj][kk][ll],
						   dxi_disp_x[ii][jj][kk][ll], deta_disp_x[ii][jj][kk][ll],
						   dxi_disp_y[ii][jj][kk][ll], deta_disp_y[ii][jj][kk][ll]);
					*/
				}

				ll = division_ele_eta;
				i = ii * division_ele_xi;
				j = (jj + 1) * division_ele_eta;
				for (kk = 1; kk <= division_ele_xi; kk++)
				{
					i++;
					rNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
								   cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
								   weight, order_xi, order_eta,
								   calc_xi[i], calc_eta[j],
								   &coord_x[i][j], &coord_y[i][j],
								   &dxi_x[ii][jj][kk][ll], &deta_x[ii][jj][kk][ll],
								   &dxi_y[ii][jj][kk][ll], &deta_y[ii][jj][kk][ll]);
					rNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
								   disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
								   weight, order_xi, order_eta,
								   calc_xi[i], calc_eta[j],
								   &disp_x[i][j], &disp_y[i][j],
								   &dxi_disp_x[ii][jj][kk][ll], &deta_disp_x[ii][jj][kk][ll],
								   &dxi_disp_y[ii][jj][kk][ll], &deta_disp_y[ii][jj][kk][ll]);
					/*
					printf("[%d][%d] [%d][%d][%d][%d]"
						   "% 1.4e % 1.4e % 1.4e % 1.4e"
						   "% 1.4e % 1.4e % 1.4e % 1.4e\n",
						   i, j, ii, jj, kk, ll,
						   dxi_x[ii][jj][kk][ll], deta_x[ii][jj][kk][ll],
						   dxi_y[ii][jj][kk][ll], deta_y[ii][jj][kk][ll],
						   dxi_disp_x[ii][jj][kk][ll], deta_disp_x[ii][jj][kk][ll],
						   dxi_disp_y[ii][jj][kk][ll], deta_disp_y[ii][jj][kk][ll]);
					*/
				}

				kk = division_ele_xi;
				ll = 0;
				i = (ii + 1) * division_ele_xi;
				j = jj * division_ele_eta;
				rlNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
								cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
								weight, order_xi, order_eta,
								calc_xi[i], calc_eta[j],
								&coord_x[i][j], &coord_y[i][j],
								&dxi_x[ii][jj][kk][ll], &deta_x[ii][jj][kk][ll],
								&dxi_y[ii][jj][kk][ll], &deta_y[ii][jj][kk][ll]);
				rlNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
								disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
								weight, order_xi, order_eta,
								calc_xi[i], calc_eta[j],
								&disp_x[i][j], &disp_y[i][j],
								&dxi_disp_x[ii][jj][kk][ll], &deta_disp_x[ii][jj][kk][ll],
								&dxi_disp_y[ii][jj][kk][ll], &deta_disp_y[ii][jj][kk][ll]);
				/*
					printf("[%d][%d] [%d][%d][%d][%d]"
						   "% 1.4e % 1.4e % 1.4e % 1.4e"
						   "% 1.4e % 1.4e % 1.4e % 1.4e\n",
						   i, j, ii, jj, kk, ll,
						   dxi_x[ii][jj][kk][ll], deta_x[ii][jj][kk][ll],
						   dxi_y[ii][jj][kk][ll], deta_y[ii][jj][kk][ll],
						   dxi_disp_x[ii][jj][kk][ll], deta_disp_x[ii][jj][kk][ll],
						   dxi_disp_y[ii][jj][kk][ll], deta_disp_y[ii][jj][kk][ll]);
				*/

				kk = 0;
				ll = division_ele_eta;
				i = ii * division_ele_xi;
				j = (jj + 1) * division_ele_eta;
				lrNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
								cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
								weight, order_xi, order_eta,
								calc_xi[i], calc_eta[j],
								&coord_x[i][j], &coord_y[i][j],
								&dxi_x[ii][jj][kk][ll], &deta_x[ii][jj][kk][ll],
								&dxi_y[ii][jj][kk][ll], &deta_y[ii][jj][kk][ll]);
				lrNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
								disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
								weight, order_xi, order_eta,
								calc_xi[i], calc_eta[j],
								&disp_x[i][j], &disp_y[i][j],
								&dxi_disp_x[ii][jj][kk][ll], &deta_disp_x[ii][jj][kk][ll],
								&dxi_disp_y[ii][jj][kk][ll], &deta_disp_y[ii][jj][kk][ll]);
				/*
					printf("[%d][%d] [%d][%d][%d][%d]"
						   "% 1.4e % 1.4e % 1.4e % 1.4e"
						   "% 1.4e % 1.4e % 1.4e % 1.4e\n",
						   i, j, ii, jj, kk, ll,
						   dxi_x[ii][jj][kk][ll], deta_x[ii][jj][kk][ll],
						   dxi_y[ii][jj][kk][ll], deta_y[ii][jj][kk][ll],
						   dxi_disp_x[ii][jj][kk][ll], deta_disp_x[ii][jj][kk][ll],
						   dxi_disp_y[ii][jj][kk][ll], deta_disp_y[ii][jj][kk][ll]);
				printf("\n");
				*/
			}
		}

		/*
		for (ii = 0; ii < element_n_xi; ii++) {
			for (jj = 0; jj < element_n_eta; jj++) {
				for (kk = 0; kk <= division_ele_xi; kk++) {
					for (ll = 0; ll <= division_ele_eta; ll++) {
						printf("[%d][%d][%d][%d]"
							   "% 1.4e % 1.4e % 1.4e % 1.4e"
							   "% 1.4e % 1.4e % 1.4e % 1.4e\n",
							   ii, jj, kk, ll,
							   dxi_x[ii][jj][kk][ll], deta_x[ii][jj][kk][ll],
							   dxi_y[ii][jj][kk][ll], deta_y[ii][jj][kk][ll],
							   dxi_disp_x[ii][jj][kk][ll], deta_disp_x[ii][jj][kk][ll],
							   dxi_disp_y[ii][jj][kk][ll], deta_disp_y[ii][jj][kk][ll]);
					}
				}
			}
		}
		printf("\n");
		*/

		//ひずみ計算
		printf("Start Calculation Strain\n\n");
		for (i = 0; i < element_n_xi; i++)
		{
			for (j = 0; j < element_n_eta; j++)
			{
				temp1 = dtilda_xi[i];
				temp2 = dtilda_eta[j];
				for (k = 0; k < division_ele_xi + 1; k++)
				{
					for (l = 0; l < division_ele_eta + 1; l++)
					{
						temp_matrix[0][0] = dxi_x[i][j][k][l] * temp1;
						temp_matrix[0][1] = dxi_y[i][j][k][l] * temp1;
						temp_matrix[1][0] = deta_x[i][j][k][l] * temp2;
						temp_matrix[1][1] = deta_y[i][j][k][l] * temp2;

						InverseMatrix_2D(temp_matrix);

						strain_xx[i][j][k][l] = temp_matrix[0][0] * temp1 * dxi_disp_x[i][j][k][l] + temp_matrix[0][1] * temp2 * deta_disp_x[i][j][k][l];
						strain_yy[i][j][k][l] = temp_matrix[1][0] * temp1 * dxi_disp_y[i][j][k][l] + temp_matrix[1][1] * temp2 * deta_disp_y[i][j][k][l];
						strain_xy[i][j][k][l] = temp_matrix[1][0] * temp1 * dxi_disp_x[i][j][k][l] + temp_matrix[1][1] * temp2 * deta_disp_x[i][j][k][l] + temp_matrix[0][0] * temp1 * dxi_disp_y[i][j][k][l] + temp_matrix[0][1] * temp2 * deta_disp_y[i][j][k][l];

						// printf("[%d][%d][%d][%d]\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\n",
						// 	   i, j, k, l, temp1, temp2,
						// 	   strain_xx[i][j][k][l], strain_yy[i][j][k][l], strain_xy[i][j][k][l]);
					}
				}
			}
			// printf("\n");
		}
		printf("End Calculation Strain\n\n");

		// Dマトリクスの計算
		double D_matrix[3][3] = {{0.0}};
		if (DM == 0)
		{ //平面応力状態
			temp1 = E * (1.0 - nu * nu);
			D_matrix[0][0] = temp1;
			D_matrix[0][1] = nu * temp1;
			D_matrix[1][0] = nu * temp1;
			D_matrix[1][1] = temp1;
			D_matrix[2][2] = (1.0 - nu) / 2.0 * temp1;
		}
		else if (DM == 1)
		{ //平面ひずみ状態(2Dの場合はこっち)
			temp1 = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
			D_matrix[0][0] = temp1;
			D_matrix[0][1] = nu / (1.0 - nu) * temp1;
			D_matrix[1][0] = nu / (1.0 - nu) * temp1;
			D_matrix[1][1] = temp1;
			D_matrix[2][2] = (1.0 - 2.0 * nu) / 2.0 / (1.0 - nu) * temp1;
		}

		printf("Start Calculation Stress\n\n");
		for (i = 0; i < element_n_xi; i++)
		{
			for (j = 0; j < element_n_eta; j++)
			{
				for (k = 0; k < division_ele_xi + 1; k++)
				{
					for (l = 0; l < division_ele_eta + 1; l++)
					{
						stress_xx[i][j][k][l] = D_matrix[0][0] * strain_xx[i][j][k][l] + D_matrix[0][1] * strain_yy[i][j][k][l];
						stress_yy[i][j][k][l] = D_matrix[1][0] * strain_xx[i][j][k][l] + D_matrix[1][1] * strain_yy[i][j][k][l];
						stress_xy[i][j][k][l] = D_matrix[2][2] * strain_xy[i][j][k][l];
						// printf("[%d][%d][%d][%d]\t% 1.4e\t% 1.4e\t% 1.4e\n", i, j, k, l,
						// 	   stress_xx[i][j][k][l], stress_yy[i][j][k][l], stress_xy[i][j][k][l]);
					}
				}
				// printf("\n");
			}
		}
		printf("End Calculation Stress\n\n");
	}

	//書き込み
	fp = fopen("view.dat", "a");
	if (fields_flag)
	{
		fprintf(fp, "%d\t%d\t%d\t%d\n",
				division_n_xi, division_n_eta,
				element_n_xi, element_n_eta);
		for (i = 0; i < division_n_xi; i++)
		{
			for (j = 0; j < division_n_eta; j++)
			{
				temp1 = disp_x[i][j];
				temp2 = disp_y[i][j];
				temp3 = sqrt(temp1 * temp1 + temp2 * temp2);
				fprintf(fp, "% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\n",
						coord_x[i][j], coord_y[i][j],
						disp_x[i][j], disp_y[i][j], temp3);
			}
		}
		for (i = 0; i < element_n_xi; i++)
		{
			for (j = 0; j < element_n_eta; j++)
			{
				for (k = 0; k < division_ele_xi + 1; k++)
				{
					for (l = 0; l < division_ele_eta + 1; l++)
					{
						fprintf(fp, "% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\n",
								strain_xx[i][j][k][l], strain_yy[i][j][k][l], strain_xy[i][j][k][l],
								stress_xx[i][j][k][l], stress_yy[i][j][k][l], stress_xy[i][j][k][l]);
					}
				}
			}
		}
	}
	else
	{
		fprintf(fp, "%d\t%d\t%d\t%d\n",
				division_n_xi, division_n_eta,
				element_n_xi, element_n_eta);
		for (i = 0; i < division_n_xi; i++)
		{
			for (j = 0; j < division_n_eta; j++)
			{
				fprintf(fp, "% 1.4e\t% 1.4e\n",
						coord_x[i][j], coord_y[i][j]);
			}
		}
	}
	fclose(fp);
	// machino
	fp = fopen("view_r_theta.dat", "a");
	if (fields_flag)
	{
		fprintf(fp, "%d\t%d\t%d\t%d\n",
				division_n_xi, division_n_eta,
				element_n_xi, element_n_eta);
		for (i = 0; i < division_n_xi; i++)
		{
			for (j = 0; j < division_n_eta; j++)
			{
				temp1 = disp_x[i][j];
				temp2 = disp_y[i][j];
				temp3 = sqrt(temp1 * temp1 + temp2 * temp2);
				fprintf(fp, "% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\n",
						coord_x[i][j], coord_y[i][j],
						disp_x[i][j], disp_y[i][j], temp3);
			}
		}
		for (i = 0; i < element_n_xi; i++)
		{
			for (j = 0; j < element_n_eta; j++)
			{
				for (k = 0; k < division_ele_xi + 1; k++)
				{
					for (l = 0; l < division_ele_eta + 1; l++)
					{
						double stress_rr, stress_theta;
						double sum = stress_xx[i][j][k][l] + stress_yy[i][j][k][l];
						double dif = stress_xx[i][j][k][l] - stress_yy[i][j][k][l];
						double tau2 = stress_xy[i][j][k][l] * stress_xy[i][j][k][l];
						stress_rr = sum * 0.5 + sqrt(dif * dif + 4 * tau2) * 0.5;
						stress_theta = sum * 0.5 - sqrt(dif * dif + 4 * tau2) * 0.5;
						fprintf(fp, "% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t0.0\n",
								strain_xx[i][j][k][l], strain_yy[i][j][k][l], strain_xy[i][j][k][l],
								stress_rr, stress_theta);
					}
				}
			}
		}
	}
	else
	{
		fprintf(fp, "%d\t%d\t%d\t%d\n",
				division_n_xi, division_n_eta,
				element_n_xi, element_n_eta);
		for (i = 0; i < division_n_xi; i++)
		{
			for (j = 0; j < division_n_eta; j++)
			{
				fprintf(fp, "% 1.4e\t% 1.4e\n",
						coord_x[i][j], coord_y[i][j]);
			}
		}
	}
	fclose(fp);
	//グラフ用ファイル書き込み
	fp = fopen("disp_graph.txt", "a");
	for (i = 0; i < division_n_xi; i++)
	{
		for (j = 0; j < division_n_eta; j++)
		{
			temp1 = disp_x[i][j];
			temp2 = disp_y[i][j];
			// temp3 = sqrt(temp1 * temp1 + temp2 * temp2);
			fprintf(fp, "%d\t% 1.15e\t% 1.15e\t% 1.15e\t% 1.15e\n",
					graph_patch_n,
					coord_x[i][j], coord_y[i][j],
					disp_x[i][j], disp_y[i][j]);
		}
	}
	fclose(fp);

	fp = fopen("stress_y_graph.txt", "a");
	for (i = 0; i < element_n_xi; i++)
	{
		for (j = 0; j < element_n_eta; j++)
		{
			for (k = 0; k < division_ele_xi + 1; k++)
			{
				for (l = 0; l < division_ele_eta + 1; l++)
				{
					fprintf(fp, "%d\t% 1.15e\t% 1.15e\t% 1.15e\n",
							graph_patch_n,
							coord_x[i * division_ele_xi + k][j * division_ele_eta + l],
							coord_y[i * division_ele_xi + k][j * division_ele_eta + l],
							stress_yy[i][j][k][l]);
				}
			}
		}
	}
	fclose(fp);

	fp = fopen("stress_y_graph_0.txt", "a");
	for (i = 0; i < element_n_xi; i++)
	{
		for (j = 0; j < element_n_eta; j++)
		{
			for (k = 0; k < division_ele_xi + 1; k++)
			{
				for (l = 0; l < division_ele_eta + 1; l++)
				{
					if (coord_y[i * division_ele_xi + k][j * division_ele_eta + l] == 0.000000000000000e+00)
					{
						fprintf(fp, "%d\t% 1.15e\t% 1.15e\t% 1.15e\n",
								graph_patch_n,
								coord_x[i * division_ele_xi + k][j * division_ele_eta + l],
								coord_y[i * division_ele_xi + k][j * division_ele_eta + l],
								stress_yy[i][j][k][l]);
					}
				}
			}
		}
	}
	fclose(fp);

	fp = fopen("stress_vm_graph.txt", "a");
	for (i = 0; i < element_n_xi; i++)
	{
		for (j = 0; j < element_n_eta; j++)
		{
			for (k = 0; k < division_ele_xi + 1; k++)
			{
				for (l = 0; l < division_ele_eta + 1; l++)
				{
					double stress_vm;
					double sum = stress_xx[i][j][k][l] + stress_yy[i][j][k][l];
					double dif = stress_xx[i][j][k][l] - stress_yy[i][j][k][l];
					double tau2 = stress_xy[i][j][k][l] * stress_xy[i][j][k][l];
					double temp1, temp2;
					double temp3;
					temp1 = 0.5 * sum;
					temp2 = 0.5 * sqrt(dif * dif + 4 * tau2);
					stress_vm = sqrt(temp1 * temp1 + 3 * temp2 * temp2);
					temp3 = sqrt(calc_xi[i * division_ele_xi + k] * calc_xi[i * division_ele_xi + k] + calc_eta[j * division_ele_eta + l] * calc_eta[j * division_ele_eta + l]);
					fprintf(fp, "% 1.10e\t% 1.10e\t% 1.10e\t% 1.10e\t"
								"% 1.10e\t% 1.10e\n",
							calc_xi[i * division_ele_xi + k],
							calc_eta[j * division_ele_eta + l],
							coord_x[i * division_ele_xi + k][j * division_ele_eta + l],
							coord_y[i * division_ele_xi + k][j * division_ele_eta + l],
							stress_vm, temp3);
				}
			}
		}
	}
	fclose(fp);
}

// 重ね合わせた結果の出力
void Calculation_overlay(int order_xi_loc, int order_eta_loc,
						 int knot_n_xi_loc, int knot_n_eta_loc,
						 int cntl_p_n_xi_loc, int cntl_p_n_eta_loc,
						 double *knot_vec_xi_loc, double *knot_vec_eta_loc,
						 double *cntl_px_loc, double *cntl_py_loc,
						 double *weight_loc,
						 int order_xi_glo, int order_eta_glo,
						 int cntl_p_n_xi_glo, int cntl_p_n_eta_glo,
						 double *knot_vec_xi_glo, double *knot_vec_eta_glo,
						 double *cntl_px_glo, double *cntl_py_glo,
						 double *disp_cntl_px_glo, double *disp_cntl_py_glo,
						 double *weight_glo)
{
	int i, j, k, l;
	double temp1, temp2, temp3;
	// double temp_matrix[2][2];

	double output_xi, output_eta;
	double disp_x_glo;
	double disp_y_glo;
	double strain_xx_glo = 0;
	double strain_yy_glo = 0;
	double strain_xy_glo = 0;
	//, strain_yy_glo, strain_xy_glo;

	//計算するξ,ηの値決定と ∂ξ/∂チルダξ, ∂η/∂チルダη の計算
	double calc_xi_loc[MAX_POINTS];	 //計算するξの値local
	double calc_eta_loc[MAX_POINTS]; //計算するηの値local
	// double dtilda_xi[MAX_KNOTS];		// ∂ξ/∂チルダξ
	// double dtilda_eta[MAX_KNOTS];	// ∂η/∂チルダη
	k = 0;
	l = 0;
	for (i = 0; i < knot_n_xi_loc - 1; i++)
	{
		if (knot_vec_xi_loc[i] != knot_vec_xi_loc[i + 1])
		{
			calc_xi_loc[k] = knot_vec_xi_loc[i];
			printf("%d\t%f\n", k, calc_xi_loc[k]);
			// dtilda_xi[l] = ( knot_vec_xi_loc[i + 1] - knot_vec_xi_loc[i] ) / 2.0;
			// printf("%d\t%f\n", k, dtilda_xi[k]);
			k++;
			l++;
			if (division_ele_xi > 1)
			{
				temp1 = (knot_vec_xi_loc[i + 1] - knot_vec_xi_loc[i]) / (double)division_ele_xi;
				for (j = 1; j < division_ele_xi; j++)
				{
					calc_xi_loc[k] = calc_xi_loc[k - 1] + temp1;
					printf("%d\t%f\n", k, calc_xi_loc[k]);
					k++;
				}
			}
		}
	}
	calc_xi_loc[k] = knot_vec_xi_loc[knot_n_xi_loc - 1];
	printf("%d\t%f\n", k, calc_xi_loc[k]);
	// printf("\n");
	division_n_xi = k + 1;
	element_n_xi = l;

	k = 0;
	l = 0;
	for (i = 0; i < knot_n_eta_loc - 1; i++)
	{
		if (knot_vec_eta_loc[i] != knot_vec_eta_loc[i + 1])
		{
			calc_eta_loc[k] = knot_vec_eta_loc[i];
			// printf("%d\t%f\n", k, calc_eta[k]);
			// dtilda_eta[l] = ( knot_vec_eta_loc[i + 1] - knot_vec_eta_loc[i] ) / 2.0;
			// printf("%d\t%f\n", k, dtilda_eta[k]);
			k++;
			l++;
			if (division_ele_eta > 1)
			{
				temp1 = (knot_vec_eta_loc[i + 1] - knot_vec_eta_loc[i]) / (double)division_ele_eta;
				for (j = 1; j < division_ele_eta; j++)
				{
					calc_eta_loc[k] = calc_eta_loc[k - 1] + temp1;
					// printf("%d\t%f\n", k, calc_eta[k]);
					k++;
				}
			}
		}
	}
	calc_eta_loc[k] = knot_vec_eta_loc[knot_n_eta_loc - 1];
	// printf("%d\t%f\n", k, calc_eta[k]);
	// printf("\n");
	division_n_eta = k + 1;
	element_n_eta = l;

	if (element_n_xi > MAX_ELEMENTS)
	{
		printf("Error!!\n");
		printf("Too many elements at xi!\n"
			   "Maximum of elements is %d (Now %d)\n"
			   "\n",
			   MAX_ELEMENTS, element_n_xi);
		exit(1);
	}
	if (element_n_eta > MAX_ELEMENTS)
	{
		printf("Error!!\n");
		printf("Too many elements at eta!\n"
			   "Maximum of elements is %d (Now %d)\n"
			   "\n",
			   MAX_ELEMENTS, element_n_eta);
		exit(1);
	}

	int ii, jj, kk, ll;

	//メッシュ座標計算
	printf("Start Calculation overlay mesh\n\n");
	printf("Start Calculation overlay displpacement\n\n");
	printf("Start Calculation overlay Strain\n\n");
	for (i = 0; i < division_n_xi; i++)
	{
		ii = i / division_ele_xi;
		kk = i % division_ele_xi;
		for (j = 0; j < division_n_eta; j++)
		{
			jj = j / division_ele_eta;
			ll = j % division_ele_eta;
			lNURBS_surface(knot_vec_xi_loc, knot_vec_eta_loc,
						   cntl_px_loc, cntl_py_loc, cntl_p_n_xi_loc, cntl_p_n_eta_loc,
						   weight_loc, order_xi_loc, order_eta_loc,
						   calc_xi_loc[i], calc_eta_loc[j],
						   &coord_x[i][j], &coord_y[i][j],
						   &dxi_x[ii][jj][kk][ll], &deta_x[ii][jj][kk][ll],
						   &dxi_y[ii][jj][kk][ll], &deta_y[ii][jj][kk][ll]);
			// printf("[%d][%d] [%d][%d][%d][%d]"
			// 	   "% 1.4e % 1.4e "
			// 	   "% 1.4e % 1.4e\n",
			// 	   i, j, ii, jj, kk, ll,
			// 	   calc_xi_loc[i], calc_eta_loc[j],
			// 	   coord_x[i][j], coord_y[i][j]);

			int itr_n = CalcXiEtaByNR(coord_x[i][j], coord_y[i][j],
									  knot_vec_xi_glo, knot_vec_eta_glo,
									  cntl_px_glo, cntl_py_glo,
									  disp_cntl_px_glo, disp_cntl_py_glo,
									  cntl_p_n_xi_glo, cntl_p_n_eta_glo,
									  weight_glo, order_xi_glo, order_eta_glo,
									  &output_xi, &output_eta,
									  &disp_x_glo, &disp_y_glo,
									  &strain_xx_glo, &strain_yy_glo, &strain_xy_glo);
			if (itr_n == 0)
			{
				// printf("itr=0\n");
			}
			// printf("iteration : %d\n",itr_n);

			//ローカル内の表示点上のグローバル変位
			// printf("disp_x_glo =% 1.4e\tdisp_y_glo =% 1.4e\n", disp_x_glo, disp_y_glo);
			// printf("%1.4e\t%1.4e\n",disp_x[i][j],disp_y[i][j]);
			disp_x[i][j] += disp_x_glo;
			disp_y[i][j] += disp_y_glo;
			// printf("% 1.4e\t% 1.4e\n",disp_x[i][j],disp_y[i][j]);

			//ローカル内の表示点上のグローバルひずみ
			// printf("strain_xx_glo =% 1.4e\n"
			//        "strain_yy_glo =% 1.4e\n"
			//        "strain_xy_glo =% 1.4e\n",
			//        strain_xx_glo, strain_yy_glo, strain_xy_glo);
			strain_xx[ii][jj][kk][ll] += strain_xx_glo;
			strain_yy[ii][jj][kk][ll] += strain_yy_glo;
			strain_xy[ii][jj][kk][ll] += strain_xy_glo;
			// printf("test[%d][%d][%d][%d]\n",ii,jj,kk,ll);
			if (jj > 0 && ll == 0)
			{
				strain_xx[ii][jj - 1][kk][division_ele_eta] += strain_xx_glo;
				strain_yy[ii][jj - 1][kk][division_ele_eta] += strain_yy_glo;
				strain_xy[ii][jj - 1][kk][division_ele_eta] += strain_xy_glo;
				// printf("test[%d][%d][%d][%d]\n",ii,jj-1,kk,division_ele_eta);
			}
			if (ii > 0 && kk == 0)
			{
				strain_xx[ii - 1][jj][division_ele_xi][ll] += strain_xx_glo;
				strain_yy[ii - 1][jj][division_ele_xi][ll] += strain_yy_glo;
				strain_xy[ii - 1][jj][division_ele_xi][ll] += strain_xy_glo;
				// printf("test[%d][%d][%d][%d]\n",ii-1,jj,division_ele_xi,ll);
			}
			if (ii > 0 && jj > 0 && kk == 0 && ll == 0)
			{
				strain_xx[ii - 1][jj - 1][division_ele_xi][division_ele_eta] += strain_xx_glo;
				strain_yy[ii - 1][jj - 1][division_ele_xi][division_ele_eta] += strain_yy_glo;
				strain_xy[ii - 1][jj - 1][division_ele_xi][division_ele_eta] += strain_xy_glo;
				// printf("test[%d][%d][%d][%d]\n",ii-1,jj-1,division_ele_xi,division_ele_eta);
			}
			// printf("% 1.4e\t% 1.4e\t% 1.4e\n",
			// 		strain_xx[ii][jj][kk][ll],
			// 		strain_yy[ii][jj][kk][ll],
			// 		strain_xy[ii][jj][kk][ll]);
		}
		// printf("\n");
	}
	printf("\n");
	printf("End Calculation overlay mesh\n\n");
	printf("End Calculation overlay displpacement\n\n");
	printf("End Calculation overlay Strain\n\n");

	// Dマトリクスの計算
	double D_matrix[3][3] = {{0.0}};
	if (DM == 0)
	{ //平面応力状態
		temp1 = E * (1.0 - nu * nu);
		D_matrix[0][0] = temp1;
		D_matrix[0][1] = nu * temp1;
		D_matrix[1][0] = nu * temp1;
		D_matrix[1][1] = temp1;
		D_matrix[2][2] = (1.0 - nu) / 2.0 * temp1;
	}
	else if (DM == 1)
	{ //平面ひずみ状態(2Dの場合はこっち)
		temp1 = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
		D_matrix[0][0] = temp1;
		D_matrix[0][1] = nu / (1.0 - nu) * temp1;
		D_matrix[1][0] = nu / (1.0 - nu) * temp1;
		D_matrix[1][1] = temp1;
		D_matrix[2][2] = (1.0 - 2.0 * nu) / 2.0 / (1.0 - nu) * temp1;
	}

	printf("Start Calculation overlay Stress\n\n");
	for (i = 0; i < element_n_xi; i++)
	{
		for (j = 0; j < element_n_eta; j++)
		{
			for (k = 0; k < division_ele_xi + 1; k++)
			{
				for (l = 0; l < division_ele_eta + 1; l++)
				{
					stress_xx[i][j][k][l] = D_matrix[0][0] * strain_xx[i][j][k][l] + D_matrix[0][1] * strain_yy[i][j][k][l];
					stress_yy[i][j][k][l] = D_matrix[1][0] * strain_xx[i][j][k][l] + D_matrix[1][1] * strain_yy[i][j][k][l];
					stress_xy[i][j][k][l] = D_matrix[2][2] * strain_xy[i][j][k][l];
					// printf("[%d][%d][%d][%d]\t% 1.4e\t% 1.4e\t% 1.4e\n", i, j, k, l,
					// 	   stress_xx[i][j][k][l], stress_yy[i][j][k][l], stress_xy[i][j][k][l]);
				}
			}
			// printf("\n");
		}
	}
	printf("End Calculation overlay Stress\n\n");

	//書き込み
	fp = fopen("overlay_view.dat", "a");
	if (fields_flag)
	{
		fprintf(fp, "%d\t%d\t%d\t%d\n",
				division_n_xi, division_n_eta,
				element_n_xi, element_n_eta);
		for (i = 0; i < division_n_xi; i++)
		{
			for (j = 0; j < division_n_eta; j++)
			{
				temp1 = disp_x[i][j];
				temp2 = disp_y[i][j];
				temp3 = sqrt(temp1 * temp1 + temp2 * temp2);
				fprintf(fp, "% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\n",
						coord_x[i][j], coord_y[i][j],
						disp_x[i][j], disp_y[i][j], temp3);
			}
		}
		for (i = 0; i < element_n_xi; i++)
		{
			for (j = 0; j < element_n_eta; j++)
			{
				for (k = 0; k < division_ele_xi + 1; k++)
				{
					for (l = 0; l < division_ele_eta + 1; l++)
					{
						fprintf(fp, "% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\n",
								strain_xx[i][j][k][l], strain_yy[i][j][k][l], strain_xy[i][j][k][l],
								stress_xx[i][j][k][l], stress_yy[i][j][k][l], stress_xy[i][j][k][l]);
					}
				}
			}
		}
	}
	else
	{
		fprintf(fp, "%d\t%d\t%d\t%d\n",
				division_n_xi, division_n_eta,
				element_n_xi, element_n_eta);
		for (i = 0; i < division_n_xi; i++)
		{
			for (j = 0; j < division_n_eta; j++)
			{
				fprintf(fp, "% 1.4e\t% 1.4e\n",
						coord_x[i][j], coord_y[i][j]);
			}
		}
	}
	fclose(fp);
	// machino
	fp = fopen("overlay_view_r_theta.dat", "a");
	if (fields_flag)
	{
		fprintf(fp, "%d\t%d\t%d\t%d\n",
				division_n_xi, division_n_eta,
				element_n_xi, element_n_eta);
		for (i = 0; i < division_n_xi; i++)
		{
			for (j = 0; j < division_n_eta; j++)
			{
				temp1 = disp_x[i][j];
				temp2 = disp_y[i][j];
				temp3 = sqrt(temp1 * temp1 + temp2 * temp2);
				fprintf(fp, "% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\n",
						coord_x[i][j], coord_y[i][j],
						disp_x[i][j], disp_y[i][j], temp3);
			}
		}
		for (i = 0; i < element_n_xi; i++)
		{
			for (j = 0; j < element_n_eta; j++)
			{
				for (k = 0; k < division_ele_xi + 1; k++)
				{
					for (l = 0; l < division_ele_eta + 1; l++)
					{
						double stress_rr, stress_theta;
						double sum = stress_xx[i][j][k][l] + stress_yy[i][j][k][l];
						double dif = stress_xx[i][j][k][l] - stress_yy[i][j][k][l];
						double tau2 = stress_xy[i][j][k][l] * stress_xy[i][j][k][l];
						stress_rr = sum * 0.5 + sqrt(dif * dif + 4 * tau2) * 0.5;
						stress_theta = sum * 0.5 - sqrt(dif * dif + 4 * tau2) * 0.5;
						fprintf(fp, "% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t0.0\n",
								strain_xx[i][j][k][l], strain_yy[i][j][k][l], strain_xy[i][j][k][l],
								stress_rr, stress_theta);
					}
				}
			}
		}
	}
	else
	{
		fprintf(fp, "%d\t%d\t%d\t%d\n",
				division_n_xi, division_n_eta,
				element_n_xi, element_n_eta);
		for (i = 0; i < division_n_xi; i++)
		{
			for (j = 0; j < division_n_eta; j++)
			{
				fprintf(fp, "% 1.4e\t% 1.4e\n",
						coord_x[i][j], coord_y[i][j]);
			}
		}
	}
	fclose(fp);
	//グラフ用ファイル書き込み
	fp = fopen("over_disp_graph.txt", "a");
	for (i = 0; i < division_n_xi; i++)
	{
		for (j = 0; j < division_n_eta; j++)
		{
			temp1 = disp_x[i][j];
			temp2 = disp_y[i][j];
			// temp3 = sqrt(temp1 * temp1 + temp2 * temp2);
			fprintf(fp, "%d\t% 1.15e\t% 1.15e\t% 1.15e\t% 1.15e\n",
					graph_patch_n,
					coord_x[i][j], coord_y[i][j],
					disp_x[i][j], disp_y[i][j]);
		}
	}
	fclose(fp);

	fp = fopen("over_stress_x_graph.txt", "a");
	for (i = 0; i < element_n_xi; i++)
	{
		for (j = 0; j < element_n_eta; j++)
		{
			for (k = 0; k < division_ele_xi + 1; k++)
			{
				for (l = 0; l < division_ele_eta + 1; l++)
				{
					fprintf(fp, "% 1.10e\t% 1.10e\t% 1.10e\n",
							coord_x[i * division_ele_xi + k][j * division_ele_eta + l],
							coord_y[i * division_ele_xi + k][j * division_ele_eta + l],
							stress_xx[i][j][k][l]);
				}
			}
		}
	}
	fclose(fp);

	fp = fopen("over_stress_y_graph.txt", "a");
	for (i = 0; i < element_n_xi; i++)
	{
		for (j = 0; j < element_n_eta; j++)
		{
			for (k = 0; k < division_ele_xi + 1; k++)
			{
				for (l = 0; l < division_ele_eta + 1; l++)
				{
					fprintf(fp, "% 1.15e\t% 1.15e\t% 1.15e\n",
							coord_x[i * division_ele_xi + k][j * division_ele_eta + l],
							coord_y[i * division_ele_xi + k][j * division_ele_eta + l],
							stress_yy[i][j][k][l]);
				}
			}
		}
	}
	fclose(fp);

	fp = fopen("over_stress_y_graph_0.txt", "a");
	for (i = 0; i < element_n_xi; i++)
	{
		for (j = 0; j < element_n_eta; j++)
		{
			for (k = 0; k < division_ele_xi + 1; k++)
			{
				for (l = 0; l < division_ele_eta + 1; l++)
				{
					if (coord_y[i * division_ele_xi + k][j * division_ele_eta + l] == 0.000000000000000e+00)
					{
						fprintf(fp, "% 1.15e\t% 1.15e\t% 1.15e\n",
								coord_x[i * division_ele_xi + k][j * division_ele_eta + l],
								coord_y[i * division_ele_xi + k][j * division_ele_eta + l],
								stress_yy[i][j][k][l]);
					}
				}
			}
		}
	}
	fclose(fp);

	fp = fopen("over_stress_r_theta_graph.txt", "a");
	for (i = 0; i < element_n_xi; i++)
	{
		for (j = 0; j < element_n_eta; j++)
		{
			for (k = 0; k < division_ele_xi + 1; k++)
			{
				for (l = 0; l < division_ele_eta + 1; l++)
				{
					double stress_rr, stress_theta;
					double sum = stress_xx[i][j][k][l] + stress_yy[i][j][k][l];
					double dif = stress_xx[i][j][k][l] - stress_yy[i][j][k][l];
					double tau2 = stress_xy[i][j][k][l] * stress_xy[i][j][k][l];
					stress_rr = sum * 0.5 + sqrt(dif * dif + 4 * tau2) * 0.5;
					stress_theta = sum * 0.5 - sqrt(dif * dif + 4 * tau2) * 0.5;
					fprintf(fp, "% 1.10e\t% 1.10e\t% 1.10e\t% 1.10e\t% 1.10e\t% 1.10e\n",
							calc_xi_loc[i * division_ele_xi + k],
							calc_eta_loc[j * division_ele_eta + l],
							coord_x[i * division_ele_xi + k][j * division_ele_eta + l],
							coord_y[i * division_ele_xi + k][j * division_ele_eta + l],
							stress_rr, stress_theta);
				}
			}
		}
	}
	fclose(fp);

	fp = fopen("over_stress_vm_graph.txt", "a");
	for (i = 0; i < element_n_xi; i++)
	{
		for (j = 0; j < element_n_eta; j++)
		{
			for (k = 0; k < division_ele_xi + 1; k++)
			{
				for (l = 0; l < division_ele_eta + 1; l++)
				{
					double stress_vm;
					double sum = stress_xx[i][j][k][l] + stress_yy[i][j][k][l];
					double dif = stress_xx[i][j][k][l] - stress_yy[i][j][k][l];
					double tau2 = stress_xy[i][j][k][l] * stress_xy[i][j][k][l];
					double temp1, temp2;
					temp1 = 0.5 * sum;
					temp2 = 0.5 * sqrt(dif * dif + 4 * tau2);
					stress_vm = sqrt(temp1 * temp1 + 3 * temp2 * temp2);
					fprintf(fp, "% 1.10e\t% 1.10e\t% 1.10e\t% 1.10e\t% 1.10e\n",
							calc_xi_loc[i * division_ele_xi + k],
							calc_eta_loc[j * division_ele_eta + l],
							coord_x[i * division_ele_xi + k][j * division_ele_eta + l],
							coord_y[i * division_ele_xi + k][j * division_ele_eta + l],
							stress_vm);
				}
			}
		}
	}
	fclose(fp);
}

/*
static void Calculation_at_GP(double E, double nu)
{
	//通常IGAでのガウス点での値
	int i, j, k, e;

	Make_gauss_array(1);

	//メッシュ座標計算
	int ele_glo_n = real_Total_Element_on_mesh[0];// グローバルメッシュの要素数

	double U_temp[MAX_KIEL_SIZE];
	double B[D_MATRIX_SIZE][MAX_KIEL_SIZE];
	double J;

	for (e = 0; e < ele_glo_n; e++)
	{
		double X_temp[No_Control_point_ON_ELEMENT[Element_patch[e]]][DIMENSION];

		// printf("ele = %d\n\n", e);

		//strain_GPの初期化
		for (i = 0; i < GP_2D; i++)
		{
			for (j = 0; j < 3; j++)
			{
				strain_GP[e][i][j] = 0.0;
			}
		}

		// printf("x\ty\n");
		for (j = 0; j < No_Control_point_ON_ELEMENT[Element_patch[e]]; j++)
		{
			for (k = 0; k < DIMENSION; k++)
			{
				U_temp[j * DIMENSION + k] = Displacement[Controlpoint_of_Element[e * MAX_NO_CCpoint_ON_ELEMENT + j] * DIMENSION + k];
				X_temp[j][k] = Node_Coordinate[Controlpoint_of_Element[e * MAX_NO_CCpoint_ON_ELEMENT + j] * (DIMENSION + 1) + k];
			}
			// printf("%.15e\t%.15e\n", X_temp[j][0], X_temp[j][1]);
		}
		// printf("\n");

		for (i = 0; i < GP_2D; i++)	//ガウス点のループ
		{
			double data_result_shape[2] = {0.0, 0.0};
			double R_shape_func;

			// printf("Gxi_x\t_Gxi_y\tR_shape_func\n");
			for (j = 0; j < No_Control_point_ON_ELEMENT[Element_patch[e]]; j++)
			{
				R_shape_func = Shape_func(j, Total_Control_Point_to_mesh[Total_mesh], Gxi[i], e);
				for (k = 0; k < DIMENSION; k++)
				{
					data_result_shape[k] += R_shape_func * X_temp[j][k];
				}
			}

			//物理座標[要素番号(ローカル内で0から始まる)][ガウス点番号][DIMENSION]
			for (j = 0; j < DIMENSION; j++)
			{
				coordinate_GP[e][i][j] = data_result_shape[j];
			}
		}

		for (i = 0; i < GP_2D; i++)	//ガウス点のループ
		{
			Make_B_Matrix(e, B, Gxi[i], X_temp, &J, Total_Control_Point_to_mesh[Total_mesh]);
			for (j = 0; j < D_MATRIX_SIZE; j++)
			{
				for (k = 0; k < KIEL_SIZE; k++)
				{
					strain_GP[e][i][j] += B[j][k] * U_temp[k];
				}
			}
			Jac[e][i] = J;
		}
	}

	//Dマトリクスの計算
	double temp1;
	double D_matrix[3][3] = {{0.0}};
	if (DM == 0) { //平面応力状態
		temp1 = E * (1.0 - nu * nu);
		D_matrix[0][0] = temp1;
		D_matrix[0][1] = nu * temp1;
		D_matrix[1][0] = nu * temp1;
		D_matrix[1][1] = temp1;
		D_matrix[2][2] = (1.0 - nu) / 2.0 * temp1;
	} else if (DM == 1) { //平面ひずみ状態(2Dの場合はこっち)
		temp1 = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
		D_matrix[0][0] = temp1;
		D_matrix[0][1] = nu / (1.0 - nu) * temp1;
		D_matrix[1][0] = nu / (1.0 - nu) * temp1;
		D_matrix[1][1] = temp1;
		D_matrix[2][2] = (1.0 - 2.0 * nu) / 2.0 / (1.0 - nu) * temp1;
	}

	for (e = 0; e < ele_glo_n; e++)
	{
		for (i = 0; i < GP_2D; i++)
		{
			stress_GP[e][i][0] = D_matrix[0][0] * strain_GP[e][i][0] + D_matrix[0][1] * strain_GP[e][i][1];
			stress_GP[e][i][1] = D_matrix[1][0] * strain_GP[e][i][0] + D_matrix[1][1] * strain_GP[e][i][1];
			stress_GP[e][i][2] = D_matrix[2][2] * strain_GP[e][i][2];
		}
	}

	//座標変換
	double theta = 0.0;

	for (e = 0; e < ele_glo_n; e++)
	{
		for (i = 0; i < GP_2D; i++)
		{
			theta = atan2(coordinate_GP[e][i][1], coordinate_GP[e][i][0]);

			stress_r_theta_GP[e][i][0] = stress_GP[e][i][0] * pow(cos(theta), 2.0)
									   + stress_GP[e][i][1] * pow(sin(theta), 2.0)
									   + 2.0 * stress_GP[e][i][2] * sin(theta) * cos(theta);
			stress_r_theta_GP[e][i][1] = stress_GP[e][i][0] * pow(sin(theta), 2.0)
									   + stress_GP[e][i][1] * pow(cos(theta), 2.0)
									   - 2.0 * stress_GP[e][i][2] * cos(theta) * sin(theta);
			stress_r_theta_GP[e][i][2] = (stress_GP[e][i][1] - stress_GP[e][i][0])
									   * sin(theta) * cos(theta) + stress_GP[e][i][2]
									   * (pow(cos(theta), 2.0) - pow(sin(theta), 2.0));
		}
	}

	//厚肉円筒の理論解
	double r_t = 0.0;

	for (e = 0; e < ele_glo_n; e++)
	{
		for (i = 0; i < GP_2D; i++)
		{
			r_t = sqrt(pow(coordinate_GP[e][i][0], 2.0) + pow(coordinate_GP[e][i][1], 2.0));

			stress_theory_r_theta[e][i][0] = (pow(r_t, 2.0) - 4.0) / (pow(r_t, 2.0) * 3.0);
			stress_theory_r_theta[e][i][1] = (pow(r_t, 2.0) + 4.0) / (pow(r_t, 2.0) * 3.0);
		}
	}

	//書き込み
	fp = fopen("at_GP_overlay_data.txt", "w");
	fprintf(fp, "e\tガウス番号\tx\ty\tstress_xx\tstress_yy\tstress_r\tstress_theta\n");
	for (e = 0; e < ele_glo_n; e++)
	{
		for (i = 0; i < GP_2D; i++)
		{
			fprintf(fp, "%d\t%d\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n", e, i, coordinate_GP[e][i][0], coordinate_GP[e][i][1], stress_GP[e][i][0], stress_GP[e][i][1], stress_r_theta_GP[e][i][0], stress_r_theta_GP[e][i][1]);
		}
	}
	fclose(fp);

	fp = fopen("at_GP_overlay_for_error_norm.txt", "w");
	fprintf(fp, "e\tガウス番号\tx\ty\tstress_r-theory\tstress_theta-theory\n");
	for (e = 0; e < ele_glo_n; e++)
	{
		for (i = 0; i < GP_2D; i++)
		{
			fprintf(fp, "%d\t%d\t%.15e\t%.15e\t%.15e\t%.15e\n", e, i, coordinate_GP[e][i][0], coordinate_GP[e][i][1], stress_r_theta_GP[e][i][0] - stress_theory_r_theta[e][i][0], stress_r_theta_GP[e][i][1] - stress_theory_r_theta[e][i][1]);
		}
	}
	fclose(fp);


	//error normを計算
	//ガウス点で出したtheoryとの差の二乗を面積分
	double temp2 = 0.0, temp3 = 0.0, temp4 = 0.0, temp5 = 0.0;
	for (e = 0; e < ele_glo_n; e++)
	{
		for (i = 0; i < GP_2D; i++)
		{
			temp2 += w[i] * pow((stress_r_theta_GP[e][i][0] - stress_theory_r_theta[e][i][0]), 2.0) * Jac[e][i];
			temp3 += w[i] * pow((stress_r_theta_GP[e][i][1] - stress_theory_r_theta[e][i][1]), 2.0) * Jac[e][i];
			temp4 += w[i] * pow(stress_r_theta_GP[e][i][0], 2.0) * Jac[e][i];
			temp5 += w[i] * pow(stress_r_theta_GP[e][i][1], 2.0) * Jac[e][i];
		}
	}
	fp = fopen("at_GP_overlay_for_error_norm_surface_integral.txt", "w");
	fprintf(fp, "(stress_r-theory)^2_surface_integral\t(stress_theta-thory)^2_surface_integral\tstress_r^2_surface_integral\n");
	fprintf(fp, "%.15e\t%.15e\t%.15e\t%.15e\n", temp2, temp3, temp4, temp5);
	fclose(fp);
}


void Calculation_overlay_at_GP(double E, double nu,
							   int order_xi_glo, int order_eta_glo,
							   int knot_n_xi_glo, int knot_n_eta_glo,
							   int cntl_p_n_xi_glo, int cntl_p_n_eta_glo,
							   double *knot_vec_xi_glo, double *knot_vec_eta_glo,
							   double *cntl_px_glo, double *cntl_py_glo,
							   double *disp_cntl_px_glo, double *disp_cntl_py_glo,
							   double *weight_glo)
{
	//s-IGAでのローカル上のガウス点での重ね合わせた値
	int i, j, k, e;

	double disp_glo[DIMENSION];
	double strain_glo[D_MATRIX_SIZE];

	Make_gauss_array(1);

	double G_GP_knot[GP_2D][DIMENSION];

	//メッシュ座標計算
	int ele_glo_n = real_Total_Element_on_mesh[0];// グローバルメッシュの要素数
	int ele_loc_n = real_Total_Element_to_mesh[Total_mesh] - real_Total_Element_on_mesh[0];// ローカルメッシュの要素数

	double U_temp[MAX_KIEL_SIZE];
	double B[D_MATRIX_SIZE][MAX_KIEL_SIZE];
	double J;

	for (e = 0; e < ele_loc_n; e++)
	{
		int El_No_loc = ele_glo_n + e;
		double X_temp[No_Control_point_ON_ELEMENT[Element_patch[El_No_loc]]][DIMENSION];

		//strain_GPの初期化
		for (i = 0; i < GP_2D; i++)
		{
			for (j = 0; j < 3; j++)
			{
				strain_GP[e][i][j] = 0.0;
			}
		}

		for (j = 0; j < No_Control_point_ON_ELEMENT[Element_patch[El_No_loc]]; j++)
		{
			for (k = 0; k < DIMENSION; k++)
			{
				U_temp[j * DIMENSION + k] = Displacement[Controlpoint_of_Element[El_No_loc * MAX_NO_CCpoint_ON_ELEMENT + j] * DIMENSION + k];
				X_temp[j][k] = Node_Coordinate[Controlpoint_of_Element[El_No_loc * MAX_NO_CCpoint_ON_ELEMENT + j] * (DIMENSION + 1) + k];
			}
		}

		for (i = 0; i < GP_2D; i++)	//ガウス点のループ
		{
			double data_result_shape[2] = {0.0, 0.0};
			double R_shape_func;

			for (j = 0; j < No_Control_point_ON_ELEMENT[Element_patch[El_No_loc]]; j++)
			{
				R_shape_func = Shape_func(j, Total_Control_Point_to_mesh[Total_mesh], Gxi[i], El_No_loc);
				for (k = 0; k < DIMENSION; k++)
				{
					data_result_shape[k] += R_shape_func * X_temp[j][k];
				}
			}

			for (j = 0; j < DIMENSION; j++)
			{
				coordinate_GP[e][i][j] = data_result_shape[j];
			}
		}

		for (i = 0; i < GP_2D; i++)	//ガウス点のループ
		{
			Make_B_Matrix(El_No_loc, B, Gxi[i], X_temp, &J, Total_Control_Point_to_mesh[Total_mesh]);
			for (j = 0; j < D_MATRIX_SIZE; j++)
			{
				for (k = 0; k < KIEL_SIZE; k++)
				{
					strain_GP[e][i][j] += B[j][k] * U_temp[k];
				}
			}
			Jac[e][i] = J;
		}

		for (i = 0; i < GP_2D; i++)	//ガウス点のループ
		{
			CalcXiEtaByNR(coordinate_GP[e][i][0], coordinate_GP[e][i][1],
						  knot_vec_xi_glo, knot_vec_eta_glo,
						  cntl_px_glo, cntl_py_glo,
						  disp_cntl_px_glo, disp_cntl_py_glo,
						  cntl_p_n_xi_glo, cntl_p_n_eta_glo,
						  weight_glo, order_xi_glo, order_eta_glo,
						  &G_GP_knot[i][0], &G_GP_knot[i][1],
						  &disp_glo[0], &disp_glo[1],
						  &strain_glo[0], &strain_glo[1], &strain_glo[2]);

			//重ね合わせ
			for (j = 0; j < 3; j++) //xx, yy, xyを重ね合わせる
			{
				strain_GP[e][i][j] += strain_glo[j];
			}
		}
	}

	//Dマトリクスの計算
	double temp1;
	double D_matrix[3][3] = {{0.0}};
	if (DM == 0) { //平面応力状態
		temp1 = E * (1.0 - nu * nu);
		D_matrix[0][0] = temp1;
		D_matrix[0][1] = nu * temp1;
		D_matrix[1][0] = nu * temp1;
		D_matrix[1][1] = temp1;
		D_matrix[2][2] = (1.0 - nu) / 2.0 * temp1;
	} else if (DM == 1) { //平面ひずみ状態(2Dの場合はこっち)
		temp1 = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
		D_matrix[0][0] = temp1;
		D_matrix[0][1] = nu / (1.0 - nu) * temp1;
		D_matrix[1][0] = nu / (1.0 - nu) * temp1;
		D_matrix[1][1] = temp1;
		D_matrix[2][2] = (1.0 - 2.0 * nu) / 2.0 / (1.0 - nu) * temp1;
	}

	for (e = 0; e < ele_loc_n; e++)
	{
		for (i = 0; i < GP_2D; i++)
		{
			stress_GP[e][i][0] = D_matrix[0][0] * strain_GP[e][i][0] + D_matrix[0][1] * strain_GP[e][i][1];
			stress_GP[e][i][1] = D_matrix[1][0] * strain_GP[e][i][0] + D_matrix[1][1] * strain_GP[e][i][1];
			stress_GP[e][i][2] = D_matrix[2][2] * strain_GP[e][i][2];
		}
	}

	//座標変換
	double theta = 0.0;

	for (e = 0; e < ele_loc_n; e++)
	{
		for (i = 0; i < GP_2D; i++)
		{
			theta = atan2(coordinate_GP[e][i][1], coordinate_GP[e][i][0]);

			stress_r_theta_GP[e][i][0] = stress_GP[e][i][0] * pow(cos(theta), 2.0)
									   + stress_GP[e][i][1] * pow(sin(theta), 2.0)
									   + 2.0 * stress_GP[e][i][2] * sin(theta) * cos(theta);
			stress_r_theta_GP[e][i][1] = stress_GP[e][i][0] * pow(sin(theta), 2.0)
									   + stress_GP[e][i][1] * pow(cos(theta), 2.0)
									   - 2.0 * stress_GP[e][i][2] * cos(theta) * sin(theta);
			stress_r_theta_GP[e][i][2] = (stress_GP[e][i][1] - stress_GP[e][i][0])
									   * sin(theta) * cos(theta) + stress_GP[e][i][2]
									   * (pow(cos(theta), 2.0) - pow(sin(theta), 2.0));
		}
	}

	//円孔を有する無限平板の理論解
	double r_t = 0.0, theta_t = 0.0;

	for (e = 0; e < ele_loc_n; e++)
	{
		for (i = 0; i < GP_2D; i++)
		{
			r_t = sqrt(pow(coordinate_GP[e][i][0], 2.0) + pow(coordinate_GP[e][i][1], 2.0));
			theta_t = atan2(coordinate_GP[e][i][1], coordinate_GP[e][i][0]);

			// stress_theory_r_theta[e][i][0] = (10.0 / 2.0) * (1.0 - pow((1.0 / r_t), 2.0))
			// 							   + (10.0 / 2.0) * (1.0 + 3.0 * pow((1.0 / r_t), 4.0)
			// 							   - 4.0 * pow((1.0 / r_t), 2.0)) * cos(2.0 * theta_t);
			// stress_theory_r_theta[e][i][1] = (10.0 / 2.0) * (1.0 + pow((1.0 / r_t), 2.0))
			// 							   - (10.0 / 2.0) * (1.0 + 3.0 * pow((1.0 / r_t), 4.0))
			// 							   * cos(2.0 * theta_t);
			// stress_theory_r_theta[e][i][2] = - (10.0 / 2.0) * (1.0 - 3.0 * pow((1.0 / r_t), 4.0)
			// 							   + 2.0 * pow((1.0 / r_t), 2.0)) * sin(2.0 * theta_t);
			stress_theory_r_theta[e][i][0] = (10.0 / 2.0) * (1.0 - pow((1.0 / r_t), 2.0))
										   - (10.0 / 2.0) * (1.0 + 3.0 * pow((1.0 / r_t), 4.0)
										   - 4.0 * pow((1.0 / r_t), 2.0)) * cos(2.0 * theta_t);
			stress_theory_r_theta[e][i][1] = (10.0 / 2.0) * (1.0 + pow((1.0 / r_t), 2.0))
										   + (10.0 / 2.0) * (1.0 + 3.0 * pow((1.0 / r_t), 4.0))
										   * cos(2.0 * theta_t);
			stress_theory_r_theta[e][i][2] = (10.0 / 2.0) * (1.0 - 3.0 * pow((1.0 / r_t), 4.0)
										   + 2.0 * pow((1.0 / r_t), 2.0)) * sin(2.0 * theta_t);
		}
	}

	//書き込み
	fp = fopen("at_GP_overlay_data.txt", "w");
	fprintf(fp, "e\tガウス番号\tx\ty\tstress_xx\tstress_yy\tstress_r\tstress_theta\n");
	for (e = 0; e < ele_loc_n; e++)
	{
		for (i = 0; i < GP_2D; i++)
		{
			fprintf(fp, "%d\t%d\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n", e, i, coordinate_GP[e][i][0], coordinate_GP[e][i][1], stress_GP[e][i][0], stress_GP[e][i][1], stress_r_theta_GP[e][i][0], stress_r_theta_GP[e][i][1]);
		}
	}
	fclose(fp);

	fp = fopen("at_GP_overlay_for_error_norm.txt", "w");
	fprintf(fp, "e\tガウス番号\tx\ty\tstress_r-theory\tstress_theta-theory\n");
	for (e = 0; e < ele_loc_n; e++)
	{
		for (i = 0; i < GP_2D; i++)
		{
			fprintf(fp, "%d\t%d\t%.15e\t%.15e\t%.15e\t%.15e\n", e, i, coordinate_GP[e][i][0], coordinate_GP[e][i][1], stress_r_theta_GP[e][i][0] - stress_theory_r_theta[e][i][0], stress_r_theta_GP[e][i][1] - stress_theory_r_theta[e][i][1]);
		}
	}
	fclose(fp);


	//error normを計算
	//ガウス点で出したtheoryとの差の二乗を面積分
	double temp2 = 0.0, temp3 = 0.0, temp4 = 0.0, temp5 = 0.0, temp6 = 0.0;
	for (e = 0; e < ele_loc_n; e++)
	{
		for (i = 0; i < GP_2D; i++)
		{
			temp2 += w[i] * pow((stress_r_theta_GP[e][i][0] - stress_theory_r_theta[e][i][0]), 2.0) * Jac[e][i];
			temp3 += w[i] * pow((stress_r_theta_GP[e][i][1] - stress_theory_r_theta[e][i][1]), 2.0) * Jac[e][i];
			temp4 += w[i] * pow(stress_r_theta_GP[e][i][0], 2.0) * Jac[e][i];
			temp5 += w[i] * pow(stress_r_theta_GP[e][i][1], 2.0) * Jac[e][i];
			temp6 += w[i] * Jac[e][i];
		}
	}
	fp = fopen("at_GP_overlay_for_error_norm_surface_integral.txt", "w");
	fprintf(fp, "(stress_r-theory)^2_surface_integral\t(stress_theta-thory)^2_surface_integral\tstress_r^2_surface_integral\tstress_theta^2_surface_integral\t面積(analysis)\n");
	fprintf(fp, "%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n", temp2, temp3, temp4, temp5, temp6);
	fclose(fp);
}
*/

void K_output_svg(int ndof)
{
	// [K] = [[K^G, K^GL], [K^GL, K^L]]

	int i, j;

	char color_vec[2][10] = {"#f5f5f5", "#ee82ee"};
	// 0	whitesmoke
	// 1	violet
	// https://www.colordic.org/

	double space = 3.0, scale = 1000.0 / (((double)ndof) + 2.0 * space);

	double width = (((double)ndof) + 2.0 * space) * scale;
	double height = width;

	char str[256] = "K_matrix.svg";
	fp = fopen(str, "w");

	fprintf(fp, "<?xml version='1.0'?>\n");
	// fprintf(fp, "<svg width='%lept' height='%lept' viewBox='0 0 %le %le' style = 'background: #eee' xmlns='http://www.w3.org/2000/svg' xmlns:xlink='http://www.w3.org/1999/xlink'>\n", width, height, width, height);
	fprintf(fp, "<svg width='%le' height='%le' version='1.1' style='background: #eee' xmlns='http://www.w3.org/2000/svg' xmlns:xlink='http://www.w3.org/1999/xlink'>\n", width, height);

	double xx = space * scale;
	double yy = space * scale;
	double ww = ndof * scale;
	double hh = ndof * scale;
	fprintf(fp, "<rect x='%le' y='%le' width='%le' height='%le' fill='%s' />\n", xx, yy, ww, hh, color_vec[0]);

	// 各行の成分を抽出
	for (i = 0; i < ndof; i++)
	{
		int *K_bool = (int *)malloc(sizeof(int) * ndof); // 一行分保存する
		for (j = 0; j < ndof; j++)
		{
			K_bool[j] = 0;
		}

		for (j = 0; j < ndof; j++)
		{
			int temp_count;
			if (i <= j)
			{
				temp_count = RowCol_to_icount(i, j);
			}
			else if (i > j)
			{
				temp_count = RowCol_to_icount(j, i);
			}

			if (temp_count != -1)
			{
				K_bool[j] = 1;
			}
		}

		for (j = 0; j < ndof; j++)
		{
			double x = (((double)j) + space) * scale;
			double y = (((double)i) + space) * scale;
			// if (K_bool[j] == 0)
			// {
			// 	fprintf(fp, "<rect x='%le' y='%le' width='%le' height='%le' fill='%s' />\n", x, y, scale, scale, color_vec[0]);
			// }
			if (K_bool[j] == 1)
			{
				fprintf(fp, "<rect x='%le' y='%le' width='%le' height='%le' fill='%s' />\n", x, y, scale, scale, color_vec[1]);
			}
		}
		free(K_bool);
	}

	fprintf(fp, "</svg>");
	fclose(fp);
}

void Make_gauss_array(int select_GP)
{
	int i, j;

	if (select_GP == 0)
	{
		GP_1dir = Ng;
	}
	else if (select_GP == 1)
	{
		GP_1dir = Ng_extended;
	}

	GP_2D = GP_1dir * GP_1dir;

	if (GP_1dir == 3)
	{
		double G1 = pow((3.0 / 5.0), 0.5);
		double G_vec[3] = {-G1, 0.0, G1};
		double w1 = 8.0 / 9.0;
		double w2 = 5.0 / 9.0;
		double w_vec[3] = {w2, w1, w2};

		for (i = 0; i < GP_1dir; i++)
		{
			for (j = 0; j < GP_1dir; j++)
			{
				w[j + (GP_1dir * i)] = w_vec[i] * w_vec[j];
				Gxi[j + (GP_1dir * i)][0] = G_vec[j];
				Gxi[j + (GP_1dir * i)][1] = G_vec[i];
			}
		}
	}
	else if (GP_1dir == 4)
	{
		double A = pow((6.0 / 5.0), 0.5);
		double G1 = pow(((3.0 - 2.0 * A) / 7.0), 0.5);
		double G2 = pow(((3.0 + 2.0 * A) / 7.0), 0.5);
		double G_vec[4] = {-G2, -G1, G1, G2};
		double B = pow(30.0, 0.5);
		double w1 = (18.0 + B) / 36.0;
		double w2 = (18.0 - B) / 36.0;
		double w_vec[4] = {w2, w1, w1, w2};

		for (i = 0; i < GP_1dir; i++)
		{
			for (j = 0; j < GP_1dir; j++)
			{
				w[j + (GP_1dir * i)] = w_vec[i] * w_vec[j];
				Gxi[j + (GP_1dir * i)][0] = G_vec[j];
				Gxi[j + (GP_1dir * i)][1] = G_vec[i];
			}
		}
	}
	else if (GP_1dir == 5)
	{
		double A = pow((10.0 / 7.0), 0.5);
		double G1 = pow((5.0 - 2.0 * A), 0.5) / 3.0;
		double G2 = pow((5.0 + 2.0 * A), 0.5) / 3.0;
		double G_vec[5] = {-G2, -G1, 0.0, G1, G2};
		double B = pow(70.0, 0.5);
		double w1 = 128.0 / 225.0;
		double w2 = (322.0 + 13.0 * B) / 900.0;
		double w3 = (322.0 - 13.0 * B) / 900.0;
		double w_vec[5] = {w3, w2, w1, w2, w3};

		for (i = 0; i < GP_1dir; i++)
		{
			for (j = 0; j < GP_1dir; j++)
			{
				w[j + (GP_1dir * i)] = w_vec[i] * w_vec[j];
				Gxi[j + (GP_1dir * i)][0] = G_vec[j];
				Gxi[j + (GP_1dir * i)][1] = G_vec[i];
			}
		}
	}
	else if (GP_1dir == 10)
	{
		double G_vec[10];
		double w_vec[10];

		G_vec[0] = -0.9739065285171717;
		G_vec[1] = -0.8650633666889845;
		G_vec[2] = -0.6794095682990244;
		G_vec[3] = -0.4333953941292472;
		G_vec[4] = -0.1488743389816312;
		G_vec[5] = 0.1488743389816312;
		G_vec[6] = 0.4333953941292472;
		G_vec[7] = 0.6794095682990244;
		G_vec[8] = 0.8650633666889845;
		G_vec[9] = 0.9739065285171717;

		w_vec[0] = 0.0666713443086881;
		w_vec[1] = 0.1494513491505804;
		w_vec[2] = 0.2190863625159820;
		w_vec[3] = 0.2692667193099965;
		w_vec[4] = 0.2955242247147530;
		w_vec[5] = 0.2955242247147530;
		w_vec[6] = 0.2692667193099965;
		w_vec[7] = 0.2190863625159820;
		w_vec[8] = 0.1494513491505804;
		w_vec[9] = 0.0666713443086881;

		for (i = 0; i < GP_1dir; i++)
		{
			for (j = 0; j < GP_1dir; j++)
			{
				w[j + (GP_1dir * i)] = w_vec[i] * w_vec[j];
				Gxi[j + (GP_1dir * i)][0] = G_vec[j];
				Gxi[j + (GP_1dir * i)][1] = G_vec[i];
			}
		}
	}
}

// about J integral
