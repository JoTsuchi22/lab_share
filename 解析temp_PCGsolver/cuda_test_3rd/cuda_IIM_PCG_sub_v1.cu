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
	int *CP = (int *)malloc(sizeof(int) * temp_i * DIMENSION);
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
	free(CP);
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
			  int *Order, int *No_knot, int *Total_Knot_to_mesh, double *Node_Coordinate,
			  double *Position_Knots, int *No_Control_point_ON_ELEMENT)
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
					difference[Total_Knot_to_patch_dim[(l + Total_Patch_to_Now) * DIMENSION + j] + kkk - Order[(l + Total_Patch_to_Now) * DIMENSION + j]]
						= Position_Knots[Total_Knot_to_patch_dim[(l + Total_Patch_to_Now) * DIMENSION + j] + kkk + 1] - Position_Knots[Total_Knot_to_patch_dim[(l + Total_Patch_to_Now) * DIMENSION + j] + kkk];
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
			Setting_Dist_Load_2D(tm, iPatch, Total_Element_to_mesh[tm + 1], iCoord, val_Coord,
								 Range_Coord, type_load, Coeff_Dist_Load, Total_Knot_to_mesh,
						  		 Controlpoint_of_Element, Order, No_knot, Total_element_all_ID,
						  		 Total_Knot_to_patch_dim, Position_Knots, Equivalent_Nodal_Force,
						  		 Total_Element_on_mesh, Total_Element_to_mesh, Element_patch, ENC,
						  		 Node_Coordinate, INC, No_Control_point_ON_ELEMENT, Total_Control_Point_to_mesh,
						  		 No_Control_point);
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
void Check_coupled_Glo_Loc_element_for_Gauss(int mesh_n_over, int mesh_n_org, int *NNLOVER, int *NELOVER,
											 double *Gauss_Coordinate, double *Gauss_Coordinate_ex, double *Jac, double *Jac_ex, double *B_Matrix, double *B_Matrix_ex, double *Loc_parameter_on_Glo, double *Loc_parameter_on_Glo_ex,
											 int *real_Total_Element_to_mesh, double *Node_Coordinate, int *Total_Control_Point_to_mesh, int *Controlpoint_of_Element,
											 int *INC, int *Element_patch, int *Order, int *No_Control_point_ON_ELEMENT, double *Position_Knots,
											 int *Total_Knot_to_patch_dim, int *No_knot, int *No_Control_point,
											 double *Control_Coord_x, double *Control_Coord_y, double *Control_Weight,
											 int *real_Total_Element_on_mesh, int *real_element, int *Total_Patch_on_mesh, int *line_No_Total_element)
{
	int re, e;
	int i, j, k, m;
	int b, l, ll;
	int n_elements_over_point[POW_Ng_extended];
	int MAX_NNLOVER = 0;

	int *temp_element_n = (int *)malloc(sizeof(int) * MAX_N_ELEMENT_OVER_POINT);
	int *element_n_point = (int *)malloc(sizeof(int) * MAX_N_ELEMENT_OVER_ELEMENT);
	int *Check_coupled_No = (int *)malloc(sizeof(int) * MAX_N_ELEMENT_OVER);

	for (i = 0; i < MAX_N_ELEMENT_OVER; i++)
	{
		Check_coupled_No[i] = 0;
	}

	for (m = 0; m < 2; m++) // 最初 Ng 個のガウス点で重なりを求め, NNLOVER[e] >= 2 の e に対して, 再度 Ng_extended 個のガウス点で重なりを求める
	{
		Make_gauss_array(m);

		// グローバルパッチの Preprocessing 作成
		if (m == 0)
		{
			double AAA = 0.0;
			for (re = 0; re < real_Total_Element_on_mesh[mesh_n_org]; re++)
			{
				e = real_element[re + real_Total_Element_to_mesh[mesh_n_org]];
				Preprocessing(m, e, Gauss_Coordinate, Gauss_Coordinate_ex, B_Matrix, B_Matrix_ex,
							  Jac, Jac_ex, real_Total_Element_to_mesh, Node_Coordinate, Total_Control_Point_to_mesh,
							  Controlpoint_of_Element, INC, Element_patch, Order, No_Control_point_ON_ELEMENT,
							  Position_Knots, Total_Knot_to_patch_dim, No_knot, No_Control_point);
			}
		}

		// ローカルパッチ(mesh_n_over)各要素の頂点の物理座標算出
		for (re = 0; re < real_Total_Element_on_mesh[mesh_n_over]; re++)
		{
			int i_gg, i_ee;
			double output_para[DIMENSION];
			int Total_n_elements;

			e = real_element[re + real_Total_Element_to_mesh[mesh_n_over]];

			if (m == 0 || (m == 1 && NNLOVER[e] >= 2))
			{
				Preprocessing(m, e, Gauss_Coordinate, Gauss_Coordinate_ex, B_Matrix, B_Matrix_ex,
							  Jac, Jac_ex, real_Total_Element_to_mesh, Node_Coordinate, Total_Control_Point_to_mesh,
							  Controlpoint_of_Element, INC, Element_patch, Order, No_Control_point_ON_ELEMENT,
							  Position_Knots, Total_Knot_to_patch_dim, No_knot, No_Control_point);

				if (m == 1)
				{
					NNLOVER[e] = 0;
					for (i = 0; i < NNLOVER[e]; i++)
					{
						NELOVER[e * MAX_N_ELEMENT_OVER + i] = 0;
					}
				}

				Total_n_elements = 0;
				k = 0;
				ll = 0;

				// ローカルパッチ各要素のガウス点の物理座標のグローバルパッチでの(xi, eta)算出
				for (i = 0; i < Total_Patch_on_mesh[mesh_n_org]; i++)
				{
					// グローバルパッチ i での各方向ノットベクトル
					double *temp_Position_Knots_xi = (double *)malloc(sizeof(double) * No_knot[i * DIMENSION + 0]);
					double *temp_Position_Knots_eta = (double *)malloc(sizeof(double) * No_knot[i * DIMENSION + 1]);

					for (j = 0; j < No_knot[i * DIMENSION + 0]; j++)
					{
						temp_Position_Knots_xi[j] = Position_Knots[Total_Knot_to_patch_dim[i * DIMENSION + 0] + j];
					}
					for (j = 0; j < No_knot[i * DIMENSION + 1]; j++)
					{
						temp_Position_Knots_eta[j] = Position_Knots[Total_Knot_to_patch_dim[i * DIMENSION + 1] + j];
					}

					for (i_ee = 0; i_ee < GP_1dir; i_ee++)
					{
						for (i_gg = 0; i_gg < GP_1dir; i_gg++)
						{
							int g_n = i_ee * GP_1dir + i_gg;
							double data_result_shape[DIMENSION] = {0.0};

							for (l = 0; l < DIMENSION; l++)
							{
								if (m == 0)
								{
									data_result_shape[l] = Gauss_Coordinate[e * GP_2D * DIMENSION + g_n * DIMENSION + l];
								}
								else if (m == 1)
								{
									data_result_shape[l] = Gauss_Coordinate_ex[e * GP_2D * DIMENSION + g_n * DIMENSION + l];
								}
							}
							int itr_n = Calc_xi_eta(data_result_shape[0], data_result_shape[1],
													temp_Position_Knots_xi, temp_Position_Knots_eta,
													No_Control_point[i * DIMENSION + 0], No_Control_point[i * DIMENSION + 1], Order[i * DIMENSION + 0], Order[i * DIMENSION + 1],
													&output_para[0], &output_para[1],
													Position_Knots, Total_Knot_to_patch_dim,
													No_Control_point, Order, Control_Coord_x, Control_Coord_y, Control_Weight, No_knot);
							
							if (m == 0)
							{
								Loc_parameter_on_Glo[e * GP_2D * DIMENSION + g_n * DIMENSION + 0] = output_para[0];
								Loc_parameter_on_Glo[e * GP_2D * DIMENSION + g_n * DIMENSION + 1] = output_para[1];
							}
							else if (m == 1)
							{
								Loc_parameter_on_Glo_ex[e * GP_2D * DIMENSION + g_n * DIMENSION + 0] = output_para[0];
								Loc_parameter_on_Glo_ex[e * GP_2D * DIMENSION + g_n * DIMENSION + 1] = output_para[1];
							}

							// Newton Laphsonによって出力されたxi,etaから重なる要素を求める
							n_elements_over_point[k] = ele_check(i, output_para, No_Control_point, Position_Knots, Total_Knot_to_patch_dim, No_knot, Order, temp_element_n, line_No_Total_element);
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
					free(temp_Position_Knots_xi), free(temp_Position_Knots_eta);
				}
				// 昇順ソート
				sort(Total_n_elements, element_n_point);
				// 重複削除
				NNLOVER[e] = duplicate_delete(Total_n_elements, e, NELOVER, element_n_point); // NNLOVER: 要素 e に重なる要素の総数
			}
		}
	}

	for (re = 0; re < real_Total_Element_on_mesh[mesh_n_over]; re++)
	{
		e = real_element[re + real_Total_Element_to_mesh[mesh_n_over]];

		Check_coupled_No[NNLOVER[e]]++;

		if (MAX_NNLOVER < NNLOVER[e])
		{
			MAX_NNLOVER = NNLOVER[e];
		}
	}

	// 重なっている要素の割合
	printf("MAX_NNLOVER = %d\n", MAX_NNLOVER);
	for (i = 0; i <= MAX_NNLOVER; i++)
	{
		double Percent_Check_coupled_No = (double)Check_coupled_No[i] * 100.0 / (double)real_Total_Element_on_mesh[mesh_n_over];
		printf("Check_coupled_No[%d] = %d\t\t%.2lf %%\n", i, Check_coupled_No[i], Percent_Check_coupled_No);
	}
	
	free(element_n_point), free(temp_element_n), free(Check_coupled_No);
}


void Make_Loc_Glo(int *real_Total_Element_on_mesh, int *real_Total_Element_to_mesh, int *real_element, int *NNLOVER, int *NELOVER)
{
	int i, j, k;
	int jj;
	int e;
	int count;

	int j_n = real_Total_Element_to_mesh[Total_mesh] - real_Total_Element_on_mesh[0];

	for (i = 0; i < real_Total_Element_on_mesh[0]; i++)
	{
		e = real_element[i];
		count = 0;

		for (j = 0; j < j_n; j++)
		{
			jj = real_element[real_Total_Element_to_mesh[1] + j]; //ローカルメッシュ上のreal element番号

			if (NNLOVER[jj] > 0)
			{
				for (k = 0; k < NNLOVER[jj]; k++)
				{
					if (NELOVER[jj * MAX_N_ELEMENT_OVER + k] == e)
					{
						NELOVER[e * MAX_N_ELEMENT_OVER + count] = jj;
						count++;
					}
				}
			}
		}
		NNLOVER[e] = count;
	}
}


// Newton Raphsonによって出力されたxi,etaから重なる要素を求める
int ele_check(int patch_n, double para_coord[DIMENSION], int *No_Control_point, double *Position_Knots, int *Total_Knot_to_patch_dim, int *No_knot, int *Order, int *temp_element_n, int *line_No_Total_element)
{
	int i, j, k, l;
	int kk, ll;
	int RangeCheck_flag;					// 要素を求め終えたら立てるフラグ
	int temp_ad[DIMENSION][MAX_ORDER + 1];	// 要素の位置を求めるための値
	int No_line[DIMENSION];					// xi, etaが含まれている要素の列数
	int n = 1;

	for (j = 0; j < DIMENSION; j++)
	{
		// 初期化
		l = 0;
		No_line[j] = 0;
		for (i = 0; i < MAX_ORDER + 1; i++)
		{
			temp_ad[j][i] = 0;
		}
		RangeCheck_flag = 0;

		for (k = 0; k < No_Control_point[patch_n * DIMENSION + j] - 1; k++)
		{
			if (RangeCheck_flag == 1)
				break;

			// Local要素の頂点がGlobalパッチ内にない場合
			if (para_coord[j] < Position_Knots[Total_Knot_to_patch_dim[patch_n * DIMENSION + j] + 0] || para_coord[j] > Position_Knots[Total_Knot_to_patch_dim[patch_n * DIMENSION + j] + No_knot[patch_n * DIMENSION + j] - 1])
			{
				RangeCheck_flag++;
			}

			// Local要素の頂点がGlobal要素内部にある場合
			if (para_coord[j] < Position_Knots[Total_Knot_to_patch_dim[patch_n * DIMENSION + j] + Order[patch_n * DIMENSION + j] + k])
			{
				int kk = 0;
				for (kk = 0; kk < k + 1; kk++)
				{
					if (para_coord[j] > Position_Knots[Total_Knot_to_patch_dim[patch_n * DIMENSION + j] + Order[patch_n * DIMENSION + j] + k - kk])
					{
						temp_ad[j][l] = k - kk;
						l++;
						RangeCheck_flag++;
						break;
					}
				}
			}

			// Local要素の頂点がGlobal要素境界上にある場合
			if (para_coord[j] == Position_Knots[Total_Knot_to_patch_dim[patch_n * DIMENSION + j] + Order[patch_n * DIMENSION + j] + k])
			{
				//頂点の座標がGlobalパッチの始点上にある場合
				if (para_coord[j] == Position_Knots[Total_Knot_to_patch_dim[patch_n * DIMENSION + j] + 0])
				{
					temp_ad[j][l] = k;
					l++;
					break;
				}
				//頂点の座標がGlobalパッチの終点上にある場合
				if (para_coord[j] == Position_Knots[Total_Knot_to_patch_dim[patch_n * DIMENSION + j] + No_knot[patch_n * DIMENSION + j] - 1])
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
				for (kk = 0; kk < Order[patch_n * DIMENSION + j]; kk++)
				{
					if (para_coord[j] == Position_Knots[Total_Knot_to_patch_dim[patch_n * DIMENSION + j] + Order[patch_n * DIMENSION + j] + k + kk + 1])
					// 多重ノット(次数分ループ)
					{
						printf("C0 continuity\n");
						temp_ad[j][l] = k + kk;
						l++;
					}
					if (para_coord[j] != Position_Knots[Total_Knot_to_patch_dim[patch_n * DIMENSION + j] + Order[patch_n * DIMENSION + j] + k + kk + 1])
						break;
				}
				RangeCheck_flag++;
			}
		}
		No_line[j] = l;
		
		// 各方向のNo_lineを掛け合わせる
		n *= l;
	}

	if (DIMENSION == 2)
	{
		for (l = 0; l < No_line[1]; l++)
		{
			for (ll = 0; ll < No_line[0]; ll++)
			{
				temp_element_n[l * No_line[0] + ll] = temp_ad[0][ll] + temp_ad[1][l] * line_No_Total_element[patch_n * DIMENSION + 0];
			}
		}
	}
	return n;
}


//昇順ソート
void sort(int total, int *element_n_point)
{
	int i, j;
	int temp;

	for (i = 0; i < total; i++)
	{
		for (j = i + 1; j < total; j++)
		{
			if (element_n_point[i] > element_n_point[j])
			{
				temp = element_n_point[i];
				element_n_point[i] = element_n_point[j];
				element_n_point[j] = temp;
			}
		}
	}
}


//重複削除
int duplicate_delete(int total, int element_n, int *NELOVER, int *element_n_point)
{
	int i, j;

	j = 0;
	NELOVER[element_n * MAX_N_ELEMENT_OVER + j] = element_n_point[0];
	j++;
	for (i = 1; i < total; i++)
	{
		if (element_n_point[i] != element_n_point[i - 1])
		{
			NELOVER[element_n * MAX_N_ELEMENT_OVER + j] = element_n_point[i];
			j++;
		}
	}
	// j = 要素element_nに重なる要素の総数
	return j;
}


// Preprocessing
void Preprocessing(int m, int e, double *Gauss_Coordinate, double *Gauss_Coordinate_ex, double *B_Matrix, double *B_Matrix_ex,
				   double *Jac, double *Jac_ex, int *real_Total_Element_to_mesh, double *Node_Coordinate, int *Total_Control_Point_to_mesh,
				   int *Controlpoint_of_Element, int *INC, int *Element_patch, int *Order, int *No_Control_point_ON_ELEMENT,
				   double *Position_Knots, int *Total_Knot_to_patch_dim, int *No_knot, int *No_Control_point)
{
	double *a_matrix = (double *)malloc(sizeof(double) * POW_Ng_extended * DIMENSION * DIMENSION); // a_matrix[POW_Ng_extended * DIMENSION * DIMENSION]
	double *dSF = (double *)malloc(sizeof(double) * POW_Ng * MAX_NO_CCpoint_ON_ELEMENT * DIMENSION); // dSF[POW_Ng * MAX_NO_CCpoint_ON_ELEMENT * DIMENSION]
	double *dSF_ex = (double *)malloc(sizeof(double) * POW_Ng_extended * MAX_NO_CCpoint_ON_ELEMENT * DIMENSION); // dSF_ex[POW_Ng_extended * MAX_NO_CCpoint_ON_ELEMENT * DIMENSION]

	// ガウス点の物理座標を計算
	Make_Gauss_Coordinate(m, e, Gauss_Coordinate, Gauss_Coordinate_ex, Node_Coordinate, Total_Control_Point_to_mesh,
						  Controlpoint_of_Element, INC, Element_patch, Order, No_Control_point_ON_ELEMENT,
						  Position_Knots, Total_Knot_to_patch_dim, No_knot, No_Control_point);

	// ガウス点でのヤコビアン, Bマトリックスを計算
	Make_dSF(m, e, dSF, dSF_ex, Element_patch, No_Control_point_ON_ELEMENT,
			 Controlpoint_of_Element, Total_Control_Point_to_mesh, Node_Coordinate, INC, Order,
			 Position_Knots, Total_Knot_to_patch_dim, No_knot, No_Control_point);
	Make_Jac(m, e, Jac, Jac_ex, dSF, dSF_ex, a_matrix, Node_Coordinate, Controlpoint_of_Element, No_Control_point_ON_ELEMENT, Element_patch);
	Make_B_Matrix(m, e, B_Matrix, B_Matrix_ex, dSF, dSF_ex, a_matrix, No_Control_point_ON_ELEMENT, Element_patch);

	free(a_matrix), free(dSF), free(dSF_ex);
}


void Make_Gauss_Coordinate(int m, int e, double *Gauss_Coordinate, double *Gauss_Coordinate_ex, double *Node_Coordinate, int *Total_Control_Point_to_mesh,
						   int *Controlpoint_of_Element, int *INC, int *Element_patch, int *Order, int *No_Control_point_ON_ELEMENT,
						   double *Position_Knots, int *Total_Knot_to_patch_dim, int *No_knot, int *No_Control_point)
{
	int i, j, k;
	double temp_coordinate[DIMENSION];
	double R;

	for (i = 0; i < GP_2D; i++)
	{
		temp_coordinate[0] = Gxi[i][0];
		temp_coordinate[1] = Gxi[i][1];
		
		for (j = 0; j < No_Control_point_ON_ELEMENT[Element_patch[e]]; j++)
		{
			R = Shape_func(j, temp_coordinate, e, Node_Coordinate, Total_Control_Point_to_mesh,
						   Controlpoint_of_Element, INC, Element_patch, Order, No_Control_point_ON_ELEMENT,
						   Position_Knots, Total_Knot_to_patch_dim, No_knot, No_Control_point);

			for (k = 0; k < DIMENSION; k++)
			{
				if (m == 0)
				{
					Gauss_Coordinate[e * GP_2D * DIMENSION + i * DIMENSION + k] += R * Node_Coordinate[Controlpoint_of_Element[e * MAX_NO_CCpoint_ON_ELEMENT + j] * (DIMENSION + 1) + k];
				}
				else if (m == 1)
				{
					Gauss_Coordinate_ex[e * GP_2D * DIMENSION + i * DIMENSION + k] += R * Node_Coordinate[Controlpoint_of_Element[e * MAX_NO_CCpoint_ON_ELEMENT + j] * (DIMENSION + 1) + k];
				}
			}
		}
	}
}


void Make_dSF(int m, int e, double *dSF, double *dSF_ex, int *Element_patch, int *No_Control_point_ON_ELEMENT,
			  int *Controlpoint_of_Element, int *Total_Control_Point_to_mesh, double *Node_Coordinate, int *INC, int *Order,
			  double *Position_Knots, int *Total_Knot_to_patch_dim, int *No_knot, int *No_Control_point)
{
	int i, j, k;
	double temp_coordinate[DIMENSION];

	for (i = 0; i < GP_2D; i++)
	{
		temp_coordinate[0] = Gxi[i][0];
		temp_coordinate[1] = Gxi[i][1];

		for (j = 0; j < No_Control_point_ON_ELEMENT[Element_patch[e]]; j++)
		{
			for (k = 0; k < DIMENSION; k++)
			{
				if (m == 0)
				{
					dSF[i * MAX_NO_CCpoint_ON_ELEMENT * DIMENSION + j * DIMENSION + k] = dShape_func(j, k, temp_coordinate, e, Controlpoint_of_Element, No_Control_point_ON_ELEMENT,
				   																					 Total_Control_Point_to_mesh, Element_patch, Node_Coordinate, INC, Order, Position_Knots,
				   																					 Total_Knot_to_patch_dim, No_knot, No_Control_point);
				}
				else if (m == 1)
				{
					dSF_ex[i * MAX_NO_CCpoint_ON_ELEMENT * DIMENSION + j * DIMENSION + k] = dShape_func(j, k, temp_coordinate, e, Controlpoint_of_Element, No_Control_point_ON_ELEMENT,
				   																						Total_Control_Point_to_mesh, Element_patch, Node_Coordinate, INC, Order, Position_Knots,
				   																						Total_Knot_to_patch_dim, No_knot, No_Control_point);
				}
			}
		}
	}
}


void Make_Jac(int m, int e, double *Jac, double *Jac_ex, double *dSF, double *dSF_ex, double *a_matrix, double *Node_Coordinate, int *Controlpoint_of_Element, int *No_Control_point_ON_ELEMENT, int *Element_patch)
{
	int i, j, k, l;
	double J;
	double a[DIMENSION][DIMENSION];

	for (i = 0; i < GP_2D; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			for (k = 0; k < DIMENSION; k++)
			{
				a[j][k] = 0.0;
				for (l = 0; l < No_Control_point_ON_ELEMENT[Element_patch[e]]; l++)
				{
					if (m == 0)
					{
						a[j][k] += dSF[i * MAX_NO_CCpoint_ON_ELEMENT * DIMENSION + l * DIMENSION + k] * Node_Coordinate[Controlpoint_of_Element[e * MAX_NO_CCpoint_ON_ELEMENT + l] * (DIMENSION + 1) + j];
					}
					else if (m == 1)
					{
						a[j][k] += dSF_ex[i * MAX_NO_CCpoint_ON_ELEMENT * DIMENSION + l * DIMENSION + k] * Node_Coordinate[Controlpoint_of_Element[e * MAX_NO_CCpoint_ON_ELEMENT + l] * (DIMENSION + 1) + j];
					}
				}
			}
		}

		if (DIMENSION == 2)
		{
			J = InverseMatrix_2x2(a);
			
		}
		else if (DIMENSION == 3)
		{
			J = InverseMatrix_3x3(a);
		}

		if (m == 0)
		{
			Jac[e * GP_2D + i] = J;
		}
		else if (m == 1)
		{
			Jac_ex[e * GP_2D + i] = J;
		}

		for (j = 0; j < DIMENSION; j++)
		{
			for (k = 0; k < DIMENSION; k++)
			{
				a_matrix[i * DIMENSION * DIMENSION + j * DIMENSION + k] = a[j][k];
			}
		}

		if (J <= 0)
		{
			printf("Error, J <= 0\n");
			exit(1);
		}
	}
}


void Make_B_Matrix(int m, int e, double *B_Matrix, double *B_Matrix_ex, double *dSF, double *dSF_ex, double *a_matrix, int *No_Control_point_ON_ELEMENT, int *Element_patch)
{
	int i, j, k, l;
	double *b = (double *)malloc(sizeof(double) * DIMENSION * MAX_NO_CCpoint_ON_ELEMENT);	// b[DIMENSION][MAX_NO_CCpoint_ON_ELEMENT]

	for (i = 0; i < GP_2D; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			for (k = 0; k < No_Control_point_ON_ELEMENT[Element_patch[e]]; k++)
			{
				b[j * MAX_NO_CCpoint_ON_ELEMENT + k] = 0.0;
				for (l = 0; l < DIMENSION; l++)
				{
					if (m == 0)
					{
						b[j * MAX_NO_CCpoint_ON_ELEMENT + k] += a_matrix[i * DIMENSION * DIMENSION + l * DIMENSION + j] * dSF[i * MAX_NO_CCpoint_ON_ELEMENT * DIMENSION + k * DIMENSION + l];
					}
					else if (m == 1)
					{
						b[j * MAX_NO_CCpoint_ON_ELEMENT + k] += a_matrix[i * DIMENSION * DIMENSION + l * DIMENSION + j] * dSF_ex[i * MAX_NO_CCpoint_ON_ELEMENT * DIMENSION + k * DIMENSION + l];
					}
				}
			}
		}

		// 2次元
		if (m == 0)
		{
			for (j = 0; j < No_Control_point_ON_ELEMENT[Element_patch[e]]; j++)
			{
				B_Matrix[e * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 0 * MAX_KIEL_SIZE + (2 * j)] = b[0 * MAX_NO_CCpoint_ON_ELEMENT + j];
				B_Matrix[e * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 0 * MAX_KIEL_SIZE + (2 * j + 1)] = 0.0;
				B_Matrix[e * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 1 * MAX_KIEL_SIZE + (2 * j)] = 0.0;
				B_Matrix[e * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 1 * MAX_KIEL_SIZE + (2 * j + 1)] = b[1 * MAX_NO_CCpoint_ON_ELEMENT + j];
				B_Matrix[e * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 2 * MAX_KIEL_SIZE + (2 * j)] = b[1 * MAX_NO_CCpoint_ON_ELEMENT + j];
				B_Matrix[e * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 2 * MAX_KIEL_SIZE + (2 * j + 1)] = b[0 * MAX_NO_CCpoint_ON_ELEMENT + j];
			}
		}
		else if (m == 1)
		{
			for (j = 0; j < No_Control_point_ON_ELEMENT[Element_patch[e]]; j++)
			{
				B_Matrix_ex[e * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 0 * MAX_KIEL_SIZE + (2 * j)] = b[0 * MAX_NO_CCpoint_ON_ELEMENT + j];
				B_Matrix_ex[e * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 0 * MAX_KIEL_SIZE + (2 * j + 1)] = 0.0;
				B_Matrix_ex[e * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 1 * MAX_KIEL_SIZE + (2 * j)] = 0.0;
				B_Matrix_ex[e * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 1 * MAX_KIEL_SIZE + (2 * j + 1)] = b[1 * MAX_NO_CCpoint_ON_ELEMENT + j];
				B_Matrix_ex[e * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 2 * MAX_KIEL_SIZE + (2 * j)] = b[1 * MAX_NO_CCpoint_ON_ELEMENT + j];
				B_Matrix_ex[e * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 2 * MAX_KIEL_SIZE + (2 * j + 1)] = b[0 * MAX_NO_CCpoint_ON_ELEMENT + j];
			}
		}
	}
	free(b);
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


// K matrix
// 応力歪マトリックス
void Make_D_Matrix(double *D)
{
	int i, j;

	// 2次元
	if (DIMENSION == 2)
	{
		if (DM == 0) // 平面応力状態
		{
			double Eone = E / (1.0 - nu * nu);
			double D1[3][3] = {{Eone, nu * Eone, 0}, {nu * Eone, Eone, 0}, {0, 0, (1 - nu) / 2 * Eone}};

			for (i = 0; i < D_MATRIX_SIZE; i++)
				for (j = 0; j < D_MATRIX_SIZE; j++)
					D[i * D_MATRIX_SIZE + j] = D1[i][j];
		}
		else if (DM == 1) // 平面ひずみ状態(2Dの場合はこっち)
		{
			double Eone = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
			double D1[3][3] = {{Eone, nu / (1.0 - nu) * Eone, 0}, {nu / (1.0 - nu) * Eone, Eone, 0}, {0, 0, (1 - 2 * nu) / 2 / (1.0 - nu) * Eone}};

			for (i = 0; i < D_MATRIX_SIZE; i++)
				for (j = 0; j < D_MATRIX_SIZE; j++)
					D[i * D_MATRIX_SIZE + j] = D1[i][j];
		}
	}
}


// 拘束されている行数を省いた行列の番号の制作
void Make_Index_Dof(int *Total_Control_Point_to_mesh, int *Total_Constraint_to_mesh, int *Constraint_Node_Dir, int *Index_Dof)
{
	int i, k = 0;

	// 拘束されている自由度(Degree Of free)をERRORにする
	for (i = 0; i < Total_Constraint_to_mesh[Total_mesh]; i++)
	{
		Index_Dof[Constraint_Node_Dir[i * 2 + 0] * DIMENSION + Constraint_Node_Dir[i * 2 + 1]] = ERROR;
	}
	// ERROR以外に番号を付ける
	for (i = 0; i < Total_Control_Point_to_mesh[Total_mesh] * DIMENSION; i++)
	{
		if (Index_Dof[i] != ERROR)
		{
			Index_Dof[i] = k;
			k++;
		}
	}

	K_Whole_Size = k;
	printf("K_Whole_Size = %d\n", k);
}


void Make_K_Whole_Ptr_Col(int *Total_Element_to_mesh, int *Total_Control_Point_to_mesh, int *Total_Control_Point_To_Node,
						  int *No_Control_point_ON_ELEMENT, int *Element_patch, int *Controlpoint_of_Element, int *Node_To_Node, int *NNLOVER,
						  int *NELOVER, int *Index_Dof, int *K_Whole_Ptr, int *K_Whole_Col)
{
	int i, j, k, ii, jj;
	int N, NE, i_index, j_index;

	// 初期化
	for (i = 0; i < K_Whole_Size + 1; i++)
		K_Whole_Ptr[i] = 0;

	// 大きく分割するためのループ
	for (N = 0; N < Total_Control_Point_to_mesh[Total_mesh]; N += K_DIVISION_LENGE)
	{
		// 各節点に接する節点を取得
		for (i = 0; i < K_DIVISION_LENGE; i++)
		{
			Total_Control_Point_To_Node[i] = 0;
		}
		for (i = 0; i < Total_Element_to_mesh[Total_mesh]; i++)
		{
			for (ii = 0; ii < No_Control_point_ON_ELEMENT[Element_patch[i]]; ii++)
			{
				NE = Controlpoint_of_Element[i * MAX_NO_CCpoint_ON_ELEMENT + ii] - N;
				if (0 <= NE && NE < K_DIVISION_LENGE)
				{
					for (j = 0; j < No_Control_point_ON_ELEMENT[Element_patch[i]]; j++) //ローカル要素
					{
						// 数字がない時
						if (Total_Control_Point_To_Node[NE] == 0)
						{
							// 節点番号を取得
							Node_To_Node[NE * Total_Control_Point_to_mesh[Total_mesh] + 0] = Controlpoint_of_Element[i * MAX_NO_CCpoint_ON_ELEMENT + j];
							Total_Control_Point_To_Node[NE]++;
						}
						// 同じものがあったら
						// k > 0 以降の取得
						// kのカウント
						for (k = 0; k < Total_Control_Point_To_Node[NE]; k++)
						{
							if (Node_To_Node[NE * Total_Control_Point_to_mesh[Total_mesh] + k] == Controlpoint_of_Element[i * MAX_NO_CCpoint_ON_ELEMENT + j])
							{
								break;
							}
						}
						// 未設定のNode_To_Node取得
						if (k == Total_Control_Point_To_Node[NE])
						{
							Node_To_Node[NE * Total_Control_Point_to_mesh[Total_mesh] + k] = Controlpoint_of_Element[i * MAX_NO_CCpoint_ON_ELEMENT + j];
							Total_Control_Point_To_Node[NE]++;
						}
					}
					// 別メッシュとの重なりを考慮
					if (NNLOVER[i] > 0)
					{
						for (jj = 0; jj < NNLOVER[i]; jj++)
						{
							for (j = 0; j < No_Control_point_ON_ELEMENT[Element_patch[NELOVER[i * MAX_N_ELEMENT_OVER + jj]]]; j++) //ローカル要素
							{
								// 数字がない時
								if (Total_Control_Point_To_Node[NE] == 0)
								{
									// 節点番号を取得
									Node_To_Node[NE * Total_Control_Point_to_mesh[Total_mesh] + 0] = Controlpoint_of_Element[NELOVER[i * MAX_N_ELEMENT_OVER + jj] * MAX_NO_CCpoint_ON_ELEMENT + j];
									Total_Control_Point_To_Node[NE]++;
								}

								// 同じものがあったら
								// k > 0 以降の取得
								// kのカウント
								for (k = 0; k < Total_Control_Point_To_Node[NE]; k++)
								{
									if (Node_To_Node[NE * Total_Control_Point_to_mesh[Total_mesh] + k] == Controlpoint_of_Element[NELOVER[i * MAX_N_ELEMENT_OVER + jj] * MAX_NO_CCpoint_ON_ELEMENT + j])
									{
										break;
									}
								}
								// 未設定のNode_To_Node取得
								if (k == Total_Control_Point_To_Node[NE])
								{
									Node_To_Node[NE * Total_Control_Point_to_mesh[Total_mesh] + k] = Controlpoint_of_Element[NELOVER[i * MAX_N_ELEMENT_OVER + jj] * MAX_NO_CCpoint_ON_ELEMENT + j];
									Total_Control_Point_To_Node[NE]++;
								}
							}
						}
					}
				}
			}
		}
		// 順番に並び替える
		for (i = 0; i < K_DIVISION_LENGE; i++)
		{
			if (N + i < Total_Control_Point_to_mesh[Total_mesh])
			{
				for (j = 0; j < Total_Control_Point_To_Node[i]; j++)
				{
					int Min = Node_To_Node[i * Total_Control_Point_to_mesh[Total_mesh] + j], No = j;
					for (k = j; k < Total_Control_Point_To_Node[i]; k++)
					{
						if (Min > Node_To_Node[i * Total_Control_Point_to_mesh[Total_mesh] + k])
						{
							Min = Node_To_Node[i * Total_Control_Point_to_mesh[Total_mesh] + k];
							No = k;
						}
					}
					for (k = No; k > j; k--)
					{
						Node_To_Node[i * Total_Control_Point_to_mesh[Total_mesh] + k] = Node_To_Node[i * Total_Control_Point_to_mesh[Total_mesh] + k - 1];
					}
					Node_To_Node[i * Total_Control_Point_to_mesh[Total_mesh] + j] = Min;
				}
			}
		}

		// 節点からcol ptrを求める
		ii = 0;
		k = 0;
		for (i = 0; i < K_DIVISION_LENGE; i++)
		{
			for (ii = 0; ii < DIMENSION; ii++)
			{
				if (N + i < Total_Control_Point_to_mesh[Total_mesh])
				{
					i_index = Index_Dof[(N + i) * DIMENSION + ii];
					k = 0;
					if (i_index >= 0)
					{
						K_Whole_Ptr[i_index + 1] = K_Whole_Ptr[i_index];
						for (j = 0; j < Total_Control_Point_To_Node[i]; j++)
						{
							for (jj = 0; jj < DIMENSION; jj++)
							{
								j_index = Index_Dof[Node_To_Node[i * Total_Control_Point_to_mesh[Total_mesh] + j] * DIMENSION + jj];
								if (j_index >= 0 && j_index >= i_index)
								{
									K_Whole_Ptr[i_index + 1]++;
									K_Whole_Col[K_Whole_Ptr[i_index] + k] = j_index;
									k++;
								}
							}
						}
					}
				}
			}
		}
	}
}


// valを求める
void Make_K_Whole_Val(int *real_Total_Element_to_mesh, int *real_element, int *Element_mesh, int *NNLOVER,
					  int *No_Control_point_ON_ELEMENT, int *Element_patch, double *Jac, double *Jac_ex, double *B_Matrix, double *B_Matrix_ex,
					  double *D, int *Index_Dof, int *Controlpoint_of_Element, int *K_Whole_Ptr, int *K_Whole_Col, double *K_Whole_Val, int *NELOVER,
					  double *Loc_parameter_on_Glo, double *Loc_parameter_on_Glo_ex, double *Position_Knots, int *Total_Knot_to_patch_dim,
					  int *Order, int *ENC, int *Total_Control_Point_to_mesh, double *Node_Coordinate, int *INC, int *No_knot, int *No_Control_point)
{
	int i, j, j1, j2, k1, k2, l;
	int a, b, re;

	double *K_EL = (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * MAX_KIEL_SIZE);			// K_EL[MAX_KIEL_SIZE][MAX_KIEL_SIZE]
	double *coupled_K_EL= (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * MAX_KIEL_SIZE);		// coupled_K_EL[MAX_KIEL_SIZE][MAX_KIEL_SIZE]

	for (re = 0; re < real_Total_Element_to_mesh[Total_mesh]; re++)
	{
		i = real_element[re];
		printf("element : %d\n\n", i);

		if (Element_mesh[i] == 0 && re == 0) // 2つめの条件は効率化のため
		{
			// GP_2D = POW_Ng;
			Make_gauss_array(0);
		}
		else if (Element_mesh[i] > 0)
		{
			if (NNLOVER[i] == 1 && (NNLOVER[real_element[re - 1]] != 1 || Element_mesh[real_element[re - 1]] == 0)) // 2つめ以降の条件は効率化のため
			{
				// GP_2D = POW_Ng;
				Make_gauss_array(0);
			}
			else if (NNLOVER[i] >= 2 && (NNLOVER[real_element[re - 1]] == 1 || Element_mesh[real_element[re - 1]] == 0)) // 2つめ以降の条件は効率化のため
			{
				// GP_2D = POW_Ng_extended;
				Make_gauss_array(1);
			}
		}

		// 各要素のK_ELを求める
		Make_K_EL(i, K_EL, No_Control_point_ON_ELEMENT, Element_patch, Jac, Jac_ex, B_Matrix, B_Matrix_ex, D);

		// Valを求める
		for (j1 = 0; j1 < No_Control_point_ON_ELEMENT[Element_patch[i]]; j1++)
		{
			for (j2 = 0; j2 < DIMENSION; j2++)
			{
				a = Index_Dof[Controlpoint_of_Element[i * MAX_NO_CCpoint_ON_ELEMENT + j1] * DIMENSION + j2];
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
										K_Whole_Val[l] += K_EL[(j1 * DIMENSION + j2) * MAX_KIEL_SIZE + k1 * DIMENSION + k2];
										break;
									}
								}
							}
						}
					}
				}
			}
		}

		if (Element_mesh[i] > 0 && NNLOVER[i] > 0) // ローカルメッシュ上の要素について, 重なっている要素が存在するとき
		{
			for (j = 0; j < NNLOVER[i]; j++)
			{
				// 各要素のcoupled_K_ELを求める
				Make_coupled_K_EL(i, NELOVER[i * MAX_N_ELEMENT_OVER + j], coupled_K_EL, No_Control_point_ON_ELEMENT, 
								  Jac, Jac_ex, B_Matrix, B_Matrix_ex, D,
								  Loc_parameter_on_Glo, Loc_parameter_on_Glo_ex, Position_Knots,
								  Total_Knot_to_patch_dim, Order, ENC, Controlpoint_of_Element,
								  Total_Control_Point_to_mesh, Element_patch, Node_Coordinate, INC,
								  No_knot, No_Control_point);

				// Valを求める
				for (j1 = 0; j1 < No_Control_point_ON_ELEMENT[Element_patch[NELOVER[i * MAX_N_ELEMENT_OVER + j]]]; j1++)
				{
					for (j2 = 0; j2 < DIMENSION; j2++)
					{
						a = Index_Dof[Controlpoint_of_Element[NELOVER[i * MAX_N_ELEMENT_OVER + j] * MAX_NO_CCpoint_ON_ELEMENT + j1] * DIMENSION + j2];
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
												K_Whole_Val[l] += coupled_K_EL[(j1 * DIMENSION + j2) * MAX_KIEL_SIZE + k1 * DIMENSION + k2];
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

	free(K_EL), free(coupled_K_EL);
}


// 要素合成マトリックス
void Make_K_EL(int El_No, double *K_EL, int *No_Control_point_ON_ELEMENT, int *Element_patch,
			   double *Jac, double *Jac_ex, double *B_Matrix, double *B_Matrix_ex, double *D)
{
	int i, j, k, l;
	
	int KIEL_SIZE = No_Control_point_ON_ELEMENT[Element_patch[El_No]] * DIMENSION;

	double *B = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * MAX_KIEL_SIZE); // B[D_MATRIX_SIZE][MAX_KIEL_SIZE];
	double *K1 = (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * MAX_KIEL_SIZE); // K1[MAX_KIEL_SIZE][MAX_KIEL_SIZE];
	double J;

	for (i = 0; i < KIEL_SIZE; i++)
	{
		for (j = 0; j < KIEL_SIZE; j++)
		{
			K_EL[i * MAX_KIEL_SIZE + j] = 0.0;
		}
	}

	for (i = 0; i < GP_2D; i++)
	{
		// J, B の作成
		if (GP_2D == POW_Ng)
		{
			J = Jac[El_No * GP_2D + i];
			for (j = 0; j < D_MATRIX_SIZE; j++)
			{
				for (k = 0; k < MAX_KIEL_SIZE; k++)
				{
					B[j * MAX_KIEL_SIZE + k] = B_Matrix[El_No * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + j * MAX_KIEL_SIZE + k];
				}
			}
		}
		else if (GP_2D == POW_Ng_extended)
		{
			J = Jac_ex[El_No * GP_2D + i];
			for (j = 0; j < D_MATRIX_SIZE; j++)
			{
				for (k = 0; k < MAX_KIEL_SIZE; k++)
				{
					B[j * MAX_KIEL_SIZE + k] = B_Matrix_ex[El_No * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + j * MAX_KIEL_SIZE + k];
				}
			}
		}

		BDBJ(KIEL_SIZE, B, D, J, K1);
		for (k = 0; k < KIEL_SIZE; k++)
		{
			for (l = 0; l < KIEL_SIZE; l++)
			{
				K_EL[k * MAX_KIEL_SIZE + l] += w[i] * K1[k * MAX_KIEL_SIZE + l];
			}
		}
	}

	free(B), free(K1);
}


// 結合要素剛性マトリックス
void Make_coupled_K_EL(int El_No_loc, int El_No_glo, double *coupled_K_EL, int *No_Control_point_ON_ELEMENT, 
					   double *Jac, double *Jac_ex, double *B_Matrix, double *B_Matrix_ex, double *D,
					   double *Loc_parameter_on_Glo, double *Loc_parameter_on_Glo_ex, double *Position_Knots,
					   int *Total_Knot_to_patch_dim, int *Order, int *ENC, int *Controlpoint_of_Element,
					   int *Total_Control_Point_to_mesh, int *Element_patch, double *Node_Coordinate, int *INC,
					   int *No_knot, int *No_Control_point)
{
	int i, j, jj, k, l;

	int BDBJ_flag;
	int KIEL_SIZE = No_Control_point_ON_ELEMENT[Element_patch[El_No_glo]] * DIMENSION;

	double *B = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * MAX_KIEL_SIZE); // B[D_MATRIX_SIZE][MAX_KIEL_SIZE];
	double *BG = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * MAX_KIEL_SIZE); // BG[D_MATRIX_SIZE][MAX_KIEL_SIZE];
	double *K1 = (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * MAX_KIEL_SIZE); // K1[MAX_KIEL_SIZE][MAX_KIEL_SIZE];
	double J;

	for (i = 0; i < KIEL_SIZE; i++)
	{
		for (j = 0; j < KIEL_SIZE; j++)
		{
			coupled_K_EL[i * MAX_KIEL_SIZE + j] = 0.0;
		}
	}

	for (i = 0; i < GP_2D; i++)
	{
		// J, B, BG の作成
		double para[DIMENSION];
		double G_Gxi[DIMENSION];

		if (GP_2D == POW_Ng)
		{
			J = Jac[El_No_loc * GP_2D + i];
			for (j = 0; j < D_MATRIX_SIZE; j++)
			{
				for (k = 0; k < MAX_KIEL_SIZE; k ++)
				{
					B[j * MAX_KIEL_SIZE + k] = B_Matrix[El_No_loc * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + j * MAX_KIEL_SIZE + k];
				}
			}
			para[0] = Loc_parameter_on_Glo[El_No_loc * GP_2D * DIMENSION + i * DIMENSION + 0];
			para[1] = Loc_parameter_on_Glo[El_No_loc * GP_2D * DIMENSION + i * DIMENSION + 1];
		}
		else if (GP_2D == POW_Ng_extended)
		{
			J = Jac_ex[El_No_loc * GP_2D + i];
			for (j = 0; j < D_MATRIX_SIZE; j++)
			{
				for (k = 0; k < MAX_KIEL_SIZE; k ++)
				{
					B[j * MAX_KIEL_SIZE + k] = B_Matrix_ex[El_No_loc * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + j * MAX_KIEL_SIZE + k];
				}
			}
			para[0] = Loc_parameter_on_Glo_ex[El_No_loc * GP_2D * DIMENSION + i * DIMENSION + 0];
			para[1] = Loc_parameter_on_Glo_ex[El_No_loc * GP_2D * DIMENSION + i * DIMENSION + 1];
		}

		// 要素内外判定
		if (para[0] >= Position_Knots[Total_Knot_to_patch_dim[0 * DIMENSION + 0] + Order[0 * DIMENSION + 0] + ENC[El_No_glo * DIMENSION + 0]] &&
			para[0] <  Position_Knots[Total_Knot_to_patch_dim[0 * DIMENSION + 0] + Order[0 * DIMENSION + 0] + ENC[El_No_glo * DIMENSION + 0] + 1] &&
			para[1] >= Position_Knots[Total_Knot_to_patch_dim[0 * DIMENSION + 1] + Order[0 * DIMENSION + 1] + ENC[El_No_glo * DIMENSION + 1]] &&
			para[1] <  Position_Knots[Total_Knot_to_patch_dim[0 * DIMENSION + 1] + Order[0 * DIMENSION + 1] + ENC[El_No_glo * DIMENSION + 1] + 1])
		{
			BDBJ_flag = 1;

			// 親要素座標の算出
			G_Gxi[0] = -1.0 + 2.0 * (para[0] - Position_Knots[Total_Knot_to_patch_dim[0 * DIMENSION + 0] + Order[0 * DIMENSION + 0] + ENC[El_No_glo * DIMENSION + 0]])
						/ (Position_Knots[Total_Knot_to_patch_dim[0 * DIMENSION + 0] + Order[0 * DIMENSION + 0] + ENC[El_No_glo * DIMENSION + 0] + 1] - Position_Knots[Total_Knot_to_patch_dim[0 * DIMENSION + 0] + Order[0 * DIMENSION + 0] + ENC[El_No_glo * DIMENSION + 0]]);
			G_Gxi[1] = -1.0 + 2.0 * (para[1] - Position_Knots[Total_Knot_to_patch_dim[0 * DIMENSION + 1] + Order[0 * DIMENSION + 1] + ENC[El_No_glo * DIMENSION + 1]])
						/ (Position_Knots[Total_Knot_to_patch_dim[0 * DIMENSION + 1] + Order[0 * DIMENSION + 1] + ENC[El_No_glo * DIMENSION + 1] + 1] - Position_Knots[Total_Knot_to_patch_dim[0 * DIMENSION + 1] + Order[0 * DIMENSION + 1] + ENC[El_No_glo * DIMENSION + 1]]);
		}
		else
		{
			BDBJ_flag = 0;
		}

		// 要素内であるとき, 結合要素剛性マトリックス計算
		if (BDBJ_flag)
		{
			// 重なるグローバル要素のBマトリックス
			Make_BG_Matrix(El_No_glo, BG, G_Gxi, Controlpoint_of_Element, No_Control_point_ON_ELEMENT,
						   Total_Control_Point_to_mesh, Element_patch, Node_Coordinate, INC, Order, Position_Knots,
						   Total_Knot_to_patch_dim, No_knot, No_Control_point);
			// BGTDBLJの計算
			coupled_BDBJ(KIEL_SIZE, B, D, BG, J, K1);
			for (k = 0; k < KIEL_SIZE; k++)
			{
				for (l = 0; l < KIEL_SIZE; l++)
				{
					coupled_K_EL[k * MAX_KIEL_SIZE + l] += w[i] * K1[k * MAX_KIEL_SIZE + l];
				}
			}
		}
	}

	free(B), free(BG), free(K1);
}


// BGマトリックスを求める関数
void Make_BG_Matrix(int El_No, double *BG, double Local_coord[DIMENSION], int *Controlpoint_of_Element, int *No_Control_point_ON_ELEMENT,
					int *Total_Control_Point_to_mesh, int *Element_patch, double *Node_Coordinate, int *INC, int *Order, double *Position_Knots,
					int *Total_Knot_to_patch_dim, int *No_knot, int *No_Control_point)
{
	double a[DIMENSION][DIMENSION];
	double *b = (double *)malloc(sizeof(double) * DIMENSION * MAX_NO_CCpoint_ON_ELEMENT); // b[DIMENSION][MAX_NO_CCpoint_ON_ELEMENT]

	int i, j, k;

	for (i = 0; i < DIMENSION; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			a[i][j] = 0.0;
			for (k = 0; k < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; k++)
			{
				a[i][j] += dShape_func(k, j, Local_coord, El_No, Controlpoint_of_Element, No_Control_point_ON_ELEMENT,
									   Total_Control_Point_to_mesh, Element_patch, Node_Coordinate, INC, Order, Position_Knots,
									   Total_Knot_to_patch_dim, No_knot, No_Control_point)
						 * Node_Coordinate[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + k] * (DIMENSION + 1) + i];
			}
		}
	}

	InverseMatrix_2x2(a);

	for (i = 0; i < DIMENSION; i++)
	{
		for (j = 0; j < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; j++)
		{
			b[i * MAX_NO_CCpoint_ON_ELEMENT + j] = 0.0;
			for (k = 0; k < DIMENSION; k++)
			{
				b[i * MAX_NO_CCpoint_ON_ELEMENT + j] += a[k][i]
													  * dShape_func(j, k, Local_coord, El_No, Controlpoint_of_Element, No_Control_point_ON_ELEMENT,
																	Total_Control_Point_to_mesh, Element_patch, Node_Coordinate, INC, Order, Position_Knots,
																	Total_Knot_to_patch_dim, No_knot, No_Control_point);
			}
		}
	}

	// 2次元
	for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++)
	{
		BG[0 * MAX_KIEL_SIZE + 2 * i] = b[0 * MAX_NO_CCpoint_ON_ELEMENT + i];
		BG[0 * MAX_KIEL_SIZE + 2 * i + 1] = 0.0;
		BG[1 * MAX_KIEL_SIZE + 2 * i] = 0.0;
		BG[1 * MAX_KIEL_SIZE + 2 * i + 1] = b[1 * MAX_NO_CCpoint_ON_ELEMENT + i];
		BG[2 * MAX_KIEL_SIZE + 2 * i] = b[1 * MAX_NO_CCpoint_ON_ELEMENT + i];
		BG[2 * MAX_KIEL_SIZE + 2 * i + 1] = b[0 * MAX_NO_CCpoint_ON_ELEMENT + i];
	}

	free(b);
}


// ガウスの数値積分
void BDBJ(int KIEL_SIZE, double *B, double *D, double J, double *K_EL)
{
	int i, j, k;
	double *BD = (double *)calloc(MAX_KIEL_SIZE * D_MATRIX_SIZE, sizeof(double)); // BD[MAX_KIEL_SIZE][D_MATRIX_SIZE];

	// [B]T[D][B]Jの計算
	for (i = 0; i < KIEL_SIZE; i++)
	{
		for (j = 0; j < D_MATRIX_SIZE; j++)
		{
			for (k = 0; k < D_MATRIX_SIZE; k++)
			{
				BD[i * D_MATRIX_SIZE + j] += B[k * MAX_KIEL_SIZE + i] * D[k * D_MATRIX_SIZE + j];
			}
		}
	}

	for (i = 0; i < KIEL_SIZE; i++)
	{
		for (j = 0; j < KIEL_SIZE; j++)
		{
			K_EL[i * MAX_KIEL_SIZE + j] = 0.0;
			for (k = 0; k < D_MATRIX_SIZE; k++)
			{
				K_EL[i * MAX_KIEL_SIZE + j] += BD[i * D_MATRIX_SIZE + k] * B[k * MAX_KIEL_SIZE + j];
			}
			K_EL[i * MAX_KIEL_SIZE + j] *= J;
		}
	}

	free(BD);
}


// 結合ガウスの数値積分
void coupled_BDBJ(int KIEL_SIZE, double *B, double *D, double *BG, double J, double *K_EL)
{
	int i, j, k;
	double *BD = (double *)calloc(MAX_KIEL_SIZE * D_MATRIX_SIZE, sizeof(double)); // BD[MAX_KIEL_SIZE][D_MATRIX_SIZE];

	//[B]T[D][B]Jの計算
	for (i = 0; i < KIEL_SIZE; i++)
	{
		for (j = 0; j < D_MATRIX_SIZE; j++)
		{
			for (k = 0; k < D_MATRIX_SIZE; k++)
			{
				BD[i * D_MATRIX_SIZE + j] += BG[k * MAX_KIEL_SIZE + i] * D[k * D_MATRIX_SIZE + j];
			}
		}
	}

	for (i = 0; i < KIEL_SIZE; i++)
	{
		for (j = 0; j < KIEL_SIZE; j++)
		{
			K_EL[i * MAX_KIEL_SIZE + j] = 0.0;
			for (k = 0; k < D_MATRIX_SIZE; k++)
			{
				K_EL[i * MAX_KIEL_SIZE + j] += BD[i * D_MATRIX_SIZE + k] * B[k * MAX_KIEL_SIZE + j];
			}
			K_EL[i * MAX_KIEL_SIZE + j] *= J;
		}
	}

	free(BD);
}


// F vector
// 荷重の行列を作る
void Make_F_Vec(double *rhs_vec, int *Total_Load_to_mesh, int *Index_Dof, int *Load_Node_Dir, double *Value_of_Load)
{
	int i, index;
	for (i = 0; i < K_Whole_Size; i++)
		rhs_vec[i] = 0.0;
	
	for (i = 0; i < Total_Load_to_mesh[Total_mesh]; i++)
	{
		index = Index_Dof[Load_Node_Dir[i * 2 + 0] * DIMENSION + Load_Node_Dir[i * 2 + 1]];
		if (index >= 0)
			rhs_vec[index] += Value_of_Load[i];
	}
}


// 強制変位対策
void Make_F_Vec_disp_const(int *Total_Constraint_to_mesh, int *real_Total_Element_to_mesh, int *No_Control_point_ON_ELEMENT, int *Element_patch,
						   int *Index_Dof, int *Controlpoint_of_Element, int *real_El_No_on_mesh, int *Total_Element_to_mesh, 
						   int *Constraint_Node_Dir, double *Value_of_Constraint, double *rhs_vec, int *real_element,
						   double *Jac, double *Jac_ex, double *B_Matrix, double *B_Matrix_ex, double *D)
{
	int ie, idir, inode, jdir, jnode, kk_const;
	int ii, iii, b, bb, jj, j1, j2, ii_local, jj_local;

	int i;

	double *K_EL = (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * MAX_KIEL_SIZE);			// K_EL[MAX_KIEL_SIZE][MAX_KIEL_SIZE]

	Make_gauss_array(0);

	for (ie = 0; ie < real_Total_Element_to_mesh[Total_mesh]; ie++)
	{
		i = real_element[ie];
		int KIEL_SIZE = No_Control_point_ON_ELEMENT[Element_patch[i]] * DIMENSION;

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

		if (iii > 0)
		{
			Make_K_EL(i, K_EL, No_Control_point_ON_ELEMENT, Element_patch, Jac, Jac_ex, B_Matrix, B_Matrix_ex, D);
			for (idir = 0; idir < DIMENSION; idir++)
			{
				for (inode = 0; inode < No_Control_point_ON_ELEMENT[Element_patch[i]]; inode++)
				{
					// ii = Controlpoint_of_Element[real_El_No_on_mesh[Total_mesh * Total_Element_to_mesh[Total_mesh] + ie] * MAX_NO_CCpoint_ON_ELEMENT + inode] * DIMENSION + idir;
					ii = Controlpoint_of_Element[0 + inode] * DIMENSION + idir;
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
									jj_local = jnode * DIMENSION + jdir;
									for (kk_const = 0; kk_const < Total_Constraint_to_mesh[Total_mesh]; kk_const++)
									{
										if (Controlpoint_of_Element[i * MAX_NO_CCpoint_ON_ELEMENT + jnode] == Constraint_Node_Dir[kk_const * 2 + 0] && jdir == Constraint_Node_Dir[kk_const * 2 + 1])
										{
											rhs_vec[b] -= K_EL[ii_local * MAX_KIEL_SIZE + jj_local] * Value_of_Constraint[kk_const];
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
	free(K_EL);
}


// 分布荷重の等価節点力を足す
void Add_Equivalent_Nodal_Force_to_F_Vec(int *Total_Control_Point_to_mesh, int *Index_Dof, double *rhs_vec, double *Equivalent_Nodal_Force)
{
	int i, j, index;
	for (j = 0; j < DIMENSION; j++)
	{
		for (i = 0; i < Total_Control_Point_to_mesh[Total_mesh]; i++)
		{
			index = Index_Dof[i * DIMENSION + j];
			if (index >= 0)
			{
				rhs_vec[index] += Equivalent_Nodal_Force[i * DIMENSION + j];
			}
		}
	}
}



// PCG solver
// 前処理付共役勾配法により[K]{d}={f}を解く
void PCG_Solver(int max_itetarion, double eps, double *K_Whole_Val, int *K_Whole_Ptr, int *K_Whole_Col, double *sol_vec, double *rhs_vec, int *Total_Control_Point_on_mesh, int *Index_Dof)
{
	int i, j, k;
	int ndof = K_Whole_Size;

	double *r = (double *)malloc(sizeof(double) * ndof);
	double *p = (double *)calloc(ndof, sizeof(double));
	double *y = (double *)malloc(sizeof(double) * ndof);
	double *r2 = (double *)calloc(ndof, sizeof(double));

	// 初期化
	for (i = 0; i < ndof; i++)
		sol_vec[i] = 0.0;

	// 前処理行列作成
	double *M = (double *)calloc(MAX_NON_ZERO, sizeof(double));
	int *M_Ptr = (int *)calloc(MAX_K_WHOLE_SIZE + 1, sizeof(int));
	int *M_Col = (int *)calloc(MAX_NON_ZERO, sizeof(int));
	Make_M(M, M_Ptr, M_Col, ndof, Total_Control_Point_on_mesh, K_Whole_Val, K_Whole_Ptr, K_Whole_Col, Index_Dof);

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
					temp1 = RowCol_to_icount(i, j, K_Whole_Ptr, K_Whole_Col); // temp_array_K[i][j]
				}
				else if (i > j)
				{
					temp1 = RowCol_to_icount(j, i, K_Whole_Ptr, K_Whole_Col); // temp_array_K[i][j] = temp_array_K[j][i]
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

		// 解x, 残差rの更新
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


void Make_M(double *M, int *M_Ptr, int *M_Col, int ndof, int *Total_Control_Point_on_mesh, double *K_Whole_Val, int *K_Whole_Ptr, int *K_Whole_Col, int *Index_Dof)
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


void CG(int ndof, double *solution_vec, double *M, int *M_Ptr, int *M_Col, double *right_vec)
{
	// CG solver
	double *gg = (double *)malloc(sizeof(double) * MAX_K_WHOLE_SIZE); // gg[MAX_K_WHOLE_SIZE]
	double *dd = (double *)malloc(sizeof(double) * MAX_K_WHOLE_SIZE); // dd[MAX_K_WHOLE_SIZE]
	double *pp = (double *)malloc(sizeof(double) * MAX_K_WHOLE_SIZE); // pp[MAX_K_WHOLE_SIZE]
	double qqq, ppp, rrr;
	double alphak, betak;
	int i, ii, itr, istop;
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
	free(gg), free(dd), free(pp);
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


int RowCol_to_icount(int row, int col, int *K_Whole_Ptr, int *K_Whole_Col)
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


// tool
// 逆行列を元の行列に代入
double InverseMatrix_2x2(double M[DIMENSION][DIMENSION])
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


double InverseMatrix_3x3(double M[DIMENSION][DIMENSION])
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

	return det;
}


// Shape Function
double Shape_func(int I_No, double Local_coord[DIMENSION], int El_No, double *Node_Coordinate, int *Total_Control_Point_to_mesh,
				  int *Controlpoint_of_Element, int *INC, int *Element_patch, int *Order, int *No_Control_point_ON_ELEMENT,
				  double *Position_Knots, int *Total_Knot_to_patch_dim, int *No_knot, int *No_Control_point)
{
	int i, j;
	double R, weight_func = 0.0;

	double Position_Data_param[DIMENSION];

	double *shape_func = (double *)malloc(sizeof(double) * MAX_NO_CCpoint_ON_ELEMENT);
	double *Shape = (double *)malloc(sizeof(double) * DIMENSION * Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER); // Shape[DIMENSION][MAX_N_NODE][10]
	double *dShape = (double *)malloc(sizeof(double) * DIMENSION * Total_Control_Point_to_mesh[Total_mesh]); // dShape[DIMENSION][MAX_N_NODE]

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
			shape_func[i] *= Shape[j * (Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + INC[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + j] * MAX_ORDER + Order[Element_patch[El_No] * DIMENSION + j]];
		}
		weight_func += shape_func[i] * Node_Coordinate[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * (DIMENSION + 1) + DIMENSION];
	}

	free(Shape), free(dShape);

	if (I_No < No_Control_point_ON_ELEMENT[Element_patch[El_No]])
		R = shape_func[I_No] * Node_Coordinate[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + I_No] * (DIMENSION + 1) + DIMENSION] / weight_func;

	else
		R = ERROR;

	free(shape_func);
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
				   int *Total_Control_Point_to_mesh, int *Element_patch, double *Node_Coordinate, int *INC, int *Order, double *Position_Knots,
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
			   * dShapeFunc_from_paren(xez, El_No, INC, Controlpoint_of_Element, Position_Knots, Total_Knot_to_patch_dim, Element_patch);
		}
		else if (xez == 1)
		{
			dR = dShape_func2[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + I_No]]
			   * dShapeFunc_from_paren(xez, El_No, INC, Controlpoint_of_Element, Position_Knots, Total_Knot_to_patch_dim, Element_patch);
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

	double *shape_func = (double *)malloc(sizeof(double) * MAX_NO_CCpoint_ON_ELEMENT);
	double *Shape = (double *)malloc(sizeof(double) * DIMENSION * Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER); // Shape[DIMENSION][MAX_N_NODE][10]
	double *dShape = (double *)malloc(sizeof(double) * DIMENSION * Total_Control_Point_to_mesh[Total_mesh]); // dShape[DIMENSION][MAX_N_NODE]

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
			shape_func[i] *= Shape[j * (Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + INC[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + j] * MAX_ORDER + Order[Element_patch[El_No] * DIMENSION + j]];
		}
		weight_func += shape_func[i] * Node_Coordinate[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * (DIMENSION + 1) + DIMENSION];
	}
	for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++)
	{
		dWeight_func1 += dShape[0 * Total_Control_Point_to_mesh[Total_mesh] + INC[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + 0]] * Shape[1 * (Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + INC[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + 1] * MAX_ORDER + Order[Element_patch[El_No] * DIMENSION + 1]] * Node_Coordinate[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * (DIMENSION + 1) + DIMENSION];
		dWeight_func2 += Shape[0 * (Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + INC[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + 0] * MAX_ORDER + Order[Element_patch[El_No] * DIMENSION + 0]] * dShape[1 * Total_Control_Point_to_mesh[Total_mesh] + INC[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + 1]] * Node_Coordinate[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * (DIMENSION + 1) + DIMENSION];
	}
	for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++)
	{
		dShape_func1[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i]] = Node_Coordinate[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * (DIMENSION + 1) + DIMENSION] * (weight_func * dShape[0 * Total_Control_Point_to_mesh[Total_mesh] + INC[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + 0]] * Shape[1 * (Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + INC[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + 1] * MAX_ORDER + Order[Element_patch[El_No] * DIMENSION + 1]] - dWeight_func1 * shape_func[i]) / (weight_func * weight_func);
		dShape_func2[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i]] = Node_Coordinate[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * (DIMENSION + 1) + DIMENSION] * (weight_func * Shape[0 * (Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + INC[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + 0] * MAX_ORDER + Order[Element_patch[El_No] * DIMENSION + 0]] * dShape[1 * Total_Control_Point_to_mesh[Total_mesh] + INC[Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + 1]] - dWeight_func2 * shape_func[i]) / (weight_func * weight_func);
	}

	free(Shape), free(dShape), free(shape_func);
}


double dShapeFunc_from_paren(int j, int e, int *INC, int *Controlpoint_of_Element, double *Position_Knots, int *Total_Knot_to_patch_dim, int *Element_patch)
{
	int i;
	double dPosition_Data_param;

	i = INC[Controlpoint_of_Element[e * MAX_NO_CCpoint_ON_ELEMENT + 0] * DIMENSION + j];

	dPosition_Data_param = (Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + i + 1] - Position_Knots[Total_Knot_to_patch_dim[Element_patch[e] * DIMENSION + j] + i]) / 2.0;
	return dPosition_Data_param;
}


// Newton-Raphson法, from NURBSviewer
double BasisFunc(double *knot_vec, int knot_index, int order, double xi, double *output, double *d_output)
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


double rBasisFunc(double *knot_vec, int knot_index, int order, double xi, double *output, double *d_output)
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


double lBasisFunc(double *knot_vec, int knot_index, int cntl_p_n, int order, double xi, double *output, double *d_output)
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
				double *output_xi, double *output_eta,
				double *Position_Knots, int *Total_Knot_to_patch_dim, int *No_Control_point, int *Order,
				double *Control_Coord_x, double *Control_Coord_y, double *Control_Weight, int *No_knot)
{
	double temp_xi, temp_eta;
	double temp_x, temp_y;
	double temp_matrix[2][2];
	double temp_dxi, temp_deta;
	double temp_tol_x, temp_tol_y;

	(*output_xi) = 0;
	(*output_eta) = 0;

	int i;
	// int repeat = 1000;
	// double tol = 10e-8;
	int repeat = 100;
	double tol = 10e-14;

	double *temp_Position_Knots_xi = (double *)malloc(sizeof(double) * No_knot[0 * DIMENSION + 0]);
	double *temp_Position_Knots_eta = (double *)malloc(sizeof(double) * No_knot[0 * DIMENSION + 1]);
	for (i = 0; i < No_knot[0 * DIMENSION + 0]; i++)
	{
		temp_Position_Knots_xi[i] = Position_Knots[Total_Knot_to_patch_dim[0 * DIMENSION + 0] + i];
	}
	for (i = 0; i < No_knot[0 * DIMENSION + 1]; i++)
	{
		temp_Position_Knots_eta[i] = Position_Knots[Total_Knot_to_patch_dim[0 * DIMENSION + 1] + i];
	}

	// 初期値の設定
	temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;

	for (i = 0; i < repeat; i++)
	{
		rNURBS_surface(temp_Position_Knots_xi, temp_Position_Knots_eta,
					   Control_Coord_x, Control_Coord_y,
					   No_Control_point[0 * DIMENSION + 0], No_Control_point[0 * DIMENSION + 1],
					   Control_Weight, Order[0 * DIMENSION + 0], Order[0 * DIMENSION + 1],
					   temp_xi, temp_eta,
					   &temp_x, &temp_y,
					   &temp_matrix[0][0], &temp_matrix[0][1],
					   &temp_matrix[1][0], &temp_matrix[1][1]);

		temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;

		// 収束した場合
		if (temp_tol_x + temp_tol_y < tol)
		{
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;
			return i;
		}

		InverseMatrix_2x2(temp_matrix);

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
	}

	// 初期値の設定
	temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;

	for (i = 0; i < repeat; i++)
	{
		lNURBS_surface(temp_Position_Knots_xi, temp_Position_Knots_eta,
					   Control_Coord_x, Control_Coord_y,
					   No_Control_point[0 * DIMENSION + 0], No_Control_point[0 * DIMENSION + 1],
					   Control_Weight, Order[0 * DIMENSION + 0], Order[0 * DIMENSION + 1],
					   temp_xi, temp_eta,
					   &temp_x, &temp_y,
					   &temp_matrix[0][0], &temp_matrix[0][1],
					   &temp_matrix[1][0], &temp_matrix[1][1]);

		temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;

		// 収束した場合
		if (temp_tol_x + temp_tol_y < tol)
		{
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;
			return i;
		}

		InverseMatrix_2x2(temp_matrix);

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
	}

	// 初期値の設定
	temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;

	for (i = 0; i < repeat; i++)
	{
		rlNURBS_surface(temp_Position_Knots_xi, temp_Position_Knots_eta,
						Control_Coord_x, Control_Coord_y,
						No_Control_point[0 * DIMENSION + 0], No_Control_point[0 * DIMENSION + 1],
						Control_Weight, Order[0 * DIMENSION + 0], Order[0 * DIMENSION + 1],
						temp_xi, temp_eta,
						&temp_x, &temp_y,
						&temp_matrix[0][0], &temp_matrix[0][1],
						&temp_matrix[1][0], &temp_matrix[1][1]);

		temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;

		// 収束した場合
		if (temp_tol_x + temp_tol_y < tol)
		{
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;
			return i;
		}

		InverseMatrix_2x2(temp_matrix);

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
	}

	// 初期値の設定
	temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;

	for (i = 0; i < repeat; i++)
	{
		lrNURBS_surface(temp_Position_Knots_xi, temp_Position_Knots_eta,
						Control_Coord_x, Control_Coord_y,
						No_Control_point[0 * DIMENSION + 0], No_Control_point[0 * DIMENSION + 1],
						Control_Weight, Order[0 * DIMENSION + 0], Order[0 * DIMENSION + 1],
						temp_xi, temp_eta,
						&temp_x, &temp_y,
						&temp_matrix[0][0], &temp_matrix[0][1],
						&temp_matrix[1][0], &temp_matrix[1][1]);

		temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;

		// 収束した場合
		if (temp_tol_x + temp_tol_y < tol)
		{
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;
			return i;
		}

		InverseMatrix_2x2(temp_matrix);

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
	}

	free(temp_Position_Knots_xi), free(temp_Position_Knots_eta);

	return 0;
}


// for s_IGA overlay


// output
void K_output_svg(int *K_Whole_Ptr, int *K_Whole_Col)
{
	// [K] = [[K^G, K^GL], [K^GL, K^L]]

	int i, j;
	int ndof = K_Whole_Size;

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
				temp_count = RowCol_to_icount(i, j, K_Whole_Ptr, K_Whole_Col);
			}
			else if (i > j)
			{
				temp_count = RowCol_to_icount(j, i, K_Whole_Ptr, K_Whole_Col);
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

