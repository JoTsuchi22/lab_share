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
void Get_Input_1(int tm, char **argv, information info)
{
	char s[256];
	int temp_i;

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
	info.Total_Patch_on_mesh[tm] = temp_i;
	info.Total_Patch_to_mesh[tm + 1] = info.Total_Patch_to_mesh[tm] + temp_i;
	printf("info.Total_Patch_to_mesh[%d] = %d\n", tm, info.Total_Patch_to_mesh[tm]);

	// コントロールポイント数
	fscanf(fp, "%d", &temp_i);
	fgets(s, 256, fp);
	printf("Total_Control_Point: %d\n", temp_i);
	info.Total_Control_Point_on_mesh[tm] = temp_i;
	info.Total_Control_Point_to_mesh[tm + 1] = info.Total_Control_Point_to_mesh[tm] + temp_i;
	printf("info.Total_Control_Point_to_mesh[%d] = %d\n", tm, info.Total_Control_Point_to_mesh[tm]);

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
	info.Total_Knot_to_mesh[tm + 1] = info.Total_Knot_to_mesh[tm];
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			fscanf(fp, "%d", &temp_i);
			info.Total_Knot_to_mesh[tm + 1] += temp_i;
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
	for (i = 0; i < info.Total_Patch_on_mesh[tm]; i++)
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
	info.Total_Constraint_to_mesh[tm + 1] = info.Total_Constraint_to_mesh[tm] + Total_Constraint;
	info.Total_Load_to_mesh[tm + 1] = info.Total_Load_to_mesh[tm] + Total_Load;
	info.Total_DistributeForce_to_mesh[tm + 1] = info.Total_DistributeForce_to_mesh[tm] + Total_DistributeForce;

	printf("Total_Constraint = %d\n", Total_Constraint);
	printf("Total_Load = %d\n", Total_Load);
	printf("Total_DistributedForce = %d\n", Total_DistributeForce);

	fclose(fp);
	free(CP);
}


// ファイル読み込み2回目
void Get_Input_2(int tm, char **argv, information info)
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
			info.Order[(i + info.Total_Patch_to_mesh[tm]) * DIMENSION + j] = temp_i;
			printf("info.Order[%d] = %d\n", (i + info.Total_Patch_to_mesh[tm]) * DIMENSION + j, info.Order[(i + info.Total_Patch_to_mesh[tm]) * DIMENSION + j]);
		}
	}
	fgets(s, 256, fp);

	// ノット数
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			fscanf(fp, "%d", &temp_i);
			info.No_knot[(i + info.Total_Patch_to_mesh[tm]) * DIMENSION + j] = temp_i;
			info.Total_Knot_to_patch_dim[(i + info.Total_Patch_to_mesh[tm]) * DIMENSION + j + 1] = info.Total_Knot_to_patch_dim[(i + info.Total_Patch_to_mesh[tm]) * DIMENSION + j] + temp_i;
			printf("info.No_knot[%d] = %d\n", (i + info.Total_Patch_to_mesh[tm]) * DIMENSION + j, info.No_knot[(i + info.Total_Patch_to_mesh[tm]) * DIMENSION + j]);
		}
	}
	fgets(s, 256, fp);

	// 各パッチ各方向のコントロールポイント数
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			fscanf(fp, "%d", &temp_i);
			info.No_Control_point[(i + info.Total_Patch_to_mesh[tm]) * DIMENSION + j] = temp_i;
			printf("info.No_Control_point[%d] = %d\n", (i + info.Total_Patch_to_mesh[tm]) * DIMENSION + j, info.No_Control_point[(i + info.Total_Patch_to_mesh[tm]) * DIMENSION + j]);
		}
	}
	fgets(s, 256, fp);

	for (i = 0; i < No_Patch; i++)
	{
		info.No_Control_point_in_patch[i + info.Total_Patch_to_mesh[tm]] = 1;
	}

	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			info.No_Control_point_in_patch[i + info.Total_Patch_to_mesh[tm]] *= info.No_Control_point[(i + info.Total_Patch_to_mesh[tm]) * DIMENSION + j];
		}
	}

	for (i = 0; i < No_Patch; i++)
	{
		info.Total_Control_Point_to_patch[i + info.Total_Patch_to_mesh[tm] + 1] = info.Total_Control_Point_to_patch[i + info.Total_Patch_to_mesh[tm]] + info.No_Control_point_in_patch[i + info.Total_Patch_to_mesh[tm]];
	}

	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			if (info.No_knot[(i + info.Total_Patch_to_mesh[tm]) * DIMENSION + j] != info.No_Control_point[(i + info.Total_Patch_to_mesh[tm]) * DIMENSION + j] + info.Order[(i + info.Total_Patch_to_mesh[tm]) * DIMENSION + j] + 1)
			{
				printf("wrong relationship between the number of knot vector and the number of control_point \n");
				printf("in mesh_No.%d in patch_No.%d direction:%d\n", tm, i, j);
			}
		}
	}

	for (i = 0; i < No_Patch; i++)
	{
		printf("info.No_Control_point_in_patch[%d] = %d\t", i + info.Total_Patch_to_mesh[tm], info.No_Control_point_in_patch[i + info.Total_Patch_to_mesh[tm]]);
	}
	printf("\n");

	// パッチコネクティビティ
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < info.No_Control_point_in_patch[i + info.Total_Patch_to_mesh[tm]]; j++)
		{
			fscanf(fp, "%d", &temp_i);
			info.Patch_Control_point[info.Total_Control_Point_to_patch[i + info.Total_Patch_to_mesh[tm]] + j] = temp_i;
			if (tm > 0)
			{
				info.Patch_Control_point[info.Total_Control_Point_to_patch[i + info.Total_Patch_to_mesh[tm]] + j] += info.Total_Control_Point_to_mesh[tm];
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
			for (k = 0; k < info.No_knot[(i + info.Total_Patch_to_mesh[tm]) * DIMENSION + j]; k++)
			{
				fscanf(fp, "%lf", &temp_d);
				info.Position_Knots[info.Total_Knot_to_patch_dim[(i + info.Total_Patch_to_mesh[tm]) * DIMENSION + j] + k] = temp_d;
				printf("%le\t", info.Position_Knots[info.Total_Knot_to_patch_dim[(i + info.Total_Patch_to_mesh[tm]) * DIMENSION + j] + k]);
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
			Total_Element += (info.No_Control_point[(i + info.Total_Patch_to_mesh[tm]) * DIMENSION + 0] - info.Order[(i + info.Total_Patch_to_mesh[tm]) * DIMENSION + 0])
						   * (info.No_Control_point[(i + info.Total_Patch_to_mesh[tm]) * DIMENSION + 1] - info.Order[(i + info.Total_Patch_to_mesh[tm]) * DIMENSION + 1]);
			info.No_Control_point_ON_ELEMENT[i + info.Total_Patch_to_mesh[tm]] = (info.Order[(i + info.Total_Patch_to_mesh[tm]) * DIMENSION + 0] + 1) * (info.Order[(i + info.Total_Patch_to_mesh[tm]) * DIMENSION + 1] + 1);
		}
		else if (DIMENSION == 3)
		{
			Total_Element += (info.No_Control_point[(i + info.Total_Patch_to_mesh[tm]) * DIMENSION + 0] - info.Order[(i + info.Total_Patch_to_mesh[tm]) * DIMENSION + 0])
						   * (info.No_Control_point[(i + info.Total_Patch_to_mesh[tm]) * DIMENSION + 1] - info.Order[(i + info.Total_Patch_to_mesh[tm]) * DIMENSION + 1])
						   * (info.No_Control_point[(i + info.Total_Patch_to_mesh[tm]) * DIMENSION + 2] - info.Order[(i + info.Total_Patch_to_mesh[tm]) * DIMENSION + 2]);
			info.No_Control_point_ON_ELEMENT[i + info.Total_Patch_to_mesh[tm]] = (info.Order[(i + info.Total_Patch_to_mesh[tm]) * DIMENSION + 0] + 1) * (info.Order[(i + info.Total_Patch_to_mesh[tm]) * DIMENSION + 1] + 1) * (info.Order[(i + info.Total_Patch_to_mesh[tm]) * DIMENSION + 2] + 1);
		}
	}
	printf("Total_Element = %d\n", Total_Element);
	info.Total_Element_on_mesh[tm] = Total_Element;
	info.Total_Element_to_mesh[tm + 1] = info.Total_Element_to_mesh[tm] + Total_Element;
	printf("info.Total_Element_on_mesh[%d] = %d\n", tm, info.Total_Element_on_mesh[tm]);

	for (i = 0; i < No_Patch; i++)
	{
		printf("info.No_Control_point_ON_ELEMENT[%d] = %d\n", i + info.Total_Patch_to_mesh[tm], info.No_Control_point_ON_ELEMENT[i + info.Total_Patch_to_mesh[tm]]);
	}

	// 節点座標
	for (i = 0; i < Total_Control_Point; i++)
	{
		fscanf(fp, "%d", &temp_i);
		for (j = 0; j < DIMENSION + 1; j++)
		{
			fscanf(fp, "%lf", &temp_d);
			info.Node_Coordinate[(temp_i + info.Total_Control_Point_to_mesh[tm]) * (DIMENSION + 1) + j] = temp_d;
		}
	}
	for (i = 0; i < Total_Control_Point; i++)
	{
		for (j = 0; j < DIMENSION + 1; j++)
		{
			// コントロールポイント座標・重みの新たな配列(for s-IGA/NewtonLaphson) DIMENSION == 2 の場合のみ記述
			if (j == 0)
			{
				info.Control_Coord_x[i + info.Total_Control_Point_to_mesh[tm]] = info.Node_Coordinate[(i + info.Total_Control_Point_to_mesh[tm]) * (DIMENSION + 1) + j];
			}
			else if (j == 1)
			{
				info.Control_Coord_y[i + info.Total_Control_Point_to_mesh[tm]] = info.Node_Coordinate[(i + info.Total_Control_Point_to_mesh[tm]) * (DIMENSION + 1) + j];
			}
			else if (j == DIMENSION)
			{
				info.Control_Weight[i + info.Total_Control_Point_to_mesh[tm]] = info.Node_Coordinate[(i + info.Total_Control_Point_to_mesh[tm]) * (DIMENSION + 1) + j];
			}
		}
	}
	fgets(s, 256, fp);

	// 拘束
	for (i = 0; i < Total_Constraint; i++)
	{
		fscanf(fp, "%d %d %lf",
			   &info.Constraint_Node_Dir[(i + info.Total_Constraint_to_mesh[tm]) * 2 + 0],
			   &info.Constraint_Node_Dir[(i + info.Total_Constraint_to_mesh[tm]) * 2 + 1],
			   &info.Value_of_Constraint[i + info.Total_Constraint_to_mesh[tm]]);
		info.Constraint_Node_Dir[(i + info.Total_Constraint_to_mesh[tm]) * 2 + 0] = info.Constraint_Node_Dir[(i + info.Total_Constraint_to_mesh[tm]) * 2 + 0] + info.Total_Control_Point_to_mesh[tm];

		printf("info.Constraint_Node_Dir[%d] = %d info.Constraint_Node_Dir[%d] = %d info.Value_of_Constraint[%d] = %le \n",
			   (i + info.Total_Constraint_to_mesh[tm]) * 2 + 0, info.Constraint_Node_Dir[(i + info.Total_Constraint_to_mesh[tm]) * 2 + 0],
			   (i + info.Total_Constraint_to_mesh[tm]) * 2 + 1, info.Constraint_Node_Dir[(i + info.Total_Constraint_to_mesh[tm]) * 2 + 1],
			   i + info.Total_Constraint_to_mesh[tm], info.Value_of_Constraint[i + info.Total_Constraint_to_mesh[tm]]);
	}
	fgets(s, 256, fp);

	// 荷重
	for (i = 0; i < Total_Load; i++)
	{
		fscanf(fp, "%d %d %lf",
			   &info.Load_Node_Dir[(i + info.Total_Load_to_mesh[tm]) * 2 + 0],
			   &info.Load_Node_Dir[(i + info.Total_Load_to_mesh[tm]) * 2 + 1],
			   &info.Value_of_Load[i + info.Total_Load_to_mesh[tm]]);
		info.Load_Node_Dir[(i + info.Total_Load_to_mesh[tm]) * 2 + 0] = info.Load_Node_Dir[(i + info.Total_Load_to_mesh[tm]) * 2 + 0] + info.Total_Control_Point_to_mesh[tm];

		printf("info.Load_Node_Dir[%d]= %d info.Load_Node_Dir[%d]= %d info.Value_of_Load[%d] = %le\n",
			   (i + info.Total_Load_to_mesh[tm]) * 2 + 0, info.Load_Node_Dir[(i + info.Total_Load_to_mesh[tm]) * 2 + 0],
			   (i + info.Total_Load_to_mesh[tm]) * 2 + 1, info.Load_Node_Dir[(i + info.Total_Load_to_mesh[tm]) * 2 + 1],
			   i + info.Total_Load_to_mesh[tm], info.Value_of_Load[i + info.Total_Load_to_mesh[tm]]);
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
		info.type_load_array[i + info.Total_DistributeForce_to_mesh[tm]] = type_load;
		info.iPatch_array[i + info.Total_DistributeForce_to_mesh[tm]] = iPatch;
		info.iCoord_array[i + info.Total_DistributeForce_to_mesh[tm]] = iCoord;
		info.val_Coord_array[i + info.Total_DistributeForce_to_mesh[tm]] = val_Coord;
		info.Range_Coord_array[(i + info.Total_DistributeForce_to_mesh[tm]) * 2 + 0] = Range_Coord[0];
		info.Range_Coord_array[(i + info.Total_DistributeForce_to_mesh[tm]) * 2 + 1] = Range_Coord[1];
		info.Coeff_Dist_Load_array[(i + info.Total_DistributeForce_to_mesh[tm]) * 3 + 0] = Coeff_Dist_Load[0];
		info.Coeff_Dist_Load_array[(i + info.Total_DistributeForce_to_mesh[tm]) * 3 + 1] = Coeff_Dist_Load[1];
		info.Coeff_Dist_Load_array[(i + info.Total_DistributeForce_to_mesh[tm]) * 3 + 2] = Coeff_Dist_Load[2];
	}
	fclose(fp);
}


// INC 等の作成
void Make_INC(int tm, information info)
{
	// info.INC の計算(節点番号をξ, ηの番号で表す為の配列)
	for (tm = 0; tm < Total_mesh; tm++)
	{
		int b, B, e, h, i, j, k, l, n, p, q, x, y, ii, jj, kkk, iiloc, jjloc, r = 0;
		int type_load, iPatch, iCoord;
		double val_Coord, Range_Coord[2], Coeff_Dist_Load[3];
		int No_Patch = info.Total_Patch_on_mesh[tm];
		int Total_Patch_to_Now = info.Total_Patch_to_mesh[tm];
		int Total_Element = info.Total_Element_on_mesh[tm];
		int Total_Element_to_Now = info.Total_Element_to_mesh[tm];
		int Total_DistributeForce = info.Total_DistributeForce_to_mesh[tm + 1] - info.Total_DistributeForce_to_mesh[tm];
		if (DIMENSION == 2) // for s-IGA
		{
			e = 0;
			for (l = 0; l < No_Patch; l++)
			{
				i = 0;
				for (jj = 0; jj < info.No_Control_point[(l + Total_Patch_to_Now) * DIMENSION + 1]; jj++)
				{
					for (ii = 0; ii < info.No_Control_point[(l + Total_Patch_to_Now) * DIMENSION + 0]; ii++)
					{
						info.INC[info.Patch_Control_point[info.Total_Control_Point_to_patch[l + Total_Patch_to_Now] + i] * DIMENSION + 0] = ii;
						info.INC[info.Patch_Control_point[info.Total_Control_Point_to_patch[l + Total_Patch_to_Now] + i] * DIMENSION + 1] = jj;

						if (ii >= info.Order[(l + Total_Patch_to_Now) * DIMENSION + 0] && jj >= info.Order[(l + Total_Patch_to_Now) * DIMENSION + 1])
						{

							for (jjloc = 0; jjloc <= info.Order[(l + Total_Patch_to_Now) * DIMENSION + 1]; jjloc++)
							{
								for (iiloc = 0; iiloc <= info.Order[(l + Total_Patch_to_Now) * DIMENSION + 0]; iiloc++)
								{
									B = info.Patch_Control_point[info.Total_Control_Point_to_patch[l + Total_Patch_to_Now] + i - jjloc * info.No_Control_point[(l + Total_Patch_to_Now) * DIMENSION + 0] - iiloc];
									b = jjloc * (info.Order[(l + Total_Patch_to_Now) * DIMENSION + 0] + 1) + iiloc;
									info.Controlpoint_of_Element[(e + Total_Element_to_Now) * MAX_NO_CCpoint_ON_ELEMENT + b] = B;
								}
							}
							info.Element_patch[e + Total_Element_to_Now] = l + Total_Patch_to_Now;
							info.Element_mesh[e + Total_Element_to_Now] = tm;
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
		// 		for (kk = 0; kk < info.No_Control_point[l][2]; kk++)
		// 		{
		// 			for (jj = 0; jj < info.No_Control_point[l][1]; jj++)
		// 			{
		// 				for (ii = 0; ii < info.No_Control_point[l][0]; ii++)
		// 				{
		// 					info.INC[l][info.Patch_Control_point[l][i]][0] = ii;
		// 					info.INC[l][info.Patch_Control_point[l][i]][1] = jj;
		// 					info.INC[l][info.Patch_Control_point[l][i]][2] = kk;
		// 					if (ii >= info.Order[l][0] && jj >= info.Order[l][1] && kk >= info.Order[l][2])
		// 					{
		// 						for (kkloc = 0; kkloc < info.Order[l][2]; kkloc++)
		// 						{
		// 							for (jjloc = 0; jjloc <= info.Order[l][1]; jjloc++)
		// 							{
		// 								for (iiloc = 0; iiloc <= info.Order[l][0]; iiloc++)
		// 								{
		// 									B = info.Patch_Control_point[l][i - jjloc * info.No_Control_point[l][0] - iiloc];
		// 									b = jjloc * (info.Order[l][0] + 1) + iiloc;
		// 									info.Controlpoint_of_Element[e * MAX_NO_CCpoint_ON_ELEMENT + b] = B;
		// 								}
		// 							}
		// 						}
		// 						info.Element_patch[e] = l + Total_Patch_to_Now;
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
				info.line_No_real_element[(l + Total_Patch_to_Now) * DIMENSION + j] = 0;
			}
		}

		for (l = 0; l < No_Patch; l++)
		{
			for (j = 0; j < DIMENSION; j++)
			{
				info.line_No_Total_element[(l + Total_Patch_to_Now) * DIMENSION + j] = info.No_knot[(l + Total_Patch_to_Now) * DIMENSION + j] - 2 * info.Order[(l + Total_Patch_to_Now) * DIMENSION + j] - 1;

				for (kkk = info.Order[(l + Total_Patch_to_Now) * DIMENSION + j]; kkk < info.No_knot[(l + Total_Patch_to_Now) * DIMENSION + j] - info.Order[(l + Total_Patch_to_Now) * DIMENSION + j] - 1; kkk++)
				{
					info.difference[info.Total_Knot_to_patch_dim[(l + Total_Patch_to_Now) * DIMENSION + j] + kkk - info.Order[(l + Total_Patch_to_Now) * DIMENSION + j]]
						= info.Position_Knots[info.Total_Knot_to_patch_dim[(l + Total_Patch_to_Now) * DIMENSION + j] + kkk + 1] - info.Position_Knots[info.Total_Knot_to_patch_dim[(l + Total_Patch_to_Now) * DIMENSION + j] + kkk];
					if (info.difference[info.Total_Knot_to_patch_dim[(l + Total_Patch_to_Now) * DIMENSION + j] + kkk - info.Order[(l + Total_Patch_to_Now) * DIMENSION + j]] != 0)
					{
						info.line_No_real_element[(l + Total_Patch_to_Now) * DIMENSION + j]++;
					}
				}
			}
		}

		// 要素に行番号, 列番号をつける
		if (DIMENSION == 2)
		{
			for (h = 0; h < Total_Element; h++)
			{
				info.Total_element_all_ID[h] = 0;
			}

			i = 0;
			for (l = 0; l < No_Patch; l++)
			{
				for (y = 0; y < info.line_No_Total_element[(l + Total_Patch_to_Now) * DIMENSION + 1]; y++)
				{
					for (x = 0; x < info.line_No_Total_element[(l + Total_Patch_to_Now) * DIMENSION + 0]; x++)
					{
						info.ENC[(i + info.Total_Element_to_mesh[tm]) * DIMENSION + 0] = x;
						info.ENC[(i + info.Total_Element_to_mesh[tm]) * DIMENSION + 1] = y;
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

				for (k = 0; k < info.line_No_Total_element[(l + Total_Patch_to_Now) * DIMENSION + j]; k++)
				{
					if (info.difference[info.Total_Knot_to_patch_dim[(l + Total_Patch_to_Now) * DIMENSION + j] + k] != 0)
					{
						info.real_element_line[(l + Total_Patch_to_Now) * (info.Total_Element_to_mesh[Total_mesh] * DIMENSION) + e * DIMENSION + j] = k;
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
				for (p = 0; p < info.line_No_real_element[(info.Element_patch[n + Total_Element_to_Now]) * DIMENSION + 0]; p++)
				{
					if (info.ENC[(n + info.Total_Element_to_mesh[tm]) * DIMENSION + 0] == info.real_element_line[(info.Element_patch[n + Total_Element_to_Now]) * (info.Total_Element_to_mesh[Total_mesh] * DIMENSION) + p * DIMENSION + 0])
					{
						for (q = 0; q < info.line_No_real_element[(info.Element_patch[n + Total_Element_to_Now]) * DIMENSION + 1]; q++)
						{
							if (info.ENC[(n + info.Total_Element_to_mesh[tm]) * DIMENSION + 1] == info.real_element_line[(info.Element_patch[n + Total_Element_to_Now]) * (info.Total_Element_to_mesh[Total_mesh] * DIMENSION) + q * DIMENSION + 1])
							{
								info.Total_element_all_ID[n]++;
							}
						}
					}
				}

				// IDが1の要素に番号を振る
				if (info.Total_element_all_ID[n] == 1)
				{
					info.real_element[r + info.real_Total_Element_to_mesh[tm]] = n + Total_Element_to_Now;
					info.real_El_No_on_mesh[tm * info.Total_Element_to_mesh[Total_mesh] + r] = n + Total_Element_to_Now;
					r++;
				}
			}

			// for s-IGA real_Total_Elementの初期化
			int real_Total_Element = 0;

			for (l = 0; l < No_Patch; l++)
			{
				real_Total_Element += info.line_No_real_element[(l + Total_Patch_to_Now) * DIMENSION + 0] * info.line_No_real_element[(l + Total_Patch_to_Now) * DIMENSION + 1];
			}
			info.real_Total_Element_on_mesh[tm] = real_Total_Element;
			info.real_Total_Element_to_mesh[tm + 1] = info.real_Total_Element_to_mesh[tm] + real_Total_Element;
		}

		for (i = 0; i < Total_DistributeForce; i++)
		{
			type_load = info.type_load_array[i + info.Total_DistributeForce_to_mesh[tm]];
			iPatch = info.iPatch_array[i + info.Total_DistributeForce_to_mesh[tm]];
			iCoord = info.iCoord_array[i + info.Total_DistributeForce_to_mesh[tm]];
			val_Coord = info.val_Coord_array[i + info.Total_DistributeForce_to_mesh[tm]];
			Range_Coord[0] = info.Range_Coord_array[(i + info.Total_DistributeForce_to_mesh[tm]) * 2 + 0];
			Range_Coord[1] = info.Range_Coord_array[(i + info.Total_DistributeForce_to_mesh[tm]) * 2 + 1];
			Coeff_Dist_Load[0] = info.Coeff_Dist_Load_array[(i + info.Total_DistributeForce_to_mesh[tm]) * 3 + 0];
			Coeff_Dist_Load[1] = info.Coeff_Dist_Load_array[(i + info.Total_DistributeForce_to_mesh[tm]) * 3 + 1];
			Coeff_Dist_Load[2] = info.Coeff_Dist_Load_array[(i + info.Total_DistributeForce_to_mesh[tm]) * 3 + 2];

			printf("type_load:%d\tiPatch:%d\tiCoord:%d\tval_coord:%lf\t", type_load, iPatch, iCoord, val_Coord);
			printf("Range0:%lf\tRange1:%lf\t", Range_Coord[0], Range_Coord[1]);
			printf("Coeff0:%lf\n", Coeff_Dist_Load[0]);
			Setting_Dist_Load_2D(tm, iPatch, iCoord, val_Coord, Range_Coord, type_load, Coeff_Dist_Load, info);
		}
	}
}


// Distributed Load
void Setting_Dist_Load_2D(int mesh_n, int iPatch, int iCoord, double val_Coord, double *Range_Coord, int type_load, double *Coeff_Dist_Load, information info)
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
	int *No_Element_for_Integration = (int *)malloc(sizeof(int) * info.Total_Knot_to_mesh[Total_mesh]); // No_Element_for_Integration[MAX_N_KNOT]
 
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

	for (iii = info.Order[iPatch * DIMENSION + iCoord]; iii < info.No_knot[iPatch * DIMENSION + iCoord] - info.Order[iPatch * DIMENSION + iCoord] - 1; iii++)
	{
		double epsi = 0.00000000001;

		if (info.Position_Knots[info.Total_Knot_to_patch_dim[iPatch * DIMENSION + iCoord] + iii] - epsi <= Range_Coord[0])
			iPos[0] = iii;
		if (info.Position_Knots[info.Total_Knot_to_patch_dim[iPatch * DIMENSION + iCoord] + iii + 1] - epsi <= Range_Coord[1])
			iPos[1] = iii + 1;
	}
	iRange_ele[0] = iPos[0] - info.Order[iPatch * DIMENSION + iCoord];
	iRange_ele[1] = iPos[1] - info.Order[iPatch * DIMENSION + iCoord] - 1;
	printf("iPos[0] = %d iPos[1] = %d\n", iPos[0], iPos[1]);
	printf("iRange_ele[0] = %d  iRange_ele[1] = %d\n", iRange_ele[0], iRange_ele[1]);

	if (iPos[0] < 0 || iPos[1] < 0)
	{
		printf("Error (Stop) iPos[0] = %d iPos[1] = %d\n", iPos[0], iPos[1]);
		exit(0);
	}

	for (jjj = info.Order[iPatch * DIMENSION + jCoord]; jjj < info.No_knot[iPatch * DIMENSION + jCoord] - info.Order[iPatch * DIMENSION + jCoord] - 1; jjj++)
	{
		double epsi = 0.00000000001;

		if (info.Position_Knots[info.Total_Knot_to_patch_dim[iPatch * DIMENSION + jCoord] + jjj] - epsi <= val_Coord
			&& info.Position_Knots[info.Total_Knot_to_patch_dim[iPatch * DIMENSION + jCoord] + jjj + 1] + epsi > val_Coord)
		{
			jPos[0] = jjj;
			jPos[1] = jjj + 1;
			val_jCoord_Local = -1.0 + 2.0 * (val_Coord - info.Position_Knots[info.Total_Knot_to_patch_dim[iPatch * DIMENSION + jCoord] + jjj])
							 / (info.Position_Knots[info.Total_Knot_to_patch_dim[iPatch * DIMENSION + jCoord] + jjj + 1] - info.Position_Knots[info.Total_Knot_to_patch_dim[iPatch * DIMENSION + jCoord] + jjj]);
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
		iX = jPos[0] - info.Order[iPatch * DIMENSION + 0];
		for (iY = iPos[0] - info.Order[iPatch * DIMENSION + 1]; iY < iPos[1] - info.Order[iPatch * DIMENSION + 1]; iY++)
		{
			No_Element_for_Integration[iii] = SearchForElement(mesh_n, iPatch, iX, iY, info);
			printf("Check No_Element_for_Integration[%d] = %d\n", iii, No_Element_for_Integration[iii]);
			iii++;
		}
	}

	if (iCoord == 0)
	{
		iY = jPos[0] - info.Order[iPatch * DIMENSION + 1];
		for (iX = iPos[0] - info.Order[iPatch * DIMENSION + 0]; iX < iPos[1] - info.Order[iPatch * DIMENSION + 0]; iX++)
		{
			No_Element_for_Integration[iii] = SearchForElement(mesh_n, iPatch, iX, iY, info);
			printf("Check No_Element_for_Integration[%d] = %d\n", iii, No_Element_for_Integration[iii]);
			iii++;
		}
	}
	No_Element_For_Dist_Load = iii;

	// Book keeping finished

	for (iii = 0; iii < No_Element_For_Dist_Load; iii++)
	{
		if (info.Total_element_all_ID[No_Element_for_Integration[iii]] == 1)
		{
			iX = info.ENC[No_Element_for_Integration[iii] * DIMENSION + 0];
			iY = info.ENC[No_Element_for_Integration[iii] * DIMENSION + 1];

			for (ic = 0; ic < (info.Order[iPatch * DIMENSION + 0] + 1) * (info.Order[iPatch * DIMENSION + 1] + 1); ic++)
				iControlpoint[ic] = info.Controlpoint_of_Element[No_Element_for_Integration[iii] * MAX_NO_CCpoint_ON_ELEMENT + ic];

			for (ig = 0; ig < NNG; ig++)
			{
				double Local_Coord[2], sfc, dxyzdge[3], detJ, XiEtaCoordParen, valDistLoad;
				int icc;
				Local_Coord[jCoord] = val_jCoord_Local;
				Local_Coord[iCoord] = GaussPt[ig];
				printf("ig = %d   Local_Coord[jCoord] = %f Local_Coord[iCoord] = %f\n", ig, Local_Coord[jCoord], Local_Coord[iCoord]);

				ShapeFunc_from_paren(Position_Data_param, Local_Coord, iCoord, No_Element_for_Integration[iii], info);
				XiEtaCoordParen = Position_Data_param[iCoord];
				printf("Check  Coeff_Dist_Load[0] = %f Coeff_Dist_Load[1] = %f  Coeff_Dist_Load[2] = %f  Position_Data_param[iCoord] = %f\n", Coeff_Dist_Load[0], Coeff_Dist_Load[1], Coeff_Dist_Load[2], Position_Data_param[iCoord]);
				valDistLoad = Coeff_Dist_Load[0] + Coeff_Dist_Load[1] * XiEtaCoordParen + Coeff_Dist_Load[2] * XiEtaCoordParen * XiEtaCoordParen;

				dxyzdge[0] = 0.0;
				dxyzdge[1] = 0.0;
				dxyzdge[2] = 0.0;
				for (icc = 0; icc < (info.Order[iPatch * DIMENSION + 0] + 1) * (info.Order[iPatch * DIMENSION + 1] + 1); icc++)
				{
					dxyzdge[0] += dShape_func(icc, iCoord, Local_Coord, No_Element_for_Integration[iii], info) * info.Node_Coordinate[iControlpoint[icc] * (DIMENSION + 1) + 0];
					dxyzdge[1] += dShape_func(icc, iCoord, Local_Coord, No_Element_for_Integration[iii], info) * info.Node_Coordinate[iControlpoint[icc] * (DIMENSION + 1) + 1];
				}

				detJ = sqrt(dxyzdge[0] * dxyzdge[0] + dxyzdge[1] * dxyzdge[1]);
				printf("Check the value of detJ etc: detJ = %f dxyzdge[0] = %f dxyzdge[1] = %f\n", detJ, dxyzdge[0], dxyzdge[1]);
				if (type_load < 2)
				{
					for (ic = 0; ic < (info.Order[iPatch * DIMENSION + 0] + 1) * (info.Order[iPatch * DIMENSION + 1] + 1); ic++)
					{
						sfc = Shape_func(ic, Local_Coord, No_Element_for_Integration[iii], info);
						info.Equivalent_Nodal_Force[iControlpoint[ic] * DIMENSION + type_load] += valDistLoad * sfc * detJ * Weight[ig];
					}
				}

				if (type_load == 2)
				{
					double LoadDir[2];
					LoadDir[0] = dxyzdge[1] / detJ;
					LoadDir[1] = -dxyzdge[0] / detJ;
					for (ic = 0; ic < (info.Order[iPatch * DIMENSION + 0] + 1) * (info.Order[iPatch * DIMENSION + 1] + 1); ic++)
					{
						sfc = Shape_func(ic, Local_Coord, No_Element_for_Integration[iii], info);
						info.Equivalent_Nodal_Force[iControlpoint[ic] * DIMENSION + 0] += LoadDir[0] * valDistLoad * sfc * detJ * Weight[ig];
						info.Equivalent_Nodal_Force[iControlpoint[ic] * DIMENSION + 1] += LoadDir[1] * valDistLoad * sfc * detJ * Weight[ig];
						printf("info.Equivalent_Nodal_Force[%d][0]=%le\nEquivalent_Nodal_Force[%d][1]=%le\n",
							   iControlpoint[ic], info.Equivalent_Nodal_Force[iControlpoint[ic] * DIMENSION + 0],
							   iControlpoint[ic], info.Equivalent_Nodal_Force[iControlpoint[ic] * DIMENSION + 1]);
					}
				}
				if (type_load == 3)
				{
					double LoadDir[2];
					LoadDir[0] = dxyzdge[0] / detJ;
					LoadDir[1] = dxyzdge[1] / detJ;
					for (ic = 0; ic < (info.Order[iPatch * DIMENSION + 0] + 1) * (info.Order[iPatch * DIMENSION + 1] + 1); ic++)
					{
						sfc = Shape_func(ic, Local_Coord, No_Element_for_Integration[iii], info);
						info.Equivalent_Nodal_Force[iControlpoint[ic] * DIMENSION + 0] += LoadDir[0] * valDistLoad * sfc * detJ * Weight[ig];
						info.Equivalent_Nodal_Force[iControlpoint[ic] * DIMENSION + 1] += LoadDir[1] * valDistLoad * sfc * detJ * Weight[ig];
					}
				}
			}
		}
	}
	free(No_Element_for_Integration);
}


int SearchForElement(int mesh_n, int iPatch, int iX, int iY, information info)
{
	int iii;

	for (iii = 0; iii < info.Total_Element_on_mesh[mesh_n]; iii++)
	{
		if (info.Element_patch[iii + info.Total_Element_to_mesh[mesh_n]] == iPatch)
		{
			if (iX == info.ENC[(iii + info.Total_Element_to_mesh[mesh_n]) * DIMENSION + 0] && iY == info.ENC[(iii + info.Total_Element_to_mesh[mesh_n]) * DIMENSION + 1])
				goto loopend;
		}
	}
	loopend:

	return (iii);
}


// for s_IGA, coupled matrix を求める, 要素の重なりを要素のガウス点から求める
void Check_coupled_Glo_Loc_element_for_Gauss(int mesh_n_over, int mesh_n_org, information info)
{
	int re, e;
	int i, j, k, m;
	int l, ll;
	int n_elements_over_point[POW_Ng_extended];
	int MAX_NNLOVER = 0;

	int *temp_element_n = (int *)malloc(sizeof(int) * MAX_N_ELEMENT_OVER_POINT);
	int *element_n_point = (int *)malloc(sizeof(int) * MAX_N_ELEMENT_OVER_ELEMENT);
	int *Check_coupled_No = (int *)malloc(sizeof(int) * MAX_N_ELEMENT_OVER);

	for (i = 0; i < MAX_N_ELEMENT_OVER; i++)
	{
		Check_coupled_No[i] = 0;
	}

	for (m = 0; m < 2; m++) // 最初 Ng 個のガウス点で重なりを求め, info.NNLOVER[e] >= 2 の e に対して, 再度 Ng_extended 個のガウス点で重なりを求める
	{
		Make_gauss_array(m);

		// グローバルパッチの Preprocessing 作成
		if (m == 0)
		{
			for (re = 0; re < info.real_Total_Element_on_mesh[mesh_n_org]; re++)
			{
				e = info.real_element[re + info.real_Total_Element_to_mesh[mesh_n_org]];
				Preprocessing(m, e, info);
			}
		}

		// ローカルパッチ(mesh_n_over)各要素の頂点の物理座標算出
		for (re = 0; re < info.real_Total_Element_on_mesh[mesh_n_over]; re++)
		{
			int i_gg, i_ee;
			double output_para[DIMENSION];
			int Total_n_elements;

			e = info.real_element[re + info.real_Total_Element_to_mesh[mesh_n_over]];

			if (m == 0 || (m == 1 && info.NNLOVER[e] >= 2))
			{
				Preprocessing(m, e, info);

				if (m == 1)
				{
					info.NNLOVER[e] = 0;
					for (i = 0; i < info.NNLOVER[e]; i++)
					{
						info.NELOVER[e * MAX_N_ELEMENT_OVER + i] = 0;
					}
				}

				Total_n_elements = 0;
				k = 0;
				ll = 0;

				// ローカルパッチ各要素のガウス点の物理座標のグローバルパッチでの(xi, eta)算出
				for (i = 0; i < info.Total_Patch_on_mesh[mesh_n_org]; i++)
				{
					// グローバルパッチ i での各方向ノットベクトル
					double *temp_Position_Knots_xi = (double *)malloc(sizeof(double) * info.No_knot[i * DIMENSION + 0]);
					double *temp_Position_Knots_eta = (double *)malloc(sizeof(double) * info.No_knot[i * DIMENSION + 1]);

					for (j = 0; j < info.No_knot[i * DIMENSION + 0]; j++)
					{
						temp_Position_Knots_xi[j] = info.Position_Knots[info.Total_Knot_to_patch_dim[i * DIMENSION + 0] + j];
					}
					for (j = 0; j < info.No_knot[i * DIMENSION + 1]; j++)
					{
						temp_Position_Knots_eta[j] = info.Position_Knots[info.Total_Knot_to_patch_dim[i * DIMENSION + 1] + j];
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
									data_result_shape[l] = info.Gauss_Coordinate[e * GP_2D * DIMENSION + g_n * DIMENSION + l];
								}
								else if (m == 1)
								{
									data_result_shape[l] = info.Gauss_Coordinate_ex[e * GP_2D * DIMENSION + g_n * DIMENSION + l];
								}
							}
							int itr_n = Calc_xi_eta(data_result_shape[0], data_result_shape[1],
													temp_Position_Knots_xi, temp_Position_Knots_eta,
													info.No_Control_point[i * DIMENSION + 0], info.No_Control_point[i * DIMENSION + 1], info.Order[i * DIMENSION + 0], info.Order[i * DIMENSION + 1],
													&output_para[0], &output_para[1], info);
							
							if (m == 0)
							{
								info.Loc_parameter_on_Glo[e * GP_2D * DIMENSION + g_n * DIMENSION + 0] = output_para[0];
								info.Loc_parameter_on_Glo[e * GP_2D * DIMENSION + g_n * DIMENSION + 1] = output_para[1];
							}
							else if (m == 1)
							{
								info.Loc_parameter_on_Glo_ex[e * GP_2D * DIMENSION + g_n * DIMENSION + 0] = output_para[0];
								info.Loc_parameter_on_Glo_ex[e * GP_2D * DIMENSION + g_n * DIMENSION + 1] = output_para[1];
							}

							// Newton Laphsonによって出力されたxi,etaから重なる要素を求める
							n_elements_over_point[k] = ele_check(i, output_para, temp_element_n, info);
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
				info.NNLOVER[e] = duplicate_delete(Total_n_elements, e, element_n_point, info); // info.NNLOVER: 要素 e に重なる要素の総数
			}
		}
	}

	for (re = 0; re < info.real_Total_Element_on_mesh[mesh_n_over]; re++)
	{
		e = info.real_element[re + info.real_Total_Element_to_mesh[mesh_n_over]];

		Check_coupled_No[info.NNLOVER[e]]++;

		if (MAX_NNLOVER < info.NNLOVER[e])
		{
			MAX_NNLOVER = info.NNLOVER[e];
		}
	}

	// 重なっている要素の割合
	printf("MAX_NNLOVER = %d\n", MAX_NNLOVER);
	for (i = 0; i <= MAX_NNLOVER; i++)
	{
		double Percent_Check_coupled_No = (double)Check_coupled_No[i] * 100.0 / (double)info.real_Total_Element_on_mesh[mesh_n_over];
		printf("Check_coupled_No[%d] = %d\t\t%.2lf %%\n", i, Check_coupled_No[i], Percent_Check_coupled_No);
	}
	
	free(element_n_point), free(temp_element_n), free(Check_coupled_No);
}


void Make_Loc_Glo(information info)
{
	int i, j, k;
	int jj;
	int e;
	int count;

	int j_n = info.real_Total_Element_to_mesh[Total_mesh] - info.real_Total_Element_on_mesh[0];

	for (i = 0; i < info.real_Total_Element_on_mesh[0]; i++)
	{
		e = info.real_element[i];
		count = 0;

		for (j = 0; j < j_n; j++)
		{
			jj = info.real_element[info.real_Total_Element_to_mesh[1] + j]; //ローカルメッシュ上のreal element番号

			if (info.NNLOVER[jj] > 0)
			{
				for (k = 0; k < info.NNLOVER[jj]; k++)
				{
					if (info.NELOVER[jj * MAX_N_ELEMENT_OVER + k] == e)
					{
						info.NELOVER[e * MAX_N_ELEMENT_OVER + count] = jj;
						count++;
					}
				}
			}
		}
		info.NNLOVER[e] = count;
	}
}


// Newton Raphsonによって出力されたxi,etaから重なる要素を求める
int ele_check(int patch_n, double *para_coord, int *temp_element_n, information info)
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

		for (k = 0; k < info.No_Control_point[patch_n * DIMENSION + j] - 1; k++)
		{
			if (RangeCheck_flag == 1)
				break;

			// Local要素の頂点がGlobalパッチ内にない場合
			if (para_coord[j] < info.Position_Knots[info.Total_Knot_to_patch_dim[patch_n * DIMENSION + j] + 0] || para_coord[j] > info.Position_Knots[info.Total_Knot_to_patch_dim[patch_n * DIMENSION + j] + info.No_knot[patch_n * DIMENSION + j] - 1])
			{
				RangeCheck_flag++;
			}

			// Local要素の頂点がGlobal要素内部にある場合
			if (para_coord[j] < info.Position_Knots[info.Total_Knot_to_patch_dim[patch_n * DIMENSION + j] + info.Order[patch_n * DIMENSION + j] + k])
			{
				int kk = 0;
				for (kk = 0; kk < k + 1; kk++)
				{
					if (para_coord[j] > info.Position_Knots[info.Total_Knot_to_patch_dim[patch_n * DIMENSION + j] + info.Order[patch_n * DIMENSION + j] + k - kk])
					{
						temp_ad[j][l] = k - kk;
						l++;
						RangeCheck_flag++;
						break;
					}
				}
			}

			// Local要素の頂点がGlobal要素境界上にある場合
			if (para_coord[j] == info.Position_Knots[info.Total_Knot_to_patch_dim[patch_n * DIMENSION + j] + info.Order[patch_n * DIMENSION + j] + k])
			{
				//頂点の座標がGlobalパッチの始点上にある場合
				if (para_coord[j] == info.Position_Knots[info.Total_Knot_to_patch_dim[patch_n * DIMENSION + j] + 0])
				{
					temp_ad[j][l] = k;
					l++;
					break;
				}
				//頂点の座標がGlobalパッチの終点上にある場合
				if (para_coord[j] == info.Position_Knots[info.Total_Knot_to_patch_dim[patch_n * DIMENSION + j] + info.No_knot[patch_n * DIMENSION + j] - 1])
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
				for (kk = 0; kk < info.Order[patch_n * DIMENSION + j]; kk++)
				{
					if (para_coord[j] == info.Position_Knots[info.Total_Knot_to_patch_dim[patch_n * DIMENSION + j] + info.Order[patch_n * DIMENSION + j] + k + kk + 1])
					// 多重ノット(次数分ループ)
					{
						printf("C0 continuity\n");
						temp_ad[j][l] = k + kk;
						l++;
					}
					if (para_coord[j] != info.Position_Knots[info.Total_Knot_to_patch_dim[patch_n * DIMENSION + j] + info.Order[patch_n * DIMENSION + j] + k + kk + 1])
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
				temp_element_n[l * No_line[0] + ll] = temp_ad[0][ll] + temp_ad[1][l] * info.line_No_Total_element[patch_n * DIMENSION + 0];
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
int duplicate_delete(int total, int element_n, int *element_n_point, information info)
{
	int i, j;

	j = 0;
	info.NELOVER[element_n * MAX_N_ELEMENT_OVER + j] = element_n_point[0];
	j++;
	for (i = 1; i < total; i++)
	{
		if (element_n_point[i] != element_n_point[i - 1])
		{
			info.NELOVER[element_n * MAX_N_ELEMENT_OVER + j] = element_n_point[i];
			j++;
		}
	}
	// j = 要素element_nに重なる要素の総数
	return j;
}


// Preprocessing
void Preprocessing(int m, int e, information info)
{
	double *a_matrix = (double *)malloc(sizeof(double) * POW_Ng_extended * DIMENSION * DIMENSION);					// a_matrix[POW_Ng_extended * DIMENSION * DIMENSION]
	double *dSF = (double *)malloc(sizeof(double) * POW_Ng * MAX_NO_CCpoint_ON_ELEMENT * DIMENSION);				// dSF[POW_Ng * MAX_NO_CCpoint_ON_ELEMENT * DIMENSION]
	double *dSF_ex = (double *)malloc(sizeof(double) * POW_Ng_extended * MAX_NO_CCpoint_ON_ELEMENT * DIMENSION);	// dSF_ex[POW_Ng_extended * MAX_NO_CCpoint_ON_ELEMENT * DIMENSION]

	// ガウス点の物理座標を計算
	Make_Gauss_Coordinate(m, e, info);

	// ガウス点でのヤコビアン, Bマトリックスを計算
	Make_dSF(m, e, dSF, dSF_ex, info);
	Make_Jac(m, e, dSF, dSF_ex, a_matrix, info);
	Make_B_Matrix(m, e, dSF, dSF_ex, a_matrix, info);

	free(a_matrix), free(dSF), free(dSF_ex);
}


void Make_Gauss_Coordinate(int m, int e, information info)
{
	int i, j, k;
	double temp_coordinate[DIMENSION];
	double R;

	for (i = 0; i < GP_2D; i++)
	{
		temp_coordinate[0] = Gxi[i][0];
		temp_coordinate[1] = Gxi[i][1];
		
		for (j = 0; j < info.No_Control_point_ON_ELEMENT[info.Element_patch[e]]; j++)
		{
			R = Shape_func(j, temp_coordinate, e, info);

			for (k = 0; k < DIMENSION; k++)
			{
				if (m == 0)
				{
					info.Gauss_Coordinate[e * GP_2D * DIMENSION + i * DIMENSION + k] += R * info.Node_Coordinate[info.Controlpoint_of_Element[e * MAX_NO_CCpoint_ON_ELEMENT + j] * (DIMENSION + 1) + k];
				}
				else if (m == 1)
				{
					info.Gauss_Coordinate_ex[e * GP_2D * DIMENSION + i * DIMENSION + k] += R * info.Node_Coordinate[info.Controlpoint_of_Element[e * MAX_NO_CCpoint_ON_ELEMENT + j] * (DIMENSION + 1) + k];
				}
			}
		}
	}
}


void Make_dSF(int m, int e, double *dSF, double *dSF_ex, information info)
{
	int i, j, k;
	double temp_coordinate[DIMENSION];

	for (i = 0; i < GP_2D; i++)
	{
		temp_coordinate[0] = Gxi[i][0];
		temp_coordinate[1] = Gxi[i][1];

		for (j = 0; j < info.No_Control_point_ON_ELEMENT[info.Element_patch[e]]; j++)
		{
			for (k = 0; k < DIMENSION; k++)
			{
				if (m == 0)
				{
					dSF[i * MAX_NO_CCpoint_ON_ELEMENT * DIMENSION + j * DIMENSION + k] = dShape_func(j, k, temp_coordinate, e, info);
				}
				else if (m == 1)
				{
					dSF_ex[i * MAX_NO_CCpoint_ON_ELEMENT * DIMENSION + j * DIMENSION + k] = dShape_func(j, k, temp_coordinate, e, info);
				}
			}
		}
	}
}


void Make_Jac(int m, int e, double *dSF, double *dSF_ex, double *a_matrix, information info)
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
				for (l = 0; l < info.No_Control_point_ON_ELEMENT[info.Element_patch[e]]; l++)
				{
					if (m == 0)
					{
						a[j][k] += dSF[i * MAX_NO_CCpoint_ON_ELEMENT * DIMENSION + l * DIMENSION + k] * info.Node_Coordinate[info.Controlpoint_of_Element[e * MAX_NO_CCpoint_ON_ELEMENT + l] * (DIMENSION + 1) + j];
					}
					else if (m == 1)
					{
						a[j][k] += dSF_ex[i * MAX_NO_CCpoint_ON_ELEMENT * DIMENSION + l * DIMENSION + k] * info.Node_Coordinate[info.Controlpoint_of_Element[e * MAX_NO_CCpoint_ON_ELEMENT + l] * (DIMENSION + 1) + j];
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
			info.Jac[e * GP_2D + i] = J;
		}
		else if (m == 1)
		{
			info.Jac_ex[e * GP_2D + i] = J;
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


void Make_B_Matrix(int m, int e, double *dSF, double *dSF_ex, double *a_matrix, information info)
{
	int i, j, k, l;
	double *b = (double *)malloc(sizeof(double) * DIMENSION * MAX_NO_CCpoint_ON_ELEMENT);	// b[DIMENSION][MAX_NO_CCpoint_ON_ELEMENT]

	for (i = 0; i < GP_2D; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			for (k = 0; k < info.No_Control_point_ON_ELEMENT[info.Element_patch[e]]; k++)
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
			for (j = 0; j < info.No_Control_point_ON_ELEMENT[info.Element_patch[e]]; j++)
			{
				info.B_Matrix[e * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 0 * MAX_KIEL_SIZE + (2 * j)] = b[0 * MAX_NO_CCpoint_ON_ELEMENT + j];
				info.B_Matrix[e * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 0 * MAX_KIEL_SIZE + (2 * j + 1)] = 0.0;
				info.B_Matrix[e * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 1 * MAX_KIEL_SIZE + (2 * j)] = 0.0;
				info.B_Matrix[e * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 1 * MAX_KIEL_SIZE + (2 * j + 1)] = b[1 * MAX_NO_CCpoint_ON_ELEMENT + j];
				info.B_Matrix[e * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 2 * MAX_KIEL_SIZE + (2 * j)] = b[1 * MAX_NO_CCpoint_ON_ELEMENT + j];
				info.B_Matrix[e * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 2 * MAX_KIEL_SIZE + (2 * j + 1)] = b[0 * MAX_NO_CCpoint_ON_ELEMENT + j];
			}
		}
		else if (m == 1)
		{
			for (j = 0; j < info.No_Control_point_ON_ELEMENT[info.Element_patch[e]]; j++)
			{
				info.B_Matrix_ex[e * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 0 * MAX_KIEL_SIZE + (2 * j)] = b[0 * MAX_NO_CCpoint_ON_ELEMENT + j];
				info.B_Matrix_ex[e * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 0 * MAX_KIEL_SIZE + (2 * j + 1)] = 0.0;
				info.B_Matrix_ex[e * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 1 * MAX_KIEL_SIZE + (2 * j)] = 0.0;
				info.B_Matrix_ex[e * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 1 * MAX_KIEL_SIZE + (2 * j + 1)] = b[1 * MAX_NO_CCpoint_ON_ELEMENT + j];
				info.B_Matrix_ex[e * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 2 * MAX_KIEL_SIZE + (2 * j)] = b[1 * MAX_NO_CCpoint_ON_ELEMENT + j];
				info.B_Matrix_ex[e * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 2 * MAX_KIEL_SIZE + (2 * j + 1)] = b[0 * MAX_NO_CCpoint_ON_ELEMENT + j];
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
// D Matrix
void Make_D_Matrix(information info)
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
					info.D[i * D_MATRIX_SIZE + j] = D1[i][j];
		}
		else if (DM == 1) // 平面ひずみ状態(2Dの場合はこっち)
		{
			double Eone = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
			double D1[3][3] = {{Eone, nu / (1.0 - nu) * Eone, 0}, {nu / (1.0 - nu) * Eone, Eone, 0}, {0, 0, (1 - 2 * nu) / 2 / (1.0 - nu) * Eone}};

			for (i = 0; i < D_MATRIX_SIZE; i++)
				for (j = 0; j < D_MATRIX_SIZE; j++)
					info.D[i * D_MATRIX_SIZE + j] = D1[i][j];
		}
	}
}


// 拘束されている行数を省いた行列の番号の制作
void Make_Index_Dof(information info)
{
	int i, k = 0;

	// 拘束されている自由度(Degree Of free)をERRORにする
	for (i = 0; i < info.Total_Constraint_to_mesh[Total_mesh]; i++)
	{
		info.Index_Dof[info.Constraint_Node_Dir[i * 2 + 0] * DIMENSION + info.Constraint_Node_Dir[i * 2 + 1]] = ERROR;
	}
	// ERROR以外に番号を付ける
	for (i = 0; i < info.Total_Control_Point_to_mesh[Total_mesh] * DIMENSION; i++)
	{
		if (info.Index_Dof[i] != ERROR)
		{
			info.Index_Dof[i] = k;
			k++;
		}
	}

	K_Whole_Size = k;
	printf("K_Whole_Size=%d\n", k);
}


void Make_K_Whole_Ptr_Col(information info)
{
	int i, j, k, ii, jj;
	int N, NE, i_index, j_index;

	// 初期化
	for (i = 0; i < K_Whole_Size + 1; i++)
		info.K_Whole_Ptr[i] = 0;

	// 大きく分割するためのループ
	for (N = 0; N < info.Total_Control_Point_to_mesh[Total_mesh]; N += K_DIVISION_LENGE)
	{
		// 各節点に接する節点を取得
		for (i = 0; i < K_DIVISION_LENGE; i++)
		{
			info.Total_Control_Point_To_Node[i] = 0;
		}
		for (i = 0; i < info.Total_Element_to_mesh[Total_mesh]; i++)
		{
			for (ii = 0; ii < info.No_Control_point_ON_ELEMENT[info.Element_patch[i]]; ii++)
			{
				NE = info.Controlpoint_of_Element[i * MAX_NO_CCpoint_ON_ELEMENT + ii] - N;
				if (0 <= NE && NE < K_DIVISION_LENGE)
				{
					for (j = 0; j < info.No_Control_point_ON_ELEMENT[info.Element_patch[i]]; j++) //ローカル要素
					{
						// 数字がない時
						if (info.Total_Control_Point_To_Node[NE] == 0)
						{
							// 節点番号を取得
							info.Node_To_Node[NE * info.Total_Control_Point_to_mesh[Total_mesh] + 0] = info.Controlpoint_of_Element[i * MAX_NO_CCpoint_ON_ELEMENT + j];
							info.Total_Control_Point_To_Node[NE]++;
						}
						// 同じものがあったら
						// k > 0 以降の取得
						// kのカウント
						for (k = 0; k < info.Total_Control_Point_To_Node[NE]; k++)
						{
							if (info.Node_To_Node[NE * info.Total_Control_Point_to_mesh[Total_mesh] + k] == info.Controlpoint_of_Element[i * MAX_NO_CCpoint_ON_ELEMENT + j])
							{
								break;
							}
						}
						// 未設定のNode_To_Node取得
						if (k == info.Total_Control_Point_To_Node[NE])
						{
							info.Node_To_Node[NE * info.Total_Control_Point_to_mesh[Total_mesh] + k] = info.Controlpoint_of_Element[i * MAX_NO_CCpoint_ON_ELEMENT + j];
							info.Total_Control_Point_To_Node[NE]++;
						}
					}
					// 別メッシュとの重なりを考慮
					if (info.NNLOVER[i] > 0)
					{
						for (jj = 0; jj < info.NNLOVER[i]; jj++)
						{
							for (j = 0; j < info.No_Control_point_ON_ELEMENT[info.Element_patch[info.NELOVER[i * MAX_N_ELEMENT_OVER + jj]]]; j++) //ローカル要素
							{
								// 数字がない時
								if (info.Total_Control_Point_To_Node[NE] == 0)
								{
									// 節点番号を取得
									info.Node_To_Node[NE * info.Total_Control_Point_to_mesh[Total_mesh] + 0] = info.Controlpoint_of_Element[info.NELOVER[i * MAX_N_ELEMENT_OVER + jj] * MAX_NO_CCpoint_ON_ELEMENT + j];
									info.Total_Control_Point_To_Node[NE]++;
								}

								// 同じものがあったら
								// k > 0 以降の取得
								// kのカウント
								for (k = 0; k < info.Total_Control_Point_To_Node[NE]; k++)
								{
									if (info.Node_To_Node[NE * info.Total_Control_Point_to_mesh[Total_mesh] + k] == info.Controlpoint_of_Element[info.NELOVER[i * MAX_N_ELEMENT_OVER + jj] * MAX_NO_CCpoint_ON_ELEMENT + j])
									{
										break;
									}
								}
								// 未設定のNode_To_Node取得
								if (k == info.Total_Control_Point_To_Node[NE])
								{
									info.Node_To_Node[NE * info.Total_Control_Point_to_mesh[Total_mesh] + k] = info.Controlpoint_of_Element[info.NELOVER[i * MAX_N_ELEMENT_OVER + jj] * MAX_NO_CCpoint_ON_ELEMENT + j];
									info.Total_Control_Point_To_Node[NE]++;
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
			if (N + i < info.Total_Control_Point_to_mesh[Total_mesh])
			{
				for (j = 0; j < info.Total_Control_Point_To_Node[i]; j++)
				{
					int Min = info.Node_To_Node[i * info.Total_Control_Point_to_mesh[Total_mesh] + j], No = j;
					for (k = j; k < info.Total_Control_Point_To_Node[i]; k++)
					{
						if (Min > info.Node_To_Node[i * info.Total_Control_Point_to_mesh[Total_mesh] + k])
						{
							Min = info.Node_To_Node[i * info.Total_Control_Point_to_mesh[Total_mesh] + k];
							No = k;
						}
					}
					for (k = No; k > j; k--)
					{
						info.Node_To_Node[i * info.Total_Control_Point_to_mesh[Total_mesh] + k] = info.Node_To_Node[i * info.Total_Control_Point_to_mesh[Total_mesh] + k - 1];
					}
					info.Node_To_Node[i * info.Total_Control_Point_to_mesh[Total_mesh] + j] = Min;
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
				if (N + i < info.Total_Control_Point_to_mesh[Total_mesh])
				{
					i_index = info.Index_Dof[(N + i) * DIMENSION + ii];
					k = 0;
					if (i_index >= 0)
					{
						info.K_Whole_Ptr[i_index + 1] = info.K_Whole_Ptr[i_index];
						for (j = 0; j < info.Total_Control_Point_To_Node[i]; j++)
						{
							for (jj = 0; jj < DIMENSION; jj++)
							{
								j_index = info.Index_Dof[info.Node_To_Node[i * info.Total_Control_Point_to_mesh[Total_mesh] + j] * DIMENSION + jj];
								if (j_index >= 0 && j_index >= i_index)
								{
									info.K_Whole_Ptr[i_index + 1]++;
									info.K_Whole_Col[info.K_Whole_Ptr[i_index] + k] = j_index;
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
void Make_K_Whole_Val(information info)
{
	int i, j, j1, j2, k1, k2, l;
	int a, b, re;

	double *K_EL = (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * MAX_KIEL_SIZE);			// K_EL[MAX_KIEL_SIZE][MAX_KIEL_SIZE]
	double *coupled_K_EL= (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * MAX_KIEL_SIZE);		// coupled_K_EL[MAX_KIEL_SIZE][MAX_KIEL_SIZE]

	for (re = 0; re < info.real_Total_Element_to_mesh[Total_mesh]; re++)
	{
		i = info.real_element[re];

		if (info.Element_mesh[i] == 0 && re == 0) // 2つめの条件は効率化のため
		{
			Make_gauss_array(0);
		}
		else if (info.Element_mesh[i] > 0)
		{
			if (info.NNLOVER[i] == 1 && (info.NNLOVER[info.real_element[re - 1]] != 1 || info.Element_mesh[info.real_element[re - 1]] == 0)) // 2つめ以降の条件は効率化のため
			{
				Make_gauss_array(0);
			}
			else if (info.NNLOVER[i] >= 2 && (info.NNLOVER[info.real_element[re - 1]] == 1 || info.Element_mesh[info.real_element[re - 1]] == 0)) // 2つめ以降の条件は効率化のため
			{
				Make_gauss_array(1);
			}
		}

		// 各要素のK_ELを求める
		Make_K_EL(i, K_EL, info);

		// Valを求める
		for (j1 = 0; j1 < info.No_Control_point_ON_ELEMENT[info.Element_patch[i]]; j1++)
		{
			for (j2 = 0; j2 < DIMENSION; j2++)
			{
				a = info.Index_Dof[info.Controlpoint_of_Element[i * MAX_NO_CCpoint_ON_ELEMENT + j1] * DIMENSION + j2];
				if (a >= 0)
				{
					for (k1 = 0; k1 < info.No_Control_point_ON_ELEMENT[info.Element_patch[i]]; k1++)
					{
						for (k2 = 0; k2 < DIMENSION; k2++)
						{
							b = info.Index_Dof[info.Controlpoint_of_Element[i * MAX_NO_CCpoint_ON_ELEMENT + k1] * DIMENSION + k2];
							if (b >= 0 && b >= a)
							{
								for (l = info.K_Whole_Ptr[a]; l < info.K_Whole_Ptr[a + 1]; l++)
								{
									if (info.K_Whole_Col[l] == b)
									{
										info.K_Whole_Val[l] += K_EL[(j1 * DIMENSION + j2) * MAX_KIEL_SIZE + k1 * DIMENSION + k2];
										break;
									}
								}
							}
						}
					}
				}
			}
		}

		if (info.Element_mesh[i] > 0 && info.NNLOVER[i] > 0) // ローカルメッシュ上の要素について, 重なっている要素が存在するとき
		{
			for (j = 0; j < info.NNLOVER[i]; j++)
			{
				// 各要素のcoupled_K_ELを求める
				Make_coupled_K_EL(i, info.NELOVER[i * MAX_N_ELEMENT_OVER + j], coupled_K_EL, info);

				// Valを求める
				for (j1 = 0; j1 < info.No_Control_point_ON_ELEMENT[info.Element_patch[info.NELOVER[i * MAX_N_ELEMENT_OVER + j]]]; j1++)
				{
					for (j2 = 0; j2 < DIMENSION; j2++)
					{
						a = info.Index_Dof[info.Controlpoint_of_Element[info.NELOVER[i * MAX_N_ELEMENT_OVER + j] * MAX_NO_CCpoint_ON_ELEMENT + j1] * DIMENSION + j2];
						if (a >= 0)
						{
							for (k1 = 0; k1 < info.No_Control_point_ON_ELEMENT[info.Element_patch[i]]; k1++)
							{
								for (k2 = 0; k2 < DIMENSION; k2++)
								{
									b = info.Index_Dof[info.Controlpoint_of_Element[i * MAX_NO_CCpoint_ON_ELEMENT + k1] * DIMENSION + k2];
									if (b >= 0 && b >= a)
									{
										for (l = info.K_Whole_Ptr[a]; l < info.K_Whole_Ptr[a + 1]; l++)
										{
											if (info.K_Whole_Col[l] == b)
											{
												info.K_Whole_Val[l] += coupled_K_EL[(j1 * DIMENSION + j2) * MAX_KIEL_SIZE + k1 * DIMENSION + k2];
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
void Make_K_EL(int El_No, double *K_EL, information info)
{
	int i, j, k, l;
	
	int KIEL_SIZE = info.No_Control_point_ON_ELEMENT[info.Element_patch[El_No]] * DIMENSION;

	double *B = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * MAX_KIEL_SIZE); // B[D_MATRIX_SIZE][MAX_KIEL_SIZE];
	double *K1 = (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * MAX_KIEL_SIZE); // K1[MAX_KIEL_SIZE][MAX_KIEL_SIZE];
	double J = 0;

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
			J = info.Jac[El_No * GP_2D + i];
			for (j = 0; j < D_MATRIX_SIZE; j++)
			{
				for (k = 0; k < MAX_KIEL_SIZE; k++)
				{
					B[j * MAX_KIEL_SIZE + k] = info.B_Matrix[El_No * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + j * MAX_KIEL_SIZE + k];
				}
			}
		}
		else if (GP_2D == POW_Ng_extended)
		{
			J = info.Jac_ex[El_No * GP_2D + i];
			for (j = 0; j < D_MATRIX_SIZE; j++)
			{
				for (k = 0; k < MAX_KIEL_SIZE; k++)
				{
					B[j * MAX_KIEL_SIZE + k] = info.B_Matrix_ex[El_No * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + j * MAX_KIEL_SIZE + k];
				}
			}
		}

		BDBJ(KIEL_SIZE, B, J, K1, info);
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
void Make_coupled_K_EL(int El_No_loc, int El_No_glo, double *coupled_K_EL, information info)
{
	int i, j, k, l;

	int BDBJ_flag;
	int KIEL_SIZE = info.No_Control_point_ON_ELEMENT[info.Element_patch[El_No_glo]] * DIMENSION;

	double *B = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * MAX_KIEL_SIZE); // B[D_MATRIX_SIZE][MAX_KIEL_SIZE];
	double *BG = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * MAX_KIEL_SIZE); // BG[D_MATRIX_SIZE][MAX_KIEL_SIZE];
	double *K1 = (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * MAX_KIEL_SIZE); // K1[MAX_KIEL_SIZE][MAX_KIEL_SIZE];
	double J = 0;

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
		double para[DIMENSION] = {0.0};
		double G_Gxi[DIMENSION] = {0.0};

		if (GP_2D == POW_Ng)
		{
			J = info.Jac[El_No_loc * GP_2D + i];
			for (j = 0; j < D_MATRIX_SIZE; j++)
			{
				for (k = 0; k < MAX_KIEL_SIZE; k++)
				{
					B[j * MAX_KIEL_SIZE + k] = info.B_Matrix[El_No_loc * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + j * MAX_KIEL_SIZE + k];
				}
			}
			para[0] = info.Loc_parameter_on_Glo[El_No_loc * GP_2D * DIMENSION + i * DIMENSION + 0];
			para[1] = info.Loc_parameter_on_Glo[El_No_loc * GP_2D * DIMENSION + i * DIMENSION + 1];
		}
		else if (GP_2D == POW_Ng_extended)
		{
			J = info.Jac_ex[El_No_loc * GP_2D + i];
			for (j = 0; j < D_MATRIX_SIZE; j++)
			{
				for (k = 0; k < MAX_KIEL_SIZE; k++)
				{
					B[j * MAX_KIEL_SIZE + k] = info.B_Matrix_ex[El_No_loc * GP_2D * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + j * MAX_KIEL_SIZE + k];
				}
			}
			para[0] = info.Loc_parameter_on_Glo_ex[El_No_loc * GP_2D * DIMENSION + i * DIMENSION + 0];
			para[1] = info.Loc_parameter_on_Glo_ex[El_No_loc * GP_2D * DIMENSION + i * DIMENSION + 1];
		}

		// 要素内外判定
		if (para[0] >= info.Position_Knots[info.Total_Knot_to_patch_dim[0 * DIMENSION + 0] + info.Order[0 * DIMENSION + 0] + info.ENC[El_No_glo * DIMENSION + 0]] &&
			para[0] <  info.Position_Knots[info.Total_Knot_to_patch_dim[0 * DIMENSION + 0] + info.Order[0 * DIMENSION + 0] + info.ENC[El_No_glo * DIMENSION + 0] + 1] &&
			para[1] >= info.Position_Knots[info.Total_Knot_to_patch_dim[0 * DIMENSION + 1] + info.Order[0 * DIMENSION + 1] + info.ENC[El_No_glo * DIMENSION + 1]] &&
			para[1] <  info.Position_Knots[info.Total_Knot_to_patch_dim[0 * DIMENSION + 1] + info.Order[0 * DIMENSION + 1] + info.ENC[El_No_glo * DIMENSION + 1] + 1])
		{
			BDBJ_flag = 1;

			// 親要素座標の算出
			G_Gxi[0] = -1.0 + 2.0 * (para[0] - info.Position_Knots[info.Total_Knot_to_patch_dim[0 * DIMENSION + 0] + info.Order[0 * DIMENSION + 0] + info.ENC[El_No_glo * DIMENSION + 0]])
						/ (info.Position_Knots[info.Total_Knot_to_patch_dim[0 * DIMENSION + 0] + info.Order[0 * DIMENSION + 0] + info.ENC[El_No_glo * DIMENSION + 0] + 1] - info.Position_Knots[info.Total_Knot_to_patch_dim[0 * DIMENSION + 0] + info.Order[0 * DIMENSION + 0] + info.ENC[El_No_glo * DIMENSION + 0]]);
			G_Gxi[1] = -1.0 + 2.0 * (para[1] - info.Position_Knots[info.Total_Knot_to_patch_dim[0 * DIMENSION + 1] + info.Order[0 * DIMENSION + 1] + info.ENC[El_No_glo * DIMENSION + 1]])
						/ (info.Position_Knots[info.Total_Knot_to_patch_dim[0 * DIMENSION + 1] + info.Order[0 * DIMENSION + 1] + info.ENC[El_No_glo * DIMENSION + 1] + 1] - info.Position_Knots[info.Total_Knot_to_patch_dim[0 * DIMENSION + 1] + info.Order[0 * DIMENSION + 1] + info.ENC[El_No_glo * DIMENSION + 1]]);
		}
		else
		{
			BDBJ_flag = 0;
		}

		// 要素内であるとき, 結合要素剛性マトリックス計算
		if (BDBJ_flag)
		{
			// 重なるグローバル要素のBマトリックス
			Make_BG_Matrix(El_No_glo, BG, G_Gxi, info);
			// BGTDBLJの計算
			coupled_BDBJ(KIEL_SIZE, B, BG, J, K1, info);
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
void Make_BG_Matrix(int El_No, double *BG, double *Local_coord, information info)
{
	double a[DIMENSION][DIMENSION];
	double *b = (double *)malloc(sizeof(double) * DIMENSION * MAX_NO_CCpoint_ON_ELEMENT); // b[DIMENSION][MAX_NO_CCpoint_ON_ELEMENT]

	int i, j, k;

	for (i = 0; i < DIMENSION; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			a[i][j] = 0.0;
			for (k = 0; k < info.No_Control_point_ON_ELEMENT[info.Element_patch[El_No]]; k++)
			{
				a[i][j] += dShape_func(k, j, Local_coord, El_No, info) * info.Node_Coordinate[info.Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + k] * (DIMENSION + 1) + i];
			}
		}
	}

	InverseMatrix_2x2(a);

	for (i = 0; i < DIMENSION; i++)
	{
		for (j = 0; j < info.No_Control_point_ON_ELEMENT[info.Element_patch[El_No]]; j++)
		{
			b[i * MAX_NO_CCpoint_ON_ELEMENT + j] = 0.0;
			for (k = 0; k < DIMENSION; k++)
			{
				b[i * MAX_NO_CCpoint_ON_ELEMENT + j] += a[k][i] * dShape_func(j, k, Local_coord, El_No, info);
			}
		}
	}

	// 2次元
	for (i = 0; i < info.No_Control_point_ON_ELEMENT[info.Element_patch[El_No]]; i++)
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
void BDBJ(int KIEL_SIZE, double *B, double J, double *K_EL, information info)
{
	int i, j, k;
	double *BD = (double *)calloc(MAX_KIEL_SIZE * D_MATRIX_SIZE, sizeof(double)); // BD[MAX_KIEL_SIZE][D_MATRIX_SIZE];

	// [B]T[info.D][B]Jの計算
	for (i = 0; i < KIEL_SIZE; i++)
	{
		for (j = 0; j < D_MATRIX_SIZE; j++)
		{
			for (k = 0; k < D_MATRIX_SIZE; k++)
			{
				BD[i * D_MATRIX_SIZE + j] += B[k * MAX_KIEL_SIZE + i] * info.D[k * D_MATRIX_SIZE + j];
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
void coupled_BDBJ(int KIEL_SIZE, double *B, double *BG, double J, double *K_EL, information info)
{
	int i, j, k;
	double *BD = (double *)calloc(MAX_KIEL_SIZE * D_MATRIX_SIZE, sizeof(double)); // BD[MAX_KIEL_SIZE][D_MATRIX_SIZE];

	//[B]T[info.D][B]Jの計算
	for (i = 0; i < KIEL_SIZE; i++)
	{
		for (j = 0; j < D_MATRIX_SIZE; j++)
		{
			for (k = 0; k < D_MATRIX_SIZE; k++)
			{
				BD[i * D_MATRIX_SIZE + j] += BG[k * MAX_KIEL_SIZE + i] * info.D[k * D_MATRIX_SIZE + j];
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
void Make_F_Vec(information info)
{
	int i, index;
	for (i = 0; i < K_Whole_Size; i++)
		info.rhs_vec[i] = 0.0;
	for (i = 0; i < info.Total_Load_to_mesh[Total_mesh]; i++)
	{
		index = info.Index_Dof[info.Load_Node_Dir[i * 2 + 0] * DIMENSION + info.Load_Node_Dir[i * 2 + 1]];
		if (index >= 0)
			info.rhs_vec[index] += info.Value_of_Load[i];
	}
}


// 強制変位対策
void Make_F_Vec_disp_const(information info)
{
	int ie, idir, inode, jdir, jnode, kk_const;
	int ii, iii, b, bb, jj, ii_local, jj_local;

	int i;

	double *K_EL = (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * MAX_KIEL_SIZE);			// K_EL[MAX_KIEL_SIZE][MAX_KIEL_SIZE]

	Make_gauss_array(0);

	for (ie = 0; ie < info.real_Total_Element_to_mesh[Total_mesh]; ie++)
	{
		i = info.real_element[ie];

		iii = 0;
		for (idir = 0; idir < DIMENSION; idir++)
		{
			for (inode = 0; inode < info.No_Control_point_ON_ELEMENT[info.Element_patch[i]]; inode++)
			{
				b = info.Index_Dof[info.Controlpoint_of_Element[i * MAX_NO_CCpoint_ON_ELEMENT + inode] * DIMENSION + idir];
				if (b < 0)
					iii++;
			}
		}

		if (iii > 0)
		{
			Make_K_EL(i, K_EL, info);
			for (idir = 0; idir < DIMENSION; idir++)
			{
				for (inode = 0; inode < info.No_Control_point_ON_ELEMENT[info.Element_patch[i]]; inode++)
				{
					// ii = info.Controlpoint_of_Element[info.real_El_No_on_mesh[Total_mesh * info.Total_Element_to_mesh[Total_mesh] + ie] * MAX_NO_CCpoint_ON_ELEMENT + inode] * DIMENSION + idir;
					ii = info.Controlpoint_of_Element[info.real_element[ie] * MAX_NO_CCpoint_ON_ELEMENT + inode] * DIMENSION + idir;
					// ↑ real_El_No_on_mesh を real_element に変更
					b = info.Index_Dof[ii];
					if (b >= 0)
					{
						ii_local = inode * DIMENSION + idir;
						for (jdir = 0; jdir < DIMENSION; jdir++)
						{
							for (jnode = 0; jnode < info.No_Control_point_ON_ELEMENT[info.Element_patch[i]]; jnode++)
							{
								jj = info.Controlpoint_of_Element[i * MAX_NO_CCpoint_ON_ELEMENT + jnode] * DIMENSION + jdir;
								bb = info.Index_Dof[jj];
								if (bb < 0)
								{
									jj_local = jnode * DIMENSION + jdir;
									for (kk_const = 0; kk_const < info.Total_Constraint_to_mesh[Total_mesh]; kk_const++)
									{
										if (info.Controlpoint_of_Element[i * MAX_NO_CCpoint_ON_ELEMENT + jnode] == info.Constraint_Node_Dir[kk_const * 2 + 0] && jdir == info.Constraint_Node_Dir[kk_const * 2 + 1])
										{
											info.rhs_vec[b] -= K_EL[ii_local * MAX_KIEL_SIZE + jj_local] * info.Value_of_Constraint[kk_const];
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
void Add_Equivalent_Nodal_Force_to_F_Vec(information info)
{
	int i, j, index;
	for (j = 0; j < DIMENSION; j++)
	{
		for (i = 0; i < info.Total_Control_Point_to_mesh[Total_mesh]; i++)
		{
			index = info.Index_Dof[i * DIMENSION + j];
			if (index >= 0)
			{
				info.rhs_vec[index] += info.Equivalent_Nodal_Force[i * DIMENSION + j];
			}
		}
	}
}


// PCG solver, 前処理付共役勾配法により[K]{d}={f}を解く
void PCG_Solver(int max_itetarion, double eps, information info)
{
	int i, j, k;
	int ndof = K_Whole_Size;

	double *r = (double *)malloc(sizeof(double) * ndof);
	double *p = (double *)calloc(ndof, sizeof(double));
	double *y = (double *)malloc(sizeof(double) * ndof);
	double *r2 = (double *)calloc(ndof, sizeof(double));

	double *gg = (double *)malloc(sizeof(double) * MAX_K_WHOLE_SIZE); // gg[MAX_K_WHOLE_SIZE]
	double *dd = (double *)malloc(sizeof(double) * MAX_K_WHOLE_SIZE); // dd[MAX_K_WHOLE_SIZE]
	double *pp = (double *)malloc(sizeof(double) * MAX_K_WHOLE_SIZE); // pp[MAX_K_WHOLE_SIZE]

	// 初期化
	for (i = 0; i < ndof; i++)
		info.sol_vec[i] = 0.0;

	// 前処理行列作成
	double *M = (double *)malloc(sizeof(double) * MAX_NON_ZERO);
	int *M_Ptr = (int *)malloc(sizeof(int) * MAX_K_WHOLE_SIZE + 1);
	int *M_Col = (int *)malloc(sizeof(int) * MAX_NON_ZERO);
	Make_M(M, M_Ptr, M_Col, ndof, info);

	// 第0近似解に対する残差の計算
	double *ax = (double *)calloc(ndof, sizeof(double));
	for (i = 0; i < ndof; i++)
	{
		for (j = info.K_Whole_Ptr[i]; j < info.K_Whole_Ptr[i + 1]; j++)
		{
			ax[i] += info.K_Whole_Val[j] * info.sol_vec[info.K_Whole_Col[j]];
			if (i != info.K_Whole_Col[j])
			{
				ax[info.K_Whole_Col[j]] += info.K_Whole_Val[j] * info.sol_vec[i];
			}
		}
	}
	for (i = 0; i < ndof; i++)
	{
		r[i] = info.rhs_vec[i] - ax[i];
	}
	free(ax);

	// 第0近似解に対する残差の計算
	// for (i = 0; i < ndof; i++)
	// {
	// 	r[i] = info.rhs_vec[i];
	// }

	// p_0 = (LDL^T)^-1 r_0 の計算 <- CG法で M = [[K^G, 0], [0, K^L]] とし,p_0 = (LDL^T)^-1 r_0 = M^-1 r_0
	CG(ndof, p, M, M_Ptr, M_Col, r, gg, dd, pp);

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
					temp1 = RowCol_to_icount(i, j, info); // temp_array_K[i][j]
				}
				else if (i > j)
				{
					temp1 = RowCol_to_icount(j, i, info); // temp_array_K[i][j] = temp_array_K[j][i]
				}

				if (temp1 != -1)
				{
					temp_array_K[j] = info.K_Whole_Val[temp1];
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
			info.sol_vec[i] += alpha * p[i];
			r[i] -= alpha * y[i];
		}

		// (r*r)_(k+1)の計算
		CG(ndof, r2, M, M_Ptr, M_Col, r, gg, dd, pp);

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
			e2 += info.sol_vec[i] * info.sol_vec[i];
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
		beta = - inner_product(ndof, y, r2) / temp_scaler;

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
	free(gg), free(dd), free(pp);
}


void Make_M(double *M, int *M_Ptr, int *M_Col, int ndof, information info)
{
	int i, j;
	int ndof_glo = 0;

	// グローバルパッチのdofを求める
	for (i = 0; i < info.Total_Control_Point_on_mesh[0] * DIMENSION; i++)
	{
		if (info.Index_Dof[i] != ERROR)
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

		for (j = info.K_Whole_Ptr[i]; j < info.K_Whole_Ptr[i + 1]; j++)
		{
			if (i < ndof_glo && info.K_Whole_Col[j] < ndof_glo)
			{
				M[counter] = info.K_Whole_Val[j];
				M_Col[counter] = info.K_Whole_Col[j];
				counter++;
				M_Ptr[i + 1]++;
			}
			else if (i >= ndof_glo)
			{
				M[counter] = info.K_Whole_Val[j];
				M_Col[counter] = info.K_Whole_Col[j];
				counter++;
				M_Ptr[i + 1]++;
			}
		}
	}
}


void CG(int ndof, double *solution_vec, double *M, int *M_Ptr, int *M_Col, double *right_vec, double *gg, double *dd, double *pp)
{
	// CG solver
	// double *gg = (double *)malloc(sizeof(double) * MAX_K_WHOLE_SIZE); // gg[MAX_K_WHOLE_SIZE]
	// double *dd = (double *)malloc(sizeof(double) * MAX_K_WHOLE_SIZE); // dd[MAX_K_WHOLE_SIZE]
	// double *pp = (double *)malloc(sizeof(double) * MAX_K_WHOLE_SIZE); // pp[MAX_K_WHOLE_SIZE]
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
	// free(gg), free(dd), free(pp);
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


int RowCol_to_icount(int row, int col, information info)
{
	for (int j = info.K_Whole_Ptr[row]; j < info.K_Whole_Ptr[row + 1]; j++)
	{
		if (info.K_Whole_Col[j] == col)
		{
			return j;
		}
		else if (info.K_Whole_Col[j] > col)
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
double Shape_func(int I_No, double *Local_coord, int El_No, information info)
{
	int i, j;
	double R, weight_func = 0.0;

	double Position_Data_param[DIMENSION];

	double *shape_func = (double *)malloc(sizeof(double) * MAX_NO_CCpoint_ON_ELEMENT);
	double *Shape = (double *)malloc(sizeof(double) * DIMENSION * info.Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER); // Shape[DIMENSION][MAX_N_NODE][10]
	double *dShape = (double *)malloc(sizeof(double) * DIMENSION * info.Total_Control_Point_to_mesh[Total_mesh]); // dShape[DIMENSION][MAX_N_NODE]

	for (i = 0; i < MAX_NO_CCpoint_ON_ELEMENT; i++)
	{
		shape_func[i] = 1.0;
	}

	for (i = 0; i < info.No_Control_point_ON_ELEMENT[info.Element_patch[El_No]]; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			ShapeFunc_from_paren(Position_Data_param, Local_coord, j, El_No, info);
			ShapeFunction1D(Position_Data_param, j, El_No, Shape, dShape, info);
			shape_func[i] *= Shape[j * (info.Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info.INC[info.Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + j] * MAX_ORDER + info.Order[info.Element_patch[El_No] * DIMENSION + j]];
		}
		weight_func += shape_func[i] * info.Node_Coordinate[info.Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * (DIMENSION + 1) + DIMENSION];
	}

	free(Shape), free(dShape);

	if (I_No < info.No_Control_point_ON_ELEMENT[info.Element_patch[El_No]])
		R = shape_func[I_No] * info.Node_Coordinate[info.Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + I_No] * (DIMENSION + 1) + DIMENSION] / weight_func;

	else
		R = ERROR;

	free(shape_func);
	return R;
}


void ShapeFunc_from_paren(double *Position_Data_param, double *Local_coord, int j, int e, information info)
{
	int i = 0;

	i = info.INC[info.Controlpoint_of_Element[e * MAX_NO_CCpoint_ON_ELEMENT + 0] * DIMENSION + j];
	Position_Data_param[j] = ((info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + i + 1] - info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + i]) * Local_coord[j] + (info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + i + 1] + info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + i])) / 2.0;
}


void ShapeFunction1D(double *Position_Data_param, int j, int e, double *Shape, double *dShape, information info)
{
	int ii;
	int p;

	for (ii = 0; ii < info.No_knot[info.Element_patch[e] * DIMENSION + j]; ii++)
	{
		if (info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + ii] == info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + ii + 1])
		{
			Shape[j * (info.Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + 0] = 0.0;
		}
		else if (info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + ii] != info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + ii + 1] && info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + ii] <= Position_Data_param[j] && Position_Data_param[j] < info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + ii + 1])
		{
			Shape[j * (info.Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + 0] = 1.0;
		}
		else if (info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + ii] != info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + ii + 1] && info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + ii + 1] == info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + (info.No_knot[info.Element_patch[e] * DIMENSION + j] - 1)] && info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + ii] <= Position_Data_param[j] && Position_Data_param[j] <= info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + ii + 1])
		{
			Shape[j * (info.Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + 0] = 1.0;
		}
		else
		{
			Shape[j * (info.Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + 0] = 0.0;
		}
	}

	for (ii = 0; ii < info.No_knot[info.Element_patch[e] * DIMENSION + j]; ii++)
	{
		for (p = 1; p <= info.Order[info.Element_patch[e] * DIMENSION + j]; p++)
		{
			Shape[j * (info.Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + p] = 0.0;
		}
	}

	double left_term, right_term;
	for (p = 1; p <= info.Order[info.Element_patch[e] * DIMENSION + j]; p++)
	{
		for (ii = 0; ii < info.No_knot[info.Element_patch[e] * DIMENSION + j]; ii++)
		{
			left_term = 0.0;
			right_term = 0.0;

			if ((Position_Data_param[j] - info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + ii]) * Shape[j * (info.Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + p - 1] == 0 && info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + ii + p] - info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + ii] == 0)
				left_term = 0.0;
			else
			{
				left_term = (Position_Data_param[j] - info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + ii]) / (info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + ii + p] - info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + ii]) * Shape[j * (info.Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + p - 1];
			}
			if ((info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + ii + p + 1] - Position_Data_param[j]) * Shape[j * (info.Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + (ii + 1) * MAX_ORDER + p - 1] == 0 && info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + ii + p + 1] - info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + ii + 1] == 0)
				right_term = 0.0;
			else
			{
				right_term = (info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + ii + p + 1] - Position_Data_param[j]) / (info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + ii + p + 1] - info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + ii + 1]) * Shape[j * (info.Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + (ii + 1) * MAX_ORDER + p - 1];
			}
			Shape[j * (info.Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + p] = left_term + right_term;
		}
	}

	double dleft_term, dright_term;
	for (ii = 0; ii < info.No_Control_point[info.Element_patch[e] * DIMENSION + j] + 1; ii++)
	{
		dleft_term = 0.0;
		dright_term = 0.0;

		if (info.Order[info.Element_patch[e] * DIMENSION + j] * Shape[j * (info.Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + info.Order[info.Element_patch[e] * DIMENSION + j] - 1] == 0 && info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + ii + info.Order[info.Element_patch[e] * DIMENSION + j]] - info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + ii] == 0)
			dleft_term = 0.0;
		else
			dleft_term = info.Order[info.Element_patch[e] * DIMENSION + j] / (info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + ii + info.Order[info.Element_patch[e] * DIMENSION + j]] - info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + ii]) * Shape[j * (info.Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + info.Order[info.Element_patch[e] * DIMENSION + j] - 1];

		if (info.Order[info.Element_patch[e] * DIMENSION + j] * Shape[j * (info.Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + (ii + 1) * MAX_ORDER + info.Order[info.Element_patch[e] * DIMENSION + j] - 1] == 0 && info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + ii + info.Order[info.Element_patch[e] * DIMENSION + j] + 1] - info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + ii + 1] == 0)
			dright_term = 0.0;
		else
			dright_term = info.Order[info.Element_patch[e] * DIMENSION + j] / (info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + ii + info.Order[info.Element_patch[e] * DIMENSION + j] + 1] - info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + ii + 1]) * Shape[j * (info.Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + (ii + 1) * MAX_ORDER + info.Order[info.Element_patch[e] * DIMENSION + j] - 1];

		dShape[j * info.Total_Control_Point_to_mesh[Total_mesh] + ii] = dleft_term - dright_term;
	}
}


double dShape_func(int I_No, int xez, double *Local_coord, int El_No, information info)
{
	double dR;

	double *dShape_func1 = (double *)malloc(sizeof(double) * info.Total_Control_Point_to_mesh[Total_mesh]); // dShape_func1[MAX_N_NODE];
	double *dShape_func2 = (double *)malloc(sizeof(double) * info.Total_Control_Point_to_mesh[Total_mesh]); // dShape_func2[MAX_N_NODE];

	NURBS_deriv(Local_coord, El_No, dShape_func1, dShape_func2, info);

	if (xez != 0 && xez != 1)
		dR = ERROR;

	else if (I_No < info.No_Control_point_ON_ELEMENT[info.Element_patch[El_No]])
	{
		if (xez == 0)
		{
			dR = dShape_func1[info.Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + I_No]] * dShapeFunc_from_paren(xez, El_No, info);
		}
		else if (xez == 1)
		{
			dR = dShape_func2[info.Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + I_No]] * dShapeFunc_from_paren(xez, El_No, info);
		}
	}
	else
		dR = ERROR;

	free(dShape_func1), free(dShape_func2);
	return dR;
}


void NURBS_deriv(double *Local_coord, int El_No, double *dShape_func1, double *dShape_func2, information info)
{
	int i, j;
	double weight_func = 0.0;
	double dWeight_func1 = 0.0;
	double dWeight_func2 = 0.0;

	double Position_Data_param[DIMENSION];

	double *shape_func = (double *)malloc(sizeof(double) * MAX_NO_CCpoint_ON_ELEMENT);
	double *Shape = (double *)malloc(sizeof(double) * DIMENSION * info.Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER); // Shape[DIMENSION][MAX_N_NODE][10]
	double *dShape = (double *)malloc(sizeof(double) * DIMENSION * info.Total_Control_Point_to_mesh[Total_mesh]); // dShape[DIMENSION][MAX_N_NODE]

	for (i = 0; i < MAX_NO_CCpoint_ON_ELEMENT; i++)
	{
		shape_func[i] = 1.0;
	}

	for (i = 0; i < info.No_Control_point_ON_ELEMENT[info.Element_patch[El_No]]; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			ShapeFunc_from_paren(Position_Data_param, Local_coord, j, El_No, info);
			ShapeFunction1D(Position_Data_param, j, El_No, Shape, dShape, info);
			shape_func[i] *= Shape[j * (info.Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info.INC[info.Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + j] * MAX_ORDER + info.Order[info.Element_patch[El_No] * DIMENSION + j]];
		}
		weight_func += shape_func[i] * info.Node_Coordinate[info.Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * (DIMENSION + 1) + DIMENSION];
	}
	for (i = 0; i < info.No_Control_point_ON_ELEMENT[info.Element_patch[El_No]]; i++)
	{
		dWeight_func1 += dShape[0 * info.Total_Control_Point_to_mesh[Total_mesh] + info.INC[info.Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + 0]] * Shape[1 * (info.Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info.INC[info.Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + 1] * MAX_ORDER + info.Order[info.Element_patch[El_No] * DIMENSION + 1]] * info.Node_Coordinate[info.Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * (DIMENSION + 1) + DIMENSION];
		dWeight_func2 += Shape[0 * (info.Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info.INC[info.Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + 0] * MAX_ORDER + info.Order[info.Element_patch[El_No] * DIMENSION + 0]] * dShape[1 * info.Total_Control_Point_to_mesh[Total_mesh] + info.INC[info.Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + 1]] * info.Node_Coordinate[info.Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * (DIMENSION + 1) + DIMENSION];
	}
	for (i = 0; i < info.No_Control_point_ON_ELEMENT[info.Element_patch[El_No]]; i++)
	{
		dShape_func1[info.Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i]] = info.Node_Coordinate[info.Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * (DIMENSION + 1) + DIMENSION] * (weight_func * dShape[0 * info.Total_Control_Point_to_mesh[Total_mesh] + info.INC[info.Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + 0]] * Shape[1 * (info.Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info.INC[info.Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + 1] * MAX_ORDER + info.Order[info.Element_patch[El_No] * DIMENSION + 1]] - dWeight_func1 * shape_func[i]) / (weight_func * weight_func);
		dShape_func2[info.Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i]] = info.Node_Coordinate[info.Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * (DIMENSION + 1) + DIMENSION] * (weight_func * Shape[0 * (info.Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info.INC[info.Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + 0] * MAX_ORDER + info.Order[info.Element_patch[El_No] * DIMENSION + 0]] * dShape[1 * info.Total_Control_Point_to_mesh[Total_mesh] + info.INC[info.Controlpoint_of_Element[El_No * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + 1]] - dWeight_func2 * shape_func[i]) / (weight_func * weight_func);
	}

	free(Shape), free(dShape), free(shape_func);
}


double dShapeFunc_from_paren(int j, int e, information info)
{
	int i;
	double dPosition_Data_param;

	i = info.INC[info.Controlpoint_of_Element[e * MAX_NO_CCpoint_ON_ELEMENT + 0] * DIMENSION + j];

	dPosition_Data_param = (info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + i + 1] - info.Position_Knots[info.Total_Knot_to_patch_dim[info.Element_patch[e] * DIMENSION + j] + i]) / 2.0;
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

			// for (int temp_i = 0; temp_i < info.No_knot[0][0]; temp_i++)
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
				double *output_xi, double *output_eta, information info)
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

	double *temp_Position_Knots_xi = (double *)malloc(sizeof(double) * info.No_knot[0 * DIMENSION + 0]);
	double *temp_Position_Knots_eta = (double *)malloc(sizeof(double) * info.No_knot[0 * DIMENSION + 1]);
	for (i = 0; i < info.No_knot[0 * DIMENSION + 0]; i++)
	{
		temp_Position_Knots_xi[i] = info.Position_Knots[info.Total_Knot_to_patch_dim[0 * DIMENSION + 0] + i];
	}
	for (i = 0; i < info.No_knot[0 * DIMENSION + 1]; i++)
	{
		temp_Position_Knots_eta[i] = info.Position_Knots[info.Total_Knot_to_patch_dim[0 * DIMENSION + 1] + i];
	}

	// 初期値の設定
	temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;

	for (i = 0; i < repeat; i++)
	{
		rNURBS_surface(temp_Position_Knots_xi, temp_Position_Knots_eta,
					   info.Control_Coord_x, info.Control_Coord_y,
					   info.No_Control_point[0 * DIMENSION + 0], info.No_Control_point[0 * DIMENSION + 1],
					   info.Control_Weight, info.Order[0 * DIMENSION + 0], info.Order[0 * DIMENSION + 1],
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
					   info.Control_Coord_x, info.Control_Coord_y,
					   info.No_Control_point[0 * DIMENSION + 0], info.No_Control_point[0 * DIMENSION + 1],
					   info.Control_Weight, info.Order[0 * DIMENSION + 0], info.Order[0 * DIMENSION + 1],
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
						info.Control_Coord_x, info.Control_Coord_y,
						info.No_Control_point[0 * DIMENSION + 0], info.No_Control_point[0 * DIMENSION + 1],
						info.Control_Weight, info.Order[0 * DIMENSION + 0], info.Order[0 * DIMENSION + 1],
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
						info.Control_Coord_x, info.Control_Coord_y,
						info.No_Control_point[0 * DIMENSION + 0], info.No_Control_point[0 * DIMENSION + 1],
						info.Control_Weight, info.Order[0 * DIMENSION + 0], info.Order[0 * DIMENSION + 1],
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


// Postprocessing
void Make_Displacement(information info)
{
	int i, j;

	for(i = 0; i < info.Total_Constraint_to_mesh[Total_mesh]; i++)
    {
		info.Displacement[info.Constraint_Node_Dir[i * 2 + 0] * DIMENSION + info.Constraint_Node_Dir[i * 2 + 1]] = info.Value_of_Constraint[i];
    }

	for(i = 0; i < info.Total_Control_Point_to_mesh[Total_mesh]; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			int index = info.Index_Dof[i * DIMENSION + j];
			if (index >= 0)
				info.Displacement[i * DIMENSION + j] = info.sol_vec[index];
		}
	}
}


// for s_IGA overlay
void s_IGA_overlay(information info)
{
	int i, j, k;
	int patch_n = info.Total_Patch_to_mesh[Total_mesh];
	int n_patch_glo = info.Total_Patch_to_mesh[1];

	// 一要素の分割数
	division_ele_xi = DIVISION_ELE;
	division_ele_eta = DIVISION_ELE;

	fp = fopen("view.dat", "w");
	fprintf(fp, "%d\t%d\t%d\n", 1, division_ele_xi, division_ele_eta);
	fclose(fp);

	// 重ね合わせ結果出力のためのoverlay_view.dat作成
	fp = fopen("overlay_view.dat", "w");
	fprintf(fp, "%d\t%d\t%d\n", 1, division_ele_xi, division_ele_eta);
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

		// パッチ i での各方向ノットベクトル, Coodinate, displacement
		double *temp_Position_Knots_xi = (double *)malloc(sizeof(double) * info.No_knot[i * DIMENSION + 0]);
		double *temp_Position_Knots_eta = (double *)malloc(sizeof(double) * info.No_knot[i * DIMENSION + 1]);
		double *temp_Coord_x = (double *)malloc(sizeof(double) * info.No_Control_point_in_patch[i]);
		double *temp_Coord_y = (double *)malloc(sizeof(double) * info.No_Control_point_in_patch[i]);
		double *temp_Coord_w = (double *)malloc(sizeof(double) * info.No_Control_point_in_patch[i]);
		double *temp_disp_x = (double *)malloc(sizeof(double) * info.No_Control_point_in_patch[i]);
		double *temp_disp_y = (double *)malloc(sizeof(double) * info.No_Control_point_in_patch[i]);

		for (j = 0; j < info.No_knot[i * DIMENSION + 0]; j++)
		{
			temp_Position_Knots_xi[j] = info.Position_Knots[info.Total_Knot_to_patch_dim[i * DIMENSION + 0] + j];
		}
		for (j = 0; j < info.No_knot[i * DIMENSION + 1]; j++)
		{
			temp_Position_Knots_eta[j] = info.Position_Knots[info.Total_Knot_to_patch_dim[i * DIMENSION + 1] + j];
		}

		for (j = 0; j < info.No_Control_point_in_patch[i]; j++)
		{
			temp_Coord_x[j] = info.Control_Coord_x[info.Total_Control_Point_to_patch[i] + j];
			temp_Coord_y[j] = info.Control_Coord_y[info.Total_Control_Point_to_patch[i] + j];
			temp_Coord_w[j] = info.Control_Weight[info.Total_Control_Point_to_patch[i] + j];
			temp_disp_x[j] = info.Displacement[(info.Total_Control_Point_to_patch[i] + j) * DIMENSION + 0];
			temp_disp_y[j] = info.Displacement[(info.Total_Control_Point_to_patch[i] + j) * DIMENSION + 1];
		}

		Calculation(info.Order[i * DIMENSION + 0], info.Order[i * DIMENSION + 1], info.No_knot[i * DIMENSION + 0], info.No_knot[i * DIMENSION + 1], info.No_Control_point[i * DIMENSION + 0], info.No_Control_point[i * DIMENSION + 1],
					temp_Position_Knots_xi, temp_Position_Knots_eta, temp_Coord_x, temp_Coord_y, temp_Coord_w,
					temp_disp_x, temp_disp_y, i, info);

		if (i >= n_patch_glo) // ローカル上のパッチに対しては重合計算行う
		{
			for (j = 0; j < n_patch_glo; j++)
			{
				// パッチ j での各方向ノットベクトル, Coodinate, displacement
				double *temp_Position_Knots_xi_glo = (double *)malloc(sizeof(double) * info.No_knot[j * DIMENSION + 0]);
				double *temp_Position_Knots_eta_glo = (double *)malloc(sizeof(double) * info.No_knot[j * DIMENSION + 1]);
				double *temp_Coord_x_glo = (double *)malloc(sizeof(double) * info.No_Control_point_in_patch[j]);
				double *temp_Coord_y_glo = (double *)malloc(sizeof(double) * info.No_Control_point_in_patch[j]);
				double *temp_Coord_w_glo = (double *)malloc(sizeof(double) * info.No_Control_point_in_patch[j]);
				double *temp_disp_x_glo = (double *)malloc(sizeof(double) * info.No_Control_point_in_patch[j]);
				double *temp_disp_y_glo = (double *)malloc(sizeof(double) * info.No_Control_point_in_patch[j]);

				for (k = 0; k < info.No_knot[j * DIMENSION + 0]; k++)
				{
					temp_Position_Knots_xi_glo[k] = info.Position_Knots[info.Total_Knot_to_patch_dim[j * DIMENSION + 0] + k];
				}
				for (k = 0; k < info.No_knot[j * DIMENSION + 1]; k++)
				{
					temp_Position_Knots_eta_glo[k] = info.Position_Knots[info.Total_Knot_to_patch_dim[j * DIMENSION + 1] + k];
				}

				for (k = 0; k < info.No_Control_point_in_patch[j]; k++)
				{
					temp_Coord_x_glo[k] = info.Control_Coord_x[info.Total_Control_Point_to_patch[j] + k];
					temp_Coord_y_glo[k] = info.Control_Coord_y[info.Total_Control_Point_to_patch[j] + k];
					temp_Coord_w_glo[k] = info.Control_Weight[info.Total_Control_Point_to_patch[j] + k];
					temp_disp_x_glo[k] = info.Displacement[(info.Total_Control_Point_to_patch[j] + k) * DIMENSION + 0];
					temp_disp_y_glo[k] = info.Displacement[(info.Total_Control_Point_to_patch[j] + k) * DIMENSION + 1];
				}

				Calculation_overlay(info.Order[i * DIMENSION + 0], info.Order[i * DIMENSION + 1],
									info.No_knot[i * DIMENSION + 0], info.No_knot[i * DIMENSION + 1],
									info.No_Control_point[i * DIMENSION + 0], info.No_Control_point[i * DIMENSION + 1],
									temp_Position_Knots_xi, temp_Position_Knots_eta,
									temp_Coord_x, temp_Coord_y,
									temp_Coord_w,
									info.Order[j * DIMENSION + 0], info.Order[j * DIMENSION + 1],
									info.No_Control_point[j * DIMENSION + 0], info.No_Control_point[j * DIMENSION + 1],
									temp_Position_Knots_xi_glo, temp_Position_Knots_eta_glo,	
									temp_Coord_x_glo, temp_Coord_y_glo,
									temp_disp_x_glo, temp_disp_y_glo,
									temp_Coord_w_glo, i, info);

				free(temp_Position_Knots_xi_glo), free(temp_Position_Knots_eta_glo);
				free(temp_Coord_x_glo), free(temp_Coord_y_glo), free(temp_Coord_w_glo), free(temp_disp_x_glo), free(temp_disp_y_glo);
			}
		}

		free(temp_Position_Knots_xi), free(temp_Position_Knots_eta);
		free(temp_Coord_x), free(temp_Coord_y), free(temp_Coord_w), free(temp_disp_x), free(temp_disp_y);
	}
}


void Calculation(int order_xi, int order_eta, int knot_n_xi, int knot_n_eta, int cntl_p_n_xi, int cntl_p_n_eta,
				 double *input_knot_vec_xi, double *input_knot_vec_eta, double *cntl_px, double *cntl_py, double *weight,
				 double *disp_cntl_px, double *disp_cntl_py, int patch_num, information info)
{
	int i, j, k, l;
	int ii, jj, kk, ll;
	int e, p, ep;
	double temp1, temp2, temp3;
	double temp_matrix[2][2];

	// ξ, ηの値決定と ∂ξ/∂チルダξ, ∂η/∂チルダη
	double *calc_xi = (double *)malloc(sizeof(double) * knot_n_xi * division_ele_xi);		// calc_xi[MAX_POINTS]
	double *calc_eta = (double *)malloc(sizeof(double) * knot_n_eta * division_ele_eta);	// calc_eta[MAX_POINTS]
	double *dtilda_xi = (double *)malloc(sizeof(double) * knot_n_xi * division_ele_xi);		// dtilda_xi[MAX_KNOTS]
	double *dtilda_eta = (double *)malloc(sizeof(double) * knot_n_eta * division_ele_eta);	// dtilda_eta[MAX_KNOTS]

	double *dxi_x = (double *)malloc(sizeof(double) * info.Total_Element_to_mesh[Total_mesh] * DIVISION_ELEMENT);		// dxi_x[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]
	double *dxi_y = (double *)malloc(sizeof(double) * info.Total_Element_to_mesh[Total_mesh] * DIVISION_ELEMENT);		// dxi_y[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]
	double *deta_x = (double *)malloc(sizeof(double) * info.Total_Element_to_mesh[Total_mesh] * DIVISION_ELEMENT);		// deta_x[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]
	double *deta_y = (double *)malloc(sizeof(double) * info.Total_Element_to_mesh[Total_mesh] * DIVISION_ELEMENT);		// deta_y[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]
	double *dxi_disp_x = (double *)malloc(sizeof(double) * info.Total_Element_to_mesh[Total_mesh] * DIVISION_ELEMENT);	// dxi_disp_x[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]
	double *dxi_disp_y = (double *)malloc(sizeof(double) * info.Total_Element_to_mesh[Total_mesh] * DIVISION_ELEMENT);	// dxi_disp_y[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]
	double *deta_disp_x = (double *)malloc(sizeof(double) * info.Total_Element_to_mesh[Total_mesh] * DIVISION_ELEMENT);	// deta_disp_x[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]
	double *deta_disp_y = (double *)malloc(sizeof(double) * info.Total_Element_to_mesh[Total_mesh] * DIVISION_ELEMENT);	// deta_disp_y[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]

	k = 0;
	l = 0;
	for (i = 0; i < knot_n_xi - 1; i++)
	{
		if (input_knot_vec_xi[i] != input_knot_vec_xi[i + 1])
		{
			calc_xi[k] = input_knot_vec_xi[i];
			dtilda_xi[l] = (input_knot_vec_xi[i + 1] - input_knot_vec_xi[i]) / 2.0;
			k++;
			l++;
			if (division_ele_xi > 1)
			{
				temp1 = (input_knot_vec_xi[i + 1] - input_knot_vec_xi[i]) / (double)division_ele_xi;
				for (j = 1; j < division_ele_xi; j++)
				{
					calc_xi[k] = calc_xi[k - 1] + temp1;
					k++;
				}
			}
		}
	}
	calc_xi[k] = input_knot_vec_xi[knot_n_xi - 1];
	division_n_xi = k + 1;
	element_n_xi = l;

	k = 0;
	l = 0;
	for (i = 0; i < knot_n_eta - 1; i++)
	{
		if (input_knot_vec_eta[i] != input_knot_vec_eta[i + 1])
		{
			calc_eta[k] = input_knot_vec_eta[i];
			dtilda_eta[l] = (input_knot_vec_eta[i + 1] - input_knot_vec_eta[i]) / 2.0;
			k++;
			l++;
			if (division_ele_eta > 1)
			{
				temp1 = (input_knot_vec_eta[i + 1] - input_knot_vec_eta[i]) / (double)division_ele_eta;
				for (j = 1; j < division_ele_eta; j++)
				{
					calc_eta[k] = calc_eta[k - 1] + temp1;
					k++;
				}
			}
		}
	}
	calc_eta[k] = input_knot_vec_eta[knot_n_eta - 1];
	division_n_eta = k + 1;
	element_n_eta = l;

	// メッシュ座標計算
	for (i = 0; i < division_n_xi; i++)
	{
		ii = i / division_ele_xi;
		kk = i % division_ele_xi;
		for (j = 0; j < division_n_eta; j++)
		{
			jj = j / division_ele_eta;
			ll = j % division_ele_eta;

			e = ii * (element_n_eta + 1) + jj;
			p = kk * (division_ele_eta + 1) + ll;
			ep = e * DIVISION_ELEMENT + p;

			lNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
						   cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
						   weight, order_xi, order_eta,
						   calc_xi[i], calc_eta[j],
						   &info.coord_x[i * division_n_eta + j], &info.coord_y[i * division_n_eta + j],
						   &dxi_x[ep], &deta_x[ep],
						   &dxi_y[ep], &deta_y[ep]);
		}
	}

	// 変位計算
	for (i = 0; i < division_n_xi; i++)
	{
		ii = i / division_ele_xi;
		kk = i % division_ele_xi;
		for (j = 0; j < division_n_eta; j++)
		{
			jj = j / division_ele_eta;
			ll = j % division_ele_eta;

			e = ii * (element_n_eta + 1) + jj;
			p = kk * (division_ele_eta + 1) + ll;
			ep = e * DIVISION_ELEMENT + p;

			lNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
						   disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
						   weight, order_xi, order_eta,
						   calc_xi[i], calc_eta[j],
						   &info.disp_x[i * division_n_eta + j], &info.disp_y[i * division_n_eta + j],
						   &dxi_disp_x[ep], &deta_disp_x[ep],
						   &dxi_disp_y[ep], &deta_disp_y[ep]);
		}
	}

	// 足りない微分値計算
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

				e = ii * (element_n_eta + 1) + jj;
				p = kk * (division_ele_eta + 1) + ll;
				ep = e * DIVISION_ELEMENT + p;

				rNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
							   cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
							   weight, order_xi, order_eta,
							   calc_xi[i], calc_eta[j],
							   &info.coord_x[i * division_n_eta + j], &info.coord_y[i * division_n_eta + j],
							   &dxi_x[ep], &deta_x[ep],
							   &dxi_y[ep], &deta_y[ep]);
				rNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
							   disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
							   weight, order_xi, order_eta,
							   calc_xi[i], calc_eta[j],
							   &info.disp_x[i * division_n_eta + j], &info.disp_y[i * division_n_eta + j],
							   &dxi_disp_x[ep], &deta_disp_x[ep],
							   &dxi_disp_y[ep], &deta_disp_y[ep]);
			}

			ll = division_ele_eta;
			i = ii * division_ele_xi;
			j = (jj + 1) * division_ele_eta;
			for (kk = 1; kk <= division_ele_xi; kk++)
			{
				i++;
				
				e = ii * (element_n_eta + 1) + jj;
				p = kk * (division_ele_eta + 1) + ll;
				ep = e * DIVISION_ELEMENT + p;

				rNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
							   cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
							   weight, order_xi, order_eta,
							   calc_xi[i], calc_eta[j],
							   &info.coord_x[i * division_n_eta + j], &info.coord_y[i * division_n_eta + j],
							   &dxi_x[ep], &deta_x[ep],
							   &dxi_y[ep], &deta_y[ep]);
				rNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
							   disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
							   weight, order_xi, order_eta,
							   calc_xi[i], calc_eta[j],
							   &info.disp_x[i * division_n_eta + j], &info.disp_y[i * division_n_eta + j],
							   &dxi_disp_x[ep], &deta_disp_x[ep],
							   &dxi_disp_y[ep], &deta_disp_y[ep]);
			}

			kk = division_ele_xi;
			ll = 0;
			i = (ii + 1) * division_ele_xi;
			j = jj * division_ele_eta;

			e = ii * (element_n_eta + 1) + jj;
			p = kk * (division_ele_eta + 1) + ll;
			ep = e * DIVISION_ELEMENT + p;

			rlNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
							cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
							weight, order_xi, order_eta,
							calc_xi[i], calc_eta[j],
							&info.coord_x[i * division_n_eta + j], &info.coord_y[i * division_n_eta + j],
							&dxi_x[ep], &deta_x[ep],
							&dxi_y[ep], &deta_y[ep]);
			rlNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
							disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
							weight, order_xi, order_eta,
							calc_xi[i], calc_eta[j],
							&info.disp_x[i * division_n_eta + j], &info.disp_y[i * division_n_eta + j],
							&dxi_disp_x[ep], &deta_disp_x[ep],
							&dxi_disp_y[ep], &deta_disp_y[ep]);

			kk = 0;
			ll = division_ele_eta;
			i = ii * division_ele_xi;
			j = (jj + 1) * division_ele_eta;

			e = ii * (element_n_eta + 1) + jj;
			p = kk * (division_ele_eta + 1) + ll;
			ep = e * DIVISION_ELEMENT + p;

			lrNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
							cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
							weight, order_xi, order_eta,
							calc_xi[i], calc_eta[j],
							&info.coord_x[i * division_n_eta + j], &info.coord_y[i * division_n_eta + j],
							&dxi_x[ep], &deta_x[ep],
							&dxi_y[ep], &deta_y[ep]);
			lrNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
							disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
							weight, order_xi, order_eta,
							calc_xi[i], calc_eta[j],
							&info.disp_x[i * division_n_eta + j], &info.disp_y[i * division_n_eta + j],
							&dxi_disp_x[ep], &deta_disp_x[ep],
							&dxi_disp_y[ep], &deta_disp_y[ep]);
		}
	}

	// calc strain
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
					e = i * (element_n_eta + 1) + j;
					p = k * (division_ele_eta + 1) + l;
					ep = e * DIVISION_ELEMENT + p;

					temp_matrix[0][0] = dxi_x[ep] * temp1;
					temp_matrix[0][1] = dxi_y[ep] * temp1;
					temp_matrix[1][0] = deta_x[ep] * temp2;
					temp_matrix[1][1] = deta_y[ep] * temp2;

					InverseMatrix_2x2(temp_matrix);

					info.strain[ep * 3 + 0] = temp_matrix[0][0] * temp1 * dxi_disp_x[ep] + temp_matrix[0][1] * temp2 * deta_disp_x[ep];
					info.strain[ep * 3 + 1] = temp_matrix[1][0] * temp1 * dxi_disp_y[ep] + temp_matrix[1][1] * temp2 * deta_disp_y[ep];
					info.strain[ep * 3 + 2] = temp_matrix[1][0] * temp1 * dxi_disp_x[ep] + temp_matrix[1][1] * temp2 * deta_disp_x[ep] + temp_matrix[0][0] * temp1 * dxi_disp_y[ep] + temp_matrix[0][1] * temp2 * deta_disp_y[ep];
				}
			}
		}
	}

	// calc stress
	for (i = 0; i < element_n_xi; i++)
	{
		for (j = 0; j < element_n_eta; j++)
		{
			for (k = 0; k < division_ele_xi + 1; k++)
			{
				for (l = 0; l < division_ele_eta + 1; l++)
				{
					e = i * (element_n_eta + 1) + j;
					p = k * (division_ele_eta + 1) + l;
					ep = e * DIVISION_ELEMENT + p;

					info.stress[ep * 3 + 0] = info.D[0 * D_MATRIX_SIZE + 0] * info.strain[ep * 3 + 0] + info.D[0 * D_MATRIX_SIZE + 1] * info.strain[ep * 3 + 1];
					info.stress[ep * 3 + 1] = info.D[1 * D_MATRIX_SIZE + 0] * info.strain[ep * 3 + 0] + info.D[1 * D_MATRIX_SIZE + 1] * info.strain[ep * 3 + 1];
					info.stress[ep * 3 + 2] = info.D[2 * D_MATRIX_SIZE + 2] * info.strain[ep * 3 + 2];
				}
			}
		}
	}

	// 書き込み
	fp = fopen("view.dat", "a");
	fprintf(fp, "%d\t%d\t%d\t%d\n", division_n_xi, division_n_eta, element_n_xi, element_n_eta);
	for (i = 0; i < division_n_xi; i++)
	{
		for (j = 0; j < division_n_eta; j++)
		{
			temp1 = info.disp_x[i * division_n_eta + j];
			temp2 = info.disp_y[i * division_n_eta + j];
			temp3 = sqrt(temp1 * temp1 + temp2 * temp2);
			fprintf(fp, "% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\n", info.coord_x[i * division_n_eta + j], info.coord_y[i * division_n_eta + j], info.disp_x[i * division_n_eta + j], info.disp_y[i * division_n_eta + j], temp3);
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
					e = i * (element_n_eta + 1) + j;
					p = k * (division_ele_eta + 1) + l;
					ep = e * DIVISION_ELEMENT + p;

					fprintf(fp, "% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\n",
							info.strain[ep * 3 + 0], info.strain[ep * 3 + 1], info.strain[ep * 3 + 2],
							info.stress[ep * 3 + 0], info.stress[ep * 3 + 1], info.stress[ep * 3 + 2]);
				}
			}
		}
	}
	fclose(fp);

	// グラフ用ファイル書き込み
	fp = fopen("disp_graph.txt", "a");
	for (i = 0; i < division_n_xi; i++)
	{
		for (j = 0; j < division_n_eta; j++)
		{
			temp1 = info.disp_x[i * division_n_eta + j];
			temp2 = info.disp_y[i * division_n_eta + j];
			fprintf(fp, "%d\t% 1.15e\t% 1.15e\t% 1.15e\t% 1.15e\n",
					patch_num,
					info.coord_x[i * division_n_eta + j], info.coord_y[i * division_n_eta + j],
					info.disp_x[i * division_n_eta + j], info.disp_y[i * division_n_eta + j]);
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
					e = i * (element_n_eta + 1) + j;
					p = k * (division_ele_eta + 1) + l;
					ep = e * DIVISION_ELEMENT + p;

					fprintf(fp, "%d\t% 1.15e\t% 1.15e\t% 1.15e\n",
							patch_num,
							info.coord_x[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
							info.coord_y[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
							info.stress[ep * 3 + 1]);
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
					if (info.coord_y[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)] == 0.000000000000000e+00)
					{
						e = i * (element_n_eta + 1) + j;
						p = k * (division_ele_eta + 1) + l;
						ep = e * DIVISION_ELEMENT + p;

						fprintf(fp, "%d\t% 1.15e\t% 1.15e\t% 1.15e\n",
								patch_num,
								info.coord_x[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
								info.coord_y[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
								info.stress[ep * 3 + 1]);
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
					e = i * (element_n_eta + 1) + j;
					p = k * (division_ele_eta + 1) + l;
					ep = e * DIVISION_ELEMENT + p;

					double stress_vm;
					double sum = info.stress[ep * 3 + 0] + info.stress[ep * 3 + 1];
					double dif = info.stress[ep * 3 + 0] - info.stress[ep * 3 + 1];
					double tau2 = info.stress[ep * 3 + 2] * info.stress[ep * 3 + 2];
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
							info.coord_x[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
							info.coord_y[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
							stress_vm, temp3);
				}
			}
		}
	}
	fclose(fp);

	free(calc_xi), free(calc_eta), free(dtilda_xi), free(dtilda_eta);
	free(dxi_x), free(dxi_y), free(deta_x), free(deta_y), free(dxi_disp_x), free(dxi_disp_y), free(deta_disp_x), free(deta_disp_y);
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
						 double *weight_glo, int patch_num, information info)
{
	int i, j, k, l;
	int ii, jj, kk, ll;
	int e, p, ep;
	double temp1, temp2, temp3;

	double output_xi, output_eta;
	double disp_x_glo, disp_y_glo;
	double strain_xx_glo = 0;
	double strain_yy_glo = 0;
	double strain_xy_glo = 0;

	// 計算するξ,ηの値決定と ∂ξ/∂チルダξ, ∂η/∂チルダη の計算
	double *calc_xi_loc = (double *)malloc(sizeof(double) * knot_n_xi_loc * division_ele_xi);		// calc_xi_loc[MAX_POINTS], 計算するξの値local
	double *calc_eta_loc = (double *)malloc(sizeof(double) * knot_n_eta_loc * division_ele_eta);	// calc_eta_loc[MAX_POINTS], 計算するηの値local

	double *dxi_x = (double *)malloc(sizeof(double) * info.Total_Element_to_mesh[Total_mesh] * DIVISION_ELEMENT);		// dxi_x[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]
	double *dxi_y = (double *)malloc(sizeof(double) * info.Total_Element_to_mesh[Total_mesh] * DIVISION_ELEMENT);		// dxi_y[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]
	double *deta_x = (double *)malloc(sizeof(double) * info.Total_Element_to_mesh[Total_mesh] * DIVISION_ELEMENT);		// deta_x[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]
	double *deta_y = (double *)malloc(sizeof(double) * info.Total_Element_to_mesh[Total_mesh] * DIVISION_ELEMENT);		// deta_y[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]

	k = 0;
	l = 0;
	for (i = 0; i < knot_n_xi_loc - 1; i++)
	{
		if (knot_vec_xi_loc[i] != knot_vec_xi_loc[i + 1])
		{
			calc_xi_loc[k] = knot_vec_xi_loc[i];
			k++;
			l++;
			if (division_ele_xi > 1)
			{
				temp1 = (knot_vec_xi_loc[i + 1] - knot_vec_xi_loc[i]) / (double)division_ele_xi;
				for (j = 1; j < division_ele_xi; j++)
				{
					calc_xi_loc[k] = calc_xi_loc[k - 1] + temp1;
					k++;
				}
			}
		}
	}
	calc_xi_loc[k] = knot_vec_xi_loc[knot_n_xi_loc - 1];
	division_n_xi = k + 1;
	element_n_xi = l;

	k = 0;
	l = 0;
	for (i = 0; i < knot_n_eta_loc - 1; i++)
	{
		if (knot_vec_eta_loc[i] != knot_vec_eta_loc[i + 1])
		{
			calc_eta_loc[k] = knot_vec_eta_loc[i];

			k++;
			l++;
			if (division_ele_eta > 1)
			{
				temp1 = (knot_vec_eta_loc[i + 1] - knot_vec_eta_loc[i]) / (double)division_ele_eta;
				for (j = 1; j < division_ele_eta; j++)
				{
					calc_eta_loc[k] = calc_eta_loc[k - 1] + temp1;
					k++;
				}
			}
		}
	}
	calc_eta_loc[k] = knot_vec_eta_loc[knot_n_eta_loc - 1];
	division_n_eta = k + 1;
	element_n_eta = l;

	// メッシュ座標計算
	for (i = 0; i < division_n_xi; i++)
	{
		ii = i / division_ele_xi;
		kk = i % division_ele_xi;
		for (j = 0; j < division_n_eta; j++)
		{
			jj = j / division_ele_eta;
			ll = j % division_ele_eta;

			e = ii * (element_n_eta + 1) + jj;
			p = kk * (division_ele_eta + 1) + ll;
			ep = e * DIVISION_ELEMENT + p;

			lNURBS_surface(knot_vec_xi_loc, knot_vec_eta_loc,
						   cntl_px_loc, cntl_py_loc, cntl_p_n_xi_loc, cntl_p_n_eta_loc,
						   weight_loc, order_xi_loc, order_eta_loc,
						   calc_xi_loc[i], calc_eta_loc[j],
						   &info.coord_x[i * division_n_eta + j], &info.coord_y[i * division_n_eta + j],
						   &dxi_x[ep], &deta_x[ep],
						   &dxi_y[ep], &deta_y[ep]);

			CalcXiEtaByNR(info.coord_x[i * division_n_eta + j], info.coord_y[i * division_n_eta + j],
						  knot_vec_xi_glo, knot_vec_eta_glo,
						  cntl_px_glo, cntl_py_glo,
						  disp_cntl_px_glo, disp_cntl_py_glo,
						  cntl_p_n_xi_glo, cntl_p_n_eta_glo,
						  weight_glo, order_xi_glo, order_eta_glo,
						  &output_xi, &output_eta,
						  &disp_x_glo, &disp_y_glo,
						  &strain_xx_glo, &strain_yy_glo, &strain_xy_glo);

			// ローカル内の表示点上のグローバル変位
			info.disp_x[i * division_n_eta + j] += disp_x_glo;
			info.disp_y[i * division_n_eta + j] += disp_y_glo;

			// ローカル内の表示点上のグローバルひずみ
			info.strain[ep * 3 + 0] += strain_xx_glo;
			info.strain[ep * 3 + 1] += strain_yy_glo;
			info.strain[ep * 3 + 2] += strain_xy_glo;
			if (jj > 0 && ll == 0)
			{
				e = ii * (element_n_eta + 1) + (jj - 1);
				p = kk * (division_ele_eta + 1) + division_ele_eta;
				ep = e * DIVISION_ELEMENT + p;

				info.strain[ep * 3 + 0] += strain_xx_glo;
				info.strain[ep * 3 + 1] += strain_yy_glo;
				info.strain[ep * 3 + 2] += strain_xy_glo;
			}
			if (ii > 0 && kk == 0)
			{
				e = (ii - 1) * (element_n_eta + 1) + jj;
				p = division_ele_xi * (division_ele_eta + 1) + ll;
				ep = e * DIVISION_ELEMENT + p;

				info.strain[ep * 3 + 0] += strain_xx_glo;
				info.strain[ep * 3 + 1] += strain_yy_glo;
				info.strain[ep * 3 + 2] += strain_xy_glo;
			}
			if (ii > 0 && jj > 0 && kk == 0 && ll == 0)
			{
				e = (ii - 1) * (element_n_eta + 1) + (jj - 1);
				p = division_ele_xi * (division_ele_eta + 1) + division_ele_eta;
				ep = e * DIVISION_ELEMENT + p;

				info.strain[ep * 3 + 0] += strain_xx_glo;
				info.strain[ep * 3 + 1] += strain_yy_glo;
				info.strain[ep * 3 + 2] += strain_xy_glo;
			}
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
					e = i * (element_n_eta + 1) + j;
					p = k * (division_ele_eta + 1) + l;
					ep = e * DIVISION_ELEMENT + p;

					info.stress[ep * 3 + 0] = info.D[0 * D_MATRIX_SIZE + 0] * info.strain[ep * 3 + 0] + info.D[0 * D_MATRIX_SIZE + 1] * info.strain[ep * 3 + 1];
					info.stress[ep * 3 + 1] = info.D[1 * D_MATRIX_SIZE + 0] * info.strain[ep * 3 + 0] + info.D[1 * D_MATRIX_SIZE + 1] * info.strain[ep * 3 + 1];
					info.stress[ep * 3 + 2] = info.D[2 * D_MATRIX_SIZE + 2] * info.strain[ep * 3 + 2];
				}
			}
		}
	}

	//書き込み
	fp = fopen("overlay_view.dat", "a");
	fprintf(fp, "%d\t%d\t%d\t%d\n", division_n_xi, division_n_eta, element_n_xi, element_n_eta);
	for (i = 0; i < division_n_xi; i++)
	{
		for (j = 0; j < division_n_eta; j++)
		{
			temp1 = info.disp_x[i * division_n_eta + j];
			temp2 = info.disp_y[i * division_n_eta + j];
			temp3 = sqrt(temp1 * temp1 + temp2 * temp2);
			fprintf(fp, "% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\n", info.coord_x[i * division_n_eta + j], info.coord_y[i * division_n_eta + j], info.disp_x[i * division_n_eta + j], info.disp_y[i * division_n_eta + j], temp3);
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
					e = i * (element_n_eta + 1) + j;
					p = k * (division_ele_eta + 1) + l;
					ep = e * DIVISION_ELEMENT + p;

					fprintf(fp, "% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\n",
							info.strain[ep * 3 + 0], info.strain[ep * 3 + 1], info.strain[ep * 3 + 2],
							info.stress[ep * 3 + 0], info.stress[ep * 3 + 1], info.stress[ep * 3 + 2]);
				}
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
			temp1 = info.disp_x[i * division_n_eta + j];
			temp2 = info.disp_y[i * division_n_eta + j];
			fprintf(fp, "%d\t% 1.15e\t% 1.15e\t% 1.15e\t% 1.15e\n",
					patch_num,
					info.coord_x[i * division_n_eta + j], info.coord_y[i * division_n_eta + j],
					info.disp_x[i * division_n_eta + j], info.disp_y[i * division_n_eta + j]);
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
					e = i * (element_n_eta + 1) + j;
					p = k * (division_ele_eta + 1) + l;
					ep = e * DIVISION_ELEMENT + p;

					fprintf(fp, "% 1.10e\t% 1.10e\t% 1.10e\n",
							info.coord_x[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
							info.coord_y[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
							info.stress[ep * 3 + 0]);
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
					e = i * (element_n_eta + 1) + j;
					p = k * (division_ele_eta + 1) + l;
					ep = e * DIVISION_ELEMENT + p;

					fprintf(fp, "% 1.15e\t% 1.15e\t% 1.15e\n",
							info.coord_x[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
							info.coord_y[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
							info.stress[ep * 3 + 1]);
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
					if (info.coord_y[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)] == 0.000000000000000e+00)
					{
						e = i * (element_n_eta + 1) + j;
						p = k * (division_ele_eta + 1) + l;
						ep = e * DIVISION_ELEMENT + p;

						fprintf(fp, "% 1.15e\t% 1.15e\t% 1.15e\n",
								info.coord_x[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
								info.coord_y[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
								info.stress[ep * 3 + 1]);
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
					e = i * (element_n_eta + 1) + j;
					p = k * (division_ele_eta + 1) + l;
					ep = e * DIVISION_ELEMENT + p;

					double stress_rr, stress_theta;
					double sum = info.stress[ep * 3 + 0] + info.stress[ep * 3 + 1];
					double dif = info.stress[ep * 3 + 0] - info.stress[ep * 3 + 1];
					double tau2 = info.stress[ep * 3 + 2] * info.stress[ep * 3 + 2];
					stress_rr = sum * 0.5 + sqrt(dif * dif + 4 * tau2) * 0.5;
					stress_theta = sum * 0.5 - sqrt(dif * dif + 4 * tau2) * 0.5;
					fprintf(fp, "% 1.10e\t% 1.10e\t% 1.10e\t% 1.10e\t% 1.10e\t% 1.10e\n",
							calc_xi_loc[i * division_ele_xi + k],
							calc_eta_loc[j * division_ele_eta + l],
							info.coord_x[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
							info.coord_y[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
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
					e = i * (element_n_eta + 1) + j;
					p = k * (division_ele_eta + 1) + l;
					ep = e * DIVISION_ELEMENT + p;

					double stress_vm;
					double sum = info.stress[ep * 3 + 0] + info.stress[ep * 3 + 1];
					double dif = info.stress[ep * 3 + 0] - info.stress[ep * 3 + 1];
					double tau2 = info.stress[ep * 3 + 2] * info.stress[ep * 3 + 2];
					double temp1, temp2;
					temp1 = 0.5 * sum;
					temp2 = 0.5 * sqrt(dif * dif + 4 * tau2);
					stress_vm = sqrt(temp1 * temp1 + 3 * temp2 * temp2);
					fprintf(fp, "% 1.10e\t% 1.10e\t% 1.10e\t% 1.10e\t% 1.10e\n",
							calc_xi_loc[i * division_ele_xi + k],
							calc_eta_loc[j * division_ele_eta + l],
							info.coord_x[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
							info.coord_y[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
							stress_vm);
				}
			}
		}
	}
	fclose(fp);

	free(calc_xi_loc), free(calc_eta_loc);
	free(dxi_x), free(dxi_y), free(deta_x), free(deta_y);
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

		// 収束した場合
		if (temp_tol_x + temp_tol_y < tol)
		{
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

			for (i = 0; i < knot_n_xi; i++)
			{
				if (input_knot_vec_xi[i] < temp_xi && temp_xi <= input_knot_vec_xi[i + 1])
				{
					dtilda_xi = (input_knot_vec_xi[i + 1] - input_knot_vec_xi[i]) / 2.0;
					break;
				}
			}
			for (i = 0; i < knot_n_eta; i++)
			{
				if (input_knot_vec_eta[i] < temp_eta && temp_eta <= input_knot_vec_eta[i + 1])
				{
					dtilda_eta = (input_knot_vec_eta[i + 1] - input_knot_vec_eta[i]) / 2.0;
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

			rNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
						   disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
						   weight, order_xi, order_eta,
						   temp_xi, temp_eta,
						   &disp_x, &disp_y,
						   &dxi_disp_x, &deta_disp_x,
						   &dxi_disp_y, &deta_disp_y);

			temp_matrix2[0][0] = dxi_x * dtilda_xi;
			temp_matrix2[0][1] = dxi_y * dtilda_xi;
			temp_matrix2[1][0] = deta_x * dtilda_eta;
			temp_matrix2[1][1] = deta_y * dtilda_eta;

			InverseMatrix_2x2(temp_matrix2);

			strain_xx = temp_matrix2[0][0] * dtilda_xi * dxi_disp_x + temp_matrix2[0][1] * dtilda_eta * deta_disp_x;
			strain_yy = temp_matrix2[1][0] * dtilda_xi * dxi_disp_y + temp_matrix2[1][1] * dtilda_eta * deta_disp_y;
			strain_xy = temp_matrix2[1][0] * dtilda_xi * dxi_disp_x + temp_matrix2[1][1] * dtilda_eta * deta_disp_x + temp_matrix2[0][0] * dtilda_xi * dxi_disp_y + temp_matrix2[0][1] * dtilda_eta * deta_disp_y;

			temp = sqrt(disp_x * disp_x + disp_y * disp_y);

			(*disp_x_glo) = disp_x;
			(*disp_y_glo) = disp_y;

			(*strain_xx_glo) = strain_xx;
			(*strain_yy_glo) = strain_yy;
			(*strain_xy_glo) = strain_xy;

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

		// 収束した場合
		if (temp_tol_x + temp_tol_y < tol)
		{
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

			for (i = 0; i < knot_n_xi; i++)
			{
				if (input_knot_vec_xi[i] <= temp_xi && temp_xi < input_knot_vec_xi[i + 1])
				{
					dtilda_xi = (input_knot_vec_xi[i + 1] - input_knot_vec_xi[i]) / 2.0;
					break;
				}
			}
			for (i = 0; i < knot_n_eta; i++)
			{
				if (input_knot_vec_eta[i] <= temp_eta && temp_eta < input_knot_vec_eta[i + 1])
				{
					dtilda_eta = (input_knot_vec_eta[i + 1] - input_knot_vec_eta[i]) / 2.0;
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

			lNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
						   disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
						   weight, order_xi, order_eta,
						   temp_xi, temp_eta,
						   &disp_x, &disp_y,
						   &dxi_disp_x, &deta_disp_x,
						   &dxi_disp_y, &deta_disp_y);

			temp_matrix2[0][0] = dxi_x * dtilda_xi;
			temp_matrix2[0][1] = dxi_y * dtilda_xi;
			temp_matrix2[1][0] = deta_x * dtilda_eta;
			temp_matrix2[1][1] = deta_y * dtilda_eta;

			InverseMatrix_2x2(temp_matrix2);

			strain_xx = temp_matrix2[0][0] * dtilda_xi * dxi_disp_x + temp_matrix2[0][1] * dtilda_eta * deta_disp_x;
			strain_yy = temp_matrix2[1][0] * dtilda_xi * dxi_disp_y + temp_matrix2[1][1] * dtilda_eta * deta_disp_y;
			strain_xy = temp_matrix2[1][0] * dtilda_xi * dxi_disp_x + temp_matrix2[1][1] * dtilda_eta * deta_disp_x + temp_matrix2[0][0] * dtilda_xi * dxi_disp_y + temp_matrix2[0][1] * dtilda_eta * deta_disp_y;

			temp = sqrt(disp_x * disp_x + disp_y * disp_y);

			(*disp_x_glo) = disp_x;
			(*disp_y_glo) = disp_y;

			(*strain_xx_glo) = strain_xx;
			(*strain_yy_glo) = strain_yy;
			(*strain_xy_glo) = strain_xy;

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

		// 収束した場合
		if (temp_tol_x + temp_tol_y < tol)
		{
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

			for (i = 0; i < knot_n_xi; i++)
			{
				if (input_knot_vec_xi[i] < temp_xi && temp_xi <= input_knot_vec_xi[i + 1])
				{
					dtilda_xi = (input_knot_vec_xi[i + 1] - input_knot_vec_xi[i]) / 2.0;
					break;
				}
			}
			for (i = 0; i < knot_n_eta; i++)
			{
				if (input_knot_vec_eta[i] <= temp_eta && temp_eta < input_knot_vec_eta[i + 1])
				{
					dtilda_eta = (input_knot_vec_eta[i + 1] - input_knot_vec_eta[i]) / 2.0;
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

			rlNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
							disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
							weight, order_xi, order_eta,
							temp_xi, temp_eta,
							&disp_x, &disp_y,
							&dxi_disp_x, &deta_disp_x,
							&dxi_disp_y, &deta_disp_y);

			temp_matrix2[0][0] = dxi_x * dtilda_xi;
			temp_matrix2[0][1] = dxi_y * dtilda_xi;
			temp_matrix2[1][0] = deta_x * dtilda_eta;
			temp_matrix2[1][1] = deta_y * dtilda_eta;

			InverseMatrix_2x2(temp_matrix2);

			strain_xx = temp_matrix2[0][0] * dtilda_xi * dxi_disp_x + temp_matrix2[0][1] * dtilda_eta * deta_disp_x;
			strain_yy = temp_matrix2[1][0] * dtilda_xi * dxi_disp_y + temp_matrix2[1][1] * dtilda_eta * deta_disp_y;
			strain_xy = temp_matrix2[1][0] * dtilda_xi * dxi_disp_x + temp_matrix2[1][1] * dtilda_eta * deta_disp_x + temp_matrix2[0][0] * dtilda_xi * dxi_disp_y + temp_matrix2[0][1] * dtilda_eta * deta_disp_y;

			temp = sqrt(disp_x * disp_x + disp_y * disp_y);

			(*disp_x_glo) = disp_x;
			(*disp_y_glo) = disp_y;

			(*strain_xx_glo) = strain_xx;
			(*strain_yy_glo) = strain_yy;
			(*strain_xy_glo) = strain_xy;

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

		// 収束した場合
		if (temp_tol_x + temp_tol_y < tol)
		{
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

			for (i = 0; i < knot_n_xi; i++)
			{
				if (input_knot_vec_xi[i] <= temp_xi && temp_xi < input_knot_vec_xi[i + 1])
				{
					dtilda_xi = (input_knot_vec_xi[i + 1] - input_knot_vec_xi[i]) / 2.0;
					break;
				}
			}
			for (i = 0; i < knot_n_eta; i++)
			{
				if (input_knot_vec_eta[i] < temp_eta && temp_eta <= input_knot_vec_eta[i + 1])
				{
					dtilda_eta = (input_knot_vec_eta[i + 1] - input_knot_vec_eta[i]) / 2.0;
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

			lrNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
							disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
							weight, order_xi, order_eta,
							temp_xi, temp_eta,
							&disp_x, &disp_y,
							&dxi_disp_x, &deta_disp_x,
							&dxi_disp_y, &deta_disp_y);

			temp_matrix2[0][0] = dxi_x * dtilda_xi;
			temp_matrix2[0][1] = dxi_y * dtilda_xi;
			temp_matrix2[1][0] = deta_x * dtilda_eta;
			temp_matrix2[1][1] = deta_y * dtilda_eta;

			InverseMatrix_2x2(temp_matrix2);

			strain_xx = temp_matrix2[0][0] * dtilda_xi * dxi_disp_x + temp_matrix2[0][1] * dtilda_eta * deta_disp_x;
			strain_yy = temp_matrix2[1][0] * dtilda_xi * dxi_disp_y + temp_matrix2[1][1] * dtilda_eta * deta_disp_y;
			strain_xy = temp_matrix2[1][0] * dtilda_xi * dxi_disp_x + temp_matrix2[1][1] * dtilda_eta * deta_disp_x + temp_matrix2[0][0] * dtilda_xi * dxi_disp_y + temp_matrix2[0][1] * dtilda_eta * deta_disp_y;

			temp = sqrt(disp_x * disp_x + disp_y * disp_y);

			(*disp_x_glo) = disp_x;
			(*disp_y_glo) = disp_y;

			(*strain_xx_glo) = strain_xx;
			(*strain_yy_glo) = strain_yy;
			(*strain_xy_glo) = strain_xy;

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
	return 0;
}


// output
void output_for_viewer(information info)
{
	// 全てのローカルパッチについて1つのinput.txt作成, 重ね合わせ結果出力のため
	int i, j, k;

	int l_patch, g_patch;
	int l_cntl_p, g_cntl_p;
	int l_cnst, g_cnst;
	int l_load, g_load;
	int l_dist, g_dist;
	int l_temp;

	g_patch = info.Total_Patch_to_mesh[1];
	g_cntl_p = info.Total_Control_Point_to_mesh[1];
	g_cnst = info.Total_Constraint_to_mesh[1];
	g_load = info.Total_Load_to_mesh[1];
	g_dist = info.Total_DistributeForce_to_mesh[1];

	fp = fopen("input_local.txt", "w");
	fprintf(fp, "%.3f\t%.6f\n\n", E, nu);

	l_patch = info.Total_Patch_to_mesh[Total_mesh] - g_patch;
	fprintf(fp, "%d\n\n", l_patch);

	l_cntl_p = info.Total_Control_Point_to_mesh[Total_mesh] - g_cntl_p;
	fprintf(fp, "%d\n\n", l_cntl_p);

	for (i = 0; i < l_patch; i++)
	{
		fprintf(fp, "%d\t%d\n", info.Order[(i + g_patch) * DIMENSION + 0], info.Order[(i + g_patch) * DIMENSION + 1]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < l_patch; i++)
	{
		fprintf(fp, "%d\t%d\n", info.No_knot[(i + g_patch) * DIMENSION + 0], info.No_knot[(i + g_patch) * DIMENSION + 1]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < l_patch; i++)
	{
		fprintf(fp, "%d\t%d\n", info.No_Control_point[(i + g_patch) * DIMENSION + 0], info.No_Control_point[(i + g_patch) * DIMENSION + 1]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < l_patch; i++)
	{
		for (j = 0; j < info.No_Control_point_in_patch[i + g_patch]; j++)
		{
			l_temp = info.Patch_Control_point[info.Total_Control_Point_to_patch[i + g_patch] + j] - g_cntl_p;
			fprintf(fp, "%d\t", l_temp);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp,"\n");

	l_cnst = info.Total_Constraint_to_mesh[Total_mesh] - g_cnst;
	l_load = info.Total_Load_to_mesh[Total_mesh] - g_load;
	l_dist = info.Total_DistributeForce_to_mesh[Total_mesh] - g_dist;
	fprintf(fp,"%d\t%d\t%d\n\n", l_cnst, l_load, l_dist);

	for (i = 0; i < l_patch; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			for (k = 0; k < info.No_knot[(i + g_patch) * DIMENSION + j]; k++)
			{
				fprintf(fp, "%.5f\t", info.Position_Knots[info.Total_Knot_to_patch_dim[(i + g_patch) * DIMENSION + j] + k]);
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
			fprintf(fp, "%.10f\t", info.Node_Coordinate[(i + g_cntl_p) * (DIMENSION + 1) + j]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp,"\n");

	for (i = 0; i < l_cnst; i++)
	{
		l_temp = info.Constraint_Node_Dir[(i + g_cnst) * 2 + 0] - g_cntl_p;
		fprintf(fp, "%d\t%d\t%.10f\n", l_temp, info.Constraint_Node_Dir[(i + g_cnst) * 2 + 1], info.Value_of_Constraint[i + g_cnst]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < l_load; i++)
	{
		l_temp = info.Load_Node_Dir[(i + g_load) * 2 + 0] - g_cntl_p;
		fprintf(fp, "%d\t%d\t%.10f\n", l_temp, info.Load_Node_Dir[(i + g_load) * 2 + 1], info.Value_of_Load[i + g_load]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < l_dist; i++)
	{
		fprintf(fp, "%d\t%d\t%d\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n",
				info.type_load_array[i + g_dist],
				info.iPatch_array[i + g_dist],
				info.iCoord_array[i + g_dist],
				info.val_Coord_array[i + g_dist],
				info.Range_Coord_array[(i + g_dist) * 2 + 0],
				info.Range_Coord_array[(i + g_dist) * 2 + 1],
				info.Coeff_Dist_Load_array[(i + g_dist) * 3 + 0],
				info.Coeff_Dist_Load_array[(i + g_dist) * 3 + 1],
				info.Coeff_Dist_Load_array[(i + g_dist) * 3 + 2]);
	}
	fclose(fp);

	fp = fopen("input_for_NURBS.txt", "w");
	fprintf(fp, "%.3f\t%.6f\n\n", E, nu);

	fprintf(fp, "%d\n\n", info.Total_Patch_to_mesh[Total_mesh]);

	fprintf(fp, "%d\n\n", info.Total_Control_Point_to_mesh[Total_mesh]);

	for (i = 0; i < info.Total_Patch_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp, "%d\t%d\n", info.Order[i * DIMENSION + 0], info.Order[i * DIMENSION + 1]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < info.Total_Patch_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp, "%d\t%d\n", info.No_knot[i * DIMENSION + 0], info.No_knot[i * DIMENSION + 1]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < info.Total_Patch_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp, "%d\t%d\n", info.No_Control_point[i * DIMENSION + 0], info.No_Control_point[i * DIMENSION + 1]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < info.Total_Patch_to_mesh[Total_mesh]; i++)
	{
		for (j = 0; j < info.No_Control_point_in_patch[i]; j++)
		{
			fprintf(fp, "%d\t", info.Patch_Control_point[info.Total_Control_Point_to_patch[i] + j]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp,"\n");

	fprintf(fp,"%d\t%d\t%d\n\n",
				info.Total_Constraint_to_mesh[Total_mesh],
				info.Total_Load_to_mesh[Total_mesh],
				info.Total_DistributeForce_to_mesh[Total_mesh]);

	for (i = 0; i < info.Total_Patch_to_mesh[Total_mesh]; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			for (k = 0; k < info.No_knot[i * DIMENSION + j]; k++)
			{
				fprintf(fp, "%.5f\t", info.Position_Knots[info.Total_Knot_to_patch_dim[i * DIMENSION + j] + k]);
			}
			fprintf(fp, "\n");
		}
	}
	fprintf(fp,"\n");

	for (i = 0; i < info.Total_Control_Point_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp, "%d\t", i);
		for (j = 0; j < DIMENSION + 1; j++)
		{
			fprintf(fp, "%.10f\t", info.Node_Coordinate[i * (DIMENSION + 1) + j]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp,"\n");

	for (i = 0; i < info.Total_Constraint_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp, "%d\t%d\t%.10f\n",
				info.Constraint_Node_Dir[i * 2 + 0],
				info.Constraint_Node_Dir[i * 2 + 1],
				info.Value_of_Constraint[i]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < info.Total_Load_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp, "%d\t%d\t%.10f\n",
				info.Load_Node_Dir[i * 2 + 0],
				info.Load_Node_Dir[i * 2 + 1],
				info.Value_of_Load[i]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < info.Total_DistributeForce_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp, "%d\t%d\t%d\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n",
				info.type_load_array[i],
				info.iPatch_array[i],
				info.iCoord_array[i],
				info.val_Coord_array[i],
				info.Range_Coord_array[i * 2 + 0],
				info.Range_Coord_array[i * 2 + 1],
				info.Coeff_Dist_Load_array[i * 3 + 0],
				info.Coeff_Dist_Load_array[i * 3 + 1],
				info.Coeff_Dist_Load_array[i * 3 + 2]);
	}

	fclose(fp);

	fp = fopen("Displacement.dat", "w");
	fprintf(fp, "label=Displacement\n");
	fprintf(fp, "num_items=%d\n", info.Total_Control_Point_to_mesh[Total_mesh]);
	fprintf(fp, "\n");
	for (j = 0; j < info.Total_Control_Point_to_mesh[Total_mesh]; j++)
	{
		fprintf(fp, "%d:	%.16e %.16e ", j, info.Displacement[j * DIMENSION + 0], info.Displacement[j * DIMENSION + 1]);
		fprintf(fp, "\n");
	}
	fclose(fp);

	fp = fopen("Displacement_loc.dat", "w");
	int loc_cntl_p_n = info.Total_Control_Point_to_mesh[Total_mesh] - info.Total_Control_Point_to_mesh[1];
	int glo_cntl_p_n = info.Total_Control_Point_to_mesh[1];
	fprintf(fp, "label=Displacement\n"
				"num_items=%d\n\n", loc_cntl_p_n);
	for (i = 0; i < loc_cntl_p_n; i++)
	{
		fprintf(fp, "%d:	%.16e %.16e \n", i, info.Displacement[(i + glo_cntl_p_n) * DIMENSION + 0], info.Displacement[(i + glo_cntl_p_n) * DIMENSION + 1]);
	}
	fclose(fp);
}


void K_output_svg(information info)
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
				temp_count = RowCol_to_icount(i, j, info);
			}
			else if (i > j)
			{
				temp_count = RowCol_to_icount(j, i, info);
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





/*

// 各種値 
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
		// printf("El_No:%d,info.No_Control_point_ON_ELEMENT[%d]:%d\n", El_No,info.Element_patch[El_No],info.No_Control_point_ON_ELEMENT[info.Element_patch[El_No]]);
		for (i = 0; i < info.No_Control_point_ON_ELEMENT[info.Element_patch[e]]; i++)
		{
			for (j = 0; j < DIMENSION; j++)
			{
				U[i * DIMENSION + j] = info.Displacement[info.Controlpoint_of_Element[e * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + j];
				X[i][j] = info.Node_Coordinate[info.Controlpoint_of_Element[e * MAX_NO_CCpoint_ON_ELEMENT + i] * (DIMENSION + 1) + j];
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


void Make_Stress_2D(double E, double nu, int Total_Element, int DM)
{
	static double info.D[D_MATRIX_SIZE][D_MATRIX_SIZE];
	int e, i, j, k;
	Make_gauss_array(0);
	Make_D_Matrix_2D(info.D, E, nu, DM);

	for (e = 0; e < Total_Element; e++)
	{
		for (k = 0; k < GP_2D; k++)
			for (i = 0; i < N_STRESS; i++)
				Stress[e][k][i] = 0.0;
		for (k = 0; k < GP_2D; k++)
			for (i = 0; i < D_MATRIX_SIZE; i++)
				for (j = 0; j < D_MATRIX_SIZE; j++)
					Stress[e][k][i] += info.D[i][j] * Strain[e][k][j];
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
	for (re = 0; re < info.real_Total_Element_to_mesh[Total_mesh]; re++)
	{
		e = info.real_element[re];
		// printf("%d,info.No_Control_point_ON_ELEMENT[%d]:%d\n",e,info.Element_patch[e], info.No_Control_point_ON_ELEMENT[info.Element_patch[e]]);
		// Bマトリックスを取得
		for (i = 0; i < info.No_Control_point_ON_ELEMENT[info.Element_patch[e]]; i++)
		{
			for (j = 0; j < DIMENSION; j++)
			{
				X[i][j] = info.Node_Coordinate[info.Controlpoint_of_Element[e * MAX_NO_CCpoint_ON_ELEMENT + i] * (DIMENSION + 1) + j];
			}
		}
		for (k = 0; k < GP_2D; k++)
		{
			Make_B_Matrix(e, B, Gxi[k], X, &J);
			for (j = 0; j < D_MATRIX_SIZE; j++)
				for (l = 0; l < DIMENSION; l++)
					for (i = 0; i < info.No_Control_point_ON_ELEMENT[info.Element_patch[e]]; i++)
						ReactionForce[info.Controlpoint_of_Element[e * MAX_NO_CCpoint_ON_ELEMENT + i] * DIMENSION + l] += B[j][i * DIMENSION + l] * Stress[e][k][j] * w[k] * J;
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

*/