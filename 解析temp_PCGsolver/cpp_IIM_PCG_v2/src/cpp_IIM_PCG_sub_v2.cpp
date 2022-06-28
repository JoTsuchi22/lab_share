#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// header
#include "S_IGA_header.h"
#include "S_IGA_sub.h"

using namespace std;

// Read file 1st time
void Get_Input_1(int tm, char **argv, information *info)
{
	char s[256];
	int temp_i;

	int i, j;

	fp = fopen(argv[tm + 1], "r");

	// DIMENSION
	fscanf(fp, "%d", &temp_i);
	info->DIMENSION = temp_i;
	fgets(s, 256, fp);
	printf("DIMENSION: %d\n", info->DIMENSION);
	if (info->DIMENSION != 2 && info->DIMENSION != 3)
	{
		printf("Error, wrong DIMENSION in input data\n");
		exit(1);
	}

	// 材料定数
	fscanf(fp, "%lf %lf", &E, &nu);
	fgets(s, 256, fp);
	printf("E: %le, nu: %le\n", E, nu);

	// パッチ数
	fscanf(fp, "%d", &temp_i);
	fgets(s, 256, fp);
	int No_Patch = temp_i;
	int *CP = (int *)malloc(sizeof(int) * temp_i * info->DIMENSION);
	printf("No_Patch: %d\n", temp_i);
	info->Total_Patch_on_mesh[tm] = temp_i;
	info->Total_Patch_to_mesh[tm + 1] = info->Total_Patch_to_mesh[tm] + temp_i;
	printf("Total_Patch_to_mesh[%d] = %d\n", tm, info->Total_Patch_to_mesh[tm]);

	// コントロールポイント数
	fscanf(fp, "%d", &temp_i);
	fgets(s, 256, fp);
	printf("Total_Control_Point: %d\n", temp_i);
	info->Total_Control_Point_on_mesh[tm] = temp_i;
	info->Total_Control_Point_to_mesh[tm + 1] = info->Total_Control_Point_to_mesh[tm] + temp_i;
	printf("Total_Control_Point_to_mesh[%d] = %d\n", tm, info->Total_Control_Point_to_mesh[tm]);

	// 各方向の次数(スキップ)
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			fscanf(fp, "%d", &temp_i);
		}
	}
	fgets(s, 256, fp);

	// ノット数
	info->Total_Knot_to_mesh[tm + 1] = info->Total_Knot_to_mesh[tm];
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			fscanf(fp, "%d", &temp_i);
			info->Total_Knot_to_mesh[tm + 1] += temp_i;
		}
	}

	// 各パッチ各方向のコントロールポイント数(スキップ)
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			fscanf(fp, "%d", &temp_i);
			CP[i * info->DIMENSION + j] = temp_i;
		}

		int temp_MAX_CP = 1;
		for (j = 0; j < info->DIMENSION; j++)
		{
			temp_MAX_CP *= CP[i * info->DIMENSION + j];
		}
		MAX_CP += temp_MAX_CP;
	}
	fgets(s, 256, fp);

	// パッチコネクティビティ(スキップ)
	for (i = 0; i < info->Total_Patch_on_mesh[tm]; i++)
	{
		int temp_MAX_CP = 1;
		for (j = 0; j < info->DIMENSION; j++)
		{
			temp_MAX_CP *= CP[i * info->DIMENSION + j];
		}
		for (j = 0; j < temp_MAX_CP; j++)
		{
			fscanf(fp, "%d", &temp_i);
		}
	}
	fgets(s, 256, fp);

	// 境界条件
	int Total_Constraint, Total_Load, Total_DistributeForce;
	fscanf(fp, "%d %d %d", &Total_Constraint, &Total_Load, &Total_DistributeForce);
	info->Total_Constraint_to_mesh[tm + 1] = info->Total_Constraint_to_mesh[tm] + Total_Constraint;
	info->Total_Load_to_mesh[tm + 1] = info->Total_Load_to_mesh[tm] + Total_Load;
	info->Total_DistributeForce_to_mesh[tm + 1] = info->Total_DistributeForce_to_mesh[tm] + Total_DistributeForce;

	printf("Total_Constraint = %d\n", Total_Constraint);
	printf("Total_Load = %d\n", Total_Load);
	printf("Total_DistributedForce = %d\n", Total_DistributeForce);

	fclose(fp);
	free(CP);
}


// Read file 2nd time
void Get_Input_2(int tm, char **argv, information *info)
{
	char s[256];
	int temp_i;
	double temp_d;

	int i, j, k;

	fp = fopen(argv[tm + 1], "r");

	// DIMENSION(スキップ)
	fscanf(fp, "%d", &temp_i);
	fgets(s, 256, fp);

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
		for (j = 0; j < info->DIMENSION; j++)
		{
			fscanf(fp, "%d", &temp_i);
			info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j] = temp_i;
			if (MAX_ORDER < temp_i)
			{
				MAX_ORDER = temp_i;
			}
			printf("Order[%d] = %d\n", (i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j, info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j]);
		}
	}
	fgets(s, 256, fp);

	// ノット数
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			fscanf(fp, "%d", &temp_i);
			info->No_knot[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j] = temp_i;
			info->Total_Knot_to_patch_dim[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j + 1] = info->Total_Knot_to_patch_dim[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j] + temp_i;
			printf("No_knot[%d] = %d\n", (i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j, info->No_knot[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j]);
		}
	}
	fgets(s, 256, fp);

	// 各パッチ各方向のコントロールポイント数
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			fscanf(fp, "%d", &temp_i);
			info->No_Control_point[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j] = temp_i;
			printf("No_Control_point[%d] = %d\n", (i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j, info->No_Control_point[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j]);
		}
	}
	fgets(s, 256, fp);

	for (i = 0; i < No_Patch; i++)
	{
		info->No_Control_point_in_patch[i + info->Total_Patch_to_mesh[tm]] = 1;
	}

	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			info->No_Control_point_in_patch[i + info->Total_Patch_to_mesh[tm]] *= info->No_Control_point[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j];
		}
	}

	for (i = 0; i < No_Patch; i++)
	{
		info->Total_Control_Point_to_patch[i + info->Total_Patch_to_mesh[tm] + 1] = info->Total_Control_Point_to_patch[i + info->Total_Patch_to_mesh[tm]] + info->No_Control_point_in_patch[i + info->Total_Patch_to_mesh[tm]];
	}

	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			if (info->No_knot[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j] != info->No_Control_point[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j] + info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j] + 1)
			{
				printf("wrong relationship between the number of knot vector and the number of control_point \n");
				printf("in mesh_No.%d in patch_No.%d direction:%d\n", tm, i, j);
			}
		}
	}

	for (i = 0; i < No_Patch; i++)
	{
		printf("No_Control_point_in_patch[%d] = %d\t", i + info->Total_Patch_to_mesh[tm], info->No_Control_point_in_patch[i + info->Total_Patch_to_mesh[tm]]);
	}
	printf("\n");

	// パッチコネクティビティ
	for (i = 0; i < No_Patch; i++)
	{
		for (j = 0; j < info->No_Control_point_in_patch[i + info->Total_Patch_to_mesh[tm]]; j++)
		{
			fscanf(fp, "%d", &temp_i);
			info->Patch_Control_point[info->Total_Control_Point_to_patch[i + info->Total_Patch_to_mesh[tm]] + j] = temp_i;
			if (tm > 0)
			{
				info->Patch_Control_point[info->Total_Control_Point_to_patch[i + info->Total_Patch_to_mesh[tm]] + j] += info->Total_Control_Point_to_mesh[tm];
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
		for (j = 0; j < info->DIMENSION; j++)
		{
			for (k = 0; k < info->No_knot[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j]; k++)
			{
				fscanf(fp, "%lf", &temp_d);
				info->Position_Knots[info->Total_Knot_to_patch_dim[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j] + k] = temp_d;
				printf("%le\t", info->Position_Knots[info->Total_Knot_to_patch_dim[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + j] + k]);
			}
			printf("\n");
		}
	}
	fgets(s, 256, fp);

	int Total_Element = 0;
	for (i = 0; i < No_Patch; i++)
	{
		if (info->DIMENSION == 2)
		{
			Total_Element += (info->No_Control_point[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 0] - info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 0])
						   * (info->No_Control_point[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 1] - info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 1]);
			info->No_Control_point_ON_ELEMENT[i + info->Total_Patch_to_mesh[tm]] = (info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 0] + 1) * (info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 1] + 1);
		}
		else if (info->DIMENSION == 3)
		{
			Total_Element += (info->No_Control_point[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 0] - info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 0])
						   * (info->No_Control_point[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 1] - info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 1])
						   * (info->No_Control_point[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 2] - info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 2]);
			info->No_Control_point_ON_ELEMENT[i + info->Total_Patch_to_mesh[tm]] = (info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 0] + 1) * (info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 1] + 1) * (info->Order[(i + info->Total_Patch_to_mesh[tm]) * info->DIMENSION + 2] + 1);
		}
	}
	printf("Total_Element = %d\n", Total_Element);
	info->Total_Element_on_mesh[tm] = Total_Element;
	info->Total_Element_to_mesh[tm + 1] = info->Total_Element_to_mesh[tm] + Total_Element;
	printf("Total_Element_on_mesh[%d] = %d\n", tm, info->Total_Element_on_mesh[tm]);

	for (i = 0; i < No_Patch; i++)
	{
		printf("No_Control_point_ON_ELEMENT[%d] = %d\n", i + info->Total_Patch_to_mesh[tm], info->No_Control_point_ON_ELEMENT[i + info->Total_Patch_to_mesh[tm]]);
	}

	// 節点座標
	for (i = 0; i < Total_Control_Point; i++)
	{
		fscanf(fp, "%d", &temp_i);
		for (j = 0; j < info->DIMENSION + 1; j++)
		{
			fscanf(fp, "%lf", &temp_d);
			info->Node_Coordinate[(temp_i + info->Total_Control_Point_to_mesh[tm]) * (info->DIMENSION + 1) + j] = temp_d;
		}
	}
	for (i = 0; i < Total_Control_Point; i++)
	{
		for (j = 0; j < info->DIMENSION + 1; j++)
		{
			if (info->DIMENSION == 2)
			{
				if (j == 0)
				{
					info->Control_Coord_x[i + info->Total_Control_Point_to_mesh[tm]] = info->Node_Coordinate[(i + info->Total_Control_Point_to_mesh[tm]) * (info->DIMENSION + 1) + j];
				}
				else if (j == 1)
				{
					info->Control_Coord_y[i + info->Total_Control_Point_to_mesh[tm]] = info->Node_Coordinate[(i + info->Total_Control_Point_to_mesh[tm]) * (info->DIMENSION + 1) + j];
				}
				else if (j == info->DIMENSION)
				{
					info->Control_Weight[i + info->Total_Control_Point_to_mesh[tm]] = info->Node_Coordinate[(i + info->Total_Control_Point_to_mesh[tm]) * (info->DIMENSION + 1) + j];
				}
			}
			else if (info->DIMENSION == 3)
			{
				if (j == 0)
				{
					info->Control_Coord_x[i + info->Total_Control_Point_to_mesh[tm]] = info->Node_Coordinate[(i + info->Total_Control_Point_to_mesh[tm]) * (info->DIMENSION + 1) + j];
				}
				else if (j == 1)
				{
					info->Control_Coord_y[i + info->Total_Control_Point_to_mesh[tm]] = info->Node_Coordinate[(i + info->Total_Control_Point_to_mesh[tm]) * (info->DIMENSION + 1) + j];
				}
				else if (j == 2)
				{
					info->Control_Coord_z[i + info->Total_Control_Point_to_mesh[tm]] = info->Node_Coordinate[(i + info->Total_Control_Point_to_mesh[tm]) * (info->DIMENSION + 1) + j];
				}
				else if (j == info->DIMENSION)
				{
					info->Control_Weight[i + info->Total_Control_Point_to_mesh[tm]] = info->Node_Coordinate[(i + info->Total_Control_Point_to_mesh[tm]) * (info->DIMENSION + 1) + j];
				}
			}
		}
	}
	fgets(s, 256, fp);

	// 拘束
	printf("\nConstraint_Node\n");
	for (i = 0; i < Total_Constraint; i++)
	{
		fscanf(fp, "%d %d %lf",
			   &info->Constraint_Node_Dir[(i + info->Total_Constraint_to_mesh[tm]) * 2 + 0],
			   &info->Constraint_Node_Dir[(i + info->Total_Constraint_to_mesh[tm]) * 2 + 1],
			   &info->Value_of_Constraint[i + info->Total_Constraint_to_mesh[tm]]);
		info->Constraint_Node_Dir[(i + info->Total_Constraint_to_mesh[tm]) * 2 + 0] = info->Constraint_Node_Dir[(i + info->Total_Constraint_to_mesh[tm]) * 2 + 0] + info->Total_Control_Point_to_mesh[tm];

		printf("%d %d %le \n",
			   info->Constraint_Node_Dir[(i + info->Total_Constraint_to_mesh[tm]) * 2 + 0],
			   info->Constraint_Node_Dir[(i + info->Total_Constraint_to_mesh[tm]) * 2 + 1],
			   info->Value_of_Constraint[i + info->Total_Constraint_to_mesh[tm]]);
	}
	fgets(s, 256, fp);

	// 荷重
	for (i = 0; i < Total_Load; i++)
	{
		fscanf(fp, "%d %d %lf",
			   &info->Load_Node_Dir[(i + info->Total_Load_to_mesh[tm]) * 2 + 0],
			   &info->Load_Node_Dir[(i + info->Total_Load_to_mesh[tm]) * 2 + 1],
			   &info->Value_of_Load[i + info->Total_Load_to_mesh[tm]]);
		info->Load_Node_Dir[(i + info->Total_Load_to_mesh[tm]) * 2 + 0] = info->Load_Node_Dir[(i + info->Total_Load_to_mesh[tm]) * 2 + 0] + info->Total_Control_Point_to_mesh[tm];

		printf("Load_Node_Dir[%d]= %d Load_Node_Dir[%d]= %d Value_of_Load[%d] = %le\n",
			   (i + info->Total_Load_to_mesh[tm]) * 2 + 0, info->Load_Node_Dir[(i + info->Total_Load_to_mesh[tm]) * 2 + 0],
			   (i + info->Total_Load_to_mesh[tm]) * 2 + 1, info->Load_Node_Dir[(i + info->Total_Load_to_mesh[tm]) * 2 + 1],
			   i + info->Total_Load_to_mesh[tm], info->Value_of_Load[i + info->Total_Load_to_mesh[tm]]);
	}
	fgets(s, 256, fp);

	if (info->DIMENSION == 2)
	{
		int type_load, iPatch, iCoord;
		double val_Coord, Range_Coord[2], Coeff_Dist_Load[3];

		for (i = 0; i < Total_DistributeForce; i++)
		{
			fscanf(fp, "%d %d %d %lf %lf %lf %lf %lf %lf", &type_load, &iPatch, &iCoord, &val_Coord, &Range_Coord[0], &Range_Coord[1], &Coeff_Dist_Load[0], &Coeff_Dist_Load[1], &Coeff_Dist_Load[2]);
			printf("Distibuted load nober: %d\n", i);
			printf("type_load: %d iPatch: %d iCoord: %d val_Coord: %le Range_Coord: %le %le\nCoef_Dist_Load: %le %le %le\n",
					type_load, iPatch, iCoord, val_Coord, Range_Coord[0], Range_Coord[1], Coeff_Dist_Load[0], Coeff_Dist_Load[1], Coeff_Dist_Load[2]);

			// for S-IGA
			info->type_load_array[i + info->Total_DistributeForce_to_mesh[tm]] = type_load;
			info->iPatch_array[i + info->Total_DistributeForce_to_mesh[tm]] = iPatch;
			info->iCoord_array[i + info->Total_DistributeForce_to_mesh[tm]] = iCoord;
			info->val_Coord_array[i + info->Total_DistributeForce_to_mesh[tm]] = val_Coord;
			info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 2 + 0] = Range_Coord[0];
			info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 2 + 1] = Range_Coord[1];
			info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 3 + 0] = Coeff_Dist_Load[0];
			info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 3 + 1] = Coeff_Dist_Load[1];
			info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 3 + 2] = Coeff_Dist_Load[2];
		}
	}
	else if (info->DIMENSION == 3)
	{
		int type_load, iPatch, iCoord, jCoord;
		double val_Coord, iRange_Coord[2], jRange_Coord[2], iCoeff_Dist_Load[3], jCoeff_Dist_Load[3];

		for (i = 0; i < Total_DistributeForce; i++)
		{
			fscanf(fp, "%d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &type_load, &iPatch, &iCoord, &jCoord, &val_Coord, &iRange_Coord[0], &iRange_Coord[1], &jRange_Coord[0], &jRange_Coord[1], &iCoeff_Dist_Load[0], &iCoeff_Dist_Load[1], &iCoeff_Dist_Load[2], &jCoeff_Dist_Load[0], &jCoeff_Dist_Load[1], &jCoeff_Dist_Load[2]);
			printf("Distibuted load nober: %d\n", i);
			printf("type_load: %d iPatch: %d iCoord: %d jCoord: %d val_Coord: %le\nRange_Coord: %le %le %le %le\nCoef_Dist_Load: %le %le %le %le %le %le\n",
					type_load, iPatch, iCoord, jCoord, val_Coord,
					iRange_Coord[0], iRange_Coord[1], jRange_Coord[0], jRange_Coord[1],
					iCoeff_Dist_Load[0], iCoeff_Dist_Load[1], iCoeff_Dist_Load[2], jCoeff_Dist_Load[0], jCoeff_Dist_Load[1], jCoeff_Dist_Load[2]);

			// for S-IGA
			info->type_load_array[i + info->Total_DistributeForce_to_mesh[tm]] = type_load;
			info->iPatch_array[i + info->Total_DistributeForce_to_mesh[tm]] = iPatch;
			info->iCoord_array[i + info->Total_DistributeForce_to_mesh[tm]] = iCoord;
			info->jCoord_array[i + info->Total_DistributeForce_to_mesh[tm]] = jCoord;
			info->val_Coord_array[i + info->Total_DistributeForce_to_mesh[tm]] = val_Coord;
			info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 4 + 0] = iRange_Coord[0];
			info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 4 + 1] = iRange_Coord[1];
			info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 4 + 2] = jRange_Coord[0];
			info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 4 + 3] = jRange_Coord[1];
			info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 0] = iCoeff_Dist_Load[0];
			info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 1] = iCoeff_Dist_Load[1];
			info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 2] = iCoeff_Dist_Load[2];
			info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 3] = jCoeff_Dist_Load[0];
			info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 4] = jCoeff_Dist_Load[1];
			info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 5] = jCoeff_Dist_Load[2];
		}
	}
	fclose(fp);
}


// INC 等の作成
void Make_INC(information *info)
{
	int tm;

	// info->INC の計算(節点番号をξ, ηの番号で表す為の配列)
	for (tm = 0; tm < Total_mesh; tm++)
	{
		int b, B, e, h, i, j, k, l, n, o, p, q, x, y, z, ii, jj, kk, kkk, iiloc, jjloc, kkloc, r = 0;
		int type_load, iPatch, iCoord, jCoord;
		double val_Coord, iRange_Coord[2], jRange_Coord[2], iCoeff_Dist_Load[3], jCoeff_Dist_Load[3];
		int No_Patch = info->Total_Patch_on_mesh[tm];
		int Total_Patch_to_Now = info->Total_Patch_to_mesh[tm];
		int Total_Element = info->Total_Element_on_mesh[tm];
		int Total_Element_to_Now = info->Total_Element_to_mesh[tm];
		int Total_DistributeForce = info->Total_DistributeForce_to_mesh[tm + 1] - info->Total_DistributeForce_to_mesh[tm];

		if (info->DIMENSION == 2)
		{
			e = 0;
			for (l = 0; l < No_Patch; l++)
			{
				i = 0;
				int patch = l + Total_Patch_to_Now;
				for (jj = 0; jj < info->No_Control_point[patch * info->DIMENSION + 1]; jj++)
				{
					for (ii = 0; ii < info->No_Control_point[patch * info->DIMENSION + 0]; ii++)
					{
						info->INC[patch * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Patch_Control_point[info->Total_Control_Point_to_patch[patch] + i] * info->DIMENSION + 0] = ii;
						info->INC[patch * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Patch_Control_point[info->Total_Control_Point_to_patch[patch] + i] * info->DIMENSION + 1] = jj;

						if (ii >= info->Order[patch * info->DIMENSION + 0] && jj >= info->Order[patch * info->DIMENSION + 1])
						{
							for (jjloc = 0; jjloc <= info->Order[patch * info->DIMENSION + 1]; jjloc++)
							{
								for (iiloc = 0; iiloc <= info->Order[patch * info->DIMENSION + 0]; iiloc++)
								{
									B = info->Patch_Control_point[info->Total_Control_Point_to_patch[patch] + i - jjloc * info->No_Control_point[patch * info->DIMENSION + 0] - iiloc];
									b = jjloc * (info->Order[patch * info->DIMENSION + 0] + 1) + iiloc;
									info->Controlpoint_of_Element[(e + Total_Element_to_Now) * MAX_NO_CP_ON_ELEMENT + b] = B;
								}
							}
							info->Element_patch[e + Total_Element_to_Now] = patch;
							info->Element_mesh[e + Total_Element_to_Now] = tm;
							e++;
						}
						i++;
					}
				}
			}
		}
		else if (info->DIMENSION == 3)
		{
			e = 0;
			for (l = 0; l < No_Patch; l++)
			{
				i = 0;
				int patch = l + Total_Patch_to_Now;
				for (kk = 0; kk < info->No_Control_point[patch * info->DIMENSION + 2]; kk++)
				{
					for (jj = 0; jj < info->No_Control_point[patch * info->DIMENSION + 1]; jj++)
					{
						for (ii = 0; ii < info->No_Control_point[patch * info->DIMENSION + 0]; ii++)
						{
							info->INC[patch * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Patch_Control_point[info->Total_Control_Point_to_patch[patch] + i] * info->DIMENSION + 0] = ii;
							info->INC[patch * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Patch_Control_point[info->Total_Control_Point_to_patch[patch] + i] * info->DIMENSION + 1] = jj;
							info->INC[patch * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Patch_Control_point[info->Total_Control_Point_to_patch[patch] + i] * info->DIMENSION + 2] = kk;

							if (ii >= info->Order[patch * info->DIMENSION + 0] && jj >= info->Order[patch * info->DIMENSION + 1] && kk >= info->Order[patch * info->DIMENSION + 2])
							{
								for (kkloc = 0; kkloc <= info->Order[patch * info->DIMENSION + 2]; kkloc++)
								{
									for (jjloc = 0; jjloc <= info->Order[patch * info->DIMENSION + 1]; jjloc++)
									{
										for (iiloc = 0; iiloc <= info->Order[patch * info->DIMENSION + 0]; iiloc++)
										{
											B = info->Patch_Control_point[info->Total_Control_Point_to_patch[patch] + i - kkloc * info->No_Control_point[patch * info->DIMENSION + 0] * info->No_Control_point[patch * info->DIMENSION + 1] - jjloc * info->No_Control_point[patch * info->DIMENSION + 0] - iiloc];
											b = kkloc * (info->Order[patch * info->DIMENSION + 0] + 1) * (info->Order[patch * info->DIMENSION + 1] + 1) + jjloc * (info->Order[patch * info->DIMENSION + 0] + 1) + iiloc;
											info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + b] = B;
										}
									}
								}
								info->Element_patch[e + Total_Element_to_Now] = patch;
								info->Element_mesh[e + Total_Element_to_Now] = tm;
								e++;
							}
							i++;
						}
					}
				}
			}
		}

		// for S-IGA line_No_real_elementの初期化
		for (l = 0; l < No_Patch; l++)
		{
			int patch = l + Total_Patch_to_Now;
			for (j = 0; j < info->DIMENSION; j++)
			{
				info->line_No_real_element[patch * info->DIMENSION + j] = 0;
			}
		}

		for (l = 0; l < No_Patch; l++)
		{
			int patch = l + Total_Patch_to_Now;
			for (j = 0; j < info->DIMENSION; j++)
			{
				info->line_No_Total_element[patch * info->DIMENSION + j] = info->No_knot[patch * info->DIMENSION + j] - 2 * info->Order[patch * info->DIMENSION + j] - 1;

				for (kkk = info->Order[patch * info->DIMENSION + j]; kkk < info->No_knot[patch * info->DIMENSION + j] - info->Order[patch * info->DIMENSION + j] - 1; kkk++)
				{
					info->difference[info->Total_Knot_to_patch_dim[patch * info->DIMENSION + j] + kkk - info->Order[patch * info->DIMENSION + j]]
						= info->Position_Knots[info->Total_Knot_to_patch_dim[patch * info->DIMENSION + j] + kkk + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[patch * info->DIMENSION + j] + kkk];
					if (info->difference[info->Total_Knot_to_patch_dim[patch * info->DIMENSION + j] + kkk - info->Order[patch * info->DIMENSION + j]] != 0)
					{
						info->line_No_real_element[patch * info->DIMENSION + j]++;
					}
				}
			}
		}

		// 要素に行番号, 列番号をつける
		for (h = 0; h < Total_Element; h++)
		{
			info->Total_element_all_ID[h] = 0;
		}

		if (info->DIMENSION == 2)
		{
			i = 0;
			for (l = 0; l < No_Patch; l++)
			{
				int patch = l + Total_Patch_to_Now;
				for (y = 0; y < info->line_No_Total_element[patch * info->DIMENSION + 1]; y++)
				{
					for (x = 0; x < info->line_No_Total_element[patch * info->DIMENSION + 0]; x++)
					{
						info->ENC[(i + info->Total_Element_to_mesh[tm]) * info->DIMENSION + 0] = x;
						info->ENC[(i + info->Total_Element_to_mesh[tm]) * info->DIMENSION + 1] = y;
						i++;
					}
				}
			}
		}
		else if (info->DIMENSION == 3)
		{
			i = 0;
			for (l = 0; l < No_Patch; l++)
			{
				int patch = l + Total_Patch_to_Now;
				for (z = 0; z < info->line_No_Total_element[patch * info->DIMENSION + 2]; z++)
				{
					for (y = 0; y < info->line_No_Total_element[patch * info->DIMENSION + 1]; y++)
					{
						for (x = 0; x < info->line_No_Total_element[patch * info->DIMENSION + 0]; x++)
						{
							info->ENC[(i + info->Total_Element_to_mesh[tm]) * info->DIMENSION + 0] = x;
							info->ENC[(i + info->Total_Element_to_mesh[tm]) * info->DIMENSION + 1] = y;
							info->ENC[(i + info->Total_Element_to_mesh[tm]) * info->DIMENSION + 2] = z;
							i++;
						}
					}
				}
			}
		}

		// 必要な要素の行と列の番号を求める
		for (j = 0; j < info->DIMENSION; j++)
		{
			for (l = 0; l < No_Patch; l++)
			{
				e = 0;
				int patch = l + Total_Patch_to_Now;
				for (k = 0; k < info->line_No_Total_element[patch * info->DIMENSION + j]; k++)
				{
					if (info->difference[info->Total_Knot_to_patch_dim[patch * info->DIMENSION + j] + k] != 0)
					{
						info->real_element_line[patch * (info->Total_Element_to_mesh[Total_mesh] * info->DIMENSION) + e * info->DIMENSION + j] = k;
						e++;
					}
				}
			}
		}

		// 必要な要素列上の要素のIDを1にする
		if (info->DIMENSION == 2)
		{
			for (n = 0; n < Total_Element; n++)
			{
				int ele = n + Total_Element_to_Now;
				for (p = 0; p < info->line_No_real_element[info->Element_patch[ele] * info->DIMENSION + 0]; p++)
				{
					if (info->ENC[ele * info->DIMENSION + 0] == info->real_element_line[info->Element_patch[ele] * (info->Total_Element_to_mesh[Total_mesh] * info->DIMENSION) + p * info->DIMENSION + 0])
					{
						for (q = 0; q < info->line_No_real_element[info->Element_patch[ele] * info->DIMENSION + 1]; q++)
						{
							if (info->ENC[ele * info->DIMENSION + 1] == info->real_element_line[info->Element_patch[ele] * (info->Total_Element_to_mesh[Total_mesh] * info->DIMENSION) + q * info->DIMENSION + 1])
							{
								info->Total_element_all_ID[n]++;
							}
						}
					}
				}

				// IDが1の要素に番号を振る
				if (info->Total_element_all_ID[n] == 1)
				{
					info->real_element[r + info->real_Total_Element_to_mesh[tm]] = n + Total_Element_to_Now;
					info->real_El_No_on_mesh[tm * info->Total_Element_to_mesh[Total_mesh] + r] = n + Total_Element_to_Now;
					r++;
				}
			}

			// for S-IGA real_Total_Elementの初期化
			int real_Total_Element = 0;

			for (l = 0; l < No_Patch; l++)
			{
				int patch = l + Total_Patch_to_Now;
				real_Total_Element += info->line_No_real_element[patch * info->DIMENSION + 0] * info->line_No_real_element[patch * info->DIMENSION + 1];
			}
			info->real_Total_Element_on_mesh[tm] = real_Total_Element;
			info->real_Total_Element_to_mesh[tm + 1] = info->real_Total_Element_to_mesh[tm] + real_Total_Element;
		}
		else if (info->DIMENSION == 3)
		{
			for (n = 0; n < Total_Element; n++)
			{
				int ele = n + Total_Element_to_Now;
				for (o = 0; o < info->line_No_real_element[info->Element_patch[ele] * info->DIMENSION + 0]; o++)
				{
					if (info->ENC[ele * info->DIMENSION + 0] == info->real_element_line[info->Element_patch[ele] * (info->Total_Element_to_mesh[Total_mesh] * info->DIMENSION) + o * info->DIMENSION + 0])
					{
						for (p = 0; p < info->line_No_real_element[info->Element_patch[ele] * info->DIMENSION + 1]; p++)
						{
							if (info->ENC[ele * info->DIMENSION + 1] == info->real_element_line[info->Element_patch[ele] * (info->Total_Element_to_mesh[Total_mesh] * info->DIMENSION) + p * info->DIMENSION + 1])
							{
								for (q = 0; q < info->line_No_real_element[info->Element_patch[ele] * info->DIMENSION + 2]; q++)
								{
									if (info->ENC[ele * info->DIMENSION + 2] == info->real_element_line[info->Element_patch[ele] * (info->Total_Element_to_mesh[Total_mesh] * info->DIMENSION) + q * info->DIMENSION + 2])
									{
										info->Total_element_all_ID[n]++;
									}
								}
							}
						}
					}
				}

				// IDが1の要素に番号を振る
				if (info->Total_element_all_ID[n] == 1)
				{
					info->real_element[r + info->real_Total_Element_to_mesh[tm]] = ele;
					info->real_El_No_on_mesh[tm * info->Total_Element_to_mesh[Total_mesh] + r] = ele;
					r++;
				}
			}

			// for S-IGA real_Total_Elementの初期化
			int real_Total_Element = 0;

			for (l = 0; l < No_Patch; l++)
			{
				int patch = l + Total_Patch_to_Now;
				real_Total_Element += info->line_No_real_element[patch * info->DIMENSION + 0] * info->line_No_real_element[patch * info->DIMENSION + 1] * info->line_No_real_element[patch * info->DIMENSION + 2];
			}
			info->real_Total_Element_on_mesh[tm] = real_Total_Element;
			info->real_Total_Element_to_mesh[tm + 1] = info->real_Total_Element_to_mesh[tm] + real_Total_Element;
		}

		// 分布荷重
		for (i = 0; i < Total_DistributeForce; i++)
		{
			if (info->DIMENSION == 2)
			{
				type_load = info->type_load_array[i + info->Total_DistributeForce_to_mesh[tm]];
				iPatch = info->iPatch_array[i + info->Total_DistributeForce_to_mesh[tm]];
				iCoord = info->iCoord_array[i + info->Total_DistributeForce_to_mesh[tm]];
				val_Coord = info->val_Coord_array[i + info->Total_DistributeForce_to_mesh[tm]];
				iRange_Coord[0] = info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 2 + 0];
				iRange_Coord[1] = info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 2 + 1];
				iCoeff_Dist_Load[0] = info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 3 + 0];
				iCoeff_Dist_Load[1] = info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 3 + 1];
				iCoeff_Dist_Load[2] = info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 3 + 2];

				Setting_Dist_Load_2D(tm, iPatch, iCoord, val_Coord, iRange_Coord, type_load, iCoeff_Dist_Load, info);
			}
			else if (info->DIMENSION == 3)
			{
				type_load = info->type_load_array[i + info->Total_DistributeForce_to_mesh[tm]];
				iPatch = info->iPatch_array[i + info->Total_DistributeForce_to_mesh[tm]];
				iCoord = info->iCoord_array[i + info->Total_DistributeForce_to_mesh[tm]];
				jCoord = info->jCoord_array[i + info->Total_DistributeForce_to_mesh[tm]];
				val_Coord = info->val_Coord_array[i + info->Total_DistributeForce_to_mesh[tm]];
				iRange_Coord[0] = info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 4 + 0];
				iRange_Coord[1] = info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 4 + 1];
				jRange_Coord[0] = info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 4 + 2];
				jRange_Coord[1] = info->Range_Coord_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 4 + 3];
				iCoeff_Dist_Load[0] = info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 0];
				iCoeff_Dist_Load[1] = info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 1];
				iCoeff_Dist_Load[2] = info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 2];
				jCoeff_Dist_Load[0] = info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 3];
				jCoeff_Dist_Load[1] = info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 4];
				jCoeff_Dist_Load[2] = info->Coeff_Dist_Load_array[(i + info->Total_DistributeForce_to_mesh[tm]) * 6 + 5];

				Setting_Dist_Load_3D(tm, iPatch, iCoord, jCoord, val_Coord, iRange_Coord, jRange_Coord, type_load, iCoeff_Dist_Load, jCoeff_Dist_Load, info);
			}
		}
	}
}


// Distributed Load
void Setting_Dist_Load_2D(int mesh_n, int iPatch, int iCoord, double val_Coord, double *Range_Coord, int type_load, double *Coeff_Dist_Load, information *info)
{
	int iii, jjj;
	int N_Seg_Load_Element_iDir = 0, jCoord = 0;
	int iPos[2] = {-10000, -10000}, jPos[2] = {-10000, -10000};
	int No_Element_For_Dist_Load;
	int iX, iY;
	int ic, ig;
	double val_jCoord_Local = 0.0;
	double Position_Data_param[MAX_DIMENSION] = {0.0};

	int *No_Element_for_Integration = (int *)malloc(sizeof(int) * info->Total_Knot_to_mesh[Total_mesh]); // No_Element_for_Integration[MAX_N_KNOT]
	int *iControlpoint = (int *)malloc(sizeof(int) * MAX_NO_CP_ON_ELEMENT);

	Make_gauss_array(0, info);

	// iCoord=0: Load on Eta=Constant
	// iCoord=1: Load on Xi=Constant
	if (iCoord == 0)
	{
		jCoord = 1;
	}
	else if (iCoord == 1)
	{
		jCoord = 0;
	}

	for (iii = info->Order[iPatch * info->DIMENSION + iCoord]; iii < info->No_knot[iPatch * info->DIMENSION + iCoord] - info->Order[iPatch * info->DIMENSION + iCoord] - 1; iii++)
	{
		double epsi = 0.00000000001;

		if (info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + iCoord] + iii] - epsi <= Range_Coord[0])
			iPos[0] = iii;
		if (info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + iCoord] + iii + 1] - epsi <= Range_Coord[1])
			iPos[1] = iii + 1;
	}

	if (iPos[0] < 0 || iPos[1] < 0)
	{
		printf("Error (Stop) iPos[0] = %d iPos[1] = %d\n", iPos[0], iPos[1]);
		exit(0);
	}

	for (jjj = info->Order[iPatch * info->DIMENSION + jCoord]; jjj < info->No_knot[iPatch * info->DIMENSION + jCoord] - info->Order[iPatch * info->DIMENSION + jCoord] - 1; jjj++)
	{
		double epsi = 0.00000000001;

		if (info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + jCoord] + jjj] - epsi <= val_Coord
			&& info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + jCoord] + jjj + 1] + epsi > val_Coord)
		{
			jPos[0] = jjj;
			jPos[1] = jjj + 1;
			val_jCoord_Local = -1.0 + 2.0 * (val_Coord - info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + jCoord] + jjj])
							 / (info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + jCoord] + jjj + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + jCoord] + jjj]);
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
		iX = jPos[0] - info->Order[iPatch * info->DIMENSION + 0];
		for (iY = iPos[0] - info->Order[iPatch * info->DIMENSION + 1]; iY < iPos[1] - info->Order[iPatch * info->DIMENSION + 1]; iY++)
		{
			No_Element_for_Integration[iii] = SearchForElement_2D(mesh_n, iPatch, iX, iY, info);
			iii++;
		}
	}

	if (iCoord == 0)
	{
		iY = jPos[0] - info->Order[iPatch * info->DIMENSION + 1];
		for (iX = iPos[0] - info->Order[iPatch * info->DIMENSION + 0]; iX < iPos[1] - info->Order[iPatch * info->DIMENSION + 0]; iX++)
		{
			No_Element_for_Integration[iii] = SearchForElement_2D(mesh_n, iPatch, iX, iY, info);
			iii++;
		}
	}
	No_Element_For_Dist_Load = iii;

	// Book keeping finished

	for (iii = 0; iii < No_Element_For_Dist_Load; iii++)
	{
		if (info->Total_element_all_ID[No_Element_for_Integration[iii]] == 1)
		{
			iX = info->ENC[No_Element_for_Integration[iii] * info->DIMENSION + 0];
			iY = info->ENC[No_Element_for_Integration[iii] * info->DIMENSION + 1];

			for (ic = 0; ic < (info->Order[iPatch * info->DIMENSION + 0] + 1) * (info->Order[iPatch * info->DIMENSION + 1] + 1); ic++)
				iControlpoint[ic] = info->Controlpoint_of_Element[No_Element_for_Integration[iii] * MAX_NO_CP_ON_ELEMENT + ic];

			for (ig = 0; ig < GP_1D; ig++)
			{
				double Local_Coord[2], sfc, dxyzdge[3], detJ, XiEtaCoordParen, valDistLoad;
				int icc;
				Local_Coord[jCoord] = val_jCoord_Local;
				Local_Coord[iCoord] = Gxi_1D[ig];

				ShapeFunc_from_paren(Position_Data_param, Local_Coord, iCoord, No_Element_for_Integration[iii], info);
				XiEtaCoordParen = Position_Data_param[iCoord];
				valDistLoad = Coeff_Dist_Load[0] + Coeff_Dist_Load[1] * XiEtaCoordParen + Coeff_Dist_Load[2] * XiEtaCoordParen * XiEtaCoordParen;

				dxyzdge[0] = 0.0;
				dxyzdge[1] = 0.0;
				dxyzdge[2] = 0.0;
				for (icc = 0; icc < (info->Order[iPatch * info->DIMENSION + 0] + 1) * (info->Order[iPatch * info->DIMENSION + 1] + 1); icc++)
				{
					dxyzdge[0] += dShape_func(icc, iCoord, Local_Coord, No_Element_for_Integration[iii], info) * info->Node_Coordinate[iControlpoint[icc] * (info->DIMENSION + 1) + 0];
					dxyzdge[1] += dShape_func(icc, iCoord, Local_Coord, No_Element_for_Integration[iii], info) * info->Node_Coordinate[iControlpoint[icc] * (info->DIMENSION + 1) + 1];
				}

				detJ = sqrt(dxyzdge[0] * dxyzdge[0] + dxyzdge[1] * dxyzdge[1]);

				if (type_load < 2)
				{
					for (ic = 0; ic < (info->Order[iPatch * info->DIMENSION + 0] + 1) * (info->Order[iPatch * info->DIMENSION + 1] + 1); ic++)
					{
						sfc = Shape_func(ic, Local_Coord, No_Element_for_Integration[iii], info);
						info->Equivalent_Nodal_Force[iControlpoint[ic] * info->DIMENSION + type_load] += valDistLoad * sfc * detJ * w_1D[ig];
					}
				}
				else if (type_load == 2) // 法線方向
				{
					double LoadDir[2];
					LoadDir[0] = dxyzdge[1] / detJ;
					LoadDir[1] = -dxyzdge[0] / detJ;
					for (ic = 0; ic < (info->Order[iPatch * info->DIMENSION + 0] + 1) * (info->Order[iPatch * info->DIMENSION + 1] + 1); ic++)
					{
						sfc = Shape_func(ic, Local_Coord, No_Element_for_Integration[iii], info);
						info->Equivalent_Nodal_Force[iControlpoint[ic] * info->DIMENSION + 0] += LoadDir[0] * valDistLoad * sfc * detJ * w_1D[ig];
						info->Equivalent_Nodal_Force[iControlpoint[ic] * info->DIMENSION + 1] += LoadDir[1] * valDistLoad * sfc * detJ * w_1D[ig];
					}
				}
				else if (type_load == 3)
				{
					double LoadDir[2];
					LoadDir[0] = dxyzdge[0] / detJ;
					LoadDir[1] = dxyzdge[1] / detJ;
					for (ic = 0; ic < (info->Order[iPatch * info->DIMENSION + 0] + 1) * (info->Order[iPatch * info->DIMENSION + 1] + 1); ic++)
					{
						sfc = Shape_func(ic, Local_Coord, No_Element_for_Integration[iii], info);
						info->Equivalent_Nodal_Force[iControlpoint[ic] * info->DIMENSION + 0] += LoadDir[0] * valDistLoad * sfc * detJ * w_1D[ig];
						info->Equivalent_Nodal_Force[iControlpoint[ic] * info->DIMENSION + 1] += LoadDir[1] * valDistLoad * sfc * detJ * w_1D[ig];
					}
				}
			}
		}
	}
	free(No_Element_for_Integration), free(iControlpoint);
}


void Setting_Dist_Load_3D(int mesh_n, int iPatch, int iCoord, int jCoord, double val_Coord, double *iRange_Coord, double *jRange_Coord, int type_load, double *iCoeff_Dist_Load, double *jCoeff_Dist_Load, information *info)
{
	int iii, jjj, kkk;
	int N_Seg_Load_Element_iDir = 0, N_Seg_Load_Element_jDir = 0, kCoord = 0;
	int iPos[2] = {-10000, -10000}, jPos[2] = {-10000, -10000}, kPos[2] = {-10000, -10000};
	int No_Element_For_Dist_Load_iDir, No_Element_For_Dist_Load_jDir;
	int iX, iY, iZ;
	int ic, ig_i, ig_j;
	double val_kCoord_Local = 0.0;
	double Position_Data_param[MAX_DIMENSION] = {0.0};

	int *No_Element_for_Integration = (int *)malloc(sizeof(int) * info->Total_Knot_to_mesh[Total_mesh] * info->Total_Knot_to_mesh[Total_mesh]); // No_Element_for_Integration[MAX_N_KNOT]
	int *iControlpoint = (int *)malloc(sizeof(int) * MAX_NO_CP_ON_ELEMENT);

	Make_gauss_array(0, info);

	if (iCoord == 0 && jCoord == 1)
	{
		kCoord = 2;
	}
	else if (iCoord == 1 && jCoord == 2)
	{
		kCoord = 0;
	}
	else if (iCoord == 2 && jCoord == 0)
	{
		kCoord = 1;
	}

	for (iii = info->Order[iPatch * info->DIMENSION + iCoord]; iii < info->No_knot[iPatch * info->DIMENSION + iCoord] - info->Order[iPatch * info->DIMENSION + iCoord] - 1; iii++)
	{
		double epsi = 0.00000000001;

		if (info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + iCoord] + iii] - epsi <= iRange_Coord[0])
			iPos[0] = iii;
		if (info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + iCoord] + iii + 1] - epsi <= iRange_Coord[1])
			iPos[1] = iii + 1;
	}
	if (iPos[0] < 0 || iPos[1] < 0)
	{
		printf("Error (Stop) iPos[0] = %d iPos[1] = %d\n", iPos[0], iPos[1]);
		exit(0);
	}

	for (jjj = info->Order[iPatch * info->DIMENSION + jCoord]; jjj < info->No_knot[iPatch * info->DIMENSION + jCoord] - info->Order[iPatch * info->DIMENSION + jCoord] - 1; jjj++)
	{
		double epsi = 0.00000000001;

		if (info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + jCoord] + jjj] - epsi <= jRange_Coord[0])
			jPos[0] = jjj;
		if (info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + jCoord] + jjj + 1] - epsi <= jRange_Coord[1])
			jPos[1] = jjj + 1;
	}
	if (jPos[0] < 0 || jPos[1] < 0)
	{
		printf("Error (Stop) jPos[0] = %d jPos[1] = %d\n", jPos[0], jPos[1]);
		exit(0);
	}

	for (kkk = info->Order[iPatch * info->DIMENSION + kCoord]; kkk < info->No_knot[iPatch * info->DIMENSION + kCoord] - info->Order[iPatch * info->DIMENSION + kCoord] - 1; kkk++)
	{
		double epsi = 0.00000000001;

		if (info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + kCoord] + kkk] - epsi <= val_Coord
			&& info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + kCoord] + kkk + 1] + epsi > val_Coord)
		{
			kPos[0] = kkk;
			kPos[1] = kkk + 1;
			val_kCoord_Local = - 1.0 + 2.0 * (val_Coord - info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + kCoord] + kkk])
							 / (info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + kCoord] + kkk + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[iPatch * info->DIMENSION + kCoord] + kkk]);
		}
	}
	if (kPos[0] < 0 || kPos[1] < 0)
	{
		printf("Error (Stop) kPos[0] = %d kPos[1] = %d\n", kPos[0], kPos[1]);
		exit(0);
	}

	for (iii = iPos[0]; iii < iPos[1]; iii++)
	{
		N_Seg_Load_Element_iDir++;
	}

	for (jjj = jPos[0]; jjj < jPos[1]; jjj++)
	{
		N_Seg_Load_Element_jDir++;
	}

	iii = 0, jjj = 0;

	if (iCoord == 0 && jCoord == 1)
	{
		iZ = kPos[0] - info->Order[iPatch * info->DIMENSION + 2];
		for (iX = iPos[0] - info->Order[iPatch * info->DIMENSION + 0]; iX < iPos[1] - info->Order[iPatch * info->DIMENSION + 0]; iX++)
		{
			for (iY = jPos[0] - info->Order[iPatch * info->DIMENSION + 1]; iY < jPos[1] - info->Order[iPatch * info->DIMENSION + 1]; iY++)
			{
				No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj] = SearchForElement_3D(mesh_n, iPatch, iX, iY, iZ, info);
				jjj++;
			}
			iii++;
			jjj = 0;
		}
		for (iY = jPos[0] - info->Order[iPatch * info->DIMENSION + 1]; iY < jPos[1] - info->Order[iPatch * info->DIMENSION + 1]; iY++)
		{
			jjj++;
		}
	}
	else if (iCoord == 1 && jCoord == 2)
	{
		iii = 0, jjj = 0;
		iX = kPos[0] - info->Order[iPatch * info->DIMENSION + 0];
		for (iY = iPos[0] - info->Order[iPatch * info->DIMENSION + 1]; iY < iPos[1] - info->Order[iPatch * info->DIMENSION + 1]; iY++)
		{
			for (iZ = jPos[0] - info->Order[iPatch * info->DIMENSION + 2]; iZ < jPos[1] - info->Order[iPatch * info->DIMENSION + 2]; iZ++)
			{
				No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj] = SearchForElement_3D(mesh_n, iPatch, iX, iY, iZ, info);
				jjj++;
			}
			iii++;
			jjj = 0;
		}
		for (iZ = jPos[0] - info->Order[iPatch * info->DIMENSION + 2]; iZ < jPos[1] - info->Order[iPatch * info->DIMENSION + 2]; iZ++)
		{
			jjj++;
		}
	}
	else if (iCoord == 2 && jCoord == 0)
	{
		
		iY = kPos[0] - info->Order[iPatch * info->DIMENSION + 1];
		for (iZ = iPos[0] - info->Order[iPatch * info->DIMENSION + 2]; iZ < iPos[1] - info->Order[iPatch * info->DIMENSION + 2]; iZ++)
		{
			for (iX = jPos[0] - info->Order[iPatch * info->DIMENSION + 0]; iX < jPos[1] - info->Order[iPatch * info->DIMENSION + 0]; iX++)
			{
				No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj] = SearchForElement_3D(mesh_n, iPatch, iX, iY, iZ, info);
				jjj++;
			}
			iii++;
			jjj = 0;
		}
		for (iX = jPos[0] - info->Order[iPatch * info->DIMENSION + 0]; iX < jPos[1] - info->Order[iPatch * info->DIMENSION + 0]; iX++)
		{
			jjj++;
		}
	}

	No_Element_For_Dist_Load_iDir = iii;
	No_Element_For_Dist_Load_jDir = jjj;

	// Book keeping finished

	for (iii = 0; iii < No_Element_For_Dist_Load_iDir; iii++)
	{
		for (jjj = 0; jjj < No_Element_For_Dist_Load_jDir; jjj++)
		{
			if (info->Total_element_all_ID[No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj]] == 1)
			{
				iX = info->ENC[No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj] * info->DIMENSION + 0];
				iY = info->ENC[No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj] * info->DIMENSION + 1];
				iZ = info->ENC[No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj] * info->DIMENSION + 2];

				for (ic = 0; ic < (info->Order[iPatch * info->DIMENSION + 0] + 1) * (info->Order[iPatch * info->DIMENSION + 1] + 1) * (info->Order[iPatch * info->DIMENSION + 2] + 1); ic++)
					iControlpoint[ic] = info->Controlpoint_of_Element[No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj] * MAX_NO_CP_ON_ELEMENT + ic];

				for (ig_i = 0; ig_i < GP_1D; ig_i++)
				{
					for (ig_j = 0; ig_j < GP_1D; ig_j++)
					{
						double Local_Coord[3], sfc, dxyzdgez_i[3] = {0.0}, dxyzdgez_j[3] = {0.0}, detJ, CoordParen_iDir, CoordParen_jDir, valDistLoad;
						int icc;
						Local_Coord[kCoord] = val_kCoord_Local;
						Local_Coord[iCoord] = Gxi_1D[ig_i];
						Local_Coord[jCoord] = Gxi_1D[ig_j];

						ShapeFunc_from_paren(Position_Data_param, Local_Coord, iCoord, No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj], info);
						CoordParen_iDir = Position_Data_param[iCoord];

						ShapeFunc_from_paren(Position_Data_param, Local_Coord, jCoord, No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj], info);
						CoordParen_jDir = Position_Data_param[jCoord];

						valDistLoad = (iCoeff_Dist_Load[0] + iCoeff_Dist_Load[1] * CoordParen_iDir + iCoeff_Dist_Load[2] * CoordParen_iDir * CoordParen_iDir)
									* (jCoeff_Dist_Load[0] + jCoeff_Dist_Load[1] * CoordParen_jDir + jCoeff_Dist_Load[2] * CoordParen_jDir * CoordParen_jDir);

						for (icc = 0; icc < (info->Order[iPatch * info->DIMENSION + 0] + 1) * (info->Order[iPatch * info->DIMENSION + 1] + 1) * (info->Order[iPatch * info->DIMENSION + 2] + 1); icc++)
						{
							dxyzdgez_i[0] += dShape_func(icc, iCoord, Local_Coord, No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj], info) * info->Node_Coordinate[iControlpoint[icc] * (info->DIMENSION + 1) + 0];
							dxyzdgez_i[1] += dShape_func(icc, iCoord, Local_Coord, No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj], info) * info->Node_Coordinate[iControlpoint[icc] * (info->DIMENSION + 1) + 1];
							dxyzdgez_i[2] += dShape_func(icc, iCoord, Local_Coord, No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj], info) * info->Node_Coordinate[iControlpoint[icc] * (info->DIMENSION + 1) + 2];

							dxyzdgez_j[0] += dShape_func(icc, jCoord, Local_Coord, No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj], info) * info->Node_Coordinate[iControlpoint[icc] * (info->DIMENSION + 1) + 0];
							dxyzdgez_j[1] += dShape_func(icc, jCoord, Local_Coord, No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj], info) * info->Node_Coordinate[iControlpoint[icc] * (info->DIMENSION + 1) + 1];
							dxyzdgez_j[2] += dShape_func(icc, jCoord, Local_Coord, No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj], info) * info->Node_Coordinate[iControlpoint[icc] * (info->DIMENSION + 1) + 2];
						}

						detJ = (dxyzdgez_i[0] * dxyzdgez_j[1] - dxyzdgez_i[1] * dxyzdgez_j[0]) + (dxyzdgez_i[1] * dxyzdgez_j[2] - dxyzdgez_i[2] * dxyzdgez_j[1]) + (dxyzdgez_i[2] * dxyzdgez_j[0] - dxyzdgez_i[0] * dxyzdgez_j[2]);

						if (type_load < 3)
						{
							for (ic = 0; ic < (info->Order[iPatch * info->DIMENSION + 0] + 1) * (info->Order[iPatch * info->DIMENSION + 1] + 1) * (info->Order[iPatch * info->DIMENSION + 2] + 1); ic++)
							{
								sfc = Shape_func(ic, Local_Coord, No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj], info);
								info->Equivalent_Nodal_Force[iControlpoint[ic] * info->DIMENSION + type_load] += valDistLoad * sfc * detJ * w_1D[ig_i] * w_1D[ig_j];
							}
						}
						else if (type_load == 3) // 法線方向
						{
							double LoadDir[3];
							LoadDir[0] = (dxyzdgez_i[1] * dxyzdgez_j[2] + dxyzdgez_i[2] * dxyzdgez_j[1]) / detJ;
							LoadDir[1] = (dxyzdgez_i[2] * dxyzdgez_j[0] + dxyzdgez_i[0] * dxyzdgez_j[2]) / detJ;
							LoadDir[2] = (dxyzdgez_i[0] * dxyzdgez_j[1] + dxyzdgez_i[1] * dxyzdgez_j[0]) / detJ;
							for (ic = 0; ic < (info->Order[iPatch * info->DIMENSION + 0] + 1) * (info->Order[iPatch * info->DIMENSION + 1] + 1) * (info->Order[iPatch * info->DIMENSION + 2] + 1); ic++)
							{
								sfc = Shape_func(ic, Local_Coord, No_Element_for_Integration[iii * info->Total_Knot_to_mesh[Total_mesh] + jjj], info);
								info->Equivalent_Nodal_Force[iControlpoint[ic] * info->DIMENSION + 0] +=   LoadDir[0] * valDistLoad * sfc * detJ * w_1D[ig_i] * w_1D[ig_j];
								info->Equivalent_Nodal_Force[iControlpoint[ic] * info->DIMENSION + 1] +=   LoadDir[1] * valDistLoad * sfc * detJ * w_1D[ig_i] * w_1D[ig_j];
								info->Equivalent_Nodal_Force[iControlpoint[ic] * info->DIMENSION + 2] += - LoadDir[2] * valDistLoad * sfc * detJ * w_1D[ig_i] * w_1D[ig_j];
							}
						}
					}
				}
			}
		}
	}
	free(No_Element_for_Integration), free(iControlpoint);
}


int SearchForElement_2D(int mesh_n, int iPatch, int iX, int iY, information *info)
{
	int iii;

	for (iii = 0; iii < info->Total_Element_on_mesh[mesh_n]; iii++)
	{
		if (info->Element_patch[iii + info->Total_Element_to_mesh[mesh_n]] == iPatch)
		{
			if (iX == info->ENC[(iii + info->Total_Element_to_mesh[mesh_n]) * info->DIMENSION + 0] && iY == info->ENC[(iii + info->Total_Element_to_mesh[mesh_n]) * info->DIMENSION + 1])
				goto loopend_2D;
		}
	}
	loopend_2D:

	return (iii);
}


int SearchForElement_3D(int mesh_n, int iPatch, int iX, int iY, int iZ, information *info)
{
	int iii;

	for (iii = 0; iii < info->Total_Element_on_mesh[mesh_n]; iii++)
	{
		if (info->Element_patch[iii + info->Total_Element_to_mesh[mesh_n]] == iPatch)
		{
			if (iX == info->ENC[(iii + info->Total_Element_to_mesh[mesh_n]) * info->DIMENSION + 0] && iY == info->ENC[(iii + info->Total_Element_to_mesh[mesh_n]) * info->DIMENSION + 1] && iZ == info->ENC[(iii + info->Total_Element_to_mesh[mesh_n]) * info->DIMENSION + 2])
				goto loopend_3D;
		}
	}
	loopend_3D:

	return (iii);
}


// for IGA
void Preprocessing_IGA(information *info)
{
	int e, global_mesh_num = 0, m = 0;

	Make_gauss_array(m, info);
	for (e = 0; e < info->real_Total_Element_on_mesh[global_mesh_num]; e++)
	{
		Preprocessing(m, e, info);
	}
}


// for S_IGA, coupled matrix を求める, 要素の重なりを要素のガウス点から求める
void Check_coupled_Glo_Loc_element(int mesh_n_over, int mesh_n_org, information *info)
{
	int re, e;
	int i, j, m;
	int l, ll;
	int MAX_NNLOVER = 0;

	int *temp_element_n = (int *)malloc(sizeof(int) * MAX_N_ELEMENT_OVER_POINT);
	int *element_n_point = (int *)malloc(sizeof(int) * MAX_N_ELEMENT_OVER_POINT * MAX_POW_NG_EXTEND);
	int *Check_coupled_No = (int *)malloc(sizeof(int) * MAX_N_ELEMENT_OVER);
	int *temp_ad = (int *)malloc(sizeof(int) * info->DIMENSION * (MAX_ORDER + 1));	// 要素の位置を求めるための値

	for (i = 0; i < MAX_N_ELEMENT_OVER; i++)
	{
		Check_coupled_No[i] = 0;
	}

	for (m = 0; m < 2; m++) // 最初 NG 個のガウス点で重なりを求め, info->NNLOVER[e] >= 2 の e に対して, 再度 NG_EXTEND 個のガウス点で重なりを求める
	{
		Make_gauss_array(m, info);

		// グローバルパッチの Preprocessing 作成
		if (m == 0)
		{
			for (re = 0; re < info->real_Total_Element_on_mesh[mesh_n_org]; re++)
			{
				e = info->real_element[re + info->real_Total_Element_to_mesh[mesh_n_org]];
				Preprocessing(m, e, info);
			}
		}

		// ローカルパッチ(mesh_n_over)各要素の頂点の物理座標算出
		for (re = 0; re < info->real_Total_Element_on_mesh[mesh_n_over]; re++)
		{
			int n;
			double output_para[MAX_DIMENSION];
			int Total_n_elements;

			e = info->real_element[re + info->real_Total_Element_to_mesh[mesh_n_over]];

			if (m == 0 || (m == 1 && info->NNLOVER[e] >= 2))
			{
				Preprocessing(m, e, info);

				if (m == 1)
				{
					info->NNLOVER[e] = 0;
					for (i = 0; i < info->NNLOVER[e]; i++)
					{
						info->NELOVER[e * MAX_N_ELEMENT_OVER + i] = 0;
					}
				}

				Total_n_elements = 0;
				ll = 0;

				// ローカルパッチ各要素のガウス点の物理座標のグローバルパッチでの(xi, eta)算出
				if (info->DIMENSION == 2)
				{
					for (i = 0; i < info->Total_Patch_on_mesh[mesh_n_org]; i++)
					{
						// グローバルパッチ i での各方向ノットベクトル
						double *temp_Position_Knots_xi = (double *)malloc(sizeof(double) * info->No_knot[i * info->DIMENSION + 0]);
						double *temp_Position_Knots_eta = (double *)malloc(sizeof(double) * info->No_knot[i * info->DIMENSION + 1]);

						for (j = 0; j < info->No_knot[i * info->DIMENSION + 0]; j++)
						{
							temp_Position_Knots_xi[j] = info->Position_Knots[info->Total_Knot_to_patch_dim[i * info->DIMENSION + 0] + j];
						}
						for (j = 0; j < info->No_knot[i * info->DIMENSION + 1]; j++)
						{
							temp_Position_Knots_eta[j] = info->Position_Knots[info->Total_Knot_to_patch_dim[i * info->DIMENSION + 1] + j];
						}

						for (n = 0; n < GP_ON_ELEMENT; n++)
						{
							double data_result_shape[MAX_DIMENSION] = {0.0};

							for (l = 0; l < info->DIMENSION; l++)
							{
								if (m == 0)
								{
									data_result_shape[l] = info->Gauss_Coordinate[e * GP_ON_ELEMENT * info->DIMENSION + n * info->DIMENSION + l];
								}
								else if (m == 1)
								{
									data_result_shape[l] = info->Gauss_Coordinate_ex[e * GP_ON_ELEMENT * info->DIMENSION + n * info->DIMENSION + l];
								}
							}
							int itr_n = Calc_xi_eta(data_result_shape[0], data_result_shape[1],
													temp_Position_Knots_xi, temp_Position_Knots_eta,
													info->No_Control_point[i * info->DIMENSION + 0], info->No_Control_point[i * info->DIMENSION + 1],
													info->Order[i * info->DIMENSION + 0], info->Order[i * info->DIMENSION + 1],
													&output_para[0], &output_para[1], info);
							
							if (m == 0)
							{
								info->Loc_parameter_on_Glo[e * GP_ON_ELEMENT * info->DIMENSION + n * info->DIMENSION + 0] = output_para[0];
								info->Loc_parameter_on_Glo[e * GP_ON_ELEMENT * info->DIMENSION + n * info->DIMENSION + 1] = output_para[1];
							}
							else if (m == 1)
							{
								info->Loc_parameter_on_Glo_ex[e * GP_ON_ELEMENT * info->DIMENSION + n * info->DIMENSION + 0] = output_para[0];
								info->Loc_parameter_on_Glo_ex[e * GP_ON_ELEMENT * info->DIMENSION + n * info->DIMENSION + 1] = output_para[1];
							}

							// Newton Laphsonによって出力されたxi, etaから重なる要素を求める
							int n_elements_over_point = ele_check(i, output_para, temp_element_n, temp_ad, info);
							if (itr_n == 0) // data_result_shapeがグローバルメッシュ上にないとき
							{
								n_elements_over_point = 0;
							}
							Total_n_elements += n_elements_over_point;
							for (l = 0; l < n_elements_over_point; l++)
							{
								element_n_point[ll] = temp_element_n[l];
								ll++;
							}
						}
						free(temp_Position_Knots_xi), free(temp_Position_Knots_eta);
					}
				}
				else if (info->DIMENSION == 3)
				{
					for (i = 0; i < info->Total_Patch_on_mesh[mesh_n_org]; i++)
					{
						// グローバルパッチ i での各方向ノットベクトル
						double *temp_Position_Knots_xi = (double *)malloc(sizeof(double) * info->No_knot[i * info->DIMENSION + 0]);
						double *temp_Position_Knots_eta = (double *)malloc(sizeof(double) * info->No_knot[i * info->DIMENSION + 1]);
						double *temp_Position_Knots_zeta = (double *)malloc(sizeof(double) * info->No_knot[i * info->DIMENSION + 2]);

						for (j = 0; j < info->No_knot[i * info->DIMENSION + 0]; j++)
						{
							temp_Position_Knots_xi[j] = info->Position_Knots[info->Total_Knot_to_patch_dim[i * info->DIMENSION + 0] + j];
						}
						for (j = 0; j < info->No_knot[i * info->DIMENSION + 1]; j++)
						{
							temp_Position_Knots_eta[j] = info->Position_Knots[info->Total_Knot_to_patch_dim[i * info->DIMENSION + 1] + j];
						}
						for (j = 0; j < info->No_knot[i * info->DIMENSION + 2]; j++)
						{
							temp_Position_Knots_zeta[j] = info->Position_Knots[info->Total_Knot_to_patch_dim[i * info->DIMENSION + 2] + j];
						}

						for (n = 0; n < GP_ON_ELEMENT; n++)
						{
							double data_result_shape[MAX_DIMENSION] = {0.0};

							for (l = 0; l < info->DIMENSION; l++)
							{
								if (m == 0)
								{
									data_result_shape[l] = info->Gauss_Coordinate[e * GP_ON_ELEMENT * info->DIMENSION + n * info->DIMENSION + l];
								}
								else if (m == 1)
								{
									data_result_shape[l] = info->Gauss_Coordinate_ex[e * GP_ON_ELEMENT * info->DIMENSION + n * info->DIMENSION + l];
								}
							}

							int itr_n = Calc_xi_eta_zeta(data_result_shape[0], data_result_shape[1], data_result_shape[2],
														 temp_Position_Knots_xi, temp_Position_Knots_eta, temp_Position_Knots_zeta,
														 info->No_Control_point[i * info->DIMENSION + 0], info->No_Control_point[i * info->DIMENSION + 1], info->No_Control_point[i * info->DIMENSION + 2],
														 info->Order[i * info->DIMENSION + 0], info->Order[i * info->DIMENSION + 1], info->Order[i * info->DIMENSION + 2],
														 &output_para[0], &output_para[1], &output_para[2], info);

							if (m == 0)
							{
								info->Loc_parameter_on_Glo[e * GP_ON_ELEMENT * info->DIMENSION + n * info->DIMENSION + 0] = output_para[0];
								info->Loc_parameter_on_Glo[e * GP_ON_ELEMENT * info->DIMENSION + n * info->DIMENSION + 1] = output_para[1];
								info->Loc_parameter_on_Glo[e * GP_ON_ELEMENT * info->DIMENSION + n * info->DIMENSION + 2] = output_para[2];
							}
							else if (m == 1)
							{
								info->Loc_parameter_on_Glo_ex[e * GP_ON_ELEMENT * info->DIMENSION + n * info->DIMENSION + 0] = output_para[0];
								info->Loc_parameter_on_Glo_ex[e * GP_ON_ELEMENT * info->DIMENSION + n * info->DIMENSION + 1] = output_para[1];
								info->Loc_parameter_on_Glo_ex[e * GP_ON_ELEMENT * info->DIMENSION + n * info->DIMENSION + 2] = output_para[2];
							}

							// Newton Laphsonによって出力されたxi, etaから重なる要素を求める
							int n_elements_over_point = ele_check(i, output_para, temp_element_n, temp_ad, info);
							if (itr_n == 0) // data_result_shapeがグローバルメッシュ上にないとき
							{
								n_elements_over_point = 0;
							}
							Total_n_elements += n_elements_over_point;
							for (l = 0; l < n_elements_over_point; l++)
							{
								element_n_point[ll] = temp_element_n[l];
								ll++;
							}
						}
						free(temp_Position_Knots_xi), free(temp_Position_Knots_eta), free(temp_Position_Knots_zeta);
					}
				}
				// 昇順ソート
				sort(Total_n_elements, element_n_point);
				// 重複削除
				info->NNLOVER[e] = duplicate_delete(Total_n_elements, e, element_n_point, info); // info->NNLOVER: 要素 e に重なる要素の総数
			}
		}
	}

	for (re = 0; re < info->real_Total_Element_on_mesh[mesh_n_over]; re++)
	{
		e = info->real_element[re + info->real_Total_Element_to_mesh[mesh_n_over]];

		Check_coupled_No[info->NNLOVER[e]]++;

		if (MAX_NNLOVER < info->NNLOVER[e])
		{
			MAX_NNLOVER = info->NNLOVER[e];
		}
	}

	// 重なっている要素の割合
	printf("MAX_NNLOVER = %d\n", MAX_NNLOVER);
	for (i = 0; i <= MAX_NNLOVER; i++)
	{
		double Percent_Check_coupled_No = (double)Check_coupled_No[i] * 100.0 / (double)info->real_Total_Element_on_mesh[mesh_n_over];
		printf("Check_coupled_No[%d] = %d\t\t%.2lf %%\n", i, Check_coupled_No[i], Percent_Check_coupled_No);
	}
	
	free(temp_element_n), free(element_n_point), free(Check_coupled_No), free(temp_ad);
}


void Make_Loc_Glo(information *info)
{
	int i, j, k;
	int jj;
	int e;
	int count;

	int j_n = info->real_Total_Element_to_mesh[Total_mesh] - info->real_Total_Element_on_mesh[0];

	for (i = 0; i < info->real_Total_Element_on_mesh[0]; i++)
	{
		e = info->real_element[i];
		count = 0;

		for (j = 0; j < j_n; j++)
		{
			jj = info->real_element[info->real_Total_Element_to_mesh[1] + j]; //ローカルメッシュ上のreal element番号

			if (info->NNLOVER[jj] > 0)
			{
				for (k = 0; k < info->NNLOVER[jj]; k++)
				{
					if (info->NELOVER[jj * MAX_N_ELEMENT_OVER + k] == e)
					{
						info->NELOVER[e * MAX_N_ELEMENT_OVER + count] = jj;
						count++;
					}
				}
			}
		}
		info->NNLOVER[e] = count;
	}
}


// Newton Raphsonによって出力されたxi,etaから重なる要素を求める
int ele_check(int patch_n, double *para_coord, int *temp_element_n, int *temp_ad, information *info)
{
	int i, j, k, l, kk;
	int RangeCheck_flag;					// 要素を求め終えたら立てるフラグ
	int No_line[MAX_DIMENSION];				// xi, etaが含まれている要素の列数
	int n = 1;

	for (j = 0; j < info->DIMENSION; j++)
	{
		// 初期化
		l = 0;
		No_line[j] = 0;
		for (i = 0; i < MAX_ORDER + 1; i++)
		{
			temp_ad[j * (MAX_ORDER + 1) + i] = 0;
		}
		RangeCheck_flag = 0;

		for (k = 0; k < info->No_Control_point[patch_n * info->DIMENSION + j] - 1; k++)
		{
			if (RangeCheck_flag == 1)
				break;

			// Local要素の頂点がGlobalパッチ内にない場合
			if (para_coord[j] < info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + j] + 0] || para_coord[j] > info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + j] + info->No_knot[patch_n * info->DIMENSION + j] - 1])
			{
				RangeCheck_flag++;
			}

			// Local要素の頂点がGlobal要素内部にある場合
			if (para_coord[j] < info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + j] + info->Order[patch_n * info->DIMENSION + j] + k])
			{
				int kk = 0;
				for (kk = 0; kk < k + 1; kk++)
				{
					if (para_coord[j] > info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + j] + info->Order[patch_n * info->DIMENSION + j] + k - kk])
					{
						temp_ad[j * (MAX_ORDER + 1) + l] = k - kk;
						l++;
						RangeCheck_flag++;
						break;
					}
				}
			}

			// Local要素の頂点がGlobal要素境界上にある場合
			if (para_coord[j] == info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + j] + info->Order[patch_n * info->DIMENSION + j] + k])
			{
				//頂点の座標がGlobalパッチの始点上にある場合
				if (para_coord[j] == info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + j] + 0])
				{
					temp_ad[j * (MAX_ORDER + 1) + l] = k;
					l++;
					break;
				}
				//頂点の座標がGlobalパッチの終点上にある場合
				if (para_coord[j] == info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + j] + info->No_knot[patch_n * info->DIMENSION + j] - 1])
				{
					temp_ad[j * (MAX_ORDER + 1) + l] = k - 1;
					l++;
					break;
				}
				//頂点の座標がGlobal要素境界上にある場合
				else
				{
					temp_ad[j * (MAX_ORDER + 1) + l] = k - 1;
					l++;
					temp_ad[j * (MAX_ORDER + 1) + l] = k;
					l++;
				}
				for (kk = 0; kk < info->Order[patch_n * info->DIMENSION + j]; kk++)
				{
					if (para_coord[j] == info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + j] + info->Order[patch_n * info->DIMENSION + j] + k + kk + 1])
					// 多重ノット(次数分ループ)
					{
						printf("C0 continuity\n");
						temp_ad[j * (MAX_ORDER + 1) + l] = k + kk;
						l++;
					}
					if (para_coord[j] != info->Position_Knots[info->Total_Knot_to_patch_dim[patch_n * info->DIMENSION + j] + info->Order[patch_n * info->DIMENSION + j] + k + kk + 1])
						break;
				}
				RangeCheck_flag++;
			}
		}
		No_line[j] = l;
		
		// 各方向のNo_lineを掛け合わせる
		n *= l;
	}

	if (info->DIMENSION == 2)
	{
		for (i = 0; i < No_line[1]; i++)
		{
			for (j = 0; j < No_line[0]; j++)
			{
				temp_element_n[i * No_line[0] + j] = temp_ad[0 * (MAX_ORDER + 1) + j] + temp_ad[1 * (MAX_ORDER + 1) + i] * info->line_No_Total_element[patch_n * info->DIMENSION + 0];
			}
		}
	}
	else if (info->DIMENSION == 3)
	{
		for (i = 0; i < No_line[2]; i++)
		{
			for (j = 0; j < No_line[1]; j++)
			{
				for (k = 0; k < No_line[0]; k++)
				{
					temp_element_n[i * No_line[0] * No_line[1] + j * No_line[0] + k] = temp_ad[0 * (MAX_ORDER + 1) + k] + temp_ad[1 * (MAX_ORDER + 1) + j] * info->line_No_Total_element[patch_n * info->DIMENSION + 0] + temp_ad[2 * (MAX_ORDER + 1) + i] * info->line_No_Total_element[patch_n * info->DIMENSION + 0] * info->line_No_Total_element[patch_n * info->DIMENSION + 1];
				}
			}
		}
	}

	return n;
}


// 昇順ソート
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


// 重複削除
int duplicate_delete(int total, int element_n, int *element_n_point, information *info)
{
	int i, j;

	j = 0;
	info->NELOVER[element_n * MAX_N_ELEMENT_OVER + j] = element_n_point[0];
	j++;
	for (i = 1; i < total; i++)
	{
		if (element_n_point[i] != element_n_point[i - 1])
		{
			info->NELOVER[element_n * MAX_N_ELEMENT_OVER + j] = element_n_point[i];
			j++;
		}
	}
	// j = 要素element_nに重なる要素の総数
	return j;
}


// Preprocessing
void Preprocessing(int m, int e, information *info)
{
	// ガウス点の物理座標を計算
	Make_Gauss_Coordinate(m, e, info);

	// ガウス点でのヤコビアン, Bマトリックスを計算
	Make_dSF(m, e, info);
	Make_Jac(m, e, info);
	Make_B_Matrix(m, e, info);
}


void Make_Gauss_Coordinate(int m, int e, information *info)
{
	int i, j, k;
	double temp_coordinate[MAX_DIMENSION];
	double R;

	for (i = 0; i < GP_ON_ELEMENT; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			temp_coordinate[j] = Gxi[i][j];
		}
		
		for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; j++)
		{
			R = Shape_func(j, temp_coordinate, e, info);

			for (k = 0; k < info->DIMENSION; k++)
			{
				if (m == 0)
				{
					info->Gauss_Coordinate[e * GP_ON_ELEMENT * info->DIMENSION + i * info->DIMENSION + k] += R * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + j] * (info->DIMENSION + 1) + k];
				}
				else if (m == 1)
				{
					info->Gauss_Coordinate_ex[e * GP_ON_ELEMENT * info->DIMENSION + i * info->DIMENSION + k] += R * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + j] * (info->DIMENSION + 1) + k];
				}
			}
		}
	}
}


void Make_dSF(int m, int e, information *info)
{
	int i, j, k;
	double temp_coordinate[MAX_DIMENSION];

	for (i = 0; i < GP_ON_ELEMENT; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			temp_coordinate[j] = Gxi[i][j];
		}

		for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; j++)
		{
			for (k = 0; k < info->DIMENSION; k++)
			{
				if (m == 0)
				{
					info->dSF[i * MAX_NO_CP_ON_ELEMENT * info->DIMENSION + j * info->DIMENSION + k] = dShape_func(j, k, temp_coordinate, e, info);
				}
				else if (m == 1)
				{
					info->dSF_ex[i * MAX_NO_CP_ON_ELEMENT * info->DIMENSION + j * info->DIMENSION + k] = dShape_func(j, k, temp_coordinate, e, info);
				}
			}
		}
	}
}


void Make_Jac(int m, int e, information *info)
{
	int i, j, k, l;
	double J = 0.0;
	double a_2x2[2][2], a_3x3[3][3];

	for (i = 0; i < GP_ON_ELEMENT; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			for (k = 0; k < info->DIMENSION; k++)
			{
				if (info->DIMENSION == 2)
				{
					a_2x2[j][k] = 0.0;
					for (l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; l++)
					{
						if (m == 0)
						{
							a_2x2[j][k] += info->dSF[i * MAX_NO_CP_ON_ELEMENT * info->DIMENSION + l * info->DIMENSION + k] * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + l] * (info->DIMENSION + 1) + j];
						}
						else if (m == 1)
						{
							a_2x2[j][k] += info->dSF_ex[i * MAX_NO_CP_ON_ELEMENT * info->DIMENSION + l * info->DIMENSION + k] * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + l] * (info->DIMENSION + 1) + j];
						}
					}
				}
				else if (info->DIMENSION == 3)
				{
					a_3x3[j][k] = 0.0;
					for (l = 0; l < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; l++)
					{
						if (m == 0)
						{
							a_3x3[j][k] += info->dSF[i * MAX_NO_CP_ON_ELEMENT * info->DIMENSION + l * info->DIMENSION + k] * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + l] * (info->DIMENSION + 1) + j];
						}
						else if (m == 1)
						{
							a_3x3[j][k] += info->dSF_ex[i * MAX_NO_CP_ON_ELEMENT * info->DIMENSION + l * info->DIMENSION + k] * info->Node_Coordinate[info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + l] * (info->DIMENSION + 1) + j];
						}
					}
				}
			}
		}

		if (info->DIMENSION == 2)
		{
			J = InverseMatrix_2x2(a_2x2);

			for (j = 0; j < info->DIMENSION; j++)
			{
				for (k = 0; k < info->DIMENSION; k++)
				{
					info->a_matrix[i * info->DIMENSION * info->DIMENSION + j * info->DIMENSION + k] = a_2x2[j][k];
				}
			}
		}
		else if (info->DIMENSION == 3)
		{
			J = InverseMatrix_3x3(a_3x3);

			for (j = 0; j < info->DIMENSION; j++)
			{
				for (k = 0; k < info->DIMENSION; k++)
				{
					info->a_matrix[i * info->DIMENSION * info->DIMENSION + j * info->DIMENSION + k] = a_3x3[j][k];
				}
			}
		}

		if (m == 0)
		{
			info->Jac[e * GP_ON_ELEMENT + i] = J;
		}
		else if (m == 1)
		{
			info->Jac_ex[e * GP_ON_ELEMENT + i] = J;
		}

		if (J <= 0)
		{
			double jac;
			if (J == ERROR)
				jac = 0.0;
			else
				jac = J;

			printf("Error, J = %le\nJ must be positive value.", jac);
			exit(1);
		}
	}
}


void Make_B_Matrix(int m, int e, information *info)
{
	int i, j, k, l;
	double *b = (double *)malloc(sizeof(double) * info->DIMENSION * MAX_NO_CP_ON_ELEMENT);	// b[info->DIMENSION][MAX_NO_CP_ON_ELEMENT]

	for (i = 0; i < GP_ON_ELEMENT; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			for (k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; k++)
			{
				b[j * MAX_NO_CP_ON_ELEMENT + k] = 0.0;
				for (l = 0; l < info->DIMENSION; l++)
				{
					if (m == 0)
					{
						b[j * MAX_NO_CP_ON_ELEMENT + k] += info->a_matrix[i * info->DIMENSION * info->DIMENSION + l * info->DIMENSION + j] * info->dSF[i * MAX_NO_CP_ON_ELEMENT * info->DIMENSION + k * info->DIMENSION + l];
					}
					else if (m == 1)
					{
						b[j * MAX_NO_CP_ON_ELEMENT + k] += info->a_matrix[i * info->DIMENSION * info->DIMENSION + l * info->DIMENSION + j] * info->dSF_ex[i * MAX_NO_CP_ON_ELEMENT * info->DIMENSION + k * info->DIMENSION + l];
					}
				}
			}
		}

		if (info->DIMENSION == 2)
		{
			if (m == 0)
			{
				for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; j++)
				{
					info->B_Matrix[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 0 * MAX_KIEL_SIZE + (2 * j)]     = b[0 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 0 * MAX_KIEL_SIZE + (2 * j + 1)] = 0.0;
					
					info->B_Matrix[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 1 * MAX_KIEL_SIZE + (2 * j)]     = 0.0;
					info->B_Matrix[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 1 * MAX_KIEL_SIZE + (2 * j + 1)] = b[1 * MAX_NO_CP_ON_ELEMENT + j];

					info->B_Matrix[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 2 * MAX_KIEL_SIZE + (2 * j)]     = b[1 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 2 * MAX_KIEL_SIZE + (2 * j + 1)] = b[0 * MAX_NO_CP_ON_ELEMENT + j];
				}
			}
			else if (m == 1)
			{
				for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; j++)
				{
					info->B_Matrix_ex[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 0 * MAX_KIEL_SIZE + (2 * j)]     = b[0 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix_ex[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 0 * MAX_KIEL_SIZE + (2 * j + 1)] = 0.0;

					info->B_Matrix_ex[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 1 * MAX_KIEL_SIZE + (2 * j)]     = 0.0;
					info->B_Matrix_ex[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 1 * MAX_KIEL_SIZE + (2 * j + 1)] = b[1 * MAX_NO_CP_ON_ELEMENT + j];
					
					info->B_Matrix_ex[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 2 * MAX_KIEL_SIZE + (2 * j)]     = b[1 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix_ex[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 2 * MAX_KIEL_SIZE + (2 * j + 1)] = b[0 * MAX_NO_CP_ON_ELEMENT + j];
				}
			}
		}
		else if (info->DIMENSION == 3)
		{
			if (m == 0)
			{
				for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; j++)
				{
					info->B_Matrix[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 0 * MAX_KIEL_SIZE + (3 * j)]     = b[0 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 0 * MAX_KIEL_SIZE + (3 * j + 1)] = 0.0;
					info->B_Matrix[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 0 * MAX_KIEL_SIZE + (3 * j + 2)] = 0.0;

					info->B_Matrix[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 1 * MAX_KIEL_SIZE + (3 * j)]     = 0.0;
					info->B_Matrix[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 1 * MAX_KIEL_SIZE + (3 * j + 1)] = b[1 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 1 * MAX_KIEL_SIZE + (3 * j + 2)] = 0.0;

					info->B_Matrix[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 2 * MAX_KIEL_SIZE + (3 * j)]     = 0.0;
					info->B_Matrix[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 2 * MAX_KIEL_SIZE + (3 * j + 1)] = 0.0;
					info->B_Matrix[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 2 * MAX_KIEL_SIZE + (3 * j + 2)] = b[2 * MAX_NO_CP_ON_ELEMENT + j];

					info->B_Matrix[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 3 * MAX_KIEL_SIZE + (3 * j)]     = b[1 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 3 * MAX_KIEL_SIZE + (3 * j + 1)] = b[0 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 3 * MAX_KIEL_SIZE + (3 * j + 2)] = 0.0;

					info->B_Matrix[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 4 * MAX_KIEL_SIZE + (3 * j)]     = 0.0;
					info->B_Matrix[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 4 * MAX_KIEL_SIZE + (3 * j + 1)] = b[2 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 4 * MAX_KIEL_SIZE + (3 * j + 2)] = b[1 * MAX_NO_CP_ON_ELEMENT + j];

					info->B_Matrix[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 5 * MAX_KIEL_SIZE + (3 * j)]     = b[2 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 5 * MAX_KIEL_SIZE + (3 * j + 1)] = 0.0;
					info->B_Matrix[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 5 * MAX_KIEL_SIZE + (3 * j + 2)] = b[0 * MAX_NO_CP_ON_ELEMENT + j];
				}
			}
			else if (m == 1)
			{
				for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[e]]; j++)
				{
					info->B_Matrix_ex[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 0 * MAX_KIEL_SIZE + (3 * j)]     = b[0 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix_ex[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 0 * MAX_KIEL_SIZE + (3 * j + 1)] = 0.0;
					info->B_Matrix_ex[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 0 * MAX_KIEL_SIZE + (3 * j + 2)] = 0.0;

					info->B_Matrix_ex[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 1 * MAX_KIEL_SIZE + (3 * j)]     = 0.0;
					info->B_Matrix_ex[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 1 * MAX_KIEL_SIZE + (3 * j + 1)] = b[1 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix_ex[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 1 * MAX_KIEL_SIZE + (3 * j + 2)] = 0.0;

					info->B_Matrix_ex[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 2 * MAX_KIEL_SIZE + (3 * j)]     = 0.0;
					info->B_Matrix_ex[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 2 * MAX_KIEL_SIZE + (3 * j + 1)] = 0.0;
					info->B_Matrix_ex[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 2 * MAX_KIEL_SIZE + (3 * j + 2)] = b[2 * MAX_NO_CP_ON_ELEMENT + j];

					info->B_Matrix_ex[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 3 * MAX_KIEL_SIZE + (3 * j)]     = b[1 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix_ex[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 3 * MAX_KIEL_SIZE + (3 * j + 1)] = b[0 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix_ex[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 3 * MAX_KIEL_SIZE + (3 * j + 2)] = 0.0;

					info->B_Matrix_ex[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 4 * MAX_KIEL_SIZE + (3 * j)]     = 0.0;
					info->B_Matrix_ex[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 4 * MAX_KIEL_SIZE + (3 * j + 1)] = b[2 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix_ex[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 4 * MAX_KIEL_SIZE + (3 * j + 2)] = b[1 * MAX_NO_CP_ON_ELEMENT + j];

					info->B_Matrix_ex[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 5 * MAX_KIEL_SIZE + (3 * j)]     = b[2 * MAX_NO_CP_ON_ELEMENT + j];
					info->B_Matrix_ex[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 5 * MAX_KIEL_SIZE + (3 * j + 1)] = 0.0;
					info->B_Matrix_ex[e * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + 5 * MAX_KIEL_SIZE + (3 * j + 2)] = b[0 * MAX_NO_CP_ON_ELEMENT + j];
				}
			}
		}
	}

	free(b);
}


void Make_gauss_array(int select_GP, information *info)
{
	int i, j, k;

	if (select_GP == 0)
	{
		GP_1D = NG;
	}
	else if (select_GP == 1)
	{
		GP_1D = NG_EXTEND;
	}

	GP_2D = GP_1D * GP_1D;
	GP_3D = GP_2D * GP_1D;

	if (info->DIMENSION == 2)
	{
		GP_ON_ELEMENT = GP_2D;
	}
	else if (info->DIMENSION == 3)
	{
		GP_ON_ELEMENT = GP_3D;
	}

	if (info->DIMENSION == 2)
	{
		if (GP_1D == 3)
		{
			double G1 = pow((3.0 / 5.0), 0.5);
			double G_vec[3] = {-G1, 0.0, G1};
			double w1 = 8.0 / 9.0;
			double w2 = 5.0 / 9.0;
			double w_vec[3] = {w2, w1, w2};

			for (i = 0; i < GP_1D; i++)
			{
				for (j = 0; j < GP_1D; j++)
				{
					w[GP_1D * i + j] = w_vec[i] * w_vec[j];
					Gxi[GP_1D * i + j][0] = G_vec[j];
					Gxi[GP_1D * i + j][1] = G_vec[i];
				}
				w_1D[i] = w_vec[i];
				Gxi_1D[i] = G_vec[i];
			}
		}
		else if (GP_1D == 4)
		{
			double A = pow((6.0 / 5.0), 0.5);
			double G1 = pow(((3.0 - 2.0 * A) / 7.0), 0.5);
			double G2 = pow(((3.0 + 2.0 * A) / 7.0), 0.5);
			double G_vec[4] = {-G2, -G1, G1, G2};
			double B = pow(30.0, 0.5);
			double w1 = (18.0 + B) / 36.0;
			double w2 = (18.0 - B) / 36.0;
			double w_vec[4] = {w2, w1, w1, w2};

			for (i = 0; i < GP_1D; i++)
			{
				for (j = 0; j < GP_1D; j++)
				{
					w[GP_1D * i + j] = w_vec[i] * w_vec[j];
					Gxi[GP_1D * i + j][0] = G_vec[j];
					Gxi[GP_1D * i + j][1] = G_vec[i];
				}
				w_1D[i] = w_vec[i];
				Gxi_1D[i] = G_vec[i];
			}
		}
		else if (GP_1D == 5)
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

			for (i = 0; i < GP_1D; i++)
			{
				for (j = 0; j < GP_1D; j++)
				{
					w[GP_1D * i + j] = w_vec[i] * w_vec[j];
					Gxi[GP_1D * i + j][0] = G_vec[j];
					Gxi[GP_1D * i + j][1] = G_vec[i];
				}
				w_1D[i] = w_vec[i];
				Gxi_1D[i] = G_vec[i];
			}
		}
		else if (GP_1D == 10)
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

			for (i = 0; i < GP_1D; i++)
			{
				for (j = 0; j < GP_1D; j++)
				{
					w[GP_1D * i + j] = w_vec[i] * w_vec[j];
					Gxi[GP_1D * i + j][0] = G_vec[j];
					Gxi[GP_1D * i + j][1] = G_vec[i];
				}
				w_1D[i] = w_vec[i];
				Gxi_1D[i] = G_vec[i];
			}
		}
	}
	else if (info->DIMENSION == 3)
	{
		if (GP_1D == 3)
		{
			double G1 = pow((3.0 / 5.0), 0.5);
			double G_vec[3] = {-G1, 0.0, G1};
			double w1 = 8.0 / 9.0;
			double w2 = 5.0 / 9.0;
			double w_vec[3] = {w2, w1, w2};

			for (i = 0; i < GP_1D; i++)
			{
				for (j = 0; j < GP_1D; j++)
				{
					for (k = 0; k < GP_1D; k++)
					{
						w[i * GP_2D + j * GP_1D + k] = w_vec[i] * w_vec[j] * w_vec[k];
						Gxi[i * GP_2D + j * GP_1D + k][0] = G_vec[k];
						Gxi[i * GP_2D + j * GP_1D + k][1] = G_vec[j];
						Gxi[i * GP_2D + j * GP_1D + k][2] = G_vec[i];
					}
				}
				w_1D[i] = w_vec[i];
				Gxi_1D[i] = G_vec[i];
			}
		}
		else if (GP_1D == 4)
		{
			double A = pow((6.0 / 5.0), 0.5);
			double G1 = pow(((3.0 - 2.0 * A) / 7.0), 0.5);
			double G2 = pow(((3.0 + 2.0 * A) / 7.0), 0.5);
			double G_vec[4] = {-G2, -G1, G1, G2};
			double B = pow(30.0, 0.5);
			double w1 = (18.0 + B) / 36.0;
			double w2 = (18.0 - B) / 36.0;
			double w_vec[4] = {w2, w1, w1, w2};

			for (i = 0; i < GP_1D; i++)
			{
				for (j = 0; j < GP_1D; j++)
				{
					for (k = 0; k < GP_1D; k++)
					{
						w[i * GP_2D + j * GP_1D + k] = w_vec[i] * w_vec[j] * w_vec[k];
						Gxi[i * GP_2D + j * GP_1D + k][0] = G_vec[k];
						Gxi[i * GP_2D + j * GP_1D + k][1] = G_vec[j];
						Gxi[i * GP_2D + j * GP_1D + k][2] = G_vec[i];
					}
				}
				w_1D[i] = w_vec[i];
				Gxi_1D[i] = G_vec[i];
			}
		}
		else if (GP_1D == 5)
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

			for (i = 0; i < GP_1D; i++)
			{
				for (j = 0; j < GP_1D; j++)
				{
					for (k = 0; k < GP_1D; k++)
					{
						w[i * GP_2D + j * GP_1D + k] = w_vec[i] * w_vec[j] * w_vec[k];
						Gxi[i * GP_2D + j * GP_1D + k][0] = G_vec[k];
						Gxi[i * GP_2D + j * GP_1D + k][1] = G_vec[j];
						Gxi[i * GP_2D + j * GP_1D + k][2] = G_vec[i];
					}
				}
				w_1D[i] = w_vec[i];
				Gxi_1D[i] = G_vec[i];
			}
		}
		else if (GP_1D == 10)
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

			for (i = 0; i < GP_1D; i++)
			{
				for (j = 0; j < GP_1D; j++)
				{
					for (k = 0; k < GP_1D; k++)
					{
						w[i * GP_2D + j * GP_1D + k] = w_vec[i] * w_vec[j] * w_vec[k];
						Gxi[i * GP_2D + j * GP_1D + k][0] = G_vec[k];
						Gxi[i * GP_2D + j * GP_1D + k][1] = G_vec[j];
						Gxi[i * GP_2D + j * GP_1D + k][2] = G_vec[i];
					}
				}
				w_1D[i] = w_vec[i];
				Gxi_1D[i] = G_vec[i];
			}
		}
	}
}


// K matrix
void Make_D_Matrix(information *info)
{
	int i, j;

	if (info->DIMENSION == 2)
	{
		if (DM == 0) // 平面応力状態
		{
			double Eone = E / (1.0 - nu * nu);
			double D1[3][3] = {{Eone, nu * Eone, 0}, {nu * Eone, Eone, 0}, {0, 0, (1 - nu) * Eone / 2.0}};

			for (i = 0; i < D_MATRIX_SIZE; i++)
				for (j = 0; j < D_MATRIX_SIZE; j++)
					info->D[i * D_MATRIX_SIZE + j] = D1[i][j];
		}
		else if (DM == 1) // 平面ひずみ状態
		{
			double Eone = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
			double D1[3][3] = {{Eone, nu / (1.0 - nu) * Eone, 0}, {nu / (1.0 - nu) * Eone, Eone, 0}, {0, 0, (1 - 2 * nu) / 2.0 / (1.0 - nu) * Eone}};

			for (i = 0; i < D_MATRIX_SIZE; i++)
				for (j = 0; j < D_MATRIX_SIZE; j++)
					info->D[i * D_MATRIX_SIZE + j] = D1[i][j];
		}
	}
	else if (info->DIMENSION == 3)
	{
		double E_ii = (1.0 -  nu) / ((1.0 + nu) * (1.0 - 2.0 * nu)) * E;
		double E_ij = nu / ((1.0 + nu) * (1.0 - 2.0 * nu)) * E;
		double G = E / (2.0 * (1.0 + nu));

		double D1[6][6] = {{E_ii, E_ij, E_ij, 0.0, 0.0, 0.0},
						   {E_ij, E_ii, E_ij, 0.0, 0.0, 0.0},
						   {E_ij, E_ij, E_ii, 0.0, 0.0, 0.0},
						   { 0.0, 0.0,  0.0,  G  , 0.0, 0.0},
						   { 0.0, 0.0,  0.0,  0.0, G  , 0.0},
						   { 0.0, 0.0,  0.0,  0.0, 0.0, G  }};
		
		for (i = 0; i < D_MATRIX_SIZE; i++)
				for (j = 0; j < D_MATRIX_SIZE; j++)
					info->D[i * D_MATRIX_SIZE + j] = D1[i][j];
	}
}


// 拘束されている行数を省いた行列の番号の制作
void Make_Index_Dof(information *info)
{
	int i, k = 0;

	// 拘束されている自由度(Degree Of free)をERRORにする
	for (i = 0; i < info->Total_Constraint_to_mesh[Total_mesh]; i++)
	{
		info->Index_Dof[info->Constraint_Node_Dir[i * 2 + 0] * info->DIMENSION + info->Constraint_Node_Dir[i * 2 + 1]] = ERROR;
	}
	// ERROR以外に番号を付ける
	for (i = 0; i < info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION; i++)
	{
		if (info->Index_Dof[i] != ERROR)
		{
			info->Index_Dof[i] = k;
			k++;
		}
	}

	K_Whole_Size = k;
	printf("\nK_Whole_Size = %d\n\n", k);
}


void Make_K_Whole_Ptr_Col(information *info, int mode_select)
{
	int i, j, k, ii, jj;
	int N, NE, i_index, j_index;

	// mode_select == 0: Ptr を作成
	// mode_select == 1: Col を作成

	if (mode_select == 0)
	{
		// 初期化
		for (i = 0; i < K_Whole_Size + 1; i++)
			info->K_Whole_Ptr[i] = 0;

		// 大きく分割するためのループ
		for (N = 0; N < info->Total_Control_Point_to_mesh[Total_mesh]; N += K_DIVISION_LENGE)
		{
			// 各節点に接する節点を取得
			for (i = 0; i < K_DIVISION_LENGE; i++)
			{
				info->Total_Control_Point_To_Node[i] = 0;
			}
			for (i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
			{
				for (ii = 0; ii < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; ii++)
				{
					NE = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + ii] - N;
					if (0 <= NE && NE < K_DIVISION_LENGE)
					{
						// ローカル要素
						for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; j++)
						{
							// 数字がない時
							if (info->Total_Control_Point_To_Node[NE] == 0)
							{
								// 節点番号を取得
								info->Node_To_Node[NE * info->Total_Control_Point_to_mesh[Total_mesh] + 0] = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j];
								info->Total_Control_Point_To_Node[NE]++;
							}
							// 同じものがあったら, k > 0 以降の取得, kのカウント
							for (k = 0; k < info->Total_Control_Point_To_Node[NE]; k++)
							{
								if (info->Node_To_Node[NE * info->Total_Control_Point_to_mesh[Total_mesh] + k] == info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j])
								{
									break;
								}
							}
							// 未設定のNode_To_Node取得
							if (k == info->Total_Control_Point_To_Node[NE])
							{
								info->Node_To_Node[NE * info->Total_Control_Point_to_mesh[Total_mesh] + k] = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j];
								info->Total_Control_Point_To_Node[NE]++;
							}
						}
						// 別メッシュとの重なりを考慮
						if (info->NNLOVER[i] > 0)
						{
							for (jj = 0; jj < info->NNLOVER[i]; jj++)
							{
								// ローカル要素
								for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[info->NELOVER[i * MAX_N_ELEMENT_OVER + jj]]]; j++)
								{
									// 数字がない時
									if (info->Total_Control_Point_To_Node[NE] == 0)
									{
										// 節点番号を取得
										info->Node_To_Node[NE * info->Total_Control_Point_to_mesh[Total_mesh] + 0] = info->Controlpoint_of_Element[info->NELOVER[i * MAX_N_ELEMENT_OVER + jj] * MAX_NO_CP_ON_ELEMENT + j];
										info->Total_Control_Point_To_Node[NE]++;
									}

									// 同じものがあったら, k > 0 以降の取得, kのカウント
									for (k = 0; k < info->Total_Control_Point_To_Node[NE]; k++)
									{
										if (info->Node_To_Node[NE * info->Total_Control_Point_to_mesh[Total_mesh] + k] == info->Controlpoint_of_Element[info->NELOVER[i * MAX_N_ELEMENT_OVER + jj] * MAX_NO_CP_ON_ELEMENT + j])
										{
											break;
										}
									}
									// 未設定のNode_To_Node取得
									if (k == info->Total_Control_Point_To_Node[NE])
									{
										info->Node_To_Node[NE * info->Total_Control_Point_to_mesh[Total_mesh] + k] = info->Controlpoint_of_Element[info->NELOVER[i * MAX_N_ELEMENT_OVER + jj] * MAX_NO_CP_ON_ELEMENT + j];
										info->Total_Control_Point_To_Node[NE]++;
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
				if (N + i < info->Total_Control_Point_to_mesh[Total_mesh])
				{
					for (j = 0; j < info->Total_Control_Point_To_Node[i]; j++)
					{
						int Min = info->Node_To_Node[i * info->Total_Control_Point_to_mesh[Total_mesh] + j], No = j;
						for (k = j; k < info->Total_Control_Point_To_Node[i]; k++)
						{
							if (Min > info->Node_To_Node[i * info->Total_Control_Point_to_mesh[Total_mesh] + k])
							{
								Min = info->Node_To_Node[i * info->Total_Control_Point_to_mesh[Total_mesh] + k];
								No = k;
							}
						}
						for (k = No; k > j; k--)
						{
							info->Node_To_Node[i * info->Total_Control_Point_to_mesh[Total_mesh] + k] = info->Node_To_Node[i * info->Total_Control_Point_to_mesh[Total_mesh] + k - 1];
						}
						info->Node_To_Node[i * info->Total_Control_Point_to_mesh[Total_mesh] + j] = Min;
					}
				}
			}

			// 節点からcol ptrを求める
			ii = 0;
			k = 0;
			for (i = 0; i < K_DIVISION_LENGE; i++)
			{
				for (ii = 0; ii < info->DIMENSION; ii++)
				{
					if (N + i < info->Total_Control_Point_to_mesh[Total_mesh])
					{
						i_index = info->Index_Dof[(N + i) * info->DIMENSION + ii];
						k = 0;
						if (i_index >= 0)
						{
							info->K_Whole_Ptr[i_index + 1] = info->K_Whole_Ptr[i_index];
							for (j = 0; j < info->Total_Control_Point_To_Node[i]; j++)
							{
								for (jj = 0; jj < info->DIMENSION; jj++)
								{
									j_index = info->Index_Dof[info->Node_To_Node[i * info->Total_Control_Point_to_mesh[Total_mesh] + j] * info->DIMENSION + jj];
									if (j_index >= 0 && j_index >= i_index)
									{
										info->K_Whole_Ptr[i_index + 1]++;
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
	else if (mode_select == 1)
	{
		// 大きく分割するためのループ
		for (N = 0; N < info->Total_Control_Point_to_mesh[Total_mesh]; N += K_DIVISION_LENGE)
		{
			// 各節点に接する節点を取得
			for (i = 0; i < K_DIVISION_LENGE; i++)
			{
				info->Total_Control_Point_To_Node[i] = 0;
			}
			for (i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
			{
				for (ii = 0; ii < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; ii++)
				{
					NE = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + ii] - N;
					if (0 <= NE && NE < K_DIVISION_LENGE)
					{
						// ローカル要素
						for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; j++)
						{
							// 数字がない時
							if (info->Total_Control_Point_To_Node[NE] == 0)
							{
								// 節点番号を取得
								info->Node_To_Node[NE * info->Total_Control_Point_to_mesh[Total_mesh] + 0] = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j];
								info->Total_Control_Point_To_Node[NE]++;
							}
							// 同じものがあったら, k > 0 以降の取得, kのカウント
							for (k = 0; k < info->Total_Control_Point_To_Node[NE]; k++)
							{
								if (info->Node_To_Node[NE * info->Total_Control_Point_to_mesh[Total_mesh] + k] == info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j])
								{
									break;
								}
							}
							// 未設定のNode_To_Node取得
							if (k == info->Total_Control_Point_To_Node[NE])
							{
								info->Node_To_Node[NE * info->Total_Control_Point_to_mesh[Total_mesh] + k] = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j];
								info->Total_Control_Point_To_Node[NE]++;
							}
						}
						// 別メッシュとの重なりを考慮
						if (info->NNLOVER[i] > 0)
						{
							for (jj = 0; jj < info->NNLOVER[i]; jj++)
							{
								// ローカル要素
								for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[info->NELOVER[i * MAX_N_ELEMENT_OVER + jj]]]; j++)
								{
									// 数字がない時
									if (info->Total_Control_Point_To_Node[NE] == 0)
									{
										// 節点番号を取得
										info->Node_To_Node[NE * info->Total_Control_Point_to_mesh[Total_mesh] + 0] = info->Controlpoint_of_Element[info->NELOVER[i * MAX_N_ELEMENT_OVER + jj] * MAX_NO_CP_ON_ELEMENT + j];
										info->Total_Control_Point_To_Node[NE]++;
									}

									// 同じものがあったら, k > 0 以降の取得, kのカウント
									for (k = 0; k < info->Total_Control_Point_To_Node[NE]; k++)
									{
										if (info->Node_To_Node[NE * info->Total_Control_Point_to_mesh[Total_mesh] + k] == info->Controlpoint_of_Element[info->NELOVER[i * MAX_N_ELEMENT_OVER + jj] * MAX_NO_CP_ON_ELEMENT + j])
										{
											break;
										}
									}
									// 未設定のNode_To_Node取得
									if (k == info->Total_Control_Point_To_Node[NE])
									{
										info->Node_To_Node[NE * info->Total_Control_Point_to_mesh[Total_mesh] + k] = info->Controlpoint_of_Element[info->NELOVER[i * MAX_N_ELEMENT_OVER + jj] * MAX_NO_CP_ON_ELEMENT + j];
										info->Total_Control_Point_To_Node[NE]++;
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
				if (N + i < info->Total_Control_Point_to_mesh[Total_mesh])
				{
					for (j = 0; j < info->Total_Control_Point_To_Node[i]; j++)
					{
						int Min = info->Node_To_Node[i * info->Total_Control_Point_to_mesh[Total_mesh] + j], No = j;
						for (k = j; k < info->Total_Control_Point_To_Node[i]; k++)
						{
							if (Min > info->Node_To_Node[i * info->Total_Control_Point_to_mesh[Total_mesh] + k])
							{
								Min = info->Node_To_Node[i * info->Total_Control_Point_to_mesh[Total_mesh] + k];
								No = k;
							}
						}
						for (k = No; k > j; k--)
						{
							info->Node_To_Node[i * info->Total_Control_Point_to_mesh[Total_mesh] + k] = info->Node_To_Node[i * info->Total_Control_Point_to_mesh[Total_mesh] + k - 1];
						}
						info->Node_To_Node[i * info->Total_Control_Point_to_mesh[Total_mesh] + j] = Min;
					}
				}
			}

			// 節点からcol ptrを求める
			ii = 0;
			k = 0;
			for (i = 0; i < K_DIVISION_LENGE; i++)
			{
				for (ii = 0; ii < info->DIMENSION; ii++)
				{
					if (N + i < info->Total_Control_Point_to_mesh[Total_mesh])
					{
						i_index = info->Index_Dof[(N + i) * info->DIMENSION + ii];
						k = 0;
						if (i_index >= 0)
						{
							for (j = 0; j < info->Total_Control_Point_To_Node[i]; j++)
							{
								for (jj = 0; jj < info->DIMENSION; jj++)
								{
									j_index = info->Index_Dof[info->Node_To_Node[i * info->Total_Control_Point_to_mesh[Total_mesh] + j] * info->DIMENSION + jj];
									if (j_index >= 0 && j_index >= i_index)
									{
										info->K_Whole_Col[info->K_Whole_Ptr[i_index] + k] = j_index;
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
}


// K matrix の値を求める
void Make_K_Whole_Val(information *info)
{
	int i, j, j1, j2, k1, k2, l;
	int a, b, re;

	double *K_EL = (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * MAX_KIEL_SIZE);			// K_EL[MAX_KIEL_SIZE][MAX_KIEL_SIZE]
	double *coupled_K_EL= (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * MAX_KIEL_SIZE);		// coupled_K_EL[MAX_KIEL_SIZE][MAX_KIEL_SIZE]

	for (re = 0; re < info->real_Total_Element_to_mesh[Total_mesh]; re++)
	{
		i = info->real_element[re];

		if (Total_mesh == 1) // IGA
		{
			Make_gauss_array(0, info);
		}
		else if (Total_mesh >= 2) // S-IGA
		{
			if (info->Element_mesh[i] == 0 && re == 0)
			{
				Make_gauss_array(0, info);
			}
			else if (info->Element_mesh[i] > 0)
			{
				if (info->NNLOVER[i] == 1 && (info->NNLOVER[info->real_element[re - 1]] != 1 || info->Element_mesh[info->real_element[re - 1]] == 0))
				{
					Make_gauss_array(0, info);
				}
				else if (info->NNLOVER[i] >= 2 && (info->NNLOVER[info->real_element[re - 1]] == 1 || info->Element_mesh[info->real_element[re - 1]] == 0))
				{
					Make_gauss_array(1, info);
				}
			}
		}

		// 各要素のK_ELを求める
		Make_K_EL(i, K_EL, info);

		// Valを求める
		for (j1 = 0; j1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; j1++)
		{
			for (j2 = 0; j2 < info->DIMENSION; j2++)
			{
				a = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j1] * info->DIMENSION + j2];
				if (a >= 0)
				{
					for (k1 = 0; k1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; k1++)
					{
						for (k2 = 0; k2 < info->DIMENSION; k2++)
						{
							b = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k1] * info->DIMENSION + k2];
							if (b >= 0 && b >= a)
							{
								for (l = info->K_Whole_Ptr[a]; l < info->K_Whole_Ptr[a + 1]; l++)
								{
									if (info->K_Whole_Col[l] == b)
									{
										info->K_Whole_Val[l] += K_EL[(j1 * info->DIMENSION + j2) * MAX_KIEL_SIZE + k1 * info->DIMENSION + k2];
										break;
									}
								}
							}
						}
					}
				}
			}
		}

		// ローカルメッシュ上の要素について, 重なっている要素が存在するとき
		if (Total_mesh >= 2 && (info->Element_mesh[i] > 0 && info->NNLOVER[i] > 0))
		{
			for (j = 0; j < info->NNLOVER[i]; j++)
			{
				// 各要素のcoupled_K_ELを求める
				Make_coupled_K_EL(i, info->NELOVER[i * MAX_N_ELEMENT_OVER + j], coupled_K_EL, info);

				// Valを求める
				for (j1 = 0; j1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[info->NELOVER[i * MAX_N_ELEMENT_OVER + j]]]; j1++)
				{
					for (j2 = 0; j2 < info->DIMENSION; j2++)
					{
						a = info->Index_Dof[info->Controlpoint_of_Element[info->NELOVER[i * MAX_N_ELEMENT_OVER + j] * MAX_NO_CP_ON_ELEMENT + j1] * info->DIMENSION + j2];
						if (a >= 0)
						{
							for (k1 = 0; k1 < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; k1++)
							{
								for (k2 = 0; k2 < info->DIMENSION; k2++)
								{
									b = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + k1] * info->DIMENSION + k2];
									if (b >= 0 && b >= a)
									{
										for (l = info->K_Whole_Ptr[a]; l < info->K_Whole_Ptr[a + 1]; l++)
										{
											if (info->K_Whole_Col[l] == b)
											{
												info->K_Whole_Val[l] += coupled_K_EL[(j1 * info->DIMENSION + j2) * MAX_KIEL_SIZE + k1 * info->DIMENSION + k2];
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
void Make_K_EL(int El_No, double *K_EL, information *info)
{
	int i, j, k, l;
	
	int KIEL_SIZE = info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]] * info->DIMENSION;

	double *B = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * MAX_KIEL_SIZE);  // B[D_MATRIX_SIZE][MAX_KIEL_SIZE];
	double *K1 = (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * MAX_KIEL_SIZE); // K1[MAX_KIEL_SIZE][MAX_KIEL_SIZE];
	double J = 0.0;

	for (i = 0; i < KIEL_SIZE; i++)
	{
		for (j = 0; j < KIEL_SIZE; j++)
		{
			K_EL[i * MAX_KIEL_SIZE + j] = 0.0;
		}
	}

	for (i = 0; i < GP_ON_ELEMENT; i++)
	{
		// J, B の作成
		if (GP_1D == NG)
		{
			J = info->Jac[El_No * GP_ON_ELEMENT + i];
			for (j = 0; j < D_MATRIX_SIZE; j++)
			{
				for (k = 0; k < MAX_KIEL_SIZE; k++)
				{
					B[j * MAX_KIEL_SIZE + k] = info->B_Matrix[El_No * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + j * MAX_KIEL_SIZE + k];
				}
			}
		}
		else if (GP_1D == NG_EXTEND)
		{
			J = info->Jac_ex[El_No * GP_ON_ELEMENT + i];
			for (j = 0; j < D_MATRIX_SIZE; j++)
			{
				for (k = 0; k < MAX_KIEL_SIZE; k++)
				{
					B[j * MAX_KIEL_SIZE + k] = info->B_Matrix_ex[El_No * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + j * MAX_KIEL_SIZE + k];
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
void Make_coupled_K_EL(int El_No_loc, int El_No_glo, double *coupled_K_EL, information *info)
{
	int i, j, k, l;

	int BDBJ_flag = 0;
	int KIEL_SIZE = info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No_glo]] * info->DIMENSION;

	double *B = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * MAX_KIEL_SIZE);  // B[D_MATRIX_SIZE][MAX_KIEL_SIZE];
	double *BG = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * MAX_KIEL_SIZE); // BG[D_MATRIX_SIZE][MAX_KIEL_SIZE];
	double *K1 = (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * MAX_KIEL_SIZE); // K1[MAX_KIEL_SIZE][MAX_KIEL_SIZE];
	double J = 0.0;

	for (i = 0; i < KIEL_SIZE; i++)
	{
		for (j = 0; j < KIEL_SIZE; j++)
		{
			coupled_K_EL[i * MAX_KIEL_SIZE + j] = 0.0;
		}
	}

	for (i = 0; i < GP_ON_ELEMENT; i++)
	{
		// J, B, BG の作成
		double para[MAX_DIMENSION] = {0.0};
		double G_Gxi[MAX_DIMENSION] = {0.0};

		if (GP_1D == NG)
		{
			J = info->Jac[El_No_loc * GP_ON_ELEMENT + i];
			for (j = 0; j < D_MATRIX_SIZE; j++)
			{
				for (k = 0; k < MAX_KIEL_SIZE; k++)
				{
					B[j * MAX_KIEL_SIZE + k] = info->B_Matrix[El_No_loc * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + j * MAX_KIEL_SIZE + k];
				}
			}
			for (j = 0; j < info->DIMENSION; j++)
			{
				para[j] = info->Loc_parameter_on_Glo[El_No_loc * GP_ON_ELEMENT * info->DIMENSION + i * info->DIMENSION + j];
			}
		}
		else if (GP_1D == NG_EXTEND)
		{
			J = info->Jac_ex[El_No_loc * GP_ON_ELEMENT + i];
			for (j = 0; j < D_MATRIX_SIZE; j++)
			{
				for (k = 0; k < MAX_KIEL_SIZE; k++)
				{
					B[j * MAX_KIEL_SIZE + k] = info->B_Matrix_ex[El_No_loc * GP_ON_ELEMENT * D_MATRIX_SIZE * MAX_KIEL_SIZE + i * D_MATRIX_SIZE * MAX_KIEL_SIZE + j * MAX_KIEL_SIZE + k];
				}
			}
			for (j = 0; j < info->DIMENSION; j++)
			{
				para[j] = info->Loc_parameter_on_Glo_ex[El_No_loc * GP_ON_ELEMENT * info->DIMENSION + i * info->DIMENSION + j];
			}
		}

		// 要素内外判定
		if (info->DIMENSION == 2)
		{
			if (para[0] >= info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 0] + info->Order[0 * info->DIMENSION + 0] + info->ENC[El_No_glo * info->DIMENSION + 0]] &&
				para[0] <  info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 0] + info->Order[0 * info->DIMENSION + 0] + info->ENC[El_No_glo * info->DIMENSION + 0] + 1] &&
				para[1] >= info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 1] + info->Order[0 * info->DIMENSION + 1] + info->ENC[El_No_glo * info->DIMENSION + 1]] &&
				para[1] <  info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 1] + info->Order[0 * info->DIMENSION + 1] + info->ENC[El_No_glo * info->DIMENSION + 1] + 1])
			{
				BDBJ_flag = 1;

				// 親要素座標の算出
				G_Gxi[0] = - 1.0 + 2.0 * (para[0] - info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 0] + info->Order[0 * info->DIMENSION + 0] + info->ENC[El_No_glo * info->DIMENSION + 0]])
						 / (info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 0] + info->Order[0 * info->DIMENSION + 0] + info->ENC[El_No_glo * info->DIMENSION + 0] + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 0] + info->Order[0 * info->DIMENSION + 0] + info->ENC[El_No_glo * info->DIMENSION + 0]]);
				G_Gxi[1] = - 1.0 + 2.0 * (para[1] - info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 1] + info->Order[0 * info->DIMENSION + 1] + info->ENC[El_No_glo * info->DIMENSION + 1]])
						 / (info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 1] + info->Order[0 * info->DIMENSION + 1] + info->ENC[El_No_glo * info->DIMENSION + 1] + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 1] + info->Order[0 * info->DIMENSION + 1] + info->ENC[El_No_glo * info->DIMENSION + 1]]);
			}
			else
			{
				BDBJ_flag = 0;
			}
		}
		else if (info->DIMENSION == 3)
		{
			if (para[0] >= info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 0] + info->Order[0 * info->DIMENSION + 0] + info->ENC[El_No_glo * info->DIMENSION + 0]] &&
				para[0] <  info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 0] + info->Order[0 * info->DIMENSION + 0] + info->ENC[El_No_glo * info->DIMENSION + 0] + 1] &&
				para[1] >= info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 1] + info->Order[0 * info->DIMENSION + 1] + info->ENC[El_No_glo * info->DIMENSION + 1]] &&
				para[1] <  info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 1] + info->Order[0 * info->DIMENSION + 1] + info->ENC[El_No_glo * info->DIMENSION + 1] + 1] &&
				para[2] >= info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 2] + info->Order[0 * info->DIMENSION + 2] + info->ENC[El_No_glo * info->DIMENSION + 2]] &&
				para[2] <  info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 2] + info->Order[0 * info->DIMENSION + 2] + info->ENC[El_No_glo * info->DIMENSION + 2] + 1])
			{
				BDBJ_flag = 1;

				// 親要素座標の算出
				G_Gxi[0] = - 1.0 + 2.0 * (para[0] - info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 0] + info->Order[0 * info->DIMENSION + 0] + info->ENC[El_No_glo * info->DIMENSION + 0]])
						 / (info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 0] + info->Order[0 * info->DIMENSION + 0] + info->ENC[El_No_glo * info->DIMENSION + 0] + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 0] + info->Order[0 * info->DIMENSION + 0] + info->ENC[El_No_glo * info->DIMENSION + 0]]);
				G_Gxi[1] = - 1.0 + 2.0 * (para[1] - info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 1] + info->Order[0 * info->DIMENSION + 1] + info->ENC[El_No_glo * info->DIMENSION + 1]])
						 / (info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 1] + info->Order[0 * info->DIMENSION + 1] + info->ENC[El_No_glo * info->DIMENSION + 1] + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 1] + info->Order[0 * info->DIMENSION + 1] + info->ENC[El_No_glo * info->DIMENSION + 1]]);
				G_Gxi[2] = - 1.0 + 2.0 * (para[2] - info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 2] + info->Order[0 * info->DIMENSION + 2] + info->ENC[El_No_glo * info->DIMENSION + 2]])
						 / (info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 2] + info->Order[0 * info->DIMENSION + 2] + info->ENC[El_No_glo * info->DIMENSION + 2] + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 2] + info->Order[0 * info->DIMENSION + 2] + info->ENC[El_No_glo * info->DIMENSION + 2]]);
			}
			else
			{
				BDBJ_flag = 0;
			}
		}

		// 要素内であるとき, 結合要素剛性マトリックス計算
		if (BDBJ_flag)
		{
			// 重なるグローバル要素のBマトリックス
			Make_B_Matrix_anypoint(El_No_glo, BG, G_Gxi, info);
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


// BGマトリックスを求める
void Make_B_Matrix_anypoint(int El_No, double *B, double *Local_coord, information *info)
{
	double a_2x2[2][2], a_3x3[3][3];
	double *b = (double *)malloc(sizeof(double) * info->DIMENSION * MAX_NO_CP_ON_ELEMENT); // b[info->DIMENSION][MAX_NO_CP_ON_ELEMENT]

	int i, j, k;

	if (info->DIMENSION == 2)
	{
		for (i = 0; i < info->DIMENSION; i++)
		{
			for (j = 0; j < info->DIMENSION; j++)
			{
				a_2x2[i][j] = 0.0;
				for (k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; k++)
				{
					a_2x2[i][j] += dShape_func(k, j, Local_coord, El_No, info) * info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i];
				}
			}
		}

		InverseMatrix_2x2(a_2x2);

		for (i = 0; i < info->DIMENSION; i++)
		{
			for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; j++)
			{
				b[i * MAX_NO_CP_ON_ELEMENT + j] = 0.0;
				for (k = 0; k < info->DIMENSION; k++)
				{
					b[i * MAX_NO_CP_ON_ELEMENT + j] += a_2x2[k][i] * dShape_func(j, k, Local_coord, El_No, info);
				}
			}
		}

		for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; i++)
		{
			B[0 * MAX_KIEL_SIZE + 2 * i]     = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			B[0 * MAX_KIEL_SIZE + 2 * i + 1] = 0.0;

			B[1 * MAX_KIEL_SIZE + 2 * i]     = 0.0;
			B[1 * MAX_KIEL_SIZE + 2 * i + 1] = b[1 * MAX_NO_CP_ON_ELEMENT + i];

			B[2 * MAX_KIEL_SIZE + 2 * i]     = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			B[2 * MAX_KIEL_SIZE + 2 * i + 1] = b[0 * MAX_NO_CP_ON_ELEMENT + i];
		}
	}
	else if (info->DIMENSION == 3)
	{
		for (i = 0; i < info->DIMENSION; i++)
		{
			for (j = 0; j < info->DIMENSION; j++)
			{
				a_3x3[i][j] = 0.0;
				for (k = 0; k < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; k++)
				{
					a_3x3[i][j] += dShape_func(k, j, Local_coord, El_No, info) * info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + k] * (info->DIMENSION + 1) + i];
				}
			}
		}

		InverseMatrix_3x3(a_3x3);

		for (i = 0; i < info->DIMENSION; i++)
		{
			for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; j++)
			{
				b[i * MAX_NO_CP_ON_ELEMENT + j] = 0.0;
				for (k = 0; k < info->DIMENSION; k++)
				{
					b[i * MAX_NO_CP_ON_ELEMENT + j] += a_3x3[k][i] * dShape_func(j, k, Local_coord, El_No, info);
				}
			}
		}

		for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; i++)
		{
			B[0 * MAX_KIEL_SIZE + 3 * i]     = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			B[0 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			B[0 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;

			B[1 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			B[1 * MAX_KIEL_SIZE + 3 * i + 1] = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			B[1 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;

			B[2 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			B[2 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			B[2 * MAX_KIEL_SIZE + 3 * i + 2] = b[2 * MAX_NO_CP_ON_ELEMENT + i];

			B[3 * MAX_KIEL_SIZE + 3 * i]     = b[1 * MAX_NO_CP_ON_ELEMENT + i];
			B[3 * MAX_KIEL_SIZE + 3 * i + 1] = b[0 * MAX_NO_CP_ON_ELEMENT + i];
			B[3 * MAX_KIEL_SIZE + 3 * i + 2] = 0.0;

			B[4 * MAX_KIEL_SIZE + 3 * i]     = 0.0;
			B[4 * MAX_KIEL_SIZE + 3 * i + 1] = b[2 * MAX_NO_CP_ON_ELEMENT + i];
			B[4 * MAX_KIEL_SIZE + 3 * i + 2] = b[1 * MAX_NO_CP_ON_ELEMENT + i];

			B[5 * MAX_KIEL_SIZE + 3 * i]     = b[2 * MAX_NO_CP_ON_ELEMENT + i];
			B[5 * MAX_KIEL_SIZE + 3 * i + 1] = 0.0;
			B[5 * MAX_KIEL_SIZE + 3 * i + 2] = b[0 * MAX_NO_CP_ON_ELEMENT + i];
		}
	}

	free(b);
}


void BDBJ(int KIEL_SIZE, double *B, double J, double *K_EL, information *info)
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
				BD[i * D_MATRIX_SIZE + j] += B[k * MAX_KIEL_SIZE + i] * info->D[k * D_MATRIX_SIZE + j];
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


void coupled_BDBJ(int KIEL_SIZE, double *B, double *BG, double J, double *K_EL, information *info)
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
				BD[i * D_MATRIX_SIZE + j] += BG[k * MAX_KIEL_SIZE + i] * info->D[k * D_MATRIX_SIZE + j];
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
void Make_F_Vec(information *info)
{
	int i, index;
	for (i = 0; i < K_Whole_Size; i++)
	{
		info->rhs_vec[i] = 0.0;
	}
	for (i = 0; i < info->Total_Load_to_mesh[Total_mesh]; i++)
	{
		index = info->Index_Dof[info->Load_Node_Dir[i * 2 + 0] * info->DIMENSION + info->Load_Node_Dir[i * 2 + 1]];
		if (index >= 0)
		{
			info->rhs_vec[index] += info->Value_of_Load[i];
		}
	}
}


// 強制変位対策
void Make_F_Vec_disp_const(information *info)
{
	int i, ie, idir, inode, jdir, jnode, kk_const;
	int ii, iii, b, bb, jj, ii_local, jj_local;

	double *K_EL = (double *)malloc(sizeof(double) * MAX_KIEL_SIZE * MAX_KIEL_SIZE);	// K_EL[MAX_KIEL_SIZE][MAX_KIEL_SIZE]

	Make_gauss_array(0, info);

	for (ie = 0; ie < info->real_Total_Element_to_mesh[Total_mesh]; ie++)
	{
		i = info->real_element[ie];

		iii = 0;
		for (idir = 0; idir < info->DIMENSION; idir++)
		{
			for (inode = 0; inode < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; inode++)
			{
				b = info->Index_Dof[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + inode] * info->DIMENSION + idir];
				if (b < 0)
					iii++;
			}
		}

		if (iii > 0)
		{
			Make_K_EL(i, K_EL, info);
			for (idir = 0; idir < info->DIMENSION; idir++)
			{
				for (inode = 0; inode < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; inode++)
				{
					// ii = info->Controlpoint_of_Element[info->real_El_No_on_mesh[Total_mesh * info->Total_Element_to_mesh[Total_mesh] + ie] * MAX_NO_CP_ON_ELEMENT + inode] * info->DIMENSION + idir;
					ii = info->Controlpoint_of_Element[info->real_element[ie] * MAX_NO_CP_ON_ELEMENT + inode] * info->DIMENSION + idir;
					// ↑ real_El_No_on_mesh を real_element に変更
					b = info->Index_Dof[ii];
					if (b >= 0)
					{
						ii_local = inode * info->DIMENSION + idir;
						for (jdir = 0; jdir < info->DIMENSION; jdir++)
						{
							for (jnode = 0; jnode < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; jnode++)
							{
								jj = info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + jnode] * info->DIMENSION + jdir;
								bb = info->Index_Dof[jj];
								if (bb < 0)
								{
									jj_local = jnode * info->DIMENSION + jdir;
									for (kk_const = 0; kk_const < info->Total_Constraint_to_mesh[Total_mesh]; kk_const++)
									{
										if (info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + jnode] == info->Constraint_Node_Dir[kk_const * 2 + 0] && jdir == info->Constraint_Node_Dir[kk_const * 2 + 1])
										{
											info->rhs_vec[b] -= K_EL[ii_local * MAX_KIEL_SIZE + jj_local] * info->Value_of_Constraint[kk_const];
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
void Add_Equivalent_Nodal_Force_to_F_Vec(information *info)
{
	int i, j, index;
	for (j = 0; j < info->DIMENSION; j++)
	{
		for (i = 0; i < info->Total_Control_Point_to_mesh[Total_mesh]; i++)
		{
			index = info->Index_Dof[i * info->DIMENSION + j];
			if (index >= 0)
			{
				info->rhs_vec[index] += info->Equivalent_Nodal_Force[i * info->DIMENSION + j];
			}
		}
	}
}


// PCG solver, 前処理付共役勾配法(PCG)により[K]{d}={f}を解く, 対角スケーリングを行ったCG法を用いる
void PCG_Solver(int max_itetarion, double eps, information *info)
{
	int i, j, k;
	int ndof = K_Whole_Size;

	double *r = (double *)malloc(sizeof(double) * ndof);
	double *p = (double *)calloc(ndof, sizeof(double));
	double *y = (double *)malloc(sizeof(double) * ndof);
	double *r2 = (double *)calloc(ndof, sizeof(double));
	double *temp_array_K = (double *)malloc(sizeof(double) * ndof);

	double *gg = (double *)malloc(sizeof(double) * ndof);
	double *dd = (double *)malloc(sizeof(double) * ndof);
	double *pp = (double *)malloc(sizeof(double) * ndof);

	double *temp_r = (double *)malloc(sizeof(double) * ndof);

	// 初期化
	for (i = 0; i < ndof; i++)
		info->sol_vec[i] = 0.0;

	// 前処理行列作成
	double *M = (double *)malloc(sizeof(double) * info->K_Whole_Ptr[ndof]);
	int *M_Ptr = (int *)malloc(sizeof(int) * (ndof + 1));
	int *M_Col = (int *)malloc(sizeof(int) * info->K_Whole_Ptr[ndof]);
	double *M_diag = (double *)malloc(sizeof(double) * ndof);
	Make_M(M, M_Ptr, M_Col, M_diag, ndof, info);

	// 第0近似解に対する残差の計算
	// double *ax = (double *)calloc(ndof, sizeof(double));
	// for (i = 0; i < ndof; i++)
	// {
	// 	for (j = info->K_Whole_Ptr[i]; j < info->K_Whole_Ptr[i + 1]; j++)
	// 	{
	// 		ax[i] += info->K_Whole_Val[j] * info->sol_vec[info->K_Whole_Col[j]];
	// 		if (i != info->K_Whole_Col[j])
	// 		{
	// 			ax[info->K_Whole_Col[j]] += info->K_Whole_Val[j] * info->sol_vec[i];
	// 		}
	// 	}
	// }
	// for (i = 0; i < ndof; i++)
	// {
	// 	r[i] = info->rhs_vec[i] - ax[i];
	// }
	// free(ax);

	// 第0近似解に対する残差の計算
	for (i = 0; i < ndof; i++)
	{
		r[i] = info->rhs_vec[i];
	}

	// p_0 = (LDL^T)^-1 r_0 の計算 <- CG法で M = [[K^G, 0], [0, K^L]] とし, p_0 = (LDL^T)^-1 r_0 = M^-1 r_0
	CG(ndof, p, M, M_Ptr, M_Col, M_diag, r, gg, dd, pp, temp_r);

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
			for (j = 0; j < ndof; j++)
			{
				temp_array_K[j] = 0.0;
			}
			
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
					temp_array_K[j] = info->K_Whole_Val[temp1];
				}
			}
			y[i] = inner_product(ndof, temp_array_K, p);
		}

		// alpha = r*r/(P*AP)の計算
		double temp_scaler = inner_product(ndof, p, y);
		alpha = rr0 / temp_scaler;
		// printf("alpha %le\n", alpha);

		// 解x, 残差rの更新
		for (i = 0; i < ndof; i++)
		{
			info->sol_vec[i] += alpha * p[i];
			r[i] -= alpha * y[i];
		}

		// (r*r)_(k+1)の計算
		CG(ndof, r2, M, M_Ptr, M_Col, M_diag, r, gg, dd, pp, temp_r);

		// rr1 = inner_product(ndof, r, r2); // 旧
		// rr1 = inner_product(ndof, y, r2); // 新
		// printf("rr1 %le\n", rr1);

		// 収束判定 (||r|| <= eps)
		// double rr1 = inner_product(ndof, y, r2);
		// e = sqrt(fabs(rr1));
		// if (e < eps)
		// {
		//     k++;
		//     break;
		// }

		// 収束判定 (CG法と同じ)
		double e1 = 0.0, e2 = 0.0;
		for (i = 0; i < ndof; i++)
		{
			e1 += p[i] * p[i];
			e2 += info->sol_vec[i] * info->sol_vec[i];
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

		printf("PCG itr %d\t", k);
		printf("eps %le", e);
		// if (rr1 < 0)
		// {
		// 	printf("\t rr1 < 0");
		// }
	}

	int max_itr_result = k;
	double eps_result = e;

	printf("\nndof = %d\n", ndof);
	printf("itr_result = %d\n", max_itr_result);
	printf("eps_result = %.15e\n", eps_result);

	free(r), free(p), free(y), free(r2), free(temp_array_K);
	free(M), free(M_Ptr), free(M_Col), free(M_diag);
	free(gg), free(dd), free(pp), free(temp_r);
}


void Make_M(double *M, int *M_Ptr, int *M_Col, double *M_diag, int ndof, information *info)
{
	int i, j;
	int ndof_glo = 0;

	// グローバルパッチのdofを求める
	for (i = 0; i < info->Total_Control_Point_on_mesh[0] * info->DIMENSION; i++)
	{
		if (info->Index_Dof[i] != ERROR)
		{
			ndof_glo++;
		}
	}
	printf("ndof		%d\n", ndof);
	printf("ndof_glo	%d\n", ndof_glo);
	printf("ndof_loc	%d\n", ndof - ndof_glo);

	int counter = 0;

	// M = [[K^G, 0], [0, K^L]] と M_diag を作成
	M_Ptr[0] = 0;
	for (i = 0; i < ndof; i++)
	{
		M_Ptr[i + 1] = M_Ptr[i];

		for (j = info->K_Whole_Ptr[i]; j < info->K_Whole_Ptr[i + 1]; j++)
		{
			if (i < ndof_glo && info->K_Whole_Col[j] < ndof_glo)
			{
				M[counter] = info->K_Whole_Val[j];
				M_Col[counter] = info->K_Whole_Col[j];
				counter++;
				M_Ptr[i + 1]++;
			}
			else if (i >= ndof_glo)
			{
				M[counter] = info->K_Whole_Val[j];
				M_Col[counter] = info->K_Whole_Col[j];
				counter++;
				M_Ptr[i + 1]++;
			}
		}
	}

	// diag scaling preprocess
	M_diag[0] = 1.0 / sqrt(M[0]);
	for (i = 1; i < ndof; i++)
	{
		M_diag[i] = 1.0 / sqrt(M[M_Ptr[i]]);
	}

	counter = 0;
	for (i = 0; i < ndof; i++)
	{
		for (j = M_Ptr[i]; j < M_Ptr[i + 1]; j++)
		{
			M[counter] = M[counter] * M_diag[i] * M_diag[M_Col[j]];
			counter++;
		}
	}
}


void CG(int ndof, double *solution_vec, double *M, int *M_Ptr, int *M_Col, double *M_diag, double *right_vec, double *gg, double *dd, double *pp, double *temp_r)
{
	double qqq, ppp, rrr;
	double alphak, betak;
	int i, j, ii, itr, istop = 0;
	int max_itr = ndof;
	double eps = 1.0e-13;
	double rrr3 = 0.0;

	// diag scaling preprocess
	for (i = 0; i < ndof; i++)
	{
		temp_r[i] = right_vec[i] * M_diag[i];
	}

	// CG solver (diag scaling)
	for (i = 0; i < ndof; i++)
	{
		solution_vec[i] = 0.0;
	}
	M_mat_vec_crs(M, M_Ptr, M_Col, dd, solution_vec, ndof);
	for (i = 0; i < ndof; i++)
	{
		gg[i] = temp_r[i] - dd[i];
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

		// M_check_conv_CG
		double rrr1 = 0.0, rrr2 = 0.0;
		for (j = 0; j < ndof; j++)
		{
			rrr1 += pp[j] * pp[j];
			rrr2 += solution_vec[j] * solution_vec[j];
		}
		rrr3 = fabs(alphak) * sqrt(rrr1 / rrr2);
		if (rrr3 < eps)
		{
			istop = 1;
		}

		if (istop == 1)
			break;
	}

	// diag scaling postprocess
	for (i = 0; i < ndof; i++)
	{
		solution_vec[i] = solution_vec[i] * M_diag[i];
	}

	printf("\t\tCG itr %d\teps %le\n", itr, rrr3);
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


int RowCol_to_icount(int row, int col, information *info)
{
	for (int j = info->K_Whole_Ptr[row]; j < info->K_Whole_Ptr[row + 1]; j++)
	{
		if (info->K_Whole_Col[j] == col)
		{
			return j;
		}
		else if (info->K_Whole_Col[j] > col)
		{
			return -1;
		}
	}
	return -1;
}


// tool
double InverseMatrix_2x2(double Matrix[2][2])
{
	int i, j;
	double a[2][2];
	double det = Matrix[0][0] * Matrix[1][1] - Matrix[0][1] * Matrix[1][0];

	if (det == 0)
		return ERROR;

	for (i = 0; i < 2; i++)
	{
		for (j = 0; j < 2; j++)
			a[i][j] = Matrix[i][j];
	}
	Matrix[0][0] = a[1][1] / det;
	Matrix[0][1] = a[0][1] * (-1) / det;
	Matrix[1][0] = a[1][0] * (-1) / det;
	Matrix[1][1] = a[0][0] / det;

	return det;
}


double InverseMatrix_3x3(double Matrix[3][3])
{
	int i, j;
	double a[3][3];
	double det = Matrix[0][0] * Matrix[1][1] * Matrix[2][2] + Matrix[1][0] * Matrix[2][1] * Matrix[0][2] + Matrix[2][0] * Matrix[0][1] * Matrix[1][2] - Matrix[0][0] * Matrix[2][1] * Matrix[1][2] - Matrix[2][0] * Matrix[1][1] * Matrix[0][2] - Matrix[1][0] * Matrix[0][1] * Matrix[2][2];

	if (det == 0)
		return ERROR;

	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
			a[i][j] = Matrix[i][j];
	}
	Matrix[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) / det;
	Matrix[0][1] = (a[0][2] * a[2][1] - a[0][1] * a[2][2]) / det;
	Matrix[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) / det;
	Matrix[1][0] = (a[1][2] * a[2][0] - a[1][0] * a[2][2]) / det;
	Matrix[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]) / det;
	Matrix[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) / det;
	Matrix[2][0] = (a[1][0] * a[2][1] - a[1][1] * a[2][0]) / det;
	Matrix[2][1] = (a[0][1] * a[2][0] - a[0][0] * a[2][1]) / det;
	Matrix[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / det;

	return det;
}


// Shape Function
double Shape_func(int I_No, double *Local_coord, int El_No, information *info)
{
	int i, j;
	double R, weight_func = 0.0;

	double Position_Data_param[MAX_DIMENSION];

	double *shape_func = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT);
	double *Shape = (double *)malloc(sizeof(double) * info->DIMENSION * info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER); // Shape[info->DIMENSION][MAX_N_NODE][10]
	double *dShape = (double *)malloc(sizeof(double) * info->DIMENSION * info->Total_Control_Point_to_mesh[Total_mesh]); // dShape[info->DIMENSION][MAX_N_NODE]

	for (i = 0; i < MAX_NO_CP_ON_ELEMENT; i++)
	{
		shape_func[i] = 1.0;
	}

	for (j = 0; j < info->DIMENSION; j++)
	{
		ShapeFunc_from_paren(Position_Data_param, Local_coord, j, El_No, info);
		ShapeFunction1D(Position_Data_param, j, El_No, Shape, dShape, info);
		for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; i++)
		{
			shape_func[i] *= Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + j] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + j]];
		}
	}
	for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; i++)
	{
		weight_func += shape_func[i] * info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION];
	}
	free(Shape), free(dShape);

	if (I_No < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]])
		R = shape_func[I_No] * info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + I_No] * (info->DIMENSION + 1) + info->DIMENSION] / weight_func;

	else
		R = ERROR;

	free(shape_func);
	return R;
}


void ShapeFunc_from_paren(double *Position_Data_param, double *Local_coord, int j, int e, information *info)
{
	int i = info->INC[info->Element_patch[e] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + 0] * info->DIMENSION + j];
	Position_Data_param[j] = ((info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + i + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + i]) * Local_coord[j] + (info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + i + 1] + info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + i])) / 2.0;
}


void ShapeFunction1D(double *Position_Data_param, int j, int e, double *Shape, double *dShape, information *info)
{
	int ii;
	int p;

	for (ii = 0; ii < info->No_knot[info->Element_patch[e] * info->DIMENSION + j]; ii++)
	{
		if (info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii] == info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + 1])
		{
			Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + 0] = 0.0;
		}
		else if (info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii] != info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + 1] && info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii] <= Position_Data_param[j] && Position_Data_param[j] < info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + 1])
		{
			Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + 0] = 1.0;
		}
		else if (info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii] != info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + 1] && info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + 1] == info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + (info->No_knot[info->Element_patch[e] * info->DIMENSION + j] - 1)] && info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii] <= Position_Data_param[j] && Position_Data_param[j] <= info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + 1])
		{
			Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + 0] = 1.0;
		}
		else
		{
			Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + 0] = 0.0;
		}
	}

	for (ii = 0; ii < info->No_knot[info->Element_patch[e] * info->DIMENSION + j]; ii++)
	{
		for (p = 1; p <= info->Order[info->Element_patch[e] * info->DIMENSION + j]; p++)
		{
			Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + p] = 0.0;
		}
	}

	double left_term, right_term;
	for (p = 1; p <= info->Order[info->Element_patch[e] * info->DIMENSION + j]; p++)
	{
		for (ii = 0; ii < info->No_knot[info->Element_patch[e] * info->DIMENSION + j]; ii++)
		{
			left_term = 0.0;
			right_term = 0.0;

			if ((Position_Data_param[j] - info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii]) * Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + p - 1] == 0 && info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + p] - info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii] == 0)
				left_term = 0.0;
			else
			{
				left_term = (Position_Data_param[j] - info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii]) / (info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + p] - info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii]) * Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + p - 1];
			}
			if ((info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + p + 1] - Position_Data_param[j]) * Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + (ii + 1) * MAX_ORDER + p - 1] == 0 && info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + p + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + 1] == 0)
				right_term = 0.0;
			else
			{
				right_term = (info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + p + 1] - Position_Data_param[j]) / (info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + p + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + 1]) * Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + (ii + 1) * MAX_ORDER + p - 1];
			}
			Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + p] = left_term + right_term;
		}
	}

	double dleft_term, dright_term;
	for (ii = 0; ii < info->No_Control_point[info->Element_patch[e] * info->DIMENSION + j] + 1; ii++)
	{
		dleft_term = 0.0;
		dright_term = 0.0;

		if (info->Order[info->Element_patch[e] * info->DIMENSION + j] * Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + info->Order[info->Element_patch[e] * info->DIMENSION + j] - 1] == 0 && info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + info->Order[info->Element_patch[e] * info->DIMENSION + j]] - info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii] == 0)
			dleft_term = 0.0;
		else
			dleft_term = info->Order[info->Element_patch[e] * info->DIMENSION + j] / (info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + info->Order[info->Element_patch[e] * info->DIMENSION + j]] - info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii]) * Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + ii * MAX_ORDER + info->Order[info->Element_patch[e] * info->DIMENSION + j] - 1];

		if (info->Order[info->Element_patch[e] * info->DIMENSION + j] * Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + (ii + 1) * MAX_ORDER + info->Order[info->Element_patch[e] * info->DIMENSION + j] - 1] == 0 && info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + info->Order[info->Element_patch[e] * info->DIMENSION + j] + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + 1] == 0)
			dright_term = 0.0;
		else
			dright_term = info->Order[info->Element_patch[e] * info->DIMENSION + j] / (info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + info->Order[info->Element_patch[e] * info->DIMENSION + j] + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + ii + 1]) * Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + (ii + 1) * MAX_ORDER + info->Order[info->Element_patch[e] * info->DIMENSION + j] - 1];

		dShape[j * info->Total_Control_Point_to_mesh[Total_mesh] + ii] = dleft_term - dright_term;
	}
}


double dShape_func(int I_No, int xez, double *Local_coord, int El_No, information *info)
{
	double dR = 0.0;

	if (info->DIMENSION == 2)
	{
		double *dShape_func1 = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT); // dShape_func1[MAX_N_NODE];
		double *dShape_func2 = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT); // dShape_func2[MAX_N_NODE];

		NURBS_deriv_2D(Local_coord, El_No, dShape_func1, dShape_func2, info);

		if (xez != 0 && xez != 1)
		{
			dR = ERROR;
		}
		else if (I_No < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]])
		{
			if (xez == 0)
			{
				dR = dShape_func1[I_No] * dShapeFunc_from_paren(xez, El_No, info);
			}
			else if (xez == 1)
			{
				dR = dShape_func2[I_No] * dShapeFunc_from_paren(xez, El_No, info);
			}
		}
		else
		{
			dR = ERROR;
		}

		free(dShape_func1), free(dShape_func2);
	}
	else if (info->DIMENSION == 3)
	{
		double *dShape_func1 = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT); // dShape_func1[MAX_N_NODE];
		double *dShape_func2 = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT); // dShape_func2[MAX_N_NODE];
		double *dShape_func3 = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT); // dShape_func3[MAX_N_NODE];

		NURBS_deriv_3D(Local_coord, El_No, dShape_func1, dShape_func2, dShape_func3, info);

		if (xez != 0 && xez != 1 && xez != 2)
		{
			dR = ERROR;
		}
		else if (I_No < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]])
		{
			if (xez == 0)
			{
				dR = dShape_func1[I_No] * dShapeFunc_from_paren(xez, El_No, info);
			}
			else if (xez == 1)
			{
				dR = dShape_func2[I_No] * dShapeFunc_from_paren(xez, El_No, info);
			}
			else if (xez == 2)
			{
				dR = dShape_func3[I_No] * dShapeFunc_from_paren(xez, El_No, info);
			}
		}
		else
		{
			dR = ERROR;
		}

		free(dShape_func1), free(dShape_func2), free(dShape_func3);
	}

	return dR;
}


void NURBS_deriv_2D(double *Local_coord, int El_No, double *dShape_func1, double *dShape_func2, information *info)
{
	int i, j;
	double weight_func = 0.0;
	double dWeight_func1 = 0.0;
	double dWeight_func2 = 0.0;

	double Position_Data_param[MAX_DIMENSION] = {0.0};

	double *shape_func = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT);
	double *Shape = (double *)malloc(sizeof(double) * info->DIMENSION * info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER); // Shape[info->DIMENSION][MAX_N_NODE][10]
	double *dShape = (double *)malloc(sizeof(double) * info->DIMENSION * info->Total_Control_Point_to_mesh[Total_mesh]);			// dShape[info->DIMENSION][MAX_N_NODE]

	for (i = 0; i < MAX_NO_CP_ON_ELEMENT; i++)
	{
		shape_func[i] = 1.0;
	}

	for (j = 0; j < info->DIMENSION; j++)
	{
		ShapeFunc_from_paren(Position_Data_param, Local_coord, j, El_No, info);
		ShapeFunction1D(Position_Data_param, j, El_No, Shape, dShape, info);
		for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; i++)
		{
			shape_func[i] *= Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + j] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + j]];
		}
	}
	for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; i++)
	{
		weight_func += shape_func[i] * info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION];
	}
	for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; i++)
	{
		dWeight_func1 += dShape[0 * info->Total_Control_Point_to_mesh[Total_mesh] + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 0]] * Shape[1 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 1] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 1]] * info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION];
		dWeight_func2 += Shape[0 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 0] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 0]] * dShape[1 * info->Total_Control_Point_to_mesh[Total_mesh] + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 1]] * info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION];
	}
	for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; i++)
	{
		dShape_func1[i] = info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION] * (weight_func * dShape[0 * info->Total_Control_Point_to_mesh[Total_mesh] + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 0]] * Shape[1 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 1] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 1]] - dWeight_func1 * shape_func[i]) / (weight_func * weight_func);
		dShape_func2[i] = info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION] * (weight_func * Shape[0 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 0] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 0]] * dShape[1 * info->Total_Control_Point_to_mesh[Total_mesh] + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 1]] - dWeight_func2 * shape_func[i]) / (weight_func * weight_func);
	}

	free(Shape), free(dShape), free(shape_func);
}


void NURBS_deriv_3D(double *Local_coord, int El_No, double *dShape_func1, double *dShape_func2, double *dShape_func3, information *info)
{
	int i, j;
	double weight_func = 0.0;
	double dWeight_func1 = 0.0;
	double dWeight_func2 = 0.0;
	double dWeight_func3 = 0.0;

	double Position_Data_param[MAX_DIMENSION] = {0.0};

	double *shape_func = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT);
	double *Shape = (double *)malloc(sizeof(double) * info->DIMENSION * info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER); // Shape[info->DIMENSION][MAX_N_NODE][10]
	double *dShape = (double *)malloc(sizeof(double) * info->DIMENSION * info->Total_Control_Point_to_mesh[Total_mesh]);			// dShape[info->DIMENSION][MAX_N_NODE]

	for (i = 0; i < MAX_NO_CP_ON_ELEMENT; i++)
	{
		shape_func[i] = 1.0;
	}

	for (j = 0; j < info->DIMENSION; j++)
	{
		ShapeFunc_from_paren(Position_Data_param, Local_coord, j, El_No, info);
		ShapeFunction1D(Position_Data_param, j, El_No, Shape, dShape, info);
		for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; i++)
		{
			shape_func[i] *= Shape[j * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + j] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + j]];
		}
	}
	for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; i++)
	{
		weight_func += shape_func[i] * info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION];
	}
	for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; i++)
	{
		dWeight_func1 += dShape[0 * info->Total_Control_Point_to_mesh[Total_mesh] + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 0]]
					   * Shape[1 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 1] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 1]]
					   * Shape[2 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 2] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 2]]
					   * info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION];

		dWeight_func2 += Shape[0 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 0] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 0]]
					   * dShape[1 * info->Total_Control_Point_to_mesh[Total_mesh] + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 1]]
					   * Shape[2 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 2] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 2]]
					   * info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION];

		dWeight_func3 += Shape[0 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 0] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 0]]
					   * Shape[1 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 1] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 1]]
					   * dShape[2 * info->Total_Control_Point_to_mesh[Total_mesh] + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 2]]
					   * info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION];
	}
	for (i = 0; i < info->No_Control_point_ON_ELEMENT[info->Element_patch[El_No]]; i++)
	{
		dShape_func1[i]
			= info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION]
			* (weight_func
			* dShape[0 * info->Total_Control_Point_to_mesh[Total_mesh] + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 0]]
			* Shape[1 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 1] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 1]]
			* Shape[2 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 2] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 2]]
			- dWeight_func1	* shape_func[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i]])
			/ (weight_func * weight_func);

		dShape_func2[i]
			= info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION]
			* (weight_func
			* Shape[0 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 0] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 0]]
			* dShape[1 * info->Total_Control_Point_to_mesh[Total_mesh] + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 1]]
			* Shape[2 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 2] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 1]]
			- dWeight_func2	* shape_func[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i]])
			/ (weight_func * weight_func);

		dShape_func3[i]
			= info->Node_Coordinate[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * (info->DIMENSION + 1) + info->DIMENSION]
			* (weight_func
			* Shape[0 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 0] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 0]]
			* Shape[1 * (info->Total_Control_Point_to_mesh[Total_mesh] * MAX_ORDER) + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 1] * MAX_ORDER + info->Order[info->Element_patch[El_No] * info->DIMENSION + 1]]
			* dShape[2 * info->Total_Control_Point_to_mesh[Total_mesh] + info->INC[info->Element_patch[El_No] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i] * info->DIMENSION + 2]]
			- dWeight_func3 * shape_func[info->Controlpoint_of_Element[El_No * MAX_NO_CP_ON_ELEMENT + i]])
			/ (weight_func * weight_func);
	}

	free(Shape), free(dShape), free(shape_func);
}


double dShapeFunc_from_paren(int j, int e, information *info)
{
	int i = info->INC[info->Element_patch[e] * info->Total_Control_Point_to_mesh[Total_mesh] * info->DIMENSION + info->Controlpoint_of_Element[e * MAX_NO_CP_ON_ELEMENT + 0] * info->DIMENSION + j];;
	double dPosition_Data_param = (info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + i + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[info->Element_patch[e] * info->DIMENSION + j] + i]) / 2.0;
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

			// for (int temp_i = 0; temp_i < info->No_knot[0][0]; temp_i++)
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
					 double *output_dxi_x, double *output_data_x,
					 double *output_dxi_y, double *output_data_y)
{
	int i, j, temp_index;
	double temp1, temp2, temp3;
	double molecule_x, molecule_y;
	double dxi_molecule_x, dxi_molecule_y;
	double data_molecule_x, data_molecule_y;
	double denominator, dxi_denominator, data_denominator;
	double temp_output_xi, temp_output_eta;
	double temp_d_output_xi, temp_d_output_eta;
	molecule_x = 0.0;
	molecule_y = 0.0;
	denominator = 0.0;
	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_denominator = 0.0;
	data_molecule_x = 0.0;
	data_molecule_y = 0.0;
	data_denominator = 0.0;

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
			data_molecule_x += temp3 * cntl_px[temp_index];
			data_molecule_y += temp3 * cntl_py[temp_index];
			data_denominator += temp3;
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;

	temp1 = denominator * denominator;
	(*output_dxi_x) = (dxi_molecule_x * denominator - molecule_x * dxi_denominator) / temp1;
	(*output_dxi_y) = (dxi_molecule_y * denominator - molecule_y * dxi_denominator) / temp1;
	(*output_data_x) = (data_molecule_x * denominator - molecule_x * data_denominator) / temp1;
	(*output_data_y) = (data_molecule_y * denominator - molecule_y * data_denominator) / temp1;
	return denominator;
}


double rNURBS_surface(double *input_knot_vec_xi, double *input_knot_vec_eta,
					  double *cntl_px, double *cntl_py,
					  int cntl_p_n_xi, int cntl_p_n_eta,
					  double *weight, int order_xi, int order_eta,
					  double xi, double eta,
					  double *output_x, double *output_y,
					  double *output_dxi_x, double *output_data_x,
					  double *output_dxi_y, double *output_data_y)
{
	int i, j, temp_index;
	double temp1, temp2, temp3;
	double molecule_x, molecule_y;
	double dxi_molecule_x, dxi_molecule_y;
	double data_molecule_x, data_molecule_y;
	double denominator, dxi_denominator, data_denominator;
	double temp_output_xi, temp_output_eta;
	double temp_d_output_xi, temp_d_output_eta;
	molecule_x = 0.0;
	molecule_y = 0.0;
	denominator = 0.0;
	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_denominator = 0.0;
	data_molecule_x = 0.0;
	data_molecule_y = 0.0;
	data_denominator = 0.0;

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
			data_molecule_x += temp3 * cntl_px[temp_index];
			data_molecule_y += temp3 * cntl_py[temp_index];
			data_denominator += temp3;
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;

	temp1 = denominator * denominator;
	(*output_dxi_x) = (dxi_molecule_x * denominator - molecule_x * dxi_denominator) / temp1;
	(*output_dxi_y) = (dxi_molecule_y * denominator - molecule_y * dxi_denominator) / temp1;
	(*output_data_x) = (data_molecule_x * denominator - molecule_x * data_denominator) / temp1;
	(*output_data_y) = (data_molecule_y * denominator - molecule_y * data_denominator) / temp1;
	return denominator;
}


double lNURBS_surface(double *input_knot_vec_xi, double *input_knot_vec_eta,
					  double *cntl_px, double *cntl_py,
					  int cntl_p_n_xi, int cntl_p_n_eta,
					  double *weight, int order_xi, int order_eta,
					  double xi, double eta,
					  double *output_x, double *output_y,
					  double *output_dxi_x, double *output_data_x,
					  double *output_dxi_y, double *output_data_y)
{
	int i, j, temp_index;
	double temp1, temp2, temp3;
	double molecule_x, molecule_y;
	double dxi_molecule_x, dxi_molecule_y;
	double data_molecule_x, data_molecule_y;
	double denominator, dxi_denominator, data_denominator;
	double temp_output_xi, temp_output_eta;
	double temp_d_output_xi, temp_d_output_eta;
	molecule_x = 0.0;
	molecule_y = 0.0;
	denominator = 0.0;
	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_denominator = 0.0;
	data_molecule_x = 0.0;
	data_molecule_y = 0.0;
	data_denominator = 0.0;

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
			data_molecule_x += temp3 * cntl_px[temp_index];
			data_molecule_y += temp3 * cntl_py[temp_index];
			data_denominator += temp3;
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;

	temp1 = denominator * denominator;
	(*output_dxi_x) = (dxi_molecule_x * denominator - molecule_x * dxi_denominator) / temp1;
	(*output_dxi_y) = (dxi_molecule_y * denominator - molecule_y * dxi_denominator) / temp1;
	(*output_data_x) = (data_molecule_x * denominator - molecule_x * data_denominator) / temp1;
	(*output_data_y) = (data_molecule_y * denominator - molecule_y * data_denominator) / temp1;
	return denominator;
}


double rlNURBS_surface(double *input_knot_vec_xi, double *input_knot_vec_eta,
					   double *cntl_px, double *cntl_py,
					   int cntl_p_n_xi, int cntl_p_n_eta,
					   double *weight, int order_xi, int order_eta,
					   double xi, double eta,
					   double *output_x, double *output_y,
					   double *output_dxi_x, double *output_data_x,
					   double *output_dxi_y, double *output_data_y)
{
	int i, j, temp_index;
	double temp1, temp2, temp3;
	double molecule_x, molecule_y;
	double dxi_molecule_x, dxi_molecule_y;
	double data_molecule_x, data_molecule_y;
	double denominator, dxi_denominator, data_denominator;
	double temp_output_xi, temp_output_eta;
	double temp_d_output_xi, temp_d_output_eta;
	molecule_x = 0.0;
	molecule_y = 0.0;
	denominator = 0.0;
	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_denominator = 0.0;
	data_molecule_x = 0.0;
	data_molecule_y = 0.0;
	data_denominator = 0.0;

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
			data_molecule_x += temp3 * cntl_px[temp_index];
			data_molecule_y += temp3 * cntl_py[temp_index];
			data_denominator += temp3;
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;

	temp1 = denominator * denominator;
	(*output_dxi_x) = (dxi_molecule_x * denominator - molecule_x * dxi_denominator) / temp1;
	(*output_dxi_y) = (dxi_molecule_y * denominator - molecule_y * dxi_denominator) / temp1;
	(*output_data_x) = (data_molecule_x * denominator - molecule_x * data_denominator) / temp1;
	(*output_data_y) = (data_molecule_y * denominator - molecule_y * data_denominator) / temp1;
	return denominator;
}


double lrNURBS_surface(double *input_knot_vec_xi, double *input_knot_vec_eta,
					   double *cntl_px, double *cntl_py,
					   int cntl_p_n_xi, int cntl_p_n_eta,
					   double *weight, int order_xi, int order_eta,
					   double xi, double eta,
					   double *output_x, double *output_y,
					   double *output_dxi_x, double *output_data_x,
					   double *output_dxi_y, double *output_data_y)
{
	int i, j, temp_index;
	double temp1, temp2, temp3;
	double molecule_x, molecule_y;
	double dxi_molecule_x, dxi_molecule_y;
	double data_molecule_x, data_molecule_y;
	double denominator, dxi_denominator, data_denominator;
	double temp_output_xi, temp_output_eta;
	double temp_d_output_xi, temp_d_output_eta;
	molecule_x = 0.0;
	molecule_y = 0.0;
	denominator = 0.0;
	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_denominator = 0.0;
	data_molecule_x = 0.0;
	data_molecule_y = 0.0;
	data_denominator = 0.0;

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
			data_molecule_x += temp3 * cntl_px[temp_index];
			data_molecule_y += temp3 * cntl_py[temp_index];
			data_denominator += temp3;
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;

	temp1 = denominator * denominator;
	(*output_dxi_x) = (dxi_molecule_x * denominator - molecule_x * dxi_denominator) / temp1;
	(*output_dxi_y) = (dxi_molecule_y * denominator - molecule_y * dxi_denominator) / temp1;
	(*output_data_x) = (data_molecule_x * denominator - molecule_x * data_denominator) / temp1;
	(*output_data_y) = (data_molecule_y * denominator - molecule_y * data_denominator) / temp1;
	return denominator;
}


double rrrNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta,
                       double *cntl_px, double *cntl_py, double *cntl_pz,
                       int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta,
                       double *weight, int order_xi, int order_eta, int order_zeta,
                       double xi, double eta, double zeta,
                       double *output_x, double *output_y, double *output_z,
                       double *output_dxi_x, double *output_data_x, double *output_dzeta_x,
                       double *output_dxi_y, double *output_data_y, double *output_dzeta_y,
                       double *output_dxi_z, double *output_data_z, double *output_dzeta_z)
{
	int i, j, k, temp_index = 0;
	double temp1, temp2, temp3, temp4;
	double molecule_x, molecule_y, molecule_z;
	double dxi_molecule_x, dxi_molecule_y, dxi_molecule_z;
	double data_molecule_x, data_molecule_y, data_molecule_z;
	double dzeta_molecule_x, dzeta_molecule_y, dzeta_molecule_z;
	double denominator, dxi_denominator, data_denominator, dzeta_denominator;
	double temp_output_xi, temp_output_eta, temp_output_zeta;
	double temp_d_output_xi, temp_d_output_eta, temp_d_output_zeta;

	molecule_x = 0.0;
	molecule_y = 0.0;
	molecule_z = 0.0;
	denominator = 0.0;

	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_molecule_z = 0.0;
	dxi_denominator = 0.0;

	data_molecule_x = 0.0;
	data_molecule_y = 0.0;
	data_molecule_z = 0.0;
	data_denominator = 0.0;

	dzeta_molecule_x = 0.0;
	dzeta_molecule_y = 0.0;
	dzeta_molecule_z = 0.0;
	dzeta_denominator = 0.0;

	int index_min_xi   = 0;
	int index_max_xi   = cntl_p_n_xi - 1;
	int index_min_eta  = 0;
	int index_max_eta  = cntl_p_n_eta - 1;
	int index_min_zeta = 0;
	int index_max_zeta = cntl_p_n_zeta - 1;

	// xi
	for (i = 0; i < cntl_p_n_xi; i++)
	{
		if (knot_vec_xi[i + 1] >= xi)
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}

	if (index_min_xi < 0)
	{
		index_min_xi = 0;
		index_max_xi = order_xi + 1;
	}

	if (index_max_xi > cntl_p_n_xi)
	{
		index_min_xi = cntl_p_n_xi - order_xi - 1;
		index_max_xi = cntl_p_n_xi;
	}

	// eta
	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if (knot_vec_eta[i + 1] >= eta)
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0)
	{
		index_min_eta = 0;
		index_max_eta = order_eta + 1;
	}

	if (index_max_eta > cntl_p_n_eta)
	{
		index_min_eta = cntl_p_n_eta - order_eta - 1;
		index_max_eta = cntl_p_n_eta;
	}

	// zeta
	for (i = 0; i < cntl_p_n_zeta; i++)
	{
		if (knot_vec_zeta[i + 1] >= zeta)
		{
			index_min_zeta = i - order_zeta;
			index_max_zeta = i + 1;
			break;
		}
	}
	if (index_min_zeta < 0)
	{
		index_min_zeta = 0;
		index_max_zeta = order_zeta + 1;
	}

	if (index_max_zeta > cntl_p_n_zeta)
	{
		index_min_zeta = cntl_p_n_zeta - order_zeta - 1;
		index_max_zeta = cntl_p_n_zeta;
	}

	for (i = index_min_xi; i < index_max_xi; i++)
	{
		rBasisFunc(knot_vec_xi, i, order_xi, xi, &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j < index_max_eta; j++)
		{
			rBasisFunc(knot_vec_eta, j, order_eta, eta, &temp_output_eta, &temp_d_output_eta);
			for (k = index_min_zeta; k < index_max_zeta; k++)
			{
				rBasisFunc(knot_vec_zeta, k, order_zeta, zeta, &temp_output_zeta, &temp_d_output_zeta);

				temp_index = i + (j * cntl_p_n_xi) + (k * cntl_p_n_xi * cntl_p_n_eta);

				temp1 = temp_output_xi   * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp2 = temp_d_output_xi * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp3 = temp_output_xi   * temp_d_output_eta * temp_output_zeta   * weight[temp_index];
				temp4 = temp_output_xi   * temp_output_eta   * temp_d_output_zeta * weight[temp_index];

				molecule_x += temp1 * cntl_px[temp_index];
				molecule_y += temp1 * cntl_py[temp_index];
				molecule_z += temp1 * cntl_pz[temp_index];
				denominator += temp1;

				dxi_molecule_x += temp2 * cntl_px[temp_index];
				dxi_molecule_y += temp2 * cntl_py[temp_index];
				dxi_molecule_z += temp2 * cntl_pz[temp_index];
				dxi_denominator += temp2;

				data_molecule_x += temp3 * cntl_px[temp_index];
				data_molecule_y += temp3 * cntl_py[temp_index];
				data_molecule_z += temp3 * cntl_pz[temp_index];
				data_denominator += temp3;

				dzeta_molecule_x += temp4 * cntl_px[temp_index];
				dzeta_molecule_y += temp4 * cntl_py[temp_index];
				dzeta_molecule_z += temp4 * cntl_pz[temp_index];
				dzeta_denominator += temp4;
			}
		}
	}

	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;
	(*output_z) = molecule_z / denominator;
	(*output_dxi_x)   = weight[temp_index] * (dxi_molecule_x   * denominator - molecule_x * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_y)   = weight[temp_index] * (dxi_molecule_y   * denominator - molecule_y * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_z)   = weight[temp_index] * (dxi_molecule_z   * denominator - molecule_z * dxi_denominator)   / (denominator * denominator);
	(*output_data_x)  = weight[temp_index] * (data_molecule_x  * denominator - molecule_x * data_denominator)  / (denominator * denominator);
	(*output_data_y)  = weight[temp_index] * (data_molecule_y  * denominator - molecule_y * data_denominator)  / (denominator * denominator);
	(*output_data_z)  = weight[temp_index] * (data_molecule_z  * denominator - molecule_z * data_denominator)  / (denominator * denominator);
	(*output_dzeta_x) = weight[temp_index] * (dzeta_molecule_x * denominator - molecule_x * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_y) = weight[temp_index] * (dzeta_molecule_y * denominator - molecule_y * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_z) = weight[temp_index] * (dzeta_molecule_z * denominator - molecule_z * dzeta_denominator) / (denominator * denominator);

	return denominator;
}


double lrrNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta,
                       double *cntl_px, double *cntl_py, double *cntl_pz,
                       int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta,
                       double *weight, int order_xi, int order_eta, int order_zeta,
                       double xi, double eta, double zeta,
                       double *output_x, double *output_y, double *output_z,
                       double *output_dxi_x, double *output_data_x, double *output_dzeta_x,
                       double *output_dxi_y, double *output_data_y, double *output_dzeta_y,
                       double *output_dxi_z, double *output_data_z, double *output_dzeta_z)
{
	int i, j, k, temp_index = 0;
	double temp1, temp2, temp3, temp4;
	double molecule_x, molecule_y, molecule_z;
	double dxi_molecule_x, dxi_molecule_y, dxi_molecule_z;
	double data_molecule_x, data_molecule_y, data_molecule_z;
	double dzeta_molecule_x, dzeta_molecule_y, dzeta_molecule_z;
	double denominator, dxi_denominator, data_denominator, dzeta_denominator;
	double temp_output_xi, temp_output_eta, temp_output_zeta;
	double temp_d_output_xi, temp_d_output_eta, temp_d_output_zeta;

	molecule_x = 0.0;
	molecule_y = 0.0;
	molecule_z = 0.0;
	denominator = 0.0;

	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_molecule_z = 0.0;
	dxi_denominator = 0.0;

	data_molecule_x = 0.0;
	data_molecule_y = 0.0;
	data_molecule_z = 0.0;
	data_denominator = 0.0;

	dzeta_molecule_x = 0.0;
	dzeta_molecule_y = 0.0;
	dzeta_molecule_z = 0.0;
	dzeta_denominator = 0.0;

	int index_min_xi   = 0;
	int index_max_xi   = cntl_p_n_xi - 1;
	int index_min_eta  = 0;
	int index_max_eta  = cntl_p_n_eta - 1;
	int index_min_zeta = 0;
	int index_max_zeta = cntl_p_n_zeta - 1;

	// xi
	for (i = 0; i < cntl_p_n_xi; i++)
	{
		if (knot_vec_xi[i + 1] >= xi)
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}

	if (index_min_xi < 0)
	{
		index_min_xi = 0;
		index_max_xi = order_xi + 1;
	}

	if (index_max_xi > cntl_p_n_xi)
	{
		index_min_xi = cntl_p_n_xi - order_xi - 1;
		index_max_xi = cntl_p_n_xi;
	}

	// eta
	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if (knot_vec_eta[i + 1] >= eta)
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0)
	{
		index_min_eta = 0;
		index_max_eta = order_eta + 1;
	}

	if (index_max_eta > cntl_p_n_eta)
	{
		index_min_eta = cntl_p_n_eta - order_eta - 1;
		index_max_eta = cntl_p_n_eta;
	}

	// zeta
	for (i = 0; i < cntl_p_n_zeta; i++)
	{
		if (knot_vec_zeta[i + 1] >= zeta)
		{
			index_min_zeta = i - order_zeta;
			index_max_zeta = i + 1;
			break;
		}
	}
	if (index_min_zeta < 0)
	{
		index_min_zeta = 0;
		index_max_zeta = order_zeta + 1;
	}

	if (index_max_zeta > cntl_p_n_zeta)
	{
		index_min_zeta = cntl_p_n_zeta - order_zeta - 1;
		index_max_zeta = cntl_p_n_zeta;
	}

	for (i = index_min_xi; i < index_max_xi; i++)
	{
		lBasisFunc(knot_vec_xi, i, cntl_p_n_xi, order_xi, xi, &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j < index_max_eta; j++)
		{
			rBasisFunc(knot_vec_eta, j, order_eta, eta, &temp_output_eta, &temp_d_output_eta);
			for (k = index_min_zeta; k < index_max_zeta; k++)
			{
				rBasisFunc(knot_vec_zeta, k, order_zeta, zeta, &temp_output_zeta, &temp_d_output_zeta);

				temp_index = i + j * cntl_p_n_xi + k * cntl_p_n_xi * cntl_p_n_eta;

				temp1 = temp_output_xi   * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp2 = temp_d_output_xi * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp3 = temp_output_xi   * temp_d_output_eta * temp_output_zeta   * weight[temp_index];
				temp4 = temp_output_xi   * temp_output_eta   * temp_d_output_zeta * weight[temp_index];

				molecule_x += temp1 * cntl_px[temp_index];
				molecule_y += temp1 * cntl_py[temp_index];
				molecule_z += temp1 * cntl_pz[temp_index];
				denominator += temp1;

				dxi_molecule_x += temp2 * cntl_px[temp_index];
				dxi_molecule_y += temp2 * cntl_py[temp_index];
				dxi_molecule_z += temp2 * cntl_pz[temp_index];
				dxi_denominator += temp2;

				data_molecule_x += temp3 * cntl_px[temp_index];
				data_molecule_y += temp3 * cntl_py[temp_index];
				data_molecule_z += temp3 * cntl_pz[temp_index];
				data_denominator += temp3;

				dzeta_molecule_x += temp4 * cntl_px[temp_index];
				dzeta_molecule_y += temp4 * cntl_py[temp_index];
				dzeta_molecule_z += temp4 * cntl_pz[temp_index];
				dzeta_denominator += temp4;
			}
		}
	}

	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;
	(*output_z) = molecule_z / denominator;
	(*output_dxi_x)   = weight[temp_index] * (dxi_molecule_x   * denominator - molecule_x * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_y)   = weight[temp_index] * (dxi_molecule_y   * denominator - molecule_y * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_z)   = weight[temp_index] * (dxi_molecule_z   * denominator - molecule_z * dxi_denominator)   / (denominator * denominator);
	(*output_data_x)  = weight[temp_index] * (data_molecule_x  * denominator - molecule_x * data_denominator)  / (denominator * denominator);
	(*output_data_y)  = weight[temp_index] * (data_molecule_y  * denominator - molecule_y * data_denominator)  / (denominator * denominator);
	(*output_data_z)  = weight[temp_index] * (data_molecule_z  * denominator - molecule_z * data_denominator)  / (denominator * denominator);
	(*output_dzeta_x) = weight[temp_index] * (dzeta_molecule_x * denominator - molecule_x * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_y) = weight[temp_index] * (dzeta_molecule_y * denominator - molecule_y * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_z) = weight[temp_index] * (dzeta_molecule_z * denominator - molecule_z * dzeta_denominator) / (denominator * denominator);

	return denominator;
}


double rlrNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta,
                       double *cntl_px, double *cntl_py, double *cntl_pz,
                       int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta,
                       double *weight, int order_xi, int order_eta, int order_zeta,
                       double xi, double eta, double zeta,
                       double *output_x, double *output_y, double *output_z,
                       double *output_dxi_x, double *output_data_x, double *output_dzeta_x,
                       double *output_dxi_y, double *output_data_y, double *output_dzeta_y,
                       double *output_dxi_z, double *output_data_z, double *output_dzeta_z)
{
	int i, j, k, temp_index = 0;
	double temp1, temp2, temp3, temp4;
	double molecule_x, molecule_y, molecule_z;
	double dxi_molecule_x, dxi_molecule_y, dxi_molecule_z;
	double data_molecule_x, data_molecule_y, data_molecule_z;
	double dzeta_molecule_x, dzeta_molecule_y, dzeta_molecule_z;
	double denominator, dxi_denominator, data_denominator, dzeta_denominator;
	double temp_output_xi, temp_output_eta, temp_output_zeta;
	double temp_d_output_xi, temp_d_output_eta, temp_d_output_zeta;

	molecule_x = 0.0;
	molecule_y = 0.0;
	molecule_z = 0.0;
	denominator = 0.0;

	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_molecule_z = 0.0;
	dxi_denominator = 0.0;

	data_molecule_x = 0.0;
	data_molecule_y = 0.0;
	data_molecule_z = 0.0;
	data_denominator = 0.0;

	dzeta_molecule_x = 0.0;
	dzeta_molecule_y = 0.0;
	dzeta_molecule_z = 0.0;
	dzeta_denominator = 0.0;

	int index_min_xi   = 0;
	int index_max_xi   = cntl_p_n_xi - 1;
	int index_min_eta  = 0;
	int index_max_eta  = cntl_p_n_eta - 1;
	int index_min_zeta = 0;
	int index_max_zeta = cntl_p_n_zeta - 1;

	// xi
	for (i = 0; i < cntl_p_n_xi; i++)
	{
		if (knot_vec_xi[i + 1] >= xi)
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}

	if (index_min_xi < 0)
	{
		index_min_xi = 0;
		index_max_xi = order_xi + 1;
	}

	if (index_max_xi > cntl_p_n_xi)
	{
		index_min_xi = cntl_p_n_xi - order_xi - 1;
		index_max_xi = cntl_p_n_xi;
	}

	// eta
	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if (knot_vec_eta[i + 1] >= eta)
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0)
	{
		index_min_eta = 0;
		index_max_eta = order_eta + 1;
	}

	if (index_max_eta > cntl_p_n_eta)
	{
		index_min_eta = cntl_p_n_eta - order_eta - 1;
		index_max_eta = cntl_p_n_eta;
	}

	// zeta
	for (i = 0; i < cntl_p_n_zeta; i++)
	{
		if (knot_vec_zeta[i + 1] >= zeta)
		{
			index_min_zeta = i - order_zeta;
			index_max_zeta = i + 1;
			break;
		}
	}
	if (index_min_zeta < 0)
	{
		index_min_zeta = 0;
		index_max_zeta = order_zeta + 1;
	}

	if (index_max_zeta > cntl_p_n_zeta)
	{
		index_min_zeta = cntl_p_n_zeta - order_zeta - 1;
		index_max_zeta = cntl_p_n_zeta;
	}

	for (i = index_min_xi; i < index_max_xi; i++)
	{
		rBasisFunc(knot_vec_xi, i, order_xi, xi, &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j < index_max_eta; j++)
		{
			lBasisFunc(knot_vec_eta, j, cntl_p_n_eta, order_eta, eta, &temp_output_eta, &temp_d_output_eta);
			for (k = index_min_zeta; k < index_max_zeta; k++)
			{
				rBasisFunc(knot_vec_zeta, k, order_zeta, zeta, &temp_output_zeta, &temp_d_output_zeta);

				temp_index = i + j * cntl_p_n_xi + k * cntl_p_n_xi * cntl_p_n_eta;

				temp1 = temp_output_xi   * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp2 = temp_d_output_xi * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp3 = temp_output_xi   * temp_d_output_eta * temp_output_zeta   * weight[temp_index];
				temp4 = temp_output_xi   * temp_output_eta   * temp_d_output_zeta * weight[temp_index];

				molecule_x += temp1 * cntl_px[temp_index];
				molecule_y += temp1 * cntl_py[temp_index];
				molecule_z += temp1 * cntl_pz[temp_index];
				denominator += temp1;

				dxi_molecule_x += temp2 * cntl_px[temp_index];
				dxi_molecule_y += temp2 * cntl_py[temp_index];
				dxi_molecule_z += temp2 * cntl_pz[temp_index];
				dxi_denominator += temp2;

				data_molecule_x += temp3 * cntl_px[temp_index];
				data_molecule_y += temp3 * cntl_py[temp_index];
				data_molecule_z += temp3 * cntl_pz[temp_index];
				data_denominator += temp3;

				dzeta_molecule_x += temp4 * cntl_px[temp_index];
				dzeta_molecule_y += temp4 * cntl_py[temp_index];
				dzeta_molecule_z += temp4 * cntl_pz[temp_index];
				dzeta_denominator += temp4;
			}
		}
	}

	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;
	(*output_z) = molecule_z / denominator;
	(*output_dxi_x)   = weight[temp_index] * (dxi_molecule_x   * denominator - molecule_x * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_y)   = weight[temp_index] * (dxi_molecule_y   * denominator - molecule_y * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_z)   = weight[temp_index] * (dxi_molecule_z   * denominator - molecule_z * dxi_denominator)   / (denominator * denominator);
	(*output_data_x)  = weight[temp_index] * (data_molecule_x  * denominator - molecule_x * data_denominator)  / (denominator * denominator);
	(*output_data_y)  = weight[temp_index] * (data_molecule_y  * denominator - molecule_y * data_denominator)  / (denominator * denominator);
	(*output_data_z)  = weight[temp_index] * (data_molecule_z  * denominator - molecule_z * data_denominator)  / (denominator * denominator);
	(*output_dzeta_x) = weight[temp_index] * (dzeta_molecule_x * denominator - molecule_x * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_y) = weight[temp_index] * (dzeta_molecule_y * denominator - molecule_y * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_z) = weight[temp_index] * (dzeta_molecule_z * denominator - molecule_z * dzeta_denominator) / (denominator * denominator);

	return denominator;
}


double rrlNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta,
                       double *cntl_px, double *cntl_py, double *cntl_pz,
                       int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta,
                       double *weight, int order_xi, int order_eta, int order_zeta,
                       double xi, double eta, double zeta,
                       double *output_x, double *output_y, double *output_z,
                       double *output_dxi_x, double *output_data_x, double *output_dzeta_x,
                       double *output_dxi_y, double *output_data_y, double *output_dzeta_y,
                       double *output_dxi_z, double *output_data_z, double *output_dzeta_z)
{
	int i, j, k, temp_index = 0;
	double temp1, temp2, temp3, temp4;
	double molecule_x, molecule_y, molecule_z;
	double dxi_molecule_x, dxi_molecule_y, dxi_molecule_z;
	double data_molecule_x, data_molecule_y, data_molecule_z;
	double dzeta_molecule_x, dzeta_molecule_y, dzeta_molecule_z;
	double denominator, dxi_denominator, data_denominator, dzeta_denominator;
	double temp_output_xi, temp_output_eta, temp_output_zeta;
	double temp_d_output_xi, temp_d_output_eta, temp_d_output_zeta;

	molecule_x = 0.0;
	molecule_y = 0.0;
	molecule_z = 0.0;
	denominator = 0.0;

	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_molecule_z = 0.0;
	dxi_denominator = 0.0;

	data_molecule_x = 0.0;
	data_molecule_y = 0.0;
	data_molecule_z = 0.0;
	data_denominator = 0.0;

	dzeta_molecule_x = 0.0;
	dzeta_molecule_y = 0.0;
	dzeta_molecule_z = 0.0;
	dzeta_denominator = 0.0;

	int index_min_xi   = 0;
	int index_max_xi   = cntl_p_n_xi - 1;
	int index_min_eta  = 0;
	int index_max_eta  = cntl_p_n_eta - 1;
	int index_min_zeta = 0;
	int index_max_zeta = cntl_p_n_zeta - 1;

	// xi
	for (i = 0; i < cntl_p_n_xi; i++)
	{
		if (knot_vec_xi[i + 1] >= xi)
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}

	if (index_min_xi < 0)
	{
		index_min_xi = 0;
		index_max_xi = order_xi + 1;
	}

	if (index_max_xi > cntl_p_n_xi)
	{
		index_min_xi = cntl_p_n_xi - order_xi - 1;
		index_max_xi = cntl_p_n_xi;
	}

	// eta
	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if (knot_vec_eta[i + 1] >= eta)
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0)
	{
		index_min_eta = 0;
		index_max_eta = order_eta + 1;
	}

	if (index_max_eta > cntl_p_n_eta)
	{
		index_min_eta = cntl_p_n_eta - order_eta - 1;
		index_max_eta = cntl_p_n_eta;
	}

	// zeta
	for (i = 0; i < cntl_p_n_zeta; i++)
	{
		if (knot_vec_zeta[i + 1] >= zeta)
		{
			index_min_zeta = i - order_zeta;
			index_max_zeta = i + 1;
			break;
		}
	}
	if (index_min_zeta < 0)
	{
		index_min_zeta = 0;
		index_max_zeta = order_zeta + 1;
	}

	if (index_max_zeta > cntl_p_n_zeta)
	{
		index_min_zeta = cntl_p_n_zeta - order_zeta - 1;
		index_max_zeta = cntl_p_n_zeta;
	}

	for (i = index_min_xi; i < index_max_xi; i++)
	{
		rBasisFunc(knot_vec_xi, i, order_xi, xi, &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j < index_max_eta; j++)
		{
			rBasisFunc(knot_vec_eta, j, order_eta, eta, &temp_output_eta, &temp_d_output_eta);
			for (k = index_min_zeta; k < index_max_zeta; k++)
			{
				lBasisFunc(knot_vec_zeta, k, cntl_p_n_zeta, order_zeta, zeta, &temp_output_zeta, &temp_d_output_zeta);

				temp_index = i + j * cntl_p_n_xi + k * cntl_p_n_xi * cntl_p_n_eta;

				temp1 = temp_output_xi   * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp2 = temp_d_output_xi * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp3 = temp_output_xi   * temp_d_output_eta * temp_output_zeta   * weight[temp_index];
				temp4 = temp_output_xi   * temp_output_eta   * temp_d_output_zeta * weight[temp_index];

				molecule_x += temp1 * cntl_px[temp_index];
				molecule_y += temp1 * cntl_py[temp_index];
				molecule_z += temp1 * cntl_pz[temp_index];
				denominator += temp1;

				dxi_molecule_x += temp2 * cntl_px[temp_index];
				dxi_molecule_y += temp2 * cntl_py[temp_index];
				dxi_molecule_z += temp2 * cntl_pz[temp_index];
				dxi_denominator += temp2;

				data_molecule_x += temp3 * cntl_px[temp_index];
				data_molecule_y += temp3 * cntl_py[temp_index];
				data_molecule_z += temp3 * cntl_pz[temp_index];
				data_denominator += temp3;

				dzeta_molecule_x += temp4 * cntl_px[temp_index];
				dzeta_molecule_y += temp4 * cntl_py[temp_index];
				dzeta_molecule_z += temp4 * cntl_pz[temp_index];
				dzeta_denominator += temp4;
			}
		}
	}

	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;
	(*output_z) = molecule_z / denominator;
	(*output_dxi_x)   = weight[temp_index] * (dxi_molecule_x   * denominator - molecule_x * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_y)   = weight[temp_index] * (dxi_molecule_y   * denominator - molecule_y * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_z)   = weight[temp_index] * (dxi_molecule_z   * denominator - molecule_z * dxi_denominator)   / (denominator * denominator);
	(*output_data_x)  = weight[temp_index] * (data_molecule_x  * denominator - molecule_x * data_denominator)  / (denominator * denominator);
	(*output_data_y)  = weight[temp_index] * (data_molecule_y  * denominator - molecule_y * data_denominator)  / (denominator * denominator);
	(*output_data_z)  = weight[temp_index] * (data_molecule_z  * denominator - molecule_z * data_denominator)  / (denominator * denominator);
	(*output_dzeta_x) = weight[temp_index] * (dzeta_molecule_x * denominator - molecule_x * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_y) = weight[temp_index] * (dzeta_molecule_y * denominator - molecule_y * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_z) = weight[temp_index] * (dzeta_molecule_z * denominator - molecule_z * dzeta_denominator) / (denominator * denominator);

	return denominator;
}


double llrNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta,
                       double *cntl_px, double *cntl_py, double *cntl_pz,
                       int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta,
                       double *weight, int order_xi, int order_eta, int order_zeta,
                       double xi, double eta, double zeta,
                       double *output_x, double *output_y, double *output_z,
                       double *output_dxi_x, double *output_data_x, double *output_dzeta_x,
                       double *output_dxi_y, double *output_data_y, double *output_dzeta_y,
                       double *output_dxi_z, double *output_data_z, double *output_dzeta_z)
{
	int i, j, k, temp_index = 0;
	double temp1, temp2, temp3, temp4;
	double molecule_x, molecule_y, molecule_z;
	double dxi_molecule_x, dxi_molecule_y, dxi_molecule_z;
	double data_molecule_x, data_molecule_y, data_molecule_z;
	double dzeta_molecule_x, dzeta_molecule_y, dzeta_molecule_z;
	double denominator, dxi_denominator, data_denominator, dzeta_denominator;
	double temp_output_xi, temp_output_eta, temp_output_zeta;
	double temp_d_output_xi, temp_d_output_eta, temp_d_output_zeta;

	molecule_x = 0.0;
	molecule_y = 0.0;
	molecule_z = 0.0;
	denominator = 0.0;

	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_molecule_z = 0.0;
	dxi_denominator = 0.0;

	data_molecule_x = 0.0;
	data_molecule_y = 0.0;
	data_molecule_z = 0.0;
	data_denominator = 0.0;

	dzeta_molecule_x = 0.0;
	dzeta_molecule_y = 0.0;
	dzeta_molecule_z = 0.0;
	dzeta_denominator = 0.0;

	int index_min_xi   = 0;
	int index_max_xi   = cntl_p_n_xi - 1;
	int index_min_eta  = 0;
	int index_max_eta  = cntl_p_n_eta - 1;
	int index_min_zeta = 0;
	int index_max_zeta = cntl_p_n_zeta - 1;

	// xi
	for (i = 0; i < cntl_p_n_xi; i++)
	{
		if (knot_vec_xi[i + 1] >= xi)
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}

	if (index_min_xi < 0)
	{
		index_min_xi = 0;
		index_max_xi = order_xi + 1;
	}

	if (index_max_xi > cntl_p_n_xi)
	{
		index_min_xi = cntl_p_n_xi - order_xi - 1;
		index_max_xi = cntl_p_n_xi;
	}

	// eta
	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if (knot_vec_eta[i + 1] >= eta)
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0)
	{
		index_min_eta = 0;
		index_max_eta = order_eta + 1;
	}

	if (index_max_eta > cntl_p_n_eta)
	{
		index_min_eta = cntl_p_n_eta - order_eta - 1;
		index_max_eta = cntl_p_n_eta;
	}

	// zeta
	for (i = 0; i < cntl_p_n_zeta; i++)
	{
		if (knot_vec_zeta[i + 1] >= zeta)
		{
			index_min_zeta = i - order_zeta;
			index_max_zeta = i + 1;
			break;
		}
	}
	if (index_min_zeta < 0)
	{
		index_min_zeta = 0;
		index_max_zeta = order_zeta + 1;
	}

	if (index_max_zeta > cntl_p_n_zeta)
	{
		index_min_zeta = cntl_p_n_zeta - order_zeta - 1;
		index_max_zeta = cntl_p_n_zeta;
	}

	for (i = index_min_xi; i < index_max_xi; i++)
	{
		lBasisFunc(knot_vec_xi, i, cntl_p_n_xi, order_xi, xi, &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j < index_max_eta; j++)
		{
			lBasisFunc(knot_vec_eta, j, cntl_p_n_eta, order_eta, eta, &temp_output_eta, &temp_d_output_eta);
			for (k = index_min_zeta; k < index_max_zeta; k++)
			{
				rBasisFunc(knot_vec_zeta, k, order_zeta, zeta, &temp_output_zeta, &temp_d_output_zeta);

				temp_index = i + j * cntl_p_n_xi + k * cntl_p_n_xi * cntl_p_n_eta;

				temp1 = temp_output_xi   * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp2 = temp_d_output_xi * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp3 = temp_output_xi   * temp_d_output_eta * temp_output_zeta   * weight[temp_index];
				temp4 = temp_output_xi   * temp_output_eta   * temp_d_output_zeta * weight[temp_index];

				molecule_x += temp1 * cntl_px[temp_index];
				molecule_y += temp1 * cntl_py[temp_index];
				molecule_z += temp1 * cntl_pz[temp_index];
				denominator += temp1;

				dxi_molecule_x += temp2 * cntl_px[temp_index];
				dxi_molecule_y += temp2 * cntl_py[temp_index];
				dxi_molecule_z += temp2 * cntl_pz[temp_index];
				dxi_denominator += temp2;

				data_molecule_x += temp3 * cntl_px[temp_index];
				data_molecule_y += temp3 * cntl_py[temp_index];
				data_molecule_z += temp3 * cntl_pz[temp_index];
				data_denominator += temp3;

				dzeta_molecule_x += temp4 * cntl_px[temp_index];
				dzeta_molecule_y += temp4 * cntl_py[temp_index];
				dzeta_molecule_z += temp4 * cntl_pz[temp_index];
				dzeta_denominator += temp4;
			}
		}
	}

	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;
	(*output_z) = molecule_z / denominator;
	(*output_dxi_x)   = weight[temp_index] * (dxi_molecule_x   * denominator - molecule_x * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_y)   = weight[temp_index] * (dxi_molecule_y   * denominator - molecule_y * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_z)   = weight[temp_index] * (dxi_molecule_z   * denominator - molecule_z * dxi_denominator)   / (denominator * denominator);
	(*output_data_x)  = weight[temp_index] * (data_molecule_x  * denominator - molecule_x * data_denominator)  / (denominator * denominator);
	(*output_data_y)  = weight[temp_index] * (data_molecule_y  * denominator - molecule_y * data_denominator)  / (denominator * denominator);
	(*output_data_z)  = weight[temp_index] * (data_molecule_z  * denominator - molecule_z * data_denominator)  / (denominator * denominator);
	(*output_dzeta_x) = weight[temp_index] * (dzeta_molecule_x * denominator - molecule_x * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_y) = weight[temp_index] * (dzeta_molecule_y * denominator - molecule_y * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_z) = weight[temp_index] * (dzeta_molecule_z * denominator - molecule_z * dzeta_denominator) / (denominator * denominator);

	return denominator;
}


double lrlNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta,
                       double *cntl_px, double *cntl_py, double *cntl_pz,
                       int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta,
                       double *weight, int order_xi, int order_eta, int order_zeta,
                       double xi, double eta, double zeta,
                       double *output_x, double *output_y, double *output_z,
                       double *output_dxi_x, double *output_data_x, double *output_dzeta_x,
                       double *output_dxi_y, double *output_data_y, double *output_dzeta_y,
                       double *output_dxi_z, double *output_data_z, double *output_dzeta_z)
{
	int i, j, k, temp_index = 0;
	double temp1, temp2, temp3, temp4;
	double molecule_x, molecule_y, molecule_z;
	double dxi_molecule_x, dxi_molecule_y, dxi_molecule_z;
	double data_molecule_x, data_molecule_y, data_molecule_z;
	double dzeta_molecule_x, dzeta_molecule_y, dzeta_molecule_z;
	double denominator, dxi_denominator, data_denominator, dzeta_denominator;
	double temp_output_xi, temp_output_eta, temp_output_zeta;
	double temp_d_output_xi, temp_d_output_eta, temp_d_output_zeta;

	molecule_x = 0.0;
	molecule_y = 0.0;
	molecule_z = 0.0;
	denominator = 0.0;

	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_molecule_z = 0.0;
	dxi_denominator = 0.0;

	data_molecule_x = 0.0;
	data_molecule_y = 0.0;
	data_molecule_z = 0.0;
	data_denominator = 0.0;

	dzeta_molecule_x = 0.0;
	dzeta_molecule_y = 0.0;
	dzeta_molecule_z = 0.0;
	dzeta_denominator = 0.0;

	int index_min_xi   = 0;
	int index_max_xi   = cntl_p_n_xi - 1;
	int index_min_eta  = 0;
	int index_max_eta  = cntl_p_n_eta - 1;
	int index_min_zeta = 0;
	int index_max_zeta = cntl_p_n_zeta - 1;

	// xi
	for (i = 0; i < cntl_p_n_xi; i++)
	{
		if (knot_vec_xi[i + 1] >= xi)
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}

	if (index_min_xi < 0)
	{
		index_min_xi = 0;
		index_max_xi = order_xi + 1;
	}

	if (index_max_xi > cntl_p_n_xi)
	{
		index_min_xi = cntl_p_n_xi - order_xi - 1;
		index_max_xi = cntl_p_n_xi;
	}

	// eta
	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if (knot_vec_eta[i + 1] >= eta)
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0)
	{
		index_min_eta = 0;
		index_max_eta = order_eta + 1;
	}

	if (index_max_eta > cntl_p_n_eta)
	{
		index_min_eta = cntl_p_n_eta - order_eta - 1;
		index_max_eta = cntl_p_n_eta;
	}

	// zeta
	for (i = 0; i < cntl_p_n_zeta; i++)
	{
		if (knot_vec_zeta[i + 1] >= zeta)
		{
			index_min_zeta = i - order_zeta;
			index_max_zeta = i + 1;
			break;
		}
	}
	if (index_min_zeta < 0)
	{
		index_min_zeta = 0;
		index_max_zeta = order_zeta + 1;
	}

	if (index_max_zeta > cntl_p_n_zeta)
	{
		index_min_zeta = cntl_p_n_zeta - order_zeta - 1;
		index_max_zeta = cntl_p_n_zeta;
	}

	for (i = index_min_xi; i < index_max_xi; i++)
	{
		lBasisFunc(knot_vec_xi, i, cntl_p_n_xi, order_xi, xi, &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j < index_max_eta; j++)
		{
			rBasisFunc(knot_vec_eta, j, order_eta, eta, &temp_output_eta, &temp_d_output_eta);
			for (k = index_min_zeta; k < index_max_zeta; k++)
			{
				lBasisFunc(knot_vec_zeta, k, cntl_p_n_zeta, order_zeta, zeta, &temp_output_zeta, &temp_d_output_zeta);

				temp_index = i + j * cntl_p_n_xi + k * cntl_p_n_xi * cntl_p_n_eta;

				temp1 = temp_output_xi   * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp2 = temp_d_output_xi * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp3 = temp_output_xi   * temp_d_output_eta * temp_output_zeta   * weight[temp_index];
				temp4 = temp_output_xi   * temp_output_eta   * temp_d_output_zeta * weight[temp_index];

				molecule_x += temp1 * cntl_px[temp_index];
				molecule_y += temp1 * cntl_py[temp_index];
				molecule_z += temp1 * cntl_pz[temp_index];
				denominator += temp1;

				dxi_molecule_x += temp2 * cntl_px[temp_index];
				dxi_molecule_y += temp2 * cntl_py[temp_index];
				dxi_molecule_z += temp2 * cntl_pz[temp_index];
				dxi_denominator += temp2;

				data_molecule_x += temp3 * cntl_px[temp_index];
				data_molecule_y += temp3 * cntl_py[temp_index];
				data_molecule_z += temp3 * cntl_pz[temp_index];
				data_denominator += temp3;

				dzeta_molecule_x += temp4 * cntl_px[temp_index];
				dzeta_molecule_y += temp4 * cntl_py[temp_index];
				dzeta_molecule_z += temp4 * cntl_pz[temp_index];
				dzeta_denominator += temp4;
			}
		}
	}

	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;
	(*output_z) = molecule_z / denominator;
	(*output_dxi_x)   = weight[temp_index] * (dxi_molecule_x   * denominator - molecule_x * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_y)   = weight[temp_index] * (dxi_molecule_y   * denominator - molecule_y * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_z)   = weight[temp_index] * (dxi_molecule_z   * denominator - molecule_z * dxi_denominator)   / (denominator * denominator);
	(*output_data_x)  = weight[temp_index] * (data_molecule_x  * denominator - molecule_x * data_denominator)  / (denominator * denominator);
	(*output_data_y)  = weight[temp_index] * (data_molecule_y  * denominator - molecule_y * data_denominator)  / (denominator * denominator);
	(*output_data_z)  = weight[temp_index] * (data_molecule_z  * denominator - molecule_z * data_denominator)  / (denominator * denominator);
	(*output_dzeta_x) = weight[temp_index] * (dzeta_molecule_x * denominator - molecule_x * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_y) = weight[temp_index] * (dzeta_molecule_y * denominator - molecule_y * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_z) = weight[temp_index] * (dzeta_molecule_z * denominator - molecule_z * dzeta_denominator) / (denominator * denominator);

	return denominator;
}


double rllNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta,
                       double *cntl_px, double *cntl_py, double *cntl_pz,
                       int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta,
                       double *weight, int order_xi, int order_eta, int order_zeta,
                       double xi, double eta, double zeta,
                       double *output_x, double *output_y, double *output_z,
                       double *output_dxi_x, double *output_data_x, double *output_dzeta_x,
                       double *output_dxi_y, double *output_data_y, double *output_dzeta_y,
                       double *output_dxi_z, double *output_data_z, double *output_dzeta_z)
{
	int i, j, k, temp_index = 0;
	double temp1, temp2, temp3, temp4;
	double molecule_x, molecule_y, molecule_z;
	double dxi_molecule_x, dxi_molecule_y, dxi_molecule_z;
	double data_molecule_x, data_molecule_y, data_molecule_z;
	double dzeta_molecule_x, dzeta_molecule_y, dzeta_molecule_z;
	double denominator, dxi_denominator, data_denominator, dzeta_denominator;
	double temp_output_xi, temp_output_eta, temp_output_zeta;
	double temp_d_output_xi, temp_d_output_eta, temp_d_output_zeta;

	molecule_x = 0.0;
	molecule_y = 0.0;
	molecule_z = 0.0;
	denominator = 0.0;

	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_molecule_z = 0.0;
	dxi_denominator = 0.0;

	data_molecule_x = 0.0;
	data_molecule_y = 0.0;
	data_molecule_z = 0.0;
	data_denominator = 0.0;

	dzeta_molecule_x = 0.0;
	dzeta_molecule_y = 0.0;
	dzeta_molecule_z = 0.0;
	dzeta_denominator = 0.0;

	int index_min_xi   = 0;
	int index_max_xi   = cntl_p_n_xi - 1;
	int index_min_eta  = 0;
	int index_max_eta  = cntl_p_n_eta - 1;
	int index_min_zeta = 0;
	int index_max_zeta = cntl_p_n_zeta - 1;

	// xi
	for (i = 0; i < cntl_p_n_xi; i++)
	{
		if (knot_vec_xi[i + 1] >= xi)
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}

	if (index_min_xi < 0)
	{
		index_min_xi = 0;
		index_max_xi = order_xi + 1;
	}

	if (index_max_xi > cntl_p_n_xi)
	{
		index_min_xi = cntl_p_n_xi - order_xi - 1;
		index_max_xi = cntl_p_n_xi;
	}

	// eta
	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if (knot_vec_eta[i + 1] >= eta)
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0)
	{
		index_min_eta = 0;
		index_max_eta = order_eta + 1;
	}

	if (index_max_eta > cntl_p_n_eta)
	{
		index_min_eta = cntl_p_n_eta - order_eta - 1;
		index_max_eta = cntl_p_n_eta;
	}

	// zeta
	for (i = 0; i < cntl_p_n_zeta; i++)
	{
		if (knot_vec_zeta[i + 1] >= zeta)
		{
			index_min_zeta = i - order_zeta;
			index_max_zeta = i + 1;
			break;
		}
	}
	if (index_min_zeta < 0)
	{
		index_min_zeta = 0;
		index_max_zeta = order_zeta + 1;
	}

	if (index_max_zeta > cntl_p_n_zeta)
	{
		index_min_zeta = cntl_p_n_zeta - order_zeta - 1;
		index_max_zeta = cntl_p_n_zeta;
	}

	for (i = index_min_xi; i < index_max_xi; i++)
	{
		rBasisFunc(knot_vec_xi, i, order_xi, xi, &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j < index_max_eta; j++)
		{
			lBasisFunc(knot_vec_eta, j, cntl_p_n_eta, order_eta, eta, &temp_output_eta, &temp_d_output_eta);
			for (k = index_min_zeta; k < index_max_zeta; k++)
			{
				lBasisFunc(knot_vec_zeta, k, cntl_p_n_zeta, order_zeta, zeta, &temp_output_zeta, &temp_d_output_zeta);

				temp_index = i + j * cntl_p_n_xi + k * cntl_p_n_xi * cntl_p_n_eta;

				temp1 = temp_output_xi   * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp2 = temp_d_output_xi * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp3 = temp_output_xi   * temp_d_output_eta * temp_output_zeta   * weight[temp_index];
				temp4 = temp_output_xi   * temp_output_eta   * temp_d_output_zeta * weight[temp_index];

				molecule_x += temp1 * cntl_px[temp_index];
				molecule_y += temp1 * cntl_py[temp_index];
				molecule_z += temp1 * cntl_pz[temp_index];
				denominator += temp1;

				dxi_molecule_x += temp2 * cntl_px[temp_index];
				dxi_molecule_y += temp2 * cntl_py[temp_index];
				dxi_molecule_z += temp2 * cntl_pz[temp_index];
				dxi_denominator += temp2;

				data_molecule_x += temp3 * cntl_px[temp_index];
				data_molecule_y += temp3 * cntl_py[temp_index];
				data_molecule_z += temp3 * cntl_pz[temp_index];
				data_denominator += temp3;

				dzeta_molecule_x += temp4 * cntl_px[temp_index];
				dzeta_molecule_y += temp4 * cntl_py[temp_index];
				dzeta_molecule_z += temp4 * cntl_pz[temp_index];
				dzeta_denominator += temp4;
			}
		}
	}

	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;
	(*output_z) = molecule_z / denominator;
	(*output_dxi_x)   = weight[temp_index] * (dxi_molecule_x   * denominator - molecule_x * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_y)   = weight[temp_index] * (dxi_molecule_y   * denominator - molecule_y * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_z)   = weight[temp_index] * (dxi_molecule_z   * denominator - molecule_z * dxi_denominator)   / (denominator * denominator);
	(*output_data_x)  = weight[temp_index] * (data_molecule_x  * denominator - molecule_x * data_denominator)  / (denominator * denominator);
	(*output_data_y)  = weight[temp_index] * (data_molecule_y  * denominator - molecule_y * data_denominator)  / (denominator * denominator);
	(*output_data_z)  = weight[temp_index] * (data_molecule_z  * denominator - molecule_z * data_denominator)  / (denominator * denominator);
	(*output_dzeta_x) = weight[temp_index] * (dzeta_molecule_x * denominator - molecule_x * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_y) = weight[temp_index] * (dzeta_molecule_y * denominator - molecule_y * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_z) = weight[temp_index] * (dzeta_molecule_z * denominator - molecule_z * dzeta_denominator) / (denominator * denominator);

	return denominator;
}


double lllNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta,
                       double *cntl_px, double *cntl_py, double *cntl_pz,
                       int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta,
                       double *weight, int order_xi, int order_eta, int order_zeta,
                       double xi, double eta, double zeta,
                       double *output_x, double *output_y, double *output_z,
                       double *output_dxi_x, double *output_data_x, double *output_dzeta_x,
                       double *output_dxi_y, double *output_data_y, double *output_dzeta_y,
                       double *output_dxi_z, double *output_data_z, double *output_dzeta_z)
{
	int i, j, k, temp_index = 0;
	double temp1, temp2, temp3, temp4;
	double molecule_x, molecule_y, molecule_z;
	double dxi_molecule_x, dxi_molecule_y, dxi_molecule_z;
	double data_molecule_x, data_molecule_y, data_molecule_z;
	double dzeta_molecule_x, dzeta_molecule_y, dzeta_molecule_z;
	double denominator, dxi_denominator, data_denominator, dzeta_denominator;
	double temp_output_xi, temp_output_eta, temp_output_zeta;
	double temp_d_output_xi, temp_d_output_eta, temp_d_output_zeta;

	molecule_x = 0.0;
	molecule_y = 0.0;
	molecule_z = 0.0;
	denominator = 0.0;

	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_molecule_z = 0.0;
	dxi_denominator = 0.0;

	data_molecule_x = 0.0;
	data_molecule_y = 0.0;
	data_molecule_z = 0.0;
	data_denominator = 0.0;

	dzeta_molecule_x = 0.0;
	dzeta_molecule_y = 0.0;
	dzeta_molecule_z = 0.0;
	dzeta_denominator = 0.0;

	int index_min_xi   = 0;
	int index_max_xi   = cntl_p_n_xi - 1;
	int index_min_eta  = 0;
	int index_max_eta  = cntl_p_n_eta - 1;
	int index_min_zeta = 0;
	int index_max_zeta = cntl_p_n_zeta - 1;

	// xi
	for (i = 0; i < cntl_p_n_xi; i++)
	{
		if (knot_vec_xi[i + 1] >= xi)
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}

	if (index_min_xi < 0)
	{
		index_min_xi = 0;
		index_max_xi = order_xi + 1;
	}

	if (index_max_xi > cntl_p_n_xi)
	{
		index_min_xi = cntl_p_n_xi - order_xi - 1;
		index_max_xi = cntl_p_n_xi;
	}

	// eta
	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if (knot_vec_eta[i + 1] >= eta)
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0)
	{
		index_min_eta = 0;
		index_max_eta = order_eta + 1;
	}

	if (index_max_eta > cntl_p_n_eta)
	{
		index_min_eta = cntl_p_n_eta - order_eta - 1;
		index_max_eta = cntl_p_n_eta;
	}

	// zeta
	for (i = 0; i < cntl_p_n_zeta; i++)
	{
		if (knot_vec_zeta[i + 1] >= zeta)
		{
			index_min_zeta = i - order_zeta;
			index_max_zeta = i + 1;
			break;
		}
	}
	if (index_min_zeta < 0)
	{
		index_min_zeta = 0;
		index_max_zeta = order_zeta + 1;
	}

	if (index_max_zeta > cntl_p_n_zeta)
	{
		index_min_zeta = cntl_p_n_zeta - order_zeta - 1;
		index_max_zeta = cntl_p_n_zeta;
	}

	for (i = index_min_xi; i < index_max_xi; i++)
	{
		lBasisFunc(knot_vec_xi, i, cntl_p_n_xi, order_xi, xi, &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j < index_max_eta; j++)
		{
			lBasisFunc(knot_vec_eta, j, cntl_p_n_eta, order_eta, eta, &temp_output_eta, &temp_d_output_eta);
			for (k = index_min_zeta; k < index_max_zeta; k++)
			{
				lBasisFunc(knot_vec_zeta, k, cntl_p_n_zeta, order_zeta, zeta, &temp_output_zeta, &temp_d_output_zeta);

				temp_index = i + j * cntl_p_n_xi + k * cntl_p_n_xi * cntl_p_n_eta;

				temp1 = temp_output_xi   * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp2 = temp_d_output_xi * temp_output_eta   * temp_output_zeta   * weight[temp_index];
				temp3 = temp_output_xi   * temp_d_output_eta * temp_output_zeta   * weight[temp_index];
				temp4 = temp_output_xi   * temp_output_eta   * temp_d_output_zeta * weight[temp_index];

				molecule_x += temp1 * cntl_px[temp_index];
				molecule_y += temp1 * cntl_py[temp_index];
				molecule_z += temp1 * cntl_pz[temp_index];
				denominator += temp1;

				dxi_molecule_x += temp2 * cntl_px[temp_index];
				dxi_molecule_y += temp2 * cntl_py[temp_index];
				dxi_molecule_z += temp2 * cntl_pz[temp_index];
				dxi_denominator += temp2;

				data_molecule_x += temp3 * cntl_px[temp_index];
				data_molecule_y += temp3 * cntl_py[temp_index];
				data_molecule_z += temp3 * cntl_pz[temp_index];
				data_denominator += temp3;

				dzeta_molecule_x += temp4 * cntl_px[temp_index];
				dzeta_molecule_y += temp4 * cntl_py[temp_index];
				dzeta_molecule_z += temp4 * cntl_pz[temp_index];
				dzeta_denominator += temp4;
			}
		}
	}

	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;
	(*output_z) = molecule_z / denominator;
	(*output_dxi_x)   = weight[temp_index] * (dxi_molecule_x   * denominator - molecule_x * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_y)   = weight[temp_index] * (dxi_molecule_y   * denominator - molecule_y * dxi_denominator)   / (denominator * denominator);
	(*output_dxi_z)   = weight[temp_index] * (dxi_molecule_z   * denominator - molecule_z * dxi_denominator)   / (denominator * denominator);
	(*output_data_x)  = weight[temp_index] * (data_molecule_x  * denominator - molecule_x * data_denominator)  / (denominator * denominator);
	(*output_data_y)  = weight[temp_index] * (data_molecule_y  * denominator - molecule_y * data_denominator)  / (denominator * denominator);
	(*output_data_z)  = weight[temp_index] * (data_molecule_z  * denominator - molecule_z * data_denominator)  / (denominator * denominator);
	(*output_dzeta_x) = weight[temp_index] * (dzeta_molecule_x * denominator - molecule_x * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_y) = weight[temp_index] * (dzeta_molecule_y * denominator - molecule_y * dzeta_denominator) / (denominator * denominator);
	(*output_dzeta_z) = weight[temp_index] * (dzeta_molecule_z * denominator - molecule_z * dzeta_denominator) / (denominator * denominator);

	return denominator;
}


// 算出したローカルパッチ各要素の頂点の物理座標のグローバルパッチでの(xi,eta)算出
int Calc_xi_eta(double px, double py,
				double *input_knot_vec_xi, double *input_knot_vec_eta,
				int cntl_p_n_xi, int cntl_p_n_eta, int order_xi, int order_eta,
				double *output_xi, double *output_eta, information *info)
{
	double temp_xi, temp_eta;
	double temp_x, temp_y;
	double temp_matrix[2][2];
	double temp_dxi, temp_data;
	double temp_tol_x, temp_tol_y;

	(*output_xi) = 0;
	(*output_eta) = 0;

	int i;
	// int repeat = 1000;
	// double tol = 10e-8;
	int repeat = 10;
	double tol = 10e-14;

	double *temp_Position_Knots_xi = (double *)malloc(sizeof(double) * info->No_knot[0 * info->DIMENSION + 0]);
	double *temp_Position_Knots_eta = (double *)malloc(sizeof(double) * info->No_knot[0 * info->DIMENSION + 1]);
	for (i = 0; i < info->No_knot[0 * info->DIMENSION + 0]; i++)
	{
		temp_Position_Knots_xi[i] = info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 0] + i];
	}
	for (i = 0; i < info->No_knot[0 * info->DIMENSION + 1]; i++)
	{
		temp_Position_Knots_eta[i] = info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 1] + i];
	}

	// 初期値の設定
	temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;

	for (i = 0; i < repeat; i++)
	{
		rNURBS_surface(temp_Position_Knots_xi, temp_Position_Knots_eta,
					   info->Control_Coord_x, info->Control_Coord_y,
					   info->No_Control_point[0 * info->DIMENSION + 0], info->No_Control_point[0 * info->DIMENSION + 1],
					   info->Control_Weight, info->Order[0 * info->DIMENSION + 0], info->Order[0 * info->DIMENSION + 1],
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
		temp_data = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_data;
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
					   info->Control_Coord_x, info->Control_Coord_y,
					   info->No_Control_point[0 * info->DIMENSION + 0], info->No_Control_point[0 * info->DIMENSION + 1],
					   info->Control_Weight, info->Order[0 * info->DIMENSION + 0], info->Order[0 * info->DIMENSION + 1],
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
		temp_data = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_data;
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
						info->Control_Coord_x, info->Control_Coord_y,
						info->No_Control_point[0 * info->DIMENSION + 0], info->No_Control_point[0 * info->DIMENSION + 1],
						info->Control_Weight, info->Order[0 * info->DIMENSION + 0], info->Order[0 * info->DIMENSION + 1],
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
		temp_data = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_data;
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
						info->Control_Coord_x, info->Control_Coord_y,
						info->No_Control_point[0 * info->DIMENSION + 0], info->No_Control_point[0 * info->DIMENSION + 1],
						info->Control_Weight, info->Order[0 * info->DIMENSION + 0], info->Order[0 * info->DIMENSION + 1],
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
		temp_data = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_data;
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


int Calc_xi_eta_zeta(double px, double py, double pz,
				     double *input_knot_vec_xi, double *input_knot_vec_eta, double *input_knot_vec_zeta,
				     int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta,
                     int order_xi, int order_eta, int order_zeta,
				     double *output_xi, double *output_eta, double *output_zeta, information *info)
{
	double temp_xi, temp_eta, temp_zeta;
	double temp_x, temp_y, temp_z;
	double temp_matrix[3][3];
	double temp_dxi, temp_data, temp_dzeta;
    double temp_tol_x, temp_tol_y, temp_tol_z;

	(*output_xi)   = 0.0;
	(*output_eta)  = 0.0;
	(*output_zeta) = 0.0;

	int i, j;
	// int repeat  = 5;
	// int repeat2 = 5;
	// double tol = 10e-22;
	int repeat  = 10;
	int repeat2 = 5;
	double tol = 10e-14;
	double coef = 0.0;

	double *temp_Position_Knots_xi = (double *)malloc(sizeof(double) * info->No_knot[0 * info->DIMENSION + 0]);
	double *temp_Position_Knots_eta = (double *)malloc(sizeof(double) * info->No_knot[0 * info->DIMENSION + 1]);
	double *temp_Position_Knots_zeta = (double *)malloc(sizeof(double) * info->No_knot[0 * info->DIMENSION + 2]);
	for (i = 0; i < info->No_knot[0 * info->DIMENSION + 0]; i++)
	{
		temp_Position_Knots_xi[i] = info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 0] + i];
	}
	for (i = 0; i < info->No_knot[0 * info->DIMENSION + 1]; i++)
	{
		temp_Position_Knots_eta[i] = info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 1] + i];
	}
	for (i = 0; i < info->No_knot[0 * info->DIMENSION + 2]; i++)
	{
		temp_Position_Knots_eta[i] = info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 2] + i];
	}

	for (j = 0; j < repeat2; j++)
	{
		// 初期値の設定
		if (j == 0)
		{
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= 0.5;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= 0.5;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= 0.5;
		}
		else
		{
			coef = j / repeat2;
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= coef;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= coef;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= coef;
		}

		for (i = 0; i < repeat; i++)
		{
			rrrNURBS_volume(temp_Position_Knots_xi, temp_Position_Knots_eta, temp_Position_Knots_zeta,
							info->Control_Coord_x, info->Control_Coord_y, info->Control_Coord_z,
							info->No_Control_point[0 * info->DIMENSION + 0], info->No_Control_point[0 * info->DIMENSION + 1], info->No_Control_point[0 * info->DIMENSION + 2],
							info->Control_Weight, info->Order[0 * info->DIMENSION + 0], info->Order[0 * info->DIMENSION + 1], info->Order[0 * info->DIMENSION + 2],
							temp_xi, temp_eta, temp_zeta,
							&temp_x, &temp_y, &temp_z,
							&temp_matrix[0][0], &temp_matrix[0][1], &temp_matrix[0][2],
							&temp_matrix[1][0], &temp_matrix[1][1], &temp_matrix[1][2],
							&temp_matrix[2][0], &temp_matrix[2][1], &temp_matrix[2][2]);

			temp_tol_x = px - temp_x;
			temp_tol_x *= temp_tol_x;
			temp_tol_y = py - temp_y;
			temp_tol_y *= temp_tol_y;
			temp_tol_z = pz - temp_z;
			temp_tol_z *= temp_tol_z;

			// 収束した場合
			if (temp_tol_x < tol && temp_tol_y < tol && temp_tol_z < tol && i != 0)
			{
				(*output_xi)   = temp_xi;
				(*output_eta)  = temp_eta;
				(*output_zeta) = temp_zeta;
				return i;
			}

			InverseMatrix_3x3(temp_matrix);

			temp_dxi   = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y) + temp_matrix[0][2] * (pz - temp_z);
			temp_data  = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y) + temp_matrix[1][2] * (pz - temp_z);
			temp_dzeta = temp_matrix[2][0] * (px - temp_x) + temp_matrix[2][1] * (py - temp_y) + temp_matrix[2][2] * (pz - temp_z);
			temp_xi   = temp_xi + temp_dxi;
			temp_eta  = temp_eta + temp_data;
			temp_zeta = temp_zeta + temp_dzeta;
			if (temp_xi < input_knot_vec_xi[0])
				temp_xi = input_knot_vec_xi[0];
			if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
				temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
			if (temp_eta < input_knot_vec_eta[0])
				temp_eta = input_knot_vec_eta[0];
			if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
				temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];
			if (temp_zeta < input_knot_vec_zeta[0])
				temp_zeta = input_knot_vec_zeta[0];
			if (temp_zeta > input_knot_vec_zeta[cntl_p_n_zeta + order_zeta])
				temp_zeta = input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
		}
	}

	for (j = 0; j < repeat2; j++)
	{
		// 初期値の設定
		if (j == 0)
		{
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= 0.5;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= 0.5;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= 0.5;
		}
		else
		{
			coef = j / repeat2;
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= coef;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= coef;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= coef;
		}

		for (i = 0; i < repeat; i++)
		{
			lrrNURBS_volume(temp_Position_Knots_xi, temp_Position_Knots_eta, temp_Position_Knots_zeta,
							info->Control_Coord_x, info->Control_Coord_y, info->Control_Coord_z,
							info->No_Control_point[0 * info->DIMENSION + 0], info->No_Control_point[0 * info->DIMENSION + 1], info->No_Control_point[0 * info->DIMENSION + 2],
							info->Control_Weight, info->Order[0 * info->DIMENSION + 0], info->Order[0 * info->DIMENSION + 1], info->Order[0 * info->DIMENSION + 2],
							temp_xi, temp_eta, temp_zeta,
							&temp_x, &temp_y, &temp_z,
							&temp_matrix[0][0], &temp_matrix[0][1], &temp_matrix[0][2],
							&temp_matrix[1][0], &temp_matrix[1][1], &temp_matrix[1][2],
							&temp_matrix[2][0], &temp_matrix[2][1], &temp_matrix[2][2]);

			temp_tol_x = px - temp_x;
			temp_tol_x *= temp_tol_x;
			temp_tol_y = py - temp_y;
			temp_tol_y *= temp_tol_y;
			temp_tol_z = pz - temp_z;
			temp_tol_z *= temp_tol_z;

			// 収束した場合
			if (temp_tol_x < tol && temp_tol_y < tol && temp_tol_z < tol && i != 0)
			{
				(*output_xi)   = temp_xi;
				(*output_eta)  = temp_eta;
				(*output_zeta) = temp_zeta;
				return i;
			}

			InverseMatrix_3x3(temp_matrix);

			temp_dxi   = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y) + temp_matrix[0][2] * (pz - temp_z);
			temp_data  = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y) + temp_matrix[1][2] * (pz - temp_z);
			temp_dzeta = temp_matrix[2][0] * (px - temp_x) + temp_matrix[2][1] * (py - temp_y) + temp_matrix[2][2] * (pz - temp_z);
			temp_xi   = temp_xi + temp_dxi;
			temp_eta  = temp_eta + temp_data;
			temp_zeta = temp_zeta + temp_dzeta;
			if (temp_xi < input_knot_vec_xi[0])
				temp_xi = input_knot_vec_xi[0];
			if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
				temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
			if (temp_eta < input_knot_vec_eta[0])
				temp_eta = input_knot_vec_eta[0];
			if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
				temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];
			if (temp_zeta < input_knot_vec_zeta[0])
				temp_zeta = input_knot_vec_zeta[0];
			if (temp_zeta > input_knot_vec_zeta[cntl_p_n_zeta + order_zeta])
				temp_zeta = input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
		}
	}

	for (j = 0; j < repeat2; j++)
	{
		// 初期値の設定
		if (j == 0)
		{
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= 0.5;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= 0.5;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= 0.5;
		}
		else
		{
			coef = j / repeat2;
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= coef;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= coef;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= coef;
		}

		for (i = 0; i < repeat; i++)
		{
			rlrNURBS_volume(temp_Position_Knots_xi, temp_Position_Knots_eta, temp_Position_Knots_zeta,
							info->Control_Coord_x, info->Control_Coord_y, info->Control_Coord_z,
							info->No_Control_point[0 * info->DIMENSION + 0], info->No_Control_point[0 * info->DIMENSION + 1], info->No_Control_point[0 * info->DIMENSION + 2],
							info->Control_Weight, info->Order[0 * info->DIMENSION + 0], info->Order[0 * info->DIMENSION + 1], info->Order[0 * info->DIMENSION + 2],
							temp_xi, temp_eta, temp_zeta,
							&temp_x, &temp_y, &temp_z,
							&temp_matrix[0][0], &temp_matrix[0][1], &temp_matrix[0][2],
							&temp_matrix[1][0], &temp_matrix[1][1], &temp_matrix[1][2],
							&temp_matrix[2][0], &temp_matrix[2][1], &temp_matrix[2][2]);

			temp_tol_x = px - temp_x;
			temp_tol_x *= temp_tol_x;
			temp_tol_y = py - temp_y;
			temp_tol_y *= temp_tol_y;
			temp_tol_z = pz - temp_z;
			temp_tol_z *= temp_tol_z;

			// 収束した場合
			if (temp_tol_x < tol && temp_tol_y < tol && temp_tol_z < tol && i != 0)
			{
				(*output_xi)   = temp_xi;
				(*output_eta)  = temp_eta;
				(*output_zeta) = temp_zeta;
				return i;
			}

			InverseMatrix_3x3(temp_matrix);

			temp_dxi   = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y) + temp_matrix[0][2] * (pz - temp_z);
			temp_data  = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y) + temp_matrix[1][2] * (pz - temp_z);
			temp_dzeta = temp_matrix[2][0] * (px - temp_x) + temp_matrix[2][1] * (py - temp_y) + temp_matrix[2][2] * (pz - temp_z);
			temp_xi   = temp_xi + temp_dxi;
			temp_eta  = temp_eta + temp_data;
			temp_zeta = temp_zeta + temp_dzeta;
			if (temp_xi < input_knot_vec_xi[0])
				temp_xi = input_knot_vec_xi[0];
			if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
				temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
			if (temp_eta < input_knot_vec_eta[0])
				temp_eta = input_knot_vec_eta[0];
			if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
				temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];
			if (temp_zeta < input_knot_vec_zeta[0])
				temp_zeta = input_knot_vec_zeta[0];
			if (temp_zeta > input_knot_vec_zeta[cntl_p_n_zeta + order_zeta])
				temp_zeta = input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
		}
	}

	for (j = 0; j < repeat2; j++)
	{
		// 初期値の設定
		if (j == 0)
		{
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= 0.5;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= 0.5;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= 0.5;
		}
		else
		{
			coef = j / repeat2;
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= coef;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= coef;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= coef;
		}

		for (i = 0; i < repeat; i++)
		{
			rrlNURBS_volume(temp_Position_Knots_xi, temp_Position_Knots_eta, temp_Position_Knots_zeta,
							info->Control_Coord_x, info->Control_Coord_y, info->Control_Coord_z,
							info->No_Control_point[0 * info->DIMENSION + 0], info->No_Control_point[0 * info->DIMENSION + 1], info->No_Control_point[0 * info->DIMENSION + 2],
							info->Control_Weight, info->Order[0 * info->DIMENSION + 0], info->Order[0 * info->DIMENSION + 1], info->Order[0 * info->DIMENSION + 2],
							temp_xi, temp_eta, temp_zeta,
							&temp_x, &temp_y, &temp_z,
							&temp_matrix[0][0], &temp_matrix[0][1], &temp_matrix[0][2],
							&temp_matrix[1][0], &temp_matrix[1][1], &temp_matrix[1][2],
							&temp_matrix[2][0], &temp_matrix[2][1], &temp_matrix[2][2]);

			temp_tol_x = px - temp_x;
			temp_tol_x *= temp_tol_x;
			temp_tol_y = py - temp_y;
			temp_tol_y *= temp_tol_y;
			temp_tol_z = pz - temp_z;
			temp_tol_z *= temp_tol_z;

			// 収束した場合
			if (temp_tol_x < tol && temp_tol_y < tol && temp_tol_z < tol && i != 0)
			{
				(*output_xi)   = temp_xi;
				(*output_eta)  = temp_eta;
				(*output_zeta) = temp_zeta;
				return i;
			}

			InverseMatrix_3x3(temp_matrix);

			temp_dxi   = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y) + temp_matrix[0][2] * (pz - temp_z);
			temp_data  = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y) + temp_matrix[1][2] * (pz - temp_z);
			temp_dzeta = temp_matrix[2][0] * (px - temp_x) + temp_matrix[2][1] * (py - temp_y) + temp_matrix[2][2] * (pz - temp_z);
			temp_xi   = temp_xi + temp_dxi;
			temp_eta  = temp_eta + temp_data;
			temp_zeta = temp_zeta + temp_dzeta;
			if (temp_xi < input_knot_vec_xi[0])
				temp_xi = input_knot_vec_xi[0];
			if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
				temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
			if (temp_eta < input_knot_vec_eta[0])
				temp_eta = input_knot_vec_eta[0];
			if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
				temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];
			if (temp_zeta < input_knot_vec_zeta[0])
				temp_zeta = input_knot_vec_zeta[0];
			if (temp_zeta > input_knot_vec_zeta[cntl_p_n_zeta + order_zeta])
				temp_zeta = input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
		}
	}

	for (j = 0; j < repeat2; j++)
	{
		// 初期値の設定
		if (j == 0)
		{
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= 0.5;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= 0.5;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= 0.5;
		}
		else
		{
			coef = j / repeat2;
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= coef;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= coef;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= coef;
		}

		for (i = 0; i < repeat; i++)
		{
			llrNURBS_volume(temp_Position_Knots_xi, temp_Position_Knots_eta, temp_Position_Knots_zeta,
							info->Control_Coord_x, info->Control_Coord_y, info->Control_Coord_z,
							info->No_Control_point[0 * info->DIMENSION + 0], info->No_Control_point[0 * info->DIMENSION + 1], info->No_Control_point[0 * info->DIMENSION + 2],
							info->Control_Weight, info->Order[0 * info->DIMENSION + 0], info->Order[0 * info->DIMENSION + 1], info->Order[0 * info->DIMENSION + 2],
							temp_xi, temp_eta, temp_zeta,
							&temp_x, &temp_y, &temp_z,
							&temp_matrix[0][0], &temp_matrix[0][1], &temp_matrix[0][2],
							&temp_matrix[1][0], &temp_matrix[1][1], &temp_matrix[1][2],
							&temp_matrix[2][0], &temp_matrix[2][1], &temp_matrix[2][2]);

			temp_tol_x = px - temp_x;
			temp_tol_x *= temp_tol_x;
			temp_tol_y = py - temp_y;
			temp_tol_y *= temp_tol_y;
			temp_tol_z = pz - temp_z;
			temp_tol_z *= temp_tol_z;

			// 収束した場合
			if (temp_tol_x < tol && temp_tol_y < tol && temp_tol_z < tol && i != 0)
			{
				(*output_xi)   = temp_xi;
				(*output_eta)  = temp_eta;
				(*output_zeta) = temp_zeta;
				return i;
			}

			InverseMatrix_3x3(temp_matrix);

			temp_dxi   = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y) + temp_matrix[0][2] * (pz - temp_z);
			temp_data  = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y) + temp_matrix[1][2] * (pz - temp_z);
			temp_dzeta = temp_matrix[2][0] * (px - temp_x) + temp_matrix[2][1] * (py - temp_y) + temp_matrix[2][2] * (pz - temp_z);
			temp_xi   = temp_xi + temp_dxi;
			temp_eta  = temp_eta + temp_data;
			temp_zeta = temp_zeta + temp_dzeta;
			if (temp_xi < input_knot_vec_xi[0])
				temp_xi = input_knot_vec_xi[0];
			if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
				temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
			if (temp_eta < input_knot_vec_eta[0])
				temp_eta = input_knot_vec_eta[0];
			if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
				temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];
			if (temp_zeta < input_knot_vec_zeta[0])
				temp_zeta = input_knot_vec_zeta[0];
			if (temp_zeta > input_knot_vec_zeta[cntl_p_n_zeta + order_zeta])
				temp_zeta = input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
		}
	}

	for (j = 0; j < repeat2; j++)
	{
		// 初期値の設定
		if (j == 0)
		{
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= 0.5;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= 0.5;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= 0.5;
		}
		else
		{
			coef = j / repeat2;
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= coef;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= coef;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= coef;
		}

		for (i = 0; i < repeat; i++)
		{
			lrlNURBS_volume(temp_Position_Knots_xi, temp_Position_Knots_eta, temp_Position_Knots_zeta,
							info->Control_Coord_x, info->Control_Coord_y, info->Control_Coord_z,
							info->No_Control_point[0 * info->DIMENSION + 0], info->No_Control_point[0 * info->DIMENSION + 1], info->No_Control_point[0 * info->DIMENSION + 2],
							info->Control_Weight, info->Order[0 * info->DIMENSION + 0], info->Order[0 * info->DIMENSION + 1], info->Order[0 * info->DIMENSION + 2],
							temp_xi, temp_eta, temp_zeta,
							&temp_x, &temp_y, &temp_z,
							&temp_matrix[0][0], &temp_matrix[0][1], &temp_matrix[0][2],
							&temp_matrix[1][0], &temp_matrix[1][1], &temp_matrix[1][2],
							&temp_matrix[2][0], &temp_matrix[2][1], &temp_matrix[2][2]);

			temp_tol_x = px - temp_x;
			temp_tol_x *= temp_tol_x;
			temp_tol_y = py - temp_y;
			temp_tol_y *= temp_tol_y;
			temp_tol_z = pz - temp_z;
			temp_tol_z *= temp_tol_z;

			// 収束した場合
			if (temp_tol_x < tol && temp_tol_y < tol && temp_tol_z < tol && i != 0)
			{
				(*output_xi)   = temp_xi;
				(*output_eta)  = temp_eta;
				(*output_zeta) = temp_zeta;
				return i;
			}

			InverseMatrix_3x3(temp_matrix);

			temp_dxi   = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y) + temp_matrix[0][2] * (pz - temp_z);
			temp_data  = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y) + temp_matrix[1][2] * (pz - temp_z);
			temp_dzeta = temp_matrix[2][0] * (px - temp_x) + temp_matrix[2][1] * (py - temp_y) + temp_matrix[2][2] * (pz - temp_z);
			temp_xi   = temp_xi + temp_dxi;
			temp_eta  = temp_eta + temp_data;
			temp_zeta = temp_zeta + temp_dzeta;
			if (temp_xi < input_knot_vec_xi[0])
				temp_xi = input_knot_vec_xi[0];
			if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
				temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
			if (temp_eta < input_knot_vec_eta[0])
				temp_eta = input_knot_vec_eta[0];
			if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
				temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];
			if (temp_zeta < input_knot_vec_zeta[0])
				temp_zeta = input_knot_vec_zeta[0];
			if (temp_zeta > input_knot_vec_zeta[cntl_p_n_zeta + order_zeta])
				temp_zeta = input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
		}
	}

	for (j = 0; j < repeat2; j++)
	{
		// 初期値の設定
		if (j == 0)
		{
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= 0.5;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= 0.5;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= 0.5;
		}
		else
		{
			coef = j / repeat2;
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= coef;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= coef;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= coef;
		}

		for (i = 0; i < repeat; i++)
		{
			rllNURBS_volume(temp_Position_Knots_xi, temp_Position_Knots_eta, temp_Position_Knots_zeta,
							info->Control_Coord_x, info->Control_Coord_y, info->Control_Coord_z,
							info->No_Control_point[0 * info->DIMENSION + 0], info->No_Control_point[0 * info->DIMENSION + 1], info->No_Control_point[0 * info->DIMENSION + 2],
							info->Control_Weight, info->Order[0 * info->DIMENSION + 0], info->Order[0 * info->DIMENSION + 1], info->Order[0 * info->DIMENSION + 2],
							temp_xi, temp_eta, temp_zeta,
							&temp_x, &temp_y, &temp_z,
							&temp_matrix[0][0], &temp_matrix[0][1], &temp_matrix[0][2],
							&temp_matrix[1][0], &temp_matrix[1][1], &temp_matrix[1][2],
							&temp_matrix[2][0], &temp_matrix[2][1], &temp_matrix[2][2]);

			temp_tol_x = px - temp_x;
			temp_tol_x *= temp_tol_x;
			temp_tol_y = py - temp_y;
			temp_tol_y *= temp_tol_y;
			temp_tol_z = pz - temp_z;
			temp_tol_z *= temp_tol_z;

			// 収束した場合
			if (temp_tol_x < tol && temp_tol_y < tol && temp_tol_z < tol && i != 0)
			{
				(*output_xi)   = temp_xi;
				(*output_eta)  = temp_eta;
				(*output_zeta) = temp_zeta;
				return i;
			}

			InverseMatrix_3x3(temp_matrix);

			temp_dxi   = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y) + temp_matrix[0][2] * (pz - temp_z);
			temp_data  = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y) + temp_matrix[1][2] * (pz - temp_z);
			temp_dzeta = temp_matrix[2][0] * (px - temp_x) + temp_matrix[2][1] * (py - temp_y) + temp_matrix[2][2] * (pz - temp_z);
			temp_xi   = temp_xi + temp_dxi;
			temp_eta  = temp_eta + temp_data;
			temp_zeta = temp_zeta + temp_dzeta;
			if (temp_xi < input_knot_vec_xi[0])
				temp_xi = input_knot_vec_xi[0];
			if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
				temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
			if (temp_eta < input_knot_vec_eta[0])
				temp_eta = input_knot_vec_eta[0];
			if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
				temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];
			if (temp_zeta < input_knot_vec_zeta[0])
				temp_zeta = input_knot_vec_zeta[0];
			if (temp_zeta > input_knot_vec_zeta[cntl_p_n_zeta + order_zeta])
				temp_zeta = input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
		}
	}

	for (j = 0; j < repeat2; j++)
	{
		// 初期値の設定
		if (j == 0)
		{
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= 0.5;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= 0.5;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= 0.5;
		}
		else
		{
			coef = j / repeat2;
			temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= coef;
			temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= coef;
			temp_zeta = input_knot_vec_zeta[0] + input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= coef;
		}

		for (i = 0; i < repeat; i++)
		{
			lllNURBS_volume(temp_Position_Knots_xi, temp_Position_Knots_eta, temp_Position_Knots_zeta,
							info->Control_Coord_x, info->Control_Coord_y, info->Control_Coord_z,
							info->No_Control_point[0 * info->DIMENSION + 0], info->No_Control_point[0 * info->DIMENSION + 1], info->No_Control_point[0 * info->DIMENSION + 2],
							info->Control_Weight, info->Order[0 * info->DIMENSION + 0], info->Order[0 * info->DIMENSION + 1], info->Order[0 * info->DIMENSION + 2],
							temp_xi, temp_eta, temp_zeta,
							&temp_x, &temp_y, &temp_z,
							&temp_matrix[0][0], &temp_matrix[0][1], &temp_matrix[0][2],
							&temp_matrix[1][0], &temp_matrix[1][1], &temp_matrix[1][2],
							&temp_matrix[2][0], &temp_matrix[2][1], &temp_matrix[2][2]);

			temp_tol_x = px - temp_x;
			temp_tol_x *= temp_tol_x;
			temp_tol_y = py - temp_y;
			temp_tol_y *= temp_tol_y;
			temp_tol_z = pz - temp_z;
			temp_tol_z *= temp_tol_z;

			// 収束した場合
			if (temp_tol_x < tol && temp_tol_y < tol && temp_tol_z < tol && i != 0)
			{
				(*output_xi)   = temp_xi;
				(*output_eta)  = temp_eta;
				(*output_zeta) = temp_zeta;
				return i;
			}

			InverseMatrix_3x3(temp_matrix);

			temp_dxi   = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y) + temp_matrix[0][2] * (pz - temp_z);
			temp_data  = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y) + temp_matrix[1][2] * (pz - temp_z);
			temp_dzeta = temp_matrix[2][0] * (px - temp_x) + temp_matrix[2][1] * (py - temp_y) + temp_matrix[2][2] * (pz - temp_z);
			temp_xi   = temp_xi + temp_dxi;
			temp_eta  = temp_eta + temp_data;
			temp_zeta = temp_zeta + temp_dzeta;
			if (temp_xi < input_knot_vec_xi[0])
				temp_xi = input_knot_vec_xi[0];
			if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
				temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
			if (temp_eta < input_knot_vec_eta[0])
				temp_eta = input_knot_vec_eta[0];
			if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
				temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];
			if (temp_zeta < input_knot_vec_zeta[0])
				temp_zeta = input_knot_vec_zeta[0];
			if (temp_zeta > input_knot_vec_zeta[cntl_p_n_zeta + order_zeta])
				temp_zeta = input_knot_vec_zeta[cntl_p_n_zeta + order_zeta];
		}
	}

	free(temp_Position_Knots_xi), free(temp_Position_Knots_eta), free(temp_Position_Knots_zeta);

	return 0;
}


// Postprocessing
void Make_Displacement(information *info)
{
	int i, j;

	for (i = 0; i < info->Total_Constraint_to_mesh[Total_mesh]; i++)
    {
		info->Displacement[info->Constraint_Node_Dir[i * 2 + 0] * info->DIMENSION + info->Constraint_Node_Dir[i * 2 + 1]] = info->Value_of_Constraint[i];
    }

	for (i = 0; i < info->Total_Control_Point_to_mesh[Total_mesh]; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			int index = info->Index_Dof[i * info->DIMENSION + j];
			if (index >= 0)
				info->Displacement[i * info->DIMENSION + j] = info->sol_vec[index];
		}
	}
}


void Make_Strain(information *info)
{
	int i, j, k, l;
	int KIEL_SIZE;
	double *temp_disp = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT * info->DIMENSION);

	for (i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
	{
		KIEL_SIZE = info->No_Control_point_ON_ELEMENT[info->Element_patch[i]] * info->DIMENSION;

		for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; j++)
		{
			for (k = 0; k < info->DIMENSION; k++)
			{
				temp_disp[j * info->DIMENSION + k] = info->Displacement[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + j] * info->DIMENSION + k];
			}
		}

		for (j = 0; j < GP_ON_ELEMENT; j++)
		{
			for (k = 0; k < D_MATRIX_SIZE; k++)
			{
				for (l = 0; l < KIEL_SIZE; l++)
				{
					info->Strain_at_GP[i * GP_ON_ELEMENT * N_STRAIN + j * N_STRAIN + k] += info->B_Matrix[i * MAX_POW_NG * D_MATRIX_SIZE * MAX_KIEL_SIZE + j * D_MATRIX_SIZE * MAX_KIEL_SIZE + k * MAX_KIEL_SIZE + l] * temp_disp[l];
				}
			}
		}
	}

	free(temp_disp);
}


void Make_Stress(information *info)
{
	int i, j, k, l;

	for (i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
	{
		for (j = 0; j < GP_ON_ELEMENT; j++)
		{
			for (k = 0; k < D_MATRIX_SIZE; k++)
			{
				for (l = 0; l < D_MATRIX_SIZE; l++)
				{
					info->Stress_at_GP[i * GP_ON_ELEMENT * N_STRESS + j * N_STRESS + k] += info->D[k * D_MATRIX_SIZE + l] * info->Strain_at_GP[i * GP_ON_ELEMENT * N_STRAIN + j * N_STRAIN + l];
				}
			}
		}
	}
}


void Make_Parameter_z(information *info)
{
	int i, j;

	// 平面応力状態
	if (DM == 0)
	{
		for (i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
		{
			for (j = 0; j < GP_ON_ELEMENT; j++)
			{
				info->Strain_at_GP[i * GP_ON_ELEMENT * N_STRAIN + j * N_STRAIN + 3] = - 1.0 * nu / E * (info->Stress_at_GP[i * GP_ON_ELEMENT * N_STRESS + j * N_STRESS + 0] + info->Stress_at_GP[i * GP_ON_ELEMENT * N_STRESS + j * N_STRESS + 1]);
			}
		}
	}
	// 平面ひずみ状態
	else if (DM == 1)
	{
		for (i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
		{
			for (j = 0; j < GP_ON_ELEMENT; j++)
			{
				info->Stress_at_GP[i * GP_ON_ELEMENT * N_STRESS + j * N_STRESS + 3] = E * nu / (1.0 + nu) / (1.0 - 2.0 * nu) * (info->Strain_at_GP[i * GP_ON_ELEMENT * N_STRAIN + j * N_STRAIN + 0] + info->Strain_at_GP[i * GP_ON_ELEMENT * N_STRAIN + j * N_STRAIN + 1]);
			}
		}
	}
}


void Make_ReactionForce(information *info)
{
	int i, j, k, l, m;

	for (i = 0; i < info->Total_Element_to_mesh[Total_mesh]; i++)
	{
		for (j = 0; j < GP_ON_ELEMENT; j++)
		{
			for (k = 0; k < D_MATRIX_SIZE; k++)
			{
				for (l = 0; l < info->DIMENSION; l++)
				{
					for (m = 0; m < info->No_Control_point_ON_ELEMENT[info->Element_patch[i]]; m++)
					{
						info->ReactionForce[info->Controlpoint_of_Element[i * MAX_NO_CP_ON_ELEMENT + m] * info->DIMENSION + l] += info->B_Matrix[i * MAX_POW_NG * D_MATRIX_SIZE * MAX_KIEL_SIZE + j * D_MATRIX_SIZE * MAX_KIEL_SIZE + k * MAX_KIEL_SIZE + m * info->DIMENSION + l] * info->Stress_at_GP[i * GP_ON_ELEMENT * N_STRESS + j * N_STRESS + k] * w[j] * info->Jac[i * MAX_POW_NG + j];
					}
				}
			}
		}
	}
}



// for S-IGA overlay
void S_IGA_overlay(information *info)
{
	int i, j, k;
	int n_patch_glo = info->Total_Patch_to_mesh[1];

	// 一要素の分割数
	division_ele_xi = DIVISION_ELE;
	division_ele_eta = DIVISION_ELE;

	// output
	fp = fopen("view.dat", "w");
	fprintf(fp, "%d\t%d\t%d\n", 1, division_ele_xi, division_ele_eta);
	fclose(fp);

	fp = fopen("overlay_view.dat", "w");
	fprintf(fp, "%d\t%d\t%d\n", 1, division_ele_xi, division_ele_eta);
	fclose(fp);

	fp = fopen("disp.txt", "w");
	fprintf(fp, "patch_n\tx\ty\tdisp_x\tdisp_y\n");
	fclose(fp);

	fp = fopen("stress_y.txt", "w");
	fprintf(fp, "patch_n\tx\ty\tstress_yy\n");
	fclose(fp);

	fp = fopen("stress_y_0.txt", "w");
	fprintf(fp, "patch_n\tx\ty\tstress_yy\n");
	fclose(fp);

	fp = fopen("stress_vm.txt", "w");
	fprintf(fp, "xi\teta\tx\ty\tstress_vm\n");
	fclose(fp);

	fp = fopen("overlay_disp.txt", "w");
	fprintf(fp, "patch_n\tx\ty\tdisp_x\tdisp_y\n");
	fclose(fp);

	fp = fopen("overlay_stress_x.txt", "w");
	fprintf(fp, "x\ty\tstress_xx\n");
	fclose(fp);

	fp = fopen("overlay_stress_y.txt", "w");
	fprintf(fp, "x\ty\tstress_yy\n");
	fclose(fp);

	fp = fopen("overlay_stress_y_0.txt", "w");
	fprintf(fp, "x\ty\tstress_yy\n");
	fclose(fp);

	fp = fopen("overlay_stress_r_theta.txt", "w");
	fprintf(fp, "xi\teta\tx\ty\tstress_r\tstress_theta\n");
	fclose(fp);

	fp = fopen("overlay_stress_vm.txt", "w");
	fprintf(fp, "xi\teta\tx\ty\tstress_vm\n");
	fclose(fp);

	// 重ね合わせ
	for (i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
	{
		fp = fopen("stress_vm.txt", "a");
		fprintf(fp, "\npatch_n;%d\n\n", i);
		fclose(fp);

		// パッチ i での各方向ノットベクトル, Coodinate, displacement
		double *temp_Position_Knots_xi = (double *)malloc(sizeof(double) * info->No_knot[i * info->DIMENSION + 0]);
		double *temp_Position_Knots_eta = (double *)malloc(sizeof(double) * info->No_knot[i * info->DIMENSION + 1]);
		double *temp_Coord_x = (double *)malloc(sizeof(double) * info->No_Control_point_in_patch[i]);
		double *temp_Coord_y = (double *)malloc(sizeof(double) * info->No_Control_point_in_patch[i]);
		double *temp_Coord_w = (double *)malloc(sizeof(double) * info->No_Control_point_in_patch[i]);
		double *temp_disp_x = (double *)malloc(sizeof(double) * info->No_Control_point_in_patch[i]);
		double *temp_disp_y = (double *)malloc(sizeof(double) * info->No_Control_point_in_patch[i]);

		for (j = 0; j < info->No_knot[i * info->DIMENSION + 0]; j++)
		{
			temp_Position_Knots_xi[j] = info->Position_Knots[info->Total_Knot_to_patch_dim[i * info->DIMENSION + 0] + j];
		}
		for (j = 0; j < info->No_knot[i * info->DIMENSION + 1]; j++)
		{
			temp_Position_Knots_eta[j] = info->Position_Knots[info->Total_Knot_to_patch_dim[i * info->DIMENSION + 1] + j];
		}

		for (j = 0; j < info->No_Control_point_in_patch[i]; j++)
		{
			temp_Coord_x[j] = info->Control_Coord_x[info->Patch_Control_point[info->Total_Control_Point_to_patch[i] + j]];
			temp_Coord_y[j] = info->Control_Coord_y[info->Patch_Control_point[info->Total_Control_Point_to_patch[i] + j]];
			temp_Coord_w[j] = info->Control_Weight[info->Patch_Control_point[info->Total_Control_Point_to_patch[i] + j]];
			temp_disp_x[j] = info->Displacement[info->Patch_Control_point[info->Total_Control_Point_to_patch[i] + j] * info->DIMENSION + 0];
			temp_disp_y[j] = info->Displacement[info->Patch_Control_point[info->Total_Control_Point_to_patch[i] + j] * info->DIMENSION + 1];
		}

		Calculation(info->Order[i * info->DIMENSION + 0], info->Order[i * info->DIMENSION + 1], info->No_knot[i * info->DIMENSION + 0], info->No_knot[i * info->DIMENSION + 1], info->No_Control_point[i * info->DIMENSION + 0], info->No_Control_point[i * info->DIMENSION + 1],
					temp_Position_Knots_xi, temp_Position_Knots_eta, temp_Coord_x, temp_Coord_y, temp_Coord_w,
					temp_disp_x, temp_disp_y, i, info);

		if (i >= n_patch_glo) // ローカル上のパッチに対しては重合計算行う
		{
			for (j = 0; j < n_patch_glo; j++)
			{
				// パッチ j での各方向ノットベクトル, Coodinate, displacement
				double *temp_Position_Knots_xi_glo = (double *)malloc(sizeof(double) * info->No_knot[j * info->DIMENSION + 0]);
				double *temp_Position_Knots_eta_glo = (double *)malloc(sizeof(double) * info->No_knot[j * info->DIMENSION + 1]);
				double *temp_Coord_x_glo = (double *)malloc(sizeof(double) * info->No_Control_point_in_patch[j]);
				double *temp_Coord_y_glo = (double *)malloc(sizeof(double) * info->No_Control_point_in_patch[j]);
				double *temp_Coord_w_glo = (double *)malloc(sizeof(double) * info->No_Control_point_in_patch[j]);
				double *temp_disp_x_glo = (double *)malloc(sizeof(double) * info->No_Control_point_in_patch[j]);
				double *temp_disp_y_glo = (double *)malloc(sizeof(double) * info->No_Control_point_in_patch[j]);

				for (k = 0; k < info->No_knot[j * info->DIMENSION + 0]; k++)
				{
					temp_Position_Knots_xi_glo[k] = info->Position_Knots[info->Total_Knot_to_patch_dim[j * info->DIMENSION + 0] + k];
				}
				for (k = 0; k < info->No_knot[j * info->DIMENSION + 1]; k++)
				{
					temp_Position_Knots_eta_glo[k] = info->Position_Knots[info->Total_Knot_to_patch_dim[j * info->DIMENSION + 1] + k];
				}

				for (k = 0; k < info->No_Control_point_in_patch[j]; k++)
				{
					temp_Coord_x_glo[k] = info->Control_Coord_x[info->Patch_Control_point[info->Total_Control_Point_to_patch[j] + k]];
					temp_Coord_y_glo[k] = info->Control_Coord_y[info->Patch_Control_point[info->Total_Control_Point_to_patch[j] + k]];
					temp_Coord_w_glo[k] = info->Control_Weight[info->Patch_Control_point[info->Total_Control_Point_to_patch[j] + k]];
					temp_disp_x_glo[k] = info->Displacement[(info->Patch_Control_point[info->Total_Control_Point_to_patch[j] + k]) * info->DIMENSION + 0];
					temp_disp_y_glo[k] = info->Displacement[(info->Patch_Control_point[info->Total_Control_Point_to_patch[j] + k]) * info->DIMENSION + 1];
				}

				Calculation_overlay(info->Order[i * info->DIMENSION + 0], info->Order[i * info->DIMENSION + 1],
									info->No_knot[i * info->DIMENSION + 0], info->No_knot[i * info->DIMENSION + 1],
									info->No_Control_point[i * info->DIMENSION + 0], info->No_Control_point[i * info->DIMENSION + 1],
									temp_Position_Knots_xi, temp_Position_Knots_eta,
									temp_Coord_x, temp_Coord_y, temp_Coord_w,
									info->Order[j * info->DIMENSION + 0], info->Order[j * info->DIMENSION + 1],
									info->No_Control_point[j * info->DIMENSION + 0], info->No_Control_point[j * info->DIMENSION + 1],
									temp_Position_Knots_xi_glo, temp_Position_Knots_eta_glo,	
									temp_Coord_x_glo, temp_Coord_y_glo,
									temp_disp_x_glo, temp_disp_y_glo, temp_Coord_w_glo,
									i, info);

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
				 double *disp_cntl_px, double *disp_cntl_py, int patch_num, information *info)
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

	double *dxi_x = (double *)malloc(sizeof(double) * info->Total_Element_to_mesh[Total_mesh] * DIVISION_ELEMENT);		// dxi_x[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]
	double *dxi_y = (double *)malloc(sizeof(double) * info->Total_Element_to_mesh[Total_mesh] * DIVISION_ELEMENT);		// dxi_y[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]
	double *data_x = (double *)malloc(sizeof(double) * info->Total_Element_to_mesh[Total_mesh] * DIVISION_ELEMENT);		// data_x[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]
	double *data_y = (double *)malloc(sizeof(double) * info->Total_Element_to_mesh[Total_mesh] * DIVISION_ELEMENT);		// data_y[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]
	double *dxi_disp_x = (double *)malloc(sizeof(double) * info->Total_Element_to_mesh[Total_mesh] * DIVISION_ELEMENT);	// dxi_disp_x[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]
	double *dxi_disp_y = (double *)malloc(sizeof(double) * info->Total_Element_to_mesh[Total_mesh] * DIVISION_ELEMENT);	// dxi_disp_y[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]
	double *data_disp_x = (double *)malloc(sizeof(double) * info->Total_Element_to_mesh[Total_mesh] * DIVISION_ELEMENT);	// data_disp_x[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]
	double *data_disp_y = (double *)malloc(sizeof(double) * info->Total_Element_to_mesh[Total_mesh] * DIVISION_ELEMENT);	// data_disp_y[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]

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
						   &info->coord_x[i * division_n_eta + j], &info->coord_y[i * division_n_eta + j],
						   &dxi_x[ep], &data_x[ep],
						   &dxi_y[ep], &data_y[ep]);
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
						   &info->disp_x[i * division_n_eta + j], &info->disp_y[i * division_n_eta + j],
						   &dxi_disp_x[ep], &data_disp_x[ep],
						   &dxi_disp_y[ep], &data_disp_y[ep]);
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
							   &info->coord_x[i * division_n_eta + j], &info->coord_y[i * division_n_eta + j],
							   &dxi_x[ep], &data_x[ep],
							   &dxi_y[ep], &data_y[ep]);
				rNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
							   disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
							   weight, order_xi, order_eta,
							   calc_xi[i], calc_eta[j],
							   &info->disp_x[i * division_n_eta + j], &info->disp_y[i * division_n_eta + j],
							   &dxi_disp_x[ep], &data_disp_x[ep],
							   &dxi_disp_y[ep], &data_disp_y[ep]);
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
							   &info->coord_x[i * division_n_eta + j], &info->coord_y[i * division_n_eta + j],
							   &dxi_x[ep], &data_x[ep],
							   &dxi_y[ep], &data_y[ep]);
				rNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
							   disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
							   weight, order_xi, order_eta,
							   calc_xi[i], calc_eta[j],
							   &info->disp_x[i * division_n_eta + j], &info->disp_y[i * division_n_eta + j],
							   &dxi_disp_x[ep], &data_disp_x[ep],
							   &dxi_disp_y[ep], &data_disp_y[ep]);
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
							&info->coord_x[i * division_n_eta + j], &info->coord_y[i * division_n_eta + j],
							&dxi_x[ep], &data_x[ep],
							&dxi_y[ep], &data_y[ep]);
			rlNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
							disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
							weight, order_xi, order_eta,
							calc_xi[i], calc_eta[j],
							&info->disp_x[i * division_n_eta + j], &info->disp_y[i * division_n_eta + j],
							&dxi_disp_x[ep], &data_disp_x[ep],
							&dxi_disp_y[ep], &data_disp_y[ep]);

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
							&info->coord_x[i * division_n_eta + j], &info->coord_y[i * division_n_eta + j],
							&dxi_x[ep], &data_x[ep],
							&dxi_y[ep], &data_y[ep]);
			lrNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
							disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
							weight, order_xi, order_eta,
							calc_xi[i], calc_eta[j],
							&info->disp_x[i * division_n_eta + j], &info->disp_y[i * division_n_eta + j],
							&dxi_disp_x[ep], &data_disp_x[ep],
							&dxi_disp_y[ep], &data_disp_y[ep]);
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
					temp_matrix[1][0] = data_x[ep] * temp2;
					temp_matrix[1][1] = data_y[ep] * temp2;

					InverseMatrix_2x2(temp_matrix);

					info->strain[ep * 3 + 0] = temp_matrix[0][0] * temp1 * dxi_disp_x[ep] + temp_matrix[0][1] * temp2 * data_disp_x[ep];
					info->strain[ep * 3 + 1] = temp_matrix[1][0] * temp1 * dxi_disp_y[ep] + temp_matrix[1][1] * temp2 * data_disp_y[ep];
					info->strain[ep * 3 + 2] = temp_matrix[1][0] * temp1 * dxi_disp_x[ep] + temp_matrix[1][1] * temp2 * data_disp_x[ep] + temp_matrix[0][0] * temp1 * dxi_disp_y[ep] + temp_matrix[0][1] * temp2 * data_disp_y[ep];
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

					info->stress[ep * 3 + 0] = info->D[0 * D_MATRIX_SIZE + 0] * info->strain[ep * 3 + 0] + info->D[0 * D_MATRIX_SIZE + 1] * info->strain[ep * 3 + 1];
					info->stress[ep * 3 + 1] = info->D[1 * D_MATRIX_SIZE + 0] * info->strain[ep * 3 + 0] + info->D[1 * D_MATRIX_SIZE + 1] * info->strain[ep * 3 + 1];
					info->stress[ep * 3 + 2] = info->D[2 * D_MATRIX_SIZE + 2] * info->strain[ep * 3 + 2];
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
			temp1 = info->disp_x[i * division_n_eta + j];
			temp2 = info->disp_y[i * division_n_eta + j];
			temp3 = sqrt(temp1 * temp1 + temp2 * temp2);
			fprintf(fp, "% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\n", info->coord_x[i * division_n_eta + j], info->coord_y[i * division_n_eta + j], info->disp_x[i * division_n_eta + j], info->disp_y[i * division_n_eta + j], temp3);
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
							info->strain[ep * 3 + 0], info->strain[ep * 3 + 1], info->strain[ep * 3 + 2],
							info->stress[ep * 3 + 0], info->stress[ep * 3 + 1], info->stress[ep * 3 + 2]);
				}
			}
		}
	}
	fclose(fp);

	// グラフ用ファイル書き込み
	fp = fopen("disp.txt", "a");
	for (i = 0; i < division_n_xi; i++)
	{
		for (j = 0; j < division_n_eta; j++)
		{
			temp1 = info->disp_x[i * division_n_eta + j];
			temp2 = info->disp_y[i * division_n_eta + j];
			fprintf(fp, "%d\t% 1le\t% 1le\t% 1le\t% 1le\n",
					patch_num,
					info->coord_x[i * division_n_eta + j], info->coord_y[i * division_n_eta + j],
					info->disp_x[i * division_n_eta + j], info->disp_y[i * division_n_eta + j]);
		}
	}
	fclose(fp);

	fp = fopen("stress_y.txt", "a");
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

					fprintf(fp, "%d\t% 1le\t% 1le\t% 1le\n",
							patch_num,
							info->coord_x[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
							info->coord_y[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
							info->stress[ep * 3 + 1]);
				}
			}
		}
	}
	fclose(fp);

	fp = fopen("stress_y_0.txt", "a");
	for (i = 0; i < element_n_xi; i++)
	{
		for (j = 0; j < element_n_eta; j++)
		{
			for (k = 0; k < division_ele_xi + 1; k++)
			{
				for (l = 0; l < division_ele_eta + 1; l++)
				{
					if (info->coord_y[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)] == 0.000000000000000e+00)
					{
						e = i * (element_n_eta + 1) + j;
						p = k * (division_ele_eta + 1) + l;
						ep = e * DIVISION_ELEMENT + p;

						fprintf(fp, "%d\t% 1le\t% 1le\t% 1le\n",
								patch_num,
								info->coord_x[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
								info->coord_y[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
								info->stress[ep * 3 + 1]);
					}
				}
			}
		}
	}
	fclose(fp);

	fp = fopen("stress_vm.txt", "a");
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
					double sum = info->stress[ep * 3 + 0] + info->stress[ep * 3 + 1];
					double dif = info->stress[ep * 3 + 0] - info->stress[ep * 3 + 1];
					double tau2 = info->stress[ep * 3 + 2] * info->stress[ep * 3 + 2];
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
							info->coord_x[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
							info->coord_y[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
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
						 double *cntl_px_loc, double *cntl_py_loc, double *weight_loc,
						 int order_xi_glo, int order_eta_glo,
						 int cntl_p_n_xi_glo, int cntl_p_n_eta_glo,
						 double *knot_vec_xi_glo, double *knot_vec_eta_glo,
						 double *cntl_px_glo, double *cntl_py_glo,
						 double *disp_cntl_px_glo, double *disp_cntl_py_glo, double *weight_glo,
						 int patch_num, information *info)
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

	double *dxi_x = (double *)malloc(sizeof(double) * info->Total_Element_to_mesh[Total_mesh] * DIVISION_ELEMENT);		// dxi_x[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]
	double *dxi_y = (double *)malloc(sizeof(double) * info->Total_Element_to_mesh[Total_mesh] * DIVISION_ELEMENT);		// dxi_y[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]
	double *data_x = (double *)malloc(sizeof(double) * info->Total_Element_to_mesh[Total_mesh] * DIVISION_ELEMENT);		// data_x[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]
	double *data_y = (double *)malloc(sizeof(double) * info->Total_Element_to_mesh[Total_mesh] * DIVISION_ELEMENT);		// data_y[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]

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
						   &info->coord_x[i * division_n_eta + j], &info->coord_y[i * division_n_eta + j],
						   &dxi_x[ep], &data_x[ep],
						   &dxi_y[ep], &data_y[ep]);

			CalcXiEtaByNR(info->coord_x[i * division_n_eta + j], info->coord_y[i * division_n_eta + j],
						  knot_vec_xi_glo, knot_vec_eta_glo,
						  cntl_px_glo, cntl_py_glo,
						  disp_cntl_px_glo, disp_cntl_py_glo,
						  cntl_p_n_xi_glo, cntl_p_n_eta_glo,
						  weight_glo, order_xi_glo, order_eta_glo,
						  &output_xi, &output_eta,
						  &disp_x_glo, &disp_y_glo,
						  &strain_xx_glo, &strain_yy_glo, &strain_xy_glo);

			// ローカル内の表示点上のグローバル変位
			info->disp_x[i * division_n_eta + j] += disp_x_glo;
			info->disp_y[i * division_n_eta + j] += disp_y_glo;

			// ローカル内の表示点上のグローバルひずみ
			info->strain[ep * 3 + 0] += strain_xx_glo;
			info->strain[ep * 3 + 1] += strain_yy_glo;
			info->strain[ep * 3 + 2] += strain_xy_glo;
			if (jj > 0 && ll == 0)
			{
				e = ii * (element_n_eta + 1) + (jj - 1);
				p = kk * (division_ele_eta + 1) + division_ele_eta;
				ep = e * DIVISION_ELEMENT + p;

				info->strain[ep * 3 + 0] += strain_xx_glo;
				info->strain[ep * 3 + 1] += strain_yy_glo;
				info->strain[ep * 3 + 2] += strain_xy_glo;
			}
			if (ii > 0 && kk == 0)
			{
				e = (ii - 1) * (element_n_eta + 1) + jj;
				p = division_ele_xi * (division_ele_eta + 1) + ll;
				ep = e * DIVISION_ELEMENT + p;

				info->strain[ep * 3 + 0] += strain_xx_glo;
				info->strain[ep * 3 + 1] += strain_yy_glo;
				info->strain[ep * 3 + 2] += strain_xy_glo;
			}
			if (ii > 0 && jj > 0 && kk == 0 && ll == 0)
			{
				e = (ii - 1) * (element_n_eta + 1) + (jj - 1);
				p = division_ele_xi * (division_ele_eta + 1) + division_ele_eta;
				ep = e * DIVISION_ELEMENT + p;

				info->strain[ep * 3 + 0] += strain_xx_glo;
				info->strain[ep * 3 + 1] += strain_yy_glo;
				info->strain[ep * 3 + 2] += strain_xy_glo;
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

					info->stress[ep * 3 + 0] = info->D[0 * D_MATRIX_SIZE + 0] * info->strain[ep * 3 + 0] + info->D[0 * D_MATRIX_SIZE + 1] * info->strain[ep * 3 + 1];
					info->stress[ep * 3 + 1] = info->D[1 * D_MATRIX_SIZE + 0] * info->strain[ep * 3 + 0] + info->D[1 * D_MATRIX_SIZE + 1] * info->strain[ep * 3 + 1];
					info->stress[ep * 3 + 2] = info->D[2 * D_MATRIX_SIZE + 2] * info->strain[ep * 3 + 2];
				}
			}
		}
	}

	// 書き込み
	fp = fopen("overlay_view.dat", "a");
	fprintf(fp, "%d\t%d\t%d\t%d\n", division_n_xi, division_n_eta, element_n_xi, element_n_eta);
	for (i = 0; i < division_n_xi; i++)
	{
		for (j = 0; j < division_n_eta; j++)
		{
			temp1 = info->disp_x[i * division_n_eta + j];
			temp2 = info->disp_y[i * division_n_eta + j];
			temp3 = sqrt(temp1 * temp1 + temp2 * temp2);
			fprintf(fp, "% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\n", info->coord_x[i * division_n_eta + j], info->coord_y[i * division_n_eta + j], info->disp_x[i * division_n_eta + j], info->disp_y[i * division_n_eta + j], temp3);
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
							info->strain[ep * 3 + 0], info->strain[ep * 3 + 1], info->strain[ep * 3 + 2],
							info->stress[ep * 3 + 0], info->stress[ep * 3 + 1], info->stress[ep * 3 + 2]);
				}
			}
		}
	}
	fclose(fp);

	// グラフ用ファイル書き込み
	fp = fopen("overlay_disp.txt", "a");
	for (i = 0; i < division_n_xi; i++)
	{
		for (j = 0; j < division_n_eta; j++)
		{
			temp1 = info->disp_x[i * division_n_eta + j];
			temp2 = info->disp_y[i * division_n_eta + j];
			fprintf(fp, "%d\t% 1le\t% 1le\t% 1le\t% 1le\n",
					patch_num,
					info->coord_x[i * division_n_eta + j], info->coord_y[i * division_n_eta + j],
					info->disp_x[i * division_n_eta + j], info->disp_y[i * division_n_eta + j]);
		}
	}
	fclose(fp);

	fp = fopen("overlay_stress_x.txt", "a");
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
							info->coord_x[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
							info->coord_y[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
							info->stress[ep * 3 + 0]);
				}
			}
		}
	}
	fclose(fp);

	fp = fopen("overlay_stress_y.txt", "a");
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

					fprintf(fp, "% 1le\t% 1le\t% 1le\n",
							info->coord_x[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
							info->coord_y[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
							info->stress[ep * 3 + 1]);
				}
			}
		}
	}
	fclose(fp);

	fp = fopen("overlay_stress_y_0.txt", "a");
	for (i = 0; i < element_n_xi; i++)
	{
		for (j = 0; j < element_n_eta; j++)
		{
			for (k = 0; k < division_ele_xi + 1; k++)
			{
				for (l = 0; l < division_ele_eta + 1; l++)
				{
					if (info->coord_y[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)] == 0.000000000000000e+00)
					{
						e = i * (element_n_eta + 1) + j;
						p = k * (division_ele_eta + 1) + l;
						ep = e * DIVISION_ELEMENT + p;

						fprintf(fp, "% 1le\t% 1le\t% 1le\n",
								info->coord_x[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
								info->coord_y[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
								info->stress[ep * 3 + 1]);
					}
				}
			}
		}
	}
	fclose(fp);

	fp = fopen("overlay_stress_r_theta.txt", "a");
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
					double sum = info->stress[ep * 3 + 0] + info->stress[ep * 3 + 1];
					double dif = info->stress[ep * 3 + 0] - info->stress[ep * 3 + 1];
					double tau2 = info->stress[ep * 3 + 2] * info->stress[ep * 3 + 2];
					stress_rr = sum * 0.5 + sqrt(dif * dif + 4 * tau2) * 0.5;
					stress_theta = sum * 0.5 - sqrt(dif * dif + 4 * tau2) * 0.5;
					fprintf(fp, "% 1.10e\t% 1.10e\t% 1.10e\t% 1.10e\t% 1.10e\t% 1.10e\n",
							calc_xi_loc[i * division_ele_xi + k],
							calc_eta_loc[j * division_ele_eta + l],
							info->coord_x[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
							info->coord_y[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
							stress_rr, stress_theta);
				}
			}
		}
	}
	fclose(fp);

	fp = fopen("overlay_stress_vm.txt", "a");
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
					double sum = info->stress[ep * 3 + 0] + info->stress[ep * 3 + 1];
					double dif = info->stress[ep * 3 + 0] - info->stress[ep * 3 + 1];
					double tau2 = info->stress[ep * 3 + 2] * info->stress[ep * 3 + 2];
					double temp1, temp2;
					temp1 = 0.5 * sum;
					temp2 = 0.5 * sqrt(dif * dif + 4 * tau2);
					stress_vm = sqrt(temp1 * temp1 + 3 * temp2 * temp2);
					fprintf(fp, "% 1.10e\t% 1.10e\t% 1.10e\t% 1.10e\t% 1.10e\n",
							calc_xi_loc[i * division_ele_xi + k],
							calc_eta_loc[j * division_ele_eta + l],
							info->coord_x[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
							info->coord_y[(i * division_ele_xi + k) * division_n_eta + (j * division_ele_eta + l)],
							stress_vm);
				}
			}
		}
	}
	fclose(fp);

	free(calc_xi_loc), free(calc_eta_loc);
	free(dxi_x), free(dxi_y), free(data_x), free(data_y);
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
	double temp_dxi, temp_data;
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
			double dxi_x, data_x, dxi_y, data_y;
			double dxi_disp_x, data_disp_x, dxi_disp_y, data_disp_y;
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
						   &dxi_x, &data_x,
						   &dxi_y, &data_y);

			rNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
						   disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
						   weight, order_xi, order_eta,
						   temp_xi, temp_eta,
						   &disp_x, &disp_y,
						   &dxi_disp_x, &data_disp_x,
						   &dxi_disp_y, &data_disp_y);

			temp_matrix2[0][0] = dxi_x * dtilda_xi;
			temp_matrix2[0][1] = dxi_y * dtilda_xi;
			temp_matrix2[1][0] = data_x * dtilda_eta;
			temp_matrix2[1][1] = data_y * dtilda_eta;

			InverseMatrix_2x2(temp_matrix2);

			strain_xx = temp_matrix2[0][0] * dtilda_xi * dxi_disp_x + temp_matrix2[0][1] * dtilda_eta * data_disp_x;
			strain_yy = temp_matrix2[1][0] * dtilda_xi * dxi_disp_y + temp_matrix2[1][1] * dtilda_eta * data_disp_y;
			strain_xy = temp_matrix2[1][0] * dtilda_xi * dxi_disp_x + temp_matrix2[1][1] * dtilda_eta * data_disp_x + temp_matrix2[0][0] * dtilda_xi * dxi_disp_y + temp_matrix2[0][1] * dtilda_eta * data_disp_y;

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
		temp_data = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_data;
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
			double dxi_x, data_x, dxi_y, data_y;
			double dxi_disp_x, data_disp_x, dxi_disp_y, data_disp_y;
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
						   &dxi_x, &data_x,
						   &dxi_y, &data_y);

			lNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
						   disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
						   weight, order_xi, order_eta,
						   temp_xi, temp_eta,
						   &disp_x, &disp_y,
						   &dxi_disp_x, &data_disp_x,
						   &dxi_disp_y, &data_disp_y);

			temp_matrix2[0][0] = dxi_x * dtilda_xi;
			temp_matrix2[0][1] = dxi_y * dtilda_xi;
			temp_matrix2[1][0] = data_x * dtilda_eta;
			temp_matrix2[1][1] = data_y * dtilda_eta;

			InverseMatrix_2x2(temp_matrix2);

			strain_xx = temp_matrix2[0][0] * dtilda_xi * dxi_disp_x + temp_matrix2[0][1] * dtilda_eta * data_disp_x;
			strain_yy = temp_matrix2[1][0] * dtilda_xi * dxi_disp_y + temp_matrix2[1][1] * dtilda_eta * data_disp_y;
			strain_xy = temp_matrix2[1][0] * dtilda_xi * dxi_disp_x + temp_matrix2[1][1] * dtilda_eta * data_disp_x + temp_matrix2[0][0] * dtilda_xi * dxi_disp_y + temp_matrix2[0][1] * dtilda_eta * data_disp_y;

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
		temp_data = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_data;
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
			double dxi_x, data_x, dxi_y, data_y;
			double dxi_disp_x, data_disp_x, dxi_disp_y, data_disp_y;
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
							&dxi_x, &data_x,
							&dxi_y, &data_y);

			rlNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
							disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
							weight, order_xi, order_eta,
							temp_xi, temp_eta,
							&disp_x, &disp_y,
							&dxi_disp_x, &data_disp_x,
							&dxi_disp_y, &data_disp_y);

			temp_matrix2[0][0] = dxi_x * dtilda_xi;
			temp_matrix2[0][1] = dxi_y * dtilda_xi;
			temp_matrix2[1][0] = data_x * dtilda_eta;
			temp_matrix2[1][1] = data_y * dtilda_eta;

			InverseMatrix_2x2(temp_matrix2);

			strain_xx = temp_matrix2[0][0] * dtilda_xi * dxi_disp_x + temp_matrix2[0][1] * dtilda_eta * data_disp_x;
			strain_yy = temp_matrix2[1][0] * dtilda_xi * dxi_disp_y + temp_matrix2[1][1] * dtilda_eta * data_disp_y;
			strain_xy = temp_matrix2[1][0] * dtilda_xi * dxi_disp_x + temp_matrix2[1][1] * dtilda_eta * data_disp_x + temp_matrix2[0][0] * dtilda_xi * dxi_disp_y + temp_matrix2[0][1] * dtilda_eta * data_disp_y;

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
		temp_data = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_data;
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
			double dxi_x, data_x, dxi_y, data_y;
			double dxi_disp_x, data_disp_x, dxi_disp_y, data_disp_y;
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
							&dxi_x, &data_x,
							&dxi_y, &data_y);

			lrNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
							disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
							weight, order_xi, order_eta,
							temp_xi, temp_eta,
							&disp_x, &disp_y,
							&dxi_disp_x, &data_disp_x,
							&dxi_disp_y, &data_disp_y);

			temp_matrix2[0][0] = dxi_x * dtilda_xi;
			temp_matrix2[0][1] = dxi_y * dtilda_xi;
			temp_matrix2[1][0] = data_x * dtilda_eta;
			temp_matrix2[1][1] = data_y * dtilda_eta;

			InverseMatrix_2x2(temp_matrix2);

			strain_xx = temp_matrix2[0][0] * dtilda_xi * dxi_disp_x + temp_matrix2[0][1] * dtilda_eta * data_disp_x;
			strain_yy = temp_matrix2[1][0] * dtilda_xi * dxi_disp_y + temp_matrix2[1][1] * dtilda_eta * data_disp_y;
			strain_xy = temp_matrix2[1][0] * dtilda_xi * dxi_disp_x + temp_matrix2[1][1] * dtilda_eta * data_disp_x + temp_matrix2[0][0] * dtilda_xi * dxi_disp_y + temp_matrix2[0][1] * dtilda_eta * data_disp_y;

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
		temp_data = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_data;
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
void output_for_viewer(information *info)
{
	// 全てのローカルパッチについて1つのinput.txt作成, 重ね合わせ結果出力のため
	int i, j, k;

	int l_patch, g_patch;
	int l_cntl_p, g_cntl_p;
	int l_cnst, g_cnst;
	int l_load, g_load;
	int l_dist, g_dist;
	int l_temp;

	g_patch = info->Total_Patch_to_mesh[1];
	g_cntl_p = info->Total_Control_Point_to_mesh[1];
	g_cnst = info->Total_Constraint_to_mesh[1];
	g_load = info->Total_Load_to_mesh[1];
	g_dist = info->Total_DistributeForce_to_mesh[1];

	fp = fopen("input_local.txt", "w");
	fprintf(fp, "%.3f\t%.6f\n\n", E, nu);

	l_patch = info->Total_Patch_to_mesh[Total_mesh] - g_patch;
	fprintf(fp, "%d\n\n", l_patch);

	l_cntl_p = info->Total_Control_Point_to_mesh[Total_mesh] - g_cntl_p;
	fprintf(fp, "%d\n\n", l_cntl_p);

	for (i = 0; i < l_patch; i++)
	{
		fprintf(fp, "%d\t%d\n", info->Order[(i + g_patch) * info->DIMENSION + 0], info->Order[(i + g_patch) * info->DIMENSION + 1]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < l_patch; i++)
	{
		fprintf(fp, "%d\t%d\n", info->No_knot[(i + g_patch) * info->DIMENSION + 0], info->No_knot[(i + g_patch) * info->DIMENSION + 1]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < l_patch; i++)
	{
		fprintf(fp, "%d\t%d\n", info->No_Control_point[(i + g_patch) * info->DIMENSION + 0], info->No_Control_point[(i + g_patch) * info->DIMENSION + 1]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < l_patch; i++)
	{
		for (j = 0; j < info->No_Control_point_in_patch[i + g_patch]; j++)
		{
			l_temp = info->Patch_Control_point[info->Total_Control_Point_to_patch[i + g_patch] + j] - g_cntl_p;
			fprintf(fp, "%d\t", l_temp);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp,"\n");

	l_cnst = info->Total_Constraint_to_mesh[Total_mesh] - g_cnst;
	l_load = info->Total_Load_to_mesh[Total_mesh] - g_load;
	l_dist = info->Total_DistributeForce_to_mesh[Total_mesh] - g_dist;
	fprintf(fp,"%d\t%d\t%d\n\n", l_cnst, l_load, l_dist);

	for (i = 0; i < l_patch; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			for (k = 0; k < info->No_knot[(i + g_patch) * info->DIMENSION + j]; k++)
			{
				fprintf(fp, "%.5f\t", info->Position_Knots[info->Total_Knot_to_patch_dim[(i + g_patch) * info->DIMENSION + j] + k]);
			}
			fprintf(fp, "\n");
		}
	}
	fprintf(fp,"\n");

	for (i = 0; i < l_cntl_p; i++)
	{
		fprintf(fp, "%d\t", i);
		for (j = 0; j < info->DIMENSION + 1; j++)
		{
			fprintf(fp, "%.10f\t", info->Node_Coordinate[(i + g_cntl_p) * (info->DIMENSION + 1) + j]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp,"\n");

	for (i = 0; i < l_cnst; i++)
	{
		l_temp = info->Constraint_Node_Dir[(i + g_cnst) * 2 + 0] - g_cntl_p;
		fprintf(fp, "%d\t%d\t%.10f\n", l_temp, info->Constraint_Node_Dir[(i + g_cnst) * 2 + 1], info->Value_of_Constraint[i + g_cnst]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < l_load; i++)
	{
		l_temp = info->Load_Node_Dir[(i + g_load) * 2 + 0] - g_cntl_p;
		fprintf(fp, "%d\t%d\t%.10f\n", l_temp, info->Load_Node_Dir[(i + g_load) * 2 + 1], info->Value_of_Load[i + g_load]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < l_dist; i++)
	{
		fprintf(fp, "%d\t%d\t%d\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n",
				info->type_load_array[i + g_dist],
				info->iPatch_array[i + g_dist],
				info->iCoord_array[i + g_dist],
				info->val_Coord_array[i + g_dist],
				info->Range_Coord_array[(i + g_dist) * 2 + 0],
				info->Range_Coord_array[(i + g_dist) * 2 + 1],
				info->Coeff_Dist_Load_array[(i + g_dist) * 3 + 0],
				info->Coeff_Dist_Load_array[(i + g_dist) * 3 + 1],
				info->Coeff_Dist_Load_array[(i + g_dist) * 3 + 2]);
	}
	fclose(fp);

	fp = fopen("input_for_NURBS.txt", "w");
	fprintf(fp, "%.3f\t%.6f\n\n", E, nu);

	fprintf(fp, "%d\n\n", info->Total_Patch_to_mesh[Total_mesh]);

	fprintf(fp, "%d\n\n", info->Total_Control_Point_to_mesh[Total_mesh]);

	for (i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp, "%d\t%d\n", info->Order[i * info->DIMENSION + 0], info->Order[i * info->DIMENSION + 1]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp, "%d\t%d\n", info->No_knot[i * info->DIMENSION + 0], info->No_knot[i * info->DIMENSION + 1]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp, "%d\t%d\n", info->No_Control_point[i * info->DIMENSION + 0], info->No_Control_point[i * info->DIMENSION + 1]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
	{
		for (j = 0; j < info->No_Control_point_in_patch[i]; j++)
		{
			fprintf(fp, "%d\t", info->Patch_Control_point[info->Total_Control_Point_to_patch[i] + j]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp,"\n");

	fprintf(fp,"%d\t%d\t%d\n\n",
				info->Total_Constraint_to_mesh[Total_mesh],
				info->Total_Load_to_mesh[Total_mesh],
				info->Total_DistributeForce_to_mesh[Total_mesh]);

	for (i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
	{
		for (j = 0; j < info->DIMENSION; j++)
		{
			for (k = 0; k < info->No_knot[i * info->DIMENSION + j]; k++)
			{
				fprintf(fp, "%.5f\t", info->Position_Knots[info->Total_Knot_to_patch_dim[i * info->DIMENSION + j] + k]);
			}
			fprintf(fp, "\n");
		}
	}
	fprintf(fp,"\n");

	for (i = 0; i < info->Total_Control_Point_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp, "%d\t", i);
		for (j = 0; j < info->DIMENSION + 1; j++)
		{
			fprintf(fp, "%.10f\t", info->Node_Coordinate[i * (info->DIMENSION + 1) + j]);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp,"\n");

	for (i = 0; i < info->Total_Constraint_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp, "%d\t%d\t%.10f\n",
				info->Constraint_Node_Dir[i * 2 + 0],
				info->Constraint_Node_Dir[i * 2 + 1],
				info->Value_of_Constraint[i]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < info->Total_Load_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp, "%d\t%d\t%.10f\n",
				info->Load_Node_Dir[i * 2 + 0],
				info->Load_Node_Dir[i * 2 + 1],
				info->Value_of_Load[i]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < info->Total_DistributeForce_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp, "%d\t%d\t%d\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n",
				info->type_load_array[i],
				info->iPatch_array[i],
				info->iCoord_array[i],
				info->val_Coord_array[i],
				info->Range_Coord_array[i * 2 + 0],
				info->Range_Coord_array[i * 2 + 1],
				info->Coeff_Dist_Load_array[i * 3 + 0],
				info->Coeff_Dist_Load_array[i * 3 + 1],
				info->Coeff_Dist_Load_array[i * 3 + 2]);
	}

	fclose(fp);

	fp = fopen("Displacement.dat", "w");
	fprintf(fp, "label=Displacement\n");
	fprintf(fp, "num_items=%d\n", info->Total_Control_Point_to_mesh[Total_mesh]);
	fprintf(fp, "\n");
	for (j = 0; j < info->Total_Control_Point_to_mesh[Total_mesh]; j++)
	{
		fprintf(fp, "%d:	%.16e %.16e ", j, info->Displacement[j * info->DIMENSION + 0], info->Displacement[j * info->DIMENSION + 1]);
		fprintf(fp, "\n");
	}
	fclose(fp);

	fp = fopen("Displacement_loc.dat", "w");
	int loc_cntl_p_n = info->Total_Control_Point_to_mesh[Total_mesh] - info->Total_Control_Point_to_mesh[1];
	int glo_cntl_p_n = info->Total_Control_Point_to_mesh[1];
	fprintf(fp, "label=Displacement\n"
				"num_items=%d\n\n", loc_cntl_p_n);
	for (i = 0; i < loc_cntl_p_n; i++)
	{
		fprintf(fp, "%d:	%.16e %.16e \n", i, info->Displacement[(i + glo_cntl_p_n) * info->DIMENSION + 0], info->Displacement[(i + glo_cntl_p_n) * info->DIMENSION + 1]);
	}
	fclose(fp);
}


void IGA_view(information *info)
{
	int i, j;

	// 一要素の分割数
	division_ele_xi = DIVISION_ELE;
	division_ele_eta = DIVISION_ELE;

	// output
	fp = fopen("view.dat", "w");
	fprintf(fp, "%d\t%d\t%d\n", 1, division_ele_xi, division_ele_eta);
	fclose(fp);

	fp = fopen("disp.txt", "w");
	fprintf(fp, "patch_n\tx\ty\tdisp_x\tdisp_y\n");
	fclose(fp);

	fp = fopen("stress_y.txt", "w");
	fprintf(fp, "patch_n\tx\ty\tstress_yy\n");
	fclose(fp);

	fp = fopen("stress_y_0.txt", "w");
	fprintf(fp, "patch_n\tx\ty\tstress_yy\n");
	fclose(fp);

	fp = fopen("stress_vm.txt", "w");
	fprintf(fp, "xi\teta\tx\ty\tstress_vm\n");
	fclose(fp);

	for (i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
	{
		fp = fopen("stress_vm.txt", "a");
		fprintf(fp, "\npatch_n;%d\n\n", i);
		fclose(fp);

		// パッチ i での各方向ノットベクトル, Coodinate, displacement
		double *temp_Position_Knots_xi = (double *)malloc(sizeof(double) * info->No_knot[i * info->DIMENSION + 0]);
		double *temp_Position_Knots_eta = (double *)malloc(sizeof(double) * info->No_knot[i * info->DIMENSION + 1]);
		double *temp_Coord_x = (double *)malloc(sizeof(double) * info->No_Control_point_in_patch[i]);
		double *temp_Coord_y = (double *)malloc(sizeof(double) * info->No_Control_point_in_patch[i]);
		double *temp_Coord_w = (double *)malloc(sizeof(double) * info->No_Control_point_in_patch[i]);
		double *temp_disp_x = (double *)malloc(sizeof(double) * info->No_Control_point_in_patch[i]);
		double *temp_disp_y = (double *)malloc(sizeof(double) * info->No_Control_point_in_patch[i]);

		for (j = 0; j < info->No_knot[i * info->DIMENSION + 0]; j++)
		{
			temp_Position_Knots_xi[j] = info->Position_Knots[info->Total_Knot_to_patch_dim[i * info->DIMENSION + 0] + j];
		}
		for (j = 0; j < info->No_knot[i * info->DIMENSION + 1]; j++)
		{
			temp_Position_Knots_eta[j] = info->Position_Knots[info->Total_Knot_to_patch_dim[i * info->DIMENSION + 1] + j];
		}

		for (j = 0; j < info->No_Control_point_in_patch[i]; j++)
		{
			temp_Coord_x[j] = info->Control_Coord_x[info->Patch_Control_point[info->Total_Control_Point_to_patch[i] + j]];
			temp_Coord_y[j] = info->Control_Coord_y[info->Patch_Control_point[info->Total_Control_Point_to_patch[i] + j]];
			temp_Coord_w[j] = info->Control_Weight[info->Patch_Control_point[info->Total_Control_Point_to_patch[i] + j]];
			temp_disp_x[j] = info->Displacement[info->Patch_Control_point[info->Total_Control_Point_to_patch[i] + j] * info->DIMENSION + 0];
			temp_disp_y[j] = info->Displacement[info->Patch_Control_point[info->Total_Control_Point_to_patch[i] + j] * info->DIMENSION + 1];
		}

		Calculation(info->Order[i * info->DIMENSION + 0], info->Order[i * info->DIMENSION + 1], info->No_knot[i * info->DIMENSION + 0], info->No_knot[i * info->DIMENSION + 1], info->No_Control_point[i * info->DIMENSION + 0], info->No_Control_point[i * info->DIMENSION + 1],
					temp_Position_Knots_xi, temp_Position_Knots_eta, temp_Coord_x, temp_Coord_y, temp_Coord_w,
					temp_disp_x, temp_disp_y, i, info);
	}
}


void K_output_svg(information *info)
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


// paraview
void output_for_paraview(information *info)
{
	int i, j;
	int point_on_element = 0;

	int num = info->Total_Control_Point_to_mesh[Total_mesh];
	int digit = 0;
	while(num != 0)
	{
		num = num / 10;
		digit++;
  	}

	Make_connectivity(info);
	Make_info_for_viewer(info);

	// global patch
	char str_glo[256] = "global_patch.xmf";
	fp = fopen(str_glo, "w");

	fprintf(fp, "<?xml version=\"1.0\" ?>\n");
	fprintf(fp, "<Xdmf Version=\"2.0\">\n");
	fprintf(fp, "  <Domain>\n");
	fprintf(fp, "    <Grid Name=\"ien\">\n");

	// TopologyType
	if (info->DIMENSION == 2)
	{
		point_on_element = 9;
		fprintf(fp, "      <Topology TopologyType=\"Quadrilateral_9\" NumberOfElements=\"%d\">\n", info->Total_Element_to_mesh[1]);
	}
	else if (info->DIMENSION == 3)
	{
		point_on_element = 27;
		fprintf(fp, "      <Topology TopologyType=\"HEXAHEDRON_27\" NumberOfElements=\"%d\">\n", info->Total_Element_to_mesh[1]);
	}
	fprintf(fp, "        <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/DataItem[@Name=&quot;ien&quot;]</DataItem>\n");
	fprintf(fp, "      </Topology>\n");

	// Geometry
	if (info->DIMENSION == 2)
	{
		fprintf(fp, "      <Geometry GeometryType=\"XY\">\n");
	}
	else if (info->DIMENSION == 3)
	{
		fprintf(fp, "      <Geometry GeometryType=\"XYZ\">\n");
	}
	fprintf(fp, "        <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/DataItem[@Name=&quot;xyz&quot;]</DataItem>\n");
	fprintf(fp, "      </Geometry>\n");
	if (info->DIMENSION == 2)
	{
		fprintf(fp, "      <Geometry GeometryType=\"XY\">\n");
	}
	else if (info->DIMENSION == 3)
	{
		fprintf(fp, "      <Geometry GeometryType=\"XYZ\">\n");
	}
	fprintf(fp, "        <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/DataItem[@Name=&quot;xyz&quot;]</DataItem>\n");
	fprintf(fp, "      </Geometry>\n");

	// Attribute stress
	fprintf(fp, "      <Attribute Name=\"stress\" AttributeType=\"Vector\" Dimensions=\"%d 3\">\n", Total_connectivity_glo);
	fprintf(fp, "        <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/DataItem[@Name=&quot;stress&quot;]</DataItem>\n");
	fprintf(fp, "      </Attribute>\n");

	// Attribute strain
	fprintf(fp, "      <Attribute Name=\"strain\" AttributeType=\"Vector\" Dimensions=\"%d 3\">\n", Total_connectivity_glo);
	fprintf(fp, "        <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/DataItem[@Name=&quot;strain&quot;]</DataItem>\n");
	fprintf(fp, "      </Attribute>\n");

	// xyz
	fprintf(fp, "      <DataItem Name=\"xyz\" Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\" Endian=\"Big\">\n", Total_connectivity_glo, info->DIMENSION);
	if (info->DIMENSION == 2)
	{
		for (i = 0; i < Total_connectivity_glo; i++)
		{
			fprintf(fp, "\t\t\t\t%le %le\n", info->Connectivity_coord[i * info->DIMENSION + 0], info->Connectivity_coord[i * info->DIMENSION + 1]);
		}
	}
	else if (info->DIMENSION == 3)
	{
		for (i = 0; i < Total_connectivity_glo; i++)
		{
			fprintf(fp, "\t\t\t\t%le %le %le\n", info->Connectivity_coord[i * info->DIMENSION + 0], info->Connectivity_coord[i * info->DIMENSION + 1], info->Connectivity_coord[i * info->DIMENSION + 2]);
		}
	}
	fprintf(fp, "      </DataItem>\n");

	// ien
	fprintf(fp, "      <DataItem Name=\"ien\" Dimensions=\"%d %d\" NumberType=\"Int\" Format=\"XML\" Endian=\"Big\">\n", info->Total_Element_to_mesh[1], point_on_element);
	for (i = 0; i < info->Total_Element_to_mesh[1]; i++)
	{
		fprintf(fp, "\t\t\t\t");
		if (info->DIMENSION == 2)
		{
			for (j = 0; j < point_on_element; j++)
			{
				fprintf(fp, "%*d ", -digit, info->Connectivity[i * point_on_element + j]);
			}
		}
		else if (info->DIMENSION == 3)
		{
			for (j = 0; j < point_on_element; j++)
			{
				fprintf(fp, "%*d ", -digit, info->Connectivity[i * point_on_element + j]);
			}
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "      </DataItem>\n");

	// stress
	fprintf(fp, "      <DataItem Name=\"stress\" Dimensions=\"%d 3\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\" Endian=\"Big\">\n", Total_connectivity_glo);
	if (info->DIMENSION == 2)
	{
		for (i = 0; i < Total_connectivity_glo; i++)
		{
			fprintf(fp, "\t\t\t\t%le %le %le\n", info->stress_at_connectivity[i * N_STRESS + 0], info->stress_at_connectivity[i * N_STRESS + 1], info->stress_at_connectivity[i * N_STRESS + 3]);
		}
	}
	else if (info->DIMENSION == 3)
	{
		for (i = 0; i < Total_connectivity_glo; i++)
		{
			fprintf(fp, "\t\t\t\t%le %le %le\n", info->stress_at_connectivity[i * N_STRESS + 0], info->stress_at_connectivity[i * N_STRESS + 1], info->stress_at_connectivity[i * N_STRESS + 2]);
		}
	}
	fprintf(fp, "      </DataItem>\n");

	// strain
	fprintf(fp, "      <DataItem Name=\"strain\" Dimensions=\"%d 3\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\" Endian=\"Big\">\n", Total_connectivity_glo);
	if (info->DIMENSION == 2)
	{
		for (i = 0; i < Total_connectivity_glo; i++)
		{
			fprintf(fp, "\t\t\t\t%le %le %le\n", info->strain_at_connectivity[i * N_STRAIN + 0], info->strain_at_connectivity[i * N_STRAIN + 1], info->strain_at_connectivity[i * N_STRAIN + 3]);
		}
	}
	else if (info->DIMENSION == 3)
	{
		for (i = 0; i < Total_connectivity_glo; i++)
		{
			fprintf(fp, "\t\t\t\t%le %le %le\n", info->strain_at_connectivity[i * N_STRAIN + 0], info->strain_at_connectivity[i * N_STRAIN + 1], info->strain_at_connectivity[i * N_STRAIN + 2]);
		}
	}
	fprintf(fp, "      </DataItem>\n");
	fprintf(fp, "    </Grid>\n");
	fprintf(fp, "  </Domain>\n");
	fprintf(fp, "</Xdmf>\n");

	fclose(fp);

	// local patch (overlay)
	if (Total_mesh >= 2)
	{
		char str_loc[256] = "local_patch.xmf";
		fp = fopen(str_loc, "w");

		fprintf(fp, "<?xml version=\"1.0\" ?>\n");
		fprintf(fp, "<Xdmf Version=\"2.0\">\n");
		fprintf(fp, "  <Domain>\n");
		fprintf(fp, "    <Grid Name=\"ien\">\n");

		// TopologyType
		if (info->DIMENSION == 2)
		{
			point_on_element = 9;
			fprintf(fp, "      <Topology TopologyType=\"Quadrilateral_9\" NumberOfElements=\"%d\">\n", info->Total_Element_to_mesh[2] - info->Total_Element_to_mesh[1]);
		}
		else if (info->DIMENSION == 3)
		{
			point_on_element = 27;
			fprintf(fp, "      <Topology TopologyType=\"HEXAHEDRON_27\" NumberOfElements=\"%d\">\n", info->Total_Element_to_mesh[2] - info->Total_Element_to_mesh[1]);
		}
		fprintf(fp, "        <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/DataItem[@Name=&quot;ien&quot;]</DataItem>\n");
		fprintf(fp, "      </Topology>\n");

		// Geometry
		if (info->DIMENSION == 2)
		{
			fprintf(fp, "      <Geometry GeometryType=\"XY\">\n");
		}
		else if (info->DIMENSION == 3)
		{
			fprintf(fp, "      <Geometry GeometryType=\"XYZ\">\n");
		}
		fprintf(fp, "        <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/DataItem[@Name=&quot;xyz&quot;]</DataItem>\n");
		fprintf(fp, "      </Geometry>\n");

		// Attribute stress
		fprintf(fp, "      <Attribute Name=\"stress\" AttributeType=\"Vector\" Dimensions=\"%d 3\">\n", Total_connectivity - Total_connectivity_glo);
		fprintf(fp, "        <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/DataItem[@Name=&quot;stress&quot;]</DataItem>\n");
		fprintf(fp, "      </Attribute>\n");

		// Attribute strain
		fprintf(fp, "      <Attribute Name=\"strain\" AttributeType=\"Vector\" Dimensions=\"%d 3\">\n", Total_connectivity - Total_connectivity_glo);
		fprintf(fp, "        <DataItem Reference=\"XML\">/Xdmf/Domain/Grid/DataItem[@Name=&quot;strain&quot;]</DataItem>\n");
		fprintf(fp, "      </Attribute>\n");

		// xyz
		fprintf(fp, "      <DataItem Name=\"xyz\" Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"8\" Format=\"XML\" Endian=\"Big\">\n", Total_connectivity - Total_connectivity_glo, info->DIMENSION);
		if (info->DIMENSION == 2)
		{
			for (i = Total_connectivity_glo; i < Total_connectivity; i++)
			{
				fprintf(fp, "\t\t\t\t%le %le\n", info->Connectivity_coord[i * info->DIMENSION + 0], info->Connectivity_coord[i * info->DIMENSION + 1]);
			}
		}
		else if (info->DIMENSION == 3)
		{
			for (i = Total_connectivity_glo; i < Total_connectivity; i++)
			{
				fprintf(fp, "\t\t\t\t%le %le %le\n", info->Connectivity_coord[i * info->DIMENSION + 0], info->Connectivity_coord[i * info->DIMENSION + 1], info->Connectivity_coord[i * info->DIMENSION + 2]);
			}
		}
		fprintf(fp, "      </DataItem>\n");

		// ien
		fprintf(fp, "      <DataItem Name=\"ien\" Dimensions=\"%d %d\" NumberType=\"Int\" Format=\"XML\" Endian=\"Big\">\n", info->Total_Element_to_mesh[2] - info->Total_Element_to_mesh[1], point_on_element);
		for (i = info->Total_Element_to_mesh[1]; i < info->Total_Element_to_mesh[2]; i++)
		{
			fprintf(fp, "\t\t\t\t");
			if (info->DIMENSION == 2)
			{
				for (j = 0; j < point_on_element; j++)
				{
					fprintf(fp, "%*d ", -digit, info->Connectivity[i * point_on_element + j] - Total_connectivity_glo);
				}
			}
			else if (info->DIMENSION == 3)
			{
				for (j = 0; j < point_on_element; j++)
				{
					fprintf(fp, "%*d ", -digit, info->Connectivity[i * point_on_element + j] - Total_connectivity_glo);
				}
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "      </DataItem>\n");

		// stress
		fprintf(fp, "      <DataItem Name=\"stress\" Dimensions=\"%d 3\" NumberType=\"Float\" Format=\"XML\" Endian=\"Big\">\n", Total_connectivity - Total_connectivity_glo);
		if (info->DIMENSION == 2)
		{
			for (i = Total_connectivity_glo; i < Total_connectivity; i++)
			{
				fprintf(fp, "\t\t\t\t%le %le %le\n", info->stress_at_connectivity[i * N_STRESS + 0], info->stress_at_connectivity[i * N_STRESS + 1], info->stress_at_connectivity[i * N_STRESS + 3]);
			}
		}
		else if (info->DIMENSION == 3)
		{
			for (i = Total_connectivity_glo; i < Total_connectivity; i++)
			{
				fprintf(fp, "\t\t\t\t%le %le %le\n", info->stress_at_connectivity[i * N_STRESS + 0], info->stress_at_connectivity[i * N_STRESS + 1], info->stress_at_connectivity[i * N_STRESS + 2]);
			}
		}
		fprintf(fp, "      </DataItem>\n");

		// strain
		fprintf(fp, "      <DataItem Name=\"strain\" Dimensions=\"%d 3\" NumberType=\"Float\" Format=\"XML\" Endian=\"Big\">\n", Total_connectivity - Total_connectivity_glo);
		if (info->DIMENSION == 2)
		{
			for (i = Total_connectivity_glo; i < Total_connectivity; i++)
			{
				fprintf(fp, "\t\t\t\t%le %le %le\n", info->strain_at_connectivity[i * N_STRAIN + 0], info->strain_at_connectivity[i * N_STRAIN + 1], info->strain_at_connectivity[i * N_STRAIN + 3]);
			}
		}
		else if (info->DIMENSION == 3)
		{
			for (i = Total_connectivity_glo; i < Total_connectivity; i++)
			{
				fprintf(fp, "\t\t\t\t%le %le %le\n", info->strain_at_connectivity[i * N_STRAIN + 0], info->strain_at_connectivity[i * N_STRAIN + 1], info->strain_at_connectivity[i * N_STRAIN + 2]);
			}
		}
		fprintf(fp, "      </DataItem>\n");
		fprintf(fp, "    </Grid>\n");
		fprintf(fp, "  </Domain>\n");
		fprintf(fp, "</Xdmf>\n");

		fclose(fp);
	}
}


void Make_connectivity(information *info)
{
	int i, j, k, l, m;

	// Make Patch_check
	int Patch_check_counter = 0;
	if (info->DIMENSION == 2)
    {
		for (i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
		{
			int CP_counter = info->Total_Control_Point_to_patch[i];

			// 辺0 点0
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter];
			Patch_check_counter++;
			// 辺0 点1
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] - 1)];
			Patch_check_counter++;
			// 辺1 点0
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] - 1)];
			Patch_check_counter++;
			// 辺1 点1
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter];
			Patch_check_counter++;
			// 辺2 点0
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] - 1)];
			Patch_check_counter++;
			// 辺2 点1
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] * info->No_Control_point[i * info->DIMENSION + 1] - 1)];
			Patch_check_counter++;
			// 辺3 点0
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] * info->No_Control_point[i * info->DIMENSION + 1] - 1)];
			Patch_check_counter++;
			// 辺3 点1
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] - 1)];
			Patch_check_counter++;
			// 辺4 点0
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] * (info->No_Control_point[i * info->DIMENSION + 1] - 1))];
			Patch_check_counter++;
			// 辺4 点1
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] * info->No_Control_point[i * info->DIMENSION + 1] - 1)];
			Patch_check_counter++;
			// 辺5 点0
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] * info->No_Control_point[i * info->DIMENSION + 1] - 1)];
			Patch_check_counter++;
			// 辺5 点1
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] * (info->No_Control_point[i * info->DIMENSION + 1] - 1))];
			Patch_check_counter++;
			// 辺6 点0
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter];
			Patch_check_counter++;
			// 辺6 点1
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] * (info->No_Control_point[i * info->DIMENSION + 1] - 1))];
			Patch_check_counter++;
			// 辺7 点0
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] * (info->No_Control_point[i * info->DIMENSION + 1] - 1))];
			Patch_check_counter++;
			// 辺7 点1
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter];
			Patch_check_counter++;
		}
	}
	else if (info->DIMENSION == 3)
	{
		for (i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
		{
			int CP_counter = info->Total_Control_Point_to_patch[i];
			int temp = info->No_Control_point[i * info->DIMENSION] * info->No_Control_point[i * info->DIMENSION + 1] * (info->No_Control_point[i * info->DIMENSION + 2] - 1);
			int temp_Patch_check_counter;
			int temp_Patch_check[8];

			// 面0 点0
			temp_Patch_check[0] = Patch_check_counter;
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter];
			Patch_check_counter++;
			// 面0 点1
			temp_Patch_check[1] = Patch_check_counter;
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] - 1)];
			Patch_check_counter++;
			// 面0 点2
			temp_Patch_check[2] = Patch_check_counter;
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] * info->No_Control_point[i * info->DIMENSION + 1] - 1)];
			Patch_check_counter++;
			// 面0 点3
			temp_Patch_check[3] = Patch_check_counter;
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + (info->No_Control_point[i * info->DIMENSION] * (info->No_Control_point[i * info->DIMENSION + 1] - 1))];
			Patch_check_counter++;
			// 面1 点0
			temp_Patch_check[4] = Patch_check_counter;
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + temp];
			Patch_check_counter++;
			// 面1 点1
			temp_Patch_check[5] = Patch_check_counter;
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + temp + (info->No_Control_point[i * info->DIMENSION] - 1)];
			Patch_check_counter++;
			// 面1 点2
			temp_Patch_check[6] = Patch_check_counter;
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + temp + (info->No_Control_point[i * info->DIMENSION] * info->No_Control_point[i * info->DIMENSION + 1] - 1)];
			Patch_check_counter++;
			// 面1 点3
			temp_Patch_check[7] = Patch_check_counter;
			info->Patch_check[Patch_check_counter] = info->Patch_Control_point[CP_counter + temp + (info->No_Control_point[i * info->DIMENSION] * (info->No_Control_point[i * info->DIMENSION + 1] - 1))];
			Patch_check_counter++;
			// 面2 点0
			temp_Patch_check_counter = temp_Patch_check[0];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面2 点1
			temp_Patch_check_counter = temp_Patch_check[3];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面2 点2
			temp_Patch_check_counter = temp_Patch_check[7];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面2 点3
			temp_Patch_check_counter = temp_Patch_check[4];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面3 点0
			temp_Patch_check_counter = temp_Patch_check[0];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面3 点1
			temp_Patch_check_counter = temp_Patch_check[1];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面3 点2
			temp_Patch_check_counter = temp_Patch_check[5];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面3 点3
			temp_Patch_check_counter = temp_Patch_check[4];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面4 点0
			temp_Patch_check_counter = temp_Patch_check[1];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面4 点1
			temp_Patch_check_counter = temp_Patch_check[2];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面4 点2
			temp_Patch_check_counter = temp_Patch_check[6];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面4 点3
			temp_Patch_check_counter = temp_Patch_check[5];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面5 点0
			temp_Patch_check_counter = temp_Patch_check[3];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面5 点1
			temp_Patch_check_counter = temp_Patch_check[2];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面5 点2
			temp_Patch_check_counter = temp_Patch_check[6];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
			// 面5 点3
			temp_Patch_check_counter = temp_Patch_check[7];
			info->Patch_check[Patch_check_counter] = info->Patch_check[temp_Patch_check_counter];
			Patch_check_counter++;
		}
	}

	// Check
	int patch_start = 0;
	if (info->DIMENSION == 2)
	{
		for (i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
		{
			if (i == info->Total_Patch_to_mesh[1])
			{
				patch_start = i;
			}

			for (j = patch_start; j < i; j++)
			{
				int own = i * 16;
				int opp = j * 16;
				for (k = 0; k < 4; k++)
				{
					int kk = 2 * k;
					for (l = 0; l < 8; l++)
					{
						// 辺が一致している場合 Face_Edge_info を True
						if (info->Patch_check[own + kk * 2] == info->Patch_check[opp + l * 2] && info->Patch_check[own + kk * 2 + 1] == info->Patch_check[opp + l * 2 + 1])
						{
							info->Face_Edge_info[i * 32 + k * 8 + l] = 1;
							info->Opponent_patch_num[i * 4 + k] = j;
							goto end_loop_a;
						}
					}
				}
				end_loop_a:;
			}
		}
	}
	else if (info->DIMENSION == 3)
	{
		int check_a[4], check_b[4];
		for (i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
		{
			if (i == info->Total_Patch_to_mesh[1])
			{
				patch_start = i;
			}

			for (j = patch_start; j < i; j++)
			{
				int own = i * 24;
				int opp = j * 24;
				for (k = 0; k < 6; k++)
				{
					int kk = k * 4;
					for (m = 0; m < 4; m++)
					{
						check_a[m] = info->Patch_check[own + kk + m];
					}
					sort(4, check_a);

					for (l = 0; l < 6; l++)
					{
						int ll = l * 4;
						for (m = 0; m < 4; m++)
						{
							check_b[m] = info->Patch_check[opp + ll + m];
						}
						sort(4, check_b);

						// 面が一致している場合 Face_Edge_info を Mode 番号 (m) に
						if (check_a[0] == check_b[0] && check_a[1] == check_b[1] && check_a[2] == check_b[2] && check_a[3] == check_b[3])
						{
							for (m = 0; m < 4; m++)
							{
								if (info->Patch_check[own + kk] == info->Patch_check[opp + ll + m])
								{
									info->Face_Edge_info[i * 36 + k * 6 + l] = m;
									info->Opponent_patch_num[i * 6 + k] = j;
									goto end_loop_b;
								}
							}
						}
					}
				}
				end_loop_b:;
			}
		}
	}

	// Make connectivity
	int counter = 0;
	int point_counter = 0;
	int ele_counter = 0;
	if (info->DIMENSION == 2)
	{
		int p_x[9] = {0, 2, 2, 0, 1, 2, 1, 0, 1};
		int p_y[9] = {0, 0, 2, 2, 0, 1, 2, 1, 1};
		for (i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
		{
			int own = 0;
			for (j = 0; j < i; j++)
			{
				int a, b;
				a = info->line_No_Total_element[j * info->DIMENSION + 0] * 2 + 1;
				b = info->line_No_Total_element[j * info->DIMENSION + 1] * 2 + 1;
				own += 2 * (a + b);
			}

			int Edge[4] = {0};
    		int Edge_counter = i * 32;
			int temp_n = 0;

			int ii = info->line_No_Total_element[i * info->DIMENSION + 0] * 2 + 1;
			int jj = info->line_No_Total_element[i * info->DIMENSION + 1] * 2 + 1;

			// 重なっている辺の Patch_array を作成
			for (j = 0; j < 4; j++)
			{
				// j == 0 ではなにもしない
				if (j == 1)
				{
					own += ii;
				}
				else if (j == 2)
				{
					own += jj;
				}
				else if (j == 3)
				{
					own += ii;
				}

				for (k = 0; k < 8; k++)
				{
					if (info->Face_Edge_info[Edge_counter + j * 8 + k] == 1)
					{
						int opp = 0;
						for (l = 0; l < info->Opponent_patch_num[i * 4 + j]; l++)
						{
							int a, b;
							a = info->line_No_Total_element[l * info->DIMENSION + 0] * 2 + 1;
							b = info->line_No_Total_element[l * info->DIMENSION + 1] * 2 + 1;
							opp += 2 * (a + b);
						}

						Edge[j] = 1;

						int p = k / 2;
						int q = k % 2;

						int temp_ii = info->line_No_Total_element[info->Opponent_patch_num[i * 4 + j] * info->DIMENSION + 0] * 2 + 1;
						int temp_jj = info->line_No_Total_element[info->Opponent_patch_num[i * 4 + j] * info->DIMENSION + 1] * 2 + 1;

						if (p == 0)
						{
							temp_n = temp_ii;
						}
						else if (p == 1)
						{
							temp_n = temp_jj;
							opp += temp_ii;
						}
						else if (p == 2)
						{
							temp_n = temp_ii;
							opp += temp_ii + temp_jj;
						}
						else if (p == 3)
						{
							temp_n = temp_jj;
							opp += 2 * temp_ii + temp_jj;
						}

						if (q == 0)
						{
							for (l = 0; l < temp_n; l++)
							{
								info->Patch_array[own + l] = info->Patch_array[opp + l];
							}
							break;
						}
						else if (q == 1)
						{
							for (l = 0; l < temp_n; l++)
							{
								info->Patch_array[own + l] = info->Patch_array[opp + (temp_n - 1) - l];
							}
							break;
						}
					}
				}
			}

			// コネクティビティを作成
			int x, y;
			int point;

			own = 0;
			for (j = 0; j < i; j++)
			{
				int a, b;
				a = info->line_No_Total_element[j * info->DIMENSION + 0] * 2 + 1;
				b = info->line_No_Total_element[j * info->DIMENSION + 1] * 2 + 1;
				own += 2 * (a + b);
			}

			int e_x_max = info->line_No_Total_element[i * info->DIMENSION + 0];
			int e_y_max = info->line_No_Total_element[i * info->DIMENSION + 1];
			for (y = 0; y < e_y_max; y++)
			{
				for (x = 0; x < e_x_max; x++)
				{
					for (point = 0; point < 9; point++)
					{
						// point 0
						if (point == 0)
						{
							if (x == 0 && Edge[3] == 1)
							{
								int eta = 2 * y + p_y[point];
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Patch_array[own + 2 * ii + jj + eta];
							}
							else if (y == 0 && Edge[0] == 1)
							{
								int xi = 2 * x + p_x[point];
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Patch_array[own + xi];
							}
							else if (x > 0)
							{
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Connectivity[point_counter + y * e_x_max * 9 + (x - 1) * 9 + 1];
							}
							else if (y > 0)
							{
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Connectivity[point_counter + (y - 1) * e_x_max * 9 + x * 9 + 3];
							}
							else
							{
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = counter;
								info->Connectivity_ele[counter] = ele_counter + y * e_x_max + x;
								info->Connectivity_point[counter] = point;
								counter++;
							}
						}
						// point 1
						else if (point == 1)
						{
							if (x == e_x_max - 1 && Edge[1] == 1)
							{
								int eta = 2 * y + p_y[point];
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Patch_array[own + ii + eta];
							}
							else if (y == 0 && Edge[0] == 1)
							{
								int xi = 2 * x + p_x[point];
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Patch_array[own + xi];
							}
							else if (y > 0)
							{
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Connectivity[point_counter + (y - 1) * e_x_max * 9 + x * 9 + 2];
							}
							else
							{
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = counter;
								info->Connectivity_ele[counter] = ele_counter + y * e_x_max + x;
								info->Connectivity_point[counter] = point;
								counter++;
							}
						}
						// point 2
						else if (point == 2)
						{
							if (x == e_x_max - 1 && Edge[1] == 1)
							{
								int eta = 2 * y + p_y[point];
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Patch_array[own + ii + eta];
							}
							else if (y == e_y_max - 1 && Edge[2] == 1)
							{
								int xi = 2 * x + p_x[point];
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Patch_array[own + ii + jj + xi];
							}
							else
							{
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = counter;
								info->Connectivity_ele[counter] = ele_counter + y * e_x_max + x;
								info->Connectivity_point[counter] = point;
								counter++;
							}
						}
						// point 3
						else if (point == 3)
						{
							if (x == 0 && Edge[3] == 1)
							{
								int eta = 2 * y + p_y[point];
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Patch_array[own + 2 * ii + jj + eta];
							}
							else if (y == e_y_max - 1 && Edge[2] == 1)
							{
								int xi = 2 * x + p_x[point];
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Patch_array[own + ii + jj + xi];
							}
							else if (x > 0)
							{
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Connectivity[point_counter + y * e_x_max * 9 + (x - 1) * 9 + 2];
							}
							else
							{
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = counter;
								info->Connectivity_ele[counter] = ele_counter + y * e_x_max + x;
								info->Connectivity_point[counter] = point;
								counter++;
							}
						}
						// point 4
						else if (point == 4)
						{
							if (y == 0 && Edge[0] == 1)
							{
								int xi = 2 * x + p_x[point];
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Patch_array[own + xi];
							}
							else if (y > 0)
							{
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Connectivity[point_counter + (y - 1) * e_x_max * 9 + x * 9 + 6];
							}
							else
							{
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = counter;
								info->Connectivity_ele[counter] = ele_counter + y * e_x_max + x;
								info->Connectivity_point[counter] = point;
								counter++;
							}
						}
						// point 5
						else if (point == 5)
						{
							if (x == e_x_max - 1 && Edge[1] == 1)
							{
								int eta = 2 * y + p_y[point];
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Patch_array[own + ii + eta];
							}
							else
							{
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = counter;
								info->Connectivity_ele[counter] = ele_counter + y * e_x_max + x;
								info->Connectivity_point[counter] = point;
								counter++;
							}
						}
						// point 6
						else if (point == 6)
						{
							if (y == e_y_max - 1 && Edge[2] == 1)
							{
								int xi = 2 * x + p_x[point];
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Patch_array[own + ii + jj + xi];
							}
							else
							{
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = counter;
								info->Connectivity_ele[counter] = ele_counter + y * e_x_max + x;
								info->Connectivity_point[counter] = point;
								counter++;
							}
						}
						// point 7
						else if (point == 7)
						{
							if (x == 0 && Edge[3] == 1)
							{
								int eta = 2 * y + p_y[point];
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Patch_array[own + 2 * ii + jj + eta];
							}
							else if (x > 0)
							{
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = info->Connectivity[point_counter + y * e_x_max * 9 + (x - 1) * 9 + 5];
							}
							else
							{
								info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = counter;
								info->Connectivity_ele[counter] = ele_counter + y * e_x_max + x;
								info->Connectivity_point[counter] = point;
								counter++;
							}
						}
						// point 8
						else if (point == 8)
						{
							info->Connectivity[point_counter + y * e_x_max * 9 + x * 9 + point] = counter;
							info->Connectivity_ele[counter] = ele_counter + y * e_x_max + x;
							info->Connectivity_point[counter] = point;
							counter++;
						}
					}
				}
			}

			// Patch_array の作ってない分を作成
			int xi, eta;
			for (eta = 0; eta < jj; eta++)
			{
				for (xi = 0; xi < ii; xi++)
				{
					int e_x, e_y;

					if (eta == 0 && Edge[0] == 0)
					{
						Search_ele_point_2D(xi, eta, e_x_max, e_y_max, p_x, p_y, &e_x, &e_y, &point);
						info->Patch_array[own + xi] = info->Connectivity[point_counter + e_y * e_x_max * 9 + e_x * 9 + point];
					}
					if (eta == jj - 1 && Edge[2] == 0)
					{
						Search_ele_point_2D(xi, eta, e_x_max, e_y_max, p_x, p_y, &e_x, &e_y, &point);
						info->Patch_array[own + ii + jj + xi] = info->Connectivity[point_counter + e_y * e_x_max * 9 + e_x * 9 + point];
					}
					if (xi == 0 && Edge[3] == 0)
					{
						Search_ele_point_2D(xi, eta, e_x_max, e_y_max, p_x, p_y, &e_x, &e_y, &point);
						info->Patch_array[own + 2 * ii + jj + eta] = info->Connectivity[point_counter + e_y * e_x_max * 9 + e_x * 9 + point];
					}
					if (xi == ii - 1 && Edge[1] == 0)
					{
						Search_ele_point_2D(xi, eta, e_x_max, e_y_max, p_x, p_y, &e_x, &e_y, &point);
						info->Patch_array[own + ii + eta] = info->Connectivity[point_counter + e_y * e_x_max * 9 + e_x * 9 + point];
					}
				}
			}
			point_counter += e_x_max * e_y_max * 9;
			ele_counter += e_x_max * e_y_max;

			if (i == info->Total_Patch_to_mesh[1] - 1)
			{
				Total_connectivity_glo = counter;
			}
		}
		Total_connectivity = counter;
	}
	else if (info->DIMENSION == 3)
	{
		int p_x[27] = {0, 2, 2, 0, 0, 2, 2, 0, 1, 2, 1, 0, 1, 2, 1, 0, 0, 2, 2, 0, 0, 2, 1, 1, 1, 1, 1};
		int p_y[27] = {0, 0, 2, 2, 0, 0, 2, 2, 0, 1, 2, 1, 0, 1, 2, 1, 0, 0, 2, 2, 1, 1, 0, 2, 1, 1, 1};
		int p_z[27] = {0, 0, 0, 0, 2, 2, 2, 2, 0, 0, 0, 0, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 0, 2, 1};
		for (i = 0; i < info->Total_Patch_to_mesh[Total_mesh]; i++)
		{
			int own = 0;
			for (j = 0; j < i; j++)
			{
				int a, b, c;
				a = info->line_No_Total_element[j * info->DIMENSION + 0] * 2 + 1;
				b = info->line_No_Total_element[j * info->DIMENSION + 1] * 2 + 1;
				c = info->line_No_Total_element[j * info->DIMENSION + 2] * 2 + 1;
				own += 2 * (a * b + b * c + a * c);
			}

			int Face[6] = {0};
			int Face_counter = i * 36;
			int own_a = 0, own_b = 0;

			int ii = info->line_No_Total_element[i * info->DIMENSION + 0] * 2 + 1;
			int jj = info->line_No_Total_element[i * info->DIMENSION + 1] * 2 + 1;
			int kk = info->line_No_Total_element[i * info->DIMENSION + 2] * 2 + 1;

			// 重なっている面の Patch_array 配列を作成
			for (j = 0; j < 6; j++)
			{
				if (j == 0)
				{
					own_a = ii;
					own_b = jj;
				}
				else if (j == 1)
				{
					own += ii * jj;
					own_a = ii;
					own_b = jj;
				}
				else if (j == 2)
				{
					own += ii * jj;
					own_a = jj;
					own_b = kk;
				}
				else if (j == 3)
				{
					own += jj * kk;
					own_a = ii;
					own_b = kk;
				}
				else if (j == 4)
				{
					own += ii * kk;
					own_a = jj;
					own_b = kk;
				}
				else if (j == 5)
				{
					own += jj * kk;
					own_a = ii;
					own_b = kk;
				}

				for (k = 0; k < 6; k++)
				{
					if (info->Face_Edge_info[Face_counter + j * 6 + k] >= 0)
					{
						int opp = 0;
						for (l = 0; l < info->Opponent_patch_num[i * 6 + j]; l++)
						{
							int a, b, c;
							a = info->line_No_Total_element[l * info->DIMENSION + 0] * 2 + 1;
							b = info->line_No_Total_element[l * info->DIMENSION + 1] * 2 + 1;
							c = info->line_No_Total_element[l * info->DIMENSION + 2] * 2 + 1;
							opp += 2 * (a * b + b * c + a * c);
						}

						Face[j] = 1;

						int temp_ii = info->line_No_Total_element[info->Opponent_patch_num[i * 6 + j] * info->DIMENSION + 0] * 2 + 1;
						int temp_jj = info->line_No_Total_element[info->Opponent_patch_num[i * 6 + j] * info->DIMENSION + 1] * 2 + 1;
						int temp_kk = info->line_No_Total_element[info->Opponent_patch_num[i * 6 + j] * info->DIMENSION + 2] * 2 + 1;

						for (l = 1; l <= k; l++)
						{
							if (l == 1)
							{
								opp += temp_ii * temp_jj;
							}
							else if (l == 2)
							{
								opp += temp_ii * temp_jj;
							}
							else if (l == 3)
							{
								opp += temp_jj * temp_kk;
							}
							else if (l == 4)
							{
								opp += temp_ii * temp_kk;
							}
							else if (l == 5)
							{
								opp += temp_jj * temp_kk;
							}
						}
						
						if (info->Face_Edge_info[Face_counter + j * 6 + k] == 0)
						{
							for (l = 0; l < own_b; l++)
							{
								for (m = 0; m < own_a; m++)
								{
									info->Patch_array[own + l * own_a + m] = info->Patch_array[opp + l * own_a + m];
								}
							}
							break;
						}
						else if (info->Face_Edge_info[Face_counter + j * 6 + k]  == 1)
						{
							for (l = 0; l < own_b; l++)
							{
								for (m = 0; m < own_a; m++)
								{
									info->Patch_array[own + l * own_a + m] = info->Patch_array[opp + l * own_a + ((own_a - 1) - m)];
								}
							}
							break;
						}
						else if (info->Face_Edge_info[Face_counter + j * 6 + k]  == 2)
						{
							for (l = 0; l < own_b; l++)
							{
								for (m = 0; m < own_a; m++)
								{
									info->Patch_array[own + l * own_a + m] = info->Patch_array[opp + ((own_b - 1) - l) * own_a + ((own_a - 1) - m)];
								}
							}
							break;
						}
						else if (info->Face_Edge_info[Face_counter + j * 6 + k]  == 3)
						{
							for (l = 0; l < own_b; l++)
							{
								for (m = 0; m < own_a; m++)
								{
									info->Patch_array[own + l * own_a + m] = info->Patch_array[opp + ((own_b - 1) - l) * own_a + m];
								}
							}
							break;
						}
					}
				}
			}

			// コネクティビティを作成
			int x, y, z;
			int point;

			own = 0;
			for (j = 0; j < i; j++)
			{
				int a, b, c;
				a = info->line_No_Total_element[j * info->DIMENSION + 0] * 2 + 1;
				b = info->line_No_Total_element[j * info->DIMENSION + 1] * 2 + 1;
				c = info->line_No_Total_element[j * info->DIMENSION + 2] * 2 + 1;
				own += 2 * (a * b + b * c + a * c);
			}

			int temp1, temp2, temp3, temp4, temp5;
			temp1 = ii * jj;
			temp2 = temp1 + ii * jj;
			temp3 = temp2 + jj * kk;
			temp4 = temp3 + ii * kk;
			temp5 = temp4 + jj * kk;

			int e_x_max = info->line_No_Total_element[i * info->DIMENSION + 0];
			int e_y_max = info->line_No_Total_element[i * info->DIMENSION + 1];
			int e_z_max = info->line_No_Total_element[i * info->DIMENSION + 2];
			for (z = 0; z < e_z_max; z++)
			{
				for (y = 0; y < e_y_max; y++)
				{
					for (x = 0; x < e_x_max; x++)
					{
						for (point = 0; point < 27; point++)
						{
							int temp = point_counter + z * e_y_max * e_x_max * 27 + y * e_x_max * 27 + x * 27 + point;

							// point 0
							if (point == 0)
							{
								if (x == 0 && Face[2] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp2 + zeta * jj + eta];
								}
								else if (y == 0 && Face[3] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp3 + zeta * ii + xi];
								}
								else if (z == 0 && Face[0] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + eta * ii + xi];
								}
								else if (x > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + y * e_x_max * 27 + (x - 1) * 27 + 1];
								}
								else if (y > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + (y - 1) * e_x_max * 27 + x * 27 + 3];
								}
								else if (z > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + (z - 1) * e_y_max * e_x_max * 27 + y * e_x_max * 27 + x * 27 + 4];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 1
							else if (point == 1)
							{
								if (x == e_x_max - 1 && Face[4] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp4 + zeta * jj + eta];
								}
								else if (y == 0 && Face[3] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp3 + zeta * ii + xi];
								}
								else if (z == 0 && Face[0] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + eta * ii + xi];
								}
								else if (y > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + (y - 1) * e_x_max * 27 + x * 27 + 2];
								}
								else if (z > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + (z - 1) * e_y_max * e_x_max * 27 + y * e_x_max * 27 + x * 27 + 5];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 2
							else if (point == 2)
							{
								if (x == e_x_max - 1 && Face[4] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp4 + zeta * jj + eta];
								}
								else if (y == e_y_max - 1 && Face[5] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp5 + zeta * ii + xi];
								}
								else if (z == 0 && Face[0] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + eta * ii + xi];
								}
								else if (z > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + (z - 1) * e_y_max * e_x_max * 27 + y * e_x_max * 27 + x * 27 + 6];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 3
							else if (point == 3)
							{
								if (x == 0 && Face[2] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp2 + zeta * jj + eta];
								}
								else if (y == e_y_max - 1 && Face[5] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp5 + zeta * ii + xi];
								}
								else if (z == 0 && Face[0] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + eta * ii + xi];
								}
								else if (x > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + y * e_x_max * 27 + (x - 1) * 27 + 2];
								}
								else if (z > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + (z - 1) * e_y_max * e_x_max * 27 + y * e_x_max * 27 + x * 27 + 7];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 4
							else if (point == 4)
							{
								if (x == 0 && Face[2] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp2 + zeta * jj + eta];
								}
								else if (y == 0 && Face[3] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp3 + zeta * ii + xi];
								}
								else if (z == e_z_max - 1 && Face[1] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + temp1 + eta * ii + xi];
								}
								else if (x > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + y * e_x_max * 27 + (x - 1) * 27 + 5];
								}
								else if (y > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + (y - 1) * e_x_max * 27 + x * 27 + 7];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 5
							else if (point == 5)
							{
								if (x == e_x_max - 1 && Face[4] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp4 + zeta * jj + eta];
								}
								else if (y == 0 && Face[3] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp3 + zeta * ii + xi];
								}
								else if (z == e_z_max - 1 && Face[1] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + temp1 + eta * ii + xi];
								}
								else if (y > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + (y - 1) * e_x_max * 27 + x * 27 + 6];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 6
							else if (point == 6)
							{
								if (x == e_x_max - 1 && Face[4] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp4 + zeta * jj + eta];
								}
								else if (y == e_y_max - 1 && Face[5] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp5 + zeta * ii + xi];
								}
								else if (z == e_z_max - 1 && Face[1] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + temp1 + eta * ii + xi];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 7
							else if (point == 7)
							{
								if (x == 0 && Face[2] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp2 + zeta * jj + eta];
								}
								else if (y == e_y_max - 1 && Face[5] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp5 + zeta * ii + xi];
								}
								else if (z == e_z_max - 1 && Face[1] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + temp1 + eta * ii + xi];
								}
								else if (x > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + y * e_x_max * 27 + (x - 1) * 27 + 6];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 8
							else if (point == 8)
							{
								if (y == 0 && Face[3] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp3 + zeta * ii + xi];
								}
								else if (z == 0 && Face[0] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + eta * ii + xi];
								}
								else if (y > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + (y - 1) * e_x_max * 27 + x * 27 + 10];
								}
								else if (z > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + (z - 1) * e_y_max * e_x_max * 27 + y * e_x_max * 27 + x * 27 + 12];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 9
							else if (point == 9)
							{
								if (x == e_x_max - 1 && Face[4] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp4 + zeta * jj + eta];
								}
								else if (z == 0 && Face[0] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + eta * ii + xi];
								}
								else if (z > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + (z - 1) * e_y_max * e_x_max * 27 + y * e_x_max * 27 + x * 27 + 13];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 10
							else if (point == 10)
							{
								if (y == e_y_max - 1 && Face[5] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp5 + zeta * ii + xi];
								}
								else if (z == 0 && Face[0] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + eta * ii + xi];
								}
								else if (z > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + (z - 1) * e_y_max * e_x_max * 27 + y * e_x_max * 27 + x * 27 + 14];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 11
							else if (point == 11)
							{
								if (x == 0 && Face[2] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp2 + zeta * jj + eta];
								}
								else if (z == 0 && Face[0] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + eta * ii + xi];
								}
								else if (x > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + y * e_x_max * 27 + (x - 1) * 27 + 9];
								}
								else if (z > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + (z - 1) * e_y_max * e_x_max * 27 + y * e_x_max * 27 + x * 27 + 15];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 12
							else if (point == 12)
							{
								if (y == 0 && Face[3] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp3 + zeta * ii + xi];
								}
								else if (z == e_z_max - 1 && Face[1] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + temp1 + eta * ii + xi];
								}
								else if (y > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + (y - 1) * e_x_max * 27 + x * 27 + 14];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 13
							else if (point == 13)
							{
								if (x == e_x_max - 1 && Face[4] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp4 + zeta * jj + eta];
								}
								else if (z == e_z_max - 1 && Face[1] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + temp1 + eta * ii + xi];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 14
							else if (point == 14)
							{
								if (y == e_y_max - 1 && Face[5] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp5 + zeta * ii + xi];
								}
								else if (z == e_z_max - 1 && Face[1] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + temp1 + eta * ii + xi];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 15
							else if (point == 15)
							{
								if (x == 0 && Face[2] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp2 + zeta * jj + eta];
								}
								else if (z == e_z_max - 1 && Face[1] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + temp1 + eta * ii + xi];
								}
								else if (x > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + y * e_x_max * 27 + (x - 1) * 27 + 13];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 16
							else if (point == 16)
							{
								if (x == 0 && Face[2] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp2 + zeta * jj + eta];
								}
								else if (y == 0 && Face[3] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp3 + zeta * ii + xi];
								}
								else if (x > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + y * e_x_max * 27 + (x - 1) * 27 + 17];
								}
								else if (y > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + (y - 1) * e_x_max * 27 + x * 27 + 19];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 17
							else if (point == 17)
							{
								if (x == e_x_max - 1 && Face[4] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp4 + zeta * jj + eta];
								}
								else if (y == 0 && Face[3] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp3 + zeta * ii + xi];
								}
								else if (y > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + (y - 1) * e_x_max * 27 + x * 27 + 18];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 18
							else if (point == 18)
							{
								if (x == e_x_max - 1 && Face[4] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp4 + zeta * jj + eta];
								}
								else if (y == e_y_max - 1 && Face[5] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp5 + zeta * ii + xi];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 19
							else if (point == 19)
							{
								if (x == 0 && Face[2] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp2 + zeta * jj + eta];
								}
								else if (y == e_y_max - 1 && Face[5] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp5 + zeta * ii + xi];
								}
								else if (x > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + y * e_x_max * 27 + (x - 1) * 27 + 18];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 20
							else if (point == 20)
							{
								if (x == 0 && Face[2] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp2 + zeta * jj + eta];
								}
								else if (x > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + y * e_x_max * 27 + (x - 1) * 27 + 21];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 21
							else if (point == 21)
							{
								if (x == e_x_max - 1 && Face[4] == 1)
								{
									int eta = 2 * y + p_y[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp4 + zeta * jj + eta];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 22
							else if (point == 22)
							{
								if (y == 0 && Face[3] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp3 + zeta * ii + xi];
								}
								else if (y > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + z * e_y_max * e_x_max * 27 + (y - 1) * e_x_max * 27 + x * 27 + 23];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 23
							else if (point == 23)
							{
								if (y == e_y_max - 1 && Face[5] == 1)
								{
									int xi = 2 * x + p_x[point];
									int zeta = 2 * z + p_z[point];
									info->Connectivity[temp] = info->Patch_array[own + temp5 + zeta * ii + xi];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 24
							else if (point == 24)
							{
								if (z == 0 && Face[0] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + eta * ii + xi];
								}
								else if (z > 0)
								{
									info->Connectivity[temp] = info->Connectivity[point_counter + (z - 1) * e_y_max * e_x_max * 27 + y * e_x_max * 27 + x * 27 + 25];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 25
							else if (point == 25)
							{
								if (z == e_z_max - 1 && Face[1] == 1)
								{
									int xi = 2 * x + p_x[point];
									int eta = 2 * y + p_y[point];
									info->Connectivity[temp] = info->Patch_array[own + temp1 + eta * ii + xi];
								}
								else
								{
									info->Connectivity[temp] = counter;
									info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
									info->Connectivity_point[counter] = point;
									counter++;
								}
							}
							// point 26
							else if (point == 26)
							{
								info->Connectivity[temp] = counter;
								info->Connectivity_ele[counter] = ele_counter + z * e_y_max * e_x_max + y * e_x_max + x;
								info->Connectivity_point[counter] = point;
								counter++;
							}
						}
					}
				}
			}

			// Patch_array の作ってない分を作成
			int xi, eta, zeta;
			for (zeta = 0; zeta < kk; zeta++)
			{
				for (eta = 0; eta < jj; eta++)
				{
					for (xi = 0; xi < ii; xi++)
					{
						int e_x, e_y, e_z;

						if (zeta == 0 && Face[0] == 0)
						{
							Search_ele_point_3D(xi, eta, zeta, e_x_max, e_y_max, e_z_max, p_x, p_y, p_z, &e_x, &e_y, &e_z, &point);
							info->Patch_array[own + eta * ii + xi] = info->Connectivity[point_counter + e_z * e_y_max * e_x_max * 27 + e_y * e_x_max * 27 + e_x * 27 + point];
						}
						if (zeta == kk - 1 && Face[1] == 0)
						{
							Search_ele_point_3D(xi, eta, zeta, e_x_max, e_y_max, e_z_max, p_x, p_y, p_z, &e_x, &e_y, &e_z, &point);
							info->Patch_array[own + temp1 + eta * ii + xi] = info->Connectivity[point_counter + e_z * e_y_max * e_x_max * 27 + e_y * e_x_max * 27 + e_x * 27 + point];
						}
						if (eta == 0 && Face[3] == 0)
						{
							Search_ele_point_3D(xi, eta, zeta, e_x_max, e_y_max, e_z_max, p_x, p_y, p_z, &e_x, &e_y, &e_z, &point);
							info->Patch_array[own + temp3 + zeta * ii + xi] = info->Connectivity[point_counter + e_z * e_y_max * e_x_max * 27 + e_y * e_x_max * 27 + e_x * 27 + point];
						}
						if (eta == jj - 1 && Face[5] == 0)
						{
							Search_ele_point_3D(xi, eta, zeta, e_x_max, e_y_max, e_z_max, p_x, p_y, p_z, &e_x, &e_y, &e_z, &point);
							info->Patch_array[own + temp5 + zeta * ii + xi] = info->Connectivity[point_counter + e_z * e_y_max * e_x_max * 27 + e_y * e_x_max * 27 + e_x * 27 + point];
						}
						if (xi == 0 && Face[2] == 0)
						{
							Search_ele_point_3D(xi, eta, zeta, e_x_max, e_y_max, e_z_max, p_x, p_y, p_z, &e_x, &e_y, &e_z, &point);
							info->Patch_array[own + temp2 + zeta * jj + eta] = info->Connectivity[point_counter + e_z * e_y_max * e_x_max * 27 + e_y * e_x_max * 27 + e_x * 27 + point];
						}
						if (xi == ii - 1 && Face[4] == 0)
						{
							Search_ele_point_3D(xi, eta, zeta, e_x_max, e_y_max, e_z_max, p_x, p_y, p_z, &e_x, &e_y, &e_z, &point);
							info->Patch_array[own + temp4 + zeta * jj + eta] = info->Connectivity[point_counter + e_z * e_y_max * e_x_max * 27 + e_y * e_x_max * 27 + e_x * 27 + point];
						}
					}
				}
			}
			point_counter += e_x_max * e_y_max * e_z_max * 27;
			ele_counter += e_x_max * e_y_max * e_z_max;

			if (i == info->Total_Patch_to_mesh[1] - 1)
			{
				Total_connectivity_glo = counter;
			}
		}
		Total_connectivity = counter;
	}
}


void Search_ele_point_2D(int xi, int eta, int e_x_max, int e_y_max, int *p_x, int *p_y, int *e_x, int *e_y, int *point)
{
	int i;

	*e_x = xi / 2;
	*e_y = eta / 2;
	int p_num_x = xi % 2;
	int p_num_y = eta % 2;

	// 例外
	if (*e_x == e_x_max)
	{
		(*e_x) -= 1;
		p_num_x = 2;
	}
	if (*e_y == e_y_max)
	{
		(*e_y) -= 1;
		p_num_y = 2;
	}

	// ポイント番号を探索
	for (i = 0; i < 9 ; i++)
	{
		if (p_x[i] == p_num_x && p_y[i] == p_num_y)
		{
			*point = i;
			return;
		}
	}
}


void Search_ele_point_3D(int xi, int eta, int zeta, int e_x_max, int e_y_max, int e_z_max, int *p_x, int *p_y, int *p_z, int *e_x, int *e_y, int *e_z, int *point)
{
	int i;

	*e_x = xi / 2;
	*e_y = eta / 2;
	*e_z = zeta / 2;
	int p_num_x = xi % 2;
	int p_num_y = eta % 2;
	int p_num_z = zeta % 2;

	// 例外
	if (*e_x == e_x_max)
	{
		(*e_x) -= 1;
		p_num_x = 2;
	}
	if (*e_y == e_y_max)
	{
		(*e_y) -= 1;
		p_num_y = 2;
	}
	if (*e_z == e_z_max)
	{
		(*e_z) -= 1;
		p_num_z = 2;
	}

	// ポイント番号を探索
	for (i = 0; i < 27 ; i++)
	{
		if (p_x[i] == p_num_x && p_y[i] == p_num_y && p_z[i] == p_num_z)
		{
			*point = i;
			return;
		}
	}
}


void Make_info_for_viewer(information *info)
{
	int i, j, k;
	int point_on_element = pow(3, info->DIMENSION);
	double *point_array = (double *)malloc(sizeof(double) * point_on_element * info->DIMENSION);

	// Make point in element
	int counter = 0;
	if (info->DIMENSION == 2)
	{
		// 0, 2, 8, 6, 1, 5, 7, 3, 4

		// point 0
		point_array[counter] = -1.0;		point_array[counter + 1] = -1.0;		counter += 2;
		// point 1
		point_array[counter] =  1.0;		point_array[counter + 1] = -1.0;		counter += 2;
		// point 2
		point_array[counter] =  1.0;		point_array[counter + 1] =  1.0;		counter += 2;
		// point 3
		point_array[counter] = -1.0;		point_array[counter + 1] =  1.0;		counter += 2;
		// point 4
		point_array[counter] =  0.0;		point_array[counter + 1] = -1.0;		counter += 2;
		// point 5
		point_array[counter] =  1.0;		point_array[counter + 1] =  0.0;		counter += 2;
		// point 6
		point_array[counter] =  0.0;		point_array[counter + 1] =  1.0;		counter += 2;
		// point 7
		point_array[counter] = -1.0;		point_array[counter + 1] =  0.0;		counter += 2;
		// point 8
		point_array[counter] =  0.0;		point_array[counter + 1] =  0.0;
	}
	else if (info->DIMENSION == 3)
	{
		// 0, 2, 8, 6, 18, 20, 26, 24, 1, 5, 7, 3, 19, 23, 25, 21, 9, 11, 17, 15, 12, 14, 10, 16, 4, 22, 13

		// point 0
		point_array[counter] = -1.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 1
		point_array[counter] =  1.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 2
		point_array[counter] =  1.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 3
		point_array[counter] = -1.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 4
		point_array[counter] = -1.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 5
		point_array[counter] =  1.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 6
		point_array[counter] =  1.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 7
		point_array[counter] = -1.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 8
		point_array[counter] =  0.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 9
		point_array[counter] =  1.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 10
		point_array[counter] =  0.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 11
		point_array[counter] = -1.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 12
		point_array[counter] =  0.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 13
		point_array[counter] =  1.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 14
		point_array[counter] =  0.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 15
		point_array[counter] = -1.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 16
		point_array[counter] = -1.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 17
		point_array[counter] =  1.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 18
		point_array[counter] =  1.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 19
		point_array[counter] = -1.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 20
		point_array[counter] = -1.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 21
		point_array[counter] =  1.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 22
		point_array[counter] =  0.0;		point_array[counter + 1] = -1.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 23
		point_array[counter] =  0.0;		point_array[counter + 1] =  1.0;		point_array[counter + 2] =  0.0;		counter += 3;
		// point 24
		point_array[counter] =  0.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] = -1.0;		counter += 3;
		// point 25
		point_array[counter] =  0.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] =  1.0;		counter += 3;
		// point 26
		point_array[counter] =  0.0;		point_array[counter + 1] =  0.0;		point_array[counter + 2] =  0.0;
	}

	// calculation physical coordinate, B matrix, strain, stress
	double temp_point[MAX_DIMENSION];
	double temp_point_glo[MAX_DIMENSION];
	double temp_para_glo[MAX_DIMENSION];

	double *B = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * MAX_KIEL_SIZE);
	double *BG = (double *)malloc(sizeof(double) * D_MATRIX_SIZE * MAX_KIEL_SIZE);
	double *temp_disp = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT * info->DIMENSION);
	double *temp_disp_glo = (double *)malloc(sizeof(double) * MAX_NO_CP_ON_ELEMENT * info->DIMENSION);

	int *temp_element_n = (int *)malloc(sizeof(int) * MAX_N_ELEMENT_OVER_POINT);
	int *temp_ad = (int *)malloc(sizeof(int) * info->DIMENSION * (MAX_ORDER + 1));

	for (i = 0; i < Total_connectivity; i++)
	{
		int element = info->Connectivity_ele[i], element_glo = 0;
		int point = info->Connectivity_point[i];

		// make temp_point
		for (j = 0; j < info->DIMENSION; j++)
		{
			temp_point[j] = point_array[point * info->DIMENSION + j];
		}

		// B matrix
		Make_B_Matrix_anypoint(element, B, temp_point, info);

		for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[element]]; j++)
		{
			double R = Shape_func(j, temp_point, element, info);

			// make physical coordinate, disp_at_connectivity, temp_disp
			for (k = 0; k < info->DIMENSION; k++)
			{
				double d = info->Displacement[info->Controlpoint_of_Element[element * MAX_NO_CP_ON_ELEMENT + j] * info->DIMENSION + k];
				info->Connectivity_coord[i * info->DIMENSION + k] += R * info->Node_Coordinate[info->Controlpoint_of_Element[element * MAX_NO_CP_ON_ELEMENT + j] * (info->DIMENSION + 1) + k];
				info->disp_at_connectivity[i * info->DIMENSION + k] += R * d;
				temp_disp[j * info->DIMENSION + k] = d;
			}
		}

		// overlay
		if (i >= Total_connectivity_glo)
		{
			// make temp_point_glo
			if (info->DIMENSION == 2)
			{
				double *temp_Position_Knots_xi = (double *)malloc(sizeof(double) * info->No_knot[0 * info->DIMENSION + 0]);
				double *temp_Position_Knots_eta = (double *)malloc(sizeof(double) * info->No_knot[0 * info->DIMENSION + 1]);
				for (j = 0; j < info->No_knot[0 * info->DIMENSION + 0]; j++)
				{
					temp_Position_Knots_xi[j] = info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 0] + j];
				}
				for (j = 0; j < info->No_knot[0 * info->DIMENSION + 1]; j++)
				{
					temp_Position_Knots_eta[j] = info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 1] + j];
				}

				Calc_xi_eta(info->Connectivity_coord[i * info->DIMENSION + 0], info->Connectivity_coord[i * info->DIMENSION + 1],
							temp_Position_Knots_xi, temp_Position_Knots_eta,
							info->No_Control_point[0 * info->DIMENSION + 0], info->No_Control_point[0 * info->DIMENSION + 1],
							info->Order[0 * info->DIMENSION + 0], info->Order[0 * info->DIMENSION + 1],
							&temp_para_glo[0], &temp_para_glo[1], info);

				ele_check(0, temp_para_glo, temp_element_n, temp_ad, info);
				element_glo = temp_element_n[0];

				// 親要素座標の算出
				temp_point_glo[0] = - 1.0 + 2.0 * (temp_para_glo[0] - info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 0] + info->Order[0 * info->DIMENSION + 0] + info->ENC[element_glo * info->DIMENSION + 0]])
								  / (info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 0] + info->Order[0 * info->DIMENSION + 0] + info->ENC[element_glo * info->DIMENSION + 0] + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 0] + info->Order[0 * info->DIMENSION + 0] + info->ENC[element_glo * info->DIMENSION + 0]]);
				temp_point_glo[1] = - 1.0 + 2.0 * (temp_para_glo[1] - info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 1] + info->Order[0 * info->DIMENSION + 1] + info->ENC[element_glo * info->DIMENSION + 1]])
								  / (info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 1] + info->Order[0 * info->DIMENSION + 1] + info->ENC[element_glo * info->DIMENSION + 1] + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 1] + info->Order[0 * info->DIMENSION + 1] + info->ENC[element_glo * info->DIMENSION + 1]]);

			}
			else if (info->DIMENSION == 3)
			{
				double *temp_Position_Knots_xi = (double *)malloc(sizeof(double) * info->No_knot[0 * info->DIMENSION + 0]);
				double *temp_Position_Knots_eta = (double *)malloc(sizeof(double) * info->No_knot[0 * info->DIMENSION + 1]);
				double *temp_Position_Knots_zeta = (double *)malloc(sizeof(double) * info->No_knot[0 * info->DIMENSION + 2]);

				for (j = 0; j < info->No_knot[0 * info->DIMENSION + 0]; j++)
				{
					temp_Position_Knots_xi[j] = info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 0] + j];
				}
				for (j = 0; j < info->No_knot[0 * info->DIMENSION + 1]; j++)
				{
					temp_Position_Knots_eta[j] = info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 1] + j];
				}
				for (j = 0; j < info->No_knot[0 * info->DIMENSION + 2]; j++)
				{
					temp_Position_Knots_zeta[j] = info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 2] + j];
				}

				Calc_xi_eta_zeta(info->Connectivity_coord[i * info->DIMENSION + 0], info->Connectivity_coord[i * info->DIMENSION + 1], info->Connectivity_coord[i * info->DIMENSION + 2],
								 temp_Position_Knots_xi, temp_Position_Knots_eta, temp_Position_Knots_zeta,
								 info->No_Control_point[0 * info->DIMENSION + 0], info->No_Control_point[0 * info->DIMENSION + 1], info->No_Control_point[0 * info->DIMENSION + 2],
								 info->Order[0 * info->DIMENSION + 0], info->Order[0 * info->DIMENSION + 1], info->Order[0 * info->DIMENSION + 2],
								 &temp_para_glo[0], &temp_para_glo[1], &temp_para_glo[2], info);

				ele_check(0, temp_para_glo, temp_element_n, temp_ad, info);
				element_glo = temp_element_n[0];

				// 親要素座標の算出
				temp_point_glo[0] = - 1.0 + 2.0 * (temp_para_glo[0] - info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 0] + info->Order[0 * info->DIMENSION + 0] + info->ENC[element_glo * info->DIMENSION + 0]])
						 		  / (info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 0] + info->Order[0 * info->DIMENSION + 0] + info->ENC[element_glo * info->DIMENSION + 0] + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 0] + info->Order[0 * info->DIMENSION + 0] + info->ENC[element_glo * info->DIMENSION + 0]]);
				temp_point_glo[1] = - 1.0 + 2.0 * (temp_para_glo[1] - info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 1] + info->Order[0 * info->DIMENSION + 1] + info->ENC[element_glo * info->DIMENSION + 1]])
								  / (info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 1] + info->Order[0 * info->DIMENSION + 1] + info->ENC[element_glo * info->DIMENSION + 1] + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 1] + info->Order[0 * info->DIMENSION + 1] + info->ENC[element_glo * info->DIMENSION + 1]]);
				temp_point_glo[2] = - 1.0 + 2.0 * (temp_para_glo[2] - info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 2] + info->Order[0 * info->DIMENSION + 2] + info->ENC[element_glo * info->DIMENSION + 2]])
								  / (info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 2] + info->Order[0 * info->DIMENSION + 2] + info->ENC[element_glo * info->DIMENSION + 2] + 1] - info->Position_Knots[info->Total_Knot_to_patch_dim[0 * info->DIMENSION + 2] + info->Order[0 * info->DIMENSION + 2] + info->ENC[element_glo * info->DIMENSION + 2]]);
			}

			// BG matrix
			Make_B_Matrix_anypoint(element_glo, BG, temp_point_glo, info);

			for (j = 0; j < info->No_Control_point_ON_ELEMENT[info->Element_patch[element_glo]]; j++)
			{
				double R_glo = Shape_func(j, temp_point_glo, element_glo, info);

				// overlay displacement, make temp_disp_glo
				for (k = 0; k < info->DIMENSION; k++)
				{
					double d_glo = info->Displacement[info->Controlpoint_of_Element[element_glo * MAX_NO_CP_ON_ELEMENT + j] * info->DIMENSION + k];
					info->disp_at_connectivity[i * info->DIMENSION + k] += R_glo * d_glo;
					temp_disp_glo[j * info->DIMENSION + k] = d_glo;
				}
			}
		}

		// strain
		int KIEL_SIZE = info->No_Control_point_ON_ELEMENT[info->Element_patch[element]] * info->DIMENSION;
		for (j = 0; j < D_MATRIX_SIZE; j++)
		{
			for (k = 0; k < KIEL_SIZE; k++)
			{
				info->strain_at_connectivity[i * N_STRAIN + j] += B[j * MAX_KIEL_SIZE + k] * temp_disp[k];
			}
		}

		// overlay strain
		if (i >= Total_connectivity_glo)
		{
			int KIEL_SIZE_glo = info->No_Control_point_ON_ELEMENT[info->Element_patch[element_glo]] * info->DIMENSION;
			for (j = 0; j < D_MATRIX_SIZE; j++)
			{
				for (k = 0; k < KIEL_SIZE_glo; k++)
				{
					info->strain_at_connectivity[i * N_STRAIN + j] += BG[j * MAX_KIEL_SIZE + k] * temp_disp_glo[k];
				}
			}
		}

		// stress
		for (j = 0; j < D_MATRIX_SIZE; j++)
		{
			for (k = 0; k < D_MATRIX_SIZE; k++)
			{
				info->stress_at_connectivity[i * N_STRESS + j] += info->D[j * D_MATRIX_SIZE + k] * info->strain_at_connectivity[i * N_STRAIN + k];
			}
		}

		// if DIMENSION == 2, make strain zz
		if (info->DIMENSION == 2)
		{
			// 平面応力状態
			if (DM == 0)
			{
				info->strain_at_connectivity[i * N_STRAIN + 3] = - 1.0 * nu / E * (info->stress_at_connectivity[i * N_STRESS + 0] + info->stress_at_connectivity[i * N_STRESS + 1]);
			}
			// 平面ひずみ状態
			else if (DM == 1)
			{
				info->stress_at_connectivity[i * N_STRESS + 3] = E * nu / (1.0 + nu) / (1.0 - 2.0 * nu) * (info->strain_at_connectivity[i * N_STRAIN + 0] + info->strain_at_connectivity[i * N_STRAIN + 1]);
			}
		}

		// von Mises stress

	}

	free(point_array), free(B), free(BG), free(temp_disp), free(temp_disp_glo), free(temp_element_n), free(temp_ad);
}