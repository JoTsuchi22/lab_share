#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#include "../NURBS_view/constant.h"
#include "NURBS_calc.h"

int stress_type_flag = 1;			//0:平面応力 1:平面ひずみ
double E;							//ヤング率(GPa)
double nu;							//ポアソン比(-)
int patch_n;						//パッチ数
int cntl_p_n;						//コントロールポイント数
int order_xi[MAX_PATCHES];		//ξ基底関数の次数(p)
int order_eta[MAX_PATCHES];		//η基底関数の次数(p)
int knot_n_xi[MAX_PATCHES];		//ξノットベクトルの数(n+p+1)
int knot_n_eta[MAX_PATCHES];		//ηノットベクトルの数(n+p+1)
int cntl_p_n_xi[MAX_PATCHES];	//ξ方向コントロールポイント数(n)
int cntl_p_n_eta[MAX_PATCHES];	//η方向コントロールポイント数(n)

double knot_vec_xi[MAX_PATCHES][MAX_KNOTS];		//ξノットベクトル
double knot_vec_eta[MAX_PATCHES][MAX_KNOTS];	//ηノットベクトル
double cntl_px[MAX_PATCHES][MAX_CNRL_P];		//コントロールポイントx座標
double cntl_py[MAX_PATCHES][MAX_CNRL_P];		//コントロールポイントy座標
double disp_cntl_px[MAX_PATCHES][MAX_CNRL_P];	//コントロールポイント上のx方向変位
double disp_cntl_py[MAX_PATCHES][MAX_CNRL_P];	//コントロールポイント上のy方向変位
double weight[MAX_PATCHES][MAX_CNRL_P];			//重み

double coord_x[MAX_POINTS][MAX_POINTS];		//メッシュx座標
double coord_y[MAX_POINTS][MAX_POINTS];		//メッシュy座標
double dxi_x[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];		// ∂x/∂ξ
double dxi_y[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];		// ∂y/∂ξ
double deta_x[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];		// ∂x/∂η
double deta_y[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];		// ∂y/∂η

double disp_x[MAX_POINTS][MAX_POINTS];			//x方向変位
double disp_y[MAX_POINTS][MAX_POINTS];			//y方向変位
double dxi_disp_x[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	// ∂u/∂ξ
double dxi_disp_y[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	// ∂v/∂ξ
double deta_disp_x[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	// ∂u/∂η
double deta_disp_y[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	// ∂v/∂η

double strain_xx[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//x方向ひずみ
double strain_yy[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//y方向ひずみ
double strain_xy[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//剪断ひずみ

double stress_xx[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//x方向垂直応力
double stress_yy[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//y方向垂直応力
double stress_xy[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//剪断応力

int fields_flag;
int division_ele_xi;		//ξ方向の一要素あたりの分割数
int division_ele_eta;		//η方向の一要素あたりの分割数
int division_n_xi;		//ξ方向の表示する点の数
int division_n_eta;		//η方向の表示する点の数
int element_n_xi;			//ξ方向要素数
int element_n_eta;		//η方向要素数

int temp_index[MAX_PATCHES][MAX_CNRL_P];
double temp_cntl_px[MAX_CNRL_P];
double temp_cntl_py[MAX_CNRL_P];
double temp_weight[MAX_CNRL_P];
double temp_disp_x[MAX_CNRL_P];
double temp_disp_y[MAX_CNRL_P];

FILE *fp;	// FILE型構造体

void ReadFile (char *filename_input, char *filename_disp) {
	int i, j;
	double temp1, temp2, temp3;
	int temp_int;

	printf("Start Reading input\n\n");
	fp = fopen(filename_input, "r");

	fscanf(fp, "%lf%lf", &E, &nu);
	fscanf(fp, "\n");
	printf("E nu: % 1.4e % 1.4e\n", E, nu);

	fscanf(fp, "%d", &patch_n);
	fscanf(fp, "\n");
	printf("patches: %d \n", patch_n);
	if (patch_n > MAX_PATCHES) {
		printf("Error!!\n");
		printf("Too many patches!\n"
			   "Maximum of patches is %d (Now %d)\n"
			   "\n", MAX_PATCHES, patch_n);
		exit(1);
	}

	fscanf(fp, "%d", &cntl_p_n);
	fscanf(fp, "\n");
	printf("total control points:%d \n", cntl_p_n);
	if (cntl_p_n > MAX_CNRL_P) {
		printf("Error!!\n");
		printf("Too many control points!\n"
			   "Maximum of control points is %d (Now %d)\n"
			   "\n", MAX_CNRL_P, cntl_p_n);
		exit(1);
	}

	for (i = 0; i < patch_n; i++) {
		fscanf(fp, "%d%d", &order_xi[i], &order_eta[i]);
		fscanf(fp, "\n");
		printf("order %d: %d %d\n", i, order_xi[i], order_eta[i]);
		if (order_xi[i] > MAX_ORDER) {
			printf("Error!!\n");
			printf("Order too big at xi!\n"
				   "Maximum of order is %d (Now %d at patch %d)\n"
				   "\n", MAX_ORDER, order_xi[i], i);
			exit(1);
		}
		if (order_eta[i] > MAX_ORDER) {
			printf("Error!!\n");
			printf("Order too big at eta!\n"
				   "Maximum of order is %d (Now %d at patch %d)\n"
				   "\n", MAX_ORDER, order_eta[i], i);
			exit(1);
		}
	}

	for (i = 0; i < patch_n; i++) {
		fscanf(fp, "%d%d", &knot_n_xi[i], &knot_n_eta[i]);
		fscanf(fp, "\n");
		printf("knots %d: %d %d\n", i, knot_n_xi[i], knot_n_eta[i]);
		if (knot_n_xi[i] > MAX_KNOTS) {
			printf("Error!!\n");
			printf("Knot vector too long at xi!\n"
				   "Maximum of knot vector is %d (Now %d at patch %d)\n"
				   "\n", MAX_KNOTS, knot_n_xi[i], i);
			exit(1);
		}
		if (knot_n_eta[i] > MAX_KNOTS) {
			printf("Error!!\n");
			printf("Knot vector too long at eta!\n"
				   "Maximum of knot vector is %d (Now %d at patch %d)\n"
				   "\n", MAX_KNOTS, knot_n_eta[i], i);
			exit(1);
		}
	}

	for (i = 0; i < patch_n; i++) {
		fscanf(fp, "%d%d", &cntl_p_n_xi[i], &cntl_p_n_eta[i]);
		printf("control points %d: %d %d\n",
			   i, cntl_p_n_xi[i], cntl_p_n_eta[i]);
		fscanf(fp, "\n");
	}
	printf("\n");


	for (i = 0; i < patch_n; i++) {
		for (j = 0; j < cntl_p_n_xi[i] * cntl_p_n_eta[i]; j++) {
			fscanf(fp, "%d", &temp_index[i][j]);
			printf("%d ", temp_index[i][j]);
		}
		fscanf(fp, "\n");
		printf("\n");
	}
	printf("\n");

	fscanf(fp, "%lf%lf%lf", &temp1, &temp2, &temp3);
	fscanf(fp, "\n");

	for (i = 0; i < patch_n; i++) {
		for (j = 0; j < knot_n_xi[i]; j++) {
			fscanf(fp, "%le", &knot_vec_xi[i][j]);
			printf("%f\t", knot_vec_xi[i][j]);
		}
		printf("\n");
		for (j = 0; j < knot_n_eta[i]; j++) {
			fscanf(fp, "%le", &knot_vec_eta[i][j]);
			printf("%f\t", knot_vec_eta[i][j]);
		}
		printf("\n");
	}
	printf("\n");

	for (i = 0; i < cntl_p_n; i++) {
		fscanf(fp, "%d%le%le%le",
			   &temp_int,
			   &temp_cntl_px[i], &temp_cntl_py[i], &temp_weight[i]);
		printf("%d\t%f\t%f\t%f\n",
			   temp_int,
			   temp_cntl_px[i], temp_cntl_py[i], temp_weight[i]);
	}
	printf("\n");

	for (i = 0; i < patch_n; i++) {
		for (j = 0; j < cntl_p_n_xi[i] * cntl_p_n_eta[i]; j++) {
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

	if (fields_flag) {
		printf("Start Reading displacement\n\n");
		fp = fopen(filename_disp, "r");
		char buff[256];

		fscanf(fp, "%s", buff);
		fscanf(fp, "%s", buff);

		for (i = 0; i < cntl_p_n; i++) {
			fscanf(fp, "%d:%le%le",
				   &temp_int, &temp_disp_x[i], &temp_disp_y[i]);
			printf("%d\t%1.6e\t%1.6e\n",
				   temp_int, temp_disp_x[i], temp_disp_y[i]);
		}
		printf("\n");

		for (i = 0; i < patch_n; i++) {
			for (j = 0; j < cntl_p_n_xi[i] * cntl_p_n_eta[i]; j++) {
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
}

static void Calculation (char *filename_output,
						 int order_xi, int order_eta,
						 int knot_n_xi, int knot_n_eta,
						 int cntl_p_n_xi, int cntl_p_n_eta,
						 double *knot_vec_xi, double *knot_vec_eta,
						 double *cntl_px, double *cntl_py,
						 double *disp_cntl_px, double *disp_cntl_py,
						 double *weight) {
	int i, j, k, l;
	double temp1, temp2, temp3;
	double temp_matrix[2][2];

	//計算するξ,ηの値決定と ∂ξ/∂チルダξ, ∂η/∂チルダη の計算
	double calc_xi[MAX_POINTS];		//計算するξの値
	double calc_eta[MAX_POINTS];		//計算するηの値
	double dtilda_xi[MAX_KNOTS];		// ∂ξ/∂チルダξ
	double dtilda_eta[MAX_KNOTS];	// ∂η/∂チルダη
	k = 0;
	l = 0;
	for (i = 0; i < knot_n_xi - 1; i++) {
		if ( knot_vec_xi[i] != knot_vec_xi[i + 1] ) {
			calc_xi[k] = knot_vec_xi[i];
			//printf("%d\t%f\n", k, calc_xi[k]);
			dtilda_xi[l] = ( knot_vec_xi[i + 1] - knot_vec_xi[i] ) / 2.0;
			//printf("%d\t%f\n", k, dtilda_xi[k]);
			k++;
			l++;
			if (division_ele_xi > 1) {
				temp1 = (knot_vec_xi[i + 1] - knot_vec_xi[i])
						/ (double)division_ele_xi;
				for (j = 1; j < division_ele_xi; j++) {
					calc_xi[k] = calc_xi[k - 1] + temp1;
					//printf("%d\t%f\n", k, calc_xi[k]);
					k++;
				}
			}
		}
	}
	calc_xi[k] = knot_vec_xi[knot_n_xi - 1];
	//printf("%d\t%f\n", k, calc_xi[k]);
	//printf("\n");
	division_n_xi = k + 1;
	element_n_xi = l;

	k = 0;
	l = 0;
	for (i = 0; i < knot_n_eta - 1; i++) {
		if ( knot_vec_eta[i] != knot_vec_eta[i + 1] ) {
			calc_eta[k] = knot_vec_eta[i];
			//printf("%d\t%f\n", k, calc_eta[k]);
			dtilda_eta[l] = ( knot_vec_eta[i + 1] - knot_vec_eta[i] ) / 2.0;
			//printf("%d\t%f\n", k, dtilda_eta[k]);
			k++;
			l++;
			if (division_ele_eta > 1) {
				temp1 = (knot_vec_eta[i + 1] - knot_vec_eta[i])
						/ (double)division_ele_eta;
				for (j = 1; j < division_ele_eta; j++) {
					calc_eta[k] = calc_eta[k - 1] + temp1;
					//printf("%d\t%f\n", k, calc_eta[k]);
					k++;
				}
			}
		}
	}
	calc_eta[k] = knot_vec_eta[knot_n_eta - 1];
	//printf("%d\t%f\n", k, calc_eta[k]);
	//printf("\n");
	division_n_eta = k + 1;
	element_n_eta = l;

	if (element_n_xi > MAX_ELEMENTS) {
		printf("Error!!\n");
		printf("Too many elements at xi!\n"
			   "Maximum of elements is %d (Now %d)\n"
			   "\n", MAX_ELEMENTS, element_n_xi);
		exit(1);
	}
	if (element_n_eta > MAX_ELEMENTS) {
		printf("Error!!\n");
		printf("Too many elements at eta!\n"
			   "Maximum of elements is %d (Now %d)\n"
			   "\n", MAX_ELEMENTS, element_n_eta);
		exit(1);
	}

	int ii, jj, kk, ll;

	//メッシュ座標計算
	printf("Start Calculation mesh\n\n");
	for (i = 0; i < division_n_xi; i++) {
		ii = i / division_ele_xi;
		kk = i % division_ele_xi;
		for (j = 0; j < division_n_eta; j++) {
			jj = j / division_ele_eta;
			ll = j % division_ele_eta;
			lNURBS_surface(knot_vec_xi, knot_vec_eta,
						   cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
						   weight, order_xi, order_eta,
						   calc_xi[i], calc_eta[j],
						   &coord_x[i][j], &coord_y[i][j],
						   &dxi_x[ii][jj][kk][ll], &deta_x[ii][jj][kk][ll],
						   &dxi_y[ii][jj][kk][ll], &deta_y[ii][jj][kk][ll]);
			printf("[%d][%d] [%d][%d][%d][%d]"
				   "% 1.4e % 1.4e "
				   "% 1.4e % 1.4e\n",
				   i, j, ii, jj, kk, ll,
				   calc_xi[i], calc_eta[j],
				   disp_x[i][j], disp_y[i][j]);
		}
		//printf("\n");
	}
	printf("\n");
	printf("End Calculation mesh\n\n");

	if (fields_flag) {
		//変位計算
		printf("Start Calculation displpacement\n\n");
		for (i = 0; i < division_n_xi; i++) {
			ii = i / division_ele_xi;
			kk = i % division_ele_xi;
			for (j = 0; j < division_n_eta; j++) {
				jj = j / division_ele_eta;
				ll = j % division_ele_eta;
				lNURBS_surface(knot_vec_xi, knot_vec_eta,
							   disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
							   weight, order_xi, order_eta,
							   calc_xi[i], calc_eta[j],
							   &disp_x[i][j], &disp_y[i][j],
							   &dxi_disp_x[ii][jj][kk][ll], &deta_disp_x[ii][jj][kk][ll],
							   &dxi_disp_y[ii][jj][kk][ll], &deta_disp_y[ii][jj][kk][ll]);
				printf("[%d][%d] [%d][%d][%d][%d]"
					   "% 1.4e % 1.4e "
					   "% 1.4e % 1.4e\n",
					   i, j, ii, jj, kk, ll,
					   calc_xi[i], calc_eta[j],
					   disp_x[i][j], disp_y[i][j]);
			}
			//printf("\n");
		}
		printf("\n");
		printf("End Calculation displpacement\n\n");

		//足りない微分値計算
		for (ii = 0; ii < element_n_xi; ii++) {
			for (jj = 0; jj < element_n_eta; jj++) {
				kk = division_ele_xi;
				i = (ii + 1) * division_ele_xi;
				j = jj * division_ele_eta;
				for (ll = 1; ll < division_ele_eta; ll++) {
					j++;
					rNURBS_surface(knot_vec_xi, knot_vec_eta,
								   cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
								   weight, order_xi, order_eta,
								   calc_xi[i], calc_eta[j],
								   &coord_x[i][j], &coord_y[i][j],
								   &dxi_x[ii][jj][kk][ll], &deta_x[ii][jj][kk][ll],
								   &dxi_y[ii][jj][kk][ll], &deta_y[ii][jj][kk][ll]);
					rNURBS_surface(knot_vec_xi, knot_vec_eta,
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
				for (kk = 1; kk <= division_ele_xi; kk++) {
					i++;
					rNURBS_surface(knot_vec_xi, knot_vec_eta,
								   cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
								   weight, order_xi, order_eta,
								   calc_xi[i], calc_eta[j],
								   &coord_x[i][j], &coord_y[i][j],
								   &dxi_x[ii][jj][kk][ll], &deta_x[ii][jj][kk][ll],
								   &dxi_y[ii][jj][kk][ll], &deta_y[ii][jj][kk][ll]);
					rNURBS_surface(knot_vec_xi, knot_vec_eta,
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
				rlNURBS_surface(knot_vec_xi, knot_vec_eta,
								cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
								weight, order_xi, order_eta,
								calc_xi[i], calc_eta[j],
								&coord_x[i][j], &coord_y[i][j],
								&dxi_x[ii][jj][kk][ll], &deta_x[ii][jj][kk][ll],
								&dxi_y[ii][jj][kk][ll], &deta_y[ii][jj][kk][ll]);
				rlNURBS_surface(knot_vec_xi, knot_vec_eta,
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
				lrNURBS_surface(knot_vec_xi, knot_vec_eta,
								cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
								weight, order_xi, order_eta,
								calc_xi[i], calc_eta[j],
								&coord_x[i][j], &coord_y[i][j],
								&dxi_x[ii][jj][kk][ll], &deta_x[ii][jj][kk][ll],
								&dxi_y[ii][jj][kk][ll], &deta_y[ii][jj][kk][ll]);
				lrNURBS_surface(knot_vec_xi, knot_vec_eta,
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
		for (i = 0; i < element_n_xi; i++) {
			for (j = 0; j < element_n_eta; j++) {
				temp1 = dtilda_xi[i];
				temp2 = dtilda_eta[j];
				for (k = 0; k < division_ele_xi + 1; k++) {
					for (l = 0; l < division_ele_eta + 1; l++) {
						temp_matrix[0][0] = dxi_x[i][j][k][l] * temp1;
						temp_matrix[0][1] = dxi_y[i][j][k][l] * temp1;
						temp_matrix[1][0] = deta_x[i][j][k][l] * temp2;
						temp_matrix[1][1] = deta_y[i][j][k][l] * temp2;

						InverseMatrix_2D(temp_matrix);

						strain_xx[i][j][k][l] = temp_matrix[0][0] * temp1
												* dxi_disp_x[i][j][k][l]
												+ temp_matrix[0][1] * temp2
												* deta_disp_x[i][j][k][l];
						strain_yy[i][j][k][l] = temp_matrix[1][0] * temp1
												* dxi_disp_y[i][j][k][l]
												+ temp_matrix[1][1] * temp2
												* deta_disp_y[i][j][k][l];
						strain_xy[i][j][k][l] = temp_matrix[1][0] * temp1
												* dxi_disp_x[i][j][k][l]
												+ temp_matrix[1][1] * temp2
												* deta_disp_x[i][j][k][l]
												+ temp_matrix[0][0] * temp1
												* dxi_disp_y[i][j][k][l]
												+ temp_matrix[0][1] * temp2
												* deta_disp_y[i][j][k][l];

						printf("[%d][%d][%d][%d]\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\n",
							   i, j, k, l, temp1, temp2,
							   strain_xx[i][j][k][l], strain_yy[i][j][k][l], strain_xy[i][j][k][l]);
					}
				}
			}
			printf("\n");
		}
		printf("End Calculation Strain\n\n");

		//Dマトリクスの計算
		double D_matrix[3][3] = {{0.0}};
		if (stress_type_flag == 0) { //平面応力状態
			temp1 = E * (1.0 - nu * nu);
			D_matrix[0][0] = temp1;
			D_matrix[0][1] = nu * temp1;
			D_matrix[1][0] = nu * temp1;
			D_matrix[1][1] = temp1;
			D_matrix[2][2] = (1.0 - nu) / 2.0 * temp1;
		} else if (stress_type_flag == 1) { //平面ひずみ状態(2Dの場合はこっち)
			temp1 = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
			D_matrix[0][0] = temp1;
			D_matrix[0][1] = nu / (1.0 - nu) * temp1;
			D_matrix[1][0] = nu / (1.0 - nu) * temp1;
			D_matrix[1][1] = temp1;
			D_matrix[2][2] = (1.0 - 2.0 * nu) / 2.0 / (1.0 - nu) * temp1;
		}

		printf("Start Calculation Stress\n\n");
		for (i = 0; i < element_n_xi; i++) {
			for (j = 0; j < element_n_eta; j++) {
				for (k = 0; k < division_ele_xi + 1; k++) {
					for (l = 0; l < division_ele_eta + 1; l++) {
						stress_xx[i][j][k][l] = D_matrix[0][0] * strain_xx[i][j][k][l]
												+ D_matrix[0][1] * strain_yy[i][j][k][l];
						stress_yy[i][j][k][l] = D_matrix[1][0] * strain_xx[i][j][k][l]
												+ D_matrix[1][1] * strain_yy[i][j][k][l];
						stress_xy[i][j][k][l] = D_matrix[2][2] * strain_xy[i][j][k][l];
						printf("[%d][%d][%d][%d]\t% 1.4e\t% 1.4e\t% 1.4e\n", i, j, k, l,
							   stress_xx[i][j][k][l], stress_yy[i][j][k][l], stress_xy[i][j][k][l]);
					}
				}
				printf("\n");
			}
		}
		printf("End Calculation Stress\n\n");
	}

	//書き込み
	fp = fopen(filename_output, "a");
	if (fields_flag) {
		fprintf(fp, "%d\t%d\t%d\t%d\n",
				division_n_xi, division_n_eta,
				element_n_xi, element_n_eta);
		for (i = 0; i < division_n_xi; i++) {
			for (j = 0; j < division_n_eta; j++) {
				temp1 = disp_x[i][j];
				temp2 = disp_y[i][j];
				temp3 = sqrt(temp1 * temp1 + temp2 * temp2);
				fprintf(fp, "% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\n",
						coord_x[i][j], coord_y[i][j],
						disp_x[i][j], disp_y[i][j], temp3);
			}
		}
		for (i = 0; i < element_n_xi; i++) {
			for (j = 0; j < element_n_eta; j++) {
				for (k = 0; k < division_ele_xi + 1; k++) {
					for (l = 0; l < division_ele_eta + 1; l++) {
						fprintf(fp, "% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\n",
								strain_xx[i][j][k][l], strain_yy[i][j][k][l], strain_xy[i][j][k][l],
								stress_xx[i][j][k][l], stress_yy[i][j][k][l], stress_xy[i][j][k][l]);
					}
				}
			}
		}
	} else {
		fprintf(fp, "%d\t%d\t%d\t%d\n",
				division_n_xi, division_n_eta,
				element_n_xi, element_n_eta);
		for (i = 0; i < division_n_xi; i++) {
			for (j = 0; j < division_n_eta; j++) {
				fprintf(fp, "% 1.4e\t% 1.4e\n",
						coord_x[i][j], coord_y[i][j]);
			}
		}
	}
	fclose(fp);
}

int main(int argc, char *argv[]) {
	int i;

	if (argc != 5 && argc != 6) {
		printf("Error!!\n");
		printf("Usage:\n"
			   "./NURBS_input input.txt "
			   "view.dat div_xi div_eta\n"
			   "or\n"
			   "./NURBS_input input.txt Displacement.dat "
			   "view.dat div_xi div_eta\n"
			   "\n");
		exit(1);
	}
	fp = fopen(argv[1], "r");
	if (fp == NULL) {
		printf("Error!!\n");
		printf("%s file not open!\n", argv[1]);
		exit(1);
	}
	if (argc == 6) {
		fp = fopen(argv[2], "r");
		if (fp == NULL) {
			printf("Error!!\n");
			printf("%s file not open!\n", argv[2]);
			exit(1);
		}
	}
	fclose(fp);

	division_ele_xi  = atoi(argv[argc - 2]);
	division_ele_eta = atoi(argv[argc - 1]);

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

	fields_flag = argc - 5;
	printf("fields flag:%d\n", fields_flag);

	ReadFile (argv[1], argv[2]);
	fp = fopen(argv[argc - 3], "w");
	fprintf(fp, "%d\t%d\t%d\n",
			fields_flag, division_ele_xi, division_ele_eta);
	fclose(fp);
	for (i = 0; i < patch_n; i++) {
		printf("----------Start calculation at patch %d----------\n\n", i);
		Calculation(argv[argc - 3],
					order_xi[i], order_eta[i],
					knot_n_xi[i], knot_n_eta[i],
					cntl_p_n_xi[i], cntl_p_n_eta[i],
					knot_vec_xi[i], knot_vec_eta[i],
					cntl_px[i], cntl_py[i],
					disp_cntl_px[i], disp_cntl_py[i],
					weight[i]);
		printf("-----------End calculation at patch %d-----------\n\n", i);
	}
	return 0;
}

/*
	fp = fopen("Virtual_Crack_Extension.dat", "w");
	fprintf(fp, "label=Virtual_Crack_Extension\n");
	fprintf(fp, "num_items=%d\n", Total_Control_Point);
	fprintf( fp, "\n");
	for (j = 0; j < Total_Control_Point; j++ ) {
		fprintf(fp, "%d:	%le %le ", j,
				Virtual_Crack_Extension_Ct_Pt[j][0],
				Virtual_Crack_Extension_Ct_Pt[j][1]);
		fprintf(fp, "\n");
	}
	fclose(fp);
	*/