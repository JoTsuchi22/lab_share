#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#include "../NURBS_view/constant.h"

#define DIMENSION 2
#define MAX_N_NODE 110000

int stress_type_flag = 0;			//0:平面応力 1:平面ひずみ
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

double coord_x[MAX_POINTS][MAX_POINTS];		//変形前x座標
double coord_y[MAX_POINTS][MAX_POINTS];		//変形前y座標
double dxi_x[MAX_POINTS][MAX_POINTS];		// ∂x/∂ξ
double dxi_y[MAX_POINTS][MAX_POINTS];		// ∂y/∂ξ
double deta_x[MAX_POINTS][MAX_POINTS];		// ∂x/∂η
double deta_y[MAX_POINTS][MAX_POINTS];		// ∂y/∂η

double disp_x[MAX_POINTS][MAX_POINTS];			//x方向変位
double disp_y[MAX_POINTS][MAX_POINTS];			//y方向変位
double dxi_disp_x[MAX_POINTS][MAX_POINTS];		// ∂u/∂ξ
double dxi_disp_y[MAX_POINTS][MAX_POINTS];		// ∂v/∂ξ
double deta_disp_x[MAX_POINTS][MAX_POINTS];		// ∂u/∂η
double deta_disp_y[MAX_POINTS][MAX_POINTS];		// ∂v/∂η

double strain_xx[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//x方向ひずみ
double strain_yy[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//y方向ひずみ
double strain_xy[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//剪断ひずみ

double stress_xx[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//x方向応力
double stress_yy[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//y方向応力
double stress_xy[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//剪断応力

int fields_flag;
int division_ele_xi;		//ξ方向の一要素あたりの分割数
int division_ele_eta;		//η方向の一要素あたりの分割数
int division_n_xi;		//ξ方向の表示する点の数
int division_n_eta;		//η方向の表示する点の数
int element_n_xi;			//ξ方向要素数
int element_n_eta;		//η方向要素数

FILE *fp;	// FILE型構造体

void ReadFile (char *filename_input) {
	int i, j;
	double temp1, temp2, temp3;
	int temp_int;

	printf("Start Reading input\n\n");
	fp = fopen(filename_input, "r");

	fscanf(fp, "%lf%lf", &E, &nu);
	fscanf(fp, "\n");
	printf("E nu: % 1.8e % 1.8e\n", E, nu);

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

	int temp_index[MAX_PATCHES][MAX_CNRL_P];
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
		for (j = 0; j < knot_n_eta[i]; j++) {
			fscanf(fp, "%le", &knot_vec_eta[i][j]);
			printf("%f\t", knot_vec_eta[i][j]);
		}
		printf("\n");
	}
	printf("\n");


	double temp_cntl_px[MAX_CNRL_P * MAX_PATCHES];
	double temp_cntl_py[MAX_CNRL_P * MAX_PATCHES];
	double temp_weight[MAX_CNRL_P * MAX_PATCHES];
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
}



int main(int argc, char *argv[]) {
	int Location_Crack_Tip_Patch;
	int ii, jj, kk, j;
	int Num_Non_Zero_Pt;
	int kk_control_pt;
	double Location_Local_Coordinates[DIMENSION];
	double Virtual_Crack_Extension_Ct_Pt[MAX_N_NODE][DIMENSION];
	double DeltaA;

	division_ele_xi  = atoi(argv[2]);
	division_ele_eta = atoi(argv[3]);
	ReadFile(argv[1]);

	for (jj = 0; jj < DIMENSION; jj++) {
		for (ii = 0; ii < cntl_p_n; ii++) {
			Virtual_Crack_Extension_Ct_Pt[ii][jj] = 0.0;
		}
	}

	fp = fopen("J-int.inp", "r");
	fscanf(fp, "%d %lf %lf %lf",
	       &Location_Crack_Tip_Patch,
	       &Location_Local_Coordinates[0],
	       &Location_Local_Coordinates[1], &DeltaA);
	fscanf(fp, "%d", &Num_Non_Zero_Pt);
	printf("%d\n", Num_Non_Zero_Pt);
	for (kk = 0; kk < Num_Non_Zero_Pt; kk++) {
		fscanf(fp, "%d", &kk_control_pt);
		for (jj = 0; jj < DIMENSION; jj++)
			fscanf(fp, "%lf", &Virtual_Crack_Extension_Ct_Pt[kk_control_pt][jj]);
		printf("kk = %d  kk_control_pt = %d  coord values %f %f\n",
		       kk, kk_control_pt,
		       Virtual_Crack_Extension_Ct_Pt[kk_control_pt][0],
		       Virtual_Crack_Extension_Ct_Pt[kk_control_pt][1]);
	}
	fclose(fp);

	fp = fopen("Virtual_Crack_Extension.dat", "w");
	fprintf(fp, "label=Virtual_Crack_Extension\n");
	fprintf(fp, "num_items=%d\n", cntl_p_n);
	fprintf( fp, "\n");
	for (j = 0; j < cntl_p_n; j++ ) {
		fprintf(fp, "%d:	%le %le ",
		        j,
		        Virtual_Crack_Extension_Ct_Pt[j][0],
		        Virtual_Crack_Extension_Ct_Pt[j][1]);
		fprintf(fp, "\n");
	}
	fclose(fp);
	return 0;
}