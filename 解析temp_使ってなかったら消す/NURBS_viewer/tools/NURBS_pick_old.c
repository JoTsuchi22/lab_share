#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#define MAX_ELEMENTS 10000	//最大要素数
#define MAX_DIVISION 20		//一要素一方向あたりの最大分割数
#define MAX_POINTS 2000		//一方向あたりの最大点数
#define MAX_CNRL_P 5000		//最大コントロールポイント数(n)
#define MAX_KNOTS 300			//ノットベクトルの最大長さ(n+p+1)
#define MAX_ORDER 3			//最大次数(p)

int order_xi;					//ξ基底関数の次数(p)
int order_eta;					//η基底関数の次数(p)
int knot_n_xi;					//ξノットベクトルの数(n+p+1)
int knot_n_eta;				//ηノットベクトルの数(n+p+1)
int cntl_p_n_xi;				//ξ方向コントロールポイント数(n)
int cntl_p_n_eta;				//η方向コントロールポイント数(n)

double knot_vec_xi[MAX_KNOTS];		//ξノットベクトル
double knot_vec_eta[MAX_KNOTS];		//ηノットベクトル
double cntl_px[MAX_CNRL_P];			//コントロールポイントx座標
double cntl_py[MAX_CNRL_P];			//コントロールポイントy座標
double disp_cntl_px[MAX_CNRL_P];		//コントロールポイント上のx方向変位
double disp_cntl_py[MAX_CNRL_P];		//コントロールポイント上のy方向変位
double weight[MAX_CNRL_P];			//重み

double coord_x[MAX_POINTS][MAX_POINTS];			//変形前x座標
double coord_y[MAX_POINTS][MAX_POINTS];			//変形前y座標

double disp_x[MAX_POINTS][MAX_POINTS];			//x方向変位
double disp_y[MAX_POINTS][MAX_POINTS];			//y方向変位
double disp_val[MAX_POINTS][MAX_POINTS];		//変位の絶対値

double strain_xx[MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//x方向ひずみ
double strain_yy[MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//y方向ひずみ
double strain_xy[MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//剪断ひずみ

double stress_xx[MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//x方向応力
double stress_yy[MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//y方向応力
double stress_xy[MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//剪断応力

double strain_rr[MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//r方向ひずみ
double strain_tt[MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//θ方向ひずみ
double strain_rt[MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//極座標剪断ひずみ

double stress_rr[MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//r方向応力
double stress_tt[MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//θ方向応力
double stress_rt[MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//極座標剪断応力

int division_ele_xi;		//ξ方向の一要素あたりの分割数
int division_ele_eta;		//η方向の一要素あたりの分割数
int division_n_xi;		//ξ方向の表示する点の数
int division_n_eta;		//η方向の表示する点の数
int element_n_xi;			//ξ方向要素数
int element_n_eta;		//η方向要素数

FILE *fp; // FILE型構造体

/* 読み込み */
static void ReadFile (void) {
	char fname[] = "input.txt";
	int i, j, k, temp_int;
	double temp1, temp2, temp3;

	fp = fopen(fname, "r"); // ファイルを開く。失敗するとNULLを返す。
	if (fp == NULL) {
		printf("%s file not open!\n", fname);
	}

	printf("Start Reading input\n");

	fscanf(fp, "%d", &temp_int);
	fscanf(fp, "\n");

	fscanf(fp, "%d%d", &order_xi, &order_eta);
	fscanf(fp, "\n");
	printf("order:%d %d\n", order_xi, order_eta);

	fscanf(fp, "%d%d", &knot_n_xi, &knot_n_eta);
	fscanf(fp, "\n");
	printf("knots:%d %d\n", knot_n_xi, knot_n_eta);

	fscanf(fp, "%d%d", &cntl_p_n_xi, &cntl_p_n_eta);
	printf("control points:%d %d\n", cntl_p_n_xi, cntl_p_n_eta);
	fscanf(fp, "\n");

	fscanf(fp, "%lf%lf%lf", &temp1, &temp2, &temp3);
	fscanf(fp, "\n");

	for (i = 0; i < knot_n_xi; i++) {
		fscanf(fp, "%le", &knot_vec_xi[i]);
		printf("% 1.8e\t", knot_vec_xi[i]);
	}
	printf("\n");

	for (i = 0; i < knot_n_eta; i++) {
		fscanf(fp, "%le", &knot_vec_eta[i]);
		printf("% 1.8e\t", knot_vec_eta[i]);
	}
	printf("\n\n");

	fscanf(fp, "\n");

	for (i = 0; i < cntl_p_n_xi * cntl_p_n_eta; i++) {
		fscanf(fp, "%d%le%le%le",
		       &temp_int, &cntl_px[i], &cntl_py[i], &weight[i]);
		printf("%d\t% 1.8e\t% 1.8e\t% 1.8e\n",
		       temp_int, cntl_px[i], cntl_py[i], weight[i]);
	}
	fclose(fp); // ファイルを閉じる
	printf("End Reading input\n\n");

	fp = fopen("points.dat", "r");
	if (fp == NULL) {
		printf("%s points file not open!\n", "points.dat");
	}

	printf("Start Reading points\n");

	fscanf(fp, "%d%d%d%d%d%d",
	       &division_ele_xi, &division_ele_eta,
	       &division_n_xi, &division_n_eta,
	       &element_n_xi, &element_n_eta);
	printf("%d\t%d\t%d\t%d\t%d\t%d\n",
	       division_ele_xi, division_ele_eta,
	       division_n_xi, division_n_eta,
	       element_n_xi, element_n_eta);

	for (i = 0; i < division_n_xi; i++) {
		for (j = 0; j < division_n_eta; j++) {
			fscanf(fp, "%le%le%le%le%le",
			       &coord_x[i][j], &coord_y[i][j],
			       &disp_x[i][j], &disp_y[i][j], &disp_val[i][j]);
			printf("%d\t%d\t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\n",
			       i, j, coord_x[i][j], coord_y[i][j],
			       disp_x[i][j], disp_y[i][j], disp_val[i][j]);
		}
	}

	for (i = 0; i < element_n_xi * element_n_eta; i++) {
		for (j = 0; j < division_ele_xi + 1; j++) {
			for (k = 0; k < division_ele_eta + 1; k++) {
				fscanf(fp, "%le%le%le%le%le%le\n",
				       &strain_xx[i][j][k], &strain_yy[i][j][k],
				       &strain_xy[i][j][k],
				       &stress_xx[i][j][k], &stress_yy[i][j][k],
				       &stress_xy[i][j][k]);
				printf("%d\t%d\t%d\t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\n",
				       i, j, k,
				       strain_xx[i][j][k], strain_yy[i][j][k],
				       strain_xy[i][j][k],
				       stress_xx[i][j][k], stress_yy[i][j][k],
				       stress_xy[i][j][k]);
			}
		}
	}
	printf("End Reading points\n\n");
	fclose(fp);
}

static void PickDisp (int size1, double value[size1][size1]) {
	int i, j;
	double x, y, scalar, xi, eta;
	double temp1, temp2;
	temp1 = (knot_vec_xi[knot_n_xi - 1] - knot_vec_xi[0])
	        / (double)(division_n_xi - 1);
	temp2 = (knot_vec_eta[knot_n_eta - 1] - knot_vec_eta[0])
	        / (double)(division_n_eta - 1);

	fp = fopen("pick_disp.dat", "w");
	for (i = 0; i < division_n_xi; i++) {
		//for (i = 0; i < division_n_xi; i++) {
		for (j = 0; j < 1; j++) {
			//for (j = 0; j < division_n_eta; j++) {
			xi = temp1 * (double)i + knot_vec_xi[0];
			eta = temp2 * (double)j + knot_vec_eta[0];
			x = coord_x[i][j] + disp_x[i][j];
			y = coord_y[i][j] + disp_y[i][j];
			scalar = value[i][j];
			printf("% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\n",
			       xi, eta, x, y, scalar);
			fprintf(fp, "% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\n",
			        xi, eta, x, y, scalar);
		}
	}
	fclose(fp);
}

static void PickStrainStress (int coord, int num, int size1, int size2,
                              double value1[size1][size2 + 1][size2 + 1],
                              double value2[size1][size2 + 1][size2 + 1],
                              double value3[size1][size2 + 1][size2 + 1]) {
	//coord 0:ξ固定 1:η固定
	int i, j = 0, k, l, temp_index_xi = 0, temp_index_eta = 0;
	double x, y, scalar1, scalar2, scalar3;

	fp = fopen("pick_disp.dat", "w");
	if (coord == 0) {
		j = num % division_ele_xi;
		i = num / division_ele_xi;
		for (l = 0; l < element_n_eta; l++) {
			temp_index_xi = (i % element_n_xi) * division_ele_xi;
			temp_index_eta = (i / element_n_xi) * division_ele_eta;
			for (k = 0; k < division_ele_eta; k++) {
				x = coord_x[temp_index_xi + j][temp_index_eta + k];
				y = coord_y[temp_index_xi + j][temp_index_eta + k];
				scalar1 = value1[i][j][k];
				scalar2 = value2[i][j][k];
				scalar3 = value3[i][j][k];
				printf("t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\n",
				       x, y, scalar1, scalar2, scalar3);
				fprintf(fp, "% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\n",
				        x, y, scalar1, scalar2, scalar3);
			}
			i += element_n_xi;
		}
		x = coord_x[temp_index_xi + j][temp_index_eta + division_ele_eta];
		y = coord_y[temp_index_xi + j][temp_index_eta + division_ele_eta];
		scalar1 = value1[i - 1][j][k];
		scalar2 = value2[i - 1][j][k];
		scalar3 = value3[i - 1][j][k];
		printf("% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\n",
		       x, y, scalar1, scalar2, scalar3);
		fprintf(fp, "% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\n",
		        x, y, scalar1, scalar2, scalar3);
		printf("\n");
	} else if (coord == 1) {
		if (num != element_n_eta * division_ele_eta) {
			k = num / division_ele_eta;
			i = k * element_n_xi;
			for (l = 0; l < element_n_xi; l++) {
				temp_index_xi = (i % element_n_xi) * division_ele_xi;
				temp_index_eta = (i / element_n_xi) * division_ele_eta;
				for (j = 0; j < division_ele_eta; j++) {
					x = coord_x[temp_index_xi + j][temp_index_eta + k];
					y = coord_y[temp_index_xi + j][temp_index_eta + k];
					scalar1 = value1[i][j][k];
					scalar2 = value2[i][j][k];
					scalar3 = value3[i][j][k];
					printf("% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\n",
					       x, y, scalar1, scalar2, scalar3);
					fprintf(fp, "% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\n",
					        x, y, scalar1, scalar2, scalar3);
				}
				i += 1;
			}
			x = coord_x[temp_index_xi + division_ele_eta][temp_index_eta + k];
			y = coord_y[temp_index_xi + division_ele_eta][temp_index_eta + k];
			scalar1 = value1[i - 1][j][k];
			scalar2 = value2[i - 1][j][k];
			scalar3 = value3[i - 1][j][k];
			printf("% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\n",
			       x, y, scalar1, scalar2, scalar3);
			fprintf(fp, "% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\n",
			        x, y, scalar1, scalar2, scalar3);
			printf("\n");
		} else {
			k = division_ele_eta;
			i = element_n_xi * (element_n_eta - 1);
			for (l = 0; l < element_n_xi; l++) {
				temp_index_xi = (i % element_n_xi) * division_ele_xi;
				temp_index_eta = (i / element_n_xi) * division_ele_eta;
				for (j = 0; j < division_ele_eta; j++) {
					x = coord_x[temp_index_xi + j][temp_index_eta + k];
					y = coord_y[temp_index_xi + j][temp_index_eta + k];
					scalar1 = value1[i][j][k];
					scalar2 = value2[i][j][k];
					scalar3 = value3[i][j][k];
					printf("% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\n",
					       x, y, scalar1, scalar2, scalar3);
					fprintf(fp, "% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\n",
					        x, y, scalar1, scalar2, scalar3);
				}
				i += 1;
			}
			x = coord_x[temp_index_xi + division_ele_eta][temp_index_eta + k];
			y = coord_y[temp_index_xi + division_ele_eta][temp_index_eta + k];
			scalar1 = value1[i - 1][j][k];
			scalar2 = value2[i - 1][j][k];
			scalar3 = value3[i - 1][j][k];
			printf("% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\n",
			       x, y, scalar1, scalar2, scalar3);
			fprintf(fp, "% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\n",
			        x, y, scalar1, scalar2, scalar3);
			printf("\n");
		}
	} else {
		printf("error!!!\n");
	}
	fclose(fp);
}


static void CalcRadiusStrainStress (void) {
	int i, j, k, temp_index_xi, temp_index_eta;
	double temp_x, temp_y;
	double temp1, temp2, temp3;
	double temp_cos, temp_sin;

	for (i = 0; i < element_n_xi * element_n_eta; i++) {
		temp_index_xi = (i % element_n_xi) * division_ele_xi;
		temp_index_eta = (i / element_n_xi) * division_ele_eta;
		for (j = 0; j <= division_ele_xi; j++) {
			for (k = 0; k <= division_ele_eta; k++) {
				temp_x = coord_x[temp_index_xi + j][temp_index_eta + k];
				temp_y = coord_y[temp_index_xi + j][temp_index_eta + k];
				temp1 = sqrt(temp_x * temp_x + temp_y * temp_y);
				temp_cos = temp_x / temp1;
				temp_sin = temp_y / temp1;
				temp1 = temp_cos * temp_cos;
				temp2 = temp_sin * temp_sin;
				temp3 = temp_cos * temp_sin;
				strain_rr[i][j][k] = strain_xx[i][j][k] * temp1
				                     + strain_yy[i][j][k] * temp2
				                     + 2.0 * strain_xy[i][j][k] * temp3;
				strain_tt[i][j][k] = strain_xx[i][j][k] * temp2
				                     + strain_yy[i][j][k] * temp1
				                     - 2.0 * strain_xy[i][j][k] * temp3;
				strain_rt[i][j][k] = - strain_xx[i][j][k] * temp3
				                     + strain_yy[i][j][k] * temp3
				                     + strain_xy[i][j][k] * (temp1 - temp2);
				printf("% 1.8e\t% 1.8e\t% 1.8e\t",
				       strain_rr[i][j][k], strain_tt[i][j][k], strain_rt[i][j][k]);

				stress_rr[i][j][k] = stress_xx[i][j][k] * temp1
				                     + stress_yy[i][j][k] * temp2
				                     + 2.0 * stress_xy[i][j][k] * temp3;
				stress_tt[i][j][k] = stress_xx[i][j][k] * temp2
				                     + stress_yy[i][j][k] * temp1
				                     - 2.0 * stress_xy[i][j][k] * temp3;
				stress_rt[i][j][k] = - stress_xx[i][j][k] * temp3
				                     + stress_yy[i][j][k] * temp3
				                     + stress_xy[i][j][k] * (temp1 - temp2);
				printf("% 1.8e\t% 1.8e\t% 1.8e\n",
				       stress_rr[i][j][k], stress_tt[i][j][k], stress_rt[i][j][k]);
			}
		}
	}
	printf("\n");

	fp = fopen("points.dat", "w");
	fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\n",
	        division_ele_xi, division_ele_eta,
	        division_n_xi, division_n_eta,
	        element_n_xi, element_n_eta);
	for (i = 0; i < division_n_xi; i++) {
		for (j = 0; j < division_n_eta; j++) {
			fprintf(fp, "%f\t%f\t% 1.8e\t% 1.8e\t% 1.8e\n",
			        coord_x[i][j], coord_y[i][j],
			        disp_x[i][j], disp_y[i][j], disp_val[i][j]);
		}
	}
	for (i = 0; i < element_n_xi * element_n_eta; i++) {
		for (j = 0; j <= division_ele_xi; j++) {
			for (k = 0; k <= division_ele_eta; k++) {
				fprintf(fp, "% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\n",
				        strain_rr[i][j][k], strain_tt[i][j][k], strain_rt[i][j][k],
				        stress_rr[i][j][k], stress_tt[i][j][k], stress_rt[i][j][k]);
			}
		}
	}
	fclose(fp);
}

int main(int argc, char *argv[]) {
	ReadFile();

	if (argc != 2) {
		exit(1);
	}

	int i = atoi(argv[2]);
	if (i < 0 || i > 3) {
		exit(1);
	}

	switch (i) {
	case 0:
		CalcRadiusStrainStress();
		break;
	case 1:
		PickDisp (MAX_POINTS, disp_val);
		break;
	case 2:
		PickStrainStress (1, 0, MAX_ELEMENTS, MAX_DIVISION,
		                  stress_xx, stress_yy, stress_xy);
		break;
	}
	return 0;
}