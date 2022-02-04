#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "../NURBS_view/constant.h"

double E;						//ヤング率(GPa)
double nu;						//ポアソン比(-)

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

double coord_x[MAX_PATCHES][MAX_POINTS][MAX_POINTS];		//変形前x座標
double coord_y[MAX_PATCHES][MAX_POINTS][MAX_POINTS];		//変形前y座標

double disp_x[MAX_PATCHES][MAX_POINTS][MAX_POINTS];		//x方向変位
double disp_y[MAX_PATCHES][MAX_POINTS][MAX_POINTS];		//y方向変位
double disp_val[MAX_PATCHES][MAX_POINTS][MAX_POINTS];		//変位の絶対値

double strain_xx[MAX_PATCHES][MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//x方向ひずみ
double strain_yy[MAX_PATCHES][MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//y方向ひずみ
double strain_xy[MAX_PATCHES][MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//剪断ひずみ

double stress_xx[MAX_PATCHES][MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//x方向応力
double stress_yy[MAX_PATCHES][MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//y方向応力
double stress_xy[MAX_PATCHES][MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//剪断応力

double strain_rr[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//x方向ひずみ
double strain_tt[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//y方向ひずみ
double strain_rt[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//剪断ひずみ

double stress_rr[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//x方向応力
double stress_tt[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//y方向応力
double stress_rt[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//剪断応力

int fields_flag;
int division_ele_xi;					//ξ方向の一要素あたりの分割数
int division_ele_eta;					//η方向の一要素あたりの分割数
int division_n_xi[MAX_PATCHES];		//ξ方向の表示する点の数
int division_n_eta[MAX_PATCHES];		//η方向の表示する点の数
int element_n_xi[MAX_PATCHES];		//ξ方向要素数
int element_n_eta[MAX_PATCHES];		//η方向要素数

int temp_index[MAX_PATCHES][MAX_CNRL_P];
double temp_cntl_px[MAX_CNRL_P];
double temp_cntl_py[MAX_CNRL_P];
double temp_weight[MAX_CNRL_P];
double temp_disp_x[MAX_CNRL_P];
double temp_disp_y[MAX_CNRL_P];

FILE *fp;	// FILE型構造体

void ReadFile (char *filename_input, char *filename_view, char *filename_disp) {
	int i, j, k, l, m, temp_int;
	double temp1, temp2, temp3;

	fp = fopen(filename_input, "r");
	printf("Start Reading input\n");

	fscanf(fp, "%lf%lf", &E, &nu);
	fscanf(fp, "\n");
	printf("E nu:% 1.8e % 1.8e\n", E, nu);

	fscanf(fp, "%d", &patch_n);
	fscanf(fp, "\n");
	printf("patches:%d \n", patch_n);

	fscanf(fp, "%d", &cntl_p_n);
	fscanf(fp, "\n");
	printf("control points:%d \n", cntl_p_n);

	for (i = 0; i < patch_n; i++) {
		fscanf(fp, "%d%d", &order_xi[i], &order_eta[i]);
		fscanf(fp, "\n");
		printf("order%d:%d %d\n", i, order_xi[i], order_eta[i]);
	}


	for (i = 0; i < patch_n; i++) {
		fscanf(fp, "%d%d", &knot_n_xi[i], &knot_n_eta[i]);
		fscanf(fp, "\n");
		printf("knots%d:%d %d\n", i, knot_n_xi[i], knot_n_eta[i]);
	}

	for (i = 0; i < patch_n; i++) {
		fscanf(fp, "%d%d", &cntl_p_n_xi[i], &cntl_p_n_eta[i]);
		printf("control points%d (xi, eta):%d %d\n",
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


	printf("Start Reading points\n");
	fp = fopen(filename_view, "r");
	fscanf(fp, "%d\t%d\t%d\n",
	       &fields_flag, &division_ele_xi, &division_ele_eta);
	if (fields_flag) {
		for (i = 0; i < patch_n; i++) {
			fscanf(fp, "%d%d%d%d",
			       &division_n_xi[i], &division_n_eta[i],
			       &element_n_xi[i], &element_n_eta[i]);
			printf("%d\t%d\t%d\t%d\n",
			       division_n_xi[i], division_n_eta[i],
			       element_n_xi[i], element_n_eta[i]);
			for (j = 0; j < division_n_xi[i]; j++) {
				for (k = 0; k < division_n_eta[i]; k++) {
					fscanf(fp, "%le%le%le%le%le",
					       &coord_x[i][j][k], &coord_y[i][j][k],
					       &disp_x[i][j][k], &disp_y[i][j][k], &disp_val[i][j][k]);
					printf("%d\t%d\t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\n",
					       i, j, coord_x[i][j][k], coord_y[i][j][k],
					       disp_x[i][j][k], disp_y[i][j][k], disp_val[i][j][k]);
				}
			}

			for (j = 0; j < element_n_xi[i]; j++) {
				for (k = 0; k < element_n_eta[i]; k++) {
					for (l = 0; l < division_ele_xi + 1; l++) {
						for (m = 0; m < division_ele_eta + 1; m++) {
							fscanf(fp, "%le%le%le%le%le%le\n",
							       &strain_xx[i][j][k][l][m], &strain_yy[i][j][k][l][m],
							       &strain_xy[i][j][k][l][m],
							       &stress_xx[i][j][k][l][m], &stress_yy[i][j][k][l][m],
							       &stress_xy[i][j][k][l][m]);
							printf("[%d][%d][%d][%d][%d]\t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\n",
							       i, j, k, l, m,
							       strain_xx[i][j][k][l][m], strain_yy[i][j][k][l][m],
							       strain_xy[i][j][k][l][m],
							       stress_xx[i][j][k][l][m], stress_yy[i][j][k][l][m],
							       stress_xy[i][j][k][l][m]);
						}
					}
				}
			}
		}
	} else {
		for (i = 0; i < patch_n; i++) {
			fscanf(fp, "%d%d%d%d",
			       &division_n_xi[i], &division_n_eta[i],
			       &element_n_xi[i], &element_n_eta[i]);
			printf("%d\t%d\t%d\t%d\n",
			       division_n_xi[i], division_n_eta[i],
			       element_n_xi[i], element_n_eta[i]);
			for (j = 0; j < division_n_xi[i]; j++) {
				for (k = 0; k < division_n_eta[i]; k++) {
					fscanf(fp, "%le%le",
					       &coord_x[i][j][k], &coord_y[i][j][k]);
					printf("%d\t%d\t% 1.8e\t% 1.8e\n",
					       i, j, coord_x[i][j][k], coord_y[i][j][k]);
				}
			}
		}
	}
	fclose(fp);
	printf("End Reading points\n\n");

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

void CalcRadiusStrainStress(char *filename_output,
                            int size1, int size2, int size3,
                            int division_n_xi, int division_n_eta,
                            int element_n_xi, int element_n_eta,
                            double coord_x[size1][size1],
                            double coord_y[size1][size1],
                            double disp_x[size1][size1],
                            double disp_y[size1][size1],
                            double disp_val[size1][size1],
                            double strain_xx[size2][size2][size3 + 1][size3 + 1],
                            double strain_yy[size2][size2][size3 + 1][size3 + 1],
                            double strain_xy[size2][size2][size3 + 1][size3 + 1],
                            double stress_xx[size2][size2][size3 + 1][size3 + 1],
                            double stress_yy[size2][size2][size3 + 1][size3 + 1],
                            double stress_xy[size2][size2][size3 + 1][size3 + 1]) {
	int i, j, k, l;
	double temp_x, temp_y;
	double temp1;
	for (i = 0; i < element_n_xi; i++) {
		for (j = 0; j < element_n_eta; j++) {
			for (k = 0; k < division_ele_xi + 1; k++) {
				for (l = 0; l < division_ele_eta + 1; l++) {
					temp_x = coord_x[i * division_ele_xi + k][j * division_ele_eta + l];
					temp_y = coord_y[i * division_ele_xi + k][j * division_ele_eta + l];
					temp1 = temp_x * temp_x + temp_y * temp_y;

					stress_rr[i][j][k][l] = -4.0 / 3.0 / temp1 + 1.0 / 3.0;
					stress_tt[i][j][k][l] =  4.0 / 3.0 / temp1 + 1.0 / 3.0;
					stress_rt[i][j][k][l] = 0.0;

					strain_rr[i][j][k][l] = 0.91 / 206000.0 * (stress_rr[i][j][k][l]
					                        - 3.0 / 7.0 * stress_tt[i][j][k][l]);
					strain_tt[i][j][k][l] = 0.91 / 206000.0 * (stress_tt[i][j][k][l]
					                        - 3.0 / 7.0 * stress_rr[i][j][k][l]);
					strain_rt[i][j][k][l] = 0.0;

					disp_val[i * division_ele_xi + k][j * division_ele_eta + l]
					    = strain_tt[i][j][k][l] * sqrt(temp1);

					printf("% 1.8e\t% 1.8e\t% 1.8e\t",
					       strain_rr[i][j][k][l], strain_tt[i][j][k][l], strain_rt[i][j][k][l]);
					printf("% 1.8e\t% 1.8e\t% 1.8e\n",
					       stress_rr[i][j][k][l], stress_tt[i][j][k][l], stress_rt[i][j][k][l]);
				}
			}
		}
	}

	fp = fopen(filename_output, "a");
	fprintf(fp, "%d\t%d\t%d\t%d\n",
	        division_n_xi, division_n_eta,
	        element_n_xi, element_n_eta);
	for (i = 0; i < division_n_xi; i++) {
		for (j = 0; j < division_n_eta; j++) {
			temp1 = disp_val[i][j];
			fprintf(fp, "% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\n",
			        coord_x[i][j], coord_y[i][j],
			        disp_x[i][j], disp_y[i][j], temp1);
		}
	}
	for (i = 0; i < element_n_xi; i++) {
		for (j = 0; j < element_n_eta; j++) {
			for (k = 0; k < division_ele_xi + 1; k++) {
				for (l = 0; l < division_ele_eta + 1; l++) {
					fprintf(fp, "% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\n",
					        strain_rr[i][j][k][l], strain_tt[i][j][k][l], strain_rt[i][j][k][l],
					        stress_rr[i][j][k][l], stress_tt[i][j][k][l], stress_rt[i][j][k][l]);
				}
			}
		}
	}
	fclose(fp);
}

int main(int argc, char *argv[]) {
	if (argc != 5) {
		exit(1);
	}

	ReadFile(argv[1], argv[2], argv[3]);

	fp = fopen(argv[4], "w");
	fprintf(fp, "%d\t%d\t%d\n",
	        fields_flag, division_ele_xi, division_ele_eta);
	fclose(fp);

	int i;
	for (i = 0; i < patch_n; i++) {
		printf("----------Start calculation at patch %d----------\n\n", i);
		CalcRadiusStrainStress(argv[4],
		                       MAX_POINTS, MAX_ELEMENTS, MAX_DIVISION,
		                       division_n_xi[i], division_n_eta[i],
		                       element_n_xi[i], element_n_eta[i],
		                       coord_x[i], coord_y[i],
		                       disp_x[i], disp_y[i],  disp_val[i],
		                       strain_xx[i], strain_yy[i], strain_xy[i],
		                       stress_xx[i], stress_yy[i], stress_xy[i]);
		printf("-----------End calculation at patch %d-----------\n\n", i);
	}
	return 0;
}