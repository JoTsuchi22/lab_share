#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#include "NURBS_calc.h"

/////変更後はmake cleanしてからmakeする/////////////////////////////////////////////////////
///////最大値////////////////////////////////////////////////////////////////////////////////
#define MAX_ORDER 2										//最大次数(p)
#define MAX_CNRL_P 5000									//最大コントロールポイント数(n)
#define MAX_KNOTS (MAX_CNRL_P + MAX_ORDER + 1)			//ノットベクトルの最大長さ(n+p+1)
///////各方向での最大値//////////////////////////////////////////////////////////////////////
#define MAX_ELEMENTS 100									//最大要素数
#define MAX_DIVISION 100									//一要素あたりの最大分割数
#define MAX_POINTS (MAX_ELEMENTS	 * MAX_DIVISION + 1)	//最大点数
/////////////////////////////////////////////////////////////////////////////////////////////

double E;			//ヤング率(GPa)
double nu;			//ポアソン比(-)
int order;			//ξ基底関数の次数(p)
int knot_n;		//ξノットベクトルの数(n+p+1)
int cntl_p_n;		//ξ方向コントロールポイント数(n)

int division_ele;	//ξ方向の一要素あたりの分割数
int point_n;		//ξ方向の表示する点の数
int element_n;		//ξ方向要素数

static double knot_vec[MAX_KNOTS];		//ξノットベクトル
static double cntl_px[MAX_CNRL_P];		//コントロールポイントx座標
static double cntl_py[MAX_CNRL_P];		//コントロールポイントy座標
static double weight[MAX_CNRL_P];		//重み


static double output[MAX_CNRL_P][MAX_POINTS];
static double doutput[MAX_CNRL_P][MAX_POINTS];
static double output_x[MAX_POINTS];
static double output_y[MAX_POINTS];
static double calc_xi[MAX_POINTS];		//計算するξの値

FILE *fp;	// FILE型構造体

void ReadFile (char *filename_input) {
	int j;

	printf("Start Reading input\n\n");
	fp = fopen(filename_input, "r");

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

	fscanf(fp, "%d", &order);
	fscanf(fp, "\n");
	printf("order: %d\n", order);
	if (order > MAX_ORDER) {
		printf("Error!!\n");
		printf("Order too big at xi!\n"
		       "Maximum of order is %d (Now %d)\n"
		       "\n", MAX_ORDER, order);
		exit(1);
	}

	fscanf(fp, "%d", &knot_n);
	fscanf(fp, "\n");
	printf("knots: %d\n", knot_n);
	if (knot_n > MAX_KNOTS) {
		printf("Error!!\n");
		printf("Knot vector too long at xi!\n"
		       "Maximum of knot vector is %d (Now %d)\n"
		       "\n", MAX_KNOTS, knot_n);
		exit(1);
	}

	for (j = 0; j < knot_n; j++) {
		fscanf(fp, "%le", &knot_vec[j]);
		printf("%f\t", knot_vec[j]);
	}
	printf("\n\n");

	for (j = 0; j < cntl_p_n; j++) {
		fscanf(fp, "%d%le%le%le",
		       &j, &cntl_px[j], &cntl_py[j], &weight[j]);
		printf("%d\t%f\t%f\t%f\n",
		       j, cntl_px[j], cntl_py[j], weight[j]);
	}
	printf("\n");
	fclose(fp);
	printf("End Reading input\n\n");
}

static void Calculation (int order, int knot_n,
                         int cntl_p_n, double *knot_vec,
                         double *cntl_px, double *cntl_py,
                         double *weight) {
	int i, j, k, l;
	double temp1;

	//計算するξ,ηの値決定
	k = 0;
	l = 0;
	for (i = 0; i < knot_n - 1; i++) {
		if ( knot_vec[i] != knot_vec[i + 1] ) {
			calc_xi[k] = knot_vec[i];
			//printf("%d\t%f\n", k, calc_xi[k]);
			k++;
			l++;
			if (division_ele > 1) {
				temp1 = (knot_vec[i + 1] - knot_vec[i])
				        / (double)division_ele;
				for (j = 1; j < division_ele; j++) {
					calc_xi[k] = calc_xi[k - 1] + temp1;
					//printf("%d\t%f\n", k, calc_xi[k]);
					k++;
				}
			}
		}
	}
	calc_xi[k] = knot_vec[knot_n - 1];
	//printf("%d\t%f\n", k, calc_xi[k]);
	//printf("\n");
	point_n = k + 1;
	element_n = l;

	if (element_n > MAX_ELEMENTS) {
		printf("Error!!\n");
		printf("Too many elements at xi!\n"
		       "Maximum of elements is %d (Now %d)\n"
		       "\n", MAX_ELEMENTS, element_n);
		exit(1);
	}

	double temp_double;

	for (i = 0; i < cntl_p_n; i++) {
		for (j = 0; j < point_n; j++) {
			rBasisFunc(knot_vec, i, order, calc_xi[j],
			          &output[i][j], &doutput[i][j]);
			//lBasisFunc(knot_vec, i,
			//           cntl_p_n, order, calc_xi[j],
			//           &output[i][j], &doutput[i][j]);
		}
	}

	for (j = 0; j < point_n; j++) {
		NURBS_line(knot_vec, cntl_px, weight,
		           cntl_p_n, order, calc_xi[j],
		           &output_x[j], &temp_double);
		NURBS_line(knot_vec, cntl_py, weight,
		           cntl_p_n, order, calc_xi[j],
		           &output_y[j], &temp_double);
	}

//書き込み
	fp = fopen("BasisFunc.dat", "w");
	for (i = 0; i < cntl_p_n; i++) {
		for (j = 0; j < point_n; j++) {
			fprintf(fp, "% 1.8e % 1.8e % 1.8e\n",
			        calc_xi[j], output[i][j], doutput[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	fp = fopen("NURBSline.dat", "w");
	for (j = 0; j < point_n; j++) {
		fprintf(fp, "% 1.8e % 1.8e % 1.8e\n",
		        calc_xi[j], output_x[j], output_y[j]);
	}
	fclose(fp);
}

int main(int argc, char *argv[]) {
	if (argc != 3) {
		printf("Error!!\n");
		printf("Usage: ./NURBS_graph input.txt xi_div\n");
		printf("argc = %d\n", argc);
		exit(1);
	}

	division_ele  = atoi(argv[2]);

	if (division_ele > MAX_DIVISION) {
		printf("Error!!\n");
		printf("Too many Divsion at xi!\n"
		       "Maximum of division is %d (Now %d)\n"
		       "\n", MAX_DIVISION, division_ele);
		exit(1);
	}

	ReadFile (argv[1]);
	Calculation(order, knot_n,
	            cntl_p_n, knot_vec,
	            cntl_px, cntl_py,
	            weight);
	return 0;
}
