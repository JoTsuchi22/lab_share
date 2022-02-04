#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[]) {
	int xi  = atoi(argv[1]); //ξ方向分割数
	int eta = atoi(argv[2]); //η方向分割数
	double r1 = 1.0; //内径
	double r2 = 2.0; //外形

	int i, j;
	double comp_xi;
	double comp_eta;

	//ヤング率　ポアソン比
	printf("206000  0.3\n\n");

	//パッチ数
	printf("1\n");
	printf("\n");

	//全体のコントロールポイント数
	printf("%d\n", (2 * xi + 1) * (2 * eta + 1));
	printf("\n");

	//次元
	printf("2\t2\n");
	printf("\n");

	//ノットベクトル数
	printf("%d\t%d\n", 2 * xi + 4, 2 * eta + 4);
	printf("\n");

	//コントロールポイント数
	printf("%d\t%d\n", 2 * xi + 1, 2 * eta + 1);
	printf("\n");

	//使用するコントロール番号
	for (i = 0; i < (2 * xi + 1) * (2 * eta + 1); i++) {
		printf("%d ", i);
	}
	printf("\n\n");

	//境界条件数
	printf("%d\t0\t1\n", 4 * xi + 2);
	printf("\n");

	//ノットベクトル
	printf("0\t0\t0\t");
	for (i = 0; i < xi - 1; i++) {
		comp_xi = ((double)i + 1.0) / (double)xi;
		printf("%.15f\t%.15f\t", comp_xi, comp_xi);
	}
	printf("1\t1\t1\n");
	printf("0\t0\t0\t");
	for (i = 0; i < eta - 1; i++) {
		comp_eta = ((double)i + 1.0) / (double)eta;
		printf("%.15f\t%.15f\t", comp_eta, comp_eta);
	}
	printf("1\t1\t1\n");
	printf("\n");

//コントロールポイント座標
	double coor_xi;
	double coor_eta;
	double r;
	double theta = M_PI / (2.0 * (double)eta);
	double theta_calc;

	int n = 0;

	for (i = 0; i < eta; i++) {
		for (j = 0; j < 2 * xi + 1; j++) {
			r = (r2 - r1) * (double)j / (double)(2 * xi) + r1;
			theta_calc = (double)i * theta;
			coor_xi =  r * cos(theta_calc);
			coor_eta = r * sin(theta_calc);
			printf("%d\t%.15f\t%.15f\t1\n", n, coor_xi, coor_eta);
			n++;
		}
		for (j = 0; j < 2 * xi + 1; j++) {
			r = (r2 - r1) * (double)j / (double)(2 * xi) + r1;
			theta_calc = ((double)i + 0.5) * theta;
			coor_xi =  r / cos(theta / 2.0) * cos(theta_calc);
			coor_eta = r / cos(theta / 2.0) * sin(theta_calc);
			printf("%d\t%.15f\t%.15f\t%.15f\n",
			       n, coor_xi, coor_eta, cos(theta / 2.0));
			n++;
		}
	}
	//theta_calc=π/2
	for (j = 0; j < 2 * xi + 1; j++) {
		r = (r2 - r1) * (double)j / (double)(2 * xi) + r1;
		theta_calc = (double)eta * theta;
		coor_xi =  r * cos(theta_calc);
		coor_eta = r * sin(theta_calc);
		printf("%d\t%.15f\t%.15f\t1\n", n, coor_xi, coor_eta);
		n++;
	}
	printf("\n");

	//変位固定
	for (i = 0; i < 2 * xi + 1; i++) {
		printf("%d\t1\t0\n", i);
	}
	for (i = 2 * eta * (2 * xi + 1); i < (2 * eta + 1) * (2 * xi + 1); i++) {
		printf("%d\t0\t0\n", i);
	}
	printf("\n");

	//集中荷重(なし)

	//エッジ荷重
	printf("2\t0\t1\t0\t0\t1\t1\t0\t0\n");

	return 0;
}
