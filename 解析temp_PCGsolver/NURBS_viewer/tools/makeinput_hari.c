#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
	int xi  = atoi(argv[1]);   //ξ方向分割数
	int eta = atoi(argv[2]);   //η方向分割数
	double l = 1.0;	//梁の長さ
	double d = 2.0;	//梁の幅

	int i, j;
	double comp_xi;
	double comp_eta;

	//ヤング率　ポアソン比
	printf("10  0.4\n\n");

	//パッチ数
	printf("1\n");
	printf("\n");

	//全体のコントロールポイント数
	printf("%d\n", (2 * xi + 1) * (2 * eta + 1));
	printf("\n");

	//次元
	printf("2   2\n");
	printf("\n");

	//ノットベクトル数
	printf("%d   %d\n", 2 * xi + 4, 2 * eta + 4);
	printf("\n");

	//コントロールポイント数
	printf("%d   %d\n", 2 * xi + 1, 2 * eta + 1);
	printf("\n");

	//使用するコントロールポイント番号
	for (i = 0; i < (2 * xi + 1) * (2 * eta + 1); i++) {
		printf("%d ", i);
	}
	printf("\n\n");

	//境界条件数
	//printf("%d   0   1\n", 4 * eta + 2);
	printf("%d   0   1\n", 2 * eta + 1 + 2 * xi + 1);
	printf("\n");

	//ノットベクトル数
	printf("0   0   0   ");
	for (i = 0; i < xi - 1; i++) {
		comp_xi = ((double)i + 1.0) / (double)xi;
		printf("%f   %f   ", comp_xi, comp_xi);
	}
	printf("1   1   1\n");

	printf("0   0   0   ");
	for (i = 0; i < eta - 1; i++) {
		comp_eta = ((double)i + 1.0) / (double)eta;
		printf("%f %f ", comp_eta, comp_eta);
	}
	printf("1   1   1\n");
	printf("\n");

	//コントロールポイント座標
	for (i = 0; i < 2 * eta + 1; i++) {
		comp_eta  = d * (double)i / (2.0 * (double)eta);
		for (j = 0; j < 2 * xi + 1; j++) {
			comp_xi =  l * (double)j / (2.0 * (double)xi);
			printf("%d\t%f\t%f\t1\n", (2 * xi + 1)*i + j, comp_xi, comp_eta);
		}
	}
	printf("\n");

	//変位固定
	for (i = 0; i < 2 * eta + 1; i++) {
		printf("%d   0   0\n", i * (2 * xi + 1));
		//printf("%d   1   0\n", i * (2 * xi + 1));
	}
	for (i = 0; i < 2 * xi + 1; i++) {
		//printf("%d   0   0\n", i * (2 * xi + 1));
		printf("%d   1   0\n", i);
	}
	printf("\n");

	//集中荷重
	//printf("440   1   -10   0   0\n");

	//エッジ荷重
	//printf("1   0   0   1   0   1   -10   0   0\n");
	printf("0   0   1   1.0   0.0   1.0   1.0   0.0   0.0\n");

	return 0;
}
