#ifndef S_IGA_MAIN_H
#define S_IGA_MAIN_H

// gauss array
int GP_1D;                                    // 1方向のガウス点数
int GP_2D;                                    // 2次元のガウス点数
int GP_3D;                                    // 3次元のガウス点数
int GP_ON_ELEMENT;                            // 要素内のガウス点数
double Gxi_1D[MAX_POW_NG];                    // 1次元のガウス点
double w_1D[MAX_POW_NG];                      // 1次元のガウス点での重み
double Gxi[MAX_POW_NG_EXTEND][MAX_DIMENSION]; // ガウス点
double w[MAX_POW_NG_EXTEND];                  // ガウス点での重み

// for viewer
int division_ele_xi;  // ξ方向の一要素あたりの分割数
int division_ele_eta; // η方向の一要素あたりの分割数
int division_n_xi;    // ξ方向の表示する点の数
int division_n_eta;   // η方向の表示する点の数
int element_n_xi;     // ξ方向要素数
int element_n_eta;    // η方向要素数

// for s-IGA
double E;             // ヤング率
double nu;            // ポアソン比
int Total_mesh;

int MAX_ORDER = 0; // 基底関数の次数の最大値 + 1
int MAX_CP = 0;
int MAX_NO_CP_ON_ELEMENT;
int MAX_KIEL_SIZE;

int D_MATRIX_SIZE;
int N_STRAIN;
int N_STRESS;
int MAX_K_WHOLE_SIZE;
int K_Whole_Size;

int Total_connectivity;
int Total_connectivity_glo;
int DIVISION_ELEMENT;

// file pointer
FILE *fp;

#endif