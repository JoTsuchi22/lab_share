#ifndef S_IGA_SUB_H
#define S_IGA_SUB_H

// gauss array
extern int GP_1D;                                    // 1方向のガウス点数
extern int GP_2D;                                    // 2次元のガウス点数
extern int GP_3D;                                    // 3次元のガウス点数
extern int GP_ON_ELEMENT;                            // 要素内のガウス点数
extern double Gxi_1D[MAX_POW_NG];                    // 1次元のガウス点
extern double w_1D[MAX_POW_NG];                      // 1次元のガウス点での重み
extern double Gxi[MAX_POW_NG_EXTEND][MAX_DIMENSION]; // ガウス点
extern double w[MAX_POW_NG_EXTEND];                  // ガウス点での重み

// for viewer
extern int division_ele_xi;  // ξ方向の一要素あたりの分割数
extern int division_ele_eta; // η方向の一要素あたりの分割数
extern int division_n_xi;    // ξ方向の表示する点の数
extern int division_n_eta;   // η方向の表示する点の数
extern int element_n_xi;     // ξ方向要素数
extern int element_n_eta;    // η方向要素数

// for s-IGA
extern double E;             // ヤング率
extern double nu;            // ポアソン比
extern int Total_mesh;

extern int MAX_ORDER; // 基底関数の次数の最大値 + 1
extern int MAX_CP;
extern int MAX_NO_CP_ON_ELEMENT;
extern int MAX_KIEL_SIZE;

extern int D_MATRIX_SIZE;
extern int N_STRAIN;
extern int N_STRESS;
extern int MAX_K_WHOLE_SIZE;
extern int K_Whole_Size;

extern int Total_connectivity;
extern int Total_connectivity_glo;
extern int DIVISION_ELEMENT;

// file pointer
extern FILE *fp;

#endif