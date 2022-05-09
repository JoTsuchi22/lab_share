#ifndef S_IGA_MAIN_H
#define S_IGA_MAIN_H

// gauss array
int GP_1dir;                            // 1方向のガウス点数
int GP_2D;                              // 2次元のガウス点数
double Gxi[POW_Ng_extended][DIMENSION]; // ガウス点
double w[POW_Ng_extended];              // ガウス点での重み

// for s-IGA
int Total_mesh;

double E;                      // ヤング率(GPa)
double nu;                     // ポアソン比(-)

int D_MATRIX_SIZE = 0;
int MAX_K_WHOLE_SIZE;
int K_Whole_Size;
int MAX_NON_ZERO;

// file pointer
FILE *fp;

#endif