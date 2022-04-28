#ifndef S_IGA_SUB_H
#define S_IGA_SUB_H

// gauss array
extern int GP_1dir;                            // 1方向のガウス点数
extern int GP_2D;                              // 2次元のガウス点数
extern double Gxi[POW_Ng_extended][DIMENSION]; // ガウス点
extern double w[POW_Ng_extended];              // ガウス点での重み

// for s-IGA
extern int Total_mesh;

//重ね合わせの結果
extern double E;                      // ヤング率(GPa)
extern double nu;                     // ポアソン比(-)
extern FILE *fp;

#endif