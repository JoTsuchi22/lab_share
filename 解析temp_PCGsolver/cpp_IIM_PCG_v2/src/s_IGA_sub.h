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

// extern double Strain[MAX_N_ELEMENT][POW_NG_EXTEND][N_STRAIN];
// extern double Strain_glo[MAX_N_ELEMENT][POW_NG_EXTEND][N_STRAIN];
// extern double Strain_overlay[MAX_N_ELEMENT][POW_NG_EXTEND][N_STRAIN];
// extern double Strain_overlay_loc[MAX_N_ELEMENT][POW_NG_EXTEND][N_STRAIN];
// extern double Strain_aux_mode1[MAX_N_ELEMENT][POW_NG_EXTEND][N_STRAIN];
// extern double Strain_aux_mode2[MAX_N_ELEMENT][POW_NG_EXTEND][N_STRAIN];
// extern double Strain_aux_mode1_local[MAX_N_ELEMENT][POW_NG_EXTEND][N_STRAIN];
// extern double Stress[MAX_N_ELEMENT][POW_NG_EXTEND][N_STRESS];
// extern double Stress_glo[MAX_N_ELEMENT][POW_NG_EXTEND][N_STRESS];
// extern double Stress_overlay[MAX_N_ELEMENT][POW_NG_EXTEND][N_STRESS];
// extern double Stress_overlay_loc[MAX_N_ELEMENT][POW_NG_EXTEND][N_STRESS];
// extern double Stress_aux_mode1_local[MAX_N_ELEMENT][POW_NG_EXTEND][N_STRESS];
// extern double Stress_aux_mode1[MAX_N_ELEMENT][POW_NG_EXTEND][N_STRESS];
// extern double Stress_aux_mode2_local[MAX_N_ELEMENT][POW_NG_EXTEND][N_STRESS];
// extern double Stress_aux_mode2[MAX_N_ELEMENT][POW_NG_EXTEND][N_STRESS];
// extern double StrainEnergyDensity[MAX_N_ELEMENT][POW_NG_EXTEND];
// extern double StrainEnergyDensity_overlay[MAX_N_ELEMENT][POW_NG_EXTEND];
// extern double StrainEnergyDensity_aux_mode1[MAX_N_ELEMENT][POW_NG_EXTEND];
// extern double StrainEnergyDensity_aux_mode2[MAX_N_ELEMENT][POW_NG_EXTEND];
// extern double StrainEnergyDensity_aux_mode1_loc[MAX_N_ELEMENT][POW_NG_EXTEND];
// extern double StrainEnergyDensity_aux_only_mode1[MAX_N_ELEMENT][POW_NG_EXTEND];
// extern double StrainEnergyDensity_aux_only_mode2[MAX_N_ELEMENT][POW_NG_EXTEND];
// extern double StrainEnergyDensity_aux_only_mode1_local[MAX_N_ELEMENT][POW_NG_EXTEND];
// extern double Disp_grad[MAX_N_ELEMENT][POW_NG_EXTEND][DIMENSION * DIMENSION]; // Disp_grad[MAX_N_ELEMENT][POW_NG_EXTEND][0] = 𝜕u1/𝜕x1  Disp_grad[MAX_N_ELEMENT][POW_NG_EXTEND][1] = 𝜕u1/𝜕x2  Disp_grad[MAX_N_ELEMENT][POW_NG_EXTEND][2] = 𝜕u2/𝜕x1 Disp_grad[MAX_N_ELEMENT][POW_NG_EXTEND][3] = 𝜕u2/𝜕x2
// extern double Disp_grad_glo[MAX_N_ELEMENT][POW_NG_EXTEND][DIMENSION * DIMENSION];
// extern double Disp_grad_overlay[MAX_N_ELEMENT][POW_NG_EXTEND][DIMENSION * DIMENSION];
// extern double Disp_grad_aux_mode1_local[MAX_N_ELEMENT][POW_NG_EXTEND][DIMENSION * DIMENSION];
// extern double Disp_grad_aux_mode1[MAX_N_ELEMENT][POW_NG_EXTEND][DIMENSION * DIMENSION];
// extern double Disp_grad_aux_mode2_local[MAX_N_ELEMENT][POW_NG_EXTEND][DIMENSION * DIMENSION];
// extern double Disp_grad_aux_mode2[MAX_N_ELEMENT][POW_NG_EXTEND][DIMENSION * DIMENSION];
// extern double ReactionForce[MAX_K_WHOLE_SIZE];

extern int division_ele_xi;  // ξ方向の一要素あたりの分割数
extern int division_ele_eta; // η方向の一要素あたりの分割数
extern int division_n_xi;    // ξ方向の表示する点の数
extern int division_n_eta;   // η方向の表示する点の数
extern int element_n_xi;     // ξ方向要素数
extern int element_n_eta;    // η方向要素数

// for s-IGA
extern double E;                      // ヤング率(GPa)
extern double nu;                     // ポアソン比(-)
extern int Total_mesh;

extern int MAX_ORDER;          // 基底関数の次数の最大値 + 1
extern int MAX_CP;
extern int MAX_NO_CP_ON_ELEMENT;  // new
extern int MAX_KIEL_SIZE;              // new

extern int MAX_ORDER;          // 基底関数の次数の最大値 + 1
extern int D_MATRIX_SIZE;
extern int MAX_K_WHOLE_SIZE;
extern int K_Whole_Size;
extern int MAX_NON_ZERO;
extern int DIVISION_ELEMENT;

// file pointer
extern FILE *fp;

#endif