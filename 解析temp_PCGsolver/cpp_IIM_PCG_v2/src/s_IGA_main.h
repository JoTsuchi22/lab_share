#ifndef S_IGA_MAIN_H
#define S_IGA_MAIN_H

// gauss array
int GP_1dir;                            // 1方向のガウス点数
int GP_2D;                              // 2次元のガウス点数
double Gxi[POW_NG_EXTEND][DIMENSION]; // ガウス点
double w[POW_NG_EXTEND];              // ガウス点での重み

// double Strain[MAX_N_ELEMENT][POW_NG_EXTEND][N_STRAIN];
// double Strain_glo[MAX_N_ELEMENT][POW_NG_EXTEND][N_STRAIN];
// double Strain_overlay[MAX_N_ELEMENT][POW_NG_EXTEND][N_STRAIN];
// double Strain_overlay_loc[MAX_N_ELEMENT][POW_NG_EXTEND][N_STRAIN];
// double Strain_aux_mode1[MAX_N_ELEMENT][POW_NG_EXTEND][N_STRAIN];
// double Strain_aux_mode2[MAX_N_ELEMENT][POW_NG_EXTEND][N_STRAIN];
// double Strain_aux_mode1_local[MAX_N_ELEMENT][POW_NG_EXTEND][N_STRAIN];
// double Stress[MAX_N_ELEMENT][POW_NG_EXTEND][N_STRESS];
// double Stress_glo[MAX_N_ELEMENT][POW_NG_EXTEND][N_STRESS];
// double Stress_overlay[MAX_N_ELEMENT][POW_NG_EXTEND][N_STRESS];
// double Stress_overlay_loc[MAX_N_ELEMENT][POW_NG_EXTEND][N_STRESS];
// double Stress_aux_mode1_local[MAX_N_ELEMENT][POW_NG_EXTEND][N_STRESS];
// double Stress_aux_mode1[MAX_N_ELEMENT][POW_NG_EXTEND][N_STRESS];
// double Stress_aux_mode2_local[MAX_N_ELEMENT][POW_NG_EXTEND][N_STRESS];
// double Stress_aux_mode2[MAX_N_ELEMENT][POW_NG_EXTEND][N_STRESS];
// double StrainEnergyDensity[MAX_N_ELEMENT][POW_NG_EXTEND];
// double StrainEnergyDensity_overlay[MAX_N_ELEMENT][POW_NG_EXTEND];
// double StrainEnergyDensity_aux_mode1[MAX_N_ELEMENT][POW_NG_EXTEND];
// double StrainEnergyDensity_aux_mode2[MAX_N_ELEMENT][POW_NG_EXTEND];
// double StrainEnergyDensity_aux_mode1_loc[MAX_N_ELEMENT][POW_NG_EXTEND];
// double StrainEnergyDensity_aux_only_mode1[MAX_N_ELEMENT][POW_NG_EXTEND];
// double StrainEnergyDensity_aux_only_mode2[MAX_N_ELEMENT][POW_NG_EXTEND];
// double StrainEnergyDensity_aux_only_mode1_local[MAX_N_ELEMENT][POW_NG_EXTEND];
// double Disp_grad[MAX_N_ELEMENT][POW_NG_EXTEND][DIMENSION * DIMENSION]; // Disp_grad[MAX_N_ELEMENT][POW_NG_EXTEND][0] = 𝜕u1/𝜕x1  Disp_grad[MAX_N_ELEMENT][POW_NG_EXTEND][1] = 𝜕u1/𝜕x2  Disp_grad[MAX_N_ELEMENT][POW_NG_EXTEND][2] = 𝜕u2/𝜕x1 Disp_grad[MAX_N_ELEMENT][POW_NG_EXTEND][3] = 𝜕u2/𝜕x2
// double Disp_grad_glo[MAX_N_ELEMENT][POW_NG_EXTEND][DIMENSION * DIMENSION];
// double Disp_grad_overlay[MAX_N_ELEMENT][POW_NG_EXTEND][DIMENSION * DIMENSION];
// double Disp_grad_aux_mode1_local[MAX_N_ELEMENT][POW_NG_EXTEND][DIMENSION * DIMENSION];
// double Disp_grad_aux_mode1[MAX_N_ELEMENT][POW_NG_EXTEND][DIMENSION * DIMENSION];
// double Disp_grad_aux_mode2_local[MAX_N_ELEMENT][POW_NG_EXTEND][DIMENSION * DIMENSION];
// double Disp_grad_aux_mode2[MAX_N_ELEMENT][POW_NG_EXTEND][DIMENSION * DIMENSION];
// double ReactionForce[MAX_K_WHOLE_SIZE];

int division_ele_xi;  // ξ方向の一要素あたりの分割数
int division_ele_eta; // η方向の一要素あたりの分割数
int division_n_xi;    // ξ方向の表示する点の数
int division_n_eta;   // η方向の表示する点の数
int element_n_xi;     // ξ方向要素数
int element_n_eta;    // η方向要素数

// for s-IGA
double E;                      // ヤング率(GPa)
double nu;                     // ポアソン比(-)
int Total_mesh;

int D_MATRIX_SIZE;
int MAX_K_WHOLE_SIZE;
int K_Whole_Size;
int MAX_NON_ZERO;
int DIVISION_ELEMENT;

// file pointer
FILE *fp;

#endif