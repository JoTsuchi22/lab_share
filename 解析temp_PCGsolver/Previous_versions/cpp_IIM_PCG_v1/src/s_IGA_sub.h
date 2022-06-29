#ifndef S_IGA_SUB_H
#define S_IGA_SUB_H

// gauss array
extern int GP_1dir;                            // 1方向のガウス点数
extern int GP_2D;                              // 2次元のガウス点数
extern double Gxi[POW_Ng_extended][DIMENSION]; // ガウス点
extern double w[POW_Ng_extended];              // ガウス点での重み

// extern double Strain[MAX_N_ELEMENT][POW_Ng_extended][N_STRAIN];
// extern double Strain_glo[MAX_N_ELEMENT][POW_Ng_extended][N_STRAIN];
// extern double Strain_overlay[MAX_N_ELEMENT][POW_Ng_extended][N_STRAIN];
// extern double Strain_overlay_loc[MAX_N_ELEMENT][POW_Ng_extended][N_STRAIN];
// extern double Strain_aux_mode1[MAX_N_ELEMENT][POW_Ng_extended][N_STRAIN];
// extern double Strain_aux_mode2[MAX_N_ELEMENT][POW_Ng_extended][N_STRAIN];
// extern double Strain_aux_mode1_local[MAX_N_ELEMENT][POW_Ng_extended][N_STRAIN];
// extern double Stress[MAX_N_ELEMENT][POW_Ng_extended][N_STRESS];
// extern double Stress_glo[MAX_N_ELEMENT][POW_Ng_extended][N_STRESS];
// extern double Stress_overlay[MAX_N_ELEMENT][POW_Ng_extended][N_STRESS];
// extern double Stress_overlay_loc[MAX_N_ELEMENT][POW_Ng_extended][N_STRESS];
// extern double Stress_aux_mode1_local[MAX_N_ELEMENT][POW_Ng_extended][N_STRESS];
// extern double Stress_aux_mode1[MAX_N_ELEMENT][POW_Ng_extended][N_STRESS];
// extern double Stress_aux_mode2_local[MAX_N_ELEMENT][POW_Ng_extended][N_STRESS];
// extern double Stress_aux_mode2[MAX_N_ELEMENT][POW_Ng_extended][N_STRESS];
// extern double StrainEnergyDensity[MAX_N_ELEMENT][POW_Ng_extended];
// extern double StrainEnergyDensity_overlay[MAX_N_ELEMENT][POW_Ng_extended];
// extern double StrainEnergyDensity_aux_mode1[MAX_N_ELEMENT][POW_Ng_extended];
// extern double StrainEnergyDensity_aux_mode2[MAX_N_ELEMENT][POW_Ng_extended];
// extern double StrainEnergyDensity_aux_mode1_loc[MAX_N_ELEMENT][POW_Ng_extended];
// extern double StrainEnergyDensity_aux_only_mode1[MAX_N_ELEMENT][POW_Ng_extended];
// extern double StrainEnergyDensity_aux_only_mode2[MAX_N_ELEMENT][POW_Ng_extended];
// extern double StrainEnergyDensity_aux_only_mode1_local[MAX_N_ELEMENT][POW_Ng_extended];
// extern double Disp_grad[MAX_N_ELEMENT][POW_Ng_extended][DIMENSION * DIMENSION]; // Disp_grad[MAX_N_ELEMENT][POW_Ng_extended][0] = 𝜕u1/𝜕x1  Disp_grad[MAX_N_ELEMENT][POW_Ng_extended][1] = 𝜕u1/𝜕x2  Disp_grad[MAX_N_ELEMENT][POW_Ng_extended][2] = 𝜕u2/𝜕x1 Disp_grad[MAX_N_ELEMENT][POW_Ng_extended][3] = 𝜕u2/𝜕x2
// extern double Disp_grad_glo[MAX_N_ELEMENT][POW_Ng_extended][DIMENSION * DIMENSION];
// extern double Disp_grad_overlay[MAX_N_ELEMENT][POW_Ng_extended][DIMENSION * DIMENSION];
// extern double Disp_grad_aux_mode1_local[MAX_N_ELEMENT][POW_Ng_extended][DIMENSION * DIMENSION];
// extern double Disp_grad_aux_mode1[MAX_N_ELEMENT][POW_Ng_extended][DIMENSION * DIMENSION];
// extern double Disp_grad_aux_mode2_local[MAX_N_ELEMENT][POW_Ng_extended][DIMENSION * DIMENSION];
// extern double Disp_grad_aux_mode2[MAX_N_ELEMENT][POW_Ng_extended][DIMENSION * DIMENSION];
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

extern int D_MATRIX_SIZE;
extern int MAX_K_WHOLE_SIZE;
extern int K_Whole_Size;
extern int MAX_NON_ZERO;
extern int DIVISION_ELEMENT;

// file pointer
extern FILE *fp;

#endif