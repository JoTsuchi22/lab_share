#ifndef S_IGA_SUB_H
#define S_IGA_SUB_H

// gauss array
extern int GP_1dir;                            // 1方向のガウス点数
extern int GP_2D;                              // 2次元のガウス点数
extern double Gxi[POW_Ng_extended][DIMENSION]; // ガウス点
extern double w[POW_Ng_extended];              // ガウス点での重み



extern double Displacement[MAX_K_WHOLE_SIZE];
extern double Strain[MAX_N_ELEMENT][POW_Ng_extended][N_STRAIN];
extern double Strain_glo[MAX_N_ELEMENT][POW_Ng_extended][N_STRAIN];
extern double Strain_overlay[MAX_N_ELEMENT][POW_Ng_extended][N_STRAIN];
extern double Strain_overlay_loc[MAX_N_ELEMENT][POW_Ng_extended][N_STRAIN];
extern double Strain_aux_mode1[MAX_N_ELEMENT][POW_Ng_extended][N_STRAIN];
extern double Strain_aux_mode2[MAX_N_ELEMENT][POW_Ng_extended][N_STRAIN];
extern double Strain_aux_mode1_local[MAX_N_ELEMENT][POW_Ng_extended][N_STRAIN];
extern double Stress[MAX_N_ELEMENT][POW_Ng_extended][N_STRESS];
extern double Stress_glo[MAX_N_ELEMENT][POW_Ng_extended][N_STRESS];
extern double Stress_overlay[MAX_N_ELEMENT][POW_Ng_extended][N_STRESS];
extern double Stress_overlay_loc[MAX_N_ELEMENT][POW_Ng_extended][N_STRESS];
extern double Stress_aux_mode1_local[MAX_N_ELEMENT][POW_Ng_extended][N_STRESS];
extern double Stress_aux_mode1[MAX_N_ELEMENT][POW_Ng_extended][N_STRESS];
extern double Stress_aux_mode2_local[MAX_N_ELEMENT][POW_Ng_extended][N_STRESS];
extern double Stress_aux_mode2[MAX_N_ELEMENT][POW_Ng_extended][N_STRESS];
extern double StrainEnergyDensity[MAX_N_ELEMENT][POW_Ng_extended];
extern double StrainEnergyDensity_overlay[MAX_N_ELEMENT][POW_Ng_extended];
extern double StrainEnergyDensity_aux_mode1[MAX_N_ELEMENT][POW_Ng_extended];
extern double StrainEnergyDensity_aux_mode2[MAX_N_ELEMENT][POW_Ng_extended];
extern double StrainEnergyDensity_aux_mode1_loc[MAX_N_ELEMENT][POW_Ng_extended];
extern double StrainEnergyDensity_aux_only_mode1[MAX_N_ELEMENT][POW_Ng_extended];
extern double StrainEnergyDensity_aux_only_mode2[MAX_N_ELEMENT][POW_Ng_extended];
extern double StrainEnergyDensity_aux_only_mode1_local[MAX_N_ELEMENT][POW_Ng_extended];
extern double Disp_grad[MAX_N_ELEMENT][POW_Ng_extended][DIMENSION * DIMENSION]; // Disp_grad[MAX_N_ELEMENT][POW_Ng_extended][0] = 𝜕u1/𝜕x1  Disp_grad[MAX_N_ELEMENT][POW_Ng_extended][1] = 𝜕u1/𝜕x2  Disp_grad[MAX_N_ELEMENT][POW_Ng_extended][2] = 𝜕u2/𝜕x1 Disp_grad[MAX_N_ELEMENT][POW_Ng_extended][3] = 𝜕u2/𝜕x2
extern double Disp_grad_glo[MAX_N_ELEMENT][POW_Ng_extended][DIMENSION * DIMENSION];
extern double Disp_grad_overlay[MAX_N_ELEMENT][POW_Ng_extended][DIMENSION * DIMENSION];
extern double Disp_grad_aux_mode1_local[MAX_N_ELEMENT][POW_Ng_extended][DIMENSION * DIMENSION];
extern double Disp_grad_aux_mode1[MAX_N_ELEMENT][POW_Ng_extended][DIMENSION * DIMENSION];
extern double Disp_grad_aux_mode2_local[MAX_N_ELEMENT][POW_Ng_extended][DIMENSION * DIMENSION];
extern double Disp_grad_aux_mode2[MAX_N_ELEMENT][POW_Ng_extended][DIMENSION * DIMENSION];
extern double ReactionForce[MAX_K_WHOLE_SIZE];


extern int Constraint_ID[MAX_N_NODE * DIMENSION];


// extern int Check_BDBJ_flag[MAX_N_ELEMENT];
// extern int Total_BDBJ_flag;
// extern int Same_BDBJ_flag[POW_Ng_extended];


//重ね合わせの結果
extern int patch_n;                   // パッチ数
extern int cntl_p_n;                  // コントロールポイント数
extern int order_xi[MAX_PATCHES];     // ξ基底関数の次数(p)
extern int order_eta[MAX_PATCHES];    // η基底関数の次数(p)
extern int knot_n_xi[MAX_PATCHES];    // ξノットベクトルの数(n+p+1)
extern int knot_n_eta[MAX_PATCHES];   // ηノットベクトルの数(n+p+1)
extern int cntl_p_n_xi[MAX_PATCHES];  // ξ方向コントロールポイント数(n)
extern int cntl_p_n_eta[MAX_PATCHES]; // η方向コントロールポイント数(n)

extern double knot_vec_xi[MAX_PATCHES][MAX_KNOTS];   // ξノットベクトル
extern double knot_vec_eta[MAX_PATCHES][MAX_KNOTS];  // ηノットベクトル
extern double cntl_px[MAX_PATCHES][MAX_CNRL_P];      // コントロールポイントx座標
extern double cntl_py[MAX_PATCHES][MAX_CNRL_P];      // コントロールポイントy座標
extern double disp_cntl_px[MAX_PATCHES][MAX_CNRL_P]; // コントロールポイント上のx方向変位
extern double disp_cntl_py[MAX_PATCHES][MAX_CNRL_P]; // コントロールポイント上のy方向変位
extern double weight[MAX_PATCHES][MAX_CNRL_P];       // 重み

extern double output_xi_loc[MAX_ELEMENTS][Ng];
extern double output_eta_loc[MAX_ELEMENTS][Ng];
extern double coord_x[MAX_POINTS][MAX_POINTS];                                              // メッシュx座標
extern double coord_y[MAX_POINTS][MAX_POINTS];                                              // メッシュy座標
extern double coord_x_gauss[MAX_POINTS][MAX_POINTS];                                        // メッシュx座標 for gauss
extern double coord_y_gauss[MAX_POINTS][MAX_POINTS];                                        // メッシュy座標 for gauss
extern double dxi_x[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];        // ∂x/∂ξ
extern double dxi_y[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];        // ∂y/∂ξ
extern double deta_x[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];       // ∂x/∂η
extern double deta_y[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];       // ∂y/∂η
extern double dxi_x_gauss[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];  // ∂x/∂ξ for Gauss
extern double dxi_y_gauss[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];  // ∂y/∂ξ for Gauss
extern double deta_x_gauss[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]; // ∂x/∂η for Gauss
extern double deta_y_gauss[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]; // ∂y/∂η for Gauss

extern double disp_px_glo[MAX_ELEMENTS][MAX_N_Controlpoint_in_Patch];
extern double disp_py_glo[MAX_ELEMENTS][MAX_N_Controlpoint_in_Patch];
extern double disp_x[MAX_POINTS][MAX_POINTS];                                              // x方向変位
extern double disp_y[MAX_POINTS][MAX_POINTS];                                              // y方向変位
extern double disp_x_glo_gauss[MAX_ELEMENTS * Ng][MAX_ELEMENTS * Ng];                      // x方向変位 for Gauss
extern double disp_y_glo_gauss[MAX_ELEMENTS * Ng][MAX_ELEMENTS * Ng];                      // y方向変位 for Gauss
extern double dxi_disp_x[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];  // ∂u/∂ξ
extern double dxi_disp_y[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];  // ∂v/∂ξ
extern double deta_disp_x[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]; // ∂u/∂η
extern double deta_disp_y[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]; // ∂v/∂η

extern double strain_xx[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];       // x方向ひずみ
extern double strain_yy[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];       // y方向ひずみ
extern double strain_xy[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];       // 剪断ひずみ
extern double strain_xx_gauss[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]; // x方向ひずみ for Gauss
extern double strain_yy_gauss[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]; // y方向ひずみ for Gauss
extern double strain_xy_gauss[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]; // 剪断ひずみ for Gauss

extern double stress_xx[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];       // x方向垂直応力
extern double stress_yy[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];       // y方向垂直応力
extern double stress_xy[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];       // 剪断応力
extern double stress_xx_gauss[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]; // x方向垂直応力
extern double stress_yy_gauss[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]; // y方向垂直応力
extern double stress_xy_gauss[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1]; // 剪断応力

extern int fields_flag;  // s-IGAのためのNURBS_inputでは変位データは必ず読み込ませる
extern int division_ele_xi;  // ξ方向の一要素あたりの分割数
extern int division_ele_eta; // η方向の一要素あたりの分割数
extern int division_n_xi;    // ξ方向の表示する点の数
extern int division_n_eta;   // η方向の表示する点の数
extern int element_n_xi;     // ξ方向要素数
extern int element_n_eta;    // η方向要素数

extern int temp_index[MAX_PATCHES][MAX_CNRL_P];
extern double temp_cntl_px[MAX_CNRL_P];
extern double temp_cntl_py[MAX_CNRL_P];
extern double temp_weight[MAX_CNRL_P];
extern double temp_disp_x[MAX_CNRL_P];
extern double temp_disp_y[MAX_CNRL_P];

// for s-IGA
extern int n_patch_glo;  // グローバルメッシュ上のパッチ数
extern int n_patch_loc;  // ローカルメッシュ上のパッチ数
extern int glo_cntl_p_n; // グローバルメッシュ上のコントロールポイント数
extern int loc_cntl_p_n; // ローカルメッシュ上のコントロールポイント数

// for graph
extern int graph_patch_n; // グラフ作成用出力ファイル内のパッチ番号

// for s-IGA
extern double E;                      // ヤング率(GPa)
extern double nu;                     // ポアソン比(-)
extern int Total_mesh;


// extern int KIEL_SIZE;                           //要素分割マトリックスの大きさ
extern int D_MATRIX_SIZE;
extern int MAX_K_WHOLE_SIZE;
extern int K_Whole_Size;
extern int MAX_NON_ZERO;

// file pointer
extern FILE *fp;

#endif