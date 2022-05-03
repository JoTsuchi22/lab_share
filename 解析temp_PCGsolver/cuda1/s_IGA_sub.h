#ifndef S_IGA_SUB_H
#define S_IGA_SUB_H

// gauss array
extern int GP_1dir;                            // 1方向のガウス点数
extern int GP_2D;                              // 2次元のガウス点数
extern double Gxi[POW_Ng_extended][DIMENSION]; // ガウス点
extern double w[POW_Ng_extended];              // ガウス点での重み

extern int KIEL_SIZE;                           //要素分割マトリックスの大きさ
// extern int Controlpoint_of_Element[MAX_N_ELEMENT][MAX_NO_CCpoint_ON_ELEMENT];
// extern double Equivalent_Nodal_Force[MAX_N_NODE][DIMENSION]; // Equivalent nodal forces arising from the distributed load
extern int K_Whole_Ptr[MAX_K_WHOLE_SIZE + 1], K_Whole_Col[MAX_NON_ZERO];
extern double K_Whole_Val[MAX_NON_ZERO];
extern int Index_Dof[MAX_K_WHOLE_SIZE];
// extern int INC[MAX_N_PATCH][MAX_N_ELEMENT][DIMENSION];
// extern int Adress_Controlpoint[MAX_N_PATCH][1000][1000]; // INCの配列をいじったものAdress_Controlpoint[ξ][η]；コントールポイント番号、任意のパッチ上でξ方向[]番目、η方向[]番目のコントロールポイント番号を示す
extern double element_coordinate_Nopoint[MAX_N_ELEMENT][DIMENSION];
extern double Gausspoint_coordinates[MAX_N_ELEMENT][POW_Ng_extended][DIMENSION];
extern int same_point_in_Element[MAX_N_NODE];
// extern int Element_patch[MAX_N_ELEMENT];                                 // 要素がどのパッチに属しているか示す配列(要素番号は1つのモデルで通し番号)

extern int Node_To_Node[K_DIVISION_LENGE][10000], Total_Control_Point_To_Node[K_DIVISION_LENGE]; //ある節点に関係する節点番号s

extern double sol_vec[MAX_K_WHOLE_SIZE];
extern double rhs_vec[MAX_K_WHOLE_SIZE];
extern double diag_scaling[MAX_K_WHOLE_SIZE];

// extern double Shape[DIMENSION][MAX_N_NODE][10];
// extern double shape_func[MAX_N_NODE];
// extern double dShape_func1[MAX_N_NODE];
// extern double dShape_func2[MAX_N_NODE];
// extern double dShape[DIMENSION][MAX_N_NODE];
// extern double Position_Data_param[DIMENSION];

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

// extern int real_Total_Element;                                       // ゼロエレメントを除いた要素数
// extern double difference[MAX_N_PATCH][MAX_N_KNOT][DIMENSION];        // 隣り合うノットベクトルの差
// extern int ENC[MAX_N_PATCH][MAX_N_ELEMENT][DIMENSION];               // ENC[パッチ][全ての要素][0, 1] = x, y方向の何番目の要素か
// extern int real_element[MAX_N_ELEMENT];                              // ゼロエレメントではない要素の番号
// extern int Total_element_all_ID[MAX_N_ELEMENT];                      // ゼロエレメントではない要素 = 1, ゼロエレメント = 0
// extern int line_No_Total_element[MAX_N_PATCH][DIMENSION];            // ゼロエレメントを含むすべての要素列の数
// extern int line_No_real_element[MAX_N_PATCH][DIMENSION];             // ゼロエレメントではない要素列の数
// extern int real_element_line[MAX_N_PATCH][MAX_N_ELEMENT][DIMENSION]; // ゼロエレメントではない要素列

extern int No_points_for_colored_points;
extern int No_points_for_new_zarusoba;

// for s-IGA
extern int Total_mesh;

// extern int Element_mesh[MAX_N_ELEMENT]; // 要素がどのメッシュ内にあるかを示す配列
extern int Patch_mesh[MAX_N_PATCH];     // パッチがどのメッシュ内にあるかを示す配列

// extern int Total_Patch_on_mesh[MAX_N_MESH];     // 各メッシュ上のパッチ数
// extern int Total_Patch_to_mesh[MAX_N_MESH + 1]; // メッシュ[]までのパッチ数(メッシュ[]内のパッチ数は含まない)
// extern int Total_Patch_to_Now;                  // 現メッシュまでのパッチ数(現メッシュのパッチ数は含まない)

// extern int Total_Control_Point_on_mesh[MAX_N_MESH];     // 各メッシュ上のコントロールポイント数
// extern int Total_Control_Point_to_mesh[MAX_N_MESH + 1]; // メッシュ[]までのコントロールポイント数(メッシュ[]内のコントロールポイント数は含まない)
// extern int Total_Control_Point_to_Now;                  // 現メッシュまでのコントロールポイント数(現メッシュのコントロールポイント数は含まない)

// extern int Total_Element_on_mesh[MAX_N_MESH];
// extern int Total_Element_to_mesh[MAX_N_MESH + 1];
// extern int Total_Element_to_Now;

// extern int Total_Constraint_all_mesh;
// extern int Total_Constraint_on_mesh[MAX_N_MESH];
// extern int Total_Constraint_to_mesh[MAX_N_MESH + 1];
// extern int Total_Load_on_mesh[MAX_N_MESH];
// extern int Total_Load_to_mesh[MAX_N_MESH + 1];
// extern int Total_DistributeForce_on_mesh[MAX_N_MESH];
// extern int Total_DistributeForce_to_mesh[MAX_N_MESH + 1];
// extern double Node_Coordinate[MAX_N_NODE][DIMENSION + 1];
// extern int Order[MAX_N_PATCH][DIMENSION];
// extern int No_knot[MAX_N_PATCH][DIMENSION];
// extern int No_Control_point[MAX_N_PATCH][DIMENSION];
// extern int Patch_controlpoint[MAX_N_PATCH][MAX_N_Controlpoint_in_Patch]; // バッチとコントロールポイント番号の要素コネクティビティ
// extern int No_Controlpoint_in_patch[MAX_N_PATCH];
// extern int No_Control_point_ON_ELEMENT[10000];

// extern double Control_Coord[DIMENSION][MAX_N_NODE];
// extern double Control_Weight[MAX_N_NODE];

// extern double Position_Knots[MAX_N_PATCH][DIMENSION][MAX_N_KNOT];
// extern int Constraint_Node_Dir_on_mesh[MAX_N_MESH][MAX_N_CONSTRAINT][2];
// extern double Value_of_Constraint_on_mesh[MAX_N_MESH][MAX_N_CONSTRAINT];

// extern int iPatch_array[MAX_N_DISTRIBUTE_FORCE], iCoord_array[MAX_N_DISTRIBUTE_FORCE], type_load_array[MAX_N_DISTRIBUTE_FORCE];
// extern double val_Coord_array[MAX_N_DISTRIBUTE_FORCE], Range_Coord_array[MAX_N_DISTRIBUTE_FORCE][2], Coeff_Dist_Load_array[MAX_N_DISTRIBUTE_FORCE][3];

// extern int real_Total_Element_to_Now;

// extern int El_No_on_mesh[MAX_N_MESH][MAX_N_ELEMENT]; // メッシュ内でのコントロールポイント配列
extern int Constraint_ID[MAX_N_NODE * DIMENSION];

// extern int real_Total_Element_on_mesh[MAX_N_MESH];
// extern int real_Total_Element_to_mesh[MAX_N_MESH + 1];
// extern int real_El_No_on_mesh[MAX_N_MESH][MAX_N_ELEMENT];

extern int temp_element_n[MAX_N_ELEMENT_OVER_POINT];
extern int element_n_point[MAX_N_ELEMENT_OVER_ELEMENT];
extern int NNLOVER[MAX_N_ELEMENT];
extern int NELOVER[MAX_N_ELEMENT][MAX_N_ELEMENT_OVER];
extern int Check_BDBJ_flag[MAX_N_ELEMENT];
extern int Total_BDBJ_flag;
extern int Same_BDBJ_flag[POW_Ng_extended];


// for Interaction integral
extern double T[DIMENSION][DIMENSION];
extern double K_mode1;
extern double K_mode2;
extern double J_integral_value_aux_mode1;
extern double J_integral_value_aux_mode2;

//重ね合わせの結果
extern double E;                      // ヤング率(GPa)
extern double nu;                     // ポアソン比(-)
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

// for GP info
//  extern double coordinate_GP[MAX_ELEMENTS*MAX_ELEMENTS][POW_Ng_extended][DIMENSION];
//  extern double strain_GP[MAX_ELEMENTS*MAX_ELEMENTS][POW_Ng_extended][3];
//  extern double stress_GP[MAX_ELEMENTS*MAX_ELEMENTS][POW_Ng_extended][3];
//  extern double stress_r_theta_GP[MAX_ELEMENTS*MAX_ELEMENTS][POW_Ng_extended][3];
//  extern double stress_theory_r_theta[MAX_ELEMENTS*MAX_ELEMENTS][POW_Ng_extended][3];
//  extern double Jac[MAX_ELEMENTS*MAX_ELEMENTS][POW_Ng_extended];

// 解析条件パラメータの設定
// extern int DM;                   // 平面応力状態:DM=0	平面ひずみ状態:DM=1
// extern int check_over_parameter; // 要素の重なりの判定(要素の物体上の端点:0 ガウス点:1)

extern int n_patch;

extern FILE *fp;

#endif