#ifndef S_IGA_HEADER_H
#define S_IGA_HEADER_H

#define SKIP_S_IGA 2 // 重ね合わせとJ積分を行う 0, 重ね合わせをスキップしてJ積分を行う 1, J積分を行わない 2
#define DM           // 平面応力状態:DM=0	平面ひずみ状態:DM=1
#define ERROR -999
#define PI  3.14159265359
#define MAX_NO_CCpoint_ON_ELEMENT 16						// 分割節点数
#define DIMENSION 2											// 次元数
#define MAX_KIEL_SIZE MAX_NO_CCpoint_ON_ELEMENT * DIMENSION	// 要素分割マトリックスの大きさ
#define Ng 4												// Gauss-Legendreの足す回数
#define POW_Ng Ng * Ng										// NgのDIMENSION乗の計算
#define Ng_extended 10										// Gauss-Legendreの足す回数
#define POW_Ng_extended Ng_extended * Ng_extended			// NgのDIMENSION乗の計算
#define D_MATRIX_SIZE 3										// 応力歪マトリックスの大きさ（2次元:3 3次元:6）
#define K_DIVISION_LENGE 10 	// 全体剛性マトリックスのcol&ptrを制作時に分ける節点数
#define EPS 1.0e-10				// 連立1次方程式の残差
#define N_STRAIN 4
#define N_STRESS 4
// 各種最大配置可能数
// #define MAX_N_KNOT 1000
// #define MAX_N_ELEMENT 12000
// #define MAX_N_NODE 110000
// #define MAX_N_LOAD 100000
// #define MAX_N_CONSTRAINT 100000
// #define MAX_K_WHOLE_SIZE MAX_N_NODE * DIMENSION
// #define MAX_NON_ZERO 10000000
// #define MAX_N_PATCH 100
// #define MAX_N_Controlpoint_in_Patch 12000
// #define MAX_N_ORDER	5
// #define MAX_N_DISTRIBUTE_FORCE 100

// for s-IGA
#define GAUSS_1DIR	Ng_extended						// 重なり判定のための一方向ガウス点数
#define NO_GAUSS_PT		GAUSS_1DIR * GAUSS_1DIR		// 重なり判定のためのガウス点総数
#define MAX_N_POINT_OVER	GAUSS_1DIR * GAUSS_1DIR	// 要素重なり判定に用いるローカルメッシュ上1要素内の点数
#define MAX_N_MESH  10								// 重合IGAを行うモデルの総数（ローカルメッシュ+1）
#define MAX_N_ELEMENT_OVER	1000					// グローバルメッシュ内の1要素に重なる最大要素数
#define MAX_N_ELEMENT_OVER_POINT	5				// ローカル要素内の1点に重なるグローバル要素
#define MAX_N_ELEMENT_OVER_ELEMENT	MAX_N_ELEMENT_OVER_POINT * MAX_N_POINT_OVER		// ローカルメッシュ内の1要素に重なる最大要素数

// 重ね合わせの結果
#define DBL_MAX          1.7976931348623158e+308 // max value
#define DIVISION_ELE_XI 10
#define DIVISION_ELE_ETA 10
// 最大値
#define MAX_PATCHES MAX_N_PATCH							//最大パッチ数
#define MAX_ORDER MAX_N_ORDER							//最大次数(p)
#define MAX_CNRL_P MAX_N_Controlpoint_in_Patch			//最大コントロールポイント数(n)
// 各パッチでの最大値
#define MAX_KNOTS (MAX_CNRL_P + MAX_ORDER + 1)			//ノットベクトルの最大長さ(n+p+1)
// 各パッチ、各方向での最大値
#define MAX_ELEMENTS 100								//最大要素数
#define MAX_DIVISION 10									//一要素あたりの最大分割数
#define MAX_POINTS (MAX_ELEMENTS * MAX_DIVISION + 1)	//最大点数

void Get_Input_1(int tm, int *Total_Knot_to_mesh,
				 int *Total_Patch_on_mesh, int *Total_Patch_to_mesh,
				 int *Total_Control_Point_on_mesh, int *Total_Control_Point_to_mesh,
				 int *Total_Constraint_to_mesh, int *Total_Load_to_mesh, int *Total_DistributeForce_to_mesh,
				 char **argv);
void Get_Input_2(int tm, int *Total_Knot_to_mesh,
				 int *Total_Patch_to_mesh, int *Total_Control_Point_to_mesh,
				 int *Total_Element_on_mesh, int *Total_Element_to_mesh,
				 int *Total_Constraint_to_mesh, int *Total_Load_to_mesh, int *Total_DistributeForce_to_mesh,
				 int *Order, int *No_knot, int *No_Control_point, int *No_Controlpoint_in_patch,
				 int *Patch_controlpoint, double *Position_Knots, int *No_Control_point_ON_ELEMENT,
				 double *Node_Coordinate, double *Control_Coord_x, double *Control_Coord_y, double *Control_Weight,
				 int *Constraint_Node_Dir, double *Value_of_Constraint,
				 int *Load_Node_Dir, double *Value_of_Load,
				 int *type_load_array, int *iPatch_array, int *iCoord_array,
				 double *val_Coord_array, double *Range_Coord_array, double *Coeff_Dist_Load_array,
                 int *Total_Control_Point_to_patch, int *Total_Knot_to_patch_dim,
				 char **argv);
void Make_INC(int tm, int *Total_Patch_on_mesh, int *Total_Patch_to_mesh,
			  int *Total_Element_on_mesh, int *Total_Element_to_mesh,
			  int *Total_Control_Point_on_mesh, int *Total_Control_Point_to_mesh,
			  int *Total_DistributeForce_to_mesh,
			  int *INC, int *Patch_controlpoint, int *Total_Control_Point_to_patch,
			  int *No_Control_point, int *Controlpoint_of_Element,
			  int *Element_patch, int *Element_mesh,
			  double *difference, int *Total_Knot_to_patch_dim,
			  int *Total_element_all_ID, int *ENC,
			  int *real_element_line, int *line_No_real_element,
			  int *line_No_Total_element, int *real_element,
			  int *real_El_No_on_mesh, int *real_Total_Element_on_mesh,
			  int *real_Total_Element_to_mesh, double *Equivalent_Nodal_Force,
			  int *type_load_array, int *iPatch_array, int *iCoord_array,
			  double *val_Coord_array, double *Range_Coord_array, double *Coeff_Dist_Load_array,
			  int *Order, int *No_knot);
void Make_gauss_array(int select_GP);
int Make_K_EL(int El_No, double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], double K_EL[MAX_KIEL_SIZE][MAX_KIEL_SIZE], double E, double nu, int DM);
int Make_coupled_K_EL(int El_No_loc, int El_No_glo,
					  double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION],
					  double XG[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION],
					  double K_EL[MAX_KIEL_SIZE][MAX_KIEL_SIZE],
					  double E, double nu, int DM);
int Make_Displacement_grad_glo(int El_No_loc, int El_No_glo,
					  		   double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION],
					  		   double XG[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION]);


// 全体剛性マトリックス
int Make_Index_Dof(int Total_Control_Point, int Total_Constraint, int Constraint_Node_Dir[MAX_N_CONSTRAINT][2]);
void Make_K_Whole_Ptr_Col(int Total_Element, int Total_Control_Point, int K_Whole_Size);
void Make_K_Whole_Val(double E, double nu, int Total_Element, int DM);
void Make_Displacement_grad_glo_check(int Total_Element);
// for s-IGA
void Check_coupled_Glo_Loc_element_for_end(double element_loc[DIMENSION], int mesh_n_over, int mesh_n_org);
void Check_coupled_Glo_Loc_element_for_Gauss(double element_loc[DIMENSION], int mesh_n_over, int mesh_n_org);
void Make_Loc_Glo();
//連立1次方程式
void Make_F_Vec(int Total_Load, int Load_Node_Dir[MAX_N_LOAD][2], double Value_of_Load[MAX_N_LOAD], int K_Whole_Size);
void Make_F_Vec_disp_const(int Mesh_No, int Total_Constraint, int Constraint_Node_Dir[MAX_N_CONSTRAINT][2], double Value_of_Constraint[MAX_N_CONSTRAINT], double E, double nu, int DM);
void mat_vec_crs(double vec_result[], double vec[], const int ndof);
double inner_product(int ndof, double vec1[], double vec2[]);
int check_conv_CG(int ndof, double alphak, double pp[], double eps, int itr);
void Diag_Scaling_CG_pre(int ndof, int flag_operation);
void CG_Solver(int ndof, int max_itr, double eps, int flag_ini_val);
//PCG solver
int RowCol_to_icount(int row, int col);
void PCG_Solver(int ndof, int max_itr, double eps);
void Make_M(double *M, int *M_Ptr, int *M_Col, int ndof);
void M_mat_vec_crs(double *M, int *M_Ptr, int *M_Col, double vec_result[], double vec[], const int ndof);
int M_check_conv_CG(int ndof, double alphak, double pp[], double eps, double *solution_vec);
void CG(int ndof, double *solution_vec, double *M, int *M_Ptr, int *M_Col, double *right_vec);
//各種値
void Make_Strain(int Total_Element);
void Make_Stress_2D(double E, double nu, int Total_Element, int DM);
void Make_Stress_2D_glo(double E, double nu, int Total_Element, int DM);
void Make_StrainEnergyDensity_2D();
void Make_Displacement_grad(int El_No);
void Make_StrainEnergyDensity_2D_overlay();
void Make_ReactionForce(int Total_Control_Point);
void Make_Parameter_z(int Total_Element, double E, double nu, int DM);
void Make_Parameter_z_overlay(int Total_Element, double E, double nu, int DM);
//分布荷重
void Force_dis(int Distriction_Force[DIMENSION][3], double Val_Distribute_Force[DIMENSION], double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], double Fe[DIMENSION]);
//NURBSの計算
void element_coordinate(int Total_Element);
void calculate_Controlpoint_using_NURBS(double element[DIMENSION], int Total_Element);
void calculate_extendmesh_using_NURBS(double element_emsh[DIMENSION]);
void Gausspoint_coordinate(int Total_Element);
int Jacobian(int El_No, double a[DIMENSION][DIMENSION], double Local_coord[DIMENSION], double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION]);
int Make_B_Matrix(int El_No, double B[D_MATRIX_SIZE][MAX_KIEL_SIZE], double Local_coord[DIMENSION], double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], double *J);

// J積分
void J_Integral_Input_Data(int Total_Control_Point, int *Location_Crack_Tip_Patch, double Location_Local_Coordinates[DIMENSION], double Virtual_Crack_Extension_Ct_Pt[MAX_N_NODE][DIMENSION], double *DeltaA);
double J_Integral_Computation(int Total_Control_Point, double Virtual_Crack_Extension_Ct_Pt[MAX_N_NODE][DIMENSION], double DeltaA);
double J_Integral_Computation_Interaction(int Total_Control_Point, double Location_Local_Coordinates[DIMENSION], double Virtual_Crack_Extension_Ct_Pt[MAX_N_NODE][DIMENSION], double DeltaA, double E, double nu, int DM);
// Shape Function
double Shape_func (int I_No, double Local_coord[DIMENSION],int El_No);
// Distributed Load
int SerchForElement(int mesh_n, int iPatch, int Total_Element, int iX, int iY);
void Setting_Dist_Load_2D(int mesh_n, int iPatch, int Total_Element, int iCoord, double val_Coord, double Range_Coord[2], int type_load, double Coeff_Dist_Load[3]);
void Add_Equivalent_Nodal_Forec_to_F_Vec(int Total_Control_Point);

//重ね合わせの結果
void GetLocData();
void ReadFile();
int CalcXiEtaByNR(double px, double py,
                  double *input_knot_vec_xi, double *input_knot_vec_eta,
                  double *cntl_px, double *cntl_py,
                  double *disp_cntl_px, double *disp_cntl_py,
                  int cntl_p_n_xi, int cntl_p_n_eta,
                  double *weight, int order_xi, int order_eta,
                  double *output_xi, double *output_eta,
				  double *disp_x_glo, double *disp_y_glo,
                  double *strain_xx_glo, double *strain_yy_glo, double *strain_xy_glo);
void Calculation(int order_xi, int order_eta,
                 int knot_n_xi, int knot_n_eta,
                 int cntl_p_n_xi, int cntl_p_n_eta,
                 double *input_knot_vec_xi, double *input_knot_vec_eta,
                 double *cntl_px, double *cntl_py,
                 double *disp_cntl_px, double *disp_cntl_py,
                 double *weight);
void Calculation_overlay(int order_xi_loc, int order_eta_loc,
                         int knot_n_xi_loc, int knot_n_eta_loc,
                         int cntl_p_n_xi_loc, int cntl_p_n_eta_loc,
                         double *knot_vec_xi_loc, double *knot_vec_eta_loc,
                         double *cntl_px_loc, double *cntl_py_loc,
                         double *weight_loc,
                         int order_xi_glo, int order_eta_glo,
                         int cntl_p_n_xi_glo, int cntl_p_n_eta_glo,
                         double *knot_vec_xi_glo, double *knot_vec_eta_glo,
                         double *cntl_px_glo, double *cntl_py_glo,
                         double *disp_cntl_px_glo, double *disp_cntl_py_glo,
                         double *weight_glo);
// void Calculation_overlay_at_GP(double E, double nu,
// 							   int order_xi_glo, int order_eta_glo,
// 							   int knot_n_xi_glo, int knot_n_eta_glo,
// 							   int cntl_p_n_xi_glo, int cntl_p_n_eta_glo,
// 							   double *knot_vec_xi_glo, double *knot_vec_eta_glo,
// 							   double *cntl_px_glo, double *cntl_py_glo,
// 							   double *disp_cntl_px_glo, double *disp_cntl_py_glo,
// 							   double *weight_glo);
// static void Calculation_at_GP(double E, double nu);

void K_output_svg(int ndof);

// 要素合成マトリックス
int Make_b_grad_Matrix(int El_No, double b_grad[DIMENSION * DIMENSION][2 * MAX_NO_CCpoint_ON_ELEMENT], double Local_coord[DIMENSION], double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], double *J);
int Make_D_Matrix_2D(double D[D_MATRIX_SIZE][D_MATRIX_SIZE], double E, double nu, int DM);

// Interaction integral
void Make_auxiliary_mode1(int e, double E, double nu, int DM, double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], double crack_front_coordinates_x, double crack_front_coordinates_y);
void Make_auxiliary_mode2(int e, double E, double nu, int DM, double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], double crack_front_coordinates_x, double crack_front_coordinates_y);
void Coordinate_Transform(double T[DIMENSION][DIMENSION], double S[DIMENSION][DIMENSION], double result[DIMENSION][DIMENSION]);

#endif