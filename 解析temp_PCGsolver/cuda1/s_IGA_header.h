#ifndef S_IGA_HEADER_H
#define S_IGA_HEADER_H

#define SKIP_S_IGA 2 // 重ね合わせとJ積分を行う 0, 重ね合わせをスキップしてJ積分を行う 1, J積分を行わない 2
#define DM 1         // 平面応力状態:DM=0	平面ひずみ状態:DM=1
#define ERROR -999
#define PI  3.14159265359
#define MAX_NO_CCpoint_ON_ELEMENT 16						// 分割節点数
#define DIMENSION 2											// 次元数
#define MAX_ORDER 4											// 基底関数の次数の最大値 + 1
#define MAX_KIEL_SIZE MAX_NO_CCpoint_ON_ELEMENT * DIMENSION	// 要素分割マトリックスの大きさ
#define Ng 4												// Gauss-Legendreの足す回数
#define POW_Ng Ng * Ng										// NgのDIMENSION乗の計算
#define Ng_extended 10										// Gauss-Legendreの足す回数
#define POW_Ng_extended Ng_extended * Ng_extended			// NgのDIMENSION乗の計算
#define K_DIVISION_LENGE 10 	// 全体剛性マトリックスのcol&ptrを制作時に分ける節点数
#define EPS 1.0e-10				// 連立1次方程式の残差
#define N_STRAIN 4
#define N_STRESS 4
#define MAX_N_ELEMENT_OVER	100  					// グローバルメッシュ内の1要素に重なる最大要素数
#define MAX_N_ELEMENT_OVER_POINT	5				// ローカル要素内の1点に重なるグローバル要素
#define MAX_N_ELEMENT_OVER_ELEMENT	MAX_N_ELEMENT_OVER_POINT * POW_Ng_extended		// ローカルメッシュ内の1要素に重なる最大要素数
// #define D_MATRIX_SIZE 3										// 応力歪マトリックスの大きさ（2次元:3 3次元:6）
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
// #define GAUSS_1DIR	Ng_extended						// 重なり判定のための一方向ガウス点数
// #define NO_GAUSS_PT		GAUSS_1DIR * GAUSS_1DIR		// 重なり判定のためのガウス点総数
// #define MAX_N_POINT_OVER	GAUSS_1DIR * GAUSS_1DIR	// 要素重なり判定に用いるローカルメッシュ上1要素内の点数
// #define MAX_N_MESH  10								// 重合IGAを行うモデルの総数（ローカルメッシュ+1）

// 重ね合わせの結果
#define DBL_MAX          1.7976931348623158e+308 // max value
#define DIVISION_ELE_XI 10
#define DIVISION_ELE_ETA 10
#define MAX_PATCHES MAX_N_PATCH							//最大パッチ数
// #define MAX_ORDER MAX_ORDER							//最大次数(p)
// #define MAX_CNRL_P MAX_N_Controlpoint_in_Patch			//最大コントロールポイント数(n)
// #define MAX_KNOTS (MAX_CNRL_P + MAX_ORDER + 1)			//ノットベクトルの最大長さ(n+p+1)
#define MAX_ELEMENTS 100								//最大要素数
#define MAX_DIVISION 10									//一要素あたりの最大分割数
#define MAX_POINTS (MAX_ELEMENTS * MAX_DIVISION + 1)	//最大点数

// Get input file data
void Get_Input_1(int tm, int *Total_Knot_to_mesh,
				 int *Total_Patch_on_mesh, int *Total_Patch_to_mesh,
				 int *Total_Control_Point_on_mesh, int *Total_Control_Point_to_mesh,
				 int *Total_Constraint_to_mesh, int *Total_Load_to_mesh, int *Total_DistributeForce_to_mesh,
				 char **argv);
void Get_Input_2(int tm, int *Total_Knot_to_mesh,
				 int *Total_Patch_to_mesh, int *Total_Control_Point_to_mesh,
				 int *Total_Element_on_mesh, int *Total_Element_to_mesh,
				 int *Total_Constraint_to_mesh, int *Total_Load_to_mesh, int *Total_DistributeForce_to_mesh,
				 int *Order, int *No_knot, int *No_Control_point, int *No_Control_point_in_patch,
				 int *Patch_Control_point, double *Position_Knots, int *No_Control_point_ON_ELEMENT,
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
			  int *INC, int *Patch_Control_point, int *Total_Control_Point_to_patch,
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
			  int *Order, int *No_knot, int *Total_Knot_to_mesh, double *Node_Coordinate,
			  double *Position_Knots, int *No_Control_point_ON_ELEMENT);
// Distributed Load
void Setting_Dist_Load_2D(int mesh_n, int iPatch, int Total_Element, int iCoord, double val_Coord,
						  double Range_Coord[2], int type_load, double Coeff_Dist_Load[3], int *Total_Knot_to_mesh,
						  int *Controlpoint_of_Element, int *Order, int *No_knot, int *Total_element_all_ID,
						  int *Total_Knot_to_patch_dim, double *Position_Knots, double *Equivalent_Nodal_Force,
						  int *Total_Element_on_mesh, int *Total_Element_to_mesh, int *Element_patch, int *ENC,
						  double *Node_Coordinate, int *INC, int *No_Control_point_ON_ELEMENT, int *Total_Control_Point_to_mesh,
						  int *No_Control_point);
int SearchForElement(int mesh_n, int iPatch, int iX, int iY, int *Total_Element_on_mesh, int *Total_Element_to_mesh, int *Element_patch, int *ENC);
// for s_IGA
void Check_coupled_Glo_Loc_element_for_Gauss(int mesh_n_over, int mesh_n_org, int *NNLOVER, int *NELOVER,
											 double *Gauss_Coordinate, double *Gauss_Coordinate_ex, double *Jac, double *Jac_ex, double *B_Matrix, double *B_Matrix_ex, double *Loc_parameter_on_Glo, double *Loc_parameter_on_Glo_ex,
											 int *real_Total_Element_to_mesh, double *Node_Coordinate, int *Total_Control_Point_to_mesh, int *Controlpoint_of_Element,
											 int *INC, int *Element_patch, int *Order, int *No_Control_point_ON_ELEMENT, double *Position_Knots,
											 int *Total_Knot_to_patch_dim, int *No_knot, int *No_Control_point,
											 double *Control_Coord_x, double *Control_Coord_y, double *Control_Weight,
											 int *real_Total_Element_on_mesh, int *real_element, int *Total_Patch_on_mesh, int *line_No_Total_element);
void Make_Loc_Glo(int *real_Total_Element_on_mesh, int *real_Total_Element_to_mesh, int *real_element, int *NNLOVER, int *NELOVER);
int ele_check(int patch_n, double para_coord[DIMENSION], int *No_Control_point, double *Position_Knots, int *Total_Knot_to_patch_dim, int *No_knot, int *Order, int *temp_element_n, int *line_No_Total_element);
void sort(int total, int *element_n_point);
int duplicate_delete(int total, int element_n, int *NELOVER, int *element_n_point);
// Preprocessing
void Preprocessing(int m, int e, double *Gauss_Coordinate, double *Gauss_Coordinate_ex, double *B_Matrix, double *B_Matrix_ex,
				   double *Jac, double *Jac_ex, int *real_Total_Element_to_mesh, double *Node_Coordinate, int *Total_Control_Point_to_mesh,
				   int *Controlpoint_of_Element, int *INC, int *Element_patch, int *Order, int *No_Control_point_ON_ELEMENT,
				   double *Position_Knots, int *Total_Knot_to_patch_dim, int *No_knot, int *No_Control_point);
void Make_Gauss_Coordinate(int m, int e, double *Gauss_Coordinate, double *Gauss_Coordinate_ex, double *Node_Coordinate, int *Total_Control_Point_to_mesh,
						   int *Controlpoint_of_Element, int *INC, int *Element_patch, int *Order, int *No_Control_point_ON_ELEMENT,
						   double *Position_Knots, int *Total_Knot_to_patch_dim, int *No_knot, int *No_Control_point);
void Make_dSF(int m, int e, double *dSF, double *dSF_ex, int *Element_patch, int *No_Control_point_ON_ELEMENT,
			  int *Controlpoint_of_Element, int *Total_Control_Point_to_mesh, double *Node_Coordinate, int *INC, int *Order,
			  double *Position_Knots, int *Total_Knot_to_patch_dim, int *No_knot, int *No_Control_point);
void Make_Jac(int m, int e, double *Jac, double *Jac_ex, double *dSF, double *dSF_ex, double *a_matrix, double *Node_Coordinate, int *Controlpoint_of_Element, int *No_Control_point_ON_ELEMENT, int *Element_patch);
void Make_B_Matrix(int m, int e, double *B_Matrix, double *B_Matrix_ex, double *dSF, double *dSF_ex, double *a_matrix, int *No_Control_point_ON_ELEMENT, int *Element_patch);
void Make_gauss_array(int select_GP);
// K matrix
void Make_D_Matrix(double *D);
void Make_Index_Dof(int *Total_Control_Point_to_mesh, int *Total_Constraint_to_mesh, int *Constraint_Node_Dir, int *Index_Dof);
void Make_K_Whole_Ptr_Col(int *K_Whole_Ptr, int *Total_Element_to_mesh, int *Total_Control_Point_to_mesh, int *Total_Control_Point_To_Node,
						  int *No_Control_point_ON_ELEMENT, int *Element_patch, int *Controlpoint_of_Element, int *Node_To_Node, int *NNLOVER,
						  int *NELOVER, int *Index_Dof, int *K_Whole_Ptr, int *K_Whole_Col);
void Make_K_Whole_Val(double E, double nu, int Total_Element, int DM);
void Make_K_EL(int El_No, double *K_EL, int *No_Control_point_ON_ELEMENT, int *Element_patch,
			   double *Jac, double *Jac_ex, double *B_Matrix, double *B_Matrix_ex, double *D);
void Make_BG_Matrix(int El_No, double B[D_MATRIX_SIZE][MAX_KIEL_SIZE], double Local_coord[DIMENSION], double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], double *J);
void Make_coupled_K_EL();
void BDBJ(int *KIEL_SIZE, double *B, double *D, double J, double *K_EL)
void coupled_BDBJ();
// F vector
void Make_F_Vec(int Total_Load, int Load_Node_Dir[MAX_N_LOAD][2], double Value_of_Load[MAX_N_LOAD], int K_Whole_Size);
void Make_F_Vec_disp_const(int Mesh_No, int Total_Constraint, int Constraint_Node_Dir[MAX_N_CONSTRAINT][2], double Value_of_Constraint[MAX_N_CONSTRAINT], double E, double nu, int DM);
void Add_Equivalent_Nodal_Forec_to_F_Vec(int Total_Control_Point);
// tool
double InverseMatrix_2x2(double M[DIMENSION][DIMENSION]);
double InverseMatrix_3x3(double M[DIMENSION][DIMENSION]);
// Shape Function
double Shape_func(int I_No, double Local_coord[DIMENSION], int El_No, double *Node_Coordinate, int *Total_Control_Point_to_mesh,
				  int *Controlpoint_of_Element, int *INC, int *Element_patch, int *Order, int *No_Control_point_ON_ELEMENT,
				  double *Position_Knots, int *Total_Knot_to_patch_dim, int *No_knot, int *No_Control_point);
void ShapeFunc_from_paren(double Position_Data_param[DIMENSION], double Local_coord[DIMENSION], int j, int e, int *INC,
						  double *Position_Knots, int *Total_Knot_to_patch_dim, int *Controlpoint_of_Element, int *Element_patch);
void ShapeFunction1D(double Position_Data_param[DIMENSION], int j, int e, double *Shape, double *dShape, int *No_knot, int *Total_Control_Point_to_mesh,
					 double *Position_Knots, int *Total_Knot_to_patch_dim,int *Element_patch, int *Order, int *No_Control_point);
double dShape_func(int I_No, int xez, double Local_coord[DIMENSION], int El_No, int *Controlpoint_of_Element, int *No_Control_point_ON_ELEMENT,
				   int *Total_Control_Point_to_mesh, int *Element_patch, double *Node_Coordinate, int *INC, int *Order, double *Position_Knots,
				   int *Total_Knot_to_patch_dim, int *No_knot, int *No_Control_point);
void NURBS_deriv(double Local_coord[DIMENSION], int El_No, double *Node_Coordinate, int *Total_Control_Point_to_mesh,
				 int *Controlpoint_of_Element, int *INC, int *Element_patch, int *Order, int *No_Control_point_ON_ELEMENT,
				 double *Position_Knots, int *Total_Knot_to_patch_dim, int *No_knot, int *No_Control_point,
				 double *dShape_func1, double *dShape_func2);
double dShapeFunc_from_paren(int j, int e, int *INC, int *Controlpoint_of_Element, double *Position_Knots, int *Total_Knot_to_patch_dim, int *Element_patch);
// CG solver
void mat_vec_crs(double *vec_result, double *vec, const int ndof);
double inner_product(int ndof, double *vec1, double *vec2);
int check_conv_CG(int ndof, double alphak, double *pp, double eps, int itr);
void Diag_Scaling_CG_pre(int ndof, int flag_operation);
void CG_Solver(int ndof, int max_itr, double eps, int flag_ini_val);
// PCG solver
int RowCol_to_icount(int row, int col);
void PCG_Solver(int ndof, int max_itr, double eps);
void Make_M(double *M, int *M_Ptr, int *M_Col, int ndof);
void M_mat_vec_crs(double *M, int *M_Ptr, int *M_Col, double *vec_result, double *vec, const int ndof);
int M_check_conv_CG(int ndof, double alphak, double *pp, double eps, double *solution_vec);
void CG(int ndof, double *solution_vec, double *M, int *M_Ptr, int *M_Col, double *right_vec);
// Newton-Raphson
double BasisFunc(double *knot_vec, int knot_index, int order, double xi,
				 double *output, double *d_output);
double rBasisFunc(double *knot_vec, int knot_index,
				  int order, double xi,
				  double *output, double *d_output);
double lBasisFunc(double *knot_vec, int knot_index,
				  int cntl_p_n, int order, double xi,
				  double *output, double *d_output);
double NURBS_surface(double *input_knot_vec_xi, double *input_knot_vec_eta,
					 double *cntl_px, double *cntl_py,
					 int cntl_p_n_xi, int cntl_p_n_eta,
					 double *weight, int order_xi, int order_eta,
					 double xi, double eta,
					 double *output_x, double *output_y,
					 double *output_dxi_x, double *output_deta_x,
					 double *output_dxi_y, double *output_deta_y);
double rNURBS_surface(double *input_knot_vec_xi, double *input_knot_vec_eta,
					  double *cntl_px, double *cntl_py,
					  int cntl_p_n_xi, int cntl_p_n_eta,
					  double *weight, int order_xi, int order_eta,
					  double xi, double eta,
					  double *output_x, double *output_y,
					  double *output_dxi_x, double *output_deta_x,
					  double *output_dxi_y, double *output_deta_y);
double lNURBS_surface(double *input_knot_vec_xi, double *input_knot_vec_eta,
					  double *cntl_px, double *cntl_py,
					  int cntl_p_n_xi, int cntl_p_n_eta,
					  double *weight, int order_xi, int order_eta,
					  double xi, double eta,
					  double *output_x, double *output_y,
					  double *output_dxi_x, double *output_deta_x,
					  double *output_dxi_y, double *output_deta_y);
double rlNURBS_surface(double *input_knot_vec_xi, double *input_knot_vec_eta,
					   double *cntl_px, double *cntl_py,
					   int cntl_p_n_xi, int cntl_p_n_eta,
					   double *weight, int order_xi, int order_eta,
					   double xi, double eta,
					   double *output_x, double *output_y,
					   double *output_dxi_x, double *output_deta_x,
					   double *output_dxi_y, double *output_deta_y);
double lrNURBS_surface(double *input_knot_vec_xi, double *input_knot_vec_eta,
					   double *cntl_px, double *cntl_py,
					   int cntl_p_n_xi, int cntl_p_n_eta,
					   double *weight, int order_xi, int order_eta,
					   double xi, double eta,
					   double *output_x, double *output_y,
					   double *output_dxi_x, double *output_deta_x,
					   double *output_dxi_y, double *output_deta_y);
int Calc_xi_eta(double px, double py,
				double *input_knot_vec_xi, double *input_knot_vec_eta,
				int cntl_p_n_xi, int cntl_p_n_eta, int order_xi, int order_eta,
				double *output_xi, double *output_eta,
				double *Position_Knots, int *Total_Knot_to_patch_dim, int *No_Control_point, int *Order,
				double *Control_Coord_x, double *Control_Coord_y, double *Control_Weight, int *No_knot);
// for s_IGA overlay
void s_IGA_overlay();
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
// output
void output_file();
void K_output_svg(int ndof);













// J積分


// 仮











//各種値
void Make_Strain(int Total_Element);
void Make_Stress_2D(double E, double nu, int Total_Element, int DM);
void Make_ReactionForce(int Total_Control_Point);
void Make_Parameter_z(int Total_Element, double E, double nu, int DM);
void Make_Parameter_z_overlay(int Total_Element, double E, double nu, int DM);

// J積分


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

#endif