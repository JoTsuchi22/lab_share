#ifndef S_IGA_HEADER_H
#define S_IGA_HEADER_H

#define SKIP_S_IGA 2 // 重ね合わせとJ積分を行う 0, 重ね合わせをスキップしてJ積分を行う 1, J積分を行わない 2
#define OUTPUT_SVG 1 // SVG 出力を行わない 0, 行う 1
#define DM 1         // 平面応力状態:DM=0	平面ひずみ状態:DM=1
#define ERROR -999
#define MAX_NO_CCpoint_ON_ELEMENT 16						// 分割節点数
#define DIMENSION 2											// 次元数
#define MAX_ORDER 4											// 基底関数の次数の最大値 + 1
#define MAX_KIEL_SIZE MAX_NO_CCpoint_ON_ELEMENT * DIMENSION	// 要素分割マトリックスの大きさ
#define NG 4												// Gauss-Legendreの積分点数
#define POW_NG NG * NG										// NGのDIMENSION乗の計算
#define NG_EXTEND 10										// Gauss-Legendreの積分点数
#define POW_NG_EXTEND NG_EXTEND * NG_EXTEND			        // NGのDIMENSION乗の計算
#define K_DIVISION_LENGE 10 	                            // 全体剛性マトリックスのcol&ptrを制作時に分ける節点数
#define EPS 1.0e-10				                            // 連立1次方程式の残差
#define N_STRAIN 4
#define N_STRESS 4
#define MAX_N_ELEMENT_OVER 100  				            // グローバルメッシュ内の1要素に重なる最大要素数
#define MAX_N_ELEMENT_OVER_POINT 5				            // ローカル要素内の1点に重なるグローバル要素
#define MAX_N_ELEMENT_OVER_ELEMENT MAX_N_ELEMENT_OVER_POINT * POW_NG_EXTEND // ローカルメッシュ内の1要素に重なる最大要素数
#define DIVISION_ELE 10                                     // 一要素あたりの分割数
#define DBL_MAX 1.7976931348623158e+308                     // max value

// stract array ptr def
struct information {
    int *Total_Knot_to_mesh;
    int *Total_Patch_on_mesh;
    int *Total_Patch_to_mesh;
    int *Total_Control_Point_on_mesh;
    int *Total_Control_Point_to_mesh;
    int *Total_Element_on_mesh;
    int *Total_Element_to_mesh;
    int *real_Total_Element_on_mesh;
    int *real_Total_Element_to_mesh;
    int *Total_Load_to_mesh;
    int *Total_Constraint_to_mesh;
    int *Total_DistributeForce_to_mesh;

    int *Order;
    int *No_knot;
    int *Total_Control_Point_to_patch;
    int *Total_Knot_to_patch_dim;
    double *Position_Knots;
    int *No_Control_point;
    int *No_Control_point_in_patch;
    int *Patch_Control_point;
    int *No_Control_point_ON_ELEMENT;
    double *Node_Coordinate;
    double *Control_Coord_x;
    double *Control_Coord_y;
    double *Control_Weight;
    int *Constraint_Node_Dir;
    double *Value_of_Constraint;
    int *Load_Node_Dir;
    double *Value_of_Load;
    int *iPatch_array;
    int *iCoord_array;
    int *type_load_array;
    double *val_Coord_array;
    double *Range_Coord_array;
    double *Coeff_Dist_Load_array;

    int *INC;
    int *Controlpoint_of_Element;
    int *Element_patch;
    int *Element_mesh;
    int *line_No_real_element;
    int *line_No_Total_element;
    double *difference;
    int *Total_element_all_ID;
    int *ENC;
    int *real_element_line;
    int *real_element;
    int *real_El_No_on_mesh;
    double *Equivalent_Nodal_Force;

    int *NNLOVER;
    int *NELOVER;
    double *Gauss_Coordinate;
    double *Gauss_Coordinate_ex;
    double *Jac;
    double *Jac_ex;
    double *B_Matrix;
    double *B_Matrix_ex;
    double *Loc_parameter_on_Glo;
    double *Loc_parameter_on_Glo_ex;

    double *D;
    int *Node_To_Node;
    int *Total_Control_Point_To_Node;
    int *Index_Dof;
    int *K_Whole_Ptr;
    int *K_Whole_Col;
    double *K_Whole_Val;

    double *sol_vec;
    double *rhs_vec;

    double *Displacement;

	double *coord_x;
	double *coord_y;
	double *disp_x;
	double *disp_y;
	double *strain;
	double *stress;
};

// Get input file data
void Get_Input_1(int tm, char **argv, information *info);
void Get_Input_2(int tm, char **argv, information *info);
void Make_INC(int tm, information *info);
// Distributed Load
void Setting_Dist_Load_2D(int mesh_n, int iPatch, int iCoord, double val_Coord, double *Range_Coord, int type_load, double *Coeff_Dist_Load, information *info);
int SearchForElement(int mesh_n, int iPatch, int iX, int iY, information *info);
// for s_IGA
void Check_coupled_Glo_Loc_element_for_Gauss(int mesh_n_over, int mesh_n_org, information *info);
void Make_Loc_Glo(information *info);
int ele_check(int patch_n, double *para_coord, int *temp_element_n, information *info);
void sort(int total, int *element_n_point);
int duplicate_delete(int total, int element_n, int *element_n_point, information *info);
// Preprocessing
void Preprocessing(int m, int e, information *info);
void Make_Gauss_Coordinate(int m, int e, information *info);
void Make_dSF(int m, int e, double *dSF, double *dSF_ex, information *info);
void Make_Jac(int m, int e, double *dSF, double *dSF_ex, double *a_matrix, information *info);
void Make_B_Matrix(int m, int e, double *dSF, double *dSF_ex, double *a_matrix, information *info);
void Make_gauss_array(int select_GP);
// K matrix
void Make_D_Matrix(information *info);
void Make_Index_Dof(information *info);
void Make_K_Whole_Ptr_Col(information *info);
void Make_K_Whole_Val(information *info);
void Make_K_EL(int El_No, double *K_EL, information *info);
void Make_coupled_K_EL(int El_No_loc, int El_No_glo, double *coupled_K_EL, information *info);
void Make_BG_Matrix(int El_No, double *BG, double *Local_coord, information *info);
void BDBJ(int KIEL_SIZE, double *B, double J, double *K_EL, information *info);
void coupled_BDBJ(int KIEL_SIZE, double *B, double *BG, double J, double *K_EL, information *info);
// F vector
void Make_F_Vec(information *info);
void Make_F_Vec_disp_const(information *info);
void Add_Equivalent_Nodal_Force_to_F_Vec(information *info);
// PCG solver
void PCG_Solver(int max_itetarion, double eps, information *info);
void Make_M(double *M, int *M_Ptr, int *M_Col, int ndof, information *info);
void CG(int ndof, double *solution_vec, double *M, int *M_Ptr, int *M_Col, double *right_vec, double *gg, double *dd, double *pp);
void M_mat_vec_crs(double *M, int *M_Ptr, int *M_Col, double *vec_result, double *vec, const int ndof);
double inner_product(int ndof, double *vec1, double *vec2);
int M_check_conv_CG(int ndof, double alphak, double *pp, double eps, double *solution_vec);
int RowCol_to_icount(int row, int col, information *info);
// tool
double InverseMatrix_2x2(double M[DIMENSION][DIMENSION]);
double InverseMatrix_3x3(double M[DIMENSION][DIMENSION]);
// Shape Function
double Shape_func(int I_No, double *Local_coord, int El_No, information *info);
void ShapeFunc_from_paren(double *Position_Data_param, double *Local_coord, int j, int e, information *info);
void ShapeFunction1D(double *Position_Data_param, int j, int e, double *Shape, double *dShape, information *info);
double dShape_func(int I_No, int xez, double *Local_coord, int El_No, information *info);
void NURBS_deriv(double *Local_coord, int El_No, double *dShape_func1, double *dShape_func2, information *info);
double dShapeFunc_from_paren(int j, int e, information *info);
// Newton-Raphson
double BasisFunc(double *knot_vec, int knot_index, int order, double xi, double *output, double *d_output);
double rBasisFunc(double *knot_vec, int knot_index, int order, double xi, double *output, double *d_output);
double lBasisFunc(double *knot_vec, int knot_index, int cntl_p_n, int order, double xi, double *output, double *d_output);
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
				double *output_xi, double *output_eta, information *info);
// Postprocessing
void Make_Displacement(information *info);
// for s_IGA overlay
void s_IGA_overlay(information *info);
void Calculation(int order_xi, int order_eta, int knot_n_xi, int knot_n_eta, int cntl_p_n_xi, int cntl_p_n_eta,
				 double *input_knot_vec_xi, double *input_knot_vec_eta, double *cntl_px, double *cntl_py, double *weight,
				 double *disp_cntl_px, double *disp_cntl_py, int patch_num, information *info);
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
						 double *weight_glo, int patch_num, information *info);
int CalcXiEtaByNR(double px, double py,
				  double *input_knot_vec_xi, double *input_knot_vec_eta,
				  double *cntl_px, double *cntl_py,
				  double *disp_cntl_px, double *disp_cntl_py,
				  int cntl_p_n_xi, int cntl_p_n_eta,
				  double *weight, int order_xi, int order_eta,
				  double *output_xi, double *output_eta,
				  double *disp_x_glo, double *disp_y_glo,
				  double *strain_xx_glo, double *strain_yy_glo, double *strain_xy_glo);
// output
void output_for_viewer(information *info);
void K_output_svg(information *info);




// J integral



// 各種値
// void Make_Strain(int Total_Element);
// void Make_Stress_2D(double E, double nu, int Total_Element, int DM);
// void Make_ReactionForce(int Total_Control_Point);
// void Make_Parameter_z(int Total_Element, double E, double nu, int DM);
// void Make_Parameter_z_overlay(int Total_Element, double E, double nu, int DM);


#endif