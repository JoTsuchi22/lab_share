/**************************************
s-IGA
multipch

Newton_Raphson法
(要素の重なり判定については2通りの判定方法
check_over_parameter=0:要素の端点
or
check_over_parameter=1:ガウス点
結局ガウス点で判定するcheck_over_parameter=1が正しい
--->明らかに要素の端点のみによる要素判定で判定可能であり、
	計算時間を少しでも短縮したいときはcheck_over_parameter=0にする)

NURBSviewerのための入力データ出力(input_for_NURBS.txt)
s_NURBS_viewerのための入力データ出力(input_local.txt)
ガウス積分の積分点数：4(一部10)

仮定条件
	ローカルメッシュ同士は原則被りなしと仮定
	グローバルパッチはシングルパッチ

Input
inputファイルを複数読み込む
2つ目以降がローカルメッシュ

Output
mkdir checkAns
	  colored_point
	  Gauss_stress
	  mesh_net
	  new_zarusoba
	  NURBS
	  shapefunc

**************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#define SKIP_S_IGA 2 // 重ね合わせとJ積分を行う 0, 重ね合わせをスキップしてJ積分を行う 1, J積分を行わない 2

#define ERROR -999
#define PI  3.14159265359

#define MAX_NO_CCpoint_ON_ELEMENT 16						//分割節点数
#define DIMENSION 2											//次元数
#define MAX_KIEL_SIZE MAX_NO_CCpoint_ON_ELEMENT * DIMENSION	//要素分割マトリックスの大きさ
#define Ng 4												//Gauss-Legendreの足す回数
#define POW_Ng Ng * Ng										//NgのDIMENSION乗の計算
#define Ng_extended 10										//Gauss-Legendreの足す回数
#define POW_Ng_extended Ng_extended * Ng_extended			//NgのDIMENSION乗の計算
#define D_MATRIX_SIZE 3										//応力歪マトリックスの大きさ（2次元:3 3次元:6）

#define K_DIVISION_LENGE 10 	//全体剛性マトリックスのcol&ptrを制作時に分ける節点数
#define EPS 1.0e-10				//連立1次方程式の残差
#define N_STRAIN 4
#define N_STRESS 4
//各種最大配置可能数
#define MAX_N_KNOT 1000
#define MAX_N_ELEMENT 14000
// #define MAX_N_ELEMENT 1000
#define MAX_N_NODE 110000
#define MAX_N_LOAD 100000
#define MAX_N_CONSTRAINT 100000
#define MAX_K_WHOLE_SIZE MAX_N_NODE * DIMENSION
#define MAX_NON_ZERO 10000000
#define MAX_N_PATCH 10
#define MAX_N_Controlpoint_in_Patch 12000
#define MAX_N_ORDER	5

#define MAX_N_DISTRIBUTE_FORCE 100
// #define DISTRIBUTE_FORCE_Ng 3

//for s-IGA
#define GAUSS_1DIR	Ng_extended						//重なり判定のための一方向ガウス点数
#define NO_GAUSS_PT		GAUSS_1DIR * GAUSS_1DIR		//重なり判定のためのガウス点総数
#define MAX_N_POINT_OVER	GAUSS_1DIR * GAUSS_1DIR	//要素重なり判定に用いるローカルメッシュ上1要素内の点数
#define MAX_N_MESH  10								//重合IGAを行うモデルの総数（ローカルメッシュ+1）
#define MAX_N_ELEMENT_OVER	1000					//グローバルメッシュ内の1要素に重なる最大要素数
#define MAX_N_ELEMENT_OVER_POINT	5				//ローカル要素内の1点に重なるグローバル要素
#define MAX_N_ELEMENT_OVER_ELEMENT	MAX_N_ELEMENT_OVER_POINT * MAX_N_POINT_OVER		//ローカルメッシュ内の1要素に重なる最大要素数


//重ね合わせの結果

#define DBL_MAX          1.7976931348623158e+308 // max value

#define DIVISION_ELE_XI 10
#define DIVISION_ELE_ETA 10

/////変更後はmake cleanしてからmakeする/////////////////////////////////////////////////////
///////最大値////////////////////////////////////////////////////////////////////////////////
#define MAX_PATCHES MAX_N_PATCH							//最大パッチ数
// #define MAX_PATCHES 15								//最大パッチ数
#define MAX_ORDER MAX_N_ORDER							//最大次数(p)
// #define MAX_ORDER 3									//最大次数(p)
#define MAX_CNRL_P MAX_N_Controlpoint_in_Patch			//最大コントロールポイント数(n)
// #define MAX_CNRL_P 2500								//最大コントロールポイント数(n)
///////各パッチでの最大値/////////////////////////////////////////////////////////////////////
#define MAX_KNOTS (MAX_CNRL_P + MAX_ORDER + 1)			//ノットベクトルの最大長さ(n+p+1)
// #define MAX_KNOTS MAX_N_KNOT							//ノットベクトルの最大長さ(n+p+1)
///////各パッチ、各方向での最大値//////////////////////////////////////////////////////////////
// #define MAX_ELEMENTS MAX_N_ELEMENT					//最大要素数
#define MAX_ELEMENTS 100								//最大要素数
#define MAX_DIVISION 10									//一要素あたりの最大分割数
#define MAX_POINTS (MAX_ELEMENTS * MAX_DIVISION + 1)	//最大点数
/////////////////////////////////////////////////////////////////////////////////////////////

//void Force_Dis( int Total_DistributeForce, int DistributeForce[MAX_N_DISTRIBUTE_FORCE][3], double Val_DistributeForce[MAX_N_DISTRIBUTE_FORCE],int *Total_Load,int Load_Node_Dir[MAX_N_LOAD][2],double Value_of_Load[MAX_N_LOAD],int Total_Control_Point, int El_No, int *Total_Element );
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
void Get_InputData(int tm,
                   double *E, double *nu, int *Total_Element, int *Total_Control_Point,
				   int *Total_Load, int *No_Patch, int Load_Node_Dir[MAX_N_LOAD][2], double Value_of_Load[MAX_N_LOAD],
				   int *Total_Constraint, int Constraint_Node_Dir[MAX_N_CONSTRAINT][2], double Value_of_Constraint[MAX_N_CONSTRAINT],
				   int *Total_DistributeForce, char *argv[]);

//全体剛性マトリックス
int Make_Index_Dof(int Total_Control_Point, int Total_Constraint, int Constraint_Node_Dir[MAX_N_CONSTRAINT][2]);
void Make_K_Whole_Ptr_Col(int Total_Element, int Total_Control_Point, int K_Whole_Size);
void Make_K_Whole_Val(double E, double nu, int Total_Element, int DM);
void Make_Displacement_grad_glo_check(int Total_Element);
//for s-IGA
void Check_coupled_Glo_Loc_element_for_end(double element_loc[DIMENSION], int mesh_n_over, int mesh_n_org);
void Check_coupled_Glo_Loc_element_for_Gauss(double element_loc[DIMENSION], int mesh_n_over, int mesh_n_org);
void Make_Loc_Glo();
//void Make_K_Whole_overlaid(int K_Whole_overlaid_Size);
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
// void Make_Displacement_grad_overlay(double E, double nu, int Total_Element , int El_No, int Total_Control_Point);
void Make_ReactionForce(int Total_Control_Point);
void Make_Parameter_z(int Total_Element, double E, double nu, int DM);
void Make_Parameter_z_overlay(int Total_Element, double E, double nu, int DM);
//分布荷重
void Force_dis(int Distriction_Force[DIMENSION][3], double Val_Distribute_Force[DIMENSION], double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], double Fe[DIMENSION]);
//void Make_Output( int Total_Control_Point, int Total_Element );
//NURBSの計算
void element_coordinate(int Total_Element);
void calculate_Controlpoint_using_NURBS(double element[DIMENSION], int Total_Element);
void calculate_extendmesh_using_NURBS(double element_emsh[DIMENSION]);
void Gausspoint_coordinate(int Total_Element);
int Jacobian(int El_No, double a[DIMENSION][DIMENSION], double Local_coord[DIMENSION], double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION]);
int Make_B_Matrix(int El_No, double B[D_MATRIX_SIZE][MAX_KIEL_SIZE], double Local_coord[DIMENSION], double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], double *J);

/*///J積分
int Make_B_x_Matrix_Quad_4(double B_x[DIMENSION][KIEL_SIZE], double Local_coord[DIMENSION], double X[No_Control_point_ON_ELEMENT][DIMENSION], double *J );
void Make_Strain_x_Quad_4(double E, double nu, int Total_Element);
void Make_EMT(double E, double nu, int Total_Element);*/
void J_Integral_Input_Data(int Total_Control_Point, int *Location_Crack_Tip_Patch, double Location_Local_Coordinates[DIMENSION], double Virtual_Crack_Extension_Ct_Pt[MAX_N_NODE][DIMENSION], double *DeltaA);
double J_Integral_Computation(int Total_Control_Point, double Virtual_Crack_Extension_Ct_Pt[MAX_N_NODE][DIMENSION], double DeltaA);
double J_Integral_Computation_Interaction(int Total_Control_Point, double Location_Local_Coordinates[DIMENSION], double Virtual_Crack_Extension_Ct_Pt[MAX_N_NODE][DIMENSION], double DeltaA, double E, double nu, int DM);
// double J_Integral_Computation6by6(int Total_Control_Point, int Total_Element,int Location_Crack_Tip_Patch, double Location_Local_Coordinates[DIMENSION],double Virtual_Crack_Extension_Ct_Pt[MAX_N_NODE][DIMENSION], double DeltaA, int El_No,  double E, double nu, int DM1);

/* Shape Function */
double Shape_func (int I_No, double Local_coord[DIMENSION],int El_No);

/* Distributed Load */

int SerchForElement(int mesh_n, int iPatch, int Total_Element, int iX, int iY);

void Setting_Dist_Load_2D(int mesh_n, int iPatch, int Total_Element, int iCoord, double val_Coord,
						  double Range_Coord[2], int type_load, double Coeff_Dist_Load[3]);

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
static void Calculation(int order_xi, int order_eta,
						int knot_n_xi, int knot_n_eta,
						int cntl_p_n_xi, int cntl_p_n_eta,
						double *input_knot_vec_xi, double *input_knot_vec_eta,
						double *cntl_px, double *cntl_py,
						double *disp_cntl_px, double *disp_cntl_py,
						double *weight);
static void Calculation_overlay(int order_xi_loc, int order_eta_loc,
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

//要素合成マトリックス
int Make_B_Matrix(int El_No, double B[D_MATRIX_SIZE][MAX_KIEL_SIZE], double Local_coord[DIMENSION], double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], double *J);
int Make_b_grad_Matrix(int El_No, double b_grad[DIMENSION * DIMENSION][2 * MAX_NO_CCpoint_ON_ELEMENT], double Local_coord[DIMENSION], double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], double *J);
int Make_D_Matrix_2D(double D[D_MATRIX_SIZE][D_MATRIX_SIZE], double E, double nu, int DM);

//Interaction integral
void Make_auxiliary_mode1(int e, double E, double nu, int DM, double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], double crack_front_coordinates_x, double crack_front_coordinates_y);
void Make_auxiliary_mode2(int e, double E, double nu, int DM, double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], double crack_front_coordinates_x, double crack_front_coordinates_y);

//gauss array
static int GP_1dir;						//1方向のガウス点数
static int GP_2D;						//2次元のガウス点数
static double Gxi[POW_Ng_extended][DIMENSION];	//ガウス点
static double w[POW_Ng_extended];				//ガウス点での重み

//static int DIMENSION;
static int KIEL_SIZE; //要素分割マトリックスの大きさ

static int Controlpoint_of_Element[MAX_N_ELEMENT][MAX_NO_CCpoint_ON_ELEMENT];
static double Node_Coordinate[MAX_N_NODE][DIMENSION + 1];
static double Equivalent_Nodal_Force[MAX_N_NODE][DIMENSION]; // Equivalent nodal forces arising from the distributed load
static int K_Whole_Ptr[MAX_K_WHOLE_SIZE + 1], K_Whole_Col[MAX_NON_ZERO];
//static int K_Whole_Ptr[MAX_N_MESH][MAX_K_WHOLE_SIZE + 1], K_Whole_Col[MAX_N_MESH][MAX_NON_ZERO];
static double K_Whole_Val[MAX_NON_ZERO];
//static double K_Whole_Val[MAX_N_MESH][MAX_NON_ZERO];
static int Index_Dof[MAX_K_WHOLE_SIZE];
static int INC[MAX_N_PATCH][MAX_N_ELEMENT][DIMENSION];
static int Adress_Controlpoint[MAX_N_PATCH][1000][1000]; //INCの配列をいじったものAdress_Controlpoint[ξ][η]；コントールポイント番号、任意のパッチ上でξ方向[]番目、η方向[]番目のコントロールポイント番号を示す
static int Order[MAX_N_PATCH][DIMENSION];
static int No_knot[MAX_N_PATCH][DIMENSION];
static int No_Control_point[MAX_N_PATCH][DIMENSION];
static double element_coordinate_Nopoint[MAX_N_ELEMENT][DIMENSION];
static double Gausspoint_coordinates[MAX_N_ELEMENT][POW_Ng_extended][DIMENSION];
//static int same_point[100];
static int same_point_in_Element[MAX_N_NODE];
static int Patch_controlpoint[MAX_N_PATCH][MAX_N_Controlpoint_in_Patch]; //バッチとコントロールポイント番号の要素コネクティビティ
static int Element_patch[MAX_N_ELEMENT];								 //要素がどのパッチに属しているか示す配列(要素番号は1つのモデルで通し番号)
static int No_Controlpoint_in_patch[MAX_N_PATCH];
static int No_Control_point_ON_ELEMENT[10000];

static int Node_To_Node[K_DIVISION_LENGE][10000], Total_Control_Point_To_Node[K_DIVISION_LENGE]; //ある節点に関係する節点番号s
//static int col_N[10][1000];

// static double b_grad_glo[MAX_ELEMENTS][POW_Ng][DIMENSION * DIMENSION][2 * POW_Ng];

static double sol_vec[MAX_K_WHOLE_SIZE];
static double rhs_vec[MAX_K_WHOLE_SIZE];
static double diag_scaling[MAX_K_WHOLE_SIZE];

static double Shape[DIMENSION][MAX_N_NODE][10];
static double shape_func[MAX_N_NODE];
static double dShape_func1[MAX_N_NODE];
static double dShape_func2[MAX_N_NODE];
static double dShape[DIMENSION][MAX_N_NODE];
static double Position_Knots[MAX_N_PATCH][DIMENSION][MAX_N_KNOT];
static double Position_Data_param[DIMENSION];

static double Displacement[MAX_K_WHOLE_SIZE];
// static double Disp_overlay[MAX_K_WHOLE_SIZE];
static double Strain[MAX_N_ELEMENT][POW_Ng_extended][N_STRAIN];
static double Strain_glo[MAX_N_ELEMENT][POW_Ng_extended][N_STRAIN];
static double Strain_overlay[MAX_N_ELEMENT][POW_Ng_extended][N_STRAIN];
static double Strain_aux_mode1[MAX_N_ELEMENT][POW_Ng_extended][N_STRAIN];
static double Strain_aux_mode2[MAX_N_ELEMENT][POW_Ng_extended][N_STRAIN];
static double Stress[MAX_N_ELEMENT][POW_Ng_extended][N_STRESS];
static double Stress_glo[MAX_N_ELEMENT][POW_Ng_extended][N_STRESS];
static double Stress_overlay[MAX_N_ELEMENT][POW_Ng_extended][N_STRESS];
static double Stress_aux_mode1_local[MAX_N_ELEMENT][POW_Ng_extended][N_STRESS];
static double Stress_aux_mode1[MAX_N_ELEMENT][POW_Ng_extended][N_STRESS];
static double Stress_aux_mode2_local[MAX_N_ELEMENT][POW_Ng_extended][N_STRESS];
static double Stress_aux_mode2[MAX_N_ELEMENT][POW_Ng_extended][N_STRESS];
static double StrainEnergyDensity[MAX_N_ELEMENT][POW_Ng_extended];
static double StrainEnergyDensity_overlay[MAX_N_ELEMENT][POW_Ng_extended];
static double StrainEnergyDensity_aux_mode1[MAX_N_ELEMENT][POW_Ng_extended];
static double StrainEnergyDensity_aux_mode2[MAX_N_ELEMENT][POW_Ng_extended];
static double Disp_grad[MAX_N_ELEMENT][POW_Ng_extended][DIMENSION * DIMENSION]; //Disp_grad[MAX_N_ELEMENT][POW_Ng_extended][0] = 𝜕u1/𝜕x1  Disp_grad[MAX_N_ELEMENT][POW_Ng_extended][1] = 𝜕u1/𝜕x2  Disp_grad[MAX_N_ELEMENT][POW_Ng_extended][2] = 𝜕u2/𝜕x1 Disp_grad[MAX_N_ELEMENT][POW_Ng_extended][3] = 𝜕u2/𝜕x2
static double Disp_grad_glo[MAX_N_ELEMENT][POW_Ng_extended][DIMENSION * DIMENSION];
static double Disp_grad_overlay[MAX_N_ELEMENT][POW_Ng_extended][DIMENSION * DIMENSION];
static double Disp_grad_aux_mode1_local[MAX_N_ELEMENT][POW_Ng_extended][DIMENSION * DIMENSION];
static double Disp_grad_aux_mode1[MAX_N_ELEMENT][POW_Ng_extended][DIMENSION * DIMENSION];
static double Disp_grad_aux_mode2_local[MAX_N_ELEMENT][POW_Ng_extended][DIMENSION * DIMENSION];
static double Disp_grad_aux_mode2[MAX_N_ELEMENT][POW_Ng_extended][DIMENSION * DIMENSION];
static double ReactionForce[MAX_K_WHOLE_SIZE];

static double difference[MAX_N_PATCH][MAX_N_KNOT][DIMENSION];		  /*隣り合うノットベクトルの差*/
static int ENC[MAX_N_PATCH][MAX_N_ELEMENT][DIMENSION];				  /*ENC[パッチ][全ての要素][0,1]=x,y方向の何番目の要素か*/
static int real_Total_Element;							  /*ゼロエレメントを除いた要素数*/
static int real_element[MAX_N_ELEMENT];					  /*ゼロエレメントではない要素の番号*/
static int Total_element_all_ID[MAX_N_ELEMENT];			  /*ゼロエレメントではない要素＝１、ゼロエレメント＝０*/
static int line_No_Total_element[MAX_N_PATCH][DIMENSION]; /*ゼロエレメントを含むすべての要素列の数*/
static int line_No_real_element[MAX_N_PATCH][DIMENSION];  /*ゼロエレメントではない要素列の数*/
static int real_element_line[MAX_N_PATCH][MAX_N_ELEMENT][DIMENSION];   /*ゼロエレメントではない要素列*/

static int No_points_for_colored_points; /*zarusobaで点に色付ける時の全ての点の数*/
/*
static double data_result_shape_x[10000];
static double data_result_shape_y[10000];
static double data_result_disp_x[10000];
static double data_result_disp_y[10000];
*/

static int No_points_for_new_zarusoba; /*zarusobaで点に色付ける時の全ての点の数*/
/*
static double data_result_shape_x_for_new_zarusoba[10000];
static double data_result_shape_y_for_new_zarusoba[10000];
static double data_result_disp_x_for_new_zarusoba[10000];
static double data_result_disp_y_for_new_zarusoba[10000];
*/

//static double Strain_x[MAX_N_ELEMENT][POW_Ng_extended][N_STRAIN];

//for s-IGA
static int Total_mesh;

static int Element_mesh[MAX_N_ELEMENT]; //要素がどのメッシュ内にあるかを示す配列
static int Patch_mesh[MAX_N_PATCH]; //パッチがどのメッシュ内にあるかを示す配列

static int Total_Patch_on_mesh[MAX_N_MESH];	//各メッシュ上のパッチ数
static int Total_Patch_to_mesh[MAX_N_MESH + 1];	//メッシュ[]までのパッチ数（メッシュ[]内のパッチ数は含まない）
static int Total_Patch_to_Now;				//現メッシュまでのパッチ数（現メッシュのパッチ数は含まない）

static int Total_Control_Point_on_mesh[MAX_N_MESH];	//各メッシュ上のコントロールポイント数
static int Total_Control_Point_to_mesh[MAX_N_MESH + 1];	//メッシュ[]までのコントロールポイント数（メッシュ[]内のコントロールポイント数は含まない）
static int Total_Control_Point_to_Now;				//現メッシュまでのコントロールポイント数（現メッシュのコントロールポイント数は含まない）

static int Total_Element_on_mesh[MAX_N_MESH];
static int Total_Element_to_mesh[MAX_N_MESH + 1];
static int Total_Element_to_Now;
//static int Total_Element_to_patch[MAX_N_PATCH];
static int El_No_on_mesh[MAX_N_MESH][MAX_N_ELEMENT];   //メッシュ内でのコントロールポイント配列

static int Total_Constraint_all_mesh;
static int Total_Constraint_on_mesh[MAX_N_MESH];
static int Total_Constraint_to_mesh[MAX_N_MESH + 1];
static int Total_Load_on_mesh[MAX_N_MESH];
static int Total_Load_to_mesh[MAX_N_MESH + 1];
static int Total_DistributeForce_on_mesh[MAX_N_MESH];
static int Total_DistributeForce_to_mesh[MAX_N_MESH + 1];
static int Constraint_ID[MAX_N_NODE * DIMENSION];
static int Constraint_Node_Dir_on_mesh[MAX_N_MESH][MAX_N_CONSTRAINT][2];
static double Value_of_Constraint_on_mesh[MAX_N_MESH][MAX_N_CONSTRAINT];

static double Control_Coord[DIMENSION][MAX_N_NODE];
static double Control_Weight[MAX_N_NODE];

static int real_Total_Element_on_mesh[MAX_N_MESH];
static int real_Total_Element_to_mesh[MAX_N_MESH + 1];
static int real_Total_Element_to_Now;
static int real_El_No_on_mesh[MAX_N_MESH][MAX_N_ELEMENT];

static int temp_element_n[MAX_N_ELEMENT_OVER_POINT];
static int element_n_point[MAX_N_ELEMENT_OVER_ELEMENT];
static int NNLOVER[MAX_N_ELEMENT];
//static int NELOVER[MAX_N_ELEMENT][MAX_N_ELEMENT_OVER_ELEMENT];
static int NELOVER[MAX_N_ELEMENT][MAX_N_ELEMENT_OVER];
static int Check_BDBJ_flag[MAX_N_ELEMENT];
static int Total_BDBJ_flag;
static int Same_BDBJ_flag[POW_Ng_extended];

static int iPatch_array[MAX_N_DISTRIBUTE_FORCE], iCoord_array[MAX_N_DISTRIBUTE_FORCE], type_load_array[MAX_N_DISTRIBUTE_FORCE];
static double val_Coord_array[MAX_N_DISTRIBUTE_FORCE], Range_Coord_array[MAX_N_DISTRIBUTE_FORCE][2], Coeff_Dist_Load_array[MAX_N_DISTRIBUTE_FORCE][3];

//for test
//static int shape_check_frag;

//for Interaction integral
static double T[DIMENSION][DIMENSION];
static double K_mode1;
static double K_mode2;
static double J_integral_value_aux_mode1;
static double J_integral_value_aux_mode2;

//重ね合わせの結果
double E;							//ヤング率(GPa)
double nu;							//ポアソン比(-)
int patch_n;						//パッチ数
int cntl_p_n;						//コントロールポイント数
int order_xi[MAX_PATCHES];		//ξ基底関数の次数(p)
int order_eta[MAX_PATCHES];		//η基底関数の次数(p)
int knot_n_xi[MAX_PATCHES];		//ξノットベクトルの数(n+p+1)
int knot_n_eta[MAX_PATCHES];		//ηノットベクトルの数(n+p+1)
int cntl_p_n_xi[MAX_PATCHES];	//ξ方向コントロールポイント数(n)
int cntl_p_n_eta[MAX_PATCHES];	//η方向コントロールポイント数(n)

double knot_vec_xi[MAX_PATCHES][MAX_KNOTS];		//ξノットベクトル
double knot_vec_eta[MAX_PATCHES][MAX_KNOTS];	//ηノットベクトル
double cntl_px[MAX_PATCHES][MAX_CNRL_P];		//コントロールポイントx座標
double cntl_py[MAX_PATCHES][MAX_CNRL_P];		//コントロールポイントy座標
double disp_cntl_px[MAX_PATCHES][MAX_CNRL_P];	//コントロールポイント上のx方向変位
double disp_cntl_py[MAX_PATCHES][MAX_CNRL_P];	//コントロールポイント上のy方向変位
double weight[MAX_PATCHES][MAX_CNRL_P];			//重み

double output_xi_loc[MAX_ELEMENTS][Ng];
double output_eta_loc[MAX_ELEMENTS][Ng];
double coord_x[MAX_POINTS][MAX_POINTS];		//メッシュx座標
double coord_y[MAX_POINTS][MAX_POINTS];		//メッシュy座標
double coord_x_gauss[MAX_POINTS][MAX_POINTS];		//メッシュx座標 for gauss
double coord_y_gauss[MAX_POINTS][MAX_POINTS];		//メッシュy座標 for gauss
double dxi_x[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];		// ∂x/∂ξ
double dxi_y[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];		// ∂y/∂ξ
double deta_x[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];		// ∂x/∂η
double deta_y[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];		// ∂y/∂η
double dxi_x_gauss[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];			// ∂x/∂ξ for Gauss
double dxi_y_gauss[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];			// ∂y/∂ξ for Gauss
double deta_x_gauss[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];		// ∂x/∂η for Gauss
double deta_y_gauss[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];		// ∂y/∂η for Gauss

double disp_px_glo[MAX_ELEMENTS][MAX_N_Controlpoint_in_Patch];
double disp_py_glo[MAX_ELEMENTS][MAX_N_Controlpoint_in_Patch];
double disp_x[MAX_POINTS][MAX_POINTS];			//x方向変位
double disp_y[MAX_POINTS][MAX_POINTS];			//y方向変位
double disp_x_glo_gauss[MAX_ELEMENTS * Ng][MAX_ELEMENTS * Ng];			//x方向変位 for Gauss
double disp_y_glo_gauss[MAX_ELEMENTS * Ng][MAX_ELEMENTS * Ng];			//y方向変位 for Gauss
double dxi_disp_x[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	// ∂u/∂ξ
double dxi_disp_y[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	// ∂v/∂ξ
double deta_disp_x[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	// ∂u/∂η
double deta_disp_y[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	// ∂v/∂η

double strain_xx[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//x方向ひずみ
double strain_yy[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//y方向ひずみ
double strain_xy[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//剪断ひずみ
double strain_xx_gauss[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//x方向ひずみ for Gauss
double strain_yy_gauss[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//y方向ひずみ for Gauss
double strain_xy_gauss[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//剪断ひずみ for Gauss

double stress_xx[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//x方向垂直応力
double stress_yy[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//y方向垂直応力
double stress_xy[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//剪断応力
double stress_xx_gauss[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//x方向垂直応力
double stress_yy_gauss[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//y方向垂直応力
double stress_xy_gauss[MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//剪断応力

int fields_flag = 1;		//s-IGAのためのNURBS_inputでは変位データは必ず読み込ませる
int division_ele_xi;		//ξ方向の一要素あたりの分割数
int division_ele_eta;		//η方向の一要素あたりの分割数
int division_n_xi;		//ξ方向の表示する点の数
int division_n_eta;		//η方向の表示する点の数
int element_n_xi;			//ξ方向要素数
int element_n_eta;		//η方向要素数

int temp_index[MAX_PATCHES][MAX_CNRL_P];
double temp_cntl_px[MAX_CNRL_P];
double temp_cntl_py[MAX_CNRL_P];
double temp_weight[MAX_CNRL_P];
double temp_disp_x[MAX_CNRL_P];
double temp_disp_y[MAX_CNRL_P];

//for s-IGA
int n_patch_glo;	//グローバルメッシュ上のパッチ数
int n_patch_loc;	//ローカルメッシュ上のパッチ数
int glo_cntl_p_n;	//グローバルメッシュ上のコントロールポイント数
int loc_cntl_p_n;	//ローカルメッシュ上のコントロールポイント数

//for graph
int graph_patch_n;	//グラフ作成用出力ファイル内のパッチ番号

//for GP info
// static double coordinate_GP[MAX_ELEMENTS*MAX_ELEMENTS][POW_Ng_extended][DIMENSION];
// static double strain_GP[MAX_ELEMENTS*MAX_ELEMENTS][POW_Ng_extended][3];
// static double stress_GP[MAX_ELEMENTS*MAX_ELEMENTS][POW_Ng_extended][3];
// static double stress_r_theta_GP[MAX_ELEMENTS*MAX_ELEMENTS][POW_Ng_extended][3];
// static double stress_theory_r_theta[MAX_ELEMENTS*MAX_ELEMENTS][POW_Ng_extended][3];
// static double Jac[MAX_ELEMENTS*MAX_ELEMENTS][POW_Ng_extended];

//解析条件パラメータの設定
static int DM = 1;		//平面応力状態:DM=0	平面ひずみ状態:DM=1
static int check_over_parameter = 1;	//要素の重なりの判定(要素の物体上の端点:0 ガウス点:1)

FILE *fp;

int main(int argc, char *argv[])
{
	clock_t start, end, t1;

	int i, j, k;
	//int l;
	int re;
	//int p;
	//int q, r;
	int Total_Element;
	int Total_Control_Point;
	int No_Patch = 0;
	//int Total_net = 0;
	static int Total_Load = 0, Load_Node_Dir[MAX_N_LOAD][2];
	static double Value_of_Load[MAX_N_LOAD];
	static int Total_Constraint = 0, Constraint_Node_Dir[MAX_N_CONSTRAINT][2];
	static double Value_of_Constraint[MAX_N_CONSTRAINT];
	static int Total_DistributeForce = 0;
	int K_Whole_Size = 0;
	int El_No = 0;
	// static double element[DIMENSION];
	//static double element_emsh[DIMENSION];
    static double element_loc[DIMENSION];
	static double E, nu;
	static int max_itr;

    //for s-IGA
    int tm;
    //int K_Whole_overlaid_Size = 0;

	// 複数ファイル読み込みのためコメントアウト
	//for s-IGA
    //引数の個数確認
	if (argc <= 1)
	{
		printf("Argument is missing\n");
	}
    if (argc == 2)  /*通常IGA：input file 1つ*/
    {
        printf("IGA carried out.(No local mesh)\n");
    }
    if (argc >= 3)  /*s-IGA：input file 複数*/
    {
        printf("s-IGA carried out.(%d local meshes)\n", argc - 2);
    }
    Total_mesh = argc - 1;

	start = clock();
    //ローカルメッシュのためにinput file読み込みのループ
    for ( tm = 0 ; tm < Total_mesh ; tm++ )
    {
	    Get_InputData(tm,&E, &nu, &Total_Element, &Total_Control_Point, &Total_Load, &No_Patch, Load_Node_Dir, Value_of_Load, &Total_Constraint, Constraint_Node_Dir, Value_of_Constraint, &Total_DistributeForce, argv);
	    //printf("Finish Get_InputData\n");
        if(tm==0 && Total_mesh>1)
        {
            printf("Finish Get_InputData(Global mesh:%s)\n", argv[1]);
        }
        if(tm>0)
        {
            printf("Finish Get_inputData(Local mesh No.[%d]:%s)\n", tm, argv[tm+1]);
        }
	    printf("Total Element=%d Node=%d Constraint=%d Load=%d\n", Total_Element, Total_Control_Point, Total_Constraint, Total_Load);
	    printf("E;%le nu;%le\n", E, nu);
		//printf("Finish Get_InputData\n");
	}	//ローカルメッシュのためにinput fileの読み込みのループ終了
	//グローバルパッチとローカルパッチの重なり
	//ローカルメッシュ同士は被りなしと仮定

	printf("\ncheck_over_parameter;%d\n\n", check_over_parameter);
	for (i = 1; i < Total_mesh; i++)
	{
		printf("mesh_n_org;0\tmesh_n_over;%d\n", i);
		//NNLOVER[over_ele][]=org_eleの算出
		if (check_over_parameter == 0)
		{
			Check_coupled_Glo_Loc_element_for_end(element_loc, i, 0);
		}
		if (check_over_parameter == 1)
		{
			Check_coupled_Glo_Loc_element_for_Gauss(element_loc, i, 0);
		}
		//NNLOVER[org_ele][]=over_eleの算出
		Make_Loc_Glo();
	}
	/*
	for (i = 0; i < Total_mesh - 1; i++)	//mesh_n_org
	{
		for (j = i + 1; j < Total_mesh; j++)	//mesh_n_over
		{
			printf("mesh_n_org;%d\tmesh_n_over;%d\n",i,j);
			if (check_over_parameter == 0)
			{
				Check_coupled_Glo_Loc_element_for_end(element_loc,
											  		  Total_Element,
											  		  Total_Control_Point,
											  		  j,i);
			}
			if (check_over_parameter == 1)
			{
				Check_coupled_Glo_Loc_element_for_Gauss(element_loc,
														Total_Element,
														Total_Control_Point,
														j,i);
			}
		}
	}
	*/

    /////////////////////全体剛性マトリックスの制作////////////////////////////
	K_Whole_Size = Make_Index_Dof(Total_Control_Point_to_mesh[Total_mesh],
								  Total_Constraint_to_mesh[Total_mesh],
								  Constraint_Node_Dir);
    //K_Whole_Size_array[tm]=K_Whole_Size;
    /*
	for(i=0;i<tm;i++)
    {
        K_Whole_Size_array_to_Now += K_Whole_Size_array[i];
    }
	*/
    printf("K_Whole_Size=%d\n",K_Whole_Size);
    Make_K_Whole_Ptr_Col(Total_Element_to_mesh[Total_mesh], Total_Control_Point_to_mesh[Total_mesh], K_Whole_Size);
    Make_K_Whole_Val(E, nu, real_Total_Element_to_mesh[Total_mesh], DM);
    printf("Finish Make_K_Whole\n");
	/*
	if(tm==0 && Total_mesh>1)
    {
        printf("Finish Make_K_Whole(Global mesh:%s)\n",argv[1]);
    }
    if(tm>0)
    {
        printf("Finish Make_K_Whole(Local mesh No.[%d]:%s)\n",tm,argv[1+tm]);
    }
	*/
	//for s-IGA　複数メッシュのループ内に移動
	//荷重ベクトルの算出部分
	for (i = 0; i < Total_Load_to_mesh[Total_mesh]; i++)
	{
		// printf("Value_of_Load;%11.10e\n", Value_of_Load[i]);
		printf("Value_of_Load;%.20e\n", Value_of_Load[i]);
	}
	printf("pp\n");
	//max_itr = K_Whole_Size;
	//printf("K_Whole_Size:%d\n", K_Whole_Size);
	Make_F_Vec(Total_Load_to_mesh[Total_mesh],
			   Load_Node_Dir, Value_of_Load, K_Whole_Size);
	Make_F_Vec_disp_const(tm, Total_Constraint_to_mesh[Total_mesh],
						  Constraint_Node_Dir, Value_of_Constraint,
						  E, nu, DM);
	Add_Equivalent_Nodal_Forec_to_F_Vec(Total_Control_Point_to_mesh[Total_mesh]);
    //K_Whole_overlaid_Size += K_Whole_Size_array[tm];
    //printf("K_Whole_Size_array[%d]=%d\n",tm,K_Whole_Size_array[tm]);
    //}   //ローカルメッシュのためにinput fileの読み込みのループ終了
    //printf("K_Whole_overlaid_Size=%d\n",K_Whole_overlaid_Size);
    printf("Finish Make_K_Whole\n");
	/*
	//グローバルパッチとローカルパッチの重なり
	for (i = 0; i < Total_mesh - 1; i++)	//mesh_n_org
	{
		for (j = i + 1; j < Total_mesh; j++)	//mesh_n_over
		{
			printf("mesh_n_org;%d\tmesh_n_over;%d\n",i,j);
			Check_coupled_Glo_Loc_element(element_loc, Total_Element, Total_Control_Point,j,i);
		}
	}
    //Check_coupled_Glo_Loc_element(element_loc, Total_Element, Total_Control_Point,1,0);
    */
	//Make_K_Whole_overlaid(K_Whole_overlaid_Size);
    //printf("Finish Make_K_Whole_overlaid\n");
	//printf("check\n");
    //printf("K_Whole_overlaid_Size=%d\n",K_Whole_overlaid_Size);
	///////////////連立一次方程式/////////////////////////////////////////
	//printf("DistributeForce=%d\n", DistributeForce);
	//printf("Val_DistributeForce=%d\n", Val_DistributeForce);
	//printf("Load_Node_Dir=%s\n", );
	/*		Force_Dis(Total_DistributeForce, DistributeForce, Val_DistributeForce, &Total_Load, Load_Node_Dir, Value_of_Load, Total_Control_Point,El_No,&Total_Element);*/

	//s-IGAのため複数メッシュのループ内に移動
	/*
	for (i = 0; i < Total_Load; i++)
		printf("%11.10e\n", Value_of_Load[i]);
	printf("pp");
	max_itr = K_Whole_Size;
	printf("K_Whole_Size:%d\n", K_Whole_Size);
	Make_F_Vec(Total_Load, Load_Node_Dir, Value_of_Load, K_Whole_Size);
	Make_F_Vec_disp_const(Total_Constraint, Constraint_Node_Dir, Value_of_Constraint, Total_Element, E, nu, DM, Total_Control_Point);
	Add_Equivalent_Nodal_Forec_to_F_Vec(Total_Control_Point);
	*/

	// Kマトリックスのsvg出力
	// K_output_svg(K_Whole_Size);

    //反復回数の設定
    max_itr = K_Whole_Size;
	// max_itr = K_Whole_Size * 5;

	// CG法
	// Diag_Scaling_CG_pre(K_Whole_Size, 0);
    // printf("Finish 1st Diag_Scaling_CG_Pre\n");
	// CG_Solver(K_Whole_Size, max_itr, EPS, 0);
	// Diag_Scaling_CG_pre(K_Whole_Size, 1);
	// printf("Finish CG_Solver\n");

	// PCG法
    printf("\nStart PCG solver\n");
	PCG_Solver(K_Whole_Size, max_itr, EPS);
	printf("Finish PCG solver\n\n");

	/////////////変位と歪と応力//////////////////////////////////////
    //for s-IGA
    //Total_Control_Point_to_Now += Total_Control_Point_on_mesh[Total_mesh-1];
    //printf("Total_Control_Point_to_Now;%d\n",Total_Control_Point_to_Now);
    /*
	for(i = 0 ; i < Total_mesh ; i++)
    {
        //printf("Total_Control_Point_to_mesh[%d]=%d\n",i,Total_Control_Point_to_mesh[i]);
        for(j = 0 ; j < Total_Constraint_to_mesh[Total_mesh] ; j++)
        {
            Constraint_ID[(Constraint_Node_Dir_on_mesh[i][j][0] + Total_Control_Point_to_mesh[i+1]) * DIMENSION + Constraint_Node_Dir_on_mesh[i][j][1]] = 1;
            Displacement[(Constraint_Node_Dir_on_mesh[i][j][0] + Total_Control_Point_to_mesh[i+1]) * DIMENSION + Constraint_Node_Dir_on_mesh[i][j][1]]
                = Value_of_Constraint_on_mesh[i][j];
            //printf("Constraint_Node_Dir_on_mesh[i][j][0] + Total_Control_Point_to_mesh[i+1]=%d\n",Constraint_Node_Dir_on_mesh[i][j][0] + Total_Control_Point_to_mesh[i+1]);
        }
        //printf("i=%d\n",i);
    }*/
    for(i = 0; i < Total_Constraint_to_mesh[Total_mesh]; i++)
    {
        Constraint_ID[Constraint_Node_Dir[i][0] * DIMENSION + Constraint_Node_Dir[i][1]] = 1;
        Displacement[Constraint_Node_Dir[i][0] * DIMENSION + Constraint_Node_Dir[i][1]]
            = Value_of_Constraint[i];
        //printf("Constraint_Node_Dir_on_mesh[i][j][0] + Total_Control_Point_to_mesh[i+1]=%d\n",Constraint_Node_Dir_on_mesh[i][j][0] + Total_Control_Point_to_mesh[i+1]);
    }

	for(i = 0; i < Total_Control_Point_to_mesh[Total_mesh]; i++)
	{
		for(j = 0; j < DIMENSION; j++)
		{

            int index = Index_Dof[i * DIMENSION + j];
			if (index >= 0)
				Displacement[i * DIMENSION + j] = sol_vec[index];

            // printf("%d\t%le\t", i * DIMENSION + j, Displacement[i * DIMENSION + j]);
		}
        // printf("\n");
	}
	printf("Finish Make_Displacement\n");
	end = clock();
	printf("Analysis time:%.2f[s]\n",(double)(end-start)/CLOCKS_PER_SEC);
	// Make_Strain(E, nu, Total_Element_to_mesh[Total_mesh],
	// 			El_No, Total_Control_Point);
	// printf("Finish Make_Strain\n");
	// Make_Stress_2D(E, nu, Total_Element, DM);
	// printf("Finish Make_Stress\n");
	// Make_ReactionForce(Total_Element, Total_Control_Point, El_No);
	// printf("Finish Make_ReactionForce\n");
	// puts("sol_vec");
	// Make_Parameter_z(Total_Element, E, nu, DM);
	// printf("Finish Make_Parameter_z\n");
	Make_Strain(real_Total_Element_to_mesh[Total_mesh]);
	printf("Finish Make_Strain\n");
	Make_Stress_2D(E, nu, real_Total_Element_to_mesh[Total_mesh], DM);
	printf("Finish Make_Stress\n");
	Make_ReactionForce(Total_Control_Point_to_mesh[Total_mesh]);
	printf("Finish Make_ReactionForce\n");
	puts("sol_vec");
	Make_Parameter_z(real_Total_Element_to_mesh[Total_mesh], E, nu, DM);
	printf("Finish Make_Parameter_z\n");
	//Make_Output( Total_Control_Point, Total_Element );
	//printf("Finish Make_Output\n" );

	//Make_Strain_x_Quad_4( E, nu, Total_Element);
	//printf("Finish Make_Strain_x_Quad_4\n");
	//Make_EMT(E, nu, Total_Element);
	//printf("Finish Make_EMT\n");
	/*

	for(i = 0; i <K_Whole_Size; i++ ){
		printf("%le ",sol_vec[i]);
		printf("\n");
	}
	puts("DIMENSION");
	for(j = 0; j < Total_Control_Point; j++ ){
		for(i = 0; i < DIMENSION; i++ )
			printf("%le ",Displacement[j*DIMENSION+i]);
		printf("\n");
	}
	printf("\n");
	puts("Strain");
	for( k = 0; k < POW_Ng; k++){
		for( i = 0; i < Total_Element; i++ ){
			for( j = 0; j < N_STRAIN; j++ )
				printf("%le ",Strain[i][k][j]);
			printf("\n");
		}
	}
	printf("\n");
	puts("Stress");
	for( k = 0; k < POW_Ng; k++){
		for( i = 0; i < Total_Element; i++ ){
			for( j = 0; j < N_STRESS; j++ )
				printf("%le ",Stress[i][k][j]);
			printf("\n");
		}
	}
	printf("\n");
	puts("ReactionForce");
	for(j = 0; j < Total_Control_Point; j++ ){
		for(i = 0; i < DIMENSION; i++ )
			printf("%le ",ReactionForce[j*DIMENSION+i]);
		printf("\n");
	}

	*/
	//全てのローカルパッチについて1つのinput.txt作成
	//(重ね合わせ結果出力のため)
	int l_patch, g_patch;
	int l_cntl_p, g_cntl_p;
	int l_cnst, g_cnst;
	int l_load, g_load;
	int l_dist, g_dist;
	int l_temp;

	g_patch = Total_Patch_to_mesh[1];
	g_cntl_p = Total_Control_Point_to_mesh[1];
	g_cnst = Total_Constraint_to_mesh[1];
	g_load = Total_Load_to_mesh[1];
	g_dist = Total_DistributeForce_to_mesh[1];

	fp = fopen("input_local.txt", "w");
	fprintf(fp, "%.3f\t%.6f\n\n", E, nu);

	l_patch = Total_Patch_to_mesh[Total_mesh] - g_patch;
	fprintf(fp, "%d\n\n", l_patch);

	l_cntl_p = Total_Control_Point_to_mesh[Total_mesh] - g_cntl_p;
	fprintf(fp, "%d\n\n", l_cntl_p);

	for (i = 0; i < l_patch; i++)
	{
		fprintf(fp, "%d\t%d\n",
					Order[i + g_patch][0], Order[i + g_patch][1]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < l_patch; i++)
	{
		fprintf(fp, "%d\t%d\n",
					No_knot[i + g_patch][0], No_knot[i + g_patch][1]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < l_patch; i++)
	{
		fprintf(fp, "%d\t%d\n",
					No_Control_point[i + g_patch][0],
					No_Control_point[i + g_patch][1]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < l_patch; i++)
	{
		for (j = 0; j < No_Controlpoint_in_patch[i + g_patch]; j++)
		{
			l_temp = Patch_controlpoint[i + g_patch][j] - g_cntl_p;
			fprintf(fp,"%d\t", l_temp);
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n");

	l_cnst = Total_Constraint_to_mesh[Total_mesh] - g_cnst;
	l_load = Total_Load_to_mesh[Total_mesh] - g_load;
	l_dist = Total_DistributeForce_to_mesh[Total_mesh] - g_dist;
	fprintf(fp,"%d\t%d\t%d\n\n",
				l_cnst, l_load,	l_dist);

	for (i = 0; i < l_patch; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			for (k = 0; k < No_knot[i + g_patch][j]; k++)
			{
				fprintf(fp,"%.5f\t",Position_Knots[i + g_patch][j][k]);
			}
			fprintf(fp,"\n");
		}
	}
	fprintf(fp,"\n");

	for (i = 0; i < l_cntl_p; i++)
	{
		fprintf(fp,"%d\t",i);
		for (j = 0; j < DIMENSION+1; j++)
		{
			fprintf(fp,"%.10f\t",Node_Coordinate[i + g_cntl_p][j]);
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n");

	for (i = 0; i < l_cnst; i++)
	{
		l_temp = Constraint_Node_Dir[i + g_cnst][0] - g_cntl_p;
		fprintf(fp,"%d\t%d\t%.10f\n",
					l_temp,
					Constraint_Node_Dir[i + g_cnst][1],
					Value_of_Constraint[i + g_cnst]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < l_load; i++)
	{
		l_temp = Load_Node_Dir[i + g_load][0] - g_cntl_p;
		fprintf(fp,"%d\t%d\t%.10f\n",
					l_temp,
					Load_Node_Dir[i + g_load][1],
					Value_of_Load[i + g_load]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < l_dist; i++)
	{
		fprintf(fp,"%d\t%d\t%d\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n",
					type_load_array[i + g_dist],
					iPatch_array[i + g_dist],
					iCoord_array[i + g_dist],
					val_Coord_array[i + g_dist],
					Range_Coord_array[i + g_dist][0],
					Range_Coord_array[i + g_dist][1],
					Coeff_Dist_Load_array[i + g_dist][0],
					Coeff_Dist_Load_array[i + g_dist][1],
					Coeff_Dist_Load_array[i + g_dist][2]);
	}
	fclose(fp);

	fp = fopen("input_for_NURBS.txt", "w");
	fprintf(fp, "%.3f\t%.6f\n\n", E, nu);

	fprintf(fp, "%d\n\n", Total_Patch_to_mesh[Total_mesh]);

	fprintf(fp, "%d\n\n", Total_Control_Point_to_mesh[Total_mesh]);

	for (i = 0; i < Total_Patch_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp, "%d\t%d\n", Order[i][0], Order[i][1]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < Total_Patch_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp, "%d\t%d\n", No_knot[i][0], No_knot[i][1]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < Total_Patch_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp, "%d\t%d\n", No_Control_point[i][0], No_Control_point[i][1]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < Total_Patch_to_mesh[Total_mesh]; i++)
	{
		for (j = 0; j < No_Controlpoint_in_patch[i]; j++)
		{
			fprintf(fp,"%d\t",Patch_controlpoint[i][j]);
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n");

	fprintf(fp,"%d\t%d\t%d\n\n",
				Total_Constraint_to_mesh[Total_mesh],
				Total_Load_to_mesh[Total_mesh],
				Total_DistributeForce_to_mesh[Total_mesh]);

	for (i = 0; i < Total_Patch_to_mesh[Total_mesh]; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			for (k = 0; k < No_knot[i][j]; k++)
			{
				fprintf(fp,"%.5f\t",Position_Knots[i][j][k]);
			}
			fprintf(fp,"\n");
		}
	}
	fprintf(fp,"\n");

	for (i = 0; i < Total_Control_Point_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp,"%d\t",i);
		for (j = 0; j < DIMENSION+1; j++)
		{
			fprintf(fp,"%.10f\t",Node_Coordinate[i][j]);
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"\n");

	for (i = 0; i < Total_Constraint_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp,"%d\t%d\t%.10f\n",
					Constraint_Node_Dir[i][0],
					Constraint_Node_Dir[i][1],
					Value_of_Constraint[i]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < Total_Load_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp,"%d\t%d\t%.10f\n",
					Load_Node_Dir[i][0],
					Load_Node_Dir[i][1],
					Value_of_Load[i]);
	}
	fprintf(fp,"\n");

	for (i = 0; i < Total_DistributeForce_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp,"%d\t%d\t%d\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n",
					type_load_array[i],
					iPatch_array[i],
					iCoord_array[i],
					val_Coord_array[i],
					Range_Coord_array[i][0],
					Range_Coord_array[i][1],
					Coeff_Dist_Load_array[i][0],
					Coeff_Dist_Load_array[i][1],
					Coeff_Dist_Load_array[i][2]);
	}

	fclose(fp);

	fp = fopen("coord_data.txt","w");
	for (i = 0; i < Total_Control_Point_to_mesh[Total_mesh]; i++)
	{
		fprintf(fp,"%d\t",i);
		for (j = 0; j < DIMENSION; j++)
		{
			fprintf(fp,"%.10e\t",Node_Coordinate[i][j]);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);

	// fp = fopen("checkAns/checkAns.txt", "w");
	// for (i = 0; i < K_Whole_Size; i++)
	// {
	// 	fprintf(fp, "%le ", sol_vec[i]);
	// 	fprintf(fp, "\n");
	// }
	// fprintf(fp, "\n\n\nDisplacement\n");
	// for (j = 0; j < Total_Control_Point_to_mesh[Total_mesh]; j++)
	// {
	// 	fprintf(fp, "%d\t", j);
	// 	for (i = 0; i < DIMENSION; i++)
	// 		fprintf(fp, "%.13e\t ", Displacement[j * DIMENSION + i]);
	// 	fprintf(fp, "\n");
	// }
	// fprintf(fp, "\n\n\nStrain\n");
	// //for( i = 0; i < Total_Element; i++ ){
	// // for (re = 0; re < real_Total_Element; re++)
	// for (re = 0; re < real_Total_Element_to_mesh[Total_mesh]; re++)
	// {
	// 	i = real_element[re];
	// 	for (k = 0; k < POW_Ng; k++)
	// 	{
	// 		fprintf(fp, "%d\t%d\t%d\t", re, i, k);
	// 		for (j = 0; j < N_STRAIN; j++)
	// 			fprintf(fp, "%.13e\t", Strain[i][k][j]);
	// 		fprintf(fp, "\n");
	// 	}
	// }
	// /*	fprintf(fp,"\n\n\n\n\n\n\n\n\n\n\nStrain_x\n");
	// for( k = 0; k < POW_Ng; k++){
	// 	for( i = 0; i < Total_Element; i++ ){
	// 		for( j = 0; j < DIMENSION; j++ )
	// 			fprintf(fp,"%d	%d	%le\t ",k,i,Strain_x[i][k][j]);
	// 		fprintf(fp,"\n");
	// 	}
	// }*/

	// fprintf(fp, "\n\n\n\n\n\n\n\n\n\n\nStress\n");
	// //for( i = 0; i < Total_Element; i++ ){
	// // for (re = 0; re < real_Total_Element; re++)
	// for (re = 0; re < real_Total_Element_to_mesh[Total_mesh]; re++)
	// {
	// 	i = real_element[re];
	// 	for (k = 0; k < POW_Ng; k++)
	// 	{
	// 		fprintf(fp, "%d\t%d\t", i, k);
	// 		for (j = 0; j < N_STRESS; j++)
	// 			fprintf(fp, "%.13e\t", Stress[i][k][j]);
	// 		fprintf(fp, "\n");
	// 	}
	// }
	// fprintf(fp, "\n\n\nReaction Force\n");
	// // for (j = 0; j < Total_Control_Point; j++)
	// for (j = 0; j < Total_Control_Point_to_mesh[Total_mesh]; j++)
	// {
	// 	for (i = 0; i < DIMENSION; i++)
	// 		fprintf(fp, "%.13e\t ", ReactionForce[j * DIMENSION + i]);
	// 	fprintf(fp, "\n");
	// }
	// fclose(fp);

	fp = fopen("Displacement.dat", "w");
	fprintf(fp, "label=Displacement\n");
	fprintf(fp, "num_items=%d\n", Total_Control_Point_to_mesh[Total_mesh]);
	fprintf(fp, "\n");
	for (j = 0; j < Total_Control_Point_to_mesh[Total_mesh]; j++)
	{
		// fprintf(fp, "%d:	%le %le ", j, Displacement[j * DIMENSION + 0], Displacement[j * DIMENSION + 1]);
		fprintf(fp, "%d:	%.16e %.16e ", j, Displacement[j * DIMENSION + 0], Displacement[j * DIMENSION + 1]);
		fprintf(fp, "\n");
	}
	fclose(fp);

/*
	fp = fopen("Strain@IntegrationPoint.dat", "w");
	fprintf(fp, "label=Strain@IntegrationPoint\n");
	fprintf(fp, "num_items=%d\n", real_Total_Element);
	fprintf(fp, "\n");
	//for( i = 0; i < Total_Element; i++ ){
	for (re = 0; re < real_Total_Element; re++)
	{
		i = real_element[re];
		fprintf(fp, "%d:\t", i);
		for (k = 0; k < POW_Ng; k++)
		{
			fprintf(fp, "%13e\t%13e\t%13e\t%13e\t", Strain[i][k][0], Strain[i][k][1], Strain[i][k][2], Strain[i][k][3]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	fp = fopen("Stress@IntegrationPoint.dat", "w");
	fprintf(fp, "label=Stress@IntegrationPoint\n");
	fprintf(fp, "num_items=%d\n", real_Total_Element);
	fprintf(fp, "\n");

	//for( i = 0; i < Total_Element; i++ ){
	for (re = 0; re < real_Total_Element; re++)
	{
		i = real_element[re];
		fprintf(fp, "%d\t:", i);
		for (k = 0; k < POW_Ng; k++)
		{
			fprintf(fp, "%13e\t%13e\t%13e\t%13e\t", Stress[i][k][0], Stress[i][k][1], Stress[i][k][2], Stress[i][k][3]);
		}
		fprintf(fp, "\n");
	}

	fclose(fp);

	fp = fopen("Disp_grad@IntegrationPoint.dat", "w");
	fprintf(fp,"label=Disp_grad@IntegrationPoint\n");
	fprintf(fp, "num_items=%d\n", real_Total_Element_to_mesh[Total_mesh]);
	fprintf( fp, "\n");

	//for( i = 0; i < Total_Element; i++ ){
	for(re = 0; re < real_Total_Element_to_mesh[Total_mesh]; re++){
	i = real_element[re];
		fprintf(fp,"%d\t:\t",i );
		for( k = 0; k < POW_Ng; k++){
				fprintf(fp,"%13e\t%13e\t%13e\t%13e\t",Disp_grad[i][k][0],Disp_grad[i][k][1],Disp_grad[i][k][2],Disp_grad[i][k][3]);
			}fprintf(fp,"\n");
		}

	fclose(fp);

	fp = fopen("ReactionForce.dat", "w");
	fprintf(fp, "label=ReactionForce\n");
	fprintf(fp, "num_items=%d\n", Total_Control_Point);
	fprintf(fp, "\n");
	for (j = 0; j < Total_Control_Point; j++)
	{
		fprintf(fp, "%d:	%13e %13e ", j, ReactionForce[j * DIMENSION + 0], ReactionForce[j * DIMENSION + 1]);
		fprintf(fp, "\n");
	}

	fclose(fp);

	fp = fopen("StrainEnergyDensity@IntegrationPoint.dat", "w");
	fprintf(fp,"label=StrainEnergyDensity@IntegrationPoint\n");
	fprintf(fp, "num_items=%d\n", real_Total_Element_to_mesh[Total_mesh]);
	fprintf( fp, "\n");
	//for( i = 0; i < Total_Element; i++ ){
	for(re = 0; re < real_Total_Element_to_mesh[Total_mesh]; re++){
		i=real_element[re];
		fprintf(fp,"%d\t:\t", i);
		for( k = 0; k < POW_Ng; k++){
			fprintf(fp,"%13e\t",StrainEnergyDensity[i][k]);
			}fprintf(fp,"\n");
		}
		
	fclose(fp);

	fp = fopen("mesh.msh", "w");

	fprintf(fp, "%d\n", real_Total_Element);

	//for( i = 0; i < Total_Element; i ++ ){
	for (re = 0; re < real_Total_Element; re++)
	{
		i = real_element[re];
		fprintf(fp, "%d %d %d %d %d %d %d %d %d", Controlpoint_of_Element[i][8], Controlpoint_of_Element[i][6], Controlpoint_of_Element[i][0], Controlpoint_of_Element[i][2], Controlpoint_of_Element[i][7], Controlpoint_of_Element[i][3], Controlpoint_of_Element[i][1], Controlpoint_of_Element[i][5], Controlpoint_of_Element[i][4]);
		//fprintf(fp,"%d %d %d %d %d %d %d %d %d",Controlpoint_of_Element[i][0],Controlpoint_of_Element[i][2],Controlpoint_of_Element[i][8],Controlpoint_of_Element[i][6],Controlpoint_of_Element[i][1],Controlpoint_of_Element[i][5],Controlpoint_of_Element[i][7],Controlpoint_of_Element[i][3],Controlpoint_of_Element[i][4] );
		fprintf(fp, "\n");
	}
	fprintf(fp, "%d\n", Total_Control_Point);
	for (i = 0; i < Total_Control_Point; i++)
	{
		fprintf(fp, "%lf %lf\n", Node_Coordinate[i][0], Node_Coordinate[i][1]);
	}
	fprintf(fp, "\n");

	fclose(fp);

	fp = fopen("mesh_net/control_net.msh", "w");
	for (l = 0; l < No_Patch; l++)
	{
		Total_net += (No_Control_point[l][0] - 1) * (No_Control_point[l][1] - 1);
	}

	fprintf(fp, "%d\n", Total_net);
	for (l = 0; l < No_Patch; l++)
	{
		for (j = 0; j < No_Control_point[l][1] - 1; j++)
		{
			for (i = 0; i < No_Control_point[l][0] - 1; i++)
			{
				fprintf(fp, "%d %d %d %d\n", Adress_Controlpoint[l][i][j], Adress_Controlpoint[l][i + 1][j], Adress_Controlpoint[l][i + 1][j + 1], Adress_Controlpoint[l][i][j + 1]);
			}
		}
	}

	fprintf(fp, "%d\n", Total_Control_Point);

	for (i = 0; i < Total_Control_Point; i++)
	{
		fprintf(fp, "%lf %lf\n", Node_Coordinate[i][0], Node_Coordinate[i][1]);
	}
	fclose(fp);

	fp = fopen("mesh_net/element.msh", "w");

	fprintf(fp, "%d\n", Total_Element);
	element_coordinate(Total_Element, Total_Control_Point);
	for (i = 0; i < Total_Element * 9; i += 9)
	{
		fprintf(fp, "%d %d %d %d %d %d %d %d %d\n", same_point_in_Element[i], same_point_in_Element[i + 1], same_point_in_Element[i + 2], same_point_in_Element[i + 3], same_point_in_Element[i + 4],
				same_point_in_Element[i + 5], same_point_in_Element[i + 6], same_point_in_Element[i + 7], same_point_in_Element[i + 8]);
	}
	fprintf(fp, "%d\n", Total_Element * 9);
	for (i = 0; i < Total_Element * 9; i++)
	{
		fprintf(fp, "%lf %lf\n", element_coordinate_Nopoint[i][0], element_coordinate_Nopoint[i][1]);
	}
	fclose(fp);
	*/

	//////////////////////
	///////NURBS//////////
	//////////////////////
	// fp = fopen("NURBS/NURBS_disp.dat", "w");
	// fprintf(fp, "%d\t", real_Total_Element_to_mesh[Total_mesh]);
	// fclose(fp);

	// calculate_Controlpoint_using_NURBS(element,
	// 								   real_Total_Element_to_mesh[Total_mesh],
	// 								   Total_Control_Point_to_mesh[Total_mesh]);

	// //fclose(fp);

	// fp = fopen("NURBS/control_point.dat", "w");
	// for (j = 0; j < Total_Control_Point_to_mesh[Total_mesh]; j++)
	// {
	// 	fprintf(fp, "%d:  %le  %le	 %le %le ", j, Node_Coordinate[j][0], Node_Coordinate[j][1], Displacement[j * DIMENSION + 0], Displacement[j * DIMENSION + 1]);
	// 	fprintf(fp, "\n");
	// }
	// fclose(fp);

	/*

	fp = fopen("new_zarusoba/control_point.dat", "w");
	fprintf(fp, "%d\n", Total_Control_Point);
	for (j = 0; j < Total_Control_Point; j++)
	{
		fprintf(fp, "%le  %le\n", Node_Coordinate[j][0], Node_Coordinate[j][1]);
	}
	fclose(fp);

	calculate_extendmesh_using_NURBS(element_emsh, Total_Element, Total_Control_Point);

	fp = fopen("new_zarusoba/extended_mesh.emsh", "w");
	fprintf(fp, "%d\n", real_Total_Element * 10);

	q = 0;
	r = 0;
	for (i = Total_Control_Point; i < (Total_Control_Point + No_points_for_new_zarusoba - 2); i = i + 2)
	{
		if (q != 10 + 11 * r)
		{
			fprintf(fp, "%d %d %d %d\n", i, i + 1, i + 3, i + 2);
			printf("q:%d\n", q);
		}
		if (q == 10 + 11 * r)
		{
			r++;
		}
		printf("r:%d\n", r);
		q++;
	}

	for (l = 0; l < No_Patch; l++)
	{
		Total_net += (No_Control_point[l][0] - 1) * (No_Control_point[l][1] - 1);
	}

	fprintf(fp, "%d\n", Total_net * 2);
	for (l = 0; l < No_Patch; l++)
	{
		for (j = 0; j < No_Control_point[l][1] - 1; j++)
		{
			for (i = 0; i < No_Control_point[l][0] - 1; i++)
			{
				fprintf(fp, "%d %d\n", Adress_Controlpoint[l][i][j], Adress_Controlpoint[l][i + 1][j]);
				fprintf(fp, "%d %d\n", Adress_Controlpoint[l][i + 1][j], Adress_Controlpoint[l][i + 1][j + 1]);
				fprintf(fp, "%d %d\n", Adress_Controlpoint[l][i + 1][j + 1], Adress_Controlpoint[l][i][j + 1]);
				fprintf(fp, "%d %d\n", Adress_Controlpoint[l][i][j + 1], Adress_Controlpoint[l][i][j]);
			}
		}
	}

	fprintf(fp, "%d\n", Total_Control_Point + No_points_for_new_zarusoba);

	for (i = 0; i < Total_Control_Point; i++)
	{
		fprintf(fp, "%lf %lf\n", Node_Coordinate[i][0], Node_Coordinate[i][1]);
	}

	for (i = Total_Control_Point; i < Total_Control_Point + No_points_for_new_zarusoba; i++)
	{
		fprintf(fp, "%lf %lf\n", data_result_shape_x_for_new_zarusoba[i], data_result_shape_y_for_new_zarusoba[i]);
	}

	fclose(fp);

	fp = fopen("new_zarusoba/radial_displacement.dat", "w");
	fprintf(fp, "label=Radial Displacement\nnum_items=%d\n\n", Total_Control_Point + No_points_for_new_zarusoba);
	for (i = 0; i < Total_Control_Point; i++)
	{
		fprintf(fp, "%d:0.0\n", i);
	}
	for (i = Total_Control_Point; i < (Total_Control_Point + No_points_for_new_zarusoba); i++)
		fprintf(fp, "%d:%.16e\n", i, pow((data_result_disp_x_for_new_zarusoba[i] * data_result_disp_x_for_new_zarusoba[i] + data_result_disp_y_for_new_zarusoba[i] * data_result_disp_y_for_new_zarusoba[i]), 0.5));
	fclose(fp);
	*/

	/*printf("Finish Make_Displacement\n");
	Make_Strain_Quad_4( E, nu, Total_Element, El_No,Total_Control_Point);
	printf("Finish Make_Strain\n");
	Make_Stress_2D( E, nu, Total_Element, DM);
	printf("Finish Make_Stress\n");
	Make_ReactionForce_Quad_4( Total_Element, Total_Control_Point, El_No);
	printf("Finish Make_ReactionForce\n");
	puts("sol_vec");
	Make_Parameter_z( Total_Element, E, nu, DM);
	printf("Finish Make_Parameter_z\n");
	Make_Output( Total_Control_Point, Total_Element );
	printf("Finish Make_Output\n" );*/

/*
	fp = fopen("colored_point/NURBS_points.txt", "w");

	fprintf(fp, "%d\n", No_points_for_colored_points);
	for (p = 0; p < No_points_for_colored_points; p++)
	{
		fprintf(fp, "%-.13lf  %-.13lf\n", data_result_shape_x[p], data_result_shape_y[p]);
	}

	fclose(fp);

	fp = fopen("colored_point/NURBS_disp_x.dat", "w");

	fprintf(fp, "label=Displacement\nnum_items=%d\n\n", No_points_for_colored_points);
	for (p = 0; p < No_points_for_colored_points; p++)
	{
		fprintf(fp, "%d:%-.13le\n", p, data_result_disp_x[p]);
	}

	fclose(fp);

	fp = fopen("colored_point/NURBS_disp_y.dat", "w");

	fprintf(fp, "label=Displacement\nnum_items=%d\n\n", No_points_for_colored_points);
	for (p = 0; p < No_points_for_colored_points; p++)
	{
		fprintf(fp, "%d:%-.13le\n", p, data_result_disp_y[p]);
	}

	fclose(fp);

	fp = fopen("colored_point/NURBS_disp_radius.dat", "w");

	fprintf(fp, "label=Displacement\nnum_items=%d\n\n", No_points_for_colored_points);
	for (p = 0; p < No_points_for_colored_points; p++)
	{
		fprintf(fp, "%d:%-.13le\n", p, pow((data_result_disp_x[p] * data_result_disp_x[p] + data_result_disp_y[p] * data_result_disp_y[p]), 0.5));
	}

	fclose(fp);
*/
	// fp = fopen("Gauss_stress/Gausspoint_coordinates.dat", "w");
	// fprintf(fp, "%d\n", real_Total_Element_to_mesh[Total_mesh] * POW_Ng);
	// Gausspoint_coordinate(Total_Element_to_mesh[Total_mesh],
	// 					  Total_Control_Point_to_mesh[Total_mesh]);
	// //for (i = 0; i < Total_Element; i++)
	// for (re = 0; re < real_Total_Element_to_mesh[Total_mesh]; re++)
	// {
	// 	i = real_element[re];
	// 	for (j = 0; j < POW_Ng; j++)
	// 	{
	// 		fprintf(fp, "%.16e %.16e\n", Gausspoint_coordinates[i][j][0], Gausspoint_coordinates[i][j][1]);
	// 	}
	// }
	// fclose(fp);

	// fp = fopen("Gauss_stress/Gauss_stress_x.dat", "w");
	// k = 0;
	// fprintf(fp, "label=Stress x\nnum_items=%d\n\n", real_Total_Element_to_mesh[Total_mesh] * POW_Ng);
	// //for (i = 0; i < Total_Element; i++)
	// for (re = 0; re < real_Total_Element_to_mesh[Total_mesh]; re++)
	// {
	// 	i = real_element[re];
	// 	for (j = 0; j < POW_Ng; j++)
	// 	{
	// 		fprintf(fp, "%d:%.16e\n", k, Stress[i][j][0]);
	// 		k++;
	// 	}
	// }

	// fclose(fp);

	// fp = fopen("Gauss_stress/Gauss_stress_y.dat", "w");
	// k = 0;
	// fprintf(fp, "label=Stress y\nnum_items=%d\n\n", real_Total_Element_to_mesh[Total_mesh] * POW_Ng);
	// //for (i = 0; i < Total_Element; i++)
	// for (re = 0; re < real_Total_Element_to_mesh[Total_mesh]; re++)
	// {
	// 	i = real_element[re];
	// 	for (j = 0; j < POW_Ng; j++)
	// 	{
	// 		fprintf(fp, "%d:%.16e\n", k, Stress[i][j][1]);
	// 		k++;
	// 	}
	// }

	// fclose(fp);

	// fp = fopen("Gauss_stress/Gauss_stress_xy.dat", "w");

	// fprintf(fp, "label=Stress xy\nnum_items=%d\n\n", real_Total_Element_to_mesh[Total_mesh] * POW_Ng);
	// k = 0;
	// //for (i = 0; i < Total_Element; i++)
	// for (re = 0; re < real_Total_Element_to_mesh[Total_mesh]; re++)
	// {
	// 	i = real_element[re];
	// 	for (j = 0; j < POW_Ng; j++)
	// 	{
	// 		fprintf(fp, "%d:%.16e\n", k, Stress[i][j][2]);
	// 		k++;
	// 	}
	// }

	// fclose(fp);

	if (Total_mesh == 1)
	{
		printf("start GP calc\n");
		// Calculation_at_GP(E, nu);
		printf("end GP calc\n");
	}


	// 重ね合わせをスキップしてJ積分を行う
	if (SKIP_S_IGA != 1)
	{
		//重ね合わせの結果
		division_ele_xi  = DIVISION_ELE_XI;
		division_ele_eta = DIVISION_ELE_ETA;

		if (division_ele_xi > MAX_DIVISION) {
			printf("Error!!\n");
			printf("Too many Divsion at xi!\n"
				"Maximum of division is %d (Now %d)\n"
				"\n", MAX_DIVISION, division_ele_xi);
			exit(1);
		}
		if (division_ele_eta > MAX_DIVISION) {
			printf("Error!!\n");
			printf("Too many Divsion at eta!\n"
				"Maximum of division is %d (Now %d)\n"
				"\n", MAX_DIVISION, division_ele_eta);
			exit(1);
		}

		//ローカルメッシュの情報取得
		GetLocData();

		int patch_n_loc = 0, patch_n_glo = 0;	//パッチ番号

		ReadFile();
		fp = fopen("view.dat", "w");
		fprintf(fp, "%d\t%d\t%d\n",
				fields_flag, division_ele_xi, division_ele_eta);
		fclose(fp);
		//machino
		fp = fopen("view_r_theta.dat", "w");
		fprintf(fp, "%d\t%d\t%d\n",
				fields_flag, division_ele_xi, division_ele_eta);
		fclose(fp);

		n_patch_glo = patch_n - n_patch_loc;

		//for s-IGA
		//重ね合わせ結果出力のためのoverlay_view.dat作成
		fp = fopen("overlay_view.dat", "w");
		fprintf(fp, "%d\t%d\t%d\n",
					fields_flag, division_ele_xi, division_ele_eta);
		fclose(fp);
		//machino
		fp = fopen("overlay_view_r_theta.dat", "w");
		fprintf(fp, "%d\t%d\t%d\n",
					fields_flag, division_ele_xi, division_ele_eta);
		fclose(fp);

		//グラフ作成のための出力
		fp = fopen("disp_graph.txt", "w");
		fprintf(fp, "patch_n\tx\ty\tdisp_x\tdisp_y\n");
		fclose(fp);

		fp = fopen("stress_y_graph.txt", "w");
		fprintf(fp, "patch_n\tx\ty\tstress_yy\n");
		fclose(fp);

		fp = fopen("stress_y_graph_0.txt", "w");
		fprintf(fp, "patch_n\tx\ty\tstress_yy\n");
		fclose(fp);

		fp = fopen("stress_vm_graph.txt", "w");
		fprintf(fp, "xi\teta\tx\ty\tstress_vm\n");
		fclose(fp);

		fp = fopen("over_disp_graph.txt", "w");
		fprintf(fp, "patch_n\tx\ty\tdisp_x\tdisp_y\n");
		fclose(fp);

		fp = fopen("over_stress_x_graph.txt", "w");
		fprintf(fp, "x\ty\tstress_xx\n");
		fclose(fp);

		fp = fopen("over_stress_y_graph.txt", "w");
		fprintf(fp, "x\ty\tstress_yy\n");
		fclose(fp);

		fp = fopen("over_stress_y_graph_0.txt", "w");
		fprintf(fp, "x\ty\tstress_yy\n");
		fclose(fp);

		fp = fopen("over_stress_r_theta_graph.txt", "w");
		fprintf(fp, "xi\teta\tx\ty\tstress_r\tstress_theta\n");
		fclose(fp);

		fp = fopen("over_stress_vm_graph.txt", "w");
		fprintf(fp, "xi\teta\tx\ty\tstress_vm\n");
		fclose(fp);

		// S-IGAでローカル上のガウス点で重ね合わせた値をデータ整理する場合に使う
		if (Total_mesh >= 2)
		{
			printf("start GP calc\n");

			// patch_n_glo = 0;
			// Calculation_overlay_at_GP(E, nu,
			// 						  order_xi[patch_n_glo],order_eta[patch_n_glo],
			// 						  knot_n_xi[patch_n_glo], knot_n_eta[patch_n_glo],
			// 						  cntl_p_n_xi[patch_n_glo], cntl_p_n_eta[patch_n_glo],
			// 						  knot_vec_xi[patch_n_glo], knot_vec_eta[patch_n_glo],
			// 						  cntl_px[patch_n_glo], cntl_py[patch_n_glo],
			// 						  disp_cntl_px[patch_n_glo], disp_cntl_py[patch_n_glo],
			// 						  weight[patch_n_glo]);

			printf("end GP calc\n");
		}

		// 重ね合わせ
		for (i = 0; i < patch_n; i++) {

			fp = fopen("stress_vm_graph.txt", "a");
			fprintf(fp, "\npatch_n;%d\n\n",i);
			fclose(fp);

			graph_patch_n = i;

			printf("----------Start calculation at patch %d----------\n\n", i);
			Calculation(order_xi[i], order_eta[i],
						knot_n_xi[i], knot_n_eta[i],
						cntl_p_n_xi[i], cntl_p_n_eta[i],
						knot_vec_xi[i], knot_vec_eta[i],
						cntl_px[i], cntl_py[i],
						disp_cntl_px[i], disp_cntl_py[i],
						weight[i]);
			printf("-----------End calculation at patch %d-----------\n\n", i);

			if (i >= n_patch_glo)	//ローカル上のパッチに対しては重合計算行う
			{
				patch_n_loc = i;
				printf("----------Start overlay calculation at patch %d in LOCAL patch----------\n\n", i);
				for (j = 0; j < n_patch_glo; j++)
				{
					patch_n_glo = j;
					//printf("patch_n_loc: %d \tpatch_n_glo: %d\n",
					//		  patch_n_loc, patch_n_glo);
					Calculation_overlay(order_xi[patch_n_loc],order_eta[patch_n_loc],
										knot_n_xi[patch_n_loc], knot_n_eta[patch_n_loc],
										cntl_p_n_xi[patch_n_loc], cntl_p_n_eta[patch_n_loc],
										knot_vec_xi[patch_n_loc], knot_vec_eta[patch_n_loc],
										cntl_px[patch_n_loc], cntl_py[patch_n_loc],
										weight[patch_n_loc],
										order_xi[patch_n_glo],order_eta[patch_n_glo],
										cntl_p_n_xi[patch_n_glo], cntl_p_n_eta[patch_n_glo],
										knot_vec_xi[patch_n_glo], knot_vec_eta[patch_n_glo],
										cntl_px[patch_n_glo], cntl_py[patch_n_glo],
										disp_cntl_px[patch_n_glo], disp_cntl_py[patch_n_glo],
										weight[patch_n_glo]);
				}
			}
		}
		printf("End S-IGA\n\n");
	}
	else if (SKIP_S_IGA == 1)
	{
		printf("SKIP_S_IGA = 1, skip S-IGA\n");
	}

	if (SKIP_S_IGA == 2) {
		printf("SKIP_S_IGA = 2, exit before J integration\n");
		exit(0);
	}

	// 重ね合わせをスキップした場合ここから
	printf("Start J Integration Mixed Mode\n\n");

	Make_Displacement_grad(El_No);
	printf("Finish Make_Displacement_grad\n");
	Make_StrainEnergyDensity_2D();
	printf("Finish Make_StrainEnergyDensity\n");

	//globalでのDisp_gradを求める
	int e, N;
	Make_Displacement_grad_glo_check(real_Total_Element_to_mesh[Total_mesh]);
	
	Make_gauss_array(0);

	for(re = 0; re < real_Total_Element_to_mesh[Total_mesh]; re++){
        e = real_element[re];
		for(N = 0; N < GP_2D; N++){
			//Disp_gradをlocalとglobalで重ね合わせる
			for(i = 0; i <  DIMENSION * DIMENSION; i++){
				Disp_grad_overlay[e][N][i] = Disp_grad[e][N][i] + Disp_grad_glo[e][N][i];
				// printf("Disp_grad[%d][%d][%d] = %.10e\tDisp_grad_glo[%d][%d][%d] = %.10e\n", e, N, i, Disp_grad[e][N][i], e, N, i, Disp_grad_glo[e][N][i]);
				printf("Disp_grad_overlay[%d][%d][%d] = %.10e\n", e, N, i, Disp_grad_overlay[e][N][i]);
			}
			//Strainをlocalとglobalで重ね合わせる
			for (i = 0; i < D_MATRIX_SIZE; i++){
				Strain_overlay[e][N][i] = Strain[e][N][i] + Strain_glo[e][N][i];
				printf("Strain_overlay[%d][%d][%d] = %.10e\n", e, N, i, Strain_overlay[e][N][i]);
			}
		}
	}
	Make_Stress_2D_glo(E, nu, Total_Element_to_mesh[Total_mesh], DM);

	Make_gauss_array(0);

	for(re = 0; re < real_Total_Element_to_mesh[Total_mesh]; re++){
		e = real_element[re];
		for(N = 0; N < GP_2D; N++){
			for (i = 0; i < D_MATRIX_SIZE; i++){
				Stress_overlay[e][N][i] = Stress[e][N][i] + Stress_glo[e][N][i];
				// printf("Stress_overlay[%d][%d][%d] = %.10e\n", e, N, i, Stress_overlay[e][N][i]);
			}
		}
	}

	// Make_Displacement_grad_overlay(E, nu, Total_Element_to_mesh[Total_mesh] , El_No, Total_Control_Point_to_mesh[Total_mesh]);
	// printf("Finish Make_Displacement_grad\n");
	Make_StrainEnergyDensity_2D_overlay();
	printf("Finish Make_StrainEnergyDensity\n");
	Make_Parameter_z_overlay(Total_Element_to_mesh[Total_mesh], E, nu, DM);
	printf("Finish Make_Parameter_z_overlay\n");

	/* For the J-integral Evaluation */
	int Location_Crack_Tip_Patch;
	double Location_Local_Coordinates[DIMENSION];
	double Virtual_Crack_Extension_Ct_Pt[MAX_N_NODE][DIMENSION];
	double /*J_value,*/ DeltaA;
	J_Integral_Input_Data(Total_Control_Point_to_mesh[Total_mesh],&Location_Crack_Tip_Patch,Location_Local_Coordinates,Virtual_Crack_Extension_Ct_Pt, &DeltaA);

	printf("Start J_Integral_Computation\n");
	// J_value = J_Integral_Computation(Total_Control_Point_to_mesh[Total_mesh], Total_Element_to_mesh[Total_mesh], Location_Crack_Tip_Patch,  Location_Local_Coordinates, Virtual_Crack_Extension_Ct_Pt, DeltaA, El_No);

	double /*Kfact1,*/ Kfact1_ref, Kfact2_ref, ref_Sigma = 1.0, ref_a = 2.5, ref_W = 50.0, Pi = 4.0 * atan(1.0), rad = Pi / 4.0;
	printf("Pi = %.15e E = %.15e  nu = %.15e\n", Pi, E, nu);
	printf("o = %1.8e\ta = %1.8e\tW = %1.8e\n", ref_Sigma, ref_a, ref_W);
	Kfact1_ref = ref_Sigma * sin(rad) * sin(rad) * sqrt(Pi * ref_a);
	Kfact2_ref = ref_Sigma * sin(rad) * cos(rad) * sqrt(Pi * ref_a);

	// Kfact1 = sqrt(2.0 * J_value * E / (1 - nu*nu));

	J_Integral_Computation_Interaction(Total_Control_Point_to_mesh[Total_mesh], Location_Local_Coordinates, Virtual_Crack_Extension_Ct_Pt, DeltaA, E, nu, DM);

	double K_aux_mode1, K_aux_mode2;
	K_aux_mode1 = sqrt(J_integral_value_aux_mode1 * E / (1.0 - nu * nu));
	K_aux_mode2 = sqrt(J_integral_value_aux_mode2 * E / (1.0 - nu * nu));

	//for interaction integral（モデルの対称性によりかける倍率を変える。 1/4モデル → * 2.0、　フルモデル → * 1.0）
	printf("J_value_aux_mode1 = %1.10e\n", J_integral_value_aux_mode1);
	printf("J_value_aux_mode2 = %1.10e\n", J_integral_value_aux_mode2);

	printf("K_aux_mode1 = %.15e\n", K_aux_mode1);
	printf("K_aux_mode2 = %.15e\n", K_aux_mode2);

	printf("**** J-Integral Value 3 by 3 Gauss Integration ***\n");
	printf("K_mode1 = %.15e\n", K_mode1);
	printf("K_mode2 = %.15e\n", K_mode2);
	// printf("J-integral value = %20.10e\n",J_value);
	// printf("K_IGA = %.15e K_Ref = %.15e\n",Kfact1, Kfact1_ref);

	//誤差評価
	printf("Kfact1_ref = %.15e\n", Kfact1_ref);
	printf("Kfact2_ref = %.15e\n", Kfact2_ref);
	// printf("K_error(J_integral) = %.15e\n",(Kfact1 - Kfact1_ref) / Kfact1_ref * 100);
	printf("K_error(I_integral_mode1) = %.15e\n", (K_mode1 - Kfact1_ref) / Kfact1_ref * 100.0);
	printf("K_error(I_integral_mode2) = %.15e\n", (K_mode2 - Kfact2_ref) / Kfact2_ref * 100.0);
	printf("****                                           ***\n");

    //     J_value = J_Integral_Computation6by6(Total_Control_Point_to_mesh[Total_mesh], Total_Element_to_mesh[Total_mesh],Location_Crack_Tip_Patch,  Location_Local_Coordinates, Virtual_Crack_Extension_Ct_Pt, DeltaA, El_No, E, nu, DM);

    //     printf("Pi = %.15e E = %.15e  nu = %.15e\n",Pi,E,nu);
    //     Kfact1_ref = ref_Sigma * sqrt(Pi * ref_a) * sqrt(1.0/cos(alpha_ref * Pi / 2.0))* (1.0 - 0.025 * pow(alpha_ref, 2.0) + 0.06 * pow(alpha_ref, 4.0));

    //     Kfact1 = sqrt(2.0* J_value * E / (1 - nu*nu));
	
	// printf("**** J-Integral Value 6 by 6 Gauss Integration ***\n");
	// printf("J-integral value = %20.10e\n",J_value);
    //     printf("K_IGA(6by6 Integration) = %.15e K_Ref = %.15e\n",Kfact1, Kfact1_ref);
	// printf("****                                           ***\n");

	t1 = clock();
	printf("Total calculation time:%.2f[s]\n",(double)(t1-start)/CLOCKS_PER_SEC);

	return 0;
}


//ファイルからデータをもらう
void J_Integral_Input_Data(int Total_Control_Point, int *Location_Crack_Tip_Patch, double Location_Local_Coordinates[DIMENSION], double Virtual_Crack_Extension_Ct_Pt[MAX_N_NODE][DIMENSION], double *DeltaA)
{
	// FILE *fp;
	int ii, jj, kk;
	int Num_Non_Zero_Pt;
	int kk_control_pt;

	fp = fopen("J-int.inp","r");

	for(jj = 0; jj < DIMENSION; jj++) 
		for(ii = 0; ii < Total_Control_Point; ii++) 
			Virtual_Crack_Extension_Ct_Pt[ii][jj] = 0.0;

	fscanf(fp,"%d %lf %lf %lf", Location_Crack_Tip_Patch, &Location_Local_Coordinates[0], &Location_Local_Coordinates[1], DeltaA); 
	printf("J-integral input: ID of patch (crack tip) %d Normalize coords %f %f %f\n", *Location_Crack_Tip_Patch,Location_Local_Coordinates[0], Location_Local_Coordinates[1], *DeltaA);
	fscanf(fp, "%d", &Num_Non_Zero_Pt);
	printf("%d\n", Num_Non_Zero_Pt);
	for(kk=0; kk < Num_Non_Zero_Pt; kk++){
		fscanf(fp, "%d", &kk_control_pt); 
			for(jj=0; jj < DIMENSION; jj++) 
				fscanf(fp, "%lf", &Virtual_Crack_Extension_Ct_Pt[kk_control_pt][jj]);
		printf("kk = %d  kk_control_pt = %d  cood values %f %f\n", kk, kk_control_pt, Virtual_Crack_Extension_Ct_Pt[kk_control_pt][0], Virtual_Crack_Extension_Ct_Pt[kk_control_pt][1]);
	}
	fclose(fp);
}


double J_Integral_Computation(int Total_Control_Point, double Virtual_Crack_Extension_Ct_Pt[MAX_N_NODE][DIMENSION], double DeltaA)
{
	int ii, jj, kk; 
	int e, i, j, re, N;
	double J;
	double J_integral_value = 0.0;
	double q_func_grad[DIMENSION][DIMENSION], EMT[DIMENSION][DIMENSION];
	double Virtual_Crack_Extension_El_Nodes[MAX_NO_CCpoint_ON_ELEMENT * DIMENSION];
	double GXY[9][2], qXY[2];
	double check_int_r = 0.0;

	Make_gauss_array(0);

	/*
	ii = 0;
	for(i=0; i < 6; i++) 
		for(j=0; j < 6; j++){
			G2D36[ii][0] = G1D6[j]; G2D36[ii][2] = G1D6[i];
			w2D36[ii] = w1D6[j] * w1D6[i];
			ii++;
			}
	*/			


	for( re = 0; re < real_Total_Element_to_mesh[Total_mesh]; re++ ) {
		e = real_element[re];
		/* printf("re = %d  e = %d\n",re,e);*/
		double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION];
		double b_grad[DIMENSION * DIMENSION][2 * MAX_NO_CCpoint_ON_ELEMENT];
		for (ii = 0; ii < DIMENSION; ii++) for (jj = 0; jj < DIMENSION; jj++) q_func_grad[ii][jj] = 0.0;

		//printf("%d,No_Control_point_ON_ELEMENT[%d]:%d\n",e,Element_patch[e], No_Control_point_ON_ELEMENT[Element_patch[e]]);
		//Bマトリックスを取得
		for( i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[e]]; i++ ){
			for( j = 0; j < DIMENSION; j++ ){
				X[i][j] = Node_Coordinate[ Controlpoint_of_Element[e][i] ][j];
				Virtual_Crack_Extension_El_Nodes[i * DIMENSION + j] = Virtual_Crack_Extension_Ct_Pt[Controlpoint_of_Element[e][i]][j];
			}	
			//printf("3^3-2check re = %d  Virtual_Crack_Extension_El_Nodes[2 * %d] = %lf   Virtual_Crack_Extension_Ct_Pt = %lf\n",re, i, Virtual_Crack_Extension_El_Nodes[2 * i], Virtual_Crack_Extension_Ct_Pt[Controlpoint_of_Element[e][i]][0]);
			//printf("Node_Coordinate[%d][%d]:%le\n",Controlpoint_of_Element[e][i],j, Node_Coordinate[ Controlpoint_of_Element[e][i] ][j]);
			//printf("X[%d][%d]:%le\n",Controlpoint_of_Element[e][i],j,X[i][j] );
		}
		double zero_thresh = 1.0e-12;
		int zero_flag = 0;
		for(kk = 0; kk < 2 * No_Control_point_ON_ELEMENT[Element_patch[e]]; kk++)
			if(fabs(Virtual_Crack_Extension_El_Nodes[kk]) < zero_thresh) zero_flag++;
		if(zero_flag < 2 * No_Control_point_ON_ELEMENT[Element_patch[e]]) {
			printf("chekc2 re = %d  e = %d\n", re, e);
			for(ii = 0; ii < No_Control_point_ON_ELEMENT[Element_patch[e]]; ii++) 
				printf("Virtual Crack Extension %.15e %.15e\n",Virtual_Crack_Extension_El_Nodes[ii * DIMENSION],Virtual_Crack_Extension_El_Nodes[ii * DIMENSION + 1]);
			double CheckArea = 0.0, CheckArea_Norm=0.0;
			for( N = 0; N < GP_2D; N++ ) {
				GXY[N][0] = 0.0; GXY[N][1] = 0.0; qXY[0] = 0.0; qXY[1] = 0.0;

				for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[e]]; i++) {
					double R_shape_func = Shape_func(i, Gxi[N],e);
					for (j = 0; j < DIMENSION; j++) {
						GXY[N][j] += R_shape_func * X[i][j];
						qXY[j] += R_shape_func * Virtual_Crack_Extension_El_Nodes[i * DIMENSION + j];
					}
				}
				//printf("3^3chekc3 re = %d  e = %d N = %d  qXY = %lf\n",re,e,N, qXY[0]);
				for (ii=0; ii < DIMENSION; ii++) for (jj=0; jj < DIMENSION; jj++) {
					q_func_grad[ii][jj] = 0.0;
					EMT[ii][jj] = 0.0;
					}
				printf("chekc4 re = %d  e = %d N = %d\n",re,e,N);
				Make_b_grad_Matrix( e, b_grad, Gxi[N], X, &J);
				printf("check5 re = %d  e = %d N = %d J = %.15e\n",re,e,N,J);
				for( ii = 0; ii < DIMENSION; ii++)
					for(jj = 0; jj < DIMENSION; jj++)			
						for( j = 0; j < 2 * No_Control_point_ON_ELEMENT[Element_patch[e]]; j++ )
							q_func_grad[ii][jj] += b_grad[ii*DIMENSION + jj][j] * Virtual_Crack_Extension_El_Nodes[j];

				printf("chekc6 re = %d  e = %d N = %d\n", re, e, N);
				/* Setting up the Energy Momentum Tensor */	
				printf("check 9 %.15e  %.15e  %.15e  %.15e  %.15e\n",StrainEnergyDensity_overlay[e][N],Stress_overlay[e][N][0],Disp_grad_overlay[e][N][0],Stress_overlay[e][N][2],Disp_grad_overlay[e][N][2]);
				EMT[0][0] =  StrainEnergyDensity_overlay[e][N] - (Stress_overlay[e][N][0] * Disp_grad_overlay[e][N][0] + Stress_overlay[e][N][2] * Disp_grad_overlay[e][N][2]);
				printf("G E: %.15e %.15e\n", Gxi[N][0], Gxi[N][1]);
				printf("EMT[0][0] = %.15e\n",EMT[0][0]);
				/* W - s11 * du1/dx1 - s12 * du2/dx1 */
				EMT[0][1] =                            - (Stress_overlay[e][N][0] * Disp_grad_overlay[e][N][1] + Stress_overlay[e][N][2] * Disp_grad_overlay[e][N][3]);
				printf("EMT[0][1] = %.15e\n",EMT[0][1]);
				/* - s11 * du1/dx2 - s22 * du2/dx2 */
				EMT[1][0] =                            - (Stress_overlay[e][N][2] * Disp_grad_overlay[e][N][0] + Stress_overlay[e][N][1] * Disp_grad_overlay[e][N][2]);
				printf("EMT[1][0] = %.15e\n",EMT[1][0]);
				/* - s21 * du1/dx1 - s22 * du2/dx1 */
				EMT[1][1] =  StrainEnergyDensity_overlay[e][N] - (Stress_overlay[e][N][2] * Disp_grad_overlay[e][N][1] + Stress_overlay[e][N][1] * Disp_grad_overlay[e][N][3]);
				printf("EMT[1][1] = %.15e\n",EMT[1][1]);
				/* W - s21 * du1/dx2 - s22 * du2/dx2 */

				printf("chekc8 re = %d  e = %d N = %d\n", re, e, N);

				printf("check strain energry e N %d %d %.15e\n", e, N, StrainEnergyDensity_overlay[e][N]);
				printf("check Stresses e N %d %d %.15e %.15e %.15e\n", e, N, Stress_overlay[e][N][0], Stress_overlay[e][N][1], Stress_overlay[e][N][2]);
				printf("check Strain e N %d %d %.15e %.15e %.15e\n", e, N, Strain_overlay[e][N][0], Strain_overlay[e][N][1], Strain_overlay[e][N][2]);
				printf("check Dispgrad e N %d %d %.15e %.15e %.15e %.15e\n", e, N, Disp_grad_overlay[e][N][0], Disp_grad_overlay[e][N][1], Disp_grad_overlay[e][N][2], Disp_grad_overlay[e][N][3]);

				printf("check q_func_grad e N %d %d %.15e %.15e %.15e %.15e\n",e, N, q_func_grad[0][0], q_func_grad[0][1],q_func_grad[1][0], q_func_grad[1][1]);
				/* q_func_grad[0][0] = -(GXY[N][0] - 25.0) / sqrt((GXY[N][0] - 25.0)*(GXY[N][0] - 25.0)+GXY[N][1]*GXY[N][1])/5.0;
								q_func_grad[0][1] = -(GXY[N][1] - 0.0) / sqrt((GXY[N][0] - 25.0)*(GXY[N][0] - 25.0)+GXY[N][1]*GXY[N][1])/5.0;

								printf("check q_func_grad2 e N %d %d %.15e %.15e %.15e %.15e\n",e, N, q_func_grad[0][0], q_func_grad[0][1],q_func_grad[1][0], q_func_grad[1][1]); */
				double rrrr = sqrt((GXY[N][0] - 25.0)*(GXY[N][0] - 25.0) +(GXY[N][1] - 0.0)*(GXY[N][1] - 0.0));
				J_integral_value -= (EMT[0][0] * q_func_grad[0][0] + EMT[1][0] * q_func_grad[0][1]) * w[N] * J;
				/* EMT11 * dq1/dq1 + EMT21 * dq1/dx2 */
				J_integral_value -= (EMT[0][1] * q_func_grad[1][0] + EMT[1][1] * q_func_grad[1][1]) * w[N] * J; 
				/* EMT12 * dq2/dq1 + EMT22 * dq2/dx2 */ 
				CheckArea += J*w[N];
				CheckArea_Norm += w[N];
				double NPP = 2.0;
				check_int_r += 1.0 / rrrr / pow(5.0, NPP-1) * pow(5.0 - rrrr, NPP-1) * J * w[N];
				/*check_int_r += 1.0/rrrr * J * w[N];*/

				printf("check7 re = %d  e = %d N = %d J = %.15e\n",re,e,N,J);
			}
			printf("Area = %.15e\n", CheckArea);
			printf("Area Norm = %.15e\n", CheckArea_Norm);
			/* printf("check_int_r = %.15e 1/ 5*Pi = %.15e\n",check_int_r,4.0*atan(1.0)*5.0/5.0);*/
			printf("check_int_r = %.15e 5*Pi = %.15e\n", check_int_r, 4.0 * atan(1.0) * 5.0 / 2);
		}
	}
	printf("check_int_r = %.15e  5*Pi = %.15e\n", check_int_r, 4.0 * atan(1.0) * 5.0 / 2);
	fp = fopen("Virtual_Crack_Extension.dat", "w");
	fprintf(fp,"label=Virtual_Crack_Extension\n");
	fprintf(fp, "num_items=%d\n", Total_Control_Point);
	fprintf( fp, "\n");
	for(j = 0; j < Total_Control_Point; j++) {
		fprintf(fp, "%d:	%le %le ", j, Virtual_Crack_Extension_El_Nodes[j * DIMENSION + 0], Virtual_Crack_Extension_El_Nodes[j * DIMENSION + 1]);
		fprintf(fp, "\n");
		}
	fclose(fp);
	return(J_integral_value / DeltaA);
}


double J_Integral_Computation_Interaction(int Total_Control_Point, double Location_Local_Coordinates[DIMENSION], double Virtual_Crack_Extension_Ct_Pt[MAX_N_NODE][DIMENSION], double DeltaA, double E, double nu, int DM)
{
	int ii, jj, kk; 
	int e, i, j, re, N;
	double unit_basis_local[DIMENSION] = {0.0};
	double r_tip = sqrt(Location_Local_Coordinates[0] * Location_Local_Coordinates[0] + Location_Local_Coordinates[1] * Location_Local_Coordinates[1]);
	double J;
	double K_1 = 0.0;
	double K_2 = 0.0;
	double q_func_grad[DIMENSION][DIMENSION], q_func_grad_local[DIMENSION][DIMENSION], EMT_mode1[DIMENSION][DIMENSION], EMT_mode2[DIMENSION][DIMENSION], EMT_aux_mode1[DIMENSION][DIMENSION], EMT_aux_mode2[DIMENSION][DIMENSION];
	double Virtual_Crack_Extension_El_Nodes[MAX_NO_CCpoint_ON_ELEMENT * DIMENSION];
	double GXY[9][2], qXY[2];
	double check_int_r = 0.0;

	Make_gauss_array(0);

	//x'-y'（き裂先端）座標における単位基底ベクトル
	unit_basis_local[0] = Location_Local_Coordinates[0] / r_tip;
	unit_basis_local[1] = Location_Local_Coordinates[1] / r_tip;
	printf("unit_basis[0] : % 1.8e\n", unit_basis_local[0]);
	printf("unit_basis[1] : % 1.8e\n", unit_basis_local[1]);

	//2 × 2 の行列として局所座標系を登録
	T[0][0] = unit_basis_local[0];   T[0][1] = unit_basis_local[1];
	T[1][0] = -unit_basis_local[1];  T[1][1] = unit_basis_local[0];
	printf("T[0][0] = %1.8e\tT[0][1] = %1.8e\tT[1][0] = %1.8e\tT[1][1] = %1.8e\n", T[0][0], T[0][1], T[1][0], T[1][1]);

	for(re = 0; re < real_Total_Element_to_mesh[Total_mesh]; re++){
		e = real_element[re];
		/* printf("re = %d  e = %d\n",re,e);*/
		double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION];
		double b_grad[DIMENSION * DIMENSION][2 * MAX_NO_CCpoint_ON_ELEMENT];
		for(ii=0; ii < DIMENSION; ii++)
			for(jj=0; jj < DIMENSION; jj++){
				q_func_grad_local[ii][jj] = 0.0;
				q_func_grad[ii][jj] = 0.0;
			}

		//printf("%d,No_Control_point_ON_ELEMENT[%d]:%d\n",e,Element_patch[e], No_Control_point_ON_ELEMENT[Element_patch[e]]);
		//Bマトリックスを取得
		for(i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[e]]; i++){
			for( j = 0; j < DIMENSION; j++ ){
				X[i][j] = Node_Coordinate[Controlpoint_of_Element[e][i]][j];
				Virtual_Crack_Extension_El_Nodes[i * DIMENSION + j] = Virtual_Crack_Extension_Ct_Pt[Controlpoint_of_Element[e][i]][j];
			}	
			//printf("3^3-2check re = %d  Virtual_Crack_Extension_El_Nodes[2 * %d] = %lf   Virtual_Crack_Extension_Ct_Pt = %lf\n",re, i, Virtual_Crack_Extension_El_Nodes[2 * i], Virtual_Crack_Extension_Ct_Pt[Controlpoint_of_Element[e][i]][0]);
			//printf("Node_Coordinate[%d][%d]:%le\n",Controlpoint_of_Element[e][i],j, Node_Coordinate[ Controlpoint_of_Element[e][i] ][j]);
			//printf("X[%d][%d]:%le\n",Controlpoint_of_Element[e][i],j,X[i][j]);
		}

		Make_auxiliary_mode1(e, E, nu, DM, X, Location_Local_Coordinates[0], Location_Local_Coordinates[1]);
		Make_auxiliary_mode2(e, E, nu, DM, X, Location_Local_Coordinates[0], Location_Local_Coordinates[1]);

		Make_gauss_array(0);

		double zero_thresh = 1.0e-12;
		int zero_flag = 0;
		for(kk = 0; kk < 2 * No_Control_point_ON_ELEMENT[Element_patch[e]]; kk++)
			if(fabs(Virtual_Crack_Extension_El_Nodes[kk]) < zero_thresh) zero_flag++;
		if(zero_flag < 2 * No_Control_point_ON_ELEMENT[Element_patch[e]]){
			printf("chekc2 re = %d  e = %d\n", re, e);
			for(ii = 0; ii < No_Control_point_ON_ELEMENT[Element_patch[e]]; ii++) 
				printf("Virtual Crack Extension %.15e %.15e\n", Virtual_Crack_Extension_El_Nodes[ii * DIMENSION], Virtual_Crack_Extension_El_Nodes[ii * DIMENSION + 1]);
			double CheckArea = 0.0, CheckArea_Norm = 0.0;
			for(N = 0; N < GP_2D; N++){
				GXY[N][0] = 0.0; GXY[N][1] = 0.0; qXY[0] = 0.0; qXY[1] = 0.0;
				for(i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[e]]; i++){
					double R_shape_func = Shape_func(i, Gxi[N],e);
					for(j = 0; j < DIMENSION; j++){
						GXY[N][j] += R_shape_func * X[i][j];
						qXY[j] += R_shape_func * Virtual_Crack_Extension_El_Nodes[i * DIMENSION + j];
					}
				}
				//printf("3^3chekc3 re = %d  e = %d N = %d  qXY = %lf\n",re,e,N, qXY[0]);
				for(ii = 0; ii < DIMENSION; ii++)
					for(jj = 0; jj < DIMENSION; jj++){
						q_func_grad_local[ii][jj] = 0.0;
						q_func_grad[ii][jj] = 0.0;
						EMT_mode1[ii][jj] = 0.0;
						EMT_mode2[ii][jj] = 0.0;
						EMT_aux_mode1[ii][jj] = 0.0;
						EMT_aux_mode2[ii][jj] = 0.0;
					}
				// printf("chekc4 re = %d  e = %d N = %d\n",re,e,N);
				Make_b_grad_Matrix(e, b_grad, Gxi[N], X, &J);
				printf("check5 re = %d  e = %d N = %d J = %.15e\n", re, e, N, J);
				for(ii = 0; ii < DIMENSION; ii++)
					for(jj = 0; jj < DIMENSION; jj++)	
						for(j = 0; j < 2 * No_Control_point_ON_ELEMENT[Element_patch[e]]; j++)
							q_func_grad_local[ii][jj] += b_grad[ii * DIMENSION + jj][j] * Virtual_Crack_Extension_El_Nodes[j];
				
				//q関数をx'-y'座標からx-y座標に変換する
				q_func_grad[0][0] = T[0][0] * T[0][0] * q_func_grad_local[0][0] + T[0][0] * T[1][0] * (q_func_grad_local[0][1] + q_func_grad_local[1][0]) + T[1][0] * T[1][0] * q_func_grad_local[1][1];
				q_func_grad[0][1] = T[0][0] * T[0][1] * q_func_grad_local[0][0] + T[0][0] * T[1][1] * q_func_grad_local[0][1] + T[0][1] * T[1][0] * q_func_grad_local[1][0] + T[1][0] * T[1][1] * q_func_grad_local[1][1];
				q_func_grad[1][0] = T[0][0] * T[0][1] * q_func_grad_local[0][0] + T[0][0] * T[1][1] * q_func_grad_local[1][0] + T[0][1] * T[1][0] * q_func_grad_local[0][1] + T[1][0] * T[1][1] * q_func_grad_local[1][1];
				q_func_grad[1][1] = T[0][1] * T[0][1] * q_func_grad_local[0][0] + T[0][1] * T[1][1] * (q_func_grad_local[0][1] + q_func_grad_local[1][0]) + T[1][1] * T[1][1] * q_func_grad_local[1][1];
				printf("q_func_grad[0][0] = %1.10e\n", q_func_grad[0][0]);
				printf("q_func_grad[0][1] = %1.10e\n", q_func_grad[0][1]);
				printf("q_func_grad[1][0] = %1.10e\n", q_func_grad[1][0]);
				printf("q_func_grad[1][1] = %1.10e\n", q_func_grad[1][1]);

				// printf("chekc6 re = %d  e = %d N = %d\n",re,e,N);

				/* Setting up the Energy Momentum Tensor */	
				//モード1のEMTの計算
				printf("check 9 %.15e  %.15e  %.15e  %.15e  %.15e\n", 
						StrainEnergyDensity_aux_mode1[e][N], Stress_aux_mode1[e][N][0], Disp_grad_aux_mode1[e][N][0], Stress_aux_mode1[e][N][2], Disp_grad_aux_mode1[e][N][2]);
				EMT_mode1[0][0] = StrainEnergyDensity_aux_mode1[e][N] - (Stress_overlay[e][N][0] * Disp_grad_aux_mode1[e][N][0] + Stress_aux_mode1[e][N][0] * Disp_grad_overlay[e][N][0] + Stress_overlay[e][N][2] * Disp_grad_aux_mode1[e][N][2] + Stress_aux_mode1[e][N][2] * Disp_grad_overlay[e][N][2]);
				printf("G E: %.15e %.15e\n", Gxi[N][0], Gxi[N][1]);
				printf("EMT_mode1[0][0] = %.15e\n", EMT_mode1[0][0]);
				/* W - s11 * du1/dx1 - s12 * du2/dx1 */
				EMT_mode1[0][1] =                            - (Stress_overlay[e][N][0] * Disp_grad_aux_mode1[e][N][1] + Stress_aux_mode1[e][N][0] * Disp_grad_overlay[e][N][1] + Stress_overlay[e][N][2] * Disp_grad_aux_mode1[e][N][3] + Stress_aux_mode1[e][N][2] * Disp_grad_overlay[e][N][3]);
				printf("EMT_mode1[0][1] = %.15e\n", EMT_mode1[0][1]);
				/* - s11 * du1/dx2 - s22 * du2/dx2 */
				EMT_mode1[1][0] =                            - (Stress_overlay[e][N][2] * Disp_grad_aux_mode1[e][N][0] + Stress_aux_mode1[e][N][2] * Disp_grad_overlay[e][N][0] + Stress_overlay[e][N][1] * Disp_grad_aux_mode1[e][N][2] + Stress_aux_mode1[e][N][1] * Disp_grad_overlay[e][N][2]);
				printf("EMT_mode1[1][0] = %.15e\n", EMT_mode1[1][0]);
				/* - s21 * du1/dx1 - s22 * du2/dx1 */
				EMT_mode1[1][1] = StrainEnergyDensity_aux_mode1[e][N] - (Stress_overlay[e][N][2] * Disp_grad_aux_mode1[e][N][1] + Stress_aux_mode1[e][N][2] * Disp_grad_overlay[e][N][1] + Stress_overlay[e][N][1] * Disp_grad_aux_mode1[e][N][3] + Stress_aux_mode1[e][N][1] * Disp_grad_overlay[e][N][3]);
				printf("EMT_mode1[1][1] = %.15e\n", EMT_mode1[1][1]);
				/* W - s21 * du1/dx2 - s22 * du2/dx2 */

				//モード2のEMTの計算
				// printf("check 9 %.15e  %.15e  %.15e  %.15e  %.15e\n",StrainEnergyDensity_aux_mode2[e][N],Stress_aux_mode2[e][N][0],Disp_grad_aux_mode2[e][N][0],Stress_aux_mode2[e][N][2],Disp_grad_aux_mode2[e][N][2]);
				EMT_mode2[0][0] = StrainEnergyDensity_aux_mode2[e][N] - (Stress_overlay[e][N][0] * Disp_grad_aux_mode2[e][N][0] + Stress_aux_mode2[e][N][0] * Disp_grad_overlay[e][N][0] + Stress_overlay[e][N][2] * Disp_grad_aux_mode2[e][N][2] + Stress_aux_mode2[e][N][2] * Disp_grad_overlay[e][N][2]);
				printf("EMT_mode2[0][0] = %.15e\n", EMT_mode2[0][0]);
				/* W - s11 * du1/dx1 - s12 * du2/dx1 */
				EMT_mode2[0][1] =                            - (Stress_overlay[e][N][0] * Disp_grad_aux_mode2[e][N][1] + Stress_aux_mode2[e][N][0] * Disp_grad_overlay[e][N][1] + Stress_overlay[e][N][2] * Disp_grad_aux_mode2[e][N][3] + Stress_aux_mode2[e][N][2] * Disp_grad_overlay[e][N][3]);
				printf("EMT_mode2[0][1] = %.15e\n", EMT_mode2[0][1]);
				/* - s11 * du1/dx2 - s22 * du2/dx2 */
				EMT_mode2[1][0] =                            - (Stress_overlay[e][N][2] * Disp_grad_aux_mode2[e][N][0] + Stress_aux_mode2[e][N][2] * Disp_grad_overlay[e][N][0] + Stress_overlay[e][N][1] * Disp_grad_aux_mode2[e][N][2] + Stress_aux_mode2[e][N][1] * Disp_grad_overlay[e][N][2]);
				printf("EMT_mode2[1][0] = %.15e\n", EMT_mode2[1][0]);
				/* - s21 * du1/dx1 - s22 * du2/dx1 */
				EMT_mode2[1][1] = StrainEnergyDensity_aux_mode2[e][N] - (Stress_overlay[e][N][2] * Disp_grad_aux_mode2[e][N][1] + Stress_aux_mode2[e][N][2] * Disp_grad_overlay[e][N][1] + Stress_overlay[e][N][1] * Disp_grad_aux_mode2[e][N][3] + Stress_aux_mode2[e][N][1] * Disp_grad_overlay[e][N][3]);
				printf("EMT_mode2[1][1] = %.15e\n", EMT_mode1[1][1]);
				/* W - s21 * du1/dx2 - s22 * du2/dx2 */

				//補助場のみでのモードIにおけるEMT
				EMT_aux_mode1[0][0] =  StrainEnergyDensity_aux_mode1[e][N] - (Stress_aux_mode1[e][N][0] * Disp_grad_aux_mode1[e][N][0] + Stress_aux_mode1[e][N][2] * Disp_grad_aux_mode1[e][N][2]);
				printf("EMT_aux_mode1[0][0] = %.15e\n",EMT_aux_mode1[0][0]);
				/* W - s11 * du1/dx1 - s12 * du2/dx1 */
				EMT_aux_mode1[0][1] =                            - (Stress_aux_mode1[e][N][0] * Disp_grad_aux_mode1[e][N][1] + Stress_aux_mode1[e][N][2] * Disp_grad_aux_mode1[e][N][3]);
				printf("EMT_aux_mode1[0][1] = %.15e\n",EMT_aux_mode1[0][1]);
				/* - s11 * du1/dx2 - s22 * du2/dx2 */
				EMT_aux_mode1[1][0] =                            - (Stress_aux_mode1[e][N][2] * Disp_grad_aux_mode1[e][N][0] + Stress_aux_mode1[e][N][1] * Disp_grad_aux_mode1[e][N][2]);
				printf("EMT_aux_mode1[1][0] = %.15e\n",EMT_aux_mode1[1][0]);
				/* - s21 * du1/dx1 - s22 * du2/dx1 */
				EMT_aux_mode1[1][1] =  StrainEnergyDensity_aux_mode1[e][N] - (Stress_aux_mode1[e][N][2] * Disp_grad_aux_mode1[e][N][1] + Stress_aux_mode1[e][N][1] * Disp_grad_aux_mode1[e][N][3]);
				printf("EMT_aux_mode1[1][1] = %.15e\n",EMT_aux_mode1[1][1]);
				/* W - s21 * du1/dx2 - s22 * du2/dx2 */

				//補助場のみでのモードIIにおけるEMT
				EMT_aux_mode2[0][0] =  StrainEnergyDensity_aux_mode2[e][N] - (Stress_aux_mode2[e][N][0] * Disp_grad_aux_mode2[e][N][0] + Stress_aux_mode2[e][N][2] * Disp_grad_aux_mode2[e][N][2]);
				printf("EMT_aux_mode2[0][0] = %.15e\n",EMT_aux_mode2[0][0]);
				/* W - s11 * du1/dx1 - s12 * du2/dx1 */
				EMT_aux_mode2[0][1] =                            - (Stress_aux_mode2[e][N][0] * Disp_grad_aux_mode2[e][N][1] + Stress_aux_mode2[e][N][2] * Disp_grad_aux_mode2[e][N][3]);
				printf("EMT_aux_mode2[0][1] = %.15e\n",EMT_aux_mode2[0][1]);
				/* - s11 * du1/dx2 - s22 * du2/dx2 */
				EMT_aux_mode2[1][0] =                            - (Stress_aux_mode2[e][N][2] * Disp_grad_aux_mode2[e][N][0] + Stress_aux_mode2[e][N][1] * Disp_grad_aux_mode2[e][N][2]);
				printf("EMT_aux_mode2[1][0] = %.15e\n",EMT_aux_mode2[1][0]);
				/* - s21 * du1/dx1 - s22 * du2/dx1 */
				EMT_aux_mode2[1][1] =  StrainEnergyDensity_aux_mode2[e][N] - (Stress_aux_mode2[e][N][2] * Disp_grad_aux_mode2[e][N][1] + Stress_aux_mode2[e][N][1] * Disp_grad_aux_mode2[e][N][3]);
				printf("EMT_aux_mode2[1][1] = %.15e\n",EMT_aux_mode2[1][1]);
				/* W - s21 * du1/dx2 - s22 * du2/dx2 */

				// printf("chekc8 re = %d  e = %d N = %d\n",re,e,N);

				// printf("check strain energy e N %d %d %.15e\n",e, N, StrainEnergyDensity_aux_mode1[e][N]);
				// printf("check Stresses e N %d %d %.15e %.15e %.15e\n",e, N, Stress_overlay[e][N][0], Stress_overlay[e][N][1], Stress_overlay[e][N][2]);
				// printf("check Strain e N %d %d %.15e %.15e %.15e\n",e, N, Strain_overlay[e][N][0], Strain_overlay[e][N][1], Strain_overlay[e][N][2]);
				// printf("check Disp_grad e N %d %d %.15e %.15e %.15e %.15e\n",e, N, Disp_grad_overlay[e][N][0], Disp_grad_overlay[e][N][1], Disp_grad_overlay[e][N][2], Disp_grad_overlay[e][N][3]);

				printf("check q_func_grad e N %d %d %.15e %.15e %.15e %.15e\n", e, N, q_func_grad[0][0], q_func_grad[0][1], q_func_grad[1][0], q_func_grad[1][1]);
				/* q_func_grad[0][0] = -(GXY[N][0] - 25.0) / sqrt((GXY[N][0] - 25.0)*(GXY[N][0] - 25.0)+GXY[N][1]*GXY[N][1])/5.0;
								q_func_grad[0][1] = -(GXY[N][1] - 0.0) / sqrt((GXY[N][0] - 25.0)*(GXY[N][0] - 25.0)+GXY[N][1]*GXY[N][1])/5.0;

					printf("check q_func_grad2 e N %d %d %.15e %.15e %.15e %.15e\n",e, N, q_func_grad[0][0], q_func_grad[0][1],q_func_grad[1][0], q_func_grad[1][1]); */
				double rrrr = sqrt((GXY[N][0] - 25.0) * (GXY[N][0] - 25.0) + (GXY[N][1] - 0.0) * (GXY[N][1] - 0.0));
				
				//モード1
				K_1 -= (EMT_mode1[0][0] * q_func_grad[0][0] + EMT_mode1[1][0] * q_func_grad[0][1]) * w[N] * J;
				K_1 -= (EMT_mode1[0][1] * q_func_grad[1][0] + EMT_mode1[1][1] * q_func_grad[1][1]) * w[N] * J; 

				//モード2
				K_2 -= (EMT_mode2[0][0] * q_func_grad[0][0] + EMT_mode2[1][0] * q_func_grad[0][1]) * w[N] * J;
				K_2 -= (EMT_mode2[0][1] * q_func_grad[1][0] + EMT_mode2[1][1] * q_func_grad[1][1]) * w[N] * J; 

				//補助場のモードI
				J_integral_value_aux_mode1 -= (EMT_aux_mode1[0][0] * q_func_grad[0][0] + EMT_aux_mode1[1][0] * q_func_grad[0][1]) * w[N] * J;
				J_integral_value_aux_mode1 -= (EMT_aux_mode1[0][1] * q_func_grad[1][0] + EMT_aux_mode1[1][1] * q_func_grad[1][1]) * w[N] * J;

				//補助場のモードII
				J_integral_value_aux_mode2 -= (EMT_aux_mode2[0][0] * q_func_grad[0][0] + EMT_aux_mode2[1][0] * q_func_grad[0][1]) * w[N] * J;
				J_integral_value_aux_mode2 -= (EMT_aux_mode2[0][1] * q_func_grad[1][0] + EMT_aux_mode2[1][1] * q_func_grad[1][1]) * w[N] * J;

				CheckArea += J * w[N];
				CheckArea_Norm += w[N];
				double NPP = 2.0;
				check_int_r += 1.0 / rrrr / pow(5.0, NPP-1) * pow(5.0 - rrrr, NPP-1) * J * w[N];
				/*check_int_r += 1.0/rrrr * J * w[N];*/
			}
			printf("Area = %.15e\n", CheckArea);
			printf("Area Norm = %.15e\n", CheckArea_Norm);
			/* printf("check_int_r = %.15e 1/ 5*Pi = %.15e\n",check_int_r,4.0*atan(1.0)*5.0/5.0);*/
			printf("check_int_r = %.15e 5*Pi = %.15e\n", check_int_r, 4.0 * atan(1.0) * 5.0 / 2);
		}
	}
	printf("check_int_r = %.15e  5*Pi = %.15e\n", check_int_r, 4.0 * atan(1.0) * 5.0 / 2);
	fp = fopen("Virtual_Crack_Extension.dat", "w");
	fprintf(fp, "label=Virtual_Crack_Extension\n");
	fprintf(fp, "num_items=%d\n", Total_Control_Point);
	fprintf( fp, "\n");
	for(j = 0; j < Total_Control_Point; j++ ){
		fprintf(fp, "%d:	%le %le ", j, Virtual_Crack_Extension_El_Nodes[j * DIMENSION + 0], Virtual_Crack_Extension_El_Nodes[j * DIMENSION + 1]);
		fprintf(fp, "\n");
	}
	fclose(fp); 

	if(DM == 0){
		K_mode1 = K_1 / DeltaA / (2.0 * (1.0 / E));
		K_mode2 = K_2 / DeltaA / (2.0 * (1.0 / E));
	}else if(DM == 1){
		K_mode1 = K_1 / DeltaA / (2.0 * ((1.0 - nu * nu) / E));
		K_mode2 = K_2 / DeltaA / (2.0 * ((1.0 - nu * nu) / E));
	}else{
		return ERROR;
	}
	return 0;
}

/*
double J_Integral_Computation6by6(int Total_Control_Point, int Total_Element,int Location_Crack_Tip_Patch, double Location_Local_Coordinates[DIMENSION], double Virtual_Crack_Extension_Ct_Pt[MAX_N_NODE][DIMENSION], double DeltaA, int El_No, double E, double nu, int DM)
{
	int ii, jj, kk; 
	int e, i, j, k, re, N;
	double B[D_MATRIX_SIZE][MAX_KIEL_SIZE], J;

	double G1D6[6] = {-0.932469514203152, -0.661209386466265, -0.238619186083197, 0.238619186083197, 0.661209386466265, 0.932469514203152};
	double w1D6[6] = {0.171324492379170, 0.360761573048139, 0.467913934572691, 0.467913934572691, 0.360761573048139, 0.171324492379170};
	double G2D36org[36][2], w2D36org[36]; 
	//double G = pow(0.6,0.5);
	//double Gxi[POW_Ng][DIMENSION] = { {-G,-G},{0.0,-G},{G,-G},{-G,0.0},{0.0,0.0},{G,0.0},{-G,G},{0.0,G},{G,G} };
	double J_integral_value=0.0;
	double q_func_grad[DIMENSION][DIMENSION], EMT[DIMENSION][DIMENSION];
	double Virtual_Crack_Extension_El_Nodes[MAX_NO_CCpoint_ON_ELEMENT*DIMENSION];
	double Strain_Gauss[3], Stress2D_Gauss[4], StrainEnergyDensity_Gauss,Disp_grad_Gauss[4];
	double U[MAX_KIEL_SIZE];
	double D[D_MATRIX_SIZE][D_MATRIX_SIZE];
	double GXY[36*25][2], qXY[2];
	double GXY_Corner[8][2]={{-1.0, -1.0},{0.0, -1.0}, {1.0, -1.0},{1.0, 0.0}, {1.0, 1.0},{0.0, 1.0}, {-1.0, 1.0}, {-1.0, 0.0}};
	double check_int_r = 0.0;
	int idiv, jdiv;
	int Idiv = 1, Jdiv = 1;
	double G2D36[36*5*5][2], w2D36[36*5*5];

	ii = 0;
	for(i=0; i < 6; i++) 
		for(j=0; j < 6; j++){
			G2D36org[ii][0] = G1D6[j]; G2D36org[ii][1] = G1D6[i];
			w2D36org[ii] = w1D6[j] * w1D6[i];
			ii++;
			}

	ii=0;
	for(idiv=0; idiv<Idiv; idiv++)
		for(jdiv=0; jdiv<Jdiv; jdiv++)
			for(jj=0; jj<36; jj++){
			G2D36[ii][0] = -1.0 + 1.0/Jdiv + 2.0/Jdiv*jdiv + G2D36org[jj][0] / Jdiv;
			G2D36[ii][1] = -1.0 + 1.0/Idiv + 2.0/Idiv*idiv + G2D36org[jj][1] / Idiv;
			w2D36[ii] = w2D36org[jj]/(Idiv*Jdiv);
			ii++;
			}
			


	Make_D_Matrix_2D( D, E, nu , DM);			

	for( re = 0; re < real_Total_Element_to_mesh[Total_mesh]; re++ ){
		e=real_element[re];
		// printf("re = %d  e = %d\n",re,e);
		double  X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION];
		double b_grad[DIMENSION * DIMENSION][2 * MAX_NO_CCpoint_ON_ELEMENT];
		double XY_Corner[8][2] = {{0.0, 0.0},{0.0, 0.0},{0.0, 0.0},{0.0, 0.0},{0.0, 0.0},{0.0, 0.0},{0.0, 0.0},{0.0, 0.0}};
		double VC_value[8][2] =  {{0.0, 0.0},{0.0, 0.0},{0.0, 0.0},{0.0, 0.0},{0.0, 0.0},{0.0, 0.0},{0.0, 0.0},{0.0, 0.0}}; 
		for (ii=0; ii < DIMENSION; ii++) for (jj=0; jj < DIMENSION; jj++) q_func_grad[ii][jj] = 0.0;

		//printf("%d,No_Control_point_ON_ELEMENT[%d]:%d\n",e,Element_patch[e], No_Control_point_ON_ELEMENT[Element_patch[e]]);
		//Bマトリックスを取得
		for( i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[e]]; i++ ){
			for( j = 0; j < DIMENSION; j++ ){
				X[i][j] = Node_Coordinate[ Controlpoint_of_Element[e][i] ][j];
				U[ i*DIMENSION +j ] = Displacement[ Controlpoint_of_Element[e][i]*DIMENSION + j ];
				Virtual_Crack_Extension_El_Nodes[i * DIMENSION + j] = Virtual_Crack_Extension_Ct_Pt[Controlpoint_of_Element[e][i]][j];
				}	
				//printf("Node_Coordinate[%d][%d]:%le\n",Controlpoint_of_Element[e][i],j, Node_Coordinate[ Controlpoint_of_Element[e][i] ][j]);
				//printf("X[%d][%d]:%le\n",Controlpoint_of_Element[e][i],j,X[i][j] );
			}
			double zero_thresh = 1.0e-12;
			int zero_flag = 0;
				for(kk = 0; kk < 2 * No_Control_point_ON_ELEMENT[Element_patch[e]]; kk++)
					if(fabs(Virtual_Crack_Extension_El_Nodes[kk]) < zero_thresh) zero_flag++;
						if(zero_flag < 2 * No_Control_point_ON_ELEMENT[Element_patch[e]]){

							for(k=0; k<8; k++)
								for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[e]]; i++){
									double R_shape_func = Shape_func(i,Total_Control_Point, GXY_Corner[k],e);
									for (j = 0; j < DIMENSION; j++){ 
										XY_Corner[k][j] += R_shape_func * X[i][j];
										VC_value[k][j] += R_shape_func * Virtual_Crack_Extension_El_Nodes[i * DIMENSION + j]; 
									}
								}

		printf("*** Coordinates and Virtual Crack Extension at 4 corners of the patch ***\n");
			for(k = 0; k < 8; k++) printf("Corner Coords V-Crack Extension e k %d %d %.15e %.15e %.15e %.15e\n",e,k,XY_Corner[k][0], XY_Corner[k][1],VC_value[k][0],VC_value[k][1]);                                        
		printf("***                                                                   ***\n");

			printf("chekc2 re = %d  e = %d\n", re, e);
		for(ii = 0; ii < No_Control_point_ON_ELEMENT[Element_patch[e]]; ii++) 
			printf("Virtual Crack Extension %.15e %.15e\n",Virtual_Crack_Extension_El_Nodes[ii * DIMENSION],Virtual_Crack_Extension_El_Nodes[ii * DIMENSION + 1]);
		double CheckArea = 0.0, CheckAreaNorm=0.0;
		for( N = 0; N < 36*Idiv*Jdiv; N++ ){
		// computeing the coordinates of the Gauss Point
		GXY[N][0] = 0.0; GXY[N][1] = 0.0; qXY[0] = 0.0; qXY[1] = 0.0;

		for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[e]]; i++){
				double R_shape_func = Shape_func(i, Total_Control_Point, G2D36[N], e);
				for (j = 0; j < DIMENSION; j++){
					GXY[N][j] += R_shape_func * X[i][j];
					qXY[j] += R_shape_func * Virtual_Crack_Extension_El_Nodes[i * DIMENSION + j];
			}
		}

		printf("6^6chekc3 re = %d  e = %d N = %d  qXY = %lf\n",re,e,N, qXY[0]);
		for (ii=0; ii < DIMENSION; ii++) for (jj=0; jj < DIMENSION; jj++){
			q_func_grad[ii][jj] = 0.0;
			EMT[ii][jj] = 0.0;
			}
		printf("chekc4 re = %d  e = %d N = %d\n",re,e,N);
		Make_b_grad_Matrix( e, b_grad, G2D36[N], X ,&J , Total_Control_Point);
		Disp_grad_Gauss[0] = 0.0; Disp_grad_Gauss[1] = 0.0;Disp_grad_Gauss[2] = 0.0;Disp_grad_Gauss[3] = 0.0;
		printf("check5 re = %d  e = %d N = %d J = %.15e\n",re,e,N,J);
			for( ii = 0; ii < DIMENSION; ii++)
				for(jj = 0; jj < DIMENSION; jj++)			
					for( j = 0; j < 2 * No_Control_point_ON_ELEMENT[Element_patch[e]]; j++ ){
						q_func_grad[ii][jj] += b_grad[ii*DIMENSION + jj][j] * Virtual_Crack_Extension_El_Nodes[j] ;
						Disp_grad_Gauss[ii*DIMENSION + jj] += b_grad[ii*DIMENSION + jj][j] *  U[j];}

		Make_B_Matrix( e, B, G2D36[N], X ,&J , Total_Control_Point);
			Strain_Gauss[0]=0.0; Strain_Gauss[1]=0.0;Strain_Gauss[2]=0.0;	
			Stress2D_Gauss[0] = 0.0; Stress2D_Gauss[1] = 0.0; Stress2D_Gauss[2] = 0.0; Stress2D_Gauss[3] = 0.0;
				for( i = 0; i < D_MATRIX_SIZE; i++ )
									for( j = 0; j < KIEL_SIZE; j++ )
												Strain_Gauss[i] += B[i][j] * U[j] ;
			
			for( i = 0; i < D_MATRIX_SIZE; i ++ )
									for( j = 0; j < D_MATRIX_SIZE; j++ )
											Stress2D_Gauss[i] += D[i][j] * Strain_Gauss[j];
			StrainEnergyDensity_Gauss = 0.5 * (Stress2D_Gauss[0] * Strain_Gauss[0] + Stress2D_Gauss[1] * Strain_Gauss[1] + Stress2D_Gauss[2] * Strain_Gauss[2]); 
		

			// printf("check 9 %.15e  %.15e  %.15e  %.15e  %.15e\n",StrainEnergyDensity_overlay[e][N],Stress_overlay[e][N][0],Disp_grad_overlay[e][N][0],Stress_overlay[e][N][2],Disp_grad_overlay[e][N][2]);
			// printf("chekc6 re = %d  e = %d N = %d\n",re,e,N);
			// //Setting up the Energy Momentum Tensor
			// printf("check 9 %.15e  %.15e  %.15e  %.15e  %.15e\n",StrainEnergyDensity_overlay[e][N],Stress_overlay[e][N][0],Disp_grad_overlay[e][N][0],Stress_overlay[e][N][2],Disp_grad_overlay[e][N][2]);
			EMT[0][0] =  StrainEnergyDensity_Gauss - (Stress2D_Gauss[0] * Disp_grad_Gauss[0] + Stress2D_Gauss[2] * Disp_grad_Gauss[2]);
			// printf("G E: %.15e %.15e\n",Gxi[N][0], Gxi[N][1]);
			printf("EMT[0][0] = %.15e\n",EMT[0][0]);
			// W - s11 * du1/dx1 - s12 * du2/dx1
			EMT[0][1] =                            - (Stress2D_Gauss[0] * Disp_grad_Gauss[1] + Stress2D_Gauss[2] * Disp_grad_Gauss[3]);
			printf("EMT[0][1] = %.15e\n",EMT[0][1]);
			//  - s11 * du1/dx2 - s22 * du2/dx2
			EMT[1][0] =                            - (Stress2D_Gauss[2] * Disp_grad_Gauss[0] + Stress2D_Gauss[1] * Disp_grad_Gauss[2]);
			printf("EMT[1][0] = %.15e\n",EMT[1][0]);
			//  - s21 * du1/dx1 - s22 * du2/dx1
			EMT[1][1] =  StrainEnergyDensity_Gauss - (Stress2D_Gauss[2] * Disp_grad_Gauss[1] + Stress2D_Gauss[1] * Disp_grad_Gauss[3]);
			printf("EMT[1][1] = %.15e\n",EMT[1][1]);
			//  W - s21 * du1/dx2 - s22 * du2/dx2

			printf("check Coordinates q-values G-E e N %d %d %.15e %.15e %.15e %.15e %.15e %.15e\n",e,N,GXY[N][0],GXY[N][1],qXY[0],qXY[1],G2D36[N][0],G2D36[N][1]);
			printf("check strain energry e N %d %d %.15e\n",e, N, StrainEnergyDensity_Gauss);
			printf("check Stresses e N %d %d %.15e %.15e %.15e\n",e, N, Stress2D_Gauss[0], Stress2D_Gauss[1], Stress2D_Gauss[2]);
			printf("check Strain e N %d %d %.15e %.15e %.15e\n",e, N, Strain_Gauss[0], Strain_Gauss[1], Strain_Gauss[2]);
			printf("check Dispgrad e N %d %d %.15e %.15e %.15e %.15e\n",e, N, Disp_grad_Gauss[0], Disp_grad_Gauss[1], Disp_grad_Gauss[2], Disp_grad_Gauss[3]);

			printf("check q_func_grad e N %d %d %.15e %.15e %.15e %.15e\n",e, N, q_func_grad[0][0], q_func_grad[0][1],q_func_grad[1][0], q_func_grad[1][1]);
			
			// q_func_grad[0][0] = -(GXY[N][0] - 25.0) / sqrt((GXY[N][0] - 25.0)*(GXY[N][0] - 25.0)+GXY[N][1]*GXY[N][1])/5.0;
			// q_func_grad[0][1] = -(GXY[N][1] - 0.0) / sqrt((GXY[N][0] - 25.0)*(GXY[N][0] - 25.0)+GXY[N][1]*GXY[N][1])/5.0;

			double rrrr = sqrt((GXY[N][0] - 25.0)*(GXY[N][0] - 25.0) +(GXY[N][1] - 0.0)*(GXY[N][1] - 0.0));
			// double rrr0 = 5.0;
			// 				q_func_grad[0][0] = 2.0 * (rrrr - rrr0) / (rrr0*rrr0) * (GXY[N][0] - 25.0)/rrrr;
			// 				q_func_grad[0][1] = 2.0 * (rrrr - rrr0) / (rrr0*rrr0) * (GXY[N][1] - 0.0)/rrrr;

			// printf("check q_func_grad2 e N %d %d %.15e %.15e %.15e %.15e\n",e, N, q_func_grad[0][0], q_func_grad[0][1],q_func_grad[1][0], q_func_grad[1][1]);

			J_integral_value -= (EMT[0][0] * q_func_grad[0][0] + EMT[1][0] * q_func_grad[0][1]) * w2D36[N] * J;
			//  EMT11 * dq1/dq1 + EMT21 * dq1/dx2
			J_integral_value -= (EMT[0][1] * q_func_grad[1][0] + EMT[1][1] * q_func_grad[1][1]) * w2D36[N] * J; 
			//  EMT12 * dq2/dq1 + EMT22 * dq2/dx2
			CheckArea += J*w2D36[N];
			CheckAreaNorm += w2D36[N];
			printf("check rrrr = %.15e\n",rrrr);
			double NPP=2.0;
			//  check_int_r += 1.0/rrrr * J * w2D36[N];
			check_int_r += 1.0/rrrr / pow(5.0, NPP-1) * pow(5.0 - rrrr, NPP-1) * J * w2D36[N]; 

			printf("check7 re = %d  e = %d N = %d J = %.15e\n",re,e,N,J);
			}
			printf("Area = %.15e\n",CheckArea);
			printf("Area Norm = %.15e\n",CheckAreaNorm);
			//  printf("check_int_r = %.15e  Pi = %.15e\n",check_int_r,4.0*atan(1.0)*5.0/5.0);
			printf("check_int_r = %.15e  Pi*5/NPP = %.15e\n",check_int_r,4.0*atan(1.0)*5.0/2);
		}

		}
		printf("check_int_r = %.15e  5*Pi/NPP = %.15e\n",check_int_r,4.0*atan(1.0)*5.0/2);
	return(J_integral_value/DeltaA);
}
*/

//ファイルからデータをもらう
void Get_InputData(int tm,
                   double *E, double *nu, int *Total_Element, int *Total_Control_Point,
				   int *Total_Load, int *No_Patch, int Load_Node_Dir[MAX_N_LOAD][2], double Value_of_Load[MAX_N_LOAD],
				   int *Total_Constraint, int Constraint_Node_Dir[MAX_N_CONSTRAINT][2], double Value_of_Constraint[MAX_N_CONSTRAINT],
				   int *Total_DistributeForce, char *argv[])
{
	int i, j, k, l, iii;
	int n, p, q, h, x, y;
	char s[256];
	int ii, jj, kk, kkk;
	int e, b, B;
	int iiloc, jjloc, kkloc;
	int r = 0;
	/* for the distributed loads*/

	if ((fp = fopen(argv[tm+1], "r")) == NULL)
		printf("file open error!!\n");
	//次元数
	/*fscanf(fp,"%d",&DIMENSION );
	printf("DIMENSION=%d\n",DIMENSION);
	fgets(s, 256, fp);*/
	//材料定数
	fscanf(fp, "%le %le", &*E, &*nu);
	fgets(s, 256, fp);
	printf("E:%le nu:%le\n", *E, *nu);
    //E[tm]=*E;
    //nu[tm]=*nu;

	//パッチ数
	fscanf(fp, "%d", &*No_Patch);
	fgets(s, 256, fp);
	printf("No_Patch:%d\n", *No_Patch);
    Total_Patch_on_mesh[tm] = *No_Patch;
	Total_Patch_to_mesh[tm+1] = Total_Patch_to_mesh[tm] + *No_Patch;
	//Total_Patch_to_mesh[tm+1] += *No_Patch;
	printf("Total_Patch_to_mesh[%d]=%d\n",
			tm,Total_Patch_to_mesh[tm]);
	Total_Patch_to_Now = Total_Patch_to_mesh[tm];
    printf("Total_Patch_to_Now=%d\n",Total_Patch_to_Now);

	//コントロールポイント数
	fscanf(fp, "%d", &*Total_Control_Point);
	fgets(s, 256, fp);
	printf("Total_Control_Point:%d\n", *Total_Control_Point);
    Total_Control_Point_on_mesh[tm] = *Total_Control_Point;
	Total_Control_Point_to_mesh[tm+1] = Total_Control_Point_to_mesh[tm] + *Total_Control_Point;
    //Total_Control_Point_to_mesh[tm+1] += *Total_Control_Point;
    Total_Control_Point_to_Now = Total_Control_Point_to_mesh[tm];
    printf("Total_Control_Point_to_mesh[%d]=%d\n",
            tm,Total_Control_Point_to_mesh[tm]);
	printf("Total_Control_Point_to_mesh[%d]=%d\n",
            tm+1,Total_Control_Point_to_mesh[tm+1]);
    printf("Total_Control_Point_to_Now=%d\n",Total_Control_Point_to_Now);

	//ξη方向の各次数
	for (l = 0; l < *No_Patch; l++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
            fscanf(fp, "%d", &Order[l+Total_Patch_to_mesh[tm]][j]);
		    printf("Order[%d][%d]=%d\n", l+Total_Patch_to_mesh[tm], j
                                       , Order[l+Total_Patch_to_mesh[tm]][j]);
		}
	}

	fgets(s, 256, fp);
	//ノット数
	for (l = 0; l < *No_Patch; l++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
            fscanf(fp, "%d", &No_knot[l+Total_Patch_to_mesh[tm]][j]);
		    printf("No_knot[%d][%d]=%d\n", l+Total_Patch_to_mesh[tm], j
                                         , No_knot[l+Total_Patch_to_mesh[tm]][j]);
		}
	}

	fgets(s, 256, fp);
	//各パッチ各方向のコントロールポイント数
	for (l = 0; l < *No_Patch; l++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
    		fscanf(fp, "%d", &No_Control_point[l+Total_Patch_to_mesh[tm]][j]);
			printf("No_Control_point[%d][%d]:%d\n", l+Total_Patch_to_mesh[tm], j
                                                  , No_Control_point[l+Total_Patch_to_mesh[tm]][j]);
		}
	}

	fgets(s, 256, fp);

	for (l = 0; l < *No_Patch; l++)
	{
		No_Controlpoint_in_patch[l+Total_Patch_to_mesh[tm]] = 1.0;
	}

	for (l = 0; l < *No_Patch; l++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			No_Controlpoint_in_patch[l+Total_Patch_to_mesh[tm]]
                    *= No_Control_point[l+Total_Patch_to_mesh[tm]][j];
		}
	}

	for (l = 0; l < *No_Patch; l++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			if (No_knot[l+Total_Patch_to_mesh[tm]][j] != No_Control_point[l+Total_Patch_to_mesh[tm]][j] + Order[l+Total_Patch_to_mesh[tm]][j] + 1)
			{
				printf("wrong relationship between the number of knot vector and the number of control_point \n");
				printf("in mesh_No.%d in patch_No.%d direction:%d\n",tm, l, j);
			}
		}
	}

	for (l = 0; l < *No_Patch; l++)
	{
		printf("No_Controlpoint_in_patch[%d]:%d\t", l+Total_Patch_to_mesh[tm]
                                                  , No_Controlpoint_in_patch[l+Total_Patch_to_mesh[tm]]);
	}
	printf("\n");

	for (l = 0; l < *No_Patch; l++)
	{
		//printf("l;%d\n",l);
		for (i = 0; i < No_Controlpoint_in_patch[l+Total_Patch_to_mesh[tm]]; i++)
		{
			//printf("i:%d\n",i );
			fscanf(fp, "%d", &Patch_controlpoint[l+Total_Patch_to_mesh[tm]][i]);
			if (tm > 0)
			{
				Patch_controlpoint[l+Total_Patch_to_mesh[tm]][i] += Total_Control_Point_to_mesh[tm];
			}
			//printf("Patch_controlpoint[%d][%d]=%d\n",
			//		l+Total_Patch_to_mesh[tm], i,
			//		Patch_controlpoint[l+Total_Patch_to_mesh[tm]][i]);
		}
	}

	/*for (l = 0; l < *No_Patch; l++) {
			//printf("l:%d\n",l);
			for ( i = 0; i < No_Controlpoint_in_patch[l]; i++) {
				//printf("i:%d\n", i);
				printf("Patch_controlpoint[%d][%d]:%d\n",l,i,Patch_controlpoint[l][i]);
			}printf("\n");
		}*/

	fscanf(fp, "%d %d %d", Total_Constraint, Total_Load, Total_DistributeForce);
	Total_Constraint_on_mesh[tm] = *Total_Constraint;
	Total_Constraint_to_mesh[tm+1] = Total_Constraint_to_mesh[tm]
								   + *Total_Constraint;
    if( tm == Total_mesh-1 )
    {
        for(i = 0 ; i < Total_mesh ; i++ )
        {
            Total_Constraint_all_mesh += Total_Constraint_on_mesh[i];
        }
        printf("Total_Constraint_all_mesh;%d\n",Total_Constraint_all_mesh);
    }
    printf("Total_Constraint;%d\n", *Total_Constraint);
    printf("Total_Constraint_on_mesh[%d];%d\n",tm, Total_Constraint_on_mesh[tm]);
    printf("Total_Constraint_to_mesh[%d];%d\n",tm, Total_Constraint_to_mesh[tm]);
	Total_Load_on_mesh[tm] = *Total_Load;
	Total_Load_to_mesh[tm+1] = Total_Load_to_mesh[tm] + *Total_Load;
	printf("Total_Load;%d\n", *Total_Load);
	printf("Total_DistributedForce;%d\n", *Total_DistributeForce);
	Total_DistributeForce_on_mesh[tm] = *Total_DistributeForce;
	Total_DistributeForce_to_mesh[tm+1] = Total_DistributeForce_to_mesh[tm]
										+ *Total_DistributeForce;
	fgets(s, 256, fp);

	//ノットベクトルの読み込み
	for (l = 0; l < *No_Patch; l++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			for (k = 0; k < No_knot[l+Total_Patch_to_mesh[tm]][j]; k++)
			{
				fscanf(fp, "%le", &Position_Knots[l+Total_Patch_to_mesh[tm]][j][k]);
				printf("%le\t", Position_Knots[l+Total_Patch_to_mesh[tm]][j][k]);
                //printf("Position_Knots[%d][%d][%d]=%le\n",
                //        l+Total_Patch_to_mesh[tm], j, k,
                //        Position_Knots[l+Total_Patch_to_mesh[tm]][j][k]);
			}
			printf("\n");
		}
	}

	for (l = 0; l < *No_Patch; l++)
	{
		No_Control_point_ON_ELEMENT[l+Total_Patch_to_mesh[tm]] = 1.0;
	}

	*Total_Element = 0.0;

	for (l = 0; l < *No_Patch; l++)
	{
		if (DIMENSION == 2)
		{
			*Total_Element += (No_Control_point[l+Total_Patch_to_mesh[tm]][0] - Order[l+Total_Patch_to_mesh[tm]][0])
                            * (No_Control_point[l+Total_Patch_to_mesh[tm]][1] - Order[l+Total_Patch_to_mesh[tm]][1]);
			No_Control_point_ON_ELEMENT[l+Total_Patch_to_mesh[tm]]
                = (Order[l+Total_Patch_to_mesh[tm]][0] + 1) * (Order[l+Total_Patch_to_mesh[tm]][1] + 1);
		}
		else
		{
			*Total_Element += (No_Control_point[l+Total_Patch_to_mesh[tm]][0] - Order[l+Total_Patch_to_mesh[tm]][0])
                            * (No_Control_point[l+Total_Patch_to_mesh[tm]][1] - Order[l+Total_Patch_to_mesh[tm]][1])
                            * (No_Control_point[l+Total_Patch_to_mesh[tm]][2] - Order[l+Total_Patch_to_mesh[tm]][2]);
			No_Control_point_ON_ELEMENT[l+Total_Patch_to_mesh[tm]]
                = (Order[l+Total_Patch_to_mesh[tm]][0] + 1) * (Order[l+Total_Patch_to_mesh[tm]][1] + 1) * (Order[l+Total_Patch_to_mesh[tm]][2] + 1);
		}
	}
	printf("Total_Element=%d\n", *Total_Element);
    Total_Element_on_mesh[tm] = *Total_Element;
	Total_Element_to_mesh[tm+1] = Total_Element_to_mesh[tm] + *Total_Element;
	//Total_Element_to_mesh[tm] += *Total_Element;
	Total_Element_to_Now = Total_Element_to_mesh[tm];
	printf("Total_Element_on_mesh[%d]=%d\n",
			tm , Total_Element_on_mesh[tm]);
	printf("Total_Element_to_Now=%d\n",
			Total_Element_to_Now);

	for ( l = 0; l < *No_Patch; l++)
    {
		printf("No_Control_point_ON_ELEMENT[%d]=%d\n",
               l+Total_Patch_to_mesh[tm],
			   No_Control_point_ON_ELEMENT[l+Total_Patch_to_mesh[tm]]);
	}

	//節点座標
	for (i = 0; i < *Total_Control_Point; i++)
	{
		fscanf(fp, "%d", &ii);
		for (j = 0; j < DIMENSION + 1; j++)
			fscanf(fp, "%le", &Node_Coordinate[ii+Total_Control_Point_to_mesh[tm]][j]); //Node_Coordinate[i][2]:重み
	}
	for (i = 0; i < *Total_Control_Point; i++)
	{
		for (j = 0; j < DIMENSION + 1; j++)
		{
			//コントロールポイント座標・重みの新たな配列（for s-IGA/NewtonLaphson）
			if (j < DIMENSION)
			{
				Control_Coord[j][i+Total_Control_Point_to_mesh[tm]]
					=Node_Coordinate[i+Total_Control_Point_to_mesh[tm]][j];
			}
			else if (j == DIMENSION)
			{
				Control_Weight[i+Total_Control_Point_to_mesh[tm]]
					=Node_Coordinate[i+Total_Control_Point_to_mesh[tm]][DIMENSION];
			}
			printf("Node_Coordinate[%d][%d]=%e\n", i+Total_Control_Point_to_mesh[tm], j, Node_Coordinate[i+Total_Control_Point_to_mesh[tm]][j]);
		}
	}
	fgets(s, 256, fp);


	//拘束
	for (i = 0; i < *Total_Constraint; i++)
		fscanf(fp, "%d %d %le",
					&Constraint_Node_Dir[i+Total_Constraint_to_mesh[tm]][0],
					&Constraint_Node_Dir[i+Total_Constraint_to_mesh[tm]][1],
					&Value_of_Constraint[i+Total_Constraint_to_mesh[tm]]);
	for (i = 0; i < *Total_Constraint; i++)
    {
		Constraint_Node_Dir[i+Total_Constraint_to_mesh[tm]][0]
			= Constraint_Node_Dir[i+Total_Constraint_to_mesh[tm]][0] + Total_Control_Point_to_mesh[tm];
		Constraint_Node_Dir_on_mesh[tm][i][0] = Constraint_Node_Dir[i][0];
		Constraint_Node_Dir_on_mesh[tm][i][1] = Constraint_Node_Dir[i][1];
        Value_of_Constraint_on_mesh[tm][i] = Value_of_Constraint[i];
        printf("Constraint_Node_Dir[%d][0]= %d Constraint_Node_Dir[%d][1]=%d Value_of_Constraint[%d]= %e \n",
				i+Total_Constraint_to_mesh[tm],
				Constraint_Node_Dir[i+Total_Constraint_to_mesh[tm]][0],
				i+Total_Constraint_to_mesh[tm],
				Constraint_Node_Dir[i+Total_Constraint_to_mesh[tm]][1],
				i+Total_Constraint_to_mesh[tm],
				Value_of_Constraint[i+Total_Constraint_to_mesh[tm]]);
    }
    fgets(s, 256, fp);

	//荷重
	for (i = 0; i < *Total_Load; i++)
	{
		fscanf(fp, "%d %d %le",
					&Load_Node_Dir[i+Total_Load_to_mesh[tm]][0],
					&Load_Node_Dir[i+Total_Load_to_mesh[tm]][1],
					&Value_of_Load[i+Total_Load_to_mesh[tm]]);
		Load_Node_Dir[i+Total_Load_to_mesh[tm]][0]
			= Load_Node_Dir[i+Total_Load_to_mesh[tm]][0] + Total_Control_Point_to_mesh[tm];

		printf("Load_Node_Dir[%d][0]= %d Load_Node_Dir[%d][1]= %d Value_of_Load[%d]= %e\n",
				i+Total_Load_to_mesh[tm],
				Load_Node_Dir[i+Total_Load_to_mesh[tm]][0],
				i+Total_Load_to_mesh[tm],
				Load_Node_Dir[i+Total_Load_to_mesh[tm]][1],
				i+Total_Load_to_mesh[tm],
				Value_of_Load[i+Total_Load_to_mesh[tm]]);
	}

	int iPatch, iCoord, type_load;
	double Range_Coord[2], val_Coord, Coeff_Dist_Load[3];
	//int iPatch_array[MAX_N_DISTRIBUTE_FORCE], iCoord_array[MAX_N_DISTRIBUTE_FORCE], type_load_array[MAX_N_DISTRIBUTE_FORCE];
	//double val_Coord_array[MAX_N_DISTRIBUTE_FORCE], Range_Coord_array[MAX_N_DISTRIBUTE_FORCE][2], Coeff_Dist_Load_array[MAX_N_DISTRIBUTE_FORCE][3];

	fgets(s, 256, fp);
	for (i = 0; i < *Total_DistributeForce; i++)
	{
		fscanf(fp, "%d %d %d %lf %lf %lf %lf %lf %lf", &type_load, &iPatch, &iCoord, &val_Coord, &Range_Coord[0], &Range_Coord[1], &Coeff_Dist_Load[0], &Coeff_Dist_Load[1], &Coeff_Dist_Load[2]);
		printf("Distibuted load nober: %d\n", i);
		printf("type_load: %d  iPatch: %d iCoord: %d  val_Coord: %.15e  Range_Coord: %.15e  %.15e\n Coef_Dist_Load: %.15e %.15e %.15e\n",
			   type_load, iPatch, iCoord,
			   val_Coord, Range_Coord[0], Range_Coord[1], Coeff_Dist_Load[0], Coeff_Dist_Load[1], Coeff_Dist_Load[2]);
		/*
 		type_load: Direction of distributed load: 0-x direction, 1-y direction, 2-normal to the segemet/surface
		iPatch: Patch number to which the distributed load is assigned., 0, 1, ...
		iCoord: 0: Distributed load is applied to line along Xi axis.
                        1: Distributed load is applied to line along Eta axis
		val_Coord: その時のもう片方の座標
		Range_Coord[0]: Local coordinate value at which the distributed load starts.
		Range_Coord[1]: Local coordinate value at which the distributed load ends.
		Coeff_Dist_Load[0], &Coeff_Dist_Load[1], &Coeff_Dist_Load[2]: The coefficients of distributed load value:
			Coeff_Dist_Load[0]*Xi + Coeff_Dist_Load[1]*Xi + Coeff_Dist_Load[2]*Xi^2
		or
			Coeff_Dist_Load[0]*Xi + Coeff_Dist_Load[1]*Eta + Coeff_Dist_Load[2]*Eta^2
  		*/
	  	//for s-IGA
		type_load_array[i+Total_DistributeForce_to_mesh[tm]] = type_load;
		iPatch_array[i+Total_DistributeForce_to_mesh[tm]] = iPatch;
		iCoord_array[i+Total_DistributeForce_to_mesh[tm]] = iCoord;
		val_Coord_array[i+Total_DistributeForce_to_mesh[tm]] = val_Coord;
		Range_Coord_array[i+Total_DistributeForce_to_mesh[tm]][0] = Range_Coord[0];
		Range_Coord_array[i+Total_DistributeForce_to_mesh[tm]][1] = Range_Coord[1];
		Coeff_Dist_Load_array[i+Total_DistributeForce_to_mesh[tm]][0] = Coeff_Dist_Load[0];
		Coeff_Dist_Load_array[i+Total_DistributeForce_to_mesh[tm]][1] = Coeff_Dist_Load[1];
		Coeff_Dist_Load_array[i+Total_DistributeForce_to_mesh[tm]][2] = Coeff_Dist_Load[2];
	}

	/* Setting_Dist_Load_2D(&*Total_Control_Point,&iPatch, &*Total_Element, iCoord, val_Coord,
		     Range_Coord[2],  norm_dir, type_load, Coeff_Dist_Load); */
	//自由度共有の計算
	//(同じ座標を計算して要素コネクティビィティのコントロールポイント番号を入れ替える)
	/*(2018_01_31)for (ii = 0; ii < *Total_Control_Point; ii++) {
			same_point[ii]=ii;
		}

		for (ii = 0 ; ii < *Total_Control_Point; ii++) {
			for ( jj = ii-1; jj >= 0 ; jj--) {
				if (Node_Coordinate[ii][0]== Node_Coordinate[jj][0] && Node_Coordinate[ii][1]==Node_Coordinate[jj][1]) {
					printf("同じ座標の番号ii:%d jj:%d\n",ii,jj);
					same_point[ii]=jj;
					//printf("same_point_1[%d]:%d\n",ii,same_point[ii]);
				}
			}
		}*/

	/*	for (ii = 0; ii < *Total_Control_Point; ii++) {
			printf("same_point[%d]:%d\n",ii,same_point[ii]);
		}*/
	//INC\の計算（節点番号をξ、ηの番号で表す為の配列）
	if (DIMENSION == 2) //for s-IGA
	{
		e = 0;
		for (l = 0; l < *No_Patch; l++)
		{
			i = 0;
			for (jj = 0; jj < No_Control_point[l+Total_Patch_to_Now][1]; jj++)
			{
				for (ii = 0; ii < No_Control_point[l+Total_Patch_to_Now][0]; ii++)
				{

					INC[l+Total_Patch_to_Now][Patch_controlpoint[l+Total_Patch_to_Now][i]][0] = ii;
					INC[l+Total_Patch_to_Now][Patch_controlpoint[l+Total_Patch_to_Now][i]][1] = jj;
					//printf("INC[%d][%d][0]=%d INC[%d][%d][1]=%d\n"
                    //        ,l+Total_Patch_to_Now, Patch_controlpoint[l+Total_Patch_to_Now][i]
                    //        ,INC[l+Total_Patch_to_Now][Patch_controlpoint[l+Total_Patch_to_Now][i]][0]
                    //        ,l+Total_Patch_to_Now, Patch_controlpoint[l+Total_Patch_to_Now][i]
                    //        ,INC[l+Total_Patch_to_Now][Patch_controlpoint[l+Total_Patch_to_Now][i]][1]);
					//printf("INC[%d][0]=ii;%d INC[%d][1]=jj;%d\n",i, ii, i,jj);
					//printf("INC[%d][0]=%d INC[%d][1]=%d\n",i, INC[i][0], i,INC[i][1] );
					Adress_Controlpoint[l+Total_Patch_to_Now][ii][jj] = Patch_controlpoint[l+Total_Patch_to_Now][i];
                    //printf("Adress_con[%d][%d][%d]=%d\n",
                    //        l+Total_Patch_to_Now,ii,jj,Adress_Controlpoint[l+Total_Patch_to_Now][ii][jj]);

					// if (ii >= Order[l][0] && jj >= Order[l+Total_Patch_to_Now][1])
					if (ii >= Order[l+Total_Patch_to_Now][0] && jj >= Order[l+Total_Patch_to_Now][1])
					{

						for (jjloc = 0; jjloc <= Order[l+Total_Patch_to_Now][1]; jjloc++)
						{
							for (iiloc = 0; iiloc <= Order[l+Total_Patch_to_Now][0]; iiloc++)
							{
								//printf("jjloc:%d iiloc:%d\n",jjloc,iiloc);
								B = Patch_controlpoint[l+Total_Patch_to_Now][i - jjloc * No_Control_point[l+Total_Patch_to_Now][0] - iiloc];
								b = jjloc * (Order[l+Total_Patch_to_Now][0] + 1) + iiloc;
								//printf("B=%d b=%d e=%d\n",B,b,e);
								//Controlpoint_of_Element[e][b] = B;
								Controlpoint_of_Element[e +Total_Element_to_Now][b] //= B;
                                    = B;
                                //	= B + Total_Control_Point_to_Now;
                                //printf("Controlpoint_of_Element[%d][%d]=%d\n",
                                //        e+Total_Element_to_Now,b,
								// 		Controlpoint_of_Element[e+Total_Element_to_Now][b]);
							}
						}
						//Element_patch[e] = l+Total_Patch_to_Now;
                        Element_patch[e+Total_Element_to_Now]
                            = l+Total_Patch_to_Now;
						//printf("Element_patch[%d]=%d\n",
						//		e+Total_Element_to_Now,
						//		Element_patch[e+Total_Element_to_Now]);

                        Element_mesh[e+Total_Element_to_Now] = tm;
                        //printf("Element_mesh[%d]=%d\n",
                        //        e+Total_Element_to_Now,
                        //        Element_mesh[e+Total_Element_to_Now]);

                        El_No_on_mesh[tm][e] = e+Total_Element_to_Now;
                        //printf("El_No_on_mesh[%d][%d]=%d\n",
                        //        tm,e,El_No_on_mesh[tm][e]);
						e++;
					}
					i++;
				}
			}
            Patch_mesh[l+Total_Patch_to_Now] = tm;
            //printf("Patch_mesh[%d]=%d\n",
            //        l+Total_Patch_to_Now,
            //        Patch_mesh[l+Total_Patch_to_Now]);

		}
		/*for ( i = 0; i < *Total_Control_Point; i++){
				for ( l = 0; l < *No_Patch; l++) {
					printf("INC[%d][%d][0]=%d INC[%d][%d][1]=%d\n",l,i, INC[l][i][0], l,i,INC[l][i][1] );
				}
		}*/

		/*for ( l = 0; l < *No_Patch; l++) {
			for (i = 0; i < *Total_Element; i++) {
					for (b = 0; b < No_Control_point_ON_ELEMENT[l]; b++) {
					printf("Controlpoint_of_Element[%d][%d]=%d\n",i,b,Controlpoint_of_Element[i][b]);
				}
			}
		}*/

		/*for ( i = 0; i < *Total_Element; i++) {
			printf("Element_patch[%d]:%d\n",i,Element_patch[i]);
		}*/
		/*for (i = 0; i < *Total_Element; i++) {
				for (b = 0; b < No_Control_point_ON_ELEMENT; b++) {
					printf("Controlpoint_of_Element_before[%d][%d]=%d\n",i,b,Controlpoint_of_Element[i][b]);
					Controlpoint_of_Element[i][b]=same_point[Controlpoint_of_Element[i][b]];
					printf("Controlpoint_of_Element_after[%d][%d]=%d\n",i,b,Controlpoint_of_Element[i][b]);
				}
			}*/
	}

	if (DIMENSION == 3)
	{
		e = 0;
		for (l = 0; l < *No_Patch; l++)
		{
			i = 0;
			for (kk = 0; kk < No_Control_point[l][2]; kk++)
			{
				for (jj = 0; jj < No_Control_point[l][1]; jj++)
				{
					for (ii = 0; ii < No_Control_point[l][0]; ii++)
					{

						//printf("kk=%d\n",kk );
						INC[l][Patch_controlpoint[l][i]][0] = ii;
						INC[l][Patch_controlpoint[l][i]][1] = jj;
						INC[l][Patch_controlpoint[l][i]][2] = kk;
						//printf("INC[%d][0]=%d INC[%d][1]=%d\n",i, INC[i][0], i,INC[i][1] );
						if (ii >= Order[l][0] && jj >= Order[l][1] && kk >= Order[l][2])
						{
							for (kkloc = 0; kkloc < Order[l][2]; kkloc++)
							{
								for (jjloc = 0; jjloc <= Order[l][1]; jjloc++)
								{
									for (iiloc = 0; iiloc <= Order[l][0]; iiloc++)
									{
										//printf("jjloc:%d iiloc:%d\n",jjloc,iiloc);
										B = Patch_controlpoint[l][i - jjloc * No_Control_point[l][0] - iiloc];
										b = jjloc * (Order[l][0] + 1) + iiloc;
										//printf("B=%d b=%d e=%d\n",B,b,e);
										Controlpoint_of_Element[e][b] = B;
									}
								}
							}
							Element_patch[e] = l+Total_Patch_to_Now;
							e++;
						}
						i++;
					}
				}
			}
		}
	}
	/*----------------------------------------------------------------------------------------------*/

	//#include<stdio.h>

	//#define DIMENSION               2
	//#define MAX_N_KNOT 				1000
	//#define MAX_N_ELEMENT 		    110000

	//int main()
	//{
	//static double Position_Knots[DIMENSION][MAX_N_KNOT];    /*ノットベクトル*/
	//static double difference[MAX_N_KNOT][DIMENSION];    /*隣り合うノットベクトルの差*/
	//static int ENC[MAX_N_ELEMENT][DIMENSION];   /*ENC[全ての要素][0,1]=x,y方向の何番目の要素か*/
	//int /*i,j,k,e,*/m,n,h,p,q,x,y;
	//int kk,rr;
	//int Order[DIMENSION];   /*次数*/
	//int No_knot[DIMENSION];   /*ノット数*/
	//int *Total_Element;  /*ゼロエレメントを含むすべての要素の数*/
	//int Total_element_all_ID[MAX_N_ELEMENT];    /*ゼロエレメントではない要素＝１、ゼロエレメント＝０*/
	//int line_No_Total_element[DIMENSION];   /*ゼロエレメントを含むすべての要素列の数*/
	//int line_No_real_element[DIMENSION];   /*ゼロエレメントではない要素列の数*/
	//int real_element_line[MAX_N_ELEMENT][DIMENSION];   /*ゼロエレメントではない要素列*/
	//int real_element[MAX_N_ELEMENT];    /*ゼロエレメントではない要素の番号*/
	//int r=0;

    //for s-IGA line_No_real_elementの初期化
    for (l = 0; l < *No_Patch; l++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
           line_No_real_element[l+Total_Patch_to_Now][j]=0;
        }
    }

	for (l = 0; l < *No_Patch; l++)
	{
		for (j = 0; j < DIMENSION; j++)
		{

			//printf("Order[%d]= ",j);
			//scanf("%d",&Order[j]);

			//printf("No_knot[%d]= ",j);
			//scanf("%d",&No_knot[j]);

			line_No_Total_element[l+Total_Patch_to_Now][j]
                = No_knot[l+Total_Patch_to_Now][j] - 2
                * Order[l+Total_Patch_to_Now][j] - 1;

			//for(i=0;i<No_knot[j];i++){
			//printf("Position_Knots[%d][%d]= ",j,i);
			//scanf("%lf",&Position_Knots[j][i]);
			//}

			for (kkk = Order[l+Total_Patch_to_Now][j]; kkk < No_knot[l+Total_Patch_to_Now][j] - Order[l+Total_Patch_to_Now][j] - 1; kkk++)
			{
				difference[l+Total_Patch_to_Now][kkk - Order[l+Total_Patch_to_Now][j]][j]
                    = Position_Knots[l+Total_Patch_to_Now][j][kkk + 1]
                    - Position_Knots[l+Total_Patch_to_Now][j][kkk];
				//printf("[[[%d]]] ξ[%d]-ξ[%d]=%lf\n",kkk-Order[l][j],kkk+1,kkk,difference[kkk-Order[l][j]][j]);

				if (difference[l+Total_Patch_to_Now][kkk - Order[l+Total_Patch_to_Now][j]][j] != 0)
				{
					line_No_real_element[l+Total_Patch_to_Now][j]++;
				}
			}
			//  printf("line_No_real_element[%d][%d]=%d\n",l,j,line_No_real_element[l][j]);
		}
	}

	/*要素に行番号、列番号をつける*/

	if (DIMENSION == 2)
	{
		//printf("Total_element_all= ");
		//scanf("%d",&Total_element_all);

		for (h = 0; h < *Total_Element; h++)
		{
			Total_element_all_ID[h] = 0;
			//Total_element_all_ID[h]=h;
		}

		i = 0;
		for (l = 0; l < *No_Patch; l++)
		{
			for (y = 0; y < line_No_Total_element[l+Total_Patch_to_Now][1]; y++)
			{
				for (x = 0; x < line_No_Total_element[l+Total_Patch_to_Now][0]; x++)
				{
					ENC[l+Total_Patch_to_mesh[tm]][i+Total_Element_to_mesh[tm]][0] = x;
					ENC[l+Total_Patch_to_mesh[tm]][i+Total_Element_to_mesh[tm]][1] = y;
					//printf("ENC[%d][%d][0]=%d\tENC[%d][%d][1]=%d\n",
					//		l+Total_Patch_to_mesh[tm],
					//		i+Total_Element_to_mesh[tm],
					//		ENC[l+Total_Patch_to_mesh[tm]][i+Total_Element_to_mesh[tm]][0],
					//		l+Total_Patch_to_mesh[tm],
					//		i+Total_Element_to_mesh[tm],
					//		ENC[l+Total_Patch_to_mesh[tm]][i+Total_Element_to_mesh[tm]][1]);
					i++;
				}
			}
		}

	}

	/*必要な要素の行と列の番号を求める*/

	for (j = 0; j < DIMENSION; j++)
	{
		for (l = 0; l < *No_Patch; l++)
		{
			e=0;
			for (k = 0; k < line_No_Total_element[l+Total_Patch_to_Now][j]; k++)
			{
				//printf("//%d,%d//\n",j,line_No_Total_element[j]);
				if (difference[l+Total_Patch_to_Now][k][j] != 0)
				{
					//printf("k=%d\n",k);
					real_element_line[l+Total_Patch_to_Now][e][j] = k;
                    //printf("real_element_line[%d][%d]=%d\n",e,j,real_element_line[e][j]);
					e++;
					//printf("e=%d\n",e);
				}
			}

			/*for(e=0;e<line_No_real_element[l][j];e++){
            printf("real_element_line[%d][%d]=%d\n",e,j,real_element_line[e][j]);
        }*/
		}
	}

	/*必要な要素列上の要素のIDを1にする*/

	if (DIMENSION == 2)
	{
		for (n = 0; n < *Total_Element; n++)
		{
			//for(j=0;j<DIMENSION;j++){
			for (p = 0; p < line_No_real_element[Element_patch[n+Total_Element_to_Now]][0]; p++)
			{
				if (ENC[Element_patch[n+Total_Element_to_mesh[tm]]][n+Total_Element_to_mesh[tm]][0] == real_element_line[Element_patch[n+Total_Element_to_Now]][p][0])
				{
					for (q = 0; q < line_No_real_element[Element_patch[n+Total_Element_to_Now]][1]; q++)
					{
						if (ENC[Element_patch[n+Total_Element_to_mesh[tm]]][n+Total_Element_to_mesh[tm]][1] == real_element_line[Element_patch[n+Total_Element_to_Now]][q][1])
						{
							//Total_element_all_ID[n]++;
							Total_element_all_ID[n]++;
                            //printf("Total_element_all_ID[%d]=%d\n"
                            //        ,n,Total_element_all_ID[n]);
							//break;
						}
					}
				}
			}
			//}

			//printf("Total_element_all_ID[%d]=%d\n",n,Total_element_all_ID[n]);

			/*IDが1の要素に番号を振る*/

			if (Total_element_all_ID[n] == 1)
			{
				real_element[r+real_Total_Element_to_mesh[tm]]
                    = n + Total_Element_to_Now;
				// printf("real_element[%d]=%d\n",
                //        r + real_Total_Element_to_mesh[tm],
                //        real_element[r + real_Total_Element_to_mesh[tm]]);
                real_El_No_on_mesh[tm][r] = n + Total_Element_to_Now;
                //printf("real_El_No_on_mesh[%d][%d]=%d\n",
                //        tm,r,real_El_No_on_mesh[tm][r]);
				r++;
			}
		}
        //for s-IGA real_Total_Elementの初期化
        real_Total_Element = 0;

		for (l = 0; l < *No_Patch; l++)
		{
			real_Total_Element
                += line_No_real_element[l+Total_Patch_to_Now][0]
                 * line_No_real_element[l+Total_Patch_to_Now][1];
		}
		real_Total_Element_on_mesh[tm] = real_Total_Element;
        real_Total_Element_to_mesh[tm+1] = real_Total_Element_to_mesh[tm] + real_Total_Element;
		//real_Total_Element_to_mesh[tm] += real_Total_Element;
		real_Total_Element_to_Now = real_Total_Element_to_mesh[tm];
		//printf("real_Total_Element_on_mesh[%d]=%d\n",
		//		tm, real_Total_Element_on_mesh[tm]);
		//printf("real_Total_Element_to_mesh[%d]=%d\n",
		//		tm, real_Total_Element_to_mesh[tm]);
	}

	/*
    for(rr=0;rr<real_Total_Element;rr++){
        printf("real_element[%d]=%d\n",rr,real_element[rr]);
    }
    */
	//}
	//
	/* For distributed load 2D */

	for (iii = 0; iii < *Total_Control_Point; iii++)
	{
		Equivalent_Nodal_Force[iii+Total_Control_Point_to_mesh[tm]][0] = 0.0;
		Equivalent_Nodal_Force[iii+Total_Control_Point_to_mesh[tm]][1] = 0.0;
	}

	for (i = 0; i < *Total_DistributeForce; i++)
	{

		type_load = type_load_array[i+Total_DistributeForce_to_mesh[tm]];
		iPatch = iPatch_array[i+Total_DistributeForce_to_mesh[tm]];
		iCoord = iCoord_array[i+Total_DistributeForce_to_mesh[tm]];
		val_Coord = val_Coord_array[i+Total_DistributeForce_to_mesh[tm]];
		Range_Coord[0] = Range_Coord_array[i+Total_DistributeForce_to_mesh[tm]][0];
		Range_Coord[1] = Range_Coord_array[i+Total_DistributeForce_to_mesh[tm]][1];
		Coeff_Dist_Load[0] = Coeff_Dist_Load_array[i+Total_DistributeForce_to_mesh[tm]][0];
		Coeff_Dist_Load[1] = Coeff_Dist_Load_array[i+Total_DistributeForce_to_mesh[tm]][1];
		Coeff_Dist_Load[2] = Coeff_Dist_Load_array[i+Total_DistributeForce_to_mesh[tm]][2];

		printf("type_load:%d\tiPatch:%d\tiCoord:%d\tval_coord:%lf\t",
				type_load, iPatch, iCoord, val_Coord);
		printf("Range0:%lf\tRange1:%lf\t",
				Range_Coord[0], Range_Coord[1]);
		printf("Coeff0:%lf\n",Coeff_Dist_Load[0]);
		Setting_Dist_Load_2D(tm, iPatch, Total_Element_to_mesh[tm+1], iCoord, val_Coord,
							 Range_Coord, type_load, Coeff_Dist_Load);
	}
	/*-------------------------------------------------------------------------------------*/
}


//拘束されている行数を省いた行列の番号の制作
int Make_Index_Dof(int Total_Control_Point,
				   int Total_Constraint,
				   int Constraint_Node_Dir[MAX_N_CONSTRAINT][2])
{
	int i, k = 0;

    //printf("Total_Constraint=%d\n",Total_Constraint);
    //Index_Dofの初期化(複数メッシュ読み込みのため)
    for (i = 0; i < Total_Control_Point * 2; i++)
    {
        Index_Dof[i] = 0;
        //printf("test_DOF[%d]:%d\n",i,Index_Dof[i]);
    }
    //拘束されている自由度(Degree Of free)をERRORにする
	for (i = 0; i < Total_Constraint; i++)
    {
        //printf("%d\n",Constraint_Node_Dir[i][0] * DIMENSION + Constraint_Node_Dir[i][1]);
		Index_Dof[Constraint_Node_Dir[i][0] * DIMENSION + Constraint_Node_Dir[i][1]] = ERROR;
    }
    //ERROR以外に番号を付ける
	for (i = 0; i < Total_Control_Point * DIMENSION; i++)
	{
		if (Index_Dof[i] != ERROR)
		{
			Index_Dof[i] = k;
			k++;
		}
        //printf("Index_Dof[%d]=%d\n",i,Index_Dof[i]);
	}
	printf("Max_Index_Dof=%d\n", k);
	return k;
}


void Make_K_Whole_Ptr_Col(int Total_Element,
						  int Total_Control_Point,
						  int K_Whole_Size)
{
	int i, ii, j, jj, k;
	int NE;
	int N, i_index, j_index;

	//初期化
	//for (i = 0; i < Total_Control_Point * DIMENSION; i++)
		//Total_Control_Point_To_Node[i] = 0;
	for (i = 0; i < K_Whole_Size + 1; i++)
		K_Whole_Ptr[i] = 0;

	for (N = 0; N < Total_Control_Point; N += K_DIVISION_LENGE)
	{ 	//大きく分割するためのループ
		//各節点に接する節点を取得
		for (i = 0; i < K_DIVISION_LENGE; i++)
		{
			Total_Control_Point_To_Node[i] = 0;
		}
        for (i = 0; i < Total_Element; i++)
		{
			for (ii = 0; ii < No_Control_point_ON_ELEMENT[Element_patch[i]]; ii++)
			{
				NE = Controlpoint_of_Element[i][ii] - N;
				//printf("NE=%d\n",NE);
				//printf("K_DIVISION_LENGE=%d,N=%d,NE=%d\n",K_DIVISION_LENGE,N,NE);    //K_DIVISION_LENGE=0,N=0,NE=コネクティビティ的な
				if (0 <= NE && NE < K_DIVISION_LENGE)
				{
					for (j = 0; j < No_Control_point_ON_ELEMENT[Element_patch[i]]; j++)	//ローカル要素
					{
						//printf("j=%d\n",j);
						//数字がない時
						if (Total_Control_Point_To_Node[NE] == 0)
						{
							//節点番号を取得
							Node_To_Node[NE][0] = Controlpoint_of_Element[i][j];
							Total_Control_Point_To_Node[NE]++;
							//printf("Node_To_Node[%d][0]=%d\n",NE,Node_To_Node[NE][0]);
						}
						//printf("②Total_Control_Point_To_Node[%d]=%d\n",NE,Total_Control_Point_To_Node[NE]);
						//同じものがあったら
						//k > 0 以降の取得
						//kのカウント
						for (k = 0; k < Total_Control_Point_To_Node[NE]; k++)
						{
							//printf("Node_To_Node[%d][%d]=%d\n",NE,k,Node_To_Node[NE][k]);
							//
							//printf("k_1=%d\t",k);
							if (Node_To_Node[NE][k] == Controlpoint_of_Element[i][j])
							{
								//printf("break\t");
								break;
							}
						}
						//printf("\nk_2=%d\n",k);
						//未設定のNode_To_Node取得
						if (k == Total_Control_Point_To_Node[NE])
						{
							Node_To_Node[NE][k] = Controlpoint_of_Element[i][j];
							//printf("Node_To_Node[%d][%d]=%d\n",NE,k,Node_To_Node[NE][k]);
							Total_Control_Point_To_Node[NE]++;
							//printf("③Total_Control_Point_To_Node[%d]=%d\n",NE,Total_Control_Point_To_Node[NE]);
						}
					}
					//別メッシュとの重なりを考慮
					if (NNLOVER[i] > 0)
					{
						for (jj = 0; jj < NNLOVER[i]; jj++)
						{
							for (j = 0; j < No_Control_point_ON_ELEMENT[Element_patch[NELOVER[i][jj]]]; j++)	//ローカル要素
							{
								//printf("j=%d\n",j);
								//数字がない時
								if (Total_Control_Point_To_Node[NE] == 0)
								{
									//節点番号を取得
									Node_To_Node[NE][0] = Controlpoint_of_Element[NELOVER[i][jj]][j];
									Total_Control_Point_To_Node[NE]++;
									//printf("Node_To_Node[%d][0]=%d\n",NE,Node_To_Node[NE][0]);
								}
								//printf("②Total_Control_Point_To_Node[%d]=%d\n",NE,Total_Control_Point_To_Node[NE]);
								//同じものがあったら
								//k > 0 以降の取得
								//kのカウント
								for (k = 0; k < Total_Control_Point_To_Node[NE]; k++)
								{
									//printf("Node_To_Node[%d][%d]=%d\n",NE,k,Node_To_Node[NE][k]);
									//
									//printf("k_1=%d\t",k);
									if (Node_To_Node[NE][k] == Controlpoint_of_Element[NELOVER[i][jj]][j])
									{
										//printf("break\t");
										break;
									}
								}
								//printf("\nk_2=%d\n",k);
								//未設定のNode_To_Node取得
								if (k == Total_Control_Point_To_Node[NE])
								{
									Node_To_Node[NE][k] = Controlpoint_of_Element[NELOVER[i][jj]][j];
									//printf("Node_To_Node[%d][%d]=%d\n",NE,k,Node_To_Node[NE][k]);
									Total_Control_Point_To_Node[NE]++;
									//printf("③Total_Control_Point_To_Node[%d]=%d\n",NE,Total_Control_Point_To_Node[NE]);
								}
							}
						}
					}
				}
				//printf("\n");
			}
			//printf("\n");
		}
		//順番に並び替える
		for (i = 0; i < K_DIVISION_LENGE; i++)
		{
			if (N + i < Total_Control_Point)
			{
				//printf("Node[%d] T=%d; \n",N+i, Total_Control_Point_To_Node[ i ]);
				for (j = 0; j < Total_Control_Point_To_Node[i]; j++)
				{
					int Min = Node_To_Node[i][j], No = j;
					for (k = j; k < Total_Control_Point_To_Node[i]; k++)
					{
						if (Min > Node_To_Node[i][k])
						{
							Min = Node_To_Node[i][k];
							No = k;
						}
					}
					for (k = No; k > j; k--)
					{
						Node_To_Node[i][k] = Node_To_Node[i][k - 1];
					}
					Node_To_Node[i][j] = Min;
					//				printf("%d ",Node_To_Node[i][j]);
				}
				//			printf("\n");
			}

			//並べ替えたNode_To_Node確認
			// for (j = 0; j < Total_Control_Point_To_Node[i]; j++)
			// {
			// 	printf("sort_Node_To_Node[%d][%d]=%d\n",i,j,Node_To_Node[i][j]);
			// }
		}

		//節点からcol ptrを求める
		ii = 0;
		k = 0;
		for (i = 0; i < K_DIVISION_LENGE; i++)
		{
			for (ii = 0; ii < DIMENSION; ii++)
			{
				if (N + i < Total_Control_Point)
				{
					i_index = Index_Dof[(N + i) * DIMENSION + ii];
					// printf("i = %d\n", i);
					// printf("N = %d\n", N);
					// printf("i_index = %d\n", i_index);
					k = 0;
					if (i_index >= 0)
					{
                        //K_Whole_Ptr[i_index + 1] = K_Whole_Ptr[i_index];
                        K_Whole_Ptr[i_index + 1] = K_Whole_Ptr[i_index];
                        //printf("K_Whole_Ptr[%d][%d]=%d\n",tm,i_index,K_Whole_Ptr[tm][i_index+1]);
						for (j = 0; j < Total_Control_Point_To_Node[i]; j++)
						{
							// printf("Total_Control_Point_To_Node[%d] = %d\n", i, Total_Control_Point_To_Node[i]);
							for (jj = 0; jj < DIMENSION; jj++)
							{
								j_index = Index_Dof[Node_To_Node[i][j] * DIMENSION + jj];
								if (j_index >= 0 && j_index >= i_index)
								{
									K_Whole_Ptr[i_index + 1]++;
									//col_N[N/K_DIVISION_LENGE][k] = j_index;
									K_Whole_Col[K_Whole_Ptr[i_index] + k] = j_index;
                                    //printf("K_Whole_Col[%d]=%d\n"
                                    //       ,K_Whole_Ptr[i_index]+k
                                    //       ,K_Whole_Col[K_Whole_Ptr[i_index]+k]);
									k++;
									// printf("ptr[%d]=%d,col[%d]=%d\n",i_index+1,K_Whole_Ptr[i_index+1],K_Whole_Ptr[i_index]+k,K_Whole_Col[K_Whole_Ptr[i_index]+k]);
								}
							}
						}
					}
				}
			}
		}
		//col_N[N/K_DIVISION_LENGE][ k ] = -1;
	}
    //for(i=0;i<K_Whole_Size_array[tm]+1;i++)
	/*
    for(i=0;i<K_Whole_Size+1;i++)
    {
        printf("K_Whole_Ptr[%d]=%d\n",
                i,K_Whole_Ptr[i]);
    }*/
	/*
	for( i = 0; i < K_Whole_Size+1; i++ )//printf("K_Whole_Ptr[%d]= %d\n",i,K_Whole_Ptr[i]);
	//col合成
	k = 0;
	for( N = 0; N < Total_Control_Point ; N +=K_DIVISION_LENGE ){
		for(i = 0; col_N[ N/K_DIVISION_LENGE ][i] != -1; i++ ){
			K_Whole_Col[k] = col_N[ N/K_DIVISION_LENGE ][i];
			k++;
		}
	}
    */
}


//valを求める
void Make_K_Whole_Val(double E, double nu, int Total_Element, int DM)
{
	int i, j, j1, j2, k1, k2, l;
	int a, b, re;

	for (i = 0; i < MAX_NON_ZERO; i++)
	{
		K_Whole_Val[i] = 0.0;
	}

	for (re = 0; re < Total_Element; re++)
	{
		i = real_element[re];
		Check_BDBJ_flag[i] = 0;
	}

	/*for(rr=0;rr<line_No_real_element[0]*line_No_real_element[1];rr++){
        	printf("real_element[%d]=%d\n",rr,real_element[rr]);
    	}
		printf("real_Total_element=%d\n",real_Total_Element);*/
	//re=0;
	for (re = 0; re < Total_Element; re++)
	{
		//if(i==real_element[re]){
		//printf("re=%d\n",re);
		i = real_element[re];
		// printf("El_No;i=%d\n", real_element[re]);

		if(Element_mesh[i] == 0 && re == 0)/*2つめの条件は効率化のため*/
		{
			Make_gauss_array(0);
		}
		else if(Element_mesh[i] > 0)
		{
			printf("NNLOVER[%d]:%d\tNNLOVER[%d]:%d\tElement_mesh[%d]:%d\n", i, NNLOVER[i], real_element[re - 1], NNLOVER[real_element[re - 1]], real_element[re - 1], Element_mesh[real_element[re - 1]]);
			if(NNLOVER[i] == 1 && (NNLOVER[real_element[re - 1]] != 1 || Element_mesh[real_element[re - 1]] == 0))/*2つめ以降の条件は効率化のため*/
			{
				Make_gauss_array(0);
			}
			else if(NNLOVER[i] >= 2 && (NNLOVER[real_element[re - 1]] == 1 || Element_mesh[real_element[re - 1]] == 0))/*2つめ以降の条件は効率化のため*/
			{
				Make_gauss_array(1);
			}
		}
		//printf("i= %d\tGaussPt_3D=%d\n",i ,GaussPt_3D);

		for (j = 0; j < GP_2D; j++)
		{
			Same_BDBJ_flag[j] = 0;
		}

		KIEL_SIZE = No_Control_point_ON_ELEMENT[Element_patch[i]] * DIMENSION;
		double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], K_EL[MAX_KIEL_SIZE][MAX_KIEL_SIZE];
		// printf("Total_Element=%d\tre=%d\tEl_No=%d\n", Total_Element, re, i);
		//各要素のKelを求める
		for (j1 = 0; j1 < No_Control_point_ON_ELEMENT[Element_patch[i]]; j1++)
		{
			for (j2 = 0; j2 < DIMENSION; j2++)
			{
				X[j1][j2] = Node_Coordinate[Controlpoint_of_Element[i][j1]][j2];
			}
		}

		Make_K_EL(i, X, K_EL, E, nu, DM);

		//Valを求める
		//for (j1 = 0; j1 < No_Control_point_ON_ELEMENT[Element_patch[El_No_on_mesh[tm][i]]]; j1++)
		for (j1 = 0; j1 < No_Control_point_ON_ELEMENT[Element_patch[i]]; j1++)
		{
			for (j2 = 0; j2 < DIMENSION; j2++)
			{
				a = Index_Dof[Controlpoint_of_Element[i][j1] * DIMENSION + j2];
				if (a >= 0)
				{
					//for (k1 = 0; k1 < No_Control_point_ON_ELEMENT[Element_patch[El_No_on_mesh[tm][i]]]; k1++)
					for (k1 = 0; k1 < No_Control_point_ON_ELEMENT[Element_patch[i]]; k1++)
					{
						for (k2 = 0; k2 < DIMENSION; k2++)
						{
							b = Index_Dof[Controlpoint_of_Element[i][k1] * DIMENSION + k2];
							if (b >= 0 && b >= a)
							{
								for (l = K_Whole_Ptr[a]; l < K_Whole_Ptr[a + 1]; l++)
								{
									if (K_Whole_Col[l] == b)
									{
										// printf("real_element[%d] = %d\n", re, real_element[re]);
										// printf("Controlpoint_of_Element[%d][%d] = %d\n", i, k1, Controlpoint_of_Element[i][k1]);
										// printf("l = %d\n", l);
										// printf("K_Whole_Col[%d] = %d\n", l, K_Whole_Col[l]);
										// printf("b = %d\n", b);
										K_Whole_Val[l] += K_EL[j1 * DIMENSION + j2][k1 * DIMENSION + k2];
										// printf("K_Whole_Val[%d] = %.20e\n", l, K_Whole_Val[l]);
										//printf("K_Whole_Val[%d]=%le\n",l,K_Whole_Val[l]);
										break;
									}
								}
							}
						}
					}
				}
			}
		}

		if (Element_mesh[i] > 0)	//ローカルメッシュ上の要素について
		{
			if (NNLOVER[i] > 0)		//重なっている要素が存在するとき
			{
				for (j = 0; j < NNLOVER[i]; j++)
				{
					double XG[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION];
					KIEL_SIZE = No_Control_point_ON_ELEMENT[Element_patch[NELOVER[i][j]]] * DIMENSION;
					double coupled_K_EL[MAX_KIEL_SIZE][MAX_KIEL_SIZE];
					for (j1 = 0; j1 < No_Control_point_ON_ELEMENT[Element_patch[NELOVER[i][j]]]; j1++)
					{
						for (j2 = 0; j2 < DIMENSION; j2++)
						{
							XG[j1][j2] = Node_Coordinate[Controlpoint_of_Element[NELOVER[i][j]][j1]][j2];
							//重なっている要素の物理座標取得
						}
					}
					Make_coupled_K_EL(i, NELOVER[i][j],
									  X,
									  XG,
									  coupled_K_EL,
									  E, nu, DM);
					
					Check_BDBJ_flag[i] += Total_BDBJ_flag;
					if (j == NNLOVER[i] - 1)
					{
						for (j1 = 0; j1 < GP_2D; j1++)
						{
							//printf("Same_BDBJ_flag[%d]=%d\n",j1,Same_BDBJ_flag[j1]);
							if (Same_BDBJ_flag[j1] != 1)
							{
								printf("ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR\n");
							}
						}
						printf("-------------------------Check_BDBJ_flag[%d]=%d-------------------------\n",i ,Check_BDBJ_flag[i]);
						if (Check_BDBJ_flag[i] != GP_2D)
						{
							printf("ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR\n");
						}
					}

					//Valを求める
					//for (j1 = 0; j1 < No_Control_point_ON_ELEMENT[Element_patch[El_No_on_mesh[tm][i]]]; j1++)
					for (j1 = 0; j1 < No_Control_point_ON_ELEMENT[Element_patch[NELOVER[i][j]]]; j1++)
					{
						for (j2 = 0; j2 < DIMENSION; j2++)
						{
							a = Index_Dof[Controlpoint_of_Element[NELOVER[i][j]][j1] * DIMENSION + j2];
							if (a >= 0)
							{
								//for (k1 = 0; k1 < No_Control_point_ON_ELEMENT[Element_patch[El_No_on_mesh[tm][i]]]; k1++)
								for (k1 = 0; k1 < No_Control_point_ON_ELEMENT[Element_patch[i]]; k1++)
								{
									for (k2 = 0; k2 < DIMENSION; k2++)
									{
										b = Index_Dof[Controlpoint_of_Element[i][k1] * DIMENSION + k2];
										if (b >= 0 && b >= a)
										{
											for (l = K_Whole_Ptr[a]; l < K_Whole_Ptr[a + 1]; l++)
											{
												if (K_Whole_Col[l] == b)
												{
													K_Whole_Val[l] += coupled_K_EL[j1 * DIMENSION + j2][k1 * DIMENSION + k2];
													//printf("coupled_K_Whole_Val[%d]=%le\n",l,K_Whole_Val[l]);
													break;
												}
											}
										}
									}
								}
							}
						}
					}

				}
			}
		}
		//re++;
		//}
	}

}


///////////////////////////////////////////////////////////////////////////
/////////////////////連立1次方程式の解法
/////////////////////////////////////////////////////////////////////
//分布荷重の等価節点力を足す
void Add_Equivalent_Nodal_Forec_to_F_Vec(int Total_Control_Point)
{
	int i, j, index;
	for (j = 0; j < DIMENSION; j++)
    {
		for (i = 0; i < Total_Control_Point; i++)
		{
			index = Index_Dof[i * DIMENSION + j];
			if (index >= 0)
			{
				rhs_vec[index] += Equivalent_Nodal_Force[i][j];
				//printf("i = %d index = %d rhs_vec[index] = %f\n", i, index, rhs_vec[index]);

                //for s-IGA
				/*
                //結合剛性行列の右辺の組み立て
                rhs_overlaid_vec[index+K_Whole_Size_array_to_Now]=rhs_vec[index];
		        printf("i = %d index = %d rhs_overlaid_vec[index] = %f\n"
                       , i, index+K_Whole_Size_array_to_Now
                       , rhs_overlaid_vec[index+K_Whole_Size_array_to_Now]);
				*/
			}
			//printf("i = %d  j = %d  Equivalent_Nodal_Force[i][j] = %f\n",i, j, Equivalent_Nodal_Force[i][j]);
		}
    }
}


//荷重の行列を作る
void Make_F_Vec(int Total_Load, int Load_Node_Dir[MAX_N_LOAD][2], double Value_of_Load[MAX_N_LOAD], int K_Whole_Size)
{
	int i, index;
	for (i = 0; i < K_Whole_Size; i++)
		rhs_vec[i] = 0.0;
	for (i = 0; i < Total_Load; i++)
	{
		index = Index_Dof[Load_Node_Dir[i][0] * DIMENSION + Load_Node_Dir[i][1]];
		if (index >= 0)
			rhs_vec[index] += Value_of_Load[i];
	}
}


//強制変位対策
void Make_F_Vec_disp_const(int Mesh_No, int Total_Constraint,
                           int Constraint_Node_Dir[MAX_N_CONSTRAINT][2],
                           double Value_of_Constraint[MAX_N_CONSTRAINT],
                           double E, double nu, int DM)
{
	int ie, idir, inode, jdir, jnode, kk_const;
	int ii, iii, b, bb, jj, j1, j2, ii_local, jj_local;
	int iee;

	int i;

	Make_gauss_array(0);

	// for (ie = 0; ie < real_Total_Element; ie++)
	for (ie = 0; ie < real_Total_Element_to_mesh[Total_mesh]; ie++)
	{
		i = real_element[ie];

		// double X[No_Control_point_ON_ELEMENT[Element_patch[real_El_No_on_mesh[Mesh_No][ie]]]][DIMENSION];
		
		KIEL_SIZE = No_Control_point_ON_ELEMENT[Element_patch[i]] * DIMENSION;
		double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], K_EL[MAX_KIEL_SIZE][MAX_KIEL_SIZE];
		
		iii = 0;
		for (idir = 0; idir < DIMENSION; idir++)
		{
			for (inode = 0; inode < No_Control_point_ON_ELEMENT[Element_patch[i]]; inode++)
			{
				//b = Index_Dof[Controlpoint_of_Element[ie ][inode] * DIMENSION + idir];
				b = Index_Dof[Controlpoint_of_Element[i][inode] * DIMENSION + idir];
				if (b < 0)
					iii++;
			}
		}
		//printf("iii;%d\n",iii);
		if (iii > 0)
		{
			for (j1 = 0; j1 < No_Control_point_ON_ELEMENT[Element_patch[i]]; j1++)
			{
				for (j2 = 0; j2 < DIMENSION; j2++)
				{
					//X[j1][j2] = Node_Coordinate[Controlpoint_of_Element[ie ][j1]][j2];
					X[j1][j2] = Node_Coordinate[Controlpoint_of_Element[i][j1]][j2];
				} //end for j2
			}	 //end for j1
			iee = i;
			Make_K_EL(iee, X, K_EL, E, nu, DM);
			for (idir = 0; idir < DIMENSION; idir++)
			{
				for (inode = 0; inode < No_Control_point_ON_ELEMENT[Element_patch[i]]; inode++)
				{
					//ii = Controlpoint_of_Element[ie ][inode] * DIMENSION + idir;
					ii = Controlpoint_of_Element[real_El_No_on_mesh[Mesh_No][ie]][inode] * DIMENSION + idir;
					b = Index_Dof[ii];
					if (b >= 0)
					{
						ii_local = inode * DIMENSION + idir;
						for (jdir = 0; jdir < DIMENSION; jdir++)
						{
							for (jnode = 0; jnode < No_Control_point_ON_ELEMENT[Element_patch[i]]; jnode++)
							{
								//jj = Controlpoint_of_Element[ie ][jnode] * DIMENSION + jdir;
								jj = Controlpoint_of_Element[i][jnode] * DIMENSION + jdir;
								bb = Index_Dof[jj];
								if (bb < 0)
								{
									jj_local = jnode * DIMENSION + jdir; //printf("%d,%d\n",ie,jnode);
									for (kk_const = 0; kk_const < Total_Constraint; kk_const++)
									{
										//if (Controlpoint_of_Element[ie ][jnode] == Constraint_Node_Dir[kk_const][0] && jdir == Constraint_Node_Dir[kk_const][1])
										if (Controlpoint_of_Element[i][jnode] == Constraint_Node_Dir[kk_const][0] && jdir == Constraint_Node_Dir[kk_const][1])
										{
											rhs_vec[b] -= K_EL[ii_local][jj_local] * Value_of_Constraint[kk_const]; //if(kk_const >= 28){printf("%d , %d ,%16.15e\n",ii_local, jj_local ,  K_EL[ii_local][jj_local]);}
										}																			//end if Controlpoint_of_Element[ie][jnode]
									}																				//end for kk_const
								}																					//end if bb
							}																						//end for jnode
						}																							//end for jdir
					}																								//end if b>=0
				}																									//end for inode
			}																										//end for idir
		}																											//end if iii>0
	}																												//end for ie
} //end


void mat_vec_crs(double vec_result[], double vec[], const int ndof)
{
	int i, j, icount = 0;
	/* zero clear */

	for (i = 0; i < ndof; i++)
		vec_result[i] = 0;
	for (i = 0; i < ndof; i++)
	{
		//for (j = K_Whole_Ptr[i]; j < K_Whole_Ptr[i + 1]; j++)
        for (j = K_Whole_Ptr[i]; j < K_Whole_Ptr[i + 1]; j++)
		{
			vec_result[i] += K_Whole_Val[icount] * vec[K_Whole_Col[j]];
			//vec_result[i] += K_Whole_overlaid_Val[icount] * vec[K_Whole_overlaid_Col[j]];
            //if (i != K_Whole_overlaid_Col[j])
            if (i != K_Whole_Col[j])
				vec_result[K_Whole_Col[j]] += K_Whole_Val[icount] * vec[i];
				//vec_result[K_Whole_overlaid_Col[j]] += K_Whole_overlaid_Val[icount] * vec[i];
			icount++;
		}
	}
}


double inner_product(int ndof, double vec1[], double vec2[])
{
	double rrr = 0.0;
	int i;
	for (i = 0; i < ndof; i++)
	{
		rrr += vec1[i] * vec2[i];
		//printf("vec1[%d]=%f vec2[%d]=%f\n",i,vec1[i],i,vec2[i]); /*-nan 10/23*/
	}
	return (rrr);
}


int check_conv_CG(int ndof, double alphak, double pp[], double eps, int itr)
{
	double rrr1 = 0.0, rrr2 = 0.0, rrr3;
	int i, istop = 0;
	/* Checking the convergence of the CG solver */
	/* istop =0; Not converged, istop = 1; converged */
	printf("ndof=%d alphak= %15e\t", ndof, alphak);
	for (i = 0; i < ndof; i++)
	{
		rrr1 += pp[i] * pp[i];
		rrr2 += sol_vec[i] * sol_vec[i];
		//printf("pp[%d]=%f sol_vec[%d]=%f\n",i,pp[i],i,sol_vec[i]); /*-nan 10/23*/
	}
	rrr3 = fabs(alphak) * sqrt(rrr1 / rrr2);
	printf("Iteration# = %d  residual = %15e (%15e)\n", itr, rrr3, eps);
	if (rrr3 < eps)
		istop = 1;
	/* Temporaty Oct. 10, 2017 by H.Okada */
	//if(itr < 100) istop=0;
	//printf("Iteration# = %d  residual = %15e (%15e)\n",itr, rrr3, eps);
	return (istop);
}


void Diag_Scaling_CG_pre(int ndof, int flag_operation)
{
	int i, j;
	int icount = 0;
	/* flag_opertion = 0: Preprocess to the CG solver
			A <-- Dt A D  and b <-- Dt b */
	/* flag_operation = 1: Post process to the CG solver
			b <-- Dt b  */
    printf("ndof=%d\n",ndof);
	if (flag_operation == 0)
	{
		diag_scaling[0] = 1.0 / sqrt(K_Whole_Val[0]);
        //diag_scaling[0] = 1.0 / sqrt(K_Whole_overlaid_Val[0]);
		/* diag_scaling[0] = 1.0; */
		for (i = 1; i < ndof; i++)
		{
			//printf("%d %le\n",K_Whole_Ptr[i], K_Whole_Val[K_Whole_Ptr[i]]);
			diag_scaling[i] = 1.0 / sqrt(K_Whole_Val[K_Whole_Ptr[i]]);
			printf("diag=%le\n",diag_scaling[i]);
			printf("K_Whole_Val[%d] = %.16e\n",K_Whole_Ptr[i],K_Whole_Val[K_Whole_Ptr[i]]);
            printf("sqrt=%le\n",sqrt(K_Whole_Val[K_Whole_Ptr[i]]) );
			//diag_scaling[i]
            //    = 1.0 / sqrt(K_Whole_overlaid_Val[K_Whole_overlaid_Ptr[i]]);
			/* diag_scaling[i] = 1.0; */
		}
		for (i = 0; i < ndof; i++)
		{
			for (j = K_Whole_Ptr[i]; j < K_Whole_Ptr[i + 1]; j++)
            //for (j = K_Whole_overlaid_Ptr[i]; j < K_Whole_overlaid_Ptr[i + 1]; j++)
			{
				//printf("Check scling icount=%d i=%d K_Whole_Col[%d] = %d\n",icount,i,j,K_Whole_Col[j]);
				K_Whole_Val[icount] = K_Whole_Val[icount] * diag_scaling[i] * diag_scaling[K_Whole_Col[j]];
                //K_Whole_overlaid_Val[icount]
                //    = K_Whole_overlaid_Val[icount] * diag_scaling[i]
                //    * diag_scaling[K_Whole_overlaid_Col[j]];
				//printf("K_Whole_Val = %f\n",K_Whole_Val[icount]);
				icount++;
			}
			printf("rhs_vec_before[%d]:%le diag_scaling[%d]:%le\n"
                    , i, rhs_vec[i], i, diag_scaling[i]);
			//rhs_overlaid_vec[i] = rhs_overlaid_vec[i] * diag_scaling[i];
			rhs_vec[i] = rhs_vec[i] * diag_scaling[i];
			//printf("rhs_vec[%d]:%le\n",i,rhs_vec[i]);
		}
	}
	if (flag_operation == 1)
		for (i = 0; i < ndof; i++)
		{
			//printf("solvec[%d] = %f\n",i, sol_vec[i]);
			sol_vec[i] = sol_vec[i] * diag_scaling[i];
		}

	printf("\nqq\n");
}


void CG_Solver(int ndof, int max_itr, double eps, int flag_ini_val)
{
	static double gg[MAX_K_WHOLE_SIZE], dd[MAX_K_WHOLE_SIZE], pp[MAX_K_WHOLE_SIZE];
	static double qqq, ppp, rrr;
	static double alphak, betak;
	int i;
	int itr;
	int ii, istop;

	/* Program to solve linear equations by using the CG method */
	if (flag_ini_val == 0)
		for (i = 0; i < ndof; i++)
			sol_vec[i] = 0.0;
	/* Initializing the solution vector if it were not given */
	mat_vec_crs(dd, sol_vec, ndof);
	for (i = 0; i < ndof; i++)
	{
		//gg[i] = rhs_overlaid_vec[i] - dd[i];
		gg[i] = rhs_vec[i] - dd[i];
		//printf("rhs_vec[%d]=%f dd[%d]=%f\n", i, rhs_overlaid_vec[i], i, dd[i]);
		printf("rhs_vec[%d]=%f dd[%d]=%f\n", i, rhs_vec[i], i, dd[i]);
		pp[i] = gg[i];
	}

	printf("\nrr");

	for (itr = 0; itr < max_itr; itr++)
	{
		ppp = inner_product(ndof, gg, gg);
		mat_vec_crs(dd, pp, ndof);
		rrr = inner_product(ndof, dd, pp);
		alphak = ppp / rrr;
		//printf("ppp=%f rrr=%f\n", ppp, rrr); /*ppp,rrrも-nan,10/22*/
		//printf("i=%d",i);
		for (ii = 0; ii < ndof; ii++)
		{
			sol_vec[ii] += alphak * pp[ii];
			gg[ii] -= alphak * dd[ii];
		}
		qqq = inner_product(ndof, gg, dd);
		betak = qqq / rrr;
		for (ii = 0; ii < ndof; ii++)
			pp[ii] = gg[ii] - betak * pp[ii];
		istop = check_conv_CG(ndof, alphak, pp, eps, itr);
		if (istop == 1)
			break;
	}

	printf("\nss");
}


void Make_M(double *M, int *M_Ptr, int *M_Col, int ndof)
{
	int i, j;
	int ndof_glo = 0;

	// グローバルパッチのdofを求める
	for (i = 0; i < Total_Control_Point_on_mesh[0] * DIMENSION; i++)
	{
		if (Index_Dof[i] != ERROR)
		{
			ndof_glo++;
		}
	}
	printf("ndof		%d\n", ndof);
	printf("ndof_glo	%d\n", ndof_glo);

	int counter = 0;
	
	// M = [[K^G, 0], [0, K^L]] を作成
	M_Ptr[0] = 0;
	for (i = 0; i < ndof; i++)
	{
		M_Ptr[i + 1] = M_Ptr[i];

        for (j = K_Whole_Ptr[i]; j < K_Whole_Ptr[i + 1]; j++)
		{
			if (i < ndof_glo && K_Whole_Col[j] < ndof_glo)
			{
				M[counter] = K_Whole_Val[j];
				M_Col[counter] = K_Whole_Col[j];
				counter++;
				M_Ptr[i + 1]++;
			}
			else if (i >= ndof_glo)
			{
				M[counter] = K_Whole_Val[j];
				M_Col[counter] = K_Whole_Col[j];
				counter++;
				M_Ptr[i + 1]++;
			}
        }
    }
}


void M_mat_vec_crs(double *M, int *M_Ptr, int *M_Col, double *vec_result, double *vec, const int ndof)
{
	int i, j, icount = 0;

	for (i = 0; i < ndof; i++)
		vec_result[i] = 0;
	for (i = 0; i < ndof; i++)
	{
        for (j = M_Ptr[i]; j < M_Ptr[i + 1]; j++)
		{
			vec_result[i] += M[icount] * vec[M_Col[j]];
            if (i != M_Col[j])
				vec_result[M_Col[j]] += M[icount] * vec[i];
			icount++;
		}
	}
}


int M_check_conv_CG(int ndof, double alphak, double *pp, double eps, double *solution_vec)
{
	double rrr1 = 0.0, rrr2 = 0.0, rrr3;
	int i, istop = 0;
	for (i = 0; i < ndof; i++)
	{
		rrr1 += pp[i] * pp[i];
		rrr2 += solution_vec[i] * solution_vec[i];
	}
	rrr3 = fabs(alphak) * sqrt(rrr1 / rrr2);
	if (rrr3 < eps)
		istop = 1;
	return (istop);
}


void CG(int ndof, double *solution_vec, double *M, int *M_Ptr, int *M_Col, double *right_vec)
{
	int i/*, j*/;
	// int icount = 0;

	// 対角スケーリング 前処理
	// diag_scaling[0] = 1.0 / sqrt(M[0]);
	// for (i = 1; i < ndof; i++)
	// {
	// 	diag_scaling[i] = 1.0 / sqrt(M[M_Ptr[i]]);
	// }
	// for (i = 0; i < ndof; i++)
	// {
	// 	for (j = M_Ptr[i]; j < M_Ptr[i + 1]; j++)
	// 	{
	// 		M[icount] = M[icount] * diag_scaling[i] * diag_scaling[M_Col[j]];
	// 		icount++;
	// 	}
	// 	right_vec[i] = right_vec[i] * diag_scaling[i];
	// }

	// CG solver
	static double gg[MAX_K_WHOLE_SIZE], dd[MAX_K_WHOLE_SIZE], pp[MAX_K_WHOLE_SIZE];
	static double qqq, ppp, rrr;
	static double alphak, betak;
	int itr;
	int ii, istop;
	int max_itr = ndof;
	double eps = 1.0e-13;

	for (i = 0; i < ndof; i++)
	{
		solution_vec[i] = 0.0;
	}
	M_mat_vec_crs(M, M_Ptr, M_Col, dd, solution_vec, ndof);
	for (i = 0; i < ndof; i++)
	{
		gg[i] = right_vec[i] - dd[i];
		pp[i] = gg[i];
	}
	for (itr = 0; itr < max_itr; itr++)
	{
		ppp = inner_product(ndof, gg, gg);
		M_mat_vec_crs(M, M_Ptr, M_Col, dd, pp, ndof);
		rrr = inner_product(ndof, dd, pp);
		alphak = ppp / rrr;
		for (ii = 0; ii < ndof; ii++)
		{
			solution_vec[ii] += alphak * pp[ii];
			gg[ii] -= alphak * dd[ii];
		}
		qqq = inner_product(ndof, gg, dd);
		betak = qqq / rrr;
		for (ii = 0; ii < ndof; ii++)
			pp[ii] = gg[ii] - betak * pp[ii];
		istop = M_check_conv_CG(ndof, alphak, pp, eps, solution_vec);
		if (istop == 1)
			break;
	}
	printf("\titr %d\n", itr);

	// 対角スケーリング 後処理
	// for (i = 0; i < ndof; i++)
	// {
	// 	solution_vec[i] = solution_vec[i] * diag_scaling[i];
	// }
}


int RowCol_to_icount(int row, int col)
{
	for (int j = K_Whole_Ptr[row]; j < K_Whole_Ptr[row + 1]; j++)
	{
		if (K_Whole_Col[j] == col)
		{
			return j;
		}
		else if (K_Whole_Col[j] > col)
		{
			return -1;
		}
	}
	return -1;
}


// 前処理付共役勾配法により[K]{d}={f}を解く
void PCG_Solver(int ndof, int max_itetarion, double eps)
{
	int i, j, k;

	double *r = (double *)malloc(sizeof(double) * ndof);
	double *p = (double *)calloc(ndof, sizeof(double));
	double *y = (double *)malloc(sizeof(double) * ndof);
	double *r2 = (double *)calloc(ndof, sizeof(double));

	// 初期化
	for (i = 0; i < ndof; i++)
	sol_vec[i] = 0.0;

	// 前処理行列作成
	double *M = (double *)malloc(sizeof(double) * MAX_NON_ZERO);
	int *M_Ptr = (int *)malloc(sizeof(int) * MAX_K_WHOLE_SIZE + 1);
	int *M_Col = (int *)malloc(sizeof(int) * MAX_NON_ZERO);
	Make_M(M, M_Ptr, M_Col, ndof);

	// 第0近似解に対する残差の計算
	double *ax = (double *)calloc(ndof, sizeof(double));
    for (i = 0; i < ndof; i++)
	{
        for (j = K_Whole_Ptr[i]; j < K_Whole_Ptr[i + 1]; j++)
		{
            ax[i] += K_Whole_Val[j] * sol_vec[K_Whole_Col[j]];
			if (i != K_Whole_Col[j])
			{
				ax[K_Whole_Col[j]] += K_Whole_Val[j] * sol_vec[i];
			}
        }
    }
	for (i = 0; i < ndof; i++)
	{
		r[i] = rhs_vec[i] - ax[i];
	}
	free(ax);

	// 第0近似解に対する残差の計算
	// for (i = 0; i < ndof; i++)
	// {
	// 	r[i] = rhs_vec[i];
	// }

	// p_0 = (LDL^T)^-1 r_0 の計算 <- CG法で M = [[K^G, 0], [0, K^L]] とし，p_0 = (LDL^T)^-1 r_0 = M^-1 r_0
	CG(ndof, p, M, M_Ptr, M_Col, r);


    // double rr0 = inner_product(ndof, r, p), rr1;
	double rr0;
    double alpha, beta;

	double e = 0.0;
    for(k = 0; k < max_itetarion; k++)
	{
		// rr0 の計算
		rr0 = inner_product(ndof, r, p);

        // y = AP の計算
        for(i = 0; i < ndof; i++)
		{
			double *temp_array_K = (double *)calloc(ndof, sizeof(double));
			for (j = 0; j < ndof; j++)
			{
				int temp1;
				if (i <= j)
				{
					temp1 = RowCol_to_icount(i, j); // temp_array_K[i][j]
				}
				else if (i > j)
				{
					temp1 = RowCol_to_icount(j, i); // temp_array_K[i][j] = temp_array_K[j][i]
				}

				if (temp1 != -1)
				{
					temp_array_K[j] = K_Whole_Val[temp1];
				}
			}
            y[i] = inner_product(ndof, temp_array_K, p);
			free(temp_array_K);
        }
 
        // alpha = r*r/(P*AP)の計算
		double temp_scaler = inner_product(ndof, p, y);
        alpha = rr0 / temp_scaler;
		// printf("alpha %le\n", alpha);
 
        // 解x、残差rの更新
        for(i = 0; i < ndof; i++)
		{
            sol_vec[i] += alpha * p[i];
            r[i] -= alpha * y[i];
        }

		// (r*r)_(k+1)の計算
		CG(ndof, r2, M, M_Ptr, M_Col, r);
	
	
        // rr1 = inner_product(ndof, r, r2); // 旧
		// rr1 = inner_product(ndof, y, r2); // 新
		// printf("rr1 %le\n", rr1);
 
        // 収束判定 (||r||<=eps)
		// double rr1 = inner_product(ndof, y, r2);
        // e = sqrt(fabs(rr1));
        // if(e < eps)
		// {
        //     k++;
        //     break;
        // }

		// 収束判定 (CG法と同じ)
		double e1 = 0.0, e2 = 0.0;
		for (i = 0; i < ndof; i++)
		{
			e1 += p[i] * p[i];
			e2 += sol_vec[i] * sol_vec[i];
		}
        e = fabs(alpha) * sqrt(e1 / e2);
        if(e < eps)
		{
            k++;
            break;
        }
 
        // βの計算とPの更新
        // beta = rr1 / rr0; //旧
		// beta = - rr1 / temp_scaler; // 新
		beta = - inner_product(ndof, y, r2) / temp_scaler;

        for(i = 0; i < ndof; i++)
		{
            // p[i] = r2[i] - beta * p[i];
			p[i] = r2[i] + beta * p[i];
        }
		// printf("beta %le\n", beta);
 
        // (r*r)_(k+1)を次のステップのために確保しておく
        // rr0 = rr1;
	
		printf("itr %d\t", k);
		printf("eps %.15e", e);
		// if (rr1 < 0)
		// {
		// 	printf("\t rr1 < 0");
		// }
		printf("\n");
    }
 
    int max_itr_result = k;
    double eps_result = e;

	printf("\nndof = %d\n", ndof);
	printf("itr_result = %d\n", max_itr_result);
	printf("eps_result = %.15e\n", eps_result);

	free(r), free(p), free(y), free(r2);
	free(M), free(M_Ptr), free(M_Col);
}


////////////////////////////////////////////////////////////////////////
/////////////////基底関数
////////////////////////////////////////////////////////////////////////
//IGAの基底関数
void ShapeFunction1D(double Position_Data_param[DIMENSION], int j, int e)
{

	int ii;
	int p;

	//printf("shapefuc_Position_Data_param[%d]:%le\n", j,Position_Data_param[j]);

	for (ii = 0; ii < No_knot[Element_patch[e]][j]; ii++)
	{
		if (Position_Knots[Element_patch[e]][j][ii] == Position_Knots[Element_patch[e]][j][ii + 1])
		{
			Shape[j][ii][0] = 0.0;
		}
		else if (Position_Knots[Element_patch[e]][j][ii] != Position_Knots[Element_patch[e]][j][ii + 1] && Position_Knots[Element_patch[e]][j][ii] <= Position_Data_param[j] && Position_Data_param[j] < Position_Knots[Element_patch[e]][j][ii + 1])
		{
			Shape[j][ii][0] = 1.0;
		}
		else if (Position_Knots[Element_patch[e]][j][ii] != Position_Knots[Element_patch[e]][j][ii + 1] && Position_Knots[Element_patch[e]][j][ii + 1] == Position_Knots[Element_patch[e]][j][(No_knot[Element_patch[e]][j] - 1)] && Position_Knots[Element_patch[e]][j][ii] <= Position_Data_param[j] && Position_Data_param[j] <= Position_Knots[Element_patch[e]][j][ii + 1])
		{
			Shape[j][ii][0] = 1.0;
		}
		else
			Shape[j][ii][0] = 0.0;
		//printf("Shape[%d][%d][0]=%le   ",j,ii,Shape[j][ii][0]);
	}

	for (ii = 0; ii < No_knot[Element_patch[e]][j]; ii++)
	{
		for (p = 1; p <= Order[Element_patch[e]][j]; p++)
		{
			Shape[j][ii][p] = 0.0;
		}
	}
	double left_term, right_term;
	for (p = 1; p <= Order[Element_patch[e]][j]; p++)
	{
		for (ii = 0; ii < No_knot[Element_patch[e]][j]; ii++)
		{
			left_term = 0.0;
			right_term = 0.0;

			if ((Position_Data_param[j] - Position_Knots[Element_patch[e]][j][ii]) * Shape[j][ii][p - 1] == 0 && Position_Knots[Element_patch[e]][j][ii + p] - Position_Knots[Element_patch[e]][j][ii] == 0)
				left_term = 0.0;
			else
            {
				left_term = (Position_Data_param[j] - Position_Knots[Element_patch[e]][j][ii]) / (Position_Knots[Element_patch[e]][j][ii + p] - Position_Knots[Element_patch[e]][j][ii]) * Shape[j][ii][p - 1];
                //printf("else left\tleft_term;%le\n",left_term);
            }
			if ((Position_Knots[Element_patch[e]][j][ii + p + 1] - Position_Data_param[j]) * Shape[j][ii + 1][p - 1] == 0 && Position_Knots[Element_patch[e]][j][ii + p + 1] - Position_Knots[Element_patch[e]][j][ii + 1] == 0)
				right_term = 0.0;
			else
            {
				right_term = (Position_Knots[Element_patch[e]][j][ii + p + 1] - Position_Data_param[j]) / (Position_Knots[Element_patch[e]][j][ii + p + 1] - Position_Knots[Element_patch[e]][j][ii + 1]) * Shape[j][ii + 1][p - 1];
                //printf("else right\tright_term;%le\n",right_term);
            }
			Shape[j][ii][p] = left_term + right_term;
			//printf("Shape[%d][%d][%d]=%le\n",j,ii,p,Shape[j][ii][p]);
		}
	}
	//printf("order[%d]:%d\n",j,Order[Element_patch[e]][j] );
	double dleft_term, dright_term;
	for (ii = 0; ii < No_Control_point[Element_patch[e]][j] + 1; ii++)
	{
		//printf("No_Control_point[%d]=%d\n",j,No_Control_point[j] );
		dleft_term = 0.0;
		dright_term = 0.0;

		if (Order[Element_patch[e]][j] * Shape[j][ii][Order[Element_patch[e]][j] - 1] == 0 && Position_Knots[Element_patch[e]][j][ii + Order[Element_patch[e]][j]] - Position_Knots[Element_patch[e]][j][ii] == 0)
			dleft_term = 0.0;
		else
			dleft_term = Order[Element_patch[e]][j] / (Position_Knots[Element_patch[e]][j][ii + Order[Element_patch[e]][j]] - Position_Knots[Element_patch[e]][j][ii]) * Shape[j][ii][Order[Element_patch[e]][j] - 1];
		/*printf("test_Shape_left[%d][%d][%d]=%le\n", j,ii,Order[Element_patch[e]][j]-1,Shape[j][ii][Order[Element_patch[e]][j]-1]);
		printf("Position_Knots[Element_patch[e]][%d][%d]:%le\n", j,ii+Order[Element_patch[e]][j],Position_Knots[Element_patch[e]][j][ii+Order[Element_patch[e]][j]]);
		printf("Position_Knots[Element_patch[e]][%d][%d]:%le\n", j,ii,Position_Knots[Element_patch[e]][j][ii]);
		printf("dleft_term=%f\n",dleft_term );*/

		if (Order[Element_patch[e]][j] * Shape[j][ii + 1][Order[Element_patch[e]][j] - 1] == 0 && Position_Knots[Element_patch[e]][j][ii + Order[Element_patch[e]][j] + 1] - Position_Knots[Element_patch[e]][j][ii + 1] == 0)
			dright_term = 0.0;
		else
			dright_term = Order[Element_patch[e]][j] / (Position_Knots[Element_patch[e]][j][ii + Order[Element_patch[e]][j] + 1] - Position_Knots[Element_patch[e]][j][ii + 1]) * Shape[j][ii + 1][Order[Element_patch[e]][j] - 1];
		/*printf("test_Shape_right[%d][%d][%d]=%le\n", j,ii+1,Order[Element_patch[e]][j]-1,Shape[j][ii+1][Order[Element_patch[e]][j]-1]);
		printf("Position_Knots[%d][%d]:%le\n", j,ii+Order[Element_patch[e]][j]+1,Position_Knots[j][ii+Order[Element_patch[e]][j]+1]);
		printf("Position_Knots[%d][%d]:%le\n", j,ii+1,Position_Knots[j][ii+1]);
		printf("dright_term=%f\n",dright_term );*/

		dShape[j][ii] = dleft_term - dright_term;

		//printf("PP=%d\n",PP );

		//printf("dShape[%d][%d]= %f\n",j,ii,dShape[j][ii]);
	}
}


void ShapeFunc_from_paren(double Local_coord[DIMENSION], int j, int e)
{
	int i = 0;
	//printf("Local_coord[%d]:%le\n",j,Local_coord[j]);
	i = INC[Element_patch[e]][Controlpoint_of_Element[e][0]][j];
	//printf("El_No:%d\n",e );
	//printf("i:%d\n",i);
	//printf("Position_Knots[%d][%d][%d]:%le Position_Knots[%d][%d][%d]:%le\n"
    //        ,Element_patch[e],j,i+1,Position_Knots[Element_patch[e]][j][i+1]
    //        ,Element_patch[e],j,i,Position_Knots[Element_patch[e]][j][i]);
	Position_Data_param[j] = ((Position_Knots[Element_patch[e]][j][i + 1] - Position_Knots[Element_patch[e]][j][i]) * Local_coord[j] + (Position_Knots[Element_patch[e]][j][i + 1] + Position_Knots[Element_patch[e]][j][i])) / 2;
	//Position_Data_param[j] = ((Position_Knots[Element_patch[e]][j][i + 1] - Position_Knots[Element_patch[e]][j][i]) * Local_coord[j] + (Position_Knots[Element_patch[e]][j][i + 1] + Position_Knots[Element_patch[e]][j][i])) / 2;
	//printf("Position_Data_param[%d]:%le\n", j, Position_Data_param[j]);
}


double dShapeFunc_from_paren(int j, int e)
{
	int i;
	double dPosition_Data_param;

	i = INC[Element_patch[e]][Controlpoint_of_Element[e ][0]][j];
	//printf("El_No:%d\n",e );
	//printf("i:%d\n",i);
	//printf("Position_Knots[Element_patch[e]][%d][%d]:%le Position_Knots[%d][%d]:%le\n", j,i+1,Position_Knots[j][i+1],j,i,Position_Knots[j][i]);
	dPosition_Data_param = (Position_Knots[Element_patch[e]][j][i + 1] - Position_Knots[Element_patch[e]][j][i]) / 2;
	//printf("dPosition_Data_param:%le\t",dPosition_Data_param);
	return dPosition_Data_param;
}


double Shape_func(int I_No, double Local_coord[DIMENSION], int El_No)
{

	int i, j;
	double R;
	double weight_func;
	weight_func = 0.0;
	//shape_func[]={0.0};

	// MAX_N_NODEでshape_funcを初期化していたところを修正
	// // for (i = 0; i < MAX_N_NODE; i++)
	// {
	// 	shape_func[i] = 1.0;
	// }

	// for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++)
	// {
	// 	for (j = 0; j < DIMENSION; j++)
	// 	{
	// 		ShapeFunc_from_paren(Local_coord, j, El_No);
	// 		ShapeFunction1D(Position_Data_param, j, El_No);
	// 		shape_func[Controlpoint_of_Element[El_No][i]] *= Shape[j][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][j]][Order[Element_patch[El_No]][j]];
	// 		//shape_func[Controlpoint_of_Element[El_No][i]] *= Shape[j][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][j]][Order[Element_patch[El_No]][j]]; /*基底関数*/
	// 	    //printf("%d",shape_func[0][0]);
    //         //printf("Shape[%d][%d][%d]=%le\n",
    //         //        j,
    //         //        INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][j],
    //         //        Order[Element_patch[El_No]][j],
    //         //        Shape[j][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][j]][Order[Element_patch[El_No]][j]]);
	// 	}
    //     //printf("Controlpoint_of_element[%d][%d]=%d\n",
    //     //        El_No,i,Controlpoint_of_Element[El_No][i]);
    //     //printf("shape_func=%le\n",shape_func[Controlpoint_of_Element[El_No][i]]);
	// 	//printf("weight[%d];%f\n",Controlpoint_of_Element[El_No][i],Node_Coordinate[Controlpoint_of_Element[El_No][i]][DIMENSION]);
	// 	weight_func += shape_func[Controlpoint_of_Element[El_No][i]] * Node_Coordinate[Controlpoint_of_Element[El_No][i]][DIMENSION];
	// }

	int Total_CP = 0;

	for (i = 0; i < Total_mesh; i++)
	{
		Total_CP += Total_Control_Point_on_mesh[i];
	}
	for (i = 0; i < Total_CP; i++)
	{
		shape_func[i] = 1.0;
	}

	for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			ShapeFunc_from_paren(Local_coord, j, El_No);
			ShapeFunction1D(Position_Data_param, j, El_No);
			shape_func[Controlpoint_of_Element[El_No][i]] *= Shape[j][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][j]][Order[Element_patch[El_No]][j]];
		}
		weight_func += shape_func[Controlpoint_of_Element[El_No][i]] * Node_Coordinate[Controlpoint_of_Element[El_No][i]][DIMENSION];
	}
    //printf("weight_func_inShapefunc=%le\n",weight_func);
	if (I_No < No_Control_point_ON_ELEMENT[Element_patch[El_No]])
		R = shape_func[Controlpoint_of_Element[El_No][I_No]] * Node_Coordinate[Controlpoint_of_Element[El_No][I_No]][DIMENSION] / weight_func;

	else
		R = ERROR;
    //printf("R=%le\n",R);
	return R;
}


void NURBS_deriv(double Local_coord[DIMENSION], int El_No)
{
	double weight_func;
	//double shape_func[100][50];

	double dWeight_func1;
	double dWeight_func2;

	int i, j;
	//int ii;

	//for(ii = 0; ii < NN+1; ii++){
	//printf("NdShape3[%d]= %f\n",ii,dShape3[ii]);
	//}

	//for (ii = 0; ii < NN+1; ii++)printf("NdShape1[%d]= %f\n",ii,dShape1[ii]);
	//for (jj = 0; jj < MM+1; jj++)printf("NdShape2[%d]= %f\n",jj,dShape2[jj]);
	//printf("\n");

	// //for (i = 0; i < Total_Control_Point_to_mesh[Total_mesh]; i++)	//"Total_Control_Point_to_mesh[Total_mesh]"仮
	//"Total_Control_Point_to_mesh[Total_mesh]"仮
	// MAX_N_NODEでshape_funcを初期化していたところを修正
	// for (i = 0; i < MAX_N_NODE; i++)
	// {
	// 	shape_func[i] = 1.0;
	// }

	int Total_CP = 0;

	for (i = 0; i < Total_mesh; i++)
	{
		Total_CP += Total_Control_Point_on_mesh[i];
	}
	for (i = 0; i < Total_CP; i++)
	{
		shape_func[i] = 1.0;
	}

	weight_func = 0.0;

	dWeight_func1 = 0.0;
	dWeight_func2 = 0.0;

	for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			//printf("b_before_shape_func:%lf\n",
			//		shape_func[Controlpoint_of_Element[El_No][i]]);
			ShapeFunc_from_paren(Local_coord, j, El_No);
			ShapeFunction1D(Position_Data_param, j, El_No);
			//printf("before_shape_func:%lf\n",
			//		shape_func[Controlpoint_of_Element[El_No][i]]);
			shape_func[Controlpoint_of_Element[El_No][i]] *= Shape[j][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][j]][Order[Element_patch[El_No]][j]];
			//printf("Shape[%d][%d][%d]:%le\n",
            //        j,
            //        INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][j],
            //        Order[Element_patch[El_No]][j],
            //        Shape[j][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][j]][Order[Element_patch[El_No]][j]]);
			//printf("after_shape_func:%lf\n",
			//		shape_func[Controlpoint_of_Element[El_No][i]]);
		}
		weight_func += shape_func[Controlpoint_of_Element[El_No][i]] * Node_Coordinate[Controlpoint_of_Element[El_No][i]][DIMENSION];
		//printf("Node_Coordinate[%d][%d]:%le\n", Controlpoint_of_Element[El_No][i],DIMENSION,Node_Coordinate[Controlpoint_of_Element[El_No][i]][DIMENSION]);
		//printf("shape_func[%d]:%le\n",Controlpoint_of_Element[El_No][i],shape_func[Controlpoint_of_Element[El_No][i]]);
		//printf("weight_func:%le\n", weight_func);
	}
	//printf("weight_func:%le\n", weight_func);
	//for(jj=0;jj<NN;jj++) for(kk=0;kk<MM;kk++)shape_func[jj][kk] = Shape1[jj][PP]*Shape2[kk][QQ]*Node_Coordinate[MM*jj+kk][2]/weight_func;
	for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++)
	{
		dWeight_func1
            += dShape[0][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][0]]
            * Shape[1][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][1]][Order[Element_patch[El_No]][1]]
            * Node_Coordinate[Controlpoint_of_Element[El_No][i]][DIMENSION];
		dWeight_func2
            += Shape[0][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][0]][Order[Element_patch[El_No]][0]]
            * dShape[1][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][1]]
            * Node_Coordinate[Controlpoint_of_Element[El_No][i]][DIMENSION];
	}
	//printf("dWeight_func1:%le dWeight_func2:%le\n",dWeight_func1,dWeight_func2);
	for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++)
	{
		dShape_func1[Controlpoint_of_Element[El_No][i]]
            = Node_Coordinate[Controlpoint_of_Element[El_No][i]][DIMENSION]
            * (weight_func * dShape[0][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][0]]
            * Shape[1][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][1]][Order[Element_patch[El_No]][1]] - dWeight_func1 * shape_func[Controlpoint_of_Element[El_No][i]]) / (weight_func * weight_func);
		dShape_func2[Controlpoint_of_Element[El_No][i]]
            = Node_Coordinate[Controlpoint_of_Element[El_No][i]][DIMENSION]
			* (weight_func * Shape[0][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][0]][Order[Element_patch[El_No]][0]] * dShape[1][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][1]] - dWeight_func2 * shape_func[Controlpoint_of_Element[El_No][i]])
			/ (weight_func * weight_func);
		//printf("NURBS_deriv;Controlpoint_of_Element[%d][%d]:%d\n",El_No,i,Controlpoint_of_Element[El_No][i]);
		//printf("dShape_func1[%d]:%le\n",Controlpoint_of_Element[El_No][i],dShape_func1[Controlpoint_of_Element[El_No][i]]);
		//printf("dShape_func2[%d]:%le\n",Controlpoint_of_Element[El_No][i],dShape_func2[Controlpoint_of_Element[El_No][i]]);
	}
}


double dShape_func(int I_No, int xez, double Local_coord[DIMENSION], int El_No)
{
	double dR;

	//printf("El_No=%d\n",El_No);

	NURBS_deriv(Local_coord, El_No);

	if (xez != 0 && xez != 1)
		dR = ERROR;

	else if (I_No < No_Control_point_ON_ELEMENT[Element_patch[El_No]])
	{
		if (xez == 0)
		{
			//dR = dShape_func1[Controlpoint_of_Element[El_No][I_No]] * dShapeFunc_from_paren(xez, El_No);
			dR =
                dShape_func1[Controlpoint_of_Element[El_No][I_No]]
                * dShapeFunc_from_paren(xez, El_No);
			//printf("dShape_func1[%d]:%le\n",Controlpoint_of_Element[El_No][I_No],dShape_func1[Controlpoint_of_Element[El_No][I_No]]);
		}
		else if (xez == 1)
		{
			//dR = dShape_func2[Controlpoint_of_Element[El_No][I_No]] * dShapeFunc_from_paren(xez, El_No);
			dR =
                dShape_func2[Controlpoint_of_Element[El_No][I_No]]
                * dShapeFunc_from_paren(xez, El_No);
			//printf("dShape_func2[%d]:%le\n",Controlpoint_of_Element[El_No][I_No],dShape_func2[Controlpoint_of_Element[El_No][I_No]]);
		}
		//printf("dR%lf\n",dR);
	}

	else
		dR = ERROR;
	//printf("dR:%le\n",dR);

	//printf("I_No=%d xez=%d dR=%le\t", I_No, xez, dR );
	//printf("Controlpoint_of_Element[%d][%d]:%d\t",El_No,I_No,Controlpoint_of_Element[El_No][I_No]);
	//printf("dShape_func1:%le\t",dShape_func1[Controlpoint_of_Element[El_No][I_No]]);
	//printf("dShape_func2:%le\n",dShape_func2[Controlpoint_of_Element[El_No][I_No]]);

	/*for (i = 0; i < DIMENSION; i++) {
	printf("Local_coord[%d]=%lf\n",i,Local_coord[i] );
}*/
	//printf("dR:%le\t", dR);
	return dR;
}
/*
//形状関数
double N_Quad_4(int I_No, double Local_coord[DIMENSION] )
{
	double N;
	if(I_No==0) N= (1.0+Local_coord[0])*(1.0-Local_coord[1])/4.0;
	else if(I_No==1) N = (1.0+Local_coord[0])*(1.0+Local_coord[1])/4.0;
	else if(I_No==2) N = (1.0-Local_coord[0])*(1.0+Local_coord[1])/4.0;
	else if(I_No==3) N = (1.0-Local_coord[0])*(1.0-Local_coord[1])/4.0;
	else N = ERROR;
	return N;
}

//形状関数の偏微分（I_No:節点番号 xez:偏微分の分母部分0ξ1η2ζ）
double dN_Quad_4(int I_No, double Local_coord[DIMENSION], int xez)
{
	double dN;
	if(xez!=0 && xez!=1) dN = ERROR;

	else if(I_No==0)
		{
		if( xez == 0 )	    dN = (1.0-Local_coord[1])/4.0;
		else if( xez == 1 ) dN = (1.0+Local_coord[0])*(-1)/4.0;
		}

	else if(I_No==1)
		{
		if( xez == 0 )      dN = (1.0+Local_coord[1])/4.0;
		else if( xez == 1 ) dN = (1.0+Local_coord[0])/4.0;
		}

	else if(I_No==2)
		{
		if( xez == 0 )      dN = (1.0+Local_coord[1])*(-1)/4.0;
		else if( xez == 1 ) dN = (1.0-Local_coord[0])/4.0;
		}

	else if(I_No==3)
		{
		if( xez == 0 )      dN = (1.0-Local_coord[1])*(-1)/4.0;
		else if( xez == 1 ) dN = (1.0-Local_coord[0])*(-1)/4.0;
		}

	else dN = ERROR;

	return dN;
}
*/


//逆行列を元の行列に代入
double InverseMatrix_2D(double M[2][2])
{
	int i, j;
	double a[2][2];
	double det = M[0][0] * M[1][1] - M[0][1] * M[1][0];

	if (det == 0)
		return ERROR;

	for (i = 0; i < 2; i++)
	{
		for (j = 0; j < 2; j++)
			a[i][j] = M[i][j];
	}
	M[0][0] = a[1][1] / det;
	M[0][1] = a[0][1] * (-1) / det;
	M[1][0] = a[1][0] * (-1) / det;
	M[1][1] = a[0][0] / det;
	//printf("det;%le\n", det);
	return det;
}


double InverseMatrix_3X3(double M[3][3])
{
	int i, j;
	double a[3][3];
	double det = M[0][0] * M[1][1] * M[2][2] + M[1][0] * M[2][1] * M[0][2] + M[2][0] * M[0][1] * M[1][2] - M[0][0] * M[2][1] * M[1][2] - M[2][0] * M[1][1] * M[0][2] - M[1][0] * M[0][1] * M[2][2];

	if (det == 0)
		return ERROR;

	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
			a[i][j] = M[i][j];
	}
	M[0][0] = (a[1][1] * a[2][2] - a[1][2] * a[2][1]) / det;
	M[0][1] = (a[0][2] * a[2][1] - a[0][1] * a[2][2]) / det;
	M[0][2] = (a[0][1] * a[1][2] - a[0][2] * a[1][1]) / det;
	M[1][0] = (a[1][2] * a[2][0] - a[1][0] * a[2][2]) / det;
	M[1][1] = (a[0][0] * a[2][2] - a[0][2] * a[2][0]) / det;
	M[1][2] = (a[0][2] * a[1][0] - a[0][0] * a[1][2]) / det;
	M[2][0] = (a[1][0] * a[2][1] - a[1][1] * a[2][0]) / det;
	M[2][1] = (a[0][1] * a[2][0] - a[0][0] * a[2][1]) / det;
	M[2][2] = (a[0][0] * a[1][1] - a[0][1] * a[1][0]) / det;
	//printf("det;%le\n", det);
	return det;
}


////////////////////////////////////////////////////////////////////////
/////////////////Newton-Raphson法
////////////////////////////////////////////////////////////////////////
//from NURBSviewer
double BasisFunc(double *knot_vec, int knot_index, int order, double xi,
                 double *output, double *d_output) {
	int p, j;
	double sum1 = 0.0;
	double sum2 = 0.0;
	double temp_basis[order][order];
	(*output) = 0.0;
	(*d_output) = 0.0;

	if ( knot_vec[knot_index] <= xi && knot_vec[knot_index + order + 1] >= xi ) {
		for (j = 0; j <= order; j++) {
			if ( (knot_vec[knot_index + j] <= xi) && (xi <= knot_vec[knot_index + j + 1]) ) {
				temp_basis[j][0] = 1.0;
			} else {
				temp_basis[j][0] = 0.0;
			}
		}

		if ( order > 0 ) {
			for (p = 1; p <= order; p++) {
				for (j = 0; j <= order - p; j++) {
					sum1 = 0.0;
					sum2 = 0.0;
					if ( (knot_vec[knot_index + j + p] - knot_vec[knot_index + j]) != 0.0) {
						sum1 = ( xi - knot_vec[knot_index + j] )
						       / (knot_vec[knot_index + j + p] - knot_vec[knot_index + j])
						       * temp_basis[j][p - 1];
					}
					if ( (knot_vec[knot_index + j + p + 1] - knot_vec[knot_index + j + 1]) != 0.0) {
						sum2 = (knot_vec[knot_index + j + p + 1] - xi)
						       / (knot_vec[knot_index + j + p + 1] - knot_vec[knot_index + j + 1])
						       * temp_basis[j + 1][p - 1];
					}
					temp_basis[j][p] = sum1 + sum2;
				}
			}
			sum1 = 0.0;
			sum2 = 0.0;
			if ( (knot_vec[knot_index + order] - knot_vec[knot_index]) != 0.0) {
				sum1 = order
				       / (knot_vec[knot_index + order] - knot_vec[knot_index])
				       * temp_basis[0][order - 1];
			}
			if ( (knot_vec[knot_index + order + 1] - knot_vec[knot_index + 1]) != 0.0) {
				sum2 = order
				       / (knot_vec[knot_index + order + 1] - knot_vec[knot_index + 1])
				       * temp_basis[1][order - 1];
			}
		}
		(*output) = temp_basis[0][order];
		(*d_output) = sum1 - sum2;
	}
	return (*output);
}


double rBasisFunc(double *knot_vec, int knot_index,
				  int order, double xi,
				  double *output, double *d_output)
{
	int p, j;
	double sum1 = 0.0;
	double sum2 = 0.0;
	double temp_basis[order][order];
	(*output) = 0.0;
	(*d_output) = 0.0;


	if ( knot_vec[knot_index] <= xi && xi <= knot_vec[knot_index + order + 1] ) {
		if (knot_index == 0) {
			for (j = 0; j <= order; j++) {
				if ( (knot_vec[j] <= xi) && (xi <= knot_vec[j + 1]) ) {
					temp_basis[j][0] = 1.0;
				} else {
					temp_basis[j][0] = 0.0;
				}
			}
		} else {
			for (j = 0; j <= order; j++) {
				if ( (knot_vec[knot_index + j] < xi) && (xi <= knot_vec[knot_index + j + 1]) ) {
					temp_basis[j][0] = 1.0;
				} else {
					temp_basis[j][0] = 0.0;
				}
			}
		}

		if ( order > 0 ) {
			for (p = 1; p <= order; p++) {
				for (j = 0; j <= order - p; j++) {
					sum1 = 0.0;
					sum2 = 0.0;
					if ( (knot_vec[knot_index + j + p] - knot_vec[knot_index + j]) != 0.0) {
						sum1 = ( xi - knot_vec[knot_index + j] )
							   / (knot_vec[knot_index + j + p] - knot_vec[knot_index + j])
							   * temp_basis[j][p - 1];
					}
					if ( (knot_vec[knot_index + j + p + 1] - knot_vec[knot_index + j + 1]) != 0.0) {
						sum2 = (knot_vec[knot_index + j + p + 1] - xi)
							   / (knot_vec[knot_index + j + p + 1] - knot_vec[knot_index + j + 1])
							   * temp_basis[j + 1][p - 1];
					}
					temp_basis[j][p] = sum1 + sum2;
				}
			}
			sum1 = 0.0;
			sum2 = 0.0;

			// for (int temp_i = 0; temp_i < No_knot[0][0]; temp_i++)
			// {
			// 	printf("knot_vec[%d] = %f\n", temp_i, knot_vec[temp_i]);
			// }

			if ( (knot_vec[knot_index + order] - knot_vec[knot_index]) != 0.0) {
				sum1 = order
					   / (knot_vec[knot_index + order] - knot_vec[knot_index])
					   * temp_basis[0][order - 1];
			}
			if ( (knot_vec[knot_index + order + 1] - knot_vec[knot_index + 1]) != 0.0) {
				sum2 = order
					   / (knot_vec[knot_index + order + 1] - knot_vec[knot_index + 1])
					   * temp_basis[1][order - 1];
			}
		}
		(*output) = temp_basis[0][order];
		(*d_output) = sum1 - sum2;
	}
	return (*output);
}


double lBasisFunc(double *knot_vec, int knot_index,
				  int cntl_p_n, int order, double xi,
				  double *output, double *d_output)
{
	int p, j;
	double sum1 = 0.0;
	double sum2 = 0.0;
	double temp_basis[order][order];
	(*output) = 0.0;
	(*d_output) = 0.0;


	if ( knot_vec[knot_index] <= xi && xi <= knot_vec[knot_index + order + 1] ) {
		if (knot_index == cntl_p_n - 1) {
			for (j = 0; j <= order; j++) {
				if ( (knot_vec[cntl_p_n - 1 + j] <= xi) && (xi <= knot_vec[cntl_p_n + j]) ) {
					temp_basis[j][0] = 1.0;
				} else {
					temp_basis[j][0] = 0.0;
				}
			}
		} else {
			for (j = 0; j <= order; j++) {
				if ( (knot_vec[knot_index + j] <= xi) && (xi < knot_vec[knot_index + j + 1]) ) {
					temp_basis[j][0] = 1.0;
				} else {
					temp_basis[j][0] = 0.0;
				}
			}
		}

		if ( order > 0 ) {
			for (p = 1; p <= order; p++) {
				for (j = 0; j <= order - p; j++) {
					sum1 = 0.0;
					sum2 = 0.0;
					if ( (knot_vec[knot_index + j + p] - knot_vec[knot_index + j]) != 0.0) {
						sum1 = ( xi - knot_vec[knot_index + j] )
							   / (knot_vec[knot_index + j + p] - knot_vec[knot_index + j])
							   * temp_basis[j][p - 1];
					}
					if ( (knot_vec[knot_index + j + p + 1] - knot_vec[knot_index + j + 1]) != 0.0) {
						sum2 = (knot_vec[knot_index + j + p + 1] - xi)
							   / (knot_vec[knot_index + j + p + 1] - knot_vec[knot_index + j + 1])
							   * temp_basis[j + 1][p - 1];
					}
					temp_basis[j][p] = sum1 + sum2;
				}
			}
			sum1 = 0.0;
			sum2 = 0.0;
			if ( (knot_vec[knot_index + order] - knot_vec[knot_index]) != 0.0) {
				sum1 = order
					   / (knot_vec[knot_index + order] - knot_vec[knot_index])
					   * temp_basis[0][order - 1];
			}
			if ( (knot_vec[knot_index + order + 1] - knot_vec[knot_index + 1]) != 0.0) {
				sum2 = order
					   / (knot_vec[knot_index + order + 1] - knot_vec[knot_index + 1])
					   * temp_basis[1][order - 1];
			}
		}
		(*output) = temp_basis[0][order];
		(*d_output) = sum1 - sum2;
	}
	return (*output);
}


double NURBS_surface(double *input_knot_vec_xi, double *input_knot_vec_eta,
                     double *cntl_px, double *cntl_py,
                     int cntl_p_n_xi, int cntl_p_n_eta,
                     double *weight, int order_xi, int order_eta,
                     double xi, double eta,
                     double *output_x, double *output_y,
                     double *output_dxi_x, double *output_deta_x,
                     double *output_dxi_y, double *output_deta_y) {
	int i, j, temp_index;
	double temp1, temp2, temp3;
	double molecule_x, molecule_y;
	double dxi_molecule_x, dxi_molecule_y;
	double deta_molecule_x, deta_molecule_y;
	double denominator, dxi_denominator, deta_denominator;
	double temp_output_xi, temp_output_eta;
	double temp_d_output_xi, temp_d_output_eta;
	molecule_x = 0.0;
	molecule_y = 0.0;
	denominator = 0.0;
	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_denominator = 0.0;
	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_denominator = 0.0;

	int index_min_xi = 0;
	//int index_max_xi = cntl_p_n_xi; //2020_09_12
	int index_max_xi = cntl_p_n_xi - 1; //2020_09_12
	int index_min_eta = 0;
	//int index_max_eta = cntl_p_n_eta; //2020_09_12
	int index_max_eta = cntl_p_n_eta - 1; //2020_09_12

	for (i = 0; i < cntl_p_n_xi; i++) {
		if ( input_knot_vec_xi[i + 1] > xi ) {
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}
	if (index_min_xi < 0) index_min_xi = 0;
	if (index_max_xi > cntl_p_n_xi) index_max_xi = cntl_p_n_xi; //2020_09_12

	for (i = 0; i < cntl_p_n_eta; i++) {
		if ( input_knot_vec_eta[i + 1] > eta ) {
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0) index_min_eta = 0;
	if (index_max_eta > cntl_p_n_eta) index_max_eta = cntl_p_n_eta; //2020_09_12

	for (i = index_min_xi; i <= index_max_xi; i++) {
		BasisFunc(input_knot_vec_xi, i, order_xi, xi,
		          &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j <= index_max_eta; j++) {
			BasisFunc(input_knot_vec_eta, j, order_eta, eta,
			          &temp_output_eta, &temp_d_output_eta);
			temp_index = i + j * cntl_p_n_xi;
			temp1 = temp_output_xi * temp_output_eta
			        * weight[temp_index];
			temp2 = temp_d_output_xi * temp_output_eta
			        * weight[temp_index];
			temp3 = temp_output_xi * temp_d_output_eta
			        * weight[temp_index];
			molecule_x += temp1 * cntl_px[temp_index];
			molecule_y += temp1 * cntl_py[temp_index];
			denominator += temp1;
			dxi_molecule_x += temp2 * cntl_px[temp_index];
			dxi_molecule_y += temp2 * cntl_py[temp_index];
			dxi_denominator += temp2;
			deta_molecule_x += temp3 * cntl_px[temp_index];
			deta_molecule_y += temp3 * cntl_py[temp_index];
			deta_denominator += temp3;
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;

	temp1 = denominator * denominator;
	(*output_dxi_x) = ( dxi_molecule_x * denominator
	                    - molecule_x * dxi_denominator) / temp1;
	(*output_dxi_y) = ( dxi_molecule_y * denominator
	                    - molecule_y * dxi_denominator) / temp1;
	(*output_deta_x) = ( deta_molecule_x * denominator
	                     - molecule_x * deta_denominator) / temp1;
	(*output_deta_y) = ( deta_molecule_y * denominator
	                     - molecule_y * deta_denominator) / temp1;
	return denominator;
}


double rNURBS_surface(double *input_knot_vec_xi, double *input_knot_vec_eta,
                      double *cntl_px, double *cntl_py,
                      int cntl_p_n_xi, int cntl_p_n_eta,
                      double *weight, int order_xi, int order_eta,
                      double xi, double eta,
                      double *output_x, double *output_y,
                      double *output_dxi_x, double *output_deta_x,
                      double *output_dxi_y, double *output_deta_y)
{
	int i, j, temp_index;
	double temp1, temp2, temp3;
	double molecule_x, molecule_y;
	double dxi_molecule_x, dxi_molecule_y;
	double deta_molecule_x, deta_molecule_y;
	double denominator, dxi_denominator, deta_denominator;
	double temp_output_xi, temp_output_eta;
	double temp_d_output_xi, temp_d_output_eta;
	molecule_x = 0.0;
	molecule_y = 0.0;
	denominator = 0.0;
	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_denominator = 0.0;
	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_denominator = 0.0;

	int index_min_xi = 0;
	int index_max_xi = cntl_p_n_xi - 1;
	int index_min_eta = 0;
	int index_max_eta = cntl_p_n_eta - 1;

	for (i = 0; i < cntl_p_n_xi; i++) {
		if ( input_knot_vec_xi[i + 1] >= xi ) {
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}
	if (index_min_xi < 0) index_min_xi = 0;
	if (index_max_xi > cntl_p_n_xi) index_max_xi = cntl_p_n_xi;

	for (i = 0; i < cntl_p_n_eta; i++) {
		if ( input_knot_vec_eta[i + 1] >= eta ) {
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0) index_min_eta = 0;
	if (index_max_eta > cntl_p_n_eta) index_max_eta = cntl_p_n_eta;

	for (i = index_min_xi; i <= index_max_xi; i++) {
		rBasisFunc(input_knot_vec_xi, i, order_xi, xi,
				   &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j <= index_max_eta; j++) {
			rBasisFunc(input_knot_vec_eta, j, order_eta, eta,
					   &temp_output_eta, &temp_d_output_eta);
			temp_index = i + j * cntl_p_n_xi;
			temp1 = temp_output_xi * temp_output_eta
					* weight[temp_index];
			temp2 = temp_d_output_xi * temp_output_eta
					* weight[temp_index];
			temp3 = temp_output_xi * temp_d_output_eta
					* weight[temp_index];
			molecule_x += temp1 * cntl_px[temp_index];
			molecule_y += temp1 * cntl_py[temp_index];
			denominator += temp1;
			dxi_molecule_x += temp2 * cntl_px[temp_index];
			dxi_molecule_y += temp2 * cntl_py[temp_index];
			dxi_denominator += temp2;
			deta_molecule_x += temp3 * cntl_px[temp_index];
			deta_molecule_y += temp3 * cntl_py[temp_index];
			deta_denominator += temp3;
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;

	temp1 = denominator * denominator;
	(*output_dxi_x) = ( dxi_molecule_x * denominator
					 - molecule_x * dxi_denominator) / temp1;
	(*output_dxi_y) = ( dxi_molecule_y * denominator
					 - molecule_y * dxi_denominator) / temp1;
	(*output_deta_x) = ( deta_molecule_x * denominator
					  - molecule_x * deta_denominator) / temp1;
	(*output_deta_y) = ( deta_molecule_y * denominator
					  - molecule_y * deta_denominator) / temp1;
	return denominator;
}


double lNURBS_surface(double *input_knot_vec_xi, double *input_knot_vec_eta,
                      double *cntl_px, double *cntl_py,
                      int cntl_p_n_xi, int cntl_p_n_eta,
                      double *weight, int order_xi, int order_eta,
                      double xi, double eta,
                      double *output_x, double *output_y,
                      double *output_dxi_x, double *output_deta_x,
                      double *output_dxi_y, double *output_deta_y) {
	int i, j, temp_index;
	double temp1, temp2, temp3;
	double molecule_x, molecule_y;
	double dxi_molecule_x, dxi_molecule_y;
	double deta_molecule_x, deta_molecule_y;
	double denominator, dxi_denominator, deta_denominator;
	double temp_output_xi, temp_output_eta;
	double temp_d_output_xi, temp_d_output_eta;
	molecule_x = 0.0;
	molecule_y = 0.0;
	denominator = 0.0;
	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_denominator = 0.0;
	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_denominator = 0.0;

	int index_min_xi = 0;
	int index_max_xi = cntl_p_n_xi - 1;
	int index_min_eta = 0;
	int index_max_eta = cntl_p_n_eta - 1;

	for (i = 0; i < cntl_p_n_xi; i++) {
		if ( input_knot_vec_xi[i + 1] >= xi ) {
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}
	if (index_min_xi < 0) index_min_xi = 0;
	if (index_max_xi > cntl_p_n_xi) index_max_xi = cntl_p_n_xi;

	for (i = 0; i < cntl_p_n_eta; i++) {
		if ( input_knot_vec_eta[i + 1] >= eta ) {
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0) index_min_eta = 0;
	if (index_max_eta > cntl_p_n_eta) index_max_eta = cntl_p_n_eta;

	for (i = index_min_xi; i <= index_max_xi; i++) {
		lBasisFunc(input_knot_vec_xi, i,
				   cntl_p_n_xi, order_xi, xi,
				   &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j <= index_max_eta; j++) {
			lBasisFunc(input_knot_vec_eta, j,
					   cntl_p_n_eta, order_eta, eta,
					   &temp_output_eta, &temp_d_output_eta);
			temp_index = i + j * cntl_p_n_xi;
			temp1 = temp_output_xi * temp_output_eta
					* weight[temp_index];
			temp2 = temp_d_output_xi * temp_output_eta
					* weight[temp_index];
			temp3 = temp_output_xi * temp_d_output_eta
					* weight[temp_index];
			molecule_x += temp1 * cntl_px[temp_index];
			molecule_y += temp1 * cntl_py[temp_index];
			denominator += temp1;
			dxi_molecule_x += temp2 * cntl_px[temp_index];
			dxi_molecule_y += temp2 * cntl_py[temp_index];
			dxi_denominator += temp2;
			deta_molecule_x += temp3 * cntl_px[temp_index];
			deta_molecule_y += temp3 * cntl_py[temp_index];
			deta_denominator += temp3;
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;

	temp1 = denominator * denominator;
	(*output_dxi_x) = ( dxi_molecule_x * denominator
					 - molecule_x * dxi_denominator) / temp1;
	(*output_dxi_y) = ( dxi_molecule_y * denominator
					 - molecule_y * dxi_denominator) / temp1;
	(*output_deta_x) = ( deta_molecule_x * denominator
					  - molecule_x * deta_denominator) / temp1;
	(*output_deta_y) = ( deta_molecule_y * denominator
					  - molecule_y * deta_denominator) / temp1;
	return denominator;
}


double rlNURBS_surface(double *input_knot_vec_xi, double *input_knot_vec_eta,
                       double *cntl_px, double *cntl_py,
                       int cntl_p_n_xi, int cntl_p_n_eta,
                       double *weight, int order_xi, int order_eta,
                       double xi, double eta,
                       double *output_x, double *output_y,
                       double *output_dxi_x, double *output_deta_x,
                       double *output_dxi_y, double *output_deta_y)
{
	int i, j, temp_index;
	double temp1, temp2, temp3;
	double molecule_x, molecule_y;
	double dxi_molecule_x, dxi_molecule_y;
	double deta_molecule_x, deta_molecule_y;
	double denominator, dxi_denominator, deta_denominator;
	double temp_output_xi, temp_output_eta;
	double temp_d_output_xi, temp_d_output_eta;
	molecule_x = 0.0;
	molecule_y = 0.0;
	denominator = 0.0;
	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_denominator = 0.0;
	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_denominator = 0.0;

	int index_min_xi = 0;
	int index_max_xi = cntl_p_n_xi - 1;
	int index_min_eta = 0;
	int index_max_eta = cntl_p_n_eta - 1;

	for (i = 0; i < cntl_p_n_xi; i++) {
		if ( input_knot_vec_xi[i + 1] >= xi ) {
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}
	if (index_min_xi < 0) index_min_xi = 0;
	if (index_max_xi > cntl_p_n_xi) index_max_xi = cntl_p_n_xi;

	for (i = 0; i < cntl_p_n_eta; i++) {
		if ( input_knot_vec_eta[i + 1] >= eta ) {
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0) index_min_eta = 0;
	if (index_max_eta > cntl_p_n_eta) index_max_eta = cntl_p_n_eta;

	for (i = index_min_xi; i <= index_max_xi; i++) {
		rBasisFunc(input_knot_vec_xi, i, order_xi, xi,
				   &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j <= index_max_eta; j++) {
			lBasisFunc(input_knot_vec_eta, j,
					   cntl_p_n_eta, order_eta, eta,
					   &temp_output_eta, &temp_d_output_eta);
			temp_index = i + j * cntl_p_n_xi;
			temp1 = temp_output_xi * temp_output_eta
					* weight[temp_index];
			temp2 = temp_d_output_xi * temp_output_eta
					* weight[temp_index];
			temp3 = temp_output_xi * temp_d_output_eta
					* weight[temp_index];
			molecule_x += temp1 * cntl_px[temp_index];
			molecule_y += temp1 * cntl_py[temp_index];
			denominator += temp1;
			dxi_molecule_x += temp2 * cntl_px[temp_index];
			dxi_molecule_y += temp2 * cntl_py[temp_index];
			dxi_denominator += temp2;
			deta_molecule_x += temp3 * cntl_px[temp_index];
			deta_molecule_y += temp3 * cntl_py[temp_index];
			deta_denominator += temp3;
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;

	temp1 = denominator * denominator;
	(*output_dxi_x) = ( dxi_molecule_x * denominator
					 - molecule_x * dxi_denominator) / temp1;
	(*output_dxi_y) = ( dxi_molecule_y * denominator
					 - molecule_y * dxi_denominator) / temp1;
	(*output_deta_x) = ( deta_molecule_x * denominator
					  - molecule_x * deta_denominator) / temp1;
	(*output_deta_y) = ( deta_molecule_y * denominator
					  - molecule_y * deta_denominator) / temp1;
	return denominator;
}


double lrNURBS_surface(double *input_knot_vec_xi, double *input_knot_vec_eta,
                       double *cntl_px, double *cntl_py,
                       int cntl_p_n_xi, int cntl_p_n_eta,
                       double *weight, int order_xi, int order_eta,
                       double xi, double eta,
                       double *output_x, double *output_y,
                       double *output_dxi_x, double *output_deta_x,
                       double *output_dxi_y, double *output_deta_y)
{
	int i, j, temp_index;
	double temp1, temp2, temp3;
	double molecule_x, molecule_y;
	double dxi_molecule_x, dxi_molecule_y;
	double deta_molecule_x, deta_molecule_y;
	double denominator, dxi_denominator, deta_denominator;
	double temp_output_xi, temp_output_eta;
	double temp_d_output_xi, temp_d_output_eta;
	molecule_x = 0.0;
	molecule_y = 0.0;
	denominator = 0.0;
	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_denominator = 0.0;
	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_denominator = 0.0;

	int index_min_xi = 0;
	int index_max_xi = cntl_p_n_xi - 1;
	int index_min_eta = 0;
	int index_max_eta = cntl_p_n_eta - 1;

	for (i = 0; i < cntl_p_n_xi; i++) {
		if ( input_knot_vec_xi[i + 1] >= xi ) {
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}
	if (index_min_xi < 0) index_min_xi = 0;
	if (index_max_xi > cntl_p_n_xi) index_max_xi = cntl_p_n_xi;

	for (i = 0; i < cntl_p_n_eta; i++) {
		if ( input_knot_vec_eta[i + 1] >= eta ) {
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0) index_min_eta = 0;
	if (index_max_eta > cntl_p_n_eta) index_max_eta = cntl_p_n_eta;

	for (i = index_min_xi; i <= index_max_xi; i++) {
		lBasisFunc(input_knot_vec_xi, i,
				   cntl_p_n_xi, order_xi, xi,
				   &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j <= index_max_eta; j++) {
			rBasisFunc(input_knot_vec_eta, j, order_eta, eta,
					   &temp_output_eta, &temp_d_output_eta);
			temp_index = i + j * cntl_p_n_xi;
			temp1 = temp_output_xi * temp_output_eta
					* weight[temp_index];
			temp2 = temp_d_output_xi * temp_output_eta
					* weight[temp_index];
			temp3 = temp_output_xi * temp_d_output_eta
					* weight[temp_index];
			molecule_x += temp1 * cntl_px[temp_index];
			molecule_y += temp1 * cntl_py[temp_index];
			denominator += temp1;
			dxi_molecule_x += temp2 * cntl_px[temp_index];
			dxi_molecule_y += temp2 * cntl_py[temp_index];
			dxi_denominator += temp2;
			deta_molecule_x += temp3 * cntl_px[temp_index];
			deta_molecule_y += temp3 * cntl_py[temp_index];
			deta_denominator += temp3;
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;

	temp1 = denominator * denominator;
	(*output_dxi_x) = ( dxi_molecule_x * denominator
					 - molecule_x * dxi_denominator) / temp1;
	(*output_dxi_y) = ( dxi_molecule_y * denominator
					 - molecule_y * dxi_denominator) / temp1;
	(*output_deta_x) = ( deta_molecule_x * denominator
					  - molecule_x * deta_denominator) / temp1;
	(*output_deta_y) = ( deta_molecule_y * denominator
					  - molecule_y * deta_denominator) / temp1;
	return denominator;
}


//算出したローカルパッチ各要素の頂点の物理座標のグローバルパッチでの(xi,eta)算出
int Calc_xi_eta(double px, double py,
                  double *input_knot_vec_xi, double *input_knot_vec_eta,
                  int cntl_p_n_xi, int cntl_p_n_eta, int order_xi, int order_eta,
                  double *output_xi, double *output_eta){
	double temp_xi, temp_eta;
	double temp_x, temp_y;
	double temp_matrix[2][2];
	double temp_dxi, temp_deta;
    //double temp_tol_x = DBL_MAX;
	//double temp_tol_y = DBL_MAX;
    double temp_tol_x, temp_tol_y;

	(*output_xi) = 0;
	(*output_eta) = 0;

	int i;
	// int repeat = 1000;
	// double tol = 10e-8;
	// int repeat = 10000;
	int repeat = 100;
	double tol = 10e-14;

	//初期値の設定
	temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;
	//printf("r_temp_xi_eta % 1.8e % 1.8e\n", temp_xi, temp_eta);

	for (i = 0; i < repeat; i++)
	{
		rNURBS_surface(Position_Knots[0][0], Position_Knots[0][1],
		               Control_Coord[0], Control_Coord[1],
	  	  			   No_Control_point[0][0], No_Control_point[0][1],
		               Control_Weight, Order[0][0], Order[0][1],
		               temp_xi, temp_eta,
		               &temp_x, &temp_y,
		               &temp_matrix[0][0], &temp_matrix[0][1],
		               &temp_matrix[1][0], &temp_matrix[1][1]);

        temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;

        //printf("xi_0:  % 1.8e\n", temp_xi);
		//printf("eta_0: % 1.8e\n", temp_eta);
        //printf("px: % 1.8e\n",px);
        //printf("temp_x: % 1.8e\n",temp_x);
        //printf("py: % 1.8e\n",py);
        //printf("temp_y: % 1.8e\n",temp_y);
        //printf("temp_tol_x:  % 1.8e\n", temp_tol_x);
		//printf("temp_tol_y:  % 1.8e\n", temp_tol_y);

		//収束した場合////////////////////////////////////////////////////////////////
        //if (temp_tol_x < tol && temp_tol_y < tol)
        if (temp_tol_x + temp_tol_y < tol)
		{
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;

            //printf("r_xi:  % 1.8e\n", temp_xi);
			//printf("r_eta: % 1.8e\n", temp_eta);

            return i;
        }

		InverseMatrix_2D(temp_matrix);

		temp_dxi = temp_matrix[0][0] * (px - temp_x)
		           + temp_matrix[0][1] * (py - temp_y);
		temp_deta = temp_matrix[1][0] * (px - temp_x)
		            + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_deta;
		if (temp_xi < input_knot_vec_xi[0])
			temp_xi = input_knot_vec_xi[0];
		if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
			temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
		if (temp_eta < input_knot_vec_eta[0])
			temp_eta = input_knot_vec_eta[0];
		if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
			temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];

        //printf("r_xi:  % 1.8e\n", temp_xi);
		//printf("r_eta: % 1.8e\n", temp_eta);
        //printf("i=%d\n",i);

		//double temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
		//printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
    }

	//初期値の設定
	temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;
	//printf("l_temp_xi_eta % 1.8e % 1.8e\n", temp_xi, temp_eta);

	for (i = 0; i < repeat; i++)
	{
		lNURBS_surface(Position_Knots[0][0], Position_Knots[0][1],
		               Control_Coord[0], Control_Coord[1],
	  	  			   No_Control_point[0][0], No_Control_point[0][1],
		               Control_Weight, Order[0][0], Order[0][1],
		               temp_xi, temp_eta,
		               &temp_x, &temp_y,
		               &temp_matrix[0][0], &temp_matrix[0][1],
		               &temp_matrix[1][0], &temp_matrix[1][1]);

        temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;

		//収束した場合////////////////////////////////////////////////////////////////
        //if (temp_tol_x < tol && temp_tol_y < tol)
        if (temp_tol_x + temp_tol_y < tol)
        {
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;

            //printf("l_xi:  % 1.8e\n", temp_xi);
			//printf("l_eta: % 1.8e\n", temp_eta);
            //printf("i=%d\n",i);

            return i;
        }

		InverseMatrix_2D(temp_matrix);

		temp_dxi = temp_matrix[0][0] * (px - temp_x)
		           + temp_matrix[0][1] * (py - temp_y);
		temp_deta = temp_matrix[1][0] * (px - temp_x)
		            + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_deta;
		if (temp_xi < input_knot_vec_xi[0])
			temp_xi = input_knot_vec_xi[0];
		if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
			temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
		if (temp_eta < input_knot_vec_eta[0])
			temp_eta = input_knot_vec_eta[0];
		if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
			temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];

		//double temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
		//printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
    }

	//初期値の設定
	temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;
	//printf("rl_temp_xi_eta % 1.8e % 1.8e\n", temp_xi, temp_eta);

	for (i = 0; i < repeat; i++)
	{
		rlNURBS_surface(Position_Knots[0][0], Position_Knots[0][1],
		                Control_Coord[0], Control_Coord[1],
	  	  			    No_Control_point[0][0], No_Control_point[0][1],
		                Control_Weight, Order[0][0], Order[0][1],
		                temp_xi, temp_eta,
		                &temp_x, &temp_y,
		                &temp_matrix[0][0], &temp_matrix[0][1],
		                &temp_matrix[1][0], &temp_matrix[1][1]);

        temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;

		//収束した場合////////////////////////////////////////////////////////////////
        //if (temp_tol_x < tol && temp_tol_y < tol)
        if (temp_tol_x + temp_tol_y < tol)
        {
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;

            //printf("rl_xi:  % 1.8e\n", temp_xi);
			//printf("rl_eta: % 1.8e\n", temp_eta);
            //printf("i=%d\n",i);

            return i;
        }

		InverseMatrix_2D(temp_matrix);

		temp_dxi = temp_matrix[0][0] * (px - temp_x)
		           + temp_matrix[0][1] * (py - temp_y);
		temp_deta = temp_matrix[1][0] * (px - temp_x)
		            + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_deta;
		if (temp_xi < input_knot_vec_xi[0])
			temp_xi = input_knot_vec_xi[0];
		if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
			temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
		if (temp_eta < input_knot_vec_eta[0])
			temp_eta = input_knot_vec_eta[0];
		if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
			temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];

		//double temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
		//printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
    }

	//初期値の設定
	temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;
	//printf("lr_temp_xi_eta % 1.8e % 1.8e\n", temp_xi, temp_eta);

	for (i = 0; i < repeat; i++)
	{
		lrNURBS_surface(Position_Knots[0][0], Position_Knots[0][1],
		                Control_Coord[0], Control_Coord[1],
	  	  			    No_Control_point[0][0], No_Control_point[0][1],
		                Control_Weight, Order[0][0], Order[0][1],
		                temp_xi, temp_eta,
		                &temp_x, &temp_y,
		                &temp_matrix[0][0], &temp_matrix[0][1],
		                &temp_matrix[1][0], &temp_matrix[1][1]);

        temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;

		//収束した場合////////////////////////////////////////////////////////////////
        //if (temp_tol_x < tol && temp_tol_y < tol)
        if (temp_tol_x + temp_tol_y < tol)
        {
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;

            //printf("lr_xi:  % 1.8e\n", temp_xi);
			//printf("lr_eta: % 1.8e\n", temp_eta);
            //printf("i=%d\n",i);

            return i;
        }

		InverseMatrix_2D(temp_matrix);

		temp_dxi = temp_matrix[0][0] * (px - temp_x)
		           + temp_matrix[0][1] * (py - temp_y);
		temp_deta = temp_matrix[1][0] * (px - temp_x)
		            + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_deta;
		if (temp_xi < input_knot_vec_xi[0])
			temp_xi = input_knot_vec_xi[0];
		if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
			temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
		if (temp_eta < input_knot_vec_eta[0])
			temp_eta = input_knot_vec_eta[0];
		if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
			temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];

		//double temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
		//printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
    }

	//printf("% 1.8e % 1.8e\n", temp_x, temp_y);
	return 0;
}


////////////////////////////////////////////////////////////////////////
/////////////////要素剛性マトリックス
////////////////////////////////////////////////////////////////////////
int Jacobian(int El_No, double a[DIMENSION][DIMENSION], double Local_coord[DIMENSION], double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION])
{
	int i, j, k;
	//printf("El_No_jacobi:%d\n",El_No);
	for (i = 0; i < DIMENSION; i++)
	{
		for (j = 0; j < DIMENSION; j++)
		{
			a[i][j] = 0.0;
			for (k = 0; k < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; k++)
			{
				//printf("Local_coord[%d]:%le\n",j,Local_coord[j] );
				a[i][j] += dShape_func(k, j, Local_coord, El_No) * X[k][i];
				//printf(" X[%d][%d]=%lf\t",k,i, X[k][i] );
				//printf("k=%d a[%d][%d]:%le\n",k,i,j,a[i][j]);
			}
			//printf("<<<最終a[%d][%d]:%le>>>\n",i,j,a[i][j]);
		}
	}
	/*for (i = 0; i < DIMENSION; i++) {
		for (j = 0; j < DIMENSION; j++) {
			printf("a[%d][%d]:%le\n",i,j,a[i][j]);
		}
	}*/
	return 0;
}


//Bマトリックスを求める関数
int Make_B_Matrix(int El_No,
                  double B[D_MATRIX_SIZE][MAX_KIEL_SIZE],
                  double Local_coord[DIMENSION],
				  double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION],
                  double *J)
{
	double a[DIMENSION][DIMENSION], b[DIMENSION][MAX_NO_CCpoint_ON_ELEMENT];

	int i, j, k;

	Jacobian(El_No, a, Local_coord, X);

	*J = InverseMatrix_2D(a);
	//printf("B_Matri_J:%le\n",*J);
	if (*J <= 0)
		return -999;

	for (i = 0; i < DIMENSION; i++)
	{
		for (j = 0; j < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; j++)
		{
			b[i][j] = 0.0;
			for (k = 0; k < DIMENSION; k++)
			{
				b[i][j] += a[k][i] * dShape_func(j, k, Local_coord, El_No);
			}
		}
	}

	for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++)
	{
		B[0][2 * i] = b[0][i];
		B[0][2 * i + 1] = 0.0;
		B[1][2 * i] = 0.0;
		B[1][2 * i + 1] = b[1][i];
		B[2][2 * i] = b[1][i];
		B[2][2 * i + 1] = b[0][i];
	}
	/*for (i = 0; i < D_MATRIX_SIZE; i++)for (j = 0; j < KIEL_SIZE; j++) {
		printf("B[%d][%d]_B_mat:%le\n",i,j,B[i][j]);
	}*/
	return 0;
}


//変位勾配を求めるためのマトリクスb_gradを求める関数
int Make_b_grad_Matrix(int El_No, double b_grad[DIMENSION * DIMENSION][2 * MAX_NO_CCpoint_ON_ELEMENT], double Local_coord[DIMENSION],
					   double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], double *J)
{
		double a[DIMENSION][DIMENSION], b[DIMENSION][MAX_NO_CCpoint_ON_ELEMENT];
		int i, j, k;

		// printf("El_No : %d\n", El_No);
		// printf("No_Control_point_ON_ELEMENT[Element_patch[El_No]] : %d\n", No_Control_point_ON_ELEMENT[Element_patch[El_No]]);

		for(i = 0; i < DIMENSION * DIMENSION; i++){
			for(j = 0; j < 2 * No_Control_point_ON_ELEMENT[Element_patch[El_No]]; j++){
				b_grad[i][j] = 0.0;
			}
		}

		Jacobian( El_No, a, Local_coord, X);

		*J = InverseMatrix_2D( a );
		
		if( *J <= 0 )return -999;

		for( i = 0; i < DIMENSION; i++ ){
			for( j = 0; j < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; j++ ){
				b[i][j] = 0.0;
				for( k = 0; k < DIMENSION; k++ ){
					b[i][j] += a[k][i] * dShape_func(j, k, Local_coord, El_No);
				}
			}
		}
		
		/*３次元にする場合は修正が必要*/
		for( j = 0; j < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; j++ ){
			b_grad[0][2 * j] = b[0][j];    b_grad[0][2 * j + 1] = 0.0;         //   𝜕N1/𝜕x1    0.0      𝜕N2/𝜕x1    0.0      𝜕N3/𝜕x1      0.0    𝜕N4/𝜕x1    0.0   ...
			b_grad[1][2 * j] = b[1][j];    b_grad[1][2 * j + 1] = 0.0;         //   𝜕N1/𝜕x2    0.0      𝜕N2/𝜕x2    0.0      𝜕N3/𝜕x2      0.0    𝜕N4/𝜕x2    0.0   ...
			b_grad[2][2 * j] = 0.0;        b_grad[2][2 * j + 1] = b[0][j];     //     0.0      𝜕N1/𝜕x1    0.0      𝜕N2/𝜕x1    0.0      𝜕N3/𝜕x1    0.0    𝜕N4/𝜕x1 ...
			b_grad[3][2 * j] = 0.0;        b_grad[3][2 * j + 1] = b[1][j];     //     0.0      𝜕N1/𝜕x2    0.0      𝜕N2/𝜕x2    0.0      𝜕N3/𝜕x2    0.0    𝜕N4/𝜕x2 ...
		}

		/*for (i = 0; i < DIMENSION * DIMENSION; i++){
		    for (j = 0; j < 2 * No_Control_point_ON_ELEMENT[Element_patch[El_No]]; j++) {
			    printf("b_grad[%d][%d] =  %le\t",i, j, b_grad[i][j]);
			}
			printf("\n");
		}*/
		return 0;
	}


//応力歪マトリックス
int Make_D_Matrix_2D(double D[D_MATRIX_SIZE][D_MATRIX_SIZE], double E, double nu, int DM)
{
	int i, j;

	if (DM == 0) //平面応力状態
	{
		//printf("E:%le nu:%le\n",E,nu);
		double Eone = E / (1.0 - nu * nu);
		double D1[D_MATRIX_SIZE][D_MATRIX_SIZE] = {{Eone, nu * Eone, 0}, {nu * Eone, Eone, 0}, {0, 0, (1 - nu) / 2 * Eone}};

		for (i = 0; i < D_MATRIX_SIZE; i++)
			for (j = 0; j < D_MATRIX_SIZE; j++)
				D[i][j] = D1[i][j];
	}

	else if (DM == 1) //平面ひずみ状態(2Dの場合はこっち)
	{
		//printf("E:%le nu:%le\n",E,nu);
		double Eone = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
		double D1[D_MATRIX_SIZE][D_MATRIX_SIZE] = {{Eone, nu / (1.0 - nu) * Eone, 0}, {nu / (1.0 - nu) * Eone, Eone, 0}, {0, 0, (1 - 2 * nu) / 2 / (1.0 - nu) * Eone}};

		for (i = 0; i < D_MATRIX_SIZE; i++)
			for (j = 0; j < D_MATRIX_SIZE; j++)
				D[i][j] = D1[i][j];
	}

	else
		return ERROR;

	return 0;
}


//ガウスの数値積分法の中身
int BDBJ(double B[D_MATRIX_SIZE][MAX_KIEL_SIZE], double D[D_MATRIX_SIZE][D_MATRIX_SIZE], double J, double K_EL[MAX_KIEL_SIZE][MAX_KIEL_SIZE])
{
	int i, j, k;
	double BD[MAX_KIEL_SIZE][D_MATRIX_SIZE];

	//[B]T[D][B]Jの計算
	for (i = 0; i < KIEL_SIZE; i++)
	{
		for (j = 0; j < D_MATRIX_SIZE; j++)
		{
			BD[i][j] = 0;
			for (k = 0; k < D_MATRIX_SIZE; k++)
			{
				//printf("B[%d][%d]=%le D[%d][%d]=%le\n", k,i,B[k][i],k,j,D[k][j] );
				//printf("B[%d][%d]=%le\n", k,i,B[k][i]);
				BD[i][j] += B[k][i] * D[k][j];
				//printf("BD[%d][%d]=%e\n",i,j,BD[i][j] );
			}
		}
	}
	//for( j = 0; j < D_MATRIX_SIZE; j++ )for( i = 0; i < KIEL_SIZE; i++ )printf("B[%d][%d]=%le\n",j,i,B[j][i] );
	for (i = 0; i < KIEL_SIZE; i++)
	{
		for (j = 0; j < KIEL_SIZE; j++)
		{
			K_EL[i][j] = 0;
			for (k = 0; k < D_MATRIX_SIZE; k++)
			{
				K_EL[i][j] += BD[i][k] * B[k][j];
			}
			K_EL[i][j] *= J;
		}
	}
	return 0;
}


//結合ガウスの数値積分法の中身
int coupled_BDBJ(double B[D_MATRIX_SIZE][MAX_KIEL_SIZE],
				 double D[D_MATRIX_SIZE][D_MATRIX_SIZE],
				 double BG[D_MATRIX_SIZE][MAX_KIEL_SIZE],
				 double J, double K_EL[MAX_KIEL_SIZE][MAX_KIEL_SIZE])
{
	int i, j, k;
	double BD[MAX_KIEL_SIZE][D_MATRIX_SIZE];

	//[B]T[D][B]Jの計算
	for (i = 0; i < KIEL_SIZE; i++)
	{
		for (j = 0; j < D_MATRIX_SIZE; j++)
		{
			BD[i][j] = 0;
			for (k = 0; k < D_MATRIX_SIZE; k++)
			{
				//printf("B[%d][%d]=%le D[%d][%d]=%le\n", k,i,B[k][i],k,j,D[k][j] );
				//printf("B[%d][%d]=%le\n", k,i,B[k][i]);
				BD[i][j] += BG[k][i] * D[k][j];
				//printf("BD[%d][%d]=%e\n",i,j,BD[i][j] );
			}
		}
	}
	//for( j = 0; j < D_MATRIX_SIZE; j++ )for( i = 0; i < KIEL_SIZE; i++ )printf("B[%d][%d]=%le\n",j,i,B[j][i] );
	for (i = 0; i < KIEL_SIZE; i++)
	{
		for (j = 0; j < KIEL_SIZE; j++)
		{
			K_EL[i][j] = 0;
			for (k = 0; k < D_MATRIX_SIZE; k++)
			{
				K_EL[i][j] += BD[i][k] * B[k][j];
			}
			K_EL[i][j] *= J;
		}
	}
	return 0;
}


//要素合成マトリックス
int Make_K_EL(int El_No, double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], double K_EL[MAX_KIEL_SIZE][MAX_KIEL_SIZE], double E, double nu, int DM)
{
	int i, j, k, l;

	double K1[MAX_KIEL_SIZE][MAX_KIEL_SIZE], B[D_MATRIX_SIZE][MAX_KIEL_SIZE];
	double D[D_MATRIX_SIZE][D_MATRIX_SIZE];
	double J = 0.0;
	double J_test = 0.0;

	for (i = 0; i < KIEL_SIZE; i++)
	{
		for (j = 0; j < KIEL_SIZE; j++)
		{
			K_EL[i][j] = 0.0;
		}
	}

	Make_D_Matrix_2D(D, E, nu, DM);

	// for (i = 0; i < POW_Ng; i++)
	for (i = 0; i < GP_2D; i++)
	{

		//printf("i=%d\n",i );
		Make_B_Matrix(El_No, B, Gxi[i], X, &J);

		BDBJ(B, D, J, K1);
		J_test += J;
		for (k = 0; k < KIEL_SIZE; k++)
		{
			for (l = 0; l < KIEL_SIZE; l++)
			{
				K_EL[k][l] += w[i] * K1[k][l];
			}
		} //printf("w[%d]=%f\n",i,w[i]);
	}
	//printf("El_No:%d J_test=%e\n", El_No, J_test);
	//printf("G=%f\n",G );
	/*for ( k = 0; k < KIEL_SIZE; k++) {
		for ( l = 0; l < KIEL_SIZE; l++) {
			printf("K_EL[%d][%d]:%le\n",k,l,K_EL[k][l]);
		}
	}*/

	return 0;
}


//結合要素剛性マトリックス
int Make_coupled_K_EL(int El_No_loc, int El_No_glo,
					  double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION],
					  double XG[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION],
					  double K_EL[MAX_KIEL_SIZE][MAX_KIEL_SIZE],
					  double E, double nu, int DM)
{
	int i, j, jj, k, l;
	int BDBJ_flag;

	double K1[MAX_KIEL_SIZE][MAX_KIEL_SIZE];
	double B[D_MATRIX_SIZE][MAX_KIEL_SIZE], BG[D_MATRIX_SIZE][MAX_KIEL_SIZE];
	double D[D_MATRIX_SIZE][D_MATRIX_SIZE];
	double J = 0.0;
	double J_test = 0.0;
	double G_Gxi[GP_2D][DIMENSION];	//グローバルパッチ上での親要素内座標xi_bar,eta_bar

	Total_BDBJ_flag = 0;

	//printf("El_No;%d\tcoupled_El_No;%d\n",El_No_loc, El_No_glo);
	for (i = 0; i < KIEL_SIZE; i++)
	{
		for (j = 0; j < KIEL_SIZE; j++)
		{
			K_EL[i][j] = 0.0;
		}
	}

	Make_D_Matrix_2D(D, E, nu, DM);

	for (i = 0; i < GP_2D; i++)	//ガウス点のループ(local)
	{
		// printf("gauss point number:%d\n", i);

		////ローカルガウス点がグローバル要素に含まれているかの判定
		//ローカル要素ガウス点の物理座標算出
		double data_result_shape[2] = {0.0};
		double output_xi, output_eta;
		int patch_n = 0;

		for (j = 0; j < No_Control_point_ON_ELEMENT[Element_patch[El_No_loc]]; j++)
		{
			double R_shape_func = Shape_func(j, Gxi[i], El_No_loc);

			for (jj = 0; jj < DIMENSION; jj++)
			{
				data_result_shape[jj] += R_shape_func * X[j][jj];
			}
		}

		//ローカル要素ガウス点のグローバルパッチ上のパラメータ空間座標算出
		for (j = 0; j < Total_Patch_on_mesh[0]; j++)	//グローバルメッシュ[0]上
		{
			Calc_xi_eta(data_result_shape[0], data_result_shape[1],
	            	    Position_Knots[j][0], Position_Knots[j][1],
	                	No_Control_point[j][0], No_Control_point[j][1], Order[j][0], Order[j][1],
	                	&output_xi, &output_eta);
			// printf("  x: % 1.8e\n", data_result_shape[0]);
			// printf("  y: % 1.8e\n", data_result_shape[1]);
			// printf(" xi: % 1.8e\n", output_xi);
			// printf("eta: % 1.8e\n", output_eta);
			// printf("patch_n: %d\n", j);
			patch_n = j;
		}
		//要素内外判定

		if (output_xi >= Position_Knots[patch_n][0][Order[patch_n][0]+ENC[patch_n][El_No_glo][0]]
			&&
			output_xi < Position_Knots[patch_n][0][Order[patch_n][0]+ENC[patch_n][El_No_glo][0]+1]
			&&
			output_eta >= Position_Knots[patch_n][1][Order[patch_n][1]+ENC[patch_n][El_No_glo][1]]
			&&
			output_eta < Position_Knots[patch_n][1][Order[patch_n][1]+ENC[patch_n][El_No_glo][1]+1])	//要素内であるとき
		{
			BDBJ_flag = 1;
			// printf("BDBJ_flag\n");

			//親要素座標の算出
			G_Gxi[i][0] = -1.0
					  + 2.0 * (output_xi - Position_Knots[patch_n][0][Order[patch_n][0]+ENC[patch_n][El_No_glo][0]]) /
					  (Position_Knots[patch_n][0][Order[patch_n][0]+ENC[patch_n][El_No_glo][0]+1] - Position_Knots[patch_n][0][Order[patch_n][0]+ENC[patch_n][El_No_glo][0]]);
			G_Gxi[i][1] = -1.0
					  + 2.0 * (output_eta - Position_Knots[patch_n][1][Order[patch_n][1]+ENC[patch_n][El_No_glo][1]]) /
					  (Position_Knots[patch_n][1][Order[patch_n][1]+ENC[patch_n][El_No_glo][1]+1] - Position_Knots[patch_n][1][Order[patch_n][1]+ENC[patch_n][El_No_glo][1]]);
			//printf("G_Gxi[][]=\n");
		}
		else	//要素外であるとき
		{
			BDBJ_flag = 0;
		}

		//printf("i=%d\n",i );

		////結合要素剛性マトリックス計算
		//要素内であるとき、次を計算
		if (BDBJ_flag)
		{
			Total_BDBJ_flag++;
			Same_BDBJ_flag[i]++;
			//printf("BDBJ_flag\ti=%d\n",i );
			//重なるグローバル要素のBマトリックス
			Make_B_Matrix(El_No_glo, BG, G_Gxi[i], XG, &J);
			//ローカル要素のBマトリックス
			Make_B_Matrix(El_No_loc, B, Gxi[i], X, &J);
			//BGTDBLの計算
			coupled_BDBJ(B, D, BG, J, K1);
			J_test += J;
			for (k = 0; k < KIEL_SIZE; k++)
			{
				for (l = 0; l < KIEL_SIZE; l++)
				{
					K_EL[k][l] += w[i] * K1[k][l];
				}
			} //printf("w[%d]=%f\n",i,w[i]);
		}
	}
	//printf("El_No:%d J_test=%e\n", El_No_loc, J_test);
	//printf("G=%f\n",G );
	/*for ( k = 0; k < KIEL_SIZE; k++) {
		for ( l = 0; l < KIEL_SIZE; l++) {
			printf("K_EL[%d][%d]:%le\n",k,l,K_EL[k][l]);
		}
	}*/

	// if (i == GP_2D - 1)
	// {
	// 	printf("-------------------Total_BDBJ_flag=%d-------------------\n", Total_BDBJ_flag);
	// }

	return 0;
}


///////////////////////////////////////////////////
//////////////歪と応力, ひずみエネルギ密度, 変位勾配
///////////////////////////////////////////////////
void Make_Strain(int Total_Element)
{
	static double U[MAX_KIEL_SIZE];
	double B[D_MATRIX_SIZE][MAX_KIEL_SIZE], X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], J;

	int N, e, i, j;
	//printf("Strain\n");

	Make_gauss_array(0);

	for (e = 0; e < Total_Element; e++)
	{
		//printf("\nElementNo.:%d\n",e);
		for (N = 0; N < GP_2D; N++)
			for (i = 0; i < N_STRAIN; i++)
				Strain[e][N][i] = 0.0;
		//Bマトリックスと各要素の変位を取得
		//printf("El_No:%d,No_Control_point_ON_ELEMENT[%d]:%d\n", El_No,Element_patch[El_No],No_Control_point_ON_ELEMENT[Element_patch[El_No]]);
		for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[e]]; i++)
		{
			for (j = 0; j < DIMENSION; j++)
			{
				U[i * DIMENSION + j] = Displacement[Controlpoint_of_Element[e][i] * DIMENSION + j];
				X[i][j] = Node_Coordinate[Controlpoint_of_Element[e][i]][j];
			}
		}
		//歪
		for (N = 0; N < GP_2D; N++)
		{
			//printf("N:%d\n",N);
			Make_B_Matrix(e, B, Gxi[N], X, &J);
			for (i = 0; i < D_MATRIX_SIZE; i++)
				for (j = 0; j < KIEL_SIZE; j++)
				{
					Strain[e][N][i] += B[i][j] * U[j];
					//printf("B[%d][%d]_in_strain:%le * ",i,j,B[i][j]);
					//if(e==1){
					//printf("U[%d]=%le = %le\n",j,U[j],B[i][j]*U[j]);
					//}
				}
		}
	}
}


//応力
void Make_Stress_2D(double E, double nu, int Total_Element, int DM)
{

	static double D[D_MATRIX_SIZE][D_MATRIX_SIZE];
	int e, i, j, k;
	Make_gauss_array(0);
	Make_D_Matrix_2D(D, E, nu, DM);

	for (e = 0; e < Total_Element; e++)
	{
		for (k = 0; k < GP_2D; k++)
			for (i = 0; i < N_STRESS; i++)
				Stress[e][k][i] = 0.0;
		for (k = 0; k < GP_2D; k++)
			for (i = 0; i < D_MATRIX_SIZE; i++)
				for (j = 0; j < D_MATRIX_SIZE; j++)
					Stress[e][k][i] += D[i][j] * Strain[e][k][j];
	}
}


void Make_Stress_2D_glo(double E, double nu, int Total_Element, int DM)
{
	Make_gauss_array(0);

	static double D[D_MATRIX_SIZE][D_MATRIX_SIZE];
	int e, i, j, k;
	Make_D_Matrix_2D(D, E, nu, DM);

	for (e = 0; e < Total_Element; e++)
	{
		for (k = 0; k < GP_2D; k++)
			for (i = 0; i < N_STRESS; i++)
				Stress_glo[e][k][i] = 0.0;
		for (k = 0; k < GP_2D; k++)
			for (i = 0; i < D_MATRIX_SIZE; i++)
				for (j = 0; j < D_MATRIX_SIZE; j++)
					Stress_glo[e][k][i] += D[i][j] * Strain_glo[e][k][j];
	}
}


void Make_Displacement_grad(int El_No)
{
	static double U[MAX_KIEL_SIZE];
	double b_grad[DIMENSION * DIMENSION][2 * MAX_NO_CCpoint_ON_ELEMENT]; 
	double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], J;

	Make_gauss_array(0);

	int N, e, re, i, j;

	//printf("変位勾配\n");
	for( re = 0; re < real_Total_Element_to_mesh[Total_mesh]; re++ ){
		e = real_element[re];
		//printf("\nElementNo.:%d\n",e);
		for( N = 0; N < GP_2D; N++)
			for( i = 0; i <  DIMENSION * DIMENSION; i ++ )
				Disp_grad[e][N][i] = 0.0;
		//Bマトリックスと各要素の変位を取得
		//printf("El_No:%d,No_Control_point_ON_ELEMENT[%d]:%d\n", El_No,Element_patch[El_No],No_Control_point_ON_ELEMENT[Element_patch[El_No]]);
		for( i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[e]]; i++ ){
			for( j = 0; j < DIMENSION; j++ ){
				U[ i*DIMENSION +j ] = Displacement[ Controlpoint_of_Element[e][i]*DIMENSION + j ];
				X[i][j] = Node_Coordinate[ Controlpoint_of_Element[e][i] ][j];
			}
		}
		//変位勾配
		for( N = 0; N < GP_2D; N++ ){
			//printf("N:%d\n",N);
			Make_b_grad_Matrix( e, b_grad, Gxi[N], X, &J);
			for( i = 0; i < DIMENSION * DIMENSION; i++ )
				for( j = 0; j < 2 * No_Control_point_ON_ELEMENT[Element_patch[El_No]]; j++ ){
					Disp_grad[e][N][i] += b_grad[i][j] * U[j] ;
					//printf("b_grad[%d][%d] = %lf  また　U = %lf",i,j,b_grad[i][j],U[j]);
				}
		}
	}
}


void Make_StrainEnergyDensity_2D()
{
	int re,e, k;
	//Make_D_Matrix_2D( D, E, nu ,DM);
	for( re = 0; re < real_Total_Element_to_mesh[Total_mesh]; re++ ){
		e=real_element[re];	
		for( k = 0; k < POW_Ng; k++){
			StrainEnergyDensity[e][k] = 0.0;
		}
	}

	for( re = 0; re < real_Total_Element_to_mesh[Total_mesh]; re++ ){
		e=real_element[re];
		for( k = 0; k < POW_Ng; k++){
			StrainEnergyDensity[e][k] = ( Stress[e][k][0] * Strain[e][k][0] + Stress[e][k][1] * Strain[e][k][1] + Stress[e][k][2] * Strain[e][k][2]) / 2.0;
		}
	}
}


	// void Make_Displacement_grad_overlay(double E, double nu, int Total_Element , int El_No, int Total_Control_Point){
	// 	static double U[MAX_KIEL_SIZE];
	// 	double b_grad[DIMENSION * DIMENSION][2 * No_Control_point_ON_ELEMENT[Element_patch[El_No]]]; 
	// 	double X[No_Control_point_ON_ELEMENT[Element_patch[El_No]]][DIMENSION],J;
	// 	double G = pow(0.6,0.5);
	// 	double Gxi[POW_Ng][DIMENSION] = { {-G,-G},{0.0,-G},{G,-G},{-G,0.0},{0.0,0.0},{G,0.0},{-G,G},{0.0,G},{G,G} };


	// 	int N, e, re, i, j;
	// 	//printf("変位勾配\n");
    //     	for(re = 0; re < real_Total_Element_to_mesh[Total_mesh]; re++){
    //         	e = real_element[re];
	// 		//printf("\nElementNo.:%d\n",e);
	// 		for(N = 0; N < POW_Ng; N++)
	// 			for(i = 0; i <  DIMENSION * DIMENSION; i++)
	// 				Disp_grad_overlay[e][N][i] = 0.0;
	// 		//Bマトリックスと各要素の変位を取得
	// 		//printf("El_No:%d,No_Control_point_ON_ELEMENT[%d]:%d\n", El_No,Element_patch[El_No],No_Control_point_ON_ELEMENT[Element_patch[El_No]]);
	// 		for(i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[e]]; i++){
	// 			for(j = 0; j < DIMENSION; j++){
	// 				U[i * DIMENSION + j] = Disp_overlay[Controlpoint_of_Element[e][i] * DIMENSION + j];
	// 				X[i][j] = Node_Coordinate[Controlpoint_of_Element[e][i]][j];
	// 			}
	// 		}
	// 		//変位勾配
	// 		for(N = 0; N < POW_Ng; N++){
	// 			//printf("N:%d\n",N);
	// 			Make_b_grad_Matrix(e, b_grad, Gxi[N], X ,&J , Total_Control_Point);
	// 			for(i = 0; i < DIMENSION * DIMENSION; i++)
	// 				for(j = 0; j < 2 * No_Control_point_ON_ELEMENT[Element_patch[El_No]]; j++){
	// 				Disp_grad_overlay[e][N][i] += b_grad[i][j] * U[j];
	// 				//printf("b_grad[%d][%d] = %lf  また　U = %lf",i,j,b_grad[i][j],U[j]);
	// 			}
	// 		}
	// 	}
	// }


void Make_StrainEnergyDensity_2D_overlay()
{
	int re, e, k;
	for(re = 0; re < real_Total_Element_to_mesh[Total_mesh]; re++){
		e = real_element[re];	
		for(k = 0; k < POW_Ng; k++){
			StrainEnergyDensity_overlay[e][k] = 0.0;
		}
	}

	for(re = 0; re < real_Total_Element_to_mesh[Total_mesh]; re++){
		e = real_element[re];
		for(k = 0; k < POW_Ng; k++){
			StrainEnergyDensity_overlay[e][k] = (Stress_overlay[e][k][0] * Strain_overlay[e][k][0] + Stress_overlay[e][k][1] * Strain_overlay[e][k][1] + Stress_overlay[e][k][2] * Strain_overlay[e][k][2]) / 2.0;
		}
	}
}


void Make_ReactionForce(int Total_Control_Point)
{
	int e, i, j, k, l, re;
	double B[D_MATRIX_SIZE][MAX_KIEL_SIZE], X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], J;

	Make_gauss_array(0);

	for (i = 0; i < Total_Control_Point * DIMENSION; i++)
		ReactionForce[i] = 0.0;
	//printf("ReactionForce\n");
	for (re = 0; re < real_Total_Element_to_mesh[Total_mesh]; re++)
	{
		e = real_element[re];
		//printf("%d,No_Control_point_ON_ELEMENT[%d]:%d\n",e,Element_patch[e], No_Control_point_ON_ELEMENT[Element_patch[e]]);
		//Bマトリックスを取得
		for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[e]]; i++)
		{
			for (j = 0; j < DIMENSION; j++)
			{
				X[i][j] = Node_Coordinate[Controlpoint_of_Element[e][i]][j];
				//printf("Node_Coordinate[%d][%d]:%le\n",Controlpoint_of_Element[e][i],j, Node_Coordinate[ Controlpoint_of_Element[e][i] ][j]);
				//printf("X[%d][%d]:%le\n",Controlpoint_of_Element[e][i],j,X[i][j] );
			}
		}
		for (k = 0; k < GP_2D; k++)
		{
			Make_B_Matrix(e, B, Gxi[k], X, &J);
			for (j = 0; j < D_MATRIX_SIZE; j++)
				for (l = 0; l < DIMENSION; l++)
					for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[e]]; i++)
						ReactionForce[Controlpoint_of_Element[e][i] * DIMENSION + l] += B[j][i * DIMENSION + l] * Stress[e][k][j] * w[k] * J;
			//printf("J:%le\n", J);
		}
	}
}


////////////////////////////////////////////////////////////////
//////////////////分布荷重//////////////////////////////////////
////////////////////////////////////////////////////////////////


//面番号ごとの節点番号の取り方の指定
/*void Force_Dis_NodeOfElement( int Number, int DistributeForce[MAX_N_DISTRIBUTE_FORCE][3], int ForceDis_NoE[No_Control_point_ON_ELEMENT[Element_patch[El_No]] ){
	int i;
	if( DistributeForce[Number][1] == 0 )
		for(i = 0; i < No_Control_point_ON_ELEMENT[iPatch]; i++ )
			ForceDis_NoE[i] = Controlpoint_of_Element[ DistributeForce[Number][0] ][i];
	else if( DistributeForce[Number][1] == 1){
		ForceDis_NoE[0] = Controlpoint_of_Element[ DistributeForce[Number][0] ][1];		ForceDis_NoE[1] = Controlpoint_of_Element[ DistributeForce[Number][0] ][2];
		ForceDis_NoE[2] = Controlpoint_of_Element[ DistributeForce[Number][0] ][3];		ForceDis_NoE[3] = Controlpoint_of_Element[ DistributeForce[Number][0] ][0];
	}
	else if( DistributeForce[Number][1] == 2){
		ForceDis_NoE[0] = Controlpoint_of_Element[ DistributeForce[Number][0] ][2];		ForceDis_NoE[1] = Controlpoint_of_Element[ DistributeForce[Number][0] ][3];
		ForceDis_NoE[2] = Controlpoint_of_Element[ DistributeForce[Number][0] ][0];		ForceDis_NoE[3] = Controlpoint_of_Element[ DistributeForce[Number][0] ][1];
	}
	else if( DistributeForce[Number][1] == 3){
		ForceDis_NoE[0] = Controlpoint_of_Element[ DistributeForce[Number][0] ][3];		ForceDis_NoE[1] = Controlpoint_of_Element[ DistributeForce[Number][0] ][0];
		ForceDis_NoE[2] = Controlpoint_of_Element[ DistributeForce[Number][0] ][1];		ForceDis_NoE[3] = Controlpoint_of_Element[ DistributeForce[Number][0] ][2];
	}

}
void Force_Dis( int Total_DistributeForce, int DistributeForce[MAX_N_DISTRIBUTE_FORCE][3], double Val_DistributeForce[MAX_N_DISTRIBUTE_FORCE],
					 int *Total_Load,int Load_Node_Dir[MAX_N_LOAD][2],double Value_of_Load[MAX_N_LOAD],int Total_Control_Point, int El_No, int *Total_Element ){


	int i,j, DF;
	static int ForceDistribute_Controlpoint_of_Element[MAX_NO_CCpoint_ON_ELEMENT];
	static double Out_Force_Distribute[MAX_NO_CCpoint_ON_ELEMENT];
	double a[DIMENSION][DIMENSION];
	static double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION];
	double G = pow( 0.6 , 0.5 );
	double J = 0.0,Sum=0.0;
	double w[DISTRIBUTE_FORCE_Ng] = {(5.0/9.0),(8.0/9.0),(5.0/9.0)};
	double Gxi[DISTRIBUTE_FORCE_Ng][DIMENSION] = { { 1.0, G },{1.0,0.0},{ 1.0, (-1.0)*G } };



	for( DF = 0; DF < Total_DistributeForce; DF++){
		//回転させた要素の節点座標の取得
		Force_Dis_NodeOfElement(DF, DistributeForce, ForceDistribute_Controlpoint_of_Element,El_No);

		for( i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++ ){
			for(j = 0; j < DIMENSION; j++ )
				X[i][j] = Node_Coordinate[ ForceDistribute_Controlpoint_of_Element[i] ][j];
			Out_Force_Distribute[i] = 0.0;
		}
		//各節点ごとの力の計算
		for( j = 0; j < DISTRIBUTE_FORCE_Ng; j++ ){
			Jacobian( a, Gxi[j], X , El_No,Total_Control_Point);
			J = pow(a[0][1],2)+pow(a[1][1],2);
			J = sqrt(J);
			for(i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++ )
				Out_Force_Distribute[i] +=Shape_func(i,Total_Control_Point,Gxi[j],El_No) * Val_DistributeForce[DF] * J * w[j];
		}
		//求めた値を荷重(Load)に追加
		for( i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++ ){
			if( Out_Force_Distribute[i] != 0 ){
				Load_Node_Dir[ *Total_Load ][0] = ForceDistribute_Controlpoint_of_Element[i];
				Load_Node_Dir[ *Total_Load ][1] = DistributeForce[DF][2];
				Value_of_Load[ *Total_Load ]  = Out_Force_Distribute[i];
				*Total_Load += 1;
				Sum += Out_Force_Distribute[i];
			}
		}
		printf("DistributeForce :Element= %d Factor= %d Direction= %d SumForce= %le\n",
			DistributeForce[DF][0],DistributeForce[DF][1],DistributeForce[DF][2],Sum);
		Sum=0.0;

	}

}*/


void Make_Parameter_z(int Total_Element, double E, double nu, int DM)
{
	int e, k;
	Make_gauss_array(0);

	if (DM == 0)
	{
		//Make_strain_z
		for (e = 0; e < Total_Element; e++)
			for (k = 0; k < GP_2D; k++)
				Strain[e][k][3] = 0.0;

		for (e = 0; e < Total_Element; e++)
			for (k = 0; k < GP_2D; k++)
				Strain[e][k][3] = -1.0 * nu / E * (Stress[e][k][0] + Stress[e][k][1]);
	}

	if (DM == 1)
	{
		//Make_stree_z
		for (e = 0; e < Total_Element; e++)
			for (k = 0; k < GP_2D; k++)
				Stress[e][k][3] = 0.0;

		for (e = 0; e < Total_Element; e++)
			for (k = 0; k < GP_2D; k++)
				Stress[e][k][3] = E * nu / (1.0 + nu) / (1 - 2.0 * nu) * (Strain[e][k][0] + Strain[e][k][1]);
	}
}


void Make_Parameter_z_overlay(int Total_Element, double E, double nu, int DM)
{
	int e, k;

	if (DM == 0)
	{
		//Make_strain_z
		for (e = 0; e < Total_Element; e++)
			for (k = 0; k < POW_Ng; k++)
				Strain_overlay[e][k][3] = 0.0;

		for (e = 0; e < Total_Element; e++)
			for (k = 0; k < POW_Ng; k++)
				Strain_overlay[e][k][3] = -1.0 * nu / E * (Stress_overlay[e][k][0] + Stress_overlay[e][k][1]);
	}

	if (DM == 1)
	{
		//Make_stree_z
		for (e = 0; e < Total_Element; e++)
			for (k = 0; k < POW_Ng; k++)
				Stress_overlay[e][k][3] = 0.0;

		for (e = 0; e < Total_Element; e++)
			for (k = 0; k < POW_Ng; k++)
				Stress_overlay[e][k][3] = E * nu / (1.0 + nu) / (1 - 2.0 * nu) * (Strain_overlay[e][k][0] + Strain_overlay[e][k][1]);
	}
}


void element_coordinate(int Total_Element)
{
	int i, j, k, e, l=0;
	double element_edge[9][DIMENSION] = {{-1.0, -1.0}, {1.0, -1.0}, {1.0, 1.0}, {-1.0, 1.0}, {0.0, -1.0}, {1.0, 0.0}, {0.0, 1.0}, {-1.0, 0.0}, {0.0, 0.0}};

	for (e = 0; e < Total_Element; e++)
	{
		for (k = 0; k < 9; k++)
		{
			double data_result_shape[2] = {0.0};
			for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[e]]; i++)
			{
				double R_shape_func = Shape_func(i, element_edge[k], e);
				for (j = 0; j < DIMENSION; j++)
				{
					data_result_shape[j] += R_shape_func * Node_Coordinate[Controlpoint_of_Element[e][i]][j];
				}
			}
			element_coordinate_Nopoint[l][0] = data_result_shape[0];
			element_coordinate_Nopoint[l][1] = data_result_shape[1];
			l++;
		}
	}
	for (l = 0; l < 9 * Total_Element; l++)
		same_point_in_Element[l] = l;
	for (l = 0; l < 9 * Total_Element; l++)
	{
		for (i = l - 1; i >= 0; i--)
		{
			if (element_coordinate_Nopoint[l][0] == element_coordinate_Nopoint[i][0] && element_coordinate_Nopoint[l][1] == element_coordinate_Nopoint[i][1])
			{
				//printf("同じ座標の番号l:%d i:%d\n",l,i);
				same_point_in_Element[l] = i;
			}
		}
	}
}


/*
void ourput_graph_2D(FILE *fp, int e, double element_gg, double element_ee, double data_result_shape_x, double data_result_shape_y, double data_result_disp_x, double data_result_disp_y)
{
	fp = fopen("NURBS/NURBS_disp.dat", "a");
	fprintf(fp, "%d %20.13e  %20.13e %20.13e %20.13e %20.13e %20.13e\n", e, element_gg, element_ee, data_result_shape_x, data_result_shape_y, data_result_disp_x, data_result_disp_y);
	fclose(fp);
}
*/


void calculate_Controlpoint_using_NURBS(double element[DIMENSION], int Total_Element)
{
	int e, b, j, re, i;
	//int p = 0;
	//for(e=0; e < Total_Element; e++){

	double R_shape_func;

	for (re = 0; re < Total_Element; re++)
	{
		e = real_element[re];
		//printf("\n");
		// printf("Element_No:%d\n",e );
		double element_gg = 0.0, element_ee = 0.0, element_delta;

		int i_gg, i_ee;
		// int element_ndiv = 10; 分割数
		int element_ndiv = 1;

		No_points_for_colored_points = (element_ndiv + 1) * (element_ndiv + 1) * real_Total_Element;

		element_delta = 2.0 / element_ndiv;

		for (i_ee = 0; i_ee < element_ndiv + 1; i_ee++)
			for (i_gg = 0; i_gg < element_ndiv + 1; i_gg++)
			{
				double data_result_shape[3] = {0.0};
				double data_result_disp[3] = {0.0};

				element_gg = -1.0 + element_delta * i_gg;
				element_ee = -1.0 + element_delta * i_ee;
				element[0] = element_gg;
				element[1] = element_ee;

				//printf("element_gg:%le\n",element_gg);
				//printf("element_ee:%le\n",element_ee);

				for (b = 0; b < No_Control_point_ON_ELEMENT[Element_patch[e]]; b++)
				{
					R_shape_func = Shape_func(b, element, e);
					for (j = 0; j < DIMENSION; j++)
					{
						data_result_disp[j] += R_shape_func * Displacement[Controlpoint_of_Element[e][b] * DIMENSION + j];
						data_result_shape[j] += R_shape_func * Node_Coordinate[Controlpoint_of_Element[e][b]][j];
					}
				}

				fp = fopen("NURBS/NURBS_points.dat", "a");
				fprintf(fp, "%lf	%lf\n",
							data_result_shape[0],
							data_result_shape[1]);
				fclose(fp);

				fp = fopen("shapefunc/shape_func_xi.dat", "a");
				fprintf(fp, "%lf	%lf\t", element[0], Position_Data_param[0]);
				for (i = 0; i < No_Control_point[Element_patch[e]][0]; i++)
				{
					fprintf(fp, "%lf\t", Shape[0][i][Order[Element_patch[e]][0]]);
				}
				fprintf(fp, "\n");
				fclose(fp);

				fp = fopen("shapefunc/shape_func_eta.dat", "a");
				fprintf(fp, "%lf	%lf\t", element[1], Position_Data_param[1]);
				for (i = 0; i < No_Control_point[Element_patch[e]][1]; i++)
				{
					fprintf(fp, "%lf\t", Shape[1][i][Order[Element_patch[e]][1]]);
				}
				fprintf(fp, "\n");
				fclose(fp);

				//data_result_shape_x[p]=data_result_shape[0];
				//data_result_shape_y[p]=data_result_shape[1];
				//p++;

				//printf("\n");

				//ourput_graph_2D(fp, e, element_gg, element_ee, data_result_shape[0], data_result_shape[1], data_result_disp[0], data_result_disp[1]);

				/*
				data_result_shape_x[p] = data_result_shape[0];
				data_result_shape_y[p] = data_result_shape[1];
				data_result_disp_x[p] = data_result_disp[0];
				data_result_disp_y[p] = data_result_disp[1];
				p++;
				*/

				//NURBS_points(fp,No_points_for_colored_points e,data_result_shape[0],data_result_shape[1])
			}
	}
}


void Gausspoint_coordinate(int Total_Element)
{
	int i, j, k, e;

	Make_gauss_array(0);

	for (e = 0; e < Total_Element; e++)
	{
		for (k = 0; k < GP_2D; k++)
		{
			double data_result_shape[2] = {0.0};
			for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[e]]; i++)
			{
				double R_shape_func = Shape_func(i, Gxi[k], e);
				for (j = 0; j < DIMENSION; j++)
				{
					data_result_shape[j] += R_shape_func * Node_Coordinate[Controlpoint_of_Element[e][i]][j];
				}
			}
			Gausspoint_coordinates[e][k][0] = data_result_shape[0];
			Gausspoint_coordinates[e][k][1] = data_result_shape[1];
		}
	}
}


void calculate_extendmesh_using_NURBS(double element_emsh[DIMENSION])
{
	int e, b, j, re;

	for (re = 0; re < real_Total_Element; re++)
	{
		e = real_element[re];
		//printf("\n");
		//printf("Element_No:%d\n",e );
		double element_gg = 0.0, element_ee = 0.0, element_delta_ee, element_delta_gg;

		int i_gg, i_ee, element_ndiv_ee = 10, element_ndiv_gg = 1;

		No_points_for_new_zarusoba = (element_ndiv_ee + 1) * (element_ndiv_gg + 1) * real_Total_Element;

		element_delta_ee = 2.0 / element_ndiv_ee;
		element_delta_gg = 2.0 / element_ndiv_gg;

		for (i_ee = 0; i_ee < element_ndiv_ee + 1; i_ee++)
			for (i_gg = 0; i_gg < element_ndiv_gg + 1; i_gg++)
			{
				double data_result_shape[3] = {0.0};
				double data_result_disp[3] = {0.0};

				element_gg = -1.0 + element_delta_gg * i_gg;
				element_ee = -1.0 + element_delta_ee * i_ee;
				element_emsh[0] = element_gg;
				element_emsh[1] = element_ee;

				// printf("element_gg:%le\n",element_gg);
				//printf("element_ee:%le\n",element_ee);

				for (b = 0; b < No_Control_point_ON_ELEMENT[Element_patch[e]]; b++)
				{
					double R_shape_func = Shape_func(b, element_emsh, e);
					for (j = 0; j < DIMENSION; j++)
					{
						data_result_disp[j] += R_shape_func * Displacement[Controlpoint_of_Element[e][b] * DIMENSION + j];
						data_result_shape[j] += R_shape_func * Node_Coordinate[Controlpoint_of_Element[e][b]][j];
					}
				}

				//data_result_shape_x[p]=data_result_shape[0];
				//data_result_shape_y[p]=data_result_shape[1];
				//p++;

				//printf("\n");
				/*
				data_result_shape_x_for_new_zarusoba[p] = data_result_shape[0];
				data_result_shape_y_for_new_zarusoba[p] = data_result_shape[1];
				data_result_disp_x_for_new_zarusoba[p] = data_result_disp[0];
				data_result_disp_y_for_new_zarusoba[p] = data_result_disp[1];
				p++;
				*/

				//NURBS_points(fp,No_points_for_colored_points e,data_result_shape[0],data_result_shape[1])
			}
	}
}


//Newton Raphsonによって出力されたxi,etaから重なる要素を求める
int ele_check(int patch_n, double para_coord[DIMENSION])
{
	int i;
	int j;
	int k,kk;
	int l,ll;
    int RangeCheck_flag;						//要素を求め終えたら立てるフラグ
    int temp_ad[DIMENSION][MAX_N_ORDER + 1];	//要素の位置を求めるための値
	int No_line[DIMENSION];						//xi,etaが含まれている要素の列数
	int n = 1;

	for(j = 0; j < DIMENSION; j++)
    {
		//初期化
    	l = 0;
		No_line[j] = 0;
		//printf("No_line[%d]=%d\n",j,No_line[j]);
		// temp_ad[j][MAX_N_ORDER + 1] =  0;
		for(i = 0; i < MAX_N_ORDER + 1; i++)
		{
			temp_ad[j][i] = 0;
		}
	 	RangeCheck_flag = 0;

    	for(k = 0; k < No_Control_point[patch_n][j]-1; k++)
    	{
            if (RangeCheck_flag == 1) break;
			//Local要素の頂点がGlobalパッチ内にない場合
			if (para_coord[j] < Position_Knots[patch_n][j][0] || para_coord[j] > Position_Knots[patch_n][j][No_knot[patch_n][j]-1])
			{
				//printf("no over element\n");
				RangeCheck_flag++;
			}
			//Local要素の頂点がGlobal要素内部にある場合
	        if ( para_coord[j] < Position_Knots[patch_n][j][Order[patch_n][j]+k] )
            {
                //printf("if\nPosition_Knots[%d][%d][%d]=%le\n",
                //        patch_n,j,Order[patch_n][j]+k,Position_Knots[patch_n][j][Order[patch_n][j]+k]);
                int kk = 0;
                for (kk = 0; kk < k + 1; kk++)
                {
                    //if ( RangeCheck_flag == 1 ) break;
                    if (para_coord[j] > Position_Knots[patch_n][j][Order[patch_n][j] + k - kk])
                    {
                        //printf("ifif\nPosition_Knots[%d][%d][%d]=%le\n",
                        //        patch_n,j,Order[patch_n][j]+k-kk,
						//		Position_Knots[patch_n][j][Order[patch_n][j]+k-kk]);
                        temp_ad[j][l] = k - kk;
						//printf("temp_ad[%d][%d]=%d\n",j,l,temp_ad[j][l]);
						l++;
                        RangeCheck_flag++;
                        break;
                    }
                }
            }
			//Local要素の頂点がGlobal要素境界上にある場合
			if (para_coord[j] == Position_Knots[patch_n][j][Order[patch_n][j]+k])
			{
				//頂点の座標がGlobalパッチの始点上にある場合
				if (para_coord[j] == Position_Knots[patch_n][j][0])
				{
					temp_ad[j][l] = k;
					//printf("start point\n");
					//printf("temp_ad[%d][%d]=%d\n",j,l,temp_ad[j][l]);
					l++;
					break;
				}
				//頂点の座標がGlobalパッチの終点上にある場合
				if (para_coord[j] == Position_Knots[patch_n][j][No_knot[patch_n][j]-1])
				{
					temp_ad[j][l] = k - 1;
					//printf("finish point\n");
					//printf("temp_ad[%d][%d]=%d\n",j,l,temp_ad[j][l]);
					l++;
					break;
				}
				//頂点の座標がGlobal要素境界上にある場合
				else
				{
					temp_ad[j][l] = k - 1;
					//printf("temp_ad[%d][%d]=%d\n",j,l,temp_ad[j][l]);
					//printf("data_result_shape[%d]=%le in element_line[%d] on patch[%d] on mesh[0]\n",
					//		j,phys_coord[j],k-1,patch_n);
					l++;
					temp_ad[j][l] = k;
					//printf("temp_ad[%d][%d]=%d\n",j,l,temp_ad[j][l]);
					l++;
					//break;
				}
				for(kk = 0; kk < Order[patch_n][j]; kk++)
				{
					if(para_coord[j] == Position_Knots[patch_n][j][Order[patch_n][j]+k+kk+1])
					//多重ノット（次数分ループ）
					{
						printf("C0 continuity\n");
						temp_ad[j][l] = k + kk;
						//printf("temp_ad[%d][%d]=%d\n",j,l,temp_ad[j][l]);
						l++;
					}
					if(para_coord[j] != Position_Knots[patch_n][j][Order[patch_n][j] + k + kk + 1]) break;
				}
				RangeCheck_flag++;
			}
        }
		No_line[j] = l;
        // printf("No_line[%d]=%d\n", j, No_line[j]);
		n *= l;	//各方向のNo_lineを掛け合わせる
    }
	/*
	for (j = 0; j < DIMENSION; j++)
	{
		printf("No_line[%d];%d\n",j,No_line[j]);
	}*/
	if (DIMENSION == 2)
	{
		for (l = 0; l < No_line[1]; l++)
		{
			for (ll = 0; ll < No_line[0]; ll++)
			{
				temp_element_n[l*No_line[0]+ll] = temp_ad[0][ll]
  											    + temp_ad[1][l] * line_No_Total_element[patch_n][0];
				//printf("el[%d](x,y)=(%le,%le) in element[%d] on patch[%d] on mesh[0]\n"
				//		,element_n_over
				//		,phys_coord[0],phys_coord[1]
				//		,temp_element_n[l*No_line[0]+ll],patch_n);
				//printf("temp_element_n[%d];%d\n",
				//		l*No_line[0]+ll,temp_element_n[l*No_line[0]+ll]);
			}
		}
	}
	// printf("n_element_over;%d\n", n);	//重なっている要素の総数
	return n;
}


//昇順ソート
void sort(int total)
{
	int i,j;
	int tmp;

	for (i = 0; i < total; i++)
	{
		for (j = i + 1; j < total; j++)
		{
			if (element_n_point[i] > element_n_point[j])
			{
				tmp = element_n_point[i];
				element_n_point[i] = element_n_point[j];
				element_n_point[j] = tmp;
			}
		}
	}
}


//重複削除
int duplicate_delete(int total, int element_n)
{
	int i,j;

	j = 0;
	NELOVER[element_n][j] = element_n_point[0];
	j++;
	for (i = 1; i < total; i++)
	{
		if (element_n_point[i] != element_n_point[i-1])
		{
			NELOVER[element_n][j] = element_n_point[i];
			//printf("NELOVER\n");
			j++;
		}
	}
	//j;要素element_nに重なる要素の総数
	return j;
}


//coupled matrix求める
//要素の重なりを求める(要素の端点から求める)
void Check_coupled_Glo_Loc_element_for_end(double element_loc[DIMENSION], int mesh_n_over, int mesh_n_org)
{
    int re;
    int e;
    int b,j;
	int i;
    int k;
	int l, ll;
	int n_elements_over_point[MAX_N_POINT_OVER];
	//int n_points;
	//int n_elements; // = 0;
	// int patch_n, itr_n;
	int patch_n = 0, itr_n = 0;

	int element_ndiv = 1;
	/*n_points = (element_ndiv + 1) * (element_ndiv + 1)
		     * real_Total_Element_on_mesh[mesh_n_over];*/

	/*for (i = 0; i < MAX_N_ELEMENT_OVER; i++)
	{
		temp_element_n[i] = 0;
		printf("i;%d\n",i);
	}*/
	//printf("i;%d\n",i);
	//ローカルパッチ(mesh_n_over)各要素の頂点の物理座標算出
	//from Func.:calculate_Controlpoint_using_NURBS
	for (re = 0; re < real_Total_Element_on_mesh[mesh_n_over]; re++)
	{
		e = real_element[re + real_Total_Element_to_mesh[mesh_n_over]];
		//printf("in_array_real_element:%d\n",
		//		re + real_Total_Element_to_mesh[mesh_n_over+1]);
		//printf("\n");
		//printf("Element_No:%d\n",e );
		double element_gg = 0.0, element_ee = 0.0, element_delta;

		int i_gg, i_ee;

		double output_para[DIMENSION];
		int Total_n_elements = 0;

		k = 0;
		ll = 0;
		element_delta = 2.0 / element_ndiv;

		for (i_ee = 0; i_ee < element_ndiv + 1; i_ee++)
		{
			for (i_gg = 0; i_gg < element_ndiv + 1; i_gg++)
			{
				double data_result_shape[3] = {0.0};

				element_gg = -1.0 + element_delta * i_gg;
				element_ee = -1.0 + element_delta * i_ee;
				element_loc[0] = element_gg;
				element_loc[1] = element_ee;

				//printf("element_gg:%le\n",element_gg);
				//printf("element_ee:%le\n",element_ee);

                //printf("e=%d\n",e);
                //printf("element_patch[e]=%d\n",Element_patch[e]);

				for (b = 0; b < No_Control_point_ON_ELEMENT[Element_patch[e]]; b++)
				{
					double R_shape_func = Shape_func(b, element_loc, e);
					for (j = 0; j < DIMENSION; j++)
					{
                        data_result_shape[j] += R_shape_func * Node_Coordinate[Controlpoint_of_Element[e][b]][j];
                    }
				}
                //printf("data_result_shape[0]=%le\n",data_result_shape[0]);
                //printf("data_result_shape[1]=%le\n",data_result_shape[1]);
				//算出したローカルパッチ各要素の頂点の物理座標のグローバルパッチでの(xi,eta)算出
				//from NURBSviewer/NURBS_view/clickcalc.c/func.:calcXiEtaByNR
				for (i = 0; i < Total_Patch_on_mesh[mesh_n_org]; i++)
				{
					int ii = Calc_xi_eta(data_result_shape[0], data_result_shape[1],
	            				         Position_Knots[i][0], Position_Knots[i][1],
	                			         No_Control_point[i][0], No_Control_point[i][1], Order[i][0], Order[i][1],
	                			         &output_para[0], &output_para[1]);
    			    //printf("Newton_iteration;%d\n",ii);
					//printf("patch: %d\n", i);
					//printf("  x: % 1.8e\n", data_result_shape[0]);
					//printf("  y: % 1.8e\n", data_result_shape[1]);
    			    //printf(" xi: % 1.8e\n", output_para[0]);
     	 	 	    //printf("eta: % 1.8e\n", output_para[1]);
					patch_n = i;
					itr_n = ii;
				}
				//Newton Laphsonによって出力されたxi,etaから重なる要素を求める
			    n_elements_over_point[k] = ele_check(patch_n, output_para);
				//printf("itr_n;%d\n",itr_n);
			    if (itr_n == 0)	//data_result_shapeがグローバルメッシュ上にないとき
				{
					n_elements_over_point[k] = 0;
				}
				//printf("n_elements_over_point[%d];%d\n",
				//		k,n_elements_over_point[k]);
				Total_n_elements += n_elements_over_point[k];
				//printf("Total_n_elements;%d\n",Total_n_elements);
				for (l = 0; l < n_elements_over_point[k]; l++)
				{
					element_n_point[ll] = temp_element_n[l];
					//printf("element_n_point[%d]=%d\n",
					//		ll,element_n_point[ll]);
					ll++;
				}
				k++;
                    //printf("(x,y)=(%le,%le) in element[%d] on patch[%d] on mesh[0]\n",
                    //        data_result_shape[0],data_result_shape[1],
                    //        temp_ad[0]+temp_ad[1]*line_No_Total_element[i][0],i);
				//}

			}
		}
		//printf("Total_n_elements;%d\n",Total_n_elements);
		if (Total_n_elements > 0)
		{
			//昇順ソート
			sort(Total_n_elements);
			//重複削除
			NNLOVER[e] = duplicate_delete(Total_n_elements,e);
		}
		else
		{
			NNLOVER[e] = 0;
		}
		/*
		for (i = 0; i < Total_n_elements; i++)
		{
			printf("element_n_point[%d]=%d\n",i,element_n_point[i]);
		}
		*/
		printf("NNLOVER[%d]=%d\n",e,NNLOVER[e]);

		for (i = 0; i < NNLOVER[e]; i++)
		{
			printf("NELOVER[%d][%d]=%d\n",e,i,NELOVER[e][i]);
		}

	}
}


//要素の重なりを求める(要素のガウス点から求める)
void Check_coupled_Glo_Loc_element_for_Gauss(double element_loc[DIMENSION], int mesh_n_over, int mesh_n_org)
{
    int re;
    int e;
    int i, j, k, m;
	int b;
	int l, ll;
	int n_elements_over_point[MAX_N_POINT_OVER];
	// int patch_n, itr_n;
	int patch_n = 0, itr_n = 0;

	int Check_coupled_No[MAX_N_ELEMENT_OVER];
	double Percent_Check_coupled_No;
	int MAX_NNLOVER = 0;

	//int gauss_1dir = 3;	//重なり判定のための一方向ガウス点数
	//int no_gauss_pt = gauss_1dir * gauss_1dir;	//重なり判定のためのガウス点総数

	for (i = 0; i < MAX_N_ELEMENT_OVER; i++)
	{
		Check_coupled_No[i] = 0;
	}

	for (m = 0; m < 2; m++) //最初Ng個のガウス点で重なりを求め，NNLOVER[e]>=2のeに対して，再度10個のガウス点で重なりを求める
	{
		Make_gauss_array(m);
		
		//ローカルパッチ(mesh_n_over)各要素の頂点の物理座標算出
		//from Func.:calculate_Controlpoint_using_NURBS
		for (re = 0; re < real_Total_Element_on_mesh[mesh_n_over]; re++)
		{
			e = real_element[re + real_Total_Element_to_mesh[mesh_n_over]];
			//printf("Element_No:%d\n",e );

			int i_gg, i_ee;
			int g_n;

			double output_para[DIMENSION];
			int Total_n_elements;

			if (m == 0 || (m == 1 && NNLOVER[e] >= 2))
			{
				//printf("-----------------------------------------m:%d Element_No:%d-----------------------------------------\n",m,e );
				if (m == 1)
				{
					NNLOVER[e] = 0;
					for (i = 0; i < NNLOVER[e]; i++)
					{
						NELOVER[e][i] = 0;
					}
				}

				Total_n_elements = 0;
				k = 0;
				ll = 0;

				for (i_ee = 0; i_ee < GP_1dir; i_ee++)
				{
					for (i_gg = 0; i_gg < GP_1dir; i_gg++)
					{
						double data_result_shape[3] = {0.0};

						g_n = i_ee * GP_1dir + i_gg;
						element_loc[0] = Gxi[g_n][0];
						element_loc[1] = Gxi[g_n][1];

						//printf("e=%d\n",e);
						//printf("element_patch[e]=%d\n",Element_patch[e]);

						for (b = 0; b < No_Control_point_ON_ELEMENT[Element_patch[e]]; b++)
						{
							double R_shape_func = Shape_func(b, element_loc, e);
							for (j = 0; j < DIMENSION; j++)
							{
								data_result_shape[j] += R_shape_func * Node_Coordinate[Controlpoint_of_Element[e][b]][j];
							}
						}
						
						//算出したローカルパッチ各要素の頂点の物理座標のグローバルパッチでの(xi,eta)算出
						//from NURBSviewer/NURBS_view/clickcalc.c/func.:calcXiEtaByNR
						for (i = 0; i < Total_Patch_on_mesh[mesh_n_org]; i++)
						{
							int ii = Calc_xi_eta(data_result_shape[0], data_result_shape[1],
												Position_Knots[i][0], Position_Knots[i][1],
												No_Control_point[i][0], No_Control_point[i][1], Order[i][0], Order[i][1],
												&output_para[0], &output_para[1]);
							//printf("Newton_iteration;%d\n",ii);
							//printf("patch: %d\n", i);
							//printf("  x: % 1.8e\n", data_result_shape[0]);
							//printf("  y: % 1.8e\n", data_result_shape[1]);
							//printf(" xi: % 1.8e\n", output_para[0]);
							//printf("eta: % 1.8e\n", output_para[1]);
							patch_n = i;
							itr_n = ii;
						}
						//Newton Laphsonによって出力されたxi,etaから重なる要素を求める
						n_elements_over_point[k] = ele_check(patch_n, output_para);
						//printf("itr_n;%d\n",itr_n);
						if (itr_n == 0)	//data_result_shapeがグローバルメッシュ上にないとき
						{
							n_elements_over_point[k] = 0;
						}
						// printf("n_elements_over_point[%d];%d\n",
						// 		k,n_elements_over_point[k]);
						Total_n_elements += n_elements_over_point[k];
						//printf("Total_n_elements;%d\n",Total_n_elements);
						for (l = 0; l < n_elements_over_point[k]; l++)
						{
							element_n_point[ll] = temp_element_n[l];
							//printf("element_n_point[%d]=%d\n",
							//		ll,element_n_point[ll]);
							ll++;
						}
						k++;
							//printf("(x,y)=(%le,%le) in element[%d] on patch[%d] on mesh[0]\n",
							//        data_result_shape[0],data_result_shape[1],
							//        temp_ad[0]+temp_ad[1]*line_No_Total_element[i][0],i);
						//}

					}
				}
				// printf("Total_n_elements;%d\n",Total_n_elements);

				//昇順ソート
				sort(Total_n_elements);
				//重複削除
				NNLOVER[e] = duplicate_delete(Total_n_elements, e); //NNLOVER:要素eに重なる要素の総数

				/*
				for (i = 0; i < Total_n_elements; i++)
				{
					printf("element_n_point[%d]=%d\n",i,element_n_point[i]);
				}
				*/
				// printf("NNLOVER[%d]=%d\n",e,NNLOVER[e]);

				// for (i = 0; i < NNLOVER[e]; i++)
				// {
				// 	// printf("NELOVER[%d][%d]=%d\n",e,i,NELOVER[e][i]);
				// }
			}
		}
	}

	for (re = 0; re < real_Total_Element_on_mesh[mesh_n_over]; re++)
	{
		e = real_element[re + real_Total_Element_to_mesh[mesh_n_over]];
		printf("-----------------------------------------Element_No:%d-----------------------------------------\n",e );

		Check_coupled_No[NNLOVER[e]]++;

		if (MAX_NNLOVER <  NNLOVER[e])
		{
			MAX_NNLOVER =  NNLOVER[e];
		}

		printf("NNLOVER[%d] = %d\n", e, NNLOVER[e]);

		for (i = 0; i < NNLOVER[e]; i++)
		{
			printf("\tNELOVER[%d][%d] = %d\n", e, i, NELOVER[e][i]); //要素eに重なるi番目の要素番号
		}
	}

	printf("MAX_NNLOVER = %d\n", MAX_NNLOVER);
	for (i = 0; i <= MAX_NNLOVER; i++)
	{
		Percent_Check_coupled_No = (double)Check_coupled_No[i] * 100.0 / (double)real_Total_Element_on_mesh[mesh_n_over];
		printf("Check_coupled_No[%d] = %d\t\t%3.1lf %%\n",i , Check_coupled_No[i], Percent_Check_coupled_No);
	}
	printf("---------------------------------------------------------------------------------------------------------------------------\n");
}


void Make_Loc_Glo()
{
	int i,j,k;
	int jj;
	int e;
	int j_n;
	int count;

	j_n = real_Total_Element_to_mesh[Total_mesh] - real_Total_Element_on_mesh[0];
	//printf("j_n=%d\n",j_n);
	//printf("%d\t%d\n",)
	for (i = 0; i < real_Total_Element_on_mesh[0]; i++)
	{
		e = real_element[i];
		count = 0;

		for(j = 0; j < j_n ; j++)
		{
			jj = real_element[real_Total_Element_to_mesh[1] + j];	//ローカルメッシュ上のreal element番号
			//printf("jj=%d\n",jj);
			if (NNLOVER[jj] > 0)
			{
				//printf("jj=%d\n",jj);
				for (k = 0; k < NNLOVER[jj]; k++)
				{
					if (NELOVER[jj][k] == e)
					{
						NELOVER[e][count] = jj;
						//printf("NELOVER[%d][%d]=%d\n",e,count,
						//							  NELOVER[e][count]);
						count++;
					}
				}
			}
		}
		NNLOVER[e] = count;
		//printf("NNLOVER[%d]=%d\n",e,NNLOVER[e]);
	}
}

/*
////////////////////////////////////////////////////////////////
//////////////////AVS出力///////////////////////////////////////
////////////////////////////////////////////////////////////////

//節点座標と要素の節点番号の書き込み
void AVS_inputInp_Quad_4( int Total_Element, int Total_Control_Point ){
	int i;

	fprintf(fp,"%d	%d\n", Total_Control_Point, Total_Element);	//総節点数、総要素数
	for(i = 0; i < Total_Control_Point; i++ ){
		fprintf(fp,"%d	%e	%e	%e\n", i+1,Node_Coordinate[i][0],Node_Coordinate[i][1],0.0);
														//節点番号(1〜)、ＸＹＺ座標（２次元だとZ=0と記述
	}
	for( i = 0; i < Total_Element; i++ ){
		fprintf(fp,"%d	0	quad	", i+1);					//要素番号(1〜)、材料番号(0)、要素の形
		//要素タイプごとに順番を確認のとこ
		fprintf(fp,"%d	", Controlpoint_of_Element[i][2]+1);
		fprintf(fp,"%d	", Controlpoint_of_Element[i][3]+1);
		fprintf(fp,"%d	", Controlpoint_of_Element[i][0]+1);
		fprintf(fp,"%d	", Controlpoint_of_Element[i][1]+1);
		fprintf(fp,"\n");
	}

}
//計算結果のデータの書き込み
void AVS_inputAns_2D( int Total_Control_Point, int Total_Element ){

	int i,j,k, NodeDataCom = 6;
	int ElementDataCom = 12;	//N_STRAIN + N_STRESS + 4
	double Str;

	fprintf(fp,"%d	%d\n", NodeDataCom, ElementDataCom);	//各節点に存在するデータ数、各要素に存在するデータ数

	//節点のデータの書き込み（変位ＸＹＺ方向＋自由利用分XYZ）
	fprintf(fp,"%d	", NodeDataCom);						//節点データ成分数
	for( i = 0; i < NodeDataCom; i++ )
		fprintf(fp,"1	");								//各成分の構成数
	fprintf(fp,"\n");

	fprintf(fp,"DisX,\nDisY,\nDisZ,\nOpX,\nOpY,\nOpZ,\n");		//各節点データ成分のラベル
	for( i = 0; i < Total_Control_Point; i++ ){
		fprintf(fp,"%d	%e	%le	", i+1, Displacement[i*DIMENSION], Displacement[i*DIMENSION+1]);
		fprintf(fp,"%e	%le	%le	%le\n", 0.0, 0.0, 0.0, 0.0 );
	}

	//要素データの書き込み（応力＋自由利用２、歪＋自由利用２）
	fprintf(fp,"%d	", ElementDataCom);						//要素データ成分数
	for( i = 0; i < ElementDataCom; i++ )
		fprintf(fp,"1	");									//各成分の構成数
	fprintf(fp,"\n");

	fprintf(fp,"SigXX,\nSigYY,\nSigXY,\nSigZZ,\nSigOp1,\nSigOp2,\n");		//各節点データ成分のラベル
	fprintf(fp,"IpuXX,\nIpuYY,\nIpuXY,\nIpuZZ,\nIpuOp1,\nIpuOp2,\n");		//各節点データ成分のラベル

	for( i = 0; i < Total_Element; i++ ){
		fprintf(fp,"%d	", i+1);
		for( j = 0; j < N_STRESS; j++ ){
			Str = 0.0;
			for( k = 0; k < POW_Ng; k++ )		Str += Stress[i][k][j];
			fprintf(fp,"	%e",Str / (double)(POW_Ng) );
		}
		for( ; j < 6; j++ )				fprintf(fp,"	%e", 0.0 );

		for( j = 0; j < N_STRAIN; j++ ){
			Str = 0.0;
			for( k = 0; k < POW_Ng; k++ )
				Str += Strain[i][k][j];
			fprintf(fp,"	%e",Str / (double)(POW_Ng) );
		}
		for( ; j < 6; j++ )				fprintf(fp,"	%e", 0.0 );

		fprintf(fp,"\n");
	}
}

void Make_Output( int Total_Control_Point, int Total_Element ){
	int StepMax=1, StepNo=1;

	//AVS用のinpファイルの制作
	fp = fopen( "AVS/1_1_force_120_0122.inp", "w");
	fprintf(fp,"# AVS field file\n");					//注釈文（必ず先頭に「#」）
	fprintf(fp,"%d\n", StepMax);						//ステップ数の設定
	fprintf(fp,"data\n" );								//データの繰り返しタイプ	data,geom,data_geom

	fprintf(fp,"step%d\n", StepNo);					//ステップ番号
	AVS_inputInp_Quad_4( Total_Element, Total_Control_Point);	//節点座標と要素の取得
	AVS_inputAns_2D( Total_Control_Point, Total_Element );			//各種計算結果の取得
	printf("Finish Make_AVS_Step%d\n",StepNo);

	fclose(fp);
}
*/

/*
////////////////////////////////////////////////////////////////
//////////////////J積分/////////////////////////////////////////
////////////////////////////////////////////////////////////////



//B_xマトリックスを求める関数
int Make_B_x_Matrix_Quad_4(double B_x[DIMENSION][KIEL_SIZE], double Local_coord[DIMENSION], double X[No_Control_point_ON_ELEMENT][DIMENSION], double *J ){
	double a[DIMENSION][DIMENSION], b[DIMENSION][No_Control_point_ON_ELEMENT];
	int i,j,k;

	Jacobian_Quad_4( a, Local_coord, X );

	*J = InverseMatrix_2D( a );
	if( *J <= 0 )return -999;

	for( i = 0; i < DIMENSION; i++ ){
		for( j = 0; j < No_Control_point_ON_ELEMENT; j++ ){
			b[i][j] = 0.0;
			for( k = 0; k < DIMENSION; k++ ){
				b[i][j] += a[k][i] * dN_Quad_4( j, Local_coord, k);
			}
		}
	}
	for( i = 0; i < No_Control_point_ON_ELEMENT; i++ ){
		B_x[0][2*i] = b[0][i];	B_x[0][2*i+1] = 0.0;
		B_x[1][2*i] = 0.0;		B_x[1][2*i+1] = b[0][i];
	}
	return 0;
}







void Make_Strain_x_Quad_4(double E, double nu, int Total_Element){
	static double U[KIEL_SIZE];
	static double B_x[DIMENSION][KIEL_SIZE],X[No_Control_point_ON_ELEMENT][DIMENSION],J;
	double w[POW_Ng] = {1.0,1.0,1.0,1.0};
	double G = 1/pow(3,0.5);
	double Gxi[Total_Element][POW_Ng][DIMENSION] = { {{(1+G)/2,(1-G)/2},{(1+G)/2,(1+G)/2},{(1-G)/2,(1+G)/2},{(1-G)/2,(1-G)/2}},
	 {{(1+G)/2,(1-G)/2},{(1+G)/2,(1+G)/2},{(1-G)/2,(1+G)/2},{(1-G)/2,(1-G)/2}}};

	int N,e,i,j;

	for( e = 0; e < Total_Element; e++ ){
		for( N = 0; N < POW_Ng; N++)
			for( i = 0; i < DIMENSION; i ++ )
				Strain_x[e][N][i] = 0.0;
		//Bマトリックスと各要素の変位を取得
		for( i = 0; i < No_Control_point_ON_ELEMENT; i++ ){
			for( j = 0; j < DIMENSION; j++ ){
				U[ i*DIMENSION +j ] = Displacement[ Controlpoint_of_Element[e][i]*DIMENSION + j ];
				X[i][j] = Node_Coordinate[ Controlpoint_of_Element[e][i] ][j];
			}
		}
		//歪
		for( N = 0; N < POW_Ng; N++ ){
			Make_B_x_Matrix_Quad_4( B_x, Gxi[N], X ,&J );
			for( i = 0; i < DIMENSION; i++ )
				for( j = 0; j < KIEL_SIZE; j++ )
					 Strain_x[e][N][i] += B_x[i][j] * U[j] * w[N];
		}
	}
}












//エネルギーモーメンタムテンソルを求める関数
void Make_EMT(double E, double nu, int Total_Element){

	int i, j, k;
	double W;
	double W_x[DIMENSION];

	int K_D[DIMENSION] = {1, 0};
	double W_K_D[DIMENSION];
	double P_1j[DIMENSION];

	Make_Stress_2D(E, nu, Total_Element);
	Make_Strain_Quad_4(E, nu, Total_Element);
	Make_Strain_x_Quad_4(E, nu, Total_Element);

	W = 0.0;


	for( i = 0; i < DIMENSION; i++){
		W_x[i] = 0.0;
	}


	for( k = 0; k < POW_Ng; k++ ){
		for( i = 0; i < Total_Element; i++ ){
			for( j = 0; j < DIMENSION+1; j++ ){

				W += (1.0/2.0)*Stress[i][k][j]*Strain[i][k][j];
			}
		}
	}
	printf("\nW = %lf\n", W);


	for( i = 0; i < DIMENSION; i++ ){
		W_K_D[i] = W*K_D[i];

		printf("\nW_K_D[%d] = %lf\n", i, W_K_D[i]);
	}


	for( k = 0; k < POW_Ng; k++ ){
		for( i = 0; i < Total_Element; i++ ){
			for( j = 0; j < DIMENSION; j++ ){

				if(j == 0){
					W_x[j] += Stress[i][k][j]*Strain_x[i][k][j]+Stress[i][k][j+2]*Strain_x[i][k][j+1];
				}

				else if(j == 1){
					W_x[j] += Stress[i][k][j+1]*Strain_x[i][k][j-1]+Stress[i][k][j]*Strain_x[i][k][j];
				}
			}
		}
	}


	for( j = 0; j < DIMENSION; j++ ){
		printf("\nW_x[%d] = %le\n", j, W_x[j]);
	}


	for( i = 0; i < DIMENSION; i++ ){
		P_1j[i] = W_K_D[i] - W_x[i];

		printf("\nP_1j[%d] = %lf\n", i, P_1j[i]);
	}
}
*/


int SerchForElement(int mesh_n, int iPatch, int Total_Element, int iX, int iY)
{
	int iii = 0;

	printf("iii%d\n",iii);
	printf("Total_Element%d\n",Total_Element);
	for (iii = 0; iii < Total_Element; iii++)
	{
		// printf("%d\n",iii+Total_Element_to_mesh[mesh_n]);
		// printf("ipatch%d\n",iPatch);
		if (Element_patch[iii+Total_Element_to_mesh[mesh_n]] == iPatch)
		{
			// printf("Check SerchForElement 1 iii = %d\n", iii+Total_Element_to_mesh[mesh_n]);
			// printf("ENC[iPatch][iii][0] = %d ENC[iPatch][iii][1] = %d  iX = %d  iY = %d\n",
			// 		ENC[iPatch][iii+Total_Element_to_mesh[mesh_n]][0],
			// 		ENC[iPatch][iii+Total_Element_to_mesh[mesh_n]][1],
			// 		iX, iY);
			if (iX == ENC[iPatch][iii+Total_Element_to_mesh[mesh_n]][0] && iY == ENC[iPatch][iii+Total_Element_to_mesh[mesh_n]][1])
				goto loopend;
			/* iii --; */

			// printf("Check SerchForElement 2 iii = %d\n", iii);
		}
	}
loopend:

	return (iii);
}


void Setting_Dist_Load_2D(int mesh_n, int iPatch, int Total_Element, int iCoord, double val_Coord,
						  double Range_Coord[2], int type_load, double Coeff_Dist_Load[3])
{
	int iii, jjj;
	// int iDir_Element[MAX_N_KNOT], jDir_Element;
	int N_Seg_Load_Element_iDir = 0, jCoord;
	int iRange_ele[2]/*, jRange_ele*/;
	int iPos[2] = {-10000, -10000}, jPos[2] = {-10000, -10000};
	// double Coord_Seg_Load_Element_iDir[MAX_N_KNOT][2], Coord_Seg_Load_Element_jDir[2];
	int No_Element_for_Integration[MAX_N_KNOT], No_Element_For_Dist_Load;
	int iX, iY;
	//int Element_Integration;
	int iControlpoint[MAX_NO_CCpoint_ON_ELEMENT], ic, ig, NNG = 3;
	double val_jCoord_Local = 0.0;
	double GaussPt[3], Weight[3];
	double Gg = pow(3.0 / 5.0, 0.5);

	/* type_load: 0: Dist load in x direction
 * 	              1:              y direction
 * 	              2:              normal direciton */

	GaussPt[0] = -Gg;
	GaussPt[1] = 0.0;
	GaussPt[2] = Gg;
	Weight[0] = 5.0 / 9.0;
	Weight[1] = 8.0 / 9.0;
	Weight[2] = 5.0 / 9.0;

	/* iCoord=0: Load on Eta=Constant
	   iCoord=1: Load on Xi=Constant */
	if (iCoord == 0)
		jCoord = 1;
	if (iCoord == 1)
		jCoord = 0;

	/* val_Coord: Value of Eta or Xi of the line or surface to give the distributed load */

	/* Setting elements needed to computed the distributed load */

	for (iii = Order[iPatch][iCoord]; iii < No_knot[iPatch][iCoord] - Order[iPatch][iCoord] - 1; iii++)
	{
		double epsi = 0.00000000001;
		/* iPos[0] = -10000; iPos[1] = -10000; jPos[0] = -10000; jPos[1] = -10000;*/
		printf("Check1 iii = %d\n", iii);
		printf("Check2 Position_Knots[iCoord][iii]= %f  Range_Coord[0] =%f Position_Knots[iCoord][iii+1] = %f\n", Position_Knots[iPatch][iCoord][iii], Range_Coord[0], Position_Knots[iPatch][iCoord][iii + 1]);
		/*

		if(Position_Knots[iCoord][iii]-epsi <= Range_Coord[0] &&
			Position_Knots[iCoord][iii+1]+epsi > Range_Coord[0]) iPos[0] = iii;

		if(Position_Knots[iCoord][iii]-epsi <= Range_Coord[1] &&
                        Position_Knots[iCoord][iii+1]+epsi > Range_Coord[1]) iPos[1] = iii+1;
	*/
		if (Position_Knots[iPatch][iCoord][iii] - epsi <= Range_Coord[0])
			iPos[0] = iii;
		if (Position_Knots[iPatch][iCoord][iii + 1] - epsi <= Range_Coord[1])
			iPos[1] = iii + 1;
	}
	iRange_ele[0] = iPos[0] - Order[iPatch][iCoord];
	iRange_ele[1] = iPos[1] - Order[iPatch][iCoord] - 1;
	printf("iPos[0] = %d iPos[1] = %d\n", iPos[0], iPos[1]);
	printf("iRange_ele[0] = %d  iRange_ele[1] = %d\n", iRange_ele[0], iRange_ele[1]);

	if (iPos[0] < 0 || iPos[1] < 0)
	{
		printf("Error (Stop) iPos[0] = %d iPos[1] = %d\n", iPos[0], iPos[1]);
		exit(0);
	}

	for (jjj = Order[iPatch][jCoord]; jjj < No_knot[iPatch][jCoord] - Order[iPatch][jCoord] - 1; jjj++)
	{
		double epsi = 0.00000000001;
		/* jjj=Order[jCoord]; */
		if (Position_Knots[iPatch][jCoord][jjj] - epsi <= val_Coord &&
			Position_Knots[iPatch][jCoord][jjj + 1] + epsi > val_Coord)
		{
			jPos[0] = jjj;
			jPos[1] = jjj + 1;
			val_jCoord_Local = -1.0 + 2.0 * (val_Coord - Position_Knots[iPatch][jCoord][jjj]) /
										  (Position_Knots[iPatch][jCoord][jjj + 1] - Position_Knots[iPatch][jCoord][jjj]);
		}
		//(2019_06_13)printf("Check jjj count: jjj =  %d\n",jjj);
	}
	// jRange_ele = jPos[0] - Order[iPatch][jCoord];
	//(2019_06_13)printf("jPos[0] = %d jPos[1] = %d  jRange_ele = %d val_jCoord_Local = %f\n",jPos[0],jPos[1],jRange_ele, val_jCoord_Local);

	if (jPos[0] < 0 || jPos[1] < 0)
	{
		printf("Error (Stop) jPos[0] = %d jPos[1] = %d\n", jPos[0], jPos[1]);
		exit(0);
	}

	for (iii = iPos[0]; iii < iPos[1]; iii++)
	{
		// Coord_Seg_Load_Element_iDir[iii][0] = Position_Knots[iPatch][iCoord][iii + iPos[0]];
		// Coord_Seg_Load_Element_iDir[iii][1] = Position_Knots[iPatch][iCoord][iii + iPos[0] + 1];
		// iDir_Element[N_Seg_Load_Element_iDir] = iii - Order[iPatch][iCoord];
		N_Seg_Load_Element_iDir++;
	}

	// Coord_Seg_Load_Element_jDir[0] = Position_Knots[iPatch][iCoord][jPos[0]];
	// Coord_Seg_Load_Element_jDir[1] = Position_Knots[iPatch][iCoord][jPos[1]];
	// jDir_Element = jPos[0] - Order[iPatch][iCoord];
	iii = 0;
	if (iCoord == 1)
	{
		// int iX, iY;
		iX = jPos[0] - Order[iPatch][0];
		//(2019_06_13)printf("Check iPos[0] = %d  iPos[1] = %d\n",iPos[0],iPos[1]);
		for (iY = iPos[0] - Order[iPatch][1]; iY < iPos[1] - Order[iPatch][1]; iY++)
		{
			//(2019_06_13)printf("Check iY = %d\n",iY);
			No_Element_for_Integration[iii] = SerchForElement(mesh_n, iPatch, Total_Element, iX, iY);
			printf("Check No_Element_for_Integration[%d] = %d\n"
                    ,iii,No_Element_for_Integration[iii]);//(2019_06_13)
			iii++;
		}
	}

	if (iCoord == 0)
	{
		// int iX, iY;
		iY = jPos[0] - Order[iPatch][1];
		//(2019_06_13)printf("Check iPos[0] = %d  iPos[1] = %d\n",iPos[0],iPos[1]);
		for (iX = iPos[0] - Order[iPatch][0]; iX < iPos[1] - Order[iPatch][0]; iX++)
		{
			//(2019_06_13)printf("Check iX = %d\n",iX);
			No_Element_for_Integration[iii] = SerchForElement(mesh_n, iPatch, Total_Element, iX, iY);
		    printf("Check No_Element_for_Integration[%d] = %d\n"
                   ,iii,No_Element_for_Integration[iii]);//(2019_06_13)
			iii++;
		}
	}
	No_Element_For_Dist_Load = iii;
	//(2019_06_13)printf("No_Element_For_Dist_Load = %d\n",No_Element_For_Dist_Load);

	/* Book keeping finished */

	for (iii = 0; iii < No_Element_For_Dist_Load; iii++)
	{ //B
		//(2019_06_13)printf("Check3 iii = %d\n",iii);
		//(2019_06_13)printf("Total_element_all_ID[No_Element_for_Integration[iii]] = %d\n No_Element_for_Integration[iii] = %d  iii = %d\n",
		//Total_element_all_ID[No_Element_for_Integration[iii]],No_Element_for_Integration[iii],iii);
		if (Total_element_all_ID[No_Element_for_Integration[iii]] == 1)
		{ //A
			iX = ENC[iPatch][No_Element_for_Integration[iii]][0];
			iY = ENC[iPatch][No_Element_for_Integration[iii]][1];
			//(2019_06_13)printf("iX = %d  iY = %d\n",iX, iY);

			for (ic = 0; ic < (Order[iPatch][0] + 1) * (Order[iPatch][1] + 1); ic++)
				iControlpoint[ic] = Controlpoint_of_Element[No_Element_for_Integration[iii]][ic];

			for (ig = 0; ig < NNG; ig++)
			{
				double Local_Coord[2], sfc, dxyzdge[3], detJ, XiEtaCoordParen, valDistLoad;
				int icc;
				Local_Coord[jCoord] = val_jCoord_Local;
				Local_Coord[iCoord] = GaussPt[ig];
				printf("ig = %d   Local_Coord[jCoord] = %f Local_Coord[iCoord] = %f\n", ig, Local_Coord[jCoord], Local_Coord[iCoord]);

				ShapeFunc_from_paren(Local_Coord, iCoord, No_Element_for_Integration[iii]);
				XiEtaCoordParen = Position_Data_param[iCoord];
				printf("Check  Coeff_Dist_Load[0] = %f Coeff_Dist_Load[1] = %f  Coeff_Dist_Load[2] = %f  Position_Data_param[iCoord] = %f\n", Coeff_Dist_Load[0], Coeff_Dist_Load[1], Coeff_Dist_Load[2], Position_Data_param[iCoord]);
				valDistLoad = Coeff_Dist_Load[0] + Coeff_Dist_Load[1] * XiEtaCoordParen + Coeff_Dist_Load[2] * XiEtaCoordParen * XiEtaCoordParen;

				//き裂の分布荷重を与える時のために...
                //printf("XiEtaCoordParen=%lf\n",XiEtaCoordParen);
                //double theta;
                //theta = XiEtaCoordParen*PI/2;
				//theta = XiEtaCoordParen*2*PI/line_No_real_element[0][1];
                //printf("theta=%lf\n",theta*180/PI);
                //valDistLoad = cos(theta);
                //printf("valDistLoad=%2.10lf\n",valDistLoad);

				dxyzdge[0] = 0.0;
				dxyzdge[1] = 0.0;
				dxyzdge[2] = 0.0;
				for (icc = 0; icc < (Order[iPatch][0] + 1) * (Order[iPatch][1] + 1); icc++)
				{
					//printf("dxy0:%lf\tdxy1:%lf\n",dxyzdge[0],dxyzdge[1]);
					//printf("Node%lf\tNode%lf\n",Node_Coordinate[iControlpoint[icc]][0],Node_Coordinate[iControlpoint[icc]][1]);
					dxyzdge[0] += dShape_func(icc, iCoord, Local_Coord, No_Element_for_Integration[iii]) * Node_Coordinate[iControlpoint[icc]][0];
					dxyzdge[1] += dShape_func(icc, iCoord, Local_Coord, No_Element_for_Integration[iii]) * Node_Coordinate[iControlpoint[icc]][1];
					//printf("dxy0:%lf\tdxy1:%lf\n",dxyzdge[0],dxyzdge[1]);
				}

				//printf("dxy")

				detJ = sqrt(dxyzdge[0] * dxyzdge[0] + dxyzdge[1] * dxyzdge[1]);
				//(2019_06_13)
				printf("Check the value of detJ etc: detJ = %f dxyzdge[0] = %f dxyzdge[1] = %f\n",detJ,dxyzdge[0],dxyzdge[1]);
				if (type_load < 2)
				{
					for (ic = 0; ic < (Order[iPatch][0] + 1) * (Order[iPatch][1] + 1); ic++)
					{
						//printf("Order[%d][0];%d,Order[%d][1]:%d\n",iPatch,Order[iPatch][0],iPatch,Order[iPatch][1]);
						sfc = Shape_func(ic, Local_Coord, No_Element_for_Integration[iii]);
						Equivalent_Nodal_Force[iControlpoint[ic]][type_load] += valDistLoad * sfc * detJ * Weight[ig];
						//(2019_06_13)printf("Check ic = %d sfc = %f   Weight[ig] = %f  valDistLoad = %f\n",ic,sfc,Weight[ig],valDistLoad);
						//(2019_06_13)printf("Equivalent_Nodal_Force[%d][%d]:%le\n",iControlpoint[ic],type_load,Equivalent_Nodal_Force[iControlpoint[ic]][type_load] );
					}
				}

				if (type_load == 2)
				{
					double LoadDir[2];
					LoadDir[0] = dxyzdge[1] / detJ;
					LoadDir[1] = -dxyzdge[0] / detJ;
					//printf("%lf\t%lf\n",LoadDir[0],LoadDir[1]);
					for (ic = 0; ic < (Order[iPatch][0] + 1) * (Order[iPatch][1] + 1); ic++)
					{
						sfc = Shape_func(ic, Local_Coord, No_Element_for_Integration[iii]);
						//printf("sfc=%lf\n",sfc);
						Equivalent_Nodal_Force[iControlpoint[ic]][0] +=
							LoadDir[0] * valDistLoad * sfc * detJ * Weight[ig];
						Equivalent_Nodal_Force[iControlpoint[ic]][1] +=
							LoadDir[1] * valDistLoad * sfc * detJ * Weight[ig];
						printf("Equivalent_Nodal_Force[%d][0]=%lf\nEquivalent_Nodal_Force[%d][1]=%lf\n",
								iControlpoint[ic],Equivalent_Nodal_Force[iControlpoint[ic]][0],
								iControlpoint[ic],Equivalent_Nodal_Force[iControlpoint[ic]][1]);
						//printf("LoadDir[0]*(0.5-0.5*cos(2*theta))*sfc*detJ*Weight[%d]=%lf*%lf*%lf*%lf*%lf=%lf\n",ig,LoadDir[0],(0.5-0.5*cos(2*theta)),sfc,detJ,Weight[ig],LoadDir[0] * (0.5-0.5*cos(2*theta)) * sfc * detJ * Weight[ig]);
						//printf("LoadDir[1]*(0.5-0.5*cos(2*theta))*sfc*detJ*Weight[%d]=%lf*%lf*%lf*%lf*%lf=%lf\n",ig,LoadDir[1],(0.5-0.5*cos(2*theta)),sfc,detJ,Weight[ig],LoadDir[1] * (0.5-0.5*cos(2*theta)) * sfc * detJ * Weight[ig]);
					}
				}
                if ( type_load == 3 )
                {
                    double LoadDir[2];
                    LoadDir[0] = dxyzdge[0] / detJ;
                    LoadDir[1] = dxyzdge[1] / detJ;
                    for (ic = 0; ic < (Order[iPatch][0] + 1) * (Order[iPatch][1] + 1); ic++)
                    {
                        sfc = Shape_func(ic, Local_Coord,
                                         No_Element_for_Integration[iii]);
                        Equivalent_Nodal_Force[iControlpoint[ic]][0] +=
                            LoadDir[0] * valDistLoad * sfc * detJ * Weight[ig];
                        Equivalent_Nodal_Force[iControlpoint[ic]][1] +=
                            LoadDir[1] * valDistLoad * sfc * detJ * Weight[ig];
						//printf("Equivalent_Nodal_Force[%d][0]=%lf\nEquivalent_Nodal_Force[%d][1]=%lf\n",iControlpoint[ic],Equivalent_Nodal_Force[iControlpoint[ic]][0],iControlpoint[ic],Equivalent_Nodal_Force[iControlpoint[ic]][1]);
						//printf("LoadDir[0]*sin(theta)*cos(theta)*sfc*detJ*Weight[%d]=%lf*%lf*%lf*%lf*%lf=%lf\n",ig,LoadDir[0],sin(theta)*cos(theta),sfc,detJ,Weight[ig],LoadDir[0] * sin(theta)*cos(theta) * sfc * detJ * Weight[ig]);
						//printf("LoadDir[1]*sin(theta)*cos(theta)*sfc*detJ*Weight[%d]=%lf*%lf*%lf*%lf*%lf=%lf\n",ig,LoadDir[1],sin(theta)*cos(theta),sfc,detJ,Weight[ig],LoadDir[1] * sin(theta)*cos(theta) * sfc * detJ * Weight[ig]);
                    }
                }
			}
		} //A
	}	 //B
}


//重ね合わせた結果の出力(NURBS_input_for_s-IGA)
void GetLocData()
{
	double temp;

	//必要なのはローカルのパッチ数とコントロールポイント数
	printf("Start Get Local Data\n\n");
	fp = fopen("input_local.txt", "r");

	fscanf(fp, "%lf%lf", &temp, &temp);
	fscanf(fp, "\n");
	fscanf(fp, "%d", &n_patch_loc);
	fscanf(fp, "\n");
	fscanf(fp, "%d", &loc_cntl_p_n);
	fscanf(fp, "\n");
	fclose(fp);
	printf("patches(in local):% d\n",n_patch_loc);
	printf("control points(in local):% d\n",loc_cntl_p_n);
	printf("\nFinish Get Local Data\n\n");

}


void ReadFile()
{
	int i, j;
	double temp1, temp2, temp3;
	int temp_int;
	printf("Start Reading input\n\n");
	fp = fopen("input_for_NURBS.txt", "r");

	fscanf(fp, "%lf%lf", &E, &nu);
	fscanf(fp, "\n");
	printf("E nu: % 1.4e % 1.4e\n", E, nu);

	fscanf(fp, "%d", &patch_n);
	fscanf(fp, "\n");
	printf("patches: %d \n", patch_n);
	if (patch_n > MAX_PATCHES) {
		printf("Error!!\n");
		printf("Too many patches!\n"
			   "Maximum of patches is %d (Now %d)\n"
			   "\n", MAX_PATCHES, patch_n);
		exit(1);
	}

	fscanf(fp, "%d", &cntl_p_n);
	fscanf(fp, "\n");
	printf("total control points:%d \n", cntl_p_n);
	if (cntl_p_n > MAX_CNRL_P) {
		printf("Error!!\n");
		printf("Too many control points!\n"
			   "Maximum of control points is %d (Now %d)\n"
			   "\n", MAX_CNRL_P, cntl_p_n);
		exit(1);
	}

	for (i = 0; i < patch_n; i++) {
		fscanf(fp, "%d%d", &order_xi[i], &order_eta[i]);
		fscanf(fp, "\n");
		printf("order %d: %d %d\n", i, order_xi[i], order_eta[i]);
		if (order_xi[i] > MAX_ORDER) {
			printf("Error!!\n");
			printf("Order too big at xi!\n"
				   "Maximum of order is %d (Now %d at patch %d)\n"
				   "\n", MAX_ORDER, order_xi[i], i);
			exit(1);
		}
		if (order_eta[i] > MAX_ORDER) {
			printf("Error!!\n");
			printf("Order too big at eta!\n"
				   "Maximum of order is %d (Now %d at patch %d)\n"
				   "\n", MAX_ORDER, order_eta[i], i);
			exit(1);
		}
	}

	for (i = 0; i < patch_n; i++) {
		fscanf(fp, "%d%d", &knot_n_xi[i], &knot_n_eta[i]);
		fscanf(fp, "\n");
		printf("knots %d: %d %d\n", i, knot_n_xi[i], knot_n_eta[i]);
		if (knot_n_xi[i] > MAX_KNOTS) {
			printf("Error!!\n");
			printf("Knot vector too long at xi!\n"
				   "Maximum of knot vector is %d (Now %d at patch %d)\n"
				   "\n", MAX_KNOTS, knot_n_xi[i], i);
			exit(1);
		}
		if (knot_n_eta[i] > MAX_KNOTS) {
			printf("Error!!\n");
			printf("Knot vector too long at eta!\n"
				   "Maximum of knot vector is %d (Now %d at patch %d)\n"
				   "\n", MAX_KNOTS, knot_n_eta[i], i);
			exit(1);
		}
	}

	for (i = 0; i < patch_n; i++) {
		fscanf(fp, "%d%d", &cntl_p_n_xi[i], &cntl_p_n_eta[i]);
		printf("control points %d: %d %d\n",
			   i, cntl_p_n_xi[i], cntl_p_n_eta[i]);
		fscanf(fp, "\n");
	}
	printf("\n");


	for (i = 0; i < patch_n; i++) {
		for (j = 0; j < cntl_p_n_xi[i] * cntl_p_n_eta[i]; j++) {
			fscanf(fp, "%d", &temp_index[i][j]);
			printf("%d ", temp_index[i][j]);
		}
		fscanf(fp, "\n");
		printf("\n");
	}
	printf("\n");

	fscanf(fp, "%lf%lf%lf", &temp1, &temp2, &temp3);
	fscanf(fp, "\n");

	for (i = 0; i < patch_n; i++) {
		for (j = 0; j < knot_n_xi[i]; j++) {
			fscanf(fp, "%le", &knot_vec_xi[i][j]);
			printf("%f\t", knot_vec_xi[i][j]);
		}
		printf("\n");
		for (j = 0; j < knot_n_eta[i]; j++) {
			fscanf(fp, "%le", &knot_vec_eta[i][j]);
			printf("%f\t", knot_vec_eta[i][j]);
		}
		printf("\n");
	}
	printf("\n");

	for (i = 0; i < cntl_p_n; i++) {
		fscanf(fp, "%d%le%le%le",
			   &temp_int,
			   &temp_cntl_px[i], &temp_cntl_py[i], &temp_weight[i]);
		printf("%d\t%f\t%f\t%f\n",
			   temp_int,
			   temp_cntl_px[i], temp_cntl_py[i], temp_weight[i]);
	}
	printf("\n");

	for (i = 0; i < patch_n; i++) {
		for (j = 0; j < cntl_p_n_xi[i] * cntl_p_n_eta[i]; j++) {
			cntl_px[i][j] = temp_cntl_px[temp_index[i][j]];
			cntl_py[i][j] = temp_cntl_py[temp_index[i][j]];
			weight[i][j] = temp_weight[temp_index[i][j]];
			printf("%d\t%f\t%f\t%f\n",
				   temp_index[i][j],
				   cntl_px[i][j], cntl_py[i][j], weight[i][j]);
		}
		printf("\n");
	}
	fclose(fp);
	printf("End Reading input\n\n");

	if (fields_flag) {
		printf("Start Reading displacement\n\n");
		fp = fopen("Displacement.dat", "r");
		char buff[256];

		fscanf(fp, "%s", buff);
		fscanf(fp, "%s", buff);

		for (i = 0; i < cntl_p_n; i++) {
			fscanf(fp, "%d:%le%le",
				   &temp_int, &temp_disp_x[i], &temp_disp_y[i]);
			printf("%d\t%1.6e\t%1.6e\n",
				   temp_int, temp_disp_x[i], temp_disp_y[i]);
		}
		printf("\n");

		for (i = 0; i < patch_n; i++) {
			for (j = 0; j < cntl_p_n_xi[i] * cntl_p_n_eta[i]; j++) {
				disp_cntl_px[i][j] = temp_disp_x[temp_index[i][j]];
				disp_cntl_py[i][j] = temp_disp_y[temp_index[i][j]];
				printf("%d\t%f\t%f\t%f\n",
					   temp_index[i][j], cntl_px[i][j], cntl_py[i][j], weight[i][j]);
			}
			printf("\n");
		}
		fclose(fp); // ファイルを閉じる
		printf("End Reading displpacement\n\n");
	}

	fp = fopen("Displacement_loc.dat", "w");
	glo_cntl_p_n = cntl_p_n - loc_cntl_p_n;
	fprintf(fp, "label=Displacement\n"
				"num_items=%d\n\n", loc_cntl_p_n);
	for (i = 0; i < loc_cntl_p_n; i++)
	{
		// fprintf(fp, "%d:	%le %le \n", i, temp_disp_x[i + glo_cntl_p_n], temp_disp_y[i + glo_cntl_p_n]);
		fprintf(fp, "%d:	%.16e %.16e \n", i, temp_disp_x[i + glo_cntl_p_n], temp_disp_y[i + glo_cntl_p_n]);
	}
	fclose(fp);
}


int CalcXiEtaByNR(double px, double py,
				  double *input_knot_vec_xi, double *input_knot_vec_eta,
				  double *cntl_px, double *cntl_py,
				  double *disp_cntl_px, double *disp_cntl_py,
				  int cntl_p_n_xi, int cntl_p_n_eta,
				  double *weight, int order_xi, int order_eta,
				  double *output_xi, double *output_eta,
				  double *disp_x_glo, double *disp_y_glo,
				  double *strain_xx_glo, double *strain_yy_glo, double *strain_xy_glo) {
	double temp_xi, temp_eta;
	double temp_x, temp_y;
	double temp_matrix[2][2];
	double temp_dxi, temp_deta;
	double temp_tol_x = DBL_MAX;
	double temp_tol_y = DBL_MAX;

	(*output_xi) = 0;
	(*output_eta) = 0;

	int i;
	// double tol = 10e-8;
	// int repeat = 10000;
	int repeat = 100;
	double tol = 10e-14;

	temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;
	//printf("% 1.8e % 1.8e\n", temp_xi, temp_eta);

	for (i = 0; i < repeat; i++) {
		rNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
		               cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
		               weight, order_xi, order_eta,
		               temp_xi, temp_eta,
		               &temp_x, &temp_y,
		               &temp_matrix[0][0], &temp_matrix[0][1],
		               &temp_matrix[1][0], &temp_matrix[1][1]);

		temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;

		//収束した場合////////////////////////////////////////////////////////////////
		//if (temp_tol_x < tol && temp_tol_y < tol) {
        if (temp_tol_x + temp_tol_y < tol) {
			// printf("rNURBS\n");
			// printf("repeat = %d\n", i);
			if (temp_xi == input_knot_vec_xi[0] || temp_eta == input_knot_vec_eta[0])
			{
				break;
			}
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;

			int knot_n_xi = cntl_p_n_xi + order_xi + 1;
			int knot_n_eta = cntl_p_n_eta + order_eta + 1;
			double dtilda_xi = 0.0;
			double dtilda_eta = 0.0;
			double disp_x, disp_y;
			double dxi_x, deta_x, dxi_y, deta_y;
			double dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y;
			double temp_matrix2[2][2];
			double temp;
			double strain_xx, strain_yy, strain_xy;
			// double stress_xx, stress_yy, stress_xy;

			for (i = 0; i < knot_n_xi; i++) {
				if ( input_knot_vec_xi[i] < temp_xi && temp_xi <= input_knot_vec_xi[i + 1]) {
					dtilda_xi = ( input_knot_vec_xi[i + 1] - input_knot_vec_xi[i] ) / 2.0;
					// printf("xi%f\n", dtilda_xi);
					break;
				}
			}
			for (i = 0; i < knot_n_eta; i++) {
				if ( input_knot_vec_eta[i] < temp_eta && temp_eta <= input_knot_vec_eta[i + 1]) {
					dtilda_eta = ( input_knot_vec_eta[i + 1] - input_knot_vec_eta[i] ) / 2.0;
					// printf("eta%f\n", dtilda_eta);
					break;
				}
			}

			rNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
			               cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
			               weight, order_xi, order_eta,
			               temp_xi, temp_eta,
			               &temp, &temp,
			               &dxi_x, &deta_x,
			               &dxi_y, &deta_y);
			// printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
			// 	   dxi_x, deta_x, dxi_y, deta_y);

			rNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
			               disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
			               weight, order_xi, order_eta,
			               temp_xi, temp_eta,
			               &disp_x, &disp_y,
			               &dxi_disp_x, &deta_disp_x,
			               &dxi_disp_y, &deta_disp_y);
			// printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
			// 	   dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y);

			temp_matrix2[0][0] = dxi_x * dtilda_xi;
			temp_matrix2[0][1] = dxi_y * dtilda_xi;
			temp_matrix2[1][0] = deta_x * dtilda_eta;
			temp_matrix2[1][1] = deta_y * dtilda_eta;

			InverseMatrix_2D(temp_matrix2);

			strain_xx = temp_matrix2[0][0] * dtilda_xi
			            * dxi_disp_x
			            + temp_matrix2[0][1] * dtilda_eta
			            * deta_disp_x;
			strain_yy = temp_matrix2[1][0] * dtilda_xi
			            * dxi_disp_y
			            + temp_matrix2[1][1] * dtilda_eta
			            * deta_disp_y;
			strain_xy = temp_matrix2[1][0] * dtilda_xi
			            * dxi_disp_x
			            + temp_matrix2[1][1] * dtilda_eta
			            * deta_disp_x
			            + temp_matrix2[0][0] * dtilda_xi
			            * dxi_disp_y
			            + temp_matrix2[0][1] * dtilda_eta
			            * deta_disp_y;

			// 応力の計算を行わないのでコメントアウト
			// double D_matrix[3][3] = {{0.0}};
			// int DM = 1;
			// if (DM == 0) { //平面応力状態
			// 	temp = E * (1.0 - nu * nu);
			// 	D_matrix[0][0] = temp;
			// 	D_matrix[0][1] = nu * temp;
			// 	D_matrix[1][0] = nu * temp;
			// 	D_matrix[1][1] = temp;
			// 	D_matrix[2][2] = (1.0 - nu) / 2.0 * temp;
			// } else if (DM == 1) { //平面ひずみ状態(2Dの場合はこっち)
			// 	temp = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
			// 	D_matrix[0][0] = temp;
			// 	D_matrix[0][1] = nu / (1.0 - nu) * temp;
			// 	D_matrix[1][0] = nu / (1.0 - nu) * temp;
			// 	D_matrix[1][1] = temp;
			// 	D_matrix[2][2] = (1.0 - 2.0 * nu) / 2.0 / (1.0 - nu) * temp;
			// }

			// strainしか使わないのでコメントアウト
			// stress_xx = D_matrix[0][0] * strain_xx + D_matrix[0][1] * strain_yy;
			// stress_yy = D_matrix[1][0] * strain_xx + D_matrix[1][1] * strain_yy;
			// stress_xy = D_matrix[2][2] * strain_xy;

			temp = sqrt(disp_x * disp_x + disp_y * disp_y);
			// printf("x:   % 1.8e\n", px);
			// printf("y:   % 1.8e\n", py);
			// printf("xi:  % 1.8e\n", temp_xi);
			// printf("eta: % 1.8e\n", temp_eta);
			// printf("Displacement x: % 1.8e\n", disp_x);
			// printf("Displacement y: % 1.8e\n", disp_y);
			// printf("Displacement  : % 1.8e\n", temp);
			// printf("Strain xx: % 1.8e\n", strain_xx);
			// printf("Strain yy: % 1.8e\n", strain_yy);
			// printf("Strain xy: % 1.8e\n", strain_xy);
			// printf("Stress xx: % 1.8e\n", stress_xx);
			// printf("Stress yy: % 1.8e\n", stress_yy);
			// printf("Stress xy: % 1.8e\n", stress_xy);

			(*disp_x_glo) = disp_x;
			(*disp_y_glo) = disp_y;

			(*strain_xx_glo) = strain_xx;
			(*strain_yy_glo) = strain_yy;
			(*strain_xy_glo) = strain_xy;

			return i;
		}

		InverseMatrix_2D(temp_matrix);

		temp_dxi = temp_matrix[0][0] * (px - temp_x)
		           + temp_matrix[0][1] * (py - temp_y);
		temp_deta = temp_matrix[1][0] * (px - temp_x)
		            + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_deta;
		if (temp_xi < input_knot_vec_xi[0])
			temp_xi = input_knot_vec_xi[0];
		if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
			temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
		if (temp_eta < input_knot_vec_eta[0])
			temp_eta = input_knot_vec_eta[0];
		if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
			temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];

		//temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
		//printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
	}
	temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;
	for (i = 0; i < repeat; i++) {
		lNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
		               cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
		               weight, order_xi, order_eta,
		               temp_xi, temp_eta,
		               &temp_x, &temp_y,
		               &temp_matrix[0][0], &temp_matrix[0][1],
		               &temp_matrix[1][0], &temp_matrix[1][1]);

		temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;

		//収束した場合////////////////////////////////////////////////////////////////
		//if (temp_tol_x < tol && temp_tol_y < tol) {
        if (temp_tol_x + temp_tol_y < tol) {
			// printf("lNURBS\n");
			// printf("repeat = %d\n", i);
			if (temp_xi == input_knot_vec_xi[cntl_p_n_xi + order_xi] || temp_eta == input_knot_vec_eta[cntl_p_n_eta + order_eta])
			{
				break;
			}
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;

			int knot_n_xi = cntl_p_n_xi + order_xi + 1;
			int knot_n_eta = cntl_p_n_eta + order_eta + 1;
			double dtilda_xi = 0.0;
			double dtilda_eta = 0.0;
			double disp_x, disp_y;
			double dxi_x, deta_x, dxi_y, deta_y;
			double dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y;
			double temp_matrix2[2][2];
			double temp;
			double strain_xx, strain_yy, strain_xy;
			// double stress_xx, stress_yy, stress_xy;

			for (i = 0; i < knot_n_xi; i++) {
				if ( input_knot_vec_xi[i] <= temp_xi && temp_xi < input_knot_vec_xi[i + 1]) {
					dtilda_xi = ( input_knot_vec_xi[i + 1] - input_knot_vec_xi[i] ) / 2.0;
					// printf("%f\n", dtilda_xi);
					break;
				}
			}
			for (i = 0; i < knot_n_eta; i++) {
				if ( input_knot_vec_eta[i] <= temp_eta && temp_eta < input_knot_vec_eta[i + 1]) {
					dtilda_eta = ( input_knot_vec_eta[i + 1] - input_knot_vec_eta[i] ) / 2.0;
					// printf("%f\n", dtilda_eta);
					break;
				}
			}

			lNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
			               cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
			               weight, order_xi, order_eta,
			               temp_xi, temp_eta,
			               &temp, &temp,
			               &dxi_x, &deta_x,
			               &dxi_y, &deta_y);
			// printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
			// 	   dxi_x, deta_x, dxi_y, deta_y);

			lNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
			               disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
			               weight, order_xi, order_eta,
			               temp_xi, temp_eta,
			               &disp_x, &disp_y,
			               &dxi_disp_x, &deta_disp_x,
			               &dxi_disp_y, &deta_disp_y);
			// printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
			// 	   dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y);

			temp_matrix2[0][0] = dxi_x * dtilda_xi;
			temp_matrix2[0][1] = dxi_y * dtilda_xi;
			temp_matrix2[1][0] = deta_x * dtilda_eta;
			temp_matrix2[1][1] = deta_y * dtilda_eta;

			InverseMatrix_2D(temp_matrix2);

			strain_xx = temp_matrix2[0][0] * dtilda_xi
			            * dxi_disp_x
			            + temp_matrix2[0][1] * dtilda_eta
			            * deta_disp_x;
			strain_yy = temp_matrix2[1][0] * dtilda_xi
			            * dxi_disp_y
			            + temp_matrix2[1][1] * dtilda_eta
			            * deta_disp_y;
			strain_xy = temp_matrix2[1][0] * dtilda_xi
			            * dxi_disp_x
			            + temp_matrix2[1][1] * dtilda_eta
			            * deta_disp_x
			            + temp_matrix2[0][0] * dtilda_xi
			            * dxi_disp_y
			            + temp_matrix2[0][1] * dtilda_eta
			            * deta_disp_y;

			// 応力の計算を行わないのでコメントアウト
			// double D_matrix[3][3] = {{0.0}};
			// int DM = 1;
			// if (DM == 0) { //平面応力状態
			// 	temp = E * (1.0 - nu * nu);
			// 	D_matrix[0][0] = temp;
			// 	D_matrix[0][1] = nu * temp;
			// 	D_matrix[1][0] = nu * temp;
			// 	D_matrix[1][1] = temp;
			// 	D_matrix[2][2] = (1.0 - nu) / 2.0 * temp;
			// } else if (DM == 1) { //平面ひずみ状態(2Dの場合はこっち)
			// 	temp = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
			// 	D_matrix[0][0] = temp;
			// 	D_matrix[0][1] = nu / (1.0 - nu) * temp;
			// 	D_matrix[1][0] = nu / (1.0 - nu) * temp;
			// 	D_matrix[1][1] = temp;
			// 	D_matrix[2][2] = (1.0 - 2.0 * nu) / 2.0 / (1.0 - nu) * temp;
			// }

			// strainしか使わないのでコメントアウト
			// stress_xx = D_matrix[0][0] * strain_xx + D_matrix[0][1] * strain_yy;
			// stress_yy = D_matrix[1][0] * strain_xx + D_matrix[1][1] * strain_yy;
			// stress_xy = D_matrix[2][2] * strain_xy;

			temp = sqrt(disp_x * disp_x + disp_y * disp_y);
			// printf("x:   % 1.8e\n", px);
			// printf("y:   % 1.8e\n", py);
			// printf("xi:  % 1.8e\n", temp_xi);
			// printf("eta: % 1.8e\n", temp_eta);
			// printf("Displacement x: % 1.8e\n", disp_x);
			// printf("Displacement y: % 1.8e\n", disp_y);
			// printf("Displacement  : % 1.8e\n", temp);
			// printf("Strain xx: % 1.8e\n", strain_xx);
			// printf("Strain yy: % 1.8e\n", strain_yy);
			// printf("Strain xy: % 1.8e\n", strain_xy);
			// printf("Stress xx: % 1.8e\n", stress_xx);
			// printf("Stress yy: % 1.8e\n", stress_yy);
			// printf("Stress xy: % 1.8e\n", stress_xy);

			(*disp_x_glo) = disp_x;
			(*disp_y_glo) = disp_y;

			(*strain_xx_glo) = strain_xx;
			(*strain_yy_glo) = strain_yy;
			(*strain_xy_glo) = strain_xy;

			return i;
		}

		InverseMatrix_2D(temp_matrix);

		temp_dxi = temp_matrix[0][0] * (px - temp_x)
		           + temp_matrix[0][1] * (py - temp_y);
		temp_deta = temp_matrix[1][0] * (px - temp_x)
		            + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_deta;
		if (temp_xi < input_knot_vec_xi[0])
			temp_xi = input_knot_vec_xi[0];
		if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
			temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
		if (temp_eta < input_knot_vec_eta[0])
			temp_eta = input_knot_vec_eta[0];
		if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
			temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];

		//temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
		//printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
	}
	temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;
	for (i = 0; i < repeat; i++) {
		rlNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
		               cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
		               weight, order_xi, order_eta,
		               temp_xi, temp_eta,
		               &temp_x, &temp_y,
		               &temp_matrix[0][0], &temp_matrix[0][1],
		               &temp_matrix[1][0], &temp_matrix[1][1]);

		temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;

		//収束した場合////////////////////////////////////////////////////////////////
		//if (temp_tol_x < tol && temp_tol_y < tol) {
        if (temp_tol_x + temp_tol_y < tol) {
			// printf("rlNURBS\n");
			if (temp_xi == input_knot_vec_xi[0])
			{
				break;
			}
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;

			int knot_n_xi = cntl_p_n_xi + order_xi + 1;
			int knot_n_eta = cntl_p_n_eta + order_eta + 1;
			double dtilda_xi = 0.0;
			double dtilda_eta = 0.0;
			double disp_x, disp_y;
			double dxi_x, deta_x, dxi_y, deta_y;
			double dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y;
			double temp_matrix2[2][2];
			double temp;
			double strain_xx, strain_yy, strain_xy;
			// double stress_xx, stress_yy, stress_xy;

			for (i = 0; i < knot_n_xi; i++) {
				if ( input_knot_vec_xi[i] < temp_xi && temp_xi <= input_knot_vec_xi[i + 1]) {
					dtilda_xi = ( input_knot_vec_xi[i + 1] - input_knot_vec_xi[i] ) / 2.0;
					//printf("%f\n", dtilda_xi);
					break;
				}
			}
			for (i = 0; i < knot_n_eta; i++) {
				if ( input_knot_vec_eta[i] <= temp_eta && temp_eta < input_knot_vec_eta[i + 1]) {
					dtilda_eta = ( input_knot_vec_eta[i + 1] - input_knot_vec_eta[i] ) / 2.0;
					//printf("%f\n", dtilda_eta);
					break;
				}
			}

			rlNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
			               cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
			               weight, order_xi, order_eta,
			               temp_xi, temp_eta,
			               &temp, &temp,
			               &dxi_x, &deta_x,
			               &dxi_y, &deta_y);
			//printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
			//	   dxi_x, deta_x, dxi_y, deta_y);

			rlNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
			               disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
			               weight, order_xi, order_eta,
			               temp_xi, temp_eta,
			               &disp_x, &disp_y,
			               &dxi_disp_x, &deta_disp_x,
			               &dxi_disp_y, &deta_disp_y);
			//printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
			//	   dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y);

			temp_matrix2[0][0] = dxi_x * dtilda_xi;
			temp_matrix2[0][1] = dxi_y * dtilda_xi;
			temp_matrix2[1][0] = deta_x * dtilda_eta;
			temp_matrix2[1][1] = deta_y * dtilda_eta;

			InverseMatrix_2D(temp_matrix2);

			strain_xx = temp_matrix2[0][0] * dtilda_xi
			            * dxi_disp_x
			            + temp_matrix2[0][1] * dtilda_eta
			            * deta_disp_x;
			strain_yy = temp_matrix2[1][0] * dtilda_xi
			            * dxi_disp_y
			            + temp_matrix2[1][1] * dtilda_eta
			            * deta_disp_y;
			strain_xy = temp_matrix2[1][0] * dtilda_xi
			            * dxi_disp_x
			            + temp_matrix2[1][1] * dtilda_eta
			            * deta_disp_x
			            + temp_matrix2[0][0] * dtilda_xi
			            * dxi_disp_y
			            + temp_matrix2[0][1] * dtilda_eta
			            * deta_disp_y;

			// 応力の計算を行わないのでコメントアウト
			// double D_matrix[3][3] = {{0.0}};
			// int DM = 1;
			// if (DM == 0) { //平面応力状態
			// 	temp = E * (1.0 - nu * nu);
			// 	D_matrix[0][0] = temp;
			// 	D_matrix[0][1] = nu * temp;
			// 	D_matrix[1][0] = nu * temp;
			// 	D_matrix[1][1] = temp;
			// 	D_matrix[2][2] = (1.0 - nu) / 2.0 * temp;
			// } else if (DM == 1) { //平面ひずみ状態(2Dの場合はこっち)
			// 	temp = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
			// 	D_matrix[0][0] = temp;
			// 	D_matrix[0][1] = nu / (1.0 - nu) * temp;
			// 	D_matrix[1][0] = nu / (1.0 - nu) * temp;
			// 	D_matrix[1][1] = temp;
			// 	D_matrix[2][2] = (1.0 - 2.0 * nu) / 2.0 / (1.0 - nu) * temp;
			// }

			// strainしか使わないのでコメントアウト
			// stress_xx = D_matrix[0][0] * strain_xx + D_matrix[0][1] * strain_yy;
			// stress_yy = D_matrix[1][0] * strain_xx + D_matrix[1][1] * strain_yy;
			// stress_xy = D_matrix[2][2] * strain_xy;

			temp = sqrt(disp_x * disp_x + disp_y * disp_y);
			// printf("x:   % 1.8e\n", px);
			// printf("y:   % 1.8e\n", py);
			// printf("xi:  % 1.8e\n", temp_xi);
			// printf("eta: % 1.8e\n", temp_eta);
			// printf("Displacement x: % 1.8e\n", disp_x);
			// printf("Displacement y: % 1.8e\n", disp_y);
			// printf("Displacement  : % 1.8e\n", temp);
			// printf("Strain xx: % 1.8e\n", strain_xx);
			// printf("Strain yy: % 1.8e\n", strain_yy);
			// printf("Strain xy: % 1.8e\n", strain_xy);
			// printf("Stress xx: % 1.8e\n", stress_xx);
			// printf("Stress yy: % 1.8e\n", stress_yy);
			// printf("Stress xy: % 1.8e\n", stress_xy);

			(*disp_x_glo) = disp_x;
			(*disp_y_glo) = disp_y;

			(*strain_xx_glo) = strain_xx;
			(*strain_yy_glo) = strain_yy;
			(*strain_xy_glo) = strain_xy;

			return i;
		}

		InverseMatrix_2D(temp_matrix);

		temp_dxi = temp_matrix[0][0] * (px - temp_x)
		           + temp_matrix[0][1] * (py - temp_y);
		temp_deta = temp_matrix[1][0] * (px - temp_x)
		            + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_deta;
		if (temp_xi < input_knot_vec_xi[0])
			temp_xi = input_knot_vec_xi[0];
		if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
			temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
		if (temp_eta < input_knot_vec_eta[0])
			temp_eta = input_knot_vec_eta[0];
		if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
			temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];

		//temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
		//printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
	}
	temp_xi = input_knot_vec_xi[0] + input_knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = input_knot_vec_eta[0] + input_knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;
	for (i = 0; i < repeat; i++) {
		lrNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
		               cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
		               weight, order_xi, order_eta,
		               temp_xi, temp_eta,
		               &temp_x, &temp_y,
		               &temp_matrix[0][0], &temp_matrix[0][1],
		               &temp_matrix[1][0], &temp_matrix[1][1]);

		temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;

		//収束した場合////////////////////////////////////////////////////////////////
		//if (temp_tol_x < tol && temp_tol_y < tol) {
        if (temp_tol_x + temp_tol_y < tol) {
			// printf("lrNURBS\n");
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;

			int knot_n_xi = cntl_p_n_xi + order_xi + 1;
			int knot_n_eta = cntl_p_n_eta + order_eta + 1;
			double dtilda_xi = 0.0;
			double dtilda_eta = 0.0;
			double disp_x, disp_y;
			double dxi_x, deta_x, dxi_y, deta_y;
			double dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y;
			double temp_matrix2[2][2];
			double temp;
			double strain_xx, strain_yy, strain_xy;
			// double stress_xx, stress_yy, stress_xy;

			for (i = 0; i < knot_n_xi; i++) {
				if ( input_knot_vec_xi[i] <= temp_xi && temp_xi < input_knot_vec_xi[i + 1]) {
					dtilda_xi = ( input_knot_vec_xi[i + 1] - input_knot_vec_xi[i] ) / 2.0;
					//printf("%f\n", dtilda_xi);
					break;
				}
			}
			for (i = 0; i < knot_n_eta; i++) {
				if ( input_knot_vec_eta[i] < temp_eta && temp_eta <= input_knot_vec_eta[i + 1]) {
					dtilda_eta = ( input_knot_vec_eta[i + 1] - input_knot_vec_eta[i] ) / 2.0;
					//printf("%f\n", dtilda_eta);
					break;
				}
			}

			lrNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
			               cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
			               weight, order_xi, order_eta,
			               temp_xi, temp_eta,
			               &temp, &temp,
			               &dxi_x, &deta_x,
			               &dxi_y, &deta_y);
			//printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
			//	   dxi_x, deta_x, dxi_y, deta_y);

			lrNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
			               disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
			               weight, order_xi, order_eta,
			               temp_xi, temp_eta,
			               &disp_x, &disp_y,
			               &dxi_disp_x, &deta_disp_x,
			               &dxi_disp_y, &deta_disp_y);
			//printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
			//	   dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y);

			temp_matrix2[0][0] = dxi_x * dtilda_xi;
			temp_matrix2[0][1] = dxi_y * dtilda_xi;
			temp_matrix2[1][0] = deta_x * dtilda_eta;
			temp_matrix2[1][1] = deta_y * dtilda_eta;

			InverseMatrix_2D(temp_matrix2);

			strain_xx = temp_matrix2[0][0] * dtilda_xi
			            * dxi_disp_x
			            + temp_matrix2[0][1] * dtilda_eta
			            * deta_disp_x;
			strain_yy = temp_matrix2[1][0] * dtilda_xi
			            * dxi_disp_y
			            + temp_matrix2[1][1] * dtilda_eta
			            * deta_disp_y;
			strain_xy = temp_matrix2[1][0] * dtilda_xi
			            * dxi_disp_x
			            + temp_matrix2[1][1] * dtilda_eta
			            * deta_disp_x
			            + temp_matrix2[0][0] * dtilda_xi
			            * dxi_disp_y
			            + temp_matrix2[0][1] * dtilda_eta
			            * deta_disp_y;

			// 応力の計算を行わないのでコメントアウト
			// double D_matrix[3][3] = {{0.0}};
			// int DM = 1;
			// if (DM == 0) { //平面応力状態
			// 	temp = E * (1.0 - nu * nu);
			// 	D_matrix[0][0] = temp;
			// 	D_matrix[0][1] = nu * temp;
			// 	D_matrix[1][0] = nu * temp;
			// 	D_matrix[1][1] = temp;
			// 	D_matrix[2][2] = (1.0 - nu) / 2.0 * temp;
			// } else if (DM == 1) { //平面ひずみ状態(2Dの場合はこっち)
			// 	temp = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
			// 	D_matrix[0][0] = temp;
			// 	D_matrix[0][1] = nu / (1.0 - nu) * temp;
			// 	D_matrix[1][0] = nu / (1.0 - nu) * temp;
			// 	D_matrix[1][1] = temp;
			// 	D_matrix[2][2] = (1.0 - 2.0 * nu) / 2.0 / (1.0 - nu) * temp;
			// }

			// strainしか使わないのでコメントアウト
			// stress_xx = D_matrix[0][0] * strain_xx + D_matrix[0][1] * strain_yy;
			// stress_yy = D_matrix[1][0] * strain_xx + D_matrix[1][1] * strain_yy;
			// stress_xy = D_matrix[2][2] * strain_xy;

			temp = sqrt(disp_x * disp_x + disp_y * disp_y);
			// printf("x:   % 1.8e\n", px);
			// printf("y:   % 1.8e\n", py);
			// printf("xi:  % 1.8e\n", temp_xi);
			// printf("eta: % 1.8e\n", temp_eta);
			// printf("Displacement x: % 1.8e\n", disp_x);
			// printf("Displacement y: % 1.8e\n", disp_y);
			// printf("Displacement  : % 1.8e\n", temp);
			// printf("Strain xx: % 1.8e\n", strain_xx);
			// printf("Strain yy: % 1.8e\n", strain_yy);
			// printf("Strain xy: % 1.8e\n", strain_xy);
			// printf("Stress xx: % 1.8e\n", stress_xx);
			// printf("Stress yy: % 1.8e\n", stress_yy);
			// printf("Stress xy: % 1.8e\n", stress_xy);

			(*disp_x_glo) = disp_x;
			(*disp_y_glo) = disp_y;

			(*strain_xx_glo) = strain_xx;
			(*strain_yy_glo) = strain_yy;
			(*strain_xy_glo) = strain_xy;

			return i;
		}

		InverseMatrix_2D(temp_matrix);

		temp_dxi = temp_matrix[0][0] * (px - temp_x)
		           + temp_matrix[0][1] * (py - temp_y);
		temp_deta = temp_matrix[1][0] * (px - temp_x)
		            + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_deta;
		if (temp_xi < input_knot_vec_xi[0])
			temp_xi = input_knot_vec_xi[0];
		if (temp_xi > input_knot_vec_xi[cntl_p_n_xi + order_xi])
			temp_xi = input_knot_vec_xi[cntl_p_n_xi + order_xi];
		if (temp_eta < input_knot_vec_eta[0])
			temp_eta = input_knot_vec_eta[0];
		if (temp_eta > input_knot_vec_eta[cntl_p_n_eta + order_eta])
			temp_eta = input_knot_vec_eta[cntl_p_n_eta + order_eta];

		//temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
		//printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
	}
	//printf("% 1.8e % 1.8e\n", temp_x, temp_y);
	return 0;
}


static void Calculation(int order_xi, int order_eta,
						int knot_n_xi, int knot_n_eta,
						int cntl_p_n_xi, int cntl_p_n_eta,
						double *input_knot_vec_xi, double *input_knot_vec_eta,
						double *cntl_px, double *cntl_py,
						double *disp_cntl_px, double *disp_cntl_py,
						double *weight) {
	int i, j, k, l;
	double temp1, temp2, temp3;
	double temp_matrix[2][2];

	//計算するξ,ηの値決定と ∂ξ/∂チルダξ, ∂η/∂チルダη の計算
	double calc_xi[MAX_POINTS];		//計算するξの値
	double calc_eta[MAX_POINTS];		//計算するηの値
	double dtilda_xi[MAX_KNOTS];		// ∂ξ/∂チルダξ
	double dtilda_eta[MAX_KNOTS];	// ∂η/∂チルダη
	k = 0;
	l = 0;
	for (i = 0; i < knot_n_xi - 1; i++) {
		if ( input_knot_vec_xi[i] != input_knot_vec_xi[i + 1] ) {
			calc_xi[k] = input_knot_vec_xi[i];
			printf("%d\t%f\n", k, calc_xi[k]);
			dtilda_xi[l] = ( input_knot_vec_xi[i + 1] - input_knot_vec_xi[i] ) / 2.0;
			printf("%d\t%f\n", k, dtilda_xi[k]);
			k++;
			l++;
			if (division_ele_xi > 1) {
				temp1 = (input_knot_vec_xi[i + 1] - input_knot_vec_xi[i])
						/ (double)division_ele_xi;
				for (j = 1; j < division_ele_xi; j++) {
					calc_xi[k] = calc_xi[k - 1] + temp1;
					printf("%d\t%f\n", k, calc_xi[k]);
					k++;
				}
			}
		}
	}
	calc_xi[k] = input_knot_vec_xi[knot_n_xi - 1];
	printf("%d\t%f\n", k, calc_xi[k]);
	//printf("\n");
	division_n_xi = k + 1;
	element_n_xi = l;

	k = 0;
	l = 0;
	for (i = 0; i < knot_n_eta - 1; i++) {
		if ( input_knot_vec_eta[i] != input_knot_vec_eta[i + 1] ) {
			calc_eta[k] = input_knot_vec_eta[i];
			//printf("%d\t%f\n", k, calc_eta[k]);
			dtilda_eta[l] = ( input_knot_vec_eta[i + 1] - input_knot_vec_eta[i] ) / 2.0;
			//printf("%d\t%f\n", k, dtilda_eta[k]);
			k++;
			l++;
			if (division_ele_eta > 1) {
				temp1 = (input_knot_vec_eta[i + 1] - input_knot_vec_eta[i])
						/ (double)division_ele_eta;
				for (j = 1; j < division_ele_eta; j++) {
					calc_eta[k] = calc_eta[k - 1] + temp1;
					//printf("%d\t%f\n", k, calc_eta[k]);
					k++;
				}
			}
		}
	}
	calc_eta[k] = input_knot_vec_eta[knot_n_eta - 1];
	//printf("%d\t%f\n", k, calc_eta[k]);
	//printf("\n");
	division_n_eta = k + 1;
	element_n_eta = l;

	if (element_n_xi > MAX_ELEMENTS) {
		printf("Error!!\n");
		printf("Too many elements at xi!\n"
			   "Maximum of elements is %d (Now %d)\n"
			   "\n", MAX_ELEMENTS, element_n_xi);
		exit(1);
	}
	if (element_n_eta > MAX_ELEMENTS) {
		printf("Error!!\n");
		printf("Too many elements at eta!\n"
			   "Maximum of elements is %d (Now %d)\n"
			   "\n", MAX_ELEMENTS, element_n_eta);
		exit(1);
	}

	int ii, jj, kk, ll;

	//メッシュ座標計算
	printf("Start Calculation mesh\n\n");
	for (i = 0; i < division_n_xi; i++) {
		ii = i / division_ele_xi;
		kk = i % division_ele_xi;
		for (j = 0; j < division_n_eta; j++) {
			jj = j / division_ele_eta;
			ll = j % division_ele_eta;
			lNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
						   cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
						   weight, order_xi, order_eta,
						   calc_xi[i], calc_eta[j],
						   &coord_x[i][j], &coord_y[i][j],
						   &dxi_x[ii][jj][kk][ll], &deta_x[ii][jj][kk][ll],
						   &dxi_y[ii][jj][kk][ll], &deta_y[ii][jj][kk][ll]);
			// printf("[%d][%d] [%d][%d][%d][%d]"
			// 	   "% 1.4e % 1.4e "
			// 	   "% 1.4e % 1.4e\n",
			// 	   i, j, ii, jj, kk, ll,
			// 	   calc_xi[i], calc_eta[j],
			// 	   coord_x[i][j], coord_y[i][j]);
		}
		//printf("\n");
	}
	printf("\n");
	printf("End Calculation mesh\n\n");

	if (fields_flag) {
		//変位計算
		printf("Start Calculation displpacement\n\n");
		for (i = 0; i < division_n_xi; i++) {
			ii = i / division_ele_xi;
			kk = i % division_ele_xi;
			for (j = 0; j < division_n_eta; j++) {
				jj = j / division_ele_eta;
				ll = j % division_ele_eta;
				lNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
							   disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
							   weight, order_xi, order_eta,
							   calc_xi[i], calc_eta[j],
							   &disp_x[i][j], &disp_y[i][j],
							   &dxi_disp_x[ii][jj][kk][ll], &deta_disp_x[ii][jj][kk][ll],
							   &dxi_disp_y[ii][jj][kk][ll], &deta_disp_y[ii][jj][kk][ll]);
				// printf("[%d][%d] [%d][%d][%d][%d]"
				// 	   "% 1.4e % 1.4e "
				// 	   "% 1.4e % 1.4e\n",
				// 	   i, j, ii, jj, kk, ll,
				// 	   calc_xi[i], calc_eta[j],
				// 	   disp_x[i][j], disp_y[i][j]);
			}
			//printf("\n");
		}
		printf("\n");
		printf("End Calculation displpacement\n\n");

		//足りない微分値計算
		for (ii = 0; ii < element_n_xi; ii++) {
			for (jj = 0; jj < element_n_eta; jj++) {
				kk = division_ele_xi;
				i = (ii + 1) * division_ele_xi;
				j = jj * division_ele_eta;
				for (ll = 1; ll < division_ele_eta; ll++) {
					j++;
					rNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
								   cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
								   weight, order_xi, order_eta,
								   calc_xi[i], calc_eta[j],
								   &coord_x[i][j], &coord_y[i][j],
								   &dxi_x[ii][jj][kk][ll], &deta_x[ii][jj][kk][ll],
								   &dxi_y[ii][jj][kk][ll], &deta_y[ii][jj][kk][ll]);
					rNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
								   disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
								   weight, order_xi, order_eta,
								   calc_xi[i], calc_eta[j],
								   &disp_x[i][j], &disp_y[i][j],
								   &dxi_disp_x[ii][jj][kk][ll], &deta_disp_x[ii][jj][kk][ll],
								   &dxi_disp_y[ii][jj][kk][ll], &deta_disp_y[ii][jj][kk][ll]);
					/*
					printf("[%d][%d] [%d][%d][%d][%d]"
						   "% 1.4e % 1.4e % 1.4e % 1.4e"
						   "% 1.4e % 1.4e % 1.4e % 1.4e\n",
						   i, j, ii, jj, kk, ll,
						   dxi_x[ii][jj][kk][ll], deta_x[ii][jj][kk][ll],
						   dxi_y[ii][jj][kk][ll], deta_y[ii][jj][kk][ll],
						   dxi_disp_x[ii][jj][kk][ll], deta_disp_x[ii][jj][kk][ll],
						   dxi_disp_y[ii][jj][kk][ll], deta_disp_y[ii][jj][kk][ll]);
					*/
				}

				ll = division_ele_eta;
				i = ii * division_ele_xi;
				j = (jj + 1) * division_ele_eta;
				for (kk = 1; kk <= division_ele_xi; kk++) {
					i++;
					rNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
								   cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
								   weight, order_xi, order_eta,
								   calc_xi[i], calc_eta[j],
								   &coord_x[i][j], &coord_y[i][j],
								   &dxi_x[ii][jj][kk][ll], &deta_x[ii][jj][kk][ll],
								   &dxi_y[ii][jj][kk][ll], &deta_y[ii][jj][kk][ll]);
					rNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
								   disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
								   weight, order_xi, order_eta,
								   calc_xi[i], calc_eta[j],
								   &disp_x[i][j], &disp_y[i][j],
								   &dxi_disp_x[ii][jj][kk][ll], &deta_disp_x[ii][jj][kk][ll],
								   &dxi_disp_y[ii][jj][kk][ll], &deta_disp_y[ii][jj][kk][ll]);
					/*
					printf("[%d][%d] [%d][%d][%d][%d]"
						   "% 1.4e % 1.4e % 1.4e % 1.4e"
						   "% 1.4e % 1.4e % 1.4e % 1.4e\n",
						   i, j, ii, jj, kk, ll,
						   dxi_x[ii][jj][kk][ll], deta_x[ii][jj][kk][ll],
						   dxi_y[ii][jj][kk][ll], deta_y[ii][jj][kk][ll],
						   dxi_disp_x[ii][jj][kk][ll], deta_disp_x[ii][jj][kk][ll],
						   dxi_disp_y[ii][jj][kk][ll], deta_disp_y[ii][jj][kk][ll]);
					*/
				}

				kk = division_ele_xi;
				ll = 0;
				i = (ii + 1) * division_ele_xi;
				j = jj * division_ele_eta;
				rlNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
								cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
								weight, order_xi, order_eta,
								calc_xi[i], calc_eta[j],
								&coord_x[i][j], &coord_y[i][j],
								&dxi_x[ii][jj][kk][ll], &deta_x[ii][jj][kk][ll],
								&dxi_y[ii][jj][kk][ll], &deta_y[ii][jj][kk][ll]);
				rlNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
								disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
								weight, order_xi, order_eta,
								calc_xi[i], calc_eta[j],
								&disp_x[i][j], &disp_y[i][j],
								&dxi_disp_x[ii][jj][kk][ll], &deta_disp_x[ii][jj][kk][ll],
								&dxi_disp_y[ii][jj][kk][ll], &deta_disp_y[ii][jj][kk][ll]);
				/*
					printf("[%d][%d] [%d][%d][%d][%d]"
						   "% 1.4e % 1.4e % 1.4e % 1.4e"
						   "% 1.4e % 1.4e % 1.4e % 1.4e\n",
						   i, j, ii, jj, kk, ll,
						   dxi_x[ii][jj][kk][ll], deta_x[ii][jj][kk][ll],
						   dxi_y[ii][jj][kk][ll], deta_y[ii][jj][kk][ll],
						   dxi_disp_x[ii][jj][kk][ll], deta_disp_x[ii][jj][kk][ll],
						   dxi_disp_y[ii][jj][kk][ll], deta_disp_y[ii][jj][kk][ll]);
				*/

				kk = 0;
				ll = division_ele_eta;
				i = ii * division_ele_xi;
				j = (jj + 1) * division_ele_eta;
				lrNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
								cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
								weight, order_xi, order_eta,
								calc_xi[i], calc_eta[j],
								&coord_x[i][j], &coord_y[i][j],
								&dxi_x[ii][jj][kk][ll], &deta_x[ii][jj][kk][ll],
								&dxi_y[ii][jj][kk][ll], &deta_y[ii][jj][kk][ll]);
				lrNURBS_surface(input_knot_vec_xi, input_knot_vec_eta,
								disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
								weight, order_xi, order_eta,
								calc_xi[i], calc_eta[j],
								&disp_x[i][j], &disp_y[i][j],
								&dxi_disp_x[ii][jj][kk][ll], &deta_disp_x[ii][jj][kk][ll],
								&dxi_disp_y[ii][jj][kk][ll], &deta_disp_y[ii][jj][kk][ll]);
				/*
					printf("[%d][%d] [%d][%d][%d][%d]"
						   "% 1.4e % 1.4e % 1.4e % 1.4e"
						   "% 1.4e % 1.4e % 1.4e % 1.4e\n",
						   i, j, ii, jj, kk, ll,
						   dxi_x[ii][jj][kk][ll], deta_x[ii][jj][kk][ll],
						   dxi_y[ii][jj][kk][ll], deta_y[ii][jj][kk][ll],
						   dxi_disp_x[ii][jj][kk][ll], deta_disp_x[ii][jj][kk][ll],
						   dxi_disp_y[ii][jj][kk][ll], deta_disp_y[ii][jj][kk][ll]);
				printf("\n");
				*/
			}
		}

		/*
		for (ii = 0; ii < element_n_xi; ii++) {
			for (jj = 0; jj < element_n_eta; jj++) {
				for (kk = 0; kk <= division_ele_xi; kk++) {
					for (ll = 0; ll <= division_ele_eta; ll++) {
						printf("[%d][%d][%d][%d]"
							   "% 1.4e % 1.4e % 1.4e % 1.4e"
							   "% 1.4e % 1.4e % 1.4e % 1.4e\n",
							   ii, jj, kk, ll,
							   dxi_x[ii][jj][kk][ll], deta_x[ii][jj][kk][ll],
							   dxi_y[ii][jj][kk][ll], deta_y[ii][jj][kk][ll],
							   dxi_disp_x[ii][jj][kk][ll], deta_disp_x[ii][jj][kk][ll],
							   dxi_disp_y[ii][jj][kk][ll], deta_disp_y[ii][jj][kk][ll]);
					}
				}
			}
		}
		printf("\n");
		*/

		//ひずみ計算
		printf("Start Calculation Strain\n\n");
		for (i = 0; i < element_n_xi; i++) {
			for (j = 0; j < element_n_eta; j++) {
				temp1 = dtilda_xi[i];
				temp2 = dtilda_eta[j];
				for (k = 0; k < division_ele_xi + 1; k++) {
					for (l = 0; l < division_ele_eta + 1; l++) {
						temp_matrix[0][0] = dxi_x[i][j][k][l] * temp1;
						temp_matrix[0][1] = dxi_y[i][j][k][l] * temp1;
						temp_matrix[1][0] = deta_x[i][j][k][l] * temp2;
						temp_matrix[1][1] = deta_y[i][j][k][l] * temp2;

						InverseMatrix_2D(temp_matrix);

						strain_xx[i][j][k][l] = temp_matrix[0][0] * temp1
												* dxi_disp_x[i][j][k][l]
												+ temp_matrix[0][1] * temp2
												* deta_disp_x[i][j][k][l];
						strain_yy[i][j][k][l] = temp_matrix[1][0] * temp1
												* dxi_disp_y[i][j][k][l]
												+ temp_matrix[1][1] * temp2
												* deta_disp_y[i][j][k][l];
						strain_xy[i][j][k][l] = temp_matrix[1][0] * temp1
												* dxi_disp_x[i][j][k][l]
												+ temp_matrix[1][1] * temp2
												* deta_disp_x[i][j][k][l]
												+ temp_matrix[0][0] * temp1
												* dxi_disp_y[i][j][k][l]
												+ temp_matrix[0][1] * temp2
												* deta_disp_y[i][j][k][l];

						// printf("[%d][%d][%d][%d]\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\n",
						// 	   i, j, k, l, temp1, temp2,
						// 	   strain_xx[i][j][k][l], strain_yy[i][j][k][l], strain_xy[i][j][k][l]);
					}
				}
			}
			// printf("\n");
		}
		printf("End Calculation Strain\n\n");

		//Dマトリクスの計算
		double D_matrix[3][3] = {{0.0}};
		if (DM == 0) { //平面応力状態
			temp1 = E * (1.0 - nu * nu);
			D_matrix[0][0] = temp1;
			D_matrix[0][1] = nu * temp1;
			D_matrix[1][0] = nu * temp1;
			D_matrix[1][1] = temp1;
			D_matrix[2][2] = (1.0 - nu) / 2.0 * temp1;
		} else if (DM == 1) { //平面ひずみ状態(2Dの場合はこっち)
			temp1 = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
			D_matrix[0][0] = temp1;
			D_matrix[0][1] = nu / (1.0 - nu) * temp1;
			D_matrix[1][0] = nu / (1.0 - nu) * temp1;
			D_matrix[1][1] = temp1;
			D_matrix[2][2] = (1.0 - 2.0 * nu) / 2.0 / (1.0 - nu) * temp1;
		}

		printf("Start Calculation Stress\n\n");
		for (i = 0; i < element_n_xi; i++) {
			for (j = 0; j < element_n_eta; j++) {
				for (k = 0; k < division_ele_xi + 1; k++) {
					for (l = 0; l < division_ele_eta + 1; l++) {
						stress_xx[i][j][k][l] = D_matrix[0][0] * strain_xx[i][j][k][l]
												+ D_matrix[0][1] * strain_yy[i][j][k][l];
						stress_yy[i][j][k][l] = D_matrix[1][0] * strain_xx[i][j][k][l]
												+ D_matrix[1][1] * strain_yy[i][j][k][l];
						stress_xy[i][j][k][l] = D_matrix[2][2] * strain_xy[i][j][k][l];
						// printf("[%d][%d][%d][%d]\t% 1.4e\t% 1.4e\t% 1.4e\n", i, j, k, l,
						// 	   stress_xx[i][j][k][l], stress_yy[i][j][k][l], stress_xy[i][j][k][l]);
					}
				}
				// printf("\n");
			}
		}
		printf("End Calculation Stress\n\n");
	}

	//書き込み
	fp = fopen("view.dat", "a");
	if (fields_flag) {
		fprintf(fp, "%d\t%d\t%d\t%d\n",
				division_n_xi, division_n_eta,
				element_n_xi, element_n_eta);
		for (i = 0; i < division_n_xi; i++) {
			for (j = 0; j < division_n_eta; j++) {
				temp1 = disp_x[i][j];
				temp2 = disp_y[i][j];
				temp3 = sqrt(temp1 * temp1 + temp2 * temp2);
				fprintf(fp, "% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\n",
						coord_x[i][j], coord_y[i][j],
						disp_x[i][j], disp_y[i][j], temp3);
			}
		}
		for (i = 0; i < element_n_xi; i++) {
			for (j = 0; j < element_n_eta; j++) {
				for (k = 0; k < division_ele_xi + 1; k++) {
					for (l = 0; l < division_ele_eta + 1; l++) {
						fprintf(fp, "% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\n",
								strain_xx[i][j][k][l], strain_yy[i][j][k][l], strain_xy[i][j][k][l],
								stress_xx[i][j][k][l], stress_yy[i][j][k][l], stress_xy[i][j][k][l]);
					}
				}
			}
		}
	} else {
		fprintf(fp, "%d\t%d\t%d\t%d\n",
				division_n_xi, division_n_eta,
				element_n_xi, element_n_eta);
		for (i = 0; i < division_n_xi; i++) {
			for (j = 0; j < division_n_eta; j++) {
				fprintf(fp, "% 1.4e\t% 1.4e\n",
						coord_x[i][j], coord_y[i][j]);
			}
		}
	}
	fclose(fp);
	//machino
	fp = fopen("view_r_theta.dat", "a");
	if (fields_flag) {
		fprintf(fp, "%d\t%d\t%d\t%d\n",
				division_n_xi, division_n_eta,
				element_n_xi, element_n_eta);
		for (i = 0; i < division_n_xi; i++) {
			for (j = 0; j < division_n_eta; j++) {
				temp1 = disp_x[i][j];
				temp2 = disp_y[i][j];
				temp3 = sqrt(temp1 * temp1 + temp2 * temp2);
				fprintf(fp, "% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\n",
						coord_x[i][j], coord_y[i][j],
						disp_x[i][j], disp_y[i][j], temp3);
			}
		}
		for (i = 0; i < element_n_xi; i++) {
			for (j = 0; j < element_n_eta; j++) {
				for (k = 0; k < division_ele_xi + 1; k++) {
					for (l = 0; l < division_ele_eta + 1; l++) {
					double stress_rr, stress_theta;
					double sum = stress_xx[i][j][k][l] + stress_yy[i][j][k][l];
					double dif = stress_xx[i][j][k][l] - stress_yy[i][j][k][l];
					double tau2 = stress_xy[i][j][k][l] * stress_xy[i][j][k][l];
					stress_rr = sum * 0.5
								+ sqrt(dif * dif + 4 * tau2) * 0.5;
					stress_theta = sum * 0.5
								  - sqrt(dif * dif + 4 * tau2) * 0.5;											
						fprintf(fp, "% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t0.0\n",
								strain_xx[i][j][k][l], strain_yy[i][j][k][l], strain_xy[i][j][k][l],
								stress_rr,stress_theta);
					}
				}
			}
		}
	} else {
		fprintf(fp, "%d\t%d\t%d\t%d\n",
				division_n_xi, division_n_eta,
				element_n_xi, element_n_eta);
		for (i = 0; i < division_n_xi; i++) {
			for (j = 0; j < division_n_eta; j++) {
				fprintf(fp, "% 1.4e\t% 1.4e\n",
						coord_x[i][j], coord_y[i][j]);
			}
		}
	}
	fclose(fp);
	//グラフ用ファイル書き込み
	fp = fopen("disp_graph.txt", "a");
	for (i = 0; i < division_n_xi; i++) {
		for (j = 0; j < division_n_eta; j++) {
			temp1 = disp_x[i][j];
			temp2 = disp_y[i][j];
			//temp3 = sqrt(temp1 * temp1 + temp2 * temp2);
			fprintf(fp, "%d\t% 1.15e\t% 1.15e\t% 1.15e\t% 1.15e\n",
					graph_patch_n,
					coord_x[i][j], coord_y[i][j],
					disp_x[i][j], disp_y[i][j]);
		}
	}
	fclose(fp);

	fp = fopen("stress_y_graph.txt", "a");
 	for (i = 0; i < element_n_xi; i++) {
		for (j = 0; j < element_n_eta; j++) {
			for (k = 0; k < division_ele_xi + 1; k++) {
				for (l = 0; l < division_ele_eta + 1; l++) {
				fprintf(fp, "%d\t% 1.15e\t% 1.15e\t% 1.15e\n",
							graph_patch_n,
							coord_x[i * division_ele_xi + k][j * division_ele_eta + l],
							coord_y[i * division_ele_xi + k][j * division_ele_eta + l],
							stress_yy[i][j][k][l]);
				}
			}
		}
	}
	fclose(fp);

	fp = fopen("stress_y_graph_0.txt", "a");
 	for (i = 0; i < element_n_xi; i++) {
		for (j = 0; j < element_n_eta; j++) {
			for (k = 0; k < division_ele_xi + 1; k++) {
				for (l = 0; l < division_ele_eta + 1; l++) {
                    if(coord_y[i * division_ele_xi + k][j * division_ele_eta + l] == 0.000000000000000e+00)
                    {
				        fprintf(fp, "%d\t% 1.15e\t% 1.15e\t% 1.15e\n", 
							graph_patch_n,
							coord_x[i * division_ele_xi + k][j * division_ele_eta + l], 
							coord_y[i * division_ele_xi + k][j * division_ele_eta + l], 
							stress_yy[i][j][k][l]);
                    }
				}
			}
		}
	}
	fclose(fp);

	fp = fopen("stress_vm_graph.txt", "a");
 	for (i = 0; i < element_n_xi; i++) {
		for (j = 0; j < element_n_eta; j++) {
			for (k = 0; k < division_ele_xi + 1; k++) {
				for (l = 0; l < division_ele_eta + 1; l++) {
					double stress_vm;
					double sum = stress_xx[i][j][k][l] + stress_yy[i][j][k][l];
					double dif = stress_xx[i][j][k][l] - stress_yy[i][j][k][l];
					double tau2 = stress_xy[i][j][k][l] * stress_xy[i][j][k][l];
					double temp1, temp2;
					double temp3;
					temp1 = 0.5 * sum;
					temp2 = 0.5 * sqrt(dif * dif + 4 * tau2);
					stress_vm = sqrt(temp1 * temp1 + 3 * temp2 * temp2);
					temp3 = sqrt(calc_xi[i * division_ele_xi + k] * calc_xi[i * division_ele_xi + k]
								+ calc_eta[j * division_ele_eta + l] * calc_eta[j * division_ele_eta + l]);
					fprintf(fp, "% 1.10e\t% 1.10e\t% 1.10e\t% 1.10e\t"
								"% 1.10e\t% 1.10e\n",
							calc_xi[i * division_ele_xi + k],
							calc_eta[j * division_ele_eta + l],
							coord_x[i * division_ele_xi + k][j * division_ele_eta + l],
							coord_y[i * division_ele_xi + k][j * division_ele_eta + l],
							stress_vm, temp3);
				}
			}
		}
	}
	fclose(fp);
}


//for s-IGA
//重ね合わせた結果の出力
static void Calculation_overlay(int order_xi_loc, int order_eta_loc,
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
							    double *weight_glo)
{
	int i, j, k, l;
	double temp1, temp2, temp3;
	//double temp_matrix[2][2];

    double output_xi, output_eta;
	double disp_x_glo;
	double disp_y_glo;
    double strain_xx_glo = 0;
    double strain_yy_glo = 0;
    double strain_xy_glo = 0;
	//, strain_yy_glo, strain_xy_glo;

	//計算するξ,ηの値決定と ∂ξ/∂チルダξ, ∂η/∂チルダη の計算
	double calc_xi_loc[MAX_POINTS];		//計算するξの値local
	double calc_eta_loc[MAX_POINTS];		//計算するηの値local
	//double dtilda_xi[MAX_KNOTS];		// ∂ξ/∂チルダξ
	//double dtilda_eta[MAX_KNOTS];	// ∂η/∂チルダη
	k = 0;
	l = 0;
	for (i = 0; i < knot_n_xi_loc - 1; i++) {
		if ( knot_vec_xi_loc[i] != knot_vec_xi_loc[i + 1] ) {
			calc_xi_loc[k] = knot_vec_xi_loc[i];
			printf("%d\t%f\n", k, calc_xi_loc[k]);
			//dtilda_xi[l] = ( knot_vec_xi_loc[i + 1] - knot_vec_xi_loc[i] ) / 2.0;
			//printf("%d\t%f\n", k, dtilda_xi[k]);
			k++;
			l++;
			if (division_ele_xi > 1) {
				temp1 = (knot_vec_xi_loc[i + 1] - knot_vec_xi_loc[i])
						/ (double)division_ele_xi;
				for (j = 1; j < division_ele_xi; j++) {
					calc_xi_loc[k] = calc_xi_loc[k - 1] + temp1;
					printf("%d\t%f\n", k, calc_xi_loc[k]);
					k++;
				}
			}
		}
	}
	calc_xi_loc[k] = knot_vec_xi_loc[knot_n_xi_loc - 1];
	printf("%d\t%f\n", k, calc_xi_loc[k]);
	//printf("\n");
	division_n_xi = k + 1;
	element_n_xi = l;

	k = 0;
	l = 0;
	for (i = 0; i < knot_n_eta_loc - 1; i++) {
		if ( knot_vec_eta_loc[i] != knot_vec_eta_loc[i + 1] ) {
			calc_eta_loc[k] = knot_vec_eta_loc[i];
			//printf("%d\t%f\n", k, calc_eta[k]);
			//dtilda_eta[l] = ( knot_vec_eta_loc[i + 1] - knot_vec_eta_loc[i] ) / 2.0;
			//printf("%d\t%f\n", k, dtilda_eta[k]);
			k++;
			l++;
			if (division_ele_eta > 1) {
				temp1 = (knot_vec_eta_loc[i + 1] - knot_vec_eta_loc[i])
						/ (double)division_ele_eta;
				for (j = 1; j < division_ele_eta; j++) {
					calc_eta_loc[k] = calc_eta_loc[k - 1] + temp1;
					//printf("%d\t%f\n", k, calc_eta[k]);
					k++;
				}
			}
		}
	}
	calc_eta_loc[k] = knot_vec_eta_loc[knot_n_eta_loc - 1];
	//printf("%d\t%f\n", k, calc_eta[k]);
	//printf("\n");
	division_n_eta = k + 1;
	element_n_eta = l;

	if (element_n_xi > MAX_ELEMENTS) {
		printf("Error!!\n");
		printf("Too many elements at xi!\n"
			   "Maximum of elements is %d (Now %d)\n"
			   "\n", MAX_ELEMENTS, element_n_xi);
		exit(1);
	}
	if (element_n_eta > MAX_ELEMENTS) {
		printf("Error!!\n");
		printf("Too many elements at eta!\n"
			   "Maximum of elements is %d (Now %d)\n"
			   "\n", MAX_ELEMENTS, element_n_eta);
		exit(1);
	}

	int ii, jj, kk, ll;

	//メッシュ座標計算
	printf("Start Calculation overlay mesh\n\n");
	printf("Start Calculation overlay displpacement\n\n");
	printf("Start Calculation overlay Strain\n\n");
	for (i = 0; i < division_n_xi; i++) {
		ii = i / division_ele_xi;
		kk = i % division_ele_xi;
		for (j = 0; j < division_n_eta; j++) {
			jj = j / division_ele_eta;
			ll = j % division_ele_eta;
			lNURBS_surface(knot_vec_xi_loc, knot_vec_eta_loc,
						   cntl_px_loc, cntl_py_loc, cntl_p_n_xi_loc, cntl_p_n_eta_loc,
						   weight_loc, order_xi_loc, order_eta_loc,
						   calc_xi_loc[i], calc_eta_loc[j],
						   &coord_x[i][j], &coord_y[i][j],
						   &dxi_x[ii][jj][kk][ll], &deta_x[ii][jj][kk][ll],
						   &dxi_y[ii][jj][kk][ll], &deta_y[ii][jj][kk][ll]);
			// printf("[%d][%d] [%d][%d][%d][%d]"
			// 	   "% 1.4e % 1.4e "
			// 	   "% 1.4e % 1.4e\n",
			// 	   i, j, ii, jj, kk, ll,
			// 	   calc_xi_loc[i], calc_eta_loc[j],
			// 	   coord_x[i][j], coord_y[i][j]);

            int itr_n = CalcXiEtaByNR(coord_x[i][j], coord_y[i][j],
                                      knot_vec_xi_glo, knot_vec_eta_glo,
				                      cntl_px_glo, cntl_py_glo,
			                          disp_cntl_px_glo, disp_cntl_py_glo,
			                          cntl_p_n_xi_glo, cntl_p_n_eta_glo,
			                          weight_glo, order_xi_glo, order_eta_glo,
			                          &output_xi, &output_eta,
									  &disp_x_glo, &disp_y_glo,
                                      &strain_xx_glo, &strain_yy_glo, &strain_xy_glo);
            if (itr_n == 0)
            {
                // printf("itr=0\n");
            }
			// printf("iteration : %d\n",itr_n);

			//ローカル内の表示点上のグローバル変位
			// printf("disp_x_glo =% 1.4e\tdisp_y_glo =% 1.4e\n", disp_x_glo, disp_y_glo);
			// printf("%1.4e\t%1.4e\n",disp_x[i][j],disp_y[i][j]);
			disp_x[i][j] += disp_x_glo;
			disp_y[i][j] += disp_y_glo;
			// printf("% 1.4e\t% 1.4e\n",disp_x[i][j],disp_y[i][j]);

			//ローカル内の表示点上のグローバルひずみ
            // printf("strain_xx_glo =% 1.4e\n"
            //        "strain_yy_glo =% 1.4e\n"
            //        "strain_xy_glo =% 1.4e\n",
            //        strain_xx_glo, strain_yy_glo, strain_xy_glo);
			strain_xx[ii][jj][kk][ll] += strain_xx_glo;
			strain_yy[ii][jj][kk][ll] += strain_yy_glo;
			strain_xy[ii][jj][kk][ll] += strain_xy_glo;
			// printf("test[%d][%d][%d][%d]\n",ii,jj,kk,ll);
			if (jj > 0 && ll == 0)
			{
				strain_xx[ii][jj-1][kk][division_ele_eta] += strain_xx_glo;
				strain_yy[ii][jj-1][kk][division_ele_eta] += strain_yy_glo;
				strain_xy[ii][jj-1][kk][division_ele_eta] += strain_xy_glo;
				// printf("test[%d][%d][%d][%d]\n",ii,jj-1,kk,division_ele_eta);
			}
			if (ii > 0 && kk == 0)
			{
				strain_xx[ii-1][jj][division_ele_xi][ll] += strain_xx_glo;
				strain_yy[ii-1][jj][division_ele_xi][ll] += strain_yy_glo;
				strain_xy[ii-1][jj][division_ele_xi][ll] += strain_xy_glo;
				// printf("test[%d][%d][%d][%d]\n",ii-1,jj,division_ele_xi,ll);
			}
			if (ii > 0 && jj > 0 && kk == 0 && ll == 0)
			{
				strain_xx[ii-1][jj-1][division_ele_xi][division_ele_eta] += strain_xx_glo;
				strain_yy[ii-1][jj-1][division_ele_xi][division_ele_eta] += strain_yy_glo;
				strain_xy[ii-1][jj-1][division_ele_xi][division_ele_eta] += strain_xy_glo;
				// printf("test[%d][%d][%d][%d]\n",ii-1,jj-1,division_ele_xi,division_ele_eta);
			}
			// printf("% 1.4e\t% 1.4e\t% 1.4e\n",
			// 		strain_xx[ii][jj][kk][ll],
			// 		strain_yy[ii][jj][kk][ll],
			// 		strain_xy[ii][jj][kk][ll]);

		}
		//printf("\n");
	}
	printf("\n");
	printf("End Calculation overlay mesh\n\n");
	printf("End Calculation overlay displpacement\n\n");
	printf("End Calculation overlay Strain\n\n");

	//Dマトリクスの計算
	double D_matrix[3][3] = {{0.0}};
	if (DM == 0) { //平面応力状態
		temp1 = E * (1.0 - nu * nu);
		D_matrix[0][0] = temp1;
		D_matrix[0][1] = nu * temp1;
		D_matrix[1][0] = nu * temp1;
		D_matrix[1][1] = temp1;
		D_matrix[2][2] = (1.0 - nu) / 2.0 * temp1;
	} else if (DM == 1) { //平面ひずみ状態(2Dの場合はこっち)
		temp1 = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
		D_matrix[0][0] = temp1;
		D_matrix[0][1] = nu / (1.0 - nu) * temp1;
		D_matrix[1][0] = nu / (1.0 - nu) * temp1;
		D_matrix[1][1] = temp1;
		D_matrix[2][2] = (1.0 - 2.0 * nu) / 2.0 / (1.0 - nu) * temp1;
	}

	printf("Start Calculation overlay Stress\n\n");
	for (i = 0; i < element_n_xi; i++) {
		for (j = 0; j < element_n_eta; j++) {
			for (k = 0; k < division_ele_xi + 1; k++) {
				for (l = 0; l < division_ele_eta + 1; l++) {
					stress_xx[i][j][k][l] = D_matrix[0][0] * strain_xx[i][j][k][l] + D_matrix[0][1] * strain_yy[i][j][k][l];
					stress_yy[i][j][k][l] = D_matrix[1][0] * strain_xx[i][j][k][l] + D_matrix[1][1] * strain_yy[i][j][k][l];
					stress_xy[i][j][k][l] = D_matrix[2][2] * strain_xy[i][j][k][l];
					// printf("[%d][%d][%d][%d]\t% 1.4e\t% 1.4e\t% 1.4e\n", i, j, k, l,
					// 	   stress_xx[i][j][k][l], stress_yy[i][j][k][l], stress_xy[i][j][k][l]);
				}
			}
			// printf("\n");
		}
	}
	printf("End Calculation overlay Stress\n\n");



	//書き込み
	fp = fopen("overlay_view.dat", "a");
	if (fields_flag) {
		fprintf(fp, "%d\t%d\t%d\t%d\n",
				division_n_xi, division_n_eta,
				element_n_xi, element_n_eta);
		for (i = 0; i < division_n_xi; i++) {
			for (j = 0; j < division_n_eta; j++) {
				temp1 = disp_x[i][j];
				temp2 = disp_y[i][j];
				temp3 = sqrt(temp1 * temp1 + temp2 * temp2);
				fprintf(fp, "% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\n",
						coord_x[i][j], coord_y[i][j],
						disp_x[i][j], disp_y[i][j], temp3);
			}
		}
		for (i = 0; i < element_n_xi; i++) {
			for (j = 0; j < element_n_eta; j++) {
				for (k = 0; k < division_ele_xi + 1; k++) {
					for (l = 0; l < division_ele_eta + 1; l++) {
						fprintf(fp, "% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\n",
								strain_xx[i][j][k][l], strain_yy[i][j][k][l], strain_xy[i][j][k][l],
								stress_xx[i][j][k][l], stress_yy[i][j][k][l], stress_xy[i][j][k][l]);
					}
				}
			}
		}
	} else {
		fprintf(fp, "%d\t%d\t%d\t%d\n",
				division_n_xi, division_n_eta,
				element_n_xi, element_n_eta);
		for (i = 0; i < division_n_xi; i++) {
			for (j = 0; j < division_n_eta; j++) {
				fprintf(fp, "% 1.4e\t% 1.4e\n",
						coord_x[i][j], coord_y[i][j]);
			}
		}
	}
	fclose(fp);
	//machino
	fp = fopen("overlay_view_r_theta.dat", "a");
	if (fields_flag) {
		fprintf(fp, "%d\t%d\t%d\t%d\n",
				division_n_xi, division_n_eta,
				element_n_xi, element_n_eta);
		for (i = 0; i < division_n_xi; i++) {
			for (j = 0; j < division_n_eta; j++) {
				temp1 = disp_x[i][j];
				temp2 = disp_y[i][j];
				temp3 = sqrt(temp1 * temp1 + temp2 * temp2);
				fprintf(fp, "% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\n",
						coord_x[i][j], coord_y[i][j],
						disp_x[i][j], disp_y[i][j], temp3);
			}
		}
		for (i = 0; i < element_n_xi; i++) {
			for (j = 0; j < element_n_eta; j++) {
				for (k = 0; k < division_ele_xi + 1; k++) {
					for (l = 0; l < division_ele_eta + 1; l++) {
					double stress_rr, stress_theta;
					double sum = stress_xx[i][j][k][l] + stress_yy[i][j][k][l];
					double dif = stress_xx[i][j][k][l] - stress_yy[i][j][k][l];
					double tau2 = stress_xy[i][j][k][l] * stress_xy[i][j][k][l];
					stress_rr = sum * 0.5
								+ sqrt(dif * dif + 4 * tau2) * 0.5;
					stress_theta = sum * 0.5
								  - sqrt(dif * dif + 4 * tau2) * 0.5;	
						fprintf(fp, "% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t% 1.4e\t0.0\n",
								strain_xx[i][j][k][l], strain_yy[i][j][k][l], strain_xy[i][j][k][l],
								stress_rr,stress_theta);
					}
				}
			}
		}
	} else {
		fprintf(fp, "%d\t%d\t%d\t%d\n",
				division_n_xi, division_n_eta,
				element_n_xi, element_n_eta);
		for (i = 0; i < division_n_xi; i++) {
			for (j = 0; j < division_n_eta; j++) {
				fprintf(fp, "% 1.4e\t% 1.4e\n",
						coord_x[i][j], coord_y[i][j]);
			}
		}
	}
	fclose(fp);
	//グラフ用ファイル書き込み
	fp = fopen("over_disp_graph.txt", "a");
	for (i = 0; i < division_n_xi; i++) {
		for (j = 0; j < division_n_eta; j++) {
			temp1 = disp_x[i][j];
			temp2 = disp_y[i][j];
			//temp3 = sqrt(temp1 * temp1 + temp2 * temp2);
			fprintf(fp, "%d\t% 1.15e\t% 1.15e\t% 1.15e\t% 1.15e\n",
					graph_patch_n,
					coord_x[i][j], coord_y[i][j],
					disp_x[i][j], disp_y[i][j]);
		}
	}
	fclose(fp);

	fp = fopen("over_stress_x_graph.txt", "a");
 	for (i = 0; i < element_n_xi; i++) {
		for (j = 0; j < element_n_eta; j++) {
			for (k = 0; k < division_ele_xi + 1; k++) {
				for (l = 0; l < division_ele_eta + 1; l++) {
				fprintf(fp, "% 1.10e\t% 1.10e\t% 1.10e\n",
							coord_x[i * division_ele_xi + k][j * division_ele_eta + l],
							coord_y[i * division_ele_xi + k][j * division_ele_eta + l],
							stress_xx[i][j][k][l]);
				}
			}
		}
	}
	fclose(fp);

	fp = fopen("over_stress_y_graph.txt", "a");
 	for (i = 0; i < element_n_xi; i++) {
		for (j = 0; j < element_n_eta; j++) {
			for (k = 0; k < division_ele_xi + 1; k++) {
				for (l = 0; l < division_ele_eta + 1; l++) {
				fprintf(fp, "% 1.15e\t% 1.15e\t% 1.15e\n",
							coord_x[i * division_ele_xi + k][j * division_ele_eta + l],
							coord_y[i * division_ele_xi + k][j * division_ele_eta + l],
							stress_yy[i][j][k][l]);
				}
			}
		}
	}
	fclose(fp);

	fp = fopen("over_stress_y_graph_0.txt", "a");
 	for (i = 0; i < element_n_xi; i++) {
		for (j = 0; j < element_n_eta; j++) {
			for (k = 0; k < division_ele_xi + 1; k++) {
				for (l = 0; l < division_ele_eta + 1; l++) {
                    if(coord_y[i * division_ele_xi + k][j * division_ele_eta + l] == 0.000000000000000e+00)
                    {
                        fprintf(fp, "% 1.15e\t% 1.15e\t% 1.15e\n", 
							coord_x[i * division_ele_xi + k][j * division_ele_eta + l], 
							coord_y[i * division_ele_xi + k][j * division_ele_eta + l], 
							stress_yy[i][j][k][l]);
                    }

				}
			}
		}
	}
	fclose(fp);

	fp = fopen("over_stress_r_theta_graph.txt", "a");
 	for (i = 0; i < element_n_xi; i++) {
		for (j = 0; j < element_n_eta; j++) {
			for (k = 0; k < division_ele_xi + 1; k++) {
				for (l = 0; l < division_ele_eta + 1; l++) {
					double stress_rr, stress_theta;
					double sum = stress_xx[i][j][k][l] + stress_yy[i][j][k][l];
					double dif = stress_xx[i][j][k][l] - stress_yy[i][j][k][l];
					double tau2 = stress_xy[i][j][k][l] * stress_xy[i][j][k][l];
					stress_rr = sum * 0.5
								+ sqrt(dif * dif + 4 * tau2) * 0.5;
					stress_theta = sum * 0.5
								  - sqrt(dif * dif + 4 * tau2) * 0.5;
					fprintf(fp, "% 1.10e\t% 1.10e\t% 1.10e\t% 1.10e\t% 1.10e\t% 1.10e\n",
							calc_xi_loc[i * division_ele_xi + k],
							calc_eta_loc[j * division_ele_eta + l],
							coord_x[i * division_ele_xi + k][j * division_ele_eta + l],
							coord_y[i * division_ele_xi + k][j * division_ele_eta + l],
							stress_rr, stress_theta);
				}
			}
		}
	}
	fclose(fp);

	fp = fopen("over_stress_vm_graph.txt", "a");
 	for (i = 0; i < element_n_xi; i++) {
		for (j = 0; j < element_n_eta; j++) {
			for (k = 0; k < division_ele_xi + 1; k++) {
				for (l = 0; l < division_ele_eta + 1; l++) {
					double stress_vm;
					double sum = stress_xx[i][j][k][l] + stress_yy[i][j][k][l];
					double dif = stress_xx[i][j][k][l] - stress_yy[i][j][k][l];
					double tau2 = stress_xy[i][j][k][l] * stress_xy[i][j][k][l];
					double temp1, temp2;
					temp1 = 0.5 * sum;
					temp2 = 0.5 * sqrt(dif * dif + 4 * tau2);
					stress_vm = sqrt(temp1 * temp1 + 3 * temp2 * temp2);
					fprintf(fp, "% 1.10e\t% 1.10e\t% 1.10e\t% 1.10e\t% 1.10e\n",
							calc_xi_loc[i * division_ele_xi + k],
							calc_eta_loc[j * division_ele_eta + l],
							coord_x[i * division_ele_xi + k][j * division_ele_eta + l],
							coord_y[i * division_ele_xi + k][j * division_ele_eta + l],
							stress_vm);
				}
			}
		}
	}
	fclose(fp);
}


/*
static void Calculation_at_GP(double E, double nu)
{
	//通常IGAでのガウス点での値
	int i, j, k, e;

	Make_gauss_array(1);

	//メッシュ座標計算
	int ele_glo_n = real_Total_Element_on_mesh[0];// グローバルメッシュの要素数

	double U_temp[MAX_KIEL_SIZE];
	double B[D_MATRIX_SIZE][MAX_KIEL_SIZE];
	double J;

	for (e = 0; e < ele_glo_n; e++)
	{
		double X_temp[No_Control_point_ON_ELEMENT[Element_patch[e]]][DIMENSION];

		// printf("ele = %d\n\n", e);

		//strain_GPの初期化
		for (i = 0; i < GP_2D; i++)
		{
			for (j = 0; j < 3; j++)
			{
				strain_GP[e][i][j] = 0.0;
			}
		}

		// printf("x\ty\n");
		for (j = 0; j < No_Control_point_ON_ELEMENT[Element_patch[e]]; j++)
		{
			for (k = 0; k < DIMENSION; k++)
			{
				U_temp[j * DIMENSION + k] = Displacement[Controlpoint_of_Element[e][j] * DIMENSION + k];
				X_temp[j][k] = Node_Coordinate[Controlpoint_of_Element[e][j]][k];
			}
			// printf("%.15e\t%.15e\n", X_temp[j][0], X_temp[j][1]);
		}
		// printf("\n");

		for (i = 0; i < GP_2D; i++)	//ガウス点のループ
		{
			double data_result_shape[2] = {0.0, 0.0};
			double R_shape_func;

			// printf("Gxi_x\t_Gxi_y\tR_shape_func\n");
			for (j = 0; j < No_Control_point_ON_ELEMENT[Element_patch[e]]; j++)
			{
				R_shape_func = Shape_func(j, Total_Control_Point_to_mesh[Total_mesh], Gxi[i], e);
				for (k = 0; k < DIMENSION; k++)
				{
					data_result_shape[k] += R_shape_func * X_temp[j][k];
				}
			}

			//物理座標[要素番号(ローカル内で0から始まる)][ガウス点番号][DIMENSION]
			for (j = 0; j < DIMENSION; j++)
			{
				coordinate_GP[e][i][j] = data_result_shape[j];
			}
		}

		for (i = 0; i < GP_2D; i++)	//ガウス点のループ
		{
			Make_B_Matrix(e, B, Gxi[i], X_temp, &J, Total_Control_Point_to_mesh[Total_mesh]);
			for (j = 0; j < D_MATRIX_SIZE; j++)
			{
				for (k = 0; k < KIEL_SIZE; k++)
				{
					strain_GP[e][i][j] += B[j][k] * U_temp[k];
				}
			}
			Jac[e][i] = J;
		}
	}

	//Dマトリクスの計算
	double temp1;
	double D_matrix[3][3] = {{0.0}};
	if (DM == 0) { //平面応力状態
		temp1 = E * (1.0 - nu * nu);
		D_matrix[0][0] = temp1;
		D_matrix[0][1] = nu * temp1;
		D_matrix[1][0] = nu * temp1;
		D_matrix[1][1] = temp1;
		D_matrix[2][2] = (1.0 - nu) / 2.0 * temp1;
	} else if (DM == 1) { //平面ひずみ状態(2Dの場合はこっち)
		temp1 = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
		D_matrix[0][0] = temp1;
		D_matrix[0][1] = nu / (1.0 - nu) * temp1;
		D_matrix[1][0] = nu / (1.0 - nu) * temp1;
		D_matrix[1][1] = temp1;
		D_matrix[2][2] = (1.0 - 2.0 * nu) / 2.0 / (1.0 - nu) * temp1;
	}

	for (e = 0; e < ele_glo_n; e++)
	{
		for (i = 0; i < GP_2D; i++)
		{
			stress_GP[e][i][0] = D_matrix[0][0] * strain_GP[e][i][0] + D_matrix[0][1] * strain_GP[e][i][1];
			stress_GP[e][i][1] = D_matrix[1][0] * strain_GP[e][i][0] + D_matrix[1][1] * strain_GP[e][i][1];
			stress_GP[e][i][2] = D_matrix[2][2] * strain_GP[e][i][2];
		}
	}

	//座標変換
	double theta = 0.0;

	for (e = 0; e < ele_glo_n; e++)
	{
		for (i = 0; i < GP_2D; i++)
		{
			theta = atan2(coordinate_GP[e][i][1], coordinate_GP[e][i][0]);
			
			stress_r_theta_GP[e][i][0] = stress_GP[e][i][0] * pow(cos(theta), 2.0)
									   + stress_GP[e][i][1] * pow(sin(theta), 2.0)
									   + 2.0 * stress_GP[e][i][2] * sin(theta) * cos(theta);
			stress_r_theta_GP[e][i][1] = stress_GP[e][i][0] * pow(sin(theta), 2.0)
									   + stress_GP[e][i][1] * pow(cos(theta), 2.0)
									   - 2.0 * stress_GP[e][i][2] * cos(theta) * sin(theta);
			stress_r_theta_GP[e][i][2] = (stress_GP[e][i][1] - stress_GP[e][i][0])
									   * sin(theta) * cos(theta) + stress_GP[e][i][2]
									   * (pow(cos(theta), 2.0) - pow(sin(theta), 2.0));
		}
	}

	//厚肉円筒の理論解
	double r_t = 0.0;

	for (e = 0; e < ele_glo_n; e++)
	{
		for (i = 0; i < GP_2D; i++)
		{
			r_t = sqrt(pow(coordinate_GP[e][i][0], 2.0) + pow(coordinate_GP[e][i][1], 2.0));

			stress_theory_r_theta[e][i][0] = (pow(r_t, 2.0) - 4.0) / (pow(r_t, 2.0) * 3.0);
			stress_theory_r_theta[e][i][1] = (pow(r_t, 2.0) + 4.0) / (pow(r_t, 2.0) * 3.0);
		}
	}

	//書き込み
	fp = fopen("at_GP_overlay_data.txt", "w");
	fprintf(fp, "e\tガウス番号\tx\ty\tstress_xx\tstress_yy\tstress_r\tstress_theta\n");
	for (e = 0; e < ele_glo_n; e++)
	{
		for (i = 0; i < GP_2D; i++)
		{
			fprintf(fp, "%d\t%d\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n", e, i, coordinate_GP[e][i][0], coordinate_GP[e][i][1], stress_GP[e][i][0], stress_GP[e][i][1], stress_r_theta_GP[e][i][0], stress_r_theta_GP[e][i][1]);
		}
	}
	fclose(fp);

	fp = fopen("at_GP_overlay_for_error_norm.txt", "w");
	fprintf(fp, "e\tガウス番号\tx\ty\tstress_r-theory\tstress_theta-theory\n");
	for (e = 0; e < ele_glo_n; e++)
	{
		for (i = 0; i < GP_2D; i++)
		{
			fprintf(fp, "%d\t%d\t%.15e\t%.15e\t%.15e\t%.15e\n", e, i, coordinate_GP[e][i][0], coordinate_GP[e][i][1], stress_r_theta_GP[e][i][0] - stress_theory_r_theta[e][i][0], stress_r_theta_GP[e][i][1] - stress_theory_r_theta[e][i][1]);
		}
	}
	fclose(fp);


	//error normを計算
	//ガウス点で出したtheoryとの差の二乗を面積分
	double temp2 = 0.0, temp3 = 0.0, temp4 = 0.0, temp5 = 0.0;
	for (e = 0; e < ele_glo_n; e++)
	{
		for (i = 0; i < GP_2D; i++)
		{
			temp2 += w[i] * pow((stress_r_theta_GP[e][i][0] - stress_theory_r_theta[e][i][0]), 2.0) * Jac[e][i];
			temp3 += w[i] * pow((stress_r_theta_GP[e][i][1] - stress_theory_r_theta[e][i][1]), 2.0) * Jac[e][i];
			temp4 += w[i] * pow(stress_r_theta_GP[e][i][0], 2.0) * Jac[e][i];
			temp5 += w[i] * pow(stress_r_theta_GP[e][i][1], 2.0) * Jac[e][i];
		}
	}
	fp = fopen("at_GP_overlay_for_error_norm_surface_integral.txt", "w");
	fprintf(fp, "(stress_r-theory)^2_surface_integral\t(stress_theta-thory)^2_surface_integral\tstress_r^2_surface_integral\n");
	fprintf(fp, "%.15e\t%.15e\t%.15e\t%.15e\n", temp2, temp3, temp4, temp5);
	fclose(fp);
}


void Calculation_overlay_at_GP(double E, double nu,
							   int order_xi_glo, int order_eta_glo,
							   int knot_n_xi_glo, int knot_n_eta_glo,
							   int cntl_p_n_xi_glo, int cntl_p_n_eta_glo,
							   double *knot_vec_xi_glo, double *knot_vec_eta_glo,
							   double *cntl_px_glo, double *cntl_py_glo,
							   double *disp_cntl_px_glo, double *disp_cntl_py_glo,
							   double *weight_glo)
{
	//s-IGAでのローカル上のガウス点での重ね合わせた値
	int i, j, k, e;

	double disp_glo[DIMENSION];
	double strain_glo[D_MATRIX_SIZE];

	Make_gauss_array(1);

	double G_GP_knot[GP_2D][DIMENSION];

	//メッシュ座標計算
	int ele_glo_n = real_Total_Element_on_mesh[0];// グローバルメッシュの要素数
	int ele_loc_n = real_Total_Element_to_mesh[Total_mesh] - real_Total_Element_on_mesh[0];// ローカルメッシュの要素数

	double U_temp[MAX_KIEL_SIZE];
	double B[D_MATRIX_SIZE][MAX_KIEL_SIZE];
	double J;

	for (e = 0; e < ele_loc_n; e++)
	{
		int El_No_loc = ele_glo_n + e;
		double X_temp[No_Control_point_ON_ELEMENT[Element_patch[El_No_loc]]][DIMENSION];

		//strain_GPの初期化
		for (i = 0; i < GP_2D; i++)
		{
			for (j = 0; j < 3; j++)
			{
				strain_GP[e][i][j] = 0.0;
			}
		}

		for (j = 0; j < No_Control_point_ON_ELEMENT[Element_patch[El_No_loc]]; j++)
		{
			for (k = 0; k < DIMENSION; k++)
			{
				U_temp[j * DIMENSION + k] = Displacement[Controlpoint_of_Element[El_No_loc][j] * DIMENSION + k];
				X_temp[j][k] = Node_Coordinate[Controlpoint_of_Element[El_No_loc][j]][k];
			}
		}

		for (i = 0; i < GP_2D; i++)	//ガウス点のループ
		{
			double data_result_shape[2] = {0.0, 0.0};
			double R_shape_func;

			for (j = 0; j < No_Control_point_ON_ELEMENT[Element_patch[El_No_loc]]; j++)
			{
				R_shape_func = Shape_func(j, Total_Control_Point_to_mesh[Total_mesh], Gxi[i], El_No_loc);
				for (k = 0; k < DIMENSION; k++)
				{
					data_result_shape[k] += R_shape_func * X_temp[j][k];
				}
			}

			for (j = 0; j < DIMENSION; j++)
			{
				coordinate_GP[e][i][j] = data_result_shape[j];
			}
		}

		for (i = 0; i < GP_2D; i++)	//ガウス点のループ
		{
			Make_B_Matrix(El_No_loc, B, Gxi[i], X_temp, &J, Total_Control_Point_to_mesh[Total_mesh]);
			for (j = 0; j < D_MATRIX_SIZE; j++)
			{
				for (k = 0; k < KIEL_SIZE; k++)
				{
					strain_GP[e][i][j] += B[j][k] * U_temp[k];
				}
			}
			Jac[e][i] = J;
		}

		for (i = 0; i < GP_2D; i++)	//ガウス点のループ
		{
			CalcXiEtaByNR(coordinate_GP[e][i][0], coordinate_GP[e][i][1],
						  knot_vec_xi_glo, knot_vec_eta_glo,
						  cntl_px_glo, cntl_py_glo,
						  disp_cntl_px_glo, disp_cntl_py_glo,
						  cntl_p_n_xi_glo, cntl_p_n_eta_glo,
						  weight_glo, order_xi_glo, order_eta_glo,
						  &G_GP_knot[i][0], &G_GP_knot[i][1],
						  &disp_glo[0], &disp_glo[1],
						  &strain_glo[0], &strain_glo[1], &strain_glo[2]);

			//重ね合わせ
			for (j = 0; j < 3; j++) //xx, yy, xyを重ね合わせる
			{
				strain_GP[e][i][j] += strain_glo[j];
			}
		}
	}

	//Dマトリクスの計算
	double temp1;
	double D_matrix[3][3] = {{0.0}};
	if (DM == 0) { //平面応力状態
		temp1 = E * (1.0 - nu * nu);
		D_matrix[0][0] = temp1;
		D_matrix[0][1] = nu * temp1;
		D_matrix[1][0] = nu * temp1;
		D_matrix[1][1] = temp1;
		D_matrix[2][2] = (1.0 - nu) / 2.0 * temp1;
	} else if (DM == 1) { //平面ひずみ状態(2Dの場合はこっち)
		temp1 = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
		D_matrix[0][0] = temp1;
		D_matrix[0][1] = nu / (1.0 - nu) * temp1;
		D_matrix[1][0] = nu / (1.0 - nu) * temp1;
		D_matrix[1][1] = temp1;
		D_matrix[2][2] = (1.0 - 2.0 * nu) / 2.0 / (1.0 - nu) * temp1;
	}

	for (e = 0; e < ele_loc_n; e++)
	{
		for (i = 0; i < GP_2D; i++)
		{
			stress_GP[e][i][0] = D_matrix[0][0] * strain_GP[e][i][0] + D_matrix[0][1] * strain_GP[e][i][1];
			stress_GP[e][i][1] = D_matrix[1][0] * strain_GP[e][i][0] + D_matrix[1][1] * strain_GP[e][i][1];
			stress_GP[e][i][2] = D_matrix[2][2] * strain_GP[e][i][2];
		}
	}

	//座標変換
	double theta = 0.0;

	for (e = 0; e < ele_loc_n; e++)
	{
		for (i = 0; i < GP_2D; i++)
		{
			theta = atan2(coordinate_GP[e][i][1], coordinate_GP[e][i][0]);
			
			stress_r_theta_GP[e][i][0] = stress_GP[e][i][0] * pow(cos(theta), 2.0)
									   + stress_GP[e][i][1] * pow(sin(theta), 2.0)
									   + 2.0 * stress_GP[e][i][2] * sin(theta) * cos(theta);
			stress_r_theta_GP[e][i][1] = stress_GP[e][i][0] * pow(sin(theta), 2.0)
									   + stress_GP[e][i][1] * pow(cos(theta), 2.0)
									   - 2.0 * stress_GP[e][i][2] * cos(theta) * sin(theta);
			stress_r_theta_GP[e][i][2] = (stress_GP[e][i][1] - stress_GP[e][i][0])
									   * sin(theta) * cos(theta) + stress_GP[e][i][2]
									   * (pow(cos(theta), 2.0) - pow(sin(theta), 2.0));
		}
	}

	//円孔を有する無限平板の理論解
	double r_t = 0.0, theta_t = 0.0;

	for (e = 0; e < ele_loc_n; e++)
	{
		for (i = 0; i < GP_2D; i++)
		{
			r_t = sqrt(pow(coordinate_GP[e][i][0], 2.0) + pow(coordinate_GP[e][i][1], 2.0));
			theta_t = atan2(coordinate_GP[e][i][1], coordinate_GP[e][i][0]);
			
			// stress_theory_r_theta[e][i][0] = (10.0 / 2.0) * (1.0 - pow((1.0 / r_t), 2.0)) 
			// 							   + (10.0 / 2.0) * (1.0 + 3.0 * pow((1.0 / r_t), 4.0)
			// 							   - 4.0 * pow((1.0 / r_t), 2.0)) * cos(2.0 * theta_t);
			// stress_theory_r_theta[e][i][1] = (10.0 / 2.0) * (1.0 + pow((1.0 / r_t), 2.0)) 
			// 							   - (10.0 / 2.0) * (1.0 + 3.0 * pow((1.0 / r_t), 4.0))
			// 							   * cos(2.0 * theta_t);
			// stress_theory_r_theta[e][i][2] = - (10.0 / 2.0) * (1.0 - 3.0 * pow((1.0 / r_t), 4.0)
			// 							   + 2.0 * pow((1.0 / r_t), 2.0)) * sin(2.0 * theta_t);
			stress_theory_r_theta[e][i][0] = (10.0 / 2.0) * (1.0 - pow((1.0 / r_t), 2.0)) 
										   - (10.0 / 2.0) * (1.0 + 3.0 * pow((1.0 / r_t), 4.0)
										   - 4.0 * pow((1.0 / r_t), 2.0)) * cos(2.0 * theta_t);
			stress_theory_r_theta[e][i][1] = (10.0 / 2.0) * (1.0 + pow((1.0 / r_t), 2.0)) 
										   + (10.0 / 2.0) * (1.0 + 3.0 * pow((1.0 / r_t), 4.0))
										   * cos(2.0 * theta_t);
			stress_theory_r_theta[e][i][2] = (10.0 / 2.0) * (1.0 - 3.0 * pow((1.0 / r_t), 4.0)
										   + 2.0 * pow((1.0 / r_t), 2.0)) * sin(2.0 * theta_t);
		}
	}

	//書き込み
	fp = fopen("at_GP_overlay_data.txt", "w");
	fprintf(fp, "e\tガウス番号\tx\ty\tstress_xx\tstress_yy\tstress_r\tstress_theta\n");
	for (e = 0; e < ele_loc_n; e++)
	{
		for (i = 0; i < GP_2D; i++)
		{
			fprintf(fp, "%d\t%d\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n", e, i, coordinate_GP[e][i][0], coordinate_GP[e][i][1], stress_GP[e][i][0], stress_GP[e][i][1], stress_r_theta_GP[e][i][0], stress_r_theta_GP[e][i][1]);
		}
	}
	fclose(fp);

	fp = fopen("at_GP_overlay_for_error_norm.txt", "w");
	fprintf(fp, "e\tガウス番号\tx\ty\tstress_r-theory\tstress_theta-theory\n");
	for (e = 0; e < ele_loc_n; e++)
	{
		for (i = 0; i < GP_2D; i++)
		{
			fprintf(fp, "%d\t%d\t%.15e\t%.15e\t%.15e\t%.15e\n", e, i, coordinate_GP[e][i][0], coordinate_GP[e][i][1], stress_r_theta_GP[e][i][0] - stress_theory_r_theta[e][i][0], stress_r_theta_GP[e][i][1] - stress_theory_r_theta[e][i][1]);
		}
	}
	fclose(fp);


	//error normを計算
	//ガウス点で出したtheoryとの差の二乗を面積分
	double temp2 = 0.0, temp3 = 0.0, temp4 = 0.0, temp5 = 0.0, temp6 = 0.0;
	for (e = 0; e < ele_loc_n; e++)
	{
		for (i = 0; i < GP_2D; i++)
		{
			temp2 += w[i] * pow((stress_r_theta_GP[e][i][0] - stress_theory_r_theta[e][i][0]), 2.0) * Jac[e][i];
			temp3 += w[i] * pow((stress_r_theta_GP[e][i][1] - stress_theory_r_theta[e][i][1]), 2.0) * Jac[e][i];
			temp4 += w[i] * pow(stress_r_theta_GP[e][i][0], 2.0) * Jac[e][i];
			temp5 += w[i] * pow(stress_r_theta_GP[e][i][1], 2.0) * Jac[e][i];
			temp6 += w[i] * Jac[e][i];
		}
	}
	fp = fopen("at_GP_overlay_for_error_norm_surface_integral.txt", "w");
	fprintf(fp, "(stress_r-theory)^2_surface_integral\t(stress_theta-thory)^2_surface_integral\tstress_r^2_surface_integral\tstress_theta^2_surface_integral\t面積(analysis)\n");
	fprintf(fp, "%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n", temp2, temp3, temp4, temp5, temp6);
	fclose(fp);
}
*/


void K_output_svg(int ndof)
{
	// [K] = [[K^G, K^GL], [K^GL, K^L]]

	int i, j;

	char color_vec[2][10] = {"#f5f5f5", "#ee82ee"};
	// 0	whitesmoke
	// 1	violet
	// https://www.colordic.org/

	double space = 3.0, scale = 1000.0 / (((double)ndof) + 2.0 * space);

	double width = (((double)ndof) + 2.0 * space) * scale;
	double height = width;

	char str[256] = "K_matrix.svg";
	fp = fopen(str, "w");

	fprintf(fp, "<?xml version='1.0'?>\n");
    // fprintf(fp, "<svg width='%lept' height='%lept' viewBox='0 0 %le %le' style = 'background: #eee' xmlns='http://www.w3.org/2000/svg' xmlns:xlink='http://www.w3.org/1999/xlink'>\n", width, height, width, height);
    fprintf(fp, "<svg width='%le' height='%le' version='1.1' style='background: #eee' xmlns='http://www.w3.org/2000/svg' xmlns:xlink='http://www.w3.org/1999/xlink'>\n", width, height);

	double xx = space * scale;
	double yy = space * scale;
	double ww = ndof * scale;
	double hh = ndof * scale;
	fprintf(fp, "<rect x='%le' y='%le' width='%le' height='%le' fill='%s' />\n", xx, yy, ww, hh, color_vec[0]);

	// 各行の成分を抽出
	for (i = 0; i < ndof; i++)
	{
		int K_bool[ndof]; //一行分保存する
		for (j = 0; j < ndof; j++)
		{
			K_bool[j] = 0;
		}

		for (j = 0; j < ndof; j++)
		{
			int temp_count;
			if (i <= j)
			{
				temp_count = RowCol_to_icount(i, j);
			}
			else if (i > j)
			{
				temp_count = RowCol_to_icount(j, i);
			}

			if(temp_count != -1)
			{
				K_bool[j] = 1;
			}
		}

		for (j = 0; j < ndof; j++)
		{
			double x = (((double)j) + space) * scale;
			double y = (((double)i) + space) * scale;
			// if (K_bool[j] == 0)
			// {
			// 	fprintf(fp, "<rect x='%le' y='%le' width='%le' height='%le' fill='%s' />\n", x, y, scale, scale, color_vec[0]);
			// }
			if (K_bool[j] == 1)
			{
				fprintf(fp, "<rect x='%le' y='%le' width='%le' height='%le' fill='%s' />\n", x, y, scale, scale, color_vec[1]);
			}
		}
	}

	fprintf(fp, "</svg>");
    fclose(fp);
}


void Make_gauss_array(int select_GP)
{
	int i, j;

	if (select_GP == 0)
	{
		GP_1dir = Ng;
	}
	else if (select_GP == 1)
	{
		GP_1dir = Ng_extended;
	}

	GP_2D = GP_1dir * GP_1dir;

	if (GP_1dir == 3)
	{
		double G1 = pow((3.0 / 5.0), 0.5);
		double G_vec[3] = {-G1, 0.0, G1};
		double w1 = 8.0 / 9.0;
		double w2 = 5.0 / 9.0;
		double w_vec[3] = {w2, w1, w2};

		for (i = 0; i < GP_1dir; i++)
		{
			for (j = 0; j < GP_1dir; j++)
			{
				w[j + (GP_1dir * i)] = w_vec[i] * w_vec[j];
				Gxi[j + (GP_1dir * i)][0] = G_vec[j];
				Gxi[j + (GP_1dir * i)][1] = G_vec[i];
			}
		}
	}
	else if (GP_1dir == 4)
	{
		double A = pow((6.0 / 5.0), 0.5);
		double G1 = pow(((3.0 - 2.0 * A) / 7.0), 0.5);
		double G2 = pow(((3.0 + 2.0 * A) / 7.0), 0.5);
		double G_vec[4] = {-G2, -G1, G1, G2};
		double B = pow(30.0, 0.5);
		double w1 = (18.0 + B) / 36.0;
		double w2 = (18.0 - B) / 36.0;
		double w_vec[4] = {w2, w1, w1, w2};

		for (i = 0; i < GP_1dir; i++)
		{
			for (j = 0; j < GP_1dir; j++)
			{
				w[j + (GP_1dir * i)] = w_vec[i] * w_vec[j];
				Gxi[j + (GP_1dir * i)][0] = G_vec[j];
				Gxi[j + (GP_1dir * i)][1] = G_vec[i];
			}
		}
	}
	else if (GP_1dir == 5)
	{
		double A = pow((10.0 / 7.0), 0.5);
		double G1 = pow((5.0 - 2.0 * A), 0.5) / 3.0;
		double G2 = pow((5.0 + 2.0 * A), 0.5) / 3.0;
		double G_vec[5] = {-G2, -G1, 0.0, G1, G2};
		double B = pow(70.0, 0.5);
		double w1 = 128.0 / 225.0;
		double w2 = (322.0 + 13.0 * B) / 900.0;
		double w3 = (322.0 - 13.0 * B) / 900.0;
		double w_vec[5] = {w3, w2, w1, w2, w3};

		for (i = 0; i < GP_1dir; i++)
		{
			for (j = 0; j < GP_1dir; j++)
			{
				w[j + (GP_1dir * i)] = w_vec[i] * w_vec[j];
				Gxi[j + (GP_1dir * i)][0] = G_vec[j];
				Gxi[j + (GP_1dir * i)][1] = G_vec[i];
			}
		}
	}
	else if (GP_1dir == 10)
	{
		double G_vec[10];
		double w_vec[10];

		G_vec[0]  = -0.9739065285171717;
		G_vec[1]  = -0.8650633666889845;
		G_vec[2]  = -0.6794095682990244;
		G_vec[3]  = -0.4333953941292472;
		G_vec[4]  = -0.1488743389816312;
		G_vec[5]  =  0.1488743389816312;
		G_vec[6]  =  0.4333953941292472;
		G_vec[7]  =  0.6794095682990244;
		G_vec[8]  =  0.8650633666889845;
		G_vec[9]  =  0.9739065285171717;

		w_vec[0]  = 0.0666713443086881;
		w_vec[1]  = 0.1494513491505804;
		w_vec[2]  = 0.2190863625159820;
		w_vec[3]  = 0.2692667193099965;
		w_vec[4]  = 0.2955242247147530;
		w_vec[5]  = 0.2955242247147530;
		w_vec[6]  = 0.2692667193099965;
		w_vec[7]  = 0.2190863625159820;
		w_vec[8]  = 0.1494513491505804;
		w_vec[9]  = 0.0666713443086881;

		for (i = 0; i < GP_1dir; i++)
		{
			for (j = 0; j < GP_1dir; j++)
			{
				w[j + (GP_1dir * i)] = w_vec[i] * w_vec[j];
				Gxi[j + (GP_1dir * i)][0] = G_vec[j];
				Gxi[j + (GP_1dir * i)][1] = G_vec[i];
			}
		}
	}
}


void Make_Displacement_grad_glo_check(int Total_Element)
{
	int i, j, j1, j2;
	int re;

	for (re = 0; re < Total_Element; re++)
	{
		i = real_element[re];
        // printf("El_No;i=%d\n", real_element[re]);
		double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION];
		// printf("Total_Element=%d\tre=%d\tEl_No=%d\n", Total_Element, re, i);

		for (j1 = 0; j1 < No_Control_point_ON_ELEMENT[Element_patch[i]]; j1++)
		{
			for (j2 = 0; j2 < DIMENSION; j2++)
			{
				X[j1][j2] = Node_Coordinate[Controlpoint_of_Element[i][j1]][j2];
				// printf("X[%d][%d] = %.10e\n", j1, j2, Node_Coordinate[Controlpoint_of_Element[i][j1]][j2]);
			}
		}

		if (Element_mesh[i] > 0)	//ローカルメッシュ上の要素について
		{
			if (NNLOVER[i] > 0)		//重なっている要素が存在するとき
			{
				for (j = 0; j < NNLOVER[i]; j++)
				{
					double XG[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION];
					for (j1 = 0; j1 < No_Control_point_ON_ELEMENT[Element_patch[NELOVER[i][j]]]; j1++)
					{
						for (j2 = 0; j2 < DIMENSION; j2++)
						{
							XG[j1][j2] = Node_Coordinate[Controlpoint_of_Element[NELOVER[i][j]][j1]][j2];
							//重なっている要素の物理座標取得
							// printf("XG[%d][%d] = %.10e\n", j1, j2, Node_Coordinate[Controlpoint_of_Element[NELOVER[i][j]][j1]][j2]);
						}
					}
					Make_Displacement_grad_glo(i, NELOVER[i][j], X, XG);
				}
			}
		}
	}
}


int Make_Displacement_grad_glo(int El_No_loc, int El_No_glo,
					  		   double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION],
					  		   double XG[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION])
{
	int i, j, jj, k;
	int BDBJ_flag;

	static double U[MAX_KIEL_SIZE];
	double B_glo[D_MATRIX_SIZE][MAX_KIEL_SIZE], b_grad_glo[DIMENSION * DIMENSION][2 * MAX_NO_CCpoint_ON_ELEMENT];
	// double X_glo[No_Control_point_ON_ELEMENT[Element_patch[El_No_glo]]][DIMENSION];
	double J = 0.0;

	Make_gauss_array(0);

	double G_Gxi[GP_2D][DIMENSION];	//グローバルパッチ上での親要素内座標xi_bar,eta_bar

	for(i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No_glo]]; i++){
		for(j = 0; j < DIMENSION; j++){
			U[i * DIMENSION + j] = Displacement[Controlpoint_of_Element[El_No_glo][i] * DIMENSION + j];
			// printf("Controlpoint_of_Element[%d][%d] = %d\n", El_No_glo, i, Controlpoint_of_Element[El_No_glo][i]);
			// printf("Displacement[Controlpoint_of_Element[%d][%d]] = %.10e\n", El_No_glo, i, Displacement[Controlpoint_of_Element[El_No_glo][i] * DIMENSION + j]);
		}
	}
	for (i = 0; i < GP_2D; i++)	//ガウス点のループ(local)
	{
		// printf("gauss point number:%d\n", i);
		////ローカルガウス点がグローバル要素に含まれているかの判定
		//ローカル要素ガウス点の物理座標算出
		double data_result_shape[2] = {0.0};
		double output_xi, output_eta;
		int patch_n = 0;
		for (j = 0; j < No_Control_point_ON_ELEMENT[Element_patch[El_No_loc]]; j++)
		{
			double R_shape_func = Shape_func(j, Gxi[i], El_No_loc);
			for (jj = 0; jj < DIMENSION; jj++)
			{
				data_result_shape[jj] += R_shape_func * X[j][jj];  //Node_coordinate[controlpoint_of_element[element_n_loc][j]][jj]
			}
		}

		//ローカル要素ガウス点のグローバルパッチ上のパラメータ空間座標算出
		for (j = 0; j < Total_Patch_on_mesh[0]; j++)	//グローバルメッシュ[0]上
		{
			Calc_xi_eta(data_result_shape[0], data_result_shape[1],
	            	    Position_Knots[j][0], Position_Knots[j][1],
	                	No_Control_point[j][0], No_Control_point[j][1], Order[j][0], Order[j][1],
	                	&output_xi, &output_eta);
			// printf("  x: % 1.8e\n", data_result_shape[0]);
			// printf("  y: % 1.8e\n", data_result_shape[1]);
			// printf(" xi: % 1.8e\n", output_xi);
			// printf("eta: % 1.8e\n", output_eta);
			// printf("patch_n: %d\n", j);
			patch_n = j;
		}
		//要素内外判定
		
		if (output_xi >= Position_Knots[patch_n][0][Order[patch_n][0]+ENC[patch_n][El_No_glo][0]] 
			&&
			output_xi < Position_Knots[patch_n][0][Order[patch_n][0]+ENC[patch_n][El_No_glo][0]+1]
			&&
			output_eta >= Position_Knots[patch_n][1][Order[patch_n][1]+ENC[patch_n][El_No_glo][1]]
			&&
			output_eta < Position_Knots[patch_n][1][Order[patch_n][1]+ENC[patch_n][El_No_glo][1]+1])	//要素内であるとき
		{
			BDBJ_flag = 1;
			printf("BDBJ_flag\n");

			//親要素座標の算出
			G_Gxi[i][0] = -1 
					  + 2 * (output_xi - Position_Knots[patch_n][0][Order[patch_n][0]+ENC[patch_n][El_No_glo][0]]) /
					  (Position_Knots[patch_n][0][Order[patch_n][0]+ENC[patch_n][El_No_glo][0]+1] - Position_Knots[patch_n][0][Order[patch_n][0]+ENC[patch_n][El_No_glo][0]]);
			G_Gxi[i][1] = -1 
					  + 2 * (output_eta - Position_Knots[patch_n][1][Order[patch_n][1]+ENC[patch_n][El_No_glo][1]]) /
					  (Position_Knots[patch_n][1][Order[patch_n][1]+ENC[patch_n][El_No_glo][1]+1] - Position_Knots[patch_n][1][Order[patch_n][1]+ENC[patch_n][El_No_glo][1]]);
			// printf("G_Gxi[][]=\n");
		}
		else	//要素外であるとき
		{
			BDBJ_flag = 0;
		}

		//printf("i=%d\n",i );
		
		////結合要素剛性マトリックス計算
		//要素内であるとき、次を計算
		if (BDBJ_flag)
		{
			//printf("BDBJ_flag\ti=%d\n",i );
			//重なるグローバル要素のBマトリックス
			// Make_B_Matrix(El_No_glo, BG, G_Gxi[i], XG, &J, 
			// 			  Total_Control_Point_to_mesh[Total_mesh]);
			Make_B_Matrix(El_No_glo, B_glo, G_Gxi[i], XG, &J);
			Make_b_grad_Matrix(El_No_glo, b_grad_glo, G_Gxi[i], XG, &J);
			for(j = 0; j < DIMENSION * DIMENSION; j++){
				for(k = 0; k < 2 * No_Control_point_ON_ELEMENT[Element_patch[El_No_glo]]; k++){
					Disp_grad_glo[El_No_loc][i][j] += b_grad_glo[j][k] * U[k];
					// printf("b_grad_glo[%d][%d] = %.10e\n", j, k, b_grad_glo[j][k]);
				}
				// printf("U[%d] = %.14e\tDisp_grad_glo[%d][%d][%d] = %.14e\n", k, U[k], El_No_loc, i, j, Disp_grad_glo[El_No_loc][i][j]);
			}
			for(j = 0; j < D_MATRIX_SIZE; j++){
				for(k = 0; k < KIEL_SIZE; k++){
					Strain_glo[El_No_loc][i][j] += B_glo[j][k] * U[k];
				}
				// printf("Strain_glo[%d][%d][%d] = %.10e\n", El_No_loc, i, j, Strain_glo[El_No_loc][i][j]);
			}
		}
	}

	return 0;	
}


void Make_auxiliary_mode1(int e, double E, double nu, int DM, double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], double crack_front_coordinates_x, double crack_front_coordinates_y){
	int i, j, N;
	static double Dinv[D_MATRIX_SIZE][D_MATRIX_SIZE];
	double unit_basis_local[DIMENSION] = {0.0};
	double r_tip = sqrt(crack_front_coordinates_x * crack_front_coordinates_x + crack_front_coordinates_y * crack_front_coordinates_y);
	double crack_front_coordinates_x_local, crack_front_coordinates_y_local;
	double data_result_shape_local[DIMENSION] = {0.0};
	double r, rad;
	double mu = E / (2.0 * (1.0 + nu));
	double Pi = 4.0 * atan(1.0);

	Make_gauss_array(0);

	Make_D_Matrix_2D(Dinv, E, nu, DM);

	InverseMatrix_3X3(Dinv);

	// printf("x_crackfront : % 1.8e\n", crack_front_coordinates_x);
	// printf("y_crackfront : % 1.8e\n", crack_front_coordinates_y);

	//x'-y'（き裂先端）座標における単位基底ベクトル
	unit_basis_local[0] = crack_front_coordinates_x / r_tip;
	unit_basis_local[1] = crack_front_coordinates_y / r_tip;
	// printf("unit_basis[0] : % 1.8e\n", unit_basis_local[0]);
	// printf("unit_basis[1] : % 1.8e\n", unit_basis_local[1]);
	crack_front_coordinates_x_local = unit_basis_local[0] * crack_front_coordinates_x + unit_basis_local[1] * crack_front_coordinates_y;
	crack_front_coordinates_y_local = -unit_basis_local[1] * crack_front_coordinates_x + unit_basis_local[0] * crack_front_coordinates_y;
	// printf("x_crackfront_local : % 1.8e\n", crack_front_coordinates_x_local);
    // printf("y_crackfront_local : % 1.8e\n", crack_front_coordinates_y_local);

	for(N = 0; N < GP_2D; N++){
		//各要素でガウス点での物理座標の算出する
		double data_result_shape[2] = {0.0};
		for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[e]]; i++){
			double R_shape_func = Shape_func(i, Gxi[N], e);
			for (j = 0; j < DIMENSION; j++){
				data_result_shape[j] += R_shape_func * X[i][j];
			}
		}
		printf("x_gauss : % 1.8e\n", data_result_shape[0]);
		printf("y_gauss : % 1.8e\n", data_result_shape[1]);

		data_result_shape_local[0] = unit_basis_local[0] * data_result_shape[0] + unit_basis_local[1] * data_result_shape[1];
		data_result_shape_local[1] = -unit_basis_local[1] * data_result_shape[0] + unit_basis_local[0] * data_result_shape[1];
		printf("x_gauss_local : % 1.8e\n", data_result_shape_local[0]);
		printf("y_gauss_local : % 1.8e\n", data_result_shape_local[1]);

		//ガウス点でのx'-y'座標をき裂先端での極座標に変換する
		r = sqrt((data_result_shape_local[0] - crack_front_coordinates_x_local) * (data_result_shape_local[0] - crack_front_coordinates_x_local) + (data_result_shape_local[1] - crack_front_coordinates_y_local) * (data_result_shape_local[1] - crack_front_coordinates_y_local));
		rad = atan2(data_result_shape_local[1] - crack_front_coordinates_y_local, data_result_shape_local[0] - crack_front_coordinates_x_local);
		printf("r_gauss : % 1.8e\n", r);
		printf("rad_gauss : % 1.8e\n", rad);
		printf("degree_gauss : % 1.8e\n", rad * 180.0 / Pi);

		//x'-y'座標での変位勾配を算出する
		Disp_grad_aux_mode1_local[e][N][0] = sqrt(1.0 / (2.0 * Pi * r)) * cos(rad / 2.0) * (1.0 - 2.0 * nu + sin(rad / 2.0) * sin(rad / 2.0) - 2.0 * sin(rad) * sin(rad / 2.0) * cos(rad / 2.0)) / (2.0 * mu);
		Disp_grad_aux_mode1_local[e][N][1] = sqrt(1.0 / (2.0 * Pi * r)) * sin(rad / 2.0) * (1.0 - 2.0 * nu + sin(rad / 2.0) * sin(rad / 2.0) + 2.0 * cos(rad / 2.0) * cos(rad / 2.0) * cos(rad)) / (2.0 * mu);
		Disp_grad_aux_mode1_local[e][N][2] = -sqrt(1.0 / (2.0 * Pi * r)) * sin(rad / 2.0) * (2.0 - 2.0 * nu - cos(rad / 2.0) * cos(rad / 2.0) + 2.0 * sin(rad) * sin(rad / 2.0) * cos(rad / 2.0)) / (2.0 * mu);
		Disp_grad_aux_mode1_local[e][N][3] = sqrt(1.0 / (2.0 * Pi * r)) * cos(rad / 2.0) * (2.0 - 2.0 * nu - cos(rad / 2.0) * cos(rad / 2.0) + 2.0 * sin(rad / 2.0) * sin(rad / 2.0) * cos(rad)) / (2.0 * mu);
		printf("Disp_grad_aux_mode1_local[%d][%d][0] = %1.10e\n", e, N, Disp_grad_aux_mode1_local[e][N][0]);
		printf("Disp_grad_aux_mode1_local[%d][%d][1] = %1.10e\n", e, N, Disp_grad_aux_mode1_local[e][N][1]);
		printf("Disp_grad_aux_mode1_local[%d][%d][2] = %1.10e\n", e, N, Disp_grad_aux_mode1_local[e][N][2]);
		printf("Disp_grad_aux_mode1_local[%d][%d][3] = %1.10e\n", e, N, Disp_grad_aux_mode1_local[e][N][3]);

		//変位勾配をx'-y'座標からx-y座標に変換する
		Disp_grad_aux_mode1[e][N][0] = T[0][0] * T[0][0] * Disp_grad_aux_mode1_local[e][N][0] + T[0][0] * T[1][0] * (Disp_grad_aux_mode1_local[e][N][1] + Disp_grad_aux_mode1_local[e][N][2]) + T[1][0] * T[1][0] * Disp_grad_aux_mode1_local[e][N][3];
		Disp_grad_aux_mode1[e][N][1] = T[0][0] * T[0][1] * Disp_grad_aux_mode1_local[e][N][0] + T[0][0] * T[1][1] * Disp_grad_aux_mode1_local[e][N][1] + T[0][1] * T[1][0] * Disp_grad_aux_mode1_local[e][N][2] + T[1][0] * T[1][1] * Disp_grad_aux_mode1_local[e][N][3];
		Disp_grad_aux_mode1[e][N][2] = T[0][0] * T[0][1] * Disp_grad_aux_mode1_local[e][N][0] + T[0][0] * T[1][1] * Disp_grad_aux_mode1_local[e][N][2] + T[0][1] * T[1][0] * Disp_grad_aux_mode1_local[e][N][1] + T[1][0] * T[1][1] * Disp_grad_aux_mode1_local[e][N][3];
		Disp_grad_aux_mode1[e][N][3] = T[0][1] * T[0][1] * Disp_grad_aux_mode1_local[e][N][0] + T[0][1] * T[1][1] * (Disp_grad_aux_mode1_local[e][N][1] + Disp_grad_aux_mode1_local[e][N][2]) + T[1][1] * T[1][1] * Disp_grad_aux_mode1_local[e][N][3];
		printf("Disp_grad_aux_mode1[%d][%d][0] = %1.10e\n", e, N, Disp_grad_aux_mode1[e][N][0]);
		printf("Disp_grad_aux_mode1[%d][%d][1] = %1.10e\n", e, N, Disp_grad_aux_mode1[e][N][1]);
		printf("Disp_grad_aux_mode1[%d][%d][2] = %1.10e\n", e, N, Disp_grad_aux_mode1[e][N][2]);
		printf("Disp_grad_aux_mode1[%d][%d][3] = %1.10e\n", e, N, Disp_grad_aux_mode1[e][N][3]);

		printf("strain_Disp_grad_aux_mode1[%d][%d][0] = %1.10e\n", e, N, Disp_grad_aux_mode1[e][N][0]);
		printf("strain_Disp_grad_aux_mode1[%d][%d][1] = %1.10e\n", e, N, Disp_grad_aux_mode1[e][N][3]);
		printf("strain_Disp_grad_aux_mode1[%d][%d][2] = %1.10e\n", e, N, Disp_grad_aux_mode1[e][N][1] + Disp_grad_aux_mode1[e][N][2]);

		//x'-y'座標での応力を算出する
		Stress_aux_mode1_local[e][N][0] = sqrt(1.0 / (2.0 * Pi * r)) * cos(rad / 2.0) * (1.0 - sin(rad / 2.0) * sin(3.0 * rad / 2.0));
		Stress_aux_mode1_local[e][N][1] = sqrt(1.0 / (2.0 * Pi * r)) * cos(rad / 2.0) * (1.0 + sin(rad / 2.0) * sin(3.0 * rad / 2.0));
		Stress_aux_mode1_local[e][N][2] = sqrt(1.0 / (2.0 * Pi * r)) * cos(rad / 2.0) * sin(rad / 2.0) * cos(3.0 * rad / 2.0);
		printf("Stress_aux_mode1_local[%d][%d][0] = %1.10e\n", e, N, Stress_aux_mode1_local[e][N][0]);
		printf("Stress_aux_mode1_local[%d][%d][1] = %1.10e\n", e, N, Stress_aux_mode1_local[e][N][1]);
		printf("Stress_aux_mode1_local[%d][%d][2] = %1.10e\n", e, N, Stress_aux_mode1_local[e][N][2]);

		//応力をx'-y'座標からx-y座標に変換する
		Stress_aux_mode1[e][N][0] = T[0][0] * T[0][0] * Stress_aux_mode1_local[e][N][0] + 2.0 * T[0][0] * T[1][0] * Stress_aux_mode1_local[e][N][2] + T[1][0] * T[1][0] * Stress_aux_mode1_local[e][N][1];
		Stress_aux_mode1[e][N][1] = T[0][1] * T[0][1] * Stress_aux_mode1_local[e][N][0] + 2.0 * T[0][1] * T[1][1] * Stress_aux_mode1_local[e][N][2] + T[1][1] * T[1][1] * Stress_aux_mode1_local[e][N][1];
		Stress_aux_mode1[e][N][2] = T[0][0] * T[0][1] * Stress_aux_mode1_local[e][N][0] + (T[0][0] * T[1][1] + T[0][1] * T[1][0]) * Stress_aux_mode1_local[e][N][2] + T[1][0] * T[1][1] * Stress_aux_mode1_local[e][N][1];
		printf("Stress_aux_mode1[%d][%d][0] = %1.10e\n", e, N, Stress_aux_mode1[e][N][0]);
		printf("Stress_aux_mode1[%d][%d][1] = %1.10e\n", e, N, Stress_aux_mode1[e][N][1]);
		printf("Stress_aux_mode1[%d][%d][2] = %1.10e\n", e, N, Stress_aux_mode1[e][N][2]);

		//ひずみを算出する
		for (i = 0; i < D_MATRIX_SIZE; i++){
			for (j = 0; j < D_MATRIX_SIZE; j++){
				Strain_aux_mode1[e][N][i] += Dinv[i][j] * Stress_aux_mode1[e][N][j];
			}
			printf("Strain_aux_mode1[%d][%d][%d] = %1.10e\n", e, N, i, Strain_aux_mode1[e][N][i]);
		}

		//相互ポテンシャルエネルギ密度Wを算出する
		StrainEnergyDensity_aux_mode1[e][N] = Stress_overlay[e][N][0] * Strain_aux_mode1[e][N][0] + Stress_overlay[e][N][1] * Strain_aux_mode1[e][N][1] + Stress_overlay[e][N][2] * Strain_aux_mode1[e][N][2];
		printf("StrainEnergyDensity_aux_mode1[%d][%d] = %1.10e\n", e, N, StrainEnergyDensity_aux_mode1[e][N]);
	}
}


void Make_auxiliary_mode2(int e, double E, double nu, int DM, double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], double crack_front_coordinates_x, double crack_front_coordinates_y)
{
	int i, j, N;
	static double Dinv[D_MATRIX_SIZE][D_MATRIX_SIZE];
	double unit_basis_local[DIMENSION] = {0.0};
	double r_tip = sqrt(crack_front_coordinates_x * crack_front_coordinates_x + crack_front_coordinates_y * crack_front_coordinates_y);
	double crack_front_coordinates_x_local, crack_front_coordinates_y_local;
	double data_result_shape_local[DIMENSION] = {0.0};
	double r, rad;
	double mu = E / (2.0 * (1 + nu));
	double Pi = 4.0 * atan(1.0);

	Make_gauss_array(0);

	Make_D_Matrix_2D(Dinv, E, nu, DM);

	InverseMatrix_3X3(Dinv);

	// printf("x_crackfront : % 1.8e\n", crack_front_coordinates_x);
	// printf("y_crackfront : % 1.8e\n", crack_front_coordinates_y);

	//x'-y'（き裂先端）座標における単位基底ベクトル
	unit_basis_local[0] = crack_front_coordinates_x / r_tip;
	unit_basis_local[1] = crack_front_coordinates_y / r_tip;
	// printf("unit_basis[0] : % 1.8e\n", unit_basis_local[0]);
	// printf("unit_basis[1] : % 1.8e\n", unit_basis_local[1]);
	crack_front_coordinates_x_local = unit_basis_local[0] * crack_front_coordinates_x + unit_basis_local[1] * crack_front_coordinates_y;
	crack_front_coordinates_y_local = -unit_basis_local[1] * crack_front_coordinates_x + unit_basis_local[0] * crack_front_coordinates_y;
	// printf("x_crackfront_local : % 1.8e\n", crack_front_coordinates_x_local);
    // printf("y_crackfront_local : % 1.8e\n", crack_front_coordinates_y_local);

	for(N = 0; N < GP_2D; N++){
		//各要素でガウス点での物理座標の算出する
		double data_result_shape[2] = {0.0};
		for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[e]]; i++){
			double R_shape_func = Shape_func(i, Gxi[N], e);
			for (j = 0; j < DIMENSION; j++){
				X[i][j] = Node_Coordinate[Controlpoint_of_Element[e][i]][j];
				data_result_shape[j] += R_shape_func * X[i][j];
			}
		}
		// printf("x_gauss : % 1.8e\n", data_result_shape[0]);
		// printf("y_gauss : % 1.8e\n", data_result_shape[1]);

		data_result_shape_local[0] = unit_basis_local[0] * data_result_shape[0] + unit_basis_local[1] * data_result_shape[1];
		data_result_shape_local[1] = -unit_basis_local[1] * data_result_shape[0] + unit_basis_local[0] * data_result_shape[1];
		// printf("x_gauss_local : % 1.8e\n", data_result_shape_local[0]);
		// printf("y_gauss_local : % 1.8e\n", data_result_shape_local[1]);

		//ガウス点での物理座標をき裂先端での円筒座標に変換する
		r = sqrt((data_result_shape_local[0] - crack_front_coordinates_x_local) * (data_result_shape_local[0] - crack_front_coordinates_x_local) + (data_result_shape_local[1] - crack_front_coordinates_y_local) * (data_result_shape_local[1] - crack_front_coordinates_y_local));
		rad = atan2(data_result_shape_local[1] - crack_front_coordinates_y_local, data_result_shape_local[0] - crack_front_coordinates_x_local);
		// printf("r_gauss : % 1.8e\n", r);
		// printf("rad_gauss : % 1.8e\n", rad);
		// printf("degree_gauss : % 1.8e\n", rad * 180.0 / Pi);

		//x'-y'座標での変位勾配を算出する
		Disp_grad_aux_mode2_local[e][N][0] = -sqrt(1.0 / (2.0 * Pi * r)) * sin(rad / 2.0) * (2.0 - 2.0 * nu + cos(rad / 2.0) * cos(rad / 2.0) - 2.0 * sin(rad) * sin(rad / 2.0) * cos(rad / 2.0)) / (2.0 * mu);
		Disp_grad_aux_mode2_local[e][N][1] = sqrt(1.0 / (2.0 * Pi * r)) * cos(rad / 2.0) * (2.0 - 2.0 * nu + cos(rad / 2.0) * cos(rad / 2.0) - 2.0 * sin(rad / 2.0) * sin(rad / 2.0) *  cos(rad)) / (2.0 * mu);
		Disp_grad_aux_mode2_local[e][N][2] = -sqrt(1.0 / (2.0 * Pi * r)) * cos(rad / 2.0) * (1.0 - 2.0 * nu - sin(rad / 2.0) * sin(rad / 2.0) + 2.0 * sin(rad) * sin(rad / 2.0) * cos(rad / 2.0)) / (2.0 * mu);
		Disp_grad_aux_mode2_local[e][N][3] = -sqrt(1.0 / (2.0 * Pi * r)) * sin(rad / 2.0) * (1.0 - 2.0 * nu - sin(rad / 2.0) * sin(rad / 2.0) + 2.0 * cos(rad / 2.0) * cos(rad / 2.0) * cos(rad)) / (2.0 * mu);
		// printf("Disp_grad_aux_mode2_local[%d][%d][0] = %1.10e\n", e, N, Disp_grad_aux_mode2_local[e][N][0]);
		// printf("Disp_grad_aux_mode2_local[%d][%d][1] = %1.10e\n", e, N, Disp_grad_aux_mode2_local[e][N][1]);
		// printf("Disp_grad_aux_mode2_local[%d][%d][2] = %1.10e\n", e, N, Disp_grad_aux_mode2_local[e][N][2]);
		// printf("Disp_grad_aux_mode2_local[%d][%d][3] = %1.10e\n", e, N, Disp_grad_aux_mode2_local[e][N][3]);

		//変位勾配をx'-y'座標からx-y座標に変換する
		Disp_grad_aux_mode2[e][N][0] = T[0][0] * T[0][0] * Disp_grad_aux_mode2_local[e][N][0] + T[0][0] * T[1][0] * (Disp_grad_aux_mode2_local[e][N][1] + Disp_grad_aux_mode2_local[e][N][2]) + T[1][0] * T[1][0] * Disp_grad_aux_mode2_local[e][N][3];
		Disp_grad_aux_mode2[e][N][1] = T[0][0] * T[0][1] * Disp_grad_aux_mode2_local[e][N][0] + T[0][0] * T[1][1] * Disp_grad_aux_mode2_local[e][N][1] + T[0][1] * T[1][0] * Disp_grad_aux_mode2_local[e][N][2] + T[1][0] * T[1][1] * Disp_grad_aux_mode2_local[e][N][3];
		Disp_grad_aux_mode2[e][N][2] = T[0][0] * T[0][1] * Disp_grad_aux_mode2_local[e][N][0] + T[0][0] * T[1][1] * Disp_grad_aux_mode2_local[e][N][2] + T[0][1] * T[1][0] * Disp_grad_aux_mode2_local[e][N][1] + T[1][0] * T[1][1] * Disp_grad_aux_mode2_local[e][N][3];
		Disp_grad_aux_mode2[e][N][3] = T[0][1] * T[0][1] * Disp_grad_aux_mode2_local[e][N][0] + T[0][1] * T[1][1] * (Disp_grad_aux_mode2_local[e][N][1] + Disp_grad_aux_mode2_local[e][N][2]) + T[1][1] * T[1][1] * Disp_grad_aux_mode2_local[e][N][3];
		// printf("Disp_grad_aux_mode2[%d][%d][0] = %1.10e\n", e, N, Disp_grad_aux_mode2[e][N][0]);
		// printf("Disp_grad_aux_mode2[%d][%d][1] = %1.10e\n", e, N, Disp_grad_aux_mode2[e][N][1]);
		// printf("Disp_grad_aux_mode2[%d][%d][2] = %1.10e\n", e, N, Disp_grad_aux_mode2[e][N][2]);
		// printf("Disp_grad_aux_mode2[%d][%d][3] = %1.10e\n", e, N, Disp_grad_aux_mode2[e][N][3]);

		// printf("strain_Disp_grad_aux_mode2[%d][%d][0] = %1.10e\n", e, N, Disp_grad_aux_mode2[e][N][0]);
		// printf("strain_Disp_grad_aux_mode2[%d][%d][1] = %1.10e\n", e, N, Disp_grad_aux_mode2[e][N][3]);
		// printf("strain_Disp_grad_aux_mode2[%d][%d][2] = %1.10e\n", e, N, Disp_grad_aux_mode2[e][N][1] + Disp_grad_aux_mode2[e][N][2]);

		//x'-y'座標での応力を算出する
		Stress_aux_mode2_local[e][N][0] = -sqrt(1.0 / (2.0 * Pi * r)) * sin(rad / 2.0) * (2.0 + cos(rad / 2.0) * cos(3.0 * rad / 2.0));
		Stress_aux_mode2_local[e][N][1] = sqrt(1.0 / (2.0 * Pi * r)) * sin(rad / 2.0) * cos(rad / 2.0) * cos(3.0 * rad / 2.0);
		Stress_aux_mode2_local[e][N][2] = sqrt(1.0 / (2.0 * Pi * r)) * cos(rad / 2.0) * (1.0 - sin(rad / 2.0) * sin(3.0 * rad / 2.0));
		// printf("Stress_aux_mode2_local[%d][%d][0] = %1.10e\n", e, N, Stress_aux_mode2_local[e][N][0]);
		// printf("Stress_aux_mode2_local[%d][%d][1] = %1.10e\n", e, N, Stress_aux_mode2_local[e][N][1]);
		// printf("Stress_aux_mode2_local[%d][%d][2] = %1.10e\n", e, N, Stress_aux_mode2_local[e][N][2]);

		//応力をx'-y'座標からx-y座標に変換する
		Stress_aux_mode2[e][N][0] = T[0][0] * T[0][0] * Stress_aux_mode2_local[e][N][0] + 2.0 * T[0][0] * T[1][0] * Stress_aux_mode2_local[e][N][2] + T[1][0] * T[1][0] * Stress_aux_mode2_local[e][N][1];
		Stress_aux_mode2[e][N][1] = T[0][1] * T[0][1] * Stress_aux_mode2_local[e][N][0] + 2.0 * T[0][1] * T[1][1] * Stress_aux_mode2_local[e][N][2] + T[1][1] * T[1][1] * Stress_aux_mode2_local[e][N][1];
		Stress_aux_mode2[e][N][2] = T[0][0] * T[0][1] * Stress_aux_mode2_local[e][N][0] + (T[0][0] * T[1][1] + T[0][1] * T[1][0]) * Stress_aux_mode2_local[e][N][2] + T[1][0] * T[1][1] * Stress_aux_mode2_local[e][N][1];
		// printf("Stress_aux_mode2[%d][%d][0] = %1.10e\n", e, N, Stress_aux_mode2[e][N][0]);
		// printf("Stress_aux_mode2[%d][%d][1] = %1.10e\n", e, N, Stress_aux_mode2[e][N][1]);
		// printf("Stress_aux_mode2[%d][%d][2] = %1.10e\n", e, N, Stress_aux_mode2[e][N][2]);

		//ひずみを算出する
		for (i = 0; i < D_MATRIX_SIZE; i++){
			for (j = 0; j < D_MATRIX_SIZE; j++){
				Strain_aux_mode2[e][N][i] += Dinv[i][j] * Stress_aux_mode2[e][N][j];
			}
			// printf("Strain_aux_mode2[%d][%d][%d] = %1.10e\n", e, N, i, Strain_aux_mode2[e][N][i]);
		}

		//相互ポテンシャルエネルギ密度Wを算出する
		StrainEnergyDensity_aux_mode2[e][N] = Stress_overlay[e][N][0] * Strain_aux_mode2[e][N][0] + Stress_overlay[e][N][1] * Strain_aux_mode2[e][N][1] + Stress_overlay[e][N][2] * Strain_aux_mode2[e][N][2];
		// printf("StrainEnergyDensity_aux_mode2[%d][%d] = %1.10e\n", e, N, StrainEnergyDensity_aux_mode2[e][N]);
	}
}