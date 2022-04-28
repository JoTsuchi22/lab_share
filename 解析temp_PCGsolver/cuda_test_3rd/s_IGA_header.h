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

#endif