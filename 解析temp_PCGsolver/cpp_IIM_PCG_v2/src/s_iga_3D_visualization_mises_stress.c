/***********************************

2020_07_17

IGA
multipch
9 gauss point


For output
mkdir checkAns
	Gauss_stress
	mesh_net



ele_check == を変更
***********************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define ERROR -999
#define ERR 0.0000000000001

#define PI  3.1415926535897932

#define MAX_NO_CCpoint_ON_ELEMENT 64			//分割節点数
#define DIMENSION 3 	 								   //次元数
#define MAX_KIEL_SIZE MAX_NO_CCpoint_ON_ELEMENT * DIMENSION //要素分割マトリックスの大きさ
#define MAX_Ng 10											   //Gauss-Legendreの足す回数
#define MAX_POW_Ng  MAX_Ng * MAX_Ng * MAX_Ng						   //NgのDIMENSION乗の計算
#define Ng 3											   //Gauss-Legendreの足す回数
#define POW_Ng  Ng * Ng * Ng						   //NgのDIMENSION乗の計算
#define D_MATRIX_SIZE 6								//応力歪マトリックスの大きさ（2次元:3 3次元:6）

#define K_DIVISION_LENGE 10		//全体剛性マトリックスのcol&ptrを制作時に分ける節点数
#define EPS 0.0000000001		//連立1次方程式の残差
#define N_STRAIN 6
#define N_STRESS 6
//各種最大配置可能数
#define MAX_N_KNOT 100
#define MAX_N_ELEMENT 40000 //100000
#define MAX_N_NODE 40000 //100000
#define MAX_N_LOAD 1
#define MAX_N_CONSTRAINT 100000
#define MAX_K_WHOLE_SIZE MAX_N_NODE * DIMENSION
#define MAX_NON_ZERO 50000000
#define MAX_N_PATCH 30
#define MAX_N_Controlpoint_in_Patch 30000
#define MAX_N_DISTRIBUTE_FORCE 15
#define MAX_Order_of_Shape_Function 3
#define MAX_N_ORDER 100
#define MAX_element_ndiv 61
#define MAX_N_REFINEMENT MAX_element_ndiv * MAX_element_ndiv

//for s-IGA
#define MAX_N_MESH 3 //ローカルメッシュの個数
#define MAX_N_POINT_OVER MAX_Ng * MAX_Ng * MAX_Ng	//要素重なり判定に用いるローカルメッシュ上1要素内の点数(メモ：つまりガウス点の数?)
#define MAX_N_ELEMENT_OVER	1000	//グローバルメッシュ内の1要素に重なる最大要素数
#define MAX_N_ELEMENT_OVER_POINT	10000	//ローカル要素内の1点に重なるグローバル要素
#define MAX_N_ELEMENT_OVER_ELEMENT	MAX_N_ELEMENT_OVER_POINT * MAX_N_POINT_OVER	//ローカルメッシュ内の1要素に重なる最大要素数


void Make_Gauss_points(int Select_GP); //Select_GP=1ならばインプットした点数で積分，Select_GP=0ならばNg点で積分
int Make_K_EL(int El_No, double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], double K_EL[MAX_KIEL_SIZE][MAX_KIEL_SIZE], double E, double nu, int Total_Element, int Total_Control_Point);

int Make_coupled_K_EL(int El_No_loc, int El_No_glo, double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], double XG[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], double K_EL[MAX_KIEL_SIZE][MAX_KIEL_SIZE], double E, double nu);

void Get_InputData(int I_F, double *E, double *nu, int *Total_Element, int *Total_Control_Point, int *Total_Load, int *No_Patch, int Load_Node_Dir[MAX_N_LOAD][2], double Value_of_Load[MAX_N_LOAD], int *Total_Constraint, int Constraint_Node_Dir[MAX_N_CONSTRAINT][2], double Value_of_Constraint[MAX_N_CONSTRAINT], int *Total_DistributeForce, int argc, char *argv[]);
void Get_displacment(int Total_Control_Point, const char *fileName);

void Check_coupled_Glo_Loc_element(double element_loc[DIMENSION], int mesh_n_over, int mesh_n_org);

void Make_Loc_Glo();

//全体剛性マトリックス
int Make_Index_Dof(int Total_Control_Point, int Total_Constraint, int Constraint_Node_Dir[MAX_N_CONSTRAINT][2]);
void Make_K_Whole_Ptr_Col(int Total_Element, int Total_Control_Point, int K_Whole_Size);
void Make_K_Whole_Val(double E, double nu, int Total_Element, int K_Whole_Size, int Total_Control_Point);
//連立1次方程式
void Make_F_Vec(int Total_Load, int Load_Node_Dir[MAX_N_LOAD][2], double Value_of_Load[MAX_N_LOAD], int K_Whole_Size);
void Make_F_Vec_disp_const(int Mesh_No, int Total_Constraint, int Constraint_Node_Dir[MAX_N_CONSTRAINT][2], double Value_of_Constraint[MAX_N_CONSTRAINT], int Total_Element, double E, double nu, int Total_Control_Point);
void mat_vec_crs(double vec_result[], double vec[], const int ndof);
double inner_product(int ndof, double vec1[], double vec2[]);
int check_conv_CG(int ndof, double alphak, double pp[], double eps, int itr);
void Diag_Scaling_CG_pre(int ndof, int flag_operation);
void CG_Solver(int ndof, int max_itr, double eps, int flag_ini_val);
//各種値
void Make_Strain(double E, double nu, int Total_Element, int El_No, int Total_Control_Point);
void Make_Stress(double E, double nu, int Total_Element/*, int DM*/);
void Make_ReactionForce(int Total_Element, int Total_Control_Point, int El_No);
//void Make_Parameter_z(int Total_Element, double E, double nu, /*int DM*/);
//分布荷重
//void Force_dis(int Distriction_Force[DIMENSION][3], double Val_Distribute_Force[DIMENSION], double X[MAX_NO_CCpoint_ON_ELEMENT][DIMENSION], double Fe[DIMENSION]);

//NURBSの計算
void element_coordinate(int Total_Element, int Total_Control_Point);
void calculate_Controlpoint_using_NURBS(double element[DIMENSION], double element_glo[DIMENSION], int Total_Element, int Total_Control_Point, double E, double nu);
//void calculate_extendmesh_using_NURBS(double element_emsh[DIMENSION], int Total_Element, int Total_Control_Point);
void Gausspoint_coordinate(int Total_Element, int Total_Control_Point);

/*///J積分
int Make_B_x_Matrix_Quad_4(double B_x[DIMENSION][KIEL_SIZE], double Local_coord[DIMENSION], double X[No_Control_point_ON_ELEMENT][DIMENSION], double *J );
void Make_Strain_x_Quad_4(double E, double nu, int Total_Element);
void Make_EMT(double E, double nu, int Total_Element);*/

/* Distributed Load */

int SerchForElement(int mesh_n, int iPatch, int Total_Element, int iX, int iY, int iZ);

void Setting_Dist_Load(int mesh_n, int Total_Control_Point, int iPatch, int Total_Element, int iCoord, int jCoord, double val_Coord, double iRange_Coord[2], double jRange_Coord[2], int type_load, double Coeff_Dist_Load_i[3], double Coeff_Dist_Load_j[3]);

void Add_Equivalent_Nodal_Forec_to_F_Vec(int Total_Control_Point);


//static int DIMENSION;
static double E, nu;
static int KIEL_SIZE; //要素分割マトリックスの大きさ
static int Load_Node_Dir[MAX_N_LOAD][2];
static double Value_of_Load[MAX_N_LOAD];
static int Constraint_Node_Dir[MAX_N_CONSTRAINT][2];
static double Value_of_Constraint[MAX_N_CONSTRAINT];


//static int Controlpoint_of_Element_nomarge[MAX_N_ELEMENT][MAX_NO_CCpoint_ON_ELEMENT];
static int Controlpoint_of_Element[MAX_N_ELEMENT][MAX_NO_CCpoint_ON_ELEMENT];
static double Node_Coordinate[MAX_N_NODE][DIMENSION + 1];
static double Equivalent_Nodal_Force[MAX_N_NODE][DIMENSION]; // Equivalent nodal forces arising from the distributed load
static int K_Whole_Ptr[MAX_K_WHOLE_SIZE + 1], K_Whole_Col[MAX_NON_ZERO];
static double K_Whole_Val[MAX_NON_ZERO];
static int Index_Dof[MAX_K_WHOLE_SIZE];
static int INC[MAX_N_PATCH][MAX_N_ELEMENT][DIMENSION];
//static int Adress_Controlpoint[MAX_N_PATCH][MAX_N_Controlpoint_in_Patch_one_Direction][MAX_N_Controlpoint_in_Patch_one_Direction][MAX_N_Controlpoint_in_Patch_one_Direction]; //INCの配列をいじったもの
static int Order[MAX_N_PATCH][DIMENSION];
static int No_knot[MAX_N_PATCH][DIMENSION];
static int No_Control_point[MAX_N_PATCH][DIMENSION];
static double element_coordinate_Nopoint[MAX_N_ELEMENT][DIMENSION];
static double Gausspoint_coordinates[MAX_N_ELEMENT][POW_Ng][DIMENSION];
//static int same_point[100];
static int same_point_in_Element[MAX_N_NODE];
static int Patch_controlpoint[MAX_N_PATCH][MAX_N_Controlpoint_in_Patch]; //バッチとコントロールポイント番号の要素コネクティビティ
static int Element_patch[MAX_N_ELEMENT];								 //要素がどのパッチに属しているか示す配列(要素番号は1つのモデルで通し番号)
static int No_Controlpoint_in_patch[MAX_N_PATCH];
static int No_Control_point_ON_ELEMENT[MAX_N_PATCH]; //修正，もとは10000

static int Node_To_Node[K_DIVISION_LENGE][10000], Total_Control_Point_To_Node[K_DIVISION_LENGE]; //ある節点に関係する節点番号


static double sol_vec[MAX_K_WHOLE_SIZE];
static double rhs_vec[MAX_K_WHOLE_SIZE];
static double diag_scaling[MAX_K_WHOLE_SIZE];

static double Shape[DIMENSION][MAX_N_NODE][MAX_Order_of_Shape_Function];  //MAX_N_NODE → MAX_N_KNOTでは？
static double shape_func[MAX_N_NODE];
static double dShape_func1[MAX_N_NODE];
static double dShape_func2[MAX_N_NODE];
static double dShape_func3[MAX_N_NODE];
static double dShape[DIMENSION][MAX_N_NODE];
static double Position_Knots[MAX_N_PATCH][DIMENSION][MAX_N_KNOT];
static double Position_Data_param[DIMENSION];

static double Displacement[MAX_K_WHOLE_SIZE];
static double Strain[MAX_N_ELEMENT][POW_Ng][N_STRAIN];
static double Stress[MAX_N_ELEMENT][POW_Ng][N_STRESS];
static double ReactionForce[MAX_K_WHOLE_SIZE];

static double difference[MAX_N_PATCH][MAX_N_KNOT][DIMENSION];		  /*隣り合うノットベクトルの差*/
static int ENC[MAX_N_PATCH][MAX_N_ELEMENT][DIMENSION];				  /*ENC[パッチ][全ての要素][0,1,2]=x,y,z方向の何番目の要素か*/
static int real_Total_Element;							  /*ゼロエレメントを除いた要素数*/
static int real_element[MAX_N_ELEMENT];					  /*ゼロエレメントではない要素の番号*/
static int Total_element_all_ID[MAX_N_ELEMENT];			  /*ゼロエレメントではない要素＝１、ゼロエレメント＝０*/
static int line_No_Total_element[MAX_N_PATCH][DIMENSION]; /*ゼロエレメントを含むすべての要素列の数*/
static int line_No_real_element[MAX_N_PATCH][DIMENSION];  /*ゼロエレメントではない要素列の数*/
static int real_element_line[MAX_N_PATCH][MAX_N_ELEMENT][DIMENSION];   /*ゼロエレメントではない要素列*/
static double Check_Volume;
static double J_test;

//static int No_points_for_colored_points; /*zarusobaで点に色付ける時の全ての点の数*/



//static int No_points_for_new_zarusoba; /*zarusobaで点に色付ける時の全ての点の数*/

//static double data_result_shape_x_for_new_zarusoba[10000];
//static double data_result_shape_y_for_new_zarusoba[10000];
//static double data_result_shape_z_for_new_zarusoba[10000];//修正2020/0717

//static double data_result_disp_x_for_new_zarusoba[10000];
//static double data_result_disp_y_for_new_zarusoba[10000];
//static double data_result_disp_z_for_new_zarusoba[10000];//修正2020/0717

//for s-IGA
static int No_Input_File; //インプットのファイル数

static int Total_Patch_on_mesh[MAX_N_MESH];	//各メッシュ上のパッチ数
static int Total_Patch_to_mesh[MAX_N_MESH + 1];	//メッシュ[]までのパッチ数（メッシュ[]内のパッチ数は含まない）

static int Total_Control_Point_on_mesh[MAX_N_MESH];	//各メッシュ上のコントロールポイント数
static int Total_Control_Point_to_mesh[MAX_N_MESH + 1];	//メッシュ[]までのコントロールポイント数（メッシュ[]内のコントロールポイント数は含まない）


static int Total_Element_on_mesh[MAX_N_MESH];
static int Total_Element_to_mesh[MAX_N_MESH + 1];

static int Total_Constraint_on_mesh[MAX_N_MESH];
static int Total_Constraint_to_mesh[MAX_N_MESH + 1];
static int Constraint_Node_Dir_on_mesh[MAX_N_MESH][MAX_N_CONSTRAINT][2];
static double Value_of_Constraint_on_mesh[MAX_N_MESH][MAX_N_CONSTRAINT];

static int Total_Load_on_mesh[MAX_N_MESH];
static int Total_Load_to_mesh[MAX_N_MESH + 1];

static int Total_DistributeForce_on_mesh[MAX_N_MESH];
static int Total_DistributeForce_to_mesh[MAX_N_MESH + 1];

static int Element_mesh[MAX_N_ELEMENT]; //要素がどのメッシュ内にあるかを示す配列
static int Patch_mesh[MAX_N_PATCH]; //パッチがどのメッシュ内にあるかを示す配列

static int real_Total_Element_on_mesh[MAX_N_MESH];
static int real_Total_Element_to_mesh[MAX_N_MESH + 1];
static int real_El_No_on_mesh[MAX_N_MESH][MAX_N_ELEMENT];	//メッシュ内でのコントロールポイント配列(ゼロエレメントを含まない)
static int El_No_on_mesh[MAX_N_MESH][MAX_N_ELEMENT];   //メッシュ内でのコントロールポイント配列(ゼロエレメントを含む)

static double Control_Coord[DIMENSION][MAX_N_NODE];
static double Control_Weight[MAX_N_NODE];

static int temp_element_n[MAX_N_ELEMENT_OVER_POINT];
static int element_n_point[MAX_N_ELEMENT_OVER_ELEMENT]; //ある点が重なっているグローバル要素
static int NNLOVER[MAX_N_ELEMENT];	//NNLOVER:要素eに重なる要素の総数
static int NELOVER[MAX_N_ELEMENT][MAX_N_ELEMENT_OVER];	//NELOVER:要素MAX_N_ELEMENT_OVERに重なるMAX_N_ELEMENT番目の要素番号
static int Check_BDBJ_flag[MAX_N_ELEMENT];
static int Total_BDBJ_flag;
static int Same_BDBJ_flag[MAX_POW_Ng];

//for ref
static int temp_element_n_ref[MAX_N_ELEMENT_OVER_POINT];
static int element_n_point_for_ref;

static double Stress_ref[N_STRESS];
static double Stress_ref_glo[N_STRESS];
static double Stress_ref_loc[N_STRESS];

static double Strain_ref[N_STRAIN];
static double Strain_ref_glo[N_STRAIN];
static double Strain_ref_loc[N_STRAIN];

static int element_ndiv_loc_glo[MAX_N_MESH];


static int No_GP;
static int GaussPt_1dir;
static int GaussPt_3D;

static double Gxi[MAX_POW_Ng][DIMENSION];
static double w[MAX_POW_Ng];
static double G_1d[MAX_Ng];
static double w_1d[MAX_Ng];




FILE *fp;

int main(int argc, char *argv[])
{
	clock_t start, end1, end2, end3;

	int i, j, k;
	int re;
	int Total_Element;
	int Total_Control_Point;
	int No_Patch = 0;
	int Total_net = 0;
	static int Total_Load = 0;
	static int Total_Constraint = 0;
	static int Total_DistributeForce = 0;
	int K_Whole_Size = 0;
	int El_No = 0;
	static double element[DIMENSION];
	static double element_glo[DIMENSION];
	//static double element_emsh[DIMENSION];
	static int max_itr;


	//for s-IGA
	static double element_loc[DIMENSION];
	int I_F;

	if (argc <= 1)
	{
		printf("Argument is missing\n");
	}
    if (argc == 2)  /*通常IGA：input file 1つ*/
    {
        printf("IGA carried out.(No local mesh)\n");
    }
    if (3 <= argc)  /*s-IGA：input file 複数*/
    {
        printf("s-IGA carried out.(%d local meshes)\n",argc - 2);
    }

	No_Input_File = argc - 2;

	printf("No_Input_File = %d\n", No_Input_File);

	start = clock();

	/////////////////////インプットデータの読み込み/////////////////////
    for (I_F = 0; I_F < No_Input_File; I_F++)
    {
		Total_Control_Point_on_mesh[I_F] = 0;
	}

    for (I_F = 0; I_F < No_Input_File; I_F++)
    {
		Get_InputData(I_F, &E, &nu, &Total_Element, &Total_Control_Point, &Total_Load, &No_Patch, Load_Node_Dir, Value_of_Load, &Total_Constraint, Constraint_Node_Dir, Value_of_Constraint, &Total_DistributeForce, argc, argv);
		if(I_F == 0 && No_Input_File > 1)
        {
            printf("Finish Get_InputData(Global mesh:%s)\n",argv[1]);
        }
        if(I_F > 0)
        {
            printf("Finish Get_inputData(Local mesh No.[%d]:%s)\n",I_F,argv[I_F + 1]);
        }
	}

	Get_displacment(Total_Control_Point, argv[No_Input_File + 1]);
	printf("Finish Get_displacment\n");


	/////////////////////重なる要素を見つける/////////////////////
	for (i = 1; i < No_Input_File; i++)
	{
		printf("mesh_n_org;0\tMesh_n_over;%d\n",i);

		Check_coupled_Glo_Loc_element(element_loc, i, 0);
		Make_Loc_Glo();

	}



	/////////////////////全体剛性マトリックスの制作////////////////////////////
	//K_Whole_Size = Make_Index_Dof(Total_Control_Point_to_mesh[No_Input_File], Total_Constraint_to_mesh[No_Input_File], Constraint_Node_Dir);
	K_Whole_Size = 1;
	printf("K_Whole_Size = %d\n",K_Whole_Size);


	//Make_K_Whole_Ptr_Col(Total_Element_to_mesh[No_Input_File], Total_Control_Point_to_mesh[No_Input_File], K_Whole_Size);
	//Make_K_Whole_Val(E, nu, real_Total_Element_to_mesh[No_Input_File], K_Whole_Size, Total_Control_Point_to_mesh[No_Input_File] /*,real_element[MAX_N_ELEMENT]*/);
	printf("Finish Make_K_Whole\n");



	///////////////連立一次方程式/////////////////////////////////////////
	for (i = 0; i < Total_Load_to_mesh[No_Input_File]; i++)
	{
		printf("Value_of_Load[%d] = %11.10e\n",i , Value_of_Load[i]);
	}

	//Make_F_Vec(Total_Load_to_mesh[No_Input_File], Load_Node_Dir, Value_of_Load, K_Whole_Size);
	//Make_F_Vec_disp_const(I_F, Total_Constraint_to_mesh[No_Input_File], Constraint_Node_Dir, Value_of_Constraint, Total_Element_to_mesh[No_Input_File], E, nu, Total_Control_Point_to_mesh[No_Input_File]);
	//Add_Equivalent_Nodal_Forec_to_F_Vec(Total_Control_Point_to_mesh[No_Input_File]);

	max_itr = K_Whole_Size;

	//printf("K_Whole_Size:%d\n", K_Whole_Size);

	//Diag_Scaling_CG_pre(K_Whole_Size, 0);
	//CG_Solver(K_Whole_Size, max_itr, EPS, 0);
	//Diag_Scaling_CG_pre(K_Whole_Size, 1);
	//printf("Finish CG_Solver\n");


	printf("Finish Make_Displacement\n");
	end1 = clock();
	printf("変位を求めるまでに%.2f秒かかりました\n",(double)(end1-start)/CLOCKS_PER_SEC);
	//Make_Strain(E, nu, Total_Element_to_mesh[No_Input_File], El_No, Total_Control_Point_to_mesh[No_Input_File]);
	printf("Finish Make_Strain\n");
	//Make_Stress(E, nu, Total_Element_to_mesh[No_Input_File]);
	printf("Finish Make_Stress\n");
	//Make_ReactionForce(Total_Element_to_mesh[No_Input_File], Total_Control_Point_to_mesh[No_Input_File], El_No);
	printf("Finish Make_ReactionForce\n");
	end2 = clock();
	printf("変位，ひずみ，応力，反力を求めるのに%.2f秒かかりました\n",(double)(end2-start)/CLOCKS_PER_SEC);
	puts("sol_vec");



	/////////////////////////////////////////////////////////////////////////////////////OUTPUT FILE/////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////checkAns.txt/////////////////////////////////////////////////////////
	/*fp = fopen("checkAns/checkAns.txt", "w");
	for (i = 0; i < K_Whole_Size; i++)
	{
		fprintf(fp, "%le", sol_vec[i]);
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n\n\nDisplacement\n");
	for (j = 0; j < Total_Control_Point_to_mesh[No_Input_File]; j++)
	{
		fprintf(fp, "%d\t", j);
		for (i = 0; i < DIMENSION; i++)
			fprintf(fp, "%.13e\t", Displacement[j * DIMENSION + i]);
		fprintf(fp, "\n");
	}
	fprintf(fp, "\n\n\nStrain\n");
	for (re = 0; re < real_Total_Element_to_mesh[No_Input_File]; re++)
	{
		i = real_element[re];
		for (k = 0; k < POW_Ng; k++)
		{
			fprintf(fp, "%d\t%d\t", i, k);
			for (j = 0; j < N_STRAIN; j++)
				fprintf(fp, "%.13e\t", Strain[i][k][j]);
			fprintf(fp, "\n");
		}
	}
	fprintf(fp, "\n\n\n\n\n\n\n\n\n\n\nStress\n");
	for (re = 0; re < real_Total_Element_to_mesh[No_Input_File]; re++)
	{
		i = real_element[re];
		for (k = 0; k < POW_Ng; k++)
		{
			fprintf(fp, "%d\t%d\t", i, k);
			for (j = 0; j < N_STRESS; j++)
				fprintf(fp, "%.13e\t", Stress[i][k][j]);
			fprintf(fp, "\n");
		}
	}
	fprintf(fp, "\n\n\nReaction Force\n");
	for (j = 0; j < Total_Control_Point_to_mesh[No_Input_File]; j++)
	{
		for (i = 0; i < DIMENSION; i++)
			fprintf(fp, "%.13e\t", ReactionForce[j * DIMENSION + i]);
		fprintf(fp, "\n");
	}
	fclose(fp);*/
	/////////////////////////////////////////////////////////checkAns.txt/////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////Displacement.dat/////////////////////////////////////////////////////////
	//fp = fopen("Displacement.dat", "w");
	////fprintf(fp, "label=Displacement\n");
	//fprintf(fp, "%d\n", Total_Control_Point_to_mesh[No_Input_File]);
	//fprintf(fp, "\n");
	//for (I_F = 0; I_F < No_Input_File; I_F++)
	//{
	//	for (j = Total_Control_Point_to_mesh[I_F]; j < Total_Control_Point_to_mesh[I_F + 1]; j++)
	//	{
	//		fprintf(fp, "%d\t\t%.15e\t%.15e\t%.15e", j, Displacement[j * DIMENSION], Displacement[j * DIMENSION + 1], Displacement[j * DIMENSION + 2]);
	//		fprintf(fp, "\n");
	//	}
	//	fprintf(fp, "\n\n\n");
	//}
	//fclose(fp);
	/////////////////////////////////////////////////////////Displacement.dat/////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////Strain@IntegrationPoint.dat/////////////////////////////////////////////////////////
	/*fp = fopen("Strain@IntegrationPoint.dat", "w");
	fprintf(fp, "label=Strain@IntegrationPoint\n");
	fprintf(fp, "num_items=%d\n", real_Total_Element_to_mesh[No_Input_File]);
	fprintf(fp, "\n");
	for (I_F = 0; I_F < No_Input_File; I_F++)
	{
		for (re = real_Total_Element_to_mesh[I_F]; re < real_Total_Element_to_mesh[I_F + 1]; re++)
		{
			i = real_element[re];
			fprintf(fp, "%d\t:", i);
			for (k = 0; k < POW_Ng; k++)
			{
				fprintf(fp, "%13e\t%13e\t%13e\t%13e\t%13e\t%13e\t",Strain[i][k][0], Strain[i][k][1], Strain[i][k][2], Strain[i][k][3], Strain[i][k][4], Strain[i][k][5]); //0:xx 1:yy 2:zz 3:xy 4:yz 5:zx
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n\n\n");
	}
	fclose(fp);*/
	/////////////////////////////////////////////////////////Strain@IntegrationPoint.dat/////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////Stress@IntegrationPoint.dat/////////////////////////////////////////////////////////
	/*fp = fopen("Stress@IntegrationPoint.dat", "w");
	fprintf(fp, "label=Stress@IntegrationPoint\n");
	fprintf(fp, "num_items=%d\n", real_Total_Element_to_mesh[No_Input_File]);
	fprintf(fp, "\n");
	for (I_F = 0; I_F < No_Input_File; I_F++)
	{
		for (re = real_Total_Element_to_mesh[I_F]; re < real_Total_Element_to_mesh[I_F + 1]; re++)
		{
			i = real_element[re];
			fprintf(fp, "%d\t:", i);
			for (k = 0; k < POW_Ng; k++)
			{
				fprintf(fp, "%13e\t%13e\t%13e\t%13e\t%13e\t%13e\t",Stress[i][k][0], Stress[i][k][1], Stress[i][k][2], Stress[i][k][3], Stress[i][k][4], Stress[i][k][5]); //0:xx 1:yy 2:zz 3:xy 4:yz 5:zx
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n\n\n");
	}
	fclose(fp);*/
	/////////////////////////////////////////////////////////Stress@IntegrationPoint.dat/////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////ReactionForce.dat/////////////////////////////////////////////////////////
	/*fp = fopen("ReactionForce.dat", "w");
	fprintf(fp, "label=ReactionForce\n");
	fprintf(fp, "num_items=%d\n", Total_Control_Point);
	fprintf(fp, "\n");
	for (I_F = 0; I_F < No_Input_File; I_F++)
	{
		for (j = Total_Control_Point_to_mesh[I_F]; j < Total_Control_Point_to_mesh[I_F + 1]; j++)
		{
			fprintf(fp, "%d\t:%13e\t%13e\t%13e", j, ReactionForce[j * DIMENSION], ReactionForce[j * DIMENSION + 1], ReactionForce[j * DIMENSION + 2]);
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n\n\n");
	}
	fclose(fp);*/
	/////////////////////////////////////////////////////////ReactionForce.dat/////////////////////////////////////////////////////////

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////一旦保留、可視化関連////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//	fp = fopen("mesh.msh", "w");
//
//	fprintf(fp, "%d\n", real_Total_Element /*Total_Element*/);
//
//	//for( i = 0; i < Total_Element; i ++ ){
//	for (re = 0; re < real_Total_Element; re++)
//	{
//		i = real_element[re];
//		fprintf(fp, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", Controlpoint_of_Element[i][26], Controlpoint_of_Element[i][24], //Controlpoint_of_Element[i][18], Controlpoint_of_Element[i][20], Controlpoint_of_Element[i][8], Controlpoint_of_Element[i][6], Controlpoint_of_Element[i][0], Controlpoint_of_Element[i][2], //Controlpoint_of_Element[i][17], Controlpoint_of_Element[i][15], Controlpoint_of_Element[i][9], Controlpoint_of_Element[i][11], Controlpoint_of_Element[i][25], Controlpoint_of_Element[i][21], //Controlpoint_of_Element[i][19], Controlpoint_of_Element[i][23], Controlpoint_of_Element[i][22], Controlpoint_of_Element[i][7], Controlpoint_of_Element[i][3], Controlpoint_of_Element[i][1], //Controlpoint_of_Element[i][5], Controlpoint_of_Element[i][4], Controlpoint_of_Element[i][16], Controlpoint_of_Element[i][12], Controlpoint_of_Element[i][10], Controlpoint_of_Element[i][14], //Controlpoint_of_Element[i][13]);
//	}
//
//	fprintf(fp, "%d\n", Total_Control_Point);
//	for (i = 0; i < Total_Control_Point; i++)
//	{
//		fprintf(fp, "%lf\t%lf\t%lf\n", Node_Coordinate[i][0], Node_Coordinate[i][1], Node_Coordinate[i][2]);
//	}
//	fprintf(fp, "\n");
//
//	fclose(fp);

//	fp = fopen("mesh_net/control_net.msh", "w");
//	for (l = 0; l < No_Patch; l++)
//	{
//		Total_net += (No_Control_point[l][0] - 1) * (No_Control_point[l][1] - 1) * (No_Control_point[l][2] - 1);
//		//修正
//	}
//
//	fprintf(fp, "%d\n", Total_net);
//	for (l = 0; l < No_Patch; l++)
//	{
//		for (k = 0; k < No_Control_point[l][2] - 1; k++)
//		{
//			for (j = 0; j < No_Control_point[l][1] - 1; j++)
//			{
//				for (i = 0; i < No_Control_point[l][0] - 1; i++)
//				{
//					fprintf(fp, "%d %d %d %d %d %d %d %d\n", Adress_Controlpoint[l][i][j][k], Adress_Controlpoint[l][i + 1][j][k], Adress_Controlpoint[l][i + 1][j + 1][k], Adress_Controlpoint[l][i][j + 1][k], //Adress_Controlpoint[l][i][j][k + 1], Adress_Controlpoint[l][i + 1][j][k + 1], Adress_Controlpoint[l][i + 1][j + 1][k + 1], Adress_Controlpoint[l][i][j + 1][k + 1]);	//修正したが必要,12コ必要なは//ず
//				}
//			}
//		}
//	}
//
//	fprintf(fp, "%d\n", Total_Control_Point);
//
//	for (i = 0; i < Total_Control_Point; i++)
//	{
//		fprintf(fp, "%lf %lf %lf\n", Node_Coordinate[i][0], Node_Coordinate[i][1], Node_Coordinate[i][2]);
//	}
//	fclose(fp);

//	fp = fopen("mesh_net/element.msh", "w");
//
//	fprintf(fp, "%d\n", Total_Element);
//	element_coordinate(Total_Element, Total_Control_Point);
//	for (i = 0; i < Total_Element * 27; i += 27) 	//修正が必要 おそらく*9じゃなくて*27?
//	{
//		fprintf(fp, "%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", same_point_in_Element[i], same_point_in_Element[i + 1], same_point_in_Element[i + 2], same_point_in_Element[i //+ 3], same_point_in_Element[i + 4], same_point_in_Element[i + 5], same_point_in_Element[i + 6], same_point_in_Element[i + 7], same_point_in_Element[i + 8], same_point_in_Element[i + 9], same_point_in_Element//[i + 10], same_point_in_Element[i + 11], same_point_in_Element[i + 12], same_point_in_Element[i + 13], same_point_in_Element[i + 14], same_point_in_Element[i + 15], same_point_in_Element[i + 16], //same_point_in_Element[i + 17], same_point_in_Element[i + 18], same_point_in_Element[i + 19], same_point_in_Element[i + 20], same_point_in_Element[i + 21], same_point_in_Element[i + 22], same_point_in_Element//[i + 23], same_point_in_Element[i + 24], same_point_in_Element[i + 25], same_point_in_Element[i + 26]);	//修正が必要
//	}
//	fprintf(fp, "%d\n", Total_Element * 27);	//修正が必要 おそらく*9じゃなくて*27?
//	for (i = 0; i < Total_Element * 27; i++)
//	{
//		fprintf(fp, "%lf %lf %lf\n", element_coordinate_Nopoint[i][0], element_coordinate_Nopoint[i][1], element_coordinate_Nopoint[i][2]);
//	}
//	fclose(fp);
//
//	fp = fopen("NURBS/NURBS_disp.dat", "w");
//
	printf("やってる1\n");
	calculate_Controlpoint_using_NURBS(element, element_glo, Total_Element, Total_Control_Point, E, nu);
	printf("やってる2\n");
//
//	fclose(fp);
//
//	fp = fopen("NURBS/control_point.dat", "w");
//	for (j = 0; j < Total_Control_Point; j++)
//	{
//		fprintf(fp, "%d:  %le  %le	%le  %le  %le %le ", j, Node_Coordinate[j][0], Node_Coordinate[j][1], Node_Coordinate[j][2], Displacement[j * DIMENSION + 0], Displacement[j * DIMENSION + 1], Displacement[j * //DIMENSION + 2]);
//		fprintf(fp, "\n");
//	}
//	fclose(fp);

//	fp = fopen("new_zarusoba/control_point.dat", "w");
//	fprintf(fp, "%d\n", Total_Control_Point);
//	for (j = 0; j < Total_Control_Point; j++)
//	{
//		fprintf(fp, "%le  %le  %le\n", Node_Coordinate[j][0], Node_Coordinate[j][1], Node_Coordinate[j][2]);
//	}
//	fclose(fp);
//
//	calculate_extendmesh_using_NURBS(element_emsh, Total_Element, Total_Control_Point);
//
//	fp = fopen("new_zarusoba/extended_mesh.emsh", "w");
//	fprintf(fp, "%d\n", real_Total_Element * 10); //修正が必要かも
//
//	q = 0;
//	r = 0;
//	for (i = Total_Control_Point; i < (Total_Control_Point + No_points_for_new_zarusoba - 2); i = i + 2)
//	{
//		if (q != 10 + 11 * r)
//		{
//			fprintf(fp, "%d %d %d %d\n", i, i + 1, i + 3, i + 2);
//			//printf("q:%d\n", q);
//		}
//		if (q == 10 + 11 * r)
//		{
//			r++;
//		}
//		//printf("r:%d\n", r);
//		q++;
//	}
//
//	for (l = 0; l < No_Patch; l++)
//	{
//		Total_net += (No_Control_point[l][0] - 1) * (No_Control_point[l][1] - 1) * (No_Control_point[l][2] - 1);
//		//修正した
//	}
//
//	fprintf(fp, "%d\n", Total_net * 2); //どう修正すればいいか分からん_zarusobaを確認する必要がある
//	for (l = 0; l < No_Patch; l++)
//	{
//		for (k = 0; k < No_Control_point[l][2] - 1; k++)
//		{
//			for (j = 0; j < No_Control_point[l][1] - 1; j++)
//			{
//				for (i = 0; i < No_Control_point[l][0] - 1; i++)
//				{
//					fprintf(fp, "%d %d\n", Adress_Controlpoint[l][i    ][j    ][k    ], Adress_Controlpoint[l][i + 1][j    ][k    ]); //おそらくここが12行になる
//					fprintf(fp, "%d %d\n", Adress_Controlpoint[l][i + 1][j    ][k    ], Adress_Controlpoint[l][i + 1][j + 1][k    ]);
//					fprintf(fp, "%d %d\n", Adress_Controlpoint[l][i + 1][j + 1][k    ], Adress_Controlpoint[l][i    ][j + 1][k    ]);
//					fprintf(fp, "%d %d\n", Adress_Controlpoint[l][i    ][j + 1][k    ], Adress_Controlpoint[l][i    ][j    ][k    ]);
//				}
//			}
//		}
//	}
//
//	fprintf(fp, "%d\n", Total_Control_Point + No_points_for_new_zarusoba);
//
//	for (i = 0; i < Total_Control_Point; i++)
//	{
//		fprintf(fp, "%lf\t%lf\t%lf\n", Node_Coordinate[i][0], Node_Coordinate[i][1], Node_Coordinate[i][2]);
//	}
//
//	for (i = Total_Control_Point; i < Total_Control_Point + No_points_for_new_zarusoba; i++)
//	{
//		fprintf(fp, "%lf\t%lf\t%lf\n", data_result_shape_x_for_new_zarusoba[i], data_result_shape_y_for_new_zarusoba[i], data_result_shape_z_for_new_zarusoba[i]);
//	}
//
//	fclose(fp);
	/*2020/07/30
	fp = fopen("new_zarusoba/radial_displacement.dat", "w");	//修正が必要 下9行 そもそも三次元で必要か?
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
	Make_ReactionForce_Quad_4( Total_Element, Total_Control_Point, El_No );
	printf("Finish Make_ReactionForce\n");
	puts("sol_vec");
	Make_Parameter_z( Total_Element, E, nu, DM);
	printf("Finish Make_Parameter_z\n");
*/

//	fp = fopen("colored_point/NURBS_points.txt", "w");
//
//	fprintf(fp, "%d\n", No_points_for_colored_points);
//	for (p = 0; p < No_points_for_colored_points; p++)
//	{
//		fprintf(fp, "%-.13lf  %-.13lf  %-.13lf\n", data_result_shape_x[p], data_result_shape_y[p], data_result_shape_z[p]);
//	}
//
//	fclose(fp);
//
//	fp = fopen("colored_point/NURBS_disp_x.dat", "w");
//	fprintf(fp, "label=Displacement\nnum_items=%d\n\n", No_points_for_colored_points);
//	for (p = 0; p < No_points_for_colored_points; p++)
//	{
//		fprintf(fp, "%d:%-.13le\n", p, data_result_disp_x[p]);
//	}
//	fclose(fp);
//
//	fp = fopen("colored_point/NURBS_disp_y.dat", "w");
//	fprintf(fp, "label=Displacement\nnum_items=%d\n\n", No_points_for_colored_points);
//	for (p = 0; p < No_points_for_colored_points; p++)
//	{
//		fprintf(fp, "%d:%-.13le\n", p, data_result_disp_y[p]);
//	}
//	fclose(fp);
//
//	fp = fopen("colored_point/NURBS_disp_z.dat", "w");
//	fprintf(fp, "label=Displacement\nnum_items=%d\n\n", No_points_for_colored_points);
//	for (p = 0; p < No_points_for_colored_points; p++)
//	{
//		fprintf(fp, "%d:%-.13le\n", p, data_result_disp_z[p]);
//	}
//	fclose(fp);
//
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////一旦保留、可視化関連////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//	fp = fopen("colored_point/NURBS_disp_radius.dat", "w");
//	fprintf(fp, "label=Displacement\nnum_items=%d\n\n", No_points_for_colored_points);
//	for (p = 0; p < No_points_for_colored_points; p++)
//	{
//		fprintf(fp, "%d:%-.13le\n", p, pow((data_result_disp_x[p] * data_result_disp_x[p] + data_result_disp_y[p] * data_result_disp_y[p]), 0.5));
//	}
//	fclose(fp);


	/////////////////////////////////////////////////////////Gausspoint_coordinates.dat/////////////////////////////////////////////////////////
	/*fp = fopen("Gauss_stress/Gausspoint_coordinates.dat", "w");
	fprintf(fp, "%d\n", Total_Element_to_mesh[No_Input_File] * POW_Ng);
	Gausspoint_coordinate(Total_Element_to_mesh[No_Input_File], Total_Control_Point_to_mesh[No_Input_File]);
	for (I_F = 0; I_F < No_Input_File; I_F++)
	{
		for (re = real_Total_Element_to_mesh[I_F]; re < real_Total_Element_to_mesh[I_F + 1]; re++)
		{
			i = real_element[re];
			for (j = 0; j < POW_Ng; j++)
			{
				fprintf(fp, "%d\t%.16e %.16e %.16e\n", i, Gausspoint_coordinates[i][j][0], Gausspoint_coordinates[i][j][1], Gausspoint_coordinates[i][j][2]);
			}
		}
		fprintf(fp, "\n\n\n");
	}
	fclose(fp);
	/////////////////////////////////////////////////////////Gausspoint_coordinates.dat/////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////Gauss_stress_x.dat/////////////////////////////////////////////////////////
	fp = fopen("Gauss_stress/Gauss_stress_x.dat", "w");
	k = 0;
	fprintf(fp, "label=Stress x\nnum_items=%d\n\n", real_Total_Element_to_mesh[No_Input_File] * POW_Ng);
	for (I_F = 0; I_F < No_Input_File; I_F++)
	{
		for (re = real_Total_Element_to_mesh[I_F]; re < real_Total_Element_to_mesh[I_F + 1]; re++)
		{
			i = real_element[re];
			for (j = 0; j < POW_Ng; j++)
			{
				fprintf(fp, "%d:%.16e\n", k, Stress[i][j][0]);
				k++;
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n\n\n");
	}
	fclose(fp);
	/////////////////////////////////////////////////////////Gauss_stress_x.dat/////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////Gauss_stress_y.dat/////////////////////////////////////////////////////////
	fp = fopen("Gauss_stress/Gauss_stress_y.dat", "w");
	k = 0;
	fprintf(fp, "label=Stress y\nnum_items=%d\n\n", real_Total_Element_to_mesh[No_Input_File] * POW_Ng);
	for (I_F = 0; I_F < No_Input_File; I_F++)
	{
		for (re = real_Total_Element_to_mesh[I_F]; re < real_Total_Element_to_mesh[I_F + 1]; re++)
		{
			i = real_element[re];
			for (j = 0; j < POW_Ng; j++)
			{
				fprintf(fp, "%d:%.16e\n", k, Stress[i][j][1]);
				k++;
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n\n\n");
	}
	fclose(fp);
	/////////////////////////////////////////////////////////Gauss_stress_y.dat/////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////Gauss_stress_z.dat/////////////////////////////////////////////////////////
	fp = fopen("Gauss_stress/Gauss_stress_z.dat", "w");
	k = 0;
	fprintf(fp, "label=Stress z\nnum_items=%d\n\n", real_Total_Element_to_mesh[No_Input_File] * POW_Ng);
	for (I_F = 0; I_F < No_Input_File; I_F++)
	{
		for (re = real_Total_Element_to_mesh[I_F]; re < real_Total_Element_to_mesh[I_F + 1]; re++)
		{
			i = real_element[re];
			for (j = 0; j < POW_Ng; j++)
			{
				fprintf(fp, "%d:%.16e\n", k, Stress[i][j][2]);
				k++;
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n\n\n");
	}
	fclose(fp);
	/////////////////////////////////////////////////////////Gauss_stress_z.dat/////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////Gauss_stress_xy.dat/////////////////////////////////////////////////////////
	fp = fopen("Gauss_stress/Gauss_stress_xy.dat", "w");
	k = 0;
	fprintf(fp, "label=Stress xy\nnum_items=%d\n\n", real_Total_Element_to_mesh[No_Input_File] * POW_Ng);
	for (I_F = 0; I_F < No_Input_File; I_F++)
	{
		for (re = real_Total_Element_to_mesh[I_F]; re < real_Total_Element_to_mesh[I_F + 1]; re++)
		{
			i = real_element[re];
			for (j = 0; j < POW_Ng; j++)
			{
				fprintf(fp, "%d:%.16e\n", k, Stress[i][j][3]);
				k++;
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n\n\n");
	}
	fclose(fp);
	/////////////////////////////////////////////////////////Gauss_stress_xy.dat/////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////Gauss_stress_yz.dat/////////////////////////////////////////////////////////
	fp = fopen("Gauss_stress/Gauss_stress_yz.dat", "w");
	k = 0;
	fprintf(fp, "label=Stress yz\nnum_items=%d\n\n", real_Total_Element_to_mesh[No_Input_File] * POW_Ng);
	for (I_F = 0; I_F < No_Input_File; I_F++)
	{
		for (re = real_Total_Element_to_mesh[I_F]; re < real_Total_Element_to_mesh[I_F + 1]; re++)
		{
			i = real_element[re];
			for (j = 0; j < POW_Ng; j++)
			{
				fprintf(fp, "%d:%.16e\n", k, Stress[i][j][4]);
				k++;
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n\n\n");
	}
	fclose(fp);
	/////////////////////////////////////////////////////////Gauss_stress_yz.dat/////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////Gauss_stress_zx.dat/////////////////////////////////////////////////////////
	fp = fopen("Gauss_stress/Gauss_stress_zx.dat", "w");
	k = 0;
	fprintf(fp, "label=Stress zx\nnum_items=%d\n\n", real_Total_Element_to_mesh[No_Input_File] * POW_Ng);
	for (I_F = 0; I_F < No_Input_File; I_F++)
	{
		for (re = real_Total_Element_to_mesh[I_F]; re < real_Total_Element_to_mesh[I_F + 1]; re++)
		{
			i = real_element[re];
			for (j = 0; j < POW_Ng; j++)
			{
				fprintf(fp, "%d:%.16e\n", k, Stress[i][j][5]);
				k++;
			}
			fprintf(fp, "\n");
		}
		fprintf(fp, "\n\n\n");
	}
	fclose(fp);*/
	/////////////////////////////////////////////////////////Gauss_stress_zx.dat/////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////////OUTPUT FILE/////////////////////////////////////////////////////////////////////////////////////

	end3 = clock();
	printf("解析終了するのに%.2f秒かかりました\n",(double)(end3-start)/CLOCKS_PER_SEC);

	return 0;
}

///////////////////////////////////////////////////////
//////////////全体剛性マトリックスの制作///////////////
///////////////////////////////////////////////////////

//ファイルからデータをもらう
void Get_InputData(int I_F, double *E, double *nu, int *Total_Element, int *Total_Control_Point, int *Total_Load, int *No_Patch, int Load_Node_Dir[MAX_N_LOAD][2], double Value_of_Load[MAX_N_LOAD], int *Total_Constraint, int Constraint_Node_Dir[MAX_N_CONSTRAINT][2], double Value_of_Constraint[MAX_N_CONSTRAINT], int *Total_DistributeForce, int argc, char *argv[])
{
	int i, j, k, l;
	int l_all, i_all, ii_all;
	int p, q, t, x, y, z;
	char s[256];
	int ii, jj, kk, kkk;
	int e, b, B;
	int iiloc, jjloc, kkloc;
	int r = 0;
	/* for the distributed loads*/

	if ((fp = fopen(argv[I_F + 1], "r")) == NULL)
	{
		printf("file open error!!\n");
	}

	//材料定数
	fscanf(fp, "%le %le", &*E, &*nu);
	fgets(s, 256, fp);

	//if (I_F == 0)
	//{
	//	printf("E:%le nu:%le\n", *E, *nu);
	//}

	printf("----------------------------------------------------------------------------------------No_input_File : %d----------------------------------------------------------------------------------------\n",I_F);

	//パッチ数
	fscanf(fp, "%d", &*No_Patch);
	fgets(s, 256, fp);
	printf("No_Patch:%d\n", *No_Patch);

	Total_Patch_on_mesh[I_F] = *No_Patch;
	Total_Patch_to_mesh[I_F + 1] = Total_Patch_to_mesh[I_F] + *No_Patch;
	printf("Total_Patch_to_mesh[%d] = %d\n", I_F, Total_Patch_to_mesh[I_F]);


	//コントロールポイント数
	fscanf(fp, "%d", &*Total_Control_Point);
	fgets(s, 256, fp);

    Total_Control_Point_on_mesh[I_F] = *Total_Control_Point;
	Total_Control_Point_to_mesh[I_F + 1] = Total_Control_Point_to_mesh[I_F] + *Total_Control_Point;
	printf("Total_Control_Point = %d\tTotal_Control_Point_to_mesh[%d] = %d\n", *Total_Control_Point, I_F, Total_Control_Point_to_mesh[I_F]);



	//ξηζ方向の各次数
	for (l = 0; l < *No_Patch; l++)
	{
		l_all = l + Total_Patch_to_mesh[I_F];
		for (j = 0; j < DIMENSION; j++)
		{
			fscanf(fp, "%d", &Order[l_all][j]);
			//printf("Order[%d][%d] = %d\n",l_all, j, Order[l_all][j]);
		}
	}
	fgets(s, 256, fp);


	//ノット数
	for (l = 0; l < *No_Patch; l++)
	{
		l_all = l + Total_Patch_to_mesh[I_F];

		for (j = 0; j < DIMENSION; j++)
		{
			fscanf(fp, "%d", &No_knot[l_all][j]);
			printf("No_knot[%d][%d] = %d\n",l_all, j, No_knot[l_all][j]);
		}
	}
	fgets(s, 256, fp);


	//各パッチ各方向のコントロールポイント数
	for (l = 0; l < *No_Patch; l++)
	{
		l_all = l + Total_Patch_to_mesh[I_F];

		for (j = 0; j < DIMENSION; j++)
		{
			fscanf(fp, "%d", &No_Control_point[l_all][j]);
			printf("No_Control_point[%d][%d] = %d\n", l_all, j, No_Control_point[l_all][j]);
		}
	}
	fgets(s, 256, fp);

	for (l = 0; l < *No_Patch; l++)
	{
		l_all = l + Total_Patch_to_mesh[I_F];

		No_Controlpoint_in_patch[l_all] = 1.0;
	}

	for (l = 0; l < *No_Patch; l++)
	{
		l_all = l + Total_Patch_to_mesh[I_F];

		for (j = 0; j < DIMENSION; j++)
		{
			No_Controlpoint_in_patch[l_all] *= No_Control_point[l_all][j];
		}
	}

	for (l = 0; l < *No_Patch; l++)
	{
		l_all = l + Total_Patch_to_mesh[I_F];

		for (j = 0; j < DIMENSION; j++)
		{
			if (No_knot[l_all][j] != No_Control_point[l_all][j] + Order[l_all][j] + 1)
			{
				printf("wrong relationship between the number of knot vector and the number of control_point \n");
				printf("in patch_No.%d direction:%d\n", l, j);
			}
		}
	}

	for (l = 0; l < *No_Patch; l++)
	{
		l_all = l + Total_Patch_to_mesh[I_F];

		printf("No_Controlpoint_in_patch[%d]:%d\n", l_all, No_Controlpoint_in_patch[l_all]);
	}
	//printf("\n");

	for (l = 0; l < *No_Patch; l++)
	{
		l_all = l + Total_Patch_to_mesh[I_F];

		for (i = 0; i < No_Controlpoint_in_patch[l_all]; i++)
		{
			fscanf(fp, "%d", &Patch_controlpoint[l_all][i]);

			if (I_F > 0)
			{
				Patch_controlpoint[l_all][i] += Total_Control_Point_to_mesh[I_F];
			}
		}
	}

	/*for (l = 0; l < *No_Patch; l++)
	{
		l_all = l + Total_Patch_to_mesh[I_F];

		printf("%d\t:\t", l_all);

		for (i = 0; i < No_Controlpoint_in_patch[l_all]; i++)
		{
			printf("%d\t", Patch_controlpoint[l_all][i]);
		}
		printf("\n");
	}*/


	//変位，荷重，分布荷重の個数
	fscanf(fp, "%d %d %d", Total_Constraint, Total_Load, Total_DistributeForce);
	fgets(s, 256, fp);

	Total_Constraint_on_mesh[I_F] = *Total_Constraint;
	Total_Constraint_to_mesh[I_F + 1] = Total_Constraint_to_mesh[I_F] + *Total_Constraint;

	Total_Load_on_mesh[I_F] = *Total_Load;
	Total_Load_to_mesh[I_F + 1] = Total_Load_to_mesh[I_F] + *Total_Load;

	Total_DistributeForce_on_mesh[I_F] = *Total_DistributeForce;
	Total_DistributeForce_to_mesh[I_F + 1] = Total_DistributeForce_to_mesh[I_F]+ *Total_DistributeForce;

	printf("Total_Constraint = %d\n", *Total_Constraint);
    printf("Total_Constraint_on_mesh[%d] = %d\n",I_F, Total_Constraint_on_mesh[I_F]);
    printf("Total_Constraint_to_mesh[%d] = %d\n",I_F, Total_Constraint_to_mesh[I_F]);

	printf("Total_Load = %d\n", *Total_Load);
	printf("Total_Load_on_mesh[%d] = %d\n",I_F, Total_Load_on_mesh[I_F]);
    printf("Total_Load_to_mesh[%d] = %d\n",I_F, Total_Load_to_mesh[I_F]);

	printf("Total_DistributedForce = %d\n", *Total_DistributeForce);
	printf("Total_DistributeForce_on_mesh[%d] = %d\n",I_F, Total_DistributeForce_on_mesh[I_F]);
    printf("Total_DistributeForce_to_mesh[%d] = %d\n",I_F, Total_DistributeForce_to_mesh[I_F]);


	//ノットベクトルの読み込み
	for (l = 0; l < *No_Patch; l++)
	{
		l_all = l + Total_Patch_to_mesh[I_F];

		for (j = 0; j < DIMENSION; j++)
		{
			for (k = 0; k < No_knot[l_all][j]; k++)
			{
				fscanf(fp, "%le", &Position_Knots[l_all][j][k]);
				//printf("%le\t", Position_Knots[l_all][j][k]);
			}
			//printf("\n");
		}
	}

	for (l = 0; l < *No_Patch; l++)
	{
		l_all = l + Total_Patch_to_mesh[I_F];

		No_Control_point_ON_ELEMENT[l_all] = 1.0;
	}

	*Total_Element = 0.0;

	for(l = 0; l < *No_Patch; l++)
	{
		l_all = l + Total_Patch_to_mesh[I_F];

		*Total_Element += (No_Control_point[l_all][0] - Order[l_all][0]) * (No_Control_point[l_all][1] - Order[l_all][1]) * (No_Control_point[l_all][2] - Order[l_all][2]);
		No_Control_point_ON_ELEMENT[l_all] = (Order[l_all][0] + 1) * (Order[l_all][1] + 1) * (Order[l_all][2] + 1);
	}

    Total_Element_on_mesh[I_F] = *Total_Element;
	Total_Element_to_mesh[I_F + 1] = Total_Element_to_mesh[I_F] + *Total_Element;

	printf("Total_Element=%d\n", *Total_Element);
	printf("Total_Element_on_mesh[%d] = %d\n",I_F, Total_Element_on_mesh[I_F]);
    printf("Total_Element_to_mesh[%d] = %d\n",I_F, Total_Element_to_mesh[I_F]);


	//節点座標
	for (i = 0; i < *Total_Control_Point; i++)
	{
		fscanf(fp, "%d", &ii);

		ii_all = ii + Total_Control_Point_to_mesh[I_F];

		for (j = 0; j < DIMENSION + 1; j++)
		{
			fscanf(fp, "%le", &Node_Coordinate[ii_all][j]); //Node_Coordinate[ii][3]:重み
		}
	}
	fgets(s, 256, fp);

	for (i = 0; i < *Total_Control_Point; i++)
	{
		i_all = i + Total_Control_Point_to_mesh[I_F];

		for (j = 0; j < DIMENSION + 1; j++)
		{
			//コントロールポイント座標・重みの新たな配列（for s-IGA/NewtonLaphson）
			if (j < DIMENSION)
			{
				Control_Coord[j][i_all] = Node_Coordinate[i_all][j];
			}
			else if (j == DIMENSION)
			{
				Control_Weight[i_all] = Node_Coordinate[i_all][DIMENSION];
			}
			//printf("Control_Coord[%d][%d]=%e\n", j, i_all, Control_Coord[j][i_all]);
		}
	}

	/*for (i = 0; i < *Total_Control_Point; i++)
	{
		i_all = i + Total_Control_Point_to_mesh[I_F];

		printf("Node_Coord/Weight\t:\t");

		for (j = 0; j < DIMENSION; j++)
		{
			printf("[%d][%d] = %e\t", i_all, j, Node_Coordinate[i_all][j]);
		}
		printf("[%d][3] = %e\n", i_all, Node_Coordinate[i_all][3]);
	}*/


	//拘束
	for (i = 0; i < *Total_Constraint; i++)
	{
		i_all = i + Total_Constraint_to_mesh[I_F];

		fscanf(fp, "%d %d %le", &Constraint_Node_Dir[i_all][0], &Constraint_Node_Dir[i_all][1], &Value_of_Constraint[i_all]);
	}
	fgets(s, 256, fp);

	for (i = 0; i < *Total_Constraint; i++)
	{
		i_all = i + Total_Constraint_to_mesh[I_F];

		Constraint_Node_Dir[i_all][0] = Constraint_Node_Dir[i_all][0] + Total_Control_Point_to_mesh[I_F];
		Constraint_Node_Dir_on_mesh[I_F][i][0] = Constraint_Node_Dir[i][0];
		Constraint_Node_Dir_on_mesh[I_F][i][1] = Constraint_Node_Dir[i][1];
        Value_of_Constraint_on_mesh[I_F][i] = Value_of_Constraint[i];

		//printf("Constraint_Node_Dir[%d][0] = %d\tConstraint_Node_Dir[%d][1] = %d\tValue_of_Constraint[%d] = %e \n", i_all, Constraint_Node_Dir[i_all][0], i_all, Constraint_Node_Dir[i_all][1], i_all, Value_of_Constraint[i_all]);
	}


	//荷重
	for (i = 0; i < *Total_Load; i++)
	{
		i_all = i + Total_Load_to_mesh[I_F];

		fscanf(fp, "%d %d %le", &Load_Node_Dir[i_all][0], &Load_Node_Dir[i_all][1], &Value_of_Load[i_all]);
	}
	fgets(s, 256, fp);

	for (i = 0; i < *Total_Load; i++)
	{
		i_all = i + Total_Load_to_mesh[I_F];

		Load_Node_Dir[i_all][0] = Load_Node_Dir[i_all][0] + Total_Control_Point_to_mesh[I_F];

		//printf("Load_Node_Dir[%d][0] = %d\tLoad_Node_Dir[%d][1] = %d\tValue_of_Load[%d] = %e\n", i_all, Load_Node_Dir[i_all][0], i_all, Load_Node_Dir[i_all][1], i_all, Value_of_Load[i_all]);
	}


	int iPatch, iCoord, jCoord, type_load;
	double iRange_Coord[2], jRange_Coord[2], val_Coord, Coeff_Dist_Load_i[3], Coeff_Dist_Load_j[3];
	int iPatch_array[MAX_N_DISTRIBUTE_FORCE], iCoord_array[MAX_N_DISTRIBUTE_FORCE], jCoord_array[MAX_N_DISTRIBUTE_FORCE], type_load_array[MAX_N_DISTRIBUTE_FORCE];
	double val_Coord_array[MAX_N_DISTRIBUTE_FORCE], iRange_Coord_array[MAX_N_DISTRIBUTE_FORCE][2], jRange_Coord_array[MAX_N_DISTRIBUTE_FORCE][2], Coeff_Dist_Load_i_array[MAX_N_DISTRIBUTE_FORCE][3], Coeff_Dist_Load_j_array[MAX_N_DISTRIBUTE_FORCE][3];


	//分布荷重
	for (i = 0; i < *Total_DistributeForce; i++)
	{
		fscanf(fp, "%d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &type_load, &iPatch, &iCoord,  &jCoord, &val_Coord, &iRange_Coord[0], &iRange_Coord[1], &jRange_Coord[0], &jRange_Coord[1], &Coeff_Dist_Load_i[0], &Coeff_Dist_Load_i[1], &Coeff_Dist_Load_i[2], &Coeff_Dist_Load_j[0], &Coeff_Dist_Load_j[1], &Coeff_Dist_Load_j[2]);
		//printf(" Total_DistributeForce_to_mesh[%d] = %d\n", I_F,  Total_DistributeForce_to_mesh[I_F]);
		//printf("Distibuted load number=%d :", i);
		//printf("type_load: %d  iPatch: %d iCoord: %d  jCoord: %d val_Coord: %15.7e \n iRange_Coord: %15.7e  %15.7e jRange_Coord: %15.7e  %15.7e\n Coef_Dist_Load_i: %15.7e %15.7e %15.7e\n Coef_Dist_Load_j: %15.7e %15.7e %15.7e\n", type_load, iPatch, iCoord, jCoord, val_Coord, iRange_Coord[0], iRange_Coord[1], jRange_Coord[0], jRange_Coord[1], Coeff_Dist_Load_i[0], Coeff_Dist_Load_i[1], Coeff_Dist_Load_i[2], Coeff_Dist_Load_j[0], Coeff_Dist_Load_j[1], Coeff_Dist_Load_j[2]);
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
		i_all = i + Total_DistributeForce_to_mesh[I_F];

		type_load_array[i_all] = type_load;
		iPatch_array[i_all] = iPatch + Total_Patch_to_mesh[I_F];
		iCoord_array[i_all] = iCoord;
		jCoord_array[i_all] = jCoord;
		val_Coord_array[i_all] = val_Coord;
		iRange_Coord_array[i_all][0] = iRange_Coord[0];
		iRange_Coord_array[i_all][1] = iRange_Coord[1];
		jRange_Coord_array[i_all][0] = jRange_Coord[0];
		jRange_Coord_array[i_all][1] = jRange_Coord[1];
		Coeff_Dist_Load_i_array[i_all][0] = Coeff_Dist_Load_i[0];
		Coeff_Dist_Load_i_array[i_all][1] = Coeff_Dist_Load_i[1];
		Coeff_Dist_Load_i_array[i_all][2] = Coeff_Dist_Load_i[2];
		Coeff_Dist_Load_j_array[i_all][0] = Coeff_Dist_Load_j[0];
		Coeff_Dist_Load_j_array[i_all][1] = Coeff_Dist_Load_j[1];
		Coeff_Dist_Load_j_array[i_all][2] = Coeff_Dist_Load_j[2];

		//printf("\t%d\t%d\t%d\t%d\t%lf\t\t%lf\t%lf\t\t%lf\t%lf\t\t%lf\t%lf\t%lf\t\t%lf\t%lf\t%lf\n",type_load_array[i_all], iPatch_array[i_all], iCoord_array[i_all], jCoord_array[i_all], val_Coord_array[i_all], iRange_Coord_array[i_all][0], iRange_Coord_array[i_all][1], jRange_Coord_array[i_all][0], jRange_Coord_array[i_all][1], Coeff_Dist_Load_i_array[i_all][0], Coeff_Dist_Load_i_array[i_all][1], Coeff_Dist_Load_i_array[i_all][2], Coeff_Dist_Load_j_array[i_all][0], Coeff_Dist_Load_j_array[i_all][1], Coeff_Dist_Load_j_array[i_all][2] );
	}

	fscanf(fp, "%d", &element_ndiv_loc_glo[I_F]);
	printf("element_ndiv_loc_glo[%d] = %d\n", I_F, element_ndiv_loc_glo[I_F]);

	fscanf(fp, "%d", &No_GP);
	printf("No_GP = %d\n",No_GP);

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
	//INC\の計算（節点番号をξ，η，ζの番号で表す為の配列）

	e = 0;
	for (l = 0; l < *No_Patch; l++)
	{
		//printf("No_patch(INC)=%d\n", l);
		i = 0;

		l_all = l + Total_Patch_to_mesh[I_F];//確認違うかも
		//e_all = e + Total_Element_to_mesh[I_F];

		for (kk = 0; kk < No_Control_point[l_all][2]; kk++)
		{
			for (jj = 0; jj < No_Control_point[l_all][1]; jj++)
			{
				for (ii = 0; ii < No_Control_point[l_all][0]; ii++)
				{
					//printf("kk=%d\n",kk );
					INC[l_all][Patch_controlpoint[l_all][i]][0] = ii;
					INC[l_all][Patch_controlpoint[l_all][i]][1] = jj;
					INC[l_all][Patch_controlpoint[l_all][i]][2] = kk;
					//printf("Patch_controlpoint[%d][%d]=%d\n", l, i, Patch_controlpoint[l][i]);
					//printf("INC[%d][%d][0]=%d INC[%d][%d][1]=%d INC[%d][%d][2]=%d\n", l_all, Patch_controlpoint[l_all][i], INC[l_all][Patch_controlpoint[l_all][i]][0], l_all, Patch_controlpoint[l_all][i], INC[l_all][Patch_controlpoint[l_all][i]][1], l_all, Patch_controlpoint[l_all][i], INC[l_all][Patch_controlpoint[l_all][i]][2]);

					//Adress_Controlpoint[l_all][ii][jj][kk] = Patch_controlpoint[l_all][i];

					if (ii >= Order[l_all][0] && jj >= Order[l_all][1] && kk >= Order[l_all][2])
					{
						for (kkloc = 0; kkloc <= Order[l_all][2]; kkloc++)
						{
							for (jjloc = 0; jjloc <= Order[l_all][1]; jjloc++)
							{
								for (iiloc = 0; iiloc <= Order[l_all][0]; iiloc++)
								{
									//printf("kkloc:%d jjloc:%d iiloc:%d\n",kkloc,jjloc,iiloc);
									B = Patch_controlpoint[l_all][i - kkloc * No_Control_point[l_all][0] * No_Control_point[l_all][1] - jjloc * No_Control_point[l_all][0] - iiloc];
									b = kkloc * (Order[l_all][1] + 1) * (Order[l_all][1] + 1) + jjloc * (Order[l_all][0] + 1) + iiloc;;

                                    //Controlpoint_of_Element_nomarge[e + Total_Element_to_mesh[I_F]][b] = B;
									Controlpoint_of_Element[e + Total_Element_to_mesh[I_F]][b] = B;
								}
							}
						}
						Element_patch[e + Total_Element_to_mesh[I_F]] = l_all;
						Element_mesh[e + Total_Element_to_mesh[I_F]] = I_F;
						El_No_on_mesh[I_F][e] = e + Total_Element_to_mesh[I_F]; //確認確認確認

						//printf("Element_patch[%d] = %d\tElement_mesh[%d] = %d\tEl_No_on_mesh[%d][%d] = %d\n", e + Total_Element_to_mesh[I_F], Element_patch[e + Total_Element_to_mesh[I_F]], e + Total_Element_to_mesh[I_F], Element_mesh[e + Total_Element_to_mesh[I_F]], I_F, e, El_No_on_mesh[I_F][e]);

						e++;
					}
					i++;
				}
			}
		}
		Patch_mesh[l_all] = I_F;
		printf("Patch_mesh[%d] = %d\n", l_all, I_F);
	}

	/*----------------------------------------------------------------------------------------------*/

    for (l = 0; l < *No_Patch; l++)
	{
		l_all = l + Total_Patch_to_mesh[I_F];//確認違うかも

		for (j = 0; j < DIMENSION; j++)
		{
			line_No_real_element[l_all][j] = 0;
        }
    }

	for (l = 0; l < *No_Patch; l++)
	{
		l_all = l + Total_Patch_to_mesh[I_F];//確認違うかも

		for (j = 0; j < DIMENSION; j++)
		{
			line_No_Total_element[l_all][j] = No_knot[l_all][j] - 2 * Order[l_all][j] - 1;

			for (kkk = Order[l_all][j]; kkk < No_knot[l_all][j] - Order[l_all][j] - 1; kkk++)
			{
				difference[l_all][kkk - Order[l_all][j]][j] = Position_Knots[l_all][j][kkk + 1] - Position_Knots[l_all][j][kkk];

				if (difference[l_all][kkk - Order[l_all][j]][j] != 0.0)
				{
					line_No_real_element[l_all][j]++;
				}
			}
			printf("line_No_real_element[%d][%d] = %d\n", l_all, j, line_No_real_element[l_all][j]);
		}
	}

	/*要素に行番号、列番号をつける*/
	for (i = 0; i < *Total_Element; i++)
	{
		Total_element_all_ID[i] = 0;
	}

	i = 0;
	for (l = 0; l < *No_Patch; l++)
	{
		l_all = l + Total_Patch_to_mesh[I_F];//確認違うかも

		for (z = 0; z < line_No_Total_element[l_all][2]; z++)
		{
			for (y = 0; y < line_No_Total_element[l_all][1]; y++)
			{
				for (x = 0; x < line_No_Total_element[l_all][0]; x++)
				{
					//printf("line_No_Total_element [%d][2] = %d\t[%d][1] = %d\t[%d][0] = %d\n", l_all, line_No_Total_element[l_all][2], l_all, line_No_Total_element[l_all][1], l_all, line_No_Total_element[l_all][0]);

					i_all = i + Total_Element_to_mesh[I_F];

					ENC[l_all][i_all][0] = x;
					ENC[l_all][i_all][1] = y;
					ENC[l_all][i_all][2] = z;

					//printf("ENC[%d][%d][0] = %d\t[%d][%d][1] = %d\t[%d][%d][2] = %d\n", l_all, i_all, ENC[l_all][i_all][0], l_all, i_all, ENC[l_all][i_all][1], l_all, i_all, ENC[l_all][i_all][2]);
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
            e = 0;

			l_all = l + Total_Patch_to_mesh[I_F];

			for (k = 0; k < line_No_Total_element[l_all][j]; k++)
			{
				if (difference[l_all][k][j] != 0.0)
				{
					real_element_line[l_all][e][j] = k;

					//printf("real_element_line[%d][%d][%d] = %d\n", l_all, e, j, real_element_line[l_all][e][j]);

					e++;
				}
			}
		}
	}


	/*必要な要素列上の要素のIDを1にする*/
	for (i = 0; i < *Total_Element; i++)
	{
		i_all = i + Total_Element_to_mesh[I_F]; //確認確認確認

		for (p = 0; p < line_No_real_element[Element_patch[i_all]][0]; p++)
		{
			if (ENC[Element_patch[i_all]][i_all][0] == real_element_line[Element_patch[i_all]][p][0])
			{
				for (q = 0; q < line_No_real_element[Element_patch[i_all]][1]; q++)
				{
					if (ENC[Element_patch[i_all]][i_all][1] == real_element_line[Element_patch[i_all]][q][1])
					{
						for (t = 0; t < line_No_real_element[Element_patch[i_all]][2]; t++)
						{
							if (ENC[Element_patch[i_all]][i_all][2] == real_element_line[Element_patch[i_all]][t][2])
							{
								Total_element_all_ID[i]++; //ほんとにi?i_allじゃない, 確認確認確認?
							}
						}
					}
				}
			}
		}


		/*IDが1の要素に番号を振る*/
		if (Total_element_all_ID[i] == 1)
		{
			real_element[r + real_Total_Element_to_mesh[I_F]] = i_all;
			real_El_No_on_mesh[I_F][r] = i_all;

			//printf("real_element[%d]=%d\n", r + real_Total_Element_to_mesh[I_F], real_element[r + real_Total_Element_to_mesh[I_F]]);
			r++;
		}
	}

	real_Total_Element = 0;

	for (l = 0; l < *No_Patch; l++)
	{
		l_all = l + Total_Patch_to_mesh[I_F];//確認違うかも

		real_Total_Element += line_No_real_element[l_all][0] * line_No_real_element[l_all][1] * line_No_real_element[l_all][2];
	}

	real_Total_Element_on_mesh[I_F] = real_Total_Element;
    real_Total_Element_to_mesh[I_F + 1] = real_Total_Element_to_mesh[I_F] + real_Total_Element;
    printf("%d:real_Total_Element:%d\n", I_F, real_Total_Element);
    /*int rr;
	for(rr=0;rr<real_Total_Element;rr++)
    {
        printf("real_element[%d]=%d\n",rr,real_element[rr]);
    }*/
	//}
	///

	for (i = 0; i < *Total_Control_Point; i++)
	{
		i_all = i + Total_Control_Point_to_mesh[I_F];

		Equivalent_Nodal_Force[i_all][0] = 0.0;
		Equivalent_Nodal_Force[i_all][1] = 0.0;
		Equivalent_Nodal_Force[i_all][2] = 0.0;
	}

	for (i = 0; i < *Total_DistributeForce; i++)
	{
		type_load = type_load_array[i + Total_DistributeForce_to_mesh[I_F]];
		iPatch = iPatch_array[i + Total_DistributeForce_to_mesh[I_F]];
		iCoord = iCoord_array[i + Total_DistributeForce_to_mesh[I_F]];
		jCoord = jCoord_array[i + Total_DistributeForce_to_mesh[I_F]];
		val_Coord = val_Coord_array[i + Total_DistributeForce_to_mesh[I_F]];
		iRange_Coord[0] = iRange_Coord_array[i + Total_DistributeForce_to_mesh[I_F]][0];
		iRange_Coord[1] = iRange_Coord_array[i + Total_DistributeForce_to_mesh[I_F]][1];
		jRange_Coord[0] = jRange_Coord_array[i + Total_DistributeForce_to_mesh[I_F]][0];
		jRange_Coord[1] = jRange_Coord_array[i + Total_DistributeForce_to_mesh[I_F]][1];
		Coeff_Dist_Load_i[0] = Coeff_Dist_Load_i_array[i + Total_DistributeForce_to_mesh[I_F]][0];
		Coeff_Dist_Load_i[1] = Coeff_Dist_Load_i_array[i + Total_DistributeForce_to_mesh[I_F]][1];
		Coeff_Dist_Load_i[2] = Coeff_Dist_Load_i_array[i + Total_DistributeForce_to_mesh[I_F]][2];
		Coeff_Dist_Load_j[0] = Coeff_Dist_Load_j_array[i + Total_DistributeForce_to_mesh[I_F]][0];
		Coeff_Dist_Load_j[1] = Coeff_Dist_Load_j_array[i + Total_DistributeForce_to_mesh[I_F]][1];
		Coeff_Dist_Load_j[2] = Coeff_Dist_Load_j_array[i + Total_DistributeForce_to_mesh[I_F]][2];

		//printf("iPatch = %d\n",iPatch);


		//Setting_Dist_Load(I_F, Total_Control_Point_to_mesh[I_F + 1], iPatch, Total_Element_to_mesh[I_F + 1], iCoord, jCoord, val_Coord, iRange_Coord,  jRange_Coord,type_load, Coeff_Dist_Load_i, Coeff_Dist_Load_j);
	}
	/*-------------------------------------------------------------------------------------*/
	printf("-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
}

void Get_displacment(int Total_Control_Point, const char *fileName)
{
    FILE *fp;

    int i, iii, disp_no;

    fp = fopen(fileName, "r");

    fscanf(fp, "%d", &disp_no);

    for (i = 0; i < disp_no; i++)
    {
        fscanf(fp, "%d %le %le %le", &iii, &Displacement[i * DIMENSION + 0], &Displacement[i * DIMENSION + 1], &Displacement[i * DIMENSION + 2]);
		//printf("%le", Displacement[i * DIMENSION + 2]);
    }

	fclose(fp);
}


void Make_Gauss_points(int Select_GP)
{
    int i, j, k, ii;

	if(Select_GP == 0)
	{
		GaussPt_1dir = Ng;
	}
	if(Select_GP == 1)
	{
		GaussPt_1dir = No_GP;
	}

	GaussPt_3D = GaussPt_1dir * GaussPt_1dir * GaussPt_1dir;

	if (GaussPt_1dir == 3)
	{
		G_1d[0]  = -0.7745966692414834;
		G_1d[1]  =  0.0000000000000000;
		G_1d[2]  =  0.7745966692414834;
		G_1d[3]  =  0.0;
		G_1d[4]  =  0.0;
		G_1d[5]  =  0.0;
		G_1d[6]  =  0.0;
		G_1d[7]  =  0.0;
		G_1d[8]  =  0.0;
		G_1d[9]  =  0.0;
		G_1d[10] =  0.0;
		G_1d[11] =  0.0;
		G_1d[12] =  0.0;
		G_1d[13] =  0.0;
		G_1d[14] =  0.0;
		G_1d[15] =  0.0;
		G_1d[16] =  0.0;
		G_1d[17] =  0.0;
		G_1d[18] =  0.0;
		G_1d[19] =  0.0;

		w_1d[0]  = 0.5555555555555556;
		w_1d[1]  = 0.8888888888888889;
		w_1d[2]  = 0.5555555555555556;
		w_1d[3]  = 0.0;
		w_1d[4]  = 0.0;
		w_1d[5]  = 0.0;
		w_1d[6]  = 0.0;
		w_1d[7]  = 0.0;
		w_1d[8]  = 0.0;
		w_1d[9]  = 0.0;
		w_1d[10] = 0.0;
		w_1d[11] = 0.0;
		w_1d[12] = 0.0;
		w_1d[13] = 0.0;
		w_1d[14] = 0.0;
		w_1d[15] = 0.0;
		w_1d[16] = 0.0;
		w_1d[17] = 0.0;
		w_1d[18] = 0.0;
		w_1d[19] = 0.0;
	}

	if (GaussPt_1dir == 4)
	{
		G_1d[0]  = -0.8611363115940526;
		G_1d[1]  = -0.3399810435848563;
		G_1d[2]  =  0.3399810435848563;
		G_1d[3]  =  0.8611363115940526;
		G_1d[4]  = 0.0;
		G_1d[5]  = 0.0;
		G_1d[6]  = 0.0;
		G_1d[7]  = 0.0;
		G_1d[8]  = 0.0;
		G_1d[9]  = 0.0;
		G_1d[10] = 0.0;
		G_1d[11] = 0.0;
		G_1d[12] = 0.0;
		G_1d[13] = 0.0;
		G_1d[14] = 0.0;
		G_1d[15] = 0.0;
		G_1d[16] = 0.0;
		G_1d[17] = 0.0;
		G_1d[18] = 0.0;
		G_1d[19] = 0.0;

		w_1d[0]  = 0.3478548451374537;
		w_1d[1]  = 0.6521451548625462;
		w_1d[2]  = 0.6521451548625462;
		w_1d[3]  = 0.3478548451374537;
		w_1d[4]  = 0.0;
		w_1d[5]  = 0.0;
		w_1d[6]  = 0.0;
		w_1d[7]  = 0.0;
		w_1d[8]  = 0.0;
		w_1d[9]  = 0.0;
		w_1d[10] = 0.0;
		w_1d[11] = 0.0;
		w_1d[12] = 0.0;
		w_1d[13] = 0.0;
		w_1d[14] = 0.0;
		w_1d[15] = 0.0;
		w_1d[16] = 0.0;
		w_1d[17] = 0.0;
		w_1d[18] = 0.0;
		w_1d[19] = 0.0;
	}

	if (GaussPt_1dir == 5)
	{
		G_1d[0]  = -0.9061798459386640;
		G_1d[1]  = -0.5384693101056831;
		G_1d[2]  =  0.0000000000000000;
		G_1d[3]  =  0.5384693101056831;
		G_1d[4]  =  0.9061798459386640;
		G_1d[5]  = 0.0;
		G_1d[6]  = 0.0;
		G_1d[7]  = 0.0;
		G_1d[8]  = 0.0;
		G_1d[9]  = 0.0;
		G_1d[10] = 0.0;
		G_1d[11] = 0.0;
		G_1d[12] = 0.0;
		G_1d[13] = 0.0;
		G_1d[14] = 0.0;
		G_1d[15] = 0.0;
		G_1d[16] = 0.0;
		G_1d[17] = 0.0;
		G_1d[18] = 0.0;
		G_1d[19] = 0.0;

		w_1d[0]  = 0.2369268850561894;
		w_1d[1]  = 0.4786286704993662;
		w_1d[2]  = 0.5688888888888890;
		w_1d[3]  = 0.4786286704993662;
		w_1d[4]  = 0.2369268850561894;
		w_1d[5]  = 0.0;
		w_1d[6]  = 0.0;
		w_1d[7]  = 0.0;
		w_1d[8]  = 0.0;
		w_1d[9]  = 0.0;
		w_1d[10] = 0.0;
		w_1d[11] = 0.0;
		w_1d[12] = 0.0;
		w_1d[13] = 0.0;
		w_1d[14] = 0.0;
		w_1d[15] = 0.0;
		w_1d[16] = 0.0;
		w_1d[17] = 0.0;
		w_1d[18] = 0.0;
		w_1d[19] = 0.0;
	}

	if (GaussPt_1dir == 6)
	{
		G_1d[0]  = -0.9324695142031521;
		G_1d[1]  = -0.6612093864662645;
		G_1d[2]  = -0.2386191860831969;
		G_1d[3]  =  0.2386191860831969;
		G_1d[4]  =  0.6612093864662645;
		G_1d[5]  =  0.9324695142031521;
		G_1d[6]  = 0.0;
		G_1d[7]  = 0.0;
		G_1d[8]  = 0.0;
		G_1d[9]  = 0.0;
		G_1d[10] = 0.0;
		G_1d[11] = 0.0;
		G_1d[12] = 0.0;
		G_1d[13] = 0.0;
		G_1d[14] = 0.0;
		G_1d[15] = 0.0;
		G_1d[16] = 0.0;
		G_1d[17] = 0.0;
		G_1d[18] = 0.0;
		G_1d[19] = 0.0;

		w_1d[0]  = 0.1713244923791697;
		w_1d[1]  = 0.3607615730481389;
		w_1d[2]  = 0.4679139345726914;
		w_1d[3]  = 0.4679139345726914;
		w_1d[4]  = 0.3607615730481389;
		w_1d[5]  = 0.1713244923791697;
		w_1d[6]  = 0.0;
		w_1d[7]  = 0.0;
		w_1d[8]  = 0.0;
		w_1d[9]  = 0.0;
		w_1d[10] = 0.0;
		w_1d[11] = 0.0;
		w_1d[12] = 0.0;
		w_1d[13] = 0.0;
		w_1d[14] = 0.0;
		w_1d[15] = 0.0;
		w_1d[16] = 0.0;
		w_1d[17] = 0.0;
		w_1d[18] = 0.0;
		w_1d[19] = 0.0;
	}

	if (GaussPt_1dir == 7)
	{
		G_1d[0]  = -0.9491079123427585;
		G_1d[1]  = -0.7415311855993945;
		G_1d[2]  = -0.4058451513773972;
		G_1d[3]  =  0.0000000000000000;
		G_1d[4]  =  0.4058451513773972;
		G_1d[5]  =  0.7415311855993945;
		G_1d[6]  =  0.9491079123427585;
		G_1d[7]  = 0.0;
		G_1d[8]  = 0.0;
		G_1d[9]  = 0.0;
		G_1d[10] = 0.0;
		G_1d[11] = 0.0;
		G_1d[12] = 0.0;
		G_1d[13] = 0.0;
		G_1d[14] = 0.0;
		G_1d[15] = 0.0;
		G_1d[16] = 0.0;
		G_1d[17] = 0.0;
		G_1d[18] = 0.0;
		G_1d[19] = 0.0;

		w_1d[0]  = 0.1294849661688706;
		w_1d[1]  = 0.2797053914892766;
		w_1d[2]  = 0.3818300505051183;
		w_1d[3]  = 0.4179591836734690;
		w_1d[4]  = 0.3818300505051183;
		w_1d[5]  = 0.2797053914892766;
		w_1d[6]  = 0.1294849661688706;
		w_1d[7]  = 0.0;
		w_1d[8]  = 0.0;
		w_1d[9]  = 0.0;
		w_1d[10] = 0.0;
		w_1d[11] = 0.0;
		w_1d[12] = 0.0;
		w_1d[13] = 0.0;
		w_1d[14] = 0.0;
		w_1d[15] = 0.0;
		w_1d[16] = 0.0;
		w_1d[17] = 0.0;
		w_1d[18] = 0.0;
		w_1d[19] = 0.0;
	}

	if (GaussPt_1dir == 8)
	{
		G_1d[0]  = -0.9602898564975362;
		G_1d[1]  = -0.7966664774136267;
		G_1d[2]  = -0.5255324099163290;
		G_1d[3]  = -0.1834346424956498;
		G_1d[4]  =  0.1834346424956498;
		G_1d[5]  =  0.5255324099163290;
		G_1d[6]  =  0.7966664774136267;
		G_1d[7]  =  0.9602898564975362;
		G_1d[8]  = 0.0;
		G_1d[9]  = 0.0;
		G_1d[10] = 0.0;
		G_1d[11] = 0.0;
		G_1d[12] = 0.0;
		G_1d[13] = 0.0;
		G_1d[14] = 0.0;
		G_1d[15] = 0.0;
		G_1d[16] = 0.0;
		G_1d[17] = 0.0;
		G_1d[18] = 0.0;
		G_1d[19] = 0.0;

		w_1d[0]  = 0.1012285362903767;
		w_1d[1]  = 0.2223810344533743;
		w_1d[2]  = 0.3137066458778870;
		w_1d[3]  = 0.3626837833783618;
		w_1d[4]  = 0.3626837833783618;
		w_1d[5]  = 0.3137066458778870;
		w_1d[6]  = 0.2223810344533743;
		w_1d[7]  = 0.1012285362903767;
		w_1d[8]  = 0.0;
		w_1d[9]  = 0.0;
		w_1d[10] = 0.0;
		w_1d[11] = 0.0;
		w_1d[12] = 0.0;
		w_1d[13] = 0.0;
		w_1d[14] = 0.0;
		w_1d[15] = 0.0;
		w_1d[16] = 0.0;
		w_1d[17] = 0.0;
		w_1d[18] = 0.0;
		w_1d[19] = 0.0;
	}


	if (GaussPt_1dir == 9)
	{
		G_1d[0]  = -0.9681602395076261;
		G_1d[1]  = -0.8360311073266358;
		G_1d[2]  = -0.6133714327005904;
		G_1d[3]  = -0.3242534234038089;
		G_1d[4]  = +0.0000000000000000;
		G_1d[5]  = +0.3242534234038089;
		G_1d[6]  = +0.6133714327005904;
		G_1d[7]  = +0.8360311073266358;
		G_1d[8]  = +0.9681602395076261;
		G_1d[9]  = 0.0;
		G_1d[10] = 0.0;
		G_1d[11] = 0.0;
		G_1d[12] = 0.0;
		G_1d[13] = 0.0;
		G_1d[14] = 0.0;
		G_1d[15] = 0.0;
		G_1d[16] = 0.0;
		G_1d[17] = 0.0;
		G_1d[18] = 0.0;
		G_1d[19] = 0.0;

		w_1d[0]  = 0.0812743883615747;
		w_1d[1]  = 0.1806481606948571;
		w_1d[2]  = 0.2606106964029357;
		w_1d[3]  = 0.3123470770400028;
		w_1d[4]  = 0.3302393550012597;
		w_1d[5]  = 0.3123470770400028;
		w_1d[6]  = 0.2606106964029357;
		w_1d[7]  = 0.1806481606948571;
		w_1d[8]  = 0.0812743883615747;
		w_1d[9]  = 0.0;
		w_1d[10] = 0.0;
		w_1d[11] = 0.0;
		w_1d[12] = 0.0;
		w_1d[13] = 0.0;
		w_1d[14] = 0.0;
		w_1d[15] = 0.0;
		w_1d[16] = 0.0;
		w_1d[17] = 0.0;
		w_1d[18] = 0.0;
		w_1d[19] = 0.0;
	}

	if (GaussPt_1dir == 10)
	{
		G_1d[0]  = -0.9739065285171717;
		G_1d[1]  = -0.8650633666889845;
		G_1d[2]  = -0.6794095682990244;
		G_1d[3]  = -0.4333953941292472;
		G_1d[4]  = -0.1488743389816312;
		G_1d[5]  =  0.1488743389816312;
		G_1d[6]  =  0.4333953941292472;
		G_1d[7]  =  0.6794095682990244;
		G_1d[8]  =  0.8650633666889845;
		G_1d[9]  =  0.9739065285171717;
		G_1d[10] = 0.0;
		G_1d[11] = 0.0;
		G_1d[12] = 0.0;
		G_1d[13] = 0.0;
		G_1d[14] = 0.0;
		G_1d[15] = 0.0;
		G_1d[16] = 0.0;
		G_1d[17] = 0.0;
		G_1d[18] = 0.0;
		G_1d[19] = 0.0;

		w_1d[0]  = 0.0666713443086881;
		w_1d[1]  = 0.1494513491505804;
		w_1d[2]  = 0.2190863625159820;
		w_1d[3]  = 0.2692667193099965;
		w_1d[4]  = 0.2955242247147530;
		w_1d[5]  = 0.2955242247147530;
		w_1d[6]  = 0.2692667193099965;
		w_1d[7]  = 0.2190863625159820;
		w_1d[8]  = 0.1494513491505804;
		w_1d[9]  = 0.0666713443086881;
		w_1d[10] = 0.0;
		w_1d[11] = 0.0;
		w_1d[12] = 0.0;
		w_1d[13] = 0.0;
		w_1d[14] = 0.0;
		w_1d[15] = 0.0;
		w_1d[16] = 0.0;
		w_1d[17] = 0.0;
		w_1d[18] = 0.0;
		w_1d[19] = 0.0;
	}

	if (GaussPt_1dir == 11)
	{
		G_1d[0]  = -0.9782286581460570;
		G_1d[1]  = -0.8870625997680953;
		G_1d[2]  = -0.7301520055740494;
		G_1d[3]  = -0.5190961292068118;
		G_1d[4]  = -0.2695431559523450;
		G_1d[5]  =  0.0000000000000000;
		G_1d[6]  =  0.2695431559523450;
		G_1d[7]  =  0.5190961292068118;
		G_1d[8]  =  0.7301520055740494;
		G_1d[9]  =  0.8870625997680953;
		G_1d[10] =  0.9782286581460570;
		G_1d[11] = 0.0;
		G_1d[12] = 0.0;
		G_1d[13] = 0.0;
		G_1d[14] = 0.0;
		G_1d[15] = 0.0;
		G_1d[16] = 0.0;
		G_1d[17] = 0.0;
		G_1d[18] = 0.0;
		G_1d[19] = 0.0;

		w_1d[0]  = 0.0556685671161732;
		w_1d[1]  = 0.1255803694649047;
		w_1d[2]  = 0.1862902109277344;
		w_1d[3]  = 0.2331937645919907;
		w_1d[4]  = 0.2628045445102468;
		w_1d[5]  = 0.2729250867779009;
		w_1d[6]  = 0.2628045445102468;
		w_1d[7]  = 0.2331937645919907;
		w_1d[8]  = 0.1862902109277344;
		w_1d[9]  = 0.1255803694649047;
		w_1d[10] = 0.0556685671161732;
		w_1d[11] = 0.0;
		w_1d[12] = 0.0;
		w_1d[13] = 0.0;
		w_1d[14] = 0.0;
		w_1d[15] = 0.0;
		w_1d[16] = 0.0;
		w_1d[17] = 0.0;
		w_1d[18] = 0.0;
		w_1d[19] = 0.0;
	}

	if (GaussPt_1dir == 12)
	{
		G_1d[0]  = -0.9815606342467192;
		G_1d[1]  = -0.9041172563704748;
		G_1d[2]  = -0.7699026741943047;
		G_1d[3]  = -0.5873179542866175;
		G_1d[4]  = -0.3678314989981802;
		G_1d[5]  = -0.1252334085114689;
		G_1d[6]  =  0.1252334085114689;
		G_1d[7]  =  0.3678314989981802;
		G_1d[8]  =  0.5873179542866175;
		G_1d[9]  =  0.7699026741943047;
		G_1d[10] =  0.9041172563704748;
		G_1d[11] =  0.9815606342467192;
		G_1d[12] = 0.0;
		G_1d[13] = 0.0;
		G_1d[14] = 0.0;
		G_1d[15] = 0.0;
		G_1d[16] = 0.0;
		G_1d[17] = 0.0;
		G_1d[18] = 0.0;
		G_1d[19] = 0.0;

		w_1d[0]  = 0.0471753363865120;
		w_1d[1]  = 0.1069393259953189;
		w_1d[2]  = 0.1600783285433461;
		w_1d[3]  = 0.2031674267230656;
		w_1d[4]  = 0.2334925365383546;
		w_1d[5]  = 0.2491470458134027;
		w_1d[6]  = 0.2491470458134027;
		w_1d[7]  = 0.2334925365383546;
		w_1d[8]  = 0.2031674267230656;
		w_1d[9]  = 0.1600783285433461;
		w_1d[10] = 0.1069393259953189;
		w_1d[11] = 0.0471753363865120;
		w_1d[12] = 0.0;
		w_1d[13] = 0.0;
		w_1d[14] = 0.0;
		w_1d[15] = 0.0;
		w_1d[16] = 0.0;
		w_1d[17] = 0.0;
		w_1d[18] = 0.0;
		w_1d[19] = 0.0;
	}

	if (GaussPt_1dir == 13)
	{
		G_1d[0]  = -0.9841830547185881;
		G_1d[1]  = -0.9175983992229779;
		G_1d[2]  = -0.8015780907333099;
		G_1d[3]  = -0.6423493394403402;
		G_1d[4]  = -0.4484927510364468;
		G_1d[5]  = -0.2304583159551348;
		G_1d[6]  =  0.0000000000000000;
		G_1d[7]  =  0.2304583159551348;
		G_1d[8]  =  0.4484927510364468;
		G_1d[9]  =  0.6423493394403402;
		G_1d[10] =  0.8015780907333099;
		G_1d[11] =  0.9175983992229779;
		G_1d[12] =  0.9841830547185881;
		G_1d[13] = 0.0;
		G_1d[14] = 0.0;
		G_1d[15] = 0.0;
		G_1d[16] = 0.0;
		G_1d[17] = 0.0;
		G_1d[18] = 0.0;
		G_1d[19] = 0.0;

		w_1d[0]  = 0.0404840047653159;
		w_1d[1]  = 0.0921214998377286;
		w_1d[2]  = 0.1388735102197874;
		w_1d[3]  = 0.1781459807619455;
		w_1d[4]  = 0.2078160475368886;
		w_1d[5]  = 0.2262831802628971;
		w_1d[6]  = 0.2325515532308739;
		w_1d[7]  = 0.2262831802628971;
		w_1d[8]  = 0.2078160475368886;
		w_1d[9]  = 0.1781459807619455;
		w_1d[10] = 0.1388735102197874;
		w_1d[11] = 0.0921214998377286;
		w_1d[12] = 0.0404840047653159;
		w_1d[13] = 0.0;
		w_1d[14] = 0.0;
		w_1d[15] = 0.0;
		w_1d[16] = 0.0;
		w_1d[17] = 0.0;
		w_1d[18] = 0.0;
		w_1d[19] = 0.0;
	}


	if (GaussPt_1dir == 14)
	{
		G_1d[0]  = -0.9862838086968123;
		G_1d[1]  = -0.9284348836635735;
		G_1d[2]  = -0.8272013150697650;
		G_1d[3]  = -0.6872929048116855;
		G_1d[4]  = -0.5152486363581541;
		G_1d[5]  = -0.3191123689278897;
		G_1d[6]  = -0.1080549487073437;
		G_1d[7]  =  0.1080549487073437;
		G_1d[8]  =  0.3191123689278897;
		G_1d[9]  =  0.5152486363581541;
		G_1d[10] =  0.6872929048116855;
		G_1d[11] =  0.8272013150697650;
		G_1d[12] =  0.9284348836635735;
		G_1d[13] =  0.9862838086968123;
		G_1d[14] = 0.0;
		G_1d[15] = 0.0;
		G_1d[16] = 0.0;
		G_1d[17] = 0.0;
		G_1d[18] = 0.0;
		G_1d[19] = 0.0;

		w_1d[0]  = 0.0351194603317524;
		w_1d[1]  = 0.0801580871597603;
		w_1d[2]  = 0.1215185706879030;
		w_1d[3]  = 0.1572031671581934;
		w_1d[4]  = 0.1855383974779376;
		w_1d[5]  = 0.2051984637212955;
		w_1d[6]  = 0.2152638534631577;
		w_1d[7]  = 0.2152638534631577;
		w_1d[8]  = 0.2051984637212955;
		w_1d[9]  = 0.1855383974779376;
		w_1d[10] = 0.1572031671581934;
		w_1d[11] = 0.1215185706879030;
		w_1d[12] = 0.0801580871597603;
		w_1d[13] = 0.0351194603317524;
		w_1d[14] = 0.0;
		w_1d[15] = 0.0;
		w_1d[16] = 0.0;
		w_1d[17] = 0.0;
		w_1d[18] = 0.0;
		w_1d[19] = 0.0;
	}

	if (GaussPt_1dir == 15)
	{
		G_1d[0]  = -0.9879925180204854;
		G_1d[1]  = -0.9372733924007060;
		G_1d[2]  = -0.8482065834104272;
		G_1d[3]  = -0.7244177313601701;
		G_1d[4]  = -0.5709721726085388;
		G_1d[5]  = -0.3941513470775634;
		G_1d[6]  = -0.2011940939974345;
		G_1d[7]  =  0.0000000000000000;
		G_1d[8]  =  0.2011940939974345;
		G_1d[9]  =  0.3941513470775634;
		G_1d[10] =  0.5709721726085388;
		G_1d[11] =  0.7244177313601701;
		G_1d[12] =  0.8482065834104272;
		G_1d[13] =  0.9372733924007060;
		G_1d[14] =  0.9879925180204854;
		G_1d[15] = 0.0;
		G_1d[16] = 0.0;
		G_1d[17] = 0.0;
		G_1d[18] = 0.0;
		G_1d[19] = 0.0;

		w_1d[0]  = 0.0307532419961186;
		w_1d[1]  = 0.0703660474881081;
		w_1d[2]  = 0.1071592204671718;
		w_1d[3]  = 0.1395706779261539;
		w_1d[4]  = 0.1662692058169938;
		w_1d[5]  = 0.1861610000155619;
		w_1d[6]  = 0.1984314853271112;
		w_1d[7]  = 0.2025782419255609;
		w_1d[8]  = 0.1984314853271112;
		w_1d[9]  = 0.1861610000155619;
		w_1d[10] = 0.1662692058169938;
		w_1d[11] = 0.1395706779261539;
		w_1d[12] = 0.1071592204671718;
		w_1d[13] = 0.0703660474881081;
		w_1d[14] = 0.0307532419961186;
		w_1d[15] = 0.0;
		w_1d[16] = 0.0;
		w_1d[17] = 0.0;
		w_1d[18] = 0.0;
		w_1d[19] = 0.0;
	}

	if (GaussPt_1dir == 16)
	{
		G_1d[0]  = -0.9894009349916499;
		G_1d[1]  = -0.9445750230732326;
		G_1d[2]  = -0.8656312023878318;
		G_1d[3]  = -0.7554044083550030;
		G_1d[4]  = -0.6178762444026438;
		G_1d[5]  = -0.4580167776572274;
		G_1d[6]  = -0.2816035507792589;
		G_1d[7]  = -0.0950125098376375;
		G_1d[8]  =  0.0950125098376375;
		G_1d[9]  =  0.2816035507792589;
		G_1d[10] =  0.4580167776572274;
		G_1d[11] =  0.6178762444026438;
		G_1d[12] =  0.7554044083550030;
		G_1d[13] =  0.8656312023878318;
		G_1d[14] =  0.9445750230732326;
		G_1d[15] =  0.9894009349916499;
		G_1d[16] = 0.0;
		G_1d[17] = 0.0;
		G_1d[18] = 0.0;
		G_1d[19] = 0.0;

		w_1d[0]  = 0.0271524594117540;
		w_1d[1]  = 0.0622535239386477;
		w_1d[2]  = 0.0951585116824926;
		w_1d[3]  = 0.1246289712555340;
		w_1d[4]  = 0.1495959888165768;
		w_1d[5]  = 0.1691565193950026;
		w_1d[6]  = 0.1826034150449236;
		w_1d[7]  = 0.1894506104550686;
		w_1d[8]  = 0.1894506104550686;
		w_1d[9]  = 0.1826034150449236;
		w_1d[10] = 0.1691565193950026;
		w_1d[11] = 0.1495959888165768;
		w_1d[12] = 0.1246289712555340;
		w_1d[13] = 0.0951585116824926;
		w_1d[14] = 0.0622535239386477;
		w_1d[15] = 0.0271524594117540;
		w_1d[16] = 0.0;
		w_1d[17] = 0.0;
		w_1d[18] = 0.0;
		w_1d[19] = 0.0;
	}

	if (GaussPt_1dir == 17)
	{
		G_1d[0]  = -0.9905754753144174;
		G_1d[1]  = -0.9506755217687678;
		G_1d[2]  = -0.8802391537269859;
		G_1d[3]  = -0.7815140038968014;
		G_1d[4]  = -0.6576711592166908;
		G_1d[5]  = -0.5126905370864769;
		G_1d[6]  = -0.3512317634538763;
		G_1d[7]  = -0.1784841814958479;
		G_1d[8]  =  0.0000000000000000;
		G_1d[9]  =  0.1784841814958479;
		G_1d[10] =  0.3512317634538763;
		G_1d[11] =  0.5126905370864769;
		G_1d[12] =  0.6576711592166908;
		G_1d[13] =  0.7815140038968014;
		G_1d[14] =  0.8802391537269859;
		G_1d[15] =  0.9506755217687678;
		G_1d[16] =  0.9905754753144174;
		G_1d[17] = 0.0;
		G_1d[18] = 0.0;
		G_1d[19] = 0.0;

		w_1d[0]  = 0.0241483028685495;
		w_1d[1]  = 0.0554595293739866;
		w_1d[2]  = 0.0850361483171791;
		w_1d[3]  = 0.1118838471934036;
		w_1d[4]  = 0.1351363684685252;
		w_1d[5]  = 0.1540457610768101;
		w_1d[6]  = 0.1680041021564500;
		w_1d[7]  = 0.1765627053669925;
		w_1d[8]  = 0.1794464703562065;
		w_1d[9]  = 0.1765627053669925;
		w_1d[10] = 0.1680041021564500;
		w_1d[11] = 0.1540457610768101;
		w_1d[12] = 0.1351363684685252;
		w_1d[13] = 0.1118838471934036;
		w_1d[14] = 0.0850361483171791;
		w_1d[15] = 0.0554595293739866;
		w_1d[16] = 0.0241483028685495;
		w_1d[17] = 0.0;
		w_1d[18] = 0.0;
		w_1d[19] = 0.0;
	}

	if (GaussPt_1dir == 18)
	{
		G_1d[0]  = -0.9915651684209309;
		G_1d[1]  = -0.9558239495713978;
		G_1d[2]  = -0.8926024664975557;
		G_1d[3]  = -0.8037049589725231;
		G_1d[4]  = -0.6916870430603532;
		G_1d[5]  = -0.5597708310739475;
		G_1d[6]  = -0.4117511614628426;
		G_1d[7]  = -0.2518862256915055;
		G_1d[8]  = -0.0847750130417353;
		G_1d[9]  =  0.0847750130417353;
		G_1d[10] =  0.2518862256915055;
		G_1d[11] =  0.4117511614628426;
		G_1d[12] =  0.5597708310739475;
		G_1d[13] =  0.6916870430603532;
		G_1d[14] =  0.8037049589725231;
		G_1d[15] =  0.8926024664975557;
		G_1d[16] =  0.9558239495713978;
		G_1d[17] =  0.9915651684209309;
		G_1d[18] = 0.0;
		G_1d[19] = 0.0;

		w_1d[0]  = 0.0216160135264841;
		w_1d[1]  = 0.0497145488949692;
		w_1d[2]  = 0.0764257302548892;
		w_1d[3]  = 0.1009420441062870;
		w_1d[4]  = 0.1225552067114784;
		w_1d[5]  = 0.1406429146706506;
		w_1d[6]  = 0.1546846751262652;
		w_1d[7]  = 0.1642764837458327;
		w_1d[8]  = 0.1691423829631436;
		w_1d[9]  = 0.1691423829631436;
		w_1d[10] = 0.1642764837458327;
		w_1d[11] = 0.1546846751262652;
		w_1d[12] = 0.1406429146706506;
		w_1d[13] = 0.1225552067114784;
		w_1d[14] = 0.1009420441062870;
		w_1d[15] = 0.0764257302548892;
		w_1d[16] = 0.0497145488949692;
		w_1d[17] = 0.0216160135264841;
		w_1d[18] = 0.0;
		w_1d[19] = 0.0;
	}

	if (GaussPt_1dir == 19)
	{
		G_1d[0]  = -0.9924068438435844;
		G_1d[1]  = -0.9602081521348300;
		G_1d[2]  = -0.9031559036148179;
		G_1d[3]  = -0.8227146565371428;
		G_1d[4]  = -0.7209661773352294;
		G_1d[5]  = -0.6005453046616810;
		G_1d[6]  = -0.4645707413759609;
		G_1d[7]  = -0.3165640999636298;
		G_1d[8]  = -0.1603586456402254;
		G_1d[9]  =  0.0000000000000000;
		G_1d[10] =  0.1603586456402254;
		G_1d[11] =  0.3165640999636298;
		G_1d[12] =  0.4645707413759609;
		G_1d[13] =  0.6005453046616810;
		G_1d[14] =  0.7209661773352294;
		G_1d[15] =  0.8227146565371428;
		G_1d[16] =  0.9031559036148179;
		G_1d[17] =  0.9602081521348300;
		G_1d[18] =  0.9924068438435844;
		G_1d[19] = 0.0;

		w_1d[0]  = 0.0194617882297276;
		w_1d[1]  = 0.0448142267656998;
		w_1d[2]  = 0.0690445427376411;
		w_1d[3]  = 0.0914900216224498;
		w_1d[4]  = 0.1115666455473338;
		w_1d[5]  = 0.1287539625393362;
		w_1d[6]  = 0.1426067021736064;
		w_1d[7]  = 0.1527660420658594;
		w_1d[8]  = 0.1589688433939541;
		w_1d[9]  = 0.1610544498487834;
		w_1d[10] = 0.1589688433939541;
		w_1d[11] = 0.1527660420658594;
		w_1d[12] = 0.1426067021736064;
		w_1d[13] = 0.1287539625393362;
		w_1d[14] = 0.1115666455473338;
		w_1d[15] = 0.0914900216224498;
		w_1d[16] = 0.0690445427376411;
		w_1d[17] = 0.0448142267656998;
		w_1d[18] = 0.0194617882297276;
		w_1d[19] = 0.0;
	}

	if (GaussPt_1dir == 20)
	{
		G_1d[0]  = -0.9931285991850949;
		G_1d[1]  = -0.9639719272779138;
		G_1d[2]  = -0.9122344282513258;
		G_1d[3]  = -0.8391169718222188;
		G_1d[4]  = -0.7463319064601508;
		G_1d[5]  = -0.6360536807265150;
		G_1d[6]  = -0.5108670019508271;
		G_1d[7]  = -0.3737060887154195;
		G_1d[8]  = -0.2277858511416451;
		G_1d[9]  = -0.0765265211334973;
		G_1d[10] =  0.0765265211334973;
		G_1d[11] =  0.2277858511416451;
		G_1d[12] =  0.3737060887154195;
		G_1d[13] =  0.5108670019508271;
		G_1d[14] =  0.6360536807265150;
		G_1d[15] =  0.7463319064601508;
		G_1d[16] =  0.8391169718222188;
		G_1d[17] =  0.9122344282513258;
		G_1d[18] =  0.9639719272779138;
		G_1d[19] =  0.9931285991850949;

		w_1d[0]  = 0.0176140071391533;
		w_1d[1]  = 0.0406014298003862;
		w_1d[2]  = 0.0626720483341094;
		w_1d[3]  = 0.0832767415767047;
		w_1d[4]  = 0.1019301198172403;
		w_1d[5]  = 0.1181945319615182;
		w_1d[6]  = 0.1316886384491765;
		w_1d[7]  = 0.1420961093183819;
		w_1d[8]  = 0.1491729864726037;
		w_1d[9]  = 0.1527533871307258;
		w_1d[10] = 0.1527533871307258;
		w_1d[11] = 0.1491729864726037;
		w_1d[12] = 0.1420961093183819;
		w_1d[13] = 0.1316886384491765;
		w_1d[14] = 0.1181945319615182;
		w_1d[15] = 0.1019301198172403;
		w_1d[16] = 0.0832767415767047;
		w_1d[17] = 0.0626720483341094;
		w_1d[18] = 0.0406014298003862;
		w_1d[19] = 0.0176140071391533;
	}

	/*for (i = 0; i < GaussPt_3D; i++)
    {
		Gxi[i][0] = 0.0;
		Gxi[i][1] = 0.0;
		Gxi[i][2] = 0.0;
		w[i]      = 0.0;
	}*/

    ii = 0;

    for (i = 0; i < GaussPt_1dir; i++)
    {
        for (j = 0; j < GaussPt_1dir; j++)
        {
            for (k = 0; k < GaussPt_1dir; k++)
            {
                Gxi[ii][0] = G_1d[k];
                Gxi[ii][1] = G_1d[j];
                Gxi[ii][2] = G_1d[i];

                w[ii] = w_1d[i] * w_1d[j] * w_1d[k];
				//printf("やってるGxi[%d][0]:%lf\t[1]:%lf[2]:%lf\tw:%8.7e\tw_1d:%8.7e\t%8.7e\t%8.7e\n", ii, Gxi[ii][0],Gxi[ii][1],Gxi[ii][2],w[ii],w_1d[i],w_1d[j],w_1d[k]);

                ii++;
            }
        }
    }

	//printf("GaussPt_3D = %d\n", GaussPt_3D);
}

//拘束されている行数を省いた行列の番号の制作
int Make_Index_Dof(int Total_Control_Point, int Total_Constraint, int Constraint_Node_Dir[MAX_N_CONSTRAINT][2])
{
	int i, k = 0;

	for (i = 0; i < Total_Control_Point * DIMENSION; i++)
	{
		Index_Dof[i] = 0;
	}

	//拘束されている自由度(Degree Of free)をERRORにする
	for (i = 0; i < Total_Constraint; i++)
	{
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
	}
	//printf("\n\nMax_Index_Dof = %d\n", k);
	return k;
}

void Make_K_Whole_Ptr_Col(int Total_Element, int Total_Control_Point, int K_Whole_Size)
{
	int i, ii, j, jj, k;
	int NE;
	int N, i_index, j_index;

	//初期化
	//for (i = 0; i < Total_Control_Point * DIMENSION; i++)
		//Total_Control_Point_To_Node[i] = 0;
	for (i = 0; i < K_Whole_Size + 1; i++)
	{
		K_Whole_Ptr[i] = 0;
	}

	for (N = 0; N < Total_Control_Point; N += K_DIVISION_LENGE)
	{ //大きく分割するためのループ
		//各節点に接する節点を取得
		for (i = 0; i < K_DIVISION_LENGE; i++)
		{
			Total_Control_Point_To_Node[i] = 0;
		}

		for (i = 0; i < Total_Element; i++)
		{
			//printf("M_K_W_P_C---------------------------------------------------Element_No:%d---------------------------------------------------\n",i);
			for (ii = 0; ii < No_Control_point_ON_ELEMENT[Element_patch[i]]; ii++)
			{
				NE = Controlpoint_of_Element[i][ii] - N;
				//printf("K_DIVISION_LENGE=%d,N=%d,NE=%d\n",K_DIVISION_LENGE,N,NE);    //K_DIVISION_LENGE=0,N=0,NE=コネクティビティ的な
				if (0 <= NE && NE < K_DIVISION_LENGE)
				{
					for (j = 0; j < No_Control_point_ON_ELEMENT[Element_patch[i]]; j++)
					{
						//printf("M_K_W_P_C------------------NE=%d\tii=%d\tj=%d------------------\n", NE, ii, j);
						//数字がない時
						if (Total_Control_Point_To_Node[NE] == 0)
						{
							//節点番号を取得
							Node_To_Node[NE][0] = Controlpoint_of_Element[i][j];
							Total_Control_Point_To_Node[NE]++;
							//printf("i=%d,ii=%d,j=%d,Total_Control_Point_To_Node[ NE ] == 0のif文やった\n",i,ii,j);
							//printf("Node_To_Node[%d][0]=%d\n",NE,Node_To_Node[NE][0]);
							//printf("1Total_Control_Point_To_Node[%d]=%d\n",NE,Total_Control_Point_To_Node[NE]);
						}
						//printf("2Total_Control_Point_To_Node[%d]=%d\n",NE,Total_Control_Point_To_Node[NE]);
						//同じものがあったら
						for (k = 0; k < Total_Control_Point_To_Node[NE]; k++)
						{
							//printf("M_K_W_P_C--k=%d\n", k);
							if (Node_To_Node[NE][k] == Controlpoint_of_Element[i][j])
								break;
						}
						//printf("M_K_W_P_C------------------k=%d\tTotal_Control_Point_To_Node[%d]=%d------------------\n", k, NE, Total_Control_Point_To_Node[NE]);
						if (k == Total_Control_Point_To_Node[NE])
						{
							Node_To_Node[NE][k] = Controlpoint_of_Element[i][j];
							Total_Control_Point_To_Node[NE]++;
							//printf("Node_To_Node[%d][%d]=%d\n",NE,k,Node_To_Node[NE][k]);
							//printf("3Total_Control_Point_To_Node[%d]=%d\tk=%d\n",NE,Total_Control_Point_To_Node[NE],k);
						}
						//printf("\n");
					}

					//別メッシュとの重なりを考慮
					if (NNLOVER[i] > 0)
					{
						//printf("M_K_W_P_C------------------重なりがあった------------------\n");
						for (jj = 0; jj < NNLOVER[i]; jj++)
						{
							//printf("M_K_W_P_C_NELOVER[%d][%d]=%d\n", i, jj, NELOVER[i][jj]);
							for (j = 0; j < No_Control_point_ON_ELEMENT[Element_patch[NELOVER[i][jj]]]; j++)	//ローカル要素
							{
								//printf("M_K_W_P_C------------------NELOVER[%d][%d]=%d\tNE=%d\tii=%d\tj=%d------------------\n", i, jj, NELOVER[i][jj], NE, ii, j);
								//printf("j=%d\n",j);
								//数字がない時
								/*if (Total_Control_Point_To_Node[NE] == 0)
								{
									printf("いらないのでは?\n");
									//節点番号を取得
									Node_To_Node[NE][0] = Controlpoint_of_Element[NELOVER[i][jj]][j];
									Total_Control_Point_To_Node[NE]++;
									//printf("Node_To_Node[%d][0]=%d\n",NE,Node_To_Node[NE][0]);
								}*/
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
										//printf("break\n");
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
									//printf("SSS_Node_To_Node[%d][%d]=%d\n",NE,k,Node_To_Node[NE][k]);
									//printf("SSS_3Total_Control_Point_To_Node[%d]=%d\n",NE,Total_Control_Point_To_Node[NE]);
								}
							}
						}
						//printf("\n");
					}
				}
			}
			//printf("\n\n");
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
		}

		/*for (i = 0; i < 2181; i++)
		{
			printf("i:%d\t\t",i);
			for (j = 0; j < 1000; j++)
			{
				if (j == 0 || Node_To_Node[i][j] > Node_To_Node[i][j - 1])
				{
					printf("%d\t",Node_To_Node[i][j]);
				}
			}
			printf("\n");
		}*/

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
					k = 0;
					if (i_index >= 0)
					{
						K_Whole_Ptr[i_index + 1] = K_Whole_Ptr[i_index];
						for (j = 0; j < Total_Control_Point_To_Node[i]; j++)
						{
							for (jj = 0; jj < DIMENSION; jj++)
							{
								j_index = Index_Dof[Node_To_Node[i][j] * DIMENSION + jj];
								if (j_index >= 0 && j_index >= i_index)
								{
									K_Whole_Ptr[i_index + 1]++;
									K_Whole_Col[K_Whole_Ptr[i_index] + k] = j_index;
									k++;
									//printf("ptr[%d]=%d,col[%d]=%d\n",i_index+1,K_Whole_Ptr[i_index+1],K_Whole_Ptr[i_index]+k,K_Whole_Col[K_Whole_Ptr[i_index]+k]);
								}
							}
						}
					}
				}
			}
		}
		//col_N[N/K_DIVISION_LENGE][ k ] = -1;
	}
	/*
	for( i = 0; i < K_Whole_Size+1; i++ )//printf("K_Whole_Ptr[%d]= %d\n",i,K_Whole_Ptr[i]);
	//col合成
	k = 0;
	for( N = 0; N < Total_Control_Point ; N +=K_DIVISION_LENGE ){
		for(i = 0; col_N[ N/K_DIVISION_LENGE ][i] != -1; i++ ){
			K_Whole_Col[k] = col_N[ N/K_DIVISION_LENGE ][i];
			k++;
		}
	}*/
}
//valを求める
void Make_K_Whole_Val(double E, double nu, int Total_Element, int K_Whole_Size, int Total_Control_Point /*,int real_element[MAX_N_ELEMENT]*/)
{
	int i, j, j1, j2, k1, k2, l;
	int a, b, re;
	int n = 0;

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
		i = real_element[re];

		//ガウス点を作る．
		if(Element_mesh[i] == 0 && re == 0)/*2つめの条件は効率化のため*/
		{
			Make_Gauss_points(0);
		}

		if(Element_mesh[i] > 0)
		{
			printf("NNLOVER[%d]:%d\tNNLOVER[%d]:%d\tElement_mesh[%d]:%d\n", i, NNLOVER[i], real_element[re - 1], NNLOVER[real_element[re - 1]], real_element[re - 1], Element_mesh[real_element[re - 1]]);
			if(NNLOVER[i] == 1 && (NNLOVER[real_element[re - 1]] != 1 || Element_mesh[real_element[re - 1]] == 0))/*2つめ以降の条件は効率化のため*/
			{
				Make_Gauss_points(0);
			}
			if(NNLOVER[i] >= 2 && (NNLOVER[real_element[re - 1]] == 1 || Element_mesh[real_element[re - 1]] == 0))/*2つめ以降の条件は効率化のため*/
			{
				Make_Gauss_points(1);
			}
		}
		//printf("i= %d\tGaussPt_3D=%d\n",i ,GaussPt_3D);

		for (j = 0; j < No_GP * No_GP * No_GP; j++)
		{
			Same_BDBJ_flag[j] = 0;
		}
		KIEL_SIZE = No_Control_point_ON_ELEMENT[Element_patch[i]] * DIMENSION;
		double X[No_Control_point_ON_ELEMENT[Element_patch[i]]][DIMENSION], K_EL[KIEL_SIZE][KIEL_SIZE];
		//printf("re=%d\n", re);
		//printf("El_No;i=%d\tElement_mesh[%d] = %d\n",real_element[re] ,i ,Element_mesh[i]);
		//printf("El_No=%d\n",i );
		//各要素のKelを求める
		for (j1 = 0; j1 < No_Control_point_ON_ELEMENT[Element_patch[i]]; j1++)
		{
			for (j2 = 0; j2 < DIMENSION; j2++)
			{
				X[j1][j2] = Node_Coordinate[Controlpoint_of_Element[i][j1]][j2];
			}
		}
		Make_K_EL(i, X, K_EL, E, nu, Total_Element, Total_Control_Point);
		printf("El_No:%d Element_patch:%d J_test = %e  Check_Volume = %.13e\n", re, Element_patch[re], J_test, Check_Volume);
		//printf("i= %d\tNNLOVER=%d\tGaussPt_3D=%d\n",i ,NNLOVER[i] ,GaussPt_3D);

		//Valを求める
		for (j1 = 0; j1 < No_Control_point_ON_ELEMENT[Element_patch[i]]; j1++)
		{
			for (j2 = 0; j2 < DIMENSION; j2++)
			{
				a = Index_Dof[Controlpoint_of_Element[i][j1] * DIMENSION + j2];
				if (a >= 0)
				{
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
										K_Whole_Val[l] += K_EL[j1 * DIMENSION + j2][k1 * DIMENSION + k2];
										//	printf("j1\t%d\tj2\t%d\ta\t%d\tk1\t%d\tk2\t%d\tb\t%d\t%d\t%22.17lf\t%22.17lf\n",j1,j2,a,k1,k2,b,l,K_Whole_Val[l], K_EL[j1 * DIMENSION + j2][k1 * DIMENSION + k2]);

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
				printf("---------------------------------------------------------------------------------------------Local_Element_No:%d---------------------------------------------------------------------------------------------\n",i);
				for (j = 0; j < NNLOVER[i]; j++)
				{
					printf("---------------------------------------------------Global_Element_No:%d(Loc:%d)---------------------------------------------------\n",NELOVER[i][j],i);
					double XG[No_Control_point_ON_ELEMENT[Element_patch[NELOVER[i][j]]]][DIMENSION];
					double coupled_K_EL[KIEL_SIZE][KIEL_SIZE];
					for (j1 = 0; j1 < No_Control_point_ON_ELEMENT[Element_patch[NELOVER[i][j]]]; j1++)
					{
						for (j2 = 0; j2 < DIMENSION; j2++)
						{
							XG[j1][j2] = Node_Coordinate[Controlpoint_of_Element[NELOVER[i][j]][j1]][j2];
							//重なっている要素の物理座標取得
						}
					}
					Make_coupled_K_EL(i, NELOVER[i][j], X, XG, coupled_K_EL, E, nu);

					Check_BDBJ_flag[i] += Total_BDBJ_flag;
					if (j == NNLOVER[i] - 1)
					{
						for (j1 = 0; j1 < GaussPt_3D; j1++)
						{
							//printf("Same_BDBJ_flag[%d]=%d\n",j1,Same_BDBJ_flag[j1]);
							if (Same_BDBJ_flag[j1] != 1)
							{
								printf("ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR-ERROR\n");
							}
						}
						printf("-------------------------Check_BDBJ_flag[%d]=%d-------------------------\n",i ,Check_BDBJ_flag[i]);
						if (Check_BDBJ_flag[i] != GaussPt_3D)
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
										if (b >= 0 && b > a) //イコールを削除
										{
											//printf("Glo_ele:%d\tCont_No_a:%d\ta:%d\t\tLoc_ele:%d\tCont_No_b:%d\tb:%d\n",NELOVER[i][j], Controlpoint_of_Element[NELOVER[i][j]][j1], a, i, Controlpoint_of_Element[i][k1], b);
											for (l = K_Whole_Ptr[a]; l < K_Whole_Ptr[a + 1]; l++)
											{
												if (K_Whole_Col[l] == b)
												{
													K_Whole_Val[l] += coupled_K_EL[j1 * DIMENSION + j2][k1 * DIMENSION + k2];
													//printf("j1\t%d\tj2\t%d\ta\t%d\tk1\t%d\tk2\t%d\tb\t%d\t%d\t%22.17lf\t%22.17lf\n",j1,j2,a,k1,k2,b,l,K_Whole_Val[l], coupled_K_EL[j1 * DIMENSION + j2][k1 * DIMENSION + k2]);
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
	}

	for (i = 0; i < MAX_NON_ZERO; i++)
	{
	//	for (j = 0; j < K_Whole_Size ; j++)
	//	{
	//		if (i == K_Whole_Ptr[j])
	//		{
	//			printf("j = %d\t",j);
	//		}
	//	}
		if (K_Whole_Val[i] != 0.0)
		{
	//		printf("K_val[%d]=%20.15lf\n",i, K_Whole_Val[i]);
			n++;
		}
	}

	printf("DOF : %d\tThe number of non zezo : %d\n", K_Whole_Size, n);
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
	{
		rhs_vec[i] = 0.0;
	}
	for (i = 0; i < Total_Load; i++)
	{
		index = Index_Dof[Load_Node_Dir[i][0] * DIMENSION + Load_Node_Dir[i][1]];
		if (index >= 0)
			rhs_vec[index] += Value_of_Load[i];
	}
}
//強制変位対策
void Make_F_Vec_disp_const(int Mesh_No, int Total_Constraint, int Constraint_Node_Dir[MAX_N_CONSTRAINT][2], double Value_of_Constraint[MAX_N_CONSTRAINT], int Total_Element, double E, double nu, /*int DM,*/ int Total_Control_Point)
{
	int ie, iee, idir, inode, jdir, jnode, kk_const;
	int ii, iii, b, bb, jj, j1, j2, ii_local, jj_local;

	double K_EL[KIEL_SIZE][KIEL_SIZE];

	Make_Gauss_points(0);

	for (ie = 0; ie < real_Total_Element; ie++)
	{

		double X[No_Control_point_ON_ELEMENT[Element_patch[real_El_No_on_mesh[Mesh_No][ie]]]][DIMENSION];
		iii = 0;
		for (idir = 0; idir < DIMENSION; idir++)
		{
			for (inode = 0; inode < No_Control_point_ON_ELEMENT[Element_patch[real_El_No_on_mesh[Mesh_No][ie]]]; inode++)
			{
				//b = Index_Dof[Controlpoint_of_Element[ie ][inode] * DIMENSION + idir];
				b = Index_Dof[Controlpoint_of_Element[real_El_No_on_mesh[Mesh_No][ie]][inode] * DIMENSION + idir];
				if (b < 0)
					iii++;
			}
		}
		//printf("iii;%d\n",iii);
		if (iii > 0)
		{
			for (j1 = 0; j1 < No_Control_point_ON_ELEMENT[Element_patch[real_El_No_on_mesh[Mesh_No][ie]]]; j1++)
			{
				for (j2 = 0; j2 < DIMENSION; j2++)
				{
					//X[j1][j2] = Node_Coordinate[Controlpoint_of_Element[ie ][j1]][j2];
					X[j1][j2] = Node_Coordinate[Controlpoint_of_Element[real_El_No_on_mesh[Mesh_No][ie]][j1]][j2];
				} //end for j2
			}	 //end for j1
			iee = real_El_No_on_mesh[Mesh_No][ie];
			Make_K_EL(iee, X, K_EL, E, nu, Total_Element, Total_Control_Point);
			for (idir = 0; idir < DIMENSION; idir++)
			{
				for (inode = 0; inode < No_Control_point_ON_ELEMENT[Element_patch[real_El_No_on_mesh[Mesh_No][ie]]]; inode++)
				{
					//ii = Controlpoint_of_Element[ie ][inode] * DIMENSION + idir;
					ii = Controlpoint_of_Element[real_El_No_on_mesh[Mesh_No][ie]][inode] * DIMENSION + idir;
					b = Index_Dof[ii];
					if (b >= 0)
					{
						ii_local = inode * DIMENSION + idir;
						for (jdir = 0; jdir < DIMENSION; jdir++)
						{
							for (jnode = 0; jnode < No_Control_point_ON_ELEMENT[Element_patch[real_El_No_on_mesh[Mesh_No][ie]]]; jnode++)
							{
								//jj = Controlpoint_of_Element[ie ][jnode] * DIMENSION + jdir;
								jj = Controlpoint_of_Element[real_El_No_on_mesh[Mesh_No][ie]][jnode] * DIMENSION + jdir;
								bb = Index_Dof[jj];
								if (bb < 0)
								{
									jj_local = jnode * DIMENSION + jdir; //printf("%d,%d\n",ie,jnode);
									for (kk_const = 0; kk_const < Total_Constraint; kk_const++)
									{
										//if (Controlpoint_of_Element[ie ][jnode] == Constraint_Node_Dir[kk_const][0] && jdir == Constraint_Node_Dir[kk_const][1])
										if (Controlpoint_of_Element[real_El_No_on_mesh[Mesh_No][ie]][jnode] == Constraint_Node_Dir[kk_const][0] && jdir == Constraint_Node_Dir[kk_const][1])
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
	{
		vec_result[i] = 0.0;
	}

	for (i = 0; i < ndof; i++)
	{
		for (j = K_Whole_Ptr[i]; j < K_Whole_Ptr[i + 1]; j++)
		{
			vec_result[i] += K_Whole_Val[icount] * vec[K_Whole_Col[j]];
			if (i != K_Whole_Col[j])
			{
				vec_result[K_Whole_Col[j]] += K_Whole_Val[icount] * vec[i];
			}
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
		//printf("vec1[%d]=%f vec2[%d]=%f\n",i,vec1[i],i,vec2[i]);
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
		//if(i == ndof - 1)
		//{
		//	printf("-Not_converged-Not_converged-Not_converged-Not_converged-Not_converged-Not_converged-Not_converged-Not_converged-Not_converged\n");
		//}
	}
	rrr3 = fabs(alphak) * sqrt(rrr1 / rrr2);
	printf("Iteration# = %d  residual = %15e (%15e)\n", itr, rrr3, eps);
	if (rrr3 < eps)
		istop = 1;

	return (istop);
}

void Diag_Scaling_CG_pre(int ndof, int flag_operation)
{
	int i, j;
	int icount = 0;
	/* flag_opertion = 0: Pre process to the CG solver
			A <-- Dt A D  and b <-- Dt b */
	/* flag_operation = 1: Post process to the CG solver
			b <-- Dt b  */
	if (flag_operation == 0)
	{
		diag_scaling[0] = 1.0 / sqrt(K_Whole_Val[0]);
		/* diag_scaling[0] = 1.0; */
		for (i = 1; i < ndof; i++)
		{
			//	printf("i = %d\tptr = %d val = %le\n", i, K_Whole_Ptr[i], K_Whole_Val[K_Whole_Ptr[i]]); //あとでコメントアウト

			diag_scaling[i] = 1.0 / sqrt(K_Whole_Val[K_Whole_Ptr[i]]);
			/* diag_scaling[i] = 1.0; */
		}
		for (i = 0; i < ndof; i++)
		{
			for (j = K_Whole_Ptr[i]; j < K_Whole_Ptr[i + 1]; j++)
			{
				//printf("Check scling icount=%d i=%d K_Whole_Col[%d] = %d\n",icount,i,j,K_Whole_Col[j]);
				K_Whole_Val[icount] = K_Whole_Val[icount] * diag_scaling[i] * diag_scaling[K_Whole_Col[j]];
				//printf("K_Whole_Val = %f\n",K_Whole_Val[icount]);
				icount++;
			}
			printf("rhs_vec_before[%d]:%le diag_scaling[%d]:%le\n", i, rhs_vec[i], i, diag_scaling[i]);
			rhs_vec[i] = rhs_vec[i] * diag_scaling[i];
			//printf("rhs_vec[%d]:%le\n",i,rhs_vec[i]);
		}
	}
	if (flag_operation == 1)
	{
		for (i = 0; i < ndof; i++)
		{
			//printf("solvec[%d] = %f\n",i, sol_vec[i]);
			sol_vec[i] = sol_vec[i] * diag_scaling[i];
		}
	}

	//printf("\nqq\n");
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
		gg[i] = rhs_vec[i] - dd[i];
		printf("rhs_vec[%d] = %f dd[%d] = %f\n", i, rhs_vec[i], i, dd[i]);
		pp[i] = gg[i];
	}

	//printf("\nrr");

	for (itr = 0; itr < max_itr; itr++)
	{
		ppp = inner_product(ndof, gg, gg);
		mat_vec_crs(dd, pp, ndof);
		rrr = inner_product(ndof, dd, pp);
		alphak = ppp / rrr;
		//printf("ppp=%f rrr=%f\n", ppp, rrr);
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

	//printf("\nss");
}
////////////////////////////////////////////////////////////////////////
/////////////////要素剛性マトリックス
////////////////////////////////////////////////////////////////////////
//IGAの基底関数
void ShapeFunction1D(double Position_Data_param[DIMENSION], int j, int e) //要チェック
{
	int ii;
    int p;
    //double err1, err2;



	for (ii = 0; ii < No_knot[Element_patch[e]][j]; ii++)
	{
		if (Position_Knots[Element_patch[e]][j][ii] == Position_Knots[Element_patch[e]][j][ii + 1])
		{
			Shape[j][ii][0] = 0.0;
		}
		else if (Position_Knots[Element_patch[e]][j][ii] != Position_Knots[Element_patch[e]][j][ii + 1] &&
				 (Position_Knots[Element_patch[e]][j][ii] < Position_Data_param[j] || fabs(Position_Knots[Element_patch[e]][j][ii] - Position_Data_param[j]) < ERR) &&
				 Position_Data_param[j] < Position_Knots[Element_patch[e]][j][ii + 1])
		{
			Shape[j][ii][0] = 1.0;
		}
		else if (Position_Knots[Element_patch[e]][j][ii] != Position_Knots[Element_patch[e]][j][ii + 1] &&
				 Position_Knots[Element_patch[e]][j][ii + 1] == Position_Knots[Element_patch[e]][j][(No_knot[Element_patch[e]][j] - 1)] &&
				 ( Position_Knots[Element_patch[e]][j][ii] < Position_Data_param[j] || fabs(Position_Knots[Element_patch[e]][j][ii] - Position_Data_param[j]) < ERR ) &&
				 ( Position_Data_param[j] < Position_Knots[Element_patch[e]][j][ii + 1] || fabs(Position_Knots[Element_patch[e]][j][ii + 1] - Position_Data_param[j]) < ERR ))
		{
			Shape[j][ii][0] = 1.0;
		}
		else
		{
			Shape[j][ii][0] = 0.0;
		}
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

			if (fabs((Position_Data_param[j] - Position_Knots[Element_patch[e]][j][ii]) * Shape[j][ii][p - 1]) < ERR &&
				Position_Knots[Element_patch[e]][j][ii + p] == Position_Knots[Element_patch[e]][j][ii])
			{
				left_term = 0.0;
			}
			else
            {
				left_term = (Position_Data_param[j] - Position_Knots[Element_patch[e]][j][ii]) / (Position_Knots[Element_patch[e]][j][ii + p] - Position_Knots[Element_patch[e]][j][ii]) * Shape[j][ii][p - 1];
                //printf("else left\tleft_term;%le\n",left_term);
            }
			if (fabs((Position_Knots[Element_patch[e]][j][ii + p + 1] - Position_Data_param[j]) * Shape[j][ii + 1][p - 1]) < ERR &&
				Position_Knots[Element_patch[e]][j][ii + p + 1] == Position_Knots[Element_patch[e]][j][ii + 1])
			{
				right_term = 0.0;
			}
			else
            {
				right_term = (Position_Knots[Element_patch[e]][j][ii + p + 1] - Position_Data_param[j]) / (Position_Knots[Element_patch[e]][j][ii + p + 1] - Position_Knots[Element_patch[e]][j][ii + 1]) * Shape[j][ii + 1][p - 1];
                //printf("else right\tright_term;%le\n",right_term);
            }
			Shape[j][ii][p] = left_term + right_term;
			/*if(j == 1 && p == 2)
			{
				printf("Shape[%d][%d][%d]=%le\n",j,ii,p,Shape[j][ii][p]);
			}*/
			//if (e == 21)
			//{
			//	printf("要素境界Position_Data_param[%d]:%lf\tShape[%d][%d][%d]:%lf\n",j,Position_Data_param[j],j,ii,p,Shape[j][ii][p]);
			//}
		}
		//if (e == 21)
		//{
		//	printf("\n");
		//}
	}
	//printf("order[%d]:%d\n",j,Order[Element_patch[e]][j] );
	double dleft_term, dright_term;
	for (ii = 0; ii < No_Control_point[Element_patch[e]][j] + 1; ii++)
	{
		//printf("No_Control_point[%d]=%d\n",j,No_Control_point[j] );
		dleft_term = 0.0;
		dright_term = 0.0;

		if (fabs(Order[Element_patch[e]][j] * Shape[j][ii][Order[Element_patch[e]][j] - 1]) < ERR &&
			Position_Knots[Element_patch[e]][j][ii + Order[Element_patch[e]][j]] == Position_Knots[Element_patch[e]][j][ii])
		{
			dleft_term = 0.0;
		}
		else
		{
			dleft_term = Order[Element_patch[e]][j] / (Position_Knots[Element_patch[e]][j][ii + Order[Element_patch[e]][j]] - Position_Knots[Element_patch[e]][j][ii]) * Shape[j][ii][Order[Element_patch[e]][j] - 1];
		}
		/*printf("test_Shape_left[%d][%d][%d]=%le\n", j,ii,Order[Element_patch[e]][j]-1,Shape[j][ii][Order[Element_patch[e]][j]-1]);
		printf("Position_Knots[Element_patch[e]][%d][%d]:%le\n", j,ii+Order[Element_patch[e]][j],Position_Knots[Element_patch[e]][j][ii+Order[Element_patch[e]][j]]);
		printf("Position_Knots[Element_patch[e]][%d][%d]:%le\n", j,ii,Position_Knots[Element_patch[e]][j][ii]);
		printf("dleft_term=%f\n",dleft_term );*/

		if (fabs(Order[Element_patch[e]][j] * Shape[j][ii + 1][Order[Element_patch[e]][j] - 1]) < ERR &&
			Position_Knots[Element_patch[e]][j][ii + Order[Element_patch[e]][j] + 1] == Position_Knots[Element_patch[e]][j][ii + 1])
		{
			dright_term = 0.0;
		}
		else
		{
			dright_term = Order[Element_patch[e]][j] / (Position_Knots[Element_patch[e]][j][ii + Order[Element_patch[e]][j] + 1] - Position_Knots[Element_patch[e]][j][ii + 1]) * Shape[j][ii + 1][Order[Element_patch[e]][j] - 1];
		}
		/*printf("test_Shape_right[%d][%d][%d]=%le\n", j,ii+1,Order[Element_patch[e]][j]-1,Shape[j][ii+1][Order[Element_patch[e]][j]-1]);
		printf("Position_Knots[%d][%d]:%le\n", j,ii+Order[Element_patch[e]][j]+1,Position_Knots[j][ii+Order[Element_patch[e]][j]+1]);
		printf("Position_Knots[%d][%d]:%le\n", j,ii+1,Position_Knots[j][ii+1]);
		printf("dright_term=%f\n",dright_term );*/

		dShape[j][ii] = dleft_term - dright_term;

		//if (e == 21)
		//{
		//	printf("要素境界Position_Data_param[%d]:%lf\tdShape[%d][%d]:%lf\n",j,Position_Data_param[j],j,ii,dShape[j][ii]);
		//}

		//printf("PP=%d\n",PP );

		//printf("dShape[1][%d]= %f\n",ii,dShape[1][ii]);
	}
}

void ShapeFunc_from_paren(double Local_coord[DIMENSION], int j, int e)
{

	int i = 0;

	//printf("Local_coord[%d]:%le\n",j,Local_coord[j]);
	//i = INC[Element_patch[e]][Controlpoint_of_Element[e][0]][j];
    i = INC[Element_patch[e]][Controlpoint_of_Element[e][0]][j];
	//printf("El_No:%d\n",e );
	//printf("i:%d\n",i);
	//printf("Position_Knots[%d][%d]:%le Position_Knots[%d][%d]:%le\n", j,i+1,Position_Knots[j][i+1],j,i,Position_Knots[j][i]);
	Position_Data_param[j] = ((Position_Knots[Element_patch[e]][j][i + 1] - Position_Knots[Element_patch[e]][j][i]) * Local_coord[j] + (Position_Knots[Element_patch[e]][j][i + 1] + Position_Knots[Element_patch[e]][j][i])) / 2.0;
	//Position_Data_param[j] = ((Position_Knots[Element_patch[e]][j][i + 1] - Position_Knots[Element_patch[e]][j][i]) * Local_coord[j] + (Position_Knots[Element_patch[e]][j][i + 1] + Position_Knots[Element_patch[e]][j][i])) / 2;
	//printf("Position_Data_param[%d]:%le\n",j,Position_Data_param[j]);
	if (fabs(Position_Data_param[j] - Position_Knots[Element_patch[e]][j][i + 1]) < ERR)
	{
		Position_Data_param[j] = Position_Knots[Element_patch[e]][j][i + 1];
		//printf("修正Position_Data_param[%d]=%20.19lf\ti+1\n",j, Position_Data_param[j]);
	}
	else if (fabs(Position_Data_param[j] - Position_Knots[Element_patch[e]][j][i]) < ERR)
	{
		Position_Data_param[j] = Position_Knots[Element_patch[e]][j][i];
		//printf("修正Position_Data_param[%d]=%20.19lf\ti\n",j, Position_Data_param[j]);
	}
}

double dShapeFunc_from_paren(int j, int e)
{
	int i;
	double dPosition_Data_param;

	//i = INC[Element_patch[e]][Controlpoint_of_Element[e][0]][j];
	i = INC[Element_patch[e]][Controlpoint_of_Element[e][0]][j];
	//printf("El_No:%d\n",e );
	//printf("i:%d\n",i);
	//printf("Position_Knots[%d][%d][%d]:%le Position_Knots[%d][%d][%d]:%le\n",Element_patch[e], j,i+1,Position_Knots[Element_patch[e]][j][i+1],Element_patch[e],j,i,Position_Knots[Element_patch[e]][j][i]);
	dPosition_Data_param = (Position_Knots[Element_patch[e]][j][i + 1] - Position_Knots[Element_patch[e]][j][i]) / 2.0;
	//printf("dPosition_Data_param:%le\t",dPosition_Data_param);
	return dPosition_Data_param;
}

double Shape_func(int I_No, double Local_coord[DIMENSION], int El_No)
{

	int i, j;
	int No_control_point = 0;
	double R;
	double weight_func;

	weight_func = 0.0;


	for (i = 0; i < No_Input_File; i++)
	{
		No_control_point += Total_Control_Point_on_mesh[i];
	}
	//printf("sssssssssssssssssss%dsssssssssssssssss\n", No_control_point);
	for (i = 0; i < No_control_point; i++)
	{
		shape_func[i] = 1.0;
	}

	for (j = 0; j < DIMENSION; j++)
	{
		ShapeFunc_from_paren(Local_coord, j, El_No);
		ShapeFunction1D(Position_Data_param, j, El_No);
		for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++)
		{
			shape_func[Controlpoint_of_Element[El_No][i]] *= Shape[j][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][j]][Order[Element_patch[El_No]][j]];
		}
	}

	for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++)
	{
		weight_func += shape_func[Controlpoint_of_Element[El_No][i]] * Node_Coordinate[Controlpoint_of_Element[El_No][i]][DIMENSION];
	}

	if (I_No < No_Control_point_ON_ELEMENT[Element_patch[El_No]])
	{
		R = shape_func[Controlpoint_of_Element[El_No][I_No]] * Node_Coordinate[Controlpoint_of_Element[El_No][I_No]][DIMENSION] / weight_func;
	}
	else
	{
		R = ERROR;
	}
	return R;
}

void NURBS_deriv(double Local_coord[DIMENSION], int El_No)
{
	double weight_func;

	double dWeight_func1;
	double dWeight_func2;
	double dWeight_func3;

	int i, j;
	int No_control_point = 0;

	//for(ii = 0; ii < NN+1; ii++){
	//printf("NdShape3[%d]= %f\n",ii,dShape3[ii]);
	//}

	//for (ii = 0; ii < NN+1; ii++)printf("NdShape1[%d]= %f\n",ii,dShape1[ii]);
	//for (jj = 0; jj < MM+1; jj++)printf("NdShape2[%d]= %f\n",jj,dShape2[jj]);
	//printf("\n");

	for (i = 0; i < No_Input_File; i++)
	{
		No_control_point += Total_Control_Point_on_mesh[i];
	}

	for (i = 0; i < No_control_point; i++)
	{
		shape_func[i] = 1.0;
	}

	weight_func = 0.0;

	dWeight_func1 = 0.0;
	dWeight_func2 = 0.0;
	dWeight_func3 = 0.0;


	for (j = 0; j < DIMENSION; j++)
	{
		ShapeFunc_from_paren(Local_coord, j, El_No);
		ShapeFunction1D(Position_Data_param, j, El_No);
		for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++)
		{
			shape_func[Controlpoint_of_Element[El_No][i]] *= Shape[j][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][j]][Order[Element_patch[El_No]][j]];

			//if (El_No == 21)
			//{
			//	printf("要素境界i=%d\tShape[%d][%d][%d]:%lf\n",i,j,INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][j],Order[Element_patch[El_No]][j],Shape[j][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][j]][Order[Element_patch[El_No]][j]]);
			//}
		}
	}
	//if (El_No == 21)
	//{
	//	for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++)
	//	{
	//		printf("要素境界i=%d\tshape_func[%d]=%lf\n",i,Controlpoint_of_Element[El_No][i],shape_func[Controlpoint_of_Element[El_No][i]]);
	//	}
	//}

	for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++)
	{
		weight_func += shape_func[Controlpoint_of_Element[El_No][i]] * Node_Coordinate[Controlpoint_of_Element[El_No][i]][DIMENSION];
		//printf("Node_Coordinate[%d][%d]:%le\n", Controlpoint_of_Element[El_No][i],DIMENSION,Node_Coordinate[Controlpoint_of_Element[El_No][i]][DIMENSION]);
		//printf("shape_func[%d]:%le\n",Controlpoint_of_Element[El_No][i],shape_func[Controlpoint_of_Element[El_No][i]]);
		//printf("weight_func:%le\n", weight_func);
	}

	for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++)
	{
		dWeight_func1 +=   dShape[0][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][0]]
						 *  Shape[1][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][1]][Order[Element_patch[El_No]][1]]
						 *  Shape[2][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][2]][Order[Element_patch[El_No]][2]]
						 *  Node_Coordinate[Controlpoint_of_Element[El_No][i]][DIMENSION];

		dWeight_func2 +=    Shape[0][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][0]][Order[Element_patch[El_No]][0]]
						 * dShape[1][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][1]]
						 *  Shape[2][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][2]][Order[Element_patch[El_No]][2]]
						 *  Node_Coordinate[Controlpoint_of_Element[El_No][i]][DIMENSION];

		dWeight_func3 +=    Shape[0][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][0]][Order[Element_patch[El_No]][0]]
						 *  Shape[1][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][1]][Order[Element_patch[El_No]][1]]
						 * dShape[2][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][2]]
						 *  Node_Coordinate[Controlpoint_of_Element[El_No][i]][DIMENSION];
	}
	//if (Local_coord[1] == -1.0 || Local_coord[1] == 1.0)
	//{
		//printf("dWeight_func1:%le dWeight_func2:%le dWeight_func3:%le\n",dWeight_func1,dWeight_func2,dWeight_func3);
	//}
	for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++)
	{
		dShape_func1[Controlpoint_of_Element[El_No][i]]
			= Node_Coordinate[Controlpoint_of_Element[El_No][i]][DIMENSION]
			 * ( weight_func * dShape[0][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][0]] * Shape[1][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][1]][Order[Element_patch[El_No]][1]] * Shape[2][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][2]][Order[Element_patch[El_No]][2]] - dWeight_func1 * shape_func[Controlpoint_of_Element[El_No][i]] ) / (weight_func * weight_func);

		dShape_func2[Controlpoint_of_Element[El_No][i]]
			= Node_Coordinate[Controlpoint_of_Element[El_No][i]][DIMENSION]
			 * ( weight_func * Shape[0][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][0]][Order[Element_patch[El_No]][0]] * dShape[1][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][1]] * Shape[2][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][2]][Order[Element_patch[El_No]][2]] - dWeight_func2 * shape_func[Controlpoint_of_Element[El_No][i]] ) / (weight_func * weight_func);

		dShape_func3[Controlpoint_of_Element[El_No][i]]
			= Node_Coordinate[Controlpoint_of_Element[El_No][i]][DIMENSION]
			 * ( weight_func * Shape[0][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][0]][Order[Element_patch[El_No]][0]] * Shape[1][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][1]][Order[Element_patch[El_No]][1]] * dShape[2][INC[Element_patch[El_No]][Controlpoint_of_Element[El_No][i]][2]] - dWeight_func3 * shape_func[Controlpoint_of_Element[El_No][i]] ) / (weight_func * weight_func);
		//if (El_No == 21)
		//{
		//	//printf("NURBS_deriv;Controlpoint_of_Element[%d][%d]:%d\n",El_No,i,Controlpoint_of_Element[El_No][i]);
		//	printf("要素境界i=%d\tdShape_func1[%d]=%le\t2=%le\t3=%le\n",i,Controlpoint_of_Element[El_No][i],dShape_func1[Controlpoint_of_Element[El_No][i]],dShape_func2[Controlpoint_of_Element[El_No][i]],dShape_func3[Controlpoint_of_Element[El_No][i]]);
		//	//printf("要素境界i=%d\tdShape_func2[%d]=%le\n",i,Controlpoint_of_Element[El_No][i],dShape_func2[Controlpoint_of_Element[El_No][i]]);
		//	//printf("要素境界i=%d\tdShape_func3[%d]=%le\n",i,Controlpoint_of_Element[El_No][i],dShape_func3[Controlpoint_of_Element[El_No][i]]);
		//}
	}
}

double dShape_func(int I_No, int xez, double Local_coord[DIMENSION], int El_No)
{
	double dR;

	//printf("El_No=%d\n",El_No );

	NURBS_deriv(Local_coord, El_No);

	//printf("No_Control_point_ON_ELEMENT[%d]=%d\n",Element_patch[El_No], No_Control_point_ON_ELEMENT[Element_patch[El_No]]);

	if (xez != 0 && xez != 1 && xez != 2)
	{
		dR = ERROR;
	}
	else if (I_No < No_Control_point_ON_ELEMENT[Element_patch[El_No]])
	{
		if (xez == 0)
		{
			//dR = dShape_func1[Controlpoint_of_Element[El_No][I_No]] * dShapeFunc_from_paren(xez, El_No);
			dR = dShape_func1[Controlpoint_of_Element[El_No][I_No]] * dShapeFunc_from_paren(xez, El_No);
			//printf("dShape_func1[%d]:%le\n",Controlpoint_of_Element[El_No][I_No],dShape_func1[Controlpoint_of_Element[El_No][I_No]]);
		}
		if (xez == 1)
		{
			//dR = dShape_func2[Controlpoint_of_Element[El_No][I_No]] * dShapeFunc_from_paren(xez, El_No);
			dR = dShape_func2[Controlpoint_of_Element[El_No][I_No]] * dShapeFunc_from_paren(xez, El_No);
			//printf("dShape_func2[%d]:%le\n",Controlpoint_of_Element[El_No][I_No],dShape_func2[Controlpoint_of_Element[El_No][I_No]]);
		}
		else if (xez == 2)
		{
			//dR = dShape_func3[Controlpoint_of_Element[El_No][I_No]] * dShapeFunc_from_paren(xez, El_No);
			dR = dShape_func3[Controlpoint_of_Element[El_No][I_No]] * dShapeFunc_from_paren(xez, El_No);
			//printf("dShape_func3[%d]:%le\n",Controlpoint_of_Element[El_No][I_No],dShape_func3[Controlpoint_of_Element[El_No][I_No]]);
		}
	}
	else
	{
		dR = ERROR;
	}
	//if (El_No == 21)
	//{
	//	printf("要素境界");
	//	printf("I_NO:%d\txez:%d\tdShapeFunc_from_paren:%lf\tdR:%le\n", I_No, xez, dShapeFunc_from_paren(xez, El_No), dR);
	//}
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


//逆行列を元の行列に代入
double InverseMatrix(double M[3][3])
{
	int i, j;
	double a[3][3];
	double det = M[0][0]*M[1][1]*M[2][2] + M[0][1]*M[1][2]*M[2][0] + M[0][2]*M[1][0]*M[2][1] - M[0][2]*M[1][1]*M[2][0] - M[0][0]*M[1][2]*M[2][1] - M[0][1]*M[1][0]*M[2][2];

	if (det == 0.0)
		return ERROR;

	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			a[i][j] = M[i][j];
		}
	}
	M[0][0] =        ( a[1][1]*a[2][2] - a[1][2]*a[2][1] ) / det;
	M[0][1] = (-1) * ( a[0][1]*a[2][2] - a[0][2]*a[2][1] ) / det;
	M[0][2] =        ( a[0][1]*a[1][2] - a[0][2]*a[1][1] ) / det;
	M[1][0] = (-1) * ( a[1][0]*a[2][2] - a[1][2]*a[2][0] ) / det;
	M[1][1] =        ( a[0][0]*a[2][2] - a[0][2]*a[2][0] ) / det;
	M[1][2] = (-1) * ( a[0][0]*a[1][2] - a[0][2]*a[1][0] ) / det;
	M[2][0] =        ( a[1][0]*a[2][1] - a[1][1]*a[2][0] ) / det;
	M[2][1] = (-1) * ( a[0][0]*a[2][1] - a[0][1]*a[2][0] ) / det;
	M[2][2] =        ( a[0][0]*a[1][1] - a[0][1]*a[1][0] ) / det;
	//printf("det;%le\n", det);
	return det;
}


/////////////////////////////////////////////////////////////////////////////////////////Newton-Raphson法////////////////////////////////////////////////////////////////////////
double rBasisFunc(double *knot_vec, int knot_index, int order, double xi, double *output, double *d_output)
{
	int p, j;
	double sum1 = 0.0;
	double sum2 = 0.0;
	double temp_basis[order][order];
	(*output) = 0.0;
	(*d_output) = 0.0;

	if ( knot_vec[knot_index] <= xi && xi <= knot_vec[knot_index + order + 1] )
	{
		if (knot_index == 0)
		{
			for (j = 0; j <= order; j++)
			{
				if ( (knot_vec[j] <= xi) && (xi <= knot_vec[j + 1]) )
				{
					temp_basis[j][0] = 1.0;
				} else
				{
					temp_basis[j][0] = 0.0;
				}
			}
		} else
		{
			for (j = 0; j <= order; j++)
			{
				if ( (knot_vec[knot_index + j] < xi) && (xi <= knot_vec[knot_index + j + 1]) )
				{
					temp_basis[j][0] = 1.0;
				} else
				{
					temp_basis[j][0] = 0.0;
				}
			}
		}

		if ( order > 0 )
		{
			for (p = 1; p <= order; p++)
			{
				for (j = 0; j <= order - p; j++)
				{
					sum1 = 0.0;
					sum2 = 0.0;
					if ( (knot_vec[knot_index + j + p] - knot_vec[knot_index + j]) != 0.0)
					{
						sum1 = ( xi - knot_vec[knot_index + j] ) / (knot_vec[knot_index + j + p] - knot_vec[knot_index + j]) * temp_basis[j][p - 1];
					}
					if ( (knot_vec[knot_index + j + p + 1] - knot_vec[knot_index + j + 1]) != 0.0)
					{
						sum2 = (knot_vec[knot_index + j + p + 1] - xi) / (knot_vec[knot_index + j + p + 1] - knot_vec[knot_index + j + 1]) * temp_basis[j + 1][p - 1];
					}
					temp_basis[j][p] = sum1 + sum2;
				}
			}
			sum1 = 0.0;
			sum2 = 0.0;
			if ( (knot_vec[knot_index + order] - knot_vec[knot_index]) != 0.0)
			{
				sum1 = order / (knot_vec[knot_index + order] - knot_vec[knot_index]) * temp_basis[0][order - 1];
			}
			if ( (knot_vec[knot_index + order + 1] - knot_vec[knot_index + 1]) != 0.0)
			{
				sum2 = order / (knot_vec[knot_index + order + 1] - knot_vec[knot_index + 1]) * temp_basis[1][order - 1];
			}
		}
		(*output) = temp_basis[0][order];
		(*d_output) = sum1 - sum2;
	}
	return (*output);
}

double lBasisFunc(double *knot_vec, int knot_index, int cntl_p_n, int order, double xi, double *output, double *d_output)
{
	int p, j;
	double sum1 = 0.0;
	double sum2 = 0.0;
	double temp_basis[order][order];
	(*output) = 0.0;
	(*d_output) = 0.0;


	if ( knot_vec[knot_index] <= xi && xi <= knot_vec[knot_index + order + 1] )
	{
		if (knot_index == cntl_p_n - 1)
		{
			for (j = 0; j <= order; j++)
			{
				if ( (knot_vec[cntl_p_n - 1 + j] <= xi) && (xi <= knot_vec[cntl_p_n + j]) )
				{
					temp_basis[j][0] = 1.0;
				} else
				{
					temp_basis[j][0] = 0.0;
				}
			}
		} else
		{
			for (j = 0; j <= order; j++)
			{
				if ( (knot_vec[knot_index + j] <= xi) && (xi < knot_vec[knot_index + j + 1]) )
				{
					temp_basis[j][0] = 1.0;
				} else
				{
					temp_basis[j][0] = 0.0;
				}
			}
		}

		if ( order > 0 )
		{
			for (p = 1; p <= order; p++)
			{
				for (j = 0; j <= order - p; j++)
				{
					sum1 = 0.0;
					sum2 = 0.0;
					if ( (knot_vec[knot_index + j + p] - knot_vec[knot_index + j]) != 0.0)
					{
						sum1 = ( xi - knot_vec[knot_index + j] ) / (knot_vec[knot_index + j + p] - knot_vec[knot_index + j]) * temp_basis[j][p - 1];
					}
					if ( (knot_vec[knot_index + j + p + 1] - knot_vec[knot_index + j + 1]) != 0.0)
					{
						sum2 = (knot_vec[knot_index + j + p + 1] - xi) / (knot_vec[knot_index + j + p + 1] - knot_vec[knot_index + j + 1]) * temp_basis[j + 1][p - 1];
					}
					temp_basis[j][p] = sum1 + sum2;
				}
			}
			sum1 = 0.0;
			sum2 = 0.0;
			if ( (knot_vec[knot_index + order] - knot_vec[knot_index]) != 0.0)
			{
				sum1 = order / (knot_vec[knot_index + order] - knot_vec[knot_index]) * temp_basis[0][order - 1];
			}
			if ( (knot_vec[knot_index + order + 1] - knot_vec[knot_index + 1]) != 0.0)
			{
				sum2 = order / (knot_vec[knot_index + order + 1] - knot_vec[knot_index + 1]) * temp_basis[1][order - 1];
			}
		}
		(*output) = temp_basis[0][order];
		(*d_output) = sum1 - sum2;
	}
	return (*output);
}

double rrrNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta, double *cntl_px, double *cntl_py, double *cntl_pz, int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta, double *weight, int order_xi, int order_eta, int order_zeta, double xi, double eta, double zeta, double *output_x, double *output_y, double *output_z, /*ここから*/double *output_dxi_x, double *output_deta_x, double *output_dzeta_x, double *output_dxi_y, double *output_deta_y, double *output_dzeta_y, double *output_dxi_z, double *output_deta_z, double *output_dzeta_z)
{
	int i, j, k, temp_index;
	double temp1, temp2, temp3, temp4;
	double molecule_x, molecule_y, molecule_z;
	double dxi_molecule_x, dxi_molecule_y, dxi_molecule_z;
	double deta_molecule_x, deta_molecule_y, deta_molecule_z;
	double dzeta_molecule_x, dzeta_molecule_y, dzeta_molecule_z;
	double denominator, dxi_denominator, deta_denominator, dzeta_denominator;
	double temp_output_xi, temp_output_eta, temp_output_zeta;
	double temp_d_output_xi, temp_d_output_eta, temp_d_output_zeta;

	molecule_x = 0.0;
	molecule_y = 0.0;
	molecule_z = 0.0;
	denominator = 0.0;

	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_molecule_z = 0.0;
	dxi_denominator = 0.0;

	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_molecule_z = 0.0;
	deta_denominator = 0.0;

	dzeta_molecule_x = 0.0;
	dzeta_molecule_y = 0.0;
	dzeta_molecule_z = 0.0;
	dzeta_denominator = 0.0;

	int index_min_xi   = 0;
	int index_max_xi   = cntl_p_n_xi - 1;
	int index_min_eta  = 0;
	int index_max_eta  = cntl_p_n_eta - 1;
	int index_min_zeta = 0;
	int index_max_zeta = cntl_p_n_zeta - 1;


	////////////////////////////////xi////////////////////////////////
	for (i = 0; i < cntl_p_n_xi; i++)
	{
		//printf("cntl_p_n_xi = %d\tknot_vec_xi[%d] = %lf\n",cntl_p_n_xi, i + 1, knot_vec_xi[i + 1]);
		if ( knot_vec_xi[i + 1] >= xi )
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			//printf("xi = %lf\tknot_vec_xi[%d] = %lf\n",xi ,i + 1, knot_vec_xi[i + 1]);
			break;
		}
	}

	if (index_min_xi < 0)
	{
		index_min_xi = 0;
		index_max_xi = order_xi + 1;
	}

	if (index_max_xi > cntl_p_n_xi)
	{
		index_min_xi = cntl_p_n_xi - order_xi - 1;
		index_max_xi = cntl_p_n_xi;
	}
	////////////////////////////////xi////////////////////////////////

	///////////////////////////////eta////////////////////////////////
	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if ( knot_vec_eta[i + 1] >= eta )
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			//printf("eta = %lf\tknot_vec_eta[%d] = %lf\n",eta ,i + 1, knot_vec_eta[i + 1]);
			break;
		}
	}
	if (index_min_eta < 0)
	{
		index_min_eta = 0;
		index_max_eta = order_eta + 1;
	}

	if (index_max_eta > cntl_p_n_eta)
	{
		index_min_eta = cntl_p_n_eta - order_eta - 1;
		index_max_eta = cntl_p_n_eta;
	}
	///////////////////////////////eta////////////////////////////////

	///////////////////////////////zeta///////////////////////////////
	for (i = 0; i < cntl_p_n_zeta; i++)
	{
		if ( knot_vec_zeta[i + 1] >= zeta )
		{
			index_min_zeta = i - order_zeta;
			index_max_zeta = i + 1;
			//printf("zeta = %lf\tknot_vec_zeta[%d] = %lf\n",zeta ,i + 1, knot_vec_zeta[i + 1]);
			break;
		}
	}
	if (index_min_zeta < 0)
	{
		index_min_zeta = 0;
		index_max_zeta = order_zeta + 1;
	}

	if (index_max_zeta > cntl_p_n_zeta)
	{
		index_min_zeta = cntl_p_n_zeta - order_zeta - 1;
		index_max_zeta = cntl_p_n_zeta;
	}
	///////////////////////////////zeta///////////////////////////////

	//printf("xi_min = %d\tmax = %d\t\teta_min = %d\tmax = %d\t\tzeta_min = %d\tmax = %d\tNo_eta = %d\n",index_min_xi, index_max_xi,  index_min_eta, index_max_eta, index_min_zeta, index_max_zeta, cntl_p_n_eta);

	for (i = index_min_xi; i < index_max_xi; i++)
	{
		rBasisFunc(knot_vec_xi, i, order_xi, xi, &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j < index_max_eta; j++)
		{
			rBasisFunc(knot_vec_eta, j, order_eta, eta, &temp_output_eta, &temp_d_output_eta);
			for (k = index_min_zeta; k < index_max_zeta; k++)
			{
				rBasisFunc(knot_vec_zeta, k, order_zeta, zeta, &temp_output_zeta, &temp_d_output_zeta);

				temp_index = i + (j * cntl_p_n_xi) + (k * cntl_p_n_xi * cntl_p_n_eta);
				//printf("i = %d\tj = %d\tk = %d\ttemp_index = %d\tPatch_controlpoint[0][%d] = %d\n", i, j, k, temp_index, temp_index , Patch_controlpoint[0][temp_index]);

				temp1 = temp_output_xi   * temp_output_eta   * temp_output_zeta   * weight[Patch_controlpoint[0][temp_index]];
				temp2 = temp_d_output_xi * temp_output_eta   * temp_output_zeta   * weight[Patch_controlpoint[0][temp_index]];
				temp3 = temp_output_xi   * temp_d_output_eta * temp_output_zeta   * weight[Patch_controlpoint[0][temp_index]];
				temp4 = temp_output_xi   * temp_output_eta   * temp_d_output_zeta * weight[Patch_controlpoint[0][temp_index]];

				molecule_x += temp1 * cntl_px[Patch_controlpoint[0][temp_index]];
				molecule_y += temp1 * cntl_py[Patch_controlpoint[0][temp_index]];
				molecule_z += temp1 * cntl_pz[Patch_controlpoint[0][temp_index]];
				denominator += temp1;

				//printf("temp1 = %lf\tpx = %lf\tpy = %lf\tpz = %lf\tconnec = %d\ttemp_index = %d\ti = %d\tj = %d\tk = %d\t\n",temp1, cntl_px[Patch_controlpoint[0][temp_index]], cntl_py[Patch_controlpoint[0][temp_index]], cntl_pz[Patch_controlpoint[0][temp_index]],Patch_controlpoint[0][temp_index],temp_index,i ,j, k);

				dxi_molecule_x += temp2 * cntl_px[Patch_controlpoint[0][temp_index]];
				dxi_molecule_y += temp2 * cntl_py[Patch_controlpoint[0][temp_index]];
				dxi_molecule_z += temp2 * cntl_pz[Patch_controlpoint[0][temp_index]];
				dxi_denominator += temp2;

				deta_molecule_x += temp3 * cntl_px[Patch_controlpoint[0][temp_index]];
				deta_molecule_y += temp3 * cntl_py[Patch_controlpoint[0][temp_index]];
				deta_molecule_z += temp3 * cntl_pz[Patch_controlpoint[0][temp_index]];
				deta_denominator += temp3;

				dzeta_molecule_x += temp4 * cntl_px[Patch_controlpoint[0][temp_index]];
				dzeta_molecule_y += temp4 * cntl_py[Patch_controlpoint[0][temp_index]];
				dzeta_molecule_z += temp4 * cntl_pz[Patch_controlpoint[0][temp_index]];
				dzeta_denominator += temp4;
			}
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;
	(*output_z) = molecule_z / denominator;

	//temp1 = denominator * denominator;
	(*output_dxi_x)   = weight[Patch_controlpoint[0][temp_index]] * ( dxi_molecule_x   * denominator - molecule_x * dxi_denominator)   / ( denominator * denominator ) ;
	(*output_dxi_y)   = weight[Patch_controlpoint[0][temp_index]] * ( dxi_molecule_y   * denominator - molecule_y * dxi_denominator)   / ( denominator * denominator ) ;
	(*output_dxi_z)   = weight[Patch_controlpoint[0][temp_index]] * ( dxi_molecule_z   * denominator - molecule_z * dxi_denominator)   / ( denominator * denominator ) ;
	(*output_deta_x)  = weight[Patch_controlpoint[0][temp_index]] * ( deta_molecule_x  * denominator - molecule_x * deta_denominator)  / ( denominator * denominator ) ;
	(*output_deta_y)  = weight[Patch_controlpoint[0][temp_index]] * ( deta_molecule_y  * denominator - molecule_y * deta_denominator)  / ( denominator * denominator ) ;
	(*output_deta_z)  = weight[Patch_controlpoint[0][temp_index]] * ( deta_molecule_z  * denominator - molecule_z * deta_denominator)  / ( denominator * denominator ) ;
	(*output_dzeta_x) = weight[Patch_controlpoint[0][temp_index]] * ( dzeta_molecule_x * denominator - molecule_x * dzeta_denominator) / ( denominator * denominator ) ;
	(*output_dzeta_y) = weight[Patch_controlpoint[0][temp_index]] * ( dzeta_molecule_y * denominator - molecule_y * dzeta_denominator) / ( denominator * denominator ) ;
	(*output_dzeta_z) = weight[Patch_controlpoint[0][temp_index]] * ( dzeta_molecule_z * denominator - molecule_z * dzeta_denominator) / ( denominator * denominator ) ;

	return denominator;
}

double lrrNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta, double *cntl_px, double *cntl_py, double *cntl_pz, int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta, double *weight, int order_xi, int order_eta, int order_zeta, double xi, double eta, double zeta, double *output_x, double *output_y, double *output_z, double *output_dxi_x, double *output_deta_x, double *output_dzeta_x, double *output_dxi_y, double *output_deta_y, double *output_dzeta_y, double *output_dxi_z, double *output_deta_z, double *output_dzeta_z)
{
	int i, j, k, temp_index;
	double temp1, temp2, temp3, temp4;
	double molecule_x, molecule_y, molecule_z;
	double dxi_molecule_x, dxi_molecule_y, dxi_molecule_z;
	double deta_molecule_x, deta_molecule_y, deta_molecule_z;
	double dzeta_molecule_x, dzeta_molecule_y, dzeta_molecule_z;
	double denominator, dxi_denominator, deta_denominator, dzeta_denominator;
	double temp_output_xi, temp_output_eta, temp_output_zeta;
	double temp_d_output_xi, temp_d_output_eta, temp_d_output_zeta;

	molecule_x = 0.0;
	molecule_y = 0.0;
	molecule_z = 0.0;
	denominator = 0.0;

	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_molecule_z = 0.0;
	dxi_denominator = 0.0;

	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_molecule_z = 0.0;
	deta_denominator = 0.0;

	dzeta_molecule_x = 0.0;
	dzeta_molecule_y = 0.0;
	dzeta_molecule_z = 0.0;
	dzeta_denominator = 0.0;

	int index_min_xi   = 0;
	int index_max_xi   = cntl_p_n_xi - 1;
	int index_min_eta  = 0;
	int index_max_eta  = cntl_p_n_eta - 1;
	int index_min_zeta = 0;
	int index_max_zeta = cntl_p_n_zeta - 1;

	////////////////////////////////xi////////////////////////////////
	for (i = 0; i < cntl_p_n_xi; i++)
	{
		//printf("cntl_p_n_xi = %d\tknot_vec_xi[%d] = %lf\n",cntl_p_n_xi, i + 1, knot_vec_xi[i + 1]);
		if ( knot_vec_xi[i + 1] >= xi )
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			//printf("xi = %lf\tknot_vec_xi[%d] = %lf\n",xi ,i + 1, knot_vec_xi[i + 1]);
			break;
		}
	}

	if (index_min_xi < 0)
	{
		index_min_xi = 0;
		index_max_xi = order_xi + 1;
	}

	if (index_max_xi > cntl_p_n_xi)
	{
		index_min_xi = cntl_p_n_xi - order_xi - 1;
		index_max_xi = cntl_p_n_xi;
	}
	////////////////////////////////xi////////////////////////////////

	///////////////////////////////eta////////////////////////////////
	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if ( knot_vec_eta[i + 1] >= eta )
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			//printf("eta = %lf\tknot_vec_eta[%d] = %lf\n",eta ,i + 1, knot_vec_eta[i + 1]);
			break;
		}
	}
	if (index_min_eta < 0)
	{
		index_min_eta = 0;
		index_max_eta = order_eta + 1;
	}

	if (index_max_eta > cntl_p_n_eta)
	{
		index_min_eta = cntl_p_n_eta - order_eta - 1;
		index_max_eta = cntl_p_n_eta;
	}
	///////////////////////////////eta////////////////////////////////

	///////////////////////////////zeta///////////////////////////////
	for (i = 0; i < cntl_p_n_zeta; i++)
	{
		if ( knot_vec_zeta[i + 1] >= zeta )
		{
			index_min_zeta = i - order_zeta;
			index_max_zeta = i + 1;
			//printf("zeta = %lf\tknot_vec_zeta[%d] = %lf\n",zeta ,i + 1, knot_vec_zeta[i + 1]);
			break;
		}
	}
	if (index_min_zeta < 0)
	{
		index_min_zeta = 0;
		index_max_zeta = order_zeta + 1;
	}

	if (index_max_zeta > cntl_p_n_zeta)
	{
		index_min_zeta = cntl_p_n_zeta - order_zeta - 1;
		index_max_zeta = cntl_p_n_zeta;
	}
	///////////////////////////////zeta///////////////////////////////




	for (i = index_min_xi; i < index_max_xi; i++)
	{
		lBasisFunc(knot_vec_xi, i, cntl_p_n_xi, order_xi, xi, &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j < index_max_eta; j++)
		{
			rBasisFunc(knot_vec_eta, j, order_eta, eta, &temp_output_eta, &temp_d_output_eta);
			for (k = index_min_zeta; k < index_max_zeta; k++)
			{
				rBasisFunc(knot_vec_zeta, k, order_zeta, zeta, &temp_output_zeta, &temp_d_output_zeta);

				temp_index = i + j * cntl_p_n_xi + k * cntl_p_n_xi * cntl_p_n_eta;

				temp1 = temp_output_xi   * temp_output_eta   * temp_output_zeta   * weight[Patch_controlpoint[0][temp_index]];
				temp2 = temp_d_output_xi * temp_output_eta   * temp_output_zeta   * weight[Patch_controlpoint[0][temp_index]];
				temp3 = temp_output_xi   * temp_d_output_eta * temp_output_zeta   * weight[Patch_controlpoint[0][temp_index]];
				temp4 = temp_output_xi   * temp_output_eta   * temp_d_output_zeta * weight[Patch_controlpoint[0][temp_index]];

				molecule_x += temp1 * cntl_px[Patch_controlpoint[0][temp_index]];
				molecule_y += temp1 * cntl_py[Patch_controlpoint[0][temp_index]];
				molecule_z += temp1 * cntl_pz[Patch_controlpoint[0][temp_index]];
				denominator += temp1;

				dxi_molecule_x += temp2 * cntl_px[Patch_controlpoint[0][temp_index]];
				dxi_molecule_y += temp2 * cntl_py[Patch_controlpoint[0][temp_index]];
				dxi_molecule_z += temp2 * cntl_pz[Patch_controlpoint[0][temp_index]];
				dxi_denominator += temp2;

				deta_molecule_x += temp3 * cntl_px[Patch_controlpoint[0][temp_index]];
				deta_molecule_y += temp3 * cntl_py[Patch_controlpoint[0][temp_index]];
				deta_molecule_z += temp3 * cntl_pz[Patch_controlpoint[0][temp_index]];
				deta_denominator += temp3;

				dzeta_molecule_x += temp4 * cntl_px[Patch_controlpoint[0][temp_index]];
				dzeta_molecule_y += temp4 * cntl_py[Patch_controlpoint[0][temp_index]];
				dzeta_molecule_z += temp4 * cntl_pz[Patch_controlpoint[0][temp_index]];
				dzeta_denominator += temp4;
			}
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;
	(*output_z) = molecule_z / denominator;

	//temp1 = denominator * denominator;
	(*output_dxi_x)   = weight[Patch_controlpoint[0][temp_index]] * ( dxi_molecule_x   * denominator - molecule_x * dxi_denominator)   / ( denominator * denominator ) ;
	(*output_dxi_y)   = weight[Patch_controlpoint[0][temp_index]] * ( dxi_molecule_y   * denominator - molecule_y * dxi_denominator)   / ( denominator * denominator ) ;
	(*output_dxi_z)   = weight[Patch_controlpoint[0][temp_index]] * ( dxi_molecule_z   * denominator - molecule_z * dxi_denominator)   / ( denominator * denominator ) ;
	(*output_deta_x)  = weight[Patch_controlpoint[0][temp_index]] * ( deta_molecule_x  * denominator - molecule_x * deta_denominator)  / ( denominator * denominator ) ;
	(*output_deta_y)  = weight[Patch_controlpoint[0][temp_index]] * ( deta_molecule_y  * denominator - molecule_y * deta_denominator)  / ( denominator * denominator ) ;
	(*output_deta_z)  = weight[Patch_controlpoint[0][temp_index]] * ( deta_molecule_z  * denominator - molecule_z * deta_denominator)  / ( denominator * denominator ) ;
	(*output_dzeta_x) = weight[Patch_controlpoint[0][temp_index]] * ( dzeta_molecule_x * denominator - molecule_x * dzeta_denominator) / ( denominator * denominator ) ;
	(*output_dzeta_y) = weight[Patch_controlpoint[0][temp_index]] * ( dzeta_molecule_y * denominator - molecule_y * dzeta_denominator) / ( denominator * denominator ) ;
	(*output_dzeta_z) = weight[Patch_controlpoint[0][temp_index]] * ( dzeta_molecule_z * denominator - molecule_z * dzeta_denominator) / ( denominator * denominator ) ;

	return denominator;
}

double rlrNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta, double *cntl_px, double *cntl_py, double *cntl_pz, int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta, double *weight, int order_xi, int order_eta, int order_zeta, double xi, double eta, double zeta, double *output_x, double *output_y, double *output_z, double *output_dxi_x, double *output_deta_x, double *output_dzeta_x, double *output_dxi_y, double *output_deta_y, double *output_dzeta_y, double *output_dxi_z, double *output_deta_z, double *output_dzeta_z)
{
	int i, j, k, temp_index;
	double temp1, temp2, temp3, temp4;
	double molecule_x, molecule_y, molecule_z;
	double dxi_molecule_x, dxi_molecule_y, dxi_molecule_z;
	double deta_molecule_x, deta_molecule_y, deta_molecule_z;
	double dzeta_molecule_x, dzeta_molecule_y, dzeta_molecule_z;
	double denominator, dxi_denominator, deta_denominator, dzeta_denominator;
	double temp_output_xi, temp_output_eta, temp_output_zeta;
	double temp_d_output_xi, temp_d_output_eta, temp_d_output_zeta;

	molecule_x = 0.0;
	molecule_y = 0.0;
	molecule_z = 0.0;
	denominator = 0.0;

	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_molecule_z = 0.0;
	dxi_denominator = 0.0;

	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_molecule_z = 0.0;
	deta_denominator = 0.0;

	dzeta_molecule_x = 0.0;
	dzeta_molecule_y = 0.0;
	dzeta_molecule_z = 0.0;
	dzeta_denominator = 0.0;

	int index_min_xi   = 0;
	int index_max_xi   = cntl_p_n_xi - 1;
	int index_min_eta  = 0;
	int index_max_eta  = cntl_p_n_eta - 1;
	int index_min_zeta = 0;
	int index_max_zeta = cntl_p_n_zeta - 1;

	////////////////////////////////xi////////////////////////////////
	for (i = 0; i < cntl_p_n_xi; i++)
	{
		//printf("cntl_p_n_xi = %d\tknot_vec_xi[%d] = %lf\n",cntl_p_n_xi, i + 1, knot_vec_xi[i + 1]);
		if ( knot_vec_xi[i + 1] >= xi )
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			//printf("xi = %lf\tknot_vec_xi[%d] = %lf\n",xi ,i + 1, knot_vec_xi[i + 1]);
			break;
		}
	}

	if (index_min_xi < 0)
	{
		index_min_xi = 0;
		index_max_xi = order_xi + 1;
	}

	if (index_max_xi > cntl_p_n_xi)
	{
		index_min_xi = cntl_p_n_xi - order_xi - 1;
		index_max_xi = cntl_p_n_xi;
	}
	////////////////////////////////xi////////////////////////////////

	///////////////////////////////eta////////////////////////////////
	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if ( knot_vec_eta[i + 1] >= eta )
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			//printf("eta = %lf\tknot_vec_eta[%d] = %lf\n",eta ,i + 1, knot_vec_eta[i + 1]);
			break;
		}
	}
	if (index_min_eta < 0)
	{
		index_min_eta = 0;
		index_max_eta = order_eta + 1;
	}

	if (index_max_eta > cntl_p_n_eta)
	{
		index_min_eta = cntl_p_n_eta - order_eta - 1;
		index_max_eta = cntl_p_n_eta;
	}
	///////////////////////////////eta////////////////////////////////

	///////////////////////////////zeta///////////////////////////////
	for (i = 0; i < cntl_p_n_zeta; i++)
	{
		if ( knot_vec_zeta[i + 1] >= zeta )
		{
			index_min_zeta = i - order_zeta;
			index_max_zeta = i + 1;
			//printf("zeta = %lf\tknot_vec_zeta[%d] = %lf\n",zeta ,i + 1, knot_vec_zeta[i + 1]);
			break;
		}
	}
	if (index_min_zeta < 0)
	{
		index_min_zeta = 0;
		index_max_zeta = order_zeta + 1;
	}

	if (index_max_zeta > cntl_p_n_zeta)
	{
		index_min_zeta = cntl_p_n_zeta - order_zeta - 1;
		index_max_zeta = cntl_p_n_zeta;
	}
	///////////////////////////////zeta///////////////////////////////



	for (i = index_min_xi; i < index_max_xi; i++)
	{
		rBasisFunc(knot_vec_xi, i, order_xi, xi, &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j < index_max_eta; j++)
		{
			lBasisFunc(knot_vec_eta, j, cntl_p_n_eta, order_eta, eta, &temp_output_eta, &temp_d_output_eta);
			for (k = index_min_zeta; k < index_max_zeta; k++)
			{
				rBasisFunc(knot_vec_zeta, k, order_zeta, zeta, &temp_output_zeta, &temp_d_output_zeta);

				temp_index = i + j * cntl_p_n_xi + k * cntl_p_n_xi * cntl_p_n_eta;

				temp1 = temp_output_xi   * temp_output_eta   * temp_output_zeta   * weight[Patch_controlpoint[0][temp_index]];
				temp2 = temp_d_output_xi * temp_output_eta   * temp_output_zeta   * weight[Patch_controlpoint[0][temp_index]];
				temp3 = temp_output_xi   * temp_d_output_eta * temp_output_zeta   * weight[Patch_controlpoint[0][temp_index]];
				temp4 = temp_output_xi   * temp_output_eta   * temp_d_output_zeta * weight[Patch_controlpoint[0][temp_index]];

				molecule_x += temp1 * cntl_px[Patch_controlpoint[0][temp_index]];
				molecule_y += temp1 * cntl_py[Patch_controlpoint[0][temp_index]];
				molecule_z += temp1 * cntl_pz[Patch_controlpoint[0][temp_index]];
				denominator += temp1;

				dxi_molecule_x += temp2 * cntl_px[Patch_controlpoint[0][temp_index]];
				dxi_molecule_y += temp2 * cntl_py[Patch_controlpoint[0][temp_index]];
				dxi_molecule_z += temp2 * cntl_pz[Patch_controlpoint[0][temp_index]];
				dxi_denominator += temp2;

				deta_molecule_x += temp3 * cntl_px[Patch_controlpoint[0][temp_index]];
				deta_molecule_y += temp3 * cntl_py[Patch_controlpoint[0][temp_index]];
				deta_molecule_z += temp3 * cntl_pz[Patch_controlpoint[0][temp_index]];
				deta_denominator += temp3;

				dzeta_molecule_x += temp4 * cntl_px[Patch_controlpoint[0][temp_index]];
				dzeta_molecule_y += temp4 * cntl_py[Patch_controlpoint[0][temp_index]];
				dzeta_molecule_z += temp4 * cntl_pz[Patch_controlpoint[0][temp_index]];
				dzeta_denominator += temp4;
			}
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;
	(*output_z) = molecule_z / denominator;

	//temp1 = denominator * denominator;
	(*output_dxi_x)   = weight[Patch_controlpoint[0][temp_index]] * ( dxi_molecule_x   * denominator - molecule_x * dxi_denominator)   / ( denominator * denominator ) ;
	(*output_dxi_y)   = weight[Patch_controlpoint[0][temp_index]] * ( dxi_molecule_y   * denominator - molecule_y * dxi_denominator)   / ( denominator * denominator ) ;
	(*output_dxi_z)   = weight[Patch_controlpoint[0][temp_index]] * ( dxi_molecule_z   * denominator - molecule_z * dxi_denominator)   / ( denominator * denominator ) ;
	(*output_deta_x)  = weight[Patch_controlpoint[0][temp_index]] * ( deta_molecule_x  * denominator - molecule_x * deta_denominator)  / ( denominator * denominator ) ;
	(*output_deta_y)  = weight[Patch_controlpoint[0][temp_index]] * ( deta_molecule_y  * denominator - molecule_y * deta_denominator)  / ( denominator * denominator ) ;
	(*output_deta_z)  = weight[Patch_controlpoint[0][temp_index]] * ( deta_molecule_z  * denominator - molecule_z * deta_denominator)  / ( denominator * denominator ) ;
	(*output_dzeta_x) = weight[Patch_controlpoint[0][temp_index]] * ( dzeta_molecule_x * denominator - molecule_x * dzeta_denominator) / ( denominator * denominator ) ;
	(*output_dzeta_y) = weight[Patch_controlpoint[0][temp_index]] * ( dzeta_molecule_y * denominator - molecule_y * dzeta_denominator) / ( denominator * denominator ) ;
	(*output_dzeta_z) = weight[Patch_controlpoint[0][temp_index]] * ( dzeta_molecule_z * denominator - molecule_z * dzeta_denominator) / ( denominator * denominator ) ;
	return denominator;
}

double rrlNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta, double *cntl_px, double *cntl_py, double *cntl_pz, int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta, double *weight, int order_xi, int order_eta, int order_zeta, double xi, double eta, double zeta, double *output_x, double *output_y, double *output_z, double *output_dxi_x, double *output_deta_x, double *output_dzeta_x, double *output_dxi_y, double *output_deta_y, double *output_dzeta_y, double *output_dxi_z, double *output_deta_z, double *output_dzeta_z)
{
	int i, j, k, temp_index;
	double temp1, temp2, temp3, temp4;
	double molecule_x, molecule_y, molecule_z;
	double dxi_molecule_x, dxi_molecule_y, dxi_molecule_z;
	double deta_molecule_x, deta_molecule_y, deta_molecule_z;
	double dzeta_molecule_x, dzeta_molecule_y, dzeta_molecule_z;
	double denominator, dxi_denominator, deta_denominator, dzeta_denominator;
	double temp_output_xi, temp_output_eta, temp_output_zeta;
	double temp_d_output_xi, temp_d_output_eta, temp_d_output_zeta;

	molecule_x = 0.0;
	molecule_y = 0.0;
	molecule_z = 0.0;
	denominator = 0.0;

	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_molecule_z = 0.0;
	dxi_denominator = 0.0;

	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_molecule_z = 0.0;
	deta_denominator = 0.0;

	dzeta_molecule_x = 0.0;
	dzeta_molecule_y = 0.0;
	dzeta_molecule_z = 0.0;
	dzeta_denominator = 0.0;

	int index_min_xi   = 0;
	int index_max_xi   = cntl_p_n_xi - 1;
	int index_min_eta  = 0;
	int index_max_eta  = cntl_p_n_eta - 1;
	int index_min_zeta = 0;
	int index_max_zeta = cntl_p_n_zeta - 1;

	////////////////////////////////xi////////////////////////////////
	for (i = 0; i < cntl_p_n_xi; i++)
	{
		//printf("cntl_p_n_xi = %d\tknot_vec_xi[%d] = %lf\n",cntl_p_n_xi, i + 1, knot_vec_xi[i + 1]);
		if ( knot_vec_xi[i + 1] >= xi )
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			//printf("xi = %lf\tknot_vec_xi[%d] = %lf\n",xi ,i + 1, knot_vec_xi[i + 1]);
			break;
		}
	}

	if (index_min_xi < 0)
	{
		index_min_xi = 0;
		index_max_xi = order_xi + 1;
	}

	if (index_max_xi > cntl_p_n_xi)
	{
		index_min_xi = cntl_p_n_xi - order_xi - 1;
		index_max_xi = cntl_p_n_xi;
	}
	////////////////////////////////xi////////////////////////////////

	///////////////////////////////eta////////////////////////////////
	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if ( knot_vec_eta[i + 1] >= eta )
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			//printf("eta = %lf\tknot_vec_eta[%d] = %lf\n",eta ,i + 1, knot_vec_eta[i + 1]);
			break;
		}
	}
	if (index_min_eta < 0)
	{
		index_min_eta = 0;
		index_max_eta = order_eta + 1;
	}

	if (index_max_eta > cntl_p_n_eta)
	{
		index_min_eta = cntl_p_n_eta - order_eta - 1;
		index_max_eta = cntl_p_n_eta;
	}
	///////////////////////////////eta////////////////////////////////

	///////////////////////////////zeta///////////////////////////////
	for (i = 0; i < cntl_p_n_zeta; i++)
	{
		if ( knot_vec_zeta[i + 1] >= zeta )
		{
			index_min_zeta = i - order_zeta;
			index_max_zeta = i + 1;
			//printf("zeta = %lf\tknot_vec_zeta[%d] = %lf\n",zeta ,i + 1, knot_vec_zeta[i + 1]);
			break;
		}
	}
	if (index_min_zeta < 0)
	{
		index_min_zeta = 0;
		index_max_zeta = order_zeta + 1;
	}

	if (index_max_zeta > cntl_p_n_zeta)
	{
		index_min_zeta = cntl_p_n_zeta - order_zeta - 1;
		index_max_zeta = cntl_p_n_zeta;
	}
	///////////////////////////////zeta///////////////////////////////



	for (i = index_min_xi; i < index_max_xi; i++)
	{
		rBasisFunc(knot_vec_xi, i, order_xi, xi, &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j < index_max_eta; j++)
		{
			rBasisFunc(knot_vec_eta, j, order_eta, eta, &temp_output_eta, &temp_d_output_eta);
			for (k = index_min_zeta; k < index_max_zeta; k++)
			{
				lBasisFunc(knot_vec_zeta, k, cntl_p_n_zeta, order_zeta, zeta, &temp_output_zeta, &temp_d_output_zeta);

				temp_index = i + j * cntl_p_n_xi + k * cntl_p_n_xi * cntl_p_n_eta;

				temp1 = temp_output_xi   * temp_output_eta   * temp_output_zeta   * weight[Patch_controlpoint[0][temp_index]];
				temp2 = temp_d_output_xi * temp_output_eta   * temp_output_zeta   * weight[Patch_controlpoint[0][temp_index]];
				temp3 = temp_output_xi   * temp_d_output_eta * temp_output_zeta   * weight[Patch_controlpoint[0][temp_index]];
				temp4 = temp_output_xi   * temp_output_eta   * temp_d_output_zeta * weight[Patch_controlpoint[0][temp_index]];

				molecule_x += temp1 * cntl_px[Patch_controlpoint[0][temp_index]];
				molecule_y += temp1 * cntl_py[Patch_controlpoint[0][temp_index]];
				molecule_z += temp1 * cntl_pz[Patch_controlpoint[0][temp_index]];
				denominator += temp1;

				dxi_molecule_x += temp2 * cntl_px[Patch_controlpoint[0][temp_index]];
				dxi_molecule_y += temp2 * cntl_py[Patch_controlpoint[0][temp_index]];
				dxi_molecule_z += temp2 * cntl_pz[Patch_controlpoint[0][temp_index]];
				dxi_denominator += temp2;

				deta_molecule_x += temp3 * cntl_px[Patch_controlpoint[0][temp_index]];
				deta_molecule_y += temp3 * cntl_py[Patch_controlpoint[0][temp_index]];
				deta_molecule_z += temp3 * cntl_pz[Patch_controlpoint[0][temp_index]];
				deta_denominator += temp3;

				dzeta_molecule_x += temp4 * cntl_px[Patch_controlpoint[0][temp_index]];
				dzeta_molecule_y += temp4 * cntl_py[Patch_controlpoint[0][temp_index]];
				dzeta_molecule_z += temp4 * cntl_pz[Patch_controlpoint[0][temp_index]];
				dzeta_denominator += temp4;
			}
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;
	(*output_z) = molecule_z / denominator;

	//temp1 = denominator * denominator;
	(*output_dxi_x)   = weight[Patch_controlpoint[0][temp_index]] * ( dxi_molecule_x   * denominator - molecule_x * dxi_denominator)   / ( denominator * denominator ) ;
	(*output_dxi_y)   = weight[Patch_controlpoint[0][temp_index]] * ( dxi_molecule_y   * denominator - molecule_y * dxi_denominator)   / ( denominator * denominator ) ;
	(*output_dxi_z)   = weight[Patch_controlpoint[0][temp_index]] * ( dxi_molecule_z   * denominator - molecule_z * dxi_denominator)   / ( denominator * denominator ) ;
	(*output_deta_x)  = weight[Patch_controlpoint[0][temp_index]] * ( deta_molecule_x  * denominator - molecule_x * deta_denominator)  / ( denominator * denominator ) ;
	(*output_deta_y)  = weight[Patch_controlpoint[0][temp_index]] * ( deta_molecule_y  * denominator - molecule_y * deta_denominator)  / ( denominator * denominator ) ;
	(*output_deta_z)  = weight[Patch_controlpoint[0][temp_index]] * ( deta_molecule_z  * denominator - molecule_z * deta_denominator)  / ( denominator * denominator ) ;
	(*output_dzeta_x) = weight[Patch_controlpoint[0][temp_index]] * ( dzeta_molecule_x * denominator - molecule_x * dzeta_denominator) / ( denominator * denominator ) ;
	(*output_dzeta_y) = weight[Patch_controlpoint[0][temp_index]] * ( dzeta_molecule_y * denominator - molecule_y * dzeta_denominator) / ( denominator * denominator ) ;
	(*output_dzeta_z) = weight[Patch_controlpoint[0][temp_index]] * ( dzeta_molecule_z * denominator - molecule_z * dzeta_denominator) / ( denominator * denominator ) ;

	return denominator;
}

double llrNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta, double *cntl_px, double *cntl_py, double *cntl_pz, int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta, double *weight, int order_xi, int order_eta, int order_zeta, double xi, double eta, double zeta, double *output_x, double *output_y, double *output_z, double *output_dxi_x, double *output_deta_x, double *output_dzeta_x, double *output_dxi_y, double *output_deta_y, double *output_dzeta_y, double *output_dxi_z, double *output_deta_z, double *output_dzeta_z)
{
	int i, j, k, temp_index;
	double temp1, temp2, temp3, temp4;
	double molecule_x, molecule_y, molecule_z;
	double dxi_molecule_x, dxi_molecule_y, dxi_molecule_z;
	double deta_molecule_x, deta_molecule_y, deta_molecule_z;
	double dzeta_molecule_x, dzeta_molecule_y, dzeta_molecule_z;
	double denominator, dxi_denominator, deta_denominator, dzeta_denominator;
	double temp_output_xi, temp_output_eta, temp_output_zeta;
	double temp_d_output_xi, temp_d_output_eta, temp_d_output_zeta;

	molecule_x = 0.0;
	molecule_y = 0.0;
	molecule_z = 0.0;
	denominator = 0.0;

	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_molecule_z = 0.0;
	dxi_denominator = 0.0;

	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_molecule_z = 0.0;
	deta_denominator = 0.0;

	dzeta_molecule_x = 0.0;
	dzeta_molecule_y = 0.0;
	dzeta_molecule_z = 0.0;
	dzeta_denominator = 0.0;

	int index_min_xi   = 0;
	int index_max_xi   = cntl_p_n_xi - 1;
	int index_min_eta  = 0;
	int index_max_eta  = cntl_p_n_eta - 1;
	int index_min_zeta = 0;
	int index_max_zeta = cntl_p_n_zeta - 1;

	////////////////////////////////xi////////////////////////////////
	for (i = 0; i < cntl_p_n_xi; i++)
	{
		//printf("cntl_p_n_xi = %d\tknot_vec_xi[%d] = %lf\n",cntl_p_n_xi, i + 1, knot_vec_xi[i + 1]);
		if ( knot_vec_xi[i + 1] >= xi )
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			//printf("xi = %lf\tknot_vec_xi[%d] = %lf\n",xi ,i + 1, knot_vec_xi[i + 1]);
			break;
		}
	}

	if (index_min_xi < 0)
	{
		index_min_xi = 0;
		index_max_xi = order_xi + 1;
	}

	if (index_max_xi > cntl_p_n_xi)
	{
		index_min_xi = cntl_p_n_xi - order_xi - 1;
		index_max_xi = cntl_p_n_xi;
	}
	////////////////////////////////xi////////////////////////////////

	///////////////////////////////eta////////////////////////////////
	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if ( knot_vec_eta[i + 1] >= eta )
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			//printf("eta = %lf\tknot_vec_eta[%d] = %lf\n",eta ,i + 1, knot_vec_eta[i + 1]);
			break;
		}
	}
	if (index_min_eta < 0)
	{
		index_min_eta = 0;
		index_max_eta = order_eta + 1;
	}

	if (index_max_eta > cntl_p_n_eta)
	{
		index_min_eta = cntl_p_n_eta - order_eta - 1;
		index_max_eta = cntl_p_n_eta;
	}
	///////////////////////////////eta////////////////////////////////

	///////////////////////////////zeta///////////////////////////////
	for (i = 0; i < cntl_p_n_zeta; i++)
	{
		if ( knot_vec_zeta[i + 1] >= zeta )
		{
			index_min_zeta = i - order_zeta;
			index_max_zeta = i + 1;
			//printf("zeta = %lf\tknot_vec_zeta[%d] = %lf\n",zeta ,i + 1, knot_vec_zeta[i + 1]);
			break;
		}
	}
	if (index_min_zeta < 0)
	{
		index_min_zeta = 0;
		index_max_zeta = order_zeta + 1;
	}

	if (index_max_zeta > cntl_p_n_zeta)
	{
		index_min_zeta = cntl_p_n_zeta - order_zeta - 1;
		index_max_zeta = cntl_p_n_zeta;
	}
	///////////////////////////////zeta///////////////////////////////


	for (i = index_min_xi; i < index_max_xi; i++)
	{
		lBasisFunc(knot_vec_xi, i, cntl_p_n_xi, order_xi, xi, &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j < index_max_eta; j++)
		{
			lBasisFunc(knot_vec_eta, j, cntl_p_n_eta, order_eta, eta, &temp_output_eta, &temp_d_output_eta);
			for (k = index_min_zeta; k < index_max_zeta; k++)
			{
				rBasisFunc(knot_vec_zeta, k, order_zeta, zeta, &temp_output_zeta, &temp_d_output_zeta);

				temp_index = i + j * cntl_p_n_xi + k * cntl_p_n_xi * cntl_p_n_eta;

				temp1 = temp_output_xi   * temp_output_eta   * temp_output_zeta   * weight[Patch_controlpoint[0][temp_index]];
				temp2 = temp_d_output_xi * temp_output_eta   * temp_output_zeta   * weight[Patch_controlpoint[0][temp_index]];
				temp3 = temp_output_xi   * temp_d_output_eta * temp_output_zeta   * weight[Patch_controlpoint[0][temp_index]];
				temp4 = temp_output_xi   * temp_output_eta   * temp_d_output_zeta * weight[Patch_controlpoint[0][temp_index]];

				molecule_x += temp1 * cntl_px[Patch_controlpoint[0][temp_index]];
				molecule_y += temp1 * cntl_py[Patch_controlpoint[0][temp_index]];
				molecule_z += temp1 * cntl_pz[Patch_controlpoint[0][temp_index]];
				denominator += temp1;

				dxi_molecule_x += temp2 * cntl_px[Patch_controlpoint[0][temp_index]];
				dxi_molecule_y += temp2 * cntl_py[Patch_controlpoint[0][temp_index]];
				dxi_molecule_z += temp2 * cntl_pz[Patch_controlpoint[0][temp_index]];
				dxi_denominator += temp2;

				deta_molecule_x += temp3 * cntl_px[Patch_controlpoint[0][temp_index]];
				deta_molecule_y += temp3 * cntl_py[Patch_controlpoint[0][temp_index]];
				deta_molecule_z += temp3 * cntl_pz[Patch_controlpoint[0][temp_index]];
				deta_denominator += temp3;

				dzeta_molecule_x += temp4 * cntl_px[temp_index];
				dzeta_molecule_y += temp4 * cntl_py[temp_index];
				dzeta_molecule_z += temp4 * cntl_pz[temp_index];
				dzeta_denominator += temp4;
			}
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;
	(*output_z) = molecule_z / denominator;

	//temp1 = denominator * denominator;
	(*output_dxi_x)   = weight[Patch_controlpoint[0][temp_index]] * ( dxi_molecule_x   * denominator - molecule_x * dxi_denominator)   / ( denominator * denominator ) ;
	(*output_dxi_y)   = weight[Patch_controlpoint[0][temp_index]] * ( dxi_molecule_y   * denominator - molecule_y * dxi_denominator)   / ( denominator * denominator ) ;
	(*output_dxi_z)   = weight[Patch_controlpoint[0][temp_index]] * ( dxi_molecule_z   * denominator - molecule_z * dxi_denominator)   / ( denominator * denominator ) ;
	(*output_deta_x)  = weight[Patch_controlpoint[0][temp_index]] * ( deta_molecule_x  * denominator - molecule_x * deta_denominator)  / ( denominator * denominator ) ;
	(*output_deta_y)  = weight[Patch_controlpoint[0][temp_index]] * ( deta_molecule_y  * denominator - molecule_y * deta_denominator)  / ( denominator * denominator ) ;
	(*output_deta_z)  = weight[Patch_controlpoint[0][temp_index]] * ( deta_molecule_z  * denominator - molecule_z * deta_denominator)  / ( denominator * denominator ) ;
	(*output_dzeta_x) = weight[Patch_controlpoint[0][temp_index]] * ( dzeta_molecule_x * denominator - molecule_x * dzeta_denominator) / ( denominator * denominator ) ;
	(*output_dzeta_y) = weight[Patch_controlpoint[0][temp_index]] * ( dzeta_molecule_y * denominator - molecule_y * dzeta_denominator) / ( denominator * denominator ) ;
	(*output_dzeta_z) = weight[Patch_controlpoint[0][temp_index]] * ( dzeta_molecule_z * denominator - molecule_z * dzeta_denominator) / ( denominator * denominator ) ;

	return denominator;
}

double lrlNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta, double *cntl_px, double *cntl_py, double *cntl_pz, int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta, double *weight, int order_xi, int order_eta, int order_zeta, double xi, double eta, double zeta, double *output_x, double *output_y, double *output_z, double *output_dxi_x, double *output_deta_x, double *output_dzeta_x, double *output_dxi_y, double *output_deta_y, double *output_dzeta_y, double *output_dxi_z, double *output_deta_z, double *output_dzeta_z)
{
	int i, j, k, temp_index;
	double temp1, temp2, temp3, temp4;
	double molecule_x, molecule_y, molecule_z;
	double dxi_molecule_x, dxi_molecule_y, dxi_molecule_z;
	double deta_molecule_x, deta_molecule_y, deta_molecule_z;
	double dzeta_molecule_x, dzeta_molecule_y, dzeta_molecule_z;
	double denominator, dxi_denominator, deta_denominator, dzeta_denominator;
	double temp_output_xi, temp_output_eta, temp_output_zeta;
	double temp_d_output_xi, temp_d_output_eta, temp_d_output_zeta;

	molecule_x = 0.0;
	molecule_y = 0.0;
	molecule_z = 0.0;
	denominator = 0.0;

	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_molecule_z = 0.0;
	dxi_denominator = 0.0;

	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_molecule_z = 0.0;
	deta_denominator = 0.0;

	dzeta_molecule_x = 0.0;
	dzeta_molecule_y = 0.0;
	dzeta_molecule_z = 0.0;
	dzeta_denominator = 0.0;

	int index_min_xi   = 0;
	int index_max_xi   = cntl_p_n_xi - 1;
	int index_min_eta  = 0;
	int index_max_eta  = cntl_p_n_eta - 1;
	int index_min_zeta = 0;
	int index_max_zeta = cntl_p_n_zeta - 1;

	////////////////////////////////xi////////////////////////////////
	for (i = 0; i < cntl_p_n_xi; i++)
	{
		//printf("cntl_p_n_xi = %d\tknot_vec_xi[%d] = %lf\n",cntl_p_n_xi, i + 1, knot_vec_xi[i + 1]);
		if ( knot_vec_xi[i + 1] >= xi )
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			//printf("xi = %lf\tknot_vec_xi[%d] = %lf\n",xi ,i + 1, knot_vec_xi[i + 1]);
			break;
		}
	}

	if (index_min_xi < 0)
	{
		index_min_xi = 0;
		index_max_xi = order_xi + 1;
	}

	if (index_max_xi > cntl_p_n_xi)
	{
		index_min_xi = cntl_p_n_xi - order_xi - 1;
		index_max_xi = cntl_p_n_xi;
	}
	////////////////////////////////xi////////////////////////////////

	///////////////////////////////eta////////////////////////////////
	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if ( knot_vec_eta[i + 1] >= eta )
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			//printf("eta = %lf\tknot_vec_eta[%d] = %lf\n",eta ,i + 1, knot_vec_eta[i + 1]);
			break;
		}
	}
	if (index_min_eta < 0)
	{
		index_min_eta = 0;
		index_max_eta = order_eta + 1;
	}

	if (index_max_eta > cntl_p_n_eta)
	{
		index_min_eta = cntl_p_n_eta - order_eta - 1;
		index_max_eta = cntl_p_n_eta;
	}
	///////////////////////////////eta////////////////////////////////

	///////////////////////////////zeta///////////////////////////////
	for (i = 0; i < cntl_p_n_zeta; i++)
	{
		if ( knot_vec_zeta[i + 1] >= zeta )
		{
			index_min_zeta = i - order_zeta;
			index_max_zeta = i + 1;
			//printf("zeta = %lf\tknot_vec_zeta[%d] = %lf\n",zeta ,i + 1, knot_vec_zeta[i + 1]);
			break;
		}
	}
	if (index_min_zeta < 0)
	{
		index_min_zeta = 0;
		index_max_zeta = order_zeta + 1;
	}

	if (index_max_zeta > cntl_p_n_zeta)
	{
		index_min_zeta = cntl_p_n_zeta - order_zeta - 1;
		index_max_zeta = cntl_p_n_zeta;
	}
	///////////////////////////////zeta///////////////////////////////



	for (i = index_min_xi; i < index_max_xi; i++)
	{
		lBasisFunc(knot_vec_xi, i, cntl_p_n_xi, order_xi, xi, &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j < index_max_eta; j++)
		{
			rBasisFunc(knot_vec_eta, j, order_eta, eta, &temp_output_eta, &temp_d_output_eta);
			for (k = index_min_zeta; k < index_max_zeta; k++)
			{
				lBasisFunc(knot_vec_zeta, k, cntl_p_n_zeta, order_zeta, zeta, &temp_output_zeta, &temp_d_output_zeta);

				temp_index = i + j * cntl_p_n_xi + k * cntl_p_n_xi * cntl_p_n_eta;

				temp1 = temp_output_xi   * temp_output_eta   * temp_output_zeta   * weight[Patch_controlpoint[0][temp_index]];
				temp2 = temp_d_output_xi * temp_output_eta   * temp_output_zeta   * weight[Patch_controlpoint[0][temp_index]];
				temp3 = temp_output_xi   * temp_d_output_eta * temp_output_zeta   * weight[Patch_controlpoint[0][temp_index]];
				temp4 = temp_output_xi   * temp_output_eta   * temp_d_output_zeta * weight[Patch_controlpoint[0][temp_index]];

				molecule_x += temp1 * cntl_px[Patch_controlpoint[0][temp_index]];
				molecule_y += temp1 * cntl_py[Patch_controlpoint[0][temp_index]];
				molecule_z += temp1 * cntl_pz[Patch_controlpoint[0][temp_index]];
				denominator += temp1;

				dxi_molecule_x += temp2 * cntl_px[Patch_controlpoint[0][temp_index]];
				dxi_molecule_y += temp2 * cntl_py[Patch_controlpoint[0][temp_index]];
				dxi_molecule_z += temp2 * cntl_pz[Patch_controlpoint[0][temp_index]];
				dxi_denominator += temp2;

				deta_molecule_x += temp3 * cntl_px[Patch_controlpoint[0][temp_index]];
				deta_molecule_y += temp3 * cntl_py[Patch_controlpoint[0][temp_index]];
				deta_molecule_z += temp3 * cntl_pz[Patch_controlpoint[0][temp_index]];
				deta_denominator += temp3;

				dzeta_molecule_x += temp4 * cntl_px[Patch_controlpoint[0][temp_index]];
				dzeta_molecule_y += temp4 * cntl_py[Patch_controlpoint[0][temp_index]];
				dzeta_molecule_z += temp4 * cntl_pz[Patch_controlpoint[0][temp_index]];
				dzeta_denominator += temp4;
			}
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;
	(*output_z) = molecule_z / denominator;

	//temp1 = denominator * denominator;
	(*output_dxi_x)   = weight[Patch_controlpoint[0][temp_index]] * ( dxi_molecule_x   * denominator - molecule_x * dxi_denominator)   / ( denominator * denominator ) ;
	(*output_dxi_y)   = weight[Patch_controlpoint[0][temp_index]] * ( dxi_molecule_y   * denominator - molecule_y * dxi_denominator)   / ( denominator * denominator ) ;
	(*output_dxi_z)   = weight[Patch_controlpoint[0][temp_index]] * ( dxi_molecule_z   * denominator - molecule_z * dxi_denominator)   / ( denominator * denominator ) ;
	(*output_deta_x)  = weight[Patch_controlpoint[0][temp_index]] * ( deta_molecule_x  * denominator - molecule_x * deta_denominator)  / ( denominator * denominator ) ;
	(*output_deta_y)  = weight[Patch_controlpoint[0][temp_index]] * ( deta_molecule_y  * denominator - molecule_y * deta_denominator)  / ( denominator * denominator ) ;
	(*output_deta_z)  = weight[Patch_controlpoint[0][temp_index]] * ( deta_molecule_z  * denominator - molecule_z * deta_denominator)  / ( denominator * denominator ) ;
	(*output_dzeta_x) = weight[Patch_controlpoint[0][temp_index]] * ( dzeta_molecule_x * denominator - molecule_x * dzeta_denominator) / ( denominator * denominator ) ;
	(*output_dzeta_y) = weight[Patch_controlpoint[0][temp_index]] * ( dzeta_molecule_y * denominator - molecule_y * dzeta_denominator) / ( denominator * denominator ) ;
	(*output_dzeta_z) = weight[Patch_controlpoint[0][temp_index]] * ( dzeta_molecule_z * denominator - molecule_z * dzeta_denominator) / ( denominator * denominator ) ;

	return denominator;
}

double rllNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta, double *cntl_px, double *cntl_py, double *cntl_pz, int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta, double *weight, int order_xi, int order_eta, int order_zeta, double xi, double eta, double zeta, double *output_x, double *output_y, double *output_z, double *output_dxi_x, double *output_deta_x, double *output_dzeta_x, double *output_dxi_y, double *output_deta_y, double *output_dzeta_y, double *output_dxi_z, double *output_deta_z, double *output_dzeta_z)
{
	int i, j, k, temp_index;
	double temp1, temp2, temp3, temp4;
	double molecule_x, molecule_y, molecule_z;
	double dxi_molecule_x, dxi_molecule_y, dxi_molecule_z;
	double deta_molecule_x, deta_molecule_y, deta_molecule_z;
	double dzeta_molecule_x, dzeta_molecule_y, dzeta_molecule_z;
	double denominator, dxi_denominator, deta_denominator, dzeta_denominator;
	double temp_output_xi, temp_output_eta, temp_output_zeta;
	double temp_d_output_xi, temp_d_output_eta, temp_d_output_zeta;

	molecule_x = 0.0;
	molecule_y = 0.0;
	molecule_z = 0.0;
	denominator = 0.0;

	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_molecule_z = 0.0;
	dxi_denominator = 0.0;

	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_molecule_z = 0.0;
	deta_denominator = 0.0;

	dzeta_molecule_x = 0.0;
	dzeta_molecule_y = 0.0;
	dzeta_molecule_z = 0.0;
	dzeta_denominator = 0.0;

	int index_min_xi   = 0;
	int index_max_xi   = cntl_p_n_xi - 1;
	int index_min_eta  = 0;
	int index_max_eta  = cntl_p_n_eta - 1;
	int index_min_zeta = 0;
	int index_max_zeta = cntl_p_n_zeta - 1;

	////////////////////////////////xi////////////////////////////////
	for (i = 0; i < cntl_p_n_xi; i++)
	{
		//printf("cntl_p_n_xi = %d\tknot_vec_xi[%d] = %lf\n",cntl_p_n_xi, i + 1, knot_vec_xi[i + 1]);
		if ( knot_vec_xi[i + 1] >= xi )
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			//printf("xi = %lf\tknot_vec_xi[%d] = %lf\n",xi ,i + 1, knot_vec_xi[i + 1]);
			break;
		}
	}

	if (index_min_xi < 0)
	{
		index_min_xi = 0;
		index_max_xi = order_xi + 1;
	}

	if (index_max_xi > cntl_p_n_xi)
	{
		index_min_xi = cntl_p_n_xi - order_xi - 1;
		index_max_xi = cntl_p_n_xi;
	}
	////////////////////////////////xi////////////////////////////////

	///////////////////////////////eta////////////////////////////////
	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if ( knot_vec_eta[i + 1] >= eta )
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			//printf("eta = %lf\tknot_vec_eta[%d] = %lf\n",eta ,i + 1, knot_vec_eta[i + 1]);
			break;
		}
	}
	if (index_min_eta < 0)
	{
		index_min_eta = 0;
		index_max_eta = order_eta + 1;
	}

	if (index_max_eta > cntl_p_n_eta)
	{
		index_min_eta = cntl_p_n_eta - order_eta - 1;
		index_max_eta = cntl_p_n_eta;
	}
	///////////////////////////////eta////////////////////////////////

	///////////////////////////////zeta///////////////////////////////
	for (i = 0; i < cntl_p_n_zeta; i++)
	{
		if ( knot_vec_zeta[i + 1] >= zeta )
		{
			index_min_zeta = i - order_zeta;
			index_max_zeta = i + 1;
			//printf("zeta = %lf\tknot_vec_zeta[%d] = %lf\n",zeta ,i + 1, knot_vec_zeta[i + 1]);
			break;
		}
	}
	if (index_min_zeta < 0)
	{
		index_min_zeta = 0;
		index_max_zeta = order_zeta + 1;
	}

	if (index_max_zeta > cntl_p_n_zeta)
	{
		index_min_zeta = cntl_p_n_zeta - order_zeta - 1;
		index_max_zeta = cntl_p_n_zeta;
	}
	///////////////////////////////zeta///////////////////////////////



	for (i = index_min_xi; i < index_max_xi; i++)
	{
		rBasisFunc(knot_vec_xi, i, order_xi, xi, &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j < index_max_eta; j++)
		{
			lBasisFunc(knot_vec_eta, j, cntl_p_n_eta, order_eta, eta, &temp_output_eta, &temp_d_output_eta);
			for (k = index_min_zeta; k < index_max_zeta; k++)
			{
				lBasisFunc(knot_vec_zeta, k, cntl_p_n_zeta, order_zeta, zeta, &temp_output_zeta, &temp_d_output_zeta);

				temp_index = i + j * cntl_p_n_xi + k * cntl_p_n_xi * cntl_p_n_eta;

				temp1 = temp_output_xi   * temp_output_eta   * temp_output_zeta   * weight[Patch_controlpoint[0][temp_index]];
				temp2 = temp_d_output_xi * temp_output_eta   * temp_output_zeta   * weight[Patch_controlpoint[0][temp_index]];
				temp3 = temp_output_xi   * temp_d_output_eta * temp_output_zeta   * weight[Patch_controlpoint[0][temp_index]];
				temp4 = temp_output_xi   * temp_output_eta   * temp_d_output_zeta * weight[Patch_controlpoint[0][temp_index]];

				molecule_x += temp1 * cntl_px[Patch_controlpoint[0][temp_index]];
				molecule_y += temp1 * cntl_py[Patch_controlpoint[0][temp_index]];
				molecule_z += temp1 * cntl_pz[Patch_controlpoint[0][temp_index]];
				denominator += temp1;

				dxi_molecule_x += temp2 * cntl_px[Patch_controlpoint[0][temp_index]];
				dxi_molecule_y += temp2 * cntl_py[Patch_controlpoint[0][temp_index]];
				dxi_molecule_z += temp2 * cntl_pz[Patch_controlpoint[0][temp_index]];
				dxi_denominator += temp2;

				deta_molecule_x += temp3 * cntl_px[Patch_controlpoint[0][temp_index]];
				deta_molecule_y += temp3 * cntl_py[Patch_controlpoint[0][temp_index]];
				deta_molecule_z += temp3 * cntl_pz[Patch_controlpoint[0][temp_index]];
				deta_denominator += temp3;

				dzeta_molecule_x += temp4 * cntl_px[temp_index];
				dzeta_molecule_y += temp4 * cntl_py[temp_index];
				dzeta_molecule_z += temp4 * cntl_pz[temp_index];
				dzeta_denominator += temp4;
			}
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;
	(*output_z) = molecule_z / denominator;

	//temp1 = denominator * denominator;
	(*output_dxi_x)   = weight[Patch_controlpoint[0][temp_index]] * ( dxi_molecule_x   * denominator - molecule_x * dxi_denominator)   / ( denominator * denominator ) ;
	(*output_dxi_y)   = weight[Patch_controlpoint[0][temp_index]] * ( dxi_molecule_y   * denominator - molecule_y * dxi_denominator)   / ( denominator * denominator ) ;
	(*output_dxi_z)   = weight[Patch_controlpoint[0][temp_index]] * ( dxi_molecule_z   * denominator - molecule_z * dxi_denominator)   / ( denominator * denominator ) ;
	(*output_deta_x)  = weight[Patch_controlpoint[0][temp_index]] * ( deta_molecule_x  * denominator - molecule_x * deta_denominator)  / ( denominator * denominator ) ;
	(*output_deta_y)  = weight[Patch_controlpoint[0][temp_index]] * ( deta_molecule_y  * denominator - molecule_y * deta_denominator)  / ( denominator * denominator ) ;
	(*output_deta_z)  = weight[Patch_controlpoint[0][temp_index]] * ( deta_molecule_z  * denominator - molecule_z * deta_denominator)  / ( denominator * denominator ) ;
	(*output_dzeta_x) = weight[Patch_controlpoint[0][temp_index]] * ( dzeta_molecule_x * denominator - molecule_x * dzeta_denominator) / ( denominator * denominator ) ;
	(*output_dzeta_y) = weight[Patch_controlpoint[0][temp_index]] * ( dzeta_molecule_y * denominator - molecule_y * dzeta_denominator) / ( denominator * denominator ) ;
	(*output_dzeta_z) = weight[Patch_controlpoint[0][temp_index]] * ( dzeta_molecule_z * denominator - molecule_z * dzeta_denominator) / ( denominator * denominator ) ;

	return denominator;
}

double lllNURBS_volume(double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta, double *cntl_px, double *cntl_py, double *cntl_pz, int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta, double *weight, int order_xi, int order_eta, int order_zeta, double xi, double eta, double zeta, double *output_x, double *output_y, double *output_z, double *output_dxi_x, double *output_deta_x, double *output_dzeta_x, double *output_dxi_y, double *output_deta_y, double *output_dzeta_y, double *output_dxi_z, double *output_deta_z, double *output_dzeta_z)
{
	int i, j, k, temp_index;
	double temp1, temp2, temp3, temp4;
	double molecule_x, molecule_y, molecule_z;
	double dxi_molecule_x, dxi_molecule_y, dxi_molecule_z;
	double deta_molecule_x, deta_molecule_y, deta_molecule_z;
	double dzeta_molecule_x, dzeta_molecule_y, dzeta_molecule_z;
	double denominator, dxi_denominator, deta_denominator, dzeta_denominator;
	double temp_output_xi, temp_output_eta, temp_output_zeta;
	double temp_d_output_xi, temp_d_output_eta, temp_d_output_zeta;

	molecule_x = 0.0;
	molecule_y = 0.0;
	molecule_z = 0.0;
	denominator = 0.0;

	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_molecule_z = 0.0;
	dxi_denominator = 0.0;

	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_molecule_z = 0.0;
	deta_denominator = 0.0;

	dzeta_molecule_x = 0.0;
	dzeta_molecule_y = 0.0;
	dzeta_molecule_z = 0.0;
	dzeta_denominator = 0.0;

	int index_min_xi   = 0;
	int index_max_xi   = cntl_p_n_xi - 1;
	int index_min_eta  = 0;
	int index_max_eta  = cntl_p_n_eta - 1;
	int index_min_zeta = 0;
	int index_max_zeta = cntl_p_n_zeta - 1;

	////////////////////////////////xi////////////////////////////////
	for (i = 0; i < cntl_p_n_xi; i++)
	{
		//printf("cntl_p_n_xi = %d\tknot_vec_xi[%d] = %lf\n",cntl_p_n_xi, i + 1, knot_vec_xi[i + 1]);
		if ( knot_vec_xi[i + 1] >= xi )
		{
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			//printf("xi = %lf\tknot_vec_xi[%d] = %lf\n",xi ,i + 1, knot_vec_xi[i + 1]);
			break;
		}
	}

	if (index_min_xi < 0)
	{
		index_min_xi = 0;
		index_max_xi = order_xi + 1;
	}

	if (index_max_xi > cntl_p_n_xi)
	{
		index_min_xi = cntl_p_n_xi - order_xi - 1;
		index_max_xi = cntl_p_n_xi;
	}
	////////////////////////////////xi////////////////////////////////

	///////////////////////////////eta////////////////////////////////
	for (i = 0; i < cntl_p_n_eta; i++)
	{
		if ( knot_vec_eta[i + 1] >= eta )
		{
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			//printf("eta = %lf\tknot_vec_eta[%d] = %lf\n",eta ,i + 1, knot_vec_eta[i + 1]);
			break;
		}
	}
	if (index_min_eta < 0)
	{
		index_min_eta = 0;
		index_max_eta = order_eta + 1;
	}

	if (index_max_eta > cntl_p_n_eta)
	{
		index_min_eta = cntl_p_n_eta - order_eta - 1;
		index_max_eta = cntl_p_n_eta;
	}
	///////////////////////////////eta////////////////////////////////

	///////////////////////////////zeta///////////////////////////////
	for (i = 0; i < cntl_p_n_zeta; i++)
	{
		if ( knot_vec_zeta[i + 1] >= zeta )
		{
			index_min_zeta = i - order_zeta;
			index_max_zeta = i + 1;
			//printf("zeta = %lf\tknot_vec_zeta[%d] = %lf\n",zeta ,i + 1, knot_vec_zeta[i + 1]);
			break;
		}
	}
	if (index_min_zeta < 0)
	{
		index_min_zeta = 0;
		index_max_zeta = order_zeta + 1;
	}

	if (index_max_zeta > cntl_p_n_zeta)
	{
		index_min_zeta = cntl_p_n_zeta - order_zeta - 1;
		index_max_zeta = cntl_p_n_zeta;
	}
	///////////////////////////////zeta///////////////////////////////



	for (i = index_min_xi; i < index_max_xi; i++)
	{
		lBasisFunc(knot_vec_xi, i, cntl_p_n_xi, order_xi, xi, &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j < index_max_eta; j++)
		{
			lBasisFunc(knot_vec_eta, j, cntl_p_n_eta, order_eta, eta, &temp_output_eta, &temp_d_output_eta);
			for (k = index_min_zeta; k < index_max_zeta; k++)
			{
				lBasisFunc(knot_vec_zeta, k, cntl_p_n_zeta, order_zeta, zeta, &temp_output_zeta, &temp_d_output_zeta);

				temp_index = i + j * cntl_p_n_xi + k * cntl_p_n_xi * cntl_p_n_eta;

				temp1 = temp_output_xi   * temp_output_eta   * temp_output_zeta   * weight[Patch_controlpoint[0][temp_index]];
				temp2 = temp_d_output_xi * temp_output_eta   * temp_output_zeta   * weight[Patch_controlpoint[0][temp_index]];
				temp3 = temp_output_xi   * temp_d_output_eta * temp_output_zeta   * weight[Patch_controlpoint[0][temp_index]];
				temp4 = temp_output_xi   * temp_output_eta   * temp_d_output_zeta * weight[Patch_controlpoint[0][temp_index]];

				molecule_x += temp1 * cntl_px[Patch_controlpoint[0][temp_index]];
				molecule_y += temp1 * cntl_py[Patch_controlpoint[0][temp_index]];
				molecule_z += temp1 * cntl_pz[Patch_controlpoint[0][temp_index]];
				denominator += temp1;

				dxi_molecule_x += temp2 * cntl_px[Patch_controlpoint[0][temp_index]];
				dxi_molecule_y += temp2 * cntl_py[Patch_controlpoint[0][temp_index]];
				dxi_molecule_z += temp2 * cntl_pz[Patch_controlpoint[0][temp_index]];
				dxi_denominator += temp2;

				deta_molecule_x += temp3 * cntl_px[Patch_controlpoint[0][temp_index]];
				deta_molecule_y += temp3 * cntl_py[Patch_controlpoint[0][temp_index]];
				deta_molecule_z += temp3 * cntl_pz[Patch_controlpoint[0][temp_index]];
				deta_denominator += temp3;

				dzeta_molecule_x += temp4 * cntl_px[Patch_controlpoint[0][temp_index]];
				dzeta_molecule_y += temp4 * cntl_py[Patch_controlpoint[0][temp_index]];
				dzeta_molecule_z += temp4 * cntl_pz[Patch_controlpoint[0][temp_index]];
				dzeta_denominator += temp4;
			}
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;
	(*output_z) = molecule_z / denominator;

	//temp1 = denominator * denominator;
	(*output_dxi_x)   = weight[Patch_controlpoint[0][temp_index]] * ( dxi_molecule_x   * denominator - molecule_x * dxi_denominator)   / ( denominator * denominator ) ;
	(*output_dxi_y)   = weight[Patch_controlpoint[0][temp_index]] * ( dxi_molecule_y   * denominator - molecule_y * dxi_denominator)   / ( denominator * denominator ) ;
	(*output_dxi_z)   = weight[Patch_controlpoint[0][temp_index]] * ( dxi_molecule_z   * denominator - molecule_z * dxi_denominator)   / ( denominator * denominator ) ;
	(*output_deta_x)  = weight[Patch_controlpoint[0][temp_index]] * ( deta_molecule_x  * denominator - molecule_x * deta_denominator)  / ( denominator * denominator ) ;
	(*output_deta_y)  = weight[Patch_controlpoint[0][temp_index]] * ( deta_molecule_y  * denominator - molecule_y * deta_denominator)  / ( denominator * denominator ) ;
	(*output_deta_z)  = weight[Patch_controlpoint[0][temp_index]] * ( deta_molecule_z  * denominator - molecule_z * deta_denominator)  / ( denominator * denominator ) ;
	(*output_dzeta_x) = weight[Patch_controlpoint[0][temp_index]] * ( dzeta_molecule_x * denominator - molecule_x * dzeta_denominator) / ( denominator * denominator ) ;
	(*output_dzeta_y) = weight[Patch_controlpoint[0][temp_index]] * ( dzeta_molecule_y * denominator - molecule_y * dzeta_denominator) / ( denominator * denominator ) ;
	(*output_dzeta_z) = weight[Patch_controlpoint[0][temp_index]] * ( dzeta_molecule_z * denominator - molecule_z * dzeta_denominator) / ( denominator * denominator ) ;

	return denominator;
}

//算出したローカルパッチ各要素の頂点の物理座標のグローバルパッチでの(xi,eta)算出
int Calc_xi_eta_zeta(double px, double py, double pz, double *knot_vec_xi, double *knot_vec_eta, double *knot_vec_zeta, double *cntl_px, double *cntl_py, double *cntl_pz, int cntl_p_n_xi, int cntl_p_n_eta, int cntl_p_n_zeta, double *weight, int order_xi, int order_eta, int order_zeta, double *output_xi, double *output_eta, double *output_zeta)
{
	double temp_xi, temp_eta, temp_zeta;
	double temp_x, temp_y, temp_z;
	double temp_matrix[3][3];
	double temp_dxi, temp_deta, temp_dzeta;
    double temp_tol_x, temp_tol_y, temp_tol_z;

	(*output_xi)   = 0.0;
	(*output_eta)  = 0.0;
	(*output_zeta) = 0.0;

	int i, j;
	int repeat  = 5;
	int repeat2 = 5;
	double tol = 10e-22;
	double coef = 0.0;




	/////////////////////////////////////////////////rrrNURBS_volume/////////////////////////////////////////////////
	for (j = 0; j < repeat2; j++)
	{
		//初期値の設定
		if (j == 0)
		{
			temp_xi = knot_vec_xi[0] + knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= 0.5;

			temp_eta = knot_vec_eta[0] + knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= 0.5;

			temp_zeta = knot_vec_zeta[0] + knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= 0.5;
		}
		else
		{
			coef = j / repeat2;

			temp_xi = knot_vec_xi[0] + knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= coef;

			temp_eta = knot_vec_eta[0] + knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= coef;

			temp_zeta = knot_vec_zeta[0] + knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= coef;
		}
		//printf("j = %d\txi = %lf\tcoef = %lf\n", j, temp_xi, coef);

		for (i = 0; i < repeat; i++)
		{
			//printf("-------------- i = %d--------------\n", i);
			rrrNURBS_volume(Position_Knots[0][0], Position_Knots[0][1], Position_Knots[0][2], Control_Coord[0], Control_Coord[1], Control_Coord[2], No_Control_point[0][0], No_Control_point[0][1], No_Control_point[0][2] 	, Control_Weight, Order[0][0], Order[0][1], Order[0][2], temp_xi, temp_eta, temp_zeta, &temp_x, &temp_y, &temp_z, &temp_matrix[0][0], &temp_matrix[0][1], &temp_matrix[0][2], &temp_matrix[1][0], &temp_matrix	[1][1], &temp_matrix[1][2], &temp_matrix[2][0], &temp_matrix[2][1], &temp_matrix[2][2]);

			temp_tol_x = px - temp_x;
			temp_tol_x *= temp_tol_x;

			temp_tol_y = py - temp_y;
			temp_tol_y *= temp_tol_y;

			temp_tol_z = pz - temp_z;
			temp_tol_z *= temp_tol_z;

    	    //printf("xi_0:  % 1.8e\n", temp_xi);
			//printf("eta_0: % 1.8e\n", temp_eta);
    	    //printf("px: % 1.8e\n",px);
    	    //printf("temp_x: % 1.8e\n",temp_x);
    	    //printf("py: % 1.8e\n",py);
    	    //printf("temp_y: % 1.8e\n",temp_y);
    	    //printf("temp_tol_x:  % 1.8e\n", temp_tol_x);
			//printf("temp_tol_y:  % 1.8e\n", temp_tol_y);

			//if (i == repeat - 1)
			//{
			//	printf("rrrNURBS_volume\n");
			//	printf("temp_tol_x:  % 1.8e\n", temp_tol_x);
			//	printf("temp_tol_y:  % 1.8e\n", temp_tol_y);
			//	printf("temp_tol_z:  % 1.8e\n", temp_tol_z);
			//}
			//printf("temp_tol_x:  % 1.8e\n", temp_tol_x);
			//printf("temp_tol_y:  % 1.8e\n", temp_tol_y);
			//printf("temp_tol_z:  % 1.8e\n", temp_tol_z);

			//収束した場合////////////////////////////////////////////////////////////////
    	    //if (temp_tol_x + temp_tol_y + temp_tol_z < tol && i != 0) //収束条件は考慮する必要あり
			if (temp_tol_x < tol && temp_tol_y < tol && temp_tol_z < tol && i != 0 /*&& i != 1*repeat/4 && i != 2*repeat/4 && i != 3*repeat/4*/)
			{
				(*output_xi)   = temp_xi;
				(*output_eta)  = temp_eta;
				(*output_zeta) = temp_zeta;

				//printf("j = %d\tconvergence:rrrNURBS_volume\n", j);
				//printf("temp_tol_x:  % 1.8e\n", temp_tol_x);
				//printf("temp_tol_y:  % 1.8e\n", temp_tol_y);
				//printf("temp_tol_z:  % 1.8e\n", temp_tol_z);

				return i;
			}

			InverseMatrix(temp_matrix);

			temp_dxi   = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y) + temp_matrix[0][2] * (pz - temp_z);
			temp_deta  = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y) + temp_matrix[1][2] * (pz - temp_z);
			temp_dzeta = temp_matrix[2][0] * (px - temp_x) + temp_matrix[2][1] * (py - temp_y) + temp_matrix[2][2] * (pz - temp_z);



			temp_xi   = temp_xi + temp_dxi;
			temp_eta  = temp_eta + temp_deta;
			temp_zeta = temp_zeta + temp_dzeta;

			if (temp_xi < knot_vec_xi[0])
				temp_xi = knot_vec_xi[0];
			if (temp_xi > knot_vec_xi[cntl_p_n_xi + order_xi])
				temp_xi = knot_vec_xi[cntl_p_n_xi + order_xi];

			if (temp_eta < knot_vec_eta[0])
				temp_eta = knot_vec_eta[0];
			if (temp_eta > knot_vec_eta[cntl_p_n_eta + order_eta])
				temp_eta = knot_vec_eta[cntl_p_n_eta + order_eta];

			if (temp_zeta < knot_vec_zeta[0])
				temp_zeta = knot_vec_zeta[0];
			if (temp_zeta > knot_vec_zeta[cntl_p_n_zeta + order_zeta])
				temp_zeta = knot_vec_zeta[cntl_p_n_zeta + order_zeta];

			//printf("j=%d\ti=%d\n",j,i);
			//printf("r_dxi:  % 1.8e\n", temp_dxi  );
			//printf("r_deta: % 1.8e\n", temp_deta );
			//printf("r_dzeta: % 1.8e\n",temp_dzeta);
    	    //printf("r_xi:  % 1.8e\n", temp_xi);
			//printf("r_eta: % 1.8e\n", temp_eta);
			//printf("r_zeta: % 1.8e\n", temp_zeta);

			//double temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
			//printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
		}
	}

	/////////////////////////////////////////////////lrrNURBS_volume/////////////////////////////////////////////////
	for (j = 0; j < repeat2; j++)
	{
		//初期値の設定
		if (j == 0)
		{
			temp_xi = knot_vec_xi[0] + knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= 0.5;

			temp_eta = knot_vec_eta[0] + knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= 0.5;

			temp_zeta = knot_vec_zeta[0] + knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= 0.5;
		}
		else
		{
			coef = j / repeat2;

			temp_xi = knot_vec_xi[0] + knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= coef;

			temp_eta = knot_vec_eta[0] + knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= coef;

			temp_zeta = knot_vec_zeta[0] + knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= coef;
		}
		//printf("j = %d\txi = %lf\tcoef = %lf\n", j, temp_xi, coef);

		for (i = 0; i < repeat; i++)
		{
			lrrNURBS_volume(Position_Knots[0][0], Position_Knots[0][1], Position_Knots[0][2], Control_Coord[0], Control_Coord[1], Control_Coord[2], No_Control_point[0][0], No_Control_point[0][1], No_Control_point[0][2] 	, Control_Weight, Order[0][0], Order[0][1], Order[0][2], temp_xi, temp_eta, temp_zeta, &temp_x, &temp_y, &temp_z, &temp_matrix[0][0], &temp_matrix[0][1], &temp_matrix[0][2], &temp_matrix[1][0], &temp_matrix	[1][1], &temp_matrix[1][2], &temp_matrix[2][0], &temp_matrix[2][1], &temp_matrix[2][2]);

			temp_tol_x = px - temp_x;
			temp_tol_x *= temp_tol_x;

			temp_tol_y = py - temp_y;
			temp_tol_y *= temp_tol_y;

			temp_tol_z = pz - temp_z;
			temp_tol_z *= temp_tol_z;

    	    //printf("xi_0:  % 1.8e\n", temp_xi);
			//printf("eta_0: % 1.8e\n", temp_eta);
    	    //printf("px: % 1.8e\n",px);
    	    //printf("temp_x: % 1.8e\n",temp_x);
    	    //printf("py: % 1.8e\n",py);
    	    //printf("temp_y: % 1.8e\n",temp_y);
    	    //printf("temp_tol_x:  % 1.8e\n", temp_tol_x);
			//printf("temp_tol_y:  % 1.8e\n", temp_tol_y);

			/*if (i == repeat - 1)
			{
				printf("lrrNURBS_volume\n");
				printf("temp_tol_x:  % 1.8e\n", temp_tol_x);
				printf("temp_tol_y:  % 1.8e\n", temp_tol_y);
				printf("temp_tol_z:  % 1.8e\n", temp_tol_z);
			}*/

			//収束した場合////////////////////////////////////////////////////////////////
    	    //if (temp_tol_x + temp_tol_y + temp_tol_z < tol && i != 1*repeat/4 && i != 2*repeat/4 && i != 3*repeat/4) //収束条件は考慮する必要あり
			if (temp_tol_x < tol && temp_tol_y < tol && temp_tol_z < tol && i != 0 /*&& i != 1*repeat/4 && i != 2*repeat/4 && i != 3*repeat/4*/)
			{
				(*output_xi)   = temp_xi;
				(*output_eta)  = temp_eta;
				(*output_zeta) = temp_zeta;

				//printf("j = %d\tconvergence:lrrNURBS_volume\n", j);
				//printf("temp_tol_x:  % 1.8e\n", temp_tol_x);
				//printf("temp_tol_y:  % 1.8e\n", temp_tol_y);
				//printf("temp_tol_z:  % 1.8e\n", temp_tol_z);

				return i;
			}

			InverseMatrix(temp_matrix);

			temp_dxi   = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y) + temp_matrix[0][2] * (pz - temp_z);
			temp_deta  = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y) + temp_matrix[1][2] * (pz - temp_z);
			temp_dzeta = temp_matrix[2][0] * (px - temp_x) + temp_matrix[2][1] * (py - temp_y) + temp_matrix[2][2] * (pz - temp_z);

			temp_xi   = temp_xi + temp_dxi;
			temp_eta  = temp_eta + temp_deta;
			temp_zeta = temp_zeta + temp_dzeta;

			if (temp_xi < knot_vec_xi[0])
				temp_xi = knot_vec_xi[0];
			if (temp_xi > knot_vec_xi[cntl_p_n_xi + order_xi])
				temp_xi = knot_vec_xi[cntl_p_n_xi + order_xi];

			if (temp_eta < knot_vec_eta[0])
				temp_eta = knot_vec_eta[0];
			if (temp_eta > knot_vec_eta[cntl_p_n_eta + order_eta])
				temp_eta = knot_vec_eta[cntl_p_n_eta + order_eta];

			if (temp_zeta < knot_vec_zeta[0])
				temp_zeta = knot_vec_zeta[0];
			if (temp_zeta > knot_vec_zeta[cntl_p_n_zeta + order_zeta])
				temp_zeta = knot_vec_zeta[cntl_p_n_zeta + order_zeta];

    	    //printf("r_xi:  % 1.8e\n", temp_xi);
			//printf("r_eta: % 1.8e\n", temp_eta);
    	    //printf("i=%d\n",i);

			//double temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
			//printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
		}
	}

	/////////////////////////////////////////////////rlrNURBS_volume/////////////////////////////////////////////////
	for (j = 0; j < repeat2; j++)
	{
		//初期値の設定
		if (j == 0)
		{
			temp_xi = knot_vec_xi[0] + knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= 0.5;

			temp_eta = knot_vec_eta[0] + knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= 0.5;

			temp_zeta = knot_vec_zeta[0] + knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= 0.5;
		}
		else
		{
			coef = j / repeat2;

			temp_xi = knot_vec_xi[0] + knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= coef;

			temp_eta = knot_vec_eta[0] + knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= coef;

			temp_zeta = knot_vec_zeta[0] + knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= coef;
		}
		//printf("j = %d\txi = %lf\tcoef = %lf\n", j, temp_xi, coef);

		for (i = 0; i < repeat; i++)
		{
			rlrNURBS_volume(Position_Knots[0][0], Position_Knots[0][1], Position_Knots[0][2], Control_Coord[0], Control_Coord[1], Control_Coord[2], No_Control_point[0][0], No_Control_point[0][1], No_Control_point[0][2] 	, Control_Weight, Order[0][0], Order[0][1], Order[0][2], temp_xi, temp_eta, temp_zeta, &temp_x, &temp_y, &temp_z, &temp_matrix[0][0], &temp_matrix[0][1], &temp_matrix[0][2], &temp_matrix[1][0], &temp_matrix	[1][1], &temp_matrix[1][2], &temp_matrix[2][0], &temp_matrix[2][1], &temp_matrix[2][2]);

			temp_tol_x = px - temp_x;
			temp_tol_x *= temp_tol_x;

			temp_tol_y = py - temp_y;
			temp_tol_y *= temp_tol_y;

			temp_tol_z = pz - temp_z;
			temp_tol_z *= temp_tol_z;

    	    //printf("xi_0:  % 1.8e\n", temp_xi);
			//printf("eta_0: % 1.8e\n", temp_eta);
    	    //printf("px: % 1.8e\n",px);
    	    //printf("temp_x: % 1.8e\n",temp_x);
    	    //printf("py: % 1.8e\n",py);
    	    //printf("temp_y: % 1.8e\n",temp_y);
    	    //printf("temp_tol_x:  % 1.8e\n", temp_tol_x);
			//printf("temp_tol_y:  % 1.8e\n", temp_tol_y);

			/*if (i == repeat - 1)
			{
				printf("rlrNURBS_volume\n");
				printf("temp_tol_x:  % 1.8e\n", temp_tol_x);
				printf("temp_tol_y:  % 1.8e\n", temp_tol_y);
				printf("temp_tol_z:  % 1.8e\n", temp_tol_z);
			}*/

			//収束した場合////////////////////////////////////////////////////////////////
    	    //if (temp_tol_x + temp_tol_y + temp_tol_z < tol && i != 1*repeat/4 && i != 2*repeat/4 && i != 3*repeat/4) //収束条件は考慮する必要あり
			if (temp_tol_x < tol && temp_tol_y < tol && temp_tol_z < tol && i != 0 /*&& i != 1*repeat/4 && i != 2*repeat/4 && i != 3*repeat/4*/)
			{
				(*output_xi)   = temp_xi;
				(*output_eta)  = temp_eta;
				(*output_zeta) = temp_zeta;

				//printf("j = %d\tconvergence:rlrNURBS_volume\n", j);
				//printf("temp_tol_x:  % 1.8e\n", temp_tol_x);
				//printf("temp_tol_y:  % 1.8e\n", temp_tol_y);
				//printf("temp_tol_z:  % 1.8e\n", temp_tol_z);

				return i;
			}

			InverseMatrix(temp_matrix);

			temp_dxi   = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y) + temp_matrix[0][2] * (pz - temp_z); //ここ違うかも，確認
			temp_deta  = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y) + temp_matrix[1][2] * (pz - temp_z);
			temp_dzeta = temp_matrix[2][0] * (px - temp_x) + temp_matrix[2][1] * (py - temp_y) + temp_matrix[2][2] * (pz - temp_z);

			temp_xi   = temp_xi + temp_dxi;
			temp_eta  = temp_eta + temp_deta;
			temp_zeta = temp_zeta + temp_dzeta;

			if (temp_xi < knot_vec_xi[0])
				temp_xi = knot_vec_xi[0];
			if (temp_xi > knot_vec_xi[cntl_p_n_xi + order_xi])
				temp_xi = knot_vec_xi[cntl_p_n_xi + order_xi];

			if (temp_eta < knot_vec_eta[0])
				temp_eta = knot_vec_eta[0];
			if (temp_eta > knot_vec_eta[cntl_p_n_eta + order_eta])
				temp_eta = knot_vec_eta[cntl_p_n_eta + order_eta];

			if (temp_zeta < knot_vec_zeta[0])
				temp_zeta = knot_vec_zeta[0];
			if (temp_zeta > knot_vec_zeta[cntl_p_n_zeta + order_zeta])
				temp_zeta = knot_vec_zeta[cntl_p_n_zeta + order_zeta];

    	    //printf("r_xi:  % 1.8e\n", temp_xi);
			//printf("r_eta: % 1.8e\n", temp_eta);
    	    //printf("i=%d\n",i);

			//double temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
			//printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
		}
	}

	/////////////////////////////////////////////////rrlNURBS_volume/////////////////////////////////////////////////
	for (j = 0; j < repeat2; j++)
	{
		//初期値の設定
		if (j == 0)
		{
			temp_xi = knot_vec_xi[0] + knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= 0.5;

			temp_eta = knot_vec_eta[0] + knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= 0.5;

			temp_zeta = knot_vec_zeta[0] + knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= 0.5;
		}
		else
		{
			coef = j / repeat2;

			temp_xi = knot_vec_xi[0] + knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= coef;

			temp_eta = knot_vec_eta[0] + knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= coef;

			temp_zeta = knot_vec_zeta[0] + knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= coef;
		}
		//printf("j = %d\txi = %lf\tcoef = %lf\n", j, temp_xi, coef);

		for (i = 0; i < repeat; i++)
		{
			rrlNURBS_volume(Position_Knots[0][0], Position_Knots[0][1], Position_Knots[0][2], Control_Coord[0], Control_Coord[1], Control_Coord[2], No_Control_point[0][0], No_Control_point[0][1], No_Control_point[0][2] 	, Control_Weight, Order[0][0], Order[0][1], Order[0][2], temp_xi, temp_eta, temp_zeta, &temp_x, &temp_y, &temp_z, &temp_matrix[0][0], &temp_matrix[0][1], &temp_matrix[0][2], &temp_matrix[1][0], &temp_matrix	[1][1], &temp_matrix[1][2], &temp_matrix[2][0], &temp_matrix[2][1], &temp_matrix[2][2]);

			temp_tol_x = px - temp_x;
			temp_tol_x *= temp_tol_x;

			temp_tol_y = py - temp_y;
			temp_tol_y *= temp_tol_y;

			temp_tol_z = pz - temp_z;
			temp_tol_z *= temp_tol_z;

    	    //printf("xi_0:  % 1.8e\n", temp_xi);
			//printf("eta_0: % 1.8e\n", temp_eta);
    	    //printf("px: % 1.8e\n",px);
    	    //printf("temp_x: % 1.8e\n",temp_x);
    	    //printf("py: % 1.8e\n",py);
    	    //printf("temp_y: % 1.8e\n",temp_y);
    	    //printf("temp_tol_x:  % 1.8e\n", temp_tol_x);
			//printf("temp_tol_y:  % 1.8e\n", temp_tol_y);

			/*if (i == repeat - 1)
			{
				printf("rrlNURBS_volume\n");
				printf("temp_tol_x:  % 1.8e\n", temp_tol_x);
				printf("temp_tol_y:  % 1.8e\n", temp_tol_y);
				printf("temp_tol_z:  % 1.8e\n", temp_tol_z);
			}*/

			//収束した場合////////////////////////////////////////////////////////////////
    	    //if (temp_tol_x + temp_tol_y + temp_tol_z < tol && i != 1*repeat/4 && i != 2*repeat/4 && i != 3*repeat/4) //収束条件は考慮する必要あり
			if (temp_tol_x < tol && temp_tol_y < tol && temp_tol_z < tol && i != 0 /*&& i != 1*repeat/4 && i != 2*repeat/4 && i != 3*repeat/4*/)
			{
				(*output_xi)   = temp_xi;
				(*output_eta)  = temp_eta;
				(*output_zeta) = temp_zeta;

				//printf("j = %d\tconvergence:rrlNURBS_volume\n", j);
				//printf("temp_tol_x:  % 1.8e\n", temp_tol_x);
				//printf("temp_tol_y:  % 1.8e\n", temp_tol_y);
				//printf("temp_tol_z:  % 1.8e\n", temp_tol_z);

				return i;
			}

			InverseMatrix(temp_matrix);

			temp_dxi   = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y) + temp_matrix[0][2] * (pz - temp_z); //ここ違うかも，確認
			temp_deta  = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y) + temp_matrix[1][2] * (pz - temp_z);
			temp_dzeta = temp_matrix[2][0] * (px - temp_x) + temp_matrix[2][1] * (py - temp_y) + temp_matrix[2][2] * (pz - temp_z);

			temp_xi   = temp_xi + temp_dxi;
			temp_eta  = temp_eta + temp_deta;
			temp_zeta = temp_zeta + temp_dzeta;

			if (temp_xi < knot_vec_xi[0])
				temp_xi = knot_vec_xi[0];
			if (temp_xi > knot_vec_xi[cntl_p_n_xi + order_xi])
				temp_xi = knot_vec_xi[cntl_p_n_xi + order_xi];

			if (temp_eta < knot_vec_eta[0])
				temp_eta = knot_vec_eta[0];
			if (temp_eta > knot_vec_eta[cntl_p_n_eta + order_eta])
				temp_eta = knot_vec_eta[cntl_p_n_eta + order_eta];

			if (temp_zeta < knot_vec_zeta[0])
				temp_zeta = knot_vec_zeta[0];
			if (temp_zeta > knot_vec_zeta[cntl_p_n_zeta + order_zeta])
				temp_zeta = knot_vec_zeta[cntl_p_n_zeta + order_zeta];

    	    //printf("r_xi:  % 1.8e\n", temp_xi);
			//printf("r_eta: % 1.8e\n", temp_eta);
    	    //printf("i=%d\n",i);

			//double temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
			//printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
		}
	}

	/////////////////////////////////////////////////llrNURBS_volume/////////////////////////////////////////////////
	for (j = 0; j < repeat2; j++)
	{
		//初期値の設定
		if (j == 0)
		{
			temp_xi = knot_vec_xi[0] + knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= 0.5;

			temp_eta = knot_vec_eta[0] + knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= 0.5;

			temp_zeta = knot_vec_zeta[0] + knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= 0.5;
		}
		else
		{
			coef = j / repeat2;

			temp_xi = knot_vec_xi[0] + knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= coef;

			temp_eta = knot_vec_eta[0] + knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= coef;

			temp_zeta = knot_vec_zeta[0] + knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= coef;
		}
		//printf("j = %d\txi = %lf\tcoef = %lf\n", j, temp_xi, coef);

		for (i = 0; i < repeat; i++)
		{
			llrNURBS_volume(Position_Knots[0][0], Position_Knots[0][1], Position_Knots[0][2], Control_Coord[0], Control_Coord[1], Control_Coord[2], No_Control_point[0][0], No_Control_point[0][1], No_Control_point[0][2] 	, Control_Weight, Order[0][0], Order[0][1], Order[0][2], temp_xi, temp_eta, temp_zeta, &temp_x, &temp_y, &temp_z, &temp_matrix[0][0], &temp_matrix[0][1], &temp_matrix[0][2], &temp_matrix[1][0], &temp_matrix	[1][1], &temp_matrix[1][2], &temp_matrix[2][0], &temp_matrix[2][1], &temp_matrix[2][2]);

			temp_tol_x = px - temp_x;
			temp_tol_x *= temp_tol_x;

			temp_tol_y = py - temp_y;
			temp_tol_y *= temp_tol_y;

			temp_tol_z = pz - temp_z;
			temp_tol_z *= temp_tol_z;

    	    //printf("xi_0:  % 1.8e\n", temp_xi);
			//printf("eta_0: % 1.8e\n", temp_eta);
    	    //printf("px: % 1.8e\n",px);
    	    //printf("temp_x: % 1.8e\n",temp_x);
    	    //printf("py: % 1.8e\n",py);
    	    //printf("temp_y: % 1.8e\n",temp_y);
    	    //printf("temp_tol_x:  % 1.8e\n", temp_tol_x);
			//printf("temp_tol_y:  % 1.8e\n", temp_tol_y);

			/*if (i == repeat - 1)
			{
				printf("llrNURBS_volume\n");
				printf("temp_tol_x:  % 1.8e\n", temp_tol_x);
				printf("temp_tol_y:  % 1.8e\n", temp_tol_y);
				printf("temp_tol_z:  % 1.8e\n", temp_tol_z);
			}*/

			//収束した場合////////////////////////////////////////////////////////////////
    	    //if (temp_tol_x + temp_tol_y + temp_tol_z < tol && i != 1*repeat/4 && i != 2*repeat/4 && i != 3*repeat/4) //収束条件は考慮する必要あり
			if (temp_tol_x < tol && temp_tol_y < tol && temp_tol_z < tol && i != 0 /*&& i != 1*repeat/4 && i != 2*repeat/4 && i != 3*repeat/4*/)
			{
				(*output_xi)   = temp_xi;
				(*output_eta)  = temp_eta;
				(*output_zeta) = temp_zeta;

				//printf("j = %d\tconvergence:llrNURBS_volume\n", j);
				//printf("temp_tol_x:  % 1.8e\n", temp_tol_x);
				//printf("temp_tol_y:  % 1.8e\n", temp_tol_y);
				//printf("temp_tol_z:  % 1.8e\n", temp_tol_z);

				return i;
			}

			InverseMatrix(temp_matrix);

			temp_dxi   = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y) + temp_matrix[0][2] * (pz - temp_z); //ここ違うかも，確認
			temp_deta  = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y) + temp_matrix[1][2] * (pz - temp_z);
			temp_dzeta = temp_matrix[2][0] * (px - temp_x) + temp_matrix[2][1] * (py - temp_y) + temp_matrix[2][2] * (pz - temp_z);

			temp_xi   = temp_xi + temp_dxi;
			temp_eta  = temp_eta + temp_deta;
			temp_zeta = temp_zeta + temp_dzeta;

			if (temp_xi < knot_vec_xi[0])
				temp_xi = knot_vec_xi[0];
			if (temp_xi > knot_vec_xi[cntl_p_n_xi + order_xi])
				temp_xi = knot_vec_xi[cntl_p_n_xi + order_xi];

			if (temp_eta < knot_vec_eta[0])
				temp_eta = knot_vec_eta[0];
			if (temp_eta > knot_vec_eta[cntl_p_n_eta + order_eta])
				temp_eta = knot_vec_eta[cntl_p_n_eta + order_eta];

			if (temp_zeta < knot_vec_zeta[0])
				temp_zeta = knot_vec_zeta[0];
			if (temp_zeta > knot_vec_zeta[cntl_p_n_zeta + order_zeta])
				temp_zeta = knot_vec_zeta[cntl_p_n_zeta + order_zeta];

    	    //printf("r_xi:  % 1.8e\n", temp_xi);
			//printf("r_eta: % 1.8e\n", temp_eta);
    	    //printf("i=%d\n",i);

			//double temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
			//printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
		}
	}

	/////////////////////////////////////////////////lrlNURBS_volume/////////////////////////////////////////////////
	for (j = 0; j < repeat2; j++)
	{
		//初期値の設定
		if (j == 0)
		{
			temp_xi = knot_vec_xi[0] + knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= 0.5;

			temp_eta = knot_vec_eta[0] + knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= 0.5;

			temp_zeta = knot_vec_zeta[0] + knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= 0.5;
		}
		else
		{
			coef = j / repeat2;

			temp_xi = knot_vec_xi[0] + knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= coef;

			temp_eta = knot_vec_eta[0] + knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= coef;

			temp_zeta = knot_vec_zeta[0] + knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= coef;
		}
		//printf("j = %d\txi = %lf\tcoef = %lf\n", j, temp_xi, coef);

		for (i = 0; i < repeat; i++)
		{
			lrlNURBS_volume(Position_Knots[0][0], Position_Knots[0][1], Position_Knots[0][2], Control_Coord[0], Control_Coord[1], Control_Coord[2], No_Control_point[0][0], No_Control_point[0][1], No_Control_point[0][2] 	, Control_Weight, Order[0][0], Order[0][1], Order[0][2], temp_xi, temp_eta, temp_zeta, &temp_x, &temp_y, &temp_z, &temp_matrix[0][0], &temp_matrix[0][1], &temp_matrix[0][2], &temp_matrix[1][0], &temp_matrix	[1][1], &temp_matrix[1][2], &temp_matrix[2][0], &temp_matrix[2][1], &temp_matrix[2][2]);

			temp_tol_x = px - temp_x;
			temp_tol_x *= temp_tol_x;

			temp_tol_y = py - temp_y;
			temp_tol_y *= temp_tol_y;

			temp_tol_z = pz - temp_z;
			temp_tol_z *= temp_tol_z;

    	    //printf("xi_0:  % 1.8e\n", temp_xi);
			//printf("eta_0: % 1.8e\n", temp_eta);
    	    //printf("px: % 1.8e\n",px);
    	    //printf("temp_x: % 1.8e\n",temp_x);
    	    //printf("py: % 1.8e\n",py);
    	    //printf("temp_y: % 1.8e\n",temp_y);
    	    //printf("temp_tol_x:  % 1.8e\n", temp_tol_x);
			//printf("temp_tol_y:  % 1.8e\n", temp_tol_y);

			/*if (i == repeat - 1)
			{
				printf("lrlNURBS_volume\n");
				printf("temp_tol_x:  % 1.8e\n", temp_tol_x);
				printf("temp_tol_y:  % 1.8e\n", temp_tol_y);
				printf("temp_tol_z:  % 1.8e\n", temp_tol_z);
			}*/

			//収束した場合////////////////////////////////////////////////////////////////
    	    //if (temp_tol_x + temp_tol_y + temp_tol_z < tol && i != 1*repeat/4 && i != 2*repeat/4 && i != 3*repeat/4) //収束条件は考慮する必要あり
			if (temp_tol_x < tol && temp_tol_y < tol && temp_tol_z < tol && i != 0 /*&& i != 1*repeat/4 && i != 2*repeat/4 && i != 3*repeat/4*/)
			{
				(*output_xi)   = temp_xi;
				(*output_eta)  = temp_eta;
				(*output_zeta) = temp_zeta;

				//printf("j = %d\tconvergence:lrlNURBS_volume\n", j);
				//printf("temp_tol_x:  % 1.8e\n", temp_tol_x);
				//printf("temp_tol_y:  % 1.8e\n", temp_tol_y);
				//printf("temp_tol_z:  % 1.8e\n", temp_tol_z);

				return i;
			}

			InverseMatrix(temp_matrix);

			temp_dxi   = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y) + temp_matrix[0][2] * (pz - temp_z); //ここ違うかも，確認
			temp_deta  = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y) + temp_matrix[1][2] * (pz - temp_z);
			temp_dzeta = temp_matrix[2][0] * (px - temp_x) + temp_matrix[2][1] * (py - temp_y) + temp_matrix[2][2] * (pz - temp_z);

			temp_xi   = temp_xi + temp_dxi;
			temp_eta  = temp_eta + temp_deta;
			temp_zeta = temp_zeta + temp_dzeta;

			if (temp_xi < knot_vec_xi[0])
				temp_xi = knot_vec_xi[0];
			if (temp_xi > knot_vec_xi[cntl_p_n_xi + order_xi])
				temp_xi = knot_vec_xi[cntl_p_n_xi + order_xi];

			if (temp_eta < knot_vec_eta[0])
				temp_eta = knot_vec_eta[0];
			if (temp_eta > knot_vec_eta[cntl_p_n_eta + order_eta])
				temp_eta = knot_vec_eta[cntl_p_n_eta + order_eta];

			if (temp_zeta < knot_vec_zeta[0])
				temp_zeta = knot_vec_zeta[0];
			if (temp_zeta > knot_vec_zeta[cntl_p_n_zeta + order_zeta])
				temp_zeta = knot_vec_zeta[cntl_p_n_zeta + order_zeta];

    	    //printf("r_xi:  % 1.8e\n", temp_xi);
			//printf("r_eta: % 1.8e\n", temp_eta);
    	    //printf("i=%d\n",i);

			//double temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
			//printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
		}
	}

	/////////////////////////////////////////////////rllNURBS_volume/////////////////////////////////////////////////
	for (j = 0; j < repeat2; j++)
	{
		//初期値の設定
		if (j == 0)
		{
			temp_xi = knot_vec_xi[0] + knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= 0.5;

			temp_eta = knot_vec_eta[0] + knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= 0.5;

			temp_zeta = knot_vec_zeta[0] + knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= 0.5;
		}
		else
		{
			coef = j / repeat2;

			temp_xi = knot_vec_xi[0] + knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= coef;

			temp_eta = knot_vec_eta[0] + knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= coef;

			temp_zeta = knot_vec_zeta[0] + knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= coef;
		}
		//printf("j = %d\txi = %lf\tcoef = %lf\n", j, temp_xi, coef);

		for (i = 0; i < repeat; i++)
		{
			rllNURBS_volume(Position_Knots[0][0], Position_Knots[0][1], Position_Knots[0][2], Control_Coord[0], Control_Coord[1], Control_Coord[2], No_Control_point[0][0], No_Control_point[0][1], No_Control_point[0][2] 	, Control_Weight, Order[0][0], Order[0][1], Order[0][2], temp_xi, temp_eta, temp_zeta, &temp_x, &temp_y, &temp_z, &temp_matrix[0][0], &temp_matrix[0][1], &temp_matrix[0][2], &temp_matrix[1][0], &temp_matrix	[1][1], &temp_matrix[1][2], &temp_matrix[2][0], &temp_matrix[2][1], &temp_matrix[2][2]);

			temp_tol_x = px - temp_x;
			temp_tol_x *= temp_tol_x;

			temp_tol_y = py - temp_y;
			temp_tol_y *= temp_tol_y;

			temp_tol_z = pz - temp_z;
			temp_tol_z *= temp_tol_z;

    	    //printf("xi_0:  % 1.8e\n", temp_xi);
			//printf("eta_0: % 1.8e\n", temp_eta);
    	    //printf("px: % 1.8e\n",px);
    	    //printf("temp_x: % 1.8e\n",temp_x);
    	    //printf("py: % 1.8e\n",py);
    	    //printf("temp_y: % 1.8e\n",temp_y);
    	    //printf("temp_tol_x:  % 1.8e\n", temp_tol_x);
			//printf("temp_tol_y:  % 1.8e\n", temp_tol_y);

			/*if (i == repeat - 1)
			{
				printf("rllNURBS_volume\n");
				printf("temp_tol_x:  % 1.8e\n", temp_tol_x);
				printf("temp_tol_y:  % 1.8e\n", temp_tol_y);
				printf("temp_tol_z:  % 1.8e\n", temp_tol_z);
			}*/

			//収束した場合////////////////////////////////////////////////////////////////
    	    //if (temp_tol_x + temp_tol_y + temp_tol_z < tol && i != 1*repeat/4 && i != 2*repeat/4 && i != 3*repeat/4) //収束条件は考慮する必要あり
			if (temp_tol_x < tol && temp_tol_y < tol && temp_tol_z < tol && i != 0 /*&& i != 1*repeat/4 && i != 2*repeat/4 && i != 3*repeat/4*/)
			{
				(*output_xi)   = temp_xi;
				(*output_eta)  = temp_eta;
				(*output_zeta) = temp_zeta;

				//printf("j = %d\tconvergence:rllNURBS_volume\n", j);
				//printf("temp_tol_x:  % 1.8e\n", temp_tol_x);
				//printf("temp_tol_y:  % 1.8e\n", temp_tol_y);
				//printf("temp_tol_z:  % 1.8e\n", temp_tol_z);

				return i;
			}

			InverseMatrix(temp_matrix);

			temp_dxi   = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y) + temp_matrix[0][2] * (pz - temp_z); //ここ違うかも，確認
			temp_deta  = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y) + temp_matrix[1][2] * (pz - temp_z);
			temp_dzeta = temp_matrix[2][0] * (px - temp_x) + temp_matrix[2][1] * (py - temp_y) + temp_matrix[2][2] * (pz - temp_z);

			temp_xi   = temp_xi + temp_dxi;
			temp_eta  = temp_eta + temp_deta;
			temp_zeta = temp_zeta + temp_dzeta;

			if (temp_xi < knot_vec_xi[0])
				temp_xi = knot_vec_xi[0];
			if (temp_xi > knot_vec_xi[cntl_p_n_xi + order_xi])
				temp_xi = knot_vec_xi[cntl_p_n_xi + order_xi];

			if (temp_eta < knot_vec_eta[0])
				temp_eta = knot_vec_eta[0];
			if (temp_eta > knot_vec_eta[cntl_p_n_eta + order_eta])
				temp_eta = knot_vec_eta[cntl_p_n_eta + order_eta];

			if (temp_zeta < knot_vec_zeta[0])
				temp_zeta = knot_vec_zeta[0];
			if (temp_zeta > knot_vec_zeta[cntl_p_n_zeta + order_zeta])
				temp_zeta = knot_vec_zeta[cntl_p_n_zeta + order_zeta];

    	    //printf("r_xi:  % 1.8e\n", temp_xi);
			//printf("r_eta: % 1.8e\n", temp_eta);
    	    //printf("i=%d\n",i);

			//double temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
			//printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
		}
	}

	/////////////////////////////////////////////////lllNURBS_volume/////////////////////////////////////////////////
	for (j = 0; j < repeat2; j++)
	{
		//初期値の設定
		if (j == 0)
		{
			temp_xi = knot_vec_xi[0] + knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= 0.5;

			temp_eta = knot_vec_eta[0] + knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= 0.5;

			temp_zeta = knot_vec_zeta[0] + knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= 0.5;
		}
		else
		{
			coef = j / repeat2;

			temp_xi = knot_vec_xi[0] + knot_vec_xi[cntl_p_n_xi + order_xi];
			temp_xi *= coef;

			temp_eta = knot_vec_eta[0] + knot_vec_eta[cntl_p_n_eta + order_eta];
			temp_eta *= coef;

			temp_zeta = knot_vec_zeta[0] + knot_vec_zeta[cntl_p_n_zeta + order_zeta];
			temp_zeta *= coef;
		}
		//printf("j = %d\txi = %lf\tcoef = %lf\n", j, temp_xi, coef);

		for (i = 0; i < repeat; i++)
		{
			lllNURBS_volume(Position_Knots[0][0], Position_Knots[0][1], Position_Knots[0][2], Control_Coord[0], Control_Coord[1], Control_Coord[2], No_Control_point[0][0], No_Control_point[0][1], No_Control_point[0][2] 	, Control_Weight, Order[0][0], Order[0][1], Order[0][2], temp_xi, temp_eta, temp_zeta, &temp_x, &temp_y, &temp_z, &temp_matrix[0][0], &temp_matrix[0][1], &temp_matrix[0][2], &temp_matrix[1][0], &temp_matrix	[1][1], &temp_matrix[1][2], &temp_matrix[2][0], &temp_matrix[2][1], &temp_matrix[2][2]);

			temp_tol_x = px - temp_x;
			temp_tol_x *= temp_tol_x;

			temp_tol_y = py - temp_y;
			temp_tol_y *= temp_tol_y;

			temp_tol_z = pz - temp_z;
			temp_tol_z *= temp_tol_z;

    	    //printf("xi_0:  % 1.8e\n", temp_xi);
			//printf("eta_0: % 1.8e\n", temp_eta);
    	    //printf("px: % 1.8e\n",px);
    	    //printf("temp_x: % 1.8e\n",temp_x);
    	    //printf("py: % 1.8e\n",py);
    	    //printf("temp_y: % 1.8e\n",temp_y);
    	    //printf("temp_tol_x:  % 1.8e\n", temp_tol_x);
			//printf("temp_tol_y:  % 1.8e\n", temp_tol_y);

			/*if (i == repeat - 1)
			{
				printf("lllNURBS_volume\n");
				printf("temp_tol_x:  % 1.8e\n", temp_tol_x);
				printf("temp_tol_y:  % 1.8e\n", temp_tol_y);
				printf("temp_tol_z:  % 1.8e\n", temp_tol_z);
			}*/

			//収束した場合////////////////////////////////////////////////////////////////
    	    //if (temp_tol_x + temp_tol_y + temp_tol_z < tol && i != 1*repeat/4 && i != 2*repeat/4 && i != 3*repeat/4) //収束条件は考慮する必要あり
			if (temp_tol_x < tol && temp_tol_y < tol && temp_tol_z < tol && i != 0 /*&& i != 1*repeat/4 && i != 2*repeat/4 && i != 3*repeat/4*/)
			{
				(*output_xi)   = temp_xi;
				(*output_eta)  = temp_eta;
				(*output_zeta) = temp_zeta;

				//printf("j = %d\tconvergence:lllNURBS_volume\n", j);
				//printf("temp_tol_x:  % 1.8e\n", temp_tol_x);
				//printf("temp_tol_y:  % 1.8e\n", temp_tol_y);
				//printf("temp_tol_z:  % 1.8e\n", temp_tol_z);

				return i;
			}

			InverseMatrix(temp_matrix);

			temp_dxi   = temp_matrix[0][0] * (px - temp_x) + temp_matrix[0][1] * (py - temp_y) + temp_matrix[0][2] * (pz - temp_z); //ここ違うかも，確認
			temp_deta  = temp_matrix[1][0] * (px - temp_x) + temp_matrix[1][1] * (py - temp_y) + temp_matrix[1][2] * (pz - temp_z);
			temp_dzeta = temp_matrix[2][0] * (px - temp_x) + temp_matrix[2][1] * (py - temp_y) + temp_matrix[2][2] * (pz - temp_z);

			temp_xi   = temp_xi + temp_dxi;
			temp_eta  = temp_eta + temp_deta;
			temp_zeta = temp_zeta + temp_dzeta;

			if (temp_xi < knot_vec_xi[0])
				temp_xi = knot_vec_xi[0];
			if (temp_xi > knot_vec_xi[cntl_p_n_xi + order_xi])
				temp_xi = knot_vec_xi[cntl_p_n_xi + order_xi];

			if (temp_eta < knot_vec_eta[0])
				temp_eta = knot_vec_eta[0];
			if (temp_eta > knot_vec_eta[cntl_p_n_eta + order_eta])
				temp_eta = knot_vec_eta[cntl_p_n_eta + order_eta];

			if (temp_zeta < knot_vec_zeta[0])
				temp_zeta = knot_vec_zeta[0];
			if (temp_zeta > knot_vec_zeta[cntl_p_n_zeta + order_zeta])
				temp_zeta = knot_vec_zeta[cntl_p_n_zeta + order_zeta];

    	    //printf("r_xi:  % 1.8e\n", temp_xi);
			//printf("r_eta: % 1.8e\n", temp_eta);
    	    //printf("i=%d\n",i);

			//double temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
			//printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
		}
	}
	return 0;
}

int Jacobian(int El_No, double a[DIMENSION][DIMENSION], double Local_coord[DIMENSION], double X[No_Control_point_ON_ELEMENT[Element_patch[El_No]]][DIMENSION], int Total_Control_Point)
{
	int i, j, k;
	//printf("El_No_jacobi:%d\n",El_No);
	for (i = 0; i < DIMENSION; i++)
	{
		//printf("最終a[%d]",i);
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
			//printf("[%d]:%10.9lf\t",j,a[i][j]);
		}
		//printf("\n");
	}
	/*for (i = 0; i < DIMENSION; i++) {
		for (j = 0; j < DIMENSION; j++) {
			printf("a[%d][%d]:%le\n",i,j,a[i][j]);
		}
	}*/
	return 0;
}

//Bマトリックスを求める関数
int Make_B_Matrix(int El_No, double B[D_MATRIX_SIZE][KIEL_SIZE], double Local_coord[DIMENSION], double X[No_Control_point_ON_ELEMENT[Element_patch[El_No]]][DIMENSION], double *J, int Total_Control_Point)
{
	double a[DIMENSION][DIMENSION], b[DIMENSION][No_Control_point_ON_ELEMENT[Element_patch[El_No]]];

	int i, j, k;

	Jacobian(El_No, a, Local_coord, X, Total_Control_Point);

	//printf("Local_coord[0]%lf\t[1]%lf\t[2]%lf\n",Local_coord[0],Local_coord[1],Local_coord[2]);

	*J = InverseMatrix(a);
	//if (El_No == 21)
	//{
	//	for (i = 0; i < DIMENSION; i++)
	//	{
	//		printf("要素境界");
	//		for (k = 0; k < DIMENSION; k++)
	//		{
	//			printf("a[%d][%d]=%lf\t",i,k,a[i][k]);
	//		}
	//		printf("\n");
	//	}
	//}


	//printf("B_Matri_J:%lf\n",*J);
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
		B[0][3 * i]     = b[0][i];
		B[0][3 * i + 1] = 0.0;
		B[0][3 * i + 2] = 0.0;

		B[1][3 * i]     = 0.0;
		B[1][3 * i + 1] = b[1][i];
		B[1][3 * i + 2] = 0.0;

		B[2][3 * i]     = 0.0;
		B[2][3 * i + 1] = 0.0;
		B[2][3 * i + 2] = b[2][i];

		B[3][3 * i]     = b[1][i];
		B[3][3 * i + 1] = b[0][i];
		B[3][3 * i + 2] = 0.0;

		B[4][3 * i]     = 0.0;
		B[4][3 * i + 1] = b[2][i];
		B[4][3 * i + 2] = b[1][i];

		B[5][3 * i]     = b[2][i];
		B[5][3 * i + 1] = 0.0;
		B[5][3 * i + 2] = b[0][i];
	}

	/*for (i = 0; i < D_MATRIX_SIZE; i++)for (j = 0; j < KIEL_SIZE; j++) {
		printf("B[%d][%d]=%le\n",i ,j ,B[i][j]);
	}*/
	return 0;
}
//応力歪マトリックス
int Make_D_Matrix(double D[D_MATRIX_SIZE][D_MATRIX_SIZE], double E, double nu/*, int DM*/)
{
	int i, j;

	double E_ii = ( 1.0 -  nu ) / ( ( 1.0 + nu ) * ( 1.0 - 2.0 * nu) ) * E;
	double E_ij = nu  / ( ( 1.0 + nu ) * ( 1.0 - 2.0 * nu) ) * E;
	double G    = E / ( 2.0 * ( 1.0 + nu ) );

	double D1[D_MATRIX_SIZE][D_MATRIX_SIZE] =  {{E_ii , E_ij , E_ij , 0.0 , 0.0 , 0.0},
												{E_ij , E_ii , E_ij , 0.0 , 0.0 , 0.0},
												{E_ij , E_ij , E_ii , 0.0 , 0.0 , 0.0},
												{ 0.0 ,  0.0 ,  0.0 ,   G , 0.0 , 0.0},
												{ 0.0 ,  0.0 ,  0.0 , 0.0 ,   G , 0.0},
												{ 0.0 ,  0.0 ,  0.0 , 0.0 , 0.0 ,   G}};

	for (i = 0; i < D_MATRIX_SIZE; i++)
	{
		for (j = 0; j < D_MATRIX_SIZE; j++)
		{
			D[i][j] = D1[i][j];
		}
	}

	/*for (i = 0; i < D_MATRIX_SIZE; i++)
	{
		for (j = 0; j < D_MATRIX_SIZE; j++)
		{
			printf("D[%d][%d]=%le\n", i, j, D[i][j]);
		}
	}*/
	return 0;
}

//ガウスの数値積分法の中身
int BDBJ(double B[D_MATRIX_SIZE][KIEL_SIZE], double D[D_MATRIX_SIZE][D_MATRIX_SIZE], double J, double K_EL[KIEL_SIZE][KIEL_SIZE])
{
	int i, j, k;
	double BD[KIEL_SIZE][D_MATRIX_SIZE];

	//[B]T[D][B]Jの計算
	for (i = 0; i < KIEL_SIZE; i++)
	{
		for (j = 0; j < D_MATRIX_SIZE; j++)
		{
			BD[i][j] = 0.0;
			for (k = 0; k < D_MATRIX_SIZE; k++)
			{
				//printf("B[%d][%d]=%le D[%d][%d]=%le\n", k,i,B[k][i],k,j,D[k][j] );
				//printf("B[%d][%d]=%le\n", k,i,B[k][i]);
				BD[i][j] += B[k][i] * D[k][j];
				//printf("BD[%d][%d]=%e\n",i,j,BD[i][j] ); //あとでコメントアウト
			}
		}
	}
	//for( j = 0; j < D_MATRIX_SIZE; j++ )for( i = 0; i < KIEL_SIZE; i++ )printf("B[%d][%d]=%le\n",j,i,B[j][i] );
	for (i = 0; i < KIEL_SIZE; i++)
	{
		for (j = 0; j < KIEL_SIZE; j++)
		{
			K_EL[i][j] = 0.0;
			for (k = 0; k < D_MATRIX_SIZE; k++)
			{
				K_EL[i][j] += BD[i][k] * B[k][j];
			}
			K_EL[i][j] *= J;
		}
	}
	return 0;
}

int coupled_BDBJ(double B[D_MATRIX_SIZE][KIEL_SIZE], double D[D_MATRIX_SIZE][D_MATRIX_SIZE], double BG[D_MATRIX_SIZE][KIEL_SIZE], double J, double K_EL[KIEL_SIZE][KIEL_SIZE])
{
	int i, j, k;
	double BD[KIEL_SIZE][D_MATRIX_SIZE];

	//[B]T[D][B]Jの計算
	for (i = 0; i < KIEL_SIZE; i++)
	{
		for (j = 0; j < D_MATRIX_SIZE; j++)
		{
			BD[i][j] = 0.0;
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
			K_EL[i][j] = 0.0;
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
int Make_K_EL(int El_No, double X[No_Control_point_ON_ELEMENT[Element_patch[El_No]]][DIMENSION], double K_EL[KIEL_SIZE][KIEL_SIZE], double E, double nu, int Total_Element, int Total_Control_Point)
{

	int i, j, k, l;

	double K1[KIEL_SIZE][KIEL_SIZE], B[D_MATRIX_SIZE][KIEL_SIZE];
	double D[D_MATRIX_SIZE][D_MATRIX_SIZE];
	double J = 0.0;

	Check_Volume = 0.0;
	J_test = 0.0;

	for (i = 0; i < KIEL_SIZE; i++)
	{
		for (j = 0; j < KIEL_SIZE; j++)
		{
			K_EL[i][j] = 0.0;
		}
	}
	Make_D_Matrix(D, E, nu);

	for (i = 0; i < GaussPt_3D; i++)
	{
		//printf("i=%d\tGxi%lf\n",i,Gxi[i][0]);
		//printf("Make_K_EL\t\tGxi[%d][0]:%lf\t[1]:%lf[2]:%lf\tw:%lf\n", i, Gxi[i][0],Gxi[i][1],Gxi[i][2],w[i]);
		Make_B_Matrix(El_No, B, Gxi[i], X, &J, Total_Control_Point);

		BDBJ(B, D, J, K1);
		J_test += J;
		Check_Volume += J * w[i];
		for (k = 0; k < KIEL_SIZE; k++)
		{
			for (l = 0; l < KIEL_SIZE; l++)
			{
				K_EL[k][l] += w[i] * K1[k][l];
			}
		} //printf("w[%d]=%f\n",i,w[i]);
	}
	//printf("El_No:%d J_test = %e  Check_Volume = %.13e\n", El_No, J_test, Check_Volume);
	//printf("G=%f\n",G );
	//for ( k = 0; k < KIEL_SIZE; k++) {
	//	for ( l = 0; l < KIEL_SIZE; l++) {
	//		printf("K_EL[%d][%d]:%le\n",k,l,K_EL[k][l]);
	//	}
	//}//あとでコメントアウト

	return 0;
}

//結合要素剛性マトリックス
int Make_coupled_K_EL(int El_No_loc, int El_No_glo, double X[No_Control_point_ON_ELEMENT[Element_patch[El_No_loc]]][DIMENSION], double XG[No_Control_point_ON_ELEMENT[Element_patch[El_No_glo]]][DIMENSION], double K_EL[KIEL_SIZE][KIEL_SIZE], double E, double nu)
{
	int i, ii, j, jj, k, l;
	int BDBJ_flag;

	double K1[KIEL_SIZE][KIEL_SIZE];
	double B[D_MATRIX_SIZE][KIEL_SIZE], BG[D_MATRIX_SIZE][KIEL_SIZE];
	double D[D_MATRIX_SIZE][D_MATRIX_SIZE];
	double J = 0.0;
	double J_test = 0.0;
	double G_Gxi[MAX_POW_Ng][DIMENSION];	//グローバルパッチ上での親要素内座標xi_bar,eta_bar



	Total_BDBJ_flag = 0;

	//printf("El_No_loc;%d\t#######El_No_glo;%d#######\n",El_No_loc, El_No_glo);
	for (i = 0; i < KIEL_SIZE; i++)
	{
		for (j = 0; j < KIEL_SIZE; j++)
		{
			K_EL[i][j] = 0.0;
		}
	}

	Make_D_Matrix(D, E, nu);

	for (i = 0; i < GaussPt_3D; i++)	//ガウス点のループ(local) //ここは要素によって変える必要がある
	{
		//printf("-----------------------------------------------------Gauss_point_No:%d-----------------------------------------------------\n",i);
		////ローカルガウス点がグローバル要素に含まれているかの判定
		//ローカル要素ガウス点の物理座標算出
		double data_result_shape[3] = {0.0};
		double output_xi, output_eta, output_zeta;
		int patch_n;
		//printf("Gxi[%d][0]:%lf\t[1]:%lf[2]:%lf\tw:%lf\n", i, Gxi[i][0],Gxi[i][1],Gxi[i][2],w[i]);
		for (j = 0; j < No_Control_point_ON_ELEMENT[Element_patch[El_No_loc]]; j++)
		{
			for (jj = 0; jj < DIMENSION; jj++)
			{
				data_result_shape[jj] += Shape_func(j, Gxi[i], El_No_loc) * X[j][jj];
			}
		}
		//ローカル要素ガウス点のグローバルパッチ上のパラメータ空間座標算出
		for (j = 0; j < Total_Patch_on_mesh[0]; j++)	//グローバルメッシュ[0]上
		{
			ii = Calc_xi_eta_zeta(data_result_shape[0], data_result_shape[1], data_result_shape[2], Position_Knots[j][0], Position_Knots[j][1], Position_Knots[j][2], Control_Coord[0], Control_Coord[1], Control_Coord[2], No_Control_point[j][0], No_Control_point[j][1], No_Control_point[j][2], Control_Weight, Order[j][0], Order[j][1], Order[j][2], &output_xi, &output_eta, &output_zeta);

			if (ii != 1)
			{
				printf("Gauss_point_No: %d\tGlobal patch: %d\tNewton_iteration: %d\n", i, j, ii);
				printf("\tx, y, z      : % 1.8e\t% 1.8e\t% 1.8e\n", data_result_shape[0], data_result_shape[1], data_result_shape[2]);
				printf("\txi, eta, zeta: % 1.8e\t% 1.8e\t% 1.8e\n", output_xi, output_eta, output_zeta);

				if (ii == 0)
				{
					printf("-ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR-\n");
				}
			}
			patch_n = j;
		}

		for (j = 0; j < No_knot[patch_n][0]; j++)
		{
			if(fabs(output_xi - Position_Knots[patch_n][0][j]) < ERR)
			{
				output_xi = Position_Knots[patch_n][0][j];
			}
		}

		for (j = 0; j < No_knot[patch_n][1]; j++)
		{
			if(fabs(output_eta - Position_Knots[patch_n][1][j]) < ERR)
			{
				output_eta = Position_Knots[patch_n][1][j];
			}
		}

		for (j = 0; j < No_knot[patch_n][2]; j++)
		{
			if(fabs(output_zeta - Position_Knots[patch_n][2][j]) < ERR)
			{
				output_zeta = Position_Knots[patch_n][2][j];
			}
		}

		//if(i == 3 || i == 4 || i == 5 || i == 12 || i == 13 || i == 14 || i == 21 || i == 22 || i == 23)
		//{
		//	printf("output_eta = %20.16lf\tPosition_Knots[0][1][%d] = %20.16lf\t[%d] = %20.16lf\n",output_eta, Order[patch_n][1] + ENC[patch_n][El_No_glo][1] , Position_Knots[patch_n][1][Order[patch_n][1] + ENC[patch_n][El_No_glo][1]] ,Order[patch_n][1] + ENC[patch_n][El_No_glo][1] + 1,  Position_Knots[patch_n][1][Order[patch_n][1] + ENC[patch_n][El_No_glo][1] + 1]);
		//}
		//要素内外判定
		if (   Position_Knots[patch_n][0][Order[patch_n][0] + ENC[patch_n][El_No_glo][0]] <= output_xi   && output_xi   < Position_Knots[patch_n][0][Order[patch_n][0] + ENC[patch_n][El_No_glo][0] + 1]
			&& Position_Knots[patch_n][1][Order[patch_n][1] + ENC[patch_n][El_No_glo][1]] <= output_eta  && output_eta  < Position_Knots[patch_n][1][Order[patch_n][1] + ENC[patch_n][El_No_glo][1] + 1]
			&& Position_Knots[patch_n][2][Order[patch_n][2] + ENC[patch_n][El_No_glo][2]] <= output_zeta && output_zeta < Position_Knots[patch_n][2][Order[patch_n][2] + ENC[patch_n][El_No_glo][2] + 1])	//要素内であるとき
		{//要チェック
			BDBJ_flag = 1;
			//if (Position_Knots[patch_n][1][Order[patch_n][1] + ENC[patch_n][El_No_glo][1]] == output_eta)
			//{
			//	printf("要素境界\n");
			//}
			//printf("BDBJ_flag\n");

			//親要素座標の算出
			G_Gxi[i][0] = - 1.0 + 2.0 * (output_xi   - Position_Knots[patch_n][0][Order[patch_n][0] + ENC[patch_n][El_No_glo][0]]) / (Position_Knots[patch_n][0][Order[patch_n][0] + ENC[patch_n][El_No_glo][0] + 1] - Position_Knots[patch_n][0][Order[patch_n][0] + ENC[patch_n][El_No_glo][0]]);
			G_Gxi[i][1] = - 1.0 + 2.0 * (output_eta  - Position_Knots[patch_n][1][Order[patch_n][1] + ENC[patch_n][El_No_glo][1]]) / (Position_Knots[patch_n][1][Order[patch_n][1] + ENC[patch_n][El_No_glo][1] + 1] - Position_Knots[patch_n][1][Order[patch_n][1] + ENC[patch_n][El_No_glo][1]]);
			G_Gxi[i][2] = - 1.0 + 2.0 * (output_zeta - Position_Knots[patch_n][2][Order[patch_n][2] + ENC[patch_n][El_No_glo][2]]) / (Position_Knots[patch_n][2][Order[patch_n][2] + ENC[patch_n][El_No_glo][2] + 1] - Position_Knots[patch_n][2][Order[patch_n][2] + ENC[patch_n][El_No_glo][2]]);
			//printf("G_Gxi[%d][1]=%20.16lf\n",i , G_Gxi[i][1]);
			if (fabs(1.0 - G_Gxi[i][0]) < ERR)
			{
				G_Gxi[i][0] = 1.0;
			}
			if (fabs(-1.0 - G_Gxi[i][0]) < ERR)
			{
				G_Gxi[i][0] = -1.0;
			}

			if (fabs(1.0 - G_Gxi[i][1]) < ERR)
			{
				G_Gxi[i][1] = 1.0;
			}
			if (fabs(-1.0 - G_Gxi[i][1]) < ERR)
			{
				G_Gxi[i][1] = -1.0;
			}

			if (fabs(1.0 - G_Gxi[i][2]) < ERR)
			{
				G_Gxi[i][2] = 1.0;
			}
			if (fabs(-1.0 - G_Gxi[i][2]) < ERR)
			{
				G_Gxi[i][2] = -1.0;
			}
		}
		else	//要素外であるとき
		{
			BDBJ_flag = 0;
		}
		//printf("i=%d\tBDBJ_flag=%d\n", i, BDBJ_flag);
		////結合要素剛性マトリックス計算
		//要素内であるとき、次を計算
		if (BDBJ_flag == 1)
		{
			Total_BDBJ_flag++;
			Same_BDBJ_flag[i]++;
			//if (Position_Knots[patch_n][1][Order[patch_n][1] + ENC[patch_n][El_No_glo][1]] == output_eta)
			//{
			//	printf("要素境界");
				//printf("G_Gxi[%d][0]=%lf\t[%d][1]=%lf\t[%d][2]=%lf\n",i,G_Gxi[i][0],i,G_Gxi[i][1],i,G_Gxi[i][2]);
			//}


			//for (m = 0; m < 27; m++)
			//{
			//	if (Position_Knots[patch_n][1][Order[patch_n][1] + ENC[patch_n][El_No_glo][1]] == output_eta)
			//	{
			//		printf("要素境界");
			//		printf("XG[%d][0]=%lf\t[%d][1]=%lf\t[%d][2]=%lf\n",m,XG[m][0],m,XG[m][1],m,XG[m][2]);
			//	}
			//}

			//重なるグローバル要素のBマトリックス
			Make_B_Matrix(El_No_glo, BG, G_Gxi[i], XG, &J, Total_Control_Point_to_mesh[No_Input_File]);
			//for (m = 0; m < 3; m++)
			//{
			//	printf("i=%d\tBG[%d]",i,m);
			//	for (mm = 0; mm < 81; mm++)
			//	{
			//		printf("[%d]=%15.14lf\t",mm,BG[m][mm]);
			//	}
			//	printf("\n");
			//}

			//ローカル要素のBマトリックス
			//printf("Gxi[%d][0]:%lf\t[1]:%lf[2]:%lf\tw:%lf\n", i, Gxi[i][0],Gxi[i][1],Gxi[i][2],w[i]);
			Make_B_Matrix(El_No_loc, B, Gxi[i], X, &J, Total_Control_Point_to_mesh[No_Input_File]);
			//for (m = 0; m < 3; m++)
			//{
			//	printf("i=%d\tBL[%d]",i,m);
			//	for (mm = 0; mm < 81; mm++)
			//	{
			//		printf("[%d]=%15.14lf\t",mm,B[m][mm]);
			//	}
			//	printf("\n");
			//}
			//BGTDBLの計算
			coupled_BDBJ(B, D, BG, J, K1);
			J_test += J;
			for (k = 0; k < KIEL_SIZE; k++)
			{
				for (l = 0; l < KIEL_SIZE; l++)
				{
					K_EL[k][l] += w[i] * K1[k][l];
					//printf("K_EL[%d][%d]=%lf\tw[%d]=%lf\tK_1[%d][%d]=%lf\n",k,l,K_EL[k][l],i,w[i],k,l,K1[k][l]);
				}
			} //printf("w[%d]=%f\n",i,w[i]);
		}

		if (i == GaussPt_3D - 1 )
		{
			printf("-------------------Total_BDBJ_flag=%d-------------------\n", Total_BDBJ_flag);
		}
	}
	return 0;
}


//for s-IGA
//Newton Raphsonによって出力されたxi,etaから重なる要素を求める
int ele_check(int patch_n, double para_coord[DIMENSION], double phys_coord[DIMENSION], int mesh_n, int element_n_over)
{
	int i, ii, iii;
	int j;
	int k, kk;
	int l;
    int RangeCheck_frag;						//要素を求め終えたら立てるフラグ
    int temp_ad[DIMENSION][MAX_N_ORDER + 1];	//要素の位置を求めるための値
	int No_line[DIMENSION];						//xi,eta,zetaが含まれている要素の列数
	int n = 1;

	for(j = 0; j < DIMENSION; j++)
    {
		//初期化
		l = 0;
		No_line[j] = 0;
		temp_ad[j][MAX_N_ORDER + 1] = 0;
		RangeCheck_frag = 0;

		for(k = 0; k < No_Control_point[patch_n][j] - 1; k++)
		{
            if ( RangeCheck_frag == 1 )
				break;


			//0////////////////////Local要素のガウス点がGlobalパッチ内にない場合////////////////////
			if ( para_coord[j] < Position_Knots[patch_n][j][0] - ERR || para_coord[j] > Position_Knots[patch_n][j][No_knot[patch_n][j] - 1] + ERR)
			{
				printf("no over element\n");
				RangeCheck_frag++;
			}


			//1////////////////////Local要素のガウス点がGlobal要素内部にある場合////////////////////
			if ( para_coord[j] < Position_Knots[patch_n][j][Order[patch_n][j] + k] )
            {
                //printf("if\nPosition_Knots[%d][%d][%d]=%le\n",
                //        patch_n,j,Order[patch_n][j]+k,Position_Knots[patch_n][j][Order[patch_n][j]+k]);
                int kk = 0;
                for (kk = 0; kk < k + 1; kk++)
                {
                    if ( para_coord[j] > Position_Knots[patch_n][j][Order[patch_n][j] + k - kk] )
                    {
						//printf("j=%d\t要素内部\n",j);
                        //printf("ifif\nPosition_Knots[%d][%d][%d]=%le\n", patch_n,j,Order[patch_n][j]+k-kk, Position_Knots[patch_n][j][Order[patch_n][j]+k-kk]);
                        temp_ad[j][l] = k - kk;
						//printf("temp_ad[%d][%d]=%d\n",j,l,temp_ad[j][l]);
						l++;
                        RangeCheck_frag++;
                        break;
                    }
                }
            }


			//2////////////////////Local要素のガウス点がGlobal要素境界上にある場合////////////////////要チェック
			if ( fabs(para_coord[j] - Position_Knots[patch_n][j][Order[patch_n][j] + k]) < ERR )
			{
				//printf("j=%d\t要素境界上",j);
				//2-1////////////////////ガウス点の座標がGlobalパッチの始点上にある場合////////////////////
				if ( fabs(para_coord[j] - Position_Knots[patch_n][j][0]) < ERR )
				{
					//printf("のパッチの始点上\n");
					temp_ad[j][l] = k;
					//printf("start point\n");
					//printf("temp_ad[%d][%d]=%d\n",j,l,temp_ad[j][l]);
					l++;
					break;
				}
				//2-2////////////////////ガウス点の座標がGlobalパッチの終点上にある場合////////////////////
				if ( fabs(para_coord[j] - Position_Knots[patch_n][j][No_knot[patch_n][j] - 1]) < ERR )
				{
					//printf("のパッチの終点上\n");
					temp_ad[j][l] = k - 1;
					//printf("finish point\n");
					//printf("temp_ad[%d][%d]=%d\n",j,l,temp_ad[j][l]);
					l++;
					break;
				}
				//2-3////////////////////ガウス点の座標がGlobal要素境界上にある場合////////////////////
				else
				{
					//printf("の要素間\n");
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
					if( fabs(para_coord[j] - Position_Knots[patch_n][j][Order[patch_n][j] + k + kk + 1]) < ERR )
					//多重ノット（次数分ループ）
					{
						printf("C0 continuity\n");
						temp_ad[j][l] = k + kk;
						//printf("temp_ad[%d][%d]=%d\n",j,l,temp_ad[j][l]);
						l++;
					}
					if( para_coord[j] != Position_Knots[patch_n][j][Order[patch_n][j] + k + kk + 1] )
					break;
				}
				RangeCheck_frag++;
			}
        }

		No_line[j] = l;
        //printf("No_line[%d] = %d\n", j ,No_line[j]);
		n *= l;	//各方向のNo_lineを掛け合わせる
    }

	//for (j = 0; j < DIMENSION; j++)
	//{
	//	printf("No_line[%d];%d\n",j,No_line[j]);
	//}
	for (iii = 0; iii < No_line[2]; iii++)
	{
		for (ii = 0; ii < No_line[1]; ii++)
		{
			for (i = 0; i < No_line[0]; i++)
			{
				temp_element_n[(iii * No_line[0] * No_line[1]) + (ii * No_line[0]) + i]     = temp_ad[0][i] + (temp_ad[1][ii] * line_No_Total_element[patch_n][0]) + (temp_ad[2][iii] * line_No_Total_element[patch_n][0] * line_No_Total_element[patch_n][1]);
				temp_element_n_ref[(iii * No_line[0] * No_line[1]) + (ii * No_line[0]) + i] = temp_ad[0][i] + (temp_ad[1][ii] * line_No_Total_element[patch_n][0]) + (temp_ad[2][iii] * line_No_Total_element[patch_n][0] * line_No_Total_element[patch_n][1]);
				//printf("el[%d](x,y)=(%le,%le) in element[%d] on patch[%d] on mesh[0]\n",element_n_over ,phys_coord[0],phys_coord[1] ,temp_element_n[i*No_line[0]+ii],patch_n);
				//printf("temp_element_n[%d];%d\n", (iii * No_line[0] * No_line[1]) + (ii * No_line[0]) + i, temp_element_n[(iii * No_line[0] * No_line[1]) + (ii * No_line[0]) + i]);
			}
		}
	}

	//printf("n_element_over;%d\n",n);	//ある一点に重なっている要素の総数

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
		if (element_n_point[i] != element_n_point[i - 1])
		{
			NELOVER[element_n][j] = element_n_point[i];
			//printf("NELOVER\n");
			j++;
		}
	}
	//j;要素element_nに重なる要素の総数
	return j;
}






//要素の重なりを求める(要素のガウス点から求める)
void Check_coupled_Glo_Loc_element(double element_loc[DIMENSION], int mesh_n_over, int mesh_n_org)
{
	int re;
    int e;
    int i, j ,k, m;
	int b;
	int l, ll;
	int n_elements_over_point[MAX_N_POINT_OVER]; //あるガウス点に重なっている要素の数
	int patch_n, itr_n;
	int i_gg, i_ee, i_zz;
	int g_n;
	double output_para[DIMENSION];
	int Total_n_elements;
	int Check_coupled_No[MAX_N_ELEMENT_OVER];
	double Percent_Check_coupled_No;
	int MAX_NNLOVER = 0;

	for (i = 0; i < MAX_N_ELEMENT_OVER; i++)
	{
		Check_coupled_No[i] = 0;
	}



	for (m = 0; m < 2; m++) //最初Ng個のガウス点で重なりを求め，NNLOVER[e]>=2のeに対して，再度10個のガウス点で重なりを求める．
	{
		Make_Gauss_points(m);

		for (re = 0; re < real_Total_Element_on_mesh[mesh_n_over]; re++)
		{
			e = real_element[re + real_Total_Element_to_mesh[mesh_n_over]];

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

				for (i_zz = 0; i_zz < GaussPt_1dir; i_zz++)
				{
					for (i_ee = 0; i_ee < GaussPt_1dir; i_ee++)
					{
						for (i_gg = 0; i_gg < GaussPt_1dir; i_gg++)
						{
							double data_result_shape[3] = {0.0};

							g_n = (i_zz * GaussPt_1dir * GaussPt_1dir) + (i_ee * GaussPt_1dir) + i_gg;
							//printf("-----------------------------------------------------Gauss_point_No:%d-----------------------------------------------------\n",g_n);
							element_loc[0] = Gxi[g_n][0];
							element_loc[1] = Gxi[g_n][1];
							element_loc[2] = Gxi[g_n][2];

							for (b = 0; b < No_Control_point_ON_ELEMENT[Element_patch[e]]; b++)
							{
								for (j = 0; j < DIMENSION; j++)
								{
    			                    data_result_shape[j] += Shape_func(b, element_loc, e) * Node_Coordinate[Controlpoint_of_Element[e][b]][j];
								}
								//printf("Element_patch[%d] = %d, Total_Control_Point_on_mesh[%d] = %d, e = %d, Node_Coordinate[%d][0] = %lf\t[%d][1] = %lf\t[%d][2] = %lf\n",e,Element_patch[e], mesh_n_over, Total_Control_Point_on_mesh[mesh_n_over], e,Controlpoint_of_Element[e][b], Node_Coordinate[Controlpoint_of_Element[e][b]][0],Controlpoint_of_Element[e][b], Node_Coordinate[Controlpoint_of_Element[e][b]][1],Controlpoint_of_Element[e][b], Node_Coordinate[Controlpoint_of_Element[e][b]][2]);
							}
							//printf("data_result_shape[0] = %lf\t[1] = %lf\t[2] = %lf\n",data_result_shape[0],data_result_shape[1],data_result_shape[2]);

							for (i = 0; i < Total_Patch_on_mesh[mesh_n_org]; i++)
							{
								int ii = Calc_xi_eta_zeta(data_result_shape[0], data_result_shape[1], data_result_shape[2], Position_Knots[i][0], Position_Knots[i][1], Position_Knots[i][2], Control_Coord[0], Control_Coord[1], Control_Coord[2], No_Control_point[i][0], No_Control_point[i][1], No_Control_point[i][2], Control_Weight, Order[i][0], Order[i][1], Order[i][2], &output_para[0], &output_para[1], &output_para[2]);
								if (ii != 1)
								{
									printf("-----------------------------------------m:%d\tElement_No:%d-----------------------------------------\n",m,e );
									printf("Gauss_point_No: %d\tGlobal patch: %d\tNewton_iteration: %d\n",g_n, i, ii);
									printf("\tx, y, z      : %20.16lf\t%20.16lf\t%20.16lf\n", data_result_shape[0], data_result_shape[1], data_result_shape[2]);
									printf("\txi, eta, zeta: %20.16lf\t%20.16lf\t%20.16lf\n", output_para[0], output_para[1], output_para[2]);

									if (ii == 0)
									{
										printf("-ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR-\n");
									}
								}
								patch_n = i;
								itr_n = ii;
							}

							for (j = 0; j < No_knot[patch_n][0]; j++)
							{
								if(fabs(output_para[0] - Position_Knots[patch_n][0][j]) < ERR)
								{
									output_para[0] = Position_Knots[patch_n][0][j];
									//printf("Check_coupled_element_output_para[0]修正\n");
								}
							}
							for (j = 0; j < No_knot[patch_n][1]; j++)
							{
								if(fabs(output_para[1] - Position_Knots[patch_n][1][j]) < ERR)
								{
									output_para[1] = Position_Knots[patch_n][1][j];
									//printf("Check_coupled_element_output_para[1]修正\n");
								}
							}
							for (j = 0; j < No_knot[patch_n][2]; j++)
							{
								if(fabs(output_para[2] - Position_Knots[patch_n][2][j]) < ERR)
								{
									output_para[2] = Position_Knots[patch_n][2][j];
									//printf("Check_coupled_element_output_para[2]修正\n");
								}
							}

							//Newton Laphsonによって出力されたxi,etaから重なる要素を求める
							n_elements_over_point[k] = ele_check(patch_n, output_para, data_result_shape, mesh_n_org, e);
							//printf("itr_n;%d\n",itr_n);
						    if (itr_n == 0)	//data_result_shapeがグローバルメッシュ上にないとき
							{
								n_elements_over_point[k] = 0;
							}
							//printf("n_elements_over_point[%d];%d\n", k, n_elements_over_point[k]);

							Total_n_elements += n_elements_over_point[k];

							//printf("Total_n_elements;%d\n",Total_n_elements);
							for (l = 0; l < n_elements_over_point[k]; l++)
							{
								element_n_point[ll] = temp_element_n[l];

								//printf("element_n_point[%d] = %d\n", ll, element_n_point[ll]);
								ll++;
							}

							k++;
						}
					}
				}
				
				//printf("Total_n_elements;%d\n",Total_n_elements);

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
				/*printf("m=%d\tNNLOVER[%d] = %d\n", m, e, NNLOVER[e]);

				for (i = 0; i < NNLOVER[e]; i++)
				{
					printf("\tNELOVER[%d][%d] = %d\n", e, i, NELOVER[e][i]); //要素eに重なるi番目の要素番号
				}*/
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
	int i, j, k;
	int jj;
	int e;
	int j_n;
	int count;

	j_n = real_Total_Element_to_mesh[No_Input_File] - real_Total_Element_on_mesh[0]; //ローカルメッシュのリアルな要素数
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
			if (NNLOVER[jj] > 0) //ローカルメッシュ上のjj番目の要素に重なる要素があるのならば
			{
				//printf("jj=%d\n",jj);
				for (k = 0; k < NNLOVER[jj]; k++)
				{
					if (NELOVER[jj][k] == e) //ローカルメッシュ上のjj番目の要素に重なるk番目の要素がグローバルメッシュ上のe番目の要素ならば
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
	/*for (i = 0; i < real_Total_Element_to_mesh[No_Input_File]; i++)
	{
		if (NNLOVER[i] > 0)
		{
			printf("NNLOVER[%d]=%d\n",i,NNLOVER[i]);
			for (j = 0; j < NNLOVER[i]; j++)
			{
				printf("\tNELOVER[%d][%d]=%d\n",i ,j ,NELOVER[i][j]);
			}
		}
	}*/
}


/////////////////////////////////////////////////////////////////歪と応力//////////////////////////////////////////////////
void Make_Strain(double E, double nu, int Total_Element, int El_No, int Total_Control_Point)
{
	static double U[MAX_KIEL_SIZE];
	double B[D_MATRIX_SIZE][KIEL_SIZE], X[No_Control_point_ON_ELEMENT[Element_patch[El_No]]][DIMENSION], J;

	int N, e, i, j;

	Make_Gauss_points(0);

	//printf("\nTotal_Control_Point:%d\n",Total_Control_Point);
	//printf("\nTotal_Element:%d\n",Total_Element);
	for (e = 0; e < Total_Element; e++)
	{
		//printf("\nElementNo.:%d\n",e);
		for (N = 0; N < POW_Ng; N++)
		{
			for (i = 0; i < N_STRAIN; i++)
			{
				Strain[e][N][i] = 0.0;
			}
		}
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
		for (N = 0; N < POW_Ng; N++)
		{
			//printf("N:%d\n",N);
			Make_B_Matrix(e, B, Gxi[N], X, &J, Total_Control_Point);
			for (i = 0; i < D_MATRIX_SIZE; i++)
			{
				for (j = 0; j < KIEL_SIZE; j++)
				{
					Strain[e][N][i] += B[i][j] * U[j];
				}
			}
		}
	}
}

//応力
void Make_Stress(double E, double nu, int Total_Element)
{

	static double D[D_MATRIX_SIZE][D_MATRIX_SIZE];
	int e, i, j, k;
	Make_Gauss_points(0);
	Make_D_Matrix(D, E, nu);

	for (e = 0; e < Total_Element; e++)
	{
		for (k = 0; k < POW_Ng; k++)
		{
			for (i = 0; i < N_STRESS; i++)
			{
				Stress[e][k][i] = 0.0;
			}
		}

		for (k = 0; k < POW_Ng; k++)
		{
			for (i = 0; i < D_MATRIX_SIZE; i++)
			{
				for (j = 0; j < D_MATRIX_SIZE; j++)
				{
					Stress[e][k][i] += D[i][j] * Strain[e][k][j];
				}
			}
		}
	}
}

void Make_ReactionForce(int Total_Element, int Total_Control_Point, int El_No)
{
	int e, i, j, k, l, re;
	double B[D_MATRIX_SIZE][KIEL_SIZE], X[No_Control_point_ON_ELEMENT[Element_patch[El_No]]][DIMENSION], J;

	Make_Gauss_points(0);

	for (i = 0; i < Total_Control_Point * DIMENSION; i++)
	{
		ReactionForce[i] = 0.0;
	}

	for (re = 0; re < real_Total_Element; re++)
	{
		e = real_element[re];
		//printf("%d,No_Control_point_ON_ELEMENT[%d]:%d\n",e,Element_patch[e], No_Control_point_ON_ELEMENT[Element_patch[e]]);
		//Bマトリックスを取得
		for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[e]]; i++)
		{
			for (j = 0; j < DIMENSION; j++)
			{
				X[i][j] = Node_Coordinate[Controlpoint_of_Element[e][i]][j];
			}
		}
		for (k = 0; k < POW_Ng; k++)
		{
			Make_B_Matrix(e, B, Gxi[k], X, &J, Total_Control_Point);
			for (j = 0; j < D_MATRIX_SIZE; j++)
				for (l = 0; l < DIMENSION; l++)
					for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[e]]; i++)
						ReactionForce[Controlpoint_of_Element[e][i] * DIMENSION + l] += B[j][i * DIMENSION + l] * Stress[e][k][j] * w[k] * J;
			//printf("J:%le\n", J);
		}
	}
}


void element_coordinate(int Total_Element, int Total_Control_Point)
{
	int i, j, k, e, l = 0;
	double element_edge[27][DIMENSION] =   {{-1.0, -1.0, -1.0},
											{ 1.0, -1.0, -1.0},
											{ 1.0,  1.0, -1.0},
											{-1.0,  1.0, -1.0},
											{ 0.0, -1.0, -1.0},
											{ 1.0,  0.0, -1.0},
											{ 0.0,  1.0, -1.0},
											{-1.0,  0.0, -1.0},
											{ 0.0,  0.0, -1.0},
											{-1.0, -1.0,  0.0},
											{ 1.0, -1.0,  0.0},
											{ 1.0,  1.0,  0.0},
											{-1.0,  1.0,  0.0},
											{ 0.0, -1.0,  0.0},
											{ 1.0,  0.0,  0.0},
											{ 0.0,  1.0,  0.0},
											{-1.0,  0.0,  0.0},
											{ 0.0,  0.0,  0.0},
											{-1.0, -1.0,  1.0},
											{ 1.0, -1.0,  1.0},
											{ 1.0,  1.0,  1.0},
											{-1.0,  1.0,  1.0},
											{ 0.0, -1.0,  1.0},
											{ 1.0,  0.0,  1.0},
											{ 0.0,  1.0,  1.0},
											{-1.0,  0.0,  1.0},
											{ 0.0,  0.0,  1.0}};
	//double data_result_shape[2]={0.0};

	for (e = 0; e < Total_Element; e++)
	{
		for (k = 0; k < 27; k++)
		{
			double data_result_shape[3] = {0.0};
			for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[e]]; i++)
			{
				for (j = 0; j < DIMENSION; j++)
				{
					data_result_shape[j] += Shape_func(i, element_edge[k], e) * Node_Coordinate[Controlpoint_of_Element[e][i]][j];
				}
			}
			element_coordinate_Nopoint[l][0] = data_result_shape[0];
			element_coordinate_Nopoint[l][1] = data_result_shape[1];
			element_coordinate_Nopoint[l][2] = data_result_shape[2];
			//printf("l;%d\n",l);
			/*printf("element_coordinate_Nopoint[%d][0]:%le element_coordinate_Nopoint[%d][1]:%le\n"
			,l,element_coordinate_Nopoint[l][0],l,element_coordinate_Nopoint[l][1]);*/
			l++;
			//printf("data_result_shape[0]:%le data_result_shape[1]:%le\n", data_result_shape[0],data_result_shape[1]);
		}
	}
	for (l = 0; l < 27 * Total_Element; l++)
	{
		same_point_in_Element[l] = l;
	}

	for (l = 0; l < 27 * Total_Element; l++)
	{
		for (i = l - 1; i >= 0; i--)
		{
			if (element_coordinate_Nopoint[l][0] == element_coordinate_Nopoint[i][0] && element_coordinate_Nopoint[l][1] == element_coordinate_Nopoint[i][1] && element_coordinate_Nopoint[l][2] == element_coordinate_Nopoint[i][2])
			{
				//printf("同じ座標の番号l:%d i:%d\n",l,i);
				same_point_in_Element[l] = i;
			}
		}
	}
	//printf("finish_element_coordinate");
}

void Make_Strain_refine(double E, double nu, int Total_Element , int El_No, int Total_Control_Point, double element[DIMENSION])
{
	static double U[MAX_KIEL_SIZE];
	double B[D_MATRIX_SIZE][KIEL_SIZE];
	double X[No_Control_point_ON_ELEMENT[Element_patch[El_No]]][DIMENSION];
	double J;

	int i, j;

	//グローバル要素について
	if (El_No < real_Total_Element_to_mesh[1])
	{
		for( i = 0; i < N_STRAIN; i ++ )
		{
			Strain_ref[i]     = 0.0;
			Strain_ref_glo[i] = 0.0;
			//Strain_ref_loc[i] = 0.0;
		}

		//Bマトリックスと各要素の変位を取得
		for( i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++ )
		{
			for( j = 0; j < DIMENSION; j++ )
			{
				U[ i * DIMENSION + j ] = Displacement[ Controlpoint_of_Element[El_No][i] * DIMENSION + j ];
				X[i][j] = Node_Coordinate[Controlpoint_of_Element[El_No][i]][j];
			}
		}

		Make_B_Matrix( El_No, B, element, X ,&J , Total_Control_Point);

		for( i = 0; i < D_MATRIX_SIZE; i++ )
		{
			for( j = 0; j < KIEL_SIZE; j++ )
			{
				Strain_ref_glo[i] += B[i][j] * U[j];
			}
			Strain_ref[i] = Strain_ref_glo[i];
		}
	}


	//ローカル要素について
	if (El_No >= real_Total_Element_to_mesh[1])
	{
		for( i = 0; i < N_STRAIN; i ++ )
		{
			Strain_ref[i]     = 0.0;
			//Strain_ref_glo[i] = 0.0;
			Strain_ref_loc[i] = 0.0;
		}

		//Bマトリックスと各要素の変位を取得
		for( i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[El_No]]; i++ )
		{
			for( j = 0; j < DIMENSION; j++ )
			{
				U[ i * DIMENSION + j ] = Displacement[ Controlpoint_of_Element[El_No][i] * DIMENSION + j ];
				X[i][j] = Node_Coordinate[Controlpoint_of_Element[El_No][i]][j];
			}
		}

		Make_B_Matrix( El_No, B, element, X ,&J , Total_Control_Point);

		for( i = 0; i < D_MATRIX_SIZE; i++ )
		{
			for( j = 0; j < KIEL_SIZE; j++ )
			{
				Strain_ref_loc[i] += B[i][j] * U[j];
			}
			Strain_ref[i] = Strain_ref_glo[i] + Strain_ref_loc[i];
		}
	}
}

void Make_Stress_refine( double E, double nu ,int Total_Element, int Total_Control_Point, double element[DIMENSION], int El_No)
{
	double D[D_MATRIX_SIZE][D_MATRIX_SIZE];
	int i, j;

	Make_D_Matrix(D, E, nu);

	//Make_Strain_refine(E, nu, Total_Element, El_No, Total_Control_Point, element);


	for( i = 0; i < N_STRESS + 2; i ++ )
	{
		Stress_ref[i]     = 0.0;
		Stress_ref_glo[i] = 0.0;
		Stress_ref_loc[i] = 0.0;
	}

	for( i = 0; i < N_STRESS; i ++ )
	{
		for( j = 0; j < D_MATRIX_SIZE; j++ )
		{
			Stress_ref[i] += D[i][j] * Strain_ref[j];
			Stress_ref_glo[i] += D[i][j] * Strain_ref_glo[j];
			Stress_ref_loc[i] += D[i][j] * Strain_ref_loc[j];
		}
	}
}

void ourput_graph_glo_loc(FILE *fp, int e, double element_gg, double element_ee, double element_zz, double data_result_shape[DIMENSION], double data_result_disp[DIMENSION])
{
	fp = fopen("NURBS/NURBS_Disp_glo_loc.dat", "a");
	fprintf(fp, "%d\t%20.16lf\t%20.16lf\t%20.16lf\t%20.13e\t%20.13e\t%20.13e\n", e, data_result_shape[0], data_result_shape[1], data_result_shape[2], data_result_disp[0], data_result_disp[1], data_result_disp[2]);
	fclose(fp);
}

void ourput_graph_glo(FILE *fp, int e, double element_gg, double element_ee, double element_zz, double data_result_shape[DIMENSION], double data_result_disp[DIMENSION])
{
	//fp = fopen("NURBS/NURBS_Disp_glo.dat", "a");
	//fprintf(fp, "%d\t%lf\t%lf\t%lf\t%20.16lf\t%20.16lf\t%20.16lf\t%20.16lf\t%20.13e\t%20.13e\t%20.13e\n", e, element_gg, element_ee, element_zz, data_result_shape[0], data_result_shape[1], data_result_shape[2], data_result_disp[0], data_result_disp[1], data_result_disp[2]);
	//fclose(fp);
}

void ourput_graph_loc(FILE *fp, int e, double element_gg, double element_ee, double element_zz, double data_result_shape[DIMENSION], double data_result_disp[DIMENSION])
{
	//fp = fopen("NURBS/NURBS_Disp_loc.dat", "a");
	//fprintf(fp, "%d\t%lf\t%lf\t%lf\t%20.16lf\t%20.16lf\t%20.16lf\t%20.16lf\t%20.13e\t%20.13e\t%20.13e\n", e, element_gg, element_ee, element_zz, data_result_shape[0], data_result_shape[1], data_result_shape[2],  data_result_disp[0], data_result_disp[1], data_result_disp[2]);
	//fclose(fp);
}


void output_Strain_refine_glo_loc(int e, double element_gg, double element_ee, double element_zz, double data_result_shape[DIMENSION])
{
	fp = fopen("NURBS/NURBS_Strain_glo_loc.dat","a");
	fprintf(fp, "%d\t%20.16lf\t%20.16lf\t%20.16lf\t%20.13e\t%20.13e\t%20.13e\t%20.13e\t%20.13e\t%20.13e\n", e, data_result_shape[0], data_result_shape[1], data_result_shape[2], Strain_ref[0], Strain_ref[1], Strain_ref[2], Strain_ref[3], Strain_ref[4], Strain_ref[5]);
	fclose(fp);
}

void output_Strain_refine_glo(int e, double element_gg, double element_ee, double element_zz, double data_result_shape[DIMENSION])
{
	//fp = fopen("NURBS/NURBS_Strain_glo.dat","a");
	//fprintf(fp, "%d\t%lf\t%lf\t%lf\t%20.16lf\t%20.16lf\t%20.16lf\t%20.16lf\t%20.13e\t%20.13e\t%20.13e\t%20.13e\t%20.13e\t%20.13e\n", e, element_gg, element_ee, element_zz, data_result_shape[0], data_result_shape[1], data_result_shape[2], Strain_ref_glo[0], Strain_ref_glo[1], Strain_ref_glo[2], Strain_ref_glo[3], Strain_ref_glo[4], Strain_ref_glo[5]);
	//fclose(fp);
}

void output_Strain_refine_loc(int e, double element_gg, double element_ee, double element_zz, double data_result_shape[DIMENSION])
{
	//fp = fopen("NURBS/NURBS_Strain_loc.dat","a");
	//fprintf(fp, "%d\t%lf\t%lf\t%lf\t%20.16lf\t%20.16lf\t%20.16lf\t%20.16lf\t%20.13e\t%20.13e\t%20.13e\t%20.13e\t%20.13e\t%20.13e\n", e, element_gg, element_ee, element_zz, data_result_shape[0], data_result_shape[1], data_result_shape[2], Strain_ref_loc[0], Strain_ref_loc[1], Strain_ref_loc[2], Strain_ref_loc[3], Strain_ref_loc[4], Strain_ref_loc[5]);
	//fclose(fp);
}


void output_Stress_refine_glo_loc(int e, double element_gg, double element_ee, double element_zz, double data_result_shape[DIMENSION])
{
	fp = fopen("NURBS/NURBS_Stress_glo_loc.dat","a");
	fprintf(fp, "%d\t%20.13lf\t%20.13lf\t%20.13lf\t%20.13e\t%20.13e\t%20.13e\n", e, data_result_shape[0], data_result_shape[1], data_result_shape[2], Stress_ref[0], Stress_ref[1], Stress_ref[2]);
	fclose(fp);
}

void output_Stress_refine_glo(int e, double element_gg, double element_ee, double element_zz, double data_result_shape[DIMENSION])
{
	//fp = fopen("NURBS/NURBS_Stress_glo.dat","a");
	//fprintf(fp, "%d\t%lf\t%lf\t%lf\t%20.13lf\t%20.13lf\t%20.13lf\t%20.13lf\t%20.13e\t%20.13e\t%20.13e\t%20.13e\t%20.13e\t%20.13e\n", e, element_gg, element_ee, element_zz, data_result_shape[0], data_result_shape[1], data_result_shape[2], Stress_ref_glo[0], Stress_ref_glo[1], Stress_ref_glo[2], Stress_ref_glo[3], Stress_ref_glo[4], Stress_ref_glo[5]);
	//fclose(fp);
}

void output_Stress_refine_loc(int e, double element_gg, double element_ee, double element_zz, double data_result_shape[DIMENSION])
{
	//fp = fopen("NURBS/NURBS_Stress_loc.dat","a");
	//fprintf(fp, "%d\t%lf\t%lf\t%lf\t%20.13lf\t%20.13lf\t%20.13lf\t%20.13lf\t%20.13e\t%20.13e\t%20.13e\t%20.13e\t%20.13e\t%20.13e\n", e, element_gg, element_ee, element_zz, data_result_shape[0], data_result_shape[1], data_result_shape[2], Stress_ref_loc[0], Stress_ref_loc[1], Stress_ref_loc[2], Stress_ref_loc[3], Stress_ref_loc[4], Stress_ref_loc[5]);
	//fclose(fp);
}

void calculate_Controlpoint_using_NURBS(double element[DIMENSION], double element_glo[DIMENSION], int Total_Element, int Total_Control_Point, double E, double nu)
{
	int e, b, i, j, k, re;
	int i_gg;
	double output_para_ref[DIMENSION];
	int n_elements_over_point_ref[MAX_N_REFINEMENT];  //ある点に重なっている要素の数
	int patch_n_ref;
	int itr_n_ref;


	double data_result_shape[DIMENSION];
	double data_result_disp[DIMENSION];
	double data_result_disp_loc[DIMENSION];
	double data_result_disp_glo[DIMENSION];


	printf("-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
	fp = fopen("NURBS/NURBS_Disp_glo_loc.dat", "a");
	fprintf(fp, "NO_element\txi\teta\tzeta\tx\ty\tz\tdegree\tDisp_x\tDisp_y\tDisp_z\tglo_loc\n");
	fclose(fp);

	fp = fopen("NURBS/NURBS_Disp_glo.dat", "a");
	fprintf(fp, "NO_element\txi\teta\tzeta\tx\ty\tz\tdegree\tDisp_x\tDisp_y\tDisp_z\tglo\n");
	fclose(fp);

	fp = fopen("NURBS/NURBS_Disp_loc.dat", "a");
	fprintf(fp, "NO_element\txi\teta\tzeta\tx\ty\tz\tdegree\tDisp_x\tDisp_y\tDisp_z\tloc\n");
	fclose(fp);


	fp = fopen("NURBS/NURBS_Strain_glo_loc.dat","a");
	fprintf(fp, "NO_element\txi\teta\tzeta\tx\ty\tz\tdegree\tStrain_xx\tStrain_yy\tStrain_zz\tStrain_xy\tStrain_yz\tStrain_zx\tglo_loc\n");
	fclose(fp);

	fp = fopen("NURBS/NURBS_Strain_glo.dat","a");
	fprintf(fp, "NO_element\txi\teta\tzeta\tx\ty\tz\tdegree\tStrain_xx\tStrain_yy\tStrain_zz\tStrain_xy\tStrain_yz\tStrain_zx\tglo\n");
	fclose(fp);

	fp = fopen("NURBS/NURBS_Strain_loc.dat","a");
	fprintf(fp, "NO_element\txi\teta\tzeta\tx\ty\tz\tdegree\tStrain_xx\tStrain_yy\tStrain_zz\tStrain_xy\tStrain_yz\tStrain_zx\tloc\n");
	fclose(fp);


	fp = fopen("NURBS/NURBS_Stress_glo_loc.dat","a");
	fprintf(fp, "NO_element\txi\teta\tzeta\tx\ty\tz\tStress_xx\tStress_yy\tStress_zz\tglo_loc\n");
	fclose(fp);

	fp = fopen("NURBS/0_NURBS_Stress_glo_loc_Crack_Edge.dat","a");
	fprintf(fp, "NO_element\txi\teta\tzeta\tx\ty\tz\tx_CE\ty_CE\tD_f_C_E\tdegree\tSlope\tStress_xx\tStress_yy\tStress_zz\tK\tFe\tglo_loc\n");
	fclose(fp);

	fp = fopen("NURBS/1_NURBS_Stress_glo_loc_x_dir.dat","a");
	fprintf(fp, "NO_element\txi\teta\tzeta\tx\ty\tz\tD_f_C_E\tD_f_C_E/c\tStress_xx\tStress_yy\tStress_zz\tK\tFe\tglo_loc\n");
	fclose(fp);

	fp = fopen("NURBS/2_NURBS_Stress_glo_loc_y_dir.dat","a");
	fprintf(fp, "NO_element\txi\teta\tzeta\tx\ty\tz\tD_f_C_E\tD_f_C_E/a\tStress_xx\tStress_yy\tStress_zz\tK\tFe\tglo_loc\n");
	fclose(fp);

	fp = fopen("NURBS/NURBS_Stress_glo.dat","a");
	fprintf(fp, "NO_element\txi\teta\tzeta\tx\ty\tz\tdegree\tStress_xx\tStress_yy\tStress_zz\tStress_xy\tStress_yz\tStress_zx\tglo\n");
	fclose(fp);

	fp = fopen("NURBS/NURBS_Stress_loc.dat","a");
	fprintf(fp, "NO_element\txi\teta\tzeta\tx\ty\tz\tdegree\tStress_xx\tStress_yy\tStress_zz\tStress_xy\tStress_yz\tStress_zx\tloc\n");
	fclose(fp);

	Make_Gauss_points(0);

	printf("No_Input_File: %d\n",No_Input_File);
	//ローカルパッチのRFINEMENT
	if (No_Input_File == 2)
	{
		printf("real_Total_Element_to_mesh[0]%d\n",real_Total_Element_to_mesh[0]);
		printf("real_Total_Element_to_mesh[1]%d\n",real_Total_Element_to_mesh[1]);
		printf("real_Total_Element_to_mesh[2]%d\n",real_Total_Element_to_mesh[2]);
		printf("real_Total_Element_to_mesh[3]%d\n",real_Total_Element_to_mesh[3]);
		for (re = real_Total_Element_to_mesh[1]; re < real_Total_Element_to_mesh[No_Input_File - 1]; re++)
		{
			k = 0;

			e = real_element[re];
			printf("eeeeee%d\n",e);

			for (i_gg = 0; i_gg < GaussPt_3D; i_gg++)
			{
				//printf("GaussPt_3D%d\n",GaussPt_3D);
				for (j = 0; j < DIMENSION; j++)
				{
					data_result_shape[j]     = 0.0;
					data_result_disp[j]      = 0.0;
					data_result_disp_loc[j]  = 0.0;
					data_result_disp_glo[j]  = 0.0;
				}

				for (b = 0; b < No_Control_point_ON_ELEMENT[Element_patch[e]]; b++)
				{
					//printf("b%d\n",b);
					for (j = 0; j < DIMENSION; j++)
					{
						//printf("%d %d %le\n",Controlpoint_of_Element[e][b],j,Displacement	[Controlpoint_of_Element[e][b]*DIMENSION+j]);
						//printf("%d %d %20.13le\n",Controlpoint_of_Element[e][b],j,Displacement	[Controlpoint_of_Element[e][b]*DIMENSION+j]);
						data_result_disp_loc[j]  += Shape_func(b, Gxi[i_gg], e) * Displacement[Controlpoint_of_Element[e][b] * DIMENSION + j];
						data_result_shape[j]     += Shape_func(b, Gxi[i_gg], e) * Node_Coordinate[Controlpoint_of_Element[e][b]][j];
					}
				}

				for (i = 0; i < Total_Patch_on_mesh[0]; i++)
				{
					int ii = Calc_xi_eta_zeta(data_result_shape[0], data_result_shape[1], data_result_shape[2], Position_Knots[i][0], Position_Knots[i][1], Position_Knots[i][2], Control_Coord[0], Control_Coord[1], Control_Coord[2], No_Control_point[i][0], No_Control_point[i][1], No_Control_point[i][2], Control_Weight, Order[i][0], Order[i][1], Order[i][2], &output_para_ref[0], &output_para_ref[1], &output_para_ref[2]);

					if (ii != 1)
					{
						printf("Division_No: %d\tGlobal patch: %d\tNewton_iteration: %d\n",i_gg, i, ii);
						printf("\tx, y, z      : % 1.8e\t% 1.8e\t% 1.8e\n", data_result_shape[0], data_result_shape[1], data_result_shape[2]);
						printf("\txi, eta, zeta: % 1.8e\t% 1.8e\t% 1.8e\n", output_para_ref[0], output_para_ref[1], output_para_ref[2]);

						if (ii == 0)
						{
							printf("-ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR--ERROR-\n");
						}
					}
					patch_n_ref = i;
					itr_n_ref = ii;
				}

				//Newton Laphsonによって出力されたxi,etaから重なる要素を求める
				n_elements_over_point_ref[k] = ele_check(patch_n_ref, output_para_ref, data_result_shape, patch_n_ref, e);
				//printf("itr_n;%d\n",itr_n);
				if (itr_n_ref == 0)	//data_result_shapeがグローバルメッシュ上にないとき
				{
					n_elements_over_point_ref[k] = 0;
				}

				//printf("n_elements_over_point_ref[%d];%d\n", k, n_elements_over_point_ref[k]);

				//for (i = 0; i < n_elements_over_point_ref[k]; i++)
				//{
				//	printf("temp_element_n_ref[%d]=%d\n",i , temp_element_n_ref[i]);
				//}

				element_n_point_for_ref = temp_element_n_ref[0];
				//printf("element_n_point_for_ref;%d\n",element_n_point_for_ref);

				k++;

				//ノットベクトルからパラメーター空間座標を求める．
				element_glo[0] = - 1.0 + 2.0 * (output_para_ref[0] - Position_Knots[patch_n_ref][0][Order[patch_n_ref][0] + ENC[patch_n_ref][element_n_point_for_ref][0]]) / (Position_Knots[patch_n_ref][0][Order[patch_n_ref][0] + ENC[patch_n_ref][element_n_point_for_ref][0] + 1] - Position_Knots[patch_n_ref][0][Order[patch_n_ref][0] + ENC[patch_n_ref][element_n_point_for_ref][0]]);
				element_glo[1] = - 1.0 + 2.0 * (output_para_ref[1] - Position_Knots[patch_n_ref][1][Order[patch_n_ref][1] + ENC[patch_n_ref][element_n_point_for_ref][1]]) / (Position_Knots[patch_n_ref][1][Order[patch_n_ref][1] + ENC[patch_n_ref][element_n_point_for_ref][1] + 1] - Position_Knots[patch_n_ref][1][Order[patch_n_ref][1] + ENC[patch_n_ref][element_n_point_for_ref][1]]);
				element_glo[2] = - 1.0 + 2.0 * (output_para_ref[2] - Position_Knots[patch_n_ref][2][Order[patch_n_ref][2] + ENC[patch_n_ref][element_n_point_for_ref][2]]) / (Position_Knots[patch_n_ref][2][Order[patch_n_ref][2] + ENC[patch_n_ref][element_n_point_for_ref][2] + 1] - Position_Knots[patch_n_ref][2][Order[patch_n_ref][2] + ENC[patch_n_ref][element_n_point_for_ref][2]]);

				for (b = 0; b < No_Control_point_ON_ELEMENT[Element_patch[element_n_point_for_ref]]; b++)
				{
					for (j = 0; j < DIMENSION; j++)
					{
						data_result_disp_glo[j]  += Shape_func(b, element_glo, element_n_point_for_ref) * Displacement[Controlpoint_of_Element[element_n_point_for_ref][b] * DIMENSION + j];
					}
					//printf("element_n_point_for_ref;%d\tlement_glo[0];%lf\t[1];%lf\t[2];%lf\n",element_n_point_for_ref,element_glo[0],element_glo[1],element_glo[2]);
				}
	
				data_result_disp[0] = data_result_disp_glo[0] + data_result_disp_loc[0];
				data_result_disp[1] = data_result_disp_glo[1] + data_result_disp_loc[1];
				data_result_disp[2] = data_result_disp_glo[2] + data_result_disp_loc[2];
		
				ourput_graph_glo_loc(fp, e,  Gxi[i_gg][0], Gxi[i_gg][1], Gxi[i_gg][2], data_result_shape, data_result_disp);
				ourput_graph_loc(fp, e,  Gxi[i_gg][0], Gxi[i_gg][1], Gxi[i_gg][2], data_result_shape, data_result_disp_loc);
				ourput_graph_glo(fp, element_n_point_for_ref, element_glo[0], element_glo[1], element_glo[2], data_result_shape, data_result_disp_glo);

				//Make_Strain_refine(E, nu, Total_Element, element_n_point_for_ref, Total_Control_Point, element_glo);
				Make_Strain_refine(E, nu, Total_Element, e, 					  Total_Control_Point, Gxi[i_gg]);
			
				output_Strain_refine_glo_loc(e, Gxi[i_gg][0], Gxi[i_gg][1], Gxi[i_gg][2], data_result_shape);
				output_Strain_refine_loc(e, Gxi[i_gg][0], Gxi[i_gg][1], Gxi[i_gg][2], data_result_shape);
				output_Strain_refine_glo(element_n_point_for_ref, element_glo[0], element_glo[1], element_glo[2], data_result_shape);
					
				Make_Stress_refine(E, nu, Total_Element, Total_Control_Point, Gxi[i_gg], e);
				output_Stress_refine_glo_loc(e, Gxi[i_gg][0], Gxi[i_gg][1], Gxi[i_gg][2], data_result_shape);
				output_Stress_refine_loc(e, Gxi[i_gg][0], Gxi[i_gg][1], Gxi[i_gg][2], data_result_shape);
				output_Stress_refine_glo(element_n_point_for_ref, element_glo[0], element_glo[1], element_glo[2], data_result_shape);
	
				k = 0;
			}
		}
	}








	//IGAのRFINEMENT
	if (No_Input_File == 1)
	{
		for (re = 0; re < real_Total_Element_to_mesh[1]; re++)
		{
			e = real_element[re];

			for (i_gg = 0; i_gg < GaussPt_3D; i_gg++)
			{
				for (j = 0; j < DIMENSION; j++)
				{
					data_result_shape[j] = 0.0;
					data_result_disp_glo[j]  = 0.0;
				}

				for (b = 0; b < No_Control_point_ON_ELEMENT[Element_patch[e]]; b++)
				{
					for (j = 0; j < DIMENSION; j++)
					{
						//printf("%d %d %le\n",Controlpoint_of_Element[e][b],j,Displacement	[Controlpoint_of_Element[e][b]*DIMENSION+j]);
						//printf("%d %d %20.13le\n",Controlpoint_of_Element[e][b],j,Displacement	[Controlpoint_of_Element[e][b]*DIMENSION+j]);
						data_result_disp_glo[j]  += Shape_func(b, Gxi[i_gg], e) * Displacement[Controlpoint_of_Element[e][b] * DIMENSION + j];
						data_result_shape[j]     += Shape_func(b, Gxi[i_gg], e) * Node_Coordinate[Controlpoint_of_Element[e][b]][j];
					}
				}

				data_result_disp[0] = data_result_disp_glo[0];
				data_result_disp[1] = data_result_disp_glo[1];
				data_result_disp[2] = data_result_disp_glo[2];

				ourput_graph_glo_loc(fp, e, Gxi[i_gg][0], Gxi[i_gg][1], Gxi[i_gg][2], data_result_shape, data_result_disp);
				ourput_graph_glo(fp, e, Gxi[i_gg][0], Gxi[i_gg][1], Gxi[i_gg][2], data_result_shape, data_result_disp_glo);

				Make_Strain_refine(E, nu, Total_Element, e, Total_Control_Point, Gxi[i_gg]);
				output_Strain_refine_glo_loc(e, Gxi[i_gg][0], Gxi[i_gg][1], Gxi[i_gg][2], data_result_shape);
				output_Strain_refine_glo(e, Gxi[i_gg][0], Gxi[i_gg][1], Gxi[i_gg][2], data_result_shape);

				Make_Stress_refine(E, nu, Total_Element, Total_Control_Point, Gxi[i_gg], e);
				output_Stress_refine_glo_loc(e, Gxi[i_gg][0], Gxi[i_gg][1], Gxi[i_gg][2], data_result_shape);
				output_Stress_refine_glo(e, Gxi[i_gg][0], Gxi[i_gg][1], Gxi[i_gg][2], data_result_shape);
			}
		}
	}
	//for(i = 0; i< p; i++)
	//	printf("%20.13e\t%20.13e\t%20.13e\t%20.13e\t%20.13e\t%20.13e\n", data_result_shape_x[i], data_result_shape_y[i], data_result_shape_z[i], data_result_disp_x[i], data_result_disp_y[i], data_result_disp_z[i]);
	//printf("-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
}

void Gausspoint_coordinate(int Total_Element, int Total_Control_Point)
{
	int i, j, k, e;

	Make_Gauss_points(0);

	for (e = 0; e < Total_Element; e++)
	{
		for (k = 0; k < POW_Ng; k++)
		{
			double data_result_shape[3] = {0.0};
			for (i = 0; i < No_Control_point_ON_ELEMENT[Element_patch[e]]; i++)
			{
				for (j = 0; j < DIMENSION; j++)
				{
					//printf("i:%d Gxi[%d][%d]:%le\n", i,k,j,Gxi[k][i]);
					data_result_shape[j] += Shape_func(i, Gxi[k], e) * Node_Coordinate[Controlpoint_of_Element[e][i]][j];//ここは3で良い
				}
			}
			Gausspoint_coordinates[e][k][0] = data_result_shape[0];
			Gausspoint_coordinates[e][k][1] = data_result_shape[1];
			Gausspoint_coordinates[e][k][2] = data_result_shape[2];
		}
	}
	//printf("finish_Gausspoint_coordinate");
}

//void calculate_extendmesh_using_NURBS(double element_emsh[DIMENSION], int Total_Element, int Total_Control_Point)
//{
//	int e, b, j, re;
//	int p = Total_Control_Point;
//	//for(e=0; e < Total_Element; e++){
//
//	for (re = 0; re < real_Total_Element; re++)
//	{
//		e = real_element[re];
//		double element_gg = 0.0, element_ee = 0.0, element_zz = 0.0, element_delta_gg, element_delta_ee, element_delta_zz;
//
//		int i_gg, i_ee, i_zz, element_ndiv_gg = 1, element_ndiv_ee = 1, element_ndiv_zz = 1;
//
//		No_points_for_new_zarusoba = (element_ndiv_gg + 1) * (element_ndiv_ee + 1) * (element_ndiv_zz + 1) * real_Total_Element;
//
//		element_delta_gg = 2.0 / element_ndiv_gg;
//		element_delta_ee = 2.0 / element_ndiv_ee;
//		element_delta_zz = 2.0 / element_ndiv_zz;
//
//		for (i_zz = 0; i_zz < element_ndiv_zz + 1; i_zz++)
//		{
//			for (i_ee = 0; i_ee < element_ndiv_ee + 1; i_ee++)
//			{
//				for (i_gg = 0; i_gg < element_ndiv_gg + 1; i_gg++)
//				{
//					double data_result_shape[3] = {0.0};
//					double data_result_disp[3] = {0.0};
//
//					element_gg = -1.0 + element_delta_gg * i_gg;
//					element_ee = -1.0 + element_delta_ee * i_ee;
//					element_zz = -1.0 + element_delta_zz * i_zz;
//					element_emsh[0] = element_gg;
//					element_emsh[1] = element_ee;
//					element_emsh[2] = element_zz;
//
//					// printf("element_gg:%le\n",element_gg);
//					//printf("element_ee:%le\n",element_ee);
//
//					for (b = 0; b < No_Control_point_ON_ELEMENT[Element_patch[e]]; b++)
//					{
//						for (j = 0; j < DIMENSION; j++)
//						{
//							//printf("%d %d %le\n",Controlpoint_of_Element[e][b],j,Displacement	[Controlpoint_of_Element[e][b]*DIMENSION+j]);
//							//printf("%d %d %20.13le\n",Controlpoint_of_Element[e][b],j,Displacement	[Controlpoint_of_Element[e][b]*DIMENSION+j]);
//							data_result_disp[j]  += Shape_func(b, Total_Control_Point, element_emsh, e) * Displacement[Controlpoint_of_Element[e][b] * DIMENSION + j];
//							data_result_shape[j] += Shape_func(b, Total_Control_Point, element_emsh, e) * Node_Coordinate[Controlpoint_of_Element[e][b]][j];
//						}
//					}
//
//					//data_result_shape_x[p]=data_result_shape[0];
//					//data_result_shape_y[p]=data_result_shape[1];
//					//p++;
//
//					//printf("\n");
//
//					data_result_shape_x_for_new_zarusoba[p] = data_result_shape[0];
//					data_result_shape_y_for_new_zarusoba[p] = data_result_shape[1];
//					data_result_shape_z_for_new_zarusoba[p] = data_result_shape[2];
//
//					data_result_disp_x_for_new_zarusoba[p] = data_result_disp[0];
//					data_result_disp_y_for_new_zarusoba[p] = data_result_disp[1];
//					data_result_disp_z_for_new_zarusoba[p] = data_result_disp[2];
//					p++;
//
//					//NURBS_points(fp,No_points_for_colored_points e,data_result_shape[0],	data_result_shape[1])
//				}
//			}
//		}
//	}
//	//for (j = 0; j < No_points_for_new_zarusoba; j++)
//	//{
//	//printf("%d: %lf\t%lf\t%lf\n", j, data_result_shape_x_for_new_zarusoba[j] ,data_result_shape_y_for_new_zarusoba[j], data_result_shape_z_for_new_zarusoba[j]);
//	//}
//}



int SerchForElement(int mesh_n, int iPatch, int Total_Element, int iX, int iY, int iZ)
{
	int iii;

	for (iii = 0; iii < Total_Element; iii++)
	{
		if (Element_patch[iii + Total_Element_to_mesh[mesh_n]] == iPatch)
		{
			//printf("Check SerchForElement 1 iii = %d\n", iii);
			//printf("ENC[iii][0] = %d ENC[iii][1] = %d ENC[iii][2] = %d  iX = %d  iY = %d  iZ = %d\n", ENC[iii][0], ENC[iii][1], ENC[iii][2], iX, iY, iZ);
			if (iX == ENC[iPatch][iii + Total_Element_to_mesh[mesh_n]][0] && iY == ENC[iPatch][iii + Total_Element_to_mesh[mesh_n]][1] && iZ == ENC[iPatch][iii + Total_Element_to_mesh[mesh_n]][2])
				goto loopend;
			/* iii --; */

			//printf("Check SerchForElement 2 iii = %d\n", iii);
		}
	}
loopend:

	return (iii);
}


void Setting_Dist_Load(int mesh_n, int Total_Control_Point, int iPatch, int Total_Element, int iCoord,int jCoord, double val_Coord, double iRange_Coord[2], double jRange_Coord[2], int type_load, double Coeff_Dist_Load_i[3], double Coeff_Dist_Load_j[3])
{
	int iii, jjj, kkk;
	int N_Seg_Load_Element_iDir = 0, N_Seg_Load_Element_jDir = 0, kCoord;
	int iPos[2] = {-10000, -10000}, jPos[2] = {-10000, -10000}, kPos[2] = {-10000, -10000};
	int No_Element_for_Integration[MAX_N_KNOT][MAX_N_KNOT], No_Element_For_Dist_Load_idir, No_Element_For_Dist_Load_jdir;
	int iX, iY, iZ;
	//int Element_Integration;
	int iControlpoint[MAX_NO_CCpoint_ON_ELEMENT], ic, ig_i, ig_j;
	double val_kCoord_Local;

	Make_Gauss_points(0);

		/* type_load: 0: Dist load in x direction
 * 	              1:              y direction
 * 	              2:              normal direciton */

	/* iCoord=0: Load on Eta=Constant
	   iCoord=1: Load on Xi=Constant */
	if (iCoord == 0 && jCoord == 1)
	{
		kCoord = 2;
	}
	if (iCoord == 1 && jCoord == 2)
	{
		kCoord = 0;
	}
	if (iCoord == 2 && jCoord == 0)
	{
		kCoord = 1;
	}

	/* val_Coord: Value of Eta or Xi of the line or surface to give the distributed load */

	/* Setting elements needed to computed the distributed load */

	for (iii = Order[iPatch][iCoord]; iii < No_knot[iPatch][iCoord] - Order[iPatch][iCoord] - 1; iii++)
	{
		double epsi = 0.00000000001;
		/* iPos[0] = -10000; iPos[1] = -10000; jPos[0] = -10000; jPos[1] = -10000;*/
		//printf("Check1 iii = %d\n", iii);
		//printf("Check2 Position_Knots[iCoord][iii]= %f  Range_Coord[0] =%f Position_Knots[iCoord][iii+1] = %f\n", Position_Knots[iPatch][iCoord][iii], iRange_Coord[0], Position_Knots[iPatch][iCoord][iii + 1]);
		/*

		if(Position_Knots[iCoord][iii]-epsi <= Range_Coord[0] &&
			Position_Knots[iCoord][iii+1]+epsi > Range_Coord[0]) iPos[0] = iii;

		if(Position_Knots[iCoord][iii]-epsi <= Range_Coord[1] &&
                        Position_Knots[iCoord][iii+1]+epsi > Range_Coord[1]) iPos[1] = iii+1;
	*/
		if (Position_Knots[iPatch][iCoord][iii] - epsi <= iRange_Coord[0])
			iPos[0] = iii;
		if (Position_Knots[iPatch][iCoord][iii + 1] - epsi <= iRange_Coord[1])
			iPos[1] = iii + 1;
	}

	if (iPos[0] < 0 || iPos[1] < 0)
	{
		printf("Error (Stop) iPos[0] = %d iPos[1] = %d\n", iPos[0], iPos[1]);
		exit(0);
	}







	for (jjj = Order[iPatch][jCoord]; jjj < No_knot[iPatch][jCoord] - Order[iPatch][jCoord] - 1; jjj++)
	{
		double epsi = 0.00000000001;

		//printf("Check1 jjj = %d\n", jjj);
		//printf("Check2 Position_Knots[iCoord][jjj]= %f  jRange_Coord[0] =%f Position_Knots[iCoord][jjj+1] = %f\n", Position_Knots[iPatch][iCoord][jjj], jRange_Coord[0], Position_Knots[iPatch][jCoord][jjj + 1]);
		/*

		if(Position_Knots[iCoord][iii]-epsi <= Range_Coord[0] &&
			Position_Knots[iCoord][iii+1]+epsi > Range_Coord[0]) iPos[0] = iii;

		if(Position_Knots[iCoord][iii]-epsi <= Range_Coord[1] &&
                        Position_Knots[iCoord][iii+1]+epsi > Range_Coord[1]) iPos[1] = iii+1;
	*/
		if (Position_Knots[iPatch][jCoord][jjj] - epsi <= jRange_Coord[0])
			jPos[0] = jjj;
		if (Position_Knots[iPatch][jCoord][jjj + 1] - epsi <= jRange_Coord[1])
			jPos[1] = jjj + 1;
	}

	if (jPos[0] < 0 || jPos[1] < 0)
	{
		printf("Error (Stop) jPos[0] = %d jPos[1] = %d\n", jPos[0], jPos[1]);
		exit(0);
	}









	for (kkk = Order[iPatch][kCoord]; kkk < No_knot[iPatch][kCoord] - Order[iPatch][kCoord] - 1; kkk++)
	{
		double epsi = 0.00000000001;
		/* kkk=Order[kCoord]; */
		if (Position_Knots[iPatch][kCoord][kkk] - epsi <= val_Coord && Position_Knots[iPatch][kCoord][kkk + 1] + epsi > val_Coord)
		{
			kPos[0] = kkk;
			kPos[1] = kkk + 1;
			val_kCoord_Local = -1.0 + 2.0 * (val_Coord - Position_Knots[iPatch][kCoord][kkk]) / (Position_Knots[iPatch][kCoord][kkk + 1] - Position_Knots[iPatch][kCoord][kkk]);
		}
		//(2019_06_13)printf("Check kkk count: kkk =  %d\n",kkk);
	}

	if (kPos[0] < 0 || kPos[1] < 0)
	{
		printf("Error (Stop) kPos[0] = %d kPos[1] = %d\n", kPos[0], kPos[1]);
		exit(0);
	}


	for (iii = iPos[0]; iii < iPos[1]; iii++)
	{
		N_Seg_Load_Element_iDir++;
	}

	for (jjj = jPos[0]; jjj < jPos[1]; jjj++)
	{
		N_Seg_Load_Element_jDir++;
	}

	//kDir_Element = jPos[0] - Order[iPatch][iCoord];
	iii = 0;
	jjj = 0;




	if (iCoord == 0 && jCoord == 1)
	{
		//int iX, iY, iZ;
		iZ = kPos[0] - Order[iPatch][2];
		//(2019_06_13)printf("Check iPos[0] = %d  iPos[1] = %d\n",iPos[0],iPos[1]);
		for (iX = iPos[0] - Order[iPatch][0]; iX < iPos[1] - Order[iPatch][0]; iX++)
		{
			for (iY = jPos[0] - Order[iPatch][1]; iY < jPos[1] - Order[iPatch][1]; iY++)
			{
				//(2019_06_13)printf("Check iX = %d\n",iX);
				No_Element_for_Integration[iii][jjj] = SerchForElement(mesh_n, iPatch, Total_Element, iX, iY, iZ);
				//printf("Check No_Element_for_Integration[%d][%d] = %d\n", iii, jjj, No_Element_for_Integration[iii][jjj]);
				jjj++;
			}
			iii++;
			jjj = 0;
		}
		for (iY = jPos[0] - Order[iPatch][1]; iY < jPos[1] - Order[iPatch][1]; iY++)
			jjj++;
	}

	if (iCoord == 1 && jCoord == 2)
	{
		//int iX, iY, iZ;
		iX = kPos[0] - Order[iPatch][0];
		//(2019_06_13)printf("Check iPos[0] = %d  iPos[1] = %d\n",iPos[0],iPos[1]);
		for (iY = iPos[0] - Order[iPatch][1]; iY < iPos[1] - Order[iPatch][1]; iY++) //積分する要素の数
		{
			for (iZ = jPos[0] - Order[iPatch][2]; iZ < jPos[1] - Order[iPatch][2]; iZ++) //積分する要素の数
			{
				//(2019_06_13)printf("Check iY = %d\n",iY);
				No_Element_for_Integration[iii][jjj] = SerchForElement(mesh_n, iPatch, Total_Element, iX, iY, iZ);
				//printf("Check No_Element_for_Integration[%d][%d] = %d\n", iii, jjj, No_Element_for_Integration[iii][jjj]);
				jjj++;
			}
			iii++;
			jjj = 0;
		}
		for (iZ = jPos[0] - Order[iPatch][2]; iZ < jPos[1] - Order[iPatch][2]; iZ++)
			jjj++;
	}

	if (iCoord == 2 && jCoord == 0)
	{
		//int iX, iY, iZ;
		iY = kPos[0] - Order[iPatch][1];
		//(2019_06_13)printf("Check iPos[0] = %d  iPos[1] = %d\n",iPos[0],iPos[1]);
		for (iZ = iPos[0] - Order[iPatch][2]; iZ < iPos[1] - Order[iPatch][2]; iZ++) //積分する要素の数
		{
			for (iX = jPos[0] - Order[iPatch][0]; iX < jPos[1] - Order[iPatch][0]; iX++) //積分する要素の数
			{
				//(2019_06_13)printf("Check iY = %d\n",iY);
				No_Element_for_Integration[iii][jjj] = SerchForElement(mesh_n, iPatch, Total_Element, iX, iY, iZ);
				//printf("Check No_Element_for_Integration[%d][%d] = %d\n", iii, jjj, No_Element_for_Integration[iii][jjj]);
				jjj++;
			}
			iii++;
			jjj = 0;
		}
		for (iX = jPos[0] - Order[iPatch][0]; iX < jPos[1] - Order[iPatch][0]; iX++)
			jjj++;
	}


	No_Element_For_Dist_Load_idir = iii;
	No_Element_For_Dist_Load_jdir = jjj;
	//printf("No_Element_For_Dist_Load_idir = %d   No_Element_For_Dist_Load_jdir = %d\n",No_Element_For_Dist_Load_idir ,No_Element_For_Dist_Load_jdir);

	/* Book keeping finished */




	for (iii = 0; iii < No_Element_For_Dist_Load_idir; iii++)
	{ //B
		for (jjj = 0; jjj < No_Element_For_Dist_Load_jdir; jjj++)
		{
			//(2019_06_13)printf("Check3 iii = %d\n",iii);
			//(2019_06_13)printf("Total_element_all_ID[No_Element_for_Integration[iii]] = %d\n No_Element_for_Integration[iii] = %d  iii = %d\n",
			//Total_element_all_ID[No_Element_for_Integration[iii]],No_Element_for_Integration[iii],iii);
			if (Total_element_all_ID[No_Element_for_Integration[iii][jjj]] == 1)
			{ //A
				iX = ENC[iPatch][No_Element_for_Integration[iii][jjj]][0];
				iY = ENC[iPatch][No_Element_for_Integration[iii][jjj]][1];
				iZ = ENC[iPatch][No_Element_for_Integration[iii][jjj]][2];
				//printf("iii = %d jjj = %d   iX = %d  iY = %d  iZ = %d\n",iii, jjj, iX, iY, iZ);

				for (ic = 0; ic < (Order[iPatch][0] + 1) * (Order[iPatch][1] + 1) * (Order[iPatch][2] + 1); ic++)
				{
					iControlpoint[ic] = Controlpoint_of_Element[No_Element_for_Integration[iii][jjj]][ic];
				}

				for (ig_i = 0; ig_i < Ng; ig_i++) //積分スタートここは3で良い
				{
					for (ig_j = 0; ig_j < Ng; ig_j++)
					{
						double Local_Coord[3], sfc, dxyzdgez_i[3], dxyzdgez_j[3], detJ, XiEtaZetaCoordParen_idir, XiEtaZetaCoordParen_jdir, valDistLoad;
						int icc;
						Local_Coord[kCoord] = val_kCoord_Local;
						Local_Coord[iCoord] = G_1d[ig_i];
						Local_Coord[jCoord] = G_1d[ig_j];
						//printf("ig_i = %d  ig_j = %d   Local_Coord[kCoord] = %f  Local_Coord[iCoord] = %f  Local_Coord[jCoord] = %f\n", ig_i, ig_j, Local_Coord[kCoord], Local_Coord[iCoord], Local_Coord[jCoord]);

						ShapeFunc_from_paren(Local_Coord, iCoord, No_Element_for_Integration[iii][jjj]);
						XiEtaZetaCoordParen_idir = Position_Data_param[iCoord];

						ShapeFunc_from_paren(Local_Coord, jCoord, No_Element_for_Integration[iii][jjj]);
						XiEtaZetaCoordParen_jdir = Position_Data_param[jCoord];

						//printf("Check  Coeff_Dist_Load[0] = %f Coeff_Dist_Load[1] = %f  Coeff_Dist_Load[2] = %f  Position_Data_param[iCoord] = %f\n", Coeff_Dist_Load[0], 	Coeff_Dist_Load[1], Coeff_Dist_Load[2], 	Position_Data_param[iCoord]);
						valDistLoad = (Coeff_Dist_Load_i[0] + Coeff_Dist_Load_i[1] * XiEtaZetaCoordParen_idir + Coeff_Dist_Load_i[2] * XiEtaZetaCoordParen_idir * XiEtaZetaCoordParen_idir)
									* (Coeff_Dist_Load_j[0] + Coeff_Dist_Load_j[1] * XiEtaZetaCoordParen_jdir + Coeff_Dist_Load_j[2] * XiEtaZetaCoordParen_jdir * XiEtaZetaCoordParen_jdir);

						//き裂の分布荷重を与える時のために...
        	        	//printf("XiEtaCoordParen=%lf\n",XiEtaCoordParen);
        	        	//double sita;
        	        	//sita = XiEtaCoordParen*PI/2;
						//sita = XiEtaCoordParen*2*PI/line_No_real_element[0][1];
        	        	//printf("sita=%lf\n",sita*180/PI);
        	        	//valDistLoad = cos(sita);
        	        	//printf("valDistLoad=%2.10lf\n",valDistLoad);

						dxyzdgez_i[0] = 0.0;
						dxyzdgez_i[1] = 0.0;
						dxyzdgez_i[2] = 0.0;

						dxyzdgez_j[0] = 0.0;
						dxyzdgez_j[1] = 0.0;
						dxyzdgez_j[2] = 0.0;

						for (icc = 0; icc < (Order[iPatch][0] + 1) * (Order[iPatch][1] + 1) * (Order[iPatch][2] + 1); icc++)
						{
							dxyzdgez_i[0] += dShape_func(icc, iCoord, Local_Coord, No_Element_for_Integration[iii][jjj]) * Node_Coordinate[iControlpoint[icc]][0];
							dxyzdgez_i[1] += dShape_func(icc, iCoord, Local_Coord, No_Element_for_Integration[iii][jjj]) * Node_Coordinate[iControlpoint[icc]][1];
							dxyzdgez_i[2] += dShape_func(icc, iCoord, Local_Coord, No_Element_for_Integration[iii][jjj]) * Node_Coordinate[iControlpoint[icc]][2];

							dxyzdgez_j[0] += dShape_func(icc, jCoord, Local_Coord, No_Element_for_Integration[iii][jjj]) * Node_Coordinate[iControlpoint[icc]][0];
							dxyzdgez_j[1] += dShape_func(icc, jCoord, Local_Coord, No_Element_for_Integration[iii][jjj]) * Node_Coordinate[iControlpoint[icc]][1];
							dxyzdgez_j[2] += dShape_func(icc, jCoord, Local_Coord, No_Element_for_Integration[iii][jjj]) * Node_Coordinate[iControlpoint[icc]][2];
						}

						detJ = (dxyzdgez_i[0]*dxyzdgez_j[1] - dxyzdgez_i[1]*dxyzdgez_j[0]) + (dxyzdgez_i[1]*dxyzdgez_j[2] - dxyzdgez_i[2]*dxyzdgez_j[1]) + (dxyzdgez_i[2]*dxyzdgez_j[0] - dxyzdgez_i[0]*dxyzdgez_j[2]);

						//printf("detJ = %lf\ndxyzdge_i[0] = %lf dxyzdge_i[1] = %lf  dxyzdge_i[2] = %lf\ndxyzdge_j[0] = %lf dxyzdge_j[1] = %lf  dxyzdge_j[2] = %lf\n\n",detJ, dxyzdgez_i[0], dxyzdgez_i[1], dxyzdgez_i[2], dxyzdgez_j[0], dxyzdgez_j[1], dxyzdgez_j[2]);
						if (type_load < 3)
						{
							for (ic = 0; ic < (Order[iPatch][0] + 1) * (Order[iPatch][1] + 1) * (Order[iPatch][2] + 1); ic++)
							{
								//printf("Order[%d][0];%d,Order[%d][1]:%d\n",iPatch,Order[iPatch][0],iPatch,Order[iPatch][1]);
								sfc = Shape_func(ic, Local_Coord, No_Element_for_Integration[iii][jjj]);

								Equivalent_Nodal_Force[iControlpoint[ic]][type_load] += valDistLoad * sfc * detJ * w_1d[ig_i] * w_1d[ig_j];
								//(2019_06_13)printf("Check ic = %d sfc = %f   Weight[ig] = %f  valDistLoad = %f\n",ic,sfc,Weight[ig],valDistLoad);
								//(2019_06_13)printf("Equivalent_Nodal_Force[%d][%d]:%le\n",iControlpoint[ic],type_load,Equivalent_Nodal_Force[iControlpoint[ic]]	[type_load] );
							}
						}


						///////////////////////////符号は検討が必要
						if (type_load == 3) //法線方向
						{
							double LoadDir[3];
							LoadDir[0] = (dxyzdgez_i[1] * dxyzdgez_j[2] + dxyzdgez_i[2] * dxyzdgez_j[1]) / detJ;
							LoadDir[1] = (dxyzdgez_i[2] * dxyzdgez_j[0] + dxyzdgez_i[0] * dxyzdgez_j[2]) / detJ;
							LoadDir[2] = (dxyzdgez_i[0] * dxyzdgez_j[1] + dxyzdgez_i[1] * dxyzdgez_j[0]) / detJ;
							//printf("LoadDir[0] = %lf\tLoadDir[1] = %lf\tLoadDir[2] = %lf\n",LoadDir[0], LoadDir[1], LoadDir[2]);
							for (ic = 0; ic < (Order[iPatch][0] + 1) * (Order[iPatch][1] + 1) * (Order[iPatch][2] + 1); ic++)
							{
								sfc = Shape_func(ic, Local_Coord, No_Element_for_Integration[iii][jjj]);
								Equivalent_Nodal_Force[iControlpoint[ic]][0] +=          LoadDir[0] * valDistLoad * sfc * detJ * w_1d[ig_i] * w_1d[ig_j];
								Equivalent_Nodal_Force[iControlpoint[ic]][1] +=          LoadDir[1] * valDistLoad * sfc * detJ * w_1d[ig_i] * w_1d[ig_j];
								Equivalent_Nodal_Force[iControlpoint[ic]][2] += (-1.0) * LoadDir[2] * valDistLoad * sfc * detJ * w_1d[ig_i] * w_1d[ig_j];
								//printf("Equivalent_Nodal_Force[%d][0]=%lf\t[%d][1]=%lf\t[%d][2]=%lf\n",iControlpoint[ic],Equivalent_Nodal_Force[iControlpoint[ic]][0],iControlpoint[ic],Equivalent_Nodal_Force[iControlpoint[ic]][1],iControlpoint[ic],Equivalent_Nodal_Force[iControlpoint[ic]][2]);

								//printf("LoadDir[0]*(0.5-0.5*cos(2*sita))*sfc*detJ*Weight[%d]=%lf*%lf*%lf*%lf*%lf=%lf\n",ig,LoadDir[0],(0.5-0.5*cos(2*sita)),sfc,detJ,Weight[ig],LoadDir[0] * (0.5-0.5*cos(2*sita)) * sfc * detJ * Weight[ig]);
								//printf("LoadDir[1]*(0.5-0.5*cos(2*sita))*sfc*detJ*Weight[%d]=%lf*%lf*%lf*%lf*%lf=%lf\n",ig,LoadDir[1],(0.5-0.5*cos(2*sita)),sfc,detJ,Weight[ig],LoadDir[1] * (0.5-0.5*cos(2*sita)) * sfc * detJ * Weight[ig]);
							}
						}
					}
				}
			} //A
		}
	}	 //B
}


