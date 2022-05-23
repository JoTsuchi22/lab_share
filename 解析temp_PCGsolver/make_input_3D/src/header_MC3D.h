#ifndef HEADER_MC3D_H
#define HEADER_MC3D_H

#define MERGE_DISTANCE 1.0e-13      // コントロールポイントが同じ点と判定する距離
// #define DIMENSION 2                 // 2次元
// #define MERGE_DISTANCE 1.0e-13      // コントロールポイントが同じ点と判定する距離
#define MAX_DISP_CONSTRAINT 10      // 変位指定する変位量の最大個数
#define MAX_DISP_CONSTRAINT_EDGE 10 // 変位指定する辺の最大個数
#define MAX_DISTRIBUTED_LOAD 5      // 分布荷重の最大個数

// static int disp_constraint_n[DIMENSION];
// static int disp_constraint_edge_n[DIMENSION][MAX_DISP_CONSTRAINT];
// static double disp_constraint_amount[DIMENSION][MAX_DISP_CONSTRAINT];
// static int disp_constraint[DIMENSION][MAX_DISP_CONSTRAINT][MAX_DISP_CONSTRAINT_EDGE][3];
// static int distributed_load_n;
// static double distributed_load_info[MAX_DISTRIBUTED_LOAD][9];
// static int counter = 0;
// static int KV_to_here, CP_to_here, CP_result_to_here;
// static int A_to_own, A_to_opponent, B_to_here;

FILE *fp;

static 

struct info_global {
    int DIMENSION;
    int Total_patch;
    double E_and_nu[2];
    int *disp_constraint_n;
    int *disp_constraint_edge_n;
    double *disp_constraint_amount;
    int *disp_constraint;
    int distributed_load_n;
    double *distributed_load_info;
};

struct info_each_DIMENSION {
    int CP_n;
    int Order;
    int CP_n_before;
    int Order_before;
    int OE_n;
    int knot_n;
    int KI_cp_n;
    int KI_non_uniform_n;
    double *CP;
    double *KV;
    double *temp_CP;
    double *temp_KV;
    double *insert_knot;
};

// get input data
void Get_inputdata_boundary_0(char *filename, info_global *info_glo);
void Get_inputdata_boundary_1(char *filename, info_global *info_glo);
void Get_inputdata_patch_0(char *filename, int *temp_Order, int *temp_KV_info, int *temp_CP_info);
void Get_inputdata_patch_1(char *filename, double *temp_KV, double *temp_CP, int *temp_A, double *temp_B, int *temp_KV_info, int *temp_CP_info, int num);
// make connectivity
void Check_B(int num_own, int num_opponent, double *temp_B, int *temp_Edge_info, int *temp_Opponent_patch_num);
void Make_connectivity(int num, int *temp_CP_info, int *temp_Edge_info, int *temp_Opponent_patch_num, int *temp_Connectivity, int *temp_A, double *temp_CP, double *temp_CP_result);
// output
void Output_inputdata(int *temp_Order, int *temp_KV_info, int *temp_CP_info, int *temp_Connectivity, double *temp_KV, double *temp_CP_result,
                      int *temp_Boundary_result, int *temp_length_before, int *temp_length_after, int total_disp_constraint_n);
void Output_by_Gnuplot(double *temp_CP_result);
void Output_SVG(double *temp_B, double *temp_CP_result);
// heap sort
void Sort(int n, int *temp_CP_info, int *temp_A, int *temp_Boundary, int *temp_Boundary_result, int *temp_length_before, int *temp_length_after);
void swap(int *a, int *b);
int getLeft(int parent);
int getRight(int parent);
int getParent(int child);
void addHeap(int *a, int size);
void removeHeap(int *a, int size);
void makeHeap(int *a, int num);
void heapSort(int *a, int num);
void Dedupe(int *a, int *num, int *a_new, int *num_new, int n);

#endif