#ifndef HEADER_MC3D_H
#define HEADER_MC3D_H

#define MERGE_DISTANCE 1.0e-13      // コントロールポイントが同じ点と判定する距離

static int counter = 0, KV_to_here = 0, CP_to_here = 0, B_to_here = 0, CP_result_to_here = 0;

FILE *fp;

struct information {
    int DIMENSION;
    int Total_patch;
    double E_and_nu[2];
    int *disp_constraint_n;
    int *disp_constraint_edge_n;
    double *disp_constraint_amount;
    int *disp_constraint;
    int distributed_load_n;
    double *distributed_load_info;
    int MAX_DISP_CONSTRAINT;
    int MAX_DISP_CONSTRAINT_EDGE;

    int *Oeder;
    int *KV_info;
    int *CP_info;
    double *CP;
    double *CP_result;
    int *A;
    double *B;
    int *Connectivity;
    double *KV;
    int *Face_Edge_info;
    int *Opponent_patch_num;

};

// get input data
void Get_inputdata_boundary_0(char *filename, information *info);
void Get_inputdata_boundary_1(char *filename, information *info);
void Get_inputdata_patch_0(char *filename, int *temp_Order, int *temp_KV_info, int *temp_CP_info);
void Get_inputdata_patch_1(char *filename, double *temp_KV, double *temp_CP, int *temp_A, double *temp_B, int *temp_KV_info, int *temp_CP_info, int num);
// make connectivity
void Check_B(int num_own, int num_opponent, double *temp_B, int *temp_Edge_info, int *temp_Opponent_patch_num);
void Make_connectivity(int num, int *temp_CP_info, int *temp_Edge_info, int *temp_Opponent_patch_num, int *temp_Connectivity, int *temp_A, double *temp_CP, double *temp_CP_result);
// output
void Output_inputdata(int *temp_Order, int *temp_KV_info, int *temp_CP_info, int *temp_Connectivity, double *temp_KV, double *temp_CP_result,
                      int *temp_Boundary_result, int *temp_length_before, int *temp_length_after, int total_disp_constraint_n);
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