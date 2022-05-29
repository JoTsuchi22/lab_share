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
    int *disp_constraint_face_edge_n;
    double *disp_constraint_amount;
    int *disp_constraint;
    int distributed_load_n;
    double *distributed_load_info;
    int MAX_DISP_CONSTRAINT;
    int MAX_DISP_CONSTRAINT_FACE_EDGE;

    int *Order;
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

    int *length_before;
    int *length_after;
    int *Boundary;
    int *Boundary_result;
};

// get input data
void Get_inputdata_boundary_0(char *filename, information *info);
void Get_inputdata_boundary_1(char *filename, information *info);
void Get_inputdata_patch_0(char *filename, information *info);
void Get_inputdata_patch_1(char *filename, information *info, int num);
// make connectivity
void Check_B_2D(int num_own, int num_opponent, information *info);
void Check_B_3D(int num_own, int num_opponent, information *info);
void Make_connectivity_2D(int num, information *info);
void Make_connectivity_3D(int num, information *info);
// output
void Output_inputdata(int total_disp_constraint_n, const information *info);
void Output_SVG(const information *info);
// heap sort
void Sort(int n, information *info);
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