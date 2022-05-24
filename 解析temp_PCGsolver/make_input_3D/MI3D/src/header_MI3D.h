#ifndef HEADER_MI3D_H
#define HEADER_MI3D_H

FILE *fp;

struct info_global {
    int DIMENSION;
    int Total_Control_Point;
    int mode;
    double *Weight;
    double *temp_Weight;
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

struct line_coordinate {
    double *line;
    double *new_line;
};

struct line_weight {
    double *line;
    double *new_line;
};

struct Bezier_coordinate {
    double *line;
    double *temp_line;
};

struct Bezier_weight {
    double *line;
    double *temp_line;
};

// Get input data
void Get_DIM(char *filename, info_global *info_glo);
void Get_InputData_1(char *filename, info_global *info_glo, info_each_DIMENSION *info);
void Get_InputData_2(char *filename, info_global *info_glo, info_each_DIMENSION *info);
// Knot Insertion
void KI_non_uniform(int insert_axis, int insert_knot_n, double *insert_knot_in_KI, info_global *info_glo, info_each_DIMENSION *info);
void KI_calc_knot_1D(int insert_axis, int insert_knot_n, double *insert_knot_in_KI, info_each_DIMENSION *info);
void KI_calc_T_1D(int insert_axis, int insert_knot_n, info_global *info_glo, info_each_DIMENSION *info, line_weight *w, line_coordinate *DIM);
void KI_cp(int insert_axis, info_global *info_glo, info_each_DIMENSION *info);
// Order Elevation
void OE(int elevation_axis, info_global *info_glo, info_each_DIMENSION *info);
void Calc_insert_knot_in_OE(int elevation_axis, int *insert_knot_n, double *insert_knot, int *removal_knot_n, double *removal_knot, info_each_DIMENSION *info);
void Calc_Bezier(int elevation_axis, info_global *info_glo, info_each_DIMENSION *info, line_weight *w, line_coordinate *DIM);
void Bezier_Order_Elevation(int elevation_axis, int Bezier_line, int *counter, info_global *info_glo, info_each_DIMENSION *info, line_weight *w, line_coordinate *DIM);
// Knot Removal
void KR_non_uniform(int removal_axis, int removal_knot_n, double *removal_knot, info_global *info_glo, info_each_DIMENSION *info);
void KR_calc_knot_1D(int removal_axis, int removal_knot_n, double *removal_knot, info_each_DIMENSION *info);
void KR_calc_Tinv_1D(int removal_axis, int removal_knot_n, info_global *info_glo, info_each_DIMENSION *info, line_weight *w, line_coordinate *DIM);
// Output
void OutputData(char *filename, info_global *info_glo, info_each_DIMENSION *info);

#endif