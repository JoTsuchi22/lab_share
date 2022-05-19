#ifndef HEADER_MI3D_H
#define HEADER_MI3D_H

// #define MAX_DIMENSION 3
// #define MAX_ORDER 3
// #define MAX_N_Controlpoint_each_parameter 1000
// #define MAX_N_Controlpoint_in_Patch MAX_N_Controlpoint_each_parameter * MAX_N_Controlpoint_each_parameter
// #define MAX_N_KNOT MAX_N_Controlpoint_each_parameter + MAX_ORDER + 1

// static int Dimension[MAX_N_INPUTFILE];
// static int Total_Control_Point[MAX_N_INPUTFILE];
// static int Order[MAX_N_INPUTFILE][MAX_DIMENSION];
// static int Order_before[MAX_N_INPUTFILE][MAX_DIMENSION];
// static int knot_n[MAX_N_INPUTFILE][MAX_DIMENSION];
// // static int Control_point_n[MAX_N_INPUTFILE][MAX_DIMENSION];
// // static int Control_point_n_before[MAX_N_INPUTFILE][MAX_DIMENSION];
// static double knot[MAX_N_INPUTFILE][MAX_DIMENSION][MAX_N_KNOT];
// // static double x[MAX_N_INPUTFILE][MAX_N_Controlpoint_in_Patch];
// // static double y[MAX_N_INPUTFILE][MAX_N_Controlpoint_in_Patch];
// // static double w[MAX_N_INPUTFILE][MAX_N_Controlpoint_in_Patch];
// static int OE_n[MAX_N_INPUTFILE][MAX_DIMENSION];
// // static int KI_uniform_n[MAX_N_INPUTFILE][MAX_DIMENSION];
// static int KI_cp_n[MAX_N_INPUTFILE][MAX_DIMENSION];
// static int KI_non_uniform_n[MAX_N_INPUTFILE][MAX_DIMENSION];
// static double insert_knot[MAX_N_INPUTFILE][MAX_DIMENSION][MAX_N_KNOT];
// static double insert_knot_in_KI[MAX_N_INPUTFILE][MAX_DIMENSION][MAX_N_KNOT];
// static double removal_knot[MAX_N_INPUTFILE][MAX_DIMENSION][MAX_N_KNOT];
// static int vec_length1[MAX_N_INPUTFILE][MAX_DIMENSION];
// static int vec_length2[MAX_N_INPUTFILE][MAX_DIMENSION];
// static double temp_knot1[MAX_N_KNOT];
// static double temp_knot2[MAX_N_KNOT];
// static double x_array[MAX_N_Controlpoint_each_parameter][MAX_N_Controlpoint_each_parameter];
// static double y_array[MAX_N_Controlpoint_each_parameter][MAX_N_Controlpoint_each_parameter];
// static double w_array[MAX_N_Controlpoint_each_parameter][MAX_N_Controlpoint_each_parameter];
// static double temp_x_array[MAX_N_Controlpoint_each_parameter];
// static double temp_y_array[MAX_N_Controlpoint_each_parameter];
// static double temp_w_array[MAX_N_Controlpoint_each_parameter];
// static double temp_x1[MAX_N_Controlpoint_in_Patch];
// static double temp_y1[MAX_N_Controlpoint_in_Patch];
// static double temp_w1[MAX_N_Controlpoint_in_Patch];
// static double temp_x2[MAX_N_Controlpoint_in_Patch];
// static double temp_y2[MAX_N_Controlpoint_in_Patch];
// static double temp_w2[MAX_N_Controlpoint_in_Patch];
// static int temp_Bezier_array[MAX_N_Controlpoint_each_parameter];
static int counter;
// static double Bezier_x[MAX_N_Controlpoint_each_parameter][MAX_ORDER];
// static double Bezier_y[MAX_N_Controlpoint_each_parameter][MAX_ORDER];
// static double Bezier_w[MAX_N_Controlpoint_each_parameter][MAX_ORDER];
static int number_of_Bezier_line;

FILE *fp;

struct info_global {
    int DIMENSION;
    int Total_Control_Point;
    int mode;

    double *Weight;
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

// Get input data
void Get_DIM(char *filename, info_global info_glo);
void Get_InputData_1(char *filename, information info_glo, info_each_DIMENSION info);
void Get_InputData_2(char *filename, information info_glo, info_each_DIMENSION info);
// Knot Insertion
void KI_calc_knot_1D(int insert_parameter_axis);
void KI_calc_T_1D(int insert_parameter_axis);
void KI_update_point_array(int line_number, int insert_parameter_axis);
void KI_update(int insert_parameter_axis);
void KI_reset_array();
void KI_non_uniform(int insert_parameter_axis, int KI_non_uniform, info_global info_glo, info_each_DIMENSION info);
void Calc_cp_insert_knot(int insert_parameter_axis);
void KI_cp(int insert_parameter_axis);
// Order Elevation
void OE_calc_point_array();
void Calc_insert_knot_in_OE(int elevation_parameter_axis);
void OE_define_temp_point_array(int line_number, int elevation_parameter_axis);
void Bezier_Order_Elevation(int elevation_parameter_axis, int Bezier_line_number);
void Bezier_update_point_array(int line_number, int elevation_parameter_axis);
void Bezier_update(int elevation_parameter_axis);
void Calc_Bezier(int elevation_parameter_axis, int other_axis);
void KR_calc_point_array();
void KR_calc_knot_1D(int removal_parameter_axis);
void KR_define_temp_point_array(int line_number, int removal_parameter_axis);
void KR_calc_Tinv_1D(int removal_parameter_axis);
void KR_update_point_array(int line_number, int removal_parameter_axis);
void KR_update(int removal_parameter_axis);
void KR_reset_array();
void KR_non_uniform(int removal_parameter_axis);
void OE(int elevation_parameter_axis);
// Output
void OutputData(char *filename);
// DBG
void Debug_printf(char *section);


#endif