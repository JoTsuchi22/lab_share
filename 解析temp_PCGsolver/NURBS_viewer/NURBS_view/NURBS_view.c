#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <png.h>
#include <float.h>

#include <autogl.h>		//AutoGL ”おてがる”を使うためのヘッダファイル
#include "constant.h"
#include "pick.h"
#include "view.h"
#include "plot.h"
//#include "option.h"
//#include "mesh.h"
//#include "field.h"
//#include "point.h"
#include "image.h"
//#include "animation.h"
#include "draw.h"

#include "clickcalc.h"
#include "NURBS_calc.h"

int MeshIsShown = 1;		//メッシュを表示するか
int EdgeIsShown = 0;		//表面を表示するか
int CntlNetIsShown = 0;	//コントロールネットを表示するか
int Cntl_PIsShown = 0;	//コントロールポイントを表示するか
int DefmMeshIsShown = 0;
int DefmEdgeIsShown = 0;
int ContourDispxIsShown = 1;
int ContourDispyIsShown = 0;
int ContourDispvalIsShown = 0;
int ContourStrainxxIsShown = 0;
int ContourStrainyyIsShown = 0;
int ContourStrainxyIsShown = 0;
int ContourStressxxIsShown = 0;
int ContourStressyyIsShown = 0;
int ContourStressxyIsShown = 0;
int patch_is_shown[MAX_PATCHES];

int contour_volume_type = DISP_X;
int contour_volume_type_flag = -1;
int format_flag_10 = DEFAULT_FORMAT_FLAG_10;

int picked_point_flag = 0;

///////////////////////////////////////////////////////////////////
char *filename_mesh = NULL;
int num_filenames_fields = 0;
char **filenames_fields = NULL;
char filename_visible_elements[256] = "";
char filename_invisible_elements[256] = "";
char filename_points[256] = "";
char filename_points_fields[256] = "";
char filename_restart[256] = "";
int batch_mode_flag = 0;
char filename_view_image[256] = DEFAULT_FILENAME_VIEW_IMAGE;
int window_size = DEFAULT_WINDOW_SIZE;
int division_number = DEFAULT_DIVISION_NUMBER;

int current_time_step = DEFAULT_FIRST_TIME_STEP;
int first_time_step = DEFAULT_FIRST_TIME_STEP;
int last_time_step = DEFAULT_FIRST_TIME_STEP;
int time_step_increment = DEFAULT_TIME_STEP_INCREMENT;
int time_step_loop_flag = DEFAULT_TIME_STEP_LOOP_FLAG;
int slow_down_level = DEFAULT_SLOW_DOWN_LEVEL;
char filename_all_view_images[256] = DEFAULT_FILENAME_ALL_VIEW_IMAGES;

int auto_view_position_flag = 1;
double view_size = 0.0;
double view_center_x = 0.0;
double view_center_y = 0.0;
double view_center_z = 0.0;
double view_direction_x = 0.0;
double view_direction_y = 0.0;
double view_direction_z = 0.0;
double view_up_vector_x = 0.0;
double view_up_vector_y = 0.0;
double view_up_vector_z = 0.0;
double section_a = 0.0;
double section_b = 0.0;
double section_c = 0.0;
double section_d = 0.0;

int element_type = -1;
int num_elements = 0;
int num_element_nodes = 0;
int num_element_triangles = 0;
int num_element_quadrangles = 0;
int **elements = NULL;
int num_visible_elements = 0;
int *visible_elements = NULL;
int num_visible_tetrahedrons = 0;
int num_visible_hexahedrons = 0;
int *visible_tetrahedrons = NULL;
int *visible_hexahedrons = NULL;
int num_triangles = 0;
int num_quadrangles = 0;
int **triangles = NULL;
int **quadrangles = NULL;
int num_visible_triangles = 0;
int num_visible_quadrangles = 0;
int *visible_triangles = NULL;
int *visible_quadrangles = NULL;
int num_triangle_edges = 0;
int num_quadrangle_edges = 0;
int *triangle_edges = NULL;
int *quadrangle_edges = NULL;
int num_crack_triangles = 0;
int num_crack_quadrangles = 0;
int *crack_triangles = NULL;
int *crack_quadrangles = NULL;
int num_crack_front_triangle_edges = 0;
int num_crack_front_quadrangle_edges = 0;
int *crack_front_triangle_edges = NULL;
int *crack_front_quadrangle_edges = NULL;
int num_segments = 0;
int **segments = NULL;
int num_nodes = 0;
int num_dimensions = 0;
double **nodes = NULL;

double *contour_field = NULL;
int num_nodal_fields = 0;
double **nodal_fields = NULL;
char **nodal_field_labels = NULL;
int num_elemental_fields = 0;
double **elemental_fields = NULL;
char **elemental_field_labels = NULL;
int num_integer_elemental_fields = 0;
int **integer_elemental_fields = NULL;
double *casted_integer_elemental_field = NULL;
char **integer_elemental_field_labels = NULL;

double *displacement_fields[3] = {NULL, NULL, NULL};
int num_displacement_fields = 0;
int *displacement_field_offsets = NULL;
char **displacement_field_labels = NULL;

int num_points = 0;
double **points = NULL;
int num_points_fields = 0;
double **points_fields = NULL;
char **points_field_labels = NULL;

int plot_type = DEFAULT_PLOT_TYPE;
int mesh_flag = DEFAULT_MESH_FLAG;
int edges_flag = DEFAULT_EDGES_FLAG;
int crack_mesh_flag = DEFAULT_CRACK_MESH_FLAG;
int crack_front_edges_flag = DEFAULT_CRACK_FRONT_EDGES_FLAG;
int segments_flag = DEFAULT_SEGMENTS_FLAG;
int deformed_mesh_flag = DEFAULT_DEFORMED_MESH_FLAG;
int deformed_edges_flag = DEFAULT_DEFORMED_EDGES_FLAG;
int deformed_crack_mesh_flag = DEFAULT_DEFORMED_CRACK_MESH_FLAG;
int deformed_crack_front_edges_flag = DEFAULT_DEFORMED_CRACK_FRONT_EDGES_FLAG;
int deformed_segments_flag = DEFAULT_DEFORMED_SEGMENTS_FLAG;
double deformation_scale_factor = DEFAULT_DEFORMATION_SCALE_FACTOR;
int points_type = DEFAULT_POINTS_TYPE;
char title[256] = DEFAULT_TITLE;
int coord_axes_flag = DEFAULT_COORD_AXES_FLAG;
int coord_axes_at_origin_flag = DEFAULT_COORD_AXES_AT_ORIGIN_FLAG;
int other_mode = DEFAULT_OTHER_MODE;
int marked_node = -1;
int marked_element = -1;
int marked_point = -1;
int scale_bar_flag = DEFAULT_SCALE_BAR_FLAG;
char scale_bar_format[256] = DEFAULT_SCALE_BAR_FORMAT;
double scale_bar_scale_factor = DEFAULT_SCALE_BAR_SCALE_FACTOR;
int section_flag = DEFAULT_SECTION_FLAG;
int section_type = DEFAULT_SECTION_TYPE;
char volume_color[256] = DEFAULT_VOLUME_COLOR;
char mesh_color[256] = DEFAULT_MESH_COLOR;
char edge_color[256] = DEFAULT_EDGE_COLOR;
char crack_mesh_color[256] = DEFAULT_CRACK_MESH_COLOR;
char crack_front_edge_color[256] = DEFAULT_CRACK_FRONT_EDGE_COLOR;
char segment_color[256] = DEFAULT_SEGMENT_COLOR;
char point_color[256] = DEFAULT_POINT_COLOR;
int mesh_line_width = DEFAULT_MESH_LINE_WIDTH;
int edge_line_width = DEFAULT_EDGE_LINE_WIDTH;
int crack_mesh_line_width = DEFAULT_CRACK_MESH_LINE_WIDTH;
int crack_front_edge_line_width = DEFAULT_CRACK_FRONT_EDGE_LINE_WIDTH;
int segment_line_width = DEFAULT_SEGMENT_LINE_WIDTH;
int point_size = DEFAULT_POINT_SIZE;
int contour_type = DEFAULT_CONTOUR_TYPE;
int reversed_contour_colors_flag = DEFAULT_REVERSED_CONTOUR_COLORS_FLAG;
int gray_contour_color_below_magenta_one_flag = DEFAULT_GRAY_CONTOUR_COLOR_BELOW_MAGENTA_ONE_FLAG;
int num_contour_grades = DEFAULT_NUM_CONTOUR_GRADES;
int contour_bar_flag = DEFAULT_CONTOUR_BAR_FLAG;
char contour_bar_label[256] = "";
char contour_bar_label2[256] = "";
int auto_contour_bar_format_flag = DEFAULT_AUTO_CONTOUR_BAR_FORMAT_FLAG;
char contour_bar_format[256] = DEFAULT_CONTOUR_BAR_FORMAT;
double contour_bar_scale_factor = DEFAULT_CONTOUR_BAR_SCALE_FACTOR;
int auto_contour_field_range_flag = 1;
double contour_field_max = 1.0;
double contour_field_min = 0.0;
int integer_contour_field_mod = DEFAULT_INTEGER_CONTOUR_FIELD_MOD;
int contour_field_id = -1;
int displacement_field_id = -1;
int point_mark = DEFAULT_POINT_MARK;
int points_reversed_contour_colors_flag = DEFAULT_POINTS_REVERSED_CONTOUR_COLORS_FLAG;
int points_contour_bar_flag = DEFAULT_POINTS_CONTOUR_BAR_FLAG;
char points_contour_bar_label[256] = "";
char points_contour_bar_label2[256] = "";
int points_auto_contour_bar_format_flag = DEFAULT_POINTS_AUTO_CONTOUR_BAR_FORMAT_FLAG;
char points_contour_bar_format[256] = DEFAULT_POINTS_CONTOUR_BAR_FORMAT;
double points_contour_bar_scale_factor = DEFAULT_POINTS_CONTOUR_BAR_SCALE_FACTOR;
int auto_points_field_range_flag = 1;
double points_field_max = 1.0;
double points_field_min = 0.0;
int points_field_id = -1;
/////////////////////////////////////////////////////////////////////////

double max_x = -DBL_MAX;
double min_x = DBL_MAX;
double max_y = -DBL_MAX;
double min_y = DBL_MAX;

double MaxRange_disp_x = -DBL_MAX;
double MinRange_disp_x = DBL_MAX;
double MaxRange_disp_y = -DBL_MAX;
double MinRange_disp_y = DBL_MAX;
double MaxRange_disp_val = -DBL_MAX;
double MinRange_disp_val = DBL_MAX;

double MaxRange_strain_xx = -DBL_MAX;
double MinRange_strain_xx = DBL_MAX;
double MaxRange_strain_yy = -DBL_MAX;
double MinRange_strain_yy = DBL_MAX;
double MaxRange_strain_xy = -DBL_MAX;
double MinRange_strain_xy = DBL_MAX;

double MaxRange_stress_xx = -DBL_MAX;
double MinRange_stress_xx = DBL_MAX;
double MaxRange_stress_yy = -DBL_MAX;
double MinRange_stress_yy = DBL_MAX;
double MaxRange_stress_xy = -DBL_MAX;
double MinRange_stress_xy = DBL_MAX;

double E;						//ヤング率(GPa)
double nu;						//ポアソン比(-)

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

double coord_x[MAX_PATCHES][MAX_POINTS][MAX_POINTS];		//変形前x座標
double coord_y[MAX_PATCHES][MAX_POINTS][MAX_POINTS];		//変形前y座標

double disp_x[MAX_PATCHES][MAX_POINTS][MAX_POINTS];		//x方向変位
double disp_y[MAX_PATCHES][MAX_POINTS][MAX_POINTS];		//y方向変位
double disp_val[MAX_PATCHES][MAX_POINTS][MAX_POINTS];		//変位の絶対値

double strain_xx[MAX_PATCHES][MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//x方向ひずみ
double strain_yy[MAX_PATCHES][MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//y方向ひずみ
double strain_xy[MAX_PATCHES][MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//剪断ひずみ

double stress_xx[MAX_PATCHES][MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//x方向応力
double stress_yy[MAX_PATCHES][MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//y方向応力
double stress_xy[MAX_PATCHES][MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];	//剪断応力

int fields_flag;
int division_ele_xi;					//ξ方向の一要素あたりの分割数
int division_ele_eta;					//η方向の一要素あたりの分割数
int division_n_xi[MAX_PATCHES];		//ξ方向の表示する点の数
int division_n_eta[MAX_PATCHES];		//η方向の表示する点の数
int element_n_xi[MAX_PATCHES];		//ξ方向要素数
int element_n_eta[MAX_PATCHES];		//η方向要素数

int temp_index[MAX_PATCHES][MAX_CNRL_P];
double temp_cntl_px[MAX_CNRL_P];
double temp_cntl_py[MAX_CNRL_P];
double temp_weight[MAX_CNRL_P];
double temp_disp_x[MAX_CNRL_P];
double temp_disp_y[MAX_CNRL_P];

FILE *fp; // FILE型構造体

static void destroyAndQuit(void);
static void redrawView(void);
static void setContourFieldRangeAndRedrawView(void);
//static void setPointsFieldRangeAndRedrawView(void);
static void handleEvent(void);
static void setPanelMain(void);
static void UpdateModel(void);
static void setPanelSelectPatch(void);
//static void setPanelOtherModeAndMarks(void);
static void setPanelScaleBar(void);
//static void setPanelSection(void);
//static void setPanelColorsAndSizes(void);
static void setPanelContourFieldOptions(void);
static void setPanelContourFieldSelection(void);
//static void setPanelDisplacementFieldSelection(void);
//static void setPanelPoints(void);
//static void setPanelViewPosition(void);
static void setPanelOutput(void);
//static void setPanelAnimation(void);

/* 読み込み */
void ReadFile (char *filename_input, char *filename_view, char *filename_disp) {
	int i, j, k, l, m, temp_int;
	double temp1, temp2, temp3;

	fp = fopen(filename_input, "r");
	printf("Start Reading input\n");

	fscanf(fp, "%lf%lf", &E, &nu);
	fscanf(fp, "\n");
	printf("E nu:% 1.8e % 1.8e\n", E, nu);

	fscanf(fp, "%d", &patch_n);
	fscanf(fp, "\n");
	printf("patches:%d \n", patch_n);

	fscanf(fp, "%d", &cntl_p_n);
	fscanf(fp, "\n");
	printf("control points:%d \n", cntl_p_n);

	for (i = 0; i < patch_n; i++) {
		fscanf(fp, "%d%d", &order_xi[i], &order_eta[i]);
		fscanf(fp, "\n");
		printf("order%d:%d %d\n", i, order_xi[i], order_eta[i]);
	}


	for (i = 0; i < patch_n; i++) {
		fscanf(fp, "%d%d", &knot_n_xi[i], &knot_n_eta[i]);
		fscanf(fp, "\n");
		printf("knots%d:%d %d\n", i, knot_n_xi[i], knot_n_eta[i]);
	}

	for (i = 0; i < patch_n; i++) {
		fscanf(fp, "%d%d", &cntl_p_n_xi[i], &cntl_p_n_eta[i]);
		printf("control points%d (xi, eta):%d %d\n",
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


	printf("Start Reading points\n");
	fp = fopen(filename_view, "r");
	fscanf(fp, "%d\t%d\t%d\n",
	       &fields_flag, &division_ele_xi, &division_ele_eta);
	if (fields_flag) {
		for (i = 0; i < patch_n; i++) {
			fscanf(fp, "%d%d%d%d",
			       &division_n_xi[i], &division_n_eta[i],
			       &element_n_xi[i], &element_n_eta[i]);
			printf("%d\t%d\t%d\t%d\n",
			       division_n_xi[i], division_n_eta[i],
			       element_n_xi[i], element_n_eta[i]);
			for (j = 0; j < division_n_xi[i]; j++) {
				for (k = 0; k < division_n_eta[i]; k++) {
					fscanf(fp, "%le%le%le%le%le",
					       &coord_x[i][j][k], &coord_y[i][j][k],
					       &disp_x[i][j][k], &disp_y[i][j][k], &disp_val[i][j][k]);
					printf("%d\t%d\t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\n",
					       i, j, coord_x[i][j][k], coord_y[i][j][k],
					       disp_x[i][j][k], disp_y[i][j][k], disp_val[i][j][k]);
				}
			}

			for (j = 0; j < element_n_xi[i]; j++) {
				for (k = 0; k < element_n_eta[i]; k++) {
					for (l = 0; l < division_ele_xi + 1; l++) {
						for (m = 0; m < division_ele_eta + 1; m++) {
							fscanf(fp, "%le%le%le%le%le%le\n",
							       &strain_xx[i][j][k][l][m], &strain_yy[i][j][k][l][m],
							       &strain_xy[i][j][k][l][m],
							       &stress_xx[i][j][k][l][m], &stress_yy[i][j][k][l][m],
							       &stress_xy[i][j][k][l][m]);
							printf("[%d][%d][%d][%d][%d]\t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\t% 1.8e\n",
							       i, j, k, l, m,
							       strain_xx[i][j][k][l][m], strain_yy[i][j][k][l][m],
							       strain_xy[i][j][k][l][m],
							       stress_xx[i][j][k][l][m], stress_yy[i][j][k][l][m],
							       stress_xy[i][j][k][l][m]);
						}
					}
				}
			}
		}
	} else {
		for (i = 0; i < patch_n; i++) {
			fscanf(fp, "%d%d%d%d",
			       &division_n_xi[i], &division_n_eta[i],
			       &element_n_xi[i], &element_n_eta[i]);
			printf("%d\t%d\t%d\t%d\n",
			       division_n_xi[i], division_n_eta[i],
			       element_n_xi[i], element_n_eta[i]);
			for (j = 0; j < division_n_xi[i]; j++) {
				for (k = 0; k < division_n_eta[i]; k++) {
					fscanf(fp, "%le%le",
					       &coord_x[i][j][k], &coord_y[i][j][k]);
					printf("%d\t%d\t% 1.8e\t% 1.8e\n",
					       i, j, coord_x[i][j][k], coord_y[i][j][k]);
				}
			}
		}
	}
	fclose(fp);
	printf("End Reading points\n\n");

	if (fields_flag) {
		printf("Start Reading displacement\n\n");
		fp = fopen(filename_disp, "r");
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
}

void plotContourVolume (int size1, int size2, int size3,
                        int division_n_xi, int division_n_eta,
                        int element_n_xi, int element_n_eta,
                        double coord_x[size1][size1],
                        double coord_y[size1][size1],
                        double disp_x[size1][size1],
                        double disp_y[size1][size1],
                        double disp_val[size1][size1],
                        double strain_xx[size2][size2][size3 + 1][size3 + 1],
                        double strain_yy[size2][size2][size3 + 1][size3 + 1],
                        double strain_xy[size2][size2][size3 + 1][size3 + 1],
                        double stress_xx[size2][size2][size3 + 1][size3 + 1],
                        double stress_yy[size2][size2][size3 + 1][size3 + 1],
                        double stress_xy[size2][size2][size3 + 1][size3 + 1]) {
	switch (contour_volume_type) {
	case DISP_X:
		plotDispContour(MAX_POINTS,
		                division_n_xi, division_n_eta,
		                coord_x, coord_y,
		                disp_x, disp_y,
		                disp_x);
		break;
	case DISP_Y:
		plotDispContour(MAX_POINTS,
		                division_n_xi, division_n_eta,
		                coord_x, coord_y,
		                disp_x, disp_y,
		                disp_y);
		break;
	case DISP_VAL:
		plotDispContour(MAX_POINTS,
		                division_n_xi, division_n_eta,
		                coord_x, coord_y,
		                disp_x, disp_y,
		                disp_val);
		break;
	case STRAIN_XX:
		plotStrainStressContour(MAX_POINTS, MAX_ELEMENTS, MAX_DIVISION,
		                        element_n_xi, element_n_eta,
		                        coord_x, coord_y,
		                        disp_x, disp_y,
		                        strain_xx);
		break;
	case STRAIN_YY:
		plotStrainStressContour(MAX_POINTS, MAX_ELEMENTS, MAX_DIVISION,
		                        element_n_xi, element_n_eta,
		                        coord_x, coord_y,
		                        disp_x, disp_y,
		                        strain_yy);
		break;
	case STRAIN_XY:
		plotStrainStressContour(MAX_POINTS, MAX_ELEMENTS, MAX_DIVISION,
		                        element_n_xi, element_n_eta,
		                        coord_x, coord_y,
		                        disp_x, disp_y,
		                        strain_xy);
		break;
	case STRESS_XX:
		plotStrainStressContour(MAX_POINTS, MAX_ELEMENTS, MAX_DIVISION,
		                        element_n_xi, element_n_eta,
		                        coord_x, coord_y,
		                        disp_x, disp_y,
		                        stress_xx);
		break;
	case STRESS_YY:
		plotStrainStressContour(MAX_POINTS, MAX_ELEMENTS, MAX_DIVISION,
		                        element_n_xi, element_n_eta,
		                        coord_x, coord_y,
		                        disp_x, disp_y,
		                        stress_yy);
		break;
	case STRESS_XY:
		plotStrainStressContour(MAX_POINTS, MAX_ELEMENTS, MAX_DIVISION,
		                        element_n_xi, element_n_eta,
		                        coord_x, coord_y,
		                        disp_x, disp_y,
		                        stress_xy);
		break;
	}
}

void setContourFieldRange(void) {
	switch (contour_volume_type) {
	case DISP_X:
		strcpy(contour_bar_label, "Displacement x");
		contour_field_max = MaxRange_disp_x;
		contour_field_min = MinRange_disp_x;
		break;
	case DISP_Y:
		strcpy(contour_bar_label, "Displacement y");
		contour_field_max = MaxRange_disp_y;
		contour_field_min = MinRange_disp_y;
		break;
	case DISP_VAL:
		strcpy(contour_bar_label, "Displacement");
		contour_field_max = MaxRange_disp_val;
		contour_field_min = MinRange_disp_val;
		break;
	case STRAIN_XX:
		strcpy(contour_bar_label, "Strain xx");
		contour_field_max = MaxRange_strain_xx;
		contour_field_min = MinRange_strain_xx;
		break;
	case STRAIN_YY:
		strcpy(contour_bar_label, "Strain yy");
		contour_field_max = MaxRange_strain_yy;
		contour_field_min = MinRange_strain_yy;
		break;
	case STRAIN_XY:
		strcpy(contour_bar_label, "Strain xy");
		contour_field_max = MaxRange_strain_xy;
		contour_field_min = MinRange_strain_xy;
		break;
	case STRESS_XX:
		strcpy(contour_bar_label, "Stress xx");
		contour_field_max = MaxRange_stress_xx;
		contour_field_min = MinRange_stress_xx;
		break;
	case STRESS_YY:
		strcpy(contour_bar_label, "Stress yy");
		contour_field_max = MaxRange_stress_yy;
		contour_field_min = MinRange_stress_yy;
		break;
	case STRESS_XY:
		strcpy(contour_bar_label, "Stress xy");
		contour_field_max = MaxRange_stress_xy;
		contour_field_min = MinRange_stress_xy;
		break;
	}
}

static void DrawView(int patch_index) {
	if (Cntl_PIsShown) {
		AutoGL_SetColor (0.4, 0.4, 0.5);
		plotCntl_P(cntl_p_n_xi[patch_index], cntl_p_n_eta[patch_index],
		           cntl_px[patch_index], cntl_py[patch_index]);
	}
	if (CntlNetIsShown) {
		AutoGL_SetLineWidth((double)mesh_line_width * 0.5);
		AutoGL_SetColor (0.4, 0.4, 0.5);
		plotCntlNet(cntl_p_n_xi[patch_index], cntl_p_n_eta[patch_index],
		            cntl_px[patch_index], cntl_py[patch_index]);
	}
	if (DefmMeshIsShown) {
		AutoGL_SetLineWidth((double)mesh_line_width);
		AutoGL_SetColor (0, 0, 1);
		plotDeformedMesh(MAX_POINTS,
		                 division_n_xi[patch_index], division_n_eta[patch_index],
		                 element_n_xi[patch_index], element_n_eta[patch_index],
		                 coord_x[patch_index], coord_y[patch_index],
		                 disp_x[patch_index], disp_y[patch_index]);
	}
	if (DefmEdgeIsShown) {
		AutoGL_SetLineWidth((double)mesh_line_width);
		AutoGL_SetColor (0, 0, 1);
		plotDeformedEdge(MAX_POINTS,
		                 division_n_xi[patch_index], division_n_eta[patch_index],
		                 element_n_xi[patch_index], element_n_eta[patch_index],
		                 coord_x[patch_index], coord_y[patch_index],
		                 disp_x[patch_index], disp_y[patch_index]);
	}
	if (MeshIsShown) {
		AutoGL_SetLineWidth((double)mesh_line_width);
		AutoGL_SetColor (0, 0, 0);
		plotMesh(MAX_POINTS,
		         division_n_xi[patch_index], division_n_eta[patch_index],
		         element_n_xi[patch_index], element_n_eta[patch_index],
		         coord_x[patch_index], coord_y[patch_index]);
	}
	if (EdgeIsShown) {
		AutoGL_SetLineWidth((double)mesh_line_width);
		AutoGL_SetColor (0, 0, 0);
		plotEdge(MAX_POINTS,
		         division_n_xi[patch_index], division_n_eta[patch_index],
		         element_n_xi[patch_index], element_n_eta[patch_index],
		         coord_x[patch_index], coord_y[patch_index]);
	}

	switch (plot_type) {
	case NOTHING_PLOT_TYPE:
		break;
	case VOLUME_PLOT_TYPE:
	case DEFORMED_VOLUME_PLOT_TYPE:
		AutoGL_SetColor (0.7, 0.7, 0.7);
		plotVolume(MAX_POINTS,
		           division_n_xi[patch_index], division_n_eta[patch_index],
		           coord_x[patch_index], coord_y[patch_index],
		           disp_x[patch_index], disp_y[patch_index]);
		break;
	case CONTOUR_VOLUME_PLOT_TYPE:
	case DEFORMED_CONTOUR_VOLUME_PLOT_TYPE:
		plotContourVolume(MAX_POINTS, MAX_ELEMENTS, MAX_DIVISION,
		                  division_n_xi[patch_index], division_n_eta[patch_index],
		                  element_n_xi[patch_index], element_n_eta[patch_index],
		                  coord_x[patch_index], coord_y[patch_index],
		                  disp_x[patch_index], disp_y[patch_index],  disp_val[patch_index],
		                  strain_xx[patch_index], strain_yy[patch_index], strain_xy[patch_index],
		                  stress_xx[patch_index], stress_yy[patch_index], stress_xy[patch_index]);
		break;
	}
}

static void destroyAndQuit(void) {
	AutoGL_Terminate();
	exit(EXIT_SUCCESS);
}

static void redrawView (void) {
	int i;
	if (fields_flag) {
		if (contour_volume_type_flag != contour_volume_type) {
			contour_volume_type_flag = contour_volume_type;
			setContourFieldRange();
		}
	}
	for (i = 0; i < patch_n; i++) {
		if (patch_is_shown[i]) {
			DrawView(i);
		}
	}
	if (coord_axes_flag)
		plotCoordAxes();

	if (coord_axes_at_origin_flag)
		plotCoordAxesAtOrigin();

	if (scale_bar_flag)
		plotScaleBar();
}

static void setContourFieldRangeAndRedrawView(void) {
	setContourFieldRange();
	AutoGL_DrawView();
}

static void handleEvent(void) {
	int event = AutoGL_GetViewEvent();
	int i, j;

	int x_ind, y_ind, z_ind;
	double x, y, z;
	double output_xi, output_eta;

	switch (event) {
	case AUTOGL_EVENT_POINTING_DEVICE_RELEASE:
//		switch (other_mode) {
//		case NODE_PICK_OTHER_MODE:
//			handleNodePickEvent();
//			break;
//		case ELEMENT_PICK_OTHER_MODE:
//			handleElementPickEvent();
//			break;
//		case POINT_PICK_OTHER_MODE:
		AutoGL_GetPointingDevicePositionDc(&x_ind, &y_ind, &z_ind);
		AutoGL_GetPositionOfDc(&x, &y, &z, x_ind, y_ind, z_ind);
		x /= COORDS_SCALE_FACTOR;
		y /= COORDS_SCALE_FACTOR;
		z /= COORDS_SCALE_FACTOR;
		if (min_x <= x && x <= max_x && min_y <= y && y <= max_y) {
			for (i = 0; i < patch_n; i++) {
				if (patch_is_shown[i]) {
					j = CalcXiEtaByNR(x, y,
					                  knot_vec_xi[i], knot_vec_eta[i],
					                  cntl_px[i], cntl_py[i],
					                  disp_cntl_px[i], disp_cntl_py[i],
					                  cntl_p_n_xi[i], cntl_p_n_eta[i],
					                  weight[i], order_xi[i], order_eta[i],
					                  &output_xi, &output_eta);
				}
			}
		}
		printf("\n");
		/*
		for (i = 0; i < patch_n; i++) {
			if (patch_is_shown[i]) {
				handlePointPickEvent(i, &picked_point_flag);
			}
		}
		*/
		/*
		if (picked_point_flag) {
			AutoGL_DrawView();
			printf("\n");
			picked_point_flag = 0;
		}
		*/
//			break;
//		}
		break;
	case AUTOGL_EVENT_KEY:
//		if (num_dimensions == 3)
//			AutoGL_HandleDefaultKeyEventInMode3D();
//		else
		AutoGL_HandleDefaultKeyEventInMode2D();
		break;
	}
}

static void setPanelMain (void) {
	AutoGL_AddGroup("Main");

	AutoGL_AddCallback (UpdateModel, "UpdateModel");
	AutoGL_SetLabel("update model");

	AutoGL_SetPanelInMode2D();

	AutoGL_AddInteger(&plot_type, "plot_type");
	AutoGL_SetLabel("plot type");
	AutoGL_AddIntegerItem("nothing");
	AutoGL_AddIntegerItem("volume");
	if (fields_flag) {
		AutoGL_AddIntegerItem("contour volume");
		AutoGL_AddIntegerItem("deformed volume");
		AutoGL_AddIntegerItem("deformed contour volume");
	}

	AutoGL_AddBoolean(&MeshIsShown, "MeshIsShown");
	AutoGL_SetLabel("mesh");
	AutoGL_AddBoolean(&EdgeIsShown, "EdgeIsShown");
	AutoGL_SetLabel("edges");
	if (fields_flag) {
		AutoGL_AddBoolean(&DefmMeshIsShown, "DefmMeshIsShown");
		AutoGL_SetLabel("deformed mesh");
		AutoGL_AddBoolean(&DefmEdgeIsShown, "DefmEdgeIsShown");
		AutoGL_SetLabel("deformed edges");
	}
	AutoGL_AddBoolean (&Cntl_PIsShown, "Cntl_PIsShown");
	AutoGL_SetLabel("Control Points");
	AutoGL_AddBoolean (&CntlNetIsShown, "CntlNetIsShown");
	AutoGL_SetLabel("Control Net");

	if (fields_flag) {
		AutoGL_AddReal(&deformation_scale_factor, "deformation_scale_factor");
		AutoGL_SetLabel("deformation scale factor");
	}

	AutoGL_AddBoolean(&coord_axes_flag, "coord_axes_flag");
	AutoGL_SetLabel("coordinate axes");
	AutoGL_AddBoolean(&coord_axes_at_origin_flag, "coord_axes_at_origin_flag");
	AutoGL_SetLabel("coordinate axes at origin");

	AutoGL_AddCallback(destroyAndQuit, "destroyAndQuit");
	AutoGL_SetLabel("quit");
}

static void UpdateModel (void) {
	AutoGL_DrawView ();
}

static void setPanelSelectPatch (void) {
	AutoGL_AddGroup("Patch selection");
	int i;
	char name[256];
	for (i = 0; i < patch_n; i++) {
		patch_is_shown[i] = 1;
		snprintf(name, 256, "patch number %d", i);
		AutoGL_AddBoolean (&patch_is_shown[i], name);
	}
}

static void setPanelScaleBar(void) {
	AutoGL_AddGroup("Scale bar");

	AutoGL_AddBoolean(&scale_bar_flag, "scale_bar_flag");
	AutoGL_SetLabel("scale bar");
	AutoGL_AddString(scale_bar_format, "scale_bar_format", 256);
	AutoGL_SetLabel("scale bar format");
	AutoGL_AddReal(&scale_bar_scale_factor, "scale_bar_scale_factor");
	AutoGL_SetLabel("scale bar scale factor");
}

static void setPanelContourFieldOptions (void) {
	AutoGL_AddGroup("Contour field options");

	AutoGL_AddInteger(&contour_type, "contour_type");
	AutoGL_SetLabel("contour type");
	AutoGL_AddIntegerItem("smooth");
	AutoGL_AddIntegerItem("band");
	AutoGL_AddBoolean(&reversed_contour_colors_flag, "reversed_contour_colors_flag");
	AutoGL_SetLabel("reversed contour colors");
	AutoGL_AddBoolean(&gray_contour_color_below_magenta_one_flag, "gray_contour_color_below_magenta_one_flag");
	AutoGL_SetLabel("gray contour color below magenta one");
	AutoGL_AddBoolean(&format_flag_10, "format_flag_10");
	AutoGL_SetLabel("use * 10 format");
	AutoGL_AddInteger(&num_contour_grades, "num_contour_grades");
	AutoGL_SetLabel("number of contour grades");
	//AutoGL_AddBoolean(&contour_bar_flag, "contour_bar_flag");
	//AutoGL_SetLabel("contour bar");
	AutoGL_AddString(contour_bar_label, "contour_bar_label", 256);
	AutoGL_SetLabel("contour bar label");
	AutoGL_AddBoolean(&auto_contour_bar_format_flag, "auto_contour_bar_format_flag");
	AutoGL_SetLabel("auto contour bar format");
	AutoGL_AddString(contour_bar_format, "contour_bar_format", 256);
	AutoGL_SetLabel("contour bar format");
	AutoGL_AddReal(&contour_bar_scale_factor, "contour_bar_scale_factor");
	AutoGL_SetLabel("contour bar scale factor");

	AutoGL_AddCallback(setContourFieldRangeAndRedrawView, "setContourFieldRangeAndRedrawView");
	AutoGL_SetLabel("set contour field range automatically");
	AutoGL_AddReal(&contour_field_max, "contour_field_max");
	AutoGL_SetLabel("max");
	AutoGL_AddReal(&contour_field_min, "contour_field_min");
	AutoGL_SetLabel("min");
}

static void setPanelContourFieldSelection (void) {
	AutoGL_AddGroup("Contour field selection");
	AutoGL_AddInteger(&contour_volume_type, "contour_volume_type");
	AutoGL_SetLabel("Contour field\t\t\t");
	AutoGL_AddIntegerItem("Displacement x");
	AutoGL_AddIntegerItem("Displacement y");
	AutoGL_AddIntegerItem("Displacement");
	AutoGL_AddIntegerItem("Strain xx");
	AutoGL_AddIntegerItem("Strain yy");
	AutoGL_AddIntegerItem("Strain xy");
	AutoGL_AddIntegerItem("Stress xx");
	AutoGL_AddIntegerItem("Stress yy");
	AutoGL_AddIntegerItem("Stress xy");
}

static void setPanelOutput(void) {
	AutoGL_AddGroup("Output");

	AutoGL_AddCallback(saveViewImage, "saveViewImage");
	AutoGL_SetLabel("save view image");
	AutoGL_AddString(filename_view_image, "filename_view_image", 256);
	AutoGL_SetLabel("filename");
}

void AutoGL_SetUp (int argc, char *argv[]) {
	/* ビューを準備する。 */
	int i;
	if (argc != 3 && argc != 4) {
		printf("Error!!\n");
		printf("Usage:"
		       "./NURBS_view input.txt view.dat\n"
		       "\n");
		exit(1);
	}
	for (i = 1; i < argc; i++) {
		fp = fopen(argv[i], "r");
		if (fp == NULL) {
			printf("Error!!\n");
			printf("%s file not open!\n", argv[i]);
			exit(1);
		}
		fclose(fp);
	}
	printf("Start reading files\n");
	if (argc == 3) {
		argv[3] = NULL;
	}
	ReadFile(argv[1], argv[2], argv[3]);
	printf("End reading files\n\n");

	AutoGL_ResizeView(DEFAULT_WINDOW_SIZE, DEFAULT_WINDOW_SIZE);
	calclateMinMax();
	calculateViewPositionFromMesh();
	AutoGL_SetViewSize (view_size);
	AutoGL_SetViewCenter(view_center_x,
	                     view_center_y,
	                     0.0);
	AutoGL_SetViewDirection(view_direction_x,
	                        view_direction_y,
	                        view_direction_z);
	AutoGL_SetViewUpVector(view_up_vector_x,
	                       view_up_vector_y,
	                       view_up_vector_z);

	AutoGL_SetViewRedrawCallback (redrawView);
	AutoGL_SetDefaultCallbackInMode2D (NULL);	/* 二次元アプリ用イベント処理 */
	AutoGL_EnableDragInMode2D ();				/* ドラッグ有効 */
	AutoGL_SetBackgroundColor (1, 1, 1);			/* 背景の色 */

	AutoGL_SetDefaultCallbackInMode2D(handleEvent);
	AutoGL_EnableDragInMode2D();

	setPanelMain();
	//setPanelOtherModeAndMarks();
	setPanelScaleBar();
	setPanelSelectPatch ();
	//setPanelSection();
	//setPanelColorsAndSizes();
	if (fields_flag) {
		setPanelContourFieldOptions();
		setPanelContourFieldSelection();
	}
	//setPanelPoints();
	//setPanelViewPosition();
	setPanelOutput();
	//setPanelAnimation();
}
