#include <stdio.h>
#include <string.h>
#include <math.h>
#include <autogl.h>

#include "plot.h"

#include "constant.h"
#include "draw.h"

extern int current_time_step;

extern double view_size;
extern double view_direction_x;
extern double view_direction_y;
extern double view_direction_z;
extern double view_up_vector_x;
extern double view_up_vector_y;
extern double view_up_vector_z;

extern int num_elements;
extern int num_element_triangles;
extern int num_element_quadrangles;
extern int **elements;
extern int **triangles;
extern int **quadrangles;
extern int num_segments;
extern int **segments;
extern int num_nodes;
extern int num_dimensions;
extern double **nodes;

extern double *contour_field;
extern int num_nodal_fields;
extern int num_elemental_fields;
extern int num_integer_elemental_fields;

extern double *displacement_fields[];

extern int num_points;
extern double **points;
extern int num_points_fields;
extern double **points_fields;

extern double deformation_scale_factor;
extern char title[];
extern int marked_node;
extern int marked_element;
extern int marked_point;
extern char scale_bar_format[];
extern double scale_bar_scale_factor;
extern char volume_color[];
extern char mesh_color[];
extern char point_color[];
extern int contour_type;
extern int reversed_contour_colors_flag;
extern int gray_contour_color_below_magenta_one_flag;
extern int num_contour_grades;
extern int contour_bar_flag;
extern char contour_bar_label[];
extern int auto_contour_bar_format_flag;
extern char contour_bar_format[];
extern double contour_bar_scale_factor;
extern double contour_field_max;
extern double contour_field_min;
extern int contour_field_id;
extern int points_reversed_contour_colors_flag;
extern int points_contour_bar_flag;
extern char points_contour_bar_label[];
extern int points_auto_contour_bar_format_flag;
extern char points_contour_bar_format[];
extern double points_contour_bar_scale_factor;
extern double points_field_max;
extern double points_field_min;
extern int points_field_id;

extern int format_flag_10;
extern int division_ele_xi;					//ξ方向の一要素あたりの分割数
extern int division_ele_eta;					//η方向の一要素あたりの分割数
extern int division_n_xi[MAX_PATCHES];		//ξ方向の表示する点の数
extern int division_n_eta[MAX_PATCHES];		//η方向の表示する点の数
extern int element_n_xi[MAX_PATCHES];		//ξ方向要素数
extern int element_n_eta[MAX_PATCHES];		//η方向要素数
extern int plot_type;

static void plotContourBar(int type,
                           double field_min, double field_max,
                           int num_grades,
                           char label[],
                           int auto_format_flag,
                           int format_flag_10,
                           char format[],
                           double scale_factor,
                           int reverse_flag,
                           int gray_flag);

static void plotCoordAxesAtCoords(double center_x, double center_y, double center_z);

static void setContourMap(int type,
                          int num_grades,
                          int reverse_flag,
                          int gray_flag);

/* メッシュを描画する */
void plotMesh (int size1,
               int division_n_xi, int division_n_eta,
               int element_n_xi, int element_n_eta,
               double coord_x[size1][size1],
               double coord_y[size1][size1]) {
	int i, j;
	double coord_old_x, coord_old_y, coord_new_x, coord_new_y;
	int temp_int;

	for (i = 0; i < element_n_xi + 1; i++) {
		temp_int = division_ele_xi * i;
		coord_old_x = coord_x[temp_int][0];
		coord_old_y = coord_y[temp_int][0];
		for (j = 1; j < division_n_eta; j++) {
			coord_new_x = coord_x[temp_int][j];
			coord_new_y = coord_y[temp_int][j];
			drawLine (coord_old_x, coord_old_y, 0,
			          coord_new_x, coord_new_y, 0);
			coord_old_x = coord_new_x;
			coord_old_y = coord_new_y;
		}
	}

	for (i = 0; i < element_n_eta + 1; i++) {
		temp_int = division_ele_eta * i;
		coord_old_x = coord_x[0][temp_int];
		coord_old_y = coord_y[0][temp_int];
		for (j = 1; j < division_n_xi; j++) {
			coord_new_x = coord_x[j][temp_int];
			coord_new_y = coord_y[j][temp_int];
			drawLine (coord_old_x, coord_old_y, 0,
			          coord_new_x, coord_new_y, 0);
			coord_old_x = coord_new_x;
			coord_old_y = coord_new_y;
		}
	}
}

/* 変形後のメッシュを描画する */
void plotDeformedMesh (int size1,
                       int division_n_xi, int division_n_eta,
                       int element_n_xi, int element_n_eta,
                       double coord_x[size1][size1],
                       double coord_y[size1][size1],
                       double disp_x[size1][size1],
                       double disp_y[size1][size1]) {
	int i, j;
	double coord_old_x, coord_old_y, coord_new_x, coord_new_y;
	int temp_int;

	for (i = 0; i < element_n_xi + 1; i++) {
		temp_int = division_ele_xi * i;
		coord_old_x = deformation_scale_factor * disp_x[temp_int][0]
		              + coord_x[temp_int][0];
		coord_old_y = deformation_scale_factor * disp_y[temp_int][0]
		              + coord_y[temp_int][0];
		for (j = 1; j < division_n_eta; j++) {
			coord_new_x = deformation_scale_factor * disp_x[temp_int][j]
			              + coord_x[temp_int][j];
			coord_new_y = deformation_scale_factor * disp_y[temp_int][j]
			              + coord_y[temp_int][j];
			drawLine (coord_old_x, coord_old_y, 0,
			          coord_new_x, coord_new_y, 0);
			coord_old_x = coord_new_x;
			coord_old_y = coord_new_y;
		}
	}

	for (i = 0; i < element_n_eta + 1; i++) {
		temp_int = division_ele_eta * i;
		coord_old_x = deformation_scale_factor * disp_x[0][temp_int]
		              + coord_x[0][temp_int];
		coord_old_y = deformation_scale_factor * disp_y[0][temp_int]
		              + coord_y[0][temp_int];
		for (j = 1; j < division_n_xi; j++) {
			coord_new_x = deformation_scale_factor * disp_x[j][temp_int]
			              + coord_x[j][temp_int];
			coord_new_y = deformation_scale_factor * disp_y[j][temp_int]
			              + coord_y[j][temp_int];
			drawLine (coord_old_x, coord_old_y, 0,
			          coord_new_x, coord_new_y, 0);
			coord_old_x = coord_new_x;
			coord_old_y = coord_new_y;
		}
	}
}

void plotVolume (int size1,
                 int division_n_xi, int division_n_eta,
                 double coord_x[size1][size1],
                 double coord_y[size1][size1],
                 double disp_x[size1][size1],
                 double disp_y[size1][size1]) {
	int i, j;
	double x1, x2, x3, x4, y1, y2, y3, y4;
	if (plot_type == VOLUME_PLOT_TYPE) {
		for (i = 0; i < division_n_xi - 1; i++) {
			for (j = 0; j < division_n_eta - 1; j++) {
				x1 = coord_x[i][j];
				x2 = coord_x[i + 1][j];
				x3 = coord_x[i + 1][j + 1];
				x4 = coord_x[i][j + 1];
				y1 = coord_y[i][j];
				y2 = coord_y[i + 1][j];
				y3 = coord_y[i + 1][j + 1];
				y4 = coord_y[i][j + 1];
				drawQuadrangle (x1, y1, 0,
				                x2, y2, 0,
				                x3, y3, 0,
				                x4, y4, 0);
			}
		}
	} else {
		for (i = 0; i < division_n_xi - 1; i++) {
			for (j = 0; j < division_n_eta - 1; j++) {
				x1 = coord_x[i][j]
				     + disp_x[i][j] * deformation_scale_factor;
				x2 = coord_x[i + 1][j]
				     + disp_x[i + 1][j] * deformation_scale_factor;
				x3 = coord_x[i + 1][j + 1]
				     + disp_x[i + 1][j + 1] * deformation_scale_factor;
				x4 = coord_x[i][j + 1]
				     + disp_x[i][j + 1] * deformation_scale_factor;
				y1 = coord_y[i][j]
				     + disp_y[i][j] * deformation_scale_factor;
				y2 = coord_y[i + 1][j]
				     + disp_y[i + 1][j] * deformation_scale_factor;
				y3 = coord_y[i + 1][j + 1]
				     + disp_y[i + 1][j + 1] * deformation_scale_factor;
				y4 = coord_y[i][j + 1]
				     + disp_y[i][j + 1] * deformation_scale_factor;
				drawQuadrangle (x1, y1, 0,
				                x2, y2, 0,
				                x3, y3, 0,
				                x4, y4, 0);
			}
		}
	}
}

/* 表面を描画する */
void plotEdge (int size1,
               int division_n_xi, int division_n_eta,
               int element_n_xi, int element_n_eta,
               double coord_x[size1][size1],
               double coord_y[size1][size1]) {
	int i;
	double coord_old_x, coord_old_y, coord_new_x, coord_new_y;
	int temp_int;

	coord_old_x = coord_x[0][0];
	coord_old_y = coord_y[0][0];
	for (i = 1; i < division_n_eta; i++) {
		coord_new_x = coord_x[0][i];
		coord_new_y = coord_y[0][i];
		drawLine (coord_old_x, coord_old_y, 0,
		          coord_new_x, coord_new_y, 0);
		coord_old_x = coord_new_x;
		coord_old_y = coord_new_y;
	}

	temp_int = division_ele_xi * element_n_xi;
	coord_old_x = coord_x[temp_int][0];
	coord_old_y = coord_y[temp_int][0];
	for (i = 1; i < division_n_eta; i++) {
		coord_new_x = coord_x[temp_int][i];
		coord_new_y = coord_y[temp_int][i];
		drawLine (coord_old_x, coord_old_y, 0,
		          coord_new_x, coord_new_y, 0);
		coord_old_x = coord_new_x;
		coord_old_y = coord_new_y;
	}

	coord_old_x = coord_x[0][0];
	coord_old_y = coord_y[0][0];
	for (i = 1; i < division_n_xi; i++) {
		coord_new_x = coord_x[i][0];
		coord_new_y = coord_y[i][0];
		drawLine (coord_old_x, coord_old_y, 0,
		          coord_new_x, coord_new_y, 0);
		coord_old_x = coord_new_x;
		coord_old_y = coord_new_y;
	}

	temp_int = division_ele_eta * element_n_eta;
	coord_old_x = coord_x[0][temp_int];
	coord_old_y = coord_y[0][temp_int];
	for (i = 1; i < division_n_xi; i++) {
		coord_new_x = coord_x[i][temp_int];
		coord_new_y = coord_y[i][temp_int];
		drawLine (coord_old_x, coord_old_y, 0,
		          coord_new_x, coord_new_y, 0);
		coord_old_x = coord_new_x;
		coord_old_y = coord_new_y;
	}
}

/* 変形後の表面を描画する */
void plotDeformedEdge (int size1,
                       int division_n_xi, int division_n_eta,
                       int element_n_xi, int element_n_eta,
                       double coord_x[size1][size1],
                       double coord_y[size1][size1],
                       double disp_x[size1][size1],
                       double disp_y[size1][size1]) {
	int i;
	double coord_old_x, coord_old_y, coord_new_x, coord_new_y;
	int temp_int;

	coord_old_x = coord_x[0][0]
	              + deformation_scale_factor * disp_x[0][0];
	coord_old_y = coord_y[0][0]
	              + deformation_scale_factor * disp_y[0][0];
	for (i = 1; i < division_n_eta; i++) {
		coord_new_x = coord_x[0][i]
		              + deformation_scale_factor * disp_x[0][i];
		coord_new_y = coord_y[0][i]
		              + deformation_scale_factor * disp_y[0][i];
		drawLine (coord_old_x, coord_old_y, 0,
		          coord_new_x, coord_new_y, 0);
		coord_old_x = coord_new_x;
		coord_old_y = coord_new_y;
	}

	temp_int = division_ele_xi * element_n_xi;
	coord_old_x = coord_x[temp_int][0]
	              + deformation_scale_factor * disp_x[temp_int][0];
	coord_old_y = coord_y[temp_int][0]
	              + deformation_scale_factor * disp_y[temp_int][0];
	for (i = 1; i < division_n_eta; i++) {
		coord_new_x = coord_x[temp_int][i]
		              + deformation_scale_factor * disp_x[temp_int][i];
		coord_new_y = coord_y[temp_int][i]
		              + deformation_scale_factor * disp_y[temp_int][i];
		drawLine (coord_old_x, coord_old_y, 0,
		          coord_new_x, coord_new_y, 0);
		coord_old_x = coord_new_x;
		coord_old_y = coord_new_y;
	}

	coord_old_x = coord_x[0][0]
	              + deformation_scale_factor * disp_x[0][0];
	coord_old_y = coord_y[0][0]
	              + deformation_scale_factor * disp_y[0][0];
	for (i = 1; i < division_n_xi; i++) {
		coord_new_x = coord_x[i][0]
		              + deformation_scale_factor * disp_x[i][0];
		coord_new_y = coord_y[i][0]
		              + deformation_scale_factor * disp_y[i][0];
		drawLine (coord_old_x, coord_old_y, 0,
		          coord_new_x, coord_new_y, 0);
		coord_old_x = coord_new_x;
		coord_old_y = coord_new_y;
	}

	temp_int = division_ele_eta * element_n_eta;
	coord_old_x = coord_x[0][temp_int]
	              + deformation_scale_factor * disp_x[0][temp_int];
	coord_old_y = coord_y[0][temp_int]
	              + deformation_scale_factor * disp_y[0][temp_int];
	for (i = 1; i < division_n_xi; i++) {
		coord_new_x = coord_x[i][temp_int];
		coord_new_y = coord_y[i][temp_int];
		coord_new_x = coord_x[i][temp_int]
		              + deformation_scale_factor * disp_x[i][temp_int];
		coord_new_y = coord_y[i][temp_int]
		              + deformation_scale_factor * disp_y[i][temp_int];
		drawLine (coord_old_x, coord_old_y, 0,
		          coord_new_x, coord_new_y, 0);
		coord_old_x = coord_new_x;
		coord_old_y = coord_new_y;
	}
}

/* 変位コンター図を描画する */
void plotDispContour (int size1,
                      int division_n_xi, int division_n_eta,
                      double coord_x[size1][size1],
                      double coord_y[size1][size1],
                      double disp_x[size1][size1],
                      double disp_y[size1][size1],
                      double value[size1][size1]) {
	int i, j;
	double x1, x2, x3, x4;
	double y1, y2, y3, y4;
	double scalar1, scalar2, scalar3, scalar4;

	setContourMap(contour_type,
	              num_contour_grades,
	              reversed_contour_colors_flag,
	              gray_contour_color_below_magenta_one_flag);

	plotContourBar(contour_type,
	               contour_field_min, contour_field_max,
	               num_contour_grades,
	               contour_bar_label,
	               auto_contour_bar_format_flag,
	               format_flag_10,
	               contour_bar_format,
	               contour_bar_scale_factor,
	               reversed_contour_colors_flag,
	               gray_contour_color_below_magenta_one_flag);

	if (plot_type == CONTOUR_VOLUME_PLOT_TYPE) {
		for (i = 0; i < division_n_xi - 1; i++) {
			for (j = 0; j < division_n_eta - 1; j++) {
				x1 = coord_x[i][j];
				x2 = coord_x[i + 1][j];
				x3 = coord_x[i + 1][j + 1];
				x4 = coord_x[i][j + 1];
				y1 = coord_y[i][j];
				y2 = coord_y[i + 1][j];
				y3 = coord_y[i + 1][j + 1];
				y4 = coord_y[i][j + 1];
				scalar1 = value[i][j];
				scalar2 = value[i + 1][j];
				scalar3 = value[i + 1][j + 1];
				scalar4 = value[i][j + 1];
				drawContourQuadrangle (x1, y1, 0, scalar1,
				                       x2, y2, 0, scalar2,
				                       x3, y3, 0, scalar3,
				                       x4, y4, 0, scalar4);
			}
		}
	} else {
		for (i = 0; i < division_n_xi - 1; i++) {
			for (j = 0; j < division_n_eta - 1; j++) {
				x1 = coord_x[i][j]
				     + deformation_scale_factor * disp_x[i][j];
				x2 = coord_x[i + 1][j]
				     + deformation_scale_factor * disp_x[i + 1][j];
				x3 = coord_x[i + 1][j + 1]
				     + deformation_scale_factor * disp_x[i + 1][j + 1];
				x4 = coord_x[i][j + 1]
				     + deformation_scale_factor * disp_x[i][j + 1];
				y1 = coord_y[i][j]
				     + deformation_scale_factor * disp_y[i][j];
				y2 = coord_y[i + 1][j]
				     + deformation_scale_factor * disp_y[i + 1][j];
				y3 = coord_y[i + 1][j + 1]
				     + deformation_scale_factor * disp_y[i + 1][j + 1];
				y4 = coord_y[i][j + 1]
				     + deformation_scale_factor * disp_y[i][j + 1];
				scalar1 = value[i][j];
				scalar2 = value[i + 1][j];
				scalar3 = value[i + 1][j + 1];
				scalar4 = value[i][j + 1];
				drawContourQuadrangle (x1, y1, 0, scalar1,
				                       x2, y2, 0, scalar2,
				                       x3, y3, 0, scalar3,
				                       x4, y4, 0, scalar4);
			}
		}
	}
}

/* 歪、応力コンター図を描画する */
void plotStrainStressContour (int size1, int size2, int size3,
                              int element_n_xi, int element_n_eta,
                              double coord_x[size1][size1],
                              double coord_y[size1][size1],
                              double disp_x[size1][size1],
                              double disp_y[size1][size1],
                              double value[size2][size2][size3 + 1][size3 + 1]) {
	int i, j, k, l;
	int jj, kk;
	int temp_index_xi, temp_index_eta;
	double x1, x2, x3, x4;
	double y1, y2, y3, y4;
	double scalar1, scalar2, scalar3, scalar4;

	setContourMap(contour_type,
	              num_contour_grades,
	              reversed_contour_colors_flag,
	              gray_contour_color_below_magenta_one_flag);

	plotContourBar(contour_type,
	               contour_field_min, contour_field_max,
	               num_contour_grades,
	               contour_bar_label,
	               auto_contour_bar_format_flag,
	               format_flag_10,
	               contour_bar_format,
	               contour_bar_scale_factor,
	               reversed_contour_colors_flag,
	               gray_contour_color_below_magenta_one_flag);

	if (plot_type == CONTOUR_VOLUME_PLOT_TYPE) {
		for (i = 0; i < element_n_xi; i++) {
			for (j = 0; j < element_n_eta; j++) {
				temp_index_xi = i * division_ele_xi;
				temp_index_eta = j * division_ele_eta;
				for (k = 0; k < division_ele_xi; k++) {
					for (l = 0; l < division_ele_eta; l++) {
						jj = temp_index_xi + k;
						kk = temp_index_eta + l;
						x1 = coord_x[jj][kk];
						x2 = coord_x[jj + 1][kk];
						x3 = coord_x[jj + 1][kk + 1];
						x4 = coord_x[jj][kk + 1];
						y1 = coord_y[jj][kk];
						y2 = coord_y[jj + 1][kk];
						y3 = coord_y[jj + 1][kk + 1];
						y4 = coord_y[jj][kk + 1];
						scalar1 = value[i][j][k][l];
						scalar2 = value[i][j][k + 1][l];
						scalar3 = value[i][j][k + 1][l + 1];
						scalar4 = value[i][j][k][l + 1];
						drawContourQuadrangle (x1, y1, 0, scalar1,
						                       x2, y2, 0, scalar2,
						                       x3, y3, 0, scalar3,
						                       x4, y4, 0, scalar4);
					}
				}
			}
		}
	} else {
		for (i = 0; i < element_n_xi; i++) {
			for (j = 0; j < element_n_eta; j++) {
				temp_index_xi = i * division_ele_xi;
				temp_index_eta = j * division_ele_eta;
				for (k = 0; k < division_ele_xi; k++) {
					for (l = 0; l < division_ele_eta; l++) {
						jj = temp_index_xi + k;
						kk = temp_index_eta + l;
						x1 = coord_x[jj][kk]
						     + deformation_scale_factor * disp_x[jj][kk];
						x2 = coord_x[jj + 1][kk]
						     + deformation_scale_factor * disp_x[jj + 1][kk];
						x3 = coord_x[jj + 1][kk + 1]
						     + deformation_scale_factor * disp_x[jj + 1][kk + 1];
						x4 = coord_x[jj][kk + 1]
						     + deformation_scale_factor * disp_x[jj][kk + 1];
						y1 = coord_y[jj][kk]
						     + deformation_scale_factor * disp_y[jj][kk];
						y2 = coord_y[jj + 1][kk]
						     + deformation_scale_factor * disp_y[jj + 1][kk];
						y3 = coord_y[jj + 1][kk + 1]
						     + deformation_scale_factor * disp_y[jj + 1][kk + 1];
						y4 = coord_y[jj][kk + 1]
						     + deformation_scale_factor * disp_y[jj][kk + 1];
						scalar1 = value[i][j][k][l];
						scalar2 = value[i][j][k + 1][l];
						scalar3 = value[i][j][k + 1][l + 1];
						scalar4 = value[i][j][k][l + 1];
						drawContourQuadrangle (x1, y1, 0, scalar1,
						                       x2, y2, 0, scalar2,
						                       x3, y3, 0, scalar3,
						                       x4, y4, 0, scalar4);
					}
				}
			}
		}
	}
}


/* コントロールネットを描画する */
void plotCntlNet (int cntl_p_n_xi, int cntl_p_n_eta,
                  double *cntl_px, double *cntl_py) {
	int i, j;
	double coord_old_x, coord_old_y, coord_new_x, coord_new_y;
	int temp = 0;

	for (i = 0; i < cntl_p_n_eta; i++) {
		coord_old_x = cntl_px[temp];
		coord_old_y = cntl_py[temp];
		temp += 1;
		for (j = 1; j < cntl_p_n_xi; j++) {
			coord_new_x = cntl_px[temp];
			coord_new_y = cntl_py[temp];
			drawLine (coord_old_x, coord_old_y, 0,
			          coord_new_x, coord_new_y, 0);
			coord_old_x = coord_new_x;
			coord_old_y = coord_new_y;
			temp += 1;
		}
	}

	for (i = 0; i < cntl_p_n_xi; i++) {
		temp = i;
		coord_old_x = cntl_px[temp];
		coord_old_y = cntl_py[temp];
		temp += cntl_p_n_xi;
		for (j = 1; j < cntl_p_n_eta; j++) {
			coord_new_x = cntl_px[temp];
			coord_new_y = cntl_py[temp];
			drawLine (coord_old_x, coord_old_y, 0,
			          coord_new_x, coord_new_y, 0);
			//printf("%d % 1.8e % 1.8e\t%d % 1.8e % 1.8e\n",
			//       temp - cntl_p_n_xi, coord_old_x, coord_old_y,
			//       temp, coord_new_x, coord_new_y);
			coord_old_x = coord_new_x;
			coord_old_y = coord_new_y;
			temp += cntl_p_n_xi;
		}
		//printf("\n");
	}
}

/* コントロールポイントを描画する */
void plotCntl_P (int cntl_p_n_xi, int cntl_p_n_eta,
                 double *cntl_px, double *cntl_py) {
	int i;
	for (i = 0; i < cntl_p_n_xi * cntl_p_n_eta; i++) {
		AutoGL_DrawDiamondMark2D(COORDS_SCALE_FACTOR * cntl_px[i],
		                         COORDS_SCALE_FACTOR * cntl_py[i], 0,
		                         5, 1);
		//printf("%d % 1.8e % 1.8e\n", i, cntl_px[i], cntl_py[i]);
	}
}


#if 1
/*
 * Plot contour bar
 */
static void plotContourBar(int type,
                           double field_min, double field_max,
                           int num_grades,
                           char label[],
                           int auto_format_flag,
                           int format_flag_10,
                           char format[],
                           double scale_factor,
                           int reverse_flag,
                           int gray_flag) {
	char buffer1[256], buffer2[256], buffer3[256];
	char format2[256];
	int string_grade_increment;
	int num_grades2;
	int width, height;
	double x, y, z;
	double x_lb, y_lb, z_lb;
	double x_rb, y_rb, z_rb;
	double x_lt, y_lt, z_lt;
	double x_rt, y_rt, z_rt;
	int i;

	double output_double, temp;
	int output_int;

	/* Set format */
	if (auto_format_flag) {
		double a = scale_factor * field_min;
		double b = scale_factor * field_max;

		if (format_flag_10 == 1)
			strcpy(format2, "%.2f");
		else if (a <= -1.0e+6 || b >= 1.0e+6)
			strcpy(format2, "%.2e");
		else if (b - a >= 1.0e+2 && (a <= -1.0e+2 || b >= 1.0e+2))
			strcpy(format2, "%.0f");
		else if (b - a >= 1.0e+1 && (a <= -1.0e+1 || b >= 1.0e+1))
			strcpy(format2, "%.1f");
		else if (b - a >= 1.0e+0 && (a <= -1.0e+0 || b >= 1.0e+0))
			strcpy(format2, "%.2f");
		else
			strcpy(format2, "%.2e");
	} else {
		strcpy(format2, format);
	}

	AutoGL_GetViewRangeDc(&width, &height);

	/* Draw label */
	AutoGL_SetColor(0.0, 0.0, 0.0);
	AutoGL_GetPositionOfDc(&x, &y, &z,
	                       width / 2 - 20 - 10 * strlen(label),
	                       height / 2 - 20 - 7,
	                       0);
	AutoGL_DrawString(x, y, z, label);

	if (format_flag_10 == 1) {
		height -= 7;
	}

	/* Draw strings */
	string_grade_increment = 1;
	while ((double)num_grades / (double)string_grade_increment > 10.0)
		string_grade_increment++;
	while (num_grades % string_grade_increment)
		string_grade_increment--;
	for (i = 0; i <= num_grades; i++) {
		if (i % string_grade_increment)
			continue;

		if (format_flag_10 == 0) {
			sprintf(buffer1, format2,
			        scale_factor
			        * ((double)(num_grades - i) * contour_field_min
			           + (double)i * contour_field_max) / (double)num_grades);
			AutoGL_GetPositionOfDc(
			    &x, &y, &z,
			    width / 2 - 50 - 10 * strlen(buffer1),
			    height / 2 - 240 - 7 + (int)((double)i / (double)num_grades * 200.0),
			    0);
			AutoGL_DrawString(x, y, z, buffer1);
		} else {
			temp = scale_factor
			       * ((double)(num_grades - i) * contour_field_min
			          + (double)i * contour_field_max) / (double)num_grades;
			if (temp == 0.0) {
				output_double = 0.0;
				output_int = 0.0;
			} else {
				output_int = floor(log10(fabs(temp)));
				output_double = temp * pow(10.0, -output_int);
			}
			//printf("%e % f*10^%d\n", temp, output_double, output_int);
			sprintf(buffer1, format2, output_double);
			sprintf(buffer2, " * 10");
			sprintf(buffer3, "% d", output_int);
			AutoGL_GetPositionOfDc(
			    &x, &y, &z,
			    width / 2 - 50 - 10 * (strlen(buffer1) + strlen(buffer2) + strlen(buffer3)),
			    height / 2 - 240 - 7 + (int)((double)i / (double)num_grades * 200.0),
			    0);
			AutoGL_DrawString(x, y, z, buffer1);

			AutoGL_GetPositionOfDc(
			    &x, &y, &z,
			    width / 2 - 50 - 10 * (strlen(buffer2) + strlen(buffer3)),
			    height / 2 - 240 - 7 + (int)((double)i / (double)num_grades * 200.0),
			    0);
			AutoGL_DrawString(x, y, z, buffer2);

			AutoGL_GetPositionOfDc(
			    &x, &y, &z,
			    width / 2 - 50 - 10 * strlen(buffer3),
			    height / 2 - 240 + (int)((double)i / (double)num_grades * 200.0),
			    0);
			AutoGL_DrawString(x, y, z, buffer3);
		}


	}

	/* Draw bar */
	AutoGL_SetLineWidth(2.0);
	AutoGL_SetColor(0.0, 0.0, 0.0);
	AutoGL_GetPositionOfDc(&x_lt, &y_lt, &z_lt,
	                       width / 2 - 40,
	                       height / 2 - 40,
	                       0);
	AutoGL_GetPositionOfDc(&x_rt, &y_rt, &z_rt,
	                       width / 2 - 20,
	                       height / 2 - 40,
	                       0);
	AutoGL_GetPositionOfDc(&x_lb, &y_lb, &z_lb,
	                       width / 2 - 40,
	                       height / 2 - 240,
	                       0);
	AutoGL_GetPositionOfDc(&x_rb, &y_rb, &z_rb,
	                       width / 2 - 20,
	                       height / 2 - 240,
	                       0);
	AutoGL_DrawLine(x_rt, y_rt, z_rt,
	                x_lt, y_lt, z_lt);
	AutoGL_DrawLine(x_lt, y_lt, z_lt,
	                x_lb, y_lb, z_lb);
	AutoGL_DrawLine(x_lb, y_lb, z_lb,
	                x_rb, y_rb, z_rb);
	AutoGL_DrawLine(x_rb, y_rb, z_rb,
	                x_rt, y_rt, z_rt);
	num_grades2 = (type == SMOOTH_CONTOUR_TYPE) ? 5 : num_grades;
	AutoGL_GetPositionOfDc(&x_lt, &y_lt, &z_lt,
	                       width / 2 - 40,
	                       height / 2 - 240,
	                       0);
	AutoGL_GetPositionOfDc(&x_rt, &y_rt, &z_rt,
	                       width / 2 - 20,
	                       height / 2 - 240,
	                       0);
	for (i = 0; i < num_grades2; i++) {
		double grade_min = (int)i / (double)num_grades2;
		double grade_max = (int)(i + 1) / (double)num_grades2;

		if (type == BAND_CONTOUR_TYPE)
			grade_max = grade_min;

		if (type == SMOOTH_CONTOUR_TYPE && gray_flag) {
			if (i == 0 && !reverse_flag)
				grade_min = EFFECTIVE_ZERO;
			else if (i == num_grades2 - 1 && reverse_flag)
				grade_max = 1.0 - EFFECTIVE_ZERO;
		}

		x_lb = x_lt; y_lb = y_lt; z_lb = z_lt;
		x_rb = x_rt; y_rb = y_rt; z_rb = z_rt;

		AutoGL_GetPositionOfDc(
		    &x_lt, &y_lt, &z_lt,
		    width / 2 - 40,
		    height / 2 - 240 + (int)((double)(i + 1) / (double)num_grades2 * 200.0),
		    0);
		AutoGL_GetPositionOfDc(
		    &x_rt, &y_rt, &z_rt,
		    width / 2 - 20,
		    height / 2 - 240 + (int)((double)(i + 1) / (double)num_grades2 * 200.0),
		    0);

		AutoGL_DrawSmoothContourQuadrangle(0.0, 1.0,
		                                   x_rt, y_rt, z_rt, grade_max,
		                                   x_lt, y_lt, z_lt, grade_max,
		                                   x_lb, y_lb, z_lb, grade_min,
		                                   x_rb, y_rb, z_rb, grade_min);
	}
}
#endif

#if 0
/*
 * Mark node
 */
void markNode(int deformation_flag) {
	double x, y, z;

	if (marked_node < 0 || marked_node >= num_nodes)
		return;

	if (deformation_flag && displacement_fields[0] == NULL)
		return;

	setColor(mesh_color);

	x = nodes[marked_node][0];
	y = nodes[marked_node][1];
	z = nodes[marked_node][2];

	if (deformation_flag) {
		x += deformation_scale_factor * displacement_fields[0][marked_node];
		y += deformation_scale_factor * displacement_fields[1][marked_node];
		z += deformation_scale_factor * displacement_fields[2][marked_node];
	}

	drawMark(x, y, z);
}

/*
 * Mark element
 */
void markElement(int deformation_flag) {
	int i;

	if (marked_element < 0 || marked_element >= num_elements)
		return;

	if (deformation_flag && displacement_fields[0] == NULL)
		return;

	setColor(mesh_color);

	/* Plot triangles */
	for (i = 0; i < num_element_triangles; i++) {
		int triangle_id = num_element_triangles * marked_element + i;

		int node0_id = triangles[triangle_id][0];
		int node1_id = triangles[triangle_id][1];
		int node2_id = triangles[triangle_id][2];

		double x0 = nodes[node0_id][0];
		double y0 = nodes[node0_id][1];
		double z0 = nodes[node0_id][2];
		double x1 = nodes[node1_id][0];
		double y1 = nodes[node1_id][1];
		double z1 = nodes[node1_id][2];
		double x2 = nodes[node2_id][0];
		double y2 = nodes[node2_id][1];
		double z2 = nodes[node2_id][2];

		if (deformation_flag) {
			x0 += deformation_scale_factor * displacement_fields[0][node0_id];
			y0 += deformation_scale_factor * displacement_fields[1][node0_id];
			z0 += deformation_scale_factor * displacement_fields[2][node0_id];
			x1 += deformation_scale_factor * displacement_fields[0][node1_id];
			y1 += deformation_scale_factor * displacement_fields[1][node1_id];
			z1 += deformation_scale_factor * displacement_fields[2][node1_id];
			x2 += deformation_scale_factor * displacement_fields[0][node2_id];
			y2 += deformation_scale_factor * displacement_fields[1][node2_id];
			z2 += deformation_scale_factor * displacement_fields[2][node2_id];
		}

		drawTriangle(x0, y0, z0,
		             x1, y1, z1,
		             x2, y2, z2);
	}

	/* Plot quadrangles */
	for (i = 0; i < num_element_quadrangles; i++) {
		int quadrangle_id = num_element_quadrangles * marked_element + i;

		int node0_id = quadrangles[quadrangle_id][0];
		int node1_id = quadrangles[quadrangle_id][1];
		int node2_id = quadrangles[quadrangle_id][2];
		int node3_id = quadrangles[quadrangle_id][3];

		double x0 = nodes[node0_id][0];
		double y0 = nodes[node0_id][1];
		double z0 = nodes[node0_id][2];
		double x1 = nodes[node1_id][0];
		double y1 = nodes[node1_id][1];
		double z1 = nodes[node1_id][2];
		double x2 = nodes[node2_id][0];
		double y2 = nodes[node2_id][1];
		double z2 = nodes[node2_id][2];
		double x3 = nodes[node3_id][0];
		double y3 = nodes[node3_id][1];
		double z3 = nodes[node3_id][2];

		if (deformation_flag) {
			x0 += deformation_scale_factor * displacement_fields[0][node0_id];
			y0 += deformation_scale_factor * displacement_fields[1][node0_id];
			z0 += deformation_scale_factor * displacement_fields[2][node0_id];
			x1 += deformation_scale_factor * displacement_fields[0][node1_id];
			y1 += deformation_scale_factor * displacement_fields[1][node1_id];
			z1 += deformation_scale_factor * displacement_fields[2][node1_id];
			x2 += deformation_scale_factor * displacement_fields[0][node2_id];
			y2 += deformation_scale_factor * displacement_fields[1][node2_id];
			z2 += deformation_scale_factor * displacement_fields[2][node2_id];
			x3 += deformation_scale_factor * displacement_fields[0][node3_id];
			y3 += deformation_scale_factor * displacement_fields[1][node3_id];
			z3 += deformation_scale_factor * displacement_fields[2][node3_id];
		}

		drawQuadrangle(x0, y0, z0,
		               x1, y1, z1,
		               x2, y2, z2,
		               x3, y3, z3);
	}
}
#endif

/*
 * Plot title
 */
void plotTitle(void) {
	char buffer[256];
	int width, height;
	double x, y, z;

	sprintf(buffer, title,
	        current_time_step, current_time_step, current_time_step);

	AutoGL_GetViewRangeDc(&width, &height);

	AutoGL_SetColor(0.0, 0.0, 0.0);
	AutoGL_GetPositionOfDc(&x, &y, &z,
	                       -width / 2 + 20,
	                       height / 2 - 20 - 7,
	                       0);
	AutoGL_DrawString(x, y, z, buffer);
}

/*
 * Plot coordinate axes
 */
void plotCoordAxes(void) {
	int width_dc, height_dc;
	int center_x_dc, center_y_dc;
	double center_x, center_y, center_z;

	AutoGL_GetViewRangeDc(&width_dc, &height_dc);

	center_x_dc = -width_dc / 2 + 50;
	center_y_dc = -height_dc / 2 + 50;

	AutoGL_GetPositionOfDc(&center_x, &center_y, &center_z,
	                       center_x_dc, center_y_dc, 0);

	plotCoordAxesAtCoords(center_x, center_y, center_z);
}

/*
 * Plot coordinate axes at origin
 */
void plotCoordAxesAtOrigin(void) {
	plotCoordAxesAtCoords(0.0, 0.0, 0.0);
}

/*
 * Plot coordinate axes at coordinates
 */
static void plotCoordAxesAtCoords(double center_x, double center_y, double center_z) {
	double length;
	double x, y, z;
	double sign_x, sign_y, sign_z;

	int i;

	double axis_length = 40 * AutoGL_GetViewSize() / AutoGL_GetViewSizeDc();

	AutoGL_SetLineWidth(2.0);
	AutoGL_SetColor(0.0, 0.0, 0.0);

	/* X arrow */
	if (fabs(fabs(view_direction_x) - 1.0) <= EFFECTIVE_ZERO) {
		AutoGL_SetColor(1.0, 1.0, 1.0);

		AutoGL_DrawCircle3D(center_x + 0.5 * axis_length, center_y, center_z,
		                    1.0, 0.0, 0.0,
		                    0.12 * axis_length, 2);
		AutoGL_DrawCircle3D(center_x - 0.5 * axis_length, center_y, center_z,
		                    1.0, 0.0, 0.0,
		                    0.12 * axis_length, 2);

		AutoGL_SetColor(0.0, 0.0, 0.0);

		for (i = 0; i < 8; i++)
			AutoGL_DrawLine(center_x + axis_length,
			                center_y + 0.12 * axis_length * cos(0.25 * M_PI * (double)i),
			                center_z + 0.12 * axis_length * sin(0.25 * M_PI * (double)i),
			                center_x + axis_length,
			                center_y + 0.12 * axis_length * cos(0.25 * M_PI * (double)(i + 1)),
			                center_z + 0.12 * axis_length * sin(0.25 * M_PI * (double)(i + 1)));
		for (i = 0; i < 8; i++)
			AutoGL_DrawLine(center_x - axis_length,
			                center_y + 0.12 * axis_length * cos(0.25 * M_PI * (double)i),
			                center_z + 0.12 * axis_length * sin(0.25 * M_PI * (double)i),
			                center_x - axis_length,
			                center_y + 0.12 * axis_length * cos(0.25 * M_PI * (double)(i + 1)),
			                center_z + 0.12 * axis_length * sin(0.25 * M_PI * (double)(i + 1)));

		AutoGL_DrawCircle3D(center_x + axis_length, center_y, center_z,
		                    1.0, 0.0, 0.0,
		                    0.05 * axis_length, 2);

		AutoGL_DrawLine(center_x - axis_length,
		                center_y + 0.12 * axis_length * cos(0.25 * M_PI),
		                center_z + 0.12 * axis_length * sin(0.25 * M_PI),
		                center_x - axis_length,
		                center_y + 0.12 * axis_length * cos(1.25 * M_PI),
		                center_z + 0.12 * axis_length * sin(1.25 * M_PI));
		AutoGL_DrawLine(center_x - axis_length,
		                center_y + 0.12 * axis_length * cos(0.75 * M_PI),
		                center_z + 0.12 * axis_length * sin(0.75 * M_PI),
		                center_x - axis_length,
		                center_y + 0.12 * axis_length * cos(1.75 * M_PI),
		                center_z + 0.12 * axis_length * sin(1.75 * M_PI));
	} else {
		AutoGL_DrawLine(center_x, center_y, center_z,
		                center_x + axis_length, center_y, center_z);
		length = sqrt(view_direction_y * view_direction_y
		              + view_direction_z * view_direction_z);
		AutoGL_DrawLine(center_x + axis_length,
		                center_y,
		                center_z,
		                center_x + 0.8 * axis_length,
		                center_y + 0.1 * axis_length * view_direction_z / length,
		                center_z - 0.1 * axis_length * view_direction_y / length);
		AutoGL_DrawLine(center_x + axis_length,
		                center_y,
		                center_z,
		                center_x + 0.8 * axis_length,
		                center_y - 0.1 * axis_length * view_direction_z / length,
		                center_z + 0.1 * axis_length * view_direction_y / length);
	}

	/* Y arrow */
	if (fabs(fabs(view_direction_y) - 1.0) <= EFFECTIVE_ZERO) {
		AutoGL_SetColor(1.0, 1.0, 1.0);

		AutoGL_DrawCircle3D(center_x, center_y + 0.5 * axis_length, center_z,
		                    0.0, 1.0, 0.0,
		                    0.12 * axis_length, 2);
		AutoGL_DrawCircle3D(center_x, center_y - 0.5 * axis_length, center_z,
		                    0.0, 1.0, 0.0,
		                    0.12 * axis_length, 2);

		AutoGL_SetColor(0.0, 0.0, 0.0);

		for (i = 0; i < 8; i++)
			AutoGL_DrawLine(center_x + 0.12 * axis_length * sin(0.25 * M_PI * (double)i),
			                center_y + axis_length,
			                center_z + 0.12 * axis_length * cos(0.25 * M_PI * (double)i),
			                center_x + 0.12 * axis_length * sin(0.25 * M_PI * (double)(i + 1)),
			                center_y + axis_length,
			                center_z + 0.12 * axis_length * cos(0.25 * M_PI * (double)(i + 1)));
		for (i = 0; i < 8; i++)
			AutoGL_DrawLine(center_x + 0.12 * axis_length * sin(0.25 * M_PI * (double)i),
			                center_y - axis_length,
			                center_z + 0.12 * axis_length * cos(0.25 * M_PI * (double)i),
			                center_x + 0.12 * axis_length * sin(0.25 * M_PI * (double)(i + 1)),
			                center_y - axis_length,
			                center_z + 0.12 * axis_length * cos(0.25 * M_PI * (double)(i + 1)));

		AutoGL_DrawCircle3D(center_x, center_y + axis_length, center_z,
		                    0.0, 1.0, 0.0,
		                    0.05 * axis_length, 2);

		AutoGL_DrawLine(center_x + 0.12 * axis_length * sin(0.25 * M_PI),
		                center_y - axis_length,
		                center_z + 0.12 * axis_length * cos(0.25 * M_PI),
		                center_x + 0.12 * axis_length * sin(1.25 * M_PI),
		                center_y - axis_length,
		                center_z + 0.12 * axis_length * cos(1.25 * M_PI));
		AutoGL_DrawLine(center_x + 0.12 * axis_length * sin(0.75 * M_PI),
		                center_y - axis_length,
		                center_z + 0.12 * axis_length * cos(0.75 * M_PI),
		                center_x + 0.12 * axis_length * sin(1.75 * M_PI),
		                center_y - axis_length,
		                center_z + 0.12 * axis_length * cos(1.75 * M_PI));
	} else {
		AutoGL_DrawLine(center_x, center_y, center_z,
		                center_x, center_y + axis_length, center_z);
		length = sqrt(view_direction_z * view_direction_z
		              + view_direction_x * view_direction_x);
		AutoGL_DrawLine(center_x,
		                center_y + axis_length,
		                center_z,
		                center_x - 0.1 * axis_length * view_direction_z / length,
		                center_y + 0.8 * axis_length,
		                center_z + 0.1 * axis_length * view_direction_x / length);
		AutoGL_DrawLine(center_x,
		                center_y + axis_length,
		                center_z,
		                center_x + 0.1 * axis_length * view_direction_z / length,
		                center_y + 0.8 * axis_length,
		                center_z - 0.1 * axis_length * view_direction_x / length);
	}

	/* Z arrow */
	if (num_dimensions == 3) {
		if (fabs(fabs(view_direction_z) - 1.0) <= EFFECTIVE_ZERO) {
			AutoGL_SetColor(1.0, 1.0, 1.0);

			AutoGL_DrawCircle3D(center_x, center_y, center_z + 0.5 * axis_length,
			                    0.0, 0.0, 1.0,
			                    0.12 * axis_length, 2);
			AutoGL_DrawCircle3D(center_x, center_y, center_z - 0.5 * axis_length,
			                    0.0, 0.0, 1.0,
			                    0.12 * axis_length, 2);

			AutoGL_SetColor(0.0, 0.0, 0.0);

			for (i = 0; i < 8; i++)
				AutoGL_DrawLine(center_x + 0.12 * axis_length * cos(0.25 * M_PI * (double)i),
				                center_y + 0.12 * axis_length * sin(0.25 * M_PI * (double)i),
				                center_z + axis_length,
				                center_x + 0.12 * axis_length * cos(0.25 * M_PI * (double)(i + 1)),
				                center_y + 0.12 * axis_length * sin(0.25 * M_PI * (double)(i + 1)),
				                center_z + axis_length);
			for (i = 0; i < 8; i++)
				AutoGL_DrawLine(center_x + 0.12 * axis_length * cos(0.25 * M_PI * (double)i),
				                center_y + 0.12 * axis_length * sin(0.25 * M_PI * (double)i),
				                center_z - axis_length,
				                center_x + 0.12 * axis_length * cos(0.25 * M_PI * (double)(i + 1)),
				                center_y + 0.12 * axis_length * sin(0.25 * M_PI * (double)(i + 1)),
				                center_z - axis_length);

			AutoGL_DrawCircle3D(center_x, center_y, center_z + axis_length,
			                    0.0, 0.0, 1.0,
			                    0.05 * axis_length, 2);

			AutoGL_DrawLine(center_x + 0.12 * axis_length * cos(0.25 * M_PI),
			                center_y + 0.12 * axis_length * sin(0.25 * M_PI),
			                center_z - axis_length,
			                center_x + 0.12 * axis_length * cos(1.25 * M_PI),
			                center_y + 0.12 * axis_length * sin(1.25 * M_PI),
			                center_z - axis_length);
			AutoGL_DrawLine(center_x + 0.12 * axis_length * cos(0.75 * M_PI),
			                center_y + 0.12 * axis_length * sin(0.75 * M_PI),
			                center_z - axis_length,
			                center_x + 0.12 * axis_length * cos(1.75 * M_PI),
			                center_y + 0.12 * axis_length * sin(1.75 * M_PI),
			                center_z - axis_length);
		} else {
			AutoGL_DrawLine(center_x, center_y, center_z,
			                center_x, center_y, center_z + axis_length);
			length = sqrt(view_direction_x * view_direction_x
			              + view_direction_y * view_direction_y);
			AutoGL_DrawLine(center_x,
			                center_y,
			                center_z + axis_length,
			                center_x + 0.1 * axis_length * view_direction_y / length,
			                center_y - 0.1 * axis_length * view_direction_x / length,
			                center_z + 0.8 * axis_length);
			AutoGL_DrawLine(center_x,
			                center_y,
			                center_z + axis_length,
			                center_x - 0.1 * axis_length * view_direction_y / length,
			                center_y + 0.1 * axis_length * view_direction_x / length,
			                center_z + 0.8 * axis_length);
		}

	}

	/* Strings */
	AutoGL_GetVectorProduct3D(&x, &y, &z,
	                          view_up_vector_x, view_up_vector_y, view_up_vector_z,
	                          view_direction_x, view_direction_y, view_direction_z);
	AutoGL_AddVector3D(&x, &y, &z,
	                   x, y, z,
	                   view_up_vector_x, view_up_vector_y, view_up_vector_z);
	sign_x = (x >= 0.0) ? 1.0 : -1.0;
	sign_y = (y >= 0.0) ? 1.0 : -1.0;
	sign_z = (z >= 0.0) ? 1.0 : -1.0;

	AutoGL_DrawString(center_x + 1.1 * axis_length,
	                  center_y + sign_y * 0.1 * axis_length,
	                  center_z + sign_z * 0.1 * axis_length,
	                  "X");
	AutoGL_DrawString(center_x + sign_x * 0.1 * axis_length,
	                  center_y + 1.1 * axis_length,
	                  center_z + sign_z * 0.1 * axis_length,
	                  "Y");
	if (num_dimensions == 3)
		AutoGL_DrawString(center_x + sign_x * 0.1 * axis_length,
		                  center_y + sign_y * 0.1 * axis_length,
		                  center_z + 1.1 * axis_length,
		                  "Z");
}

/*
 * Plot scale bar
 */
void plotScaleBar(void) {
	char buffer[256];
	double scale;
	int width, height, size;
	double x_lt, y_lt, z_lt;
	double x_lb, y_lb, z_lb;
	double x_rt, y_rt, z_rt;
	double x_rb, y_rb, z_rb;

	AutoGL_SetColor(0.0, 0.0, 0.0);

	AutoGL_GetViewRangeDc(&width, &height);
	size = AutoGL_GetViewSizeDc();

	/* Calculate scale */
	scale = pow(10.0, floor(log10(0.2 * view_size / COORDS_SCALE_FACTOR)));
	if (5.0 * scale <= 0.2 * view_size / COORDS_SCALE_FACTOR)
		scale *= 5.0;
	else if (2.0 * scale <= 0.2 * view_size / COORDS_SCALE_FACTOR)
		scale *= 2.0;

	/* Draw scale */
	sprintf(buffer, scale_bar_format,
	        scale_bar_scale_factor * scale);
	AutoGL_GetPositionOfDc(&x_lt, &y_lt, &z_lt,
	                       width / 2 - 20 - 10 * strlen(buffer),
	                       -height / 2 + 30 + 7,
	                       0);
	AutoGL_DrawString(x_lt, y_lt, z_lt, buffer);

	/* Draw bar */
	AutoGL_GetPositionOfDc(&x_lt, &y_lt, &z_lt,
	                       width / 2 - 20 - (int)(scale / (view_size / COORDS_SCALE_FACTOR) * (double)size),
	                       -height / 2 + 30,
	                       0);
	AutoGL_GetPositionOfDc(&x_lb, &y_lb, &z_lb,
	                       width / 2 - 20 - (int)(scale / (view_size / COORDS_SCALE_FACTOR) * (double)size),
	                       -height / 2 + 20,
	                       0);
	AutoGL_GetPositionOfDc(&x_rt, &y_rt, &z_rt,
	                       width / 2 - 20,
	                       -height / 2 + 30,
	                       0);
	AutoGL_GetPositionOfDc(&x_rb, &y_rb, &z_rb,
	                       width / 2 - 20,
	                       -height / 2 + 20,
	                       0);
	AutoGL_DrawQuadrangle(x_lt, y_lt, z_lt,
	                      x_lb, y_lb, z_lb,
	                      x_rb, y_rb, z_rb,
	                      x_rt, y_rt, z_rt);
}

/*
 * Set contour map
 */
static void setContourMap(int type,
                          int num_grades,
                          int reverse_flag,
                          int gray_flag) {
	double red, green, blue;

	getColorRgb(&red, &green, &blue, volume_color);

	AutoGL_ClearContourColor();

	if (type == BAND_CONTOUR_TYPE) {
		double min = 0.0;
		double max = 1.0 - 1.0 / (double)num_grades;

		if (reverse_flag) {
			if (gray_flag)
				max = 1.0 - 2.0 / (double)num_grades;

			AutoGL_AddContourColorOfGrade(0.0, 1.0, 0.0, 0.0);

			AutoGL_AddContourColorOfGrade(min, 1.0, 0.0, 0.0);
			if (max > min) {
				AutoGL_AddContourColorOfGrade(0.2 * max + 0.8 * min, 1.0, 1.0, 0.0);
				AutoGL_AddContourColorOfGrade(0.4 * max + 0.6 * min, 0.0, 1.0, 0.0);
				AutoGL_AddContourColorOfGrade(0.6 * max + 0.4 * min, 0.0, 1.0, 1.0);
				AutoGL_AddContourColorOfGrade(0.8 * max + 0.2 * min, 0.0, 0.0, 1.0);
				AutoGL_AddContourColorOfGrade(max,                   1.0, 0.0, 1.0);
			}

			if (gray_flag) {
				AutoGL_AddContourColorOfGrade(1.0 - 1.0 / (double)num_grades, red, green, blue);
				AutoGL_AddContourColorOfGrade(1.0, red, green, blue);
			} else {
				AutoGL_AddContourColorOfGrade(1.0, 1.0, 0.0, 1.0);
			}
		} else {
			if (gray_flag) {
				AutoGL_AddContourColorOfGrade(0.0, red, green, blue);

				min = 1.0 / (double)num_grades;
			} else {
				AutoGL_AddContourColorOfGrade(0.0, 1.0, 0.0, 1.0);
			}

			if (max > min) {
				AutoGL_AddContourColorOfGrade(min,                   1.0, 0.0, 1.0);
				AutoGL_AddContourColorOfGrade(0.2 * max + 0.8 * min, 0.0, 0.0, 1.0);
				AutoGL_AddContourColorOfGrade(0.4 * max + 0.6 * min, 0.0, 1.0, 1.0);
				AutoGL_AddContourColorOfGrade(0.6 * max + 0.4 * min, 0.0, 1.0, 0.0);
				AutoGL_AddContourColorOfGrade(0.8 * max + 0.2 * min, 1.0, 1.0, 0.0);
			}
			AutoGL_AddContourColorOfGrade(max, 1.0, 0.0, 0.0);

			AutoGL_AddContourColorOfGrade(1.0, 1.0, 0.0, 0.0);
		}
	} else {
		if (reverse_flag) {
			AutoGL_AddContourColorOfGrade(0.0, 1.0, 0.0, 0.0);
			AutoGL_AddContourColorOfGrade(0.2, 1.0, 1.0, 0.0);
			AutoGL_AddContourColorOfGrade(0.4, 0.0, 1.0, 0.0);
			AutoGL_AddContourColorOfGrade(0.6, 0.0, 1.0, 1.0);
			AutoGL_AddContourColorOfGrade(0.8, 0.0, 0.0, 1.0);

			if (gray_flag) {
				AutoGL_AddContourColorOfGrade(1.0 - EFFECTIVE_ZERO, 1.0, 0.0, 1.0);
				AutoGL_AddContourColorOfGrade(1.0, red, green, blue);
			} else {
				AutoGL_AddContourColorOfGrade(1.0, 1.0, 0.0, 1.0);
			}
		} else {
			if (gray_flag) {
				AutoGL_AddContourColorOfGrade(0.0, red, green, blue);
				AutoGL_AddContourColorOfGrade(EFFECTIVE_ZERO, 1.0, 0.0, 1.0);
			} else {
				AutoGL_AddContourColorOfGrade(0.0, 1.0, 0.0, 1.0);
			}

			AutoGL_AddContourColorOfGrade(0.2, 0.0, 0.0, 1.0);
			AutoGL_AddContourColorOfGrade(0.4, 0.0, 1.0, 1.0);
			AutoGL_AddContourColorOfGrade(0.6, 0.0, 1.0, 0.0);
			AutoGL_AddContourColorOfGrade(0.8, 1.0, 1.0, 0.0);
			AutoGL_AddContourColorOfGrade(1.0, 1.0, 0.0, 0.0);
		}
	}
}

/*
 * Set color
 */
void setColor(char color[]) {
	double red, green, blue;

	getColorRgb(&red, &green, &blue, color);

	AutoGL_SetColor(red, green, blue);
}

/*
 * Get color RGB
 */
void getColorRgb(double *r, double *g, double *b,
                 char color[]) {
	if (strcmp(color, "white") == 0) {
		*r = 1.0;
		*g = 1.0;
		*b = 1.0;
	} else if (strcmp(color, "whitesmoke") == 0) {
		*r = 245.0 / 255.0;
		*g = 245.0 / 255.0;
		*b = 245.0 / 255.0;
	} else if (strcmp(color, "lightgrey") == 0
	           || strcmp(color, "lightgray") == 0) {
		*r = 211.0 / 255.0;
		*g = 211.0 / 255.0;
		*b = 211.0 / 255.0;
	} else if (strcmp(color, "darkgray") == 0
	           || strcmp(color, "darkgrey") == 0) {
		*r = 169.0 / 255.0;
		*g = 169.0 / 255.0;
		*b = 169.0 / 255.0;
	} else if (strcmp(color, "dimgray") == 0
	           || strcmp(color, "dimgrey") == 0) {
		*r = 105.0 / 255.0;
		*g = 105.0 / 255.0;
		*b = 105.0 / 255.0;
	} else if (strcmp(color, "black") == 0) {
		*r = 0.0;
		*g = 0.0;
		*b = 0.0;
	} else if (strcmp(color, "red") == 0) {
		*r = 1.0;
		*g = 0.0;
		*b = 0.0;
	} else if (strcmp(color, "yellow") == 0) {
		*r = 1.0;
		*g = 1.0;
		*b = 0.0;
	} else if (strcmp(color, "green") == 0
	           || strcmp(color, "lime") == 0) {
		*r = 0.0;
		*g = 1.0;
		*b = 0.0;
	} else if (strcmp(color, "cyan") == 0
	           || strcmp(color, "aqua") == 0) {
		*r = 0.0;
		*g = 1.0;
		*b = 1.0;
	} else if (strcmp(color, "blue") == 0) {
		*r = 0.0;
		*g = 0.0;
		*b = 1.0;
	} else if (strcmp(color, "magenta") == 0
	           || strcmp(color, "fuchsia") == 0) {
		*r = 1.0;
		*g = 0.0;
		*b = 1.0;
	} else if (sscanf(color, "rgb(%lf,%lf,%lf)", r, g, b) != 3) {
		*r = 0.5;
		*g = 0.5;
		*b = 0.5;
	}
}