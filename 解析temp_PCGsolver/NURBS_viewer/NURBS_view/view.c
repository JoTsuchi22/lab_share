/* -*- mode: C; c-file-style: "bsd"; c-basic-offset: 4; indent-tabs-mode: nil; -*- vim: set filetype=c shiftwidth=4 softtabstop=4 expandtab: */
/*
 * Zarusoba Viewer
 *
 * (c) Yusa, Yasunori
 */

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <autogl.h>

#include "view.h"

#include "constant.h"

extern double view_size;
extern double view_center_x;
extern double view_center_y;
extern double view_center_z;
extern double view_direction_x;
extern double view_direction_y;
extern double view_direction_z;
extern double view_up_vector_x;
extern double view_up_vector_y;
extern double view_up_vector_z;
extern double section_a;
extern double section_b;
extern double section_c;
extern double section_d;

//extern int num_dimensions;

extern int patch_n;
extern int division_ele_xi;
extern int division_ele_eta;
extern int division_n_xi[MAX_PATCHES];
extern int division_n_eta[MAX_PATCHES];
extern int element_n_xi[MAX_PATCHES];
extern int element_n_eta[MAX_PATCHES];

extern double coord_x[MAX_PATCHES][MAX_POINTS][MAX_POINTS];
extern double coord_y[MAX_PATCHES][MAX_POINTS][MAX_POINTS];

extern double disp_x[MAX_PATCHES][MAX_POINTS][MAX_POINTS];
extern double disp_y[MAX_PATCHES][MAX_POINTS][MAX_POINTS];
extern double disp_val[MAX_PATCHES][MAX_POINTS][MAX_POINTS];

extern double strain_xx[MAX_PATCHES][MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];
extern double strain_yy[MAX_PATCHES][MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];
extern double strain_xy[MAX_PATCHES][MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];

extern double stress_xx[MAX_PATCHES][MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];
extern double stress_yy[MAX_PATCHES][MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];
extern double stress_xy[MAX_PATCHES][MAX_ELEMENTS][MAX_ELEMENTS][MAX_DIVISION + 1][MAX_DIVISION + 1];

extern double max_x;
extern double min_x;
extern double max_y;
extern double min_y;

extern double MaxRange_disp_x;
extern double MinRange_disp_x;
extern double MaxRange_disp_y;
extern double MinRange_disp_y;
extern double MaxRange_disp_val;
extern double MinRange_disp_val;

extern double MaxRange_strain_xx;
extern double MinRange_strain_xx;
extern double MaxRange_strain_yy;
extern double MinRange_strain_yy;
extern double MaxRange_strain_xy;
extern double MinRange_strain_xy;

extern double MaxRange_stress_xx;
extern double MinRange_stress_xx;
extern double MaxRange_stress_yy;
extern double MinRange_stress_yy;
extern double MaxRange_stress_xy;
extern double MinRange_stress_xy;

extern int fields_flag;

static void moveViewCenter(double distance);

void calclateMinMax (void) {
	int i, j, k, l, m;

	if (fields_flag) {
		for (i = 0; i < patch_n; i++) {
			for (j = 0; j < division_n_xi[i]; j++) {
				for (k = 0; k < division_n_eta[i]; k++) {
					if (max_x < coord_x[i][j][k]) max_x = coord_x[i][j][k];
					if (min_x > coord_x[i][j][k]) min_x = coord_x[i][j][k];
					if (max_y < coord_y[i][j][k]) max_y = coord_y[i][j][k];
					if (min_y > coord_y[i][j][k]) min_y = coord_y[i][j][k];

					if (MaxRange_disp_x < disp_x[i][j][k]) {
						MaxRange_disp_x = disp_x[i][j][k];
					}
					if (MinRange_disp_x > disp_x[i][j][k]) {
						MinRange_disp_x = disp_x[i][j][k];
					}
					if (MaxRange_disp_y < disp_y[i][j][k]) {
						MaxRange_disp_y = disp_y[i][j][k];
					}
					if (MinRange_disp_y > disp_y[i][j][k]) {
						MinRange_disp_y = disp_y[i][j][k];
					}
					if (MaxRange_disp_val < disp_val[i][j][k]) {
						MaxRange_disp_val = disp_val[i][j][k];
					}
					if (MinRange_disp_val > disp_val[i][j][k]) {
						MinRange_disp_val = disp_val[i][j][k];
					}
				}
			}

			for (j = 0; j < element_n_xi[i]; j++) {
				for (k = 0; k < element_n_eta[i]; k++) {
					for (l = 0; l < division_ele_xi + 1; l++) {
						for (m = 0; m < division_ele_eta + 1; m++) {
							if (MaxRange_strain_xx < strain_xx[i][j][k][l][m]) {
								MaxRange_strain_xx = strain_xx[i][j][k][l][m];
							}
							if (MinRange_strain_xx > strain_xx[i][j][k][l][m]) {
								MinRange_strain_xx = strain_xx[i][j][k][l][m];
							}
							if (MaxRange_strain_yy < strain_yy[i][j][k][l][m]) {
								MaxRange_strain_yy = strain_yy[i][j][k][l][m];
							}
							if (MinRange_strain_yy > strain_yy[i][j][k][l][m]) {
								MinRange_strain_yy = strain_yy[i][j][k][l][m];
							}
							if (MaxRange_strain_xy < strain_xy[i][j][k][l][m]) {
								MaxRange_strain_xy = strain_xy[i][j][k][l][m];
							}
							if (MinRange_strain_xy > strain_xy[i][j][k][l][m]) {
								MinRange_strain_xy = strain_xy[i][j][k][l][m];
							}

							if (MaxRange_stress_xx < stress_xx[i][j][k][l][m]) {
								MaxRange_stress_xx = stress_xx[i][j][k][l][m];
							}
							if (MinRange_stress_xx > stress_xx[i][j][k][l][m]) {
								MinRange_stress_xx = stress_xx[i][j][k][l][m];
							}
							if (MaxRange_stress_yy < stress_yy[i][j][k][l][m]) {
								MaxRange_stress_yy = stress_yy[i][j][k][l][m];
							}
							if (MinRange_stress_yy > stress_yy[i][j][k][l][m]) {
								MinRange_stress_yy = stress_yy[i][j][k][l][m];
							}
							if (MaxRange_stress_xy < stress_xy[i][j][k][l][m]) {
								MaxRange_stress_xy = stress_xy[i][j][k][l][m];
							}
							if (MinRange_stress_xy > stress_xy[i][j][k][l][m]) {
								MinRange_stress_xy = stress_xy[i][j][k][l][m];
							}
						}
					}
				}
			}
		}
	} else {
		for (i = 0; i < patch_n; i++) {
			for (j = 0; j < division_n_xi[i]; j++) {
				for (k = 0; k < division_n_eta[i]; k++) {
					if (max_x < coord_x[i][j][k]) max_x = coord_x[i][j][k];
					if (min_x > coord_x[i][j][k]) min_x = coord_x[i][j][k];
					if (max_y < coord_y[i][j][k]) max_y = coord_y[i][j][k];
					if (min_y > coord_y[i][j][k]) min_y = coord_y[i][j][k];
				}
			}
		}
	}
}

/*
 * Calculate view position from mesh
 */
void calculateViewPositionFromMesh(void) {
	view_size
	    = COORDS_SCALE_FACTOR
	      * 1.2                   /* For margin */
	      * 0.5 * AutoGL_GetVectorDistance3D(min_x, min_y, 0,
	              max_x, max_y, 0);

	view_center_x = COORDS_SCALE_FACTOR * 0.5 * (max_x + min_x);
	view_center_y = COORDS_SCALE_FACTOR * 0.5 * (max_y + min_y);

	view_direction_x = 0.0;
	view_direction_y = 0.0;
	view_direction_z = 1.0;

	view_up_vector_x = 0.0;
	view_up_vector_y = 1.0;
	view_up_vector_z = 0.0;

	section_a = 0.0;
	section_b = 1.0;
	section_c = 0.0;
	section_d = -0.5 * (max_y + min_y);
}

/*
 * Align view direction
 */
void alignViewDirection(void) {
	double x = fabs(view_direction_x);
	double y = fabs(view_direction_y);
	double z = fabs(view_direction_z);

	if (fabs(view_up_vector_x) >= 1.0) x = 0.0;
	if (fabs(view_up_vector_y) >= 1.0) y = 0.0;
	if (fabs(view_up_vector_z) >= 1.0) z = 0.0;

	if (x > y && x > z) {
		view_direction_x = (view_direction_x >= 0.0) ? 1.0 : -1.0;
		view_direction_y = 0.0;
		view_direction_z = 0.0;
	} else if (y > z) {
		view_direction_x = 0.0;
		view_direction_y = (view_direction_y >= 0.0) ? 1.0 : -1.0;
		view_direction_z = 0.0;
	} else {
		view_direction_x = 0.0;
		view_direction_y = 0.0;
		view_direction_z = (view_direction_z >= 0.0) ? 1.0 : -1.0;
	}

	AutoGL_SetViewDirection(view_direction_x,
	                        view_direction_y,
	                        view_direction_z);
	AutoGL_DrawView();
}

/*
 * Align view up vector
 */
void alignViewUpVector(void) {
	double x = fabs(view_up_vector_x);
	double y = fabs(view_up_vector_y);
	double z = fabs(view_up_vector_z);

	if (fabs(view_direction_x) >= 1.0) x = 0.0;
	if (fabs(view_direction_y) >= 1.0) y = 0.0;
	if (fabs(view_direction_z) >= 1.0) z = 0.0;

	if (x > y && x > z) {
		view_up_vector_x = (view_up_vector_x >= 0.0) ? 1.0 : -1.0;
		view_up_vector_y = 0.0;
		view_up_vector_z = 0.0;
	} else if (y > z) {
		view_up_vector_x = 0.0;
		view_up_vector_y = (view_up_vector_y >= 0.0) ? 1.0 : -1.0;
		view_up_vector_z = 0.0;
	} else {
		view_up_vector_x = 0.0;
		view_up_vector_y = 0.0;
		view_up_vector_z = (view_up_vector_z >= 0.0) ? 1.0 : -1.0;
	}

	AutoGL_SetViewUpVector(view_up_vector_x,
	                       view_up_vector_y,
	                       view_up_vector_z);
	AutoGL_DrawView();
}

/*
 * Step forward
 */
void stepForward(void) {
	moveViewCenter(-VIEW_CENTER_MOVING_DISTANCE);
}

/*
 * Step backward
 */
void stepBackward(void) {
	moveViewCenter(VIEW_CENTER_MOVING_DISTANCE);
}

/*
 * Move view center
 */
static void moveViewCenter(double distance) {
	view_center_x += distance * view_size * view_direction_x;
	view_center_y += distance * view_size * view_direction_y;
	view_center_z += distance * view_size * view_direction_z;

	AutoGL_SetViewCenter(view_center_x,
	                     view_center_y,
	                     view_center_z);
	AutoGL_DrawView();
}

/*
 * Print view position
 */
void printViewPosition(void) {
	printf("size: % .15E\n", view_size);
	printf("center x: % .15E\n", view_center_x);
	printf("center y: % .15E\n", view_center_y);
	/*
	if (num_dimensions == 3) {
		printf("center z: % .15E\n", view_center_z);
		printf("dir x: % .15E\n", view_direction_x);
		printf("dir y: % .15E\n", view_direction_y);
		printf("dir z: % .15E\n", view_direction_z);
		printf("up x: % .15E\n", view_up_vector_x);
		printf("up y: % .15E\n", view_up_vector_y);
		printf("up z: % .15E\n", view_up_vector_z);
	}
	*/
}
