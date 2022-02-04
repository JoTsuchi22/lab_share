#include <stdio.h>
#include <float.h>
#include <autogl.h>

#include "pick.h"

#include "constant.h"



/*
extern int num_element_nodes;
extern int num_element_triangles;
extern int num_element_quadrangles;
extern int **elements;
extern int num_triangles;
extern int num_quadrangles;
extern int **triangles;
extern int **quadrangles;
extern int num_visible_triangles;
extern int num_visible_quadrangles;
extern int *visible_triangles;
extern int *visible_quadrangles;
extern int num_crack_triangles;
extern int num_crack_quadrangles;
extern int *crack_triangles;
extern int *crack_quadrangles;
extern int num_segments;
extern int **segments;
extern int num_dimensions;
extern double **nodes;

extern int num_nodal_fields;
extern double **nodal_fields;
extern char **nodal_field_labels;
extern int num_elemental_fields;
extern double **elemental_fields;
extern char **elemental_field_labels;
extern int num_integer_elemental_fields;
extern int **integer_elemental_fields;
extern char **integer_elemental_field_labels;

extern double *displacement_fields[];

extern int num_points;
extern double **points;
extern int num_points_fields;
extern double **points_fields;

extern int mesh_flag;
extern int crack_mesh_flag;
extern int segments_flag;
extern int deformed_mesh_flag;
extern int deformed_crack_mesh_flag;
extern int deformed_segments_flag;
extern double deformation_scale_factor;
extern int points_type;
extern int marked_node;
extern int marked_element;
extern int marked_point;
*/
extern int point_size;

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

//static int doesHitNode(double *depth, int node_id, int displacement_flag);
//static int doesHitFace(double *depth, int face_id, int displacement_flag);
static int doesHitPoint(/*double *depth, */
    int patch_id, int point_id_x, int point_id_y);

#if 0
/*
 * Handle node pick event
 */
void handleNodePickEvent(void) {
	int picked_node = -1;
	double picked_node_depth = -DBL_MAX;
	double depth = 0.0;
	int i, j;

	if (mesh_flag) {
		for (i = 0; i < num_visible_triangles; i++) {
			int triangle_id = visible_triangles[i];

			for (j = 0; j < 6; j++) {
				int node_id = triangles[triangle_id][j];

				if (node_id == -1)
					continue;

				if (doesHitNode(&depth, node_id, 0)
				        && depth >= picked_node_depth) {
					picked_node = node_id;
					picked_node_depth = depth;
				}
			}
		}
		for (i = 0; i < num_visible_quadrangles; i++) {
			int quadrangle_id = visible_quadrangles[i];

			for (j = 0; j < 9; j++) {
				int node_id = quadrangles[quadrangle_id][j];

				if (node_id == -1)
					continue;

				if (doesHitNode(&depth, node_id, 0)
				        && depth >= picked_node_depth) {
					picked_node = node_id;
					picked_node_depth = depth;
				}
			}
		}
	}

	if (crack_mesh_flag) {
		for (i = 0; i < num_crack_triangles; i++) {
			int triangle_id = crack_triangles[i];

			for (j = 0; j < 6; j++) {
				int node_id = triangles[triangle_id][j];

				if (node_id == -1)
					continue;

				if (doesHitNode(&depth, node_id, 0)
				        && depth >= picked_node_depth) {
					picked_node = node_id;
					picked_node_depth = depth;
				}
			}
		}
		for (i = 0; i < num_crack_quadrangles; i++) {
			int quadrangle_id = crack_quadrangles[i];

			for (j = 0; j < 9; j++) {
				int node_id = quadrangles[quadrangle_id][j];

				if (node_id == -1)
					continue;

				if (doesHitNode(&depth, node_id, 0)
				        && depth >= picked_node_depth) {
					picked_node = node_id;
					picked_node_depth = depth;
				}
			}
		}
	}

	if (segments_flag) {
		for (i = 0; i < num_segments; i++)
			for (j = 0; j < 2; j++) {
				int node_id = segments[i][j];

				if (node_id == -1)
					continue;

				if (doesHitNode(&depth, node_id, 0)
				        && depth >= picked_node_depth) {
					picked_node = node_id;
					picked_node_depth = depth;
				}
			}
	}

	if (deformed_mesh_flag) {
		for (i = 0; i < num_visible_triangles; i++) {
			int triangle_id = visible_triangles[i];

			for (j = 0; j < 6; j++) {
				int node_id = triangles[triangle_id][j];

				if (node_id == -1)
					continue;

				if (doesHitNode(&depth, node_id, 1)
				        && depth >= picked_node_depth) {
					picked_node = node_id;
					picked_node_depth = depth;
				}
			}
		}
		for (i = 0; i < num_visible_quadrangles; i++) {
			int quadrangle_id = visible_quadrangles[i];

			for (j = 0; j < 9; j++) {
				int node_id = quadrangles[quadrangle_id][j];

				if (node_id == -1)
					continue;

				if (doesHitNode(&depth, node_id, 1)
				        && depth >= picked_node_depth) {
					picked_node = node_id;
					picked_node_depth = depth;
				}
			}
		}
	}

	if (deformed_crack_mesh_flag) {
		for (i = 0; i < num_crack_triangles; i++) {
			int triangle_id = crack_triangles[i];

			for (j = 0; j < 6; j++) {
				int node_id = triangles[triangle_id][j];

				if (node_id == -1)
					continue;

				if (doesHitNode(&depth, node_id, 1)
				        && depth >= picked_node_depth) {
					picked_node = node_id;
					picked_node_depth = depth;
				}
			}
		}
		for (i = 0; i < num_crack_quadrangles; i++) {
			int quadrangle_id = crack_quadrangles[i];

			for (j = 0; j < 9; j++) {
				int node_id = quadrangles[quadrangle_id][j];

				if (node_id == -1)
					continue;

				if (doesHitNode(&depth, node_id, 1)
				        && depth >= picked_node_depth) {
					picked_node = node_id;
					picked_node_depth = depth;
				}
			}
		}
	}

	if (deformed_segments_flag) {
		for (i = 0; i < num_segments; i++)
			for (j = 0; j < 2; j++) {
				int node_id = segments[i][j];

				if (node_id == -1)
					continue;

				if (doesHitNode(&depth, node_id, 1)
				        && depth >= picked_node_depth) {
					picked_node = node_id;
					picked_node_depth = depth;
				}
			}
	}

	if (picked_node != -1) {
		printf("Pick node %d\n", picked_node);

		printf("X: % .15E\n", nodes[picked_node][0]);
		printf("Y: % .15E\n", nodes[picked_node][1]);
		if (num_dimensions == 3)
			printf("Z: % .15E\n", nodes[picked_node][2]);

		for (i = 0; i < num_nodal_fields; i++)
			printf("%s: % .15E\n",
			       nodal_field_labels[i],
			       nodal_fields[i][picked_node]);
	}

	marked_node = picked_node;
	AutoGL_DrawView();
}

/*
 * Handle element pick event
 */
void handleElementPickEvent(void) {
	int picked_element = -1;
	int picked_element_face = -1;
	double picked_element_depth = -DBL_MAX;
	double depth = 0.0;
	int i;

	if (mesh_flag) {
		for (i = 0; i < num_visible_triangles; i++) {
			int triangle_id = visible_triangles[i];

			if (doesHitFace(&depth, triangle_id, 0)
			        && depth >= picked_element_depth) {
				picked_element       = triangle_id / num_element_triangles;
				picked_element_face  = triangle_id % num_element_triangles;
				picked_element_depth = depth;
			}
		}
		for (i = 0; i < num_visible_quadrangles; i++) {
			int quadrangle_id = visible_quadrangles[i];

			if (doesHitFace(&depth, quadrangle_id, 0)
			        && depth >= picked_element_depth) {
				picked_element       = quadrangle_id / num_element_quadrangles;
				picked_element_face  = quadrangle_id % num_element_quadrangles;
				picked_element_depth = depth;
			}
		}
	}

	if (crack_mesh_flag) {
		for (i = 0; i < num_crack_triangles; i++) {
			int triangle_id = crack_triangles[i];

			if (doesHitFace(&depth, triangle_id, 0)
			        && depth >= picked_element_depth) {
				picked_element       = triangle_id / num_element_triangles;
				picked_element_face  = triangle_id % num_element_triangles;
				picked_element_depth = depth;
			}
		}
		for (i = 0; i < num_crack_quadrangles; i++) {
			int quadrangle_id = crack_quadrangles[i];

			if (doesHitFace(&depth, quadrangle_id, 0)
			        && depth >= picked_element_depth) {
				picked_element       = quadrangle_id / num_element_quadrangles;
				picked_element_face  = quadrangle_id % num_element_quadrangles;
				picked_element_depth = depth;
			}
		}
	}
	if (deformed_mesh_flag) {
		for (i = 0; i < num_visible_triangles; i++) {
			int triangle_id = visible_triangles[i];

			if (doesHitFace(&depth, triangle_id, 1)
			        && depth >= picked_element_depth) {
				picked_element       = triangle_id / num_element_triangles;
				picked_element_face  = triangle_id % num_element_triangles;
				picked_element_depth = depth;
			}
		}
		for (i = 0; i < num_visible_quadrangles; i++) {
			int quadrangle_id = visible_quadrangles[i];

			if (doesHitFace(&depth, quadrangle_id, 1)
			        && depth >= picked_element_depth) {
				picked_element       = quadrangle_id / num_element_quadrangles;
				picked_element_face  = quadrangle_id % num_element_quadrangles;
				picked_element_depth = depth;
			}
		}
	}

	if (deformed_crack_mesh_flag) {
		for (i = 0; i < num_crack_triangles; i++) {
			int triangle_id = crack_triangles[i];

			if (doesHitFace(&depth, triangle_id, 1)
			        && depth >= picked_element_depth) {
				picked_element       = triangle_id / num_element_triangles;
				picked_element_face  = triangle_id % num_element_triangles;
				picked_element_depth = depth;
			}
		}
		for (i = 0; i < num_crack_quadrangles; i++) {
			int quadrangle_id = crack_quadrangles[i];

			if (doesHitFace(&depth, quadrangle_id, 1)
			        && depth >= picked_element_depth) {
				picked_element       = quadrangle_id / num_element_quadrangles;
				picked_element_face  = quadrangle_id % num_element_quadrangles;
				picked_element_depth = depth;
			}
		}
	}

	if (picked_element != -1 && picked_element_face != -1) {
		printf("Pick element %d (element face %d)\n",
		       picked_element, picked_element_face);

		printf("Element nodes:");
		for (i = 0; i < num_element_nodes; i++)
			printf(" %d", elements[picked_element][i]);
		printf("\n");

		for (i = 0; i < num_elemental_fields; i++)
			printf("%s: % .15E\n",
			       elemental_field_labels[i],
			       elemental_fields[i][picked_element]);

		for (i = 0; i < num_integer_elemental_fields; i++)
			printf("%s: %d\n",
			       integer_elemental_field_labels[i],
			       integer_elemental_fields[i][picked_element]);
	}

	marked_element = picked_element;
	AutoGL_DrawView();
}
#endif

void calculateStainStressIndexFromCoord(int patch_index, int index_xi, int index_eta,
                                        int *element_index_xi, int *element_index_eta,
                                        int *index_in_element_xi, int *index_in_element_eta,
                                        int *num_point_xi, int *num_point_eta) {
	int temp_num_point_xi = 1;
	int temp_num_point_eta = 1;
	int temp_element_n_xi = element_n_xi[patch_index];
	int temp_element_n_eta = element_n_eta[patch_index];
	int temp_index_in_element_xi, temp_index_in_element_eta;
	int temp_int1, temp_int2;

	temp_int1 = index_xi % division_ele_xi;
	temp_int2 = index_xi / division_ele_xi;
	if (temp_int1 == 0) {
		if (temp_int2 == 0) {
			(*element_index_xi) = 0;
			temp_index_in_element_xi = 0;
		} else if (temp_int2 == temp_element_n_xi) {
			(*element_index_xi) = temp_int2 - 1;
			temp_index_in_element_xi = division_ele_xi;
		} else {
			(*element_index_xi) = temp_int2 - 1;
			temp_num_point_xi = 2;
			temp_index_in_element_xi = division_ele_xi;
		}
	} else {
		(*element_index_xi) = temp_int2;
		temp_index_in_element_xi = temp_int1;
	}
	//printf("%d %d\n", temp_int1, temp_int2);

	temp_int1 = index_eta % division_ele_eta;
	temp_int2 = index_eta / division_ele_eta;
	if (temp_int1 == 0) {
		if (temp_int2 == 0) {
			(*element_index_eta) = 0;
			temp_index_in_element_eta = 0;
		} else if (temp_int2 == temp_element_n_eta) {
			(*element_index_eta) = temp_int2 - 1;
			temp_index_in_element_eta = division_ele_eta;
		} else {
			(*element_index_eta) = temp_int2 - 1;
			temp_num_point_eta = 2;
			temp_index_in_element_eta = division_ele_eta;
		}
	} else {
		(*element_index_eta) = temp_int2;
		temp_index_in_element_eta = temp_int1;
	}
	//printf("%d %d\n", temp_int1, temp_int2);

	//printf("%d %d\n", element_index_xi, (*element_index_eta));
	(*index_in_element_xi) = temp_index_in_element_xi;
	(*index_in_element_eta) = temp_index_in_element_eta;
	(*num_point_xi) = temp_num_point_xi;
	(*num_point_eta) = temp_num_point_eta;
}

/*
 * Handle point pick event
 */
void handlePointPickEvent(int patch_n, int *picked_point_flag) {
	/*
	int picked_point = -1;
	double picked_point_depth = -DBL_MAX;
	double depth = 0.0;
	int i;

		if (num_points == 0)
			return;

		if (points_type == NOTHING_POINTS_TYPE)
			return;

		for (i = 0; i < num_points; i++)
			if (doesHitPoint(&depth, i)
					&& depth >= picked_point_depth) {
				picked_point = i;
				picked_point_depth = depth;
			}

		if (picked_point != -1) {
			printf("Pick point %d\n", picked_point);

			printf("X: % .15E\n", points[picked_point][0]);
			printf("Y: % .15E\n", points[picked_point][1]);
			if (num_dimensions == 3)
				printf("Z: % .15E\n", points[picked_point][2]);

			for (i = 0; i < num_points_fields; i++)
				printf("Points field #%d: % .15E\n",
					   i, points_fields[picked_point][i]);
		}

		marked_point = picked_point;
		AutoGL_DrawView();
		*/

	int i, j, k, l, m;
	int jj, kk, ll, mm;
	int temp_element, temp_xi, temp_eta;
	int vis_point_num_xi, vis_point_num_eta;
	i = patch_n;
	for (j = 0; j < division_n_xi[i]; j++) {
		for (k = 0; k < division_n_eta[i]; k++) {
			if (doesHitPoint(i, j, k)) {
				printf("patch: %d\n", i);
				printf("  x: % 1.8e\n", coord_x[i][j][k]);
				printf("  y: % 1.8e\n", coord_y[i][j][k]);
				printf("  Displacement x: % 1.8e\n", disp_x[i][j][k]);
				printf("  Displacement y: % 1.8e\n", disp_y[i][j][k]);
				printf("  Displacement  : % 1.8e\n", disp_val[i][j][k]);
				//printf("%d %d %d\n", i, j, k);
				calculateStainStressIndexFromCoord(i, j, k,
				                                   &jj, &kk, &ll, &mm,
				                                   &vis_point_num_xi,
				                                   &vis_point_num_eta);
				//printf("%d %d %d %d %d   %d %d\n", i, jj, kk, ll, mm,
				//       vis_point_num_xi, vis_point_num_eta);
				for (l = 0; l < vis_point_num_xi; l++) {
					for (m = 0; m < vis_point_num_eta; m++) {
						temp_xi = ll;
						temp_eta = mm;
						temp_element = (jj + l) + (kk + m) * element_n_xi[i];
						if (vis_point_num_xi == 2) temp_xi = 0;
						if (vis_point_num_eta == 2) temp_eta = 0;
						printf("  Element: %d\n", temp_element);
						printf("    Strain xx: % 1.8e\n", strain_xx[i][jj + l][kk + m][temp_xi][temp_eta]);
						printf("    Strain yy: % 1.8e\n", strain_yy[i][jj + l][kk + m][temp_xi][temp_eta]);
						printf("    Strain xy: % 1.8e\n", strain_xy[i][jj + l][kk + m][temp_xi][temp_eta]);
						printf("    Stress xx: % 1.8e\n", stress_xx[i][jj + l][kk + m][temp_xi][temp_eta]);
						printf("    Stress yy: % 1.8e\n", stress_yy[i][jj + l][kk + m][temp_xi][temp_eta]);
						printf("    Stress xy: % 1.8e\n", stress_xy[i][jj + l][kk + m][temp_xi][temp_eta]);
					}
				}
				(*picked_point_flag) = 1;
				return;
			}
		}
	}
}

#if 0
/*
 * Does hit node?
 */
static int doesHitNode(double *depth, int node_id, int displacement_flag) {
	double x = nodes[node_id][0];
	double y = nodes[node_id][1];
	double z = nodes[node_id][2];

	if (displacement_flag) {
		x += deformation_scale_factor * displacement_fields[0][node_id];
		y += deformation_scale_factor * displacement_fields[1][node_id];
		z += deformation_scale_factor * displacement_fields[2][node_id];
	}

	if (num_dimensions == 3)
		return AutoGL_PointingDeviceIsHit3D(depth,
		                                    COORDS_SCALE_FACTOR * x,
		                                    COORDS_SCALE_FACTOR * y,
		                                    COORDS_SCALE_FACTOR * z,
		                                    10);
	else
		return AutoGL_PointingDeviceIsHit2D(COORDS_SCALE_FACTOR * x,
		                                    COORDS_SCALE_FACTOR * y,
		                                    10);
}

/*
 * Does hit face?
 */
static int doesHitFace(double *depth, int face_id, int displacement_flag) {
	double x, y, z;

	if (num_triangles) {
		x = (nodes[triangles[face_id][0]][0]
		     + nodes[triangles[face_id][1]][0]
		     + nodes[triangles[face_id][2]][0]) / 3.0;
		y = (nodes[triangles[face_id][0]][1]
		     + nodes[triangles[face_id][1]][1]
		     + nodes[triangles[face_id][2]][1]) / 3.0;
		z = (nodes[triangles[face_id][0]][2]
		     + nodes[triangles[face_id][1]][2]
		     + nodes[triangles[face_id][2]][2]) / 3.0;

		if (displacement_flag) {
			x += (displacement_fields[0][triangles[face_id][0]]
			      + displacement_fields[0][triangles[face_id][1]]
			      + displacement_fields[0][triangles[face_id][2]]) / 3.0;
			y += (displacement_fields[1][triangles[face_id][0]]
			      + displacement_fields[1][triangles[face_id][1]]
			      + displacement_fields[1][triangles[face_id][2]]) / 3.0;
			z += (displacement_fields[2][triangles[face_id][0]]
			      + displacement_fields[2][triangles[face_id][1]]
			      + displacement_fields[2][triangles[face_id][2]]) / 3.0;
		}
	} else if (num_quadrangles) {
		x = (nodes[quadrangles[face_id][0]][0]
		     + nodes[quadrangles[face_id][1]][0]
		     + nodes[quadrangles[face_id][2]][0]
		     + nodes[quadrangles[face_id][3]][0]) / 4.0;
		y = (nodes[quadrangles[face_id][0]][1]
		     + nodes[quadrangles[face_id][1]][1]
		     + nodes[quadrangles[face_id][2]][1]
		     + nodes[quadrangles[face_id][3]][1]) / 4.0;
		z = (nodes[quadrangles[face_id][0]][2]
		     + nodes[quadrangles[face_id][1]][2]
		     + nodes[quadrangles[face_id][2]][2]
		     + nodes[quadrangles[face_id][3]][2]) / 4.0;

		if (displacement_flag) {
			x += (displacement_fields[0][quadrangles[face_id][0]]
			      + displacement_fields[0][quadrangles[face_id][1]]
			      + displacement_fields[0][quadrangles[face_id][2]]
			      + displacement_fields[0][quadrangles[face_id][3]]) / 4.0;
			y += (displacement_fields[1][quadrangles[face_id][0]]
			      + displacement_fields[1][quadrangles[face_id][1]]
			      + displacement_fields[1][quadrangles[face_id][2]]
			      + displacement_fields[1][quadrangles[face_id][3]]) / 4.0;
			z += (displacement_fields[2][quadrangles[face_id][0]]
			      + displacement_fields[2][quadrangles[face_id][1]]
			      + displacement_fields[2][quadrangles[face_id][2]]
			      + displacement_fields[2][quadrangles[face_id][3]]) / 4.0;
		}
	} else {
		return 0;
	}

	if (num_dimensions == 3)
		return AutoGL_PointingDeviceIsHit3D(depth,
		                                    COORDS_SCALE_FACTOR * x,
		                                    COORDS_SCALE_FACTOR * y,
		                                    COORDS_SCALE_FACTOR * z,
		                                    10);
	else
		return AutoGL_PointingDeviceIsHit2D(COORDS_SCALE_FACTOR * x,
		                                    COORDS_SCALE_FACTOR * y,
		                                    10);
}
#endif

/*
 * Does hit point?
 */
static int doesHitPoint(/*double *depth, */
    int patch_id, int point_id_x, int point_id_y) {
	double x = coord_x[patch_id][point_id_x][point_id_y];
	double y = coord_y[patch_id][point_id_x][point_id_y];
	//double z = points[point_id][2];

	//if (num_dimensions == 3)
	//	return AutoGL_PointingDeviceIsHit3D(depth,
	//										COORDS_SCALE_FACTOR * x,
	//										COORDS_SCALE_FACTOR * y,
	//										COORDS_SCALE_FACTOR * z,
	//										point_size);
	//else
	return AutoGL_PointingDeviceIsHit2D(COORDS_SCALE_FACTOR * x,
	                                    COORDS_SCALE_FACTOR * y,
	                                    point_size);
}
