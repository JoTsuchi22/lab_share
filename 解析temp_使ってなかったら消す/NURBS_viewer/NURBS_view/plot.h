#ifndef PLOT_H
#define PLOT_H

/*
 * Get color RGB
 */
void getColorRgb(double *r, double *g, double *b,
				 char color[]);

void plotMesh (int size1,
			   int division_n_xi, int division_n_eta,
			   int element_n_xi, int element_n_eta,
			   double coord_x[size1][size1],
			   double coord_y[size1][size1]);

void plotDeformedMesh (int size1,
					   int division_n_xi, int division_n_eta,
					   int element_n_xi, int element_n_eta,
					   double coord_x[size1][size1],
					   double coord_y[size1][size1],
					   double disp_x[size1][size1],
					   double disp_y[size1][size1]);

void plotVolume (int size1,
				 int division_n_xi, int division_n_eta,
				 double coord_x[size1][size1],
				 double coord_y[size1][size1],
				 double disp_x[size1][size1],
				 double disp_y[size1][size1]);

void plotEdge (int size1,
			   int division_n_xi, int division_n_eta,
			   int element_n_xi, int element_n_eta,
			   double coord_x[size1][size1],
			   double coord_y[size1][size1]);

void plotDeformedEdge (int size1,
					   int division_n_xi, int division_n_eta,
					   int element_n_xi, int element_n_eta,
					   double coord_x[size1][size1],
					   double coord_y[size1][size1],
					   double disp_x[size1][size1],
					   double disp_y[size1][size1]);

void plotDispContour (int size1,
					  int division_n_xi, int division_n_eta,
					  double coord_x[size1][size1],
					  double coord_y[size1][size1],
					  double disp_x[size1][size1],
					  double disp_y[size1][size1],
					  double value[size1][size1]);

void plotStrainStressContour (int size1, int size2, int size3,
							  int element_n_xi, int element_n_eta,
							  double coord_x[size1][size1],
							  double coord_y[size1][size1],
							  double disp_x[size1][size1],
							  double disp_y[size1][size1],
							  double value[size2][size2][size3 + 1][size3 + 1]);

void plotCntlNet (int cntl_p_n_xi, int cntl_p_n_eta,
				  double *cntl_px, double *cntl_py);

void plotCntl_P (int cntl_p_n_xi, int cntl_p_n_eta,
				 double *cntl_px, double *cntl_py);

void plotTitle(void);

void plotCoordAxes(void);

void plotCoordAxesAtOrigin(void);

void plotScaleBar(void);

#endif  /* PLOT_H */
