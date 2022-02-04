#include <stdio.h>
#include <float.h>
#include <autogl.h>
#include <math.h>

#include "clickcalc.h"
#include "constant.h"
#include "NURBS_calc.h"

extern int fields_flag;
extern double E, nu;

double TriSignedArea(double x0, double y0,
                     double x1, double y1,
                     double x2, double y2) {
	double output = 0.5 * ((x1 - x0) * (y2 - y0)
	                       - (x2 - x0) * (y1 - y0));
	return output;
}

int SegmentHit(double p0_x0, double p0_y0,
               double p0_x1, double p0_y1,
               double p1_x0, double p1_y0,
               double p1_x1, double p1_y1) {
	double s0, s1, s2, s3;
	s0 = TriSignedArea( p0_x0, p0_y0,
	                    p0_x1, p0_y1,
	                    p1_x0, p1_y0);
	s1 = TriSignedArea( p0_x0, p0_y0,
	                    p0_x1, p0_y1,
	                    p1_x1, p1_y1);
	s2 = TriSignedArea( p1_x0, p1_y0,
	                    p1_x1, p1_y1,
	                    p0_x0, p0_y0);
	s3 = TriSignedArea( p1_x0, p1_y0,
	                    p1_x1, p1_y1,
	                    p0_x1, p0_y1);
	if (s0 * s1 < 0 && s2 * s3 < 0) {
		return 1;
	} else {
		return 0;
	}
}

int CalcXiEtaByNR(double px, double py,
                  double *knot_vec_xi, double *knot_vec_eta,
                  double *cntl_px, double *cntl_py,
                  double *disp_cntl_px, double *disp_cntl_py,
                  int cntl_p_n_xi, int cntl_p_n_eta,
                  double *weight, int order_xi, int order_eta,
                  double *output_xi, double *output_eta) {
	double temp_xi, temp_eta;
	double temp_x, temp_y;
	double temp_matrix[2][2];
	double temp_dxi, temp_deta;
	double temp_tol_x = DBL_MAX;
	double temp_tol_y = DBL_MAX;

	(*output_xi) = 0;
	(*output_eta) = 0;

	int i;
	double tol = 10e-8;

	temp_xi = knot_vec_xi[0] + knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = knot_vec_eta[0] + knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;
	//printf("% 1.8e % 1.8e\n", temp_xi, temp_eta);
	for (i = 0; i < 1000; i++) {
		rNURBS_surface(knot_vec_xi, knot_vec_eta,
		               cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
		               weight, order_xi, order_eta,
		               temp_xi, temp_eta,
		               &temp_x, &temp_y,
		               &temp_matrix[0][0], &temp_matrix[0][1],
		               &temp_matrix[1][0], &temp_matrix[1][1]);

		temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;
		
		//収束した場合////////////////////////////////////////////////////////////////
		if (temp_tol_x < tol && temp_tol_y < tol) {
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;

			int knot_n_xi = cntl_p_n_xi + order_xi + 1;
			int knot_n_eta = cntl_p_n_eta + order_eta + 1;
			double dtilda_xi = 0.0;
			double dtilda_eta = 0.0;
			double disp_x, disp_y;
			double dxi_x, deta_x, dxi_y, deta_y;
			double dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y;
			double temp_matrix2[2][2];
			double temp;
			double strain_xx, strain_yy, strain_xy;
			double stress_xx, stress_yy, stress_xy;

			for (i = 0; i < knot_n_xi; i++) {
				if ( knot_vec_xi[i] < temp_xi && temp_xi < knot_vec_xi[i + 1]) {
					dtilda_xi = ( knot_vec_xi[i + 1] - knot_vec_xi[i] ) / 2.0;
					//printf("%f\n", dtilda_xi);
					break;
				}
			}
			for (i = 0; i < knot_n_eta; i++) {
				if ( knot_vec_eta[i] < temp_eta && temp_eta < knot_vec_eta[i + 1]) {
					dtilda_eta = ( knot_vec_eta[i + 1] - knot_vec_eta[i] ) / 2.0;
					//printf("%f\n", dtilda_eta);
					break;
				}
			}

			rNURBS_surface(knot_vec_xi, knot_vec_eta,
			               cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
			               weight, order_xi, order_eta,
			               temp_xi, temp_eta,
			               &temp, &temp,
			               &dxi_x, &deta_x,
			               &dxi_y, &deta_y);
			//printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
			//	   dxi_x, deta_x, dxi_y, deta_y);

			rNURBS_surface(knot_vec_xi, knot_vec_eta,
			               disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
			               weight, order_xi, order_eta,
			               temp_xi, temp_eta,
			               &disp_x, &disp_y,
			               &dxi_disp_x, &deta_disp_x,
			               &dxi_disp_y, &deta_disp_y);
			//printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
			//	   dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y);

			temp_matrix2[0][0] = dxi_x * dtilda_xi;
			temp_matrix2[0][1] = dxi_y * dtilda_xi;
			temp_matrix2[1][0] = deta_x * dtilda_eta;
			temp_matrix2[1][1] = deta_y * dtilda_eta;

			InverseMatrix_2D(temp_matrix2);

			strain_xx = temp_matrix2[0][0] * dtilda_xi
			            * dxi_disp_x
			            + temp_matrix2[0][1] * dtilda_eta
			            * deta_disp_x;
			strain_yy = temp_matrix2[1][0] * dtilda_xi
			            * dxi_disp_y
			            + temp_matrix2[1][1] * dtilda_eta
			            * deta_disp_y;
			strain_xy = temp_matrix2[1][0] * dtilda_xi
			            * dxi_disp_x
			            + temp_matrix2[1][1] * dtilda_eta
			            * deta_disp_x
			            + temp_matrix2[0][0] * dtilda_xi
			            * dxi_disp_y
			            + temp_matrix2[0][1] * dtilda_eta
			            * deta_disp_y;

			double D_matrix[3][3] = {{0.0}};
			int stress_type_flag = 1;
			if (stress_type_flag == 0) { //平面応力状態
				temp = E * (1.0 - nu * nu);
				D_matrix[0][0] = temp;
				D_matrix[0][1] = nu * temp;
				D_matrix[1][0] = nu * temp;
				D_matrix[1][1] = temp;
				D_matrix[2][2] = (1.0 - nu) / 2.0 * temp;
			} else if (stress_type_flag == 1) { //平面ひずみ状態(2Dの場合はこっち)
				temp = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
				D_matrix[0][0] = temp;
				D_matrix[0][1] = nu / (1.0 - nu) * temp;
				D_matrix[1][0] = nu / (1.0 - nu) * temp;
				D_matrix[1][1] = temp;
				D_matrix[2][2] = (1.0 - 2.0 * nu) / 2.0 / (1.0 - nu) * temp;
			}

			stress_xx = D_matrix[0][0] * strain_xx
			            + D_matrix[0][1] * strain_yy;
			stress_yy = D_matrix[1][0] * strain_xx
			            + D_matrix[1][1] * strain_yy;
			stress_xy = D_matrix[2][2] * strain_xy;

			temp = sqrt(disp_x * disp_x + disp_y * disp_y);
			printf("x:   % 1.8e\n", px);
			printf("y:   % 1.8e\n", py);
			printf("xi:  % 1.8e\n", temp_xi);
			printf("eta: % 1.8e\n", temp_eta);
			printf("Displacement x: % 1.8e\n", disp_x);
			printf("Displacement y: % 1.8e\n", disp_y);
			printf("Displacement  : % 1.8e\n", temp);
			printf("Strain xx: % 1.8e\n", strain_xx);
			printf("Strain yy: % 1.8e\n", strain_yy);
			printf("Strain xy: % 1.8e\n", strain_xy);
			printf("Stress xx: % 1.8e\n", stress_xx);
			printf("Stress yy: % 1.8e\n", stress_yy);
			printf("Stress xy: % 1.8e\n", stress_xy);

			return i;
		}

		InverseMatrix_2D(temp_matrix);

		temp_dxi = temp_matrix[0][0] * (px - temp_x)
		           + temp_matrix[0][1] * (py - temp_y);
		temp_deta = temp_matrix[1][0] * (px - temp_x)
		            + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_deta;
		if (temp_xi < knot_vec_xi[0])
			temp_xi = knot_vec_xi[0];
		if (temp_xi > knot_vec_xi[cntl_p_n_xi + order_xi])
			temp_xi = knot_vec_xi[cntl_p_n_xi + order_xi];
		if (temp_eta < knot_vec_eta[0])
			temp_eta = knot_vec_eta[0];
		if (temp_eta > knot_vec_eta[cntl_p_n_eta + order_eta])
			temp_eta = knot_vec_eta[cntl_p_n_eta + order_eta];

		//temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
		//printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
	}
	temp_xi = knot_vec_xi[0] + knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = knot_vec_eta[0] + knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;
	for (i = 0; i < 1000; i++) {
		lNURBS_surface(knot_vec_xi, knot_vec_eta,
		               cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
		               weight, order_xi, order_eta,
		               temp_xi, temp_eta,
		               &temp_x, &temp_y,
		               &temp_matrix[0][0], &temp_matrix[0][1],
		               &temp_matrix[1][0], &temp_matrix[1][1]);

		temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;
		
		//収束した場合////////////////////////////////////////////////////////////////
		if (temp_tol_x < tol && temp_tol_y < tol) {
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;

			int knot_n_xi = cntl_p_n_xi + order_xi + 1;
			int knot_n_eta = cntl_p_n_eta + order_eta + 1;
			double dtilda_xi = 0.0;
			double dtilda_eta = 0.0;
			double disp_x, disp_y;
			double dxi_x, deta_x, dxi_y, deta_y;
			double dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y;
			double temp_matrix2[2][2];
			double temp;
			double strain_xx, strain_yy, strain_xy;
			double stress_xx, stress_yy, stress_xy;

			for (i = 0; i < knot_n_xi; i++) {
				if ( knot_vec_xi[i] < temp_xi && temp_xi < knot_vec_xi[i + 1]) {
					dtilda_xi = ( knot_vec_xi[i + 1] - knot_vec_xi[i] ) / 2.0;
					//printf("%f\n", dtilda_xi);
					break;
				}
			}
			for (i = 0; i < knot_n_eta; i++) {
				if ( knot_vec_eta[i] < temp_eta && temp_eta < knot_vec_eta[i + 1]) {
					dtilda_eta = ( knot_vec_eta[i + 1] - knot_vec_eta[i] ) / 2.0;
					//printf("%f\n", dtilda_eta);
					break;
				}
			}

			lNURBS_surface(knot_vec_xi, knot_vec_eta,
			               cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
			               weight, order_xi, order_eta,
			               temp_xi, temp_eta,
			               &temp, &temp,
			               &dxi_x, &deta_x,
			               &dxi_y, &deta_y);
			//printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
			//	   dxi_x, deta_x, dxi_y, deta_y);

			lNURBS_surface(knot_vec_xi, knot_vec_eta,
			               disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
			               weight, order_xi, order_eta,
			               temp_xi, temp_eta,
			               &disp_x, &disp_y,
			               &dxi_disp_x, &deta_disp_x,
			               &dxi_disp_y, &deta_disp_y);
			//printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
			//	   dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y);

			temp_matrix2[0][0] = dxi_x * dtilda_xi;
			temp_matrix2[0][1] = dxi_y * dtilda_xi;
			temp_matrix2[1][0] = deta_x * dtilda_eta;
			temp_matrix2[1][1] = deta_y * dtilda_eta;

			InverseMatrix_2D(temp_matrix2);

			strain_xx = temp_matrix2[0][0] * dtilda_xi
			            * dxi_disp_x
			            + temp_matrix2[0][1] * dtilda_eta
			            * deta_disp_x;
			strain_yy = temp_matrix2[1][0] * dtilda_xi
			            * dxi_disp_y
			            + temp_matrix2[1][1] * dtilda_eta
			            * deta_disp_y;
			strain_xy = temp_matrix2[1][0] * dtilda_xi
			            * dxi_disp_x
			            + temp_matrix2[1][1] * dtilda_eta
			            * deta_disp_x
			            + temp_matrix2[0][0] * dtilda_xi
			            * dxi_disp_y
			            + temp_matrix2[0][1] * dtilda_eta
			            * deta_disp_y;

			double D_matrix[3][3] = {{0.0}};
			int stress_type_flag = 1;
			if (stress_type_flag == 0) { //平面応力状態
				temp = E * (1.0 - nu * nu);
				D_matrix[0][0] = temp;
				D_matrix[0][1] = nu * temp;
				D_matrix[1][0] = nu * temp;
				D_matrix[1][1] = temp;
				D_matrix[2][2] = (1.0 - nu) / 2.0 * temp;
			} else if (stress_type_flag == 1) { //平面ひずみ状態(2Dの場合はこっち)
				temp = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
				D_matrix[0][0] = temp;
				D_matrix[0][1] = nu / (1.0 - nu) * temp;
				D_matrix[1][0] = nu / (1.0 - nu) * temp;
				D_matrix[1][1] = temp;
				D_matrix[2][2] = (1.0 - 2.0 * nu) / 2.0 / (1.0 - nu) * temp;
			}

			stress_xx = D_matrix[0][0] * strain_xx
			            + D_matrix[0][1] * strain_yy;
			stress_yy = D_matrix[1][0] * strain_xx
			            + D_matrix[1][1] * strain_yy;
			stress_xy = D_matrix[2][2] * strain_xy;

			temp = sqrt(disp_x * disp_x + disp_y * disp_y);
			printf("x:   % 1.8e\n", px);
			printf("y:   % 1.8e\n", py);
			printf("xi:  % 1.8e\n", temp_xi);
			printf("eta: % 1.8e\n", temp_eta);
			printf("Displacement x: % 1.8e\n", disp_x);
			printf("Displacement y: % 1.8e\n", disp_y);
			printf("Displacement  : % 1.8e\n", temp);
			printf("Strain xx: % 1.8e\n", strain_xx);
			printf("Strain yy: % 1.8e\n", strain_yy);
			printf("Strain xy: % 1.8e\n", strain_xy);
			printf("Stress xx: % 1.8e\n", stress_xx);
			printf("Stress yy: % 1.8e\n", stress_yy);
			printf("Stress xy: % 1.8e\n", stress_xy);

			return i;
		}

		InverseMatrix_2D(temp_matrix);

		temp_dxi = temp_matrix[0][0] * (px - temp_x)
		           + temp_matrix[0][1] * (py - temp_y);
		temp_deta = temp_matrix[1][0] * (px - temp_x)
		            + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_deta;
		if (temp_xi < knot_vec_xi[0])
			temp_xi = knot_vec_xi[0];
		if (temp_xi > knot_vec_xi[cntl_p_n_xi + order_xi])
			temp_xi = knot_vec_xi[cntl_p_n_xi + order_xi];
		if (temp_eta < knot_vec_eta[0])
			temp_eta = knot_vec_eta[0];
		if (temp_eta > knot_vec_eta[cntl_p_n_eta + order_eta])
			temp_eta = knot_vec_eta[cntl_p_n_eta + order_eta];

		//temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
		//printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
	}
	temp_xi = knot_vec_xi[0] + knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = knot_vec_eta[0] + knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;
	for (i = 0; i < 1000; i++) {
		rlNURBS_surface(knot_vec_xi, knot_vec_eta,
		               cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
		               weight, order_xi, order_eta,
		               temp_xi, temp_eta,
		               &temp_x, &temp_y,
		               &temp_matrix[0][0], &temp_matrix[0][1],
		               &temp_matrix[1][0], &temp_matrix[1][1]);

		temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;
		
		//収束した場合////////////////////////////////////////////////////////////////
		if (temp_tol_x < tol && temp_tol_y < tol) {
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;

			int knot_n_xi = cntl_p_n_xi + order_xi + 1;
			int knot_n_eta = cntl_p_n_eta + order_eta + 1;
			double dtilda_xi = 0.0;
			double dtilda_eta = 0.0;
			double disp_x, disp_y;
			double dxi_x, deta_x, dxi_y, deta_y;
			double dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y;
			double temp_matrix2[2][2];
			double temp;
			double strain_xx, strain_yy, strain_xy;
			double stress_xx, stress_yy, stress_xy;

			for (i = 0; i < knot_n_xi; i++) {
				if ( knot_vec_xi[i] < temp_xi && temp_xi < knot_vec_xi[i + 1]) {
					dtilda_xi = ( knot_vec_xi[i + 1] - knot_vec_xi[i] ) / 2.0;
					//printf("%f\n", dtilda_xi);
					break;
				}
			}
			for (i = 0; i < knot_n_eta; i++) {
				if ( knot_vec_eta[i] < temp_eta && temp_eta < knot_vec_eta[i + 1]) {
					dtilda_eta = ( knot_vec_eta[i + 1] - knot_vec_eta[i] ) / 2.0;
					//printf("%f\n", dtilda_eta);
					break;
				}
			}

			rlNURBS_surface(knot_vec_xi, knot_vec_eta,
			               cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
			               weight, order_xi, order_eta,
			               temp_xi, temp_eta,
			               &temp, &temp,
			               &dxi_x, &deta_x,
			               &dxi_y, &deta_y);
			//printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
			//	   dxi_x, deta_x, dxi_y, deta_y);

			rlNURBS_surface(knot_vec_xi, knot_vec_eta,
			               disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
			               weight, order_xi, order_eta,
			               temp_xi, temp_eta,
			               &disp_x, &disp_y,
			               &dxi_disp_x, &deta_disp_x,
			               &dxi_disp_y, &deta_disp_y);
			//printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
			//	   dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y);

			temp_matrix2[0][0] = dxi_x * dtilda_xi;
			temp_matrix2[0][1] = dxi_y * dtilda_xi;
			temp_matrix2[1][0] = deta_x * dtilda_eta;
			temp_matrix2[1][1] = deta_y * dtilda_eta;

			InverseMatrix_2D(temp_matrix2);

			strain_xx = temp_matrix2[0][0] * dtilda_xi
			            * dxi_disp_x
			            + temp_matrix2[0][1] * dtilda_eta
			            * deta_disp_x;
			strain_yy = temp_matrix2[1][0] * dtilda_xi
			            * dxi_disp_y
			            + temp_matrix2[1][1] * dtilda_eta
			            * deta_disp_y;
			strain_xy = temp_matrix2[1][0] * dtilda_xi
			            * dxi_disp_x
			            + temp_matrix2[1][1] * dtilda_eta
			            * deta_disp_x
			            + temp_matrix2[0][0] * dtilda_xi
			            * dxi_disp_y
			            + temp_matrix2[0][1] * dtilda_eta
			            * deta_disp_y;

			double D_matrix[3][3] = {{0.0}};
			int stress_type_flag = 1;
			if (stress_type_flag == 0) { //平面応力状態
				temp = E * (1.0 - nu * nu);
				D_matrix[0][0] = temp;
				D_matrix[0][1] = nu * temp;
				D_matrix[1][0] = nu * temp;
				D_matrix[1][1] = temp;
				D_matrix[2][2] = (1.0 - nu) / 2.0 * temp;
			} else if (stress_type_flag == 1) { //平面ひずみ状態(2Dの場合はこっち)
				temp = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
				D_matrix[0][0] = temp;
				D_matrix[0][1] = nu / (1.0 - nu) * temp;
				D_matrix[1][0] = nu / (1.0 - nu) * temp;
				D_matrix[1][1] = temp;
				D_matrix[2][2] = (1.0 - 2.0 * nu) / 2.0 / (1.0 - nu) * temp;
			}

			stress_xx = D_matrix[0][0] * strain_xx
			            + D_matrix[0][1] * strain_yy;
			stress_yy = D_matrix[1][0] * strain_xx
			            + D_matrix[1][1] * strain_yy;
			stress_xy = D_matrix[2][2] * strain_xy;

			temp = sqrt(disp_x * disp_x + disp_y * disp_y);
			printf("x:   % 1.8e\n", px);
			printf("y:   % 1.8e\n", py);
			printf("xi:  % 1.8e\n", temp_xi);
			printf("eta: % 1.8e\n", temp_eta);
			printf("Displacement x: % 1.8e\n", disp_x);
			printf("Displacement y: % 1.8e\n", disp_y);
			printf("Displacement  : % 1.8e\n", temp);
			printf("Strain xx: % 1.8e\n", strain_xx);
			printf("Strain yy: % 1.8e\n", strain_yy);
			printf("Strain xy: % 1.8e\n", strain_xy);
			printf("Stress xx: % 1.8e\n", stress_xx);
			printf("Stress yy: % 1.8e\n", stress_yy);
			printf("Stress xy: % 1.8e\n", stress_xy);

			return i;
		}

		InverseMatrix_2D(temp_matrix);

		temp_dxi = temp_matrix[0][0] * (px - temp_x)
		           + temp_matrix[0][1] * (py - temp_y);
		temp_deta = temp_matrix[1][0] * (px - temp_x)
		            + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_deta;
		if (temp_xi < knot_vec_xi[0])
			temp_xi = knot_vec_xi[0];
		if (temp_xi > knot_vec_xi[cntl_p_n_xi + order_xi])
			temp_xi = knot_vec_xi[cntl_p_n_xi + order_xi];
		if (temp_eta < knot_vec_eta[0])
			temp_eta = knot_vec_eta[0];
		if (temp_eta > knot_vec_eta[cntl_p_n_eta + order_eta])
			temp_eta = knot_vec_eta[cntl_p_n_eta + order_eta];

		//temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
		//printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
	}
	temp_xi = knot_vec_xi[0] + knot_vec_xi[cntl_p_n_xi + order_xi];
	temp_xi *= 0.5;
	temp_eta = knot_vec_eta[0] + knot_vec_eta[cntl_p_n_eta + order_eta];
	temp_eta *= 0.5;
	for (i = 0; i < 1000; i++) {
		lrNURBS_surface(knot_vec_xi, knot_vec_eta,
		               cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
		               weight, order_xi, order_eta,
		               temp_xi, temp_eta,
		               &temp_x, &temp_y,
		               &temp_matrix[0][0], &temp_matrix[0][1],
		               &temp_matrix[1][0], &temp_matrix[1][1]);

		temp_tol_x = px - temp_x;
		temp_tol_x *= temp_tol_x;
		temp_tol_y = py - temp_y;
		temp_tol_y *= temp_tol_y;
		
		//収束した場合////////////////////////////////////////////////////////////////
		if (temp_tol_x < tol && temp_tol_y < tol) {
			(*output_xi) = temp_xi;
			(*output_eta) = temp_eta;

			int knot_n_xi = cntl_p_n_xi + order_xi + 1;
			int knot_n_eta = cntl_p_n_eta + order_eta + 1;
			double dtilda_xi = 0.0;
			double dtilda_eta = 0.0;
			double disp_x, disp_y;
			double dxi_x, deta_x, dxi_y, deta_y;
			double dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y;
			double temp_matrix2[2][2];
			double temp;
			double strain_xx, strain_yy, strain_xy;
			double stress_xx, stress_yy, stress_xy;

			for (i = 0; i < knot_n_xi; i++) {
				if ( knot_vec_xi[i] < temp_xi && temp_xi < knot_vec_xi[i + 1]) {
					dtilda_xi = ( knot_vec_xi[i + 1] - knot_vec_xi[i] ) / 2.0;
					//printf("%f\n", dtilda_xi);
					break;
				}
			}
			for (i = 0; i < knot_n_eta; i++) {
				if ( knot_vec_eta[i] < temp_eta && temp_eta < knot_vec_eta[i + 1]) {
					dtilda_eta = ( knot_vec_eta[i + 1] - knot_vec_eta[i] ) / 2.0;
					//printf("%f\n", dtilda_eta);
					break;
				}
			}

			lrNURBS_surface(knot_vec_xi, knot_vec_eta,
			               cntl_px, cntl_py, cntl_p_n_xi, cntl_p_n_eta,
			               weight, order_xi, order_eta,
			               temp_xi, temp_eta,
			               &temp, &temp,
			               &dxi_x, &deta_x,
			               &dxi_y, &deta_y);
			//printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
			//	   dxi_x, deta_x, dxi_y, deta_y);

			lrNURBS_surface(knot_vec_xi, knot_vec_eta,
			               disp_cntl_px, disp_cntl_py, cntl_p_n_xi, cntl_p_n_eta,
			               weight, order_xi, order_eta,
			               temp_xi, temp_eta,
			               &disp_x, &disp_y,
			               &dxi_disp_x, &deta_disp_x,
			               &dxi_disp_y, &deta_disp_y);
			//printf("% 1.4e % 1.4e % 1.4e % 1.4e\n",
			//	   dxi_disp_x, deta_disp_x, dxi_disp_y, deta_disp_y);

			temp_matrix2[0][0] = dxi_x * dtilda_xi;
			temp_matrix2[0][1] = dxi_y * dtilda_xi;
			temp_matrix2[1][0] = deta_x * dtilda_eta;
			temp_matrix2[1][1] = deta_y * dtilda_eta;

			InverseMatrix_2D(temp_matrix2);

			strain_xx = temp_matrix2[0][0] * dtilda_xi
			            * dxi_disp_x
			            + temp_matrix2[0][1] * dtilda_eta
			            * deta_disp_x;
			strain_yy = temp_matrix2[1][0] * dtilda_xi
			            * dxi_disp_y
			            + temp_matrix2[1][1] * dtilda_eta
			            * deta_disp_y;
			strain_xy = temp_matrix2[1][0] * dtilda_xi
			            * dxi_disp_x
			            + temp_matrix2[1][1] * dtilda_eta
			            * deta_disp_x
			            + temp_matrix2[0][0] * dtilda_xi
			            * dxi_disp_y
			            + temp_matrix2[0][1] * dtilda_eta
			            * deta_disp_y;

			double D_matrix[3][3] = {{0.0}};
			int stress_type_flag = 1;
			if (stress_type_flag == 0) { //平面応力状態
				temp = E * (1.0 - nu * nu);
				D_matrix[0][0] = temp;
				D_matrix[0][1] = nu * temp;
				D_matrix[1][0] = nu * temp;
				D_matrix[1][1] = temp;
				D_matrix[2][2] = (1.0 - nu) / 2.0 * temp;
			} else if (stress_type_flag == 1) { //平面ひずみ状態(2Dの場合はこっち)
				temp = E * (1.0 - nu) / (1.0 + nu) / (1.0 - 2.0 * nu);
				D_matrix[0][0] = temp;
				D_matrix[0][1] = nu / (1.0 - nu) * temp;
				D_matrix[1][0] = nu / (1.0 - nu) * temp;
				D_matrix[1][1] = temp;
				D_matrix[2][2] = (1.0 - 2.0 * nu) / 2.0 / (1.0 - nu) * temp;
			}

			stress_xx = D_matrix[0][0] * strain_xx
			            + D_matrix[0][1] * strain_yy;
			stress_yy = D_matrix[1][0] * strain_xx
			            + D_matrix[1][1] * strain_yy;
			stress_xy = D_matrix[2][2] * strain_xy;

			temp = sqrt(disp_x * disp_x + disp_y * disp_y);
			printf("x:   % 1.8e\n", px);
			printf("y:   % 1.8e\n", py);
			printf("xi:  % 1.8e\n", temp_xi);
			printf("eta: % 1.8e\n", temp_eta);
			printf("Displacement x: % 1.8e\n", disp_x);
			printf("Displacement y: % 1.8e\n", disp_y);
			printf("Displacement  : % 1.8e\n", temp);
			printf("Strain xx: % 1.8e\n", strain_xx);
			printf("Strain yy: % 1.8e\n", strain_yy);
			printf("Strain xy: % 1.8e\n", strain_xy);
			printf("Stress xx: % 1.8e\n", stress_xx);
			printf("Stress yy: % 1.8e\n", stress_yy);
			printf("Stress xy: % 1.8e\n", stress_xy);

			return i;
		}

		InverseMatrix_2D(temp_matrix);

		temp_dxi = temp_matrix[0][0] * (px - temp_x)
		           + temp_matrix[0][1] * (py - temp_y);
		temp_deta = temp_matrix[1][0] * (px - temp_x)
		            + temp_matrix[1][1] * (py - temp_y);
		temp_xi = temp_xi + temp_dxi;
		temp_eta = temp_eta + temp_deta;
		if (temp_xi < knot_vec_xi[0])
			temp_xi = knot_vec_xi[0];
		if (temp_xi > knot_vec_xi[cntl_p_n_xi + order_xi])
			temp_xi = knot_vec_xi[cntl_p_n_xi + order_xi];
		if (temp_eta < knot_vec_eta[0])
			temp_eta = knot_vec_eta[0];
		if (temp_eta > knot_vec_eta[cntl_p_n_eta + order_eta])
			temp_eta = knot_vec_eta[cntl_p_n_eta + order_eta];

		//temp_tol = sqrt(temp_dxi * temp_dxi + temp_deta * temp_deta);
		//printf("% 1.15e % 1.15e % 1.15e\n", temp_xi, temp_eta, temp_tol);
	}
	//printf("% 1.8e % 1.8e\n", temp_x, temp_y);
	return 0;
}