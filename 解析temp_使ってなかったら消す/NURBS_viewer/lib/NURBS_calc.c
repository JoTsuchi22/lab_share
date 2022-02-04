#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <time.h>

#include "NURBS_calc.h"

double rBasisFunc(double *knot_vec, int knot_index,
                  int order, double xi,
                  double *output, double *d_output) {
	int p, j;
	double sum1 = 0.0;
	double sum2 = 0.0;
	double temp_basis[order][order];
	(*output) = 0.0;
	(*d_output) = 0.0;


	if ( knot_vec[knot_index] <= xi && xi <= knot_vec[knot_index + order + 1] ) {
		if (knot_index == 0) {
			for (j = 0; j <= order; j++) {
				if ( (knot_vec[j] <= xi) && (xi <= knot_vec[j + 1]) ) {
					temp_basis[j][0] = 1.0;
				} else {
					temp_basis[j][0] = 0.0;
				}
			}
		} else {
			for (j = 0; j <= order; j++) {
				if ( (knot_vec[knot_index + j] < xi) && (xi <= knot_vec[knot_index + j + 1]) ) {
					temp_basis[j][0] = 1.0;
				} else {
					temp_basis[j][0] = 0.0;
				}
			}
		}

		if ( order > 0 ) {
			for (p = 1; p <= order; p++) {
				for (j = 0; j <= order - p; j++) {
					sum1 = 0.0;
					sum2 = 0.0;
					if ( (knot_vec[knot_index + j + p] - knot_vec[knot_index + j]) != 0.0) {
						sum1 = ( xi - knot_vec[knot_index + j] )
						       / (knot_vec[knot_index + j + p] - knot_vec[knot_index + j])
						       * temp_basis[j][p - 1];
					}
					if ( (knot_vec[knot_index + j + p + 1] - knot_vec[knot_index + j + 1]) != 0.0) {
						sum2 = (knot_vec[knot_index + j + p + 1] - xi)
						       / (knot_vec[knot_index + j + p + 1] - knot_vec[knot_index + j + 1])
						       * temp_basis[j + 1][p - 1];
					}
					temp_basis[j][p] = sum1 + sum2;
				}
			}
			sum1 = 0.0;
			sum2 = 0.0;
			if ( (knot_vec[knot_index + order] - knot_vec[knot_index]) != 0.0) {
				sum1 = order
				       / (knot_vec[knot_index + order] - knot_vec[knot_index])
				       * temp_basis[0][order - 1];
			}
			if ( (knot_vec[knot_index + order + 1] - knot_vec[knot_index + 1]) != 0.0) {
				sum2 = order
				       / (knot_vec[knot_index + order + 1] - knot_vec[knot_index + 1])
				       * temp_basis[1][order - 1];
			}
		}
		(*output) = temp_basis[0][order];
		(*d_output) = sum1 - sum2;
	}
	return (*output);
}

double lBasisFunc(double *knot_vec, int knot_index,
                  int cntl_p_n, int order, double xi,
                  double *output, double *d_output) {
	int p, j;
	double sum1 = 0.0;
	double sum2 = 0.0;
	double temp_basis[order][order];
	(*output) = 0.0;
	(*d_output) = 0.0;


	if ( knot_vec[knot_index] <= xi && xi <= knot_vec[knot_index + order + 1] ) {
		if (knot_index == cntl_p_n - 1) {
			for (j = 0; j <= order; j++) {
				if ( (knot_vec[cntl_p_n - 1 + j] <= xi) && (xi <= knot_vec[cntl_p_n + j]) ) {
					temp_basis[j][0] = 1.0;
				} else {
					temp_basis[j][0] = 0.0;
				}
			}
		} else {
			for (j = 0; j <= order; j++) {
				if ( (knot_vec[knot_index + j] <= xi) && (xi < knot_vec[knot_index + j + 1]) ) {
					temp_basis[j][0] = 1.0;
				} else {
					temp_basis[j][0] = 0.0;
				}
			}
		}

		if ( order > 0 ) {
			for (p = 1; p <= order; p++) {
				for (j = 0; j <= order - p; j++) {
					sum1 = 0.0;
					sum2 = 0.0;
					if ( (knot_vec[knot_index + j + p] - knot_vec[knot_index + j]) != 0.0) {
						sum1 = ( xi - knot_vec[knot_index + j] )
						       / (knot_vec[knot_index + j + p] - knot_vec[knot_index + j])
						       * temp_basis[j][p - 1];
					}
					if ( (knot_vec[knot_index + j + p + 1] - knot_vec[knot_index + j + 1]) != 0.0) {
						sum2 = (knot_vec[knot_index + j + p + 1] - xi)
						       / (knot_vec[knot_index + j + p + 1] - knot_vec[knot_index + j + 1])
						       * temp_basis[j + 1][p - 1];
					}
					temp_basis[j][p] = sum1 + sum2;
				}
			}
			sum1 = 0.0;
			sum2 = 0.0;
			if ( (knot_vec[knot_index + order] - knot_vec[knot_index]) != 0.0) {
				sum1 = order
				       / (knot_vec[knot_index + order] - knot_vec[knot_index])
				       * temp_basis[0][order - 1];
			}
			if ( (knot_vec[knot_index + order + 1] - knot_vec[knot_index + 1]) != 0.0) {
				sum2 = order
				       / (knot_vec[knot_index + order + 1] - knot_vec[knot_index + 1])
				       * temp_basis[1][order - 1];
			}
		}
		(*output) = temp_basis[0][order];
		(*d_output) = sum1 - sum2;
	}
	return (*output);
}



double rNURBS_surface(double *knot_vec_xi, double *knot_vec_eta,
                      double *cntl_px, double *cntl_py,
                      int cntl_p_n_xi, int cntl_p_n_eta,
                      double *weight, int order_xi, int order_eta,
                      double xi, double eta,
                      double *output_x, double *output_y,
                      double *output_dxi_x, double *output_deta_x,
                      double *output_dxi_y, double *output_deta_y) {
	int i, j, temp_index;
	double temp1, temp2, temp3;
	double molecule_x, molecule_y;
	double dxi_molecule_x, dxi_molecule_y;
	double deta_molecule_x, deta_molecule_y;
	double denominator, dxi_denominator, deta_denominator;
	double temp_output_xi, temp_output_eta;
	double temp_d_output_xi, temp_d_output_eta;
	molecule_x = 0.0;
	molecule_y = 0.0;
	denominator = 0.0;
	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_denominator = 0.0;
	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_denominator = 0.0;

	int index_min_xi = 0;
	int index_max_xi = cntl_p_n_xi - 1;
	int index_min_eta = 0;
	int index_max_eta = cntl_p_n_eta - 1;

	for (i = 0; i < cntl_p_n_xi; i++) {
		if ( knot_vec_xi[i + 1] >= xi ) {
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}
	if (index_min_xi < 0) index_min_xi = 0;
	if (index_max_xi > cntl_p_n_xi) index_max_xi = cntl_p_n_xi;

	for (i = 0; i < cntl_p_n_eta; i++) {
		if ( knot_vec_eta[i + 1] >= eta ) {
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0) index_min_eta = 0;
	if (index_max_eta > cntl_p_n_eta) index_max_eta = cntl_p_n_eta;

	for (i = index_min_xi; i <= index_max_xi; i++) {
		rBasisFunc(knot_vec_xi, i, order_xi, xi,
		           &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j <= index_max_eta; j++) {
			rBasisFunc(knot_vec_eta, j, order_eta, eta,
			           &temp_output_eta, &temp_d_output_eta);
			temp_index = i + j * cntl_p_n_xi;
			temp1 = temp_output_xi * temp_output_eta
			        * weight[temp_index];
			temp2 = temp_d_output_xi * temp_output_eta
			        * weight[temp_index];
			temp3 = temp_output_xi * temp_d_output_eta
			        * weight[temp_index];
			molecule_x += temp1 * cntl_px[temp_index];
			molecule_y += temp1 * cntl_py[temp_index];
			denominator += temp1;
			dxi_molecule_x += temp2 * cntl_px[temp_index];
			dxi_molecule_y += temp2 * cntl_py[temp_index];
			dxi_denominator += temp2;
			deta_molecule_x += temp3 * cntl_px[temp_index];
			deta_molecule_y += temp3 * cntl_py[temp_index];
			deta_denominator += temp3;
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;

	temp1 = denominator * denominator;
	(*output_dxi_x) = ( dxi_molecule_x * denominator
	                    - molecule_x * dxi_denominator) / temp1;
	(*output_dxi_y) = ( dxi_molecule_y * denominator
	                    - molecule_y * dxi_denominator) / temp1;
	(*output_deta_x) = ( deta_molecule_x * denominator
	                     - molecule_x * deta_denominator) / temp1;
	(*output_deta_y) = ( deta_molecule_y * denominator
	                     - molecule_y * deta_denominator) / temp1;
	return denominator;
}

double lNURBS_surface(double *knot_vec_xi, double *knot_vec_eta,
                      double *cntl_px, double *cntl_py,
                      int cntl_p_n_xi, int cntl_p_n_eta,
                      double *weight, int order_xi, int order_eta,
                      double xi, double eta,
                      double *output_x, double *output_y,
                      double *output_dxi_x, double *output_deta_x,
                      double *output_dxi_y, double *output_deta_y) {
	int i, j, temp_index;
	double temp1, temp2, temp3;
	double molecule_x, molecule_y;
	double dxi_molecule_x, dxi_molecule_y;
	double deta_molecule_x, deta_molecule_y;
	double denominator, dxi_denominator, deta_denominator;
	double temp_output_xi, temp_output_eta;
	double temp_d_output_xi, temp_d_output_eta;
	molecule_x = 0.0;
	molecule_y = 0.0;
	denominator = 0.0;
	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_denominator = 0.0;
	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_denominator = 0.0;

	int index_min_xi = 0;
	int index_max_xi = cntl_p_n_xi - 1;
	int index_min_eta = 0;
	int index_max_eta = cntl_p_n_eta - 1;

	for (i = 0; i < cntl_p_n_xi; i++) {
		if ( knot_vec_xi[i + 1] >= xi ) {
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}
	if (index_min_xi < 0) index_min_xi = 0;
	if (index_max_xi > cntl_p_n_xi) index_max_xi = cntl_p_n_xi;

	for (i = 0; i < cntl_p_n_eta; i++) {
		if ( knot_vec_eta[i + 1] >= eta ) {
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0) index_min_eta = 0;
	if (index_max_eta > cntl_p_n_eta) index_max_eta = cntl_p_n_eta;

	for (i = index_min_xi; i <= index_max_xi; i++) {
		lBasisFunc(knot_vec_xi, i,
		           cntl_p_n_xi, order_xi, xi,
		           &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j <= index_max_eta; j++) {
			lBasisFunc(knot_vec_eta, j,
			           cntl_p_n_eta, order_eta, eta,
			           &temp_output_eta, &temp_d_output_eta);
			temp_index = i + j * cntl_p_n_xi;
			temp1 = temp_output_xi * temp_output_eta
			        * weight[temp_index];
			temp2 = temp_d_output_xi * temp_output_eta
			        * weight[temp_index];
			temp3 = temp_output_xi * temp_d_output_eta
			        * weight[temp_index];
			molecule_x += temp1 * cntl_px[temp_index];
			molecule_y += temp1 * cntl_py[temp_index];
			denominator += temp1;
			dxi_molecule_x += temp2 * cntl_px[temp_index];
			dxi_molecule_y += temp2 * cntl_py[temp_index];
			dxi_denominator += temp2;
			deta_molecule_x += temp3 * cntl_px[temp_index];
			deta_molecule_y += temp3 * cntl_py[temp_index];
			deta_denominator += temp3;
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;

	temp1 = denominator * denominator;
	(*output_dxi_x) = ( dxi_molecule_x * denominator
	                    - molecule_x * dxi_denominator) / temp1;
	(*output_dxi_y) = ( dxi_molecule_y * denominator
	                    - molecule_y * dxi_denominator) / temp1;
	(*output_deta_x) = ( deta_molecule_x * denominator
	                     - molecule_x * deta_denominator) / temp1;
	(*output_deta_y) = ( deta_molecule_y * denominator
	                     - molecule_y * deta_denominator) / temp1;
	return denominator;
}

double rlNURBS_surface(double *knot_vec_xi, double *knot_vec_eta,
                       double *cntl_px, double *cntl_py,
                       int cntl_p_n_xi, int cntl_p_n_eta,
                       double *weight, int order_xi, int order_eta,
                       double xi, double eta,
                       double *output_x, double *output_y,
                       double *output_dxi_x, double *output_deta_x,
                       double *output_dxi_y, double *output_deta_y) {
	int i, j, temp_index;
	double temp1, temp2, temp3;
	double molecule_x, molecule_y;
	double dxi_molecule_x, dxi_molecule_y;
	double deta_molecule_x, deta_molecule_y;
	double denominator, dxi_denominator, deta_denominator;
	double temp_output_xi, temp_output_eta;
	double temp_d_output_xi, temp_d_output_eta;
	molecule_x = 0.0;
	molecule_y = 0.0;
	denominator = 0.0;
	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_denominator = 0.0;
	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_denominator = 0.0;

	int index_min_xi = 0;
	int index_max_xi = cntl_p_n_xi - 1;
	int index_min_eta = 0;
	int index_max_eta = cntl_p_n_eta - 1;

	for (i = 0; i < cntl_p_n_xi; i++) {
		if ( knot_vec_xi[i + 1] >= xi ) {
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}
	if (index_min_xi < 0) index_min_xi = 0;
	if (index_max_xi > cntl_p_n_xi) index_max_xi = cntl_p_n_xi;

	for (i = 0; i < cntl_p_n_eta; i++) {
		if ( knot_vec_eta[i + 1] >= eta ) {
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0) index_min_eta = 0;
	if (index_max_eta > cntl_p_n_eta) index_max_eta = cntl_p_n_eta;

	for (i = index_min_xi; i <= index_max_xi; i++) {
		rBasisFunc(knot_vec_xi, i, order_xi, xi,
		           &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j <= index_max_eta; j++) {
			lBasisFunc(knot_vec_eta, j,
			           cntl_p_n_eta, order_eta, eta,
			           &temp_output_eta, &temp_d_output_eta);
			temp_index = i + j * cntl_p_n_xi;
			temp1 = temp_output_xi * temp_output_eta
			        * weight[temp_index];
			temp2 = temp_d_output_xi * temp_output_eta
			        * weight[temp_index];
			temp3 = temp_output_xi * temp_d_output_eta
			        * weight[temp_index];
			molecule_x += temp1 * cntl_px[temp_index];
			molecule_y += temp1 * cntl_py[temp_index];
			denominator += temp1;
			dxi_molecule_x += temp2 * cntl_px[temp_index];
			dxi_molecule_y += temp2 * cntl_py[temp_index];
			dxi_denominator += temp2;
			deta_molecule_x += temp3 * cntl_px[temp_index];
			deta_molecule_y += temp3 * cntl_py[temp_index];
			deta_denominator += temp3;
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;

	temp1 = denominator * denominator;
	(*output_dxi_x) = ( dxi_molecule_x * denominator
	                    - molecule_x * dxi_denominator) / temp1;
	(*output_dxi_y) = ( dxi_molecule_y * denominator
	                    - molecule_y * dxi_denominator) / temp1;
	(*output_deta_x) = ( deta_molecule_x * denominator
	                     - molecule_x * deta_denominator) / temp1;
	(*output_deta_y) = ( deta_molecule_y * denominator
	                     - molecule_y * deta_denominator) / temp1;
	return denominator;
}

double lrNURBS_surface(double *knot_vec_xi, double *knot_vec_eta,
                       double *cntl_px, double *cntl_py,
                       int cntl_p_n_xi, int cntl_p_n_eta,
                       double *weight, int order_xi, int order_eta,
                       double xi, double eta,
                       double *output_x, double *output_y,
                       double *output_dxi_x, double *output_deta_x,
                       double *output_dxi_y, double *output_deta_y) {
	int i, j, temp_index;
	double temp1, temp2, temp3;
	double molecule_x, molecule_y;
	double dxi_molecule_x, dxi_molecule_y;
	double deta_molecule_x, deta_molecule_y;
	double denominator, dxi_denominator, deta_denominator;
	double temp_output_xi, temp_output_eta;
	double temp_d_output_xi, temp_d_output_eta;
	molecule_x = 0.0;
	molecule_y = 0.0;
	denominator = 0.0;
	dxi_molecule_x = 0.0;
	dxi_molecule_y = 0.0;
	dxi_denominator = 0.0;
	deta_molecule_x = 0.0;
	deta_molecule_y = 0.0;
	deta_denominator = 0.0;

	int index_min_xi = 0;
	int index_max_xi = cntl_p_n_xi - 1;
	int index_min_eta = 0;
	int index_max_eta = cntl_p_n_eta - 1;

	for (i = 0; i < cntl_p_n_xi; i++) {
		if ( knot_vec_xi[i + 1] >= xi ) {
			index_min_xi = i - order_xi;
			index_max_xi = i + 1;
			break;
		}
	}
	if (index_min_xi < 0) index_min_xi = 0;
	if (index_max_xi > cntl_p_n_xi) index_max_xi = cntl_p_n_xi;

	for (i = 0; i < cntl_p_n_eta; i++) {
		if ( knot_vec_eta[i + 1] >= eta ) {
			index_min_eta = i - order_eta;
			index_max_eta = i + 1;
			break;
		}
	}
	if (index_min_eta < 0) index_min_eta = 0;
	if (index_max_eta > cntl_p_n_eta) index_max_eta = cntl_p_n_eta;

	for (i = index_min_xi; i <= index_max_xi; i++) {
		lBasisFunc(knot_vec_xi, i,
		           cntl_p_n_xi, order_xi, xi,
		           &temp_output_xi, &temp_d_output_xi);
		for (j = index_min_eta; j <= index_max_eta; j++) {
			rBasisFunc(knot_vec_eta, j, order_eta, eta,
			           &temp_output_eta, &temp_d_output_eta);
			temp_index = i + j * cntl_p_n_xi;
			temp1 = temp_output_xi * temp_output_eta
			        * weight[temp_index];
			temp2 = temp_d_output_xi * temp_output_eta
			        * weight[temp_index];
			temp3 = temp_output_xi * temp_d_output_eta
			        * weight[temp_index];
			molecule_x += temp1 * cntl_px[temp_index];
			molecule_y += temp1 * cntl_py[temp_index];
			denominator += temp1;
			dxi_molecule_x += temp2 * cntl_px[temp_index];
			dxi_molecule_y += temp2 * cntl_py[temp_index];
			dxi_denominator += temp2;
			deta_molecule_x += temp3 * cntl_px[temp_index];
			deta_molecule_y += temp3 * cntl_py[temp_index];
			deta_denominator += temp3;
		}
	}
	(*output_x) = molecule_x / denominator;
	(*output_y) = molecule_y / denominator;

	temp1 = denominator * denominator;
	(*output_dxi_x) = ( dxi_molecule_x * denominator
	                    - molecule_x * dxi_denominator) / temp1;
	(*output_dxi_y) = ( dxi_molecule_y * denominator
	                    - molecule_y * dxi_denominator) / temp1;
	(*output_deta_x) = ( deta_molecule_x * denominator
	                     - molecule_x * deta_denominator) / temp1;
	(*output_deta_y) = ( deta_molecule_y * denominator
	                     - molecule_y * deta_denominator) / temp1;
	return denominator;
}

double InverseMatrix_2D(double M[2][2]) {
	int i, j;
	double a[2][2];
	double det = M[0][0] * M[1][1] - M[0][1] * M[1][0];

	if (det == 0) return 9999999.0;

	for (i = 0; i < 2; i++) {
		for (j = 0; j < 2; j++) {
			a[i][j] = M[i][j];
		}
	}

	M[0][0] = a[1][1] / det;  M[0][1] = -a[0][1] / det;
	M[1][0] = -a[1][0] / det; M[1][1] = a[0][0] / det;

	return det;
}

double NURBS_line(double *knot_vec, double *cntl_p, double *weight,
                  int cntl_p_n, int order, double xi,
                  double *output, double *d_output) {
	int i;
	double temp_den, temp_d_den, temp_output, temp_d_output;
	(*output) = 0.0;
	(*d_output) = 0.0;
	temp_den = 0.0;
	temp_d_den = 0.0;

	for (i = 0; i < cntl_p_n; i++) {
		rBasisFunc(knot_vec, i, order, xi, &temp_output, &temp_d_output);
		(*output) += temp_output * cntl_p[i] * weight[i];
		(*d_output) += temp_d_output * cntl_p[i] * weight[i];
	}
	if (*output != 0.0 || *d_output != 0.0) {

		for (i = 0; i < cntl_p_n; i++) {
			rBasisFunc(knot_vec, i, order, xi, &temp_output, &temp_d_output);
			temp_den += temp_output * weight[i];
			temp_d_den += temp_d_output * weight[i];
		}
		(*output) /= temp_den;
		(*d_output) /= temp_d_den;
	}
	return 0.0;
}
