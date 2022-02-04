#ifndef CLICKCALC_H
#define CLICKCALC_H

double TriSignedArea(double x1, double y1,
                     double x2, double y2,
                     double x3, double y3);

int SegmentHit(double p0_x0, double p0_y0,
               double p0_x1, double p0_y1,
               double p1_x0, double p1_y0,
               double p1_x1, double p1_y1);

int CalcXiEtaByNR(double px, double py,
                  double *knot_vec_xi, double *knot_vec_eta,
                  double *cntl_px, double *cntl_py,
                  double *disp_cntl_px, double *disp_cntl_py,
                  int cntl_p_n_xi, int cntl_p_n_eta,
                  double *weight, int order_xi, int order_eta,
                  double *output_xi, double *output_eta);

#endif  /* CLICKCALC_H */
