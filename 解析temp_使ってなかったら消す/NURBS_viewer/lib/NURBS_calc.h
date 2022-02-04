#ifndef NURBS_CALC_H
#define NURBS_CALC_H

/* 基底関数と基底関数の右微分係数 */
double rBasisFunc(double *knot_vec, int knot_index,
                  int order, double xi,
                  double *output, double *d_output);
/*
入力
double *knot_vec ノットベクトル
int knot_index ノットベクトル指標
int order 基底関数の次数
double xi パラメータ座標

出力
double *output 基底関数の値
double *d_output 基底関数の微分係数
*/

/* 基底関数と基底関数の左微分係数 */
double lBasisFunc(double *knot_vec, int knot_index,
                  int knot_n, int order, double xi,
                  double *output, double *d_output);

/* NURBS曲面とNURBS曲面の右微分係数 */
double rNURBS_surface(double *knot_vec_xi, double *knot_vec_eta,
                      double *cntl_px, double *cntl_py,
                      int cntl_p_n_xi, int cntl_p_n_eta,
                      double *weight, int order_xi, int order_eta,
                      double xi, double eta,
                      double *output_x, double *output_y,
                      double *output_dxi_x, double *output_deta_x,
                      double *output_dxi_y, double *output_deta_y);
/*
入力
double *knot_vec_xi; ξノットベクトル
double *knot_vec_eta; ηノットベクトル
double *cntl_px; コントロールポイントx座標
double *cntl_py; コントロールポイントy座標
int cntl_p_n_xi; ξコントロールポイント数
int cntl_p_n_eta; ηコントロールポイント数
double *weight; 重み
int order_xi; ξ基底関数の次数
int order_eta; η基底関数の次数
double xi; ξパラメータ座標
double eta; ηパラメータ座標

出力
double *output_x; x座標値
double *output_y; y座標値
double *output_dxi_x; ∂x/∂ξ
double *output_deta_x; ∂x/∂η
double *output_dxi_y; ∂y/∂ξ
double *output_deta_y; ∂y/∂η
*/

/* NURBS曲面とNURBS曲面の左微分係数 */
double lNURBS_surface(double *knot_vec_xi, double *knot_vec_eta,
                      double *cntl_px, double *cntl_py,
                      int cntl_p_n_xi, int cntl_p_n_eta,
                      double *weight, int order_xi, int order_eta,
                      double xi, double eta,
                      double *output_x, double *output_y,
                      double *output_dxi_x, double *output_deta_x,
                      double *output_dxi_y, double *output_deta_y);

/* NURBS曲面とNURBS曲面のξ：右微分係数，η：左微分係数 */
double rlNURBS_surface(double *knot_vec_xi, double *knot_vec_eta,
                       double *cntl_px, double *cntl_py,
                       int cntl_p_n_xi, int cntl_p_n_eta,
                       double *weight, int order_xi, int order_eta,
                       double xi, double eta,
                       double *output_x, double *output_y,
                       double *output_dxi_x, double *output_deta_x,
                       double *output_dxi_y, double *output_deta_y);

/* NURBS曲面とNURBS曲面のξ：左微分係数，η：右微分係数 */
double lrNURBS_surface(double *knot_vec_xi, double *knot_vec_eta,
                       double *cntl_px, double *cntl_py,
                       int cntl_p_n_xi, int cntl_p_n_eta,
                       double *weight, int order_xi, int order_eta,
                       double xi, double eta,
                       double *output_x, double *output_y,
                       double *output_dxi_x, double *output_deta_x,
                       double *output_dxi_y, double *output_deta_y);

double InverseMatrix_2D(double M[2][2]);

double NURBS_line(double *knot_vec, double *cntl_p, double *weight,
                  int cntl_p_n, int order, double xi,
                  double *output, double *d_output);

#endif  /* NURBS_CALC_H */
