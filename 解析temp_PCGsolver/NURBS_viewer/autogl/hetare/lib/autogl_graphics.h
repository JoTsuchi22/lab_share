/* AutoGL */

#ifndef __AUTOGL_GRAPHICS__
#define __AUTOGL_GRAPHICS__



/* レンダリング */

/* ビューの背景色を設定する。 */
void AutoGL_SetBackgroundColor 
(double red, double green, double blue);  
/* red, green, blue : 赤、緑、青について、それぞれ0.0 - 1.0のレンジの強度。 */

/* 現在の色を設定する。 */
void AutoGL_SetColor 
(double red, double green, double blue);  
/* red, green, blue : 赤、緑、青について、それぞれ0.0 - 1.0のレンジの強度。 */

/* 線分を現在の色で描画する。*/
void AutoGL_DrawLine 
(double x0, double y0, double z0,    
 double x1, double y1, double z1); 
/* (x0, y0, z0), (x1, y1, z1) : 二点の世界座標 */

/* 三角形を現在の色で塗り潰す。 */
void AutoGL_DrawTriangle 
(double x0, double y0, double z0,  
 double x1, double y1, double z1,  
 double x2, double y2, double z2); 
/* (x0, y0, z0), (x1, y1, z1), (x2, y2, z2) : 三点の世界座標 */

/* 三角形を色補間により塗り潰す。 */
void AutoGL_DrawColorInterpolatedTriangle 
(double x0, double y0, double z0, double red0, double green0, double blue0,  
 double x1, double y1, double z1, double red1, double green1, double blue1,  
 double x2, double y2, double z2, double red2, double green2, double blue2); 
/* (x0, y0, z0), (x1, y1, z1), (x2, y2, z2) : 三点の世界座標 */
/* (red0, green0, blue0), (red1, green1, blue1), (red2, green2, blue2) : 三点の色 */

/* 四辺形を現在の色で塗り潰す。 */
void AutoGL_DrawQuadrangle 
(double x0, double y0, double z0,
 double x1, double y1, double z1,
 double x2, double y2, double z2,
 double x3, double y3, double z3);
/* (x0, y0, z0), (x1, y1, z1), (x2, y2, z2), (x3, y3, z3) : 四点の世界座標 */
/* 四点は同一平面上になくてもよい。 */

/* 文字列を現在の色で描画する。 */
void AutoGL_DrawString 
(double x, double y, double z,
 const char *str);
/* (x, y, z) : 文字列の原点の世界座標 */
/* str : 描画される文字列 */




/* ビューイングパラメータと座標変換 */

/* 世界座標系 : モデルが置かれている座標系。 */
/*   三次元の実数座標系。*/

/* 視点座標系 : 世界座標系が視野変換されたもの。 */
/*   三次元の実数座標系。*/
/*   画面の中心が(0, 0, 0)に対応する。*/
/*   右方向がX方向, 上がY方向, 奥行きの手前側にZ方向。*/

/* デバイス座標系 : 視点座標系が実際の画面サイズにスケールされたもの。三次元。*/
/*   三次元の整数座標系。*/
/*   画面の中心が(0, 0, 0)に対応する。*/
/*   右方向がX方向, 上がY方向, 奥行きの手前側にZ方向。*/
/*   X Window SystemやWindowsのデバイス座標系とは異なることに注意！！！ */

/* ビューのサイズを設定する。*/
void AutoGL_SetViewSize (double size);
/* size : 実際の画面の半分の大きさに対応する,世界座標におけるビューのサイズ。*/

/* ビューの注視点を設定する。*/
void AutoGL_SetViewCenter (double x, double y, double z);
/* (x, y, z) : 注視点の座標 */

/* ビューの視線方向を設定する。 */
void AutoGL_SetViewDirection (double x, double y, double z);
/* (x, y, z) : 視線方向ベクトル */

/* ビューの視線上向き方向を設定する。 */
void AutoGL_SetViewUpVector (double x, double y, double z);
/* (x, y, z) : 視線上向き方向ベクトル */

/* ビューのサイズを得る。 */
double AutoGL_GetViewSize (void);

/* ビューの注視点を得る。 */
void AutoGL_GetViewCenter 
(double *x_OUT, double *y_OUT, double *z_OUT);
/* (x_OUT, y_OUT, z_OUT) : 注視点の座標 */

/* ビューの視線方向を得る。 */
void AutoGL_GetViewDirection 
(double *x_OUT, double *y_OUT, double *z_OUT);
/* (x_OUT, y_OUT, z_OUT) : 視線方向ベクトル */

/* ビューの視線上向き方向を得る。 */
void AutoGL_GetViewUpVector 
(double *x_OUT, double *y_OUT, double *z_OUT);
/* (x_OUT, y_OUT, z_OUT) : 視線上向き方向ベクトル */

/* ビューのデバイス座標系での高さと幅を得る。 */
void AutoGL_GetViewRangeDc 
(int *widthDc_OUT, int *heightDc_OUT);
/* widthDc_OUT : ビューの幅 */
/* heightDc_OUT : ビューの高さ */

/* ビューのデバイス座標系でのサイズを得る。 */
int AutoGL_GetViewSizeDc (void); 

/* デバイス座標から視点座標を得る。 */
void AutoGL_GetVcOfDc 
(double *vcX_OUT, double *vcY_OUT, double *vcZ_OUT,
 int dcX, int dcY, int dcZ);
/* (dcX, dcY, dcZ) : デバイス座標 */
/* (vcX_OUT, vcY_OUT, vcZ_OUT) : 視点座標 */

/* 視点座標から世界座標を得る。 */
void AutoGL_GetPositionOfVc 
(double *x_OUT, double *y_OUT, double *z_OUT,
 double vcX, double vcY, double vcZ);
/* (vcX, vcY, vcZ) : 視点座標 */
/* (x_OUT, y_OUT, z_OUT) : 世界座標 */

/* デバイス座標から世界座標を得る。 */
void AutoGL_GetPositionOfDc 
(double *x_OUT, double *y_OUT, double *z_OUT,
 int dcX, int dcY, int dcZ);
/* (dcX, dcY, dcZ) : デバイス座標 */
/* (x_OUT, y_OUT, z_OUT) : 世界座標 */



/* イベント処理 */

/* AUTOGL_EVENT_ イベントのタイプ */
/* ポインタのボタンが押された。 */
#define AUTOGL_EVENT_POINTING_DEVICE_PRESS 0
/* ポインタがドラッグされた。 */
#define AUTOGL_EVENT_POINTING_DEVICE_DRAG 1
/* ポインタのボタンが離された。 */
#define AUTOGL_EVENT_POINTING_DEVICE_RELEASE 2

/* イベントのタイプを得る。 */
int AutoGL_GetViewEvent (void);   /* AUTOGL_EVENT_ */

/* ポインティングデバイスの位置を得る。 */
void AutoGL_GetPointingDevicePositionDc 
(int *dcX_OUT, int *dcY_OUT, int *dcZ_OUT);
/* dcX_OUT, dcY_OUT, dcZ_OUT : ポインタのデバイス座標 */



#endif  /* __AUTOGL_GRAPHICS__ */



