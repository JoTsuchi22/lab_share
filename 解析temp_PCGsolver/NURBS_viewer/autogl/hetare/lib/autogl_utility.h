/* AutoGL */

#ifndef __AUTOGL_UTILITY__
#define __AUTOGL_UTILITY__

/* ユーティリティ */



/* 描画プリミティブ */

/* コンタ-の値に対応する色を得る。 */
void AutoGL_GetContourColor 
(double *red_OUT, double *green_OUT, double *blue_OUT,
 double grade);
/* red_OUT, green_OUT, blue_OUT : 赤、緑、青について、それぞれ0.0 - 1.0のレンジの強度。 */
/* grade : 0.0 - 1.0に規格化された値。 */

/* 三角形をコンタ-で塗り潰す。 */
void AutoGL_DrawContourTriangle 
(double x0, double y0, double z0, double grade0,  
 double x1, double y1, double z1, double grade1,  
 double x2, double y2, double z2, double grade2); 
/* (x0, y0, z0), (x1, y1, z1), (x2, y2, z2) : 三点の世界座標 */
/* grade0, grade1, grade2 : 三点のスカラー値（0.0 - 1.0に規格化） */



/* 三次元ビューイングパラメータ操作 */

/* ビューイングパラメータ設定機能をパネルに準備する。 */
void AutoGL_SetPanelInView3D (void);

/* ビューを相対的に回転する。 */
void AutoGL_RotateView 
(double offsetX, double offsetY, double offsetZ);
/* offsetX, offsetY, offsetZ     around -1.0 ... 1.0 */



/* 三次元イベント処理 */

/*  AUTOGL_MODE_3D_  */
/* ユーザー定義イベント処理モード */
#define AUTOGL_MODE_3D_OTHER 0
/* 移動モード */
#define AUTOGL_MODE_3D_TRANSLATE 1
/* 回転モード（XY軸まわり） */
#define AUTOGL_MODE_3D_ROTATE_XY 2
/* 回転モード（Z軸まわり） */
#define AUTOGL_MODE_3D_ROTATE_Z 3
/* 拡大縮小モード */
#define AUTOGL_MODE_3D_SCALE 4

/* 移動回転拡大縮小機能をパネルに準備する。 */
void AutoGL_SetPanelInMode3D (void);

/* イベント処理モードを設定する。 */
void AutoGL_SetMode3D (int mode);    /*  AUTOGL_MODE_3D_  */

/* イベント処理モードを設定する。 */
void AutoGL_SetDefaultCallbackInMode3D 
(AutoGL_CallbackType otherCallback);
/* otherCallback : ユーザー定義イベント処理コールバック */
/* これはNULLでもよい。 */

/* ドラッグ中に描画を行う。 */
void AutoGL_EnableDragInMode3D (void);

/* イベント処理モードを得る。 */
int AutoGL_GetMode3D (void);    /*  AUTOGL_MODE_3D_  */

/* ポインティングデバイスがヒットしたかを調べる。 */
int AutoGL_PointingDeviceIsHit3D 
(double *parameter_OUT,
 double pX, double pY, double pZ,
 int toleranceDc);
/* parameter_OUT : 深さ情報 */
/* (pX, pY, pZ) : 世界座標 */
/* toleranceDc : デバイス座標系での幅 */



#endif  /* __AUTOGL_UTILITY__ */



