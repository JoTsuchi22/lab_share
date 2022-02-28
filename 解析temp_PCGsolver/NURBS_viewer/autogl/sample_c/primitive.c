/*  AutoGL ＜お手軽＞ のサンプルプログラム例 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/* AutoGL ”おてがる”を使うには,このヘッダファイルをインクルードする。*/
#include <autogl.h>

/* 基本三次元描画プリミティブの利用 */
/* 引数など詳しい情報は, autogl_graphics.hを参照。*/
/* なお、これはautogl.hからincludeされている。*/



/* コールバック関数群 */

/* ビューの再描画のためのコールバック関数 */
/* ビューが再表示されるごとに呼ばれる。 */
static void RedrawView (void) 
{
  AutoGL_SetColor (1, 1, 1);

  /* 原点ラベルを表示する。*/
  AutoGL_DrawString (0, 0, 0, "O");  /* 原点に原点ラベルの文字列 */

  /* X軸を描画する。*/
  AutoGL_DrawLine (0, 0, 0,
		   50, 0, 0);  /* 原点から(50, 0, 0)への線分 */
  AutoGL_DrawString (50, 0, 0, "X");  
  
  /* Y軸を描画する。*/
  AutoGL_DrawLine (0, 0, 0,
		   0, 50, 0);  /* 原点から(0, 50, 0)への線分 */
  AutoGL_DrawString (0, 50, 0, "Y");  
  
  /* Z軸を描画する。*/
  AutoGL_DrawLine (0, 0, 0,
		   0, 0, 50);  /* 原点から(0, 0, 50)への線分 */
  AutoGL_DrawString (0, 0, 50, "Z");  


  /* 三角形を描画する。 */
  AutoGL_SetColor (1, 1, 1);
  AutoGL_DrawTriangle (-20, 0, 0,
		       0, -20, 0,
		       0, 0, -20);  /* 三点を指定して三角形 */

  /* 直方体(0, 0, 0) -（10, 20, 30)の6枚の面（四辺形）を描画する。 */

  AutoGL_SetColor (0, 1, 1);
  AutoGL_DrawQuadrangle (0, 0, 0,
			 10, 0, 0,
			 10, 20, 0,
			 0, 20, 0);  /* 四点を指定して四辺形 */
  AutoGL_SetColor (1, 1, 0);
  AutoGL_DrawQuadrangle (0, 0, 30,
			 10, 0, 30,
			 10, 20, 30,
			 0, 20, 30);
  AutoGL_SetColor (1, 0, 1);
  AutoGL_DrawQuadrangle (0, 0, 0,
			 0, 20, 0,
			 0, 20, 30,
			 0, 0, 30);
  AutoGL_SetColor (0, 0, 1);
  AutoGL_DrawQuadrangle (10, 0, 0,
			 10, 20, 0,
			 10, 20, 30,
			 10, 0, 30);
  AutoGL_SetColor (0, 1, 0);
  AutoGL_DrawQuadrangle (0, 0, 0,
			 10, 0, 0,
			 10, 0, 30,
			 0, 0, 30);
  AutoGL_SetColor (1, 0, 0);
  AutoGL_DrawQuadrangle (0, 20, 0,
			 10, 20, 0,
			 10, 20, 30,
			 0, 20, 30);
}



/* 関数AutoGL_SetUpはユーザー側プログラムごとに必ず一つ用意すること。*/
/* この関数はプログラム起動直後にAutoGLライブラリから呼び出される。*/
/* このとき、argc, argvはmain関数から直接渡される。*/
/* ここで、コールバック関数や制御変数などを登録する。*/
void AutoGL_SetUp (int argc, char *argv[]) 
{
  /* ビューのみの構成例 */

  /* ビューイングパラメータ（サイズ,視線方向など）を設定する。*/
  AutoGL_SetViewSize (70);                 /* ビューの大きさ */
  AutoGL_SetViewCenter (10, 10, 10);       /* 注視点 */
  AutoGL_SetViewDirection (1, 2, 3);       /* 視線方向 */

  /* おまじない */
  AutoGL_SetViewRedrawCallback (RedrawView); 
  AutoGL_SetMode3D (AUTOGL_MODE_3D_ROTATE_XY); /* マウスで回転 */
  AutoGL_SetDefaultCallbackInMode3D (NULL);    
  AutoGL_SetPanelInMode3D ();
  AutoGL_EnableDragInMode3D ();                /* ドラッグ有効 */
}






