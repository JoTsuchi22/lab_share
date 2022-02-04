/*  AutoGL ＜お手軽＞ のサンプルプログラム例 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/* AutoGL ”おてがる”を使うには,このヘッダファイルをインクルードする。*/
#include <autogl.h>



/* ビューが再表示されるごとに呼ばれる。 */
static void RedrawView (void) 
{
  /* 座標軸を白で描画する。*/
  AutoGL_SetColor (1, 1, 1);
  AutoGL_DrawString (0, 0, 0, "O");   /* 原点に原点ラベルの文字列 */

  AutoGL_DrawLine (0, 0, 0,
		   50, 0, 0);         /* 原点から(50, 0, 0)への線分 */
  AutoGL_DrawString (50, 0, 0, "X");  /* X軸 */
  
  AutoGL_DrawLine (0, 0, 0,
		   0, 50, 0);         /* 原点から(0, 50, 0)への線分 */
  AutoGL_DrawString (0, 50, 0, "Y");  /* Y軸 */
  
  AutoGL_DrawLine (0, 0, 0,
		   0, 0, 50);         /* 原点から(0, 0, 50)への線分 */
  AutoGL_DrawString (0, 0, 50, "Z");  /* Z軸 */


  /* 三角形を赤で描画する。 */
  AutoGL_SetColor (1, 0, 0);
  AutoGL_DrawTriangle (-20, 0, 0,
		       0, -20, 0,
		       0, 0, -20);    /* 三点を指定して三角形 */

  /* 直方体(0, 0, 0) -（10, 20, 30)をシアンで描画する。 */
  AutoGL_SetColor (0, 1, 1);
  AutoGL_DrawQuadrangle (0, 0, 0,
			 10, 0, 0,
			 10, 20, 0,
			 0, 20, 0);   /* 四点を指定して四辺形 */
  AutoGL_DrawQuadrangle (0, 0, 30,
			 10, 0, 30,
			 10, 20, 30,
			 0, 20, 30);
  AutoGL_DrawQuadrangle (0, 0, 0,
			 0, 20, 0,
			 0, 20, 30,
			 0, 0, 30);
  AutoGL_DrawQuadrangle (10, 0, 0,
			 10, 20, 0,
			 10, 20, 30,
			 10, 0, 30);
  AutoGL_DrawQuadrangle (0, 0, 0,
			 10, 0, 0,
			 10, 0, 30,
			 0, 0, 30);
  AutoGL_DrawQuadrangle (0, 20, 0,
			 10, 20, 0,
			 10, 20, 30,
			 0, 20, 30);

  /* いわゆるお決まりのセリフ */
  AutoGL_SetColor (1, 1, 0);
  AutoGL_DrawString (10, 20, 30, "HelloWorld");
}

/* この関数はプログラム起動直後にAutoGLライブラリから呼び出される。*/
void AutoGL_SetUp (int argc, char *argv[]) 
{
  /* サイズや視線方向などを設定する。*/
  AutoGL_SetViewSize (70);                 /* ビューの大きさ */
  AutoGL_SetViewCenter (10, 10, 10);       /* 注視点 */
  AutoGL_SetViewDirection (1, 2, 3);       /* 視線方向 */

  /* おまじない */
  AutoGL_SetViewRedrawCallback (RedrawView); 
  AutoGL_SetMode3D (AUTOGL_MODE_3D_ROTATE_XY);
  AutoGL_SetDefaultCallbackInMode3D (NULL);    
  AutoGL_EnableDragInMode3D ();               
}






