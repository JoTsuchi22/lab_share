/*  AutoGL ＜お手軽＞ のサンプルプログラム例 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/* AutoGL ”おてがる”を使うには,このヘッダファイルをインクルードする。*/
#include <autogl.h>

/* AutoNOSim ”おたのしみ”を使うには,このヘッダファイルをインクルードする。*/
#include <autonosim.h>

/* ユーティリティの三次元描画プリミティブ（球、円柱など）の利用 */
/* 引数など詳しい情報は, autogl_utility3d.hを参照。*/
/* なお、これはautogl.hからincludeされている。*/



/* コールバック関数群 */

/* ビューの再描画のためのコールバック関数 */
/* ビューが再表示されるごとに呼ばれる。 */
static void RedrawView (void) 
{
  double radius = 10.0;
  int nDivisions = 4;
  int taskId = AutoNOSim_MyTaskId ();
  
  /* 鏡面反射をオンにする。*/
  AutoGL_TurnOnSpecular ();

  if (taskId == 0) {

    AutoGL_SetColor (1, 1, 1);
    AutoGL_DrawString (0, 0, 0, "O");
    AutoGL_DrawLine (0, 0, 0,
		     50, 0, 0);
    AutoGL_DrawString (50, 0, 0, "X");  
    AutoGL_DrawLine (0, 0, 0,
		     0, 50, 0);
    AutoGL_DrawString (0, 50, 0, "Y");  
    AutoGL_DrawLine (0, 0, 0,
		     0, 0, 50);
    AutoGL_DrawString (0, 0, 50, "Z");  

    /* 球を描画する。*/
    AutoGL_SetColor (1, 1, 0);
    AutoGL_DrawSphere3D (10, 10, 0,
			 radius, nDivisions);
  }

  if (taskId == 1) {

    /* 円柱を描画する。*/
    AutoGL_SetColor (0, 1, 1);
    AutoGL_DrawCylinder3D (-10, -20, 30,
			   10, 20, 30,
			   radius, nDivisions);
  }

  if (taskId == 2) {

    /* 円錐を描画する。*/
    AutoGL_SetColor (1, 0, 1);
    AutoGL_DrawCone3D (-10, -20, -30,
		       10, 20, -30,
		       radius, nDivisions);
  }
  
  if (taskId == 3) {

    /* 円を描画する。*/
    AutoGL_SetColor (1, 0, 0);
    AutoGL_DrawCircle3D (-10, -10, 0,
			 1, 1, 0,
			 radius, nDivisions);
  }

}



/* 関数AutoGL_SetUpはユーザー側プログラムごとに必ず一つ用意すること。*/
/* ここで、コールバック関数や制御変数などを登録する。*/
void AutoGL_SetUp (int argc, char *argv[]) 
{
  /* ビューのみの構成例 */

  /* ビューイングパラメータ（サイズ,視線方向など）を設定する。*/
  AutoGL_SetViewSize (70);                 /* ビューの大きさ */
  AutoGL_SetViewCenter (10, 10, 10);       /* 注視点 */
  AutoGL_SetViewDirection (1, 2, 3);       /* 視線方向 */

  /* 透視投影モードにする。*/
  AutoGL_SetPerspectiveViewFlag (1);

  /* おまじない */
  AutoGL_SetViewRedrawCallback (RedrawView); 
  AutoGL_SetMode3D (AUTOGL_MODE_3D_ROTATE_XY); /* マウスで回転 */
  AutoGL_SetDefaultCallbackInMode3D (NULL);
  AutoGL_SetPanelInMode3D ();
#if 0
  AutoGL_EnableDragInMode3D ();                /* ドラッグ有効 */
#endif
}






