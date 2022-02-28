/*  AutoGL ＜お手軽＞ のサンプルプログラム例 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/* AutoGL ”おてがる”を使うには,このヘッダファイルをインクルードする。*/
#include <autogl.h>

/* ユーティリティの三次元イベント処理の利用 */
/* 引数など詳しい情報は, autogl_utility3d.hを参照。*/
/* なお、これはautogl.hからincludeされている。*/




/* ビュー上に表示されるモデルを表現するための変数群 */

/* 現在選択されている頂点の番号 */
static int SelectedVertexId = -1;

/* 頂点の座標 */
static double Vertexs[8][3] = {
  { 0, 0, 0 },
  { 10, 0, 0 },
  { 10, 20, 0 },
  { 0, 20, 0 },
  { 0, 0, 30 },
  { 10, 0, 30 },
  { 10, 20, 30 },
  { 0, 20, 30 }
};



/* ビューの再描画のためのコールバック関数 */
/* ビューが再表示されるごとに呼ばれる。 */
static void RedrawView (void) 
{
  int i;
  int markSizeDc;
  int isFilled;

  AutoGL_SetColor (0, 1, 0);  

  AutoGL_DrawLine (-50, 0, 0,
		   50, 0, 0); 
  AutoGL_DrawString (50, 0, 0, "X = 5O");  
  AutoGL_DrawString (-50, 0, 0, "X = -5O");  
  
  AutoGL_DrawLine (0, -50, 0,
		   0, 50, 0); 
  AutoGL_DrawString (0, 50, 0, "Y = 5O");  
  AutoGL_DrawString (0, -50, 0, "Y = -5O");  
    
  AutoGL_DrawLine (0, 0, -50,
		   0, 0, 50); 
  AutoGL_DrawString (0, 0, 50, "Z = 5O");  
  AutoGL_DrawString (0, 0, -50, "Z = -5O");  
    
  AutoGL_SetColor (1, 1, 1);
  AutoGL_DrawString (0, 0, 0, "O"); 

  AutoGL_SetColor (1, 1, 1);
  AutoGL_DrawQuadrangle (0, 0, 0,
			 10, 0, 0,
			 10, 20, 0,
			 0, 20, 0);
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

  /* 各頂点について、それが選択されているかどうかを描画する。 */
  AutoGL_SetColor (1, 0, 0);
  for (i = 0; i < 8; i++) {
    isFilled = 0;
    if (i == SelectedVertexId) {
      isFilled = 1;
    }
    markSizeDc = 10;

    /* ダイアモンドのマークを表示する。 */
    AutoGL_DrawDiamondMark3D 
      (Vertexs[i][0], Vertexs[i][1], Vertexs[i][2], 
       markSizeDc , isFilled);
    /* markSizeDcにより、10 pixelぐらいの大きさとなる。 */
    /* isFilled == 1なら、この頂点は選択されている。 */
  }
}

/* ビューのイベント処理ためのコールバック関数 */
/* ビューにイベントが発生するごとに呼ばれる。 */
static void HandleEvent (void)
{
  int event = AutoGL_GetViewEvent ();

  /* もし、マウスクリックがあれば */
  if (event == AUTOGL_EVENT_POINTING_DEVICE_RELEASE) {
    int i;

    /* 各頂点について、マウスポインタがヒットしたかどうか調べる。 */
    SelectedVertexId = -1;
    for (i = 0; i < 8; i++) {
      int toleranceDc = 10;
      double parameter;

      /* この頂点がマウスポインタのそばにあれば */
      if (AutoGL_PointingDeviceIsHit3D 
	  (&parameter,
	   Vertexs[i][0], Vertexs[i][1], Vertexs[i][2], 
	   toleranceDc)) {
	/* 点の座標は世界座標（三次元）で与える。 */
	SelectedVertexId = i;
	break;
      }
    }
    if (SelectedVertexId != -1) {
      printf (" vertex %d is selected. \n",
	      SelectedVertexId);
    }

    /* ビューを再描画する。 */
    AutoGL_DrawView ();
    /* 登録されたビューの再描画関数（この場合はRedrawView ()）が呼ばれる。 */

  } else if (event == AUTOGL_EVENT_KEY) {
    AutoGL_HandleDefaultKeyEventInMode3D ();
  }
}



/* 関数AutoGL_SetUpはユーザー側プログラムごとに必ず一つ用意すること。*/
/* ここで、コールバック関数や制御変数などを登録する。*/
void AutoGL_SetUp (int argc, char *argv[]) 
{
  AutoGL_SetViewSize (70);

  /* ビューの再描画コールバックの登録 */
  AutoGL_SetViewRedrawCallback (RedrawView);       

  /* ビューのイベント処理コールバックの登録 */
  AutoGL_SetDefaultCallbackInMode3D (HandleEvent);
  /* イベント処理コールバック関数として上にあるHandleEvent ()を指定。 */
  /* 特別なイベント処理が要らなければ、NULLでよい。*/

  /* 初期状態では、ユーザー定義イベント処理を行う。 */
  AutoGL_SetMode3D (AUTOGL_MODE_3D_OTHER); 

#if 0
  /* イベント処理中, ドラッグを有効にする。*/
  AutoGL_EnableDragInMode3D ();
#endif

  /* パネル上に三次元アプリ用イベント処理機能を準備する。 */
  AutoGL_SetPanelInMode3D ();
}

