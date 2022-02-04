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

/* ユーティリティの二次元イベント処理の利用 */
/* 引数など詳しい情報は, autogl_utility2d.hを参照。*/
/* なお、これはautogl.hからincludeされている。*/



/* ビュー上に表示されるモデルを表現するための変数群 */

/* 現在選択されている頂点の番号 */
static int SelectedEdgeId = -1;

/* 頂点の座標 */
static double Vertexs[6][2] = {
  { 10, 0 },
  { 10, 20 },
  { 0, 20 },
  { -10, 0 },
  { 0, -20 },
  { 30, -20 }
};



/* ビューの再描画のためのコールバック関数 */
/* ビューが再表示されるごとに呼ばれる。 */
static void RedrawView (void) 
{
  int i;
  int nTasks = AutoNOSim_NTasks ();
  int taskId = AutoNOSim_MyTaskId ();

  AutoGL_SetColor (0, 1, 0);  

  AutoGL_DrawLine (-50, 0, 0,
		   50, 0, 0); 
  AutoGL_DrawString (50, 0, 0, "X = 5O");  
  AutoGL_DrawString (-50, 0, 0, "X = -5O");  
  
  AutoGL_DrawLine (0, -50, 0,
		   0, 50, 0); 
  AutoGL_DrawString (0, 50, 0, "Y = 5O");  
  AutoGL_DrawString (0, -50, 0, "Y = -5O");  
    
  AutoGL_SetColor (1, 1, 1);
  AutoGL_DrawString (0, 0, 0, "O"); 

  for (i = 0; i < 6 - 1; i++) {
    if (i % nTasks != taskId) continue;

    if (i == SelectedEdgeId) {
      AutoGL_SetColor (1, 0, 0);
    } else {
      AutoGL_SetColor (1, 1, 1);
    }
    AutoGL_DrawLine (Vertexs[i][0], Vertexs[i][1], 0,
		     Vertexs[i + 1][0], Vertexs[i + 1][1], 0);
  }
}

/* ビューのイベント処理ためのコールバック関数 */
/* ビューにイベントが発生するごとに呼ばれる。 */
static void HandleEvent (void)
{
  int event = AutoGL_GetViewEvent ();
  int nTasks = AutoNOSim_NTasks ();
  int taskId = AutoNOSim_MyTaskId ();

  /* もし、マウスクリックがあれば */
  if (event == AUTOGL_EVENT_POINTING_DEVICE_RELEASE) {
    int dcX, dcY, dcZ;
    int toleranceDc = 10;
    int i;

    AutoGL_GetPointingDevicePositionDc 
      (&dcX, &dcY, &dcZ);

    AutoGL_StartSelection (dcX, dcY, toleranceDc);
    
    for (i = 0; i < 6 - 1; i++) {
      if (i % nTasks != taskId) continue;

      AutoGL_SetSelectionId (i);
      AutoGL_DrawLine (Vertexs[i][0], Vertexs[i][1], 0,
		       Vertexs[i + 1][0], Vertexs[i + 1][1], 0);
    }

    AutoGL_EndSelection ();
    SelectedEdgeId = AutoGL_GetSelectedId ();

    fprintf (stderr, " selected edge ID is %d \n", SelectedEdgeId);

    /* ビューを再描画する。 */
    AutoGL_DrawView ();
    /* 登録されたビューの再描画関数（この場合はRedrawView ()）が呼ばれる。 */

  } else if (event == AUTOGL_EVENT_KEY) {
    int keyChar = AutoGL_GetKeyChar ();

    if (keyChar == ' ') {
      fprintf (stderr, " hit space \n");
    } else {
      AutoGL_HandleDefaultKeyEventInMode2D ();
    }
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
  AutoGL_SetDefaultCallbackInMode2D (HandleEvent);
  /* イベント処理コールバック関数として上にあるHandleEvent ()を指定。 */
  /* 特別なイベント処理が要らなければ、NULLでよい。*/

  /* 初期状態では、ユーザー定義イベント処理を行う。 */
  AutoGL_SetMode2D (AUTOGL_MODE_2D_OTHER); 

#if 0
  /* イベント処理中, ドラッグを有効にする。*/
  AutoGL_EnableDragInMode2D ();
#endif

  /* パネル上に二次元アプリ用イベント処理機能を準備する。 */
  AutoGL_SetPanelInMode2D ();
}

