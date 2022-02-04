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

/* ユーティリティの三次元イベント処理の利用 */
/* 引数など詳しい情報は, autogl_utility3d.hを参照。*/
/* なお、これはautogl.hからincludeされている。*/



/* ビュー上に表示されるモデルを表現するための変数群 */

/* 現在選択されている頂点の番号 */
static int SelectedFaceId = -1;




/* ビューの再描画のためのコールバック関数 */
/* ビューが再表示されるごとに呼ばれる。 */
static void RedrawView (void) 
{
  int taskId = AutoNOSim_MyTaskId ();

  if (taskId % 4 == 0) {
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
  }

  if (taskId % 4 == 1) {
    if (SelectedFaceId == 0) {
      AutoGL_SetColor (1, 0, 0);
    } else {
      AutoGL_SetColor (1, 1, 1);
    }
    AutoGL_DrawQuadrangle (0, 0, 0,
			   10, 0, 0,
			   10, 20, 0,
			   0, 20, 0);

    if (SelectedFaceId == 1) {
      AutoGL_SetColor (1, 0, 0);
    } else {
      AutoGL_SetColor (1, 1, 1);
    }
    AutoGL_DrawQuadrangle (0, 0, 30,
			   10, 0, 30,
			   10, 20, 30,
			   0, 20, 30);
  }
  if (taskId % 4 == 2) {
    if (SelectedFaceId == 2) {
      AutoGL_SetColor (1, 0, 0);
    } else {
      AutoGL_SetColor (1, 1, 1);
    }
    AutoGL_DrawQuadrangle (0, 0, 0,
			   0, 20, 0,
			   0, 20, 30,
			   0, 0, 30);

    if (SelectedFaceId == 3) {
      AutoGL_SetColor (1, 0, 0);
    } else {
      AutoGL_SetColor (1, 1, 1);
    }
    AutoGL_DrawQuadrangle (10, 0, 0,
			   10, 20, 0,
			   10, 20, 30,
			   10, 0, 30);
  }
  if (taskId % 4 == 3) {
    if (SelectedFaceId == 4) {
      AutoGL_SetColor (1, 0, 0);
    } else {
      AutoGL_SetColor (1, 1, 1);
    }
    AutoGL_DrawQuadrangle (0, 0, 0,
			   10, 0, 0,
			   10, 0, 30,
			   0, 0, 30);

    if (SelectedFaceId == 5) {
      AutoGL_SetColor (1, 0, 0);
    } else {
      AutoGL_SetColor (1, 1, 1);
    }
    AutoGL_DrawQuadrangle (0, 20, 0,
			   10, 20, 0,
			   10, 20, 30,
			   0, 20, 30);
  }

  if (taskId % 4 == 0) {
    AutoGL_SetColor (0, 0, 1);
    AutoGL_DrawString (5, 10, -2, "0");  
    AutoGL_DrawString (5, 10, 32, "1");  
    AutoGL_DrawString (-2, 10, 15, "2");  
    AutoGL_DrawString (12, 10, 15, "3");  
    AutoGL_DrawString (5, -2, 15, "4");  
    AutoGL_DrawString (5, 22, 15, "5");  
  }
}

/* ビューのイベント処理ためのコールバック関数 */
/* ビューにイベントが発生するごとに呼ばれる。 */
static void HandleEvent (void)
{
  int event = AutoGL_GetViewEvent ();
  int taskId = AutoNOSim_MyTaskId ();

  /* もし、マウスクリックがあれば */
  if (event == AUTOGL_EVENT_POINTING_DEVICE_RELEASE) {
    int dcX, dcY, dcZ;
    int toleranceDc = 10;

    AutoGL_GetPointingDevicePositionDc 
      (&dcX, &dcY, &dcZ);

    AutoGL_StartSelection (dcX, dcY, toleranceDc);
    
    if (taskId % 4 == 1) {
      AutoGL_SetSelectionId (0);
      AutoGL_DrawQuadrangle (0, 0, 0,
			     10, 0, 0,
			     10, 20, 0,
			     0, 20, 0);
      AutoGL_SetSelectionId (1);
      AutoGL_DrawQuadrangle (0, 0, 30,
			     10, 0, 30,
			     10, 20, 30,
			     0, 20, 30);
    }
    if (taskId % 4 == 2) {
      AutoGL_SetSelectionId (2);
      AutoGL_DrawQuadrangle (0, 0, 0,
			     0, 20, 0,
			     0, 20, 30,
			     0, 0, 30);
      AutoGL_SetSelectionId (3);
      AutoGL_DrawQuadrangle (10, 0, 0,
			     10, 20, 0,
			     10, 20, 30,
			     10, 0, 30);
    }
    if (taskId % 4 == 3) {
      AutoGL_SetSelectionId (4);
      AutoGL_DrawQuadrangle (0, 0, 0,
			     10, 0, 0,
			     10, 0, 30,
			     0, 0, 30);
      AutoGL_SetSelectionId (5);
      AutoGL_DrawQuadrangle (0, 20, 0,
			     10, 20, 0,
			     10, 20, 30,
			     0, 20, 30);
    }
      
    AutoGL_EndSelection ();
    SelectedFaceId = AutoGL_GetSelectedId ();

    fprintf (stderr, " selected face ID is %d \n", SelectedFaceId);

    /* ビューを再描画する。 */
    AutoGL_DrawView ();
    /* 登録されたビューの再描画関数（この場合はRedrawView ()）が呼ばれる。 */

  } else if (event == AUTOGL_EVENT_KEY) {
    int keyChar = AutoGL_GetKeyChar ();

    if (keyChar == ' ') {
      fprintf (stderr, " hit space \n");
    } else {
      AutoGL_HandleDefaultKeyEventInMode3D ();
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

