/*  AutoGL ＜お手軽＞ のサンプルプログラム例 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/* AutoGL ”おてがる”を使うには,このヘッダファイルをインクルードする。*/
#include <autogl.h>

/* 二次元差分法解析結果の可視化 */
/* 面コンタ-と三次元グラフを表示する。 */
/* これは、典型的な三次元アプリケーションの例である。 */

/* ユーティリティの可視化（面コンタ-）の利用 */
/* 引数など詳しい情報は, autogl_utility.hを参照。*/
/* なお、これはautogl.hからincludeされている。*/




/* ビュー上に表示されるモデルを表現するための変数群 */

/* 直交差分格子による二次元スカラー場 */

/* 格子の分割数 */
#if 0
#define CELLS 20
#else
#define CELLS 100
#endif

/* 格子上のスカラー値 */
static double GridValues[CELLS + 1][CELLS + 1];

/* コンタ-の最大・最小値 */
static double MinRange = -1.0;
static double MaxRange = 5.0;

/* 三次元グラフのZ方向の拡大率 */
static double ZScale = 5.0;

/* 格子を表示するか */
static int GridIsShown = 0;

/* コンタ-を表示するか */
static int ContourIsShown = 0;

/* 三次元グラフを表示するか */
static int Graph3DIsShown = 0;



/* 二次元スカラー場の定義 */
static double Function (double x, double y)
{
  return sin (x * M_PI) + y * y;
}

/* 二次元格子を初期化する。 */
static void InitializeGrid ()
{
  int i, j;

  /* モデル用変数を初期化する。 */
  for (i = 0; i <= CELLS; i++) {
    for (j = 0; j <= CELLS; j++) {
      double scale = 4.0 / CELLS;
      
      GridValues[i][j] 
	= Function ((i - CELLS / 2) * scale, 
		    (j - CELLS / 2) * scale);
    }
  }
}

/* 格子を描画する。 */
static void DrawGrid (void) 
{
  double cellSize = 60.0 / CELLS;
  int i, j;

  AutoGL_SetColor (1, 1, 1);

  for (i = 0; i < CELLS; i++) {
    for (j = 0; j < CELLS; j++) {
      double x00 = (i - CELLS / 2) * cellSize;
      double y00 = (j - CELLS / 2) * cellSize;
      double x10 = x00 + cellSize;
      double y10 = y00;
      double x11 = x00 + cellSize;
      double y11 = y00 + cellSize;
      double x01 = x00;
      double y01 = y00 + cellSize;

      AutoGL_DrawLine (x00, y00, 0.0,
		       x10, y10, 0.0);
      AutoGL_DrawLine (x10, y10, 0.0,
		       x11, y11, 0.0);
      AutoGL_DrawLine (x11, y11, 0.0,
		       x01, y01, 0.0);
      AutoGL_DrawLine (x01, y01, 0.0,
		       x00, y00, 0.0);
    }
  }
}

/* コンタ-を描画する。 */
static void DrawContour (void) 
{
  double cellSize = 60.0 / CELLS;
  int i, j;

  for (i = 0; i < CELLS; i++) {
    for (j = 0; j < CELLS; j++) {
      double x00 = (i - CELLS / 2) * cellSize;
      double y00 = (j - CELLS / 2) * cellSize;
      double value00 = GridValues[i + 0][j + 0];
      double grade00 = (value00 - MinRange) / (MaxRange - MinRange);
      
      double x10 = x00 + cellSize;
      double y10 = y00;
      double value10 = GridValues[i + 1][j + 0];
      double grade10 = (value10 - MinRange) / (MaxRange - MinRange);
      
      double x11 = x00 + cellSize;
      double y11 = y00 + cellSize;
      double value11 = GridValues[i + 1][j + 1];
      double grade11 = (value11 - MinRange) / (MaxRange - MinRange);
      
      double x01 = x00;
      double y01 = y00 + cellSize;
      double value01 = GridValues[i + 0][j + 1];
      double grade01 = (value01 - MinRange) / (MaxRange - MinRange);

      /* コンタ-つきの三角形を表示する。 */
      AutoGL_DrawContourTriangle (x00, y00, 0.0, grade00,
				  x10, y10, 0.0, grade10,
				  x11, y11, 0.0, grade11);
      AutoGL_DrawContourTriangle (x11, y11, 0.0, grade11,
				  x01, y01, 0.0, grade01,
				  x00, y00, 0.0, grade00);
      /* grade??は 0.0 - 1.0 に規格化されていなければならない。 */
    }
  }
}

/* 三次元グラフを描画する。 */
static void DrawGraph3D (void) 
{
  double cellSize = 60.0 / CELLS;
  int i, j;

  AutoGL_SetColor (1, 1, 1);

  for (i = 0; i < CELLS; i++) {
    for (j = 0; j < CELLS; j++) {
      double x00 = (i - CELLS / 2) * cellSize;
      double y00 = (j - CELLS / 2) * cellSize;
      double value00 = GridValues[i + 0][j + 0];
      
      double x10 = x00 + cellSize;
      double y10 = y00;
      double value10 = GridValues[i + 1][j + 0];
      
      double x11 = x00 + cellSize;
      double y11 = y00 + cellSize;
      double value11 = GridValues[i + 1][j + 1];
      
      double x01 = x00;
      double y01 = y00 + cellSize;
      double value01 = GridValues[i + 0][j + 1];

      AutoGL_DrawTriangle (x00, y00, value00 * ZScale,
			   x10, y10, value10 * ZScale,
			   x11, y11, value11 * ZScale);
      AutoGL_DrawTriangle (x11, y11, value11 * ZScale,
			   x01, y01, value01 * ZScale,
			   x00, y00, value00 * ZScale);
    }
  }

}



/* ビューの再描画のためのコールバック関数 */
/* ビューが再表示されるごとに呼ばれる。 */
static void RedrawView (void) 
{
  if (GridIsShown) {
    DrawGrid ();
  }
  if (ContourIsShown) {
    DrawContour ();
  }
  if (Graph3DIsShown) {
    DrawGraph3D ();
  }
}

/* 描画アップデートコマンドのためのコールバック関数 */
/* "UpdateModel"ボタンを選ぶと呼ばれる。 */
static void UpdateModel (void) 
{
  /* ビューを再描画する。 */
  AutoGL_DrawView ();  
  /* 登録されたビューの再描画関数（この場合はRedrawView ()）が呼ばれる。 */
}



/* 関数AutoGL_SetUpはユーザー側プログラムごとに必ず一つ用意すること。*/
/* ここで、コールバック関数や制御変数などを登録する。*/
void AutoGL_SetUp (int argc, char *argv[]) 
{
  InitializeGrid ();

  /* ビューを準備する。 */

  AutoGL_SetBackgroundColor (0.7, 0.7, 0.7);   /* ビューの背景色 */
  AutoGL_SetViewSize (80); 

  AutoGL_SetViewRedrawCallback (RedrawView);    
  AutoGL_SetDefaultCallbackInMode3D (NULL);  /* 三次元アプリ用イベント処理 */

  /* もし描画が重ければドラッグ中は描画させないほうがよい。 */
  AutoGL_EnableDragInMode3D ();

  AutoGL_SetPanelInMode3D ();       /* 三次元アプリ用イベント処理 */

  AutoGL_AddComment ();
  AutoGL_SetLabel ("Model Parameters");

  /* 格子 */
  AutoGL_AddBoolean (&GridIsShown, "GridIsShown"); 
  AutoGL_SetLabel ("grid");

  /* コンタ- */
  AutoGL_AddBoolean (&ContourIsShown, "ContourIsShown"); 
  AutoGL_SetLabel ("contour");
  AutoGL_AddReal (&MinRange, "MinRange");      
  AutoGL_SetLabel ("min");
  AutoGL_AddReal (&MaxRange, "MaxRange");      
  AutoGL_SetLabel ("max");

  /* 三次元グラフ */
  AutoGL_AddBoolean (&Graph3DIsShown, "Graph3DIsShown"); 
  AutoGL_SetLabel ("3-D graph");
  AutoGL_AddReal (&ZScale, "ZScale");      

  AutoGL_AddCallback (UpdateModel, "UpdateModel");
}


