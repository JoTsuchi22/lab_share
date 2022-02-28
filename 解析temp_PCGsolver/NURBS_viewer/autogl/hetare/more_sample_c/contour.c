/*  AutoGL ＜お手軽＞ のサンプルプログラム例 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/* AutoGL ”おてがる”を使うには,このヘッダファイルをインクルードする。*/
#include <autogl.h>

/* 二次元差分法解析結果の可視化 */
/* 面コンタ-を表示する。 */

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

/* コンタ-を描画する。 */
static void PlotContour (void) 
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



/* ビューの再描画のためのコールバック関数 */
/* ビューが再表示されるごとに呼ばれる。 */
static void RedrawView (void) 
{
  PlotContour ();
}

/* 関数AutoGL_SetUpはユーザー側プログラムごとに必ず一つ用意すること。*/
/* ここで、コールバック関数や制御変数などを登録する。*/
void AutoGL_SetUp (int argc, char *argv[]) 
{
  InitializeGrid ();

  AutoGL_SetViewSize (60);   /* ビューのサイズ */

  /* おまじない */
  AutoGL_SetViewRedrawCallback (RedrawView);    
  AutoGL_SetMode3D (AUTOGL_MODE_3D_SCALE);     /* マウスで拡大縮小 */
  AutoGL_SetDefaultCallbackInMode3D (NULL);

#if 0
  /* もし描画が重ければドラッグ中は描画させないほうがよい。 */
  AutoGL_EnableDragInMode3D ();
#endif
}


