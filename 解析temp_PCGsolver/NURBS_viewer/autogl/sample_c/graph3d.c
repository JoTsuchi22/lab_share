/*  AutoGL ＜お手軽＞ のサンプルプログラム例 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/* AutoGL ”おてがる”を使うには,このヘッダファイルをインクルードする。*/
#include <autogl.h>

/* 二次元差分法解析結果の可視化 */
/* 三次元グラフを表示する。 */

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

/* 三次元グラフのZ方向の拡大率 */
static double ZScale = 5.0;



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

/* 三次元グラフを描画する。 */
static void PlotGraph3D (void) 
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
  static int modelIsRendered = 0;

  /* まず、一度だけモデルをディスプレイリストへ描画しておく。 */
  if (!modelIsRendered) {

    /* ディスプレイリストを開き、これを空にする。 */
    AutoGL_OpenDisplayList ();
    
    PlotGraph3D ();

    /* ディスプレイリストを閉じる。 */
    AutoGL_CloseDisplayList ();

    modelIsRendered = 1;
  }

  /* 鏡面反射をオンにする。*/
  AutoGL_TurnOnSpecular ();

  /* 画面左下に三次元座標系を描画する。 */
  AutoGL_DrawCoordinateSystem3D (20);

  /* ディスプレイリストに入れた描画命令を実際に描画する。 */
  AutoGL_DrawDisplayList ();
}

/* 関数AutoGL_SetUpはユーザー側プログラムごとに必ず一つ用意すること。*/
/* ここで、コールバック関数や制御変数などを登録する。*/
void AutoGL_SetUp (int argc, char *argv[]) 
{
  InitializeGrid ();

  AutoGL_SetViewSize (60);   /* ビューのサイズ */

  /* おまじない */
  AutoGL_SetViewRedrawCallback (RedrawView);    
  AutoGL_SetMode3D (AUTOGL_MODE_3D_ROTATE_XY); /* マウスで回転 */
  AutoGL_SetDefaultCallbackInMode3D (NULL);
  AutoGL_SetPanelInMode3D ();

  /* もし描画が重ければドラッグ中は描画させないほうがよい。 */
  AutoGL_EnableDragInMode3D ();
}


