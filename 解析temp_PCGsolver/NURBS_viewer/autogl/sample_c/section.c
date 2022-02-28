/*  AutoGL ＜お手軽＞ のサンプルプログラム例 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/* AutoGL ”おてがる”を使うには,このヘッダファイルをインクルードする。*/
#include <autogl.h>





/* ビュー上に表示されるモデルを表現するための変数群 */

/* 直交差分格子による三次元スカラー場 */

/* 格子の大きさ（半分） */
#define HALF_GRID_SIZE 30.0

/* 格子の分割数 */
#define CELLS 40

/* 格子上のスカラー値 */
static double GridValues[CELLS + 1][CELLS + 1][CELLS + 1];

/* コンタ-の最大・最小値 */
static double MinRange = -1.0;
static double MaxRange = 5.0;

/* 断面の方程式 Ax+By+Cz+D=0 の係数 */
static double A = 1.0;
static double B = 2.0;
static double C = 3.0;
static double D = 40.0;



/* 三次元スカラー場の定義 */
static double Function (double x, double y, double z)
{
  return sin (x * M_PI) + y * y + z * z;
}

/* 三次元格子を初期化する。 */
static void InitializeGrid ()
{
  int i, j, k;

  /* モデル用変数を初期化する。 */
  for (i = 0; i <= CELLS; i++) {
    for (j = 0; j <= CELLS; j++) {
      for (k = 0; k <= CELLS; k++) {
	double scale = 4.0 / CELLS;
	
	GridValues[i][j][k] 
	  = Function ((i - CELLS / 2) * scale, 
		      (j - CELLS / 2) * scale,
		      (k - CELLS / 2) * scale);
      }
    }
  }
}

/* 格子の断面コンタ-を描画する。 */
static void PlotSection (void)
{
  double cellSize = HALF_GRID_SIZE * 2 / CELLS;
  int gx, gy, gz;

  for (gx = 0; gx < CELLS; gx++) {
    for (gy = 0; gy < CELLS; gy++) {
      for (gz = 0; gz < CELLS; gz++) {
	double x000 = (gx - CELLS / 2) * cellSize;
	double y000 = (gy - CELLS / 2) * cellSize;
	double z000 = (gz - CELLS / 2) * cellSize;
	double value000 = GridValues[gx][gy][gz];
	
	double x100 = x000 + cellSize;
	double y100 = y000;
	double z100 = z000;
	double value100 = GridValues[gx + 1][gy][gz];

	double x110 = x000 + cellSize;
	double y110 = y000 + cellSize;
	double z110 = z000;
	double value110 = GridValues[gx + 1][gy + 1][gz];

	double x010 = x000;
	double y010 = y000 + cellSize;
	double z010 = z000;
	double value010 = GridValues[gx][gy + 1][gz];
	
	double x001 = x000;
	double y001 = y000;
	double z001 = z000 + cellSize;
	double value001 = GridValues[gx][gy][gz + 1];

	double x101 = x000 + cellSize;
	double y101 = y000;
	double z101 = z000 + cellSize;
	double value101 = GridValues[gx + 1][gy][gz + 1];

	double x111 = x000 + cellSize;
	double y111 = y000 + cellSize;
	double z111 = z000 + cellSize;
	double value111 = GridValues[gx + 1][gy + 1][gz + 1];
	
	double x011 = x000;
	double y011 = y000 + cellSize;
	double z011 = z000 + cellSize;
	double value011 = GridValues[gx][gy + 1][gz + 1];

	/* 六面体の断面コンタ-を描画する。*/
	AutoGL_DrawSectionSmoothContourHexahedron 
	  (A, B, C, D, MinRange, MaxRange,
	   x000, y000, z000, value000,
	   x100, y100, z100, value100,
	   x110, y110, z110, value110,
	   x010, y010, z010, value010,
	   x001, y001, z001, value001,
	   x101, y101, z101, value101,
	   x111, y111, z111, value111,
	   x011, y011, z011, value011);
      }
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

    PlotSection ();

    /* ディスプレイリストを閉じる。 */
    AutoGL_CloseDisplayList ();

    modelIsRendered = 1;
  }

  /* 画面左下に三次元座標系を描画する。 */
  AutoGL_DrawCoordinateSystem3D (20);

  /* 画面右側にコンタ-のバーを描画する。 */
  AutoGL_DrawContourMap (80, MinRange, MaxRange);

  /* 格子の枠を描画する。 */
  {
    double size = HALF_GRID_SIZE;

    AutoGL_SetColor (1, 1, 1);
    AutoGL_DrawBox3D (-size, -size, -size,
		      size, size, size);
  }

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


