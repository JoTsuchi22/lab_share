/*  AutoGL ＜お手軽＞ のサンプルプログラム例 */

/* ユーティリティの可視化機能の利用 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/* AutoGL ”おてがる”を使うには,このヘッダファイルをインクルードする。*/
#include <autogl.h>

/* 三次元差分法解析結果の可視化 */
/* 等値面と断面を表示する。 */
/* これは、典型的な三次元アプリケーションの例である。 */

/* ユーティリティの可視化（等値面,断面）の利用 */
/* 引数など詳しい情報は, autogl_visualization.hを参照。*/
/* なお、これはautogl.hからincludeされている。*/



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

/* 真ん中の等値面の値 */
static double LevelValue = 2.0;
/* 幅（他の二枚の等値面は、LevelValue +- LevelWidth となる。） */
static double LevelWidth = 2.0;

/* 断面の方程式 Ax+By+Cz+D=0 の係数 */
static double A = 1.0;
static double B = 1.0;
static double C = 1.0;
static double D = 0.0;

/* 断面を表示するか */
static int SectionIsShown = 0;

/* 等値面を表示するか */
static int IsosurfaceIsShown = 0;



/* 三次元スカラー場の定義 */
static double Function (double x, double y, double z)
{
  return sin (x * M_PI) + y * y + z * z;
}

/* 三次元格子を初期化する。 */
static void InitializeGrid ()
{
  int i, j, k;

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

/* 格子の等値面を描画する。 */
static void PlotIsosurfaceOfLevel (double levelValue)
{
  double cellSize = HALF_GRID_SIZE * 2 / CELLS;
  double grade;
  double red, green, blue;
  int gx, gy, gz;

  /* 0.0 - 1.0 のグレード値に対応する色を得る。 */
  grade = (levelValue - MinRange) / (MaxRange - MinRange);
  AutoGL_GetContourColor (&red, &green, &blue,
			  grade);
  AutoGL_SetColor (red, green, blue);

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

	/* 六面体の等値面を描画する。*/
	AutoGL_DrawIsosurfaceHexahedron 
	  (levelValue, levelValue, 0,
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

/* 等値面を描画する。 */
static void PlotIsosurface (void)
{
  /* 三枚の等値面を描画する。 */
  PlotIsosurfaceOfLevel (LevelValue - LevelWidth);
  PlotIsosurfaceOfLevel (LevelValue);
  PlotIsosurfaceOfLevel (LevelValue + LevelWidth);
}



/* ビューの再描画のためのコールバック関数 */
/* ビューが再表示されるごとに呼ばれる。 */
static void RedrawView (void) 
{
  /* 画面左下に三次元座標系を描画する。 */
  AutoGL_DrawCoordinateSystem3D (20);

  if (SectionIsShown) {
    /* 画面右側にコンタ-のバーを描画する。 */
    AutoGL_DrawContourMap (80, MinRange, MaxRange);
  }

  /* 格子の枠を描画する。 */
  {
    double size = HALF_GRID_SIZE;

    AutoGL_SetColor (1, 1, 1);
    AutoGL_DrawBox3D (-size, -size, -size,
		      size, size, size);
  }

  /* 半透明表示を有効にする。 */
  AutoGL_TurnOnTransparency ();

  /* ディスプレイリストに入れた描画命令を実際に描画する。 */
  AutoGL_DrawDisplayList ();
}

/* 描画アップデートコマンドのためのコールバック関数 */
/* "UpdateModel"ボタンを選ぶと呼ばれる。 */
static void UpdateModel (void) 
{
  /* ディスプレイリストを開き、これを空にする。 */
  AutoGL_OpenDisplayList ();

  /* 以下の描画命令をディスプレイリストに順次入れていく。 */
  /* 実際の描画はここではまだ行われない。 */

  if (SectionIsShown) {
    PlotSection ();
  }

  if (IsosurfaceIsShown) {
    PlotIsosurface ();
  }

  /* ディスプレイリストを閉じる。 */
  AutoGL_CloseDisplayList ();

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

  AutoGL_SetViewSize (80); 

  AutoGL_SetViewRedrawCallback (RedrawView);    
  AutoGL_SetDefaultCallbackInMode3D (NULL);   /* 三次元アプリ用イベント処理 */

  /* もし描画が重ければドラッグ中は描画させないほうがよい。 */
  AutoGL_EnableDragInMode3D ();

#if 0
  /* パネルviewingを準備する。 */
  /* ユーザーが普段用いない比較的マイナーな機能。 */

  AutoGL_AddGroup ("viewing");    
  AutoGL_SetPanelInView3D ();        /* 三次元アプリ用ビューイング変数 */
  AutoGL_SetPanelForInterpreter ();  /* コマンドライン処理 */
  AutoGL_SetPanelForSave ();         /* 画面セーブ機能 */
#endif

  /* パネルmodelを準備する。 */
  /* ユーザーが頻繁に用いる機能。 */
  /* および、各プログラムごとに特有な機能。 */

  AutoGL_AddGroup ("model");    
  AutoGL_SetPanelInMode3D ();
  AutoGL_AddComment ();
  AutoGL_SetLabel ("Model Parameters");

  /* 断面コンタ- */
  AutoGL_AddBoolean (&SectionIsShown, "SectionIsShown"); 
  AutoGL_SetLabel ("section");
  AutoGL_AddComment ();
  AutoGL_SetLabel ("equation Ax+By+Cz+D=0");
  AutoGL_AddReal (&A, "A");      
  AutoGL_AddReal (&B, "B");      
  AutoGL_AddReal (&C, "C");      
  AutoGL_AddReal (&D, "D");      
  AutoGL_AddReal (&MinRange, "MinRange");      
  AutoGL_AddReal (&MaxRange, "MaxRange");      

  /* 等値面 */
  AutoGL_AddBoolean (&IsosurfaceIsShown, "IsosurfaceIsShown"); 
  AutoGL_SetLabel ("iso-surface");
  AutoGL_AddReal (&LevelValue, "LevelValue");      
  AutoGL_AddReal (&LevelWidth, "LevelWidth");      

  AutoGL_AddCallback (UpdateModel, "UpdateModel");
}
