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
/* 引数など詳しい情報は, autogl_visualization.hを参照。*/
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
static int NGrades = 20;

/* 三次元グラフのZ方向の拡大率 */
static double ZScale = 5.0;

/* 格子を表示するか */
static int GridIsShown = 0;

/* コンタ-を表示するか */
static int ContourIsShown = 0;

/* コンタ-のタイプ : 0 - 2 の整数 (CONTOUR_) */
#define CONTOUR_SMOOTH 0
#define CONTOUR_LINE 1
#define CONTOUR_BAND 2
static int ContourType = 0;

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
      
      double x10 = x00 + cellSize;
      double y10 = y00;
      double value10 = GridValues[i + 1][j + 0];
      
      double x11 = x00 + cellSize;
      double y11 = y00 + cellSize;
      double value11 = GridValues[i + 1][j + 1];
      
      double x01 = x00;
      double y01 = y00 + cellSize;
      double value01 = GridValues[i + 0][j + 1];

      switch (ContourType) {
      case CONTOUR_SMOOTH:
	/* スムーズコンタ-つきの四辺形を表示する。 */
	AutoGL_DrawSmoothContourQuadrangle 
	  (MinRange, MaxRange,
	   x00, y00, 0.0, value00,
	   x10, y10, 0.0, value10,
	   x11, y11, 0.0, value11,
	   x01, y01, 0.0, value01);
	break;
      case CONTOUR_LINE:
	/* 線コンタ-つきの四辺形を表示する。 */
	AutoGL_DrawLineContourQuadrangle 
	  (MinRange, MaxRange, NGrades,
	   x00, y00, 0.0, value00,
	   x10, y10, 0.0, value10,
	   x11, y11, 0.0, value11,
	   x01, y01, 0.0, value01);
	break;
      case CONTOUR_BAND:
	/* バンドコンタ-つきの四辺形を表示する。 */
	AutoGL_DrawBandContourQuadrangle 
	  (MinRange, MaxRange, NGrades,
	   x00, y00, 0.0, value00,
	   x10, y10, 0.0, value10,
	   x11, y11, 0.0, value11,
	   x01, y01, 0.0, value01);
	break;
      }
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
  /* 画面左下に三次元座標系を描画する。 */
  AutoGL_DrawCoordinateSystem3D (20);

  if (ContourIsShown) {
    /* 画面右側にコンタ-のバーを描画する。 */
    switch (ContourType) {
    case CONTOUR_SMOOTH:
      AutoGL_DrawContourMap (80, MinRange, MaxRange);
      break;
    case CONTOUR_LINE:
      AutoGL_DrawLineContourMap (80, MinRange, MaxRange, NGrades);
      break;
    case CONTOUR_BAND:
      AutoGL_DrawBandContourMap (80, MinRange, MaxRange, NGrades);
      break;
    }
  }

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

  if (GridIsShown) {
    DrawGrid ();
  }
  if (ContourIsShown) {
    DrawContour ();
  }
  if (Graph3DIsShown) {
    DrawGraph3D ();
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

  AutoGL_SetBackgroundColor (0.7, 0.7, 0.7);   /* ビューの背景色 */
  AutoGL_SetViewSize (80); 

  AutoGL_SetViewRedrawCallback (RedrawView);    
  AutoGL_SetDefaultCallbackInMode3D (NULL);  /* 三次元アプリ用イベント処理 */

  /* もし描画が重ければドラッグ中は描画させないほうがよい。 */
  AutoGL_EnableDragInMode3D ();

#if 1
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
  AutoGL_SetPanelInMode3D ();       /* 三次元アプリ用イベント処理 */

  AutoGL_AddComment ();
  AutoGL_SetLabel ("Model Parameters");

  /* 格子 */
  AutoGL_AddBoolean (&GridIsShown, "GridIsShown"); 
  AutoGL_SetLabel ("grid");

  /* コンタ- */
  AutoGL_AddBoolean (&ContourIsShown, "ContourIsShown"); 
  AutoGL_SetLabel ("contour");
  AutoGL_AddInteger (&ContourType, "ContourType");
  AutoGL_SetLabel ("type");
  AutoGL_AddIntegerItem ("smooth");    
  AutoGL_AddIntegerItem ("line");
  AutoGL_AddIntegerItem ("band");     
  AutoGL_AddReal (&MinRange, "MinRange");      
  AutoGL_SetLabel ("min");
  AutoGL_AddReal (&MaxRange, "MaxRange");      
  AutoGL_SetLabel ("max");
  AutoGL_AddInteger (&NGrades, "NGrades");      
  AutoGL_SetLabel ("num. grades");

  /* 三次元グラフ */
  AutoGL_AddBoolean (&Graph3DIsShown, "Graph3DIsShown"); 
  AutoGL_SetLabel ("3-D graph");
  AutoGL_AddReal (&ZScale, "ZScale");      

  AutoGL_AddCallback (UpdateModel, "UpdateModel");
}


