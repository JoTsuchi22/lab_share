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

/* 二次元有限要素法解析結果の可視化 */
/* 線コンタ-,面コンタ-と変形図を表示する。 */
/* これは、典型的な二次元アプリケーションの例である。 */

/* ユーティリティの可視化（線コンタ-）の利用 */
/* 引数など詳しい情報は, autogl_utility3d.hを参照。*/
/* なお、これはautogl.hからincludeされている。*/



/* ビュー上に表示されるモデルを表現するための変数群 */

/* 二次元三角形要素メッシュを用いた有限要素法解析 */

#define MAX_NODES 100
#define MAX_ELEMENTS 100

/* 節点座標 */
int NNodes = 12;
double Coord[MAX_NODES][2] = {
  {  0.0,  0.0 }, { 10.0,  0.0 }, { 20.0,  0.0 }, { 30.0,  0.0 }, 
  {  0.0, 10.0 }, { 10.0, 10.0 }, { 20.0, 10.0 }, { 30.0, 10.0 }, 
  {  0.0, 20.0 }, { 10.0, 20.0 }, { 20.0, 20.0 }, { 30.0, 20.0 }
};

/* 要素コネクティビティ */
int NElements = 12;
int Nop[MAX_ELEMENTS][3] = {
  { 0, 1, 4 }, { 1, 5, 4 }, 
  { 1, 2, 5 }, { 2, 6, 5 }, 
  { 2, 3, 6 }, { 3, 7, 6 }, 
  { 4, 5, 8 }, { 5, 9, 8 },
  { 5, 6, 9 }, { 6, 10, 9 },
  { 6, 7, 10 }, { 7, 11, 10 }
};

/* 適当なスカラー値, 温度とか, 応力とか */
double Scalar[MAX_NODES] = {
    0,  200,  400,   500,
  200,  500,  700,   900,
  400,  600,  900,  1000
};

/* 変位値 */
double Displacement[MAX_NODES][2] = {
  {  0.0,  0.0 }, {  0.1,  0.05 }, {  0.2,  0.15 }, {  0.3,  0.3 }, 
  {  0.0,  0.0 }, {  0.0,  0.05 }, {  0.0,  0.15 }, {  0.0,  0.3 }, 
  {  0.0,  0.0 }, { -0.1,  0.05 }, { -0.2,  0.15 }, { -0.3,  0.3 }
};

/* コンタ-の最大・最小値 */
double MinRange = 0.0;
double MaxRange = 1000.0;

/* 変形図の拡大率 */
double Scale = 3.0;

/* メッシュを表示するか */
int MeshIsShown = 1;

/* なにも表示しない */
#define PLOT_NONE 0
/* スカラー値の面コンタ- */
#define PLOT_CONTOUR 1
/* スカラー値の線コンタ- */
#define PLOT_LINE_CONTOUR 2
/* 変形図 */
#define PLOT_DEFORMATION 3

/* 表示されるプロットのタイプ */
int PlotType = PLOT_NONE;



/* メッシュを描画する。 */
static void PlotMesh (void) 
{
  int iElement;
  int taskId = AutoNOSim_MyTaskId ();
  int nTasks = AutoNOSim_NTasks ();

  AutoGL_SetColor (1, 1, 1);
  for (iElement = 0; iElement < NElements; iElement++) {
    int nodeId0 = Nop[iElement][0];
    int nodeId1 = Nop[iElement][1];
    int nodeId2 = Nop[iElement][2];

    if (iElement % nTasks != taskId) continue;

    AutoGL_DrawLine (Coord[nodeId0][0], Coord[nodeId0][1], 0,
		     Coord[nodeId1][0], Coord[nodeId1][1], 0);
    AutoGL_DrawLine (Coord[nodeId1][0], Coord[nodeId1][1], 0,
		     Coord[nodeId2][0], Coord[nodeId2][1], 0);
    AutoGL_DrawLine (Coord[nodeId2][0], Coord[nodeId2][1], 0,
		     Coord[nodeId0][0], Coord[nodeId0][1], 0);
  }
}

/* スカラー値の面コンタ-を描画する。 */
static void PlotContour (void) 
{
  int iElement;
  int taskId = AutoNOSim_MyTaskId ();
  int nTasks = AutoNOSim_NTasks ();

  AutoGL_DrawContourMap (80, MinRange, MaxRange);

  AutoGL_SetColor (1, 1, 1);
  for (iElement = 0; iElement < NElements; iElement++) {
    int nodeId0 = Nop[iElement][0];
    int nodeId1 = Nop[iElement][1];
    int nodeId2 = Nop[iElement][2];
    double grades[3];

    if (iElement % nTasks != taskId) continue;

    grades[0] = (Scalar[nodeId0] - MinRange) / (MaxRange - MinRange);
    grades[1] = (Scalar[nodeId1] - MinRange) / (MaxRange - MinRange);
    grades[2] = (Scalar[nodeId2] - MinRange) / (MaxRange - MinRange);

    AutoGL_DrawContourTriangle 
      (Coord[nodeId0][0], Coord[nodeId0][1], 0, grades[0],
       Coord[nodeId1][0], Coord[nodeId1][1], 0, grades[1], 
       Coord[nodeId2][0], Coord[nodeId2][1], 0, grades[2]);
  }
}

/* スカラー値の線コンタ-（等高線, 等値線）を描画する。 */
static void PlotLineContour (void) 
{
  int nGrades = 20;
  int iGrade;
  int taskId = AutoNOSim_MyTaskId ();
  int nTasks = AutoNOSim_NTasks ();

  AutoGL_DrawLineContourMap (80, MinRange, MaxRange, nGrades);

  for (iGrade = 1; iGrade < nGrades; iGrade++) {
    double grade = 1.0 / nGrades * iGrade;
    double levelValue = MinRange + (MaxRange - MinRange) * grade;
    double red, green, blue;
    int iElement;

    /* 0.0 - 1.0 のグレード値に対応する色を得る。 */
    AutoGL_GetContourColor (&red, &green, &blue,
			    grade);
    AutoGL_SetColor (red, green, blue);
  
    for (iElement = 0; iElement < NElements; iElement++) {
      int nodeId0 = Nop[iElement][0];
      int nodeId1 = Nop[iElement][1];
      int nodeId2 = Nop[iElement][2];
      double triPoints[3][3];
      double triValues[3];

      if (iElement % nTasks != taskId) continue;

      triPoints[0][0] = Coord[nodeId0][0];
      triPoints[0][1] = Coord[nodeId0][1];
      triPoints[0][2] = 0.0;
      triPoints[1][0] = Coord[nodeId1][0];
      triPoints[1][1] = Coord[nodeId1][1];
      triPoints[1][2] = 0.0;
      triPoints[2][0] = Coord[nodeId2][0];
      triPoints[2][1] = Coord[nodeId2][1];
      triPoints[2][2] = 0.0;

      triValues[0] = Scalar[nodeId0];
      triValues[1] = Scalar[nodeId1];
      triValues[2] = Scalar[nodeId2];
      
      /* この三角形内のlevelValueでの等高線を表示する。 */
      AutoGL_DrawTriangleIsoline
	(levelValue, triPoints, triValues);
    }
  }
}

/* 変形図を表示する。 */
static void PlotDeformation (void) 
{
  int iElement;
  int taskId = AutoNOSim_MyTaskId ();
  int nTasks = AutoNOSim_NTasks ();

  AutoGL_SetColor (0, 1, 1);
  for (iElement = 0; iElement < NElements; iElement++) {
    int nodeId0 = Nop[iElement][0];
    int nodeId1 = Nop[iElement][1];
    int nodeId2 = Nop[iElement][2];
    double p0X = Coord[nodeId0][0] + Displacement[nodeId0][0] * Scale;
    double p0Y = Coord[nodeId0][1] + Displacement[nodeId0][1] * Scale;
    double p1X = Coord[nodeId1][0] + Displacement[nodeId1][0] * Scale;
    double p1Y = Coord[nodeId1][1] + Displacement[nodeId1][1] * Scale;
    double p2X = Coord[nodeId2][0] + Displacement[nodeId2][0] * Scale;
    double p2Y = Coord[nodeId2][1] + Displacement[nodeId2][1] * Scale;

    if (iElement % nTasks != taskId) continue;

    AutoGL_DrawLine (p0X, p0Y, 0,
		     p1X, p1Y, 0);
    AutoGL_DrawLine (p1X, p1Y, 0,
		     p2X, p2Y, 0);
    AutoGL_DrawLine (p2X, p2Y, 0,
		     p0X, p0Y, 0);
  }
}



/* ビューの再描画のためのコールバック関数 */
/* ビューが再表示されるごとに呼ばれる。 */
static void RedrawView (void) 
{
  if (MeshIsShown) {
    PlotMesh ();
  }

  switch (PlotType) {
  case PLOT_NONE:
    break;
  case PLOT_CONTOUR:
    PlotContour ();
    break;
  case PLOT_LINE_CONTOUR:
    PlotLineContour ();
    break;
  case PLOT_DEFORMATION:
    PlotDeformation ();
    break;
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
  /* ビューを準備する。 */

  AutoGL_SetViewSize (30); 
  AutoGL_SetViewCenter (15, 10, 0); 

  AutoGL_SetViewRedrawCallback (RedrawView);   
  AutoGL_SetDefaultCallbackInMode2D (NULL);  /* 二次元アプリ用イベント処理 */
#if 0
  AutoGL_EnableDragInMode2D ();              /* ドラッグ有効 */
#endif

  /* パネルviewingを準備する。 */

#if 0
  AutoGL_AddGroup ("viewing");    
  AutoGL_SetPanelInView2D ();        /* 二次元アプリ用ビューイング変数 */
  AutoGL_SetPanelForInterpreter ();  /* コマンドライン処理 */
  AutoGL_SetPanelForSave ();         /* 画面セーブ機能 */
#endif

  /* パネルmodelを準備する。 */

  AutoGL_AddGroup ("model");    
  AutoGL_SetPanelInMode2D ();        /* 二次元アプリ用イベント処理 */

  AutoGL_AddComment ();
  AutoGL_SetLabel ("Model Parameters");
  AutoGL_AddBoolean (&MeshIsShown, "MeshIsShown");
  AutoGL_AddInteger (&PlotType, "PlotType");
  AutoGL_AddIntegerItem ("none");      
  AutoGL_AddIntegerItem ("contour"); 
  AutoGL_AddIntegerItem ("line_contour"); 
  AutoGL_AddIntegerItem ("deformation");    
  AutoGL_AddReal (&MinRange, "MinRange");
  AutoGL_AddReal (&MaxRange, "MaxRange");
  AutoGL_AddReal (&Scale, "Scale");
  AutoGL_AddCallback (UpdateModel, "UpdateModel");
}






