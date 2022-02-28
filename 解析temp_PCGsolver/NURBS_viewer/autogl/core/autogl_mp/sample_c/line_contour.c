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

/* 二次元差分法解析結果の可視化 */
/* 面コンタ-を表示する。 */

/* ユーティリティの可視化（面コンタ-）の利用 */
/* 引数など詳しい情報は, autogl_utility.hを参照。*/
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

/* コンタ-の最大・最小値 */
double MinRange = 0.0;
double MaxRange = 1000.0;

/* コンタ-の分割数 */
int NGrades = 20;



/* メッシュを初期化する。 */
static void InitializeMesh ()
{
  /* 特になにもなし。 */
}

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

/* スカラー値の線コンタ-（等高線, 等値線）を描画する。 */
static void PlotLineContour (void) 
{
  int iGrade;
  int taskId = AutoNOSim_MyTaskId ();
  int nTasks = AutoNOSim_NTasks ();

  for (iGrade = 1; iGrade < NGrades; iGrade++) {
    double grade = 1.0 / NGrades * iGrade;
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



/* ビューの再描画のためのコールバック関数 */
/* ビューが再表示されるごとに呼ばれる。 */
static void RedrawView (void) 
{
  static int modelIsRendered = 0;

  /* まず、一度だけモデルをディスプレイリストへ描画しておく。 */
  if (!modelIsRendered) {

    /* ディスプレイリストを開き、これを空にする。 */
    AutoGL_OpenDisplayList ();
    
    PlotMesh ();
    PlotLineContour ();

    /* ディスプレイリストを閉じる。 */
    AutoGL_CloseDisplayList ();

    modelIsRendered = 1;
  }

  /* 画面右側に線コンタ-のバーを描画する。 */
  AutoGL_DrawLineContourMap (80, MinRange, MaxRange, NGrades);

  /* ディスプレイリストに入れた描画命令を実際に描画する。 */
  AutoGL_DrawDisplayList ();
}

/* 関数AutoGL_SetUpはユーザー側プログラムごとに必ず一つ用意すること。*/
/* ここで、コールバック関数や制御変数などを登録する。*/
void AutoGL_SetUp (int argc, char *argv[]) 
{
  InitializeMesh ();

  AutoGL_SetViewCenter (15, 10, 0);   /* ビューの注視点 */
  AutoGL_SetViewSize (40);            /* ビューのサイズ */

  /* おまじない */
  AutoGL_SetViewRedrawCallback (RedrawView);    
  AutoGL_SetMode2D (AUTOGL_MODE_2D_SCALE);     /* マウスで拡大縮小 */
  AutoGL_SetDefaultCallbackInMode2D (NULL);    /* 二次元アプリ用イベント処理 */
  AutoGL_SetPanelInMode2D ();

#if 0
  /* もし描画が重ければドラッグ中は描画させないほうがよい。 */
  AutoGL_EnableDragInMode2D ();
#endif
}


