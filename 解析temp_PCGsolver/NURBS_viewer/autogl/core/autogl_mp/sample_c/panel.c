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

/* GUIの利用 */
/* 引数など詳しい情報は, autogl_uims.hを参照。*/
/* なお、これはautogl.hからincludeされている。*/



/* ビュー上に表示されるモデルを表現するための変数群 */

/* 原点のラベル : 文字列 */
static char OriginLabel[256];

/* 座標軸を表示するかどうか : 0または1の整数（論理値） */
static int AxisIsShown;

/* モデルの形のタイプ : 0 - 2 の整数 (SHAPE_) */
#define SHAPE_LINE 0
#define SHAPE_TRIANGLE 1
#define SHAPE_BOX 2
static int Shape;

/* 整数変数 */
static int Number;

/* 実数変数 */
static double Value;



/* コールバック関数群 */

/* 印刷コマンドのためのコールバック関数 */
/* "Print"ボタンを選ぶと呼ばれる。 */
static void Print (void) 
{
  printf (" PRINT \n");

  /* 現在の原点ラベルを表示する。 */
  printf (" OriginLabel is <<%s>> \n", OriginLabel);
  
  /* 現在座標軸が表示されているかどうかを表示する。 */
  printf (" AxisIsShown is %d \n", AxisIsShown);

  /* 現在のモデルの形のタイプを表示する。 */
  printf (" Shape is ");
  switch (Shape) {
  case SHAPE_LINE:
    printf ("SHAPE_LINE");
    break;
  case SHAPE_TRIANGLE:
    printf ("SHAPE_TRIANGLE");
    break;
  case SHAPE_BOX:
    printf ("SHAPE_BOX");
    break;
  }
  printf ("\n");

  /* 整数変数Numberの値を表示する。 */
  printf (" Number is %d \n", Number);

  /* 実数変数Valueの値を表示する。 */
  printf (" Value is %f \n", Value);
}

/* ビューの再描画のためのコールバック関数 */
/* ビューが再表示されるごとに呼ばれる。 */
static void RedrawView (void) 
{
  int taskId = AutoNOSim_MyTaskId ();

  if (taskId == 0) {

    /* もし座標軸を表示するならば、 */
    if (AxisIsShown) {
      
      /* X軸を描画する。 */
      AutoGL_SetColor (1, 0, 0); 
      AutoGL_DrawLine (0, 0, 0,
		       50, 0, 0);
      
      /* Y軸を描画する。*/
      AutoGL_SetColor (0, 1, 0); 
      AutoGL_DrawLine (0, 0, 0,
		       0, 50, 0);
      
      /* Z軸を描画する。*/
      AutoGL_SetColor (0, 0, 1); 
      AutoGL_DrawLine (0, 0, 0,
		       0, 0, 50);
    }
    
    /* 原点ラベルを表示する。 */
    AutoGL_SetColor (1, 1, 1);
    AutoGL_DrawString (0, 0, 0, OriginLabel); 
  }
  
  /* 形のタイプに応じて, */
  switch (Shape) {
  case SHAPE_LINE:
    if (taskId == 1) {

      /* 直方体の１２本の稜線を線分として描画する。 */
      AutoGL_SetColor (1, 1, 1);
      AutoGL_DrawLine (0, 0, 0,
		       10, 0, 0);
      AutoGL_DrawLine (10, 0, 0,
		       10, 20, 0);
      AutoGL_DrawLine (10, 20, 0,
		       0, 20, 0);
      AutoGL_DrawLine (0, 20, 0,
		       0, 0, 0);
      AutoGL_DrawLine (0, 0, 30,
		       10, 0, 30);
      AutoGL_DrawLine (10, 0, 30,
		       10, 20, 30);
      AutoGL_DrawLine (10, 20, 30,
		       0, 20, 30);
      AutoGL_DrawLine (0, 20, 30,
		       0, 0, 30);
      AutoGL_DrawLine (0, 0, 0,
		       0, 0, 30);
      AutoGL_DrawLine (0, 20, 0,
		       0, 20, 30);
      AutoGL_DrawLine (10, 20, 0,
		       10, 20, 30);
      AutoGL_DrawLine (10, 0, 0,
		       10, 0, 30);
    }
    break;

  case SHAPE_TRIANGLE:
    if (taskId == 2) {

      /* 原点を通る３枚の三角形を描画する。 */
      AutoGL_SetColor (1, 1, 1);
      AutoGL_DrawTriangle (0, 0, 0,
			   0, 20, 0,
			   0, 0, 30);
      AutoGL_DrawTriangle (10, 0, 0,
			   0, 0, 0, 
			   0, 0, 30);
      AutoGL_DrawTriangle (10, 0, 0,
			   0, 20, 0, 
			   0, 0, 0);
    }
    break;

  case SHAPE_BOX:
    if (taskId == 3) {

      /* 直方体(0, 0, 0) -（10, 20, 30)の6枚の面（四辺形）を描画する。 */
      AutoGL_SetColor (0, 1, 1);
      AutoGL_DrawQuadrangle (0, 0, 0,
			     10, 0, 0,
			     10, 20, 0,
			     0, 20, 0);
      AutoGL_SetColor (1, 1, 0);
      AutoGL_DrawQuadrangle (0, 0, 30,
			     10, 0, 30,
			     10, 20, 30,
			     0, 20, 30);
      AutoGL_SetColor (1, 0, 1);
      AutoGL_DrawQuadrangle (0, 0, 0,
			     0, 20, 0,
			     0, 20, 30,
			     0, 0, 30);
      AutoGL_SetColor (0, 0, 1);
      AutoGL_DrawQuadrangle (10, 0, 0,
			     10, 20, 0,
			     10, 20, 30,
			     10, 0, 30);
      AutoGL_SetColor (0, 1, 0);
      AutoGL_DrawQuadrangle (0, 0, 0,
			     10, 0, 0,
			     10, 0, 30,
			     0, 0, 30);    
      AutoGL_SetColor (1, 0, 0);
      AutoGL_DrawQuadrangle (0, 20, 0,
			     10, 20, 0,
			     10, 20, 30,
			     0, 20, 30);
    }
    break;
  }

  if (taskId == 0) {

    /* 位置(10, 20, 30)に整数変数Numberの値を表示する。 */
    char buf[256];

    sprintf (buf, "Number:%d", Number);
    AutoGL_SetColor (1, 1, 1);
    AutoGL_DrawString (10, 20, 30,
		       buf);

  }

  if (taskId == 0) {

    /* 画面の左下隅に実数変数Valueの値を表示する。 */
    int widthDc, heightDc;
    int dx, dy;
    double x, y, z;
    char buf[256];

    /* 画面の幅と高さを得る。 */
    AutoGL_GetViewRangeDc (&widthDc, &heightDc);

    /* 画面左下から50, 50の点のデバイス座標を計算する。 */
    dx = -widthDc / 2 + 50;
    dy = -heightDc / 2 + 50;

    /* デバイス座標を世界座標（モデル座標）に変換する。 */
    AutoGL_GetPositionOfDc (&x, &y, &z,
			    dx, dy, 0);

    sprintf (buf, "Value = %f", Value);
    AutoGL_SetColor (0, 0, 1);
    AutoGL_DrawString (x, y, z,
		       buf);
  }
}

/* 描画アップデートコマンドのためのコールバック関数 */
/* "Draw"ボタンを選ぶと呼ばれる。 */
static void Draw (void) 
{
  /* ビューを再描画する。 */
  AutoGL_DrawView ();
  /* 登録されたビューの再描画関数（この場合はRedrawView ()）が呼ばれる。 */
}



/* 関数AutoGL_SetUpはユーザー側プログラムごとに必ず一つ用意すること。*/
/* この関数はプログラム起動直後にAutoGLライブラリから呼び出される。*/
/* このとき、argc, argvはmain関数から直接渡される。*/
/* ここで、コールバック関数や制御変数などを登録する。*/
void AutoGL_SetUp (int argc, char *argv[]) 
{
  /* ビューおよびパネル１つの例 */



  /* モデル用変数を初期化する。 */
  strcpy (OriginLabel, "Hello, world !");
  AxisIsShown = 1;
  Shape = SHAPE_BOX;
  Number = 20;
  Value = 0.5;

  

  /* ビューを準備する。 */

  /* ビューイングパラメータ（サイズ,視線方向など）を設定する。 */
  AutoGL_SetViewSize (70);            /* ビューの大きさ */
  AutoGL_SetViewCenter (10, 10, 10);       /* 注視点 */
  AutoGL_SetViewDirection (3, 1, 2);       /* 視線方向 */
  
  AutoGL_SetViewRedrawCallback (RedrawView);
  /* 再描画コールバック関数として上にあるRedrawView ()を指定。 */

  /* おまじない */
  AutoGL_SetDefaultCallbackInMode3D (NULL);
#if 0
  AutoGL_EnableDragInMode3D ();              /* ドラッグ有効 */
#endif



  /* パネルを準備する。 */
  /* 制御変数とコールバック関数を登録する。 */

  /* おまじない */
  AutoGL_SetPanelInMode3D ();                /* 移動回転拡大縮小など */
#if 0
  AutoGL_SetPanelForInterpreter ();          /* コマンドライン処理 */
#endif

  /* コメントの登録 : ラベル */
  AutoGL_AddComment ();
  AutoGL_SetLabel ("Model Parameters");

  /* 文字列変数の登録 : テキスト入力 */
  /* フィールド OriginLabel */
  AutoGL_AddString (OriginLabel, "OriginLabel", 80);
  /* char配列、制御変数名、配列の長さ */

  /* 整数変数の登録 : テキスト入力 */
  /* フィールド Number */
  AutoGL_AddInteger (&Number, "Number");
  /* int型変数のポインタと制御変数名 */

  /* 実数変数の登録 : テキスト入力 */
  /* フィールド Value */
  AutoGL_AddReal (&Value, "Value");
  /* double型変数のポインタと制御変数名 */

  /* 論理変数の登録 : トグル入力 */
  /* トグル show axis */
  AutoGL_AddBoolean (&AxisIsShown, "AxisIsShown"); 
  /* int型変数のポインタと制御変数名 */
  /* AxisIsShownは0または1の値を取る。 */
  AutoGL_SetLabel ("show axis"); 
  /* ラベルを指定してもよい。 */
  /* デフォルトでは制御変数名がラベルとして使用される。 */

  /* 整数変数の登録 : ラジオボタン入力 */
  /* チョイス Shape */
  AutoGL_AddInteger (&Shape, "Shape");  /* int型変数 */
  /* int型変数のポインタと制御変数名 */
  /* Shapeは0から2までの値を取る。 */
  AutoGL_AddIntegerItem ("line");       /* 値0のアイテム */
  AutoGL_AddIntegerItem ("triangle");   /* 値1のアイテム */
  AutoGL_AddIntegerItem ("box");        /* 値2のアイテム */

  /* コールバック関数の登録 : ボタン */
  /* 描画アップデートボタン Draw */
  AutoGL_AddCallback (Draw, "Draw");      
  /* 関数へのポインタとコールバック関数名 */
  /* コールバック関数として上にあるDraw ()を指定。 */

  /* コールバック関数の登録 : ボタン */
  /* 印刷ボタン Print */
  AutoGL_AddCallback (Print, "Print");  
  /* 関数へのポインタとコールバック関数名 */
  /* コールバック関数として上にあるPrint ()を指定。 */
}






