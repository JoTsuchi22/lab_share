/*  AutoGL ＜お手軽＞ のサンプルプログラム例 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/* AutoGL ”おてがる”を使うには,このヘッダファイルをインクルードする。*/
#include <autogl.h>



#define X_DIVISION 100
#define Y_DIVISION 100
#define Z_DIVISION 100

static int GridIsShown = 0;
static int CellIsShown = 0;
static int ContourIsShown = 0;

static float Coord[X_DIVISION][Y_DIVISION + 1][Z_DIVISION + 1][3];
static float Grade[X_DIVISION][Y_DIVISION + 1][Z_DIVISION + 1];

static void MakeModel (void)
{
  int i, j, k;

  for (i = 0; i < X_DIVISION; i++) {
    for (j = 0; j <= Y_DIVISION; j++) {
      for (k = 0; k <= Z_DIVISION; k++) {
	Coord[i][j][k][0] = 100.0 / X_DIVISION * i;
	Coord[i][j][k][1] = 100.0 / Y_DIVISION * j;
	Coord[i][j][k][2] = 100.0 / Z_DIVISION * k;

	Grade[i][j][k] = 1.0 * (i + j + k) 
	  / (X_DIVISION + Y_DIVISION + Z_DIVISION);
      }
    }
  }
}

static void DrawGrid (void)
{
  int i, j, k;

  AutoGL_SetColor (1, 0, 0);

  for (i = 0; i < X_DIVISION; i++) {
    for (j = 0; j < Y_DIVISION; j++) {
      for (k = 0; k < Z_DIVISION; k++) {
	float x00 = Coord[i][j + 0][k + 0][0];
	float y00 = Coord[i][j + 0][k + 0][1];
	float z00 = Coord[i][j + 0][k + 0][2];
	float x10 = Coord[i][j + 1][k + 0][0];
	float y10 = Coord[i][j + 1][k + 0][1];
	float z10 = Coord[i][j + 1][k + 0][2];
	float x11 = Coord[i][j + 1][k + 1][0];
	float y11 = Coord[i][j + 1][k + 1][1];
	float z11 = Coord[i][j + 1][k + 1][2];
	float x01 = Coord[i][j + 0][k + 1][0];
	float y01 = Coord[i][j + 0][k + 1][1];
	float z01 = Coord[i][j + 0][k + 1][2];

	AutoGL_DrawLine (x00, y00, z00,
			 x10, y10, z10);
	AutoGL_DrawLine (x10, y10, z10,
			 x11, y11, z11);
	AutoGL_DrawLine (x11, y11, z11,
			 x01, y01, z01);
	AutoGL_DrawLine (x01, y01, z01,
			 x00, y00, z00);
      }
    }
  }
}

static void DrawCell (void)
{
  int i, j, k;

  AutoGL_SetColor (0, 1, 1);

  for (i = 0; i < X_DIVISION; i++) {
    for (j = 0; j < Y_DIVISION; j++) {
      for (k = 0; k < Z_DIVISION; k++) {
	float x00 = Coord[i][j + 0][k + 0][0];
	float y00 = Coord[i][j + 0][k + 0][1];
	float z00 = Coord[i][j + 0][k + 0][2];
	float x10 = Coord[i][j + 1][k + 0][0];
	float y10 = Coord[i][j + 1][k + 0][1];
	float z10 = Coord[i][j + 1][k + 0][2];
	float x11 = Coord[i][j + 1][k + 1][0];
	float y11 = Coord[i][j + 1][k + 1][1];
	float z11 = Coord[i][j + 1][k + 1][2];
	float x01 = Coord[i][j + 0][k + 1][0];
	float y01 = Coord[i][j + 0][k + 1][1];
	float z01 = Coord[i][j + 0][k + 1][2];

	AutoGL_DrawTriangle (x00, y00, z00,
			     x10, y10, z10,
			     x11, y11, z11);
	AutoGL_DrawTriangle (x11, y11, z11,
			     x01, y01, z01,
			     x00, y00, z00);
      }
    }
  }
}

static void DrawContour (void)
{
  int i, j, k;

  for (i = 0; i < X_DIVISION; i++) {
    for (j = 0; j < Y_DIVISION; j++) {
      for (k = 0; k < Z_DIVISION; k++) {
	float x00 = Coord[i][j + 0][k + 0][0];
	float y00 = Coord[i][j + 0][k + 0][1];
	float z00 = Coord[i][j + 0][k + 0][2];
	float x10 = Coord[i][j + 1][k + 0][0];
	float y10 = Coord[i][j + 1][k + 0][1];
	float z10 = Coord[i][j + 1][k + 0][2];
	float x11 = Coord[i][j + 1][k + 1][0];
	float y11 = Coord[i][j + 1][k + 1][1];
	float z11 = Coord[i][j + 1][k + 1][2];
	float x01 = Coord[i][j + 0][k + 1][0];
	float y01 = Coord[i][j + 0][k + 1][1];
	float z01 = Coord[i][j + 0][k + 1][2];

	float g00 = Grade[i][j + 0][k + 0];
	float g10 = Grade[i][j + 1][k + 0];
	float g11 = Grade[i][j + 1][k + 1];
	float g01 = Grade[i][j + 0][k + 1];

	AutoGL_DrawContourTriangle (x00, y00, z00, g00,
				    x10, y10, z10, g10,
				    x11, y11, z11, g11);
	AutoGL_DrawContourTriangle (x11, y11, z11, g11,
				    x01, y01, z01, g01,
				    x00, y00, z00, g00);
      }
    }
  }
}

static void RedrawView (void) 
{
  AutoGL_DrawCoordinateSystem3D (20);

  if (GridIsShown) {
    DrawGrid ();
  }
  if (CellIsShown) {
    DrawCell ();
  }
  if (ContourIsShown) {
    DrawContour ();
  }
}

static void RenderView (void)
{
  AutoGL_SetImageFileName ("batch000.ppm");
  
  GridIsShown = 1;
  CellIsShown = 1;
  ContourIsShown = 0;

  AutoGL_SetViewSize (80); 
  AutoGL_SetViewDirection (3, 2, 1);

  AutoGL_DrawView ();  

  AutoGL_SetImageFileName ("batch001.ppm");
  
  GridIsShown = 0;
  CellIsShown = 0;
  ContourIsShown = 1;

  AutoGL_SetViewSize (100); 
  AutoGL_SetViewDirection (1, 2, 3);

  AutoGL_DrawView ();  
}



void AutoGL_SetUp (int argc, char *argv[]) 
{
  MakeModel ();

  AutoGL_SetViewSize (100); 
  AutoGL_SetViewCenter (50, 50, 50); 
  AutoGL_SetViewDirection (3, 2, 1);

  AutoGL_SetViewRedrawCallback (RedrawView);    
  AutoGL_SetDefaultCallbackInMode3D (NULL); 

  AutoGL_ResizeView (800, 800);

  AutoGL_SetBatchProcessCallback (RenderView);
}


