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



#define X_DIVISION 100
#define Y_DIVISION 100
#define Z_DIVISION 100

static int GridIsShown = 0;
static int CellIsShown = 0;
static int ContourIsShown = 0;

static float Coord[X_DIVISION][Y_DIVISION + 1][Z_DIVISION + 1][3];
static float Grade[X_DIVISION][Y_DIVISION + 1][Z_DIVISION + 1];
static float Color[X_DIVISION][Y_DIVISION + 1][Z_DIVISION + 1][3];

static int StartI, EndI;

static void MakeModel (void)
{
  int i, j, k;
  int taskId = AutoNOSim_MyTaskId ();
  int nTasks = AutoNOSim_NTasks ();

  StartI = (X_DIVISION / nTasks) * taskId;
  EndI = (X_DIVISION / nTasks) * (taskId + 1);
  if (taskId == nTasks - 1) {
    EndI = X_DIVISION;
  }

  for (i = StartI; i < EndI; i++) {
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

  for (i = StartI; i < EndI; i++) {
    for (j = 0; j <= Y_DIVISION; j++) {
      for (k = 0; k <= Z_DIVISION; k++) {
	double grade = Grade[i][j][k];
	double red, green, blue;
	
	AutoGL_GetContourColor (&red, &green, &blue,
				grade);
	Color[i][j][k][0] = red;
	Color[i][j][k][1] = green;
	Color[i][j][k][2] = blue;
      }
    }
  }
}

static void DrawGrid (void)
{
  int i, j, k;

  AutoGL_SetColor2 (1, 0, 0);

  for (i = StartI; i < EndI; i++) {
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

	AutoGL_DrawLine2 (x00, y00, z00,
			  x10, y10, z10);
	AutoGL_DrawLine2 (x10, y10, z10,
			  x11, y11, z11);
	AutoGL_DrawLine2 (x11, y11, z11,
			  x01, y01, z01);
	AutoGL_DrawLine2 (x01, y01, z01,
			  x00, y00, z00);
      }
    }
  }
}

static void DrawCell (void)
{
  int i, j, k;

  AutoGL_SetColor2 (0, 1, 1);

  for (i = StartI; i < EndI; i++) {
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

	AutoGL_DrawTriangle2 (x00, y00, z00,
			      x10, y10, z10,
			      x11, y11, z11,
			      1.0f, 0.0f, 0.0f);
	AutoGL_DrawTriangle2 (x11, y11, z11,
			      x01, y01, z01,
			      x00, y00, z00,
			      1.0f, 0.0f, 0.0f);
      }
    }
  }
}

static void DrawContour (void)
{
  int i, j, k;

  for (i = StartI; i < EndI; i++) {
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

	float r00 = Color[i][j + 0][k + 0][0];
	float g00 = Color[i][j + 0][k + 0][1];
	float b00 = Color[i][j + 0][k + 0][2];
	float r10 = Color[i][j + 1][k + 0][0];
	float g10 = Color[i][j + 1][k + 0][1];
	float b10 = Color[i][j + 1][k + 0][2];
	float r11 = Color[i][j + 1][k + 1][0];
	float g11 = Color[i][j + 1][k + 1][1];
	float b11 = Color[i][j + 1][k + 1][2];
	float r01 = Color[i][j + 0][k + 1][0];
	float g01 = Color[i][j + 0][k + 1][1];
	float b01 = Color[i][j + 0][k + 1][2];

	AutoGL_DrawColorInterpolatedTriangle2 
	  (x00, y00, z00, r00, g00, b00,
	   x10, y10, z10, r10, g10, b10,
	   x11, y11, z11, r11, g11, b11,
	   1.0f, 0.0f, 0.0f);
	AutoGL_DrawColorInterpolatedTriangle2 
	  (x11, y11, z11, r11, g11, b11,
	   x01, y01, z01, r01, g01, b01,
	   x00, y00, z00, r00, g00, b00,
	   1.0f, 0.0f, 0.0f);
      }
    }
  }
}

static void RedrawView (void) 
{
  AutoGL_DrawCoordinateSystem3D (20);

  AutoGL_SetTwoSideFaceFlag (0);

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

static void UpdateModel (void) 
{
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

  AutoGL_SetPanelInMode3D ();

  AutoGL_AddBoolean (&GridIsShown, "GridIsShown"); 
  AutoGL_SetLabel ("grid");

  AutoGL_AddBoolean (&CellIsShown, "CellIsShown"); 
  AutoGL_SetLabel ("cell");

  AutoGL_AddBoolean (&ContourIsShown, "ContourIsShown"); 
  AutoGL_SetLabel ("contour");

  AutoGL_AddCallback (UpdateModel, "UpdateModel");

  AutoGL_ResizeView (800, 800);
}


