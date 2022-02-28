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

static void DrawSpheres (void)
{
  int i, j, k;
  int shpereMacroId = 0;

  AutoGL_OpenMacro (shpereMacroId);
  AutoGL_DrawSphere3D (0.0, 0.0, 0.0, 
		       1.0, 
		       1);
  AutoGL_CloseMacro ();

  AutoGL_SetColor (0, 1, 1);

  for (i = 0; i < X_DIVISION; i++) {
    for (j = 0; j < Y_DIVISION; j++) {
      for (k = 0; k < Z_DIVISION; k++) {
	double x = 100.0 / X_DIVISION * i;
	double y = 100.0 / Y_DIVISION * j;
	double z = 100.0 / Z_DIVISION * k;

	AutoGL_PushMatrix ();
	AutoGL_Translate (x, y, z);
	AutoGL_DrawMacro (shpereMacroId);
	AutoGL_PopMatrix ();
      }
    }
  }
}

static void RedrawView (void) 
{
  DrawSpheres ();
}

static void UpdateModel (void) 
{
  AutoGL_DrawView ();  
}



void AutoGL_SetUp (int argc, char *argv[]) 
{
  AutoGL_SetViewSize (70); 
  AutoGL_SetViewCenter (50, 50, 50); 
  AutoGL_SetViewDirection (3, 2, 1);

  AutoGL_SetViewRedrawCallback (RedrawView);    
  AutoGL_SetDefaultCallbackInMode3D (NULL); 
#if 0
  AutoGL_EnableDragInMode3D ();
#endif

  AutoGL_SetPanelInMode3D ();

  AutoGL_AddCallback (UpdateModel, "UpdateModel");
}


