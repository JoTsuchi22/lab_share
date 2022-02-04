/*  AutoGL ~~~~~ ~~~~~~~~~~~ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/* AutoGL ~~~~~~~~~~~,~~~~~~~~~~~~~~~~~~~*/
#include <autogl.h>

/* ~~~~~~~~~~~~~~ */
/* ~~~~~~-~~~~~~ */

/* ~~~~~~~~~~~~~~~~~~-~~~~ */
/* ~~~~~~~~~~, autogl_utility.h~~~~*/
/* ~~~~~~autogl.h~~include~~~~~~*/




/* ~~~~~~~~~~~~~~~~~~~~~~~~ */

/* ~~~~~~~~~~~~~~~~~ */

/* ~~~~~~ */
#if 0
#define CELLS 20
#else
#define CELLS 100
#endif

/* ~~~~~~~~~ */
static double GridValues[CELLS + 1][CELLS + 1];

/* ~~~-~~~~~~~ */
static double MinRange = -1.0;
static double MaxRange = 5.0;
static int NGrades = 20;



/* ~~~~~~~~~~~ */
static double Function (double x, double y)
{
  return sin (x * M_PI) + y * y;
}

/* ~~~~~~~~~~~~ */
static void InitializeGrid ()
{
  int i, j;

  /* ~~~~~~~~~~~~~ */
  for (i = 0; i <= CELLS; i++) {
    for (j = 0; j <= CELLS; j++) {
      double scale = 4.0 / CELLS;
      
      GridValues[i][j] 
	= Function ((i - CELLS / 2) * scale, 
		    (j - CELLS / 2) * scale);
    }
  }
}

/* ~~~-~~~~~~ */
static void PlotContour (void) 
{
  double cellSize = 60.0 / CELLS;
  int i, j;
  int iGrade;

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

      /* ~~~~~~-~~~~~~~~~~~~ */
      AutoGL_DrawBandContourQuadrangle 
	(MinRange, MaxRange, NGrades,
	 x00, y00, 0.0, value00,
	 x10, y10, 0.0, value10,
	 x11, y11, 0.0, value11,
	 x01, y01, 0.0, value01);

#if 1
      /* ~~~~~~~~~~~~~~~~~~~-~~~~~~~~ */
      AutoGL_SetColor (0, 0, 0);
      for (iGrade = 0; iGrade < NGrades; iGrade++) {
	double value = MinRange 
	  + (MaxRange - MinRange) / NGrades * iGrade;

	AutoGL_DrawLineContourQuadrangle 
	  (value, value, 0,
	   x00, y00, 0.0, value00,
	   x10, y10, 0.0, value10,
	   x11, y11, 0.0, value11,
	 x01, y01, 0.0, value01);
      }
#endif
    }
  }
}



/* ~~~~~~~~~~~~~~~~~~~ */
/* ~~~~~~~~~~~~~~~~~~ */
static void RedrawView (void) 
{
  static int modelIsRendered = 0;

  /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  if (!modelIsRendered) {

    /* ~~~~~~~~~~~~~~~~~~~~~ */
    AutoGL_OpenDisplayList ();
    
    PlotContour ();

    /* ~~~~~~~~~~~~~~ */
    AutoGL_CloseDisplayList ();

    modelIsRendered = 1;
  }

  /* ~~~~~~~~-~~~~~~~~~ */
  AutoGL_DrawBandContourMap (80, MinRange, MaxRange, NGrades);

  /* ~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  AutoGL_DrawDisplayList ();
}

/* ~~AutoGL_SetUp~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* ~~~~~~~~~~~~~~~~~~~~~~~~~*/
void AutoGL_SetUp (int argc, char *argv[]) 
{
  InitializeGrid ();

  AutoGL_SetViewSize (60);   /* ~~~~~~~ */

  /* ~~~~~ */
  AutoGL_SetViewRedrawCallback (RedrawView);    
  AutoGL_SetMode2D (AUTOGL_MODE_2D_SCALE);     /* ~~~~~~~~ */
  AutoGL_SetDefaultCallbackInMode2D (NULL);    /* ~~~~~~~~~~~~~ */
  AutoGL_SetPanelInMode2D ();

#if 0
  /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  AutoGL_EnableDragInMode2D ();
#endif

#if 0
  AutoGL_SetGrayScaleContourMap ();
#endif
#if 0
  AutoGL_SetContourMap_BCGYRM ();
#endif
#if 0
  AutoGL_SetContourMap_MBCGYR ();
#endif
#if 0
  AutoGL_SetContourMap_BCGYR ();
#endif
}


