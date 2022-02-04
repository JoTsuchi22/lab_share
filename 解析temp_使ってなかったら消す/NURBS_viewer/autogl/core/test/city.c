/*  AutoGL ＜お手軽＞ のサンプルプログラム例 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <unistd.h>

/* AutoGL ”おてがる”を使うには,このヘッダファイルをインクルードする。*/
#include <autogl.h>



#define MAX_CAR 100
int NCars;
int Car_type[MAX_CAR];
int Car_color[MAX_CAR];
double Car_centerX[MAX_CAR];
double Car_centerY[MAX_CAR];
double Car_angle[MAX_CAR];
double Car_speed[MAX_CAR];

static void Initialize (void)
{
  NCars = 0;

  Car_type[NCars] = 0;
  Car_color[NCars] = 0;
  Car_centerX[NCars] = 0.0;
  Car_centerY[NCars] = 0.0;
  Car_angle[NCars] = 180.0;
  Car_speed[NCars] = 2.5;
  NCars++;

  Car_type[NCars] = 0;
  Car_color[NCars] = 2;
  Car_centerX[NCars] = 25.0;
  Car_centerY[NCars] = 0.0;
  Car_angle[NCars] = 0.0;
  Car_speed[NCars] = 5.0;
  NCars++;

  Car_type[NCars] = 1;
  Car_color[NCars] = 2;
  Car_centerX[NCars] = 0.0;
  Car_centerY[NCars] = 50.0;
  Car_angle[NCars] = 0.0;
  Car_speed[NCars] = 2.5;
  NCars++;

  Car_type[NCars] = 1;
  Car_color[NCars] = 1;
  Car_centerX[NCars] = 50.0;
  Car_centerY[NCars] = -50.0;
  Car_angle[NCars] = 270.0;
  Car_speed[NCars] = 1.25;
  NCars++;

  Car_type[NCars] = 1;
  Car_color[NCars] = 0;
  Car_centerX[NCars] = 50.0;
  Car_centerY[NCars] = 50.0;
  Car_angle[NCars] = 0.0;
  Car_speed[NCars] = 2.0;
  NCars++;

  Car_type[NCars] = 0;
  Car_color[NCars] = 0;
  Car_centerX[NCars] = -100.0;
  Car_centerY[NCars] = 50.0;
  Car_angle[NCars] = 0.0;
  Car_speed[NCars] = 2.5;
  NCars++;

  Car_type[NCars] = 0;
  Car_color[NCars] = 1;
  Car_centerX[NCars] = -100.0;
  Car_centerY[NCars] = -50.0;
  Car_angle[NCars] = 0.0;
  Car_speed[NCars] = 1.25;
  NCars++;

  Car_type[NCars] = 0;
  Car_color[NCars] = 2;
  Car_centerX[NCars] = -50.0;
  Car_centerY[NCars] = 0.0;
  Car_angle[NCars] = 180.0;
  Car_speed[NCars] = 2.0;
  NCars++;
}

static void AdvanceCar (int carId)
{
  double centerX = Car_centerX[carId];
  double centerY = Car_centerY[carId];
  double angle = Car_angle[carId];
  double speed = Car_speed[carId];

  centerX += speed * cos (angle / 180 * M_PI);
  centerY += speed * sin (angle / 180 * M_PI);
  Car_centerX[carId] = centerX;
  Car_centerY[carId] = centerY;

  /* turn right / left on the cross point */
  if ((int)ceil (fabs (centerX)) % 50 <= 1
      && (int)ceil (fabs (centerY)) % 50 <= 1) {
    angle += 90;
    Car_angle[carId] = angle;
  }
}

static void NextStep (void)
{
  int iCar;

  for (iCar = 0; iCar < NCars; iCar++) {
    AdvanceCar (iCar);
  }
}



static void DrawGround (void)
{
  AutoGL_SetColor (0.8, 0.7, 0.5); 
  AutoGL_DrawQuadrangle 
    (-1000, -1000, 0,
     1000, -1000, 0,
     1000, 1000, 0,
     -1000, 1000, 0);
}

static void DrawRoad 
(double x0, double y0, double x1, double y1)
{
  double halfWidth = 4.0;
  double nx, ny;
  double x0_0, y0_0, x0_1, y0_1, x1_0, y1_0, x1_1, y1_1;
  
  AutoGL_GetAnyPerpendicularDirection2D 
    (&nx, &ny,
     x1 - x0, y1 - y0);
  x0_0 = x0 - nx * halfWidth;
  y0_0 = y0 - ny * halfWidth;
  x0_1 = x0 + nx * halfWidth;
  y0_1 = y0 + ny * halfWidth;
  x1_0 = x1 - nx * halfWidth;
  y1_0 = y1 - ny * halfWidth;
  x1_1 = x1 + nx * halfWidth;
  y1_1 = y1 + ny * halfWidth;

  AutoGL_SetColor (0.8, 0.8, 0.8); 
  AutoGL_DrawQuadrangle
    (x0_0, y0_0, 0.1, 
     x1_0, y1_0, 0.1,
     x1_1, y1_1, 0.1, 
     x0_1, y0_1, 0.1);
}

static void DrawRoads (void)
{
  int iRoad;
  
  for (iRoad = 0; iRoad < 10; iRoad++) {
    DrawRoad (-1000, -200 + iRoad * 50, 
	      1000, -200 + iRoad * 50);
    DrawRoad (-200 + iRoad * 50, -1000, 
	      -200 + iRoad * 50, 1000);
  }
}

static void MakeCarMacroType0 (void)
{
  AutoGL_OpenMacro (0);
  AutoGL_FillBox3D (-2, -1, 0.5, 
		    2, 1, 1.5);
  AutoGL_FillBox3D (-1, -1, 1.5, 
		    1, 1, 2.5);
  AutoGL_CloseMacro ();
}

static void MakeCarMacroType1 (void)
{
  AutoGL_OpenMacro (1);
  AutoGL_FillBox3D (-3, -1, 0.5, 
		    3, 1, 2.0);
  AutoGL_FillBox3D (2, -1, 2.0, 
		    3, 1, 3.0);
  AutoGL_CloseMacro ();
}

static void DrawCar (int carId)
{
  int macroId = Car_type[carId];
  double centerX = Car_centerX[carId];
  double centerY = Car_centerY[carId];
  double angle = Car_angle[carId];

  switch (Car_color[carId]) {
  case 0:
    AutoGL_SetColor (1.0, 0.0, 0.0); 
    break;
  case 1:
    AutoGL_SetColor (0.0, 0.0, 1.0); 
    break;
  case 2:
    AutoGL_SetColor (1.0, 1.0, 0.0); 
    break;
  default:
    assert(0);
    break;
  }

  AutoGL_PushMatrix ();
  AutoGL_Translate (centerX, centerY, 0.0);
  AutoGL_Rotate (angle, 0, 0, 1);
  AutoGL_DrawMacro (macroId);
  AutoGL_PopMatrix ();
}

static void DrawCity (void)
{
  int iCar;

  DrawGround ();
  DrawRoads ();

  for (iCar = 0; iCar < NCars; iCar++) {
    DrawCar (iCar);
  }
}



/* コールバック関数群 */

static void RedrawView (void) 
{
  static int macroIsInitialized = 0;
  
  if (!macroIsInitialized) {
    MakeCarMacroType0 ();
    MakeCarMacroType1 ();

    macroIsInitialized = 1;
  }

  DrawCity ();
}

static void IdleEvent (void)
{
  NextStep ();
  AutoGL_DrawView ();

  usleep (100000);
}

static void Animate (void) 
{
  static int idleEventIsOn = 0;

  if (idleEventIsOn) {
    AutoGL_SetIdleEventCallback (NULL);   
    idleEventIsOn = 0;
  } else {
    AutoGL_SetIdleEventCallback (IdleEvent);   
    idleEventIsOn = 1;
  }
}



/* 関数AutoGL_SetUpはユーザー側プログラムごとに必ず一つ用意すること。*/
/* ここで、コールバック関数や制御変数などを登録する。*/
void AutoGL_SetUp (int argc, char *argv[]) 
{
  Initialize ();

  AutoGL_SetViewSize (100);        
  AutoGL_SetViewCenter (0, 0, 0);   
  AutoGL_SetViewDirection (0, 0, 1);

  AutoGL_SetPerspectiveViewFlag (1);
  AutoGL_SetViewRedrawCallback (RedrawView); 
  AutoGL_SetDefaultCallbackInMode3D (NULL);
  AutoGL_EnableIdleEvent ();   
  AutoGL_EnableDragInMode3D ();               
  AutoGL_SetPanelInMode3D ();

  AutoGL_AddCallback (Animate, "Animate");      
}






