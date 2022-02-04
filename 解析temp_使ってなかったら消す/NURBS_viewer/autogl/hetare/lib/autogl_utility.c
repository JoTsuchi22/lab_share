/* autogl_utility.c : AutoGL utility */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "autogl_math.h"
#include "autogl_uims.h"
#include "autogl_graphics.h"
#include "autogl_gui.h"
#include "autogl_utility.h"



/* Drawing Primitive */

#define BLUE_RATIO    (0.0)
#define CYAN_RATIO    (1.0 / 21 * 5)
#define GREEN_RATIO   (1.0 / 21 * 8)
#define YELLOW_RATIO  (1.0 / 21 * 11)
#define RED_RATIO     (1.0 / 21 * 16)
#define MAGENTA_RATIO (1.0)
 
void AutoGL_GetContourColor 
(double *red_OUT, double *green_OUT, double *blue_OUT,
 double grade)
{
  if (grade < BLUE_RATIO) {
    *red_OUT = 0.0;  
    *green_OUT = 0.0;  
    *blue_OUT = 1.0;
  } else if (grade < CYAN_RATIO) {
    double v = (grade - BLUE_RATIO) / (CYAN_RATIO - BLUE_RATIO);
    *red_OUT = 0.0;  
    *green_OUT = v;  
    *blue_OUT = 1.0;
  } else if (grade < GREEN_RATIO) {
    double v = 1.0 - (grade - CYAN_RATIO) / (GREEN_RATIO - CYAN_RATIO);
    *red_OUT = 0.0;  
    *green_OUT = 1.0;  
    *blue_OUT = v;
  } else if (grade < YELLOW_RATIO) {
    double v = (grade - GREEN_RATIO) / (YELLOW_RATIO - GREEN_RATIO);
    *red_OUT = v;  
    *green_OUT = 1.0;  
    *blue_OUT = 0.0;
  } else if (grade < RED_RATIO) {
    double v = 1.0 - (grade - YELLOW_RATIO) / (RED_RATIO - YELLOW_RATIO);
    *red_OUT = 1.0;  
    *green_OUT = v;  
    *blue_OUT = 0.0;
  } else if (grade < MAGENTA_RATIO) {
    double v = (grade - RED_RATIO) / (MAGENTA_RATIO - RED_RATIO);
    *red_OUT = 1.0;  
    *green_OUT = 0.0;  
    *blue_OUT = v;
  } else {
    *red_OUT = 1.0;  
    *green_OUT = 0.0;  
    *blue_OUT = 1.0;
  }
}

void AutoGL_DrawContourTriangle 
(double x0, double y0, double z0, double grade0, 
 double x1, double y1, double z1, double grade1, 
 double x2, double y2, double z2, double grade2)
{
  double red0, green0, blue0;
  double red1, green1, blue1;
  double red2, green2, blue2;

  AutoGL_GetContourColor (&red0, &green0, &blue0,
			  grade0);
  AutoGL_GetContourColor (&red1, &green1, &blue1,
			  grade1);
  AutoGL_GetContourColor (&red2, &green2, &blue2,
			  grade2);
  AutoGL_DrawColorInterpolatedTriangle 
    (x0, y0, z0, red0, green0, blue0,
     x1, y1, z1, red1, green1, blue1,
     x2, y2, z2, red2, green2, blue2);
}



/* 3-D Viewing Parameter Handling */

static double Size = 100.0;

static double CenterX = 0.0;
static double CenterY = 0.0;
static double CenterZ = 0.0;

static double DirectionX = 0.0;
static double DirectionY = 0.0;
static double DirectionZ = 1.0;

static double UpVectorX = 0.0;
static double UpVectorY = 1.0;
static double UpVectorZ = 0.0;

static void UpdateView (void) 
{
  AutoGL_SetViewSize (Size);
  AutoGL_SetViewCenter (CenterX, CenterY, CenterZ);       
  AutoGL_SetViewDirection (DirectionX, DirectionY, DirectionZ);
  AutoGL_SetViewUpVector (UpVectorX, UpVectorY, UpVectorZ);
  AutoGL_DrawView ();  
}

static void GetViewParameter (void) 
{
  Size = AutoGL_GetViewSize ();
  AutoGL_GetViewCenter (&CenterX, &CenterY, &CenterZ);       
  AutoGL_GetViewDirection (&DirectionX, &DirectionY, &DirectionZ);
  AutoGL_GetViewUpVector (&UpVectorX, &UpVectorY, &UpVectorZ);

  /* NEED TO REDEFINE */
  /* update panel */
}

void AutoGL_SetPanelInView3D (void)
{
  GetViewParameter ();

  AutoGL_AddComment ();
  AutoGL_SetLabel ("VIEW");

  AutoGL_AddReal (&Size, "AutoGL_Size");               
  AutoGL_SetLabel ("size");               
  AutoGL_AddReal (&CenterX, "AutoGL_CenterX");        
  AutoGL_SetLabel ("center x");               
  AutoGL_AddReal (&CenterY, "AutoGL_CenterY");      
  AutoGL_SetLabel ("center y");               
  AutoGL_AddReal (&CenterZ, "AutoGL_CenterZ"); 
  AutoGL_SetLabel ("center z");               
  AutoGL_AddReal (&DirectionX, "AutoGL_DirectionX");  
  AutoGL_SetLabel ("dir x");               
  AutoGL_AddReal (&DirectionY, "AutoGL_DirectionY");
  AutoGL_SetLabel ("dir y");               
  AutoGL_AddReal (&DirectionZ, "AutoGL_DirectionZ");
  AutoGL_SetLabel ("dir z");               
  AutoGL_AddReal (&UpVectorX, "AutoGL_UpVectorX");   
  AutoGL_SetLabel ("up x");               
  AutoGL_AddReal (&UpVectorY, "AutoGL_UpVectorY");
  AutoGL_SetLabel ("up y");               
  AutoGL_AddReal (&UpVectorZ, "AutoGL_UpVectorZ");
  AutoGL_SetLabel ("up z");               

  AutoGL_AddCallback (UpdateView, "AutoGL_UpdateView");
  AutoGL_SetLabel ("update view");               
  AutoGL_AddCallback (GetViewParameter, "AutoGL_GetViewParameter");
  AutoGL_SetLabel ("get params");               
}

void AutoGL_RotateView 
(double offsetX, double offsetY, double offsetZ)
{
  double upX, upY, upZ;
  double xAxisX, xAxisY, xAxisZ;
  double yAxisX, yAxisY, yAxisZ;
  double zAxisX, zAxisY, zAxisZ;

  AutoGL_GetViewUpVector (&upX, &upY, &upZ);
  AutoGL_GetViewDirection (&zAxisX, &zAxisY, &zAxisZ);

  AutoGL_GetVectorProduct3D (&xAxisX, &xAxisY, &xAxisZ,
			     upX, upY, upZ,
			     zAxisX, zAxisY, zAxisZ);
  AutoGL_NormalizeVector3D (&xAxisX, &xAxisY, &xAxisZ,
			    xAxisX, xAxisY, xAxisZ);
  
  AutoGL_GetVectorProduct3D (&yAxisX, &yAxisY, &yAxisZ,
			     zAxisX, zAxisY, zAxisZ, 
			     xAxisX, xAxisY, xAxisZ);
  AutoGL_NormalizeVector3D (&yAxisX, &yAxisY, &yAxisZ,
			    yAxisX, yAxisY, yAxisZ);

  {
    double x, y, z;
    AutoGL_SubtractVector3D (&x, &y, &z,
			     zAxisX, zAxisY, zAxisZ,
			     xAxisX * offsetX + yAxisX * offsetY,
			     xAxisY * offsetX + yAxisY * offsetY,
			     xAxisZ * offsetX + yAxisZ * offsetY);
    AutoGL_SetViewDirection (x, y, z);
  }
  
  {
    double x, y, z;
    AutoGL_AddVector3D (&x, &y, &z,
			yAxisX, yAxisY, yAxisZ,
			zAxisX * offsetY + xAxisX * offsetZ, 
			zAxisY * offsetY + xAxisY * offsetZ, 
			zAxisZ * offsetY + xAxisZ * offsetZ);
    AutoGL_SetViewUpVector (x, y, z);
  }

}



/* 3-D Event Handler */

static int Mode = AUTOGL_MODE_3D_OTHER;    /*  AUTOGL_MODE_3D_  */
static AutoGL_CallbackType OtherCallback;
static int DragIsOn;

void AutoGL_EnableDragInMode3D (void)
{
  DragIsOn = 1;
}

static void Translate (void)
{
  static double oldX, oldY, oldZ;
  static double oldCenterX, oldCenterY, oldCenterZ;
  int event = AutoGL_GetViewEvent ();

  switch (event) {
  case AUTOGL_EVENT_POINTING_DEVICE_PRESS:
    {
      int dcX, dcY, dcZ;

      AutoGL_GetPointingDevicePositionDc (&dcX, &dcY, &dcZ);
      AutoGL_GetPositionOfDc (&oldX, &oldY, &oldZ,
			      dcX, dcY, dcZ);
      AutoGL_GetViewCenter (&oldCenterX, &oldCenterY, &oldCenterZ);
    }
    break;
  case AUTOGL_EVENT_POINTING_DEVICE_DRAG:
  case AUTOGL_EVENT_POINTING_DEVICE_RELEASE:
    if ((DragIsOn && event == AUTOGL_EVENT_POINTING_DEVICE_DRAG)
	|| (!DragIsOn && event == AUTOGL_EVENT_POINTING_DEVICE_RELEASE)) {
      int dcX, dcY, dcZ;
      double x, y, z;
      double offsetX, offsetY, offsetZ;
      double centerX, centerY, centerZ;

      AutoGL_GetPointingDevicePositionDc (&dcX, &dcY, &dcZ);
      AutoGL_SetViewCenter (oldCenterX, oldCenterY, oldCenterZ);
      AutoGL_GetPositionOfDc (&x, &y, &z,
			      dcX, dcY, dcZ);
      AutoGL_SubtractVector3D (&offsetX, &offsetY, &offsetZ,
			       x, y, z,
			       oldX, oldY, oldZ);
      AutoGL_SubtractVector3D (&centerX, &centerY, &centerZ,
			       oldCenterX, oldCenterY, oldCenterZ,
			       offsetX, offsetY, offsetZ);
      AutoGL_SetViewCenter (centerX, centerY, centerZ);
      AutoGL_DrawView ();
    }
    break;
  default:
    break;
  }
}

static void RotateXY (void)
{
  static int oldDcX, oldDcY, oldDcZ;
  int event = AutoGL_GetViewEvent ();

  switch (event) {
  case AUTOGL_EVENT_POINTING_DEVICE_PRESS:
    AutoGL_GetPointingDevicePositionDc (&oldDcX, &oldDcY, &oldDcZ);
    break;
  case AUTOGL_EVENT_POINTING_DEVICE_DRAG:
  case AUTOGL_EVENT_POINTING_DEVICE_RELEASE:
    if ((DragIsOn && event == AUTOGL_EVENT_POINTING_DEVICE_DRAG)
	|| (!DragIsOn && event == AUTOGL_EVENT_POINTING_DEVICE_RELEASE)) {
      int dcX, dcY, dcZ;
      double offsetX, offsetY;

      AutoGL_GetPointingDevicePositionDc (&dcX, &dcY, &dcZ);
      offsetX = (dcX - oldDcX) / 200.0;
      offsetY = (dcY - oldDcY) / 200.0;
      AutoGL_RotateView (offsetX, offsetY, 0.0);
      AutoGL_DrawView ();
      oldDcX = dcX;
      oldDcY = dcY;
      oldDcZ = dcZ;
    }
    break;
  default:
    break;
  }
}

static void RotateZ (void)
{
  static int oldDcX, oldDcY, oldDcZ;
  int event = AutoGL_GetViewEvent ();

  switch (event) {
  case AUTOGL_EVENT_POINTING_DEVICE_PRESS:
    {
      AutoGL_GetPointingDevicePositionDc (&oldDcX, &oldDcY, &oldDcZ);
    }
    break;
  case AUTOGL_EVENT_POINTING_DEVICE_DRAG:
  case AUTOGL_EVENT_POINTING_DEVICE_RELEASE:
    if ((DragIsOn && event == AUTOGL_EVENT_POINTING_DEVICE_DRAG)
	|| (!DragIsOn && event == AUTOGL_EVENT_POINTING_DEVICE_RELEASE)) {
      int dcX, dcY, dcZ;
      double offsetX, offsetY, offsetZ;
      double offset;

      AutoGL_GetPointingDevicePositionDc (&dcX, &dcY, &dcZ);
      offsetX = (dcX - oldDcX) / 200.0;
      offsetY = (dcY - oldDcY) / 200.0;
      offsetZ = (dcZ - oldDcZ) / 200.0;
      offset = offsetX;
      if (fabs (offset) < fabs (offsetY)) offset = offsetY;
      if (fabs (offset) < fabs (offsetZ)) offset = offsetZ;
      AutoGL_RotateView (0, 0, offset);
      AutoGL_DrawView ();
      oldDcX = dcX;
      oldDcY = dcY;
      oldDcZ = dcZ;
    }
    break;
  default:
    break;
  }
}

static void Scale (void)
{
  static int oldDcX, oldDcY, oldDcZ;
  static double oldViewSize;
  int event = AutoGL_GetViewEvent ();

  switch (event) {
  case AUTOGL_EVENT_POINTING_DEVICE_PRESS:
    {
      AutoGL_GetPointingDevicePositionDc (&oldDcX, &oldDcY, &oldDcZ);
      oldViewSize = AutoGL_GetViewSize ();
    }
    break;
  case AUTOGL_EVENT_POINTING_DEVICE_DRAG:
  case AUTOGL_EVENT_POINTING_DEVICE_RELEASE:
    if ((DragIsOn && event == AUTOGL_EVENT_POINTING_DEVICE_DRAG)
	|| (!DragIsOn && event == AUTOGL_EVENT_POINTING_DEVICE_RELEASE)) {
      int dcX, dcY, dcZ;
      double scaleX, scaleY, scaleZ;
      double scale;

      AutoGL_GetPointingDevicePositionDc (&dcX, &dcY, &dcZ);
      scaleX = (dcX - oldDcX) / 500.0;
      scaleY = (dcY - oldDcY) / 500.0;
      scaleZ = (dcZ - oldDcZ) / 500.0;
      scale = scaleX;
      if (fabs (scale) < fabs (scaleY)) scale = scaleY;
      if (fabs (scale) < fabs (scaleZ)) scale = scaleZ;
      scale += 1.0;
      if (scale < 0.1) scale = 0.1;
      AutoGL_SetViewSize (oldViewSize / scale);
      AutoGL_DrawView ();
    }
    break;
  default:
    break;
  }
}

static void DefaultViewEventCallback (void) 
{
  switch (Mode) {
  case AUTOGL_MODE_3D_OTHER:
    if (OtherCallback != 0) {
      (*OtherCallback) ();
    }
    break;
  case AUTOGL_MODE_3D_TRANSLATE:
    Translate ();
    break;
  case AUTOGL_MODE_3D_ROTATE_XY:
    RotateXY ();
    break;
  case AUTOGL_MODE_3D_ROTATE_Z:
    RotateZ ();
    break;
  case AUTOGL_MODE_3D_SCALE:
    Scale ();
    break;
  default:
    assert(0);
    break;
  }
}

void AutoGL_SetPanelInMode3D (void)
{
  AutoGL_AddInteger (&Mode, "AutoGL_Mode3D"); 
  AutoGL_SetLabel ("mode");
  AutoGL_AddIntegerItem ("other");              
  AutoGL_AddIntegerItem ("translate");             
  AutoGL_AddIntegerItem ("rotate_direction");         
  AutoGL_AddIntegerItem ("rotate_up_vector");              
  AutoGL_AddIntegerItem ("scale");              
}

void AutoGL_SetMode3D (int mode)    /*  AUTOGL_MODE_3D_  */
{
  Mode = mode;
}

void AutoGL_SetDefaultCallbackInMode3D 
(AutoGL_CallbackType otherCallback)
{
  OtherCallback = otherCallback;
  AutoGL_SetViewEventCallback (DefaultViewEventCallback);
}

int AutoGL_GetMode3D (void)    /*  AUTOGL_MODE_3D_  */
{
  return Mode;
}

int AutoGL_PointingDeviceIsHit3D 
(double *parameter_OUT,
 double pX, double pY, double pZ,
 int toleranceDc)
{
  int dcX, dcY, dcZ;
  double originX, originY, originZ;
  double dirX, dirY, dirZ;
  double x, y, z;
  double ptX, ptY, ptZ;
  double t;
  double distance;
  double tolerance;

  AutoGL_GetPointingDevicePositionDc (&dcX, &dcY, &dcZ);
  AutoGL_GetPositionOfDc (&originX, &originY, &originZ,
			  dcX, dcY, dcZ);
  AutoGL_GetViewDirection (&dirX, &dirY, &dirZ);
  
  AutoGL_SubtractVector3D (&x, &y, &z,
			   pX, pY, pZ,
			   originX, originY, originZ);
  t = AutoGL_GetScalarProduct3D (x, y, z,
				 dirX, dirY, dirZ);
  *parameter_OUT = t;

  AutoGL_AddVector3D (&ptX, &ptY, &ptZ,
		      originX, originY, originZ,
		      dirX * t, dirY * t, dirZ * t);
  distance = AutoGL_GetVectorDistance3D (pX, pY, pZ,
					 ptX, ptY, ptZ);
  tolerance = toleranceDc 
    * AutoGL_GetViewSize () / AutoGL_GetViewSizeDc ();
  return (distance <= tolerance);
}



