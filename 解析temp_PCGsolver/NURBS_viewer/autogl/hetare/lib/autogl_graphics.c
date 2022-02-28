/* autogl_graphics.c : graphics */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "autogl_math.h"
#include "autogl_uims.h"
#include "autogl_graphics.h"

#include "autogl_graphics_inside.h"



/* Module Variables */

static double View_size = 1.0;
static double View_centerX = 0.0;
static double View_centerY = 0.0;
static double View_centerZ = 0.0;
static double View_directionX = 0.0;
static double View_directionY = 0.0;
static double View_directionZ = 1.0;
static double View_upVectorX = 0.0;
static double View_upVectorY = 1.0;
static double View_upVectorZ = 0.0;

static int View_widthDc = 100;
static int View_heightDc = 100;

static int View_eventType;  /* AUTOGL_EVENT_ */

static int View_pointingDeviceDcX;
static int View_pointingDeviceDcY;
static int View_pointingDeviceDcZ;

/* background color of the view window */
static double View_backgroundRed;
static double View_backgroundGreen;
static double View_backgroundBlue;

/* ratio of Z direction (depth) against X and Y direction 
   in the viewing coordinate system of the clipping volume */
static const double View_depthRatio = 5.0;

/* whether lighting is on/off */
static int View_lightIsOn;

/* OpenGL list base for character string */
/* 128 OpenGL lists (No. 0 - 127) are reserved for drawing character string. */
static GLuint View_charBase;



/* View */

void AutoGL_SetViewSize (double size)
{
  View_size = size;
}

void AutoGL_SetViewCenter (double x, double y, double z)
{
  View_centerX = x;
  View_centerY = y;
  View_centerZ = z;
}

void AutoGL_SetViewDirection (double x, double y, double z)
{
  AutoGL_NormalizeVector3D 
    (&View_directionX, &View_directionY, &View_directionZ,
     x, y, z);
}

void AutoGL_SetViewUpVector (double x, double y, double z)
{
  AutoGL_NormalizeVector3D 
    (&View_upVectorX, &View_upVectorY, &View_upVectorZ,
     x, y, z);
}

void AutoGL_SetViewRangeDc (int widthDc, int heightDc)
{
  View_widthDc = (widthDc <= 0) ? 1 : widthDc;
  View_heightDc = (heightDc <= 0) ? 1 : heightDc;
}

double AutoGL_GetViewSize (void)
{
  return View_size;
}

void AutoGL_GetViewCenter 
(double *x_OUT, double *y_OUT, double *z_OUT)
{
  *x_OUT = View_centerX;
  *y_OUT = View_centerY;
  *z_OUT = View_centerZ;
}

void AutoGL_GetViewDirection 
(double *x_OUT, double *y_OUT, double *z_OUT)
{
  *x_OUT = View_directionX;
  *y_OUT = View_directionY;
  *z_OUT = View_directionZ;
}

void AutoGL_GetViewUpVector 
(double *x_OUT, double *y_OUT, double *z_OUT)
{
  *x_OUT = View_upVectorX; 
  *y_OUT = View_upVectorY;
  *z_OUT = View_upVectorZ;
}

void AutoGL_GetViewRangeDc 
(int *widthDc_OUT, int *heightDc_OUT)
{
  *widthDc_OUT = View_widthDc;
  *heightDc_OUT = View_heightDc;
}

int AutoGL_GetViewSizeDc (void)
{
  return (View_widthDc < View_heightDc) ?
    View_widthDc / 2 : View_heightDc / 2;
}

void AutoGL_GetVcOfDc 
(double *vcX_OUT, double *vcY_OUT, double *vcZ_OUT,
 int dcX, int dcY, int dcZ)
{
  double size = AutoGL_GetViewSize ();
  int sizeDc = AutoGL_GetViewSizeDc ();

  *vcX_OUT = size * dcX / sizeDc;
  *vcY_OUT = size * dcY / sizeDc;
  *vcZ_OUT = size * dcZ / sizeDc;
}

void AutoGL_GetPositionOfVc 
(double *x_OUT, double *y_OUT, double *z_OUT,
 double vcX, double vcY, double vcZ)
{
  double centerX, centerY, centerZ;
  double zAxisX, zAxisY, zAxisZ;
  double upX, upY, upZ;
  double xAxisX, xAxisY, xAxisZ;
  double yAxisX, yAxisY, yAxisZ;

  AutoGL_GetViewCenter (&centerX, &centerY, &centerZ);
  AutoGL_GetViewDirection (&zAxisX, &zAxisY, &zAxisZ);
  AutoGL_GetViewUpVector (&upX, &upY, &upZ);

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
  
  AutoGL_AddVector3D (x_OUT, y_OUT, z_OUT,
		      centerX, centerY, centerZ,
		      xAxisX * vcX + yAxisX * vcY + zAxisX * vcZ,
		      xAxisY * vcX + yAxisY * vcY + zAxisY * vcZ,
		      xAxisZ * vcX + yAxisZ * vcY + zAxisZ * vcZ);
}

void AutoGL_GetPositionOfDc 
(double *x_OUT, double *y_OUT, double *z_OUT,
 int dcX, int dcY, int dcZ)
{
  double vcX, vcY, vcZ;

  AutoGL_GetVcOfDc (&vcX, &vcY, &vcZ,
		    dcX, dcY, dcZ);
  AutoGL_GetPositionOfVc (x_OUT, y_OUT, z_OUT,
			  vcX, vcY, vcZ);
}



/* Event Handling */

void AutoGL_SetPointingDeviceEvent 
(int event,   /* AUTOGL_EVENT_ */
 int dcX, int dcY, int dcZ)
{
  View_eventType = event;
  View_pointingDeviceDcX = dcX;
  View_pointingDeviceDcY = dcY;
  View_pointingDeviceDcZ = dcZ;
}

int AutoGL_GetViewEvent (void)  /* AUTOGL_EVENT_ */
{
  return View_eventType;
}

void AutoGL_GetPointingDevicePositionDc 
(int *dcX_OUT, int *dcY_OUT, int *dcZ_OUT)
{
  *dcX_OUT = View_pointingDeviceDcX;
  *dcY_OUT = View_pointingDeviceDcY;
  *dcZ_OUT = View_pointingDeviceDcZ;
}



static void SetProjection (void)
{
  GLfloat size;
  int widthDc, heightDc;
  GLfloat halfWidth, halfHeight;

  size = AutoGL_GetViewSize ();
  AutoGL_GetViewRangeDc 
    (&widthDc, &heightDc);

  if (widthDc <= heightDc) {
    halfWidth = size;
    halfHeight = size * heightDc / widthDc;
  } else {
    halfWidth = size * widthDc / heightDc;
    halfHeight = size;
  }

  /* set orthogonal projection */
  glOrtho (-halfWidth, halfWidth, 
	   -halfHeight, halfHeight, 
	   -size * View_depthRatio, size * View_depthRatio);
} 

static void SetViewingTransformation (void)
{
  double centerX, centerY, centerZ;
  double dirX, dirY, dirZ;
  double upX, upY, upZ;
  GLfloat size;

  AutoGL_GetViewCenter 
    (&centerX, &centerY, &centerZ);
  AutoGL_GetViewDirection 
    (&dirX, &dirY, &dirZ);
  AutoGL_GetViewUpVector 
    (&upX, &upY, &upZ);
  size = AutoGL_GetViewSize ();

  /* set viewing transformation */
  gluLookAt (centerX + dirX, 
	     centerY + dirY, 
	     centerZ + dirZ,
	     centerX, centerY, centerZ,
	     upX, upY, upZ);
}

static void TurnOnLighting (void)
{
  if (!View_lightIsOn) {

    /* turn on lighting */
    glEnable (GL_LIGHTING);

    View_lightIsOn = 1;
  }
}

static void TurnOffLighting (void)
{
  if (View_lightIsOn) {

    /* turn off lighting */
    glDisable (GL_LIGHTING);

    View_lightIsOn = 0;
  }
}

void AutoGL_SetBackgroundColor 
(double red, double green, double blue)
{
  View_backgroundRed = red;
  View_backgroundGreen = green;
  View_backgroundBlue = blue;
}

void AutoGL_RealizeOpenGLWindow (void)
{
  GLfloat ambientLightColor[] = { 0.2f, 0.2f, 0.2f, 1.0f };
  GLfloat diffuseLightColor[] = { 0.8f, 0.8f, 0.8f, 1.0f };
  GLfloat lightPosition[4];

  /* set background color */
  glClearColor (View_backgroundRed, 
		View_backgroundGreen, 
		View_backgroundBlue, 
		1.0f);
  
  /* set foreground color : white */
  glColor3f (1.0f, 1.0f, 1.0f);
  
  /* turn on depth test */
  glEnable (GL_DEPTH_TEST);
  
  /* turn on face culling */
  glFrontFace (GL_CCW);
  glEnable (GL_CULL_FACE);
  
  /* turn on lighting */
  glEnable (GL_LIGHTING);
  View_lightIsOn = 1;
  
  /* set ambient light */
  glLightModelfv (GL_LIGHT_MODEL_AMBIENT, ambientLightColor);
  
  /* turn on material color tracking */
  glEnable (GL_COLOR_MATERIAL);
  
  /* specify ambient and diffuse as tracking color */
  glColorMaterial (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
  
  /* set color at light No. 0 */
  glLightfv (GL_LIGHT0, GL_DIFFUSE, diffuseLightColor);
  
  /* set position at light No. 0 : directional light */
  lightPosition[0] = 0.0f;
  lightPosition[1] = 0.0f;
  lightPosition[2] = 1.0f;
  lightPosition[3] = 0.0f;
  glLightfv (GL_LIGHT0, GL_POSITION, lightPosition);

  /* turn on light No. 0 */
  glEnable (GL_LIGHT0);

  /* allocate lists of bitmaps for OpenGL character string drawing */
  View_charBase = glGenLists (128);

  /* turn on edge contrasting */
  glEnable (GL_POLYGON_OFFSET_FILL);
  glPolygonOffset( 1.0, 1.0 );
} 

void AutoGL_StartDrawingOpenGLWindow (void)
{
  /* begin rendering OpenGL */
  
  /* enable Z buffer as writable */
  glDepthMask (GL_TRUE);
  
  /* turn on lighting */
  glEnable (GL_LIGHTING);
  View_lightIsOn = 1;
  
  /* set OpenGL viewport */
  {
    int widthDc, heightDc;

    AutoGL_GetViewRangeDc 
      (&widthDc, &heightDc);
    glViewport (0, 0, (GLsizei)widthDc, (GLsizei)heightDc);
  }

  /* reset projection matrix stack */
  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
  
  SetProjection ();

  /* set background color */
  glClearColor (View_backgroundRed, 
		View_backgroundGreen, 
		View_backgroundBlue, 
		1.0f);

  /* clear color and depth buffers */
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  /* reset model view matrix stack */
  glMatrixMode (GL_MODELVIEW);
  glLoadIdentity ();

  SetViewingTransformation ();

  {
    double dirX, dirY, dirZ;

    /* set initial normal */
    AutoGL_GetViewDirection 
      (&dirX, &dirY, &dirZ);
    glNormal3d (dirX, dirY, dirZ); 

    /* set initial color : white */
    glColor4f (1.0, 1.0, 1.0, 1.0); 
  }
}

void AutoGL_EndDrawingOpenGLWindow (void)
{
  /* end rendering OpenGL */

  /* flush drawing commands */
  glFlush ();
}

GLuint AutoGL_GetOpenGLCharBase (void)
{
  return View_charBase;
}

void AutoGL_SetColor 
(double red, double green, double blue)
{
  glColor4f (red, green, blue, 1.0f);
}

void AutoGL_DrawLine 
(double x0, double y0, double z0,
 double x1, double y1, double z1)
{
  TurnOffLighting ();

  /* draw a line segment */
  glBegin (GL_LINES);

  glVertex3d (x0, y0, z0);
  glVertex3d (x1, y1, z1);

  glEnd ();
}

void AutoGL_DrawTriangle 
(double x0, double y0, double z0,
 double x1, double y1, double z1,
 double x2, double y2, double z2)
{
  double nx, ny, nz;

  /* calculate the unit normal vector of the triangle */
  AutoGL_GetVectorProduct3D
    (&nx, &ny, &nz,
     x0 - x1, y0 - y1, z0 - z1,
     x1 - x2, y1 - y2, z1 - z2);
  AutoGL_NormalizeVector3D 
    (&nx, &ny, &nz,
     nx, ny, nz);

  TurnOnLighting ();

  /* draw two triangles (front side and back side) */

  glBegin (GL_TRIANGLES);

  glNormal3d (nx, ny, nz); 
  glVertex3d (x0, y0, z0);
  glVertex3d (x1, y1, z1);
  glVertex3d (x2, y2, z2);

  glNormal3d (-nx, -ny, -nz); 
  glVertex3d (x0, y0, z0);
  glVertex3d (x2, y2, z2);
  glVertex3d (x1, y1, z1);

  glEnd ();
}

void AutoGL_DrawColorInterpolatedTriangle
(double x0, double y0, double z0, double red0, double green0, double blue0,  
 double x1, double y1, double z1, double red1, double green1, double blue1,  
 double x2, double y2, double z2, double red2, double green2, double blue2)
{
  double nx, ny, nz;

  /* calculate the unit normal vector of the triangle */
  AutoGL_GetVectorProduct3D
    (&nx, &ny, &nz,
     x0 - x1, y0 - y1, z0 - z1,
     x1 - x2, y1 - y2, z1 - z2);
  AutoGL_NormalizeVector3D 
    (&nx, &ny, &nz,
     nx, ny, nz);

  TurnOnLighting ();

  /* draw two triangles (front side and back side) */

  glBegin (GL_TRIANGLES);

  glNormal3d (nx, ny, nz); 
  glColor3d (red0, green0, blue0);
  glVertex3d (x0, y0, z0);
  glColor3d (red1, green1, blue1);
  glVertex3d (x1, y1, z1);
  glColor3d (red2, green2, blue2);
  glVertex3d (x2, y2, z2);

  glNormal3d (-nx, -ny, -nz); 
  glColor3d (red0, green0, blue0);
  glVertex3d (x0, y0, z0);
  glColor3d (red2, green2, blue2);
  glVertex3d (x2, y2, z2);
  glColor3d (red1, green1, blue1);
  glVertex3d (x1, y1, z1);

  glEnd ();
}

void AutoGL_DrawQuadrangle 
(double x0, double y0, double z0,
 double x1, double y1, double z1,
 double x2, double y2, double z2,
 double x3, double y3, double z3)
{
  AutoGL_DrawTriangle (x0, y0, z0,
		       x1, y1, z1,
		       x2, y2, z2);
  AutoGL_DrawTriangle (x2, y2, z2,
		       x3, y3, z3,
		       x0, y0, z0);
}

void AutoGL_DrawString 
(double x, double y, double z,
 const char *str)
{
  TurnOffLighting ();

  /* draw character string at the point */

  glRasterPos3d (x, y, z);

  glPushAttrib (GL_LIST_BIT);
  glListBase (View_charBase);
  glCallLists (strlen (str), GL_UNSIGNED_BYTE, str);
  glPopAttrib ();
}

