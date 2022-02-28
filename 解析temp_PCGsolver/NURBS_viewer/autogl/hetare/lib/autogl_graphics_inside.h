/* header file for autogl_graphics.c */

#ifndef __AUTOGL_GRAPHICS_INSIDE__
#define __AUTOGL_GRAPHICS_INSIDE__



#define MAX_IMAGE_WIDTH 2048
#define MAX_IMAGE_HEIGHT 2048
#define MAX_PIXEL (MAX_IMAGE_WIDTH * MAX_IMAGE_HEIGHT)

void AutoGL_SetViewRangeDc (int widthDc, int heightDc);

void AutoGL_SetPointingDeviceEvent 
(int event,   /* AUTOGL_EVENT_ */
 int dcX, int dcY, int dcZ);

void AutoGL_SetKeyEvent (int keyChar);



/* OpenGL */
#include <GL/gl.h>
#include <GL/glu.h>

void AutoGL_RealizeOpenGLWindow (void);

void AutoGL_StartDrawingOpenGLWindow (void);

void AutoGL_EndDrawingOpenGLWindow (void);

GLuint AutoGL_GetOpenGLCharBase (void);



#endif  /* __AUTOGL_GRAPHICS_INSIDE__ */



