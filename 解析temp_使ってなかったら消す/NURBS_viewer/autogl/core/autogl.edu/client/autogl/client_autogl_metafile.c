/* client_autogl_metafile.c : client side AutoGL implementation */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/* usleep */
#include <unistd.h>

#include <autogl.h>



/* application model I/O */

#define MAX_COMMAND_LENGTH 1024
#define MAX_FILE_NAME_LENGTH 80
#define MAX_SERVER_HOST_NAME_LENGTH 80
#define MAX_SERVER_DIRECTORY_NAME_LENGTH 80

char ServerHostName[MAX_SERVER_HOST_NAME_LENGTH] 
= "who";
char ServerDirectoryName[MAX_SERVER_DIRECTORY_NAME_LENGTH] 
= "test_autogl_client";

#define VIEW_METAFILE "autogl_metafile"

#define COPY "/bin/cp"
#define DELETE "/bin/rm"
#define RENAME "/bin/mv"

#define REMOTE_COPY "/usr/bin/rcp"
#define REMOTE_EXECUTE "/usr/bin/rsh"

void GetFile (char *fileName) 
{
  char buf[MAX_COMMAND_LENGTH];

  assert(strlen (fileName) < MAX_FILE_NAME_LENGTH);
  assert(strlen (ServerHostName) < MAX_SERVER_HOST_NAME_LENGTH);
  assert(strlen (ServerDirectoryName) < MAX_SERVER_DIRECTORY_NAME_LENGTH);
  sprintf (buf, "%s %s:%s/%s .", 
	   REMOTE_COPY,
	   ServerHostName, ServerDirectoryName, fileName);
  system (buf);

  fprintf (stderr, "%s\n", buf);
}



static void RedrawView (void) 
{
  AutoGL_DrawCoordinateSystem3D (20);
  AutoGL_DrawDisplayList ();
}

static void UpdateModel (void) 
{
  GetFile (VIEW_METAFILE);
  AutoGL_ReadMetafile (VIEW_METAFILE, 0);
  AutoGL_DrawView ();  
}

void AutoGL_SetUp (int argc, char *argv[]) 
{
  GetFile (VIEW_METAFILE);
  AutoGL_ReadMetafile (VIEW_METAFILE, 1);

  AutoGL_SetViewRedrawCallback (RedrawView);
  AutoGL_SetDefaultCallbackInMode3D (NULL);    
  AutoGL_EnableDragInMode3D ();               
  AutoGL_SetPanelInMode3D ();

  AutoGL_AddCallback (UpdateModel, "UpdateModel");
}

