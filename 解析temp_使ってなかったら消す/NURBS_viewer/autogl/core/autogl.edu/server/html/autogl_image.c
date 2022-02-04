/* autogl_image.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "autogl_cgi.h"



#define MAX_IMAGE_SIZE (2048*2048*3)
static unsigned char Image[MAX_IMAGE_SIZE];

int main (int argc, char *argv[]) 
{
  FILE *fp; 
  int imageSize;
  
#if 0
  printf ("Content-type: image/gif\n\n");
#else
  printf ("Content-type: image/jpeg\n\n");
#endif
  
  fp = fopen (VIEW_IMAGE, "rb");
  assert(fp != NULL);
  
  imageSize = fread (Image, 1, MAX_IMAGE_SIZE, fp);
  fwrite (Image, 1, imageSize, stdout);
  
  fclose (fp);

  return 0;
}
