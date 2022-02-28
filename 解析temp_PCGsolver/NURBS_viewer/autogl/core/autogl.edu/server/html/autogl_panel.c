/* autogl_panel.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "autogl_cgi.h"



int main (int argc, char *argv[]) 
{
  char line[1024];
  FILE *fp; 
  
  printf ("Content-type:text/html\n\n");
  printf ("<html>\n");



  /* header part */

  printf ("<head>\n");

  /* title */
  printf ("<title>Panel</title>\n");

  printf ("</head>\n");



  /* body part : image and hidden form */

  printf ("<body>\n");

  /* form */
  printf ("<form name=\"panelForm\" method=\"POST\" action=\"./autogl_frame.cgi\" target=\"_top\">");

  fp = fopen (PANEL_HTML, "r");
  assert(fp != NULL);
  while (!feof (fp)) {
    char *ret;

    ret = fgets (line, 1024, fp);
    if (ret == NULL) break;
    fputs (line, stdout);
  }
  fclose (fp);

  printf ("<input type=\"submit\" value=\"AutoGL_Update\" ><br>\n");
  printf ("</form>\n");

  printf ("</body>\n");



  printf ("</html>\n");

  return 0;
}
