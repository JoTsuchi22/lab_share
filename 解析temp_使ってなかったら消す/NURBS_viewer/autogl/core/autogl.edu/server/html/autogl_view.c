/* autogl_view.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include "autogl_cgi.h"



int main (int argc, char *argv[]) 
{
  printf ("Content-type:text/html\n\n");
  printf ("<html> \n");



  /* header part */

  printf ("<head> \n");

  /* title */
  printf ("<title>View</title> \n");

  /* script */
  printf ("<script type=\"text/javascript\"> \n");
  printf ("function load () \n");
  printf ("{ \n");
  printf ("  isIE = document.all; \n");
  printf ("  isNS = document.getElementById; \n");
  printf ("  if (isIE) { \n");
  printf ("    document.onmouseup = mouseUpIE; \n");
  printf ("  } else if (isNS) { \n");
  printf ("    window.addEventListener (\"mouseup\", mouseUpNS, true); \n");
  printf ("  } \n");
  printf ("} \n");
  printf ("function mouseUpIE () \n");
  printf ("{ \n");
  printf ("  if (document.hiddenForm.AutoGL_ViewMouseClicked.value == 0) { \n");
  printf ("    document.hiddenForm.AutoGL_ViewMouseClicked.value = 1; \n");
  printf ("    document.hiddenForm.AutoGL_ViewMouseStartX.value = window.event.clientX - 10; \n");
  printf ("    document.hiddenForm.AutoGL_ViewMouseStartY.value = window.event.clientY - 10; \n");
  printf ("  } else { \n");
  printf ("    document.hiddenForm.AutoGL_ViewMouseClicked.value = 0; \n");
  printf ("    document.hiddenForm.AutoGL_ViewMouseEndX.value = window.event.clientX - 10; \n");
  printf ("    document.hiddenForm.AutoGL_ViewMouseEndY.value = window.event.clientY - 10; \n");
  printf ("    document.hiddenForm.submit (); \n");
  printf ("  } \n");
  printf ("} \n");
  printf ("function mouseUpNS (event) \n");
  printf ("{ \n");
  printf ("  if (document.hiddenForm.AutoGL_ViewMouseClicked.value == 0) { \n");
  printf ("    document.hiddenForm.AutoGL_ViewMouseClicked.value = 1; \n");
  printf ("    document.hiddenForm.AutoGL_ViewMouseStartX.value = event.clientX - 10; \n");
  printf ("    document.hiddenForm.AutoGL_ViewMouseStartY.value = event.clientY - 10; \n");
  printf ("  } else { \n");
  printf ("    document.hiddenForm.AutoGL_ViewMouseEndX.value = event.clientX - 10; \n");
  printf ("    document.hiddenForm.AutoGL_ViewMouseEndY.value = event.clientY - 10; \n");
  printf ("    document.hiddenForm.submit (); \n");
  printf ("  } \n");
  printf ("} \n");
  printf ("</script> \n");

  printf ("</head> \n");



  /* body part : image and hidden form */

  printf ("<body onLoad=\"load()\"> \n");

  /* image */
  printf ("<img name=\"viewImage\" src=\"autogl_image.cgi\" > \n");

  /* hidden form */
  printf ("<form name=\"hiddenForm\" method=\"POST\" action=\"./autogl_frame.cgi\" target=\"_top\" > \n");
  printf ("<input type=\"hidden\" name=\"AutoGL_ViewMouseClicked\" value=\"0\" ><br> \n");
  printf ("<input type=\"hidden\" name=\"AutoGL_ViewMouseStartX\" value=\"0\" ><br> \n");
  printf ("<input type=\"hidden\" name=\"AutoGL_ViewMouseStartY\" value=\"0\" ><br> \n");
  printf ("<input type=\"hidden\" name=\"AutoGL_ViewMouseEndX\" value=\"0\" ><br> \n");
  printf ("<input type=\"hidden\" name=\"AutoGL_ViewMouseEndY\" value=\"0\" ><br> \n");
  printf ("</form> \n");

  printf ("</body> \n");



  printf ("</html> \n");

  return 0;
}
