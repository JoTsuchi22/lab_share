/* autogl_post_command.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

/* usleep */
#include <unistd.h>



#define COPY "/bin/cp"
#define DELETE "/bin/rm"
#define RENAME "/bin/mv"
#define CHANGE_MODE "/bin/chmod og+w"

#define INFO_COMMAND_OUTPUT "autogl_info_command_output_file"
#define VALUE_COMMAND_OUTPUT "autogl_value_command_output_file"
#define POST_COMMAND "autogl_post_command_file"
#define POST_COMMAND_OUTPUT "autogl_post_command_output_file"

#if 0
#define VIEW_IMAGE "autogl_view_image.gif"
#else
#define VIEW_IMAGE "autogl_view_image.jpeg"
#endif

#define MAX_CGI_POST_LENGTH 10000

int main (int argc, char *argv[]) 
{
  static char buf[1024];
  FILE *in;
  FILE *out;
  int ch;
  int length;

  sprintf (buf, "%s_TMP",
	   POST_COMMAND);
  out = fopen (buf, "w");
  length = 0;
  while ((ch = getchar ()) != EOF) {
    fputc (ch, out);
    length++;
    if (MAX_CGI_POST_LENGTH <= length) break;
  }
  fclose (out);
  
  sprintf (buf, "%s %s_TMP %s_BAK",
	   COPY,
	   POST_COMMAND, 
	   POST_COMMAND);
  system (buf);

  sprintf (buf, "%s %s_TMP",
	   CHANGE_MODE,
	   POST_COMMAND);
  system (buf);
  sprintf (buf, "%s %s_TMP %s",
	   RENAME,
	   POST_COMMAND, 
	   POST_COMMAND);
  system (buf);



#if 0
  /* NFS environment */
  for (;;) {
    system ("/bin/ls > /dev/null");
    in = fopen (POST_COMMAND_OUTPUT, "r");
    if (in != NULL) break;
    usleep (500000);
  }
#else
  /* stand alone environment */
  for (;;) {
    in = fopen (POST_COMMAND_OUTPUT, "r");
    if (in != NULL) break;
    usleep (100000);
  }
#endif


  printf ("Content-Type: text/plain\n");
  printf ("\n");
  while ((ch = fgetc (in)) != EOF) {
    putchar (ch);
  }
  fclose (in);

  sprintf (buf, "%s %s %s",
	   COPY,
	   POST_COMMAND_OUTPUT,
	   VALUE_COMMAND_OUTPUT);
  system (buf);

  sprintf (buf, "%s %s",
	   DELETE, 
	   POST_COMMAND_OUTPUT);
  system (buf);

  return 0;
}
