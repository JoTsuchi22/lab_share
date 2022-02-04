/* autogl_frame.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <ctype.h>

/* usleep */
#include <unistd.h>

#include "autogl_cgi.h"



#define MAX_WORD_LENGTH 256



/* Reader */

FILE *ReaderFp;
int CurrentChar;

void InitializeReader (FILE *fp) 
{
  ReaderFp = fp;
  CurrentChar = fgetc (ReaderFp);
}

void TerminateReader (void) 
{
  fclose (ReaderFp);
  ReaderFp = NULL;
}

char *ReadWord (void) 
{
  static char word[MAX_WORD_LENGTH];
  int wordLength;

  while (isspace (CurrentChar)) {
    CurrentChar = fgetc (ReaderFp);
  }
  if (CurrentChar == EOF) return NULL;
  
  wordLength = 0;
  if (CurrentChar == '"') {
    CurrentChar = fgetc (ReaderFp);
    for (;;) {
      assert (CurrentChar != EOF);
      if (CurrentChar == '"') {
	CurrentChar = fgetc (ReaderFp);
	break;
      }
      
      word[wordLength++] = CurrentChar;
      CurrentChar = fgetc (ReaderFp);
    }
    word[wordLength] = 0;
  } else {
    for (;;) {
      if (isspace (CurrentChar)) break;
      if (CurrentChar == EOF) break;
      
      word[wordLength++] = CurrentChar;
      CurrentChar = fgetc (ReaderFp);
    }
    word[wordLength] = 0;
  }
  return word;
}
    


/* Panel */

typedef enum {
  PROPERTY_GROUP,
  PROPERTY_COMMENT,
  PROPERTY_CALLBACK,
  PROPERTY_BOOLEAN,
  PROPERTY_STRING,
  PROPERTY_INTEGER,
  PROPERTY_REAL
} PropertyType;

typedef enum {
  WIDGET_LABEL,
  WIDGET_BUTTON,
  WIDGET_TOGGLE,
  WIDGET_FIELD,
  WIDGET_CHOICE
} WidgetType;

#define MAX_WIDGETS 100
#define MAX_WIDGET_NAME_LENGTH 80
#define MAX_WIDGET_LABEL_LENGTH 80
#define MAX_WIDGET_VALUE_LENGTH 1024
#define MAX_WIDGET_ITEMS 10
#define MAX_WIDGET_ITEM_LABEL_LENGTH 80

int NWidgets;
WidgetType WidgetTypes[MAX_WIDGETS];
PropertyType WidgetPropertyTypes[MAX_WIDGETS];
char WidgetNames[MAX_WIDGETS][MAX_WIDGET_NAME_LENGTH];
char WidgetLabels[MAX_WIDGETS][MAX_WIDGET_LABEL_LENGTH];
int NWidgetItems[MAX_WIDGETS];
char WidgetItemLabels[MAX_WIDGETS]
[MAX_WIDGET_ITEMS][MAX_WIDGET_ITEM_LABEL_LENGTH];
char WidgetValues[MAX_WIDGETS][MAX_WIDGET_VALUE_LENGTH];



/* View */

int ViewMouseIsDragged;
int ViewMouseStartX;
int ViewMouseStartY;
int ViewMouseEndX;
int ViewMouseEndY;



void ReadInfoCommandOutputFile (FILE *fp) 
{
  InitializeReader (fp);

  NWidgets = 0;
  for (;;) {
    char *keywordPtr;
    char keyword[80];
    char *nameKeyword;
    char *name;
    char *labelKeyword;
    char *label;

    assert(NWidgets < MAX_WIDGETS);

    keywordPtr = ReadWord ();
    if (keywordPtr == NULL) break;
    strcpy (keyword, keywordPtr);
    
    nameKeyword = ReadWord ();
    assert(nameKeyword != NULL);
    assert(strcmp (nameKeyword, "name") == 0);
    
    name = ReadWord ();
    assert(name != NULL);
    strcpy (WidgetNames[NWidgets], name);
    
    labelKeyword = ReadWord ();
    assert(labelKeyword != NULL);
    assert(strcmp (labelKeyword, "label") == 0);
    
    label = ReadWord ();
    assert(label != NULL);
    strcpy (WidgetLabels[NWidgets], label);

    if (strcmp (keyword, "group") == 0) {
      WidgetTypes[NWidgets] = WIDGET_LABEL;
      WidgetPropertyTypes[NWidgets] = PROPERTY_GROUP;
      NWidgets++;

    } else if (strcmp (keyword, "comment") == 0) {
      WidgetTypes[NWidgets] = WIDGET_LABEL;
      WidgetPropertyTypes[NWidgets] = PROPERTY_COMMENT;
      NWidgets++;

    } else if (strcmp (keyword, "callback") == 0) {
      WidgetTypes[NWidgets] = WIDGET_BUTTON;
      WidgetPropertyTypes[NWidgets] = PROPERTY_CALLBACK;
      NWidgets++;

    } else if (strcmp (keyword, "boolean") == 0) {
      WidgetTypes[NWidgets] = WIDGET_TOGGLE;
      WidgetPropertyTypes[NWidgets] = PROPERTY_BOOLEAN;
      NWidgets++;

    } else if (strcmp (keyword, "string") == 0) {
      WidgetTypes[NWidgets] = WIDGET_FIELD;
      WidgetPropertyTypes[NWidgets] = PROPERTY_STRING;
      NWidgets++;

    } else if (strcmp (keyword, "integer") == 0) {
      char *nItemsKeyword;
      char *nItemsStr;
      int nItems;
      int iItem;
      
      nItemsKeyword = ReadWord ();
      assert(nItemsKeyword != NULL);
      assert(strcmp (nItemsKeyword, "nItems") == 0);
      
      nItemsStr = ReadWord ();
      assert(nItemsStr != NULL);
      nItems = atoi (nItemsStr);
      assert(nItems <= MAX_WIDGET_ITEMS);
      NWidgetItems[NWidgets] = nItems;

      if (nItems == 0) {
	WidgetTypes[NWidgets] = WIDGET_FIELD;
	WidgetPropertyTypes[NWidgets] = PROPERTY_INTEGER;

	NWidgets++;

      } else {
	WidgetTypes[NWidgets] = WIDGET_CHOICE;
	WidgetPropertyTypes[NWidgets] = PROPERTY_INTEGER;

	for (iItem = 0; iItem < nItems; iItem++) {
	  char *itemLabelKeyword;
	  char *itemLabel;
	  
	  itemLabelKeyword = ReadWord ();
	  assert(itemLabelKeyword != NULL);
	  assert(strcmp (itemLabelKeyword, "itemLabel") == 0);
	  
	  itemLabel = ReadWord ();
	  assert(itemLabel != NULL);
	  strcpy (WidgetItemLabels[NWidgets][iItem], itemLabel);
	}

	NWidgets++;
      }

    } else if (strcmp (keyword, "real") == 0) {
      WidgetTypes[NWidgets] = WIDGET_FIELD;
      WidgetPropertyTypes[NWidgets] = PROPERTY_REAL;
      NWidgets++;

    } else {
      assert(0);
    }
  }

  TerminateReader ();
}

void ReadValueCommandOutputFile (FILE *fp) 
{
  InitializeReader (fp);

  for (;;) {
    char *keywordPtr;
    char keyword[80];
    char *nameKeyword;
    char *name;
    int iWidget;
    int widgetId;
    char *valueKeyword;
    char *value;

    assert(NWidgets < MAX_WIDGETS);

    keywordPtr = ReadWord ();
    if (keywordPtr == NULL) break;
    strcpy (keyword, keywordPtr);
    
    nameKeyword = ReadWord ();
    assert(nameKeyword != NULL);
    assert(strcmp (nameKeyword, "name") == 0);
    
    name = ReadWord ();
    assert(name != NULL);
    
    if (strcmp (keyword, "group") == 0
	|| strcmp (keyword, "comment") == 0
	|| strcmp (keyword, "callback") == 0) continue;

    widgetId = -1;
    for (iWidget = 0; iWidget < NWidgets; iWidget++) {
      if (strcmp (WidgetNames[iWidget], name) == 0) {
	widgetId = iWidget;
	break;
      }
    }
    assert(0 <= widgetId);

    valueKeyword = ReadWord ();
    assert(valueKeyword != NULL);
    assert(strcmp (valueKeyword, "value") == 0);
    
    value = ReadWord ();
    assert(value != NULL);
    strcpy (WidgetValues[widgetId], value);

    if (WidgetTypes[widgetId] == WIDGET_CHOICE) {
      ReadWord ();
    }
  }

  TerminateReader ();
}

void WritePanelHtmlFile (void)
{
  FILE *fp;
  int iWidget;

  fp = fopen (PANEL_HTML, "w");
  assert(fp != NULL);

  for (iWidget = 0; iWidget < NWidgets; iWidget++) {
    char *name = WidgetNames[iWidget];
    char *label = WidgetLabels[iWidget];
    char *value = WidgetValues[iWidget];
    
    if (strlen (label) == 0) {
      label = WidgetNames[iWidget];
    }

    switch (WidgetPropertyTypes[iWidget]) {
    case PROPERTY_GROUP:
      fprintf 
	(fp, "<h2>%s</h2><br>", 
	 label);
      break;
    case PROPERTY_COMMENT:
      fprintf 
	(fp, "%s<br>", 
	 label);
      break;
    case PROPERTY_CALLBACK:
      fprintf 
	(fp, "<input type=\"submit\" name=\"%s\" value=\"%s\" ><br>\n", 
	 name, label);
      break;
    case PROPERTY_BOOLEAN:
      if (strcmp (value, "true") == 0) {
	fprintf 
	  (fp, "<input type=\"checkbox\" name=\"%s\" checked > %s<br>\n",
	   name, label);
      } else if (strcmp (value, "false") == 0) {
	fprintf 
	  (fp, "<input type=\"checkbox\" name=\"%s\" > %s<br>\n",
	   name, label);
      } else {
	assert(0);
      }
      break;
    case PROPERTY_STRING:
      {
	/* sanitization */
	char buf[MAX_WORD_LENGTH];
	int count;
	int iChar;
	
	count = 0;
	for (iChar = 0; iChar < strlen (value); iChar++) {
	  int ch = value[iChar];

	  if (ch == 0) break;
	  if (MAX_WORD_LENGTH - 10 <= count) break;
	  if (ch == '&') {
	    buf[count++] = '&';
	    buf[count++] = 'a';
	    buf[count++] = 'm';
	    buf[count++] = 'p';
	    buf[count++] = ';';
	  } else if (ch == '<') {
	    buf[count++] = '&';
	    buf[count++] = 'l';
	    buf[count++] = 't';
	    buf[count++] = ';';
	  } else if (ch == '>') {
	    buf[count++] = '&';
	    buf[count++] = 'g';
	    buf[count++] = 't';
	    buf[count++] = ';';
	  } else if (ch == '\"') {
	    buf[count++] = '&';
	    buf[count++] = 'q';
	    buf[count++] = 'u';
	    buf[count++] = 'o';
	    buf[count++] = 't';
	    buf[count++] = ';';
	  } else if (ch == '\'') {
	    buf[count++] = '&';
	    buf[count++] = '#';
	    buf[count++] = '3';
	    buf[count++] = '9';
	    buf[count++] = ';';
	  } else {
	    buf[count++] = ch;
	  }
	}
	buf[count] = 0;

	fprintf 
	  (fp,"%s <input type=\"text\" name=\"%s\" value=\"%s\" ><br>\n",
	   label, name, buf);
      }
      break;
    case PROPERTY_INTEGER:
      if (0 < NWidgetItems[iWidget]) {
	int nItems = NWidgetItems[iWidget];
	int iItem;

	fprintf 
	  (fp,"%s<br>\n", 
	   label);

	for (iItem = 0; iItem < nItems; iItem++) {
	  char *checkMessage = "";
	  int integerValue;
	  char *itemLabel = WidgetItemLabels[iWidget][iItem];

	  integerValue = atoi (value);
	  if (iItem == integerValue) {
	    checkMessage = "checked";
	  }
	  
	  fprintf 
	    (fp, "<input type=\"radio\" name=\"%s\" value=\"%s\" %s > %s<br>\n",
	     name, itemLabel, checkMessage, itemLabel);
	}
      } else {
	fprintf 
	  (fp,"%s <input type=\"text\" name=\"%s\" value=\"%s\" ><br>\n",
	   label, name, value);
      }
      break;
    case PROPERTY_REAL:
      {
	fprintf 
	  (fp,"%s <input type=\"text\" name=\"%s\" value=\"%s\" ><br>\n",
	   label, name, value);
      }
      break;
    default:
      assert(0);
      break;
    }
  }

  fclose (fp);
}

void SetViewVariable (char *name, char *valueString)
{
  int flag;
  int value;

  if (strcmp (name, "AutoGL_ViewMouseStartX") == 0) {
    flag = sscanf (valueString, "%d", 
		   &value);
    assert(flag == 1);
    ViewMouseIsDragged = 1;
    ViewMouseStartX = value;
  } else if (strcmp (name, "AutoGL_ViewMouseStartY") == 0) {
    flag = sscanf (valueString, "%d", 
		   &value);
    assert(flag == 1);
    ViewMouseIsDragged = 1;
    ViewMouseStartY = value;
  } else if (strcmp (name, "AutoGL_ViewMouseEndX") == 0) {
    flag = sscanf (valueString, "%d", 
		   &value);
    assert(flag == 1);
    ViewMouseIsDragged = 1;
    ViewMouseEndX = value;
  } else if (strcmp (name, "AutoGL_ViewMouseEndY") == 0) {
    flag = sscanf (valueString, "%d", 
		   &value);
    assert(flag == 1);
    ViewMouseIsDragged = 1;
    ViewMouseEndY = value;
  }
}

void WriteResizeCommand (FILE *out)
{
  fprintf (out, "resize %d %d\n",
	   VIEW_IMAGE_WIDTH, VIEW_IMAGE_HEIGHT);
}

void WriteViewCommand (FILE *out)
{
  fprintf (out, "drag %d %d %d %d\n",
	   ViewMouseStartX,
	   ViewMouseStartY,
	   ViewMouseEndX,
	   ViewMouseEndY);
}

void ClearAllToggles (FILE *out)
{
  int iWidget;

  for (iWidget = 0; iWidget < NWidgets; iWidget++) {
    if (WidgetPropertyTypes[iWidget] == PROPERTY_BOOLEAN) {
      fprintf (out, "set %s false\n",
	       WidgetNames[iWidget]);
    }
  }
}

void WriteSetCommand (FILE *out, char *name, char *valueString)
{
  int widgetId;
  int iWidget;

  widgetId = -1;
  for (iWidget = 0; iWidget < NWidgets; iWidget++) {
    if (strcmp (WidgetNames[iWidget], name) == 0) {
      widgetId = iWidget;
      break;
    }
  }
  assert(0 <= widgetId);
  
  switch (WidgetPropertyTypes[widgetId]) {
  case PROPERTY_GROUP:
  case PROPERTY_COMMENT:
    assert(0);
    break;
  case PROPERTY_CALLBACK:
    fprintf (out, "call %s\n",
	     name);
    break;
  case PROPERTY_BOOLEAN:
    if (strcmp (valueString, "on") == 0) {
      fprintf (out, "set %s true\n",
	       name);
    } else {
      assert(0);
    }
    break;
  case PROPERTY_STRING:
    fprintf (out, "set %s \"%s\"\n",
	     name, valueString);
    break;
  case PROPERTY_INTEGER:
  case PROPERTY_REAL:
    fprintf (out, "set %s %s\n",
	     name, valueString);
    break;
  default:
    assert(0);
    break;
  }
}

int HexCode (int ch)
{
  switch (ch) {
  case '0': return 0;
  case '1': return 1;
  case '2': return 2;
  case '3': return 3;
  case '4': return 4;
  case '5': return 5;
  case '6': return 6;
  case '7': return 7;
  case '8': return 8;
  case '9': return 9;
  case 'a': return 10;
  case 'b': return 11;
  case 'c': return 12;
  case 'd': return 13;
  case 'e': return 14;
  case 'f': return 15;
  case 'A': return 10;
  case 'B': return 11;
  case 'C': return 12;
  case 'D': return 13;
  case 'E': return 14;
  case 'F': return 15;
  default:
    break;
  }
  return 0;
}

void ReadCgiPostFile (FILE *out, int viewFlag)
{
  FILE *in;
  char word[MAX_WORD_LENGTH];
  int wordLength;
  char name[MAX_WORD_LENGTH];
  char value[MAX_WORD_LENGTH];
  int isValue;
  
  in = fopen (CGI_POST_ARGUMENT, "r");
  assert(in != NULL);
  
  isValue = 0;
  wordLength = 0;
  for (;;) {
    int ch = fgetc (in);

    if (isspace (ch)) {
      /* do nothing */

    } else if (ch == '=') {
      assert(!isValue);
      word[wordLength] = 0;
      strcpy (name, word);
      wordLength = 0;
      isValue = 1;

    } else if (ch == '&') {
      assert(isValue);
      word[wordLength] = 0;
      strcpy (value, word);
      if (viewFlag) {
	SetViewVariable (name, value);
      } else {
	WriteSetCommand (out, name, value);
      }
      wordLength = 0;
      isValue = 0;

    } else if (ch == EOF) {
      assert(isValue);
      word[wordLength] = 0;
      strcpy (value, word);
      if (viewFlag) {
	SetViewVariable (name, value);
      } else {
	WriteSetCommand (out, name, value);
      }
      break;

    } else {
      if (ch == '+') {
	ch = ' ';
      }
      if (ch == '%') {
	int ch0 = fgetc (in);
	int ch1 = fgetc (in);

	ch = HexCode (ch0) * 0x10 + HexCode (ch1);
      }
      if (wordLength < MAX_WORD_LENGTH - 1) {
	word[wordLength] = ch;
	wordLength++;
      }
    }
  }

  fclose (in);
}

void WritePostCommandFile (void)
{
  char buf[1024];
  FILE *out;
  int viewFlag;

  sprintf (buf, "%s_TMP",
	   POST_COMMAND);
  out = fopen (buf, "w");
  assert(out != NULL);
  
  ViewMouseIsDragged = 0;
  viewFlag = 1;
  ReadCgiPostFile (out, viewFlag);

  WriteResizeCommand (out);

  if (ViewMouseIsDragged) {
    WriteViewCommand (out);
  } else {
    ClearAllToggles (out);
    viewFlag = 0;
    ReadCgiPostFile (out, viewFlag);
  }

  fprintf (out, "draw\n");
  fprintf (out, "value\n");

  fclose (out);

  sprintf (buf, "%s %s_TMP %s_BAK",
	   COPY,
	   POST_COMMAND, 
	   POST_COMMAND);
  system (buf);

  sprintf (buf, "%s %s_TMP %s",
	   RENAME,
	   POST_COMMAND, 
	   POST_COMMAND);
  system (buf);
  sprintf (buf, "%s %s",
	   CHANGE_MODE,
	   POST_COMMAND);
  system (buf);
}



#define MAX_CGI_POST_ARGUMENT_LENGTH 10000

int main (int argc, char *argv[]) 
{
  FILE *fp; 
  char requestMethod[256];
  char buf[1024];

  fp = fopen (INFO_COMMAND_OUTPUT, "r");
  assert(fp != NULL);
  ReadInfoCommandOutputFile (fp);

  strcpy (requestMethod, "");
  if (getenv ("REQUEST_METHOD") != NULL) {
    strcpy (requestMethod, getenv ("REQUEST_METHOD"));
  }

  if (strcmp (requestMethod, "POST") == 0) {
    {
      int length;
      char *query;

      length = atoi (getenv ("CONTENT_LENGTH"));
      assert(length <= MAX_CGI_POST_ARGUMENT_LENGTH);
      query = (char *)malloc (length + 1);
      assert(query != NULL);
      scanf ("%s", query);
      
      fp = fopen (CGI_POST_ARGUMENT, "w");
      assert(fp != NULL);
      fprintf (fp, "%s", query);
      fclose (fp);

      free (query);
      query = NULL;
    }

    WritePostCommandFile ();

#if 0
    /* NFS environment */
    for (;;) {
      system ("/bin/ls > /dev/null");
      fp = fopen (POST_COMMAND_OUTPUT, "r");
      if (fp != NULL) break;
      usleep (500000);
    }
#else
    /* stand alone environment */
    for (;;) {
      fp = fopen (POST_COMMAND_OUTPUT, "r");
      if (fp != NULL) break;
      usleep (100000);
    }
#endif

    ReadValueCommandOutputFile (fp);

    WritePanelHtmlFile ();
    
    sprintf (buf, "%s %s %s",
	     COPY,
	     POST_COMMAND_OUTPUT,
	     VALUE_COMMAND_OUTPUT);
    system (buf);

    sprintf (buf, "%s %s",
	     DELETE,
	     POST_COMMAND_OUTPUT);
    system (buf);

  } else {
    /* GET */

    fp = fopen (VALUE_COMMAND_OUTPUT, "r");
    assert(fp != NULL);
    ReadValueCommandOutputFile (fp);

    WritePanelHtmlFile ();
  }

  printf ("Content-type:text/html\n\n");
  printf ("<html>\n");
  printf ("<head><title>AutoGL.edu - test</title></head>\n");
  printf ("<frameset cols=\"250,*\">\n");
  printf ("<frame src=\"autogl_panel.cgi\" name=\"Panel\">\n");
  printf ("<frame src=\"autogl_view.cgi\" name=\"View\">\n");
  printf ("</frameset>\n");
  printf ("</html>\n");

  return 0;
}

