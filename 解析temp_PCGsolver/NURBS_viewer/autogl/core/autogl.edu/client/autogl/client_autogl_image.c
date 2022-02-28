/* client_autogl_image.c : client side AutoGL implementation */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/* usleep */
#include <unistd.h>

#include <autogl.h>



/* forward declaration */

void ParseInfoCommandOutputFile (void);
void ParseValueCommandOutputFile (void);
void WriteSetCommandForEachWidget (FILE *fp);



/* basic */

void Error (char *message) 
{
  fprintf (stderr, "ERROR\n");
  fprintf (stderr, message);
  fprintf (stderr, "\n");
}



/* reader */

#define MAX_WORD_CHARS 256

FILE *_fp;
int _currentChar;
char _wordChars[MAX_WORD_CHARS];

void InitializeReader (FILE *fp) 
{
  _fp = fp;
  _currentChar = fgetc (_fp);
}

void TerminateReader (void) 
{
  fclose (_fp);
  _fp = NULL;
}

int IsEOF (void) 
{
  return _currentChar == EOF;
}

int IsWhiteSpace (void) 
{
  return _currentChar == ' ' 
    || _currentChar == '\n'
    || _currentChar == '\t';
}

char *ReadWord (void) 
{
  int count;

  while (IsWhiteSpace ()) {
    _currentChar = fgetc (_fp);
  }
  if (IsEOF ()) return NULL;
  
  count = 0;
  if (_currentChar == '"') {
    _currentChar = fgetc (_fp);
    for (;;) {
      if (IsEOF ()) {
	Error ("must be \"");
	break;
      }
      if (_currentChar == '"') {
	_currentChar = fgetc (_fp);
	break;
      }
      
      if (MAX_WORD_CHARS <= count) {
	Error ("word is too long");
	break;
      }
      _wordChars[count++] = (char)_currentChar;
      _currentChar = fgetc (_fp);
    }
  } else {
    for (;;) {
      if (IsWhiteSpace ()) break;
      if (IsEOF ()) break;
      
      if (MAX_WORD_CHARS <= count) {
	Error ("word is too long");
	break;
      }
      _wordChars[count++] = (char)_currentChar;
      _currentChar = fgetc (_fp);
    }
  }
  _wordChars[count++] = 0;

  return _wordChars;
}
    


/* application model I/O */

#define MAX_COMMAND_LENGTH 1024
#define MAX_FILE_NAME_LENGTH 80
#define MAX_SERVER_HOST_NAME_LENGTH 80
#define MAX_SERVER_DIRECTORY_NAME_LENGTH 80

char ServerHostName[MAX_SERVER_HOST_NAME_LENGTH] 
= "who";
char ServerDirectoryName[MAX_SERVER_DIRECTORY_NAME_LENGTH] 
= "test_autogl_client";

#define VIEW_IMAGE_PPM "result.ppm"
#define INFO_COMMAND_OUTPUT "autogl_info_command_output_file"
#define VALUE_COMMAND_OUTPUT "autogl_value_command_output_file"
#define POST_COMMAND "autogl_post_command_file"
#define POST_COMMAND_OUTPUT "autogl_post_command_output_file"

#define VIEW_IMAGE_WIDTH 500
#define VIEW_IMAGE_HEIGHT 500

int _imageWidth, _imageHeight;
unsigned char *_imageData;

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

void PutFile (char *fileName) 
{
  char buf[MAX_COMMAND_LENGTH];

  assert(strlen (fileName) < MAX_FILE_NAME_LENGTH);
  assert(strlen (ServerHostName) < MAX_SERVER_HOST_NAME_LENGTH);
  assert(strlen (ServerDirectoryName) < MAX_SERVER_DIRECTORY_NAME_LENGTH);

  sprintf (buf, "%s %s %s:%s/%s_TMP", 
	   REMOTE_COPY,
	   fileName, 
	   ServerHostName, ServerDirectoryName, fileName);
  system (buf);

  fprintf (stderr, "%s\n", buf);

  sprintf (buf, "%s %s %s:%s/%s_BAK", 
	   REMOTE_COPY,
	   fileName, 
	   ServerHostName, ServerDirectoryName, fileName);
  system (buf);

  fprintf (stderr, "%s\n", buf);

  sprintf (buf, "%s %s %s %s/%s_TMP %s/%s", 
	   REMOTE_EXECUTE,
	   ServerHostName, 
	   RENAME,
	   ServerDirectoryName, fileName, 
	   ServerDirectoryName, fileName);
  system (buf);

  fprintf (stderr, "%s\n", buf);

  sleep (3);
}

void ReadImage (void) 
{
  if (_imageData != NULL) {
    free (_imageData);
    _imageData = NULL;
  }

  GetFile (VIEW_IMAGE_PPM);

  AutoGL_ReadImageFromPPMFile 
    (VIEW_IMAGE_PPM,
     &_imageWidth, &_imageHeight, &_imageData);
}

void ReadModelLayout (void) 
{
  FILE *fp;

  GetFile (INFO_COMMAND_OUTPUT);

  fp = fopen (INFO_COMMAND_OUTPUT, "r");
  assert(fp);
  InitializeReader (fp);
  ParseInfoCommandOutputFile ();
  TerminateReader ();
}

void ReadModelValue (void) 
{
  FILE *fp;

  GetFile (VALUE_COMMAND_OUTPUT);

  fp = fopen (VALUE_COMMAND_OUTPUT, "r");
  assert(fp);
  InitializeReader (fp);
  ParseValueCommandOutputFile ();
  TerminateReader ();
}

void WriteResizeCommand (FILE *fp) 
{
  fprintf (fp, "resize %d %d\n",
	   VIEW_IMAGE_WIDTH, VIEW_IMAGE_HEIGHT);
}

void WriteSetCommand 
(FILE *fp, char *name, char *value) 
{
  fprintf (fp, "set %s %s\n",
	   name, value);
}

void WriteSetCommandWithQuotation  
(FILE *fp, char *name, char *value) 
{
  fprintf (fp, "set %s \"%s\"\n",
	   name, value);
}

void WriteCallCommand 
(FILE *fp, char *callbackName) 
{
  fprintf (fp, "call %s\n",
	   callbackName);
}

void WriteDrawCommand (FILE *fp) 
{
  fprintf (fp, "draw\n");
}

void WriteDragCommand 
(FILE *fp, int startX, int startY, int endX, int endY) 
{
  fprintf (fp, "drag %d %d %d %d\n",
	   startX, startY, endX, endY);
}

void WriteValueCommand (FILE *fp) 
{
  fprintf (fp, "value\n");
}

void WriteWriteCommand (FILE *fp) 
{
  fprintf (fp, "write\n");
}

void PostCallCommand (char *callbackName) 
{
  FILE *fp;

  fp = fopen (POST_COMMAND, "w");
  assert(fp);
  WriteResizeCommand (fp);
  WriteSetCommandForEachWidget (fp);
  WriteCallCommand (fp, callbackName);
  WriteValueCommand (fp);
  fclose (fp);
  
  PutFile (POST_COMMAND);

  GetFile (POST_COMMAND_OUTPUT);

  fp = fopen (POST_COMMAND_OUTPUT, "r");
  assert(fp);
  InitializeReader (fp);
  ParseValueCommandOutputFile ();
  TerminateReader ();
}

void PostDrawCommand (void) 
{
  FILE *fp;
  
  fp = fopen (POST_COMMAND, "w");
  assert(fp);
  WriteResizeCommand (fp);
  WriteSetCommandForEachWidget (fp);
  WriteDrawCommand (fp);
  WriteValueCommand (fp);
  fclose (fp);
  
  PutFile (POST_COMMAND);

  GetFile (POST_COMMAND_OUTPUT);

  fp = fopen (POST_COMMAND_OUTPUT, "r");
  assert(fp);
  InitializeReader (fp);
  ParseValueCommandOutputFile ();
  TerminateReader ();
}

void PostDragCommand 
(int startX, int startY, int endX, int endY) 
{
  FILE *fp;
  
  fp = fopen (POST_COMMAND, "w");
  assert(fp);
  WriteResizeCommand (fp);
  WriteSetCommandForEachWidget (fp);
  WriteDragCommand (fp, startX, startY, endX, endY);
  WriteDrawCommand (fp);
  WriteValueCommand (fp);
  fclose (fp);
  
  PutFile (POST_COMMAND);

  GetFile (POST_COMMAND_OUTPUT);

  fp = fopen (POST_COMMAND_OUTPUT, "r");
  assert(fp);
  InitializeReader (fp);
  ParseValueCommandOutputFile ();
  TerminateReader ();
}

void PostWriteCommand (void) 
{
  FILE *fp;

  fp = fopen (POST_COMMAND, "w");
  assert(fp);
  WriteWriteCommand (fp);
  fclose (fp);
  
  PutFile (POST_COMMAND);
}



/* widget */

#define MAX_KEYWORD_LENGTH 80
#define MAX_NAME_LENGTH 80
#define MAX_LABEL_LENGTH 80
#define MAX_STRING_LENGTH MAX_WORD_CHARS

#define MAX_WIDGETS 100
#define MAX_LABEL_COMPONENTS 100
#define MAX_COMPONENTS 50
#define MAX_CHOICE_ITEMS 10

#define PROPERTY_GROUP 0
#define PROPERTY_COMMENT 1
#define PROPERTY_CALLBACK 2
#define PROPERTY_BOOLEAN 3
#define PROPERTY_STRING 4
#define PROPERTY_INTEGER 5
#define PROPERTY_REAL 6

#define WIDGET_LABEL 0
#define WIDGET_BUTTON 1
#define WIDGET_TOGGLE 2
#define WIDGET_FIELD 3
#define WIDGET_CHOICE 4

int _nWidgets;
char _widgetNames[MAX_WIDGETS][MAX_NAME_LENGTH];
int _widgetTypes[MAX_WIDGETS];
int _widgetComponentIndexs[MAX_WIDGETS];
int _widgetPropertyTypes[MAX_WIDGETS];

int _nButtons;
int _buttonWidgetIds[MAX_COMPONENTS];

int _nToggles;
int _toggleBooleanVariables[MAX_COMPONENTS];

int _nFields;
char _fieldStringVariables[MAX_COMPONENTS][MAX_STRING_LENGTH];

int _nChoices;
int _choiceIntegerVariables[MAX_COMPONENTS];
char _nChoiceItems[MAX_COMPONENTS];
char _choiceItemLabels[MAX_COMPONENTS][MAX_CHOICE_ITEMS][MAX_LABEL_LENGTH];

int WidgetIdOfName (char *name) 
{
  int iWidget;

  for (iWidget = 0; iWidget < _nWidgets; iWidget++) {
    if (strcmp (_widgetNames[iWidget], name) == 0) {
      return iWidget;
    }
  }
  return -1;
}



/* label widget */

void MakeLabelWidget (char *name, char *label, int propertyType) 
{
  strcpy (_widgetNames[_nWidgets], name);
  _widgetTypes[_nWidgets] = WIDGET_LABEL;
  _widgetComponentIndexs[_nWidgets] = -1;
  _widgetPropertyTypes[_nWidgets] = propertyType;
  _nWidgets++;
  
  AutoGL_AddComment ();
  AutoGL_SetLabel (label);
}



/* button widget */

void ButtonCallback0 (void) 
{
  int buttonId = 0;
  int widgetId = _buttonWidgetIds[buttonId];
  
  PostCallCommand (_widgetNames[widgetId]);
  ReadImage ();
  AutoGL_DrawView ();
}

void ButtonCallback1 (void) 
{
  int buttonId = 1;
  int widgetId = _buttonWidgetIds[buttonId];
  
  PostCallCommand (_widgetNames[widgetId]);
  ReadImage ();
  AutoGL_DrawView ();
}

void ButtonCallback2 (void) 
{
  int buttonId = 2;
  int widgetId = _buttonWidgetIds[buttonId];
  
  PostCallCommand (_widgetNames[widgetId]);
  ReadImage ();
  AutoGL_DrawView ();
}

void ButtonCallback3 (void) 
{
  int buttonId = 3;
  int widgetId = _buttonWidgetIds[buttonId];
  
  PostCallCommand (_widgetNames[widgetId]);
  ReadImage ();
  AutoGL_DrawView ();
}

void ButtonCallback4 (void) 
{
  int buttonId = 4;
  int widgetId = _buttonWidgetIds[buttonId];
  
  PostCallCommand (_widgetNames[widgetId]);
  ReadImage ();
  AutoGL_DrawView ();
}

void ButtonCallback5 (void) 
{
  int buttonId = 5;
  int widgetId = _buttonWidgetIds[buttonId];
  
  PostCallCommand (_widgetNames[widgetId]);
  ReadImage ();
  AutoGL_DrawView ();
}

void ButtonCallback6 (void) 
{
  int buttonId = 6;
  int widgetId = _buttonWidgetIds[buttonId];
  
  PostCallCommand (_widgetNames[widgetId]);
  ReadImage ();
  AutoGL_DrawView ();
}

void ButtonCallback7 (void) 
{
  int buttonId = 7;
  int widgetId = _buttonWidgetIds[buttonId];
  
  PostCallCommand (_widgetNames[widgetId]);
  ReadImage ();
  AutoGL_DrawView ();
}

void ButtonCallback8 (void) 
{
  int buttonId = 8;
  int widgetId = _buttonWidgetIds[buttonId];
  
  PostCallCommand (_widgetNames[widgetId]);
  ReadImage ();
  AutoGL_DrawView ();
}

void ButtonCallback9 (void) 
{
  int buttonId = 9;
  int widgetId = _buttonWidgetIds[buttonId];
  
  PostCallCommand (_widgetNames[widgetId]);
  ReadImage ();
  AutoGL_DrawView ();
}

void MakeButtonWidget (char *name, char *label) 
{
  int widgetId = _nWidgets;
  int componentIndex = _nButtons;
  
  strcpy (_widgetNames[_nWidgets], name);
  _widgetTypes[_nWidgets] = WIDGET_BUTTON;
  _widgetComponentIndexs[_nWidgets] = componentIndex;
  _widgetPropertyTypes[_nWidgets] = PROPERTY_CALLBACK;
  _nWidgets++;
  
  switch (componentIndex) {
  case 0:  AutoGL_AddCallback (ButtonCallback0, name);   break;
  case 1:  AutoGL_AddCallback (ButtonCallback1, name);   break;
  case 2:  AutoGL_AddCallback (ButtonCallback2, name);   break;
  case 3:  AutoGL_AddCallback (ButtonCallback3, name);   break;
  case 4:  AutoGL_AddCallback (ButtonCallback4, name);   break;
  case 5:  AutoGL_AddCallback (ButtonCallback5, name);   break;
  case 6:  AutoGL_AddCallback (ButtonCallback6, name);   break;
  case 7:  AutoGL_AddCallback (ButtonCallback7, name);   break;
  case 8:  AutoGL_AddCallback (ButtonCallback8, name);   break;
  case 9:  AutoGL_AddCallback (ButtonCallback9, name);   break;
  default:
    assert(0);
    break;
  }
  if (0 < strlen (label)) {
    AutoGL_SetLabel (label); 
  }

  _buttonWidgetIds[_nButtons] = widgetId;
  _nButtons++;
}



/* toggle widget */

void MakeToggleWidget (char *name, char *label) 
{
  int componentIndex = _nToggles;
  
  strcpy (_widgetNames[_nWidgets], name);
  _widgetTypes[_nWidgets] = WIDGET_TOGGLE;
  _widgetComponentIndexs[_nWidgets] = componentIndex;
  _widgetPropertyTypes[_nWidgets] = PROPERTY_BOOLEAN;
  _nWidgets++;
  
  AutoGL_AddBoolean (&_toggleBooleanVariables[_nToggles], name); 
  if (0 < strlen (label)) {
    AutoGL_SetLabel (label); 
  }
  
  _nToggles++;
}



/* field widget */

#define MAX_FIELD_CHARS 16

void MakeFieldWidget 
(char *name, char *label, int propertyType) 
{
  int componentIndex = _nFields;

  strcpy (_widgetNames[_nWidgets], name);
  _widgetTypes[_nWidgets] = WIDGET_FIELD;
  _widgetComponentIndexs[_nWidgets] = componentIndex;
  _widgetPropertyTypes[_nWidgets] = propertyType;
  _nWidgets++;
  
  AutoGL_AddString (_fieldStringVariables[_nFields], name, 
		    MAX_STRING_LENGTH);
  if (0 < strlen (label)) {
    AutoGL_SetLabel (label); 
  }

  _nFields++;
}



/* choice widget */

void MakeChoiceWidget (char *name, char *label) 
{
  int componentIndex = _nChoices;
  
  strcpy (_widgetNames[_nWidgets], name);
  _widgetTypes[_nWidgets] = WIDGET_CHOICE;
  _widgetComponentIndexs[_nWidgets] = componentIndex;
  _widgetPropertyTypes[_nWidgets] = PROPERTY_INTEGER;
  _nWidgets++;
  
  AutoGL_AddInteger (&_choiceIntegerVariables[_nChoices], name);
  if (0 < strlen (label)) {
    AutoGL_SetLabel (label); 
  }

  _nChoices++;
}

void AddItemIntoChoiceWidget (char *itemLabel) 
{
  int choiceId = _nChoices - 1;
  int nItems = _nChoiceItems[choiceId];

  strcpy (_choiceItemLabels[choiceId][nItems], itemLabel);
  _nChoiceItems[choiceId]++;

  AutoGL_AddIntegerItem (itemLabel);
}



void ParseInfoCommandOutputFile (void) 
{
  for (;;) {
    char *word;
    char keyword[MAX_KEYWORD_LENGTH];
    char nameKeyword[MAX_KEYWORD_LENGTH];
    char name[MAX_NAME_LENGTH];
    char labelKeyword[MAX_KEYWORD_LENGTH];
    char label[MAX_LABEL_LENGTH];

    word = ReadWord ();
    if (word == NULL) break;
    strcpy (keyword, word);
    
    strcpy (nameKeyword, ReadWord ());
    if (strcmp (nameKeyword, "name") != 0) {
      Error ("must be name");
      return;
    }
    strcpy (name, ReadWord ());
    
    strcpy (labelKeyword, ReadWord ());
    if (strcmp (labelKeyword, "label") != 0) {
      Error ("must be label");
      return;
    }
    strcpy (label, ReadWord ());
    
    if (strcmp (keyword, "group") == 0) {
      MakeLabelWidget (name, label, PROPERTY_GROUP);
      
    } else if (strcmp (keyword, "comment") == 0) {
      MakeLabelWidget (name, label, PROPERTY_COMMENT);
      
    } else if (strcmp (keyword, "callback") == 0) {
      MakeButtonWidget (name, label);
      
    } else if (strcmp (keyword, "boolean") == 0) {
      MakeToggleWidget (name, label);
      
    } else if (strcmp (keyword, "string") == 0) {
      MakeFieldWidget (name, label, PROPERTY_STRING);
      
    } else if (strcmp (keyword, "integer") == 0) {
      char nItemsKeyword[MAX_KEYWORD_LENGTH];
      char nItemsStr[MAX_KEYWORD_LENGTH];
      int nItems;
      int iItem;

      strcpy (nItemsKeyword, ReadWord ());
      if (strcmp (nItemsKeyword, "nItems") != 0) {
	Error ("must be nItems");
	return;
      }
      strcpy (nItemsStr, ReadWord ());

      nItems = atoi (nItemsStr);
      if (nItems == 0) {
	MakeFieldWidget (name, label, PROPERTY_INTEGER);
      } else {
	MakeChoiceWidget (name, label);
      }
      for (iItem = 0; iItem < nItems; iItem++) {
	char itemLabelKeyword[MAX_KEYWORD_LENGTH];
	char itemLabel[MAX_LABEL_LENGTH];

	strcpy (itemLabelKeyword, ReadWord ());
	if (strcmp (itemLabelKeyword, "itemLabel") != 0) {
	  Error ("must be itemLabel");
	  return;
	}
	strcpy (itemLabel, ReadWord ());
	
	AddItemIntoChoiceWidget (itemLabel);
      }
      
    } else if (strcmp (keyword, "real") == 0) {
      MakeFieldWidget (name, label, PROPERTY_REAL);
      
    } else {
      Error ("unsupported widget type");
      return;
    }
  }
}

void ParseValueCommandOutputFile (void) 
{
  for (;;) {
    char *word;
    char keyword[MAX_KEYWORD_LENGTH];
    char nameKeyword[MAX_KEYWORD_LENGTH];
    char name[MAX_NAME_LENGTH];
    char valueKeyword[MAX_KEYWORD_LENGTH];
    char value[MAX_STRING_LENGTH];

    word = ReadWord ();
    if (word == NULL) break;
    strcpy (keyword, word);
    
    strcpy (nameKeyword, ReadWord ());
    if (strcmp (nameKeyword, "name") != 0) {
      Error ("must be name");
      return;
    }
    strcpy (name, ReadWord ());
    
    if (strcmp (keyword, "group") == 0
	|| strcmp (keyword, "comment") == 0
	|| strcmp (keyword, "callback") == 0) continue;
    
    strcpy (valueKeyword, ReadWord ());
    if (strcmp (valueKeyword, "value") != 0) {
      Error ("must be value");
      return;
    }
    strcpy (value, ReadWord ());

    if (strcmp (keyword, "boolean") == 0) {
      int widgetId = WidgetIdOfName (name);
      int toggleId = _widgetComponentIndexs[widgetId];

      if (strcmp (value, "true") == 0) {
	_toggleBooleanVariables[toggleId] = 1;
      } else if (strcmp (value, "false") == 0) {
	_toggleBooleanVariables[toggleId] = 0;
      } else {
	Error ("value must be true or false");
	return;
      }
      
    } else if (strcmp (keyword, "string") == 0) {
      int widgetId = WidgetIdOfName (name);
      int fieldId = _widgetComponentIndexs[widgetId];

      strcpy (_fieldStringVariables[fieldId], value);
      
    } else if (strcmp (keyword, "integer") == 0) {
      int widgetId = WidgetIdOfName (name);

      if (_widgetTypes[widgetId] == WIDGET_FIELD) {
	int fieldId = _widgetComponentIndexs[widgetId];

	strcpy (_fieldStringVariables[fieldId], value);
	
      } else if (_widgetTypes[widgetId] == WIDGET_CHOICE) {
	int choiceId = _widgetComponentIndexs[widgetId];
	int integerValue = atoi (value);
	char itemLabel[MAX_LABEL_LENGTH];

	_choiceIntegerVariables[choiceId] = integerValue;
	strcpy (itemLabel, ReadWord ());
	
      } else {
	Error ("integer must be field or choice");
	return;
      }
      
    } else if (strcmp (keyword, "real") == 0) {
      int widgetId = WidgetIdOfName (name);
      int fieldId = _widgetComponentIndexs[widgetId];

      strcpy (_fieldStringVariables[fieldId], value);
      
    } else {
      Error ("unsupported widget type");
      return;
    }
  }
}

void WriteSetCommandForEachWidget (FILE *fp) 
{
  int iWidget;

  for (iWidget = 0; iWidget < _nWidgets; iWidget++) {
    char *name = _widgetNames[iWidget];

    switch (_widgetTypes[iWidget]) {
    case WIDGET_LABEL:
    case WIDGET_BUTTON:
      break;
    case WIDGET_TOGGLE:
      {
	int toggleId = _widgetComponentIndexs[iWidget];

	if (_toggleBooleanVariables[toggleId]) {
	  WriteSetCommand (fp, name, "true");
	} else {
	  WriteSetCommand (fp, name, "false");
	}
      }
      break;
    case WIDGET_FIELD:
      {
	int fieldId = _widgetComponentIndexs[iWidget];
	char *str = _fieldStringVariables[fieldId];

	switch (_widgetPropertyTypes[iWidget]) {
	case PROPERTY_STRING:
	  WriteSetCommandWithQuotation (fp, name, str);
	  break;
	case PROPERTY_INTEGER:
	case PROPERTY_REAL:
	  WriteSetCommand (fp, name, str);
	  break;
	default:
	  Error ("illegal property type of field");
	  break;
	}
      }
      break;
    case WIDGET_CHOICE:
      {
	int choiceId = _widgetComponentIndexs[iWidget];
	int value = _choiceIntegerVariables[choiceId];

	assert(0 <= value);
	assert(value < _nChoiceItems[choiceId]);
	WriteSetCommand (fp, name, 
			 _choiceItemLabels[choiceId][value]);
      }
      break;
    default:
      Error ("illegal widget type");
      break;
    }
  }
}
    


static void RedrawView (void) 
{
  if (_imageData != NULL) {
    AutoGL_DrawImage (-100, -100, 0,
		      _imageWidth, _imageHeight, _imageData);
  }
}

static void HandleEvent (void)
{
  static int startX, startY;
  int event = AutoGL_GetViewEvent ();

  switch (event) {
  case AUTOGL_EVENT_POINTING_DEVICE_PRESS:
    {
      int dcX, dcY, dcZ;
      
      AutoGL_GetPointingDevicePositionDc (&dcX, &dcY, &dcZ);
      startX = dcX + VIEW_IMAGE_HEIGHT / 2;
      startY = VIEW_IMAGE_HEIGHT / 2 - dcY;
    }
    break;
  case AUTOGL_EVENT_POINTING_DEVICE_DRAG:
    break;
  case AUTOGL_EVENT_POINTING_DEVICE_RELEASE:
    {
      int dcX, dcY, dcZ;
      int endX, endY;
      
      AutoGL_GetPointingDevicePositionDc (&dcX, &dcY, &dcZ);
      endX = dcX + VIEW_IMAGE_HEIGHT / 2;
      endY = VIEW_IMAGE_HEIGHT / 2 - dcY;
      PostDragCommand (startX, startY, 
		       endX, endY);
      ReadImage ();
      AutoGL_DrawView ();
    }
    break;
  default:
    break;
  }
}

static void UpdateMetafile (void) 
{
  PostWriteCommand () ;
  AutoGL_DrawView ();  
}

void AutoGL_SetUp (int argc, char *argv[]) 
{
  AutoGL_SetViewSize (100);

  AutoGL_SetViewRedrawCallback (RedrawView);
  AutoGL_SetViewEventCallback (HandleEvent);

  AutoGL_AddGroup ("AutoGL:Panel");    
  ReadModelLayout ();
  ReadModelValue ();
  AutoGL_AddCallback (UpdateMetafile, "UpdateMetafile");

  PostDrawCommand ();
  ReadImage ();
}

