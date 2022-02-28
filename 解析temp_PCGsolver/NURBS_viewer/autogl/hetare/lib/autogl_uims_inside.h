/* AutoGL */

#ifndef __AUTOGL_UIMS_INSIDE__
#define __AUTOGL_UIMS_INSIDE__



#define AUTOGL_MAX_PROPERTY_ITEMS 10
#define AUTOGL_PROPERTY_NAME_LENGTH 80
#define AUTOGL_PROPERTY_LABEL_LENGTH 80

typedef enum {
  AUTOGL_PROPERTY_COMMENT, 
  AUTOGL_PROPERTY_CALLBACK, 
  AUTOGL_PROPERTY_BOOLEAN, 
  AUTOGL_PROPERTY_STRING, 
  AUTOGL_PROPERTY_INTEGER,
  AUTOGL_PROPERTY_REAL
} AutoGL_PropertyType;

typedef struct {
  AutoGL_PropertyType type;
  char name[AUTOGL_PROPERTY_NAME_LENGTH];
  char label[AUTOGL_PROPERTY_LABEL_LENGTH];
  AutoGL_CallbackType callback;
  int *booleanValue;
  char *stringValue;
  int stringLength;
  int *integerValue;
  int minIntegerValue;
  int maxIntegerValue;
  int nItems;
  char itemLabels[AUTOGL_MAX_PROPERTY_ITEMS][AUTOGL_PROPERTY_LABEL_LENGTH];
  double *realValue;
  double minRealValue;
  double maxRealValue;
} AutoGL_Property;

AutoGL_CallbackType AutoGL_GetViewRedrawCallback (void);

AutoGL_CallbackType AutoGL_GetViewEventCallback (void);

int AutoGL_GetNPropertys (void);

AutoGL_Property *AutoGL_GetProperty (int propertyId);

char *AutoGL_GetPropertyLabel 
(AutoGL_Property *property);

void AutoGL_SetPropertyString 
(AutoGL_Property *property, 
 const char buf[]);

void AutoGL_GetPropertyString 
(AutoGL_Property *property, 
 char buf_OUT[256]);

void AutoGL_CallCallback 
(AutoGL_Property *property);

void AutoGL_CallViewRedrawCallback (void);

void AutoGL_CallViewEventCallback (void);



#endif  /* __AUTOGL_UIMS_INSIDE__ */







