/* autogl_gui.c : GTK+, gtkglarea and OpenGL */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include "autogl_math.h"
#include "autogl_uims.h"
#include "autogl_graphics.h"
#include "autogl_gui.h"

#include "autogl_uims_inside.h"
#include "autogl_graphics_inside.h"

#include <gtk/gtk.h>
#include <gtkgl/gtkglarea.h>



/* Module Variables */

static GtkWidget *CurrentPanelWindow;
static GtkWidget *CurrentPanel;
static GtkWidget *ViewWindow;
static GtkWidget *TopGlArea;

typedef struct {
  GSList *slist;
  int *integerValue;
} ChoiceData;

#define MAX_CHOICES 1000
static int NChoices;
ChoiceData Choices[MAX_CHOICES];



/* Label */

static void AddLabel (AutoGL_Property *property) 
{
  GtkWidget *labelWidget;
  GtkWidget *parent;

  labelWidget = gtk_label_new 
    (AutoGL_GetPropertyLabel (property));
  parent = CurrentPanel;
  gtk_box_pack_start (GTK_BOX(parent), labelWidget, FALSE, FALSE, 0);
  gtk_widget_show (labelWidget);
}



/* Button */

static void ButtonCallback (GtkWidget *widget, gpointer data) 
{
  AutoGL_Property *property = (AutoGL_Property *)data;

  AutoGL_CallCallback (property);
}

static void AddButton (AutoGL_Property *property) 
{
  GtkWidget *button;
  GtkWidget *parent;

  button = gtk_button_new_with_label 
    (AutoGL_GetPropertyLabel (property));
  parent = CurrentPanel;
  gtk_box_pack_start (GTK_BOX(parent), button, FALSE, FALSE, 0);
  gtk_signal_connect 
    (GTK_OBJECT(button), "clicked", 
     GTK_SIGNAL_FUNC(ButtonCallback), 
     (gpointer)property);
  gtk_widget_show (button);
}



/* Toggle */

static void CheckButtonCallback (GtkWidget *widget, gpointer data) 
{
  int *flag = (int*)data;

  if (GTK_TOGGLE_BUTTON(widget)->active) {
    *flag = 1;
  } else {
    *flag = 0;
  }
}

static gint RedrawCheckButtonCallback 
(GtkWidget *widget, GdkEvent *event, gpointer data) 
{
  int *flag = (int*)data;

  if (*flag) {
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON(widget), TRUE);
  } else {
    gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON(widget), FALSE);
  }
  return TRUE;
}

static void AddToggle (AutoGL_Property *property) 
{
  GtkWidget *check;
  GtkWidget *parent;

  check = gtk_check_button_new_with_label 
    (AutoGL_GetPropertyLabel (property));
  parent = CurrentPanel;
  gtk_box_pack_start (GTK_BOX(parent), check, FALSE, FALSE, 0);
  gtk_signal_connect 
    (GTK_OBJECT(check), "toggled", 
     GTK_SIGNAL_FUNC(CheckButtonCallback), 
     (gpointer)property->booleanValue);
  gtk_signal_connect 
    (GTK_OBJECT(check), "expose_event", 
     GTK_SIGNAL_FUNC(RedrawCheckButtonCallback), 
     (gpointer)property->booleanValue);
  gtk_signal_connect 
    (GTK_OBJECT(check), "map_event", 
     GTK_SIGNAL_FUNC(RedrawCheckButtonCallback), 
     (gpointer)property->booleanValue);
  gtk_widget_show (check);
}



/* Field */

static void EntryCallback (GtkWidget *widget, gpointer data) 
{
  AutoGL_Property *property = (AutoGL_Property*)data;
  static char buf[256];

  strcpy (buf, gtk_entry_get_text (GTK_ENTRY(widget)));
  AutoGL_SetPropertyString (property, buf);
}

static gint RedrawEntryCallback 
(GtkWidget *widget, GdkEvent *event, gpointer data) 
{
  AutoGL_Property *property = (AutoGL_Property*)data;
  static char buf[256];

  AutoGL_GetPropertyString (property, 
			    buf);
  gtk_entry_set_text (GTK_ENTRY(widget), buf);
  return TRUE;
}

static void AddField (AutoGL_Property *property) 
{
  GtkWidget *hbox;
  GtkWidget *parent;
  GtkWidget *labelWidget;
  GtkWidget *entry;

  hbox = gtk_hbox_new (FALSE, 0);
  parent = CurrentPanel;
  gtk_box_pack_start (GTK_BOX(parent), hbox, FALSE, FALSE, 0);
  gtk_widget_show (hbox);

  labelWidget = gtk_label_new 
    (AutoGL_GetPropertyLabel (property));
  gtk_box_pack_start (GTK_BOX(hbox), labelWidget, FALSE, FALSE, 0);
  gtk_widget_show (labelWidget);

  entry = gtk_entry_new ();
  gtk_box_pack_start (GTK_BOX(hbox), entry, FALSE, FALSE, 0);
  gtk_signal_connect 
    (GTK_OBJECT(entry), "changed", 
     GTK_SIGNAL_FUNC(EntryCallback), (gpointer)property);
  gtk_signal_connect 
    (GTK_OBJECT(entry), "expose_event", 
     GTK_SIGNAL_FUNC(RedrawEntryCallback), (gpointer)property);
  gtk_signal_connect 
    (GTK_OBJECT(entry), "map_event", 
     GTK_SIGNAL_FUNC(RedrawEntryCallback), (gpointer)property);
  gtk_widget_show (entry);
}



/* Choice */

static void RadioCallback (GtkWidget *widget, gpointer data) 
{
  ChoiceData *choice;
  int count;
  GSList *link;
  GtkWidget *radio;

  choice = (ChoiceData*)data;
  count = g_slist_length (choice->slist);
  for (link = choice->slist; link != NULL; link = link->next) {
    count--;
    radio = (GtkWidget*)link->data;
    if (GTK_TOGGLE_BUTTON(radio)->active) {
      *choice->integerValue = count;
      break;
    }
  }
}

static gint RedrawRadioCallback 
(GtkWidget *widget, GdkEvent *event, gpointer data) 
{
  ChoiceData *choice;
  int count;
  GSList *link;
  GtkWidget *radio;

  choice = (ChoiceData*)data;
  count = g_slist_length (choice->slist);
  for (link = choice->slist; link != NULL; link = link->next) {
    count--;
    radio = (GtkWidget*)link->data;
    if (radio == widget) {
      if (count == *choice->integerValue) {
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON(radio), TRUE);
      } else {
	gtk_toggle_button_set_active (GTK_TOGGLE_BUTTON(radio), FALSE);
      }
      break;
    }
  }
  return TRUE;
}

static void AddChoice (AutoGL_Property *property) 
{
  GtkWidget *parent;
  int iItem;
  GtkWidget *choiceFrame;
  GtkWidget *choiceVBox;
  ChoiceData *choice;

  assert (NChoices < MAX_CHOICES);
  choice = &Choices[NChoices];
  NChoices++;

  choice->slist = NULL;
  choice->integerValue = property->integerValue;

  choiceFrame = gtk_frame_new (AutoGL_GetPropertyLabel (property));
  parent = CurrentPanel;
  gtk_box_pack_start (GTK_BOX(parent), choiceFrame, FALSE, FALSE, 0);
  gtk_widget_show (choiceFrame);
  
  choiceVBox = gtk_vbox_new (FALSE, 0);
  gtk_container_add (GTK_CONTAINER(choiceFrame), choiceVBox);
  gtk_widget_show (choiceVBox);

  for (iItem = 0; iItem < property->nItems; iItem++) {
    GtkWidget *radio;

    radio = gtk_radio_button_new_with_label 
      (choice->slist, property->itemLabels[iItem]);
    gtk_box_pack_start (GTK_BOX(choiceVBox), radio, FALSE, FALSE, 0);
    choice->slist = gtk_radio_button_group (GTK_RADIO_BUTTON(radio));
    gtk_signal_connect 
      (GTK_OBJECT(radio), "toggled", 
       GTK_SIGNAL_FUNC(RadioCallback), (gpointer)choice);
    gtk_signal_connect 
      (GTK_OBJECT(radio), "expose_event", 
       GTK_SIGNAL_FUNC(RedrawRadioCallback), (gpointer)choice);
    gtk_signal_connect 
      (GTK_OBJECT(radio), "map_event", 
       GTK_SIGNAL_FUNC(RedrawRadioCallback), (gpointer)choice);
    gtk_widget_show (radio);
  }
}



static void MakePanel (void)
{
  int nPropertys;
  int iProperty;

  nPropertys = AutoGL_GetNPropertys ();
  if (nPropertys == 0) return;

  CurrentPanelWindow = gtk_window_new (GTK_WINDOW_TOPLEVEL);
  gtk_window_set_title (GTK_WINDOW(CurrentPanelWindow), "Panel");

  CurrentPanel = gtk_vbox_new (FALSE, 0);
  gtk_container_add (GTK_CONTAINER(CurrentPanelWindow), CurrentPanel);
  gtk_widget_show (CurrentPanel);

  for (iProperty = 0; iProperty < nPropertys; iProperty++) {
    AutoGL_Property *property = AutoGL_GetProperty (iProperty);
    
    switch (property->type) {
    case AUTOGL_PROPERTY_COMMENT:
      AddLabel (property);
      break;
    case AUTOGL_PROPERTY_CALLBACK:
      AddButton (property);
      break;
    case AUTOGL_PROPERTY_BOOLEAN:
      AddToggle (property);
      break;
    case AUTOGL_PROPERTY_STRING:
    case AUTOGL_PROPERTY_INTEGER:
    case AUTOGL_PROPERTY_REAL:
      if (0 < property->nItems) {
	AddChoice (property);
      } else {
	AddField (property);
      }
      break;
    default:
      assert(0);
      break;
    }
  }

  gtk_widget_show (CurrentPanelWindow);
}



/* View */

static gint RealizeGlAreaCallback
(GtkWidget *widget)
{
  GdkFont* font;

  if (gtk_gl_area_make_current (GTK_GL_AREA(widget))) {
    AutoGL_RealizeOpenGLWindow ();

    /* make display list of character font */
    font = gdk_font_load ("a14");
    if (!font) {
      g_print ("Can`t load font a14 \n");
      exit (1);
    }
    gdk_gl_use_gdk_font (font, 0, 128, AutoGL_GetOpenGLCharBase ());
    gdk_font_unref (font);
  }

  return TRUE;
}

static gint ResizeGlAreaCallback
(GtkWidget *widget, GdkEventConfigure *event, gpointer data)
{
  if (gtk_gl_area_make_current (GTK_GL_AREA(widget))) {
    AutoGL_SetViewRangeDc 
      (widget->allocation.width, widget->allocation.height);
  }

  return TRUE;
}

static gint RedrawGlAreaCallback 
(GtkWidget *widget, GdkEventExpose *event, gpointer data)
{
  if (event->count > 0) {
    return TRUE;
  }

  if (gtk_gl_area_make_current(GTK_GL_AREA(widget))) {
    AutoGL_StartDrawingOpenGLWindow ();
    AutoGL_CallViewRedrawCallback ();
    AutoGL_EndDrawingOpenGLWindow ();
    gtk_gl_area_swapbuffers (GTK_GL_AREA(widget));
  }

  return TRUE;
}

static void GetDcOfMouse 
(int *dcX_OUT, int *dcY_OUT, int *dcZ_OUT,
 int mouseX, int mouseY)
{
  int widthDc, heightDc;

  AutoGL_GetViewRangeDc (&widthDc, &heightDc);
  *dcX_OUT = mouseX - widthDc / 2;
  *dcY_OUT = heightDc / 2 - mouseY;
  *dcZ_OUT = 0;
}

static gint ButtonPressGlAreaCallback 
(GtkWidget *widget, GdkEventButton *event, gpointer data)
{
  if (event->button == 1) {
    int dcX, dcY, dcZ;

    GetDcOfMouse (&dcX, &dcY, &dcZ,
		  (int)event->x, (int)event->y);
    AutoGL_SetPointingDeviceEvent 
      (AUTOGL_EVENT_POINTING_DEVICE_PRESS, 
       dcX, dcY, dcZ);
    AutoGL_CallViewEventCallback ();
  }

  return TRUE;
}

static gint MotionNotifyGlAreaCallback 
(GtkWidget *widget, GdkEventMotion *event, gpointer data)
{
  int x, y;
  GdkModifierType state;

  if (event->is_hint) {
    gdk_window_get_pointer (event->window, 
			    &x, &y, &state);
  } else {
    x = event->x;
    y = event->y;
    state = event->state;
  }

  if (event->state & GDK_BUTTON1_MASK) {
    int dcX, dcY, dcZ;

    GetDcOfMouse (&dcX, &dcY, &dcZ,
		  x, y);
    AutoGL_SetPointingDeviceEvent 
      (AUTOGL_EVENT_POINTING_DEVICE_DRAG, 
       dcX, dcY, dcZ);
    AutoGL_CallViewEventCallback ();
  }
  
  return TRUE;
}

static gint ButtonReleaseGlAreaCallback 
(GtkWidget *widget, GdkEventButton *event, gpointer data)
{
  if (event->button == 1) {
    int dcX, dcY, dcZ;

    GetDcOfMouse (&dcX, &dcY, &dcZ,
		  (int)event->x, (int)event->y);
    AutoGL_SetPointingDeviceEvent 
      (AUTOGL_EVENT_POINTING_DEVICE_RELEASE, 
       dcX, dcY, dcZ);
    AutoGL_CallViewEventCallback ();
  }

  return TRUE;
}

static void MakeView (void)
{
  int attributes[] = {
    GDK_GL_RGBA,
    GDK_GL_DOUBLEBUFFER,
    GDK_GL_DEPTH_SIZE, 1,
    GDK_GL_NONE
  };
  int sizeDc;

  if (AutoGL_GetViewRedrawCallback () == NULL) return;

  /* create a gtkgl GL area widget */
  
  /* Check if OpenGL is supported. */
  if (gdk_gl_query () == FALSE) {
    g_print ("OpenGL not supported\n");
    exit (1);
  }
  
  ViewWindow = gtk_window_new (GTK_WINDOW_TOPLEVEL);

  gtk_window_set_title (GTK_WINDOW(ViewWindow), "View");

  TopGlArea = GTK_WIDGET(gtk_gl_area_new (attributes));
  gtk_container_add (GTK_CONTAINER(ViewWindow), TopGlArea);

  sizeDc = 500;

  gtk_drawing_area_size (GTK_DRAWING_AREA(TopGlArea), 
			 sizeDc, sizeDc);
  gtk_widget_set_events (GTK_WIDGET(TopGlArea),
			 GDK_EXPOSURE_MASK |
			 GDK_BUTTON_PRESS_MASK |
			 GDK_POINTER_MOTION_MASK |
			 GDK_POINTER_MOTION_HINT_MASK |
			 GDK_BUTTON_RELEASE_MASK);
  gtk_signal_connect
    (GTK_OBJECT(TopGlArea), "realize",
     GTK_SIGNAL_FUNC(RealizeGlAreaCallback), NULL);
  gtk_signal_connect 
    (GTK_OBJECT(TopGlArea), "configure_event", 
     GTK_SIGNAL_FUNC(ResizeGlAreaCallback), NULL);

  gtk_signal_connect 
    (GTK_OBJECT(TopGlArea), "map_event", 
     GTK_SIGNAL_FUNC(RedrawGlAreaCallback), NULL);
  gtk_signal_connect 
    (GTK_OBJECT(TopGlArea), "expose_event", 
     GTK_SIGNAL_FUNC(RedrawGlAreaCallback), NULL);

  gtk_signal_connect 
    (GTK_OBJECT(TopGlArea), "button_press_event",
     GTK_SIGNAL_FUNC(ButtonPressGlAreaCallback), NULL);
  gtk_signal_connect 
    (GTK_OBJECT(TopGlArea), "motion_notify_event",
     GTK_SIGNAL_FUNC(MotionNotifyGlAreaCallback), NULL);
  gtk_signal_connect 
    (GTK_OBJECT(TopGlArea), "button_release_event",
     GTK_SIGNAL_FUNC(ButtonReleaseGlAreaCallback), NULL);

  gtk_widget_show (TopGlArea);
  gtk_widget_show (ViewWindow);
}



void AutoGL_DrawView (void)
{
  GdkRectangle updateRectangle;

  updateRectangle.x = 0;
  updateRectangle.y = 0;
  updateRectangle.width = TopGlArea->allocation.width;
  updateRectangle.height = TopGlArea->allocation.height;
  gtk_widget_draw (TopGlArea, &updateRectangle);
}



int main (int argc, char *argv[]) 
{
  /* call the user-side set-up routine */
  AutoGL_SetUp (argc, argv);

  gtk_init (&argc, &argv);

  MakePanel ();
  MakeView ();
  gtk_main ();

  /* never return */
  return 0;
}

