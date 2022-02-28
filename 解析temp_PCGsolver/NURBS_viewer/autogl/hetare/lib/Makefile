#
#         ADV_AutoGL : a simple graphics & GUI library 
#                         for making your own viewer
#           (pronounced as "O-te-ga-lu")
#



CC = gcc

# GTK+, gtkglarea and OpenGL (on Linux)
INCLUDES += `gtk-config --cflags`
LIBS += -lgtkgl -lGLU -lGL `gtk-config --libs`

CFLAGS = $(INCLUDES) -g -Wall 
#CFLAGS = $(INCLUDES) -DNDEBUG -O
#CFLAGS = $(INCLUDES) -DNDEBUG -O3

LDFLAGS = $(LIBS)

AR = ar 
ARFLAGS = crsv

.SUFFIXES: .c

.c.o:
	$(CC) -c $< $(CFLAGS)



all : libautogl.a



AUTOGL_HDRS = autogl_math.h autogl_uims.h \
	autogl_graphics.h autogl_gui.h \
	autogl_utility.h autogl.h 
AUTOGL_OBJS = autogl_math.o autogl_uims.o \
	autogl_graphics.o autogl_gui.o \
	autogl_utility.o 

autogl_math.o       : autogl_math.h
autogl_uims.o       : autogl_math.h \
	autogl_uims.h \
	autogl_uims_inside.h 
autogl_graphics.o   : autogl_math.h autogl_uims.h \
	autogl_graphics.h \
	autogl_graphics_inside.h 
autogl_gui.o           : autogl_math.h autogl_uims.h \
	autogl_graphics.h autogl_gui.h \
	autogl_uims_inside.h \
	autogl_graphics_inside.h 
autogl_utility.o    : autogl_math.h autogl_uims.h \
	autogl_graphics.h autogl_gui.h \
	autogl_utility.h

libautogl.a : $(AUTOGL_OBJS)
	$(AR) $(ARFLAGS) $@ $^



clean :
	rm -f *.o *.a
