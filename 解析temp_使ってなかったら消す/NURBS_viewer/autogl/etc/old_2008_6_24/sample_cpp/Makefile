#
#         ADV_AutoGL : a simple graphics & GUI library 
#                         for making your own viewer
#           (pronounced as "O-te-ga-lu")
#
# This "Makefile" is for GNU make only.
# You need to remove commands specific to GNU make, such as "ifeq" and "endif", 
# if you want to use any other kind of Unix-standard make.
#
#
# < Platform Selection >
#    
# The following platforms are currently supported.
#   * Linux   
#       Libraries are GTK+, gtkglarea and OpenGL.
#       OpenGL is either Mesa or XFree86 (Ver.4 or later).
#   * FreeBSD 
#       Libraries are GTK+, gtkglarea and OpenGL.
#       OpenGL is either Mesa or XFree86 (Ver.4 or later).
#       GTK+ should be Ver. 1.2, otherwize, it may conflict with gtkglarea.
#   * Unix
#       Libraries are Motif and OpenGL.
#   * Microsoft Windows and Cygwin
#       Libraries are WIN32 and OpenGL in Windows.
#       Cygwin is Ver. 2.78 or later.
#   * software emulation, off-line rendering
#       AutoGL OS-X (Off-Screen eXtension)
#       No libraries are required.
#   * software emulation on server 
#       AutoGL.edu (client-server / web application edition)
#       No libraries are required.
#
# You may also need to install additional packages, such as Mesa, 
# gtkglarea or Cygwin.
# For any other platforms, you may also need to change this "Makefile".
#
#      
#
# Choose any one from these platforms below.
#
#   Linux with GTK+, OpenGL and gtkglarea 
PLATFORM = Linux
#   FreeBSD with GTK+, OpenGL and gtkglarea
#PLATFORM = FreeBSD
#   Unix with Motif and OpenGL
#PLATFORM = Unix
#   Microsoft Windows and Cygwin
#PLATFORM = Cygwin
#   software emulation, off-line rendering
#PLATFORM = Software
#   software emulation on server
#PLATFORM = Server
#



AUTOGL_HEADER_DIR = ../lib
AUTOGL_LIBRARY_DIR = ../lib

INCLUDES = -I$(AUTOGL_HEADER_DIR)
LIBS = -L$(AUTOGL_LIBRARY_DIR) -lautogl

# Directories of header files and libraries
ifeq (${PLATFORM}, Linux)
  # GTK+, gtkglarea and OpenGL (on Linux)
  LIBS += -lgtkgl -lGLU -lGL `gtk-config --libs`
endif
ifeq (${PLATFORM}, FreeBSD)
  # GTK+, gtkglarea and OpenGL (on FreeBSD, using GTK+ Ver.1.2.X)
  LIBS += -lgtkgl -lGLU -lGL `gtk12-config --libs`
endif
ifeq (${PLATFORM}, Unix)
  # Motif and OpenGL
  LIBS += -L/usr/X11R6/lib
  LIBS += -L/usr/X11R6/LessTif/Motif1.2/lib
#  LIBS += -lGLw -lGLU -lGL -lXm -lXt -lXmu -lXext -lSM -lICE -lX11 -lm
  LIBS += -lMesaGLw -lGLU -lGL -lXm -lXt -lXmu -lXext -lSM -lICE -lX11 -lm
endif
ifeq (${PLATFORM}, Cygwin)
  # WIN32 and OpenGL (on Microsoft Windows, using Cygwin)
  LIBS += -lglu32 -lopengl32 -luser32 -lgdi32
endif
ifeq (${PLATFORM}, Software)
  # software emulation, off-line rendering
  LIBS += -lm
endif
ifeq (${PLATFORM}, Server)
  # software emulation on server
  LIBS += -lm
endif

CC = g++

CFLAGS = $(INCLUDES) -g -Wall
#CFLAGS = $(INCLUDES) -O -Wall
#CFLAGS = $(INCLUDES) -O3 -Wall

LDFLAGS = $(LIBS)

.SUFFIXES: .cc

.cc.o:
	$(CC) -c $< $(CFLAGS)



all : sample

SAMPLE = sample1 sample2
sample : $(SAMPLE)



sample1.o   : 
sample1     : sample1.o 
	$(CC) -o $@ $^ $(LDFLAGS)

sample2.o   : 
sample2     : sample2.o 
	$(CC) -o $@ $^ $(LDFLAGS)



clean :
	-rm *.o *.exe $(SAMPLE)
