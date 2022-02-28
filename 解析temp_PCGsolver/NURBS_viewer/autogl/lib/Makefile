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
#       Libraries : Gtk+, OpenGL and GLX
#         OpenGL is either Mesa or XFree86 (Ver.4 or later)
#         Gtk+ is Ver. 2.*
#   * Linux (using old Gtk+)
#       Libraries : Gtk+, OpenGL and gtkglarea
#         OpenGL is either Mesa or XFree86 (Ver.4 or later)
#         Gtk+ is Ver. 1.2.*
#         gtkglarea is Ver. 1.2.*
#   * Unix
#       Libraries : Motif and OpenGL
#   * Microsoft Windows and Cygwin
#       Libraries : WIN32 and OpenGL in Windows
#   * software emulation, off-line rendering
#       AutoGL OS-X (Off-Screen eXtension)
#       Libraries : nothing
#   * software emulation on server 
#       AutoGL.edu (client-server / web application edition)
#       Libraries : nothing
#
# AutoGL-MP is a parallel rendering edition.
# It uses MPI for inter-process communication.
# As MPI implementation, either LAM or mpich is supported.
#   * Linux
#       Libraries : Gtk+, OpenGL and GLX
#   * Linux (using old Gtk+)
#       Libraries : Gtk+, OpenGL and gtkglarea
#   * software emulation, off-line rendering
#       AutoGL OS-X (Off-Screen eXtension)
#       Libraries : nothing
#
# AutoGL-MMX (Multiple Monitor eXtension) is for tiling and 
#   stereo projection using graphics cluster.
# It uses MPI for inter-process communication.
# As MPI implementation, either LAM or mpich is supported.
#   * Linux
#       Libraries : Gtk+, OpenGL and GLX
#   * Linux (using old Gtk+)
#       Libraries : Gtk+, OpenGL and gtkglarea
#
# You may also need to install additional packages, such as Mesa, 
# gtkglarea or Cygwin.
# For any other platforms, you may also need to change this "Makefile".
#
#      
#
# Supported platforms for AutoGL
#   Choose any one from these platforms below.
#
#   Linux with Gtk+, OpenGL and GLX
PLATFORM = Linux
#   Linux with Gtk+, OpenGL and gtkglarea 
#PLATFORM = Linux_old
#   Unix with Motif and OpenGL
#PLATFORM = Unix
#   Microsoft Windows and Cygwin
#PLATFORM = Cygwin
#   software emulation, off-line rendering
#PLATFORM = Software
#   software emulation on server
#PLATFORM = Server
#
#
#
# Supported platforms for AutoGL-MP (parallel rendering)
#
#   Linux with Gtk+, OpenGL and GLX
#PLATFORM = Linux_MP
#   Linux with Gtk+, OpenGL and gtkglarea 
#PLATFORM = Linux_old_MP
#   software emulation, off-line rendering 
#PLATFORM = Software_MP
#
#
#
# Supported platforms for AutoGL-MMX (multiple monitor and stereo projection)
#   
#
#   Linux with Gtk+, OpenGL and GLX
#PLATFORM = Linux_MMX
#   Linux with Gtk+, OpenGL and gtkglarea 
#PLATFORM = Linux_old_MMX
#
#



CC = gcc

#
# Directories of header files and libraries
#

#
# for AutoGL (normal edition)
#
ifeq (${PLATFORM}, Linux)
  # Gtk+, OpenGL and GLX (on Linux)
  INCLUDES += `pkg-config --cflags gtk+-2.0`
  CONFIG = advautogl-config__Linux
endif
ifeq (${PLATFORM}, Linux_old)
  # Gtk+, OpenGL and gtkglarea (on Linux)
  INCLUDES += `gtk-config --cflags`
  CONFIG = advautogl-config__Linux_old
endif
ifeq (${PLATFORM}, Unix)
  # Motif and OpenGL
  INCLUDES += -I/usr/X11R6/include
  INCLUDES += -I/usr/X11R6/LessTif/Motif1.2/include
  CONFIG = advautogl-config__Unix
endif
ifeq (${PLATFORM}, Cygwin)
  # WIN32 and OpenGL (on Microsoft Windows, using Cygwin)
  INCLUDES += 
  CONFIG = advautogl-config__Cygwin
endif
ifeq (${PLATFORM}, Software)
  # software emulation, off-line rendering
  INCLUDES += 
  CONFIG = advautogl-config__Software
endif
ifeq (${PLATFORM}, Server)
  # software emulation on server
  INCLUDES += 
  CONFIG = advautogl-config__Server
endif

#
# for AutoGL-MP (parallel rendering)
#
ifeq (${PLATFORM}, Linux_MP)
  # Gtk+, OpenGL and GLX (on Linux)
  INCLUDES += `pkg-config --cflags gtk+-2.0`
  INCLUDES += -I../../autonosim/lib
  CC = mpicc -DAUTOGL_MP
  CONFIG = advautogl-config__Linux
endif
ifeq (${PLATFORM}, Linux_old_MP)
  # Gtk+, OpenGL and gtkglarea (on Linux)
  INCLUDES += `gtk-config --cflags`
  INCLUDES += -I../../autonosim/lib
  CC = mpicc -DAUTOGL_MP
  CONFIG = advautogl-config__Linux_old
endif
ifeq (${PLATFORM}, Software_MP)
  # software emulation, off-line rendering
  INCLUDES += 
  INCLUDES += -I../../autonosim/lib
  CC = mpicc -DAUTOGL_MP
  CONFIG = advautogl-config__Software
endif

#
# for AutoGL-MMX (multiple monitor and stereo projection)
#
ifeq (${PLATFORM}, Linux_MMX)
  # Gtk+, OpenGL and GLX (on Linux)
  INCLUDES += `pkg-config --cflags gtk+-2.0`
  INCLUDES += -I../../autonosim/lib
  CC = mpicc -DAUTOGL_MMX
  CONFIG = advautogl-config__Linux
endif
ifeq (${PLATFORM}, Linux_old_MMX)
  # Gtk+, OpenGL and gtkglarea (on Linux)
  INCLUDES += `gtk-config --cflags`
  INCLUDES += -I../../autonosim/lib
  CC = mpicc -DAUTOGL_MMX
  CONFIG = advautogl-config__Linux_old
endif



CFLAGS = $(INCLUDES) -g -Wall 
#CFLAGS = $(INCLUDES) -DNDEBUG -O
#CFLAGS = $(INCLUDES) -DNDEBUG -O3

AR = ar 
ARFLAGS = crsv

.SUFFIXES: .c

.c.o:
	$(CC) -c $< $(CFLAGS)



all : libautogl.a libautoglf.a advautogl-config



AUTOGL_HDRS = autogl_os.h autogl_math.h \
	autogl_uims.h autogl_graphics.h autogl_gui.h \
	autogl_cui.h \
	autogl_utility.h \
	autogl_contour.h autogl_vislow.h autogl_visualization.h \
	autogl_utility2d.h autogl_utility3d.h autogl.h 
AUTOGL_OBJS = autogl_os.o autogl_math.o \
	autogl_uims.o autogl_graphics.o autogl_gui.o \
	autogl_callback.o \
	autogl_cui.o \
	autogl_utility.o \
	autogl_contour.o autogl_vislow.o autogl_visualization.o \
	autogl_utility2d.o autogl_utility3d.o

#
# for AutoGL (normal edition)
#
ifeq (${PLATFORM}, Linux)
  # Gtk+, OpenGL and GLX
  AUTOGL_OBJS += autogl_graphics_opengl.o autogl_gui_gtk_glx.o
endif
ifeq (${PLATFORM}, Linux_old)
  # Gtk+, OpenGL and gtkglarea
  AUTOGL_OBJS += autogl_graphics_opengl.o autogl_gui_gtk_gtkglarea.o
endif
ifeq (${PLATFORM}, Unix)
  # Motif and OpenGL
  AUTOGL_OBJS += autogl_graphics_opengl.o autogl_gui_motif.o
endif
ifeq (${PLATFORM}, Cygwin)
  # WIN32 and OpenGL (on Microsoft Windows, using Cygwin)
  AUTOGL_OBJS += autogl_graphics_opengl.o autogl_gui_windows.o
endif
ifeq (${PLATFORM}, Software)
  # software emulation, off-line rendering
  AUTOGL_OBJS += autogl_graphics_software.o autogl_gui_software.o
endif
ifeq (${PLATFORM}, Server)
  # software emulation on server
  AUTOGL_OBJS += autogl_graphics_software.o autogl_gui_server.o
endif

#
# for AutoGL-MP (parallel rendering)
#
ifeq (${PLATFORM}, Linux_MP)
  # Gtk+, OpenGL and GLX
  AUTOGL_OBJS += autogl_graphics_opengl.o autogl_gui_gtk_glx.o
  AUTOGL_OBJS += autogl_com.o
endif
ifeq (${PLATFORM}, Linux_old_MP)
  # Gtk+, OpenGL and gtkglarea
  AUTOGL_OBJS += autogl_graphics_opengl.o autogl_gui_gtk_gtkglarea.o
  AUTOGL_OBJS += autogl_com.o
endif
ifeq (${PLATFORM}, Unix_MP)
  # Motif and OpenGL
  AUTOGL_OBJS += autogl_graphics_opengl.o autogl_gui_motif.o
  AUTOGL_OBJS += autogl_com.o
endif
ifeq (${PLATFORM}, Software_MP)
  # software emulation, off-line rendering
  AUTOGL_OBJS += autogl_graphics_software.o autogl_gui_software.o
  AUTOGL_OBJS += autogl_com.o
endif

#
# for AutoGL-MMX (multiple monitor and stereo projection)
#
ifeq (${PLATFORM}, Linux_MMX)
  # Gtk+, OpenGL and GLX
  AUTOGL_OBJS += autogl_graphics_opengl.o autogl_gui_gtk_glx.o
  AUTOGL_OBJS += autogl_com.o
endif
ifeq (${PLATFORM}, Linux_old_MMX)
  # Gtk+, OpenGL and gtkglarea 
  AUTOGL_OBJS += autogl_graphics_opengl.o autogl_gui_gtk_gtkglarea.o
  AUTOGL_OBJS += autogl_com.o
endif

autogl_os.o         : autogl_os.h 

autogl_math.o       : autogl_os.h autogl_math.h

autogl_uims.o       : autogl_os.h autogl_math.h \
	autogl_uims.h autogl_uims_inside.h 

autogl_graphics.o   : autogl_os.h autogl_math.h \
	autogl_uims.h autogl_uims_inside.h \
	autogl_graphics.h autogl_graphics_inside.h 

autogl_gui.o        : autogl_os.h autogl_math.h \
	autogl_uims.h autogl_uims_inside.h \
	autogl_graphics.h autogl_graphics_inside.h  \
	autogl_gui.h autogl_gui_inside.h 

autogl_com.o        : autogl_os.h autogl_math.h \
	autogl_uims.h autogl_uims_inside.h \
	autogl_graphics.h autogl_graphics_inside.h \
	autogl_gui.h autogl_gui_inside.h \
	autogl_com_inside.h

autogl_callback.o   : autogl_os.h autogl_math.h \
	autogl_uims.h autogl_uims_inside.h \
	autogl_graphics.h autogl_graphics_inside.h \
	autogl_gui.h autogl_gui_inside.h \
	autogl_com_inside.h \
	autogl_callback_inside.h 

autogl_cui.o        : autogl_os.h autogl_math.h \
	autogl_uims.h autogl_uims_inside.h \
	autogl_graphics.h autogl_graphics_inside.h \
	autogl_gui.h autogl_gui_inside.h \
	autogl_callback_inside.h \
	autogl_cui.h 

autogl_graphics_opengl.o   : autogl_os.h autogl_math.h \
	autogl_uims.h autogl_uims_inside.h \
	autogl_graphics.h autogl_graphics_inside.h \
	autogl_com_inside.h \
	autogl_callback_inside.h \
	autogl_graphics_opengl_inside.h
autogl_graphics_software.o   : autogl_os.h autogl_math.h \
	autogl_uims.h autogl_uims_inside.h \
	autogl_graphics.h autogl_graphics_inside.h \
	autogl_com_inside.h \
	autogl_callback_inside.h \
	autogl_graphics_software_inside.h

autogl_gui_gtk_glx.o           : autogl_os.h autogl_math.h \
	autogl_uims.h autogl_uims_inside.h \
	autogl_graphics.h autogl_graphics_inside.h \
	autogl_gui.h autogl_gui_inside.h \
	autogl_com_inside.h \
	autogl_callback_inside.h \
	autogl_graphics_opengl_inside.h 
autogl_gui_gtk_gtkglarea.o           : autogl_os.h autogl_math.h \
	autogl_uims.h autogl_uims_inside.h \
	autogl_graphics.h autogl_graphics_inside.h \
	autogl_gui.h autogl_gui_inside.h \
	autogl_com_inside.h \
	autogl_callback_inside.h \
	autogl_graphics_opengl_inside.h 
autogl_gui_motif.o         : autogl_os.h autogl_math.h \
	autogl_uims.h autogl_uims_inside.h \
	autogl_graphics.h autogl_graphics_inside.h \
	autogl_gui.h autogl_gui_inside.h \
	autogl_com_inside.h \
	autogl_callback_inside.h \
	autogl_graphics_opengl_inside.h 
autogl_gui_windows.o       : autogl_os.h autogl_math.h \
	autogl_uims.h autogl_uims_inside.h \
	autogl_graphics.h autogl_graphics_inside.h \
	autogl_gui.h autogl_gui_inside.h \
	autogl_com_inside.h \
	autogl_callback_inside.h \
	autogl_graphics_opengl_inside.h 
autogl_gui_software.o      : autogl_os.h autogl_math.h \
	autogl_uims.h autogl_uims_inside.h \
	autogl_graphics.h autogl_graphics_inside.h \
	autogl_gui.h autogl_gui_inside.h \
	autogl_com_inside.h \
	autogl_callback_inside.h \
	autogl_cui.h \
	autogl_graphics_software_inside.h 
autogl_gui_server.o        : autogl_os.h autogl_math.h \
	autogl_uims.h autogl_uims_inside.h \
	autogl_graphics.h autogl_graphics_inside.h \
	autogl_gui.h autogl_gui_inside.h \
	autogl_callback_inside.h \
	autogl_cui.h \
	autogl_graphics_software_inside.h 

autogl_utility.o    : autogl_os.h autogl_math.h \
	autogl_uims.h autogl_graphics.h autogl_gui.h autogl_cui.h \
	autogl_utility.h
autogl_contour.o    : autogl_os.h autogl_math.h \
	autogl_uims.h autogl_graphics.h autogl_gui.h autogl_cui.h \
	autogl_contour.h
autogl_vislow.o    : autogl_os.h autogl_math.h \
	autogl_uims.h autogl_graphics.h autogl_gui.h autogl_cui.h \
	autogl_contour.h autogl_vislow.h
autogl_visualization.o    : autogl_os.h autogl_math.h \
	autogl_uims.h autogl_graphics.h autogl_gui.h autogl_cui.h \
	autogl_contour.h autogl_vislow.h autogl_visualization.h
autogl_utility2d.o  : autogl_os.h autogl_math.h \
	autogl_uims.h autogl_graphics.h autogl_gui.h autogl_cui.h \
	autogl_utility2d.h
autogl_utility3d.o  : autogl_os.h autogl_math.h \
	autogl_uims.h autogl_graphics.h autogl_gui.h autogl_cui.h \
	autogl_utility3d.h

libautogl.a : $(AUTOGL_OBJS)
	$(AR) $(ARFLAGS) $@ $^



autogl_fortran.o  : autogl_os.h autogl_math.h \
	autogl_uims.h autogl_graphics.h autogl_gui.h autogl_cui.h \
	autogl_utility.h autogl_utility2d.h autogl_utility3d.h

libautoglf.a : $(AUTOGL_OBJS) autogl_fortran.o
	$(AR) $(ARFLAGS) $@ $^



advautogl-config :
	cp -f  $(CONFIG)  advautogl-config



clean :
	rm -f *.o *.a
	rm -f advautogl-config
