#!/bin/sh

#
# Modify autogl/lib to work it together with Zarusoba Viewer
#
# (c) Yusa, Yasunori
#

if [ $# -ne 1 ]
then
    echo "Usage $(basename $0) /path/to/autogl/lib"
    exit
fi

path_to_autogl_lib=$1

# Change fonts
sed -i -e 's/a14/10x20/g' \
    ${path_to_autogl_lib}/autogl_gui_gtk_glx.c \
    ${path_to_autogl_lib}/autogl_gui_gtk_gtkglarea.c \
    ${path_to_autogl_lib}/autogl_gui_motif.c

# Enlarge AUTOGL_MAX_PROPERTY_ITEMS
sed -i -e 's/^#define AUTOGL_MAX_PROPERTY_ITEMS .*$/#define AUTOGL_MAX_PROPERTY_ITEMS 100/' \
    ${path_to_autogl_lib}/autogl_uims_inside.h

# Add math and X11 libraries
sed -i -e 's/^libs="${opengl_libs} ${gtk_libs}"/libs="${opengl_libs} ${gtk_libs} -lm -lX11"/' \
    ${path_to_autogl_lib}/advautogl-config__Linux

# Add AutoGL_SetLineWidth function
grep AutoGL_SetLineWidth ${path_to_autogl_lib}/autogl_graphics.h > /dev/null 2>&1
if [ $? -ne 0 ]
then
    sed -i -e 's/^void AutoGL_SetColor *$/void AutoGL_SetLineWidth (double width);  void AutoGL_SetColor/' \
        ${path_to_autogl_lib}/autogl_graphics.h
    sed -i -e 's/^void AutoGL_SetColor2 *$/void AutoGL_SetLineWidth2 (float width);  void AutoGL_SetColor2/' \
        ${path_to_autogl_lib}/autogl_graphics.h

    sed -i -e 's/^  PRIMITIVE_COLOR,/  PRIMITIVE_LINE_WIDTH,   PRIMITIVE_COLOR,/' \
        ${path_to_autogl_lib}/autogl_graphics.c
    sed -i -e 's/^void AutoGL_SetColor/void AutoGL_SetLineWidth (double width) {   if (View_DisplayListIsOpened) {     AutoGL_INSIDE_AddPrimitiveBuffer ((float)PRIMITIVE_LINE_WIDTH);      AutoGL_INSIDE_AddPrimitiveBuffer ((float)width);   } else {     AutoGL_SetLineWidth2       ((float)width);   } } static int SetLineWidth (float buffer[], int count) {   float width = buffer[count++];    AutoGL_SetLineWidth2 (width);   return count; }  void AutoGL_SetColor/' \
        ${path_to_autogl_lib}/autogl_graphics.c
    sed -i -e 's/^    case PRIMITIVE_COLOR:/    case PRIMITIVE_LINE_WIDTH:       count = SetLineWidth (buffer, count);       break;     case PRIMITIVE_COLOR:/' \
        ${path_to_autogl_lib}/autogl_graphics.c

    sed -i -e 's/^void AutoGL_SetColor2/void AutoGL_SetLineWidth2 (float width) {   GLfloat range[2];    SetDefaultPrimitive ();    glGetFloatv(GL_ALIASED_LINE_WIDTH_RANGE, range);   if (width < range[0]) width = range[0];   if (width > range[1]) width = range[1];    glLineWidth(width); }  void AutoGL_SetColor2/' \
        ${path_to_autogl_lib}/autogl_graphics_opengl.c

    sed -i -e 's/^void AutoGL_SetColor2/void AutoGL_SetLineWidth2 (float width) {   \/* do nothing *\/ }  void AutoGL_SetColor2/' \
        ${path_to_autogl_lib}/autogl_graphics_software.c
fi
