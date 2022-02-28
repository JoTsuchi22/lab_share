AUTOGL_HEADER_DIR = ../lib
AUTOGL_LIBRARY_DIR = ../lib

INCLUDES = -I$(AUTOGL_HEADER_DIR)
LIBS = -L$(AUTOGL_LIBRARY_DIR) -lautogl

# invoke AutoGL config to add system-level libraries for GUI and graphics
LIBS += `${AUTOGL_LIBRARY_DIR}/advautogl-config`



CC = gcc

CFLAGS = $(INCLUDES) -g -Wall
#CFLAGS = $(INCLUDES) -O -Wall
#CFLAGS = $(INCLUDES) -O3 -Wall

LDFLAGS = $(LIBS)

.SUFFIXES: .c

.c.o:
	$(CC) -c $< $(CFLAGS)



all : sample

SAMPLE = primitive more_primitive \
	contour graph3d line_contour band_contour isosurface section  \
	panel grid2d grid3d fem2d \
	pick2d pick3d idle2d selection2d selection3d \
	heavy heavy_fast batch
sample : $(SAMPLE)



primitive.o : 
primitive   : primitive.o 
	$(CC) -o $@ $^ $(LDFLAGS)

more_primitive.o : 
more_primitive   : more_primitive.o 
	$(CC) -o $@ $^ $(LDFLAGS)

contour.o : 
contour   : contour.o 
	$(CC) -o $@ $^ $(LDFLAGS)

graph3d.o : 
graph3d   : graph3d.o 
	$(CC) -o $@ $^ $(LDFLAGS)

line_contour.o : 
line_contour   : line_contour.o 
	$(CC) -o $@ $^ $(LDFLAGS)

band_contour.o : 
band_contour   : band_contour.o 
	$(CC) -o $@ $^ $(LDFLAGS)

isosurface.o : 
isosurface   : isosurface.o 
	$(CC) -o $@ $^ $(LDFLAGS)

section.o : 
section   : section.o 
	$(CC) -o $@ $^ $(LDFLAGS)

panel.o : 
panel   : panel.o 
	$(CC) -o $@ $^ $(LDFLAGS)

grid2d.o : 
grid2d   : grid2d.o 
	$(CC) -o $@ $^ $(LDFLAGS)

grid3d.o : 
grid3d   : grid3d.o 
	$(CC) -o $@ $^ $(LDFLAGS)

fem2d.o : 
fem2d   : fem2d.o 
	$(CC) -o $@ $^ $(LDFLAGS)

pick2d.o : 
pick2d   : pick2d.o 
	$(CC) -o $@ $^ $(LDFLAGS)

pick3d.o : 
pick3d   : pick3d.o 
	$(CC) -o $@ $^ $(LDFLAGS)

idle2d.o : 
idle2d   : idle2d.o 
	$(CC) -o $@ $^ $(LDFLAGS)

selection2d.o : 
selection2d   : selection2d.o 
	$(CC) -o $@ $^ $(LDFLAGS)

selection3d.o : 
selection3d   : selection3d.o 
	$(CC) -o $@ $^ $(LDFLAGS)

heavy.o : 
heavy   : heavy.o 
	$(CC) -o $@ $^ $(LDFLAGS)

heavy_fast.o : 
heavy_fast   : heavy_fast.o 
	$(CC) -o $@ $^ $(LDFLAGS)

batch.o : 
batch   : batch.o 
	$(CC) -o $@ $^ $(LDFLAGS)



clean :
	rm -f *.o *.exe $(SAMPLE)
