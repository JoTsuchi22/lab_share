AUTOGL_HEADER_DIR = ../../lib
AUTOGL_LIBRARY_DIR = ../../lib

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

SAMPLE = spheres city
sample : $(SAMPLE)



spheres.o : 
spheres   : spheres.o 
	$(CC) -o $@ $^ $(LDFLAGS)

city.o : 
city   : city.o 
	$(CC) -o $@ $^ $(LDFLAGS)



clean :
	rm -f *.o *.exe $(SAMPLE)
