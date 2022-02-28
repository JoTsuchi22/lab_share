AUTOGL_DIR = ../lib

INCLUDES = -I$(AUTOGL_DIR)

LIBS = -L$(AUTOGL_DIR) -lautogl \
	-lgtkgl -lGLU -lGL `gtk-config --libs`



all : hello

hello   : hello.c
	gcc -o hello hello.c $(INCLUDES) $(LIBS)



clean :
	rm -f *.o *.exe hello
