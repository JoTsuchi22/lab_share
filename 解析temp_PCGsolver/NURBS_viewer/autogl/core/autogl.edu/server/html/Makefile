INCLUDES = 
LIBS = -lm

CC = gcc
CFLAGS = $(INCLUDES) -g -Wall
LDFLAGS = $(LIBS)

.SUFFIXES: .c

.c.o:
	$(CC) -c $< $(CFLAGS)



all : cgi_programs

CGI_PROGRAMS = autogl_frame.cgi autogl_panel.cgi autogl_view.cgi autogl_image.cgi
cgi_programs : $(CGI_PROGRAMS)



autogl_frame.o     : autogl_cgi.h
autogl_frame.cgi   : autogl_frame.o
	$(CC) -o $@ $^ $(LDFLAGS)

autogl_panel.o     : autogl_cgi.h
autogl_panel.cgi   : autogl_panel.o
	$(CC) -o $@ $^ $(LDFLAGS)

autogl_view.o     : autogl_cgi.h
autogl_view.cgi   : autogl_view.o
	$(CC) -o $@ $^ $(LDFLAGS)

autogl_image.o     : autogl_cgi.h
autogl_image.cgi   : autogl_image.o
	$(CC) -o $@ $^ $(LDFLAGS)



clean : 
	-rm *.o *.exe $(CGI_PROGRAMS)
