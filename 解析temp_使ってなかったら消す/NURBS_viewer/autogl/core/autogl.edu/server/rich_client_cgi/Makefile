INCLUDES = 
LIBS = -lm

CC = gcc
CFLAGS = $(INCLUDES) -g -Wall
LDFLAGS = $(LIBS)

.SUFFIXES: .c

.c.o:
	$(CC) -c $< $(CFLAGS)



all : cgi_programs

CGI_PROGRAMS = autogl_post_command.cgi
cgi_programs : $(CGI_PROGRAMS)



autogl_post_command.o     : 
autogl_post_command.cgi   : autogl_post_command.o
	$(CC) -o $@ $^ $(LDFLAGS)



clean : 
	-rm *.o *.exe $(CGI_PROGRAMS)
