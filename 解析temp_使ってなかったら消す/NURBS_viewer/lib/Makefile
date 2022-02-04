CC = gcc

CFLAGS   = -O2 -std=c99
LDFLAGS  =
INCLUDES =
LIBS     = -lm

CFLAGS += -Wall -Wextra -pedantic
DEPEND = depend.d

.SUFFIXES: .c .o
.PHONY: all clean depend

all: NURBS_calc.o

NURBS_calc.o : NURBS_calc.c
	$(CC) -c $(CFLAGS) $(INCLUDES) $<

clean:
	rm -f NURBS_calc.o $(OBJS) $(DEPEND)

depend:
	$(CC) -MM $(addsuffix .c,$(basename $(OBJS))) > $(DEPEND)

-include $(DEPEND)
