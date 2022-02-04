CC = gcc

CFLAGS   = -std=c99 -O2
LDFLAGS  =
INCLUDES = -I../lib
LIBS     = -L../lib -lm

ifeq ($(CC), gcc)
	CFLAGS += -Wall -Wextra -pedantic
endif

ifeq ($(CC), clang)
	CFLAGS += -Weverything
endif

TARGET = NURBS_input
OBJS   = NURBS_input.o ../lib/NURBS_calc.o
DEPEND = depend.d

.SUFFIXES: .c .o
.PHONY: all clean check-syntax depend

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDES) $<

clean:
	rm -f $(TARGET) $(TARGET).exe $(OBJS) $(DEPEND)

check-syntax:
	$(CC) -fsyntax-only $(CFLAGS) $(INCLUDES) $(CHK_SOURCES)

depend:
	$(CC) -MM $(addsuffix .c,$(basename $(OBJS))) > $(DEPEND)

-include $(DEPEND)
