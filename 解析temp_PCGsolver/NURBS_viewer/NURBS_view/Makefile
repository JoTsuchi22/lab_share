#ADVENTURE = ${HOME}/ADVENTURE

CC = gcc
# CC = clang

CFLAGS   = -O2 -std=c99
LDFLAGS  = -mcmodel=large
INCLUDES =
LIBS     = -lm -lpng

ifeq ($(CC), gcc)
	CFLAGS += -Wall -Wextra -pedantic
endif

ifeq ($(CC), clang)
	CFLAGS += -Weverything
endif

#ifneq ($(ADVENTURE),)
#	CFLAGS   += -DWITH_ADVIO
#	INCLUDES += `$(ADVENTURE)/bin/advsys-config --cflags`
#	LIBS     += `$(ADVENTURE)/bin/advsys-config --libs docio`
#endif

INCLUDES += -I../autogl/lib -I../lib
LIBS     += -L../autogl/lib -lautogl `../autogl/lib/advautogl-config`


TARGET = NURBS_view
OBJS   = NURBS_view.o plot.o draw.o view.o image.o pick.o
OBJS   += clickcalc.o ../lib/NURBS_calc.o
#OBJS   = main.o pick.o view.o plot.o draw.o
#OBJS  += option.o
#OBJS  += mesh.o face.o edge.o crack.o
#OBJS  += field.o
#OBJS  += point.o
#OBJS  += scalar.o
#OBJS  += image.o
#OBJS  += animation.o
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
