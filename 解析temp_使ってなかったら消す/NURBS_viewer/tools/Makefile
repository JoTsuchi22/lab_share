CC = gcc

CFLAGS   = -Wall -O2
LDFLAGS  =
INCLUDES = -lm
LIBS     =

all: cfiles

CFILES = NURBS_pick NURBS_pick2 makeinput_hari makeinput_enkan convertJintToDisp

cfiles: $(CFILES)

NURBS_pick.o:
NURBS_pick: NURBS_pick.o
	$(CC) -o $@ $^ $(INCLUDES)

NURBS_pick2.o:
NURBS_pick2: NURBS_pick2.o
	$(CC) -o $@ $^ $(INCLUDES)

makeinput_hari.o:
makeinput_hari: makeinput_hari.o
	$(CC) -o $@ $^ $(INCLUDES)

makeinput_enkan.o:
makeinput_enkan: makeinput_enkan.o
	$(CC) -o $@ $^ $(INCLUDES)

convertJintToDisp.o:
convertJintToDisp: convertJintToDisp.o
	$(CC) -o $@ $^ $(INCLUDES)

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDES) $<

clean:
	rm -f $(TARGET1) $(OBJS) && rm *.o $(CFILES) $(DEPEND)

depend:
	$(CC) -MM $(addsuffix .c,$(basename $(OBJS))) > $(DEPEND)

-include $(DEPEND)
