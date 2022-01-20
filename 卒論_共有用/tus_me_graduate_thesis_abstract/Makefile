TEX    = platex
DVIPDF = dvipdfmx

TEXFLAGS    = -kanji=utf8 -file-line-error -halt-on-error -interaction=nonstopmode
DVIPDFFLAGS = -f ptex-ipaex.map

SOURCES = main.tex
TARGET  = abstract_7518074.pdf

JOBNAME = $(basename $(TARGET))

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(JOBNAME).dvi
	$(DVIPDF) $(DVIPDFFLAGS) $<

$(JOBNAME).dvi: $(SOURCES)
	$(TEX) $(TEXFLAGS) -jobname=$(JOBNAME) $<
	while grep '^LaTeX Warning: Label(s) may have changed.' $(JOBNAME).log > /dev/null; \
	do \
		$(TEX) $(TEXFLAGS) -jobname=$(JOBNAME) $<; \
	done

clean:
	rm -f $(TARGET) $(JOBNAME).dvi
	rm -f $(JOBNAME).log $(JOBNAME).aux
	rm -f $(addsuffix .aux,$(basename $(SOURCES)))
