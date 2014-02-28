.PHONY: all clean

PDFIFY = pdflatex
TEXIFY = latex

TEXSOURCES = $(wildcard *.tex)
PDF = $(patsubst %.tex,%.pdf, $(TEXSOURCES) )
TOC = $(patsubst %.tex,%.toc, $(TEXSOURCES) )

all: $(PDF)

$(PDF): $(TEXSOURCES)
	$(PDFIFY) $(TEXSOURCES)  


clean:
	rm -f *~
	rm -f *.out
	rm -f *.aux
	rm -f *.log
	rm -f *.dvi
	rm -f *.nav
	rm -f *.snm
	rm -f *.toc
