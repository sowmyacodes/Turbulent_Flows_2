# Makefile for LATEX compiling

# Base filename of LATEX file (without .tex)
FNAME=MatRANS

all:
	latex $(FNAME)
	bibtex $(FNAME)
	latex $(FNAME)
	latex $(FNAME)

bib: 	
	latex $(FNAME)
	bibtex $(FNAME)

ps:
	latex $(FNAME)
	latex $(FNAME)
	dvips $(FNAME).dvi -Ppdf -G0 -o $(FNAME).ps
	gv $(FNAME).ps

dvi:
	latex $(FNAME)
	latex $(FNAME)
	xdvi $(FNAME)

pdf:

	latex $(FNAME)
	latex $(FNAME)
	dvips $(FNAME).dvi -Ppdf -G0 -o $(FNAME).ps

	ps2pdf -sPAPERSIZE=a4 \
	-dMaxSubsetPct=100 -dCompatibilityLevel=1.2 \
	-dSubsetFonts=true -dEmbedAllFonts=true \
	-dAutoFilterColorImages=false \
	-dAutoFilterGrayImages=false \
	-dColorImageFilter=/FlateEncode \
	-dGrayImageFilter=/FlateEncode \
	-dModoImageFilter=/FlateEncode \
	$(FNAME).ps $(FNAME).pdf
	acroread $(FNAME).pdf

clean:
#	rm -f *.pdf
	rm -f *.ps
	rm -f *.dvi
	rm -f *.log
	rm -f *.bbl
	rm -f *.blg
	rm -f *.aux
	rm -f *~