# Build standalone TikZ figures for Chapter 12 into ../images/.
#
#   make            build every ../images/figNN.pdf that is out of date
#   make svg        also build ../images/figNN.svg from the PDFs
#   make fig12_1_3  build a single figure (PDF)
#   make clean      remove LaTeX aux files (keeps the PDFs/SVGs)
#   make -j4        build in parallel
#
# Each figNN.tex is a self-contained standalone document whose ONLY
# dependency is diagrams.sty (at the repo root); TEXINPUTS is pointed
# there so the source files can live in this subdirectory.
#
# Switch the LaTeX engine with e.g.  make ENGINE=xelatex  to match the
# main document's fonts exactly (pdflatex is the faster default).

ENGINE      ?= pdflatex
LATEXFLAGS  ?= -interaction=nonstopmode -halt-on-error
IMAGES      := ../images

export TEXINPUTS := ..:$(TEXINPUTS)

SRC := $(wildcard fig*.tex)
PDF := $(patsubst %.tex,$(IMAGES)/%.pdf,$(SRC))
SVG := $(patsubst %.tex,$(IMAGES)/%.svg,$(SRC))

.PHONY: all svg clean cleanpdf
all: $(PDF)
svg: $(SVG)

# Compile a figure and move the cropped PDF into ../images/.
$(IMAGES)/%.pdf: %.tex | $(IMAGES)
	$(ENGINE) $(LATEXFLAGS) $<
	@mv $*.pdf $@

# Vector SVG (needs pdf2svg).
$(IMAGES)/%.svg: $(IMAGES)/%.pdf
	pdf2svg $< $@

$(IMAGES):
	mkdir -p $(IMAGES)

# Allow "make fig12_1_3" as a shorthand for its PDF.
fig12_%: $(IMAGES)/fig12_%.pdf ;

clean:
	rm -f *.aux *.log *.out

cleanpdf:
	rm -f $(IMAGES)/fig12_*.pdf $(IMAGES)/fig12_*.svg
