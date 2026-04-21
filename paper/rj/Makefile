# Build the ivcheck R Journal paper to PDF.
#
# rticles::rjournal_article() regenerates RJwrapper.tex on every
# render, which clobbers the "Research article" patch. This Makefile
# runs the render, re-applies the patch, then compiles RJwrapper.tex
# directly with xelatex so the patched wrapper drives the final PDF.
#
# Usage:
#   make            # build RJwrapper.pdf from paper.Rmd
#   make clean      # remove LaTeX aux files
#   make distclean  # also remove the PDF

PANDOC := /Applications/quarto/bin/tools

.PHONY: all clean distclean

all: RJwrapper.pdf

paper.tex: paper.Rmd header.tex
	RSTUDIO_PANDOC=$(PANDOC) Rscript -e \
	  'rmarkdown::render("paper.Rmd", output_format = rticles::rjournal_article(), quiet = TRUE)'
	# Re-apply the "Research article" patch that rticles clobbers.
	sed -i '' 's/Contributed research article/Research article/' RJwrapper.tex

RJwrapper.pdf: paper.tex RJwrapper.tex RJournal.sty RJreferences.bib
	xelatex -interaction=nonstopmode RJwrapper.tex
	bibtex RJwrapper
	xelatex -interaction=nonstopmode RJwrapper.tex
	xelatex -interaction=nonstopmode RJwrapper.tex

clean:
	rm -f RJwrapper.aux RJwrapper.log RJwrapper.out RJwrapper.bbl \
	      RJwrapper.blg RJwrapper.toc paper.tex paper.knit.md

distclean: clean
	rm -f RJwrapper.pdf
