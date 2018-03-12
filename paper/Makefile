PAPER = precision
TEX   = $(wildcard *.tex)
BIB   = $(PAPER).bib
FIG   = $(wildcard figures/*.pdf figures/*.png graphs/*.pdf graphs/*.png)

.PHONY: all clean

$(PAPER).pdf: $(TEX) $(BIB) $(FIG) cleanall
	pdflatex $(PAPER)
	pdflatex $(PAPER)
	bibtex $(PAPER)
	pdflatex $(PAPER)
	pdflatex $(PAPER)

clean:
	rm -f *.ilg *.aux *.log *.dvi *.idx *.toc *.lof *.lot
	rm -f *.blg *.bbl *~

cleanall: clean
	rm -f $(PAPER).out $(PAPER).pdf

update:
	git pull
