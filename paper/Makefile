PAPER = precision
TEX   = $(wildcard *.tex)
BIB   = 99-$(PAPER).bib
FIG   = $(wildcard figures/*.pdf figures/*.png graphs/*.pdf graphs/*.png)

.PHONY: all clean

$(PAPER).pdf: $(TEX) $(BIB) $(FIG) cleanall
	latex $(PAPER)
	latex $(PAPER)
	bibtex $(PAPER)
	latex $(PAPER)
	latex $(PAPER)
	dvipdfmx $(PAPER).dvi
	rm $(PAPER).dvi

clean:
	rm -f *.ilg *.aux *.log *.dvi *.idx *.toc *.lof *.lot
	rm -f *.blg *.bbl *~

cleanall: clean
	rm -f $(PAPER).out $(PAPER).pdf

update:
	git pull
