PANDOC=$(shell which pandoc)
MANUSCRIPT=manuscript
BIB=Aging.bib

$(MANUSCRIPT).docx: $(MANUSCRIPT).md $(BIB) Makefile

%.docx: %.md 
	$(PANDOC) -t docx \
		--bibliography=$(BIB) \
		$< -o $@

html2pdf:
	RDIR="R.Aging/" ; \
	pandoc -f html -t latex +RTS -K64m -RTS -o $(RDIR)Histones.pdf $(RDIR)Histones.html ; \
	pandoc -f html -t latex +RTS -K64m -RTS -o $(RDIR)TFBSs.pdf $(RDIR)TFBSs.html ; \
	pandoc -f html -t latex +RTS -K64m -RTS -o $(RDIR)ChromStates.pdf $(RDIR)ChromStates.html && \
	s  -q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -sOutputFile=figures-tables/suppl_res_S1.pdf $(RDIR)Histones.pdf $(RDIR)TFBSs.pdf $(RDIR)ChromStates.pdf