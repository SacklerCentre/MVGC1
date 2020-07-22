NAME     = mvgc_schema
TEX      = $(NAME).tex
TEXFILES = $(TEX) lcbmaths.sty
AUXFILES = $(patsubst %.tex,%.aux,$(TEX))
PDF      = $(NAME).pdf
CLEAN    = $(AUXFILES) $(NAME).log $(NAME).end $(NAME).out $(PDF)
TEXFLAGS = --shell-escape --file-line-error --halt-on-error --interaction nonstopmode
LATEX    = pdflatex $(TEXFLAGS)

.PHONY: all again clean diag

all: $(PDF)

again:
	$(LATEX) $(TEX)

clean:
	rm -f $(CLEAN)

diag:
	@echo "Nothing to see here, please move along"

$(PDF): $(TEXFILES)
	$(LATEX) $(TEX)
