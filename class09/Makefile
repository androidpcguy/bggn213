.PHONY=clean

all: Lecture09.html

%.html: %.Rmd
	@echo "Compiling R Markdown"
	@Rscript -e 'rmarkdown::render("$<")'
