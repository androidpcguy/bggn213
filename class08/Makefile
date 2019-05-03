.PHONY=clean

all: Lecture08.html

%.html: %.Rmd
	@echo "Compiling R Markdown"
	@Rscript -e 'rmarkdown::render("$<")'
