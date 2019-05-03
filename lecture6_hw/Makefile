.PHONY=clean

all: Lecture6.Rmd
	@Rscript -e 'library(rmarkdown); rmarkdown::render("Lecture6.Rmd", "html_document")'

clean:
	@rm -f Lecture6.html
