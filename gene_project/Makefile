# Author: Akshara Balachandra
# Date: 05/15/19
# Description: Makefile for markdown to pdf conversion of find a gene project

.PHONY=clean

all: project.md
	@echo "Compiling markdown..."
	@pandoc --template latex.tex project.md -o project.pdf

clean:
	@rm -f project.pdf

new: clean all
