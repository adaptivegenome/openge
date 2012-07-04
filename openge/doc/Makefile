TEX       := texi2pdf
TEX_ARGS  := -c

TEX_FILES := $(shell find . -name "*.tex")
PDF_FILES := $(patsubst %.tex,%.pdf,$(TEX_FILES))


default:
	@echo "make what? Available targets are:"
	@echo "  pdf     -- Generate pdfs from tex"
	@echo "  clean   -- Remove the generated pdfs"
.PHONY: default


pdf: $(PDF_FILES)

.PHONY: pdf


%.pdf: %.tex
	$(TEX) $(TEX_ARGS) -o $@ $^

clean:
	rm -f $(PDF_FILES)
