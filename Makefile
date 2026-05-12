.PHONY: env-annotation env-variants annotate variants clean help

help:
	@echo "Available commands:"
	@echo "  make env-annotation   — create conda environment for genome annotation"
	@echo "  make env-variants     — create conda environment for variant pipeline"
	@echo "  make annotate         — run genome annotation (annotation/)"
	@echo "  make variants         — run variant pipeline (variants/)"
	@echo "  make clean            — remove outputs/"

env-annotation:
	conda env create -f annotation/environment.yml

env-variants:
	conda env create -f variants/env/environment.yml

annotate:
	cd annotation && bash run_annotation.sh

variants:
	cd variants && bash workflow/run_all.sh

clean:
	rm -rf annotation/outputs variants/outputs
