.PHONY: env-annotation env-variants annotate variants clean help

help:
	@echo "Доступные команды:"
	@echo "  make env-annotation   — создать conda-окружение для аннотации генома"
	@echo "  make env-variants     — создать conda-окружение для вариантного пайплайна"
	@echo "  make annotate         — запустить аннотацию генома (annotation/)"
	@echo "  make variants         — запустить вариантный пайплайн (variants/)"
	@echo "  make clean            — удалить outputs/"

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
