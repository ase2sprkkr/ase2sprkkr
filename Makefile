all: install

install:
	./install.sh

clean:
	rm -rf ./dist

test:
	src/ase2sprkkr/run_tests

doc: doc-gather doc-build doc-readme

doc-gather: doc-clean
	sphinx-apidoc -feM -o sphinx/auto src/ase2sprkkr */test/

doc-clean:
	rm -rf sphinx/auto/*
	rm -rf docs/*
	rm -rf docs/.??*

doc-debug:
	sphinx-build -P sphinx docs/

doc-build:
	sphinx-build -j auto sphinx docs/
	cp -r sphinx/_root/* sphinx/_root/.??* docs/

doc-readme:
	cd sphinx; pandoc README.rst -o ../README.md

build: | build_clean
	python3 -m build

build_clean:
	rm -rf dist/*

publish: build_clean build pip

pip: | build
	twine upload dist/*

anaconda:
	~/anaconda3/bin/anaconda login && ~/anaconda3/bin/anaconda upload "`ls ~/anaconda3/conda-bld/noarch/ase2sprkkr-* | tail -n 1`"

conda:
	rm -rf ~/anaconda3/conda-bld/src_cache/ase2sprkkr-*
	PWD="`pwd`" conda build .
