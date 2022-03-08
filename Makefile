all:
	./install.sh

clean:
	rm -rf ./dist

test:
	src/ase2sprkkr/run_tests

doc: doc-gather doc-md doc-build

doc-gather: doc-clean
	sphinx-apidoc -feM -o sphinx/auto src/ase2sprkkr */test/

doc-clean:
	rm -rf sphinx/auto/*
	rm -rf docs/gen/*

doc-build:
	sphinx-build sphinx docs/gen

doc-md: README.md
	md2html README.md -f docs/README.md.html
	sed 's/<link .*<\/link>//' docs/README.md.html > sphinx/_gen/README.md.html
	sed 's#./gen/#./#' docs/index.html > sphinx/_gen/documentation.html

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
