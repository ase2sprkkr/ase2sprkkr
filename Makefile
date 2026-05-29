all: install

install:
	pip install .

editable_install:
	pip install -C cmake.build-type=Debug --no-deps --editable .

ei: editable_install

clean:
	rm -rf ./dist

test:
	src/ase2sprkkr/run_tests

doc: doc-gather doc-build doc-readme

doc-gather: doc-clean
	(cd doc_src ; sphinx-apidoc -feM -o ./auto ../src/ase2sprkkr \
	"../src/ase2sprkkr/*/test" \
	"../src/ase2sprkkr/*/test/*" \
	)

doc-clean:
	rm -rf doc_src/auto/*
	rm -rf docs/*
	rm -rf docs/.??*

doc-debug:
	(cd docs;  sphinx-build -P doc_src .)

doc-build:
	sphinx-build -j auto doc_src docs/
	cp -r doc_src/_root/* doc_src/_root/.??* docs/

doc-readme:
	cd doc_src; pandoc README.rst -o ../README.md

package: | package_clean
	python -m build --sdist

package_clean:
	rm -rf dist/*

publish: build_clean build pip

pip: | package
	twine upload --username ase2sprkkr dist/*

anaconda:
	~/anaconda3/bin/anaconda login && ~/anaconda3/bin/anaconda upload "`ls ~/anaconda3/conda-bld/noarch/ase2sprkkr-* | tail -n 1`"

conda:
	rm -rf ~/anaconda3/conda-bld/src_cache/ase2sprkkr-*
	PWD="`pwd`" conda build .

checkout-docs:
	rm -rf docs
	git checkout docs
