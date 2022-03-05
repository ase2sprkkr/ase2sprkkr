all:
	./install.sh

clean:
	rm -rf ./dist

test:
	src/ase2sprkkr/run_tests

doc: doc-gather doc-build

doc-gather: doc-clean doc-md
	sphinx-apidoc -fe -o sphinx/auto src/ase2sprkkr

doc-clean:
	rm -rf sphinx/auto/*

doc-build:
	sphinx-build sphinx docs/gen
	sed -i 's#./gen/#./#' docs/gen/_static/documentation.html

doc-md: README.md
	md2html README.md -f docs/README.md.html
	sed -r 's/###(.*)/\1\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^/' README.md > sphinx/README.md

build: | build_clean
	python3 -m build

build_clean:
	rm -rf dist/*

publish: build_clean build pip

pip: | build
	twine upload dist/*

conda: | build
	conda build .
