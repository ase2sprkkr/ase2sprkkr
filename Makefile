all:
	./install.sh

clean:
	rm -rf ./dist

test:
	src/run_tests

docs: doc-gather doc-build

doc-gather: doc-clean doc-md
	sphinx-apidoc -fe -o sphinx/auto src/ase2sprkkr

doc-clean:
	rm -rf sphinx/auto/*

doc-build:
	sphinx-build sphinx doc/gen

doc-md: README.md
	md2html README.md -f doc/README.md.html
	sed -r 's/###(.*)/\1\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^/' README.md > sphinx/README.md