Introduction
============

This package is the Python interface for the SPR-KKR package for electronic structure computation,
build upon Atomic Simulation Environment (ASE) framework.

**Usefull links:**
* [https://software.pan-data.eu/software/111/spr-kkr](https://software.pan-data.eu/software/111/spr-kkr) - the documentation of the SPR-KKR
* [https://wiki.fysik.dtu.dk/ase/](https://wiki.fysik.dtu.dk/ase/) - the documentation of the ASE
* [https://ase2sprkkr.github.io/ase2sprkkr/](https://ase2sprkkr.github.io/ase2sprkkr/) - the online version of the documentation
* [https://github.com/ase2sprkkr/ase2sprkkr/](https://github.com/ase2sprkkr/ase2sprkkr/) - GitHub repository of the package
* [https://www.ebert.cup.uni-muenchen.de/index.php/en/repository/func-download/251/chk,b2f3ab5f57c7629207b121be0d31a38d/no_html,1/lang,en-gb/](https://www.ebert.cup.uni-muenchen.de/index.php/en/repository/func-download/251/chk,b2f3ab5f57c7629207b121be0d31a38d/no_html,1/lang,en-gb/) - the SPR-KKR manual

Installation of the package using package managers
==================================================

The simplest way how to install and use the package is
to install them using package managers: either pip
```Bash
pip install ase2sprkkr
```
or conda
```Bash
conda install -c ase2sprkkr ase2sprkkr
```

### Further notes

In some systems, the pip utility for python3 is called
pip3. If it is not installed, you can install it using the linux
distribution package manager, e.g. in Debian/Ubuntu
```Bash
apt install pip3
```
or
```Bash
zypper install pip
```
in OpenSUSE

For the conda installation instructions, see the Anaconda documentation
[https://docs.anaconda.com/anaconda/install/linux/](https://docs.anaconda.com/anaconda/install/linux/)
however, for the users unexperienced with conda,
the (simpler) pip way is recommended.


Install the packages from GIT (and/or source codes)
==========================================================
If you do not want to use public package managers as pip or conda,
or you want to contribute to development, you can use GIT to obtain
the package sources.

Requirements
------------
- Python >= 3.7
- SPR-KKR (not checked by the installer)
- Python packages: ase, mendeleev, spglib, pyparsing
- Git (to obtain the sources)

Obtaining the package using GIT
--------------------------------
```Bash
git clone https://github.com/ase2sprkkr/ase2sprkkr.git
git checkout origin/release
```

The first line fetches the code of the package. The second one
checks out the recommended production version of the code.

If you want to obtain the current version of the (earlier-downloaded)
code, run
```Bash
git fetch
git checkout origin/release
```
Alternatively, you can checkout master branch
```Bash
git checkout origin/master
```
to obtain the bleeding edge version.


Using the package (without installing the pip/conda packages)
---------------------------------------------------------------

You can either just add the src directory to your PYTHONPATH, or you
can build and install the package, as it is described below.


Installation of the package from the sources
----------------------------------------------
To install the package, you have to build the "wheel package" from
the sources and install it

### Building the wheel (installation) package
If you do not have the wheel package built, you can do it
with the following steps.

```Bash
python3 -m pip install --upgrade build
python3 -m build
```

The first line installs the tool to build the package
(it is possible that you have it already installed).
The second one builds the package.

### Installing the package
To install the package (either system-wide or in an active
virtual environment), you can run

```Bash
pip install `ls ./dist/ase2sprkkr-*.whl | sort | tail -n 1`
```

One step install from the sources
------------------------------------

To do all the stuff (after cloning the GIT repository) in one step,
you can run
```Bash
make
```

To clean up the source directory after installing the package,
you can run
```Bash
make clean
```

Documentation of the package
=============================
The documentation is published online at
[https://ase2sprkkr.github.io/ase2sprkkr/](https://ase2sprkkr.github.io/ase2sprkkr/)

If you are using Git cloned sources, you can run
```
open docs/index.html
```
to see the (offline version of the) documentation.
The documentation contains parts, that are generated from the docstrings in the
source code. You can regenerate these by
```Bash
make docs
```
Sphinx and md2html python packages (installable e.g. using pip) are needed for
the build.

How to contribute or to report a bug
======================================
Please feel free to make a pull-request or post an issue at:
[https://github.com/ase2sprkkr/ase2sprkkr](https://github.com/ase2sprkkr/ase2sprkkr)
