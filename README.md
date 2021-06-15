Introduction
============

This package is the Python interface for SPR-KKR package for electronic structure computation,
build upon Atomic Simulation Environment (ASE) framework.

See
 - https://software.pan-data.eu/software/111/spr-kkr for documentation of the SPR-KRR and
 - https://wiki.fysik.dtu.dk/ase/

Installation of the devel package
=================================

Requirements
------------

- Python >= 3.8

Optional, to obtain and install the package
- Git (to obtain the sources)
- Pip (to install the package)


Obtaining the package
---------------------

```Bash
git clone github.com:ramntczcu/ase2sprkkr.git
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

Using the package
-----------------
You can either just add the src directory to your PYTHONPATH, or you
can build and install the package, as it is described below.

Building the wheel (instalation) package
----------------------------------------
If you do not have the wheel package built, you can do it
with the following steps.

```Bash
python3 -m pip install --upgrade build
python3 -m build
```
The first line installs the tool to build the package
(it is possible that you have it already installed).
The second one build the package

Installation the package
------------------------
To install the package (either system-wide, or in an active
virtual environment), run
```Bash
pip install `ls ./dist/ase2sprkkr-*.whl | sort | tail -n 1`
```
Note, that in some systems, the pip utility for python3 is called
pip3. If it is not installed, you can install it using the system
package manager, e.g. in Debian/Ubuntu
```Bash
apt install pip3
```
or
```Bash
zypper install pip
```
in OpenSUSE

One step install
----------------
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
