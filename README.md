ASE2SPRKKR
==========

ASE2SPRKKR package provide an interface that allow use of the [SPR-KKR
package](https://www.ebert.cup.uni-muenchen.de/index.php/en/software-en/13-sprkkr)
to electronic structure calculation within [Atomic Simulation
Environment](https://wiki.fysik.dtu.dk/ase/) (abbreviated as ASE)
\-\--Python tool that integrates the various tools for electronic
structure calculation.

Usefull links
-------------

-   [SPR-KKR](https://software.pan-data.eu/software/111/spr-kkr)
-   [ASE - Atomic Simulation
    Environment](https://wiki.fysik.dtu.dk/ase/)
-   [Online version of this
    documentation](https://ase2sprkkr.github.io/ase2sprkkr/)
-   [GitHub repository of the ASE2SPRKKR
    package](https://github.com/ase2sprkkr/ase2sprkkr/)
-   [SPR-KKR
    manual](https://www.ebert.cup.uni-muenchen.de/index.php/en/repository/func-download/251/chk,b2f3ab5f57c7629207b121be0d31a38d/no_html,1/lang,en-gb/)

Installation of the package using package managers
--------------------------------------------------

The simplest way how to install and use the package is to install them
using package managers: either pip

``` {.bash}
pip install ase2sprkkr
```

or conda

``` {.bash}
conda install -c ase2sprkkr ase2sprkkr
```

### Further notes

In some systems, the pip utility for python3 is called pip3. If it is
not installed, you can install it using the linux distribution package
manager, e.g. in Debian/Ubuntu

``` {.bash}
apt install pip3
```

or

``` {.bash}
zypper install pip
```

in OpenSUSE

For the conda installation instructions, see the Anaconda documentation
<https://docs.anaconda.com/anaconda/install/linux/> however, for the
users unexperienced with conda, the (simpler) pip way is recommended.

Install the packages from GIT (and/or source codes)
---------------------------------------------------

If you do not want to use public package managers as pip or conda, or
you want to contribute to development, you can use GIT to obtain the
package sources.

### Requirements

-   Python \>= 3.7
-   SPR-KKR (not checked by the installer)
-   Python packages: ase, mendeleev, spglib, pyparsing
-   Git (to obtain the sources)

### Obtaining the package using GIT

``` {.bash}
git clone https://github.com/ase2sprkkr/ase2sprkkr.git
git checkout origin/release
```

The first line fetches the code of the package. The second one checks
out the recommended production version of the code.

If you want to obtain the current version of the (earlier-downloaded)
code, run

``` {.bash}
git fetch
git checkout origin/release
```

Alternatively, you can checkout master branch

``` {.bash}
git checkout origin/master
```

to obtain the bleeding edge version.

### Using the package (without installing the pip/conda packages)

You can either just add the src directory to your PYTHONPATH, or you can
build and install the package, as it is described below.

### Installation of the package from the sources

To install the package, you have to build the "wheel package" from the
sources and install it

#### Building the wheel (installation) package

If you do not have the wheel package built, you can do it with the
following steps.

``` {.bash}
python3 -m pip install --upgrade build
python3 -m build
```

The first line installs the tool to build the package (it is possible
that you have it already installed). The second one builds the package.

#### Installing the package

To install the package (either system-wide or in an active virtual
environment), you can run

``` {.bash}
pip install `ls ./dist/ase2sprkkr-*.whl | sort | tail -n 1`
```

### One step install from the sources

To do all the stuff (after cloning the GIT repository) in one step, you
can run

``` {.bash}
make
```

To clean up the source directory after installing the package, you can
run

``` {.bash}
make clean
```

Documentation of the package
----------------------------

The documentation is published online at
<https://ase2sprkkr.github.io/ase2sprkkr/>

If you are using Git cloned sources, you can run

    open docs/index.html

to see the (offline version of the) documentation. The documentation
contains parts, that are generated from the docstrings in the source
code. You can regenerate these by

``` {.bash}
make docs
```

Sphinx and sphinx-toolbox python packages (installable e.g. using pip)
and pandoc (for generating README.md) are needed for the build.

However, the official build of Sphinx miss some attributes when it is
used to build the documentation. So, till the pullrequest that corrects
the Sphinx behavior will be merged into Sphinx, please use the following
fork for building the documentation.
<https://github.com/lokik/sphinx.git>

How to contribute or to report a bug
------------------------------------

Please feel free to make a pull-request or post an issue at:
<https://github.com/ase2sprkkr/ase2sprkkr>
