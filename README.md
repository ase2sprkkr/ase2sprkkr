# ASE2SPRKKR

ASE2SPRKKR package provide an interface that allow use of the [SPR-KKR
package](https://www.ebert.cup.uni-muenchen.de/index.php/en/software-en/13-sprkkr)
to electronic structure calculation within [Atomic Simulation
Environment](https://wiki.fysik.dtu.dk/ase/) (abbreviated as ASE)
\-\--Python tool that integrates the various tools for electronic
structure calculation.

## Usefull links

-   [SPR-KKR](https://ebert.cup.uni-muenchen.de/index.php/en/software-en/13-sprkkr)
-   [ASE - Atomic Simulation
    Environment](https://wiki.fysik.dtu.dk/ase/)
-   [Online version of this
    documentation](https://ase2sprkkr.github.io/ase2sprkkr/)
-   [GitHub repository of the ASE2SPRKKR
    package](https://github.com/ase2sprkkr/ase2sprkkr/)
-   [SPR-KKR
    manual](https://www.ebert.cup.uni-muenchen.de/index.php/en/repository/func-startdown/251/lang,en-gb/)

## Installation of the package using package managers

The simplest way how to install and use the package is to install them
using package managers: either pip

``` bash
pip install ase2sprkkr
```

or conda

``` bash
conda install -c ase2sprkkr ase2sprkkr
```

I reccomend to install beta version: there can be some bugs, but mostly
it has more bugs repaired than introduced, moreover, you can enjoy new
properties. The beta versions are available only through `pip`:

``` bash
pip install --pre ase2sprkkr
```

To use bleading edge sources (the newest features, but you risk to
encounter bugs), you can install the packages from github:

``` bash
pip install git+https://github.com/ase2sprkkr/ase2sprkkr.git@develop
```

### Further notes

In some systems, the `pip` utility for `python3` is called `pip3`. It
may be possible, that `pip` is installed, but it is not in `PATH`. In
such case, the pip utility can be runned using `python -m pip` or
`python3 -m pip3`.

If `pip` is not installed, you can install it using the linux
distribution package manager, e.g. in Debian/Ubuntu

``` bash
apt install pip3
```

or

``` bash
zypper install pip
```

For the conda installation instructions, see the Anaconda documentation
<https://docs.anaconda.com/anaconda/install/linux/> however, for the
users unexperienced with conda, the (simpler) pip way is recommended.

## Install the packages from GIT (and/or source codes)

If you do not want to use public package managers as pip or conda, or
you want to contribute to development, you can use GIT to obtain the
package sources.

### Requirements

-   Python \>= 3.8
-   SPR-KKR (not checked by the installer)
-   Python packages: see the the pypoject.toml file
-   Git (to obtain the sources)

### Obtaining the package using GIT

``` bash
git clone https://github.com/ase2sprkkr/ase2sprkkr.git
git checkout release
```

The first line fetches the code of the package. The second one checks
out the stable (production) version of the code.

If you want to obtain the current version of the (earlier-downloaded)
code, run

``` bash
git fetch
git checkout release
```

Alternatively, you can checkout `master` branch

``` bash
git checkout master
```

to obtain a newer (not thorougly tested yet) version or `develop`\`

``` bash
git checkout develop
```

to obtain the bleeding edge version (feel free to try it, test it and
report the bugs).

### Using the package (without installing the pip/conda packages)

You can install the package from the obtained sources using

``` bash
pip install .
```

Or, if you want to develop ase2sprkkr, it is better idea to do an
[editable]{.title-ref} installation, where the package will see the
changes made.

``` bash
pip install --no-build-isolation --editable .
```

You can add `--no-deps` switch for a faster rebuild.

#### The known issues of the editable build

1\. Unfortunatelly, pip is not able to install the build requirements.
So if they are not installed, the previous command will fail. You can
either install the dependencies manually, copying them out from
pyproject.toml:

``` bash
pip install scikit-build cmake ninja cython
```

or (which will be usefull if I add a build requirement and forgot to
update this help) run the following snippet to install all the
build-dependencies mentioned in pyproject.toml:

``` bash
pip install tomllib
pip install --no-build-isolation $(python -c "import tomllib; print(' '.join(tomllib.load(open('pyproject.toml','rb'))['build-system']['requires']))")
```

2\. The limitation of the editable install is, that it won\'t see newly
created files automatically: you need run the command above again to
make it notice it.

3\. If the build process fail, try to remove the `build` directory
created by the previous build (if it exists). Mostly, it happens if the
`--no-build-isolation` switch is ommited.

## Documentation of the package

The documentation is published online at
<https://ase2sprkkr.github.io/ase2sprkkr/>

If you are using Git cloned sources, you can run

    open docs/index.html

to see the (offline version of the) documentation. The documentation
contains parts, that are generated from the docstrings in the source
code. You can regenerate these by

``` bash
make docs
```

Sphinx and sphinx-toolbox python packages (installable e.g. using pip)
and pandoc (for generating README.md) are needed for the build.

However, the official build of Sphinx miss some attributes when it is
used to build the documentation. So, till the pullrequest that corrects
the Sphinx behavior will be merged into Sphinx, please use the following
fork for building the documentation.
<https://github.com/lokik/sphinx.git>

## How to contribute or to report a bug

Please feel free to make a pull-request or post an issue at:
<https://github.com/ase2sprkkr/ase2sprkkr>
