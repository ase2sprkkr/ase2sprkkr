Installation of the package using package managers
==================================================

The simplest way how to install and use the package is to install them
using package managers: either pip

.. code:: bash

   pip install ase2sprkkr

or conda

.. code:: bash

   conda install -c ase2sprkkr ase2sprkkr

Further notes
--------------

In some systems, the ``pip`` utility for ``python3`` is called ``pip3``.
If it is not installed, you can install it using the linux distribution package
manager, e.g. in Debian/Ubuntu

.. code:: bash

   apt install pip3

or

.. code:: bash

   zypper install pip

in OpenSUSE

For the conda installation instructions, see the Anaconda documentation
https://docs.anaconda.com/anaconda/install/linux/ however, for the users
unexperienced with conda, the (simpler) pip way is recommended.

Install the packages from GIT (and/or source codes)
===================================================

If you do not want to use public package managers as pip or conda, or
you want to contribute to development, you can use GIT to obtain the
package sources.

Requirements
------------

-  Python >= 3.7
-  SPR-KKR (not checked by the installer)
-  Python packages: see the the setup.cfg
-  Git (to obtain the sources)

Obtaining the package using GIT
-------------------------------

.. code:: bash

   git clone https://github.com/ase2sprkkr/ase2sprkkr.git
   git checkout release

The first line fetches the code of the package. The second one checks
out the recommended production version of the code.

If you want to obtain the current version of the (earlier-downloaded)
code, run

.. code:: bash

   git fetch
   git checkout release

Alternatively, you can checkout ``master`` or ``develop`` branch

.. code:: bash

   git checkout master

to obtain a newer (not thorougly tested yet) version or

.. code:: bash

   git checkout develop

to obtain the bleeding edge version (feel free to try it, test it
and report the bugs).

Using the package (without installing the pip/conda packages)
-------------------------------------------------------------

You can either just add the src directory to your PYTHONPATH, or you can
build and install the package, as it is described below.

Installation of the package from the sources
--------------------------------------------

To install the package, the simplest way is to use pip

.. code:: bash

   python3 -m pip install .

Maybe, you will have to replace ``python3`` with ``python``.
For an editable install, please run

.. code:: bash

   python3 setup.py develop --user

and ignore some deprecation warning. Editable install is aimed for developers:
in this type of install, only link to the current directory will be added to
your local ``site-packages``, which allows you to use the changesyou make to
the source code.


Documentation of the package
============================

The documentation is published online at
`https://ase2sprkkr.github.io/ase2sprkkr/ <https://ase2sprkkr.github.io/ase2sprkkr/>`__

If you are using Git cloned sources, you can run

::

   open docs/index.html

to see the (offline version of the) documentation. The documentation
contains parts, that are generated from the docstrings in the source
code. You can regenerate these by

.. code:: bash

   make docs

Sphinx and sphinx-toolbox python packages (installable e.g. using pip) and pandoc
(for generating README.md) are needed for the build.

However, the official build of Sphinx miss some attributes when it is used
to build the documentation. So, till the pullrequest that corrects the Sphinx
behavior will be merged into Sphinx, please use the following fork for
building the documentation.
`https://github.com/lokik/sphinx.git <https://github.com/lokik/sphinx.git>`__


How to contribute or to report a bug
====================================

Please feel free to make a pull-request or post an issue at:
`https://github.com/ase2sprkkr/ase2sprkkr <https://github.com/ase2sprkkr/ase2sprkkr>`__
