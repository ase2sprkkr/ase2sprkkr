#!/bin/bash
CURRENT=
FORCE=

while getopts "cf" opt; do
  case ${opt} in
		f )  FORCE="--force-reinstall --no-deps";;
    c )  CURRENT=1 ;;
    \? ) echo "Install the ase2sprkkr package. Use -c to install the current version of the
		           code (do not the git checkout release command)"
				 exit
      ;;
  esac
done


cd "$(dirname "$0")"
if [[ -z "$CURRENT" ]] ; then
	echo "Checking out the current version of the code"
	if ! git fetch ; then
		echo "Warning - either git is not installed, or the package source
		is not a git repository, the origin of the repository is not available.
		Attemp to get the latest version of the package failed!!!"
	fi
	if ! git checkout origin/release ; then
		echo "Warning - either git is not installed, or the package source
		is not a git repository. I will install the current (possibly development?)
		version of the code!!!"
	fi
fi

#Guessing the python executable name
if [[ -z "$PYTHON" ]] ; then
	if [[ -z "`whereis python3`" ]] ; then
		PYTHON=python
	else
		PYTHON=python3
	fi
fi

if ! echo "import build" | $PYTHON  2> /dev/null ; then
	echo "Install the build package for python to build the package."
	$PYTHON -m pip install --upgrade build
fi

echo "Building the wheel package"
$PYTHON -m build

echo "Installing the package"
$PYTHON -m pip install $FORCE `ls ./dist/ase2sprkkr-*.whl | sort | tail -n 1`
