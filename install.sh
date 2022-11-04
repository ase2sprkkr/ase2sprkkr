#!/bin/bash
CURRENT=1
FORCE=1

while getopts "cf" opt; do
  case ${opt} in
		n )  FORCE= ;;
    c )  CURRENT= ;;
    \? ) echo "Install the ase2sprkkr package. Parameters:
							 -c  check out the release version before installing
							 -n  install, only if it is not already installed "
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
	COMMIT="`git symbolic-ref --short HEAD 2> /dev/null ||  git rev-parse HEAD`"
	echo "Remebering your commit: $COMMIT, checking out origin/release"
	if ! git checkout origin/release ; then
		echo "Warning - either git is not installed, the package source
is not a git repository or the repository is "dirty" and thus it
is not possible to checkout. I will install the current (possibly development?)
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
$PYTHON -m pip install `ls ./dist/ase2sprkkr-*.whl | sort | tail -n 1`
if [[ $FORCE ]] ; then
	$PYTHON -m pip install --force-reinstall --no-deps `ls ./dist/ase2sprkkr-*.whl | sort | tail -n 1`
fi

if [[ -z "$CURRENT" ]] ; then
	 if [[ "`git rev-parse $COMMIT`" != "`git rev-parse HEAD`" ]] ; then
		  echo "Version from the online release have been installed!!!
Use -c parameter to install the current version"
	 fi
	 if [[ $COMMIT != "`git rev-parse HEAD`" ]] ; then
			echo "Restoring the original commit $COMMIT"
			git checkout $COMMIT
	 fi
fi
