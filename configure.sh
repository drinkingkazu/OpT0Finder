#!/bin/bash

# clean up previously set env
if [[ -z $FORCE_FMATCH_BASEDIR ]]; then
    where="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    export FMATCH_BASEDIR=${where}
else
    export FMATCH_BASEDIR=$FORCE_FMATCH_BASEDIR
fi

# set the build dir
unset FMATCH_BUILDDIR
if [[ -z $FMATCH_BUILDDIR ]]; then
    export FMATCH_BUILDDIR=$FMATCH_BASEDIR/build
fi

# Check python version compatibility:
export FMATCH_PYTHON_CONFIG=python-config
FMATCH_PYVERSION=0
export FMATCH_PYTHON=`which python`
if [ `command -v python` ]; then
    FMATCH_PYVERSION=$($FMATCH_PYTHON -c "import sys; print(sys.version_info.major)")
else
    export FMATCH_PYTHON=`which python3`
    FMATCH_PYVERSION=$($FMATCH_PYTHON -c "import sys; print(sys.version_info.major)")
fi
if [[ $FMATCH_PYVERSION -gt 2 ]]
then
    minor=$(python3 -c "import sys; print(sys.version_info.minor)")
    export FMATCH_PYTHON_CONFIG=python${FMATCH_PYVERSION}.${minor}-config
fi

export FMATCH_DIR=$FMATCH_BASEDIR/flashmatch/
export FMATCH_LIBDIR=$FMATCH_BUILDDIR/lib
export FMATCH_INCDIR=$FMATCH_BUILDDIR/include
export FMATCH_BINDIR=$FMATCH_BUILDDIR/bin
export FMATCH_DATADIR=$MATCH_BASEDIR/dat
export FMATCH_INCLUDES="-I${FMATCH_INCDIR} `${FMATCH_PYTHON_CONFIG} --includes` "
export FMATCH_LIBS="-L${FMATCH_LIBDIR} -lflashmatch "

# Abort if ROOT not installed. Let's check rootcint for this.
if [[ -z `command -v rootcint` ]]; then
    echo
    echo Looks like you do not have ROOT installed.
    echo You cannot use LArLite w/o ROOT!
    echo Aborting.
    echo
    return 1;
fi

# Check Numpy
export FMATCH_NUMPY=`$FMATCH_PYTHON $FMATCH_BASEDIR/bin/check_numpy`

# warning for missing support
missing=""
if [ $FMATCH_NUMPY -eq 0 ]; then
    missing+=" Numpy"
else
    FMATCH_INCLUDES="${FMATCH_INCLUDES} -I`$FMATCH_PYTHON -c\"import numpy; print(numpy.get_include())\"`"
    FMATCH_LIBS="-L`${FMATCH_PYTHON_CONFIG} --prefix`/lib/ `${FMATCH_PYTHON_CONFIG} --ldflags` ${FMATCH_LIBS}"
fi
if [[ $missing ]]; then
    printf "\033[93mWarning\033[00m ... missing$missing support. Build without them.\n";
fi

echo
printf "\033[93mFMATCH\033[00m FYI shell env. may useful for external packages:\n"
printf "    \033[95mFMATCH_INCDIR\033[00m   = $FMATCH_INCDIR\n"
printf "    \033[95mFMATCH_LIBDIR\033[00m   = $FMATCH_LIBDIR\n"
printf "    \033[95mFMATCH_BUILDDIR\033[00m = $FMATCH_BUILDDIR\n"

export PATH=$FMATCH_BASEDIR/bin:$FMATCH_BINDIR:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$FMATCH_LIBDIR
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$FMATCH_LIBDIR

mkdir -p $FMATCH_BUILDDIR;
mkdir -p $FMATCH_LIBDIR;
mkdir -p $FMATCH_BINDIR;

export LD_LIBRARY_PATH=$FMATCH_LIBDIR:$LD_LIBRARY_PATH
export PYTHONPATH=$FMATCH_BASEDIR/python:$PYTHONPATH

export FMATCH_CXX=clang++
if [ -z `command -v $FMATCH_CXX` ]; then
    export FMATCH_CXX=g++
    if [ -z `command -v $FMATCH_CXX` ]; then
        echo
        echo Looks like you do not have neither clang or g++!
        echo
        return 1;
    fi
fi

echo
echo "Finish configuration. To build, type:"
echo "> make "
echo

