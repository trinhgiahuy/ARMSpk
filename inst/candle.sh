#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR

source `cat $ROOTDIR/conf/intel.cfg` intel64
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F77=ifort
export I_MPI_F90=ifort
alias ar=`which xiar`
alias ld=`which xild`

BM="CANDLE"
VERSION="1dbcbd6c89d8655bcb5f721d1082c8d1d2ad3c3a"
if [ ! -f $ROOTDIR/dep/anaconda2/bin/anaconda ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	curl -o Anaconda2-5.1.0-Linux-x86_64.sh https://repo.continuum.io/archive/Anaconda2-5.1.0-Linux-x86_64.sh
	chmod u+x ./Anaconda2-5.1.0-Linux-x86_64.sh
	./Anaconda2-5.1.0-Linux-x86_64.sh -b -p $ROOTDIR/dep/anaconda2
	export PATH=$ROOTDIR/dep/anaconda2/bin:$PATH
	conda config --set changeps1 False
	conda install -y -c anaconda hdf5
	conda install -y -c anaconda theano
	conda install -y -c conda-forge keras
	conda install -y -c conda-forge opencv
	conda install -y -c conda-forge tqdm
	conda install -y -c intel --override-channels python=2 pip numpy
	conda install -y -c intel --override-channels intelpython2_core
	conda remove -y blaze
	conda install -y -c intel --override-channels tensorflow
	cd $ROOTDIR
fi

