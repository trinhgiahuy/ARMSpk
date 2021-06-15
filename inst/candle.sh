#!/bin/bash
exit 1 #ignore in this study

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/intel.cfg
source $INTEL_PACKAGE intel64 > /dev/null 2>&1
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F77=ifort
export I_MPI_F90=ifort
alias ar=`which xiar`
alias ld=`which xild`
export ADVISOR_2018_DIR=${ADVISOR_2019_DIR}

BM="CANDLE"
VERSION="ea14ed86d3e612f56383c56a6cff6f77210f7412"
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
	conda config --add channels intel
	conda create -y -n idp intelpython2_core python=2
	source activate idp
	conda install -y -c intel hdf5=1.10.2 theano=1.0.2 keras=2.2.4 opencv=3.4.1 tqdm=4.32.1 pip=9.0.3 numpy=1.16.2 tensorflow=1.13.1 pandas=0.24.1 scikit-learn=0.20.3
	cd $ROOTDIR
fi

