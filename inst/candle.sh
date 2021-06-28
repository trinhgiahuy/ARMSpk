#!/bin/bash
exit 1 #ignore in this study

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/inst/_common.sh
load_compiler_env "$1"

BM="CANDLE"
VERSION="ea14ed86d3e612f56383c56a6cff6f77210f7412"
if [[ "$2" = *"rebuild"* ]]; then rm -rf $BM .git/modules/$BM; git submodule update --init $BM; fi
if [ ! -f $ROOTDIR/dep/anaconda2/bin/anaconda ]; then
	cd $ROOTDIR/$BM/
	if ! [[ "$(git rev-parse --abbrev-ref HEAD)" = *"precision"* ]]; then git checkout -b precision ${VERSION}; fi
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

