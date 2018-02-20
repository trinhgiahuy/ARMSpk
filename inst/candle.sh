#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F77=ifort
export I_MPI_F90=ifort

# compile CANDLE
if [ ! -f $ROOTDIR/dep/anaconda2/bin/anaconda ]; then
	cd $ROOTDIR/CANDLE/
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

