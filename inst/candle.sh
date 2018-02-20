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
	conda install -y -c intel python
	conda install -y -c intel intelpython2_core
	conda install -y -c intel tensorflow
	conda install -y -c intel mkl
	conda install -y -c anaconda hdf5
	conda install -y -c anaconda theano
	conda install -y -c conda-forge keras=2
	conda install -y -c conda-forge opencv
	conda install -y -c conda-forge tqdm
	cd $ROOTDIR
fi

