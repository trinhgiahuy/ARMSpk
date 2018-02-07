#!/bin/bash

cd "$( dirname "${BASH_SOURCE[0]}" )"

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F77=ifort
export I_MPI_F90=ifort
ulimit -s unlimited
ulimit -n 4096

# compile AMG -> comes w/ 2 problems
if [ ! -f ./AMG/test/amg ]; then
	cd ./AMG/
	make
	make clean
	cd ../
fi

# compile CANDLE -> comes w/ 7 problems
if [ ! -f $HOME/anaconda2/bin/anaconda ]; then
	cd ./CANDLE/
	curl -o Anaconda2-4.3.1-Linux-x86_64.sh https://repo.continuum.io/archive/Anaconda2-4.3.1-Linux-x86_64.sh
	chmod u+x ./Anaconda2-4.3.1-Linux-x86_64.sh
	./Anaconda2-4.3.1-Linux-x86_64.sh -b
	export PATH=$HOME/anaconda2/bin:$PATH
	conda install -y -c conda-forge tensorflow
	conda install -y -c anaconda hdf5=1.8.17
	conda install -y -c anaconda theano
	conda install -y -c conda-forge keras=2
	conda install -y -c conda-forge opencv
	conda install -y -c conda-forge tqdm
	conda update -y -c conda-forge numpy
	cd ../
fi


