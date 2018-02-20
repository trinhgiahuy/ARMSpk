#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR

source /opt/intel/parallel_studio_xe_2018.1.038/bin/psxevars.sh intel64
export I_MPI_CC=icc
export I_MPI_CXX=icpc
export I_MPI_F77=ifort
export I_MPI_F90=ifort

# compile NGS Analyzer
if [ ! -f $ROOTDIR/NGSAnalyzer/bin/workflow ]; then
	cd $ROOTDIR/NGSAnalyzer
	sed -i -e 's/^N_THREADS=1/N_THREADS=$OMP_NUM_THREADS/' ./workflow/workflow_01.sh
	sed -i -e 's/CC.*=.*gcc/CC=icc/g' -e 's/CXX.*=.*g++/CXX=icpc/g' ./makefile.x86_64_gcc
	sed -i -e 's/CC.*=.*gcc/CC=icc/g' -e 's/CXX.*=.*g++/CXX=icpc/g' ./SNP_indel_caller/Makefile
	# bwa w/ intel breaks
	#sed -i -e 's/CC.*=.*gcc/CC=icc/g' -e 's/CXX.*=.*g++/CXX=icpc/g' ./bwa-0.5.9rc1_kei/bwt_gen/Makefile
	sed -i -e 's/^BWA_OPT.*= CC=$(CC) CXX=$(CXX)/BWA_OPT = CC=gcc CXX=g++/g' ./makefile.x86_64_gcc
	sed -i -e 's/CC.*=.*gcc/CC=icc/g' -e 's/CXX.*=.*g++/CXX=icpc/g' ./samtools-0.1.8_kei/examples/Makefile
	sed -i -e 's/CC.*=.*gcc/CC=icc/g' -e 's/CXX.*=.*g++/CXX=icpc/g' ./samtools-0.1.8_kei/misc/Makefile
	make -f makefile.x86_64_gcc
	# we also need to get an input data set
	if [ ! -f $ROOTDIR/NGSAnalyzer/ngsa_mini_input/reference.fa ]; then
		echo "Creating NGSAnalyzer input may take 1h or more (stay tuned) ..."
		mkdir -p ./ngsa_mini_input/work/../bwa_db; cd ./ngsa_mini_input
		echo "  (downloading and processing reference genome)"
		bash $ROOTDIR/NGSAnalyzer/bin/download_reference.sh ./work
		../bin/samtools faidx reference.fa
		cd ./work
		../../bin/bwa index -a bwtsw -p ../bwa_db/reference.fa ../reference.fa
		bash $ROOTDIR/NGSAnalyzer/bin/download_contig.sh
		echo "  (downloading and processing pseudo-genome data)"
		cd $ROOTDIR/NGSAnalyzer/ngsa_mini_input
		wget http://mt.aics.riken.jp/hpci-miniapp/ngsa-data/ngsa-dummy.tar.gz
		tar zxf ngsa-dummy.tar.gz
		echo "... done"
	fi
	cd $ROOTDIR
fi

