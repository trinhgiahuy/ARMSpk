#!/bin/bash

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

source $ROOTDIR/dep/spack/share/spack/setup-env.sh
spack load openmpi@3.1.6%intel@19.0.1.144
export OMPI_CC=$I_MPI_CC
export OMPI_CXX=$I_MPI_CXX
export OMPI_F77=$I_MPI_F77
export OMPI_FC=$I_MPI_F90

BM="NGSAnalyzer"
VERSION="694b38eed8a4c09160045895a1bf86fcb35e85a3"
if [ ! -f $ROOTDIR/$BM/bin/workflow ]; then
	cd $ROOTDIR/$BM/
	git checkout -b precision ${VERSION}
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	# bwa w/ intel breaks
	sed -i -e 's/LDLIBS=-lm/LDLIBS=-static -static-intel -qopenmp-link=static -lm/' ./SNP_indel_caller/Makefile
	sed -i -e 's/-L${ADVISOR/-static -static-intel -qopenmp-link=static -L${ADVISOR/' ./workflow/Makefile
	sed -i -e 's/$(CC) $(CFL/icc $(CFL/' -e 's/-lm /-static -static-intel -qopenmp-link=static -lm /g' ./bwa-0.5.9rc1_kei/Makefile
	sed -i -e 's/ -lm / -static -static-intel -qopenmp-link=static -lm /g' -e 's/-lcurses/-lcurses -ltinfo/' ./samtools-0.1.8_kei/Makefile
	sed -i -e 's/$@ $(OBJS)/$@ $(OBJS) -static -static-intel -qopenmp-link=static/g' ./splitSam2Contig2/Makefile
	for x in `find ./workflow -name 'workflow*.sh'`; do sed -i -e 's/MPI_LOCALRANKID/PMIX_RANK/g' $x; done
	make -f makefile.x86_64_intel
	# we also need to get an input data set
	if [ ! -f $ROOTDIR/$BM/ngsa_mini_input/reference.fa ]; then
		echo "Creating $BM input may take 1h or more (stay tuned) ..."
		mkdir -p ./ngsa_mini_input/work/../bwa_db; cd ./ngsa_mini_input
		echo "  (downloading and processing reference genome)"
		bash $ROOTDIR/$BM/bin/download_reference.sh ./work
		../bin/samtools faidx reference.fa
		cd ./work
		../../bin/bwa index -a bwtsw -p ../bwa_db/reference.fa ../reference.fa
		bash $ROOTDIR/$BM/bin/download_contig.sh
		echo "  (downloading and processing pseudo-genome data)"
		cd $ROOTDIR/$BM/ngsa_mini_input
		wget http://mt.r-ccs.riken.jp/hpci-miniapp/ngsa-data/ngsa-dummy.tar.gz
		tar zxf ngsa-dummy.tar.gz
		echo "... done"
	fi
	cd $ROOTDIR
fi

