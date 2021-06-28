#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/inst/_common.sh
load_compiler_env "$1"

if [[ "$1" = *"fuji"* ]] || [[ "$1" = *"gem5"* ]] || [[ "$1" = *"llvm12"* ]]; then
	echo "lol, no thanks not touching this one"; exit 1
fi

BM="NGSAnalyzer"
VERSION="694b38eed8a4c09160045895a1bf86fcb35e85a3"
if [[ "$2" = *"rebuild"* ]]; then rm -rf $BM .git/modules/$BM; git submodule update --init $BM; fi
if [ ! -f $ROOTDIR/$BM/bin/workflow ]; then
	cd $ROOTDIR/$BM/
	if ! [[ "$(git rev-parse --abbrev-ref HEAD)" = *"precision"* ]]; then git checkout -b precision ${VERSION}; fi
	git apply --check $ROOTDIR/patches/*1-${BM}*.patch
	if [ "x$?" = "x0" ]; then git am --ignore-whitespace < $ROOTDIR/patches/*1-${BM}*.patch; fi
	instrument_kernel "$1" $ROOTDIR/$BM/
	# bwa w/ intel breaks
	if [[ "$1" = *"intel"* ]]; then
		sed -i -e 's/LDLIBS=-lm/LDLIBS=-static -static-intel -qopenmp-link=static -lm/' ./SNP_indel_caller/Makefile
		sed -i -e 's/-L${ADVISOR/-static -static-intel -qopenmp-link=static -L${ADVISOR/' ./workflow/Makefile
		sed -i -e 's/$(CC) $(CFL/icc $(CFL/' -e 's/-lm /-static -static-intel -qopenmp-link=static -lm /g' ./bwa-0.5.9rc1_kei/Makefile
		sed -i -e 's/ -lm / -static -static-intel -qopenmp-link=static -lm /g' -e 's/-lcurses/-lcurses -ltinfo/' ./samtools-0.1.8_kei/Makefile
		sed -i -e 's/$@ $(OBJS)/$@ $(OBJS) -static -static-intel -qopenmp-link=static/g' ./splitSam2Contig2/Makefile
	elif [[ "$1" = *"gnu"* ]]; then
		sed -i -e 's/=icc/=gcc/' -e 's/-ipo -xHost/-march=native -flto/g' -e "s/LDLIBS=-lm/LDLIBS=-flto ${MAYBESTATIC} -lm/" ./SNP_indel_caller/Makefile
		sed -i -e 's/=icc/=gcc/' -e 's/=icpc/=g++/' -e 's/-ipo -xHost/-march=native/g' ./makefile.x86_64_intel
		sed -i -e 's/=icc/=gcc/' -e 's/=icpc/=g++/' -e 's/-ipo -xHost/-march=native -m64/g' -e 's/-lcurses/-lcurses -ltinfo/' ./samtools-0.1.8_kei/Makefile
		sed -i -e 's/=icc/=gcc/' -e 's/=icpc/=g++/' -e 's/-ipo -xHost/-march=native -m64/g' -e 's/-lcurses/-lcurses -ltinfo/' ./samtools-0.1.8_kei/misc/Makefile
		sed -i -e "s/ -lm / -flto ${MAYBESTATIC} -lm /g" ./samtools-0.1.8_kei/Makefile
		sed -i -e "s/\$@ \$(OBJS)/\$@ \$(OBJS) -flto ${MAYBESTATIC}/g" ./splitSam2Contig2/Makefile
		sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -flto ${MAYBESTATIC}#g" ./workflow/Makefile
		sed -i -e "s/-lm /-flto ${MAYBESTATIC} -lm /g" ./bwa-0.5.9rc1_kei/Makefile
		for FILE in `/usr/bin/grep 'inline void bwt_2occ' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/inline void bwt_2occ/void bwt_2occ/' $FILE; done
		for FILE in `/usr/bin/grep 'inline void bwtl_2occ' -r | cut -d':' -f1 | sort -u`; do sed -i -e 's/inline void bwtl_2occ/void bwtl_2occ/' $FILE; done
	fi
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
		URL="http://mt.r-ccs.riken.jp/hpci-miniapp/ngsa-data/ngsa-dummy.tar.gz"; DEP=$(basename $URL)
		if [ ! -f $ROOTDIR/dep/${DEP} ]; then if ! wget ${URL} -O $ROOTDIR/dep/${DEP}; then echo "ERR: download failed for ${URL}"; exit 1; fi; fi; tar xzf $ROOTDIR/dep/${DEP}
		echo "... done"
	fi
	cd $ROOTDIR
fi

