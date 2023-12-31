#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )/../" && pwd )"
cd $ROOTDIR
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/env.cfg
load_compiler_env "$1"

BM="HPL"
if [[ "$2" = *"rebuild"* ]]; then rm -rf $BM .git/modules/$BM; git submodule update --init $BM; fi
if [ ! -f $ROOTDIR/$BM/bin/Linux_Intel64/xhpl ]; then
	mkdir -p $ROOTDIR/$BM/
	cd $ROOTDIR/$BM/
	URL="http://www.netlib.org/benchmark/hpl/hpl-2.2.tar.gz"; DEP=$(basename $URL)
	if [ ! -f $ROOTDIR/dep/${DEP} ]; then if ! wget ${URL} -O $ROOTDIR/dep/${DEP}; then echo "ERR: download failed for ${URL}"; exit 1; fi; fi
	tar xzf $ROOTDIR/dep/${DEP} -C $ROOTDIR/$BM --strip-components 1
	patch -p1 --forward < $ROOTDIR/patches/*1-${BM}*.patch
	instrument_kernel "$1" $ROOTDIR/$BM/
	if [[ "$1" = *"intel"* ]]; then
		sed -i -e 's/mpiicc/mpicc/' -e 's/ -L${ADVISOR/-static -static-intel -qopenmp-link=static -L${ADVISOR/' -e "s#\$(HOME)/hpl#$ROOTDIR/$BM#g" ./Make.Linux_Intel64
	elif [[ "$1" = *"gnu"* ]]; then
		if [ -n "$FJMPI" ]; then sed -i -e 's/mpiicc/mpifcc/' -e "s/-ipo -xHost/-march=native -fno-lto ${MAYBESTATIC}/g" ./Make.Linux_Intel64;
		else                     sed -i -e 's/mpiicc/mpicc/' -e "s/-ipo -xHost/-march=native -m64 -flto ${MAYBESTATIC}/g" ./Make.Linux_Intel64; fi
		if [ -n "$MKLROOT" ]; then
			sed -i -e 's# -I${ADVISOR_2018_DIR}/include##g' -e 's# -L${ADVISOR_2018_DIR}/lib64 -littnotify##g' -e 's/libmkl_intel_thread.a/libmkl_gnu_thread.a/g' -e 's/-lpthread/-lgomp -lpthread -lm/g' -e "s#\$(HOME)/hpl#$ROOTDIR/$BM#g" -e 's/ -ansi-alias -i-static//g' -e 's/ -nocompchk//g' -e 's/ -mt_mpi//g' ./Make.Linux_Intel64
		elif [ -n "$FJBLAS" ]; then
			sed -i -e "s# -I\${ADVISOR_2018_DIR}/include# -I$(dirname `which fcc`)/../include#g" -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -L$(readlink -f $(dirname $(which mpifcc))/../lib64)#g" -e "s#\$(HOME)/hpl#$ROOTDIR/$BM#g" -e 's/ -ansi-alias -i-static//g' -e 's/ -nocompchk//g' -e 's/ -mt_mpi//g' -e 's/ -I$(LAinc)//g' -e "s# \$(LAlib)# $(readlink -f $(dirname $(which mpifcc))/../lib64)/fjhpctag.o $(readlink -f $(dirname $(which mpifcc))/../lib64)/fjlang08.o $(readlink -f $(dirname $(which mpifcc))/../lib64)/fjomp.o -lfj90rt2 -lssl2mtexsve -lssl2mtsve -lfj90i -lfj90fmt_sve -lfj90f -lfjsrcinfo -lfj90rt -lfjprofcore -lfjprofomp -lelf -lm#g" ./Make.Linux_Intel64
		fi
	elif [[ "$1" = *"fujitrad"* ]]; then
		sed -i -e 's/mpiicc/mpifcc/' -e 's/-ipo -xHost/-Kfast,openmp,ocl,largepage/g' -e 's/ -z noexecstack -z relro -z now//g' -e 's/ -Wall//g' -e "s# -I\${ADVISOR_2018_DIR}/include# -I$(readlink -f $(dirname $(which mpifcc))/../include/mpi/fujitsu)#g" -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify##g" -e "s#\$(HOME)/hpl#$ROOTDIR/$BM#g" -e 's/ -ansi-alias -i-static//g' -e 's/ -nocompchk//g' -e 's/ -mt_mpi//g' -e 's/ -I$(LAinc)//g' -e "s# \$(LAlib)# -SSL2BLAMP#g" ./Make.Linux_Intel64
	elif [[ "$1" = *"fujiclang"* ]]; then
		sed -i -e 's/mpiicc/mpifcc/' -e 's/-ipo -xHost/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-largepage/g' -e 's/ -z noexecstack -z relro -z now//g' -e 's/ -Wall//g' -e "s# -I\${ADVISOR_2018_DIR}/include# -I$(readlink -f $(dirname $(which mpifcc))/../include/mpi/fujitsu)#g" -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify##g" -e "s#\$(HOME)/hpl#$ROOTDIR/$BM#g" -e 's/ -ansi-alias -i-static//g' -e 's/ -nocompchk//g' -e 's/ -mt_mpi//g' -e 's/ -I$(LAinc)//g' -e "s# \$(LAlib)# -SSL2BLAMP#g" ./Make.Linux_Intel64
	elif [[ "$1" = *"gem5"* ]]; then
		sed -i -e 's/mpiicc/fcc/' -e 's/-ipo -xHost/-Nclang -Ofast -mcpu=a64fx+sve -fopenmp -ffj-ocl -ffj-no-largepage -fno-lto/g' -e 's/ -z noexecstack -z relro -z now//g' -e 's/ -Wall//g' -e "s# -I\${ADVISOR_2018_DIR}/include# -I$ROOTDIR/dep/mpistub/include/mpistub#g" -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -Wl,-rpath=$ROOTDIR/dep/mpistub/lib/mpistub -L$ROOTDIR/dep/mpistub/lib/mpistub -lmpi#g" -e "s#\$(HOME)/hpl#$ROOTDIR/$BM#g" -e 's/ -ansi-alias -i-static//g' -e 's/ -nocompchk//g' -e 's/ -mt_mpi//g' -e 's/ -I$(LAinc)//g' -e "s# \$(LAlib)# -SSL2BLAMP#g" ./Make.Linux_Intel64
	elif [[ "$1" = *"llvm12"* ]]; then
		sed -i -e 's/mpiicc/mpifcc/' -e 's/-ipo -xHost/-Ofast -ffast-math -mcpu=a64fx -mtune=a64fx -fopenmp -mllvm -polly -mllvm -polly-vectorizer=polly -flto=full/g' -e 's# -I${ADVISOR_2018_DIR}/include##g' -e "s# -L\${ADVISOR_2018_DIR}/lib64 -littnotify# -fuse-ld=lld -L$(readlink -f $(dirname $(which mpifcc))/../lib64) -Wl,-rpath=$(readlink -f $(dirname $(which clang))/../lib)#g" -e "s#\$(HOME)/hpl#$ROOTDIR/$BM#g" -e 's/ -ansi-alias -i-static//g' -e 's/ -nocompchk//g' -e 's/ -mt_mpi//g' -e 's/ -I$(LAinc)//g' -e "s# \$(LAlib)# $(readlink -f $(dirname $(which mpifcc))/../lib64)/fjhpctag.o $(readlink -f $(dirname $(which mpifcc))/../lib64)/fjlang08.o $(readlink -f $(dirname $(which mpifcc))/../lib64)/fjomp.o -lfj90rt2 -lssl2mtexsve -lssl2mtsve -lfj90i -lfj90fmt_sve -lfj90f -lfjsrcinfo -lfj90rt -lfjprofcore -lfjprofomp -lm#g" ./Make.Linux_Intel64
	fi
	make arch=Linux_Intel64
	cd $ROOTDIR
fi

