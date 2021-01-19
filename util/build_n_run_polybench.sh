#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/.."
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/intel.cfg
source $INTEL_PACKAGE intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
alias ar=`which xiar`
alias ld=`which xild`
export ADVISOR_2018_DIR=${ADVISOR_2019_DIR}

if [ ! -f ./polybench-c-4.2.1-beta.tar.gz ]; then
	wget https://downloads.sourceforge.net/project/polybench/polybench-c-4.2.1-beta.tar.gz
fi
if [ ! -d ./polybench-c-4.2.1-beta ]; then
	tar xzf ./polybench-c-4.2.1-beta.tar.gz
fi
cd ./polybench-c-4.2.1-beta

BMs="datamining/correlation/correlation.c|LARGE
datamining/covariance/covariance.c|LARGE
linear-algebra/blas/gemm/gemm.c|LARGE
linear-algebra/blas/gemver/gemver.c|LARGE
linear-algebra/blas/gesummv/gesummv.c|LARGE
linear-algebra/blas/symm/symm.c|LARGE
linear-algebra/blas/syr2k/syr2k.c|LARGE
linear-algebra/blas/syrk/syrk.c|LARGE
linear-algebra/blas/trmm/trmm.c|LARGE
linear-algebra/kernels/2mm/2mm.c|LARGE
linear-algebra/kernels/3mm/3mm.c|LARGE
linear-algebra/kernels/atax/atax.c|LARGE
linear-algebra/kernels/bicg/bicg.c|LARGE
linear-algebra/kernels/doitgen/doitgen.c|LARGE
linear-algebra/kernels/mvt/mvt.c|LARGE
linear-algebra/solvers/cholesky/cholesky.c|LARGE
linear-algebra/solvers/durbin/durbin.c|LARGE
linear-algebra/solvers/gramschmidt/gramschmidt.c|LARGE
linear-algebra/solvers/lu/lu.c|LARGE
linear-algebra/solvers/ludcmp/ludcmp.c|LARGE
linear-algebra/solvers/trisolv/trisolv.c|LARGE
medley/deriche/deriche.c|LARGE
medley/floyd-warshall/floyd-warshall.c|MEDIUM
medley/nussinov/nussinov.c|LARGE
stencils/adi/adi.c|LARGE
stencils/fdtd-2d/fdtd-2d.c|LARGE
stencils/heat-3d/heat-3d.c|LARGE
stencils/jacobi-1d/jacobi-1d.c|LARGE
stencils/jacobi-2d/jacobi-2d.c|LARGE
stencils/seidel-2d/seidel-2d.c|MEDIUM"

## patch timing and SSC marker
#if patch --dry-run -s -f -p1 < ../patches/0001-polybench-c-4.2.1-beta.patch; then
#	patch -p1 < ../patches/0001-polybench-c-4.2.1-beta.patch
#fi
#
#for BMconf in ${BMs}; do
#	BM="`echo ${BMconf} | cut -d '|' -f1`"
#	DS="`echo ${BMconf} | cut -d '|' -f2`"
#	rm -f `basename ${BM}`.exe
#	icc -O3 -xHost -static -static-intel -I${ADVISOR_2018_DIR}/include \
#		-I./utilities -I`dirname ${BM}` utilities/polybench.c ${BM} \
#		-D${DS}_DATASET -DPOLYBENCH_TIME \
#		-o `basename ${BM}`.exe -L${ADVISOR_2018_DIR}/lib64 -littnotify -L/usr/lib64 -lm
#done
#
#for BMconf in ${BMs}; do
#	BM="`echo ${BMconf} | cut -d '|' -f1`"
#	echo -ne "`basename ${BM} | cut -d '.' -f1` " |tee -a real.log
#	numactl -l -C 1 ./`basename ${BM}`.exe > /dev/null
#	t_min="999999.0"
#	for c in `seq 1 10`; do
#		t_curr=`numactl -l -C 1 ./\`basename ${BM}\`.exe | awk -F 'kernel: ' '{print \$2}' | cut -d ' ' -f1`
#		if (( $(echo "${t_curr} < ${t_min}" | bc -l) )); then t_min=${t_curr}; fi
#	done
#	echo "${t_min}" |tee -a real.log
#done
#
#for BMconf in ${BMs}; do
#	BM="`echo ${BMconf} | cut -d '|' -f1`"
#	OUT="`basename ${BM}`_sde"; mkdir -p ${OUT}
#	echo "SDE:" `basename ${BM}`
#	../../dep/sde-external-8.35.0-2019-03-11-lin/sde64 \
#		-sse-sde -disasm_att 1 -dcfg 1 -dcfg:write_bb 1 \
#		-align_checker_prefetch 0 -align_correct 0 -emu_fast 1 \
#		-start_ssc_mark 111:repeat -stop_ssc_mark 222:repeat \
#		-bdw  -- ./`basename ${BM}`.exe
#	mv dcfg-out.* ${OUT}
#done
#
#for BMconf in ${BMs}; do
#	BM="`echo ${BMconf} | cut -d '|' -f1`"
#	D="`basename ${BM}`_sde"
#	echo "get blocks:" `basename ${BM}`
#	../parse_basic_blocks.py \
#		-j ${D}/dcfg-out.dcfg.json.bz2 \
#		-b ${D}/dcfg-out.bb.txt.bz2 \
#		-s ${D}/dcfg-out.asm.b &
#done
#wait
#
source $ROOTDIR/dep/spack/share/spack/setup-env.sh
spack load llvm

export IACADIR=$ROOTDIR/dep/iaca-lin64
export IACAINCL=$IACADIR
export PATH=$IACADIR:$PATH

#for BMconf in ${BMs}; do
#	BM="`echo ${BMconf} | cut -d '|' -f1`"
#	D="`basename ${BM}`_sde"
#	echo "analyze blocks:" `basename ${BM}`
#	../parse_basic_blocks.py \
#		-j ${D}/dcfg-out.dcfg.json.bz2 \
#		-b ${D}/dcfg-out.bb.txt.bz2 \
#		-l ${D}/dcfg-out.asm.b > ${D}/parser.log 2>&1 &
#done
#wait
#
#for BMconf in ${BMs}; do
#	BM="`echo ${BMconf} | cut -d '|' -f1`"
#	echo `basename ${BM}` |tee -a sim.log
#	/bin/grep 'Total CPU\|Converted to' `basename ${BM}`_sde/parser.log |tee -a sim.log
#done
for BMconf in ${BMs}; do
	BM="`echo ${BMconf} | cut -d '|' -f1`"
	echo -ne "`basename ${BM}` "
	echo -ne "`/bin/grep 'LLVM: Total CPU' \`basename ${BM}\`_sde/parser.log | awk -F 'ID 0 : ' '{print $2}'` "
	echo -ne "`/bin/grep 'IACA: Total CPU' \`basename ${BM}\`_sde/parser.log | awk -F 'ID 0 : ' '{print $2}'` "
	echo     "`/bin/grep 'OSACA: Total CPU' \`basename ${BM}\`_sde/parser.log | awk -F 'ID 0 : ' '{print $2}'`"
done
