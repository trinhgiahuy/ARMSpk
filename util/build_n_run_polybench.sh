#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/.."
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/intel.cfg
source $INTEL_PACKAGE intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
source $ROOTDIR/dep/spack/share/spack/setup-env.sh
spack load gcc; spack load llvm

export IACADIR=$ROOTDIR/dep/iaca-lin64
export IACAINCL=$IACADIR
export PATH=$IACADIR:$PATH

if [ ! -f ./polybench-c-3.2.tar.gz ]; then
	wget http://web.cse.ohio-state.edu/\~pouchet.2/software/polybench/download/polybench-c-3.2.tar.gz
fi
if [ ! -d ./polybench-c-3.2 ]; then
	tar xzf polybench-c-3.2.tar.gz
fi
cd ./polybench-c-3.2

BMs="linear-algebra/kernels/2mm/2mm.c
linear-algebra/kernels/3mm/3mm.c
linear-algebra/kernels/atax/atax.c
linear-algebra/kernels/bicg/bicg.c
linear-algebra/kernels/cholesky/cholesky.c
linear-algebra/kernels/doitgen/doitgen.c
linear-algebra/kernels/gemm/gemm.c
linear-algebra/kernels/gemver/gemver.c
linear-algebra/kernels/gesummv/gesummv.c
linear-algebra/kernels/mvt/mvt.c
linear-algebra/kernels/symm/symm.c
linear-algebra/kernels/syr2k/syr2k.c
linear-algebra/kernels/syrk/syrk.c
linear-algebra/kernels/trisolv/trisolv.c
linear-algebra/kernels/trmm/trmm.c
linear-algebra/solvers/durbin/durbin.c
linear-algebra/solvers/dynprog/dynprog.c
linear-algebra/solvers/gramschmidt/gramschmidt.c
linear-algebra/solvers/lu/lu.c
linear-algebra/solvers/ludcmp/ludcmp.c
datamining/correlation/correlation.c
datamining/covariance/covariance.c
medley/floyd-warshall/floyd-warshall.c
medley/reg_detect/reg_detect.c
stencils/adi/adi.c
stencils/fdtd-2d/fdtd-2d.c
stencils/fdtd-apml/fdtd-apml.c
stencils/jacobi-1d-imper/jacobi-1d-imper.c
stencils/jacobi-2d-imper/jacobi-2d-imper.c
stencils/seidel-2d/seidel-2d.c"

for BM in ${BMs}; do
	gcc -static -O3 -march=native \
		-I./utilities -I`dirname ${BM}` utilities/polybench.c ${BM} \
		-DLARGE_DATASET -o `basename ${BM}`.exe -L/usr/lib64 -lm
done

for BM in ${BMs}; do
	echo -ne "`basename ${BM}`\t" |tee -a real.log
	numactl -l -C 1 ./`basename ${BM}`.exe
	t_min="999999.0"
	for c in `seq 1 10`; do
		ts=`date +%s%N`
		numactl -l -C 1 ./`basename ${BM}`.exe
		te=`date +%s%N`
		t_curr=`echo "(${te}-${ts})/10^9" | bc -l`
		if (( $(echo "${t_curr} < ${t_min}" | bc -l) )); then t_min=${t_curr}; fi
	done
	echo "(${te}-${ts})/10^9" | bc -l |tee -a real.log
done

for BM in ${BMs}; do
	OUT="`basename ${BM}`_sde"; mkdir -p ${OUT}
	echo "SDE:" `basename ${BM}`
	../../dep/sde-external-8.35.0-2019-03-11-lin/sde64 \
		-sse-sde -disasm_att 1 -dcfg 1 -dcfg:write_bb 1 \
		-align_checker_prefetch 0 -align_correct 0 -emu_fast 1 \
		-bdw -- ./`basename ${BM}`.exe
	mv dcfg-out.* ${OUT}
done

for BM in ${BMs}; do
	D="`basename ${BM}`_sde"
	echo "get blocks:" `basename ${BM}`
	../parse_basic_blocks.py \
		-j ${D}/dcfg-out.dcfg.json.bz2 \
		-b ${D}/dcfg-out.bb.txt.bz2 \
		-s ${D}/dcfg-out.asm.b &
done
wait

for BM in ${BMs}; do
	D="`basename ${BM}`_sde"
	echo "analyze blocks:" `basename ${BM}`
	../parse_basic_blocks.py \
		-j ${D}/dcfg-out.dcfg.json.bz2 \
		-b ${D}/dcfg-out.bb.txt.bz2 \
		-l ${D}/dcfg-out.asm.b > ${D}/parser.log 2>&1 &
done
wait

for BM in ${BMs}; do
	echo `basename ${BM}` |tee -a sim.log
	/bin/grep 'Total CPU\|Converted to' `basename ${BM}`_sde/parser.log |tee -a sim.log
done

