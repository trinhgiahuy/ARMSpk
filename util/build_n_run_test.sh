#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/.."
source $ROOTDIR/conf/host.cfg
source $ROOTDIR/conf/intel.cfg
source $INTEL_PACKAGE intel64 > /dev/null 2>&1
ulimit -s unlimited
ulimit -n 4096
MPIEXECOPT="-host `hostname` -genv I_MPI_FABRICS=shm:ofi -genv FI_PROVIDER=sockets"

if [ ! -e ../dep/sde-external-8.35.0-2019-03-11-lin/sde64 ]; then
	cd ../dep/
	tar xjf ~/sde-external-8.35.0-2019-03-11-lin.tar.bz2
	cd -
fi

set -x

TEST=sdetestC
OUT=${TEST}_out
mkdir -p ${OUT}
gcc ./${TEST}.c -o ${OUT}/${TEST}
../dep/sde-external-8.35.0-2019-03-11-lin/sde64 \
	-sse-sde -disasm_att 1 -dcfg 1 -dcfg:write_bb 1 \
	-align_checker_prefetch 0 -align_correct 0 -emu_fast 1 \
	-bdw -- \
	./${OUT}/${TEST}
mv dcfg-out.* ${OUT}

TEST=sdetestF
OUT=${TEST}_out
mkdir -p ${OUT}
gfortran ./${TEST}.f -o ${OUT}/${TEST}
../dep/sde-external-8.35.0-2019-03-11-lin/sde64 \
	-sse-sde -disasm_att 1 -dcfg 1 -dcfg:write_bb 1 \
	-align_checker_prefetch 0 -align_correct 0 -emu_fast 1 \
	-bdw -- \
	./${OUT}/${TEST}
mv dcfg-out.* ${OUT}

TEST=sdetestMPI
OUT=${TEST}_out
mkdir -p ${OUT}
mpiicc ./${TEST}.c -o ${OUT}/${TEST}
# mpiexec results will show #THREADS +2 (because mpi uses 2 additional?!)
# so here its no omp => 3 threads in result
mpiexec ${MPIEXECOPT} -np 4 bash -c "
../dep/sde-external-8.35.0-2019-03-11-lin/sde64 \
	-sse-sde -disasm_att 1 -dcfg 1 -dcfg:write_bb 1 \
	-dcfg:out_base_name dcfg-out.rank-\"\$MPI_LOCALRANKID\" \
	-align_checker_prefetch 0 -align_correct 0 -emu_fast 1 \
	-bdw -- \
	./${OUT}/${TEST}"
mv dcfg-out.* ${OUT}

TEST=sdetestOMPMPI
OUT=${TEST}_out
mkdir -p ${OUT}
mpiicc ./${TEST}.c -o ${OUT}/${TEST} -fopenmp
# mpiexec results will show #THREADS +2 (because mpi uses 2 additional?!)
# so here its 3 omp => 5 threads in result
mpiexec -genv OMP_NUM_THREADS=3 -genv OMP_SCHEDULE=static ${MPIEXECOPT} -np 4 bash -c "
../dep/sde-external-8.35.0-2019-03-11-lin/sde64 \
	-sse-sde -disasm_att 1 -dcfg 1 -dcfg:write_bb 1 \
	-dcfg:out_base_name dcfg-out.rank-\"\$MPI_LOCALRANKID\" \
	-align_checker_prefetch 0 -align_correct 0 -emu_fast 1 \
	-bdw -- \
	./${OUT}/${TEST}"
mv dcfg-out.* ${OUT}

