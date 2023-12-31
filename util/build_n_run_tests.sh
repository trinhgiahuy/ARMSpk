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
icc ./${TEST}.c -o ${OUT}/${TEST}
../dep/sde-external-8.35.0-2019-03-11-lin/sde64 \
	-sse-sde -disasm_att 1 -dcfg 1 -dcfg:write_bb 1 \
	-align_checker_prefetch 0 -align_correct 0 -emu_fast 1 \
	-bdw -- \
	./${OUT}/${TEST}
mv dcfg-out.* ${OUT}

TEST=sdetest2C
OUT=${TEST}_out_L1
mkdir -p ${OUT}
icpc -xHOST -O3 -DMIN_SIZE=$((8*1024)) -DMAX_SIZE=$((16*1024)) ./${TEST}.cpp -o ${OUT}/${TEST}
../dep/sde-external-8.35.0-2019-03-11-lin/sde64 \
	-sse-sde -disasm_att 1 -dcfg 1 -dcfg:write_bb 1 \
	-align_checker_prefetch 0 -align_correct 0 -emu_fast 1 \
	-bdw -- \
	./${OUT}/${TEST}
mv dcfg-out.* ${OUT}
OUT=${TEST}_out_L2
mkdir -p ${OUT}
icpc -xHOST -O3 -DMIN_SIZE=$((128*1024)) -DMAX_SIZE=$((256*1024)) ./${TEST}.cpp -o ${OUT}/${TEST}
../dep/sde-external-8.35.0-2019-03-11-lin/sde64 \
	-sse-sde -disasm_att 1 -dcfg 1 -dcfg:write_bb 1 \
	-align_checker_prefetch 0 -align_correct 0 -emu_fast 1 \
	-bdw -- \
	./${OUT}/${TEST}
mv dcfg-out.* ${OUT}
OUT=${TEST}_out_L3
mkdir -p ${OUT}
icpc -xHOST -O3 -DMIN_SIZE=$((2048*1024)) -DMAX_SIZE=$((4096*1024)) ./${TEST}.cpp -o ${OUT}/${TEST}
../dep/sde-external-8.35.0-2019-03-11-lin/sde64 \
	-sse-sde -disasm_att 1 -dcfg 1 -dcfg:write_bb 1 \
	-align_checker_prefetch 0 -align_correct 0 -emu_fast 1 \
	-bdw -- \
	./${OUT}/${TEST}
mv dcfg-out.* ${OUT}
OUT=${TEST}_out_DRAM
mkdir -p ${OUT}
icpc -xHOST -O3 -DMIN_SIZE=$((64*1024*1024)) -DMAX_SIZE=$((128*1024*1024)) ./${TEST}.cpp -o ${OUT}/${TEST}
../dep/sde-external-8.35.0-2019-03-11-lin/sde64 \
	-sse-sde -disasm_att 1 -dcfg 1 -dcfg:write_bb 1 \
	-align_checker_prefetch 0 -align_correct 0 -emu_fast 1 \
	-bdw -- \
	./${OUT}/${TEST}
mv dcfg-out.* ${OUT}

TEST=sdetest3C
OUT=${TEST}_out
mkdir -p ${OUT}
icc ./${TEST}.c -o ${OUT}/${TEST}
../dep/sde-external-8.35.0-2019-03-11-lin/sde64 \
	-sse-sde -disasm_att 1 -dcfg 1 -dcfg:write_bb 1 \
	-align_checker_prefetch 0 -align_correct 0 -emu_fast 1 \
	-bdw -- \
	./${OUT}/${TEST} 64 "A big black dog."
mv dcfg-out.* ${OUT}

TEST=sdetest4C
OUT=${TEST}_out
mkdir -p ${OUT}
#doesnt work with icc (or gcc & > -O1)
gcc -O1 ./${TEST}.c -o ${OUT}/${TEST}
../dep/sde-external-8.35.0-2019-03-11-lin/sde64 \
	-sse-sde -disasm_att 1 -dcfg 1 -dcfg:write_bb 1 \
	-align_checker_prefetch 0 -align_correct 0 -emu_fast 1 \
	-bdw -- \
	./${OUT}/${TEST} 1000000
mv dcfg-out.* ${OUT}

TEST=sdetestF
OUT=${TEST}_out
mkdir -p ${OUT}
ifort ./${TEST}.f -o ${OUT}/${TEST}
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

TEST=sdetest2MPI
OUT=${TEST}_out
mkdir -p ${OUT}
mpiicc -I${ADVISOR_2018_DIR}/include ./${TEST}.c -o ${OUT}/${TEST}
mpiexec ${MPIEXECOPT} -np 4 bash -c "
../dep/sde-external-8.35.0-2019-03-11-lin/sde64 \
	-sse-sde -disasm_att 1 -dcfg 1 -dcfg:write_bb 1 \
	-dcfg:out_base_name dcfg-out.rank-\"\$MPI_LOCALRANKID\" \
	-align_checker_prefetch 0 -align_correct 0 -emu_fast 1 \
	-bdw -- \
	./${OUT}/${TEST}"
mv dcfg-out.* ${OUT}

TEST=sdetest3MPI
OUT=${TEST}_out
mkdir -p ${OUT}
../dep/mpi-wrap-gen/wrap.py -f -g -c mpiicc \
	-o ${OUT}/disable-mpi-in-sde.c \
	./disable-mpi-in-sde.w
mpiicc -O2 -Wall -fPIC ${OUT}/disable-mpi-in-sde.c -shared \
	-o ${OUT}/disable-mpi-in-sde.so
mpiicc -I${ADVISOR_2018_DIR}/include ./${TEST}.c -o ${OUT}/${TEST} #-L${ADVISOR_2018_DIR}/lib64 -littnotify
mpiexec ${MPIEXECOPT} -genv LD_PRELOAD=${OUT}/disable-mpi-in-sde.so -np 4 bash -c "
../dep/sde-external-8.35.0-2019-03-11-lin/sde64 \
	-sse-sde -disasm_att 1 -dcfg 1 -dcfg:write_bb 1 \
	-dcfg:out_base_name dcfg-out.rank-\"\$MPI_LOCALRANKID\" \
	-align_checker_prefetch 0 -align_correct 0 -emu_fast 1 \
	-start_ssc_mark 11dead11 -stop_ssc_mark 22dead22 \
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

TEST=sdetestOMP
OUT=${TEST}_out
mkdir -p ${OUT}
gcc -static -O3 -march=native ./${TEST}.c -o ${OUT}/${TEST} -fopenmp -lm
for OMP in 4 16 64 256 1024; do
	mkdir -p ${OUT}/${OMP}
	OMP_NUM_THREADS=${OMP}\
		../dep/sde-external-8.35.0-2019-03-11-lin/sde64 \
		-sse-sde -disasm_att 1 -dcfg 1 -dcfg:write_bb 1 \
		-align_checker_prefetch 0 -align_correct 0 -emu_fast 1 \
		-bdw -- \
		./${OUT}/${TEST}
	mv dcfg-out.* ${OUT}/${OMP}
done
for OMP in 4 16 64 256 1024; do
	D=sdetestOMP_out/${OMP};
	./parse_basic_blocks.py -j $D/dcfg-out.dcfg.json.bz2 -b $D/dcfg-out.bb.txt.bz2 \
		-s $D/dcfg-out.asm.b > /dev/null
	./parse_basic_blocks.py -j $D/dcfg-out.dcfg.json.bz2 -b $D/dcfg-out.bb.txt.bz2 \
		-l $D/dcfg-out.asm.b > $D/parser.log;
done

TEST=SdeTestJava
OUT=${TEST}_out
mkdir -p ${OUT}
icc -c -fPIC -I${ADVISOR_2019_DIR}/include -I. -I$JAVA_HOME/include -I$JAVA_HOME/include/linux ssc.c
icc -shared -Wl,-soname,libssc.so ssc.o -o libssc.so ${ADVISOR_2019_DIR}/lib64/libittnotify.a
source $ROOTDIR/dep/spack/share/spack/setup-env.sh
spack load openjdk
javac ./SdeTestJava.java
javah -classpath . SdeTestJava
LD_LIBRARY_PATH=`pwd`:$LD_LIBRARY_PATH ../dep/sde-external-8.35.0-2019-03-11-lin/sde64 \
	-sse-sde -disasm_att 1 -dcfg 1 -dcfg:write_bb 1 \
	-align_checker_prefetch 0 -align_correct 0 -emu_fast 1 \
	-bdw -start_ssc_mark 111:repeat -stop_ssc_mark 222:repeat -- \
	`which java` -classpath . SdeTestJava
mv dcfg-out.* ${OUT}/
rm -f ssc.o libssc.so SdeTestJava.h SdeTestJava.class
