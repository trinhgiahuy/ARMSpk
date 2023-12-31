#!/bin/bash

ROOTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/../../"
cd $ROOTDIR

source $ROOTDIR/conf/host.cfg
ulimit -s unlimited
ulimit -n 4096

export HOSTNAME="kiev0" # we don't care where gem runs
GEM5="$ROOTDIR/dep/gem5_riken/build/ARM/gem5.opt"
GEM5SE="$ROOTDIR/dep/gem5_riken/configs/example/se.py"

if [ $1 -eq 1 ]; then
	ARCHCONF="--cpu-type=O3_ARM_PostK_3 --caches --l2_size=16MB --mem_bus_width=64 --mem_resp_width=128 --mem-size=32GB"    # traditional
elif [ $1 -eq 2 ]; then
	ARCHCONF="--cpu-type=O3_ARM_PostK_3 --caches --l2_size=256MB --mem_bus_width=96 --mem_resp_width=192 --mem-size=32GB"   # aggressive
else echo "err: missing which arch conf to run"; exit 1; fi

# special core pinning? XXX: no, we run all simul on a node
#if [ -n $2 ] && which numactl >/dev/null 2>&1; then PIN="numactl -C $2"; else PIN=""; fi

# ============================ PolyBench ====================================
source conf/polybench.sh
DEFLOG="$ROOTDIR/log/`hostname -s`/gem5run/polybench"
mkdir -p $DEFLOG
cd $APPDIR
for BEST in $BESTCONF; do
	for BMconf in $BINARYS; do
		NumMPI=1
		NumOMP=1	#XXX: no OMP in PolyBench: if [ $1 -eq 1 ]; then NumOMP="20"; elif [ $1 -eq 2 ]; then NumOMP="32"; fi
		BINARY="`echo ${BMconf} | tr '|' '_'`"
		BName="`basename $(echo ${BMconf} | cut -d'|' -f1)`"
		echo -e "OMP_NUM_THREADS=$NumOMP\nOMP_NUM_PARALELL=$NumOMP\nFLIB_FASTOMP=FALSE\nFLIB_CNTL_BARRIER_ERR=FALSE" > ./${BName}-omp${NumOMP}.txt
		LOG="${DEFLOG}/${BName}/conf${1}.log"; mkdir -p ${LOG}_stat
		echo "$GEM5 -d ${LOG}_stat $GEM5SE -c $BINARY -o \"$INPUT\" -n $NumOMP -e ./${BName}-omp${NumOMP}.txt $ARCHCONF" >> $LOG 2>&1
		$GEM5 -d ${LOG}_stat $GEM5SE -c $BINARY -o "$INPUT" -n $NumOMP -e ./${BName}-omp${NumOMP}.txt $ARCHCONF >> $LOG 2>&1 &
	done
	wait
done
cd $ROOTDIR
