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

# special core pinning?
if [ -n $2 ] && which numactl >/dev/null 2>&1; then PIN="numactl -C $2"; else PIN=""; fi

# ============================ BabelStream ====================================
source conf/babelstream.sh
DEFLOG="$ROOTDIR/log/`hostname -s`/gem5run/babelstream"
mkdir -p $DEFLOG
cd $APPDIR
DEFINPUT=$INPUT
for BEST in 1; do
	for BINARY in $BINARYS; do
		NumMPI=1
		if [ $1 -eq 1 ]; then NumOMP="20"; elif [ $1 -eq 2 ]; then NumOMP="32"; fi
		echo -e "OMP_NUM_THREADS=$NumOMP\nOMP_NUM_PARALELL=$NumOMP\nFLIB_FASTOMP=FALSE\nFLIB_CNTL_BARRIER_ERR=FALSE" > ./omp${NumOMP}.txt
		S="`echo $BINARY | cut -d '_' -f2`"
		BINARY="`echo $BINARY | cut -d '_' -f1`"
		LOG="${DEFLOG}/conf${1}.${S}gb.log"
		S=$((S*1024*1024*1024/8))
		INPUT="`echo $DEFINPUT | sed -e \"s/SIZE/$S/\"`"
		echo "$PIN $GEM5 -d `dirname $LOG` $GEM5SE -c $BINARY -o \"$INPUT\" -n $NumOMP -e ./omp${NumOMP}.txt $ARCHCONF" >> $LOG 2>&1
		$PIN $GEM5 -d `dirname $LOG` $GEM5SE -c $BINARY -o "$INPUT" -n $NumOMP -e ./omp${NumOMP}.txt $ARCHCONF >> $LOG 2>&1
	done
done
cd $ROOTDIR
