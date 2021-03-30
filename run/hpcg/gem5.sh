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

# ============================ HPCG ===========================================
source conf/hpcg.sh
DEFINPUT=$INPUT
LOG="$ROOTDIR/log/`hostname -s`/gem5run/hpcg/conf${1}.log"
mkdir -p `dirname $LOG`
cd $APPDIR
for BEST in $BESTCONF; do
	NumMPI=1
	if [ $1 -eq 1 ]; then NumOMP="20"; elif [ $1 -eq 2 ]; then NumOMP="32"; fi
	echo -e "OMP_NUM_THREADS=$NumOMP\nOMP_NUM_PARALELL=$NumOMP\nFLIB_FASTOMP=FALSE\nFLIB_CNTL_BARRIER_ERR=FALSE" > ./omp${NumOMP}.txt
	#XXX: npx|y|z always 1 for 1==NumMPI: test to identify hpcg's internal dimensions
	X=$MAXXYZ; Y=$MAXXYZ; Z=$MAXXYZ
	#rm -f hpcg_log_* n*.yaml
	#rm -f hpcg20*.txt HPCG-Benchmark_3*.txt
	#mpiexec $MPIEXECOPT -n $NumMPI $BINARY -n 1 > /dev/null 2>&1
	#if [ ! "x$?" = "x0" ]; then continue; fi
	#if [ -f n*.yaml ]; then
	#	X=$(($MAXXYZ / `grep 'npx:' n*.yaml | awk -F 'npx:' '{print $2}'`))
	#	Y=$(($MAXXYZ / `grep 'npy:' n*.yaml | awk -F 'npy:' '{print $2}'`))
	#	Z=$(($MAXXYZ / `grep 'npz:' n*.yaml | awk -F 'npz:' '{print $2}'`))
	#elif [ -f HPCG-Benchmark_3*.txt ]; then
	#	#non-intel version needs to be div8 https://github.com/hpcg-benchmark/hpcg/issues/47
	#	X=$((($MAXXYZ / `grep 'npx=' HPCG-Benchmark_3*.txt | awk -F 'npx=' '{print $2}'` / 8) * 8))
	#	Y=$((($MAXXYZ / `grep 'npy=' HPCG-Benchmark_3*.txt | awk -F 'npy=' '{print $2}'` / 8) * 8))
	#	Z=$((($MAXXYZ / `grep 'npz=' HPCG-Benchmark_3*.txt | awk -F 'npz=' '{print $2}'` / 8) * 8))
	#else continue; fi
	rm -f hpcg_log_* n*.yaml
	rm -f hpcg20*.txt HPCG-Benchmark_3*.txt
	INPUT="`echo $DEFINPUT | sed -e \"s/NX/$X/\" -e \"s/NY/$Y/\" -e \"s/NZ/$Z/\"`"
	echo "=== gem5 run ===" >> $LOG 2>&1
	echo "$PIN $GEM5 -d `dirname $LOG` $GEM5SE -c $BINARY -o \"$INPUT\" -n $NumOMP -e ./omp${NumOMP}.txt $ARCHCONF" >> $LOG 2>&1
	$PIN $GEM5 -d `dirname $LOG` $GEM5SE -c $BINARY -o "$INPUT" -n $NumOMP -e ./omp${NumOMP}.txt $ARCHCONF >> $LOG 2>&1
	cat hpcg_log_* >> $LOG 2>&1
	cat n*.yaml >> $LOG 2>&1
	rm -f hpcg_log_* n*.yaml
	cat hpcg20*.txt >> $LOG 2>&1
	cat HPCG-Benchmark_3*.txt >> $LOG 2>&1
	rm -f hpcg20*.txt HPCG-Benchmark_3*.txt
done
cd $ROOTDIR
